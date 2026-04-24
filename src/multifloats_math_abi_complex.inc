// Complex DD C-ABI shims. Every entry point below marshals
// complex64x2_t ↔ std::complex<multifloats::float64x2> and dispatches to
// the corresponding C++ implementation:
//   • transcendentals + abs/arg/proj → std:: specializations in
//     multifloats_math_complex_std.inc (pulled into multifloats_math.cc
//     inside `namespace std`)
//   • Kind-D overloads (log2/log10/log1p/expm1/sinpi/cospi) →
//     `namespace multifloats` overloads in
//     multifloats_math_abi_complex_kernels.inc
//   • trivial accessors (conj/real/imag) → std::complex member or
//     generic template
//
// The compensated complex DD multiply and divide keep their full bodies
// here: std::complex's `operator*` / `operator/` templates are not
// specializable, so the cancellation-gated paths are only reachable via
// this C-ABI surface. Anon-namespace helpers `dd_cross_diff` and
// `dd_x2y2m1` at file scope in multifloats_math.cc are used by those
// kernels and by the std:: log/asinh/atanh bodies.

// ---- Complex DD arithmetic ------------------------------------------------
// +, − are componentwise; std::complex's built-in operators already do the
// right thing, but the wrappers marshal explicitly to avoid accidentally
// constructing a std::complex from complex64x2_t.

complex64x2_t cadddd(complex64x2_t a, complex64x2_t b) {
  return cc_to(cc_from(a) + cc_from(b));
}

complex64x2_t csubdd(complex64x2_t a, complex64x2_t b) {
  return cc_to(cc_from(a) - cc_from(b));
}

// 4-mul form with cancellation-gated compensation. Each component is
// computed the cheap way first (two DD muls + a DD add/sub); if the
// leading-limb ratio |R|/max(|p|,|q|) drops below 2⁻⁴ the DD sub/add
// shed cancellation bits below DD precision, so re-run that component
// through `dd_cross_diff` to recover the full DD answer. The detector
// also fires on the classic exact-zero witness a = (1, ε), b = (-1, ε)
// pinned by `test_complex_mul_cancellation` in test/test.cc: Re sees
// no cancellation so the cheap path keeps the full-DD -1 - ε²; Im
// hi-limbs cancel exactly, the compensated path returns exact 0.
//
// Im = ar·bi + ai·br is a *sum*; route it through dd_cross_diff by
// negating one operand (DD-limb negation is exact, no precision cost).
//
// Karatsuba (3 muls, 5 adds, Im = (a+b)(c+d) − ac − bd) was tried and
// rejected in Fortran: ~2× slower and less precise — the subtraction
// loses bits when |ac+bd| ≫ |Im|.
complex64x2_t cmuldd(complex64x2_t a, complex64x2_t b) {
  float64x2 ar = from(a.re), ai = from(a.im);
  float64x2 br = from(b.re), bi = from(b.im);

  float64x2 p_re = ar * br;
  float64x2 q_re = ai * bi;
  float64x2 re_dd = p_re - q_re;

  float64x2 p_im = ar * bi;
  float64x2 q_im = ai * br;
  float64x2 im_dd = p_im + q_im;

  // 2⁻⁴ caps worst-case rel-err of the cheap path at ~16 ulp_dd while
  // keeping the compensated branch off the hot path on random data
  // (|R|/max < 1/16 on ~6% of uniform-quadrant inputs). Tighter 2⁻¹
  // triggered on ~25% and more than halved cdd_mul bench speedup;
  // looser 2⁻⁸ let the 100+-ulp regime through. See the `cmuldd`
  // landmark in doc/developer/INTERNALS.md §6.
  constexpr double kCancelThresh = 0x1p-4;

  double abs_re   = std::fabs(re_dd._limbs[0]);
  double abs_p_re = std::fabs(p_re._limbs[0]);
  double abs_q_re = std::fabs(q_re._limbs[0]);
  double max_re   = abs_p_re > abs_q_re ? abs_p_re : abs_q_re;

  double abs_im   = std::fabs(im_dd._limbs[0]);
  double abs_p_im = std::fabs(p_im._limbs[0]);
  double abs_q_im = std::fabs(q_im._limbs[0]);
  double max_im   = abs_p_im > abs_q_im ? abs_p_im : abs_q_im;

  float64x2 re = re_dd, im = im_dd;
  if (__builtin_expect(max_re != 0.0 && abs_re < kCancelThresh * max_re, 0)) {
    re = dd_cross_diff(ar, br, ai, bi);
  }
  if (__builtin_expect(max_im != 0.0 && abs_im < kCancelThresh * max_im, 0)) {
    // ar·bi + ai·br = ar·bi − (−ai)·br. Negating both limbs is exact.
    float64x2 neg_ai(-ai._limbs[0], -ai._limbs[1]);
    im = dd_cross_diff(ar, bi, neg_ai, br);
  }
  return { to(re), to(im) };
}

// Complex DD divide: (a + bi) / (c + di) with compensated numerators.
//   re = (a·c + b·d) / (c² + d²)         — no cancellation (sum)
//   im = (b·c − a·d) / (c² + d²)         — cancels when bc ≈ ad
// The imaginary numerator uses dd_cross_diff to recover the full DD
// precision lost by the naïve (DD_mul)−(DD_mul) subtraction.
//
// Overflow guard: the unscaled c²+d² overflows once max(|c|,|d|) >
// √DBL_MAX ≈ 1.34e154, well inside DD's representable range. Pre-scale
// the denominator by an exact power of 2 (chosen so |max| ∈ [0.5, 1)),
// compute the formula, then unscale the result. Power-of-2 scaling is
// exact at every limb (no precision lost), and the post-scale denom is
// bounded by 2 so no overflow can occur.
complex64x2_t cdivdd(complex64x2_t a, complex64x2_t b) {
  float64x2 ar = from(a.re), ai = from(a.im);
  float64x2 br = from(b.re), bi = from(b.im);

  // Fast path: b is in the range where br²+bi² won't overflow/underflow.
  // The square of √DBL_MAX (~1.34e154) is the boundary; allow a generous
  // margin and skip the ldexp/frexp pre/post-scaling for typical inputs.
  double abs_br = std::fabs(br._limbs[0]);
  double abs_bi = std::fabs(bi._limbs[0]);
  double max_b = abs_br > abs_bi ? abs_br : abs_bi;
  if (max_b < 1e150 && max_b > 1e-150) {
    float64x2 denom = br * br + bi * bi;
    float64x2 num_re = ar * br + ai * bi;
    float64x2 num_im = dd_cross_diff(ai, br, ar, bi);
    return { to(num_re / denom), to(num_im / denom) };
  }
  // Slow path: scale b by an exact power of 2 so |max element| ∈ [0.5, 1)
  // and br²+bi² is bounded by 2. Power-of-2 scaling preserves all DD limbs
  // exactly; we unscale the result by the same power at the end.
  int e = 0;
  if (max_b > 0.0 && std::isfinite(max_b)) {
    std::frexp(max_b, &e);
    br = float64x2(std::ldexp(br._limbs[0], -e),
                   std::ldexp(br._limbs[1], -e));
    bi = float64x2(std::ldexp(bi._limbs[0], -e),
                   std::ldexp(bi._limbs[1], -e));
  }
  float64x2 denom = br * br + bi * bi;
  float64x2 num_re = ar * br + ai * bi;
  float64x2 num_im = dd_cross_diff(ai, br, ar, bi);
  float64x2 res_re = num_re / denom;
  float64x2 res_im = num_im / denom;
  if (e != 0) {
    res_re = float64x2(std::ldexp(res_re._limbs[0], -e),
                       std::ldexp(res_re._limbs[1], -e));
    res_im = float64x2(std::ldexp(res_im._limbs[0], -e),
                       std::ldexp(res_im._limbs[1], -e));
  }
  return { to(res_re), to(res_im) };
}

// ---- Transcendentals → std:: specializations ------------------------------

complex64x2_t cexpdd  (complex64x2_t z) { return cc_to(std::exp  (cc_from(z))); }
complex64x2_t clogdd  (complex64x2_t z) { return cc_to(std::log  (cc_from(z))); }
complex64x2_t csqrtdd (complex64x2_t z) { return cc_to(std::sqrt (cc_from(z))); }
complex64x2_t cpowdd  (complex64x2_t z, complex64x2_t w) {
  return cc_to(std::pow(cc_from(z), cc_from(w)));
}

complex64x2_t csindd  (complex64x2_t z) { return cc_to(std::sin  (cc_from(z))); }
complex64x2_t ccosdd  (complex64x2_t z) { return cc_to(std::cos  (cc_from(z))); }
complex64x2_t ctandd  (complex64x2_t z) { return cc_to(std::tan  (cc_from(z))); }
complex64x2_t casindd (complex64x2_t z) { return cc_to(std::asin (cc_from(z))); }
complex64x2_t cacosdd (complex64x2_t z) { return cc_to(std::acos (cc_from(z))); }
complex64x2_t catandd (complex64x2_t z) { return cc_to(std::atan (cc_from(z))); }

complex64x2_t csinhdd (complex64x2_t z) { return cc_to(std::sinh (cc_from(z))); }
complex64x2_t ccoshdd (complex64x2_t z) { return cc_to(std::cosh (cc_from(z))); }
complex64x2_t ctanhdd (complex64x2_t z) { return cc_to(std::tanh (cc_from(z))); }
complex64x2_t casinhdd(complex64x2_t z) { return cc_to(std::asinh(cc_from(z))); }
complex64x2_t cacoshdd(complex64x2_t z) { return cc_to(std::acosh(cc_from(z))); }
complex64x2_t catanhdd(complex64x2_t z) { return cc_to(std::atanh(cc_from(z))); }

// ---- Kind-D overloads in namespace multifloats ---------------------------

complex64x2_t cexpm1dd (complex64x2_t z) { return cc_to(multifloats::expm1(cc_from(z))); }
complex64x2_t clog2dd  (complex64x2_t z) { return cc_to(multifloats::log2 (cc_from(z))); }
complex64x2_t clog10dd (complex64x2_t z) { return cc_to(multifloats::log10(cc_from(z))); }
complex64x2_t clog1pdd (complex64x2_t z) { return cc_to(multifloats::log1p(cc_from(z))); }
complex64x2_t csinpidd (complex64x2_t z) { return cc_to(multifloats::sinpi(cc_from(z))); }
complex64x2_t ccospidd (complex64x2_t z) { return cc_to(multifloats::cospi(cc_from(z))); }

// ---- Value accessors -----------------------------------------------------

float64x2_t   cabsdd  (complex64x2_t z) { return to(std::abs (cc_from(z))); }
float64x2_t   cargdd  (complex64x2_t z) { return to(std::arg (cc_from(z))); }
complex64x2_t cprojdd (complex64x2_t z) { return cc_to(std::proj(cc_from(z))); }

// Conjugate: std::conj on the generic template is already componentwise
// exact DD negation — each limb flips sign, no roundoff.
complex64x2_t conjdd  (complex64x2_t z) { return cc_to(std::conj(cc_from(z))); }

// real/imag: trivial accessors. Returned directly to skip the round-trip
// through std::complex construction.
float64x2_t   crealdd (complex64x2_t z) { return z.re; }
float64x2_t   cimagdd (complex64x2_t z) { return z.im; }

// ---- Matrix multiply (column-major) --------------------------------------
//
// Inner accumulator uses Neumaier-style compensation on the DD-product stream:
//   p, e = two_prod(ah, bh)          (exact)
//   cross = fma(ah, bl, al*bh)       (DD-mul cross term, ignoring bl*bl)
//   (s_hi, s_lo) += (p + cross + e)  via two_sum on s_hi
// A final two_sum renormalizes s_hi/s_lo into the DD result; if
// `renorm_interval > 0`, intermediate two_sums also fire every
// that-many reductions to keep s_lo bounded at large k.
//
// Loop order is AXPY / gaxpy: outer p (shared dim), inner i (output row)
// so A is read as contiguous columns and the output accumulators are
// m-way data-parallel. The row dimension is register-blocked to MR=8
// via a compile-time-sized panel template; any m is handled by
// MR panels + a 1..7 tail. mm calls the same mv dispatch per output
// column.
