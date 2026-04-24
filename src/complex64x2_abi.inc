// Complex DD C-ABI shims. `complex64x2` is an alias for
// `std::complex<multifloats::float64x2>` in C++, so each entry point calls
// the corresponding C++ implementation directly:
//   • transcendentals + abs/arg/proj → std:: specializations in
//     complex64x2_std.inc (pulled into multifloats_math.cc
//     inside `namespace std`)
//   • Kind-D overloads (log2/log10/log1p/expm1/sinpi/cospi) →
//     `namespace multifloats` overloads in
//     complex64x2_abi_kernels.inc
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
// right thing on std::complex<float64x2>.

complex64x2 cadddd(complex64x2 a, complex64x2 b) { return a + b; }
complex64x2 csubdd(complex64x2 a, complex64x2 b) { return a - b; }

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
complex64x2 cmuldd(complex64x2 a, complex64x2 b) {
  float64x2 ar = a.real(), ai = a.imag();
  float64x2 br = b.real(), bi = b.imag();

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

  double abs_re   = std::fabs(re_dd.limbs[0]);
  double abs_p_re = std::fabs(p_re.limbs[0]);
  double abs_q_re = std::fabs(q_re.limbs[0]);
  double max_re   = abs_p_re > abs_q_re ? abs_p_re : abs_q_re;

  double abs_im   = std::fabs(im_dd.limbs[0]);
  double abs_p_im = std::fabs(p_im.limbs[0]);
  double abs_q_im = std::fabs(q_im.limbs[0]);
  double max_im   = abs_p_im > abs_q_im ? abs_p_im : abs_q_im;

  float64x2 re = re_dd, im = im_dd;
  if (__builtin_expect(max_re != 0.0 && abs_re < kCancelThresh * max_re, 0)) {
    re = dd_cross_diff(ar, br, ai, bi);
  }
  if (__builtin_expect(max_im != 0.0 && abs_im < kCancelThresh * max_im, 0)) {
    // ar·bi + ai·br = ar·bi − (−ai)·br. DD negation is exact (componentwise).
    im = dd_cross_diff(ar, bi, -ai, br);
  }
  return { re, im };
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
complex64x2 cdivdd(complex64x2 a, complex64x2 b) {
  float64x2 ar = a.real(), ai = a.imag();
  float64x2 br = b.real(), bi = b.imag();

  // Fast path: b is in the range where br²+bi² won't overflow/underflow.
  // The square of √DBL_MAX (~1.34e154) is the boundary; allow a generous
  // margin and skip the ldexp/frexp pre/post-scaling for typical inputs.
  double abs_br = std::fabs(br.limbs[0]);
  double abs_bi = std::fabs(bi.limbs[0]);
  double max_b = abs_br > abs_bi ? abs_br : abs_bi;
  if (max_b < 1e150 && max_b > 1e-150) {
    float64x2 denom = br * br + bi * bi;
    float64x2 num_re = ar * br + ai * bi;
    float64x2 num_im = dd_cross_diff(ai, br, ar, bi);
    return { num_re / denom, num_im / denom };
  }
  // Slow path: scale b by an exact power of 2 so |max element| ∈ [0.5, 1)
  // and br²+bi² is bounded by 2. Power-of-2 scaling preserves all DD limbs
  // exactly; we unscale the result by the same power at the end.
  int e = 0;
  if (max_b > 0.0 && std::isfinite(max_b)) {
    std::frexp(max_b, &e);
    br = float64x2(std::ldexp(br.limbs[0], -e),
                   std::ldexp(br.limbs[1], -e));
    bi = float64x2(std::ldexp(bi.limbs[0], -e),
                   std::ldexp(bi.limbs[1], -e));
  }
  float64x2 denom = br * br + bi * bi;
  float64x2 num_re = ar * br + ai * bi;
  float64x2 num_im = dd_cross_diff(ai, br, ar, bi);
  float64x2 res_re = num_re / denom;
  float64x2 res_im = num_im / denom;
  if (e != 0) {
    res_re = float64x2(std::ldexp(res_re.limbs[0], -e),
                       std::ldexp(res_re.limbs[1], -e));
    res_im = float64x2(std::ldexp(res_im.limbs[0], -e),
                       std::ldexp(res_im.limbs[1], -e));
  }
  return { res_re, res_im };
}

// ---- Transcendentals → std:: specializations ------------------------------

complex64x2 cexpdd  (complex64x2 z) { return std::exp (z); }
complex64x2 clogdd  (complex64x2 z) { return std::log (z); }
complex64x2 csqrtdd (complex64x2 z) { return std::sqrt(z); }
complex64x2 cpowdd  (complex64x2 z, complex64x2 w) { return std::pow(z, w); }

complex64x2 csindd  (complex64x2 z) { return std::sin (z); }
complex64x2 ccosdd  (complex64x2 z) { return std::cos (z); }
complex64x2 ctandd  (complex64x2 z) { return std::tan (z); }
complex64x2 casindd (complex64x2 z) { return std::asin(z); }
complex64x2 cacosdd (complex64x2 z) { return std::acos(z); }
complex64x2 catandd (complex64x2 z) { return std::atan(z); }

complex64x2 csinhdd (complex64x2 z) { return std::sinh (z); }
complex64x2 ccoshdd (complex64x2 z) { return std::cosh (z); }
complex64x2 ctanhdd (complex64x2 z) { return std::tanh (z); }
complex64x2 casinhdd(complex64x2 z) { return std::asinh(z); }
complex64x2 cacoshdd(complex64x2 z) { return std::acosh(z); }
complex64x2 catanhdd(complex64x2 z) { return std::atanh(z); }

// ---- Kind-D overloads in namespace multifloats ---------------------------

complex64x2 cexpm1dd (complex64x2 z) { return multifloats::expm1(z); }
complex64x2 clog2dd  (complex64x2 z) { return multifloats::log2 (z); }
complex64x2 clog10dd (complex64x2 z) { return multifloats::log10(z); }
complex64x2 clog1pdd (complex64x2 z) { return multifloats::log1p(z); }
complex64x2 csinpidd (complex64x2 z) { return multifloats::sinpi(z); }
complex64x2 ccospidd (complex64x2 z) { return multifloats::cospi(z); }

// ---- Value accessors -----------------------------------------------------

float64x2   cabsdd  (complex64x2 z) { return std::abs(z); }
float64x2   cargdd  (complex64x2 z) { return std::arg(z); }
complex64x2 cprojdd (complex64x2 z) { return std::proj(z); }

// Conjugate: std::conj on the generic template is already componentwise
// exact DD negation — each limb flips sign, no roundoff.
complex64x2 conjdd  (complex64x2 z) { return std::conj(z); }

// real/imag: trivial accessors. Returned directly to skip the round-trip
// through std::complex construction.
float64x2   crealdd (complex64x2 z) { return z.real(); }
float64x2   cimagdd (complex64x2 z) { return z.imag(); }

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
