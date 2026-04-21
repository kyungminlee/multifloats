// ---- Complex DD arithmetic & transcendentals -----------------------------
//
// Branch cuts follow C99 Annex G (matching libquadmath cexpq/clogq/csqrtq).
// Where the real-axis formula needs both sin(y) and cos(y) (or both sinh
// and cosh), these call the fused sincos_full / sinhcosh_full kernels
// so one range reduction + Taylor pair covers both outputs. Functions that
// are built on top of cdd_sqrt / cdd_log (asin/acos/...) call the exported
// cdd_* symbols directly; within a single TU the compiler inlines those.
//
// The text "2 fused" vs "4 separate" in speed comments below is measured
// against libstdc++'s generic <complex> template path (see commit notes).

// ---- Complex DD arithmetic ------------------------------------------------
// Basic (+, -, *, /) on complex DD. These used to be inlined in the Fortran
// `cdd_arith_body` macro; moved here so the Fortran elemental interface is
// a thin bind(c) wrapper and C++ owns the algorithm. cdivdd's imaginary
// numerator uses a compensated cross-difference that sums the full 14-term
// exact expansion of a*b − c*d through a triple-double accumulator, so the
// result stays DD-accurate even under heavy |b·c| ≈ |a·d| cancellation.

complex64x2_t cadddd(complex64x2_t a, complex64x2_t b) {
  float64x2 ar = from(a.re), ai = from(a.im);
  float64x2 br = from(b.re), bi = from(b.im);
  return { to(ar + br), to(ai + bi) };
}

complex64x2_t csubdd(complex64x2_t a, complex64x2_t b) {
  float64x2 ar = from(a.re), ai = from(a.im);
  float64x2 br = from(b.re), bi = from(b.im);
  return { to(ar - br), to(ai - bi) };
}

// a*b − c*d in DD with compensated cross-difference: captures the
// sub-ulp residuals dropped by a normal DD mul so the subtraction
// stays full-DD even under catastrophic cancellation (|a*b| ≈ |c*d|).
// Each DD mul would discard ~1 DD-ulp of its product; when the two
// products agree to ~N bits, that residual becomes the 2^(-106+N)
// relative error of the result. Here we keep the residuals and fold
// them into a two_sum chain.
__attribute__((noinline, cold))
static float64x2 dd_cross_diff(float64x2 a, float64x2 b,
                               float64x2 c, float64x2 d) {
  double a0 = a._limbs[0], a1 = a._limbs[1];
  double b0 = b._limbs[0], b1 = b._limbs[1];
  double c0 = c._limbs[0], c1 = c._limbs[1];
  double d0 = d._limbs[0], d1 = d._limbs[1];
  // Full exact expansion of (a*b − c*d) with DD inputs into 14 scalar
  // doubles, accumulated in a triple-double (TD) so the leading 53-bit
  // bucket survives even when a*b and c*d don't cancel (and the trailing
  // buckets survive when they do). Each scalar two_prod via FMA is
  // exact; only the two a1*b1 / c1*d1 tails are rounded once (below
  // DD-ulp, so invisible in the DD output).
  double ab00_h = a0 * b0, ab00_l = std::fma(a0, b0, -ab00_h);
  double cd00_h = c0 * d0, cd00_l = std::fma(c0, d0, -cd00_h);
  double ab01_h = a0 * b1, ab01_l = std::fma(a0, b1, -ab01_h);
  double ab10_h = a1 * b0, ab10_l = std::fma(a1, b0, -ab10_h);
  double cd01_h = c0 * d1, cd01_l = std::fma(c0, d1, -cd01_h);
  double cd10_h = c1 * d0, cd10_l = std::fma(c1, d0, -cd10_h);
  double ab11 = a1 * b1;
  double cd11 = c1 * d1;
  double T0 = 0, T1 = 0, T2 = 0;
  #define TSUM(x) do { \
    double _x = (x), _s, _bb, _e; \
    _s = T0 + _x; _bb = _s - T0; _e = (T0 - (_s - _bb)) + (_x - _bb); \
    T0 = _s; _x = _e; \
    _s = T1 + _x; _bb = _s - T1; _e = (T1 - (_s - _bb)) + (_x - _bb); \
    T1 = _s; T2 += _e; \
  } while (0)
  // Largest-magnitude terms first so cancellation residuals land in T1/T2.
  TSUM(ab00_h);  TSUM(-cd00_h);
  TSUM(ab00_l);  TSUM(-cd00_l);
  TSUM(ab01_h);  TSUM(ab10_h);  TSUM(-cd01_h); TSUM(-cd10_h);
  TSUM(ab01_l);  TSUM(ab10_l);  TSUM(-cd01_l); TSUM(-cd10_l);
  TSUM(ab11);    TSUM(-cd11);
  #undef TSUM
  // Renormalize TD → DD: two_sum T1,T2 → T1; two_sum T0,T1 → (T0, T1);
  // fold T2 into T1; final fast_two_sum.
  { double s = T1 + T2, bb = s - T1; T2 = (T1 - (s - bb)) + (T2 - bb); T1 = s; }
  { double s = T0 + T1, bb = s - T0; T1 = (T0 - (s - bb)) + (T1 - bb); T0 = s; }
  T1 += T2;
  double hi = T0 + T1;
  double lo = T1 - (hi - T0);
  return float64x2(hi, lo);
}

// x*x + y*y − 1 in DD without cancellation. Each DD square expands to
// five scalar-double terms (two_prod + 2·x0·x1 + x1*x1); the −1 and the
// y terms feed the same TD accumulator so 0.5·log1p(x²+y²−1) stays
// full-DD even when |z|² ≈ 1. Mirrors libquadmath's __quadmath_x2y2m1q
// but on DD inputs.
static inline float64x2 dd_x2y2m1(float64x2 x, float64x2 y) {
  double x0 = x._limbs[0], x1 = x._limbs[1];
  double y0 = y._limbs[0], y1 = y._limbs[1];
  double xx00_h = x0 * x0, xx00_l = std::fma(x0, x0, -xx00_h);
  double yy00_h = y0 * y0, yy00_l = std::fma(y0, y0, -yy00_h);
  double twox0 = x0 + x0, twoy0 = y0 + y0;          // exact (×2)
  double xx01_h = twox0 * x1, xx01_l = std::fma(twox0, x1, -xx01_h);
  double yy01_h = twoy0 * y1, yy01_l = std::fma(twoy0, y1, -yy01_h);
  double xx11 = x1 * x1;
  double yy11 = y1 * y1;
  double T0 = 0, T1 = 0, T2 = 0;
  #define TSUM(v) do { \
    double _x = (v), _s, _bb, _e; \
    _s = T0 + _x; _bb = _s - T0; _e = (T0 - (_s - _bb)) + (_x - _bb); \
    T0 = _s; _x = _e; \
    _s = T1 + _x; _bb = _s - T1; _e = (T1 - (_s - _bb)) + (_x - _bb); \
    T1 = _s; T2 += _e; \
  } while (0)
  // Largest-magnitude terms first; the −1 goes in with the leading
  // two_prod hi-limbs so its cancellation residual lands in T1/T2.
  TSUM(xx00_h);  TSUM(yy00_h);  TSUM(-1.0);
  TSUM(xx01_h);  TSUM(yy01_h);
  TSUM(xx00_l);  TSUM(yy00_l);
  TSUM(xx01_l);  TSUM(yy01_l);
  TSUM(xx11);    TSUM(yy11);
  #undef TSUM
  { double s = T1 + T2, bb = s - T1; T2 = (T1 - (s - bb)) + (T2 - bb); T1 = s; }
  { double s = T0 + T1, bb = s - T0; T1 = (T0 - (s - bb)) + (T1 - bb); T0 = s; }
  T1 += T2;
  double hi = T0 + T1;
  double lo = T1 - (hi - T0);
  return float64x2(hi, lo);
}

// 4-mul form with cancellation-gated compensation. Each component is
// computed the cheap way first (two DD muls + a DD add/sub); if the
// leading-limb ratio |R|/max(|p|,|q|) drops below 2⁻¹ the DD sub/add
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
  // __builtin_expect marks the compensated branch cold so the compiler
  // keeps the hot path (no compensation) straight-line and doesn't
  // speculatively schedule dd_cross_diff ops into it.
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
complex64x2_t cdivdd(complex64x2_t a, complex64x2_t b) {
  float64x2 ar = from(a.re), ai = from(a.im);
  float64x2 br = from(b.re), bi = from(b.im);
  float64x2 denom = br * br + bi * bi;
  float64x2 num_re = ar * br + ai * bi;
  float64x2 num_im = dd_cross_diff(ai, br, ar, bi);
  return { to(num_re / denom), to(num_im / denom) };
}

// exp(a+bi) = e^a · (cos b + i·sin b). 2 transcendentals (exp + sincos).
complex64x2_t cexpdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 ea = exp_full(a);
  float64x2 s, c;
  sincos_full(b, s, c);
  return { to(ea * c), to(ea * s) };
}

// log(a+bi) = log(|z|) + i·atan2(b, a). Two paths for log(|z|):
//   • |z| near 1 (big ∈ [0.5, 2]):  log|z| = 0.5·log1p(a²+b²−1) via
//     dd_x2y2m1 — the general formula loses bits when log(big) and
//     0.5·log(1+ratio²) cancel (both ≈ 0 ± ε).
//   • otherwise: log(big) + 0.5·log(1 + (small/big)²), which avoids
//     hypot overflow when |z| exceeds DBL_MAX yet log|z| is finite.
// atan2 handles the negative-real-axis cut and signed-zero propagation.
complex64x2_t clogdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 ax = multifloats::abs(a);
  float64x2 ay = multifloats::abs(b);
  float64x2 big   = (ax > ay) ? ax : ay;
  float64x2 small = (ax > ay) ? ay : ax;
  float64x2 r;
  if (big._limbs[0] == 0.0) {
    r._limbs[0] = -std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
  } else if (big._limbs[0] >= 0.5 && big._limbs[0] <= 2.0) {
    float64x2 u = dd_x2y2m1(big, small);
    r = float64x2(0.5, 0.0) * log1p_full(u);
  } else {
    float64x2 ratio = small / big;
    float64x2 one = float64x2(1.0, 0.0);
    float64x2 half = float64x2(0.5, 0.0);
    r = log_full(big) + half * log_full(one + ratio * ratio);
  }
  float64x2 phi = atan2_full(b, a);
  return { to(r), to(phi) };
}

// log10(z) = log(z) · (1/ln 10).
complex64x2_t clog10dd(complex64x2_t z) {
  static const float64x2 inv_ln10 = float64x2(0x1.bcb7b1526e50ep-2,
                                      0x1.95355baaafad3p-57);
  complex64x2_t l = clogdd(z);
  return { to(from(l.re) * inv_ln10), to(from(l.im) * inv_ln10) };
}

// log2(z) = log(z) · (1/ln 2).
complex64x2_t clog2dd(complex64x2_t z) {
  static const float64x2 inv_ln2 =
      float64x2(1.4426950408889634, 2.0355273740931033e-17);
  complex64x2_t l = clogdd(z);
  return { to(from(l.re) * inv_ln2), to(from(l.im) * inv_ln2) };
}

// log1p(z) = log(1 + z) without losing precision for |z| << 1. Uses
// log1p on the real scalar form `|1+z|² - 1 = 2a + a² + b²`, which is
// computed without the 1 ± 1 cancellation. The imaginary part is
// `atan2(b, 1 + a)` — atan2 handles the branch cut and signed-zero
// propagation directly. Matches libquadmath clog1pq.
complex64x2_t clog1pdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 u = (a + a) + a * a + b * b;          // 2a + a² + b²
  float64x2 half(0.5, 0.0);
  float64x2 r = half * log1p_full(u);
  float64x2 phi = atan2_full(b, float64x2(1.0, 0.0) + a);
  return { to(r), to(phi) };
}

// expm1(z) = e^z - 1 without losing precision for |z| << 1. Write
//   Re = (e^a - 1) cos(b) + (cos(b) - 1)
//      = expm1(a)·cos(b) − 2·sin²(b/2)
//   Im = e^a · sin(b)
// The half-angle form of (cos(b) − 1) keeps the real part full-DD for
// small |b| where the direct `cos(b) − 1` cancels. The Im part calls
// exp_full(a) directly rather than reconstructing (1 + expm1(a)): for
// large negative a, expm1(a) = exp(a) − 1 is pinned at DD-ulp ~2⁻¹⁰⁵
// absolute, so `1 + expm1(a)` delivers only double precision relative
// to |exp(a)| ≪ 1. Cost: one extra scalar exp.
//
// Re worst case on the DD-only half-angle form grows ~linearly with
// fuzz_mpfr sample count (184 ulp at N=10^5, 1078 ulp at N=10^7). The
// shipped TD Re path below lifts that ceiling; see
// doc/developer/TRIPLE_DOUBLE.md for why the arithmetic alone is not
// enough and the constants (log2e, exp2 coefs/table, sin/cos Taylor,
// π/2 split) must be TD too.
// expm1(z) = e^z - 1. The half-angle DD form
//   Re = expm1(a)·cos(b) − 2·sin²(b/2),  Im = e^a·sin(b)
// gives full DD precision *except* on the cancellation surface
// cos(b)·e^a = 1, where the two Re terms agree to many bits and the
// DD abs-err floor ~2⁻¹⁰⁵·max_term bleeds through (1078 DD ulp at
// fuzz N=10^7). We always compute the DD form; if the cancellation
// ratio |Re|/max_term is tight enough to let that floor dominate, we
// recompute Re through a TD-precision pipeline (exp_full_td ×
// sincos_full_td − 1) rounded back to DD. Threshold picked below.
// Im stays DD either way — no cancellation structure.
complex64x2_t cexpm1dd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 ea_m1 = expm1_full(a);
  float64x2 ea = exp_full(a);
  float64x2 s, c;
  sincos_full(b, s, c);
  float64x2 shalf, chalf;
  sincos_full(float64x2(0.5, 0.0) * b, shalf, chalf);
  float64x2 cos_m1 = -((shalf * shalf) + (shalf * shalf));
  float64x2 term1 = ea_m1 * c;
  float64x2 re_dd = term1 + cos_m1;
  float64x2 im = ea * s;

  // Cancellation detector. |Re_DD|/max_term measures how many bits
  // the two Re terms agreed on before the subtraction — each bit of
  // cancellation costs ~1 DD ulp of DD-path rel-err. Threshold tuned
  // on fuzz_mpfr (N=10^7): 2⁻¹ caps worst-case DD rel-err at ~4.3
  // DD ulp (under the plan's ≤5-ulp ideal) while keeping the TD
  // branch rare enough that amortized runtime lands ~8% above the
  // old pure-DD kernel. Looser (2⁻⁴ → ~28 ulp, ~3% overhead) and
  // tighter (no split → 0.25 ulp, ~78% overhead) both tested.
  double abs_re    = std::fabs(re_dd._limbs[0]);
  double abs_term1 = std::fabs(term1._limbs[0]);
  double abs_cosm1 = std::fabs(cos_m1._limbs[0]);
  double abs_max   = abs_term1 > abs_cosm1 ? abs_term1 : abs_cosm1;
  constexpr double kCancelThresh = 0x1p-1;
  float64x2 re;
  if (abs_re >= kCancelThresh * abs_max || abs_max == 0.0) {
    re = re_dd;
  } else {
    float64x3 ea_td = exp_full_td(a);
    float64x3 s_td, c_td;
    sincos_full_td(b, s_td, c_td);
    float64x3 prod_td = td_mul_td(ea_td, c_td);
    float64x3 diff_td = td_sub_double(prod_td, 1.0);
    re = td_to_dd(diff_td);
  }
  return { to(re), to(im) };
}

// sin(π·z) = sin(π·a) cosh(π·b) + i cos(π·a) sinh(π·b). Using sinpi/cospi
// on the real part keeps full DD precision for half-integer a (where
// sin(π·a) is exactly 0 or ±1 but sin_full(π·a) loses bits to the range
// reduction on π·a). The imaginary scaling π·b is a direct DD multiply.
complex64x2_t csinpidd(complex64x2_t z) {
  static const float64x2 pi_dd(pi_dd_hi, pi_dd_lo);
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 s = sinpi_full(a);
  float64x2 c = cospi_full(a);
  float64x2 pib = pi_dd * b;
  float64x2 sh, ch;
  sinhcosh_full(pib, sh, ch);
  return { to(s * ch), to(c * sh) };
}

// cos(π·z) = cos(π·a) cosh(π·b) − i sin(π·a) sinh(π·b).
complex64x2_t ccospidd(complex64x2_t z) {
  static const float64x2 pi_dd(pi_dd_hi, pi_dd_lo);
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 s = sinpi_full(a);
  float64x2 c = cospi_full(a);
  float64x2 pib = pi_dd * b;
  float64x2 sh, ch;
  sinhcosh_full(pib, sh, ch);
  return { to(c * ch), to(-(s * sh)) };
}

// pow(z, w) = exp(w · log(z)). C99 G.6.4.1; 0^w folds to 0.
complex64x2_t cpowdd(complex64x2_t z, complex64x2_t w) {
  // True zero requires *both* limbs of *both* components to be zero. The
  // previous hi-only test misclassified subnormal DDs of the form
  // (hi=+0, lo=ε) as zero and returned 0 instead of exp(w · log(ε)).
  if (z.re.hi == 0.0 && z.re.lo == 0.0 &&
      z.im.hi == 0.0 && z.im.lo == 0.0) {
    return { {0.0, 0.0}, {0.0, 0.0} };
  }
  // Renormalize each component if its hi limb is zero but lo is not —
  // otherwise clogdd's log2 short-circuit would see hi == 0 and return
  // (−∞, …), propagating NaNs back through exp(w · log(z)). Same pattern
  // as the log1p subnormal fix (audit C5).
  if (z.re.hi == 0.0 && z.re.lo != 0.0) { z.re.hi = z.re.lo; z.re.lo = 0.0; }
  if (z.im.hi == 0.0 && z.im.lo != 0.0) { z.im.hi = z.im.lo; z.im.lo = 0.0; }
  complex64x2_t l = clogdd(z);
  float64x2 lr = from(l.re), li = from(l.im);
  float64x2 wr = from(w.re), wi = from(w.im);
  // w * log(z) = (wr·lr − wi·li) + i·(wr·li + wi·lr)
  complex64x2_t p = { to(wr * lr - wi * li), to(wr * li + wi * lr) };
  return cexpdd(p);
}

// sqrt(z). Principal branch; cut on negative real axis, continuous above.
// Uses Moshier's 2·Re·Im = Im identity to avoid cancellation in mod ± a.
complex64x2_t csqrtdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  if (a._limbs[0] == 0.0 && b._limbs[0] == 0.0)
    return { to(float64x2(0.0)), z.im };  // preserves signed zero of imag
  float64x2 mod = multifloats::hypot(a, b);
  float64x2 half = float64x2(0.5, 0.0);
  if (a._limbs[0] >= 0.0) {
    float64x2 r = multifloats::sqrt((mod + a) * half);
    float64x2 i = b / (r + r);
    return { to(r), to(i) };
  } else {
    float64x2 s = multifloats::sqrt((mod - a) * half);
    float64x2 r = multifloats::abs(b) / (s + s);
    // Use signbit so b = -0 picks the −s branch (IEEE: -0 < 0 is false,
    // so `< 0.0` mishandles the negative-real-axis branch cut per C99 G.6.4.2).
    float64x2 i = std::signbit(b._limbs[0]) ? -s : s;
    return { to(r), to(i) };
  }
}

// sin(a+bi) = sin(a)·cosh(b) + i·cos(a)·sinh(b). 2 fused transcendentals.
complex64x2_t csindd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 sa, ca, sb, cb;
  sincos_full(a, sa, ca);
  sinhcosh_full(b, sb, cb);
  return { to(sa * cb), to(ca * sb) };
}

// cos(a+bi) = cos(a)·cosh(b) − i·sin(a)·sinh(b). 2 fused transcendentals.
complex64x2_t ccosdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 sa, ca, sb, cb;
  sincos_full(a, sa, ca);
  sinhcosh_full(b, sb, cb);
  return { to(ca * cb), to(-(sa * sb)) };
}

// tan(a+bi) = (sin(a)·cos(a) + i·sinh(b)·cosh(b)) / (cos(a)² + sinh(b)²).
// libquadmath ctanq formula; 2 fused transcendentals + real div
// (vs 8 for generic sin(z)/cos(z)).
complex64x2_t ctandd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 sa, ca, sb, cb;
  sincos_full(a, sa, ca);
  sinhcosh_full(b, sb, cb);
  float64x2 den = ca * ca + sb * sb;
  // Pole at a = π/2+kπ with b = 0 ⇒ den = 0 + 0 = 0. C99 Annex G.6.3.1
  // says ctan at a real-axis pole is (NaN, NaN) — matches libquadmath ctanq.
  if (den._limbs[0] == 0.0) {
    double nan = std::numeric_limits<double>::quiet_NaN();
    return { {nan, 0.0}, {nan, 0.0} };
  }
  // Overflow in sinh(b) for |b| ≳ 710 pushes den to +∞. The mathematical
  // limit is tan(a + i·b) → i·sign(b) as |b| → ∞; bypass the division so
  // the imaginary part doesn't come out as ∞/∞ = NaN.
  if (!std::isfinite(den._limbs[0])) {
    float64x2 im_sign(std::copysign(1.0, b._limbs[0]), 0.0);
    return { to(float64x2(0.0, 0.0)), to(im_sign) };
  }
  return { to((sa * ca) / den), to((sb * cb) / den) };
}

// asin(z) = −i · log(i·z + sqrt(1 − z²)). Cuts on real axis for |a|>1.
complex64x2_t casindd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  // 1 − z² : Re = 1 − a² + b², Im = −2ab
  complex64x2_t one_mz2 = { to(float64x2(1.0) - a * a + b * b), to(-(a * b + a * b)) };
  complex64x2_t root = csqrtdd(one_mz2);
  float64x2 rr = from(root.re), ri = from(root.im);
  // i·z + root : (−b + rr) + i·(a + ri)
  complex64x2_t arg = { to(-b + rr), to(a + ri) };
  complex64x2_t l = clogdd(arg);
  // −i · (lr + i·li) = li − i·lr
  complex64x2_t r = { l.im, to(-from(l.re)) };
  // Branch-cut fixup for real-axis inputs with |Re z| > 1 (C99 G.6.2.2):
  // the sign of Im(result) must mirror the sign of Im(input), including
  // ±0. DD multiply doesn't preserve signed zero (2 * -0 → +0), so the
  // sign is lost through the (a*b + a*b) chain. If the result ended up
  // with the wrong sign, negate the WHOLE DD (both limbs) — a simple
  // copysign-per-limb would corrupt the pair because a typical DD lo
  // has the opposite sign from hi (representing the rounding residual).
  // Use a DD-aware sign: for (hi=±0, lo≠0) the hi-limb signbit does not
  // reflect the value's actual sign, so fall through to the lo limb.
  // Route to the lo limb only when hi is an exact zero AND lo carries a
  // nonzero value. For (hi=±0, lo=0) we keep signbit(hi) so ±0 is preserved.
  auto dd_signbit = [](double hi, double lo) {
    return (hi == 0.0 && lo != 0.0) ? std::signbit(lo) : std::signbit(hi);
  };
  if (dd_signbit(r.im.hi, r.im.lo) != dd_signbit(b._limbs[0], b._limbs[1])) {
    r.im.hi = -r.im.hi;
    r.im.lo = -r.im.lo;
  }
  return r;
}

// acos(z) = π/2 − asin(z). Uses DD-precision π/2; the generic <complex>
// template truncates `(_Tp)1.5707963…L` to double on arm64-macOS, losing
// the DD low limb — specializing here fixes that correctness bug.
complex64x2_t cacosdd(complex64x2_t z) {
  static const float64x2 half_pi = float64x2(0x1.921fb54442d18p+0,
                                     0x1.1a62633145c07p-54);
  complex64x2_t s = casindd(z);
  return { to(half_pi - from(s.re)), to(-from(s.im)) };
}

// atan(z) = −i · atanh(i·z). Delegating to catanhdd reuses its Kahan-
// style Im-part (log1p(4·a/den) gated on f ≥ 0.5) and its unit-circle
// denominator picker (dd_x2y2m1 near |z| ≈ 1). The old direct form
// `0.25·log(num/den)` lost ~7 decimals of Im precision when num ≈ den
// (|b| small, |a| moderate); the identity form inherits full DD.
// i·z = (−b + a·i), and DD negation is exact, so the rotation is free.
complex64x2_t catandd(complex64x2_t z) {
  complex64x2_t iz = { to(-from(z.im)), z.re };
  complex64x2_t w = catanhdd(iz);
  return { w.im, to(-from(w.re)) };
}

// sinh(a+bi) = sinh(a)·cos(b) + i·cosh(a)·sin(b). 2 fused transcendentals.
complex64x2_t csinhdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 sa, ca, sb, cb;
  sinhcosh_full(a, sa, ca);
  sincos_full(b, sb, cb);
  return { to(sa * cb), to(ca * sb) };
}

// cosh(a+bi) = cosh(a)·cos(b) + i·sinh(a)·sin(b). 2 fused transcendentals.
complex64x2_t ccoshdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 sa, ca, sb, cb;
  sinhcosh_full(a, sa, ca);
  sincos_full(b, sb, cb);
  return { to(ca * cb), to(sa * sb) };
}

// tanh(a+bi) = (sinh(a)·cosh(a) + i·sin(b)·cos(b)) / (sinh(a)² + cos(b)²).
// libquadmath ctanhq formula.
complex64x2_t ctanhdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 sa, ca, sb, cb;
  sinhcosh_full(a, sa, ca);
  sincos_full(b, sb, cb);
  float64x2 den = sa * sa + cb * cb;
  return { to((sa * ca) / den), to((sb * cb) / den) };
}

// asinh(z). Branch structure ported from libquadmath's
// __quadmath_kernel_casinhq: operates on |Re(z)|, |Im(z)|, restores
// signs at the end. The textbook log(z + sqrt(z² + 1)) form collapses
// to ~double precision once z sits near the unit circle, because
// arg = z + sqrt(1 + z²) and 1 carry the same leading digits — DD's
// 106 bits can't resolve log(|arg|) when |arg|² − 1 ≈ 0. Each branch
// below rearranges the formula so the small result is computed from
// non-cancelling pieces, summed through log1p.
complex64x2_t casinhdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 rx = multifloats::abs(a);
  float64x2 ix = multifloats::abs(b);
  // DD eps = 2⁻¹⁰⁵ ≈ 2.46e-32; constants mirror the __float128 source.
  static constexpr double dd_inv_eps    = 0x1.0p105;
  static constexpr double dd_eps_over_8 = 0x1.0p-108;
  static constexpr double dd_eps_sq     = 0x1.0p-210;
  float64x2 one(1.0, 0.0), half(0.5, 0.0), two(2.0, 0.0);
  float64x2 re, im;

  if (rx._limbs[0] >= dd_inv_eps || ix._limbs[0] >= dd_inv_eps) {
    // |z| huge: z + sqrt(1+z²) ≈ 2z; avoid overflow from squaring.
    complex64x2_t l = clogdd({ to(rx), to(ix) });
    re = from(l.re) + float64x2(ln_2_hi, ln_2_lo);
    im = from(l.im);
  } else if (rx._limbs[0] >= 0.5 && ix._limbs[0] < dd_eps_over_8) {
    float64x2 s = multifloats::hypot(one, rx);
    re = log_full(rx + s);
    im = atan2_full(ix, s);
  } else if (rx._limbs[0] < dd_eps_over_8 && ix._limbs[0] >= 1.5) {
    float64x2 s = multifloats::sqrt((ix + one) * (ix - one));
    re = log_full(ix + s);
    im = atan2_full(s, rx);
  } else if (ix._limbs[0] > 1.0 && ix._limbs[0] < 1.5 && rx._limbs[0] < 0.5) {
    float64x2 ix2m1 = (ix + one) * (ix - one);
    if (rx._limbs[0] < dd_eps_sq) {
      float64x2 s = multifloats::sqrt(ix2m1);
      re = half * log1p_full(two * (ix2m1 + ix * s));
      im = atan2_full(s, rx);
    } else {
      float64x2 rx2 = rx * rx;
      float64x2 f = rx2 * (two + rx2 + two * ix * ix);
      float64x2 d = multifloats::sqrt(ix2m1 * ix2m1 + f);
      float64x2 dp = d + ix2m1;
      float64x2 dm = f / dp;
      float64x2 r1 = multifloats::sqrt((dm + rx2) * half);
      float64x2 r2 = rx * ix / r1;
      re = half * log1p_full(rx2 + dp + two * (rx * r1 + ix * r2));
      im = atan2_full(ix + r2, rx + r1);
    }
  } else if (ix._limbs[0] == 1.0 && ix._limbs[1] == 0.0 && rx._limbs[0] < 0.5) {
    if (rx._limbs[0] < dd_eps_over_8) {
      float64x2 sr = multifloats::sqrt(rx);
      re = half * log1p_full(two * (rx + sr));
      im = atan2_full(one, sr);
    } else {
      float64x2 rx2 = rx * rx;
      float64x2 d = rx * multifloats::sqrt(float64x2(4.0, 0.0) + rx2);
      float64x2 s1 = multifloats::sqrt((d + rx2) * half);
      float64x2 s2 = multifloats::sqrt((d - rx2) * half);
      re = half * log1p_full(rx2 + d + two * (rx * s1 + s2));
      im = atan2_full(one + s2, rx + s1);
    }
  } else if (ix._limbs[0] < 1.0 && rx._limbs[0] < 0.5) {
    if (ix._limbs[0] >= 0x1.0p-105) {
      if (rx._limbs[0] < dd_eps_sq) {
        float64x2 onemix2 = (one + ix) * (one - ix);
        float64x2 s = multifloats::sqrt(onemix2);
        re = half * log1p_full(two * rx / s);
        im = atan2_full(ix, s);
      } else {
        float64x2 onemix2 = (one + ix) * (one - ix);
        float64x2 rx2 = rx * rx;
        float64x2 f = rx2 * (two + rx2 + two * ix * ix);
        float64x2 d = multifloats::sqrt(onemix2 * onemix2 + f);
        float64x2 dp = d + onemix2;
        float64x2 dm = f / dp;
        float64x2 r1 = multifloats::sqrt((dp + rx2) * half);
        float64x2 r2 = rx * ix / r1;
        re = half * log1p_full(rx2 + dm + two * (rx * r1 + ix * r2));
        im = atan2_full(ix + r2, rx + r1);
      }
    } else {
      float64x2 s = multifloats::hypot(one, rx);
      re = half * log1p_full(two * rx * (rx + s));
      im = atan2_full(ix, s);
    }
  } else {
    // General case: rx ≥ 0.5 with ix moderate. Cancellation in
    // 1+z² and log|arg| is bounded here because |arg|² is well away
    // from 1 — the textbook form carries full DD.
    float64x2 yr = (rx - ix) * (rx + ix) + one;
    float64x2 yi = two * rx * ix;
    complex64x2_t root = csqrtdd({ to(yr), to(yi) });
    complex64x2_t arg = { to(rx + from(root.re)), to(ix + from(root.im)) };
    complex64x2_t l = clogdd(arg);
    re = from(l.re);
    im = from(l.im);
  }

  auto dd_signbit = [](double hi, double lo) {
    return (hi == 0.0 && lo != 0.0) ? std::signbit(lo) : std::signbit(hi);
  };
  if (dd_signbit(a._limbs[0], a._limbs[1])) {
    re._limbs[0] = -re._limbs[0];
    re._limbs[1] = -re._limbs[1];
  }
  if (dd_signbit(b._limbs[0], b._limbs[1])) {
    im._limbs[0] = -im._limbs[0];
    im._limbs[1] = -im._limbs[1];
  }
  return { to(re), to(im) };
}

// acosh(z). Kahan's formula: 2·log(sqrt((z+1)/2) + sqrt((z−1)/2)).
// More stable near z≈1 than the naive log(z + sqrt(z²−1)) which suffers
// cancellation in z²−1.
complex64x2_t cacoshdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 half = float64x2(0.5, 0.0);
  complex64x2_t zp = { to((a + float64x2(1.0)) * half), to(b * half) };
  complex64x2_t zm = { to((a - float64x2(1.0)) * half), to(b * half) };
  complex64x2_t s1 = csqrtdd(zp);
  complex64x2_t s2 = csqrtdd(zm);
  complex64x2_t sum = { to(from(s1.re) + from(s2.re)), to(from(s1.im) + from(s2.im)) };
  complex64x2_t l = clogdd(sum);
  return { to(from(l.re) + from(l.re)), to(from(l.im) + from(l.im)) };
}

// atanh(z). Ported from libquadmath catanhq:
//   Re = 0.25 · log((1+a)² + b²)/((1−a)² + b²))
//   Im = 0.5 · atan2(2b, 1 − a² − b²)
// The Re-part uses log1p((4a)/den) when |num/den − 1| is small so the
// small answer survives cancellation in log(num/den). The Im-part's
// denominator picks between (1−|·|)(1+|·|)−(·)² and −dd_x2y2m1 based on
// which region is near the unit circle, avoiding the direct `1 − a² − b²`
// catastrophic subtraction near |z| = 1.
complex64x2_t catanhdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 one(1.0, 0.0), half(0.5, 0.0), quarter(0.25, 0.0);
  // Branch-point singularities at z = ±1 + 0i: C99 G.6.2.4 specifies
  //   catanh(+1 + 0i) = +∞ + 0i
  //   catanh(−1 + 0i) = −∞ + 0i
  // Without the short-circuit, num/den = 4/0 = inf in DD; the subsequent
  // 0.25 * inf DD multiply emits fma(0.25, inf, -inf) = NaN, producing
  // a wrong NaN-result at this exact branch point.
  if (b._limbs[0] == 0.0 && b._limbs[1] == 0.0 &&
      std::abs(a._limbs[0]) == 1.0 && a._limbs[1] == 0.0) {
    double inf = std::numeric_limits<double>::infinity();
    complex64x2_t r;
    r.re.hi = std::copysign(inf, a._limbs[0]);
    r.re.lo = 0.0;
    r.im.hi = 0.0;
    r.im.lo = 0.0;
    return r;
  }

  static constexpr double dd_eps_sq  = 0x1.0p-210;
  static constexpr double dd_eps_hlf = 0x1.0p-106;

  float64x2 absa = multifloats::abs(a);
  float64x2 absb = multifloats::abs(b);
  float64x2 re_res;

  if (absa._limbs[0] == 1.0 && absa._limbs[1] == 0.0 &&
      absb._limbs[0] < dd_eps_sq) {
    // Near-branch-point: catanh(±1 + ε·i) ≈ sign(a) · 0.5 · (ln2 − log|ε|).
    float64x2 ln2(ln_2_hi, ln_2_lo);
    float64x2 t = half * (ln2 - log_full(absb));
    re_res = (a._limbs[0] < 0.0) ? -t : t;
  } else {
    // i2: drop b² when |b| below eps² so the DD round-off of b² doesn't
    // contaminate num/den. Matches the libquadmath guard.
    float64x2 i2 = (absb._limbs[0] >= dd_eps_sq) ? (b * b)
                                                 : float64x2(0.0, 0.0);
    float64x2 ap1 = one + a;
    float64x2 am1 = one - a;
    float64x2 num = i2 + ap1 * ap1;
    float64x2 den = i2 + am1 * am1;
    float64x2 f = num / den;
    if (f._limbs[0] < 0.5) {
      re_res = quarter * log_full(f);
    } else {
      // 4a exact in DD (power-of-two scale); computing as (num−den)/den
      // would cancel when a is small — instead feed 4a/den into log1p.
      float64x2 four_a = (a + a) + (a + a);
      re_res = quarter * log1p_full(four_a / den);
    }
  }

  // Imag-part denominator picks the non-cancelling form per region.
  float64x2 big   = (absa._limbs[0] > absb._limbs[0]) ? absa : absb;
  float64x2 small_= (absa._limbs[0] > absb._limbs[0]) ? absb : absa;
  float64x2 den_im;
  if (small_._limbs[0] < dd_eps_hlf) {
    den_im = (one - big) * (one + big);
  } else if (big._limbs[0] >= 1.0) {
    den_im = (one - big) * (one + big) - small_ * small_;
  } else if (big._limbs[0] >= 0.75 || small_._limbs[0] >= 0.5) {
    // Near unit circle: 1 − big² − small² = −dd_x2y2m1(big, small).
    float64x2 t = dd_x2y2m1(big, small_);
    den_im._limbs[0] = -t._limbs[0];
    den_im._limbs[1] = -t._limbs[1];
  } else {
    den_im = (one - big) * (one + big) - small_ * small_;
  }
  float64x2 im_res = half * atan2_full(b + b, den_im);

  return { to(re_res), to(im_res) };
}

// |z| via the scale-aware path shared with clogdd: avoids overflow when
// |z| exceeds DBL_MAX even though hypot(a,b) itself would saturate.
float64x2_t cabsdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  return to(multifloats::hypot(a, b));
}

// arg(z) = atan2(im, re), with the standard branch cut on the negative
// real axis; signed zero in im selects ±pi.
float64x2_t cargdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  return to(atan2_full(b, a));
}

// Riemann-sphere projection (C99 7.3.9.4): if either component is
// infinite the result is (+inf, copysign(0, im)); otherwise z is
// returned unchanged. NaN on a finite part is left alone.
complex64x2_t cprojdd(complex64x2_t z) {
  if (std::isinf(z.re.hi) || std::isinf(z.im.hi)) {
    float64x2_t inf = { std::numeric_limits<double>::infinity(), 0.0 };
    float64x2_t zim = { std::copysign(0.0, z.im.hi), 0.0 };
    return { inf, zim };
  }
  return z;
}

// Conjugate: flip the sign of every limb of the imaginary part so the
// DD value negates exactly (no roundoff, preserves the low limb).
complex64x2_t conjdd(complex64x2_t z) {
  return { z.re, { -z.im.hi, -z.im.lo } };
}

// Real/imag accessors. Trivial but exposed as linkable symbols for
// parity with libquadmath's crealq / cimagq.
float64x2_t crealdd(complex64x2_t z) { return z.re; }
float64x2_t cimagdd(complex64x2_t z) { return z.im; }

// Matrix multiply (column-major).
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
