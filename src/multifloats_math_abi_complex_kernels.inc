// Complex overloads for the Kind-D parity gap functions — operations the
// C ABI exposes via `c*dd` kernels but which have no std::complex free
// function to specialize (expm1, log2, log1p) or whose name is outside
// the standard library (sinpi, cospi).
//
// log10 has a std::complex free-function template in C++17; placing it
// as a std:: specialization in multifloats_math_complex_std.inc lets
// `std::log10(z)` route to the DD path. The other five live here.
//
// Bodies mirror the original extern "C" c*dd implementations in
// multifloats_math_abi_complex.inc; the C shim wrappers over there now
// call these directly.

// log10(z) = log(z) · (1/ln 10). Routes through std::log (the DD complex-
// log kernel) so the DD precision profile matches the scalar path.
std::complex<float64x2> log10(std::complex<float64x2> const &z) {
  static const float64x2 inv_ln10 = float64x2(0x1.bcb7b1526e50ep-2,
                                              0x1.95355baaafad3p-57);
  std::complex<float64x2> l = std::log(z);
  return {l.real() * inv_ln10, l.imag() * inv_ln10};
}

// log2(z) = log(z) · (1/ln 2).
std::complex<float64x2> log2(std::complex<float64x2> const &z) {
  static const float64x2 inv_ln2 =
      float64x2(1.4426950408889634, 2.0355273740931033e-17);
  std::complex<float64x2> l = std::log(z);
  return {l.real() * inv_ln2, l.imag() * inv_ln2};
}

// log1p(z) = log(1 + z) without losing precision for |z| << 1. Uses
// log1p on the real scalar form |1+z|² - 1 = 2a + a² + b², which is
// computed without the 1 ± 1 cancellation. The imaginary part is
// atan2(b, 1 + a) — atan2 handles the branch cut and signed-zero
// propagation directly. Matches libquadmath clog1pq.
std::complex<float64x2> log1p(std::complex<float64x2> const &z) {
  float64x2 a = z.real(), b = z.imag();
  float64x2 u = (a + a) + a * a + b * b;          // 2a + a² + b²
  float64x2 half(0.5, 0.0);
  float64x2 r = half * log1p(u);
  float64x2 phi = atan2(b, float64x2(1.0, 0.0) + a);
  return {r, phi};
}

// expm1(z) = e^z - 1. The half-angle DD form
//   Re = expm1(a)·cos(b) − 2·sin²(b/2),  Im = e^a·sin(b)
// gives full DD precision *except* on the cancellation surface
// cos(b)·e^a = 1, where the two Re terms agree to many bits and the
// DD abs-err floor ~2⁻¹⁰⁵·max_term bleeds through (1078 DD ulp at
// fuzz N=10^7). We always compute the DD form; if the cancellation
// ratio |Re|/max_term is tight enough to let that floor dominate, we
// recompute Re through a TD-precision pipeline (exp_full_td ×
// sincos_td − 1) rounded back to DD. Im stays DD either way — no
// cancellation structure.
std::complex<float64x2> expm1(std::complex<float64x2> const &z) {
  float64x2 a = z.real(), b = z.imag();
  float64x2 ea_m1 = expm1(a);
  float64x2 ea = exp(a);
  float64x2 s, c;
  sincos(b, s, c);
  float64x2 shalf, chalf;
  sincos(float64x2(0.5, 0.0) * b, shalf, chalf);
  float64x2 cos_m1 = -((shalf * shalf) + (shalf * shalf));
  float64x2 term1 = ea_m1 * c;
  float64x2 re_dd = term1 + cos_m1;
  float64x2 im = ea * s;

  double abs_re    = std::fabs(re_dd._limbs[0]);
  double abs_term1 = std::fabs(term1._limbs[0]);
  double abs_cosm1 = std::fabs(cos_m1._limbs[0]);
  double abs_max   = abs_term1 > abs_cosm1 ? abs_term1 : abs_cosm1;
  constexpr double kCancelThresh = 0x1p-1;
  float64x2 re;
  if (abs_re >= kCancelThresh * abs_max || abs_max == 0.0) {
    re = re_dd;
  } else {
    detail::float64x3 ea_td = exp_full_td(a);
    detail::float64x3 s_td, c_td;
    sincos_td(b, s_td, c_td);
    detail::float64x3 prod_td = detail::td_mul_td(ea_td, c_td);
    detail::float64x3 diff_td = detail::td_sub_double(prod_td, 1.0);
    re = detail::td_to_dd(diff_td);
  }
  return {re, im};
}

// Compute (π·b) to triple-double precision, returning the leading DD pair
// `pib_dd` and a residual `pib_lo` (a single double). Used by complex
// sinpi/cospi to avoid the linear-in-|π·b| relative-error amplification
// when sinh/cosh are evaluated at a DD argument: a plain `pi_dd * b` carries
// ~ulp_dd of relative input error which sinh blows up to ~|π·b|·ulp_dd of
// output relative error (≈ 88 DD ulp for |b|≈30 in fuzz). Capturing the
// 159-bit residual and applying the linear correction `sh += ch·pib_lo`
// recovers full DD precision.
inline void pi_times_dd_td(float64x2 const &b, float64x2 &pib_dd, double &pib_lo) {
  float64x3 pi_td{pi_dd_hi, pi_dd_lo, pi_dd_lo2};
  float64x3 t = td_mul_dd(pi_td, b);
  pib_dd = float64x2(t._limbs[0], t._limbs[1]);
  pib_lo = t._limbs[2];
}

// sin(π·z) = sin(π·a) cosh(π·b) + i cos(π·a) sinh(π·b).
std::complex<float64x2> sinpi(std::complex<float64x2> const &z) {
  float64x2 a = z.real(), b = z.imag();
  float64x2 s = sinpi(a);
  float64x2 c = cospi(a);
  float64x2 pib_dd; double pib_lo;
  pi_times_dd_td(b, pib_dd, pib_lo);
  float64x2 sh, ch;
  sinhcosh(pib_dd, sh, ch);
  // Linear correction for the bits beyond DD: sinh(x+ε) ≈ sinh + cosh·ε,
  // cosh(x+ε) ≈ cosh + sinh·ε. ε = pib_lo is small enough that higher-order
  // terms are below DD precision.
  float64x2 sh_new = sh + ch * float64x2(pib_lo);
  float64x2 ch_new = ch + sh * float64x2(pib_lo);
  return {s * ch_new, c * sh_new};
}

// cos(π·z) = cos(π·a) cosh(π·b) − i sin(π·a) sinh(π·b).
std::complex<float64x2> cospi(std::complex<float64x2> const &z) {
  float64x2 a = z.real(), b = z.imag();
  float64x2 s = sinpi(a);
  float64x2 c = cospi(a);
  float64x2 pib_dd; double pib_lo;
  pi_times_dd_td(b, pib_dd, pib_lo);
  float64x2 sh, ch;
  sinhcosh(pib_dd, sh, ch);
  float64x2 sh_new = sh + ch * float64x2(pib_lo);
  float64x2 ch_new = ch + sh * float64x2(pib_lo);
  return {c * ch_new, -(s * sh_new)};
}
