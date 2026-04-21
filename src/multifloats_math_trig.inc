// ---- sinpi / cospi / tanpi / sin / cos / tan -------------------------------
// Family of π-scaled trig: *pi(x) = fn(π·x). Range-reduce 2·|x| to the
// nearest integer n, leaving |rx| ≤ 0.25, so |π·rx| ≤ π/4 matches the
// sin/cos eval range. Dispatch by n mod 4 for the quadrant.
//
// Large-argument cap: for |x| ≥ 2^52 a double cannot resolve the fractional
// part of 2x (ulp ≥ 1), so n mod 4 is meaningless and the result is
// mathematically indeterminate. We also need the cap to stop 2·x from
// overflowing to ±∞, which would make the subsequent static_cast<long long>
// undefined behavior. Return NaN.
namespace {
constexpr double kPiTrigArgumentCap = 0x1p52;
inline bool pi_trig_arg_too_large(double a) { return !(a < kPiTrigArgumentCap); }
inline float64x2 pi_trig_nan() {
  float64x2 r;
  r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
  r._limbs[1] = 0.0;
  return r;
}
// sin/cos/tan argument cap: through |x| ≤ 2^54, reduce_pi_half stays
// sub-1-ulp; at |x| ≥ 2^55 the n·π/2 reduction walks off because
// nearbyint(x·2/π) loses integer precision (ulp(x) ≥ 1). For normalized
// DD with |hi| > 2^54 the lo limb (|lo| ≤ 0.5 ulp(hi), ≤ 2) carries no
// sub-integer information, so the input itself is effectively an
// integer — the DD library cannot do better than what std::sin(double)
// returns via Payne–Hanek, which is outside this library's scope. Match
// the sinpi cap and return NaN.
constexpr double kTrigArgumentCap = 0x1p55;
inline bool trig_arg_too_large(double a) {
  return !(std::fabs(a) < kTrigArgumentCap);
}
}

float64x2 sinpi_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) return pi_trig_nan();
  bool sign = x._limbs[0] < 0.0;
  float64x2 ax = sign ? -x : x;
  if (pi_trig_arg_too_large(ax._limbs[0])) return pi_trig_nan();
  double n_float = std::nearbyint(2.0 * ax._limbs[0]);
  float64x2 rx = ax + float64x2(-0.5 * n_float, 0.0);
  float64x2 pi_dd = float64x2(pi_dd_hi, pi_dd_lo);
  float64x2 pr = pi_dd * rx;
  long long n_mod = static_cast<long long>(n_float) & 3LL;
  float64x2 res;
  switch (n_mod) {
  case 0: res = sin_eval(pr); break;
  case 1: res = cos_eval(pr); break;
  case 2: res = -sin_eval(pr); break;
  default: res = -cos_eval(pr); break;
  }
  return sign ? -res : res;
}

// cospi is even; no sign flip on the result.
float64x2 cospi_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) return pi_trig_nan();
  float64x2 ax = (x._limbs[0] < 0.0) ? -x : x;
  if (pi_trig_arg_too_large(ax._limbs[0])) return pi_trig_nan();
  double n_float = std::nearbyint(2.0 * ax._limbs[0]);
  float64x2 rx = ax + float64x2(-0.5 * n_float, 0.0);
  float64x2 pi_dd = float64x2(pi_dd_hi, pi_dd_lo);
  float64x2 pr = pi_dd * rx;
  long long n_mod = static_cast<long long>(n_float) & 3LL;
  switch (n_mod) {
  case 0: return cos_eval(pr);
  case 1: return -sin_eval(pr);
  case 2: return -cos_eval(pr);
  default: return sin_eval(pr);
  }
}

// tanpi: fuses sin/cos eval so both come from one Taylor pair. Poles at
// half-integers fall out naturally (sin(0)=+0 ⇒ ±c/0 = ±inf).
float64x2 tanpi_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) return pi_trig_nan();
  bool sign = x._limbs[0] < 0.0;
  float64x2 ax = sign ? -x : x;
  if (pi_trig_arg_too_large(ax._limbs[0])) return pi_trig_nan();
  double n_float = std::nearbyint(2.0 * ax._limbs[0]);
  float64x2 rx = ax + float64x2(-0.5 * n_float, 0.0);
  float64x2 pi_dd = float64x2(pi_dd_hi, pi_dd_lo);
  float64x2 pr = pi_dd * rx;
  float64x2 s, c;
  sincos_eval(pr, s, c);
  long long n_mod = static_cast<long long>(n_float) & 3LL;
  // n_mod 0, 2: tan(n·π/2 + π·rx) = tan(π·rx) = s/c
  // n_mod 1, 3: tan(n·π/2 + π·rx) = -cot(π·rx) = -c/s
  float64x2 res = (n_mod & 1LL) ? -c / s : s / c;
  return sign ? -res : res;
}

// Cody-Waite range reduction: n = nearest(x·2/π), r = x − n·(π/2) using a
// 3-part π/2 constant (three back-to-back doubles). Each n·pi_half_cwK is
// an exact DD product via FMA, then subtracted from r via Shewchuk/Bailey
// IEEE-style DD-add that carries both lo-limb contributions through a
// two_sum pair — error ≤ 4·eps²·|r|. The earlier `rl += e - pl` inline
// summed three same-magnitude terms in plain double, injecting ~eps·|pl|
// per step (~800 DD-ulp floor at x near multiples of π/2 when |n·pi_half|
// ≫ |r|, e.g. sin(47.15) near 15π).
//
// Guards (non-finite, zero) from the public DD operator- are elided: the
// caller (sincos_full etc.) pre-checks both conditions, so the kernel
// runs on a finite non-zero accumulator.
void reduce_pi_half(float64x2 const &x, float64x2 &r, int &n_mod4) {
  using multifloats::detail::two_sum;
  using multifloats::detail::fast_two_sum;
  double n_float = std::nearbyint(x._limbs[0] * 0.6366197723675814);
  double rh = x._limbs[0];
  double rl = x._limbs[1];
  double ph, pl;
  double sh, sl, tl, tll;
  // One IEEE DD-subtract step: (rh, rl) ← (rh, rl) - (ph, pl).
  auto dd_sub_step = [&](double ph, double pl) {
    two_sum(rh, -ph, sh, sl);        // sh + sl = rh - ph exact
    two_sum(rl, -pl, tl, tll);       // tl + tll = rl - pl exact
    sl += tl;                        // fold lo-limb hi part (1 rounding)
    fast_two_sum(sh, sl, rh, rl);    // renormalize
    rl += tll;                       // fold lo-limb tail (1 rounding)
    fast_two_sum(rh, rl, rh, rl);    // final renormalize
  };
  // Step 1: subtract n·cw1 (exact DD via FMA).
  ph = n_float * pi_half_cw1;
  pl = std::fma(n_float, pi_half_cw1, -ph);
  dd_sub_step(ph, pl);
  // Step 2: subtract n·cw2.
  ph = n_float * pi_half_cw2;
  pl = std::fma(n_float, pi_half_cw2, -ph);
  dd_sub_step(ph, pl);
  // Step 3: cw3 small enough that n·cw3 fits one double exactly (pl ≡ 0).
  ph = n_float * pi_half_cw3;
  dd_sub_step(ph, 0.0);
  r._limbs[0] = rh;
  r._limbs[1] = rl;
  long long nn = static_cast<long long>(std::fabs(n_float)) & 3LL;
  if (n_float < 0.0) nn = (4 - nn) & 3;
  n_mod4 = static_cast<int>(nn);
}

// Taylor kernels for |x| ≤ π/8. sin evaluates as x·poly(x²) so the low
// limb of x is preserved losslessly.
float64x2 sin_kernel(float64x2 const &x) {
  float64x2 x2 = x * x;
  // 13 coefs → degree 12 → Estrin case in neval.
  return neval(x2, sin_taylor_hi, sin_taylor_lo, 12) * x;
}
float64x2 cos_kernel(float64x2 const &x) {
  float64x2 x2 = x * x;
  // 13 coefs → degree 12 → Estrin case in neval.
  return neval(x2, cos_taylor_hi, cos_taylor_lo, 12);
}

// Evaluate sin(r)/cos(r) for |r| ≤ π/4. For |r| > π/8 shift by π/4 and use
// the angle-addition identity — this halves the polynomial argument range
// (x² ≤ 0.154 vs 0.616) so the 13-term Taylor hits full DD at the boundary.
float64x2 sin_eval(float64x2 const &r) {
  constexpr double pi8 = 0.392699081698724;
  if (std::fabs(r._limbs[0]) <= pi8) return sin_kernel(r);
  float64x2 pi4_dd = float64x2(pi_quarter_hi, pi_quarter_lo);
  float64x2 inv_sqrt2 = float64x2(inv_sqrt2_hi, inv_sqrt2_lo);
  bool pos = r._limbs[0] > 0.0;
  float64x2 rp = pos ? (r - pi4_dd) : (r + pi4_dd);
  float64x2 sk = sin_kernel(rp);
  float64x2 ck = cos_kernel(rp);
  return pos ? (sk + ck) * inv_sqrt2 : (sk - ck) * inv_sqrt2;
}

float64x2 cos_eval(float64x2 const &r) {
  constexpr double pi8 = 0.392699081698724;
  if (std::fabs(r._limbs[0]) <= pi8) return cos_kernel(r);
  float64x2 pi4_dd = float64x2(pi_quarter_hi, pi_quarter_lo);
  float64x2 inv_sqrt2 = float64x2(inv_sqrt2_hi, inv_sqrt2_lo);
  bool pos = r._limbs[0] > 0.0;
  float64x2 rp = pos ? (r - pi4_dd) : (r + pi4_dd);
  float64x2 sk = sin_kernel(rp);
  float64x2 ck = cos_kernel(rp);
  return pos ? (ck - sk) * inv_sqrt2 : (ck + sk) * inv_sqrt2;
}

// Fused sin/cos of |r| ≤ π/4 — shares the π/4 shift path and its two
// Taylor kernels between both outputs. Callers that need both s=sin(r)
// and c=cos(r) save one 13-term Horner (~40% of the eval cost).
void sincos_eval(float64x2 const &r, float64x2 &s, float64x2 &c) {
  constexpr double pi8 = 0.392699081698724;
  if (std::fabs(r._limbs[0]) <= pi8) {
    s = sin_kernel(r);
    c = cos_kernel(r);
    return;
  }
  float64x2 pi4_dd = float64x2(pi_quarter_hi, pi_quarter_lo);
  float64x2 inv_sqrt2 = float64x2(inv_sqrt2_hi, inv_sqrt2_lo);
  bool pos = r._limbs[0] > 0.0;
  float64x2 rp = pos ? (r - pi4_dd) : (r + pi4_dd);
  float64x2 sk = sin_kernel(rp);
  float64x2 ck = cos_kernel(rp);
  if (pos) {
    s = (sk + ck) * inv_sqrt2;
    c = (ck - sk) * inv_sqrt2;
  } else {
    s = (sk - ck) * inv_sqrt2;
    c = (ck + sk) * inv_sqrt2;
  }
}

// Fused sin/cos of the full-range argument — one call to reduce_pi_half
// plus one sincos_eval, with the per-quadrant sign/swap applied last.
void sincos_full(float64x2 const &x, float64x2 &s, float64x2 &c) {
  if (!std::isfinite(x._limbs[0]) || trig_arg_too_large(x._limbs[0])) {
    float64x2 nan;
    nan._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    nan._limbs[1] = 0.0;
    s = nan;
    c = nan;
    return;
  }
  float64x2 r;
  int q;
  reduce_pi_half(x, r, q);
  float64x2 ss, cc;
  sincos_eval(r, ss, cc);
  switch (q) {
    case 0: s = ss;  c = cc;  break;
    case 1: s = cc;  c = -ss; break;
    case 2: s = -ss; c = -cc; break;
    default: s = -cc; c = ss; break;
  }
}

float64x2 sin_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0]) || trig_arg_too_large(x._limbs[0])) return pi_trig_nan();
  float64x2 r;
  int q;
  reduce_pi_half(x, r, q);
  switch (q) {
  case 0: return sin_eval(r);
  case 1: return cos_eval(r);
  case 2: return -sin_eval(r);
  default: return -cos_eval(r);
  }
}

float64x2 cos_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0]) || trig_arg_too_large(x._limbs[0])) return pi_trig_nan();
  float64x2 r;
  int q;
  reduce_pi_half(x, r, q);
  switch (q) {
  case 0: return cos_eval(r);
  case 1: return -sin_eval(r);
  case 2: return -cos_eval(r);
  default: return sin_eval(r);
  }
}

float64x2 tan_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0]) || trig_arg_too_large(x._limbs[0])) return pi_trig_nan();
  float64x2 r;
  int q;
  reduce_pi_half(x, r, q);
  float64x2 s, c;
  sincos_eval(r, s, c);
  switch (q) {
    case 0:
    case 2:  return s / c;
    default: return -c / s;  // q == 1 or 3
  }
}

// ---- TD-output sincos variants --------------------------------------------
// Mirror of reduce_pi_half / sincos_full that carries ~159 bits through
// the range reduction and Taylor evaluation. Every constant (π/2 split
// cw1..cw4, pi_quarter, inv_sqrt2, sin_taylor, cos_taylor) is lifted to
// TD precision in dd_constants.hh. See doc/developer/PLAN-cexpm1-td-
// internals.md Phase B for motivation — the mock-TD-cos prototype
// stalled at the DD-cos floor of 2⁻¹⁰⁵, which this kernel removes.

void reduce_pi_half_td(float64x2 const &x, float64x3 &r, int &n_mod4) {
  double n_float = std::nearbyint(x._limbs[0] * 0.6366197723675814);
  // TSUM accumulator: start from x (DD) and subtract n·π/2 in four
  // Cody-Waite stages. Every n·cwk is captured as an exact FMA pair so
  // |n| can grow up to kTrigArgumentCap (2^55) without losing the
  // sub-double residual. Using the same FMA for cw1 is the fix for a
  // bug where |n·cw1| exceeded the 33-bit significand's exact range at
  // |n| > 2^20; the DD path hides this via dd_sub_step's FMA, but the
  // single-limb `-(n_float*cw1)` shortcut here does not.
  double T0 = 0.0, T1 = 0.0, T2 = 0.0;
  tsum3(T0, T1, T2, x._limbs[0]);
  tsum3(T0, T1, T2, x._limbs[1]);
  auto sub_fma_pair = [&](double cw) {
    double ph = n_float * cw;
    double pl = std::fma(n_float, cw, -ph);
    tsum3(T0, T1, T2, -ph);
    tsum3(T0, T1, T2, -pl);
  };
  sub_fma_pair(pi_half_cw1);
  sub_fma_pair(pi_half_cw2);
  sub_fma_pair(pi_half_cw3);
  sub_fma_pair(pi_half_cw4);
  tsum3_finalize(T0, T1, T2);
  r = {T0, T1, T2};
  long long nn = static_cast<long long>(std::fabs(n_float)) & 3LL;
  if (n_float < 0.0) nn = (4 - nn) & 3;
  n_mod4 = static_cast<int>(nn);
}

// Taylor kernels at TD precision, |x| ≤ π/8. sin keeps the factored
// sin(x)/x form so the low limbs of x survive the outer multiply.
float64x3 sin_kernel_td(float64x3 const &x) {
  float64x3 x2 = td_mul_td(x, x);
  float64x3 q = {sin_taylor_hi[12], sin_taylor_lo[12], sin_taylor_lo2[12]};
  for (int j = 11; j >= 0; --j) {
    q = td_mul_td(q, x2);
    float64x3 c = {sin_taylor_hi[j], sin_taylor_lo[j], sin_taylor_lo2[j]};
    q = td_add_td(q, c);
  }
  return td_mul_td(q, x);
}
float64x3 cos_kernel_td(float64x3 const &x) {
  float64x3 x2 = td_mul_td(x, x);
  float64x3 q = {cos_taylor_hi[12], cos_taylor_lo[12], cos_taylor_lo2[12]};
  for (int j = 11; j >= 0; --j) {
    q = td_mul_td(q, x2);
    float64x3 c = {cos_taylor_hi[j], cos_taylor_lo[j], cos_taylor_lo2[j]};
    q = td_add_td(q, c);
  }
  return q;
}

// Fused TD sin/cos for |r| ≤ π/4 — mirrors sincos_eval with TD constants
// throughout. For |r| > π/8, shift by π/4 so the Taylor arg stays inside
// [-π/8, π/8] and the 13-term polynomial hits full TD at the edge.
void sincos_eval_td(float64x3 const &r, float64x3 &s, float64x3 &c) {
  constexpr double pi8 = 0.392699081698724;
  if (std::fabs(r._limbs[0]) <= pi8) {
    s = sin_kernel_td(r);
    c = cos_kernel_td(r);
    return;
  }
  float64x3 pi4_td = {pi_quarter_hi, pi_quarter_lo, pi_quarter_lo2};
  float64x3 inv_sqrt2_td = {inv_sqrt2_hi, inv_sqrt2_lo, inv_sqrt2_lo2};
  bool pos = r._limbs[0] > 0.0;
  float64x3 rp = pos ? td_add_td(r, td_negate(pi4_td))
                     : td_add_td(r, pi4_td);
  float64x3 sk = sin_kernel_td(rp);
  float64x3 ck = cos_kernel_td(rp);
  if (pos) {
    s = td_mul_td(td_add_td(sk, ck), inv_sqrt2_td);
    c = td_mul_td(td_add_td(ck, td_negate(sk)), inv_sqrt2_td);
  } else {
    s = td_mul_td(td_add_td(sk, td_negate(ck)), inv_sqrt2_td);
    c = td_mul_td(td_add_td(ck, sk), inv_sqrt2_td);
  }
}

void sincos_full_td(float64x2 const &x, float64x3 &s, float64x3 &c) {
  if (!std::isfinite(x._limbs[0]) || trig_arg_too_large(x._limbs[0])) {
    double nan = std::numeric_limits<double>::quiet_NaN();
    s = {nan, 0.0, 0.0};
    c = {nan, 0.0, 0.0};
    return;
  }
  float64x3 r;
  int q;
  reduce_pi_half_td(x, r, q);
  float64x3 ss, cc;
  sincos_eval_td(r, ss, cc);
  switch (q) {
    case 0: s = ss;              c = cc;              break;
    case 1: s = cc;              c = td_negate(ss);   break;
    case 2: s = td_negate(ss);   c = td_negate(cc);   break;
    default: s = td_negate(cc);  c = ss;              break;
  }
}
