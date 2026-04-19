// ---- sinpi / cospi / tanpi / sin / cos / tan -------------------------------
// Family of π-scaled trig: *pi(x) = fn(π·x). Range-reduce 2·|x| to the
// nearest integer n, leaving |rx| ≤ 0.25, so |π·rx| ≤ π/4 matches the
// sin/cos eval range. Dispatch by n mod 4 for the quadrant.
float64x2 sinpi_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  bool sign = x._limbs[0] < 0.0;
  float64x2 ax = sign ? -x : x;
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
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  float64x2 ax = (x._limbs[0] < 0.0) ? -x : x;
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
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  bool sign = x._limbs[0] < 0.0;
  float64x2 ax = sign ? -x : x;
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
// an exact DD product via FMA, then subtracted from r as full DD; the
// accumulator carries the full ~106 bits DD can hold — the critical
// difference from the previous sinpi(x·inv_pi) path, which lost the
// integer part of x/π for large |x|.
//
// Inlined as three (two_sum + low-limb-merge + fast_two_sum) steps rather
// than three full DD operator- calls: operator- does two independent
// two_sums (one per limb) then two more fast_two_sums, and guards
// non-finite/zero inputs. Here we know the accumulator stays finite and
// never reaches the zero-zero branch, so we can fold the low-limb into
// one fast_two_sum per step (~11 FLOPs instead of ~20).
void reduce_pi_half(float64x2 const &x, float64x2 &r, int &n_mod4) {
  using multifloats::detail::two_sum;
  using multifloats::detail::fast_two_sum;
  double n_float = std::nearbyint(x._limbs[0] * 0.6366197723675814);
  double rh = x._limbs[0];
  double rl = x._limbs[1];
  double t, e, ph, pl;
  // Step 1: subtract n·cw1 (exact DD via FMA) from (rh, rl).
  ph = n_float * pi_half_cw1;
  pl = std::fma(n_float, pi_half_cw1, -ph);
  two_sum(rh, -ph, t, e);
  rl += e - pl;
  fast_two_sum(t, rl, rh, rl);
  // Step 2: subtract n·cw2.
  ph = n_float * pi_half_cw2;
  pl = std::fma(n_float, pi_half_cw2, -ph);
  two_sum(rh, -ph, t, e);
  rl += e - pl;
  fast_two_sum(t, rl, rh, rl);
  // Step 3: cw3 is single-double (pl ≡ 0), so skip the FMA split.
  ph = n_float * pi_half_cw3;
  two_sum(rh, -ph, t, e);
  rl += e;
  fast_two_sum(t, rl, rh, rl);
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
  float64x2 pi4_dd = float64x2(0.7853981633974483, 3.061616997868383e-17);
  float64x2 inv_sqrt2 = float64x2(0.7071067811865476, -4.833646656726457e-17);
  bool pos = r._limbs[0] > 0.0;
  float64x2 rp = pos ? (r - pi4_dd) : (r + pi4_dd);
  float64x2 sk = sin_kernel(rp);
  float64x2 ck = cos_kernel(rp);
  return pos ? (sk + ck) * inv_sqrt2 : (sk - ck) * inv_sqrt2;
}

float64x2 cos_eval(float64x2 const &r) {
  constexpr double pi8 = 0.392699081698724;
  if (std::fabs(r._limbs[0]) <= pi8) return cos_kernel(r);
  float64x2 pi4_dd = float64x2(0.7853981633974483, 3.061616997868383e-17);
  float64x2 inv_sqrt2 = float64x2(0.7071067811865476, -4.833646656726457e-17);
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
  float64x2 pi4_dd = float64x2(0.7853981633974483, 3.061616997868383e-17);
  float64x2 inv_sqrt2 = float64x2(0.7071067811865476, -4.833646656726457e-17);
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
  if (!std::isfinite(x._limbs[0])) {
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
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
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
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
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
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
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
