#include "multifloats.hh"
#include "dd_constants.hh"
#include <cstdint>
#include <cstring>
#include <vector>

// All DD math kernels have internal linkage (anonymous namespace).
// Only the extern "C" `*dd` wrappers at the bottom are exported.
namespace {

using namespace multifloats::detail;  // neval, deval, horner
using multifloats::float64x2;

// Forward declarations for internal cross-references
float64x2 exp_full(float64x2 const &x);
float64x2 exp2_full(float64x2 const &x);
float64x2 expm1_full(float64x2 const &x);
float64x2 log_full(float64x2 const &x);
float64x2 log1p_full(float64x2 const &x);
float64x2 sin_full(float64x2 const &x);
float64x2 cos_full(float64x2 const &x);
float64x2 sin_eval(float64x2 const &r);
float64x2 cos_eval(float64x2 const &r);
void sincos_eval(float64x2 const &r, float64x2 &s, float64x2 &c);
void sincos_full(float64x2 const &x, float64x2 &s, float64x2 &c);
void sinhcosh_full(float64x2 const &x, float64x2 &s, float64x2 &c);
float64x2 sinpi_full(float64x2 const &x);
float64x2 erfc_full(float64x2 const &x);
float64x2 lgamma_positive(float64x2 const &x);
float64x2 bessel_j0_full(float64x2 const &x);
float64x2 bessel_j1_full(float64x2 const &x);

// ---- exp2 / exp / log2 / log / log10 ---------------------------------------
float64x2 exp2_kernel(float64x2 const &x) {
  // n = nearest integer to x.hi, y = (x - n)/8
  double n_float = std::nearbyint(x._limbs[0]);
  float64x2 y = x + float64x2(-n_float, 0.0);
  y._limbs[0] = std::ldexp(y._limbs[0], -3);
  y._limbs[1] = std::ldexp(y._limbs[1], -3);
  // poly(y) then cube via three squarings (to undo the /8)
  float64x2 p = horner(y, exp2_coefs_hi, exp2_coefs_lo, 14);
  p = p * p;
  p = p * p;
  p = p * p;
  // multiply by 2^n via two ldexps (split to keep intermediates in range)
  int n = static_cast<int>(n_float);
  int half_n = n / 2;
  float64x2 r;
  r._limbs[0] = std::ldexp(p._limbs[0], half_n);
  r._limbs[1] = std::ldexp(p._limbs[1], half_n);
  r._limbs[0] = std::ldexp(r._limbs[0], n - half_n);
  r._limbs[1] = std::ldexp(r._limbs[1], n - half_n);
  return r;
}

float64x2 exp2_full(float64x2 const &x) {
  float64x2 r;
  if (!std::isfinite(x._limbs[0])) {
    if (x._limbs[0] > 0.0) r._limbs[0] = x._limbs[0]; // +inf
    else if (x._limbs[0] < 0.0) r._limbs[0] = 0.0;    // -inf → 0
    else r._limbs[0] = x._limbs[0];                    // NaN
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] < exp2_min) return float64x2();
  if (x._limbs[0] > exp2_max) {
    r._limbs[0] = std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
    return r;
  }
  return exp2_kernel(x);
}

float64x2 exp_full(float64x2 const &x) {
  return exp2_full(x * float64x2(log2_e_hi, log2_e_lo));
}

float64x2 log2_kernel(float64x2 const &x) {
  double hi = x._limbs[0];
  if (hi > 15.0 / 16.0 && hi < 17.0 / 16.0) {
    // direct path: t = (x - 1) / (x + 1)
    float64x2 num = x - float64x2(1.0);
    float64x2 den = x + float64x2(1.0);
    float64x2 t = num / den;
    float64x2 t_sq = t * t;
    float64x2 p = horner(t_sq, log2_wide_hi, log2_wide_lo, 9);
    return t * p;
  }
  // table path
  int e = std::ilogb(hi); // matches Julia/IEEE convention (mantissa in [1,2))
  float64x2 m;
  m._limbs[0] = std::ldexp(x._limbs[0], -e);
  m._limbs[1] = std::ldexp(x._limbs[1], -e);
  int idx = static_cast<int>((m._limbs[0] - 1.0) * 32.0);
  if (idx < 0) idx = 0;
  if (idx > 31) idx = 31;
  float64x2 c = float64x2(log2_centers[idx], 0.0);
  float64x2 v = float64x2(log2_values_hi[idx], log2_values_lo[idx]);
  float64x2 num = m - c;
  float64x2 den = m + c;
  float64x2 t = num / den;
  float64x2 t_sq = t * t;
  float64x2 p = horner(t_sq, log2_narrow_hi, log2_narrow_lo, 7);
  return float64x2(static_cast<double>(e), 0.0) + v + t * p;
}

float64x2 log2_full(float64x2 const &x) {
  float64x2 r;
  if (std::isnan(x._limbs[0]) || x._limbs[0] < 0.0) {
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] == 0.0) {
    r._limbs[0] = -std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
    return r;
  }
  if (!std::isfinite(x._limbs[0])) {
    r._limbs[0] = std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
    return r;
  }
  return log2_kernel(x);
}

float64x2 log_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::log(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  // log2(0) = -inf. DD multiply (-inf, 0) * (ln_2_hi, ln_2_lo) would refine
  // into (-inf, NaN) because the two_prod FMA hits -inf + inf; short-circuit.
  float64x2 l2 = log2_full(x);
  if (!std::isfinite(l2._limbs[0])) return l2;
  return l2 * float64x2(ln_2_hi, ln_2_lo);
}

float64x2 log10_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::log10(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  float64x2 l2 = log2_full(x);
  if (!std::isfinite(l2._limbs[0])) return l2;
  return l2 * float64x2(log10_2_hi, log10_2_lo);
}

// ---- expm1 / log1p ---------------------------------------------------------
// Both functions have a narrow cancellation regime near their input zero
// (exp(x)-1 for small x, log(1+x) for small x). Outside that regime the
// naive composition is already DD-accurate, so we gate on |x| and fall
// through to exp_full / log_full there. Inside, a dedicated Taylor /
// atanh-narrow polynomial preserves every bit.

float64x2 expm1_full(float64x2 const &x) {
  // Non-finite: expm1(+inf) = +inf, expm1(-inf) = -1, expm1(NaN) = NaN.
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    if (x._limbs[0] > 0.0)      r._limbs[0] = x._limbs[0];       // +inf
    else if (x._limbs[0] < 0.0) r._limbs[0] = -1.0;               // -inf
    else                        r._limbs[0] = x._limbs[0];       // NaN
    r._limbs[1] = 0.0;
    return r;
  }
  // For |x| >= 0.5, exp(x) moves far enough from 1 that the subtraction
  // retains ~DD precision; defer to the existing kernel (which also
  // handles overflow / subnormal tails).
  if (std::abs(x._limbs[0]) >= 0.5) {
    return exp_full(x) - float64x2(1.0);
  }
  // Small |x|: direct Taylor. expm1(x) = x * Σ x^k / (k+1)!.
  float64x2 p = horner(x, expm1_taylor_hi, expm1_taylor_lo, 25);
  return x * p;
}

float64x2 log1p_full(float64x2 const &x) {
  // Non-finite: log1p(+inf) = +inf, log1p(NaN) = NaN. log1p(-inf) reaches
  // the domain check below via (1 + -inf = -inf).
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::log1p(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  // Domain: 1 + x >= 0 with equality only for x == -1 (→ -inf).
  float64x2 arg = float64x2(1.0) + x;
  if (arg._limbs[0] <= 0.0) {
    float64x2 r;
    r._limbs[1] = 0.0;
    if (arg._limbs[0] == 0.0 && arg._limbs[1] == 0.0) {
      r._limbs[0] = -std::numeric_limits<double>::infinity();
    } else if (arg._limbs[0] < 0.0 ||
               (arg._limbs[0] == 0.0 && arg._limbs[1] < 0.0)) {
      r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    } else {
      // arg is a positive subnormal with hi == 0: fall through.
      return log_full(arg);
    }
    return r;
  }
  // For |x| >= 0.25, log(1+x) is well-conditioned; reuse log_full.
  if (std::abs(x._limbs[0]) >= 0.25) {
    return log_full(arg);
  }
  // Small |x|: log1p(x) = 2·atanh(t) with t = x/(2+x). Baking the 2
  // into the coefficient table (log1p_taylor[k] = 2/(2k+1)) gives the
  // full value in a single t-multiply.
  float64x2 t = x / (float64x2(2.0) + x);
  float64x2 t2 = t * t;
  float64x2 p = horner(t2, log1p_taylor_hi, log1p_taylor_lo, 18);
  return t * p;
}

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
  return horner(x2, sin_taylor_hi, sin_taylor_lo, 13) * x;
}
float64x2 cos_kernel(float64x2 const &x) {
  float64x2 x2 = x * x;
  return horner(x2, cos_taylor_hi, cos_taylor_lo, 13);
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

// ---- sinh / cosh / tanh ----------------------------------------------------
float64x2 sinh_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::sinh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) < 0.1) {
    // 9-term Taylor: x * (1 + x²/6 + x⁴/120 + ...)
    return x * horner(x * x, sinh_taylor_hi, sinh_taylor_lo, 9);
  }
  float64x2 e = exp_full(x);
  float64x2 ei = exp_full(-x);
  return (e - ei) * float64x2(0.5, 0.0);
}

float64x2 cosh_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::cosh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  float64x2 e = exp_full(x);
  float64x2 ei = exp_full(-x);
  return (e + ei) * float64x2(0.5, 0.0);
}

// Fused sinh/cosh. Shares the two exp evaluations (exp(x) and exp(-x)) that
// sinh_full and cosh_full would do independently, halving the exp cost
// for call sites that need both (e.g. complex sin/cos/sinh/cosh).
void sinhcosh_full(float64x2 const &x, float64x2 &s, float64x2 &c) {
  if (!std::isfinite(x._limbs[0])) {
    s._limbs[0] = std::sinh(x._limbs[0]); s._limbs[1] = 0.0;
    c._limbs[0] = std::cosh(x._limbs[0]); c._limbs[1] = 0.0;
    return;
  }
  if (std::abs(x._limbs[0]) < 0.1) {
    // Small-|x| Taylor: sinh via shared coefficients, cosh via 1 + x²·Q(x²).
    // The cosh Taylor here mirrors the structure of sinh_full's path.
    float64x2 x2 = x * x;
    s = x * horner(x2, sinh_taylor_hi, sinh_taylor_lo, 9);
    // cosh(x) = 1 + x²/2 + x⁴/24 + … — reuse separate kernel to keep code small.
    c = cosh_full(x);
    return;
  }
  float64x2 e  = exp_full(x);
  float64x2 ei = exp_full(-x);
  float64x2 half = float64x2(0.5, 0.0);
  s = (e - ei) * half;
  c = (e + ei) * half;
}

float64x2 tanh_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::tanh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  // For |x| > 37, exp(-2|x|) < 2^-106, so tanh(x) = ±1 to full DD.
  if (x._limbs[0] > 37.0) return float64x2(1.0);
  if (x._limbs[0] < -37.0) return float64x2(-1.0);
  if (std::abs(x._limbs[0]) < 0.5) {
    return sinh_full(x) / cosh_full(x);
  }
  // tanh(|x|) = (1 - em2) / (1 + em2) where em2 = exp(-2|x|) is small
  // for |x| > 0.5 — avoids the 1 − small cancellation of 1 − 2/(e^(2|x|)+1).
  bool sign = x._limbs[0] < 0.0;
  float64x2 neg_two_ax;
  neg_two_ax._limbs[0] = sign ? (2.0 * x._limbs[0]) : (-2.0 * x._limbs[0]);
  neg_two_ax._limbs[1] = sign ? (2.0 * x._limbs[1]) : (-2.0 * x._limbs[1]);
  float64x2 em2 = exp_full(neg_two_ax);
  float64x2 res = (float64x2(1.0) - em2) / (float64x2(1.0) + em2);
  return sign ? -res : res;
}

// ---- pow -------------------------------------------------------------------
float64x2 pow_full(float64x2 const &x, float64x2 const &y) {
  if (x._limbs[0] == 0.0 && y._limbs[0] == 0.0) return float64x2(1.0);
  return exp_full(y * log_full(x));
}

// ---- atan (table lookup + rational polynomial, from libquadmath atanq.c) ---
// arctan(t) = t + t^3 * P(t^2) / Q(t^2),  |t| <= 0.09375
// Argument reduced via 87-entry table: arctan(x) = atantbl[k] + arctan(t).
// Table branch (|x| < 32/3): k = floor(8|x|+0.25), t = (x - k/8)/(1 + x*k/8);
//   worst-case |t| = 1/(16*(1 + (k/8)^2)) stays well below 0.09375.
// Large branch (|x| >= 32/3): k = 86 (pi/2), t = -1/|x|; cutover placed at
//   32/3 so worst-case |t| = 3/32 = 0.09375 hits the fitted edge exactly.
float64x2 atan_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::atan(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  double ax_d = std::abs(x._limbs[0]);
  // atan(x) = x - x³/3 + …; for |x| < sqrt(3)·2⁻⁵³ ≈ 1.85e-16 the
  // cubic term is below 0.5 DD ulp, so identity is correct to DD.
  if (ax_d < 1.85e-16) return x;

  bool neg = x._limbs[0] < 0.0;
  float64x2 ax = neg ? -x : x;

  int k;
  float64x2 t;
  if (ax_d >= 32.0 / 3.0) {
    k = 86;
    t = float64x2(-1.0) / ax;
  } else {
    k = static_cast<int>(8.0 * ax_d + 0.25);
    double u_d = 0.125 * k;
    float64x2 u = float64x2(u_d, 0.0);  // exact in double
    t = (ax - u) / (float64x2(1.0) + ax * u);
  }

  // arctan(t) = t + t^3 * P(t^2) / Q(t^2)
  float64x2 t2 = t * t;
  float64x2 p = neval(t2, atan_P_hi, atan_P_lo, 4);
  float64x2 q = deval(t2, atan_Q_hi, atan_Q_lo, 4);
  float64x2 res = float64x2(atan_table_hi[k], atan_table_lo[k])
           + t + t * t2 * p / q;

  return neg ? -res : res;
}

// ---- asin (piecewise rational, from libquadmath asinq.c) -------------------
// Region 1: |x| < 0.5      — asin(x) = x + x*x^2*P(x^2)/Q(x^2)
// Region 2: 0.5 <= |x| < 0.625 — centered at 0.5625
// Region 3: 0.625 <= |x| < 1   — half-angle: asin(x) = pi/2 - 2*asin(sqrt((1-x)/2))
float64x2 asin_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::asin(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  double ax_d = std::abs(x._limbs[0]);
  if (ax_d > 1.0) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  // asin(x) = x + x³/6 + …; for |x| < sqrt(3)·2⁻⁵³ ≈ 1.85e-16 the
  // cubic term is below 0.5 DD ulp, so identity is correct to DD.
  if (ax_d < 1.85e-16) return x;

  bool neg = x._limbs[0] < 0.0;
  float64x2 ax = neg ? -x : x;
  float64x2 res;

  if (ax_d < 0.5) {
    // Region 1: asin(x) = x + x*x^2*P(x^2)/Q(x^2)
    float64x2 t = ax * ax;
    float64x2 p = neval(t, asin_pS_hi, asin_pS_lo, 9);
    float64x2 q = deval(t, asin_qS_hi, asin_qS_lo, 8);
    res = ax + ax * t * p / q;
  } else if (ax_d < 0.625) {
    // Region 2: centered at 0.5625
    float64x2 t = ax - float64x2(0.5625);
    float64x2 p = t * neval(t, asin_rS_hi, asin_rS_lo, 10);
    float64x2 q = deval(t, asin_sS_hi, asin_sS_lo, 9);
    res = float64x2(asinr5625_hi, asinr5625_lo) + p / q;
  } else {
    // Region 3: half-angle identity
    float64x2 w = float64x2(1.0) - ax;
    float64x2 t = w * float64x2(0.5);
    float64x2 s = sqrt(t);
    float64x2 p = neval(t, asin_pS_hi, asin_pS_lo, 9);
    float64x2 q = deval(t, asin_qS_hi, asin_qS_lo, 8);
    float64x2 w2 = t * p / q;
    res = float64x2(pi_half_cw1, pi_half_cw2) - float64x2(2.0) * (s + s * w2);
  }

  return neg ? -res : res;
}

// ---- acos (derived from asin polynomial) ------------------------------------
float64x2 acos_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::acos(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  double ax_d = std::abs(x._limbs[0]);
  if (ax_d > 1.0) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }

  float64x2 pi_half = float64x2(pi_half_cw1, pi_half_cw2);

  if (ax_d <= 0.5) {
    // acos(x) = pi/2 - asin(x) — no cancellation since result ∈ [π/3, 2π/3]
    return pi_half - asin_full(x);
  } else if (x._limbs[0] > 0.0) {
    // x > 0.5: acos(x) = 2*asin(sqrt((1-x)/2))
    float64x2 t = (float64x2(1.0) - x) * float64x2(0.5);
    float64x2 s = sqrt(t);
    float64x2 p = neval(t, asin_pS_hi, asin_pS_lo, 9);
    float64x2 q = deval(t, asin_qS_hi, asin_qS_lo, 8);
    return float64x2(2.0) * (s + s * t * p / q);
  } else {
    // x < -0.5: acos(x) = pi - 2*asin(sqrt((1+x)/2))
    float64x2 t = (float64x2(1.0) + x) * float64x2(0.5);
    float64x2 s = sqrt(t);
    float64x2 p = neval(t, asin_pS_hi, asin_pS_lo, 9);
    float64x2 q = deval(t, asin_qS_hi, asin_qS_lo, 8);
    return float64x2(pi_dd_hi, pi_dd_lo) - float64x2(2.0) * (s + s * t * p / q);
  }
}

// ---- atan2 ------------------------------------------------------------------
float64x2 atan2_full(float64x2 const &y, float64x2 const &x) {
  if (!std::isfinite(x._limbs[0]) || !std::isfinite(y._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::atan2(y._limbs[0], x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] == 0.0 && y._limbs[0] == 0.0) {
    float64x2 r;
    r._limbs[0] = std::atan2(y._limbs[0], x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  float64x2 pi_dd = float64x2(pi_dd_hi, pi_dd_lo);
  float64x2 pi_half_dd = float64x2(pi_half_cw1, pi_half_cw2);
  float64x2 res;
  if (std::abs(x._limbs[0]) >= std::abs(y._limbs[0])) {
    res = atan_full(y / x);
    if (x._limbs[0] < 0.0) {
      if (y._limbs[0] >= 0.0) res = res + pi_dd;
      else res = res - pi_dd;
    }
  } else {
    res = atan_full(x / y);
    if (y._limbs[0] > 0.0) res = pi_half_dd - res;
    else res = -pi_half_dd - res;
  }
  return res;
}

// ---- inverse π-scaled trig -------------------------------------------------
// fn_pi(x) = fn(x) / π via the natural path. The DD multiply by inv_pi
// costs one DD mul (~3 FLOPs) after the main eval.
float64x2 asinpi_full(float64x2 const &x) {
  return asin_full(x) * float64x2(inv_pi_hi, inv_pi_lo);
}
float64x2 acospi_full(float64x2 const &x) {
  return acos_full(x) * float64x2(inv_pi_hi, inv_pi_lo);
}
// atanpi: at ±inf, atan_full bounces through std::atan and returns π/2 to
// only dp precision — then inv_pi multiplication can't recover the missing
// low limb. Short-circuit so atanpi(±inf) = ±0.5 exactly.
float64x2 atanpi_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    if (std::isnan(x._limbs[0])) {
      r._limbs[0] = x._limbs[0];
    } else {
      r._limbs[0] = std::signbit(x._limbs[0]) ? -0.5 : 0.5;
    }
    r._limbs[1] = 0.0;
    return r;
  }
  return atan_full(x) * float64x2(inv_pi_hi, inv_pi_lo);
}
float64x2 atan2pi_full(float64x2 const &y, float64x2 const &x) {
  return atan2_full(y, x) * float64x2(inv_pi_hi, inv_pi_lo);
}

// ---- asinh / acosh / atanh -------------------------------------------------
float64x2 asinh_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::asinh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) < 0.01) {
    return x * horner(x * x, asinh_taylor_hi, asinh_taylor_lo, 15);
  }
  bool sign = x._limbs[0] < 0.0;
  float64x2 ax = sign ? -x : x;
  if (ax._limbs[0] > 1e150) {
    // log(2|x|) asymptotic to avoid x² overflow
    float64x2 r = log_full(ax) + float64x2(ln_2_hi, ln_2_lo);
    return sign ? -r : r;
  }
  float64x2 root = sqrt(ax * ax + float64x2(1.0));
  float64x2 res = log_full(ax + root);
  return sign ? -res : res;
}

float64x2 acosh_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::acosh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] < 1.0) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] > 1e150) {
    return log_full(x) + float64x2(ln_2_hi, ln_2_lo);
  }
  float64x2 root = sqrt(x * x - float64x2(1.0));
  return log_full(x + root);
}

float64x2 atanh_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::atanh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) > 1.0) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) < 0.01) {
    return x * horner(x * x, atanh_taylor_hi, atanh_taylor_lo, 15);
  }
  float64x2 num = float64x2(1.0) + x;
  float64x2 den = float64x2(1.0) - x;
  return float64x2(0.5, 0.0) * log_full(num / den);
}

// ---- erf / erfc (piecewise rational, ported from libquadmath erfq.c) -------
float64x2 erf_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::erf(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  bool neg = x._limbs[0] < 0.0;
  float64x2 ax = neg ? -x : x;

  if (ax._limbs[0] >= 1.0) {
    if (ax._limbs[0] >= 16.0) return neg ? float64x2(-1.0) : float64x2(1.0);
    float64x2 res = float64x2(1.0) - erfc_full(ax);
    return neg ? -res : res;
  }

  float64x2 y;
  if (ax._limbs[0] < 0.875) {
    if (ax._limbs[0] < 1e-18)
      return x + float64x2(erf_efx_hi, erf_efx_lo) * x;
    float64x2 z = ax * ax;
    y = ax + ax * neval(z, erf_TN1_hi, erf_TN1_lo, 8) /
                  deval(z, erf_TD1_hi, erf_TD1_lo, 8);
  } else {
    float64x2 a = ax - float64x2(1.0);
    y = float64x2(erf_const_hi, erf_const_lo) +
        neval(a, erf_TN2_hi, erf_TN2_lo, 8) /
        deval(a, erf_TD2_hi, erf_TD2_lo, 8);
  }
  return neg ? -y : y;
}

float64x2 erfc_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::erfc(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  bool neg = x._limbs[0] < 0.0;
  float64x2 ax = neg ? -x : x;

  if (ax._limbs[0] < 0.25) {
    return float64x2(1.0) - erf_full(x);
  }

  if (ax._limbs[0] > 107.0) {
    return neg ? float64x2(2.0) : float64x2();
  }

  float64x2 res;
  if (ax._limbs[0] < 1.25) {
    // 8 sub-intervals of width 0.125 in [0.25, 1.25)
    int i = static_cast<int>(8.0 * ax._limbs[0]);
    float64x2 z = ax - float64x2(erfc_x0[i - 2]);
    float64x2 y;
    switch (i) {
    case 2: y = neval(z, erfc_RN13_hi, erfc_RN13_lo, 8) / deval(z, erfc_RD13_hi, erfc_RD13_lo, 7); break;
    case 3: y = neval(z, erfc_RN14_hi, erfc_RN14_lo, 8) / deval(z, erfc_RD14_hi, erfc_RD14_lo, 7); break;
    case 4: y = neval(z, erfc_RN15_hi, erfc_RN15_lo, 8) / deval(z, erfc_RD15_hi, erfc_RD15_lo, 7); break;
    case 5: y = neval(z, erfc_RN16_hi, erfc_RN16_lo, 8) / deval(z, erfc_RD16_hi, erfc_RD16_lo, 7); break;
    case 6: y = neval(z, erfc_RN17_hi, erfc_RN17_lo, 8) / deval(z, erfc_RD17_hi, erfc_RD17_lo, 7); break;
    case 7: y = neval(z, erfc_RN18_hi, erfc_RN18_lo, 8) / deval(z, erfc_RD18_hi, erfc_RD18_lo, 7); break;
    case 8: y = neval(z, erfc_RN19_hi, erfc_RN19_lo, 8) / deval(z, erfc_RD19_hi, erfc_RD19_lo, 7); break;
    default: y = neval(z, erfc_RN20_hi, erfc_RN20_lo, 8) / deval(z, erfc_RD20_hi, erfc_RD20_lo, 7); break;
    }
    int ci = i - 2;
    res = y * z + float64x2(erfc_Cb_hi[ci], erfc_Cb_lo[ci]) +
          float64x2(erfc_Ca_hi[ci], erfc_Ca_lo[ci]);
  } else {
    // Asymptotic: erfc(x) = (1/x)*exp(-x^2 - 0.5625 + R(1/x^2))
    if (neg && ax._limbs[0] >= 9.0) return float64x2(2.0);

    float64x2 x2 = ax * ax;
    float64x2 z = float64x2(1.0) / x2;
    int i = static_cast<int>(8.0 / ax._limbs[0]);
    float64x2 p;
    switch (i) {
    case 0: p = neval(z, erfc_AN1_hi, erfc_AN1_lo, 9) / deval(z, erfc_AD1_hi, erfc_AD1_lo, 8); break;
    case 1: p = neval(z, erfc_AN2_hi, erfc_AN2_lo, 11) / deval(z, erfc_AD2_hi, erfc_AD2_lo, 10); break;
    case 2: p = neval(z, erfc_AN3_hi, erfc_AN3_lo, 11) / deval(z, erfc_AD3_hi, erfc_AD3_lo, 10); break;
    case 3: p = neval(z, erfc_AN4_hi, erfc_AN4_lo, 10) / deval(z, erfc_AD4_hi, erfc_AD4_lo, 10); break;
    case 4: p = neval(z, erfc_AN5_hi, erfc_AN5_lo, 10) / deval(z, erfc_AD5_hi, erfc_AD5_lo, 9); break;
    case 5: p = neval(z, erfc_AN6_hi, erfc_AN6_lo, 9) / deval(z, erfc_AD6_hi, erfc_AD6_lo, 9); break;
    case 6: p = neval(z, erfc_AN7_hi, erfc_AN7_lo, 9) / deval(z, erfc_AD7_hi, erfc_AD7_lo, 9); break;
    default: p = neval(z, erfc_AN8_hi, erfc_AN8_lo, 9) / deval(z, erfc_AD8_hi, erfc_AD8_lo, 8); break;
    }

    // Split x for accurate exp(-x^2): truncate hi limb by zeroing the
    // low 35 mantissa bits (keeps ~17 high bits), so s^2 fits exactly
    // in a double (2·18 = 36 bits ≤ 53). std::memcpy sidesteps strict-
    // aliasing concerns; the compiler folds it to a register move.
    uint64_t u;
    std::memcpy(&u, &ax._limbs[0], sizeof(u));
    u &= 0xFFFFFFF800000000ULL;
    double s_hi;
    std::memcpy(&s_hi, &u, sizeof(s_hi));
    float64x2 s = float64x2(s_hi, 0.0);

    float64x2 e1 = exp_full(-(s * s) - float64x2(0.5625));
    float64x2 diff_sq = (s - ax) * (s + ax);
    float64x2 e2 = exp_full(diff_sq + p);
    res = (e1 * e2) / ax;
  }

  return neg ? float64x2(2.0) - res : res;
}

// erfc_scaled(x) = exp(x²) · erfc(x).
// For |x| >= 1.25, substitute erfc(x) = (1/x)·exp(-x² - 0.5625 + R(1/x²))
// so the x² terms cancel: erfc_scaled(|x|) = (1/|x|) · exp(R - 0.5625).
// This avoids the exp(x²)·(tiny) overflow/underflow pair that the direct
// product would hit for large |x|. For the negative branch,
//   erfc_scaled(-|x|) = 2·exp(x²) - erfc_scaled(|x|)
// and exp(x²) itself can overflow; guard that.
// For |x| < 1.25, both factors are moderate and we compute the product
// directly.
float64x2 erfc_scaled_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    if (x._limbs[0] > 0.0)      r._limbs[0] = 0.0;                                            // erfc_scaled(+inf) = 0
    else if (x._limbs[0] < 0.0) r._limbs[0] = std::numeric_limits<double>::infinity();        // erfc_scaled(-inf) = +inf
    else                         r._limbs[0] = x._limbs[0];                                    // NaN
    r._limbs[1] = 0.0;
    return r;
  }
  bool neg = x._limbs[0] < 0.0;
  float64x2 ax = neg ? -x : x;
  if (ax._limbs[0] < 1.25) {
    return exp_full(x * x) * erfc_full(x);
  }
  // Asymptotic branch.
  float64x2 x2 = ax * ax;
  float64x2 z = float64x2(1.0) / x2;
  int i = static_cast<int>(8.0 / ax._limbs[0]);
  if (i > 7) i = 7;
  float64x2 p;
  switch (i) {
  case 0: p = neval(z, erfc_AN1_hi, erfc_AN1_lo, 9)  / deval(z, erfc_AD1_hi, erfc_AD1_lo, 8);  break;
  case 1: p = neval(z, erfc_AN2_hi, erfc_AN2_lo, 11) / deval(z, erfc_AD2_hi, erfc_AD2_lo, 10); break;
  case 2: p = neval(z, erfc_AN3_hi, erfc_AN3_lo, 11) / deval(z, erfc_AD3_hi, erfc_AD3_lo, 10); break;
  case 3: p = neval(z, erfc_AN4_hi, erfc_AN4_lo, 10) / deval(z, erfc_AD4_hi, erfc_AD4_lo, 10); break;
  case 4: p = neval(z, erfc_AN5_hi, erfc_AN5_lo, 10) / deval(z, erfc_AD5_hi, erfc_AD5_lo, 9);  break;
  case 5: p = neval(z, erfc_AN6_hi, erfc_AN6_lo, 9)  / deval(z, erfc_AD6_hi, erfc_AD6_lo, 9);  break;
  case 6: p = neval(z, erfc_AN7_hi, erfc_AN7_lo, 9)  / deval(z, erfc_AD7_hi, erfc_AD7_lo, 9);  break;
  default:p = neval(z, erfc_AN8_hi, erfc_AN8_lo, 9)  / deval(z, erfc_AD8_hi, erfc_AD8_lo, 8);  break;
  }
  float64x2 res = exp_full(p - float64x2(0.5625)) / ax;
  if (neg) {
    float64x2 exp_x2 = exp_full(x2);
    if (!std::isfinite(exp_x2._limbs[0])) {
      float64x2 r;
      r._limbs[0] = std::numeric_limits<double>::infinity();
      r._limbs[1] = 0.0;
      return r;
    }
    res = float64x2(2.0) * exp_x2 - res;
  }
  return res;
}

// ---- tgamma / lgamma -------------------------------------------------------
// Native DD Stirling asymptotic with shift recurrence and reflection.
// Every operation runs in double-double so no libm floor is imposed on
// the result:
//
//   log Γ(x) = (x - 1/2)·log(x) - x + (1/2)·log(2π)
//              + Σ_{k=1..13} c_k / x^{2k-1}
// where c_k = B_{2k} / (2k·(2k-1)). With a shift target of x ≥ 25 the
// k=13 remainder is ~c_13/x^25 ≈ 2e-32 — comfortably below full DD.
//
// Small-x paths: reflection for x < 1/2, upward shift recurrence
// otherwise, accumulating the product x·(x+1)·...·(x+N-1) in DD and
// subtracting a single log_full(prod).

// Stirling/gamma constants are in dd_constants.hh (via #include above).

// Evaluate the Stirling series for log Γ(x). Assumes x is well into the
// asymptotic range (x ≥ ~20) — callers shift first.
float64x2 lgamma_stirling(float64x2 const &x) {
  float64x2 inv_x = float64x2(1.0) / x;
  float64x2 inv_x2 = inv_x * inv_x;
  float64x2 poly =
      horner(inv_x2, stirling_coefs_hi, stirling_coefs_lo, 13);
  float64x2 corr = poly * inv_x;
  float64x2 logx = log_full(x);
  float64x2 xm_half = x - float64x2(0.5);
  return xm_half * logx - x +
         float64x2(half_log_2pi_hi, half_log_2pi_lo) + corr;
}

// Exact DD subtraction: z = x - c, where c is exact in double precision.
static inline float64x2 sub_exact(float64x2 const &x, double c) {
  float64x2 z;
  z._limbs[0] = x._limbs[0] - c;
  z._limbs[1] = x._limbs[1] + ((x._limbs[0] - z._limbs[0]) - c);
  return z;
}

// Piecewise rational approximation for lgamma(x), 0.5 <= x <= 13.5.
// Mirrors lgammaq.c structure: P(z)/Q(z) with z = x - center.
float64x2 lgamma_rational(float64x2 const &x) {
  int nn = (int)std::nearbyint(x._limbs[0]);
  float64x2 z, p, q, c;
  switch (nn) {
  case 0: // x ≈ 0.5: lgamma(x+1) via near-minimum, then - log(x)
    z = x + float64x2(1.0) - float64x2(lgam_x0_hi, lgam_x0_lo);
    p = neval(z, lgam_RN1r5_hi, lgam_RN1r5_lo, 8);
    q = deval(z, lgam_RD1r5_hi, lgam_RD1r5_lo, 8);
    return z * z * (p / q) + float64x2(lgam_y0_hi, lgam_y0_lo) - log_full(x);
  case 1:
    if (x._limbs[0] < 0.875) {
      if (x._limbs[0] <= 0.625) {
        z = x + float64x2(1.0) - float64x2(lgam_x0_hi, lgam_x0_lo);
        p = neval(z, lgam_RN1r5_hi, lgam_RN1r5_lo, 8);
        q = deval(z, lgam_RD1r5_hi, lgam_RD1r5_lo, 8);
        return z * z * (p / q) + float64x2(lgam_y0_hi, lgam_y0_lo) - log_full(x);
      }
      // (0.625, 0.875): z = x - 0.75, lgamma(x+1) = lgam1r75 + z*P/Q
      z = sub_exact(x, 0.75);
      p = neval(z, lgam_RN1r75_hi, lgam_RN1r75_lo, 8);
      q = deval(z, lgam_RD1r75_hi, lgam_RD1r75_lo, 8);
      return z * (p / q) + float64x2(lgam1r75_hi, lgam1r75_lo) - log_full(x);
    }
    if (x._limbs[0] < 1.0) {
      z = x - float64x2(1.0);
      p = neval(z, lgam_RNr9_hi, lgam_RNr9_lo, 8);
      q = deval(z, lgam_RDr9_hi, lgam_RDr9_lo, 8);
      return z * (p / q);
    }
    if (x._limbs[0] <= 1.125) {
      z = x - float64x2(1.0);
      p = neval(z, lgam_RN1_hi, lgam_RN1_lo, 8);
      q = deval(z, lgam_RD1_hi, lgam_RD1_lo, 7);
      return z * (p / q);
    }
    if (x._limbs[0] <= 1.375) {
      z = sub_exact(x, 1.25);
      p = neval(z, lgam_RN1r25_hi, lgam_RN1r25_lo, 9);
      q = deval(z, lgam_RD1r25_hi, lgam_RD1r25_lo, 8);
      return z * (p / q) + float64x2(lgam1r25_hi, lgam1r25_lo);
    }
    // [1.375, 1.5]: near minimum x0
    z = x - float64x2(lgam_x0_hi, lgam_x0_lo);
    p = neval(z, lgam_RN1r5_hi, lgam_RN1r5_lo, 8);
    q = deval(z, lgam_RD1r5_hi, lgam_RD1r5_lo, 8);
    return z * z * (p / q) + float64x2(lgam_y0_hi, lgam_y0_lo);
  case 2:
    if (x._limbs[0] < 1.625) {
      z = x - float64x2(lgam_x0_hi, lgam_x0_lo);
      p = neval(z, lgam_RN1r5_hi, lgam_RN1r5_lo, 8);
      q = deval(z, lgam_RD1r5_hi, lgam_RD1r5_lo, 8);
      return z * z * (p / q) + float64x2(lgam_y0_hi, lgam_y0_lo);
    }
    if (x._limbs[0] < 1.875) {
      z = sub_exact(x, 1.75);
      p = neval(z, lgam_RN1r75_hi, lgam_RN1r75_lo, 8);
      q = deval(z, lgam_RD1r75_hi, lgam_RD1r75_lo, 8);
      return z * (p / q) + float64x2(lgam1r75_hi, lgam1r75_lo);
    }
    if (x._limbs[0] < 2.375) {
      z = sub_exact(x, 2.0);
      p = neval(z, lgam_RN2_hi, lgam_RN2_lo, 9);
      q = deval(z, lgam_RD2_hi, lgam_RD2_lo, 9);
      return z * (p / q);
    }
    z = sub_exact(x, 2.5);
    p = neval(z, lgam_RN2r5_hi, lgam_RN2r5_lo, 8);
    q = deval(z, lgam_RD2r5_hi, lgam_RD2r5_lo, 8);
    return z * (p / q) + float64x2(lgam2r5_hi, lgam2r5_lo);
  case 3:
    if (x._limbs[0] < 2.75) {
      z = sub_exact(x, 2.5);
      p = neval(z, lgam_RN2r5_hi, lgam_RN2r5_lo, 8);
      q = deval(z, lgam_RD2r5_hi, lgam_RD2r5_lo, 8);
      return z * (p / q) + float64x2(lgam2r5_hi, lgam2r5_lo);
    }
    z = sub_exact(x, 3.0);
    p = neval(z, lgam_RN3_hi, lgam_RN3_lo, 9);
    q = deval(z, lgam_RD3_hi, lgam_RD3_lo, 9);
    return z * (p / q) + float64x2(lgam3_hi, lgam3_lo);
  case 4:
    z = sub_exact(x, 4.0);
    p = neval(z, lgam_RN4_hi, lgam_RN4_lo, 9);
    q = deval(z, lgam_RD4_hi, lgam_RD4_lo, 9);
    return z * (p / q) + float64x2(lgam4_hi, lgam4_lo);
  case 5:
    z = sub_exact(x, 5.0);
    p = neval(z, lgam_RN5_hi, lgam_RN5_lo, 9);
    q = deval(z, lgam_RD5_hi, lgam_RD5_lo, 8);
    return z * (p / q) + float64x2(lgam5_hi, lgam5_lo);
  case 6:
    z = sub_exact(x, 6.0);
    p = neval(z, lgam_RN6_hi, lgam_RN6_lo, 8);
    q = deval(z, lgam_RD6_hi, lgam_RD6_lo, 8);
    return z * (p / q) + float64x2(lgam6_hi, lgam6_lo);
  case 7:
    z = sub_exact(x, 7.0);
    p = neval(z, lgam_RN7_hi, lgam_RN7_lo, 8);
    q = deval(z, lgam_RD7_hi, lgam_RD7_lo, 7);
    return z * (p / q) + float64x2(lgam7_hi, lgam7_lo);
  case 8:
    z = sub_exact(x, 8.0);
    p = neval(z, lgam_RN8_hi, lgam_RN8_lo, 8);
    q = deval(z, lgam_RD8_hi, lgam_RD8_lo, 7);
    return z * (p / q) + float64x2(lgam8_hi, lgam8_lo);
  case 9:
    z = sub_exact(x, 9.0);
    p = neval(z, lgam_RN9_hi, lgam_RN9_lo, 7);
    q = deval(z, lgam_RD9_hi, lgam_RD9_lo, 7);
    return z * (p / q) + float64x2(lgam9_hi, lgam9_lo);
  case 10:
    z = sub_exact(x, 10.0);
    p = neval(z, lgam_RN10_hi, lgam_RN10_lo, 7);
    q = deval(z, lgam_RD10_hi, lgam_RD10_lo, 7);
    return z * (p / q) + float64x2(lgam10_hi, lgam10_lo);
  case 11:
    z = sub_exact(x, 11.0);
    p = neval(z, lgam_RN11_hi, lgam_RN11_lo, 7);
    q = deval(z, lgam_RD11_hi, lgam_RD11_lo, 6);
    return z * (p / q) + float64x2(lgam11_hi, lgam11_lo);
  case 12:
    z = sub_exact(x, 12.0);
    p = neval(z, lgam_RN12_hi, lgam_RN12_lo, 7);
    q = deval(z, lgam_RD12_hi, lgam_RD12_lo, 6);
    return z * (p / q) + float64x2(lgam12_hi, lgam12_lo);
  case 13:
    z = sub_exact(x, 13.0);
    p = neval(z, lgam_RN13_hi, lgam_RN13_lo, 7);
    q = deval(z, lgam_RD13_hi, lgam_RD13_lo, 6);
    return z * (p / q) + float64x2(lgam13_hi, lgam13_lo);
  default:
    return lgamma_stirling(x);
  }
}

// Compute lgamma for positive x >= 0.5.
float64x2 lgamma_positive(float64x2 const &x) {
  if (x._limbs[0] >= 13.5)
    return lgamma_stirling(x);
  return lgamma_rational(x);
}

float64x2 lgamma_full(float64x2 const &x) {
  float64x2 r;
  double hi = x._limbs[0];
  if (std::isnan(hi)) {
    r._limbs[0] = hi;
    r._limbs[1] = 0.0;
    return r;
  }
  if (!std::isfinite(hi)) {
    r._limbs[0] = std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
    return r;
  }
  // Non-positive integer pole (only exactly, i.e. lo == 0).
  if (hi <= 0.0 && hi == std::nearbyint(hi) && x._limbs[1] == 0.0) {
    r._limbs[0] = std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
    return r;
  }
  if (hi < 0.5) {
    // Reflection: log|Γ(x)| = log π - log|sin(πx)| - log|Γ(1-x)|
    float64x2 s = multifloats::abs(sinpi_full(x));
    float64x2 one_minus_x = float64x2(1.0) - x;
    float64x2 lgam_1mx = lgamma_positive(one_minus_x);
    return float64x2(log_pi_hi, log_pi_lo) - log_full(s) - lgam_1mx;
  }
  return lgamma_positive(x);
}

float64x2 tgamma_full(float64x2 const &x) {
  float64x2 r;
  double hi = x._limbs[0];
  if (std::isnan(hi)) {
    r._limbs[0] = hi;
    r._limbs[1] = 0.0;
    return r;
  }
  if (!std::isfinite(hi)) {
    if (hi > 0.0) {
      r._limbs[0] = hi;
      r._limbs[1] = 0.0;
      return r;
    }
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  if (hi <= 0.0 && hi == std::nearbyint(hi) && x._limbs[1] == 0.0) {
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  // Overflow: Γ(171.625) ≈ 2^1024. Beyond that the result overflows.
  if (hi > 171.624) {
    r._limbs[0] = std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
    return r;
  }
  if (hi < 0.5) {
    // Reflection: Γ(x) = π / (sin(πx) · Γ(1-x))
    // For very negative x, Γ(1-x) overflows and the result underflows; let
    // the DD division produce the underflowed / signed zero.
    float64x2 s = sinpi_full(x);
    float64x2 one_minus_x = float64x2(1.0) - x;
    // Use lgamma path for the (1-x) branch so we avoid a double tgamma
    // recursion and keep the large factor inside a log.
    float64x2 lg = lgamma_positive(one_minus_x);
    // Γ(1-x) = exp(lg). sin(πx) may be negative → track sign separately.
    bool neg = s._limbs[0] < 0.0;
    if (neg) {
      s._limbs[0] = -s._limbs[0];
      s._limbs[1] = -s._limbs[1];
    }
    float64x2 pi_dd = float64x2(pi_dd_hi, pi_dd_lo);
    float64x2 g1mx = exp_full(lg);
    float64x2 out = pi_dd / (s * g1mx);
    if (neg) {
      out._limbs[0] = -out._limbs[0];
      out._limbs[1] = -out._limbs[1];
    }
    return out;
  }
  return exp_full(lgamma_positive(x));
}


// =============================================================================
// Bessel functions — piecewise rational approximation from libquadmath j0q/j1q
// =============================================================================

// Shared asymptotic P(z),Q(z) selection for order 0 (used by both j0 and y0).
// z = 1/x^2, xinv_d = 1/x (double).
static void bessel_pq0(double xinv_d, float64x2 const &z, float64x2 &p, float64x2 &q) {
  if (xinv_d <= 0.25) {
    if (xinv_d <= 0.125) {
      if (xinv_d <= 0.0625) {
        p = neval(z, j0_P16_IN_hi, j0_P16_IN_lo, 9) / deval(z, j0_P16_ID_hi, j0_P16_ID_lo, 9);
        q = neval(z, j0_Q16_IN_hi, j0_Q16_IN_lo, 10) / deval(z, j0_Q16_ID_hi, j0_Q16_ID_lo, 9);
      } else {
        p = neval(z, j0_P8_16N_hi, j0_P8_16N_lo, 10) / deval(z, j0_P8_16D_hi, j0_P8_16D_lo, 10);
        q = neval(z, j0_Q8_16N_hi, j0_Q8_16N_lo, 11) / deval(z, j0_Q8_16D_hi, j0_Q8_16D_lo, 11);
      }
    } else if (xinv_d <= 0.1875) {
      p = neval(z, j0_P5_8N_hi, j0_P5_8N_lo, 10) / deval(z, j0_P5_8D_hi, j0_P5_8D_lo, 9);
      q = neval(z, j0_Q5_8N_hi, j0_Q5_8N_lo, 10) / deval(z, j0_Q5_8D_hi, j0_Q5_8D_lo, 10);
    } else {
      p = neval(z, j0_P4_5N_hi, j0_P4_5N_lo, 9) / deval(z, j0_P4_5D_hi, j0_P4_5D_lo, 9);
      q = neval(z, j0_Q4_5N_hi, j0_Q4_5N_lo, 10) / deval(z, j0_Q4_5D_hi, j0_Q4_5D_lo, 9);
    }
  } else {
    if (xinv_d <= 0.375) {
      if (xinv_d <= 0.3125) {
        p = neval(z, j0_P3r2_4N_hi, j0_P3r2_4N_lo, 9) / deval(z, j0_P3r2_4D_hi, j0_P3r2_4D_lo, 9);
        q = neval(z, j0_Q3r2_4N_hi, j0_Q3r2_4N_lo, 10) / deval(z, j0_Q3r2_4D_hi, j0_Q3r2_4D_lo, 9);
      } else {
        p = neval(z, j0_P2r7_3r2N_hi, j0_P2r7_3r2N_lo, 9) / deval(z, j0_P2r7_3r2D_hi, j0_P2r7_3r2D_lo, 8);
        q = neval(z, j0_Q2r7_3r2N_hi, j0_Q2r7_3r2N_lo, 9) / deval(z, j0_Q2r7_3r2D_hi, j0_Q2r7_3r2D_lo, 9);
      }
    } else if (xinv_d <= 0.4375) {
      p = neval(z, j0_P2r3_2r7N_hi, j0_P2r3_2r7N_lo, 9) / deval(z, j0_P2r3_2r7D_hi, j0_P2r3_2r7D_lo, 8);
      q = neval(z, j0_Q2r3_2r7N_hi, j0_Q2r3_2r7N_lo, 9) / deval(z, j0_Q2r3_2r7D_hi, j0_Q2r3_2r7D_lo, 8);
    } else {
      p = neval(z, j0_P2_2r3N_hi, j0_P2_2r3N_lo, 8) / deval(z, j0_P2_2r3D_hi, j0_P2_2r3D_lo, 8);
      q = neval(z, j0_Q2_2r3N_hi, j0_Q2_2r3N_lo, 9) / deval(z, j0_Q2_2r3D_hi, j0_Q2_2r3D_lo, 8);
    }
  }
}

// Shared asymptotic P(z),Q(z) selection for order 1 (used by both j1 and y1).
static void bessel_pq1(double xinv_d, float64x2 const &z, float64x2 &p, float64x2 &q) {
  if (xinv_d <= 0.25) {
    if (xinv_d <= 0.125) {
      if (xinv_d <= 0.0625) {
        p = neval(z, j1_P16_IN_hi, j1_P16_IN_lo, 9) / deval(z, j1_P16_ID_hi, j1_P16_ID_lo, 9);
        q = neval(z, j1_Q16_IN_hi, j1_Q16_IN_lo, 10) / deval(z, j1_Q16_ID_hi, j1_Q16_ID_lo, 9);
      } else {
        p = neval(z, j1_P8_16N_hi, j1_P8_16N_lo, 11) / deval(z, j1_P8_16D_hi, j1_P8_16D_lo, 10);
        q = neval(z, j1_Q8_16N_hi, j1_Q8_16N_lo, 11) / deval(z, j1_Q8_16D_hi, j1_Q8_16D_lo, 11);
      }
    } else if (xinv_d <= 0.1875) {
      p = neval(z, j1_P5_8N_hi, j1_P5_8N_lo, 10) / deval(z, j1_P5_8D_hi, j1_P5_8D_lo, 10);
      q = neval(z, j1_Q5_8N_hi, j1_Q5_8N_lo, 10) / deval(z, j1_Q5_8D_hi, j1_Q5_8D_lo, 10);
    } else {
      p = neval(z, j1_P4_5N_hi, j1_P4_5N_lo, 10) / deval(z, j1_P4_5D_hi, j1_P4_5D_lo, 9);
      q = neval(z, j1_Q4_5N_hi, j1_Q4_5N_lo, 10) / deval(z, j1_Q4_5D_hi, j1_Q4_5D_lo, 9);
    }
  } else {
    if (xinv_d <= 0.375) {
      if (xinv_d <= 0.3125) {
        p = neval(z, j1_P3r2_4N_hi, j1_P3r2_4N_lo, 9) / deval(z, j1_P3r2_4D_hi, j1_P3r2_4D_lo, 9);
        q = neval(z, j1_Q3r2_4N_hi, j1_Q3r2_4N_lo, 9) / deval(z, j1_Q3r2_4D_hi, j1_Q3r2_4D_lo, 9);
      } else {
        p = neval(z, j1_P2r7_3r2N_hi, j1_P2r7_3r2N_lo, 9) / deval(z, j1_P2r7_3r2D_hi, j1_P2r7_3r2D_lo, 8);
        q = neval(z, j1_Q2r7_3r2N_hi, j1_Q2r7_3r2N_lo, 9) / deval(z, j1_Q2r7_3r2D_hi, j1_Q2r7_3r2D_lo, 9);
      }
    } else if (xinv_d <= 0.4375) {
      p = neval(z, j1_P2r3_2r7N_hi, j1_P2r3_2r7N_lo, 9) / deval(z, j1_P2r3_2r7D_hi, j1_P2r3_2r7D_lo, 8);
      q = neval(z, j1_Q2r3_2r7N_hi, j1_Q2r3_2r7N_lo, 9) / deval(z, j1_Q2r3_2r7D_hi, j1_Q2r3_2r7D_lo, 8);
    } else {
      p = neval(z, j1_P2_2r3N_hi, j1_P2_2r3N_lo, 8) / deval(z, j1_P2_2r3D_hi, j1_P2_2r3D_lo, 8);
      q = neval(z, j1_Q2_2r3N_hi, j1_Q2_2r3N_lo, 9) / deval(z, j1_Q2_2r3D_hi, j1_Q2_2r3D_lo, 8);
    }
  }
}

float64x2 bessel_j0_full(float64x2 const &x) {
  float64x2 ax = x;
  if (ax._limbs[0] < 0.0) { ax._limbs[0] = -ax._limbs[0]; ax._limbs[1] = -ax._limbs[1]; }
  double xx = ax._limbs[0];
  if (xx == 0.0) return float64x2(1.0);
  if (!std::isfinite(xx)) return float64x2(0.0);

  if (xx <= 2.0) {
    float64x2 z = ax * ax;
    float64x2 r = z * z * neval(z, j0_J0_2N_hi, j0_J0_2N_lo, 6) /
                      deval(z, j0_J0_2D_hi, j0_J0_2D_lo, 6);
    return r - z * float64x2(0.25) + float64x2(1.0);
  }

  float64x2 angle = ax - float64x2(pi_quarter_hi, pi_quarter_lo);
  float64x2 s, c;
  sincos_full(angle, s, c);
  float64x2 xinv = float64x2(1.0) / ax;
  float64x2 z = xinv * xinv;
  float64x2 p, q;
  bessel_pq0(1.0 / xx, z, p, q);
  p = float64x2(1.0) + z * p;
  q = (z * q - float64x2(0.125)) * xinv;
  float64x2 tpi = float64x2(two_over_pi_hi, two_over_pi_lo);
  return multifloats::sqrt(tpi / ax) * (p * c - q * s);
}

float64x2 bessel_j1_full(float64x2 const &x) {
  // j1 is odd, so (+0, -eps) must be detected as negative — checking
  // only the hi limb would miss that case.
  bool neg = x._limbs[0] < 0.0 || (x._limbs[0] == 0.0 && x._limbs[1] < 0.0);
  float64x2 ax = neg ? -x : x;
  double xx = ax._limbs[0];
  if (xx == 0.0) return float64x2(0.0);
  if (!std::isfinite(xx)) return float64x2(0.0);

  float64x2 res;
  if (xx <= 2.0) {
    float64x2 z = ax * ax;
    res = ax * float64x2(0.5) + ax * z * neval(z, j1_J1_2N_hi, j1_J1_2N_lo, 6) /
                                     deval(z, j1_J1_2D_hi, j1_J1_2D_lo, 6);
  } else {
    float64x2 angle = ax - float64x2(three_pi_quarter_hi, three_pi_quarter_lo);
    float64x2 s, c;
    sincos_full(angle, s, c);
    float64x2 xinv = float64x2(1.0) / ax;
    float64x2 z = xinv * xinv;
    float64x2 p, q;
    bessel_pq1(1.0 / xx, z, p, q);
    p = float64x2(1.0) + z * p;
    q = (z * q + float64x2(0.375)) * xinv;
    float64x2 tpi = float64x2(two_over_pi_hi, two_over_pi_lo);
    res = multifloats::sqrt(tpi / ax) * (p * c - q * s);
  }
  if (neg) { res._limbs[0] = -res._limbs[0]; res._limbs[1] = -res._limbs[1]; }
  return res;
}

float64x2 bessel_y0_full(float64x2 const &x) {
  double xx = x._limbs[0];
  if (xx <= 0.0 || !std::isfinite(xx)) {
    float64x2 r;
    r._limbs[0] = (xx == 0.0) ? -std::numeric_limits<double>::infinity()
                               : std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }

  if (xx <= 1e-17) {
    float64x2 tpi = float64x2(two_over_pi_hi, two_over_pi_lo);
    return float64x2(bessel_u0_hi, bessel_u0_lo) + tpi * log_full(x);
  }

  if (xx <= 2.0) {
    float64x2 z = x * x;
    float64x2 p = neval(z, j0_Y0_2N_hi, j0_Y0_2N_lo, 7) /
             deval(z, j0_Y0_2D_hi, j0_Y0_2D_lo, 7);
    float64x2 tpi = float64x2(two_over_pi_hi, two_over_pi_lo);
    return tpi * log_full(x) * bessel_j0_full(x) + p;
  }

  float64x2 angle = x - float64x2(pi_quarter_hi, pi_quarter_lo);
  float64x2 s, c;
  sincos_full(angle, s, c);
  float64x2 xinv = float64x2(1.0) / x;
  float64x2 z = xinv * xinv;
  float64x2 p, q;
  bessel_pq0(1.0 / xx, z, p, q);
  p = float64x2(1.0) + z * p;
  q = (z * q - float64x2(0.125)) * xinv;
  float64x2 tpi = float64x2(two_over_pi_hi, two_over_pi_lo);
  return multifloats::sqrt(tpi / x) * (p * s + q * c);
}

float64x2 bessel_y1_full(float64x2 const &x) {
  double xx = x._limbs[0];
  if (xx <= 0.0 || !std::isfinite(xx)) {
    float64x2 r;
    r._limbs[0] = (xx == 0.0) ? -std::numeric_limits<double>::infinity()
                               : std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }

  if (xx <= 1e-30) {
    float64x2 tpi = float64x2(two_over_pi_hi, two_over_pi_lo);
    return float64x2(0.0) - tpi / x;
  }

  if (xx <= 2.0) {
    float64x2 z = x * x;
    float64x2 p = x * neval(z, j1_Y1_2N_hi, j1_Y1_2N_lo, 7) /
                 deval(z, j1_Y1_2D_hi, j1_Y1_2D_lo, 7);
    float64x2 tpi = float64x2(two_over_pi_hi, two_over_pi_lo);
    return tpi * (log_full(x) * bessel_j1_full(x) - float64x2(1.0) / x) + p;
  }

  float64x2 angle = x - float64x2(three_pi_quarter_hi, three_pi_quarter_lo);
  float64x2 s, c;
  sincos_full(angle, s, c);
  float64x2 xinv = float64x2(1.0) / x;
  float64x2 z = xinv * xinv;
  float64x2 p, q;
  bessel_pq1(1.0 / xx, z, p, q);
  p = float64x2(1.0) + z * p;
  q = (z * q + float64x2(0.375)) * xinv;
  float64x2 tpi = float64x2(two_over_pi_hi, two_over_pi_lo);
  return multifloats::sqrt(tpi / x) * (p * s + q * c);
}

// ---- Integer-order Bessel Jn / Yn (recurrence from j0/j1 / y0/y1) ----------
//
// Jn: forward recurrence  J_{i+1}(x) = (2i/x) J_i(x) − J_{i-1}(x)  is stable
// for i ≤ |x|.  For i > |x| it blows up, so we use Miller's algorithm —
// backward recurrence from a sufficiently high start index (determined by a
// continued-fraction convergence criterion), with periodic rescaling to avoid
// overflow of the unnormalized sequence, then normalized against J_0 or J_1
// (whichever has larger magnitude).  Mirrors Numerical Recipes §6.5.
//
// Yn: forward recurrence is always stable (|Y_n| grows with n for fixed x),
// so no Miller's algorithm is needed.

static float64x2 bessel_jn_full(int n, float64x2 x) {
  // Track sign flips from the identities
  //   J_{-n}(x) = (-1)^n J_n(x)     and     J_n(-x) = (-1)^n J_n(x).
  int parity = 0;
  if (n < 0) { n = -n; parity ^= (n & 1); }
  if (n == 0) return bessel_j0_full(x);
  if (n == 1) {
    float64x2 r = bessel_j1_full(x);
    return parity ? -r : r;
  }
  if (x._limbs[0] < 0.0) { x = -x; parity ^= (n & 1); }

  double xx = x._limbs[0];
  if (xx == 0.0 || !std::isfinite(xx)) return float64x2(0.0);

  float64x2 result;
  if (static_cast<double>(n) <= xx) {
    // Forward recurrence.
    float64x2 a = bessel_j0_full(x);
    float64x2 b = bessel_j1_full(x);
    for (int i = 1; i < n; ++i) {
      float64x2 temp = b;
      b = float64x2(static_cast<double>(2 * i)) / x * b - a;
      a = temp;
    }
    result = b;
  } else {
    // Continued-fraction start: iterate until q1 crosses 1e17. The index k at
    // which that happens is the safe extra depth for Miller's backward sweep.
    float64x2 w = float64x2(static_cast<double>(2 * n) / xx);
    float64x2 h = float64x2(2.0 / xx);
    float64x2 q0 = w;
    float64x2 z = w + h;
    float64x2 q1 = w * z - float64x2(1.0);
    int k = 1;
    while (q1._limbs[0] < 1e17) {
      ++k;
      z = z + h;
      float64x2 tmp = z * q1 - q0;
      q0 = q1;
      q1 = tmp;
    }
    // Downward Lentz from index 2(n+k) to 2n.
    float64x2 t(0.0);
    for (int i = 2 * (n + k); i >= 2 * n; i -= 2) {
      t = float64x2(1.0) / (float64x2(static_cast<double>(i)) / x - t);
    }
    float64x2 a = t;
    float64x2 b(1.0);
    // Decide whether this problem size can overflow during descent.
    float64x2 v = float64x2(2.0 / xx);
    float64x2 log_arg = float64x2(static_cast<double>(n)) *
                   log_full(multifloats::abs(v * float64x2(static_cast<double>(n))));
    bool rescale = log_arg._limbs[0] >= 1.1356e4;
    for (int i = n - 1; i >= 1; --i) {
      float64x2 temp = b;
      b = float64x2(static_cast<double>(2 * i)) * b / x - a;
      a = temp;
      if (rescale && b._limbs[0] > 1e100) {
        a = a / b;
        t = t / b;
        b = float64x2(1.0);
      }
    }
    float64x2 z0 = bessel_j0_full(x);
    float64x2 z1 = bessel_j1_full(x);
    if (std::abs(z0._limbs[0]) >= std::abs(z1._limbs[0])) {
      result = t * z0 / b;
    } else {
      result = t * z1 / a;
    }
  }
  return parity ? -result : result;
}

static float64x2 bessel_yn_full(int n, float64x2 x) {
  // Y_{-n}(x) = (-1)^n Y_n(x). No x-sign identity — Y is undefined on x < 0.
  int parity = 0;
  if (n < 0) { n = -n; parity ^= (n & 1); }
  double xx = x._limbs[0];
  if (xx <= 0.0) {
    float64x2 r;
    r._limbs[0] = (xx == 0.0) ? -std::numeric_limits<double>::infinity()
                              :  std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  if (n == 0) {
    float64x2 r = bessel_y0_full(x);
    return parity ? -r : r;
  }
  if (n == 1) {
    float64x2 r = bessel_y1_full(x);
    return parity ? -r : r;
  }
  if (!std::isfinite(xx)) return float64x2(0.0);
  float64x2 a = bessel_y0_full(x);
  float64x2 b = bessel_y1_full(x);
  for (int i = 1; i < n; ++i) {
    float64x2 temp = b;
    b = float64x2(static_cast<double>(2 * i)) / x * b - a;
    a = temp;
    if (!std::isfinite(b._limbs[0])) break;
  }
  return parity ? -b : b;
}

// Range variant: compute Y_{n1}, Y_{n1+1}, …, Y_{n2}(x) into `out` via a
// single forward recurrence sweep. Cheaper than n2−n1+1 independent calls.
// Caller guarantees out[] has room for n2−n1+1 elements.
static void bessel_yn_range_full(int n1, int n2, float64x2 x, float64x2 *out) {
  int size = n2 - n1 + 1;
  if (size <= 0) return;
  float64x2 a = bessel_y0_full(x);
  float64x2 b = bessel_y1_full(x);
  if (n1 == 0) out[0] = a;
  if (n1 <= 1 && n2 >= 1) out[1 - n1] = b;
  for (int n = 2; n <= n2; ++n) {
    float64x2 temp = b;
    b = float64x2(static_cast<double>(2 * (n - 1))) / x * b - a;
    a = temp;
    if (n >= n1) out[n - n1] = b;
    if (!std::isfinite(b._limbs[0])) {
      for (int j = std::max(n - n1 + 1, 0); j < size; ++j) out[j] = b;
      return;
    }
  }
}

} // anonymous namespace

// =============================================================================
// C-ABI entry points — extern "C" functions following math.h naming convention.
// These are the canonical DD implementations; the C++ template wrappers in
// multifloats.hh and the Fortran bind(C) interfaces both call these.
// =============================================================================

// multifloats_c.h is already pulled in via multifloats.hh.

namespace {
static inline float64x2 from(float64x2_t x) { float64x2 r; r._limbs[0] = x.hi; r._limbs[1] = x.lo; return r; }
static inline float64x2_t to(float64x2 const &x) { return {x._limbs[0], x._limbs[1]}; }

// Inline matmul micro-ops used by the dispatched panel templates below.
// `always_inline` is a hint: the panel templates rely on both ends of
// the MAC fusing into one register-resident body; a stray function call
// here would spill the accumulator ladder and cost ~2× the hot loop.
__attribute__((always_inline))
static inline void mac_inl(double ah, double al, double bh, double bl,
                              double &s_hi, double &s_lo) {
  double p = ah * bh;
  double e = std::fma(ah, bh, -p);
  double cross = std::fma(ah, bl, al * bh);
  double t = s_hi + p;
  double bp = t - s_hi;
  double ap = t - bp;
  double aerr = s_hi - ap;
  double berr = p - bp;
  s_lo += (e + cross) + (aerr + berr);
  s_hi = t;
}

// Renormalize an accumulator pair. Uses full two_sum rather than
// fast_two_sum because |s_hi| >= |s_lo| is not guaranteed after
// cancellation in a long compensated dot-product.
__attribute__((always_inline))
static inline void renorm_inl(double &s_hi, double &s_lo) {
  double t = s_hi + s_lo;
  double bp = t - s_hi;
  double ap = t - bp;
  double aerr = s_hi - ap;
  double berr = s_lo - bp;
  s_lo = aerr + berr;
  s_hi = t;
}

__attribute__((always_inline))
static inline float64x2_t finalize_inl(double s_hi, double s_lo) {
  renorm_inl(s_hi, s_lo);
  return {s_hi, s_lo};
}

// AXPY-style panel µkernel. Processes a row-panel of A of compile-time
// height MR against a contiguous x / b-column, writing MR finalized DD
// outputs. `lda` is the leading dim of A (so the panel can sit in a
// larger matrix); `renorm_interval > 0` triggers a two_sum of each
// accumulator pair every `renorm_interval` reductions (keeps s_lo
// bounded and precision ~DD for large k, matching dot_product).
template <int MR>
static inline void gaxpy_mv_panel(const float64x2_t *__restrict__ a,
                                     const float64x2_t *__restrict__ x,
                                     float64x2_t *__restrict__ y, int64_t lda,
                                     int64_t k, int64_t renorm_interval) {
  double s_hi[MR] = {}, s_lo[MR] = {};
  // Fast path: if no intermediate renorm is needed, drop the chunking
  // scaffolding entirely — the p-loop is a single tight sweep.
  if (renorm_interval <= 0 || k <= renorm_interval) {
    for (int64_t p = 0; p < k; ++p) {
      const double xh = x[p].hi;
      const double xl = x[p].lo;
      const float64x2_t *__restrict__ acol = a + p * lda;
      for (int i = 0; i < MR; ++i) {
        mac_inl(acol[i].hi, acol[i].lo, xh, xl, s_hi[i], s_lo[i]);
      }
    }
    for (int i = 0; i < MR; ++i) y[i] = finalize_inl(s_hi[i], s_lo[i]);
    return;
  }
  // Chunked path: for large k, do `renorm_interval` reductions then
  // two_sum each accumulator pair to keep s_lo bounded.
  const int64_t chunk = renorm_interval;
  int64_t p0 = 0;
  while (p0 < k) {
    int64_t pend = p0 + chunk;
    if (pend > k) pend = k;
    for (int64_t p = p0; p < pend; ++p) {
      const double xh = x[p].hi;
      const double xl = x[p].lo;
      const float64x2_t *__restrict__ acol = a + p * lda;
      for (int i = 0; i < MR; ++i) {
        mac_inl(acol[i].hi, acol[i].lo, xh, xl, s_hi[i], s_lo[i]);
      }
    }
    p0 = pend;
    if (p0 < k) {
      for (int i = 0; i < MR; ++i) renorm_inl(s_hi[i], s_lo[i]);
    }
  }
  for (int i = 0; i < MR; ++i) y[i] = finalize_inl(s_hi[i], s_lo[i]);
}

// Cold tail handler: 1..7 leftover rows. Kept out-of-line so the hot
// dispatcher stays small enough for the compiler to inline panel<8>
// into matmuldd_mv — otherwise the 7 extra template instantiations
// would bloat the caller and push the accumulators out of registers.
__attribute__((noinline))
static void gaxpy_mv_tail(const float64x2_t *__restrict__ a,
                             const float64x2_t *__restrict__ x,
                             float64x2_t *__restrict__ y, int tail,
                             int64_t lda, int64_t k,
                             int64_t renorm_interval) {
  switch (tail) {
    case 1: gaxpy_mv_panel<1>(a, x, y, lda, k, renorm_interval); break;
    case 2: gaxpy_mv_panel<2>(a, x, y, lda, k, renorm_interval); break;
    case 3: gaxpy_mv_panel<3>(a, x, y, lda, k, renorm_interval); break;
    case 4: gaxpy_mv_panel<4>(a, x, y, lda, k, renorm_interval); break;
    case 5: gaxpy_mv_panel<5>(a, x, y, lda, k, renorm_interval); break;
    case 6: gaxpy_mv_panel<6>(a, x, y, lda, k, renorm_interval); break;
    case 7: gaxpy_mv_panel<7>(a, x, y, lda, k, renorm_interval); break;
  }
}

// Tile any m into MR=8 row panels plus a 1..7 tail. Keep this small so
// the hot panel<8> body inlines — accumulators need to stay in registers.
static inline void gaxpy_mv_dispatch(const float64x2_t *__restrict__ a,
                                        const float64x2_t *__restrict__ x,
                                        float64x2_t *__restrict__ y, int64_t m,
                                        int64_t k, int64_t lda,
                                        int64_t renorm_interval) {
  constexpr int MR = 8;
  int64_t i = 0;
  for (; i + MR <= m; i += MR) {
    gaxpy_mv_panel<MR>(a + i, x, y + i, lda, k, renorm_interval);
  }
  int tail = static_cast<int>(m - i);
  if (tail > 0) gaxpy_mv_tail(a + i, x, y + i, tail, lda, k, renorm_interval);
}

// MR×NR GEMM µkernel: loads one column-slice of A per p and reuses it
// across NR B-column entries. Amortizes A bandwidth vs a per-column mv
// dispatch at the cost of MR*NR accumulator pairs (stack-spilled on
// machines with < 32 FP regs, but the spill cost is paid once per p
// and dwarfed by the DD-mac FLOPs).
template <int MR, int NR>
static inline void gemm_panel(const float64x2_t *__restrict__ a,
                                 const float64x2_t *__restrict__ b,
                                 float64x2_t *__restrict__ c,
                                 int64_t lda, int64_t ldb, int64_t ldc,
                                 int64_t k, int64_t renorm_interval) {
  double s_hi[MR][NR] = {}, s_lo[MR][NR] = {};
  auto kernel_block = [&](int64_t p0, int64_t pend) {
    for (int64_t p = p0; p < pend; ++p) {
      const float64x2_t *__restrict__ acol = a + p * lda;
      double ah[MR], al[MR];
      for (int i = 0; i < MR; ++i) { ah[i] = acol[i].hi; al[i] = acol[i].lo; }
      for (int jj = 0; jj < NR; ++jj) {
        const double bh = b[jj * ldb + p].hi;
        const double bl = b[jj * ldb + p].lo;
        for (int i = 0; i < MR; ++i) {
          mac_inl(ah[i], al[i], bh, bl, s_hi[i][jj], s_lo[i][jj]);
        }
      }
    }
  };
  if (renorm_interval <= 0 || k <= renorm_interval) {
    kernel_block(0, k);
  } else {
    const int64_t chunk = renorm_interval;
    int64_t p0 = 0;
    while (p0 < k) {
      int64_t pend = p0 + chunk;
      if (pend > k) pend = k;
      kernel_block(p0, pend);
      p0 = pend;
      if (p0 < k) {
        for (int jj = 0; jj < NR; ++jj)
          for (int i = 0; i < MR; ++i)
            renorm_inl(s_hi[i][jj], s_lo[i][jj]);
      }
    }
  }
  for (int jj = 0; jj < NR; ++jj)
    for (int i = 0; i < MR; ++i)
      c[jj * ldc + i] = finalize_inl(s_hi[i][jj], s_lo[i][jj]);
}
} // anonymous namespace

extern "C" {

// Arithmetic
float64x2_t adddd(float64x2_t a, float64x2_t b) { return to(from(a) + from(b)); }
float64x2_t subdd(float64x2_t a, float64x2_t b) { return to(from(a) - from(b)); }
float64x2_t muldd(float64x2_t a, float64x2_t b) { return to(from(a) * from(b)); }
float64x2_t divdd(float64x2_t a, float64x2_t b) { return to(from(a) / from(b)); }
float64x2_t negdd(float64x2_t a)  { return to(-from(a)); }
float64x2_t fabsdd(float64x2_t a)  { return to(multifloats::abs(from(a))); }
float64x2_t sqrtdd(float64x2_t a) { return to(multifloats::sqrt(from(a))); }

// Rounding (Fortran AINT / ANINT delegate here)
float64x2_t truncdd(float64x2_t a) { return to(multifloats::trunc(from(a))); }
float64x2_t rounddd(float64x2_t a) { return to(multifloats::round(from(a))); }

// Binary
float64x2_t fmindd(float64x2_t a, float64x2_t b)     { return to(multifloats::fmin(from(a), from(b))); }
float64x2_t fmaxdd(float64x2_t a, float64x2_t b)     { return to(multifloats::fmax(from(a), from(b))); }
float64x2_t hypotdd(float64x2_t a, float64x2_t b)    { return to(multifloats::hypot(from(a), from(b))); }
float64x2_t powdd(float64x2_t a, float64x2_t b)      { return to(pow_full(from(a), from(b))); }
float64x2_t fmoddd(float64x2_t a, float64x2_t b)     { return to(multifloats::fmod(from(a), from(b))); }
float64x2_t fdimdd(float64x2_t a, float64x2_t b)     { return to(multifloats::fdim(from(a), from(b))); }
float64x2_t copysigndd(float64x2_t a, float64x2_t b) { return to(multifloats::copysign(from(a), from(b))); }
float64x2_t fmadd(float64x2_t a, float64x2_t b, float64x2_t c) { return to(multifloats::fma(from(a), from(b), from(c))); }

// Exponential / logarithmic
float64x2_t expdd(float64x2_t a)   { return to(exp_full(from(a))); }
float64x2_t exp2dd(float64x2_t a)  { return to(exp2_full(from(a))); }
float64x2_t expm1dd(float64x2_t a) { return to(expm1_full(from(a))); }
float64x2_t logdd(float64x2_t a)   { return to(log_full(from(a))); }
float64x2_t log2dd(float64x2_t a)  { return to(log2_full(from(a))); }
float64x2_t log10dd(float64x2_t a) { return to(log10_full(from(a))); }
float64x2_t log1pdd(float64x2_t a) { return to(log1p_full(from(a))); }

// Trigonometric
float64x2_t sindd(float64x2_t a)   { return to(sin_full(from(a))); }
float64x2_t cosdd(float64x2_t a)   { return to(cos_full(from(a))); }
float64x2_t tandd(float64x2_t a)   { return to(tan_full(from(a))); }
float64x2_t asindd(float64x2_t a)  { return to(asin_full(from(a))); }
float64x2_t acosdd(float64x2_t a)  { return to(acos_full(from(a))); }
float64x2_t atandd(float64x2_t a)  { return to(atan_full(from(a))); }
float64x2_t atan2dd(float64x2_t a, float64x2_t b) { return to(atan2_full(from(a), from(b))); }

// π-scaled trig: fn_pi(x) = fn(π·x) or fn(x)/π.
float64x2_t sinpidd(float64x2_t a)  { return to(sinpi_full(from(a))); }
float64x2_t cospidd(float64x2_t a)  { return to(cospi_full(from(a))); }
float64x2_t tanpidd(float64x2_t a)  { return to(tanpi_full(from(a))); }
float64x2_t asinpidd(float64x2_t a) { return to(asinpi_full(from(a))); }
float64x2_t acospidd(float64x2_t a) { return to(acospi_full(from(a))); }
float64x2_t atanpidd(float64x2_t a) { return to(atanpi_full(from(a))); }
float64x2_t atan2pidd(float64x2_t a, float64x2_t b) { return to(atan2pi_full(from(a), from(b))); }

// Hyperbolic
float64x2_t sinhdd(float64x2_t a)  { return to(sinh_full(from(a))); }
float64x2_t coshdd(float64x2_t a)  { return to(cosh_full(from(a))); }
float64x2_t tanhdd(float64x2_t a)  { return to(tanh_full(from(a))); }
float64x2_t asinhdd(float64x2_t a) { return to(asinh_full(from(a))); }
float64x2_t acoshdd(float64x2_t a) { return to(acosh_full(from(a))); }
float64x2_t atanhdd(float64x2_t a) { return to(atanh_full(from(a))); }

// Error functions
float64x2_t erfdd(float64x2_t a)   { return to(erf_full(from(a))); }
float64x2_t erfcdd(float64x2_t a)  { return to(erfc_full(from(a))); }
float64x2_t erfcxdd(float64x2_t a) { return to(erfc_scaled_full(from(a))); }

// Gamma functions
float64x2_t tgammadd(float64x2_t a) { return to(tgamma_full(from(a))); }
float64x2_t lgammadd(float64x2_t a) { return to(lgamma_full(from(a))); }

// Bessel functions (math.h naming: j0, j1, y0, y1)
float64x2_t j0dd(float64x2_t a) { return to(bessel_j0_full(from(a))); }
float64x2_t j1dd(float64x2_t a) { return to(bessel_j1_full(from(a))); }
float64x2_t y0dd(float64x2_t a) { return to(bessel_y0_full(from(a))); }
float64x2_t y1dd(float64x2_t a) { return to(bessel_y1_full(from(a))); }

// Integer-order Bessel functions (recurrence from j0/j1, y0/y1).
float64x2_t jndd(int n, float64x2_t a) { return to(bessel_jn_full(n, from(a))); }
float64x2_t yndd(int n, float64x2_t a) { return to(bessel_yn_full(n, from(a))); }
void yn_rangedd(int n1, int n2, float64x2_t a, float64x2_t *out) {
  int size = n2 - n1 + 1;
  if (size <= 0) return;
  std::vector<float64x2> buf(size);
  bessel_yn_range_full(n1, n2, from(a), buf.data());
  for (int i = 0; i < size; ++i) out[i] = to(buf[i]);
}

// Fused sincos / sinhcosh. Out-pointer style (C has no multi-value return).
void sincosdd(float64x2_t a, float64x2_t *s, float64x2_t *c) {
  float64x2 ss, cc;
  sincos_full(from(a), ss, cc);
  *s = to(ss);
  *c = to(cc);
}

void sinhcoshdd(float64x2_t a, float64x2_t *s, float64x2_t *c) {
  float64x2 ss, cc;
  sinhcosh_full(from(a), ss, cc);
  *s = to(ss);
  *c = to(cc);
}

// ---- Complex DD transcendentals ------------------------------------------
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

// exp(a+bi) = e^a · (cos b + i·sin b). 2 transcendentals (exp + sincos).
complex64x2_t cexpdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 ea = exp_full(a);
  float64x2 s, c;
  sincos_full(b, s, c);
  return { to(ea * c), to(ea * s) };
}

// log(a+bi) = log(|z|) + i·atan2(b, a). Compute log(|z|) as
//   log(big) + 0.5·log(1 + (small/big)^2)
// with big = max(|a|,|b|), small = min(|a|,|b|). This avoids the
// overflow that hypot(a,b) suffers when |z| exceeds DBL_MAX even
// though log(|z|) itself is finite. atan2 handles the negative-real-
// axis branch cut (including signed-zero propagation).
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

// pow(z, w) = exp(w · log(z)). C99 G.6.4.1; 0^w folds to 0.
complex64x2_t cpowdd(complex64x2_t z, complex64x2_t w) {
  if (z.re.hi == 0.0 && z.im.hi == 0.0) return { {0.0, 0.0}, {0.0, 0.0} };
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
    float64x2 i = (b._limbs[0] < 0.0) ? -s : s;
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
  return { l.im, to(-from(l.re)) };
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

// atan(z). libstdc++ generic form:
//   Re = 0.5 · atan2(2a, 1 − a² − b²)
//   Im = 0.25 · log((a² + (b+1)²) / (a² + (b−1)²))
// Single log-ratio rather than the naive (i/2)(log(1−iz) − log(1+iz)).
complex64x2_t catandd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 a2 = a * a;
  float64x2 x = float64x2(1.0) - a2 - b * b;
  float64x2 bp1 = b + float64x2(1.0), bm1 = b - float64x2(1.0);
  float64x2 num = a2 + bp1 * bp1;
  float64x2 den = a2 + bm1 * bm1;
  float64x2 half = float64x2(0.5, 0.0), quarter = float64x2(0.25, 0.0);
  return { to(half * atan2_full(a + a, x)),
           to(quarter * log_full(num / den)) };
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

// asinh(z) = log(z + sqrt(z² + 1)). Textbook; matches libstdc++ generic.
complex64x2_t casinhdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  // 1 + z² : Re = 1 + a² − b², Im = 2ab
  complex64x2_t one_pz2 = { to(float64x2(1.0) + a * a - b * b), to(a * b + a * b) };
  complex64x2_t root = csqrtdd(one_pz2);
  complex64x2_t arg = { to(a + from(root.re)), to(b + from(root.im)) };
  return clogdd(arg);
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

// atanh(z). libstdc++ generic form (atan with a, b swapped):
//   Re = 0.25 · log((1+a)² + b²)/((1−a)² + b²))
//   Im = 0.5 · atan2(2b, 1 − a² − b²)
complex64x2_t catanhdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 b2 = b * b;
  float64x2 x = float64x2(1.0) - b2 - a * a;
  float64x2 ap1 = a + float64x2(1.0), am1 = a - float64x2(1.0);
  float64x2 num = b2 + ap1 * ap1;
  float64x2 den = b2 + am1 * am1;
  float64x2 half = float64x2(0.5, 0.0), quarter = float64x2(0.25, 0.0);
  return { to(quarter * log_full(num / den)),
           to(half * atan2_full(b + b, x)) };
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

void matmuldd_mm(const float64x2_t *__restrict__ a, const float64x2_t *__restrict__ b,
                  float64x2_t *__restrict__ c, int64_t m, int64_t k, int64_t n,
                  int64_t renorm_interval) {
  // NR-blocked mm: the MR-row × NR-col tile loads A[:,p] once per p and
  // reuses it across NR output columns, halving A-bandwidth vs a per-
  // column mv dispatch. Row-tail (1..MR-1 rows) and column-tail
  // (1..NR-1 cols) fall back to the single-column mv panel.
  constexpr int MR = 8;
  constexpr int NR = 2;
  int64_t j = 0;
  for (; j + NR <= n; j += NR) {
    int64_t i = 0;
    for (; i + MR <= m; i += MR) {
      gemm_panel<MR, NR>(a + i, b + j * k, c + j * m + i,
                            m, k, m, k, renorm_interval);
    }
    int tail_m = static_cast<int>(m - i);
    if (tail_m > 0) {
      for (int jj = 0; jj < NR; ++jj) {
        gaxpy_mv_tail(a + i, b + (j + jj) * k, c + (j + jj) * m + i,
                         tail_m, m, k, renorm_interval);
      }
    }
  }
  for (; j < n; ++j) {
    gaxpy_mv_dispatch(a, b + j * k, c + j * m, m, k, m, renorm_interval);
  }
}

void matmuldd_mv(const float64x2_t *__restrict__ a, const float64x2_t *__restrict__ x,
                  float64x2_t *__restrict__ y, int64_t m, int64_t k,
                  int64_t renorm_interval) {
  gaxpy_mv_dispatch(a, x, y, m, k, m, renorm_interval);
}

// vm: y[j] = sum_p x[p] * B[p, j]. Column-major B makes B[:, j]
// contiguous at fixed j, so one scalar accumulator per output is optimal.
void matmuldd_vm(const float64x2_t *__restrict__ x, const float64x2_t *__restrict__ b,
                  float64x2_t *__restrict__ y, int64_t k, int64_t n,
                  int64_t renorm_interval) {
  const bool simple = (renorm_interval <= 0) || (k <= renorm_interval);
  for (int64_t j = 0; j < n; ++j) {
    double s_hi = 0.0, s_lo = 0.0;
    const float64x2_t *__restrict__ bcol = b + j * k;
    if (simple) {
      for (int64_t p = 0; p < k; ++p) {
        mac_inl(x[p].hi, x[p].lo, bcol[p].hi, bcol[p].lo, s_hi, s_lo);
      }
    } else {
      const int64_t chunk = renorm_interval;
      int64_t p0 = 0;
      while (p0 < k) {
        int64_t pend = p0 + chunk;
        if (pend > k) pend = k;
        for (int64_t p = p0; p < pend; ++p) {
          mac_inl(x[p].hi, x[p].lo, bcol[p].hi, bcol[p].lo, s_hi, s_lo);
        }
        p0 = pend;
        if (p0 < k) renorm_inl(s_hi, s_lo);
      }
    }
    y[j] = finalize_inl(s_hi, s_lo);
  }
}

// Comparison
int eqdd(float64x2_t a, float64x2_t b) { return from(a) == from(b) ? 1 : 0; }
int nedd(float64x2_t a, float64x2_t b) { return from(a) != from(b) ? 1 : 0; }
int ltdd(float64x2_t a, float64x2_t b) { return from(a) <  from(b) ? 1 : 0; }
int ledd(float64x2_t a, float64x2_t b) { return from(a) <= from(b) ? 1 : 0; }
int gtdd(float64x2_t a, float64x2_t b) { return from(a) >  from(b) ? 1 : 0; }
int gedd(float64x2_t a, float64x2_t b) { return from(a) >= from(b) ? 1 : 0; }

} // extern "C"
