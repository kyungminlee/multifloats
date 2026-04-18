#include "multifloats.hh"
#include "dd_constants.hh"
#include <cstdint>
#include <cstring>

// All DD math kernels have internal linkage (anonymous namespace).
// Only the extern "C" dd_* wrappers at the bottom are exported.
namespace {

using namespace multifloats::detail;  // dd_pair, dd_neval, dd_deval, dd_horner
using MFD2 = multifloats::MultiFloat<double, 2>;

// Forward declarations for internal cross-references
MFD2 dd_exp_full(MFD2 const &x);
MFD2 dd_exp2_full(MFD2 const &x);
MFD2 dd_log_full(MFD2 const &x);
MFD2 dd_sin_full(MFD2 const &x);
MFD2 dd_cos_full(MFD2 const &x);
MFD2 dd_sin_eval(MFD2 const &r);
MFD2 dd_cos_eval(MFD2 const &r);
void dd_sincos_eval(MFD2 const &r, MFD2 &s, MFD2 &c);
void dd_sincos_full(MFD2 const &x, MFD2 &s, MFD2 &c);
void dd_sinhcosh_full(MFD2 const &x, MFD2 &s, MFD2 &c);
MFD2 dd_sinpi_full(MFD2 const &x);
MFD2 dd_erfc_full(MFD2 const &x);
MFD2 dd_lgamma_positive(MFD2 const &x);
MFD2 dd_bessel_j0_full(MFD2 const &x);
MFD2 dd_bessel_j1_full(MFD2 const &x);

// ---- exp2 / exp / log2 / log / log10 ---------------------------------------
MFD2 dd_exp2_kernel(MFD2 const &x) {
  // n = nearest integer to x.hi, y = (x - n)/8
  double n_float = std::nearbyint(x._limbs[0]);
  MFD2 y = x + dd_pair(-n_float, 0.0);
  y._limbs[0] = std::ldexp(y._limbs[0], -3);
  y._limbs[1] = std::ldexp(y._limbs[1], -3);
  // poly(y) then cube via three squarings (to undo the /8)
  MFD2 p = dd_horner(y, exp2_coefs_hi, exp2_coefs_lo, 14);
  p = p * p;
  p = p * p;
  p = p * p;
  // multiply by 2^n via two ldexps (split to keep intermediates in range)
  int n = static_cast<int>(n_float);
  int half_n = n / 2;
  MFD2 r;
  r._limbs[0] = std::ldexp(p._limbs[0], half_n);
  r._limbs[1] = std::ldexp(p._limbs[1], half_n);
  r._limbs[0] = std::ldexp(r._limbs[0], n - half_n);
  r._limbs[1] = std::ldexp(r._limbs[1], n - half_n);
  return r;
}

MFD2 dd_exp2_full(MFD2 const &x) {
  MFD2 r;
  if (!std::isfinite(x._limbs[0])) {
    if (x._limbs[0] > 0.0) r._limbs[0] = x._limbs[0]; // +inf
    else if (x._limbs[0] < 0.0) r._limbs[0] = 0.0;    // -inf → 0
    else r._limbs[0] = x._limbs[0];                    // NaN
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] < exp2_min) return MFD2();
  if (x._limbs[0] > exp2_max) {
    r._limbs[0] = std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
    return r;
  }
  return dd_exp2_kernel(x);
}

MFD2 dd_exp_full(MFD2 const &x) {
  return dd_exp2_full(x * dd_pair(log2_e_hi, log2_e_lo));
}

MFD2 dd_log2_kernel(MFD2 const &x) {
  double hi = x._limbs[0];
  if (hi > 15.0 / 16.0 && hi < 17.0 / 16.0) {
    // direct path: t = (x - 1) / (x + 1)
    MFD2 num = x - MFD2(1.0);
    MFD2 den = x + MFD2(1.0);
    MFD2 t = num / den;
    MFD2 t_sq = t * t;
    MFD2 p = dd_horner(t_sq, log2_wide_hi, log2_wide_lo, 9);
    return t * p;
  }
  // table path
  int e = std::ilogb(hi); // matches Julia/IEEE convention (mantissa in [1,2))
  MFD2 m;
  m._limbs[0] = std::ldexp(x._limbs[0], -e);
  m._limbs[1] = std::ldexp(x._limbs[1], -e);
  int idx = static_cast<int>((m._limbs[0] - 1.0) * 32.0);
  if (idx < 0) idx = 0;
  if (idx > 31) idx = 31;
  MFD2 c = dd_pair(log2_centers[idx], 0.0);
  MFD2 v = dd_pair(log2_values_hi[idx], log2_values_lo[idx]);
  MFD2 num = m - c;
  MFD2 den = m + c;
  MFD2 t = num / den;
  MFD2 t_sq = t * t;
  MFD2 p = dd_horner(t_sq, log2_narrow_hi, log2_narrow_lo, 7);
  return dd_pair(static_cast<double>(e), 0.0) + v + t * p;
}

MFD2 dd_log2_full(MFD2 const &x) {
  MFD2 r;
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
  return dd_log2_kernel(x);
}

MFD2 dd_log_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::log(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  // log2(0) = -inf. DD multiply (-inf, 0) * (ln_2_hi, ln_2_lo) would refine
  // into (-inf, NaN) because the two_prod FMA hits -inf + inf; short-circuit.
  MFD2 l2 = dd_log2_full(x);
  if (!std::isfinite(l2._limbs[0])) return l2;
  return l2 * dd_pair(ln_2_hi, ln_2_lo);
}

MFD2 dd_log10_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::log10(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 l2 = dd_log2_full(x);
  if (!std::isfinite(l2._limbs[0])) return l2;
  return l2 * dd_pair(log10_2_hi, log10_2_lo);
}

// ---- sinpi / cospi / tanpi / sin / cos / tan -------------------------------
// Family of π-scaled trig: *pi(x) = fn(π·x). Range-reduce 2·|x| to the
// nearest integer n, leaving |rx| ≤ 0.25, so |π·rx| ≤ π/4 matches the
// sin/cos eval range. Dispatch by n mod 4 for the quadrant.
MFD2 dd_sinpi_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  bool sign = x._limbs[0] < 0.0;
  MFD2 ax = sign ? -x : x;
  double n_float = std::nearbyint(2.0 * ax._limbs[0]);
  MFD2 rx = ax + dd_pair(-0.5 * n_float, 0.0);
  MFD2 pi_dd = dd_pair(pi_dd_hi, pi_dd_lo);
  MFD2 pr = pi_dd * rx;
  long long n_mod = static_cast<long long>(n_float) & 3LL;
  MFD2 res;
  switch (n_mod) {
  case 0: res = dd_sin_eval(pr); break;
  case 1: res = dd_cos_eval(pr); break;
  case 2: res = -dd_sin_eval(pr); break;
  default: res = -dd_cos_eval(pr); break;
  }
  return sign ? -res : res;
}

// cospi is even; no sign flip on the result.
MFD2 dd_cospi_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 ax = (x._limbs[0] < 0.0) ? -x : x;
  double n_float = std::nearbyint(2.0 * ax._limbs[0]);
  MFD2 rx = ax + dd_pair(-0.5 * n_float, 0.0);
  MFD2 pi_dd = dd_pair(pi_dd_hi, pi_dd_lo);
  MFD2 pr = pi_dd * rx;
  long long n_mod = static_cast<long long>(n_float) & 3LL;
  switch (n_mod) {
  case 0: return dd_cos_eval(pr);
  case 1: return -dd_sin_eval(pr);
  case 2: return -dd_cos_eval(pr);
  default: return dd_sin_eval(pr);
  }
}

// tanpi: fuses sin/cos eval so both come from one Taylor pair. Poles at
// half-integers fall out naturally (sin(0)=+0 ⇒ ±c/0 = ±inf).
MFD2 dd_tanpi_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  bool sign = x._limbs[0] < 0.0;
  MFD2 ax = sign ? -x : x;
  double n_float = std::nearbyint(2.0 * ax._limbs[0]);
  MFD2 rx = ax + dd_pair(-0.5 * n_float, 0.0);
  MFD2 pi_dd = dd_pair(pi_dd_hi, pi_dd_lo);
  MFD2 pr = pi_dd * rx;
  MFD2 s, c;
  dd_sincos_eval(pr, s, c);
  long long n_mod = static_cast<long long>(n_float) & 3LL;
  // n_mod 0, 2: tan(n·π/2 + π·rx) = tan(π·rx) = s/c
  // n_mod 1, 3: tan(n·π/2 + π·rx) = -cot(π·rx) = -c/s
  MFD2 res = (n_mod & 1LL) ? -c / s : s / c;
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
void dd_reduce_pi_half(MFD2 const &x, MFD2 &r, int &n_mod4) {
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
MFD2 dd_sin_kernel(MFD2 const &x) {
  MFD2 x2 = x * x;
  return dd_horner(x2, sin_taylor_hi, sin_taylor_lo, 13) * x;
}
MFD2 dd_cos_kernel(MFD2 const &x) {
  MFD2 x2 = x * x;
  return dd_horner(x2, cos_taylor_hi, cos_taylor_lo, 13);
}

// Evaluate sin(r)/cos(r) for |r| ≤ π/4. For |r| > π/8 shift by π/4 and use
// the angle-addition identity — this halves the polynomial argument range
// (x² ≤ 0.154 vs 0.616) so the 13-term Taylor hits full DD at the boundary.
MFD2 dd_sin_eval(MFD2 const &r) {
  constexpr double pi8 = 0.392699081698724;
  if (std::fabs(r._limbs[0]) <= pi8) return dd_sin_kernel(r);
  MFD2 pi4_dd = dd_pair(0.7853981633974483, 3.061616997868383e-17);
  MFD2 inv_sqrt2 = dd_pair(0.7071067811865476, -4.833646656726457e-17);
  bool pos = r._limbs[0] > 0.0;
  MFD2 rp = pos ? (r - pi4_dd) : (r + pi4_dd);
  MFD2 sk = dd_sin_kernel(rp);
  MFD2 ck = dd_cos_kernel(rp);
  return pos ? (sk + ck) * inv_sqrt2 : (sk - ck) * inv_sqrt2;
}

MFD2 dd_cos_eval(MFD2 const &r) {
  constexpr double pi8 = 0.392699081698724;
  if (std::fabs(r._limbs[0]) <= pi8) return dd_cos_kernel(r);
  MFD2 pi4_dd = dd_pair(0.7853981633974483, 3.061616997868383e-17);
  MFD2 inv_sqrt2 = dd_pair(0.7071067811865476, -4.833646656726457e-17);
  bool pos = r._limbs[0] > 0.0;
  MFD2 rp = pos ? (r - pi4_dd) : (r + pi4_dd);
  MFD2 sk = dd_sin_kernel(rp);
  MFD2 ck = dd_cos_kernel(rp);
  return pos ? (ck - sk) * inv_sqrt2 : (ck + sk) * inv_sqrt2;
}

// Fused sin/cos of |r| ≤ π/4 — shares the π/4 shift path and its two
// Taylor kernels between both outputs. Callers that need both s=sin(r)
// and c=cos(r) save one 13-term Horner (~40% of the eval cost).
void dd_sincos_eval(MFD2 const &r, MFD2 &s, MFD2 &c) {
  constexpr double pi8 = 0.392699081698724;
  if (std::fabs(r._limbs[0]) <= pi8) {
    s = dd_sin_kernel(r);
    c = dd_cos_kernel(r);
    return;
  }
  MFD2 pi4_dd = dd_pair(0.7853981633974483, 3.061616997868383e-17);
  MFD2 inv_sqrt2 = dd_pair(0.7071067811865476, -4.833646656726457e-17);
  bool pos = r._limbs[0] > 0.0;
  MFD2 rp = pos ? (r - pi4_dd) : (r + pi4_dd);
  MFD2 sk = dd_sin_kernel(rp);
  MFD2 ck = dd_cos_kernel(rp);
  if (pos) {
    s = (sk + ck) * inv_sqrt2;
    c = (ck - sk) * inv_sqrt2;
  } else {
    s = (sk - ck) * inv_sqrt2;
    c = (ck + sk) * inv_sqrt2;
  }
}

// Fused sin/cos of the full-range argument — one call to dd_reduce_pi_half
// plus one dd_sincos_eval, with the per-quadrant sign/swap applied last.
void dd_sincos_full(MFD2 const &x, MFD2 &s, MFD2 &c) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 nan;
    nan._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    nan._limbs[1] = 0.0;
    s = nan;
    c = nan;
    return;
  }
  MFD2 r;
  int q;
  dd_reduce_pi_half(x, r, q);
  MFD2 ss, cc;
  dd_sincos_eval(r, ss, cc);
  switch (q) {
    case 0: s = ss;  c = cc;  break;
    case 1: s = cc;  c = -ss; break;
    case 2: s = -ss; c = -cc; break;
    default: s = -cc; c = ss; break;
  }
}

MFD2 dd_sin_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 r;
  int q;
  dd_reduce_pi_half(x, r, q);
  switch (q) {
  case 0: return dd_sin_eval(r);
  case 1: return dd_cos_eval(r);
  case 2: return -dd_sin_eval(r);
  default: return -dd_cos_eval(r);
  }
}

MFD2 dd_cos_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 r;
  int q;
  dd_reduce_pi_half(x, r, q);
  switch (q) {
  case 0: return dd_cos_eval(r);
  case 1: return -dd_sin_eval(r);
  case 2: return -dd_cos_eval(r);
  default: return dd_sin_eval(r);
  }
}

MFD2 dd_tan_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 r;
  int q;
  dd_reduce_pi_half(x, r, q);
  MFD2 s, c;
  dd_sincos_eval(r, s, c);
  switch (q) {
    case 0:
    case 2:  return s / c;
    default: return -c / s;  // q == 1 or 3
  }
}

// ---- sinh / cosh / tanh ----------------------------------------------------
MFD2 dd_sinh_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::sinh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) < 0.1) {
    // 9-term Taylor: x * (1 + x²/6 + x⁴/120 + ...)
    return x * dd_horner(x * x, sinh_taylor_hi, sinh_taylor_lo, 9);
  }
  MFD2 e = dd_exp_full(x);
  MFD2 ei = dd_exp_full(-x);
  return (e - ei) * dd_pair(0.5, 0.0);
}

MFD2 dd_cosh_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::cosh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 e = dd_exp_full(x);
  MFD2 ei = dd_exp_full(-x);
  return (e + ei) * dd_pair(0.5, 0.0);
}

// Fused sinh/cosh. Shares the two exp evaluations (exp(x) and exp(-x)) that
// dd_sinh_full and dd_cosh_full would do independently, halving the exp cost
// for call sites that need both (e.g. complex sin/cos/sinh/cosh).
void dd_sinhcosh_full(MFD2 const &x, MFD2 &s, MFD2 &c) {
  if (!std::isfinite(x._limbs[0])) {
    s._limbs[0] = std::sinh(x._limbs[0]); s._limbs[1] = 0.0;
    c._limbs[0] = std::cosh(x._limbs[0]); c._limbs[1] = 0.0;
    return;
  }
  if (std::abs(x._limbs[0]) < 0.1) {
    // Small-|x| Taylor: sinh via shared coefficients, cosh via 1 + x²·Q(x²).
    // The cosh Taylor here mirrors the structure of dd_sinh_full's path.
    MFD2 x2 = x * x;
    s = x * dd_horner(x2, sinh_taylor_hi, sinh_taylor_lo, 9);
    // cosh(x) = 1 + x²/2 + x⁴/24 + … — reuse separate kernel to keep code small.
    c = dd_cosh_full(x);
    return;
  }
  MFD2 e  = dd_exp_full(x);
  MFD2 ei = dd_exp_full(-x);
  MFD2 half = dd_pair(0.5, 0.0);
  s = (e - ei) * half;
  c = (e + ei) * half;
}

MFD2 dd_tanh_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::tanh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  // For |x| > 37, exp(-2|x|) < 2^-106, so tanh(x) = ±1 to full DD.
  if (x._limbs[0] > 37.0) return MFD2(1.0);
  if (x._limbs[0] < -37.0) return MFD2(-1.0);
  if (std::abs(x._limbs[0]) < 0.5) {
    return dd_sinh_full(x) / dd_cosh_full(x);
  }
  // tanh(|x|) = (1 - em2) / (1 + em2) where em2 = exp(-2|x|) is small
  // for |x| > 0.5 — avoids the 1 − small cancellation of 1 − 2/(e^(2|x|)+1).
  bool sign = x._limbs[0] < 0.0;
  MFD2 neg_two_ax;
  neg_two_ax._limbs[0] = sign ? (2.0 * x._limbs[0]) : (-2.0 * x._limbs[0]);
  neg_two_ax._limbs[1] = sign ? (2.0 * x._limbs[1]) : (-2.0 * x._limbs[1]);
  MFD2 em2 = dd_exp_full(neg_two_ax);
  MFD2 res = (MFD2(1.0) - em2) / (MFD2(1.0) + em2);
  return sign ? -res : res;
}

// ---- pow -------------------------------------------------------------------
MFD2 dd_pow_full(MFD2 const &x, MFD2 const &y) {
  if (x._limbs[0] == 0.0 && y._limbs[0] == 0.0) return MFD2(1.0);
  return dd_exp_full(y * dd_log_full(x));
}

// ---- atan (table lookup + rational polynomial, from libquadmath atanq.c) ---
// arctan(t) = t + t^3 * P(t^2) / Q(t^2),  |t| <= 0.09375
// Argument reduced via 84-entry table: arctan(x) = atantbl[k] + arctan(t)
// where t = (x - k/8) / (1 + x*k/8).
MFD2 dd_atan_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::atan(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  double ax_d = std::abs(x._limbs[0]);
  // atan(x) = x - x³/3 + …; for |x| < sqrt(3)·2⁻⁵³ ≈ 1.85e-16 the
  // cubic term is below 0.5 DD ulp, so identity is correct to DD.
  if (ax_d < 1.85e-16) return x;

  bool neg = x._limbs[0] < 0.0;
  MFD2 ax = neg ? -x : x;

  int k;
  MFD2 t;
  if (ax_d >= 10.25) {
    k = 83;
    t = MFD2(-1.0) / ax;
  } else {
    k = static_cast<int>(8.0 * ax_d + 0.25);
    double u_d = 0.125 * k;
    MFD2 u = dd_pair(u_d, 0.0);  // exact in double
    t = (ax - u) / (MFD2(1.0) + ax * u);
  }

  // arctan(t) = t + t^3 * P(t^2) / Q(t^2)
  MFD2 t2 = t * t;
  MFD2 p = dd_neval(t2, atan_P_hi, atan_P_lo, 4);
  MFD2 q = dd_deval(t2, atan_Q_hi, atan_Q_lo, 4);
  MFD2 res = dd_pair(atan_table_hi[k], atan_table_lo[k])
           + t + t * t2 * p / q;

  return neg ? -res : res;
}

// ---- asin (piecewise rational, from libquadmath asinq.c) -------------------
// Region 1: |x| < 0.5      — asin(x) = x + x*x^2*P(x^2)/Q(x^2)
// Region 2: 0.5 <= |x| < 0.625 — centered at 0.5625
// Region 3: 0.625 <= |x| < 1   — half-angle: asin(x) = pi/2 - 2*asin(sqrt((1-x)/2))
MFD2 dd_asin_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::asin(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  double ax_d = std::abs(x._limbs[0]);
  if (ax_d > 1.0) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  // asin(x) = x + x³/6 + …; for |x| < sqrt(3)·2⁻⁵³ ≈ 1.85e-16 the
  // cubic term is below 0.5 DD ulp, so identity is correct to DD.
  if (ax_d < 1.85e-16) return x;

  bool neg = x._limbs[0] < 0.0;
  MFD2 ax = neg ? -x : x;
  MFD2 res;

  if (ax_d < 0.5) {
    // Region 1: asin(x) = x + x*x^2*P(x^2)/Q(x^2)
    MFD2 t = ax * ax;
    MFD2 p = dd_neval(t, asin_pS_hi, asin_pS_lo, 9);
    MFD2 q = dd_deval(t, asin_qS_hi, asin_qS_lo, 8);
    res = ax + ax * t * p / q;
  } else if (ax_d < 0.625) {
    // Region 2: centered at 0.5625
    MFD2 t = ax - MFD2(0.5625);
    MFD2 p = t * dd_neval(t, asin_rS_hi, asin_rS_lo, 10);
    MFD2 q = dd_deval(t, asin_sS_hi, asin_sS_lo, 9);
    res = dd_pair(asinr5625_hi, asinr5625_lo) + p / q;
  } else {
    // Region 3: half-angle identity
    MFD2 w = MFD2(1.0) - ax;
    MFD2 t = w * MFD2(0.5);
    MFD2 s = sqrt(t);
    MFD2 p = dd_neval(t, asin_pS_hi, asin_pS_lo, 9);
    MFD2 q = dd_deval(t, asin_qS_hi, asin_qS_lo, 8);
    MFD2 w2 = t * p / q;
    res = dd_pair(pi_half_cw1, pi_half_cw2) - MFD2(2.0) * (s + s * w2);
  }

  return neg ? -res : res;
}

// ---- acos (derived from asin polynomial) ------------------------------------
MFD2 dd_acos_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::acos(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  double ax_d = std::abs(x._limbs[0]);
  if (ax_d > 1.0) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }

  MFD2 pi_half = dd_pair(pi_half_cw1, pi_half_cw2);

  if (ax_d <= 0.5) {
    // acos(x) = pi/2 - asin(x) — no cancellation since result ∈ [π/3, 2π/3]
    return pi_half - dd_asin_full(x);
  } else if (x._limbs[0] > 0.0) {
    // x > 0.5: acos(x) = 2*asin(sqrt((1-x)/2))
    MFD2 t = (MFD2(1.0) - x) * MFD2(0.5);
    MFD2 s = sqrt(t);
    MFD2 p = dd_neval(t, asin_pS_hi, asin_pS_lo, 9);
    MFD2 q = dd_deval(t, asin_qS_hi, asin_qS_lo, 8);
    return MFD2(2.0) * (s + s * t * p / q);
  } else {
    // x < -0.5: acos(x) = pi - 2*asin(sqrt((1+x)/2))
    MFD2 t = (MFD2(1.0) + x) * MFD2(0.5);
    MFD2 s = sqrt(t);
    MFD2 p = dd_neval(t, asin_pS_hi, asin_pS_lo, 9);
    MFD2 q = dd_deval(t, asin_qS_hi, asin_qS_lo, 8);
    return dd_pair(pi_dd_hi, pi_dd_lo) - MFD2(2.0) * (s + s * t * p / q);
  }
}

// ---- atan2 ------------------------------------------------------------------
MFD2 dd_atan2_full(MFD2 const &y, MFD2 const &x) {
  if (!std::isfinite(x._limbs[0]) || !std::isfinite(y._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::atan2(y._limbs[0], x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] == 0.0 && y._limbs[0] == 0.0) {
    MFD2 r;
    r._limbs[0] = std::atan2(y._limbs[0], x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 pi_dd = dd_pair(pi_dd_hi, pi_dd_lo);
  MFD2 pi_half_dd = dd_pair(pi_half_cw1, pi_half_cw2);
  MFD2 res;
  if (std::abs(x._limbs[0]) >= std::abs(y._limbs[0])) {
    res = dd_atan_full(y / x);
    if (x._limbs[0] < 0.0) {
      if (y._limbs[0] >= 0.0) res = res + pi_dd;
      else res = res - pi_dd;
    }
  } else {
    res = dd_atan_full(x / y);
    if (y._limbs[0] > 0.0) res = pi_half_dd - res;
    else res = -pi_half_dd - res;
  }
  return res;
}

// ---- inverse π-scaled trig -------------------------------------------------
// fn_pi(x) = fn(x) / π via the natural path. The DD multiply by inv_pi
// costs one DD mul (~3 FLOPs) after the main eval.
MFD2 dd_asinpi_full(MFD2 const &x) {
  return dd_asin_full(x) * dd_pair(inv_pi_hi, inv_pi_lo);
}
MFD2 dd_acospi_full(MFD2 const &x) {
  return dd_acos_full(x) * dd_pair(inv_pi_hi, inv_pi_lo);
}
// atanpi: at ±inf, dd_atan_full bounces through std::atan and returns π/2 to
// only dp precision — then inv_pi multiplication can't recover the missing
// low limb. Short-circuit so atanpi(±inf) = ±0.5 exactly.
MFD2 dd_atanpi_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    if (std::isnan(x._limbs[0])) {
      r._limbs[0] = x._limbs[0];
    } else {
      r._limbs[0] = std::signbit(x._limbs[0]) ? -0.5 : 0.5;
    }
    r._limbs[1] = 0.0;
    return r;
  }
  return dd_atan_full(x) * dd_pair(inv_pi_hi, inv_pi_lo);
}
MFD2 dd_atan2pi_full(MFD2 const &y, MFD2 const &x) {
  return dd_atan2_full(y, x) * dd_pair(inv_pi_hi, inv_pi_lo);
}

// ---- asinh / acosh / atanh -------------------------------------------------
MFD2 dd_asinh_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::asinh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) < 0.01) {
    return x * dd_horner(x * x, asinh_taylor_hi, asinh_taylor_lo, 15);
  }
  bool sign = x._limbs[0] < 0.0;
  MFD2 ax = sign ? -x : x;
  if (ax._limbs[0] > 1e150) {
    // log(2|x|) asymptotic to avoid x² overflow
    MFD2 r = dd_log_full(ax) + dd_pair(ln_2_hi, ln_2_lo);
    return sign ? -r : r;
  }
  MFD2 root = sqrt(ax * ax + MFD2(1.0));
  MFD2 res = dd_log_full(ax + root);
  return sign ? -res : res;
}

MFD2 dd_acosh_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::acosh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] < 1.0) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] > 1e150) {
    return dd_log_full(x) + dd_pair(ln_2_hi, ln_2_lo);
  }
  MFD2 root = sqrt(x * x - MFD2(1.0));
  return dd_log_full(x + root);
}

MFD2 dd_atanh_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::atanh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) > 1.0) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) < 0.01) {
    return x * dd_horner(x * x, atanh_taylor_hi, atanh_taylor_lo, 15);
  }
  MFD2 num = MFD2(1.0) + x;
  MFD2 den = MFD2(1.0) - x;
  return dd_pair(0.5, 0.0) * dd_log_full(num / den);
}

// ---- erf / erfc (piecewise rational, ported from libquadmath erfq.c) -------
MFD2 dd_erf_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::erf(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  bool neg = x._limbs[0] < 0.0;
  MFD2 ax = neg ? -x : x;

  if (ax._limbs[0] >= 1.0) {
    if (ax._limbs[0] >= 16.0) return neg ? MFD2(-1.0) : MFD2(1.0);
    MFD2 res = MFD2(1.0) - dd_erfc_full(ax);
    return neg ? -res : res;
  }

  MFD2 y;
  if (ax._limbs[0] < 0.875) {
    if (ax._limbs[0] < 1e-18)
      return x + dd_pair(erf_efx_hi, erf_efx_lo) * x;
    MFD2 z = ax * ax;
    y = ax + ax * dd_neval(z, erf_TN1_hi, erf_TN1_lo, 8) /
                  dd_deval(z, erf_TD1_hi, erf_TD1_lo, 8);
  } else {
    MFD2 a = ax - MFD2(1.0);
    y = dd_pair(erf_const_hi, erf_const_lo) +
        dd_neval(a, erf_TN2_hi, erf_TN2_lo, 8) /
        dd_deval(a, erf_TD2_hi, erf_TD2_lo, 8);
  }
  return neg ? -y : y;
}

MFD2 dd_erfc_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::erfc(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  bool neg = x._limbs[0] < 0.0;
  MFD2 ax = neg ? -x : x;

  if (ax._limbs[0] < 0.25) {
    return MFD2(1.0) - dd_erf_full(x);
  }

  if (ax._limbs[0] > 107.0) {
    return neg ? MFD2(2.0) : MFD2();
  }

  MFD2 res;
  if (ax._limbs[0] < 1.25) {
    // 8 sub-intervals of width 0.125 in [0.25, 1.25)
    int i = static_cast<int>(8.0 * ax._limbs[0]);
    MFD2 z = ax - MFD2(erfc_x0[i - 2]);
    MFD2 y;
    switch (i) {
    case 2: y = dd_neval(z, erfc_RN13_hi, erfc_RN13_lo, 8) / dd_deval(z, erfc_RD13_hi, erfc_RD13_lo, 7); break;
    case 3: y = dd_neval(z, erfc_RN14_hi, erfc_RN14_lo, 8) / dd_deval(z, erfc_RD14_hi, erfc_RD14_lo, 7); break;
    case 4: y = dd_neval(z, erfc_RN15_hi, erfc_RN15_lo, 8) / dd_deval(z, erfc_RD15_hi, erfc_RD15_lo, 7); break;
    case 5: y = dd_neval(z, erfc_RN16_hi, erfc_RN16_lo, 8) / dd_deval(z, erfc_RD16_hi, erfc_RD16_lo, 7); break;
    case 6: y = dd_neval(z, erfc_RN17_hi, erfc_RN17_lo, 8) / dd_deval(z, erfc_RD17_hi, erfc_RD17_lo, 7); break;
    case 7: y = dd_neval(z, erfc_RN18_hi, erfc_RN18_lo, 8) / dd_deval(z, erfc_RD18_hi, erfc_RD18_lo, 7); break;
    case 8: y = dd_neval(z, erfc_RN19_hi, erfc_RN19_lo, 8) / dd_deval(z, erfc_RD19_hi, erfc_RD19_lo, 7); break;
    default: y = dd_neval(z, erfc_RN20_hi, erfc_RN20_lo, 8) / dd_deval(z, erfc_RD20_hi, erfc_RD20_lo, 7); break;
    }
    int ci = i - 2;
    res = y * z + dd_pair(erfc_Cb_hi[ci], erfc_Cb_lo[ci]) +
          dd_pair(erfc_Ca_hi[ci], erfc_Ca_lo[ci]);
  } else {
    // Asymptotic: erfc(x) = (1/x)*exp(-x^2 - 0.5625 + R(1/x^2))
    if (neg && ax._limbs[0] >= 9.0) return MFD2(2.0);

    MFD2 x2 = ax * ax;
    MFD2 z = MFD2(1.0) / x2;
    int i = static_cast<int>(8.0 / ax._limbs[0]);
    MFD2 p;
    switch (i) {
    case 0: p = dd_neval(z, erfc_AN1_hi, erfc_AN1_lo, 9) / dd_deval(z, erfc_AD1_hi, erfc_AD1_lo, 8); break;
    case 1: p = dd_neval(z, erfc_AN2_hi, erfc_AN2_lo, 11) / dd_deval(z, erfc_AD2_hi, erfc_AD2_lo, 10); break;
    case 2: p = dd_neval(z, erfc_AN3_hi, erfc_AN3_lo, 11) / dd_deval(z, erfc_AD3_hi, erfc_AD3_lo, 10); break;
    case 3: p = dd_neval(z, erfc_AN4_hi, erfc_AN4_lo, 10) / dd_deval(z, erfc_AD4_hi, erfc_AD4_lo, 10); break;
    case 4: p = dd_neval(z, erfc_AN5_hi, erfc_AN5_lo, 10) / dd_deval(z, erfc_AD5_hi, erfc_AD5_lo, 9); break;
    case 5: p = dd_neval(z, erfc_AN6_hi, erfc_AN6_lo, 9) / dd_deval(z, erfc_AD6_hi, erfc_AD6_lo, 9); break;
    case 6: p = dd_neval(z, erfc_AN7_hi, erfc_AN7_lo, 9) / dd_deval(z, erfc_AD7_hi, erfc_AD7_lo, 9); break;
    default: p = dd_neval(z, erfc_AN8_hi, erfc_AN8_lo, 9) / dd_deval(z, erfc_AD8_hi, erfc_AD8_lo, 8); break;
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
    MFD2 s = dd_pair(s_hi, 0.0);

    MFD2 e1 = dd_exp_full(-(s * s) - MFD2(0.5625));
    MFD2 diff_sq = (s - ax) * (s + ax);
    MFD2 e2 = dd_exp_full(diff_sq + p);
    res = (e1 * e2) / ax;
  }

  return neg ? MFD2(2.0) - res : res;
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
MFD2 dd_erfc_scaled_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    if (x._limbs[0] > 0.0)      r._limbs[0] = 0.0;                                            // erfc_scaled(+inf) = 0
    else if (x._limbs[0] < 0.0) r._limbs[0] = std::numeric_limits<double>::infinity();        // erfc_scaled(-inf) = +inf
    else                         r._limbs[0] = x._limbs[0];                                    // NaN
    r._limbs[1] = 0.0;
    return r;
  }
  bool neg = x._limbs[0] < 0.0;
  MFD2 ax = neg ? -x : x;
  if (ax._limbs[0] < 1.25) {
    return dd_exp_full(x * x) * dd_erfc_full(x);
  }
  // Asymptotic branch.
  MFD2 x2 = ax * ax;
  MFD2 z = MFD2(1.0) / x2;
  int i = static_cast<int>(8.0 / ax._limbs[0]);
  if (i > 7) i = 7;
  MFD2 p;
  switch (i) {
  case 0: p = dd_neval(z, erfc_AN1_hi, erfc_AN1_lo, 9)  / dd_deval(z, erfc_AD1_hi, erfc_AD1_lo, 8);  break;
  case 1: p = dd_neval(z, erfc_AN2_hi, erfc_AN2_lo, 11) / dd_deval(z, erfc_AD2_hi, erfc_AD2_lo, 10); break;
  case 2: p = dd_neval(z, erfc_AN3_hi, erfc_AN3_lo, 11) / dd_deval(z, erfc_AD3_hi, erfc_AD3_lo, 10); break;
  case 3: p = dd_neval(z, erfc_AN4_hi, erfc_AN4_lo, 10) / dd_deval(z, erfc_AD4_hi, erfc_AD4_lo, 10); break;
  case 4: p = dd_neval(z, erfc_AN5_hi, erfc_AN5_lo, 10) / dd_deval(z, erfc_AD5_hi, erfc_AD5_lo, 9);  break;
  case 5: p = dd_neval(z, erfc_AN6_hi, erfc_AN6_lo, 9)  / dd_deval(z, erfc_AD6_hi, erfc_AD6_lo, 9);  break;
  case 6: p = dd_neval(z, erfc_AN7_hi, erfc_AN7_lo, 9)  / dd_deval(z, erfc_AD7_hi, erfc_AD7_lo, 9);  break;
  default:p = dd_neval(z, erfc_AN8_hi, erfc_AN8_lo, 9)  / dd_deval(z, erfc_AD8_hi, erfc_AD8_lo, 8);  break;
  }
  MFD2 res = dd_exp_full(p - MFD2(0.5625)) / ax;
  if (neg) {
    MFD2 exp_x2 = dd_exp_full(x2);
    if (!std::isfinite(exp_x2._limbs[0])) {
      MFD2 r;
      r._limbs[0] = std::numeric_limits<double>::infinity();
      r._limbs[1] = 0.0;
      return r;
    }
    res = MFD2(2.0) * exp_x2 - res;
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
// subtracting a single dd_log_full(prod).

// Stirling/gamma constants are in dd_constants.hh (via #include above).

// Evaluate the Stirling series for log Γ(x). Assumes x is well into the
// asymptotic range (x ≥ ~20) — callers shift first.
MFD2 dd_lgamma_stirling(MFD2 const &x) {
  MFD2 inv_x = MFD2(1.0) / x;
  MFD2 inv_x2 = inv_x * inv_x;
  MFD2 poly =
      dd_horner(inv_x2, stirling_coefs_hi, stirling_coefs_lo, 13);
  MFD2 corr = poly * inv_x;
  MFD2 logx = dd_log_full(x);
  MFD2 xm_half = x - MFD2(0.5);
  return xm_half * logx - x +
         dd_pair(half_log_2pi_hi, half_log_2pi_lo) + corr;
}

// Exact DD subtraction: z = x - c, where c is exact in double precision.
static inline MFD2 dd_sub_exact(MFD2 const &x, double c) {
  MFD2 z;
  z._limbs[0] = x._limbs[0] - c;
  z._limbs[1] = x._limbs[1] + ((x._limbs[0] - z._limbs[0]) - c);
  return z;
}

// Piecewise rational approximation for lgamma(x), 0.5 <= x <= 13.5.
// Mirrors lgammaq.c structure: P(z)/Q(z) with z = x - center.
MFD2 dd_lgamma_rational(MFD2 const &x) {
  int nn = (int)std::nearbyint(x._limbs[0]);
  MFD2 z, p, q, c;
  switch (nn) {
  case 0: // x ≈ 0.5: lgamma(x+1) via near-minimum, then - log(x)
    z = x + MFD2(1.0) - dd_pair(lgam_x0_hi, lgam_x0_lo);
    p = dd_neval(z, lgam_RN1r5_hi, lgam_RN1r5_lo, 8);
    q = dd_deval(z, lgam_RD1r5_hi, lgam_RD1r5_lo, 8);
    return z * z * (p / q) + dd_pair(lgam_y0_hi, lgam_y0_lo) - dd_log_full(x);
  case 1:
    if (x._limbs[0] < 0.875) {
      if (x._limbs[0] <= 0.625) {
        z = x + MFD2(1.0) - dd_pair(lgam_x0_hi, lgam_x0_lo);
        p = dd_neval(z, lgam_RN1r5_hi, lgam_RN1r5_lo, 8);
        q = dd_deval(z, lgam_RD1r5_hi, lgam_RD1r5_lo, 8);
        return z * z * (p / q) + dd_pair(lgam_y0_hi, lgam_y0_lo) - dd_log_full(x);
      }
      // (0.625, 0.875): z = x - 0.75, lgamma(x+1) = lgam1r75 + z*P/Q
      z = dd_sub_exact(x, 0.75);
      p = dd_neval(z, lgam_RN1r75_hi, lgam_RN1r75_lo, 8);
      q = dd_deval(z, lgam_RD1r75_hi, lgam_RD1r75_lo, 8);
      return z * (p / q) + dd_pair(lgam1r75_hi, lgam1r75_lo) - dd_log_full(x);
    }
    if (x._limbs[0] < 1.0) {
      z = x - MFD2(1.0);
      p = dd_neval(z, lgam_RNr9_hi, lgam_RNr9_lo, 8);
      q = dd_deval(z, lgam_RDr9_hi, lgam_RDr9_lo, 8);
      return z * (p / q);
    }
    if (x._limbs[0] <= 1.125) {
      z = x - MFD2(1.0);
      p = dd_neval(z, lgam_RN1_hi, lgam_RN1_lo, 8);
      q = dd_deval(z, lgam_RD1_hi, lgam_RD1_lo, 7);
      return z * (p / q);
    }
    if (x._limbs[0] <= 1.375) {
      z = dd_sub_exact(x, 1.25);
      p = dd_neval(z, lgam_RN1r25_hi, lgam_RN1r25_lo, 9);
      q = dd_deval(z, lgam_RD1r25_hi, lgam_RD1r25_lo, 8);
      return z * (p / q) + dd_pair(lgam1r25_hi, lgam1r25_lo);
    }
    // [1.375, 1.5]: near minimum x0
    z = x - dd_pair(lgam_x0_hi, lgam_x0_lo);
    p = dd_neval(z, lgam_RN1r5_hi, lgam_RN1r5_lo, 8);
    q = dd_deval(z, lgam_RD1r5_hi, lgam_RD1r5_lo, 8);
    return z * z * (p / q) + dd_pair(lgam_y0_hi, lgam_y0_lo);
  case 2:
    if (x._limbs[0] < 1.625) {
      z = x - dd_pair(lgam_x0_hi, lgam_x0_lo);
      p = dd_neval(z, lgam_RN1r5_hi, lgam_RN1r5_lo, 8);
      q = dd_deval(z, lgam_RD1r5_hi, lgam_RD1r5_lo, 8);
      return z * z * (p / q) + dd_pair(lgam_y0_hi, lgam_y0_lo);
    }
    if (x._limbs[0] < 1.875) {
      z = dd_sub_exact(x, 1.75);
      p = dd_neval(z, lgam_RN1r75_hi, lgam_RN1r75_lo, 8);
      q = dd_deval(z, lgam_RD1r75_hi, lgam_RD1r75_lo, 8);
      return z * (p / q) + dd_pair(lgam1r75_hi, lgam1r75_lo);
    }
    if (x._limbs[0] < 2.375) {
      z = dd_sub_exact(x, 2.0);
      p = dd_neval(z, lgam_RN2_hi, lgam_RN2_lo, 9);
      q = dd_deval(z, lgam_RD2_hi, lgam_RD2_lo, 9);
      return z * (p / q);
    }
    z = dd_sub_exact(x, 2.5);
    p = dd_neval(z, lgam_RN2r5_hi, lgam_RN2r5_lo, 8);
    q = dd_deval(z, lgam_RD2r5_hi, lgam_RD2r5_lo, 8);
    return z * (p / q) + dd_pair(lgam2r5_hi, lgam2r5_lo);
  case 3:
    if (x._limbs[0] < 2.75) {
      z = dd_sub_exact(x, 2.5);
      p = dd_neval(z, lgam_RN2r5_hi, lgam_RN2r5_lo, 8);
      q = dd_deval(z, lgam_RD2r5_hi, lgam_RD2r5_lo, 8);
      return z * (p / q) + dd_pair(lgam2r5_hi, lgam2r5_lo);
    }
    z = dd_sub_exact(x, 3.0);
    p = dd_neval(z, lgam_RN3_hi, lgam_RN3_lo, 9);
    q = dd_deval(z, lgam_RD3_hi, lgam_RD3_lo, 9);
    return z * (p / q) + dd_pair(lgam3_hi, lgam3_lo);
  case 4:
    z = dd_sub_exact(x, 4.0);
    p = dd_neval(z, lgam_RN4_hi, lgam_RN4_lo, 9);
    q = dd_deval(z, lgam_RD4_hi, lgam_RD4_lo, 9);
    return z * (p / q) + dd_pair(lgam4_hi, lgam4_lo);
  case 5:
    z = dd_sub_exact(x, 5.0);
    p = dd_neval(z, lgam_RN5_hi, lgam_RN5_lo, 9);
    q = dd_deval(z, lgam_RD5_hi, lgam_RD5_lo, 8);
    return z * (p / q) + dd_pair(lgam5_hi, lgam5_lo);
  case 6:
    z = dd_sub_exact(x, 6.0);
    p = dd_neval(z, lgam_RN6_hi, lgam_RN6_lo, 8);
    q = dd_deval(z, lgam_RD6_hi, lgam_RD6_lo, 8);
    return z * (p / q) + dd_pair(lgam6_hi, lgam6_lo);
  case 7:
    z = dd_sub_exact(x, 7.0);
    p = dd_neval(z, lgam_RN7_hi, lgam_RN7_lo, 8);
    q = dd_deval(z, lgam_RD7_hi, lgam_RD7_lo, 7);
    return z * (p / q) + dd_pair(lgam7_hi, lgam7_lo);
  case 8:
    z = dd_sub_exact(x, 8.0);
    p = dd_neval(z, lgam_RN8_hi, lgam_RN8_lo, 8);
    q = dd_deval(z, lgam_RD8_hi, lgam_RD8_lo, 7);
    return z * (p / q) + dd_pair(lgam8_hi, lgam8_lo);
  case 9:
    z = dd_sub_exact(x, 9.0);
    p = dd_neval(z, lgam_RN9_hi, lgam_RN9_lo, 7);
    q = dd_deval(z, lgam_RD9_hi, lgam_RD9_lo, 7);
    return z * (p / q) + dd_pair(lgam9_hi, lgam9_lo);
  case 10:
    z = dd_sub_exact(x, 10.0);
    p = dd_neval(z, lgam_RN10_hi, lgam_RN10_lo, 7);
    q = dd_deval(z, lgam_RD10_hi, lgam_RD10_lo, 7);
    return z * (p / q) + dd_pair(lgam10_hi, lgam10_lo);
  case 11:
    z = dd_sub_exact(x, 11.0);
    p = dd_neval(z, lgam_RN11_hi, lgam_RN11_lo, 7);
    q = dd_deval(z, lgam_RD11_hi, lgam_RD11_lo, 6);
    return z * (p / q) + dd_pair(lgam11_hi, lgam11_lo);
  case 12:
    z = dd_sub_exact(x, 12.0);
    p = dd_neval(z, lgam_RN12_hi, lgam_RN12_lo, 7);
    q = dd_deval(z, lgam_RD12_hi, lgam_RD12_lo, 6);
    return z * (p / q) + dd_pair(lgam12_hi, lgam12_lo);
  case 13:
    z = dd_sub_exact(x, 13.0);
    p = dd_neval(z, lgam_RN13_hi, lgam_RN13_lo, 7);
    q = dd_deval(z, lgam_RD13_hi, lgam_RD13_lo, 6);
    return z * (p / q) + dd_pair(lgam13_hi, lgam13_lo);
  default:
    return dd_lgamma_stirling(x);
  }
}

// Compute lgamma for positive x >= 0.5.
MFD2 dd_lgamma_positive(MFD2 const &x) {
  if (x._limbs[0] >= 13.5)
    return dd_lgamma_stirling(x);
  return dd_lgamma_rational(x);
}

MFD2 dd_lgamma_full(MFD2 const &x) {
  MFD2 r;
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
    MFD2 s = multifloats::abs(dd_sinpi_full(x));
    MFD2 one_minus_x = MFD2(1.0) - x;
    MFD2 lgam_1mx = dd_lgamma_positive(one_minus_x);
    return dd_pair(log_pi_hi, log_pi_lo) - dd_log_full(s) - lgam_1mx;
  }
  return dd_lgamma_positive(x);
}

MFD2 dd_tgamma_full(MFD2 const &x) {
  MFD2 r;
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
    MFD2 s = dd_sinpi_full(x);
    MFD2 one_minus_x = MFD2(1.0) - x;
    // Use lgamma path for the (1-x) branch so we avoid a double tgamma
    // recursion and keep the large factor inside a log.
    MFD2 lg = dd_lgamma_positive(one_minus_x);
    // Γ(1-x) = exp(lg). sin(πx) may be negative → track sign separately.
    bool neg = s._limbs[0] < 0.0;
    if (neg) {
      s._limbs[0] = -s._limbs[0];
      s._limbs[1] = -s._limbs[1];
    }
    MFD2 pi_dd = dd_pair(pi_dd_hi, pi_dd_lo);
    MFD2 g1mx = dd_exp_full(lg);
    MFD2 out = pi_dd / (s * g1mx);
    if (neg) {
      out._limbs[0] = -out._limbs[0];
      out._limbs[1] = -out._limbs[1];
    }
    return out;
  }
  return dd_exp_full(dd_lgamma_positive(x));
}


// =============================================================================
// Bessel functions — piecewise rational approximation from libquadmath j0q/j1q
// =============================================================================

// Shared asymptotic P(z),Q(z) selection for order 0 (used by both j0 and y0).
// z = 1/x^2, xinv_d = 1/x (double).
static void dd_bessel_pq0(double xinv_d, MFD2 const &z, MFD2 &p, MFD2 &q) {
  if (xinv_d <= 0.25) {
    if (xinv_d <= 0.125) {
      if (xinv_d <= 0.0625) {
        p = dd_neval(z, j0_P16_IN_hi, j0_P16_IN_lo, 9) / dd_deval(z, j0_P16_ID_hi, j0_P16_ID_lo, 9);
        q = dd_neval(z, j0_Q16_IN_hi, j0_Q16_IN_lo, 10) / dd_deval(z, j0_Q16_ID_hi, j0_Q16_ID_lo, 9);
      } else {
        p = dd_neval(z, j0_P8_16N_hi, j0_P8_16N_lo, 10) / dd_deval(z, j0_P8_16D_hi, j0_P8_16D_lo, 10);
        q = dd_neval(z, j0_Q8_16N_hi, j0_Q8_16N_lo, 11) / dd_deval(z, j0_Q8_16D_hi, j0_Q8_16D_lo, 11);
      }
    } else if (xinv_d <= 0.1875) {
      p = dd_neval(z, j0_P5_8N_hi, j0_P5_8N_lo, 10) / dd_deval(z, j0_P5_8D_hi, j0_P5_8D_lo, 9);
      q = dd_neval(z, j0_Q5_8N_hi, j0_Q5_8N_lo, 10) / dd_deval(z, j0_Q5_8D_hi, j0_Q5_8D_lo, 10);
    } else {
      p = dd_neval(z, j0_P4_5N_hi, j0_P4_5N_lo, 9) / dd_deval(z, j0_P4_5D_hi, j0_P4_5D_lo, 9);
      q = dd_neval(z, j0_Q4_5N_hi, j0_Q4_5N_lo, 10) / dd_deval(z, j0_Q4_5D_hi, j0_Q4_5D_lo, 9);
    }
  } else {
    if (xinv_d <= 0.375) {
      if (xinv_d <= 0.3125) {
        p = dd_neval(z, j0_P3r2_4N_hi, j0_P3r2_4N_lo, 9) / dd_deval(z, j0_P3r2_4D_hi, j0_P3r2_4D_lo, 9);
        q = dd_neval(z, j0_Q3r2_4N_hi, j0_Q3r2_4N_lo, 10) / dd_deval(z, j0_Q3r2_4D_hi, j0_Q3r2_4D_lo, 9);
      } else {
        p = dd_neval(z, j0_P2r7_3r2N_hi, j0_P2r7_3r2N_lo, 9) / dd_deval(z, j0_P2r7_3r2D_hi, j0_P2r7_3r2D_lo, 8);
        q = dd_neval(z, j0_Q2r7_3r2N_hi, j0_Q2r7_3r2N_lo, 9) / dd_deval(z, j0_Q2r7_3r2D_hi, j0_Q2r7_3r2D_lo, 9);
      }
    } else if (xinv_d <= 0.4375) {
      p = dd_neval(z, j0_P2r3_2r7N_hi, j0_P2r3_2r7N_lo, 9) / dd_deval(z, j0_P2r3_2r7D_hi, j0_P2r3_2r7D_lo, 8);
      q = dd_neval(z, j0_Q2r3_2r7N_hi, j0_Q2r3_2r7N_lo, 9) / dd_deval(z, j0_Q2r3_2r7D_hi, j0_Q2r3_2r7D_lo, 8);
    } else {
      p = dd_neval(z, j0_P2_2r3N_hi, j0_P2_2r3N_lo, 8) / dd_deval(z, j0_P2_2r3D_hi, j0_P2_2r3D_lo, 8);
      q = dd_neval(z, j0_Q2_2r3N_hi, j0_Q2_2r3N_lo, 9) / dd_deval(z, j0_Q2_2r3D_hi, j0_Q2_2r3D_lo, 8);
    }
  }
}

// Shared asymptotic P(z),Q(z) selection for order 1 (used by both j1 and y1).
static void dd_bessel_pq1(double xinv_d, MFD2 const &z, MFD2 &p, MFD2 &q) {
  if (xinv_d <= 0.25) {
    if (xinv_d <= 0.125) {
      if (xinv_d <= 0.0625) {
        p = dd_neval(z, j1_P16_IN_hi, j1_P16_IN_lo, 9) / dd_deval(z, j1_P16_ID_hi, j1_P16_ID_lo, 9);
        q = dd_neval(z, j1_Q16_IN_hi, j1_Q16_IN_lo, 10) / dd_deval(z, j1_Q16_ID_hi, j1_Q16_ID_lo, 9);
      } else {
        p = dd_neval(z, j1_P8_16N_hi, j1_P8_16N_lo, 11) / dd_deval(z, j1_P8_16D_hi, j1_P8_16D_lo, 10);
        q = dd_neval(z, j1_Q8_16N_hi, j1_Q8_16N_lo, 11) / dd_deval(z, j1_Q8_16D_hi, j1_Q8_16D_lo, 11);
      }
    } else if (xinv_d <= 0.1875) {
      p = dd_neval(z, j1_P5_8N_hi, j1_P5_8N_lo, 10) / dd_deval(z, j1_P5_8D_hi, j1_P5_8D_lo, 10);
      q = dd_neval(z, j1_Q5_8N_hi, j1_Q5_8N_lo, 10) / dd_deval(z, j1_Q5_8D_hi, j1_Q5_8D_lo, 10);
    } else {
      p = dd_neval(z, j1_P4_5N_hi, j1_P4_5N_lo, 10) / dd_deval(z, j1_P4_5D_hi, j1_P4_5D_lo, 9);
      q = dd_neval(z, j1_Q4_5N_hi, j1_Q4_5N_lo, 10) / dd_deval(z, j1_Q4_5D_hi, j1_Q4_5D_lo, 9);
    }
  } else {
    if (xinv_d <= 0.375) {
      if (xinv_d <= 0.3125) {
        p = dd_neval(z, j1_P3r2_4N_hi, j1_P3r2_4N_lo, 9) / dd_deval(z, j1_P3r2_4D_hi, j1_P3r2_4D_lo, 9);
        q = dd_neval(z, j1_Q3r2_4N_hi, j1_Q3r2_4N_lo, 9) / dd_deval(z, j1_Q3r2_4D_hi, j1_Q3r2_4D_lo, 9);
      } else {
        p = dd_neval(z, j1_P2r7_3r2N_hi, j1_P2r7_3r2N_lo, 9) / dd_deval(z, j1_P2r7_3r2D_hi, j1_P2r7_3r2D_lo, 8);
        q = dd_neval(z, j1_Q2r7_3r2N_hi, j1_Q2r7_3r2N_lo, 9) / dd_deval(z, j1_Q2r7_3r2D_hi, j1_Q2r7_3r2D_lo, 9);
      }
    } else if (xinv_d <= 0.4375) {
      p = dd_neval(z, j1_P2r3_2r7N_hi, j1_P2r3_2r7N_lo, 9) / dd_deval(z, j1_P2r3_2r7D_hi, j1_P2r3_2r7D_lo, 8);
      q = dd_neval(z, j1_Q2r3_2r7N_hi, j1_Q2r3_2r7N_lo, 9) / dd_deval(z, j1_Q2r3_2r7D_hi, j1_Q2r3_2r7D_lo, 8);
    } else {
      p = dd_neval(z, j1_P2_2r3N_hi, j1_P2_2r3N_lo, 8) / dd_deval(z, j1_P2_2r3D_hi, j1_P2_2r3D_lo, 8);
      q = dd_neval(z, j1_Q2_2r3N_hi, j1_Q2_2r3N_lo, 9) / dd_deval(z, j1_Q2_2r3D_hi, j1_Q2_2r3D_lo, 8);
    }
  }
}

MFD2 dd_bessel_j0_full(MFD2 const &x) {
  MFD2 ax = x;
  if (ax._limbs[0] < 0.0) { ax._limbs[0] = -ax._limbs[0]; ax._limbs[1] = -ax._limbs[1]; }
  double xx = ax._limbs[0];
  if (xx == 0.0) return MFD2(1.0);
  if (!std::isfinite(xx)) return MFD2(0.0);

  if (xx <= 2.0) {
    MFD2 z = ax * ax;
    MFD2 r = z * z * dd_neval(z, j0_J0_2N_hi, j0_J0_2N_lo, 6) /
                      dd_deval(z, j0_J0_2D_hi, j0_J0_2D_lo, 6);
    return r - z * MFD2(0.25) + MFD2(1.0);
  }

  MFD2 angle = ax - dd_pair(pi_quarter_hi, pi_quarter_lo);
  MFD2 s, c;
  dd_sincos_full(angle, s, c);
  MFD2 xinv = MFD2(1.0) / ax;
  MFD2 z = xinv * xinv;
  MFD2 p, q;
  dd_bessel_pq0(1.0 / xx, z, p, q);
  p = MFD2(1.0) + z * p;
  q = (z * q - MFD2(0.125)) * xinv;
  MFD2 tpi = dd_pair(two_over_pi_hi, two_over_pi_lo);
  return multifloats::sqrt(tpi / ax) * (p * c - q * s);
}

MFD2 dd_bessel_j1_full(MFD2 const &x) {
  // j1 is odd, so (+0, -eps) must be detected as negative — checking
  // only the hi limb would miss that case.
  bool neg = x._limbs[0] < 0.0 || (x._limbs[0] == 0.0 && x._limbs[1] < 0.0);
  MFD2 ax = neg ? -x : x;
  double xx = ax._limbs[0];
  if (xx == 0.0) return MFD2(0.0);
  if (!std::isfinite(xx)) return MFD2(0.0);

  MFD2 res;
  if (xx <= 2.0) {
    MFD2 z = ax * ax;
    res = ax * MFD2(0.5) + ax * z * dd_neval(z, j1_J1_2N_hi, j1_J1_2N_lo, 6) /
                                     dd_deval(z, j1_J1_2D_hi, j1_J1_2D_lo, 6);
  } else {
    MFD2 angle = ax - dd_pair(three_pi_quarter_hi, three_pi_quarter_lo);
    MFD2 s, c;
    dd_sincos_full(angle, s, c);
    MFD2 xinv = MFD2(1.0) / ax;
    MFD2 z = xinv * xinv;
    MFD2 p, q;
    dd_bessel_pq1(1.0 / xx, z, p, q);
    p = MFD2(1.0) + z * p;
    q = (z * q + MFD2(0.375)) * xinv;
    MFD2 tpi = dd_pair(two_over_pi_hi, two_over_pi_lo);
    res = multifloats::sqrt(tpi / ax) * (p * c - q * s);
  }
  if (neg) { res._limbs[0] = -res._limbs[0]; res._limbs[1] = -res._limbs[1]; }
  return res;
}

MFD2 dd_bessel_y0_full(MFD2 const &x) {
  double xx = x._limbs[0];
  if (xx <= 0.0 || !std::isfinite(xx)) {
    MFD2 r;
    r._limbs[0] = (xx == 0.0) ? -std::numeric_limits<double>::infinity()
                               : std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }

  if (xx <= 1e-17) {
    MFD2 tpi = dd_pair(two_over_pi_hi, two_over_pi_lo);
    return dd_pair(bessel_u0_hi, bessel_u0_lo) + tpi * dd_log_full(x);
  }

  if (xx <= 2.0) {
    MFD2 z = x * x;
    MFD2 p = dd_neval(z, j0_Y0_2N_hi, j0_Y0_2N_lo, 7) /
             dd_deval(z, j0_Y0_2D_hi, j0_Y0_2D_lo, 7);
    MFD2 tpi = dd_pair(two_over_pi_hi, two_over_pi_lo);
    return tpi * dd_log_full(x) * dd_bessel_j0_full(x) + p;
  }

  MFD2 angle = x - dd_pair(pi_quarter_hi, pi_quarter_lo);
  MFD2 s, c;
  dd_sincos_full(angle, s, c);
  MFD2 xinv = MFD2(1.0) / x;
  MFD2 z = xinv * xinv;
  MFD2 p, q;
  dd_bessel_pq0(1.0 / xx, z, p, q);
  p = MFD2(1.0) + z * p;
  q = (z * q - MFD2(0.125)) * xinv;
  MFD2 tpi = dd_pair(two_over_pi_hi, two_over_pi_lo);
  return multifloats::sqrt(tpi / x) * (p * s + q * c);
}

MFD2 dd_bessel_y1_full(MFD2 const &x) {
  double xx = x._limbs[0];
  if (xx <= 0.0 || !std::isfinite(xx)) {
    MFD2 r;
    r._limbs[0] = (xx == 0.0) ? -std::numeric_limits<double>::infinity()
                               : std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }

  if (xx <= 1e-30) {
    MFD2 tpi = dd_pair(two_over_pi_hi, two_over_pi_lo);
    return MFD2(0.0) - tpi / x;
  }

  if (xx <= 2.0) {
    MFD2 z = x * x;
    MFD2 p = x * dd_neval(z, j1_Y1_2N_hi, j1_Y1_2N_lo, 7) /
                 dd_deval(z, j1_Y1_2D_hi, j1_Y1_2D_lo, 7);
    MFD2 tpi = dd_pair(two_over_pi_hi, two_over_pi_lo);
    return tpi * (dd_log_full(x) * dd_bessel_j1_full(x) - MFD2(1.0) / x) + p;
  }

  MFD2 angle = x - dd_pair(three_pi_quarter_hi, three_pi_quarter_lo);
  MFD2 s, c;
  dd_sincos_full(angle, s, c);
  MFD2 xinv = MFD2(1.0) / x;
  MFD2 z = xinv * xinv;
  MFD2 p, q;
  dd_bessel_pq1(1.0 / xx, z, p, q);
  p = MFD2(1.0) + z * p;
  q = (z * q + MFD2(0.375)) * xinv;
  MFD2 tpi = dd_pair(two_over_pi_hi, two_over_pi_lo);
  return multifloats::sqrt(tpi / x) * (p * s + q * c);
}

} // anonymous namespace

// =============================================================================
// C-ABI entry points — extern "C" functions following math.h naming convention.
// These are the canonical DD implementations; the C++ template wrappers in
// multifloats.hh and the Fortran bind(C) interfaces both call these.
// =============================================================================

// multifloats_c.h is already pulled in via multifloats.hh.

namespace {
namespace mf = multifloats;
using MF2 = mf::MultiFloat<double, 2>;
static inline MF2 from(dd_t x) { MF2 r; r._limbs[0] = x.hi; r._limbs[1] = x.lo; return r; }
static inline dd_t to(MF2 const &x) { return {x._limbs[0], x._limbs[1]}; }

// Inline matmul micro-ops used by the dispatched panel templates below.
// `always_inline` is a hint: the panel templates rely on both ends of
// the MAC fusing into one register-resident body; a stray function call
// here would spill the accumulator ladder and cost ~2× the hot loop.
__attribute__((always_inline))
static inline void dd_mac_inl(double ah, double al, double bh, double bl,
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
static inline void dd_renorm_inl(double &s_hi, double &s_lo) {
  double t = s_hi + s_lo;
  double bp = t - s_hi;
  double ap = t - bp;
  double aerr = s_hi - ap;
  double berr = s_lo - bp;
  s_lo = aerr + berr;
  s_hi = t;
}

__attribute__((always_inline))
static inline dd_t dd_finalize_inl(double s_hi, double s_lo) {
  dd_renorm_inl(s_hi, s_lo);
  return {s_hi, s_lo};
}

// AXPY-style panel µkernel. Processes a row-panel of A of compile-time
// height MR against a contiguous x / b-column, writing MR finalized DD
// outputs. `lda` is the leading dim of A (so the panel can sit in a
// larger matrix); `renorm_interval > 0` triggers a two_sum of each
// accumulator pair every `renorm_interval` reductions (keeps s_lo
// bounded and precision ~DD for large k, matching dot_product).
template <int MR>
static inline void dd_gaxpy_mv_panel(const dd_t *__restrict__ a,
                                     const dd_t *__restrict__ x,
                                     dd_t *__restrict__ y, int64_t lda,
                                     int64_t k, int64_t renorm_interval) {
  double s_hi[MR] = {}, s_lo[MR] = {};
  // Fast path: if no intermediate renorm is needed, drop the chunking
  // scaffolding entirely — the p-loop is a single tight sweep.
  if (renorm_interval <= 0 || k <= renorm_interval) {
    for (int64_t p = 0; p < k; ++p) {
      const double xh = x[p].hi;
      const double xl = x[p].lo;
      const dd_t *__restrict__ acol = a + p * lda;
      for (int i = 0; i < MR; ++i) {
        dd_mac_inl(acol[i].hi, acol[i].lo, xh, xl, s_hi[i], s_lo[i]);
      }
    }
    for (int i = 0; i < MR; ++i) y[i] = dd_finalize_inl(s_hi[i], s_lo[i]);
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
      const dd_t *__restrict__ acol = a + p * lda;
      for (int i = 0; i < MR; ++i) {
        dd_mac_inl(acol[i].hi, acol[i].lo, xh, xl, s_hi[i], s_lo[i]);
      }
    }
    p0 = pend;
    if (p0 < k) {
      for (int i = 0; i < MR; ++i) dd_renorm_inl(s_hi[i], s_lo[i]);
    }
  }
  for (int i = 0; i < MR; ++i) y[i] = dd_finalize_inl(s_hi[i], s_lo[i]);
}

// Cold tail handler: 1..7 leftover rows. Kept out-of-line so the hot
// dispatcher stays small enough for the compiler to inline panel<8>
// into dd_matmul_mv — otherwise the 7 extra template instantiations
// would bloat the caller and push the accumulators out of registers.
__attribute__((noinline))
static void dd_gaxpy_mv_tail(const dd_t *__restrict__ a,
                             const dd_t *__restrict__ x,
                             dd_t *__restrict__ y, int tail,
                             int64_t lda, int64_t k,
                             int64_t renorm_interval) {
  switch (tail) {
    case 1: dd_gaxpy_mv_panel<1>(a, x, y, lda, k, renorm_interval); break;
    case 2: dd_gaxpy_mv_panel<2>(a, x, y, lda, k, renorm_interval); break;
    case 3: dd_gaxpy_mv_panel<3>(a, x, y, lda, k, renorm_interval); break;
    case 4: dd_gaxpy_mv_panel<4>(a, x, y, lda, k, renorm_interval); break;
    case 5: dd_gaxpy_mv_panel<5>(a, x, y, lda, k, renorm_interval); break;
    case 6: dd_gaxpy_mv_panel<6>(a, x, y, lda, k, renorm_interval); break;
    case 7: dd_gaxpy_mv_panel<7>(a, x, y, lda, k, renorm_interval); break;
  }
}

// Tile any m into MR=8 row panels plus a 1..7 tail. Keep this small so
// the hot panel<8> body inlines — accumulators need to stay in registers.
static inline void dd_gaxpy_mv_dispatch(const dd_t *__restrict__ a,
                                        const dd_t *__restrict__ x,
                                        dd_t *__restrict__ y, int64_t m,
                                        int64_t k, int64_t lda,
                                        int64_t renorm_interval) {
  constexpr int MR = 8;
  int64_t i = 0;
  for (; i + MR <= m; i += MR) {
    dd_gaxpy_mv_panel<MR>(a + i, x, y + i, lda, k, renorm_interval);
  }
  int tail = static_cast<int>(m - i);
  if (tail > 0) dd_gaxpy_mv_tail(a + i, x, y + i, tail, lda, k, renorm_interval);
}

// MR×NR GEMM µkernel: loads one column-slice of A per p and reuses it
// across NR B-column entries. Amortizes A bandwidth vs a per-column mv
// dispatch at the cost of MR*NR accumulator pairs (stack-spilled on
// machines with < 32 FP regs, but the spill cost is paid once per p
// and dwarfed by the DD-mac FLOPs).
template <int MR, int NR>
static inline void dd_gemm_panel(const dd_t *__restrict__ a,
                                 const dd_t *__restrict__ b,
                                 dd_t *__restrict__ c,
                                 int64_t lda, int64_t ldb, int64_t ldc,
                                 int64_t k, int64_t renorm_interval) {
  double s_hi[MR][NR] = {}, s_lo[MR][NR] = {};
  auto kernel_block = [&](int64_t p0, int64_t pend) {
    for (int64_t p = p0; p < pend; ++p) {
      const dd_t *__restrict__ acol = a + p * lda;
      double ah[MR], al[MR];
      for (int i = 0; i < MR; ++i) { ah[i] = acol[i].hi; al[i] = acol[i].lo; }
      for (int jj = 0; jj < NR; ++jj) {
        const double bh = b[jj * ldb + p].hi;
        const double bl = b[jj * ldb + p].lo;
        for (int i = 0; i < MR; ++i) {
          dd_mac_inl(ah[i], al[i], bh, bl, s_hi[i][jj], s_lo[i][jj]);
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
            dd_renorm_inl(s_hi[i][jj], s_lo[i][jj]);
      }
    }
  }
  for (int jj = 0; jj < NR; ++jj)
    for (int i = 0; i < MR; ++i)
      c[jj * ldc + i] = dd_finalize_inl(s_hi[i][jj], s_lo[i][jj]);
}
} // anonymous namespace

extern "C" {

// Arithmetic
dd_t dd_add(dd_t a, dd_t b) { return to(from(a) + from(b)); }
dd_t dd_sub(dd_t a, dd_t b) { return to(from(a) - from(b)); }
dd_t dd_mul(dd_t a, dd_t b) { return to(from(a) * from(b)); }
dd_t dd_div(dd_t a, dd_t b) { return to(from(a) / from(b)); }
dd_t dd_neg(dd_t a)  { return to(-from(a)); }
dd_t dd_abs(dd_t a)  { return to(multifloats::abs(from(a))); }
dd_t dd_sqrt(dd_t a) { return to(multifloats::sqrt(from(a))); }

// Rounding (Fortran AINT / ANINT delegate here)
dd_t dd_trunc(dd_t a) { return to(multifloats::trunc(from(a))); }
dd_t dd_round(dd_t a) { return to(multifloats::round(from(a))); }

// Binary
dd_t dd_fmin(dd_t a, dd_t b)     { return to(multifloats::fmin(from(a), from(b))); }
dd_t dd_fmax(dd_t a, dd_t b)     { return to(multifloats::fmax(from(a), from(b))); }
dd_t dd_hypot(dd_t a, dd_t b)    { return to(multifloats::hypot(from(a), from(b))); }
dd_t dd_pow(dd_t a, dd_t b)      { return to(dd_pow_full(from(a), from(b))); }
dd_t dd_fmod(dd_t a, dd_t b)     { return to(multifloats::fmod(from(a), from(b))); }
dd_t dd_fdim(dd_t a, dd_t b)     { return to(multifloats::fdim(from(a), from(b))); }
dd_t dd_copysign(dd_t a, dd_t b) { return to(multifloats::copysign(from(a), from(b))); }
dd_t dd_fma(dd_t a, dd_t b, dd_t c) { return to(multifloats::fma(from(a), from(b), from(c))); }

// Exponential / logarithmic
dd_t dd_exp(dd_t a)   { return to(dd_exp_full(from(a))); }
dd_t dd_exp2(dd_t a)  { return to(dd_exp2_full(from(a))); }
dd_t dd_log(dd_t a)   { return to(dd_log_full(from(a))); }
dd_t dd_log2(dd_t a)  { return to(dd_log2_full(from(a))); }
dd_t dd_log10(dd_t a) { return to(dd_log10_full(from(a))); }

// Trigonometric
dd_t dd_sin(dd_t a)   { return to(dd_sin_full(from(a))); }
dd_t dd_cos(dd_t a)   { return to(dd_cos_full(from(a))); }
dd_t dd_tan(dd_t a)   { return to(dd_tan_full(from(a))); }
dd_t dd_asin(dd_t a)  { return to(dd_asin_full(from(a))); }
dd_t dd_acos(dd_t a)  { return to(dd_acos_full(from(a))); }
dd_t dd_atan(dd_t a)  { return to(dd_atan_full(from(a))); }
dd_t dd_atan2(dd_t a, dd_t b) { return to(dd_atan2_full(from(a), from(b))); }

// π-scaled trig: fn_pi(x) = fn(π·x) or fn(x)/π.
dd_t dd_sinpi(dd_t a)  { return to(dd_sinpi_full(from(a))); }
dd_t dd_cospi(dd_t a)  { return to(dd_cospi_full(from(a))); }
dd_t dd_tanpi(dd_t a)  { return to(dd_tanpi_full(from(a))); }
dd_t dd_asinpi(dd_t a) { return to(dd_asinpi_full(from(a))); }
dd_t dd_acospi(dd_t a) { return to(dd_acospi_full(from(a))); }
dd_t dd_atanpi(dd_t a) { return to(dd_atanpi_full(from(a))); }
dd_t dd_atan2pi(dd_t a, dd_t b) { return to(dd_atan2pi_full(from(a), from(b))); }

// Hyperbolic
dd_t dd_sinh(dd_t a)  { return to(dd_sinh_full(from(a))); }
dd_t dd_cosh(dd_t a)  { return to(dd_cosh_full(from(a))); }
dd_t dd_tanh(dd_t a)  { return to(dd_tanh_full(from(a))); }
dd_t dd_asinh(dd_t a) { return to(dd_asinh_full(from(a))); }
dd_t dd_acosh(dd_t a) { return to(dd_acosh_full(from(a))); }
dd_t dd_atanh(dd_t a) { return to(dd_atanh_full(from(a))); }

// Error functions
dd_t dd_erf(dd_t a)   { return to(dd_erf_full(from(a))); }
dd_t dd_erfc(dd_t a)  { return to(dd_erfc_full(from(a))); }
dd_t dd_erfcx(dd_t a) { return to(dd_erfc_scaled_full(from(a))); }

// Gamma functions
dd_t dd_tgamma(dd_t a) { return to(dd_tgamma_full(from(a))); }
dd_t dd_lgamma(dd_t a) { return to(dd_lgamma_full(from(a))); }

// Bessel functions (math.h naming: j0, j1, y0, y1)
dd_t dd_j0(dd_t a) { return to(dd_bessel_j0_full(from(a))); }
dd_t dd_j1(dd_t a) { return to(dd_bessel_j1_full(from(a))); }
dd_t dd_y0(dd_t a) { return to(dd_bessel_y0_full(from(a))); }
dd_t dd_y1(dd_t a) { return to(dd_bessel_y1_full(from(a))); }

// Fused sincos / sinhcosh. Out-pointer style (C has no multi-value return).
void dd_sincos(dd_t a, dd_t *s, dd_t *c) {
  MF2 ss, cc;
  dd_sincos_full(from(a), ss, cc);
  *s = to(ss);
  *c = to(cc);
}

void dd_sinhcosh(dd_t a, dd_t *s, dd_t *c) {
  MF2 ss, cc;
  dd_sinhcosh_full(from(a), ss, cc);
  *s = to(ss);
  *c = to(cc);
}

// ---- Complex DD transcendentals ------------------------------------------
//
// Branch cuts follow C99 Annex G (matching libquadmath cexpq/clogq/csqrtq).
// Where the real-axis formula needs both sin(y) and cos(y) (or both sinh
// and cosh), these call the fused dd_sincos_full / dd_sinhcosh_full kernels
// so one range reduction + Taylor pair covers both outputs. Functions that
// are built on top of cx_sqrt / cx_log (asin/acos/...) call the exported
// dd_cx_* symbols directly; within a single TU the compiler inlines those.
//
// The text "2 fused" vs "4 separate" in speed comments below is measured
// against libstdc++'s generic <complex> template path (see commit notes).

// exp(a+bi) = e^a · (cos b + i·sin b). 2 transcendentals (exp + sincos).
cdd_t dd_cx_exp(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  MF2 ea = dd_exp_full(a);
  MF2 s, c;
  dd_sincos_full(b, s, c);
  return { to(ea * c), to(ea * s) };
}

// log(a+bi) = log(|z|) + i·atan2(b, a). Compute log(|z|) as
//   log(big) + 0.5·log(1 + (small/big)^2)
// with big = max(|a|,|b|), small = min(|a|,|b|). This avoids the
// overflow that hypot(a,b) suffers when |z| exceeds DBL_MAX even
// though log(|z|) itself is finite. atan2 handles the negative-real-
// axis branch cut (including signed-zero propagation).
cdd_t dd_cx_log(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  MF2 ax = multifloats::abs(a);
  MF2 ay = multifloats::abs(b);
  MF2 big   = (ax > ay) ? ax : ay;
  MF2 small = (ax > ay) ? ay : ax;
  MF2 r;
  if (big._limbs[0] == 0.0) {
    r._limbs[0] = -std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
  } else {
    MF2 ratio = small / big;
    MF2 one = dd_pair(1.0, 0.0);
    MF2 half = dd_pair(0.5, 0.0);
    r = dd_log_full(big) + half * dd_log_full(one + ratio * ratio);
  }
  MF2 phi = dd_atan2_full(b, a);
  return { to(r), to(phi) };
}

// log10(z) = log(z) · (1/ln 10).
cdd_t dd_cx_log10(cdd_t z) {
  static const MF2 inv_ln10 = dd_pair(0x1.bcb7b1526e50ep-2,
                                      0x1.95355baaafad3p-57);
  cdd_t l = dd_cx_log(z);
  return { to(from(l.re) * inv_ln10), to(from(l.im) * inv_ln10) };
}

// pow(z, w) = exp(w · log(z)). C99 G.6.4.1; 0^w folds to 0.
cdd_t dd_cx_pow(cdd_t z, cdd_t w) {
  if (z.re.hi == 0.0 && z.im.hi == 0.0) return { {0.0, 0.0}, {0.0, 0.0} };
  cdd_t l = dd_cx_log(z);
  MF2 lr = from(l.re), li = from(l.im);
  MF2 wr = from(w.re), wi = from(w.im);
  // w * log(z) = (wr·lr − wi·li) + i·(wr·li + wi·lr)
  cdd_t p = { to(wr * lr - wi * li), to(wr * li + wi * lr) };
  return dd_cx_exp(p);
}

// sqrt(z). Principal branch; cut on negative real axis, continuous above.
// Uses Moshier's 2·Re·Im = Im identity to avoid cancellation in mod ± a.
cdd_t dd_cx_sqrt(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  if (a._limbs[0] == 0.0 && b._limbs[0] == 0.0)
    return { to(MF2(0.0)), z.im };  // preserves signed zero of imag
  MF2 mod = multifloats::hypot(a, b);
  MF2 half = dd_pair(0.5, 0.0);
  if (a._limbs[0] >= 0.0) {
    MF2 r = multifloats::sqrt((mod + a) * half);
    MF2 i = b / (r + r);
    return { to(r), to(i) };
  } else {
    MF2 s = multifloats::sqrt((mod - a) * half);
    MF2 r = multifloats::abs(b) / (s + s);
    MF2 i = (b._limbs[0] < 0.0) ? -s : s;
    return { to(r), to(i) };
  }
}

// sin(a+bi) = sin(a)·cosh(b) + i·cos(a)·sinh(b). 2 fused transcendentals.
cdd_t dd_cx_sin(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  MF2 sa, ca, sb, cb;
  dd_sincos_full(a, sa, ca);
  dd_sinhcosh_full(b, sb, cb);
  return { to(sa * cb), to(ca * sb) };
}

// cos(a+bi) = cos(a)·cosh(b) − i·sin(a)·sinh(b). 2 fused transcendentals.
cdd_t dd_cx_cos(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  MF2 sa, ca, sb, cb;
  dd_sincos_full(a, sa, ca);
  dd_sinhcosh_full(b, sb, cb);
  return { to(ca * cb), to(-(sa * sb)) };
}

// tan(a+bi) = (sin(a)·cos(a) + i·sinh(b)·cosh(b)) / (cos(a)² + sinh(b)²).
// libquadmath ctanq formula; 2 fused transcendentals + real div
// (vs 8 for generic sin(z)/cos(z)).
cdd_t dd_cx_tan(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  MF2 sa, ca, sb, cb;
  dd_sincos_full(a, sa, ca);
  dd_sinhcosh_full(b, sb, cb);
  MF2 den = ca * ca + sb * sb;
  return { to((sa * ca) / den), to((sb * cb) / den) };
}

// asin(z) = −i · log(i·z + sqrt(1 − z²)). Cuts on real axis for |a|>1.
cdd_t dd_cx_asin(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  // 1 − z² : Re = 1 − a² + b², Im = −2ab
  cdd_t one_mz2 = { to(MF2(1.0) - a * a + b * b), to(-(a * b + a * b)) };
  cdd_t root = dd_cx_sqrt(one_mz2);
  MF2 rr = from(root.re), ri = from(root.im);
  // i·z + root : (−b + rr) + i·(a + ri)
  cdd_t arg = { to(-b + rr), to(a + ri) };
  cdd_t l = dd_cx_log(arg);
  // −i · (lr + i·li) = li − i·lr
  return { l.im, to(-from(l.re)) };
}

// acos(z) = π/2 − asin(z). Uses DD-precision π/2; the generic <complex>
// template truncates `(_Tp)1.5707963…L` to double on arm64-macOS, losing
// the DD low limb — specializing here fixes that correctness bug.
cdd_t dd_cx_acos(cdd_t z) {
  static const MF2 half_pi = dd_pair(0x1.921fb54442d18p+0,
                                     0x1.1a62633145c07p-54);
  cdd_t s = dd_cx_asin(z);
  return { to(half_pi - from(s.re)), to(-from(s.im)) };
}

// atan(z). libstdc++ generic form:
//   Re = 0.5 · atan2(2a, 1 − a² − b²)
//   Im = 0.25 · log((a² + (b+1)²) / (a² + (b−1)²))
// Single log-ratio rather than the naive (i/2)(log(1−iz) − log(1+iz)).
cdd_t dd_cx_atan(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  MF2 a2 = a * a;
  MF2 x = MF2(1.0) - a2 - b * b;
  MF2 bp1 = b + MF2(1.0), bm1 = b - MF2(1.0);
  MF2 num = a2 + bp1 * bp1;
  MF2 den = a2 + bm1 * bm1;
  MF2 half = dd_pair(0.5, 0.0), quarter = dd_pair(0.25, 0.0);
  return { to(half * dd_atan2_full(a + a, x)),
           to(quarter * dd_log_full(num / den)) };
}

// sinh(a+bi) = sinh(a)·cos(b) + i·cosh(a)·sin(b). 2 fused transcendentals.
cdd_t dd_cx_sinh(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  MF2 sa, ca, sb, cb;
  dd_sinhcosh_full(a, sa, ca);
  dd_sincos_full(b, sb, cb);
  return { to(sa * cb), to(ca * sb) };
}

// cosh(a+bi) = cosh(a)·cos(b) + i·sinh(a)·sin(b). 2 fused transcendentals.
cdd_t dd_cx_cosh(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  MF2 sa, ca, sb, cb;
  dd_sinhcosh_full(a, sa, ca);
  dd_sincos_full(b, sb, cb);
  return { to(ca * cb), to(sa * sb) };
}

// tanh(a+bi) = (sinh(a)·cosh(a) + i·sin(b)·cos(b)) / (sinh(a)² + cos(b)²).
// libquadmath ctanhq formula.
cdd_t dd_cx_tanh(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  MF2 sa, ca, sb, cb;
  dd_sinhcosh_full(a, sa, ca);
  dd_sincos_full(b, sb, cb);
  MF2 den = sa * sa + cb * cb;
  return { to((sa * ca) / den), to((sb * cb) / den) };
}

// asinh(z) = log(z + sqrt(z² + 1)). Textbook; matches libstdc++ generic.
cdd_t dd_cx_asinh(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  // 1 + z² : Re = 1 + a² − b², Im = 2ab
  cdd_t one_pz2 = { to(MF2(1.0) + a * a - b * b), to(a * b + a * b) };
  cdd_t root = dd_cx_sqrt(one_pz2);
  cdd_t arg = { to(a + from(root.re)), to(b + from(root.im)) };
  return dd_cx_log(arg);
}

// acosh(z). Kahan's formula: 2·log(sqrt((z+1)/2) + sqrt((z−1)/2)).
// More stable near z≈1 than the naive log(z + sqrt(z²−1)) which suffers
// cancellation in z²−1.
cdd_t dd_cx_acosh(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  MF2 half = dd_pair(0.5, 0.0);
  cdd_t zp = { to((a + MF2(1.0)) * half), to(b * half) };
  cdd_t zm = { to((a - MF2(1.0)) * half), to(b * half) };
  cdd_t s1 = dd_cx_sqrt(zp);
  cdd_t s2 = dd_cx_sqrt(zm);
  cdd_t sum = { to(from(s1.re) + from(s2.re)), to(from(s1.im) + from(s2.im)) };
  cdd_t l = dd_cx_log(sum);
  return { to(from(l.re) + from(l.re)), to(from(l.im) + from(l.im)) };
}

// atanh(z). libstdc++ generic form (atan with a, b swapped):
//   Re = 0.25 · log((1+a)² + b²)/((1−a)² + b²))
//   Im = 0.5 · atan2(2b, 1 − a² − b²)
cdd_t dd_cx_atanh(cdd_t z) {
  MF2 a = from(z.re), b = from(z.im);
  MF2 b2 = b * b;
  MF2 x = MF2(1.0) - b2 - a * a;
  MF2 ap1 = a + MF2(1.0), am1 = a - MF2(1.0);
  MF2 num = b2 + ap1 * ap1;
  MF2 den = b2 + am1 * am1;
  MF2 half = dd_pair(0.5, 0.0), quarter = dd_pair(0.25, 0.0);
  return { to(quarter * dd_log_full(num / den)),
           to(half * dd_atan2_full(b + b, x)) };
}

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

void dd_matmul_mm(const dd_t *__restrict__ a, const dd_t *__restrict__ b,
                  dd_t *__restrict__ c, int64_t m, int64_t k, int64_t n,
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
      dd_gemm_panel<MR, NR>(a + i, b + j * k, c + j * m + i,
                            m, k, m, k, renorm_interval);
    }
    int tail_m = static_cast<int>(m - i);
    if (tail_m > 0) {
      for (int jj = 0; jj < NR; ++jj) {
        dd_gaxpy_mv_tail(a + i, b + (j + jj) * k, c + (j + jj) * m + i,
                         tail_m, m, k, renorm_interval);
      }
    }
  }
  for (; j < n; ++j) {
    dd_gaxpy_mv_dispatch(a, b + j * k, c + j * m, m, k, m, renorm_interval);
  }
}

void dd_matmul_mv(const dd_t *__restrict__ a, const dd_t *__restrict__ x,
                  dd_t *__restrict__ y, int64_t m, int64_t k,
                  int64_t renorm_interval) {
  dd_gaxpy_mv_dispatch(a, x, y, m, k, m, renorm_interval);
}

// vm: y[j] = sum_p x[p] * B[p, j]. Column-major B makes B[:, j]
// contiguous at fixed j, so one scalar accumulator per output is optimal.
void dd_matmul_vm(const dd_t *__restrict__ x, const dd_t *__restrict__ b,
                  dd_t *__restrict__ y, int64_t k, int64_t n,
                  int64_t renorm_interval) {
  const bool simple = (renorm_interval <= 0) || (k <= renorm_interval);
  for (int64_t j = 0; j < n; ++j) {
    double s_hi = 0.0, s_lo = 0.0;
    const dd_t *__restrict__ bcol = b + j * k;
    if (simple) {
      for (int64_t p = 0; p < k; ++p) {
        dd_mac_inl(x[p].hi, x[p].lo, bcol[p].hi, bcol[p].lo, s_hi, s_lo);
      }
    } else {
      const int64_t chunk = renorm_interval;
      int64_t p0 = 0;
      while (p0 < k) {
        int64_t pend = p0 + chunk;
        if (pend > k) pend = k;
        for (int64_t p = p0; p < pend; ++p) {
          dd_mac_inl(x[p].hi, x[p].lo, bcol[p].hi, bcol[p].lo, s_hi, s_lo);
        }
        p0 = pend;
        if (p0 < k) dd_renorm_inl(s_hi, s_lo);
      }
    }
    y[j] = dd_finalize_inl(s_hi, s_lo);
  }
}

// Comparison
int dd_eq(dd_t a, dd_t b) { return from(a) == from(b) ? 1 : 0; }
int dd_ne(dd_t a, dd_t b) { return from(a) != from(b) ? 1 : 0; }
int dd_lt(dd_t a, dd_t b) { return from(a) <  from(b) ? 1 : 0; }
int dd_le(dd_t a, dd_t b) { return from(a) <= from(b) ? 1 : 0; }
int dd_gt(dd_t a, dd_t b) { return from(a) >  from(b) ? 1 : 0; }
int dd_ge(dd_t a, dd_t b) { return from(a) >= from(b) ? 1 : 0; }

} // extern "C"
