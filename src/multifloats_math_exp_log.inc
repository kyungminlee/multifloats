// ---- exp2 / exp / log2 / log / log10 ---------------------------------------

namespace {
// exp / exp2 non-finite handling: +inf → +inf, -inf → 0, NaN → NaN, lo=0.
inline float64x2 _exp_nonfinite(double hi) {
  float64x2 r;
  if (hi > 0.0) r._limbs[0] = hi;        // +inf
  else if (hi < 0.0) r._limbs[0] = 0.0;  // -inf → 0
  else r._limbs[0] = hi;                 // NaN
  r._limbs[1] = 0.0;
  return r;
}
} // namespace

// Compute 2^(n_d + rh + rl), given |rh + rl| ≲ 0.5 and n_d an integer-valued
// double. The fractional part is further split as k/256 + s with k ∈
// [-128, 128] and |s| ≲ 1/512 from a 257-entry table of 2^(k/256). A
// degree-9 Horner on s evaluates 2^s − 1 (next term ~6e-36 ≪ ulp_dd), and
// one DD multiply against the table entry recombines. No cube, no 8×
// amplification; round-off floor is ~1 ulp_dd.
float64x2 exp2_from_reduced(double n_d, double rh, double rl) {
  if (n_d > 1100.0) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
    return r;
  }
  if (n_d < -1100.0) return float64x2();

  // Table index: k = round(rh · 256). rh − k/256 is exact (Sterbenz).
  double k_d = std::nearbyint(rh * 256.0);
  int k = static_cast<int>(k_d);
  if (k > 128) k = 128;
  if (k < -128) k = -128;
  double sh = rh - static_cast<double>(k) / 256.0;
  double sl = rl;
  fast_two_sum(sh, sl, sh, sl);
  float64x2 s_dd; s_dd._limbs[0] = sh; s_dd._limbs[1] = sl;

  // Horner on c1..c9 of exp2_coefs: q = c1 + s·(c2 + s·(… + s·c9)).
  float64x2 q; q._limbs[0] = exp2_coefs_hi[9]; q._limbs[1] = exp2_coefs_lo[9];
  for (int j = 8; j >= 1; --j) {
    q = q * s_dd;
    float64x2 c; c._limbs[0] = exp2_coefs_hi[j]; c._limbs[1] = exp2_coefs_lo[j];
    q = q + c;
  }
  // 2^s − 1 = s · q; 2^s = 1 + s·q.
  float64x2 pow2s = float64x2(1.0) + s_dd * q;

  // Mantissa = T[k+128] · 2^s.
  float64x2 T;
  T._limbs[0] = exp2_table_hi[k + 128];
  T._limbs[1] = exp2_table_lo[k + 128];
  float64x2 mantissa = T * pow2s;

  int n = static_cast<int>(n_d);
  int half_n = n / 2;
  float64x2 r;
  r._limbs[0] = std::ldexp(mantissa._limbs[0], half_n);
  r._limbs[1] = std::ldexp(mantissa._limbs[1], half_n);
  r._limbs[0] = std::ldexp(r._limbs[0], n - half_n);
  r._limbs[1] = std::ldexp(r._limbs[1], n - half_n);
  return r;
}

float64x2 exp2_kernel(float64x2 const &x) {
  double n_d = std::nearbyint(x._limbs[0]);
  double rh, rl, t;
  two_sum(x._limbs[0], -n_d, rh, t);
  rl = t + x._limbs[1];
  fast_two_sum(rh, rl, rh, rl);
  return exp2_from_reduced(n_d, rh, rl);
}

float64x2 exp2_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) return _exp_nonfinite(x._limbs[0]);
  if (x._limbs[0] < exp2_min) return float64x2();
  if (x._limbs[0] > exp2_max) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
    return r;
  }
  return exp2_kernel(x);
}

// ---- TD-output exp variants -----------------------------------------------
// Mirror of `exp2_from_reduced` / `exp_full` that operates on TD (~159 bit)
// internals throughout — every constant (log2_e, exp2_coefs, exp2_table,
// ln2_cw1..cw4) is lifted to TD precision from dd_constants.hh so the
// output genuinely carries ~158 bits rather than being capped at the
// DD-constant floor. See doc/developer/TRIPLE_DOUBLE.md §2
// for why TD arithmetic over DD constants stays capped at 2⁻¹⁰⁵.
//
// `exp_full_td(x) |> td_to_dd` matches `exp_full(x)` within one DD ulp.
// The TD pair (em, el) carries the residue used by the Re path of
// `cexpm1dd` to survive the cancellation surface cos(b)·e^a = 1.

float64x3 exp2_td_from_reduced(double n_d, float64x3 const &y) {
  if (n_d > 1100.0) {
    return {std::numeric_limits<double>::infinity(), 0.0, 0.0};
  }
  if (n_d < -1100.0) return {0.0, 0.0, 0.0};

  // Same table-index split as the DD path, but extended to TD residue.
  // y.hi − k/256 is exact (Sterbenz) provided |y.hi| ≲ 0.5 and
  // k ∈ [-128, 128].
  double k_d = std::nearbyint(y._limbs[0] * 256.0);
  int k = static_cast<int>(k_d);
  if (k > 128) k = 128;
  if (k < -128) k = -128;
  double sh = y._limbs[0] - static_cast<double>(k) / 256.0;
  float64x3 s_td = {sh, y._limbs[1], y._limbs[2]};
  // Canonicalize: the three limbs may no longer be ordered after the
  // subtraction (sh can collapse to ~y.mid if k rounded y.hi well).
  s_td = renorm3(s_td._limbs[0], s_td._limbs[1], s_td._limbs[2]);

  // TD Horner on c[1]..c[15] of exp2_coefs. Degree 15 drives the Taylor
  // truncation below ~2⁻¹⁹⁷, well under the TD floor.
  float64x3 q_td;
  q_td._limbs[0] = exp2_coefs_hi[15];
  q_td._limbs[1] = exp2_coefs_lo[15];
  q_td._limbs[2] = exp2_coefs_lo2[15];
  for (int j = 14; j >= 1; --j) {
    q_td = td_mul_td(q_td, s_td);
    float64x3 c_td = {exp2_coefs_hi[j], exp2_coefs_lo[j],
                      exp2_coefs_lo2[j]};
    q_td = td_add_td(q_td, c_td);
  }
  // 2^s − 1 = s · q; 2^s = 1 + s·q.
  float64x3 sq_td = td_mul_td(q_td, s_td);
  float64x3 pow2s_td = td_add_double(sq_td, 1.0);

  // Mantissa = T[k+128] · 2^s with the full TD table entry.
  float64x3 T_td = {exp2_table_hi[k + 128], exp2_table_lo[k + 128],
                    exp2_table_lo2[k + 128]};
  float64x3 mantissa_td = td_mul_td(pow2s_td, T_td);

  // Scale by 2^n. ldexp on each limb is exact for normal inputs; split
  // n in halves so ldexp stays inside its well-defined range at the
  // extreme ends of exp's domain.
  int n = static_cast<int>(n_d);
  int half_n = n / 2;
  float64x3 r;
  for (int i = 0; i < 3; ++i) {
    double v = std::ldexp(mantissa_td._limbs[i], half_n);
    r._limbs[i] = std::ldexp(v, n - half_n);
  }
  return r;
}

float64x3 exp_full_td(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    return {_exp_nonfinite(x._limbs[0])._limbs[0], 0.0, 0.0};
  }
  if (x._limbs[0] > exp_max) {
    return {std::numeric_limits<double>::infinity(), 0.0, 0.0};
  }
  if (x._limbs[0] < exp_min) return {0.0, 0.0, 0.0};

  double n_d = std::nearbyint(x._limbs[0] * log2_e_hi);

  // TD-precision Cody-Waite reduction: r = x − n·ln2 via the 4-part
  // split cw1..cw4. cw1 is 33-bit (n·cw1 exact in double for |n| < 2^20),
  // cw2..cw4 are 53-bit residuals totalling ~192 bits below cw1. Exact
  // products n·cwk for k≥2 are captured via FMA pairs; everything is
  // accumulated through the TSUM primitive into a TD residue.
  double T0 = 0.0, T1 = 0.0, T2 = 0.0;
  tsum3(T0, T1, T2, x._limbs[0]);
  tsum3(T0, T1, T2, x._limbs[1]);
  tsum3(T0, T1, T2, -(n_d * ln2_cw1));           // exact
  double p2h = n_d * ln2_cw2;
  double p2l = std::fma(n_d, ln2_cw2, -p2h);     // n·cw2 exact split
  tsum3(T0, T1, T2, -p2h);
  tsum3(T0, T1, T2, -p2l);
  double p3h = n_d * ln2_cw3;
  double p3l = std::fma(n_d, ln2_cw3, -p3h);
  tsum3(T0, T1, T2, -p3h);
  tsum3(T0, T1, T2, -p3l);
  tsum3(T0, T1, T2, -(n_d * ln2_cw4));           // below TD floor already
  tsum3_finalize(T0, T1, T2);
  float64x3 rr_td = {T0, T1, T2};

  // Convert r (natural-log form) to base-2: y = r · log2(e), both TD.
  // log2(e) carried to ~159 bits via log2_e_{hi,lo,lo2}.
  float64x3 log2e_td = {log2_e_hi, log2_e_lo, log2_e_lo2};
  float64x3 y_td = td_mul_td(rr_td, log2e_td);
  return exp2_td_from_reduced(n_d, y_td);
}

float64x2 exp_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) return _exp_nonfinite(x._limbs[0]);
  if (x._limbs[0] > exp_max) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] < exp_min) return float64x2();

  // Cody-Waite range reduction: r = x - n*ln2 via three-piece split.
  // n chosen from x/ln2 rounded to nearest integer; |r| <= ln2/2 ~ 0.347.
  double n_d = std::nearbyint(x._limbs[0] * log2_e_hi);

  // r = x - n*ln2 via three-piece CW subtraction, inlined as DD - double
  // steps to avoid paying for full DD-DD subtracts against zero-lo operands.
  // n_d*cw1 is exact in double (cw1 has 33-bit significand, |n_d| < 2^11).
  double rh, rl, t;
  two_sum(x._limbs[0], -(n_d * ln2_cw1), rh, t);
  rl = t + x._limbs[1];
  fast_two_sum(rh, rl, rh, rl);
  // -= n*cw2, capturing the doubled-precision product via FMA.
  double p_hi = n_d * ln2_cw2;
  double p_lo = std::fma(n_d, ln2_cw2, -p_hi);
  two_sum(rh, -p_hi, rh, t);
  rl += t - p_lo;
  fast_two_sum(rh, rl, rh, rl);
  // -= n*cw3 (single double; matters because n*cw3 ~ 1e-24 for |n| ~ 1024).
  two_sum(rh, -(n_d * ln2_cw3), rh, t);
  rl += t;
  float64x2 rr;
  fast_two_sum(rh, rl, rr._limbs[0], rr._limbs[1]);

  // Route the residue rr through exp2_from_reduced. rr is in natural-log
  // form (|rr| ≤ ln2/2 ≈ 0.347); convert to base-2 by one DD multiply
  // y_frac = rr · log2(e), |y_frac| ≤ 0.5. The integer n from the Cody-
  // Waite reduction is already the correct 2^n scale factor.
  float64x2 log2e(log2_e_hi, log2_e_lo);
  float64x2 y_frac = rr * log2e;
  return exp2_from_reduced(n_d, y_frac._limbs[0], y_frac._limbs[1]);
}

float64x2 log2_kernel(float64x2 const &x) {
  double hi = x._limbs[0];
  if (hi > 15.0 / 16.0 && hi < 17.0 / 16.0) {
    // direct path: t = (x - 1) / (x + 1)
    float64x2 num = x - float64x2(1.0);
    float64x2 den = x + float64x2(1.0);
    float64x2 t = num / den;
    float64x2 t_sq = t * t;
    // 11 coefs → degree 10 → Estrin case in neval.
    float64x2 p = neval(t_sq, log2_wide_hi, log2_wide_lo, 10);
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
  float64x2 p = neval(x, expm1_taylor_hi, expm1_taylor_lo, 24);
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
      // arg has hi == 0 but lo > 0: the "real" value is the lo limb alone
      // (a positive subnormal). Passing the pair as-is to log_full would
      // trip its leading-limb zero check and return −∞. Promote lo to the
      // leading position so log_full sees a normalized operand.
      return log_full(float64x2(arg._limbs[1], 0.0));
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
  float64x2 p = neval(t2, log1p_taylor_hi, log1p_taylor_lo, 17);
  return t * p;
}

// ---- pow -------------------------------------------------------------------
// Composition x^y = 2^(y·log2(x)). A plain DD product y·log2(x) carries
// relative error ~ulp_dd, which for large |y·log2(x)| becomes an absolute
// error that exp2 amplifies (ulp_dd·|z|·ln2 → tens of ulp_dd in the
// output). To sidestep that, factor log2(x) = e_x + log2(m) where e_x is
// an integer (exact) and m ∈ [1,2); then y·e_x is exact in triple-double
// form (two_prod + y.lo·e_x), and y·log2(m) is DD with absolute error
// ≲ ulp_dd·|y|. The fractional part of the sum — which drives 2^r — is
// thus accurate to ~ulp_dd·|y| in absolute terms regardless of |e_x|.
float64x2 pow_full(float64x2 const &x, float64x2 const &y) {
  if (x._limbs[0] == 0.0 && y._limbs[0] == 0.0) return float64x2(1.0);
  if (!std::isfinite(x._limbs[0]) || !std::isfinite(y._limbs[0]) ||
      x._limbs[0] <= 0.0) {
    return exp_full(y * log_full(x));
  }

  // x = 2^e_x · m, m ∈ [1, 2). log2_kernel(m) is DD with |Lf| < 1.
  int e_x = std::ilogb(x._limbs[0]);
  float64x2 m;
  m._limbs[0] = std::ldexp(x._limbs[0], -e_x);
  m._limbs[1] = std::ldexp(x._limbs[1], -e_x);
  float64x2 Lf = log2_kernel(m);
  // Table-path third-limb residual: log2(center) − (v_hi + v_lo) is stored
  // in log2_values_extra[idx]; ignored on the direct (near-1) path where
  // no table entry is used and the residue stays within one DD computation.
  double Lf_extra = 0.0;
  double mh = m._limbs[0];
  if (!(mh > 15.0 / 16.0 && mh < 17.0 / 16.0)) {
    int idx = static_cast<int>((mh - 1.0) * 32.0);
    if (idx < 0) idx = 0;
    if (idx > 31) idx = 31;
    Lf_extra = log2_values_extra[idx];
  }
  double e_x_d = static_cast<double>(e_x);

  // y · e_x as triple-double (exact up to the y.lo·e_x cross product).
  double AH, AM;
  two_prod(y._limbs[0], e_x_d, AH, AM);
  double AL_h, AL_l;
  two_prod(y._limbs[1], e_x_d, AL_h, AL_l);

  // y · Lf as TD: exact products y.hi·Lf.hi and y.hi·Lf.lo plus the
  // y.lo·Lf.hi cross term. Capturing the residuals as separate limbs
  // avoids the ~2e-31 rounding in DD mul's inner `p01 += p10` step on the
  // worst fuzz input (|y·Lf| ~ 17).
  double P1, P1e; two_prod(y._limbs[0], Lf._limbs[0], P1, P1e);
  double P2, P2e; two_prod(y._limbs[0], Lf._limbs[1], P2, P2e);
  double P3 = y._limbs[1] * Lf._limbs[0];  // |y.lo·Lf.hi| ~ 6e-16
  // y.hi·Lf_extra ≲ 2e-33·|y|; y.lo·Lf.lo / y.lo·Lf_extra are ~1e-50, drop.
  double PX = y._limbs[0] * Lf_extra + y._limbs[1] * Lf._limbs[1];

  // Accumulate A + B into (h, l, c). Magnitude-descending: P1 (~|y·Lf|),
  // AM, AL_h, P3 (~1e-15), then sub-ulp residuals P1e/P2/P2e/PX.
  double h = AH, l = 0.0, c = 0.0, t, u;
  two_sum(h, P1,   h, t); two_sum(l, t, l, u); c += u;
  two_sum(h, AM,   h, t); two_sum(l, t, l, u); c += u;
  two_sum(h, AL_h, h, t); two_sum(l, t, l, u); c += u;
  two_sum(h, P3,   h, t); two_sum(l, t, l, u); c += u;
  two_sum(l, P1e,  l, u); c += u;
  two_sum(l, P2,   l, u); c += u;
  two_sum(l, P2e,  l, u); c += u;
  two_sum(l, PX,   l, u); c += u;
  c += AL_l;
  // Keep (h, l, c) as triple-double through the integer-split step. Folding
  // c into l with a plain `+=` would round it away whenever |l| > ulp⁻¹·|c|,
  // which is exactly the fuzz worst case.
  double e;
  two_sum(l, c, l, e);  // l + e = old_l + c exact
  fast_two_sum(h, l, h, l);  // normalize DD pair; e is third limb.

  // Split z ≈ h + l + e into n + r with |r| ≤ 0.5. Every two_sum below is
  // exact; only the final `rl += f1 + f2` rounds, and its dropped bits
  // live at ulp(rl) ~ 1e-33, well below the DD output resolution.
  double n_d = std::nearbyint(h);
  double rh, rl, s1, s2, f1, f2;
  two_sum(h, -n_d, rh, t);  // rh + t = h - n_d exact
  two_sum(t, l, s1, f1);    // s1 + f1 = t + l exact
  two_sum(s1, e, s2, f2);   // s2 + f2 = s1 + e exact
  rl = s2 + (f1 + f2);
  fast_two_sum(rh, rl, rh, rl);

  return exp2_from_reduced(n_d, rh, rl);
}
