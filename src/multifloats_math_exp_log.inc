// ---- exp2 / exp / log2 / log / log10 ---------------------------------------
float64x2 exp2_kernel(float64x2 const &x) {
  // n = nearest integer to x.hi, y = (x - n)/8
  double n_float = std::nearbyint(x._limbs[0]);
  float64x2 y = x + float64x2(-n_float, 0.0);
  y._limbs[0] = std::ldexp(y._limbs[0], -3);
  y._limbs[1] = std::ldexp(y._limbs[1], -3);
  // poly(y) then cube via three squarings (to undo the /8).
  // 14 coefs → degree 13 → Estrin case in neval.
  float64x2 p = neval(y, exp2_coefs_hi, exp2_coefs_lo, 13);
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
    // 9 coefs → degree 8 → Estrin case in neval.
    float64x2 p = neval(t_sq, log2_wide_hi, log2_wide_lo, 8);
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
  float64x2 p = horner(t2, log1p_taylor_hi, log1p_taylor_lo, 18);
  return t * p;
}

// ---- pow -------------------------------------------------------------------
float64x2 pow_full(float64x2 const &x, float64x2 const &y) {
  if (x._limbs[0] == 0.0 && y._limbs[0] == 0.0) return float64x2(1.0);
  return exp_full(y * log_full(x));
}
