namespace {
// Non-finite DD inputs (NaN / ±Inf) are routed through the scalar std:: lib
// call, with a zero lo-limb. The six inverse-trig `*_full` functions below
// share this preamble verbatim except for the std:: callee.
template <typename ScalarFn>
inline float64x2 _nonfinite_scalar_fallback(float64x2 const &x, ScalarFn fn) {
  float64x2 r;
  r._limbs[0] = fn(x._limbs[0]);
  r._limbs[1] = 0.0;
  return r;
}
} // namespace

// ---- atan (table lookup + rational polynomial, from libquadmath atanq.c) ---
// arctan(t) = t + t^3 * P(t^2) / Q(t^2),  |t| <= 0.09375
// Argument reduced via 87-entry table: arctan(x) = atantbl[k] + arctan(t).
// Table branch (|x| < 32/3): k = floor(8|x|+0.25), t = (x - k/8)/(1 + x*k/8);
//   worst-case |t| = 1/(16*(1 + (k/8)^2)) stays well below 0.09375.
// Large branch (|x| >= 32/3): k = 86 (pi/2), t = -1/|x|; cutover placed at
//   32/3 so worst-case |t| = 3/32 = 0.09375 hits the fitted edge exactly.
float64x2 atan_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    return _nonfinite_scalar_fallback(x, [](double d) { return std::atan(d); });
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
    return _nonfinite_scalar_fallback(x, [](double d) { return std::asin(d); });
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
    return _nonfinite_scalar_fallback(x, [](double d) { return std::acos(d); });
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
    // When a DD has hi == 0 the effective sign (including ±0 and signed
    // subnormals) lives entirely in lo. Pass lo to libm's atan2 so the
    // IEEE signed-zero quadrant rules still produce the right answer.
    float64x2 r;
    r._limbs[0] = std::atan2(y._limbs[1], x._limbs[1]);
    r._limbs[1] = 0.0;
    return r;
  }
  float64x2 pi_dd = float64x2(pi_dd_hi, pi_dd_lo);
  float64x2 pi_half_dd = float64x2(pi_half_cw1, pi_half_cw2);
  float64x2 res;
  if (std::abs(x._limbs[0]) >= std::abs(y._limbs[0])) {
    res = atan_full(y / x);
    if (x._limbs[0] < 0.0) {
      // Use signbit so y = -0 picks the −π branch (IEEE: -0 >= 0 is true,
      // so `>= 0.0` mishandles the negative-real-axis branch cut).
      if (!std::signbit(y._limbs[0])) res = res + pi_dd;
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
    return _nonfinite_scalar_fallback(x, [](double d) { return std::asinh(d); });
  }
  if (std::abs(x._limbs[0]) < 0.01) {
    // 15 coefs → degree 14 → Estrin case in neval.
    return x * neval(x * x, asinh_taylor_hi, asinh_taylor_lo, 14);
  }
  bool sign = x._limbs[0] < 0.0;
  float64x2 ax = sign ? -x : x;
  if (ax._limbs[0] > 1e150) {
    // log(2|x|) asymptotic to avoid x² overflow
    float64x2 r = log_full(ax) + float64x2(ln_2_hi, ln_2_lo);
    return sign ? -r : r;
  }
  float64x2 axsq = ax * ax;
  float64x2 root = sqrt(axsq + float64x2(1.0));
  // For ax < 1, ax + root ≈ 1 → log cancels; route via log1p using
  // (root - 1) = ax² / (root + 1). For ax ≥ 1, plain log is accurate.
  float64x2 res;
  if (ax._limbs[0] < 1.0) {
    float64x2 rm1 = axsq / (root + float64x2(1.0));
    res = log1p_full(ax + rm1);
  } else {
    res = log_full(ax + root);
  }
  return sign ? -res : res;
}

float64x2 acosh_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    return _nonfinite_scalar_fallback(x, [](double d) { return std::acosh(d); });
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
    return _nonfinite_scalar_fallback(x, [](double d) { return std::atanh(d); });
  }
  if (std::abs(x._limbs[0]) > 1.0) {
    float64x2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) < 0.01) {
    // 15 coefs → degree 14 → Estrin case in neval.
    return x * neval(x * x, atanh_taylor_hi, atanh_taylor_lo, 14);
  }
  // atanh(x) = 0.5·log((1+x)/(1-x)) = 0.5·log1p(2x/(1-x)).
  // For 0.01 ≤ |x| ≪ 1 the ratio ≈ 1 and the direct log loses bits;
  // log1p of 2x/(1-x) (which stays small there) preserves them. For
  // |x| near 1 the argument to log1p becomes large and log1p routes
  // to log, so the same formula covers the whole remaining range.
  float64x2 den = float64x2(1.0) - x;
  float64x2 arg = (x + x) / den;
  return float64x2(0.5, 0.0) * log1p_full(arg);
}

