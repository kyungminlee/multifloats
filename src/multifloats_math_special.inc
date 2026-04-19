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

