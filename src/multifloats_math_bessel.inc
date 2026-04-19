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
    // Power series (j0q.c): J0(x) = 1 - x²/4 + x⁴·N(x²)/D(x²).
    float64x2 z = ax * ax;
    float64x2 r = z * z * neval(z, j0_J0_2N_hi, j0_J0_2N_lo, 6) /
                      deval(z, j0_J0_2D_hi, j0_J0_2D_lo, 6);
    return r - z * float64x2(0.25) + float64x2(1.0);
  }

  // Hankel asymptotic (x > 2):
  //   J0(x) = sqrt(2/(πx)) · (P(1/x²)·cos(x-π/4) - Q(1/x²)·sin(x-π/4))
  // P,Q picked from 8 sub-intervals via bessel_pq0 (j0q.c fits).
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
    // Power series (j1q.c): J1(x) = x/2 + x·x²·N(x²)/D(x²).
    float64x2 z = ax * ax;
    res = ax * float64x2(0.5) + ax * z * neval(z, j1_J1_2N_hi, j1_J1_2N_lo, 6) /
                                     deval(z, j1_J1_2D_hi, j1_J1_2D_lo, 6);
  } else {
    // Hankel asymptotic (x > 2), phase shifted by 3π/4:
    //   J1(x) = sqrt(2/(πx)) · (P(1/x²)·cos(x-3π/4) - Q(1/x²)·sin(x-3π/4))
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
    // Continued-fraction start: iterate until q1 exceeds the precision
    // target so J_{n+k}/J_n is small enough to truncate. For DD we need
    // ~2^106 ≈ 8e31 (vs ~2^56 ≈ 1e17 for plain double). The upper cap on k
    // prevents a runaway on pathological inputs (q1 always grows faster
    // than geometric with rate ≥ 2/x, so 1000 iterations is a vast margin).
    float64x2 w = float64x2(static_cast<double>(2 * n) / xx);
    float64x2 h = float64x2(2.0 / xx);
    float64x2 q0 = w;
    float64x2 z = w + h;
    float64x2 q1 = w * z - float64x2(1.0);
    int k = 1;
    while (q1._limbs[0] < 1e32 && k < 1000) {
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

