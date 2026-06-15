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

// sin(x − φ) and cos(x − φ) for the Hankel phase φ ∈ {π/4, 3π/4}, as DD.
//
// Forming the angle (x − φ) as a DD rounds it to ~ulp(x); for large x that
// ~x·2⁻¹⁰⁴ absolute error is amplified into a large *relative* error near the
// zeros of Y0/Y1 (where P·sin+Q·cos cancels). Instead evaluate sin x / cos x
// to triple-double via sincos_td (proper π/2 range reduction, no angle
// rounding) and apply the exact angle-subtraction identities — exact because
// cos(π/4) = sin(π/4) = 1/√2:
//   φ = π/4 : sin = (sin x − cos x)/√2,   cos = (sin x + cos x)/√2
//   φ = 3π/4: sin = −(sin x + cos x)/√2,  cos = (sin x − cos x)/√2
static void bessel_hankel_phase(float64x2 const &x, bool three_quarter,
                                float64x2 &s, float64x2 &c) {
  float64x3 sx, cx;
  sincos_td(x, sx, cx);
  float64x3 sum = td_add_td(sx, cx);                 // sin x + cos x
  float64x3 dif = td_add_td(sx, td_negate(cx));      // sin x − cos x
  float64x2 r2(inv_sqrt2_hi, inv_sqrt2_lo);
  if (!three_quarter) {
    s = td_to_dd(td_mul_dd(dif, r2));
    c = td_to_dd(td_mul_dd(sum, r2));
  } else {
    s = td_to_dd(td_mul_dd(td_negate(sum), r2));
    c = td_to_dd(td_mul_dd(dif, r2));
  }
}

// (2k/x)·b in triple-double for the Y-recurrence. The coefficient 2k/x is a DD
// division (~1e-32); near the zeros of Y_n that error is amplified into the
// result, so carry it to ~TD: c = round_dd(2k/x) plus a correction
// ccorr = (2k − c·x)/x = (2k/x − c), with c·x evaluated to TD so the residual
// is exact. (c + ccorr)·b is then accumulated in TD.
static float64x3 bessel_recur_term(float64x3 const &b, int twok,
                                   float64x2 const &x) {
  double tk = static_cast<double>(twok);
  float64x2 c = float64x2(tk) / x;
  float64x3 cx = td_mul_dd(td_from_dd(c), x);                   // c·x (TD)
  float64x2 resid = td_to_dd(td_add_double(td_negate(cx), tk)); // 2k − c·x
  float64x2 ccorr = resid / x;                                  // 2k/x − c
  return td_add_td(td_mul_dd(b, c), td_mul_dd(b, ccorr));
}

float64x2 j0(float64x2 const &x) {
  float64x2 ax = x;
  if (ax.limbs[0] < 0.0) { ax.limbs[0] = -ax.limbs[0]; ax.limbs[1] = -ax.limbs[1]; }
  double xx = ax.limbs[0];
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
  float64x2 s, c;
  bessel_hankel_phase(ax, /*three_quarter=*/false, s, c);
  float64x2 xinv = float64x2(1.0) / ax;
  float64x2 z = xinv * xinv;
  float64x2 p, q;
  bessel_pq0(1.0 / xx, z, p, q);
  p = float64x2(1.0) + z * p;
  q = (z * q - float64x2(0.125)) * xinv;
  float64x2 tpi = float64x2(two_over_pi_hi, two_over_pi_lo);
  return multifloats::sqrt(tpi / ax) * (p * c - q * s);
}

float64x2 j1(float64x2 const &x) {
  // j1 is odd, so (+0, -eps) must be detected as negative — checking
  // only the hi limb would miss that case.
  bool neg = x.limbs[0] < 0.0 || (x.limbs[0] == 0.0 && x.limbs[1] < 0.0);
  float64x2 ax = neg ? -x : x;
  double xx = ax.limbs[0];
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
    float64x2 s, c;
    bessel_hankel_phase(ax, /*three_quarter=*/true, s, c);
    float64x2 xinv = float64x2(1.0) / ax;
    float64x2 z = xinv * xinv;
    float64x2 p, q;
    bessel_pq1(1.0 / xx, z, p, q);
    p = float64x2(1.0) + z * p;
    q = (z * q + float64x2(0.375)) * xinv;
    float64x2 tpi = float64x2(two_over_pi_hi, two_over_pi_lo);
    res = multifloats::sqrt(tpi / ax) * (p * c - q * s);
  }
  if (neg) { res.limbs[0] = -res.limbs[0]; res.limbs[1] = -res.limbs[1]; }
  return res;
}

float64x2 y0(float64x2 const &x) {
  double xx = x.limbs[0];
  if (xx <= 0.0 || !std::isfinite(xx)) {
    float64x2 r;
    r.limbs[0] = (xx == 0.0) ? -std::numeric_limits<double>::infinity()
                               : std::numeric_limits<double>::quiet_NaN();
    r.limbs[1] = 0.0;
    return r;
  }

  if (xx <= 1e-17) {
    float64x2 tpi = float64x2(two_over_pi_hi, two_over_pi_lo);
    return float64x2(bessel_u0_hi, bessel_u0_lo) + tpi * log(x);
  }

  if (xx <= 2.0) {
    float64x2 z = x * x;
    float64x2 p = neval(z, j0_Y0_2N_hi, j0_Y0_2N_lo, 7) /
             deval(z, j0_Y0_2D_hi, j0_Y0_2D_lo, 7);
    float64x2 tpi = float64x2(two_over_pi_hi, two_over_pi_lo);
    return tpi * log(x) * j0(x) + p;
  }

  float64x2 s, c;
  bessel_hankel_phase(x, /*three_quarter=*/false, s, c);
  float64x2 xinv = float64x2(1.0) / x;
  float64x2 z = xinv * xinv;
  float64x2 p, q;
  bessel_pq0(1.0 / xx, z, p, q);
  p = float64x2(1.0) + z * p;
  q = (z * q - float64x2(0.125)) * xinv;
  float64x2 tpi = float64x2(two_over_pi_hi, two_over_pi_lo);
  return multifloats::sqrt(tpi / x) * (p * s + q * c);
}

float64x2 y1(float64x2 const &x) {
  double xx = x.limbs[0];
  if (xx <= 0.0 || !std::isfinite(xx)) {
    float64x2 r;
    r.limbs[0] = (xx == 0.0) ? -std::numeric_limits<double>::infinity()
                               : std::numeric_limits<double>::quiet_NaN();
    r.limbs[1] = 0.0;
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
    return tpi * (log(x) * j1(x) - float64x2(1.0) / x) + p;
  }

  float64x2 s, c;
  bessel_hankel_phase(x, /*three_quarter=*/true, s, c);
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

float64x2 jn(int n, float64x2 const &x_in) {
  // Track sign flips from the identities
  //   J_{-n}(x) = (-1)^n J_n(x)     and     J_n(-x) = (-1)^n J_n(x).
  float64x2 x = x_in;
  int parity = 0;
  if (n < 0) { n = -n; parity ^= (n & 1); }
  if (n == 0) return j0(x);
  if (n == 1) {
    float64x2 r = j1(x);
    return parity ? -r : r;
  }
  if (x.limbs[0] < 0.0) { x = -x; parity ^= (n & 1); }

  double xx = x.limbs[0];
  if (xx == 0.0 || !std::isfinite(xx)) return float64x2(0.0);

  float64x2 result;
  if (static_cast<double>(n) <= xx) {
    // Forward recurrence.
    float64x2 a = j0(x);
    float64x2 b = j1(x);
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
    while (q1.limbs[0] < 1e32 && k < 1000) {
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
                   log(multifloats::abs(v * float64x2(static_cast<double>(n))));
    bool rescale = log_arg.limbs[0] >= 1.1356e4;
    for (int i = n - 1; i >= 1; --i) {
      float64x2 temp = b;
      b = float64x2(static_cast<double>(2 * i)) * b / x - a;
      a = temp;
      if (rescale && b.limbs[0] > 1e100) {
        a = a / b;
        t = t / b;
        b = float64x2(1.0);
      }
    }
    float64x2 z0 = j0(x);
    float64x2 z1 = j1(x);
    if (std::abs(z0.limbs[0]) >= std::abs(z1.limbs[0])) {
      result = t * z0 / b;
    } else {
      result = t * z1 / a;
    }
  }
  return parity ? -result : result;
}

float64x2 yn(int n, float64x2 const &x) {
  // Y_{-n}(x) = (-1)^n Y_n(x). No x-sign identity — Y is undefined on x < 0.
  int parity = 0;
  if (n < 0) { n = -n; parity ^= (n & 1); }
  double xx = x.limbs[0];
  if (xx <= 0.0) {
    float64x2 r;
    r.limbs[0] = (xx == 0.0) ? -std::numeric_limits<double>::infinity()
                              :  std::numeric_limits<double>::quiet_NaN();
    r.limbs[1] = 0.0;
    return r;
  }
  if (n == 0) {
    float64x2 r = y0(x);
    return parity ? -r : r;
  }
  if (n == 1) {
    float64x2 r = y1(x);
    return parity ? -r : r;
  }
  if (!std::isfinite(xx)) return float64x2(0.0);
  // Forward recurrence Y_{i+1} = (2i/x)·Y_i − Y_{i-1}, carried in triple-double
  // so the per-step rounding stays below the DD floor. The recurrence is stable
  // (|Y_n| grows with n), so the residue folds away cleanly at the end.
  float64x3 a = td_from_dd(y0(x));
  float64x3 b = td_from_dd(y1(x));
  for (int i = 1; i < n; ++i) {
    float64x3 nb = td_add_td(bessel_recur_term(b, 2 * i, x), td_negate(a));
    a = b;
    b = nb;
    if (!std::isfinite(b.limbs[0])) break;
  }
  float64x2 r = td_to_dd(b);
  return parity ? -r : r;
}

// Range variant of Jn: Jn forward-recurrence is unstable for n > x, so
// each output needs its own Miller's backward-recurrence seed — no
// single-sweep speedup. Loop over the scalar kernel; the .inc-inclusion
// model means this call can be inlined by the compiler without LTO.
// Caller guarantees out[] has room for n2−n1+1 elements.
void jn(int n1, int n2, float64x2 const &x, float64x2 *out) {
  int size = n2 - n1 + 1;
  if (size <= 0) return;
  for (int i = 0; i < size; ++i) out[i] = jn(n1 + i, x);
}

// Range variant: compute Y_{n1}, Y_{n1+1}, …, Y_{n2}(x) into `out` via a
// single forward recurrence sweep. Cheaper than n2−n1+1 independent calls.
// Caller guarantees out[] has room for n2−n1+1 elements.
void yn(int n1, int n2, float64x2 const &x, float64x2 *out) {
  int size = n2 - n1 + 1;
  if (size <= 0) return;
  // Triple-double forward recurrence (see scalar yn for the rationale), folding
  // each requested order back to DD as it is produced.
  float64x3 a = td_from_dd(y0(x));
  float64x3 b = td_from_dd(y1(x));
  if (n1 == 0) out[0] = td_to_dd(a);
  if (n1 <= 1 && n2 >= 1) out[1 - n1] = td_to_dd(b);
  for (int n = 2; n <= n2; ++n) {
    float64x3 nb = td_add_td(bessel_recur_term(b, 2 * (n - 1), x), td_negate(a));
    a = b;
    b = nb;
    if (n >= n1) out[n - n1] = td_to_dd(b);
    if (!std::isfinite(b.limbs[0])) {
      float64x2 bd = td_to_dd(b);
      for (int j = std::max(n - n1 + 1, 0); j < size; ++j) out[j] = bd;
      return;
    }
  }
}

