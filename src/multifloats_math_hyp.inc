// ---- sinh / cosh / tanh ----------------------------------------------------
float64x2 sinh_full(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::sinh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) < 0.1) {
    // 9-term Taylor: x * (1 + x²/6 + x⁴/120 + ...).
    // 9 coefs → degree 8 → Estrin case in neval.
    return x * neval(x * x, sinh_taylor_hi, sinh_taylor_lo, 8);
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
    // 9 coefs → degree 8 → Estrin case in neval.
    s = x * neval(x2, sinh_taylor_hi, sinh_taylor_lo, 8);
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
