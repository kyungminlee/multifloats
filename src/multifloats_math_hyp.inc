// ---- sinh / cosh / tanh ----------------------------------------------------
float64x2 sinh(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::sinh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) < 0.5) {
    // 13-term Taylor: x * (1 + x²/6 + x⁴/120 + ... + x²²/23!).
    // 13 coefs → degree 12 → Estrin case in neval. Cutoff at |x|<0.5
    // bypasses the (e^x − e^-x) cancellation that costs ~3 ulp_dd.
    return x * neval(x * x, sinh_taylor_hi, sinh_taylor_lo, 12);
  }
  float64x2 e = exp(x);
  float64x2 ei = exp(-x);
  return (e - ei) * float64x2(0.5, 0.0);
}

float64x2 cosh(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::cosh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  float64x2 e = exp(x);
  float64x2 ei = exp(-x);
  return (e + ei) * float64x2(0.5, 0.0);
}

// Fused sinh/cosh. Shares the two exp evaluations (exp(x) and exp(-x)) that
// sinh and cosh would do independently, halving the exp cost
// for call sites that need both (e.g. complex sin/cos/sinh/cosh).
void sinhcosh(float64x2 const &x, float64x2 &s, float64x2 &c) {
  if (!std::isfinite(x._limbs[0])) {
    s._limbs[0] = std::sinh(x._limbs[0]); s._limbs[1] = 0.0;
    c._limbs[0] = std::cosh(x._limbs[0]); c._limbs[1] = 0.0;
    return;
  }
  if (std::abs(x._limbs[0]) < 0.5) {
    // Small-|x| Taylor: sinh via shared coefficients, cosh via the
    // non-cancellation-prone exp form (cosh's leading 1+1=2 doesn't cancel).
    float64x2 x2 = x * x;
    // 13 coefs → degree 12 → Estrin case in neval.
    s = x * neval(x2, sinh_taylor_hi, sinh_taylor_lo, 12);
    c = cosh(x);
    return;
  }
  float64x2 e  = exp(x);
  float64x2 ei = exp(-x);
  float64x2 half = float64x2(0.5, 0.0);
  s = (e - ei) * half;
  c = (e + ei) * half;
}

float64x2 tanh(float64x2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::tanh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  // For |x| ≥ 38, 2·exp(-2|x|) < 2⁻¹⁰⁹ ≈ 0.16 · ulp_dd(1), so tanh(x)
  // rounds to ±1 to full DD. (At |x| = 37, 2·exp(-74) ≈ 1.5 · ulp_dd(1),
  // so saturating there would lose ~1.5 DD ulp at the threshold.)
  if (x._limbs[0] >= 38.0) return float64x2(1.0);
  if (x._limbs[0] <= -38.0) return float64x2(-1.0);
  if (std::abs(x._limbs[0]) < 0.5) {
    return sinh(x) / cosh(x);
  }
  // tanh(|x|) = (1 - em2) / (1 + em2) where em2 = exp(-2|x|) is small
  // for |x| > 0.5 — avoids the 1 − small cancellation of 1 − 2/(e^(2|x|)+1).
  bool sign = x._limbs[0] < 0.0;
  float64x2 neg_two_ax;
  neg_two_ax._limbs[0] = sign ? (2.0 * x._limbs[0]) : (-2.0 * x._limbs[0]);
  neg_two_ax._limbs[1] = sign ? (2.0 * x._limbs[1]) : (-2.0 * x._limbs[1]);
  float64x2 em2 = exp(neg_two_ax);
  float64x2 res = (float64x2(1.0) - em2) / (float64x2(1.0) + em2);
  return sign ? -res : res;
}
