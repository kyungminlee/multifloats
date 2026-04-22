// DD polynomial evaluators: moved here from multifloats.hh so the public
// header doesn't need to carry ~250 lines of switch-based Estrin bodies.
// These are consumed only by the kernel .inc files compiled into
// multifloats_math.cc (same TU), so the compiler still sees the bodies
// when picking the compile-time-constant-`n` switch arm. `inline` is kept
// as a hint to the inliner; it has no ODR effect in a single TU.

// ---- Horner polynomial evaluation in DD ------------------------------------
inline float64x2 horner(float64x2 const &y, double const *hi, double const *lo,
                       int n) {
  float64x2 p(hi[n - 1], lo[n - 1]);
  for (int i = n - 2; i >= 0; --i) {
    p = p * y + float64x2(hi[i], lo[i]);
  }
  return p;
}

// ---- Estrin polynomial evaluation (neval/deval from libquadmath) ------------
// neval: evaluates P[n]*x^n + P[n-1]*x^(n-1) + ... + P[0]
// Uses Estrin's scheme for ILP when degree is known, Horner fallback otherwise.
//
// NOTE on convention: neval's `n` is the polynomial *degree*, so the
// coefficient array has n+1 entries (c[0..n]). This differs from
// `horner` above, which takes `n` as the coefficient COUNT. A common
// bug when migrating Horner → Estrin is forgetting to decrement n by 1.
//
// Tier 3 measurements (cpp_bench, relative to Horner baseline):
//   exp/exp2 (case 13), sin/cos (case 12), sinh/cosh (case 8),
//   asinh/atanh (case 14): 1.78×–1.88× speedup with no precision loss
//   outside 1 DD ULP. Extending with an x^16 case would benefit
//   expm1_taylor (degree 24) and log1p_taylor (degree 17); that
//   extension was skipped as the additional gain (~1.5× on small-|x|
//   branches only) did not justify the code size.
inline float64x2 neval(float64x2 const &x, double const *hi, double const *lo,
                      int n) {
  auto c = [&](int i) { return float64x2(hi[i], lo[i]); };
  float64x2 x2 = x * x;
  float64x2 x4 = x2 * x2;
  float64x2 x8 = x4 * x4;
  switch (n) {
  case 7: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4;
  }
  case 8: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4 + c(8) * x8;
  }
  case 9: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4 + p89 * x8;
  }
  case 10: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p810 = p89 + c(10) * x2;
    return p03 + p47 * x4 + p810 * x8;
  }
  case 11: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    return p03 + p47 * x4 + p811 * x8;
  }
  case 12: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    return p03 + p47 * x4 + (p811 + c(12) * x4) * x8;
  }
  case 13: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p1213 = c(12) + c(13) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    return p03 + p47 * x4 + (p811 + p1213 * x4) * x8;
  }
  case 14: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p1213 = c(12) + c(13) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    float64x2 p1214 = p1213 + c(14) * x2;
    return p03 + p47 * x4 + (p811 + p1214 * x4) * x8;
  }
  case 15: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p1213 = c(12) + c(13) * x;
    float64x2 p1415 = c(14) + c(15) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    float64x2 p1215 = p1213 + p1415 * x2;
    return p03 + p47 * x4 + (p811 + p1215 * x4) * x8;
  }
  case 17: {
    float64x2 x16 = x8 * x8;
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p1213 = c(12) + c(13) * x;
    float64x2 p1415 = c(14) + c(15) * x;
    float64x2 p1617 = c(16) + c(17) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    float64x2 p1215 = p1213 + p1415 * x2;
    float64x2 p0_15 = p03 + p47 * x4 + (p811 + p1215 * x4) * x8;
    return p0_15 + p1617 * x16;
  }
  case 24: {
    float64x2 x16 = x8 * x8;
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p1213 = c(12) + c(13) * x;
    float64x2 p1415 = c(14) + c(15) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    float64x2 p1215 = p1213 + p1415 * x2;
    float64x2 p0_15 = p03 + p47 * x4 + (p811 + p1215 * x4) * x8;
    float64x2 p1617 = c(16) + c(17) * x;
    float64x2 p1819 = c(18) + c(19) * x;
    float64x2 p2021 = c(20) + c(21) * x;
    float64x2 p2223 = c(22) + c(23) * x;
    float64x2 p1619 = p1617 + p1819 * x2;
    float64x2 p2023 = p2021 + p2223 * x2;
    float64x2 p16_24 = p1619 + p2023 * x4 + c(24) * x8;
    return p0_15 + p16_24 * x16;
  }
  default: {
    float64x2 y = c(n);
    for (int i = n - 1; i >= 0; --i)
      y = y * x + c(i);
    return y;
  }
  }
}

// deval: evaluates x^(n+1) + P[n]*x^n + ... + P[0] (monic leading term)
inline float64x2 deval(float64x2 const &x, double const *hi, double const *lo,
                      int n) {
  auto c = [&](int i) { return float64x2(hi[i], lo[i]); };
  float64x2 x2 = x * x;
  float64x2 x4 = x2 * x2;
  float64x2 x8 = x4 * x4;
  switch (n) {
  case 7: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4 + x8;
  }
  case 8: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4 + p89 * x8;
  }
  case 9: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p810 = p89 + x2;
    return p03 + p47 * x4 + p810 * x8;
  }
  case 10: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    return p03 + p47 * x4 + p811 * x8;
  }
  case 11: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    return p03 + p47 * x4 + (p811 + x4) * x8;
  }
  default: {
    float64x2 y = x + c(n);
    for (int i = n - 1; i >= 0; --i)
      y = y * x + c(i);
    return y;
  }
  }
}
