// Triple-double scratch primitives — out-of-line definitions matching the
// declarations in multifloats.hh (namespace multifloats::detail). Bodies
// live here so the header carries only the declaration; callers pay one
// function-call hop per TD op, which is negligible against the polynomial
// costs of the transcendental kernels that host them.
//
// Leaf helpers (two_sum, two_prod, tsum3, tsum3_finalize, three_sum) stay
// inline in the header — they're hot and tiny.

float64x3 renorm3(double h, double m, double l) {
  two_sum(m, l, m, l);
  two_sum(h, m, h, m);
  two_sum(m, l, m, l);
  return {h, m, l};
}

float64x3 td_add_double(float64x3 const &a, double d) {
  double T0 = a._limbs[0], T1 = a._limbs[1], T2 = a._limbs[2];
  tsum3(T0, T1, T2, d);
  tsum3_finalize(T0, T1, T2);
  return {T0, T1, T2};
}

float64x3 td_sub_double(float64x3 const &a, double d) {
  return td_add_double(a, -d);
}

float64x3 td_add_dd(float64x3 const &a, float64x2 const &b) {
  double T0 = 0.0, T1 = 0.0, T2 = 0.0;
  tsum3(T0, T1, T2, a._limbs[0]);
  tsum3(T0, T1, T2, b._limbs[0]);
  tsum3(T0, T1, T2, a._limbs[1]);
  tsum3(T0, T1, T2, b._limbs[1]);
  tsum3(T0, T1, T2, a._limbs[2]);
  tsum3_finalize(T0, T1, T2);
  return {T0, T1, T2};
}

float64x3 td_negate(float64x3 const &a) {
  return {-a._limbs[0], -a._limbs[1], -a._limbs[2]};
}

float64x3 td_add_td(float64x3 const &a, float64x3 const &b) {
  double T0 = 0.0, T1 = 0.0, T2 = 0.0;
  tsum3(T0, T1, T2, a._limbs[0]);
  tsum3(T0, T1, T2, b._limbs[0]);
  tsum3(T0, T1, T2, a._limbs[1]);
  tsum3(T0, T1, T2, b._limbs[1]);
  tsum3(T0, T1, T2, a._limbs[2]);
  tsum3(T0, T1, T2, b._limbs[2]);
  tsum3_finalize(T0, T1, T2);
  return {T0, T1, T2};
}

float64x3 td_mul_dd(float64x3 const &a, float64x2 const &b) {
  double a0 = a._limbs[0], a1 = a._limbs[1], a2 = a._limbs[2];
  double b0 = b._limbs[0], b1 = b._limbs[1];
  double p00h, p00l; two_prod(a0, b0, p00h, p00l);
  double p01h, p01l; two_prod(a0, b1, p01h, p01l);
  double p10h, p10l; two_prod(a1, b0, p10h, p10l);
  double p11h, p11l; two_prod(a1, b1, p11h, p11l);
  double p20h, p20l; two_prod(a2, b0, p20h, p20l);
  double p21 = a2 * b1;  // ~2^-159, below TD output resolution
  double T0 = 0.0, T1 = 0.0, T2 = 0.0;
  tsum3(T0, T1, T2, p00h);
  tsum3(T0, T1, T2, p01h); tsum3(T0, T1, T2, p10h);
  tsum3(T0, T1, T2, p00l);
  tsum3(T0, T1, T2, p11h); tsum3(T0, T1, T2, p20h);
  tsum3(T0, T1, T2, p01l); tsum3(T0, T1, T2, p10l);
  tsum3(T0, T1, T2, p21);
  tsum3(T0, T1, T2, p11l); tsum3(T0, T1, T2, p20l);
  tsum3_finalize(T0, T1, T2);
  return {T0, T1, T2};
}

float64x3 td_mul_td(float64x3 const &a, float64x3 const &b) {
  double a0 = a._limbs[0], a1 = a._limbs[1], a2 = a._limbs[2];
  double b0 = b._limbs[0], b1 = b._limbs[1], b2 = b._limbs[2];
  double p00h, p00l; two_prod(a0, b0, p00h, p00l);  // O(1)
  double p01h, p01l; two_prod(a0, b1, p01h, p01l);  // O(2^-53)
  double p10h, p10l; two_prod(a1, b0, p10h, p10l);  // O(2^-53)
  double p02h, p02l; two_prod(a0, b2, p02h, p02l);  // O(2^-106)
  double p11h, p11l; two_prod(a1, b1, p11h, p11l);  // O(2^-106)
  double p20h, p20l; two_prod(a2, b0, p20h, p20l);  // O(2^-106)
  double p12 = a1 * b2;                             // O(2^-159)
  double p21 = a2 * b1;                             // O(2^-159)
  // a2*b2 ~ 2^-212 below TD output; drop.
  double T0 = 0.0, T1 = 0.0, T2 = 0.0;
  tsum3(T0, T1, T2, p00h);
  tsum3(T0, T1, T2, p01h); tsum3(T0, T1, T2, p10h);
  tsum3(T0, T1, T2, p00l);
  tsum3(T0, T1, T2, p02h); tsum3(T0, T1, T2, p11h); tsum3(T0, T1, T2, p20h);
  tsum3(T0, T1, T2, p01l); tsum3(T0, T1, T2, p10l);
  tsum3(T0, T1, T2, p12);  tsum3(T0, T1, T2, p21);
  tsum3(T0, T1, T2, p02l); tsum3(T0, T1, T2, p11l); tsum3(T0, T1, T2, p20l);
  tsum3_finalize(T0, T1, T2);
  return {T0, T1, T2};
}

float64x2 td_to_dd(float64x3 const &a) {
  double h = a._limbs[0], m = a._limbs[1], l = a._limbs[2];
  double e;
  two_sum(m, l, m, e);      // |e| ≤ ulp(m)/2
  two_sum(h, m, h, m);      // canonical leading pair
  m += e;                   // absorb third-limb residue
  float64x2 r;
  r._limbs[0] = h;
  r._limbs[1] = m;
  return r;
}

float64x3 td_from_dd(float64x2 const &a) {
  return {a._limbs[0], a._limbs[1], 0.0};
}
