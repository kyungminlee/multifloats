// Triple-double scratch primitives — internal header. NOT part of the
// public ABI surface in multifloats.hh. Pulled in by multifloats_math.cc
// and by the direct-primitive tests in test/test.cc.
//
// Narrow-scope toolkit used by kernels whose DD output would otherwise
// suffer cancellation below the DD floor (e.g. cexpm1 Re near the
// cancellation surface cos(b)·e^a = 1). Not a general three-limb
// multifloat — just enough primitives to carry a residue through a
// single TD × TD → TD ⊖ 1 pipeline and fold back to DD at the output.
//
// See doc/developer/TRIPLE_DOUBLE.md.
#pragma once

#include "multifloats.hh"

namespace multifloats {
namespace detail {

struct float64x3 {
  double _limbs[3] = {0.0, 0.0, 0.0};
  constexpr float64x3() = default;
  constexpr float64x3(double h, double m, double l) : _limbs{h, m, l} {}
};

// a + b + c = s + t + u, exact. Two_sum variant; input ordering is not
// required. Outputs are not renormalized — callers that need the
// normalized form should follow with renorm3.
inline void three_sum(double a, double b, double c,
                      double &s, double &t, double &u) {
  double u1, v1, w;
  two_sum(a, b, u1, v1);   // a + b = u1 + v1 exact
  two_sum(u1, c, s, w);    // (a+b) + c = s + w exact
  two_sum(v1, w, t, u);    // v1 + w = t + u exact
}

// Single-step TSUM accumulator: (T0, T1, T2) += x with two sequential
// two_sums. Residuals cascade down; the last residual is absorbed into
// T2 with a plain `+=` (its magnitude sits below ulp(T1), so the rounding
// loses at most one bit of the TD's least-significant limb — well below
// the DD output resolution after td_to_dd).
inline void tsum3(double &T0, double &T1, double &T2, double x) {
  double e;
  two_sum(T0, x, T0, e);
  two_sum(T1, e, T1, e);
  T2 += e;
}

// Final renormalize of (T0, T1, T2) produced by a TSUM chain: two
// two_sum passes canonicalize the triple. Preserves the third limb.
inline void tsum3_finalize(double &T0, double &T1, double &T2) {
  double e;
  two_sum(T1, T2, T1, e); T2 = e;
  two_sum(T0, T1, T0, e); T1 = e;
}

// Renormalize three doubles to a TD with |m| ≤ ulp(h)/2, |l| ≤ ulp(m)/2.
float64x3 renorm3(double h, double m, double l);

// TD + exact double, renormalized.
float64x3 td_add_double(float64x3 const &a, double d);
float64x3 td_sub_double(float64x3 const &a, double d);

// TD + DD → TD. DD expands to two scalar doubles; TSUM accumulator keeps
// the exact sum within the 3-limb output (modulo the 3rd-limb absorb).
float64x3 td_add_dd(float64x3 const &a, float64x2 const &b);

// Exact per-limb negation. The TD sum is linear, so negating all three
// limbs negates the TD value exactly (no rounding).
float64x3 td_negate(float64x3 const &a);

// TD + TD → TD. Expands the two triples into six doubles and feeds them
// through the TSUM accumulator in a single descending pass.
float64x3 td_add_td(float64x3 const &a, float64x3 const &b);

// TD × DD → TD. Six two_prods + a scalar product (the DD has no third
// limb, so the three 9-way products involving b's missing lo are dropped;
// net cost is ~40% fewer ops than td_mul_td).
float64x3 td_mul_dd(float64x3 const &a, float64x2 const &b);

// TD × TD → TD. Exact 9-way scalar expansion via two_prod (FMA), fed into
// a TSUM accumulator in magnitude-descending order so the high bits land
// in T0 and residuals in T1/T2. The two lowest-order terms a.lo·b.mid,
// a.mid·b.lo, a.lo·b.lo sit ~2^-159 below T0 — below the TD output
// resolution; kept for symmetry (TSUM absorbs them without cost).
float64x3 td_mul_td(float64x3 const &a, float64x3 const &b);

// TD → DD: fold third limb into second via two_sum, then canonicalize the
// leading pair. Dropped residue ≤ ulp(l) ≈ 2^-159 — invisible in DD output.
float64x2 td_to_dd(float64x3 const &a);

// DD → TD: trivial zero-extension (already normalized since DD pair is).
float64x3 td_from_dd(float64x2 const &a);

} // namespace detail
} // namespace multifloats
