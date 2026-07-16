// Shared helpers for the C++ test / fuzz / bench suites.
//
// Bridges multifloats::float64x2 and __float128 (libquadmath) so tests
// can project DD results up to qp and back. Requires the compiler to
// provide __float128 (GCC) and link against libquadmath.

#pragma once

#include "multifloats.h"

#include <quadmath.h>

#include <cmath>
#include <cstdio>

namespace multifloats_test {

using q_t = __float128;

// DD → qp: exact when both limbs are finite (the DD library keeps
// |lo| <= 0.5 ULP of hi). For non-finite limbs qp inherits hi's class.
inline q_t to_q(multifloats::float64x2 const &x) {
  return (q_t)x.limbs[0] + (q_t)x.limbs[1];
}

// qp → DD: normalized via fast_two_sum (|hi| >= |lo| holds by
// construction on the finite path). For non-finite inputs both limbs
// carry the same non-finite class — this matches what the DD kernels
// themselves emit (see test_division_nonfinite).
inline multifloats::float64x2 from_q(q_t v) {
  double hi = (double)v;
  multifloats::float64x2 r;
  if (!std::isfinite(hi)) {
    r.limbs[0] = hi;
    r.limbs[1] = hi;
    return r;
  }
  double lo = (double)(v - (q_t)hi);
  double s = hi + lo;
  double err = lo - (s - hi);
  r.limbs[0] = s;
  r.limbs[1] = err;
  return r;
}

inline bool q_isnan(q_t x) { return isnanq(x); }
inline bool q_isfinite(q_t x) { return finiteq(x); }

// Relative error, with |got| used as the error when expected is zero
// (so an exact-zero expected value is not penalized when got is also
// zero, but any nonzero got is flagged at its magnitude).
inline double q_rel_err(q_t got, q_t expected) {
  q_t diff = got - expected;
  if (diff < 0) diff = -diff;
  q_t mag = expected < 0 ? -expected : expected;
  if (mag == 0) return (double)diff;
  return (double)(diff / mag);
}

// Pretty-print a __float128 to a static thread-unsafe buffer. Only
// intended for failure-diagnostic fprintf calls.
inline char const *qstr(q_t v) {
  static char buf[64];
  quadmath_snprintf(buf, sizeof(buf), "%.30Qg", v);
  return buf;
}

} // namespace multifloats_test
