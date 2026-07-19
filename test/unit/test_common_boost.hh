// Shared helpers for the boost::multiprecision::cpp_double_double
// comparison suites (test_boost_dd.cc, bench_boost_dd.cc, probe_bjn.cc).
//
// Parallels test_common.hh but bridges boost's DD backend instead of
// multifloats::float64x2, so the boost targets can project results to
// __float128 and back on the same terms as the multifloats kernels.

#pragma once

#include "multifloats.h"
#include "test_common.hh"

#include <boost/multiprecision/cpp_double_fp.hpp>

namespace multifloats_test {

// boost DD → qp: exact when both limbs are finite (same reasoning as
// to_q in test_common.hh).
inline q_t bdd_to_q(boost::multiprecision::cpp_double_double const &x) {
  auto const &r = x.backend().crep();
  return (q_t)r.first + (q_t)r.second;
}

// qp → boost DD: renormalize through the existing fast_two_sum path in
// test_common.hh so the boost limbs satisfy the same |lo| <= 0.5 ULP of
// |hi| invariant the multifloats DD representation maintains. Equivalent
// to constructing `cpp_double_double(hi) + cpp_double_double(lo)` (boost
// renormalizes on the constructor sum) but ~free.
inline boost::multiprecision::cpp_double_double bdd_from_q(q_t v) {
  multifloats::float64x2 mf = from_q(v);
  boost::multiprecision::cpp_double_double out;
  auto &r = out.backend().rep();
  r.first  = mf.limbs[0];
  r.second = mf.limbs[1];
  return out;
}

// double → boost DD: exact widening, lo limb zero.
inline boost::multiprecision::cpp_double_double bdd_from_d(double v) {
  boost::multiprecision::cpp_double_double out;
  auto &r = out.backend().rep();
  r.first  = v;
  r.second = 0.0;
  return out;
}

} // namespace multifloats_test
