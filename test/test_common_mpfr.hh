// Shared helpers for the MPFR-based precision fuzz suite.
//
// Parallels test_common.hh but upgrades the reference from 113-bit
// __float128 to arbitrary-precision mpreal. Default precision is set
// to 200 bits on first use — comfortably above both DD (~106 bits) and
// float128 (113 bits) so the MPFR value can be treated as exact for
// both.

#pragma once

#include "multifloats.h"
#include "test_common.hh"

#include <mpreal.h>

#include <cmath>
#include <cstdio>

namespace multifloats_test {

// Reference precision. 200 bits is ~4× DD, ~1.8× float128 — plenty of
// headroom for both comparisons.
inline constexpr mpfr_prec_t kMpfrPrec = 200;

inline void init_mpfr_default_prec() {
  static bool done = false;
  if (!done) {
    mpfr::mpreal::set_default_prec(kMpfrPrec);
    done = true;
  }
}

using mp_t = mpfr::mpreal;

// DD → mpreal: both limbs exact (mpreal carries more precision than the
// sum of two doubles). Non-finite inputs inherit hi's class.
inline mp_t to_mp(multifloats::float64x2 const &x) {
  init_mpfr_default_prec();
  if (!std::isfinite(x.limbs[0])) return mp_t(x.limbs[0]);
  return mp_t(x.limbs[0]) + mp_t(x.limbs[1]);
}

inline mp_t to_mp(q_t v) {
  init_mpfr_default_prec();
  // 3-double split captures all 113 qp bits exactly: each limb fits in 53
  // bits, 3*53=159 ≥ 113, and the mpreal has 200-bit precision so the sums
  // are lossless. Avoids mpfr_set_float128, which requires MPFR built with
  // --enable-float128 (not on by default in most package managers).
  double hi = (double)v;
  if (!std::isfinite(hi)) return mp_t(hi);
  double lo  = (double)(v - (q_t)hi);
  double lo2 = (double)(v - (q_t)hi - (q_t)lo);
  mp_t r(hi);
  r += mp_t(lo);
  r += mp_t(lo2);
  return r;
}

inline mp_t to_mp(double x) {
  init_mpfr_default_prec();
  return mp_t(x);
}

// Relative error in double, |got - expected| / max(|expected|, tiny).
// When expected is exact zero we return |got| so a nonzero result is
// still surfaced (matches test_common.hh::q_rel_err).
inline double mp_rel_err(mp_t const &got, mp_t const &expected) {
  if (mpfr::isnan(got) || mpfr::isnan(expected)) return 0.0;
  mp_t diff = mpfr::abs(got - expected);
  mp_t mag  = mpfr::abs(expected);
  if (mag == 0) return diff.toDouble();
  return (diff / mag).toDouble();
}

inline bool mp_isnan(mp_t const &v)    { return mpfr::isnan(v); }
inline bool mp_isfinite(mp_t const &v) { return mpfr::isfinite(v); }

} // namespace multifloats_test
