// Shared helpers for the C++ test / fuzz / bench suites.
//
// Bridges multifloats::float64x2 and __float128 (libquadmath) so tests
// can project DD results up to qp and back. Requires the compiler to
// provide __float128 (GCC) and link against libquadmath.

#pragma once

#include "multifloats.h"

#include <quadmath.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>

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

// Shared fuzz input generator — mirrors fuzz.f90's generate_pair.
//
// Single source of truth for the C++ fuzz distribution: fuzz.cc,
// test_boost_dd.cc, and ../integration/crosscheck.cc all draw from this
// struct, so an (iterations, seed) pair replays the same inputs across
// the three binaries and a failure in one can be bisected against the
// others' oracles. (fsrc's fuzz.fypp carries the Fortran translation.)
//
// Mix: 10% non-finite, 10% close numbers, 10% sum-near-zero cancellation,
// 10% near-huge, 10% near-tiny, 50% wide random 10^[-30,30]. This exposes
// cancellation, overflow, and subnormal regimes that uniform-magnitude
// inputs miss.
struct Rng {
  std::mt19937_64 engine;
  std::uniform_real_distribution<double> u01{0.0, 1.0};

  explicit Rng(uint64_t seed) : engine(seed) {}
  double u() { return u01(engine); }

  q_t pick_nonfinite(double r) {
    if (r < 0.25) return (q_t) (+1.0 / 0.0);
    if (r < 0.50) return (q_t) (-1.0 / 0.0);
    if (r < 0.70) return (q_t) (0.0 / 0.0);
    if (r < 0.85) return (q_t) (+0.0);
    return (q_t) (-0.0);
  }

  // Wide random: sign * uniform(0.5) * 10^k  with k ∈ [-30, 30].
  q_t wide(double r, double rexp) {
    int k = (int)(rexp * 60.0) - 30;
    q_t mag = powq((q_t)10.0q, (q_t)k);
    return (q_t)(r - 0.5) * mag;
  }

  // Narrow random: sign * uniform(0.5) * 10^k with k ∈ [-3, 3]. Used for
  // complex inputs — keeps re*re - im*im in c_mul / the cdivq division
  // well away from overflow and catastrophic-cancellation regimes so the
  // qp oracle stays a clean reference.
  q_t narrow(double r, double rexp) {
    int k = (int)(rexp * 6.0) - 3;
    q_t mag = powq((q_t)10.0q, (q_t)k);
    return (q_t)(r - 0.5) * mag;
  }

  void generate_pair(q_t &q1, q_t &q2) {
    double r1 = u(), r2 = u(), r3 = u(), r4 = u();
    int mode = (int)(r1 * 10.0);
    switch (mode) {
    case 0: {
      q1 = pick_nonfinite(r2);
      q2 = pick_nonfinite(r3);
      break;
    }
    case 1: {
      int k = (int)(r3 * 20.0) - 10;
      q1 = (q_t)(r2 - 0.5) * powq((q_t)10.0q, (q_t)k);
      q2 = q1 * ((q_t)1.0q + (q_t)((r4 - 0.5) * 2e-15));
      break;
    }
    case 2: {
      int k = (int)(r3 * 20.0) - 10;
      q1 = (q_t)(r2 - 0.5) * powq((q_t)10.0q, (q_t)k);
      q2 = -q1 * ((q_t)1.0q + (q_t)((r4 - 0.5) * 2e-15));
      break;
    }
    case 3: {
      // Near-huge: magnitudes in 10^[20, 30]; stays clear of DBL_MAX so
      // binary ops don't all overflow, while still exercising the
      // large-exponent regime.
      int k = (int)(r3 * 10.0) + 20;
      q1 = (q_t)(r2 - 0.5) * 2.0q * powq((q_t)10.0q, (q_t)k);
      q2 = (q_t)(r4 - 0.5) * 2.0q * powq((q_t)10.0q, (q_t)k);
      break;
    }
    case 4: {
      // Near-tiny: magnitudes in 10^[-30, -20]; kept above the
      // kSubnormalFloor=1e-290 gate so the stat recorder doesn't
      // silently discard every draw.
      int k = -((int)(r3 * 10.0) + 20);
      q1 = (q_t)(r2 - 0.5) * 2.0q * powq((q_t)10.0q, (q_t)k);
      q2 = (q_t)(r4 - 0.5) * 2.0q * powq((q_t)10.0q, (q_t)k);
      break;
    }
    default: {
      q1 = wide(r2, r3);
      q2 = wide(r4, r1);
      break;
    }
    }
  }
};

} // namespace multifloats_test
