// Property-based fuzz tests for multifloats.hh.
//
// Each random pair of MultiFloat<double, 2> inputs is exercised through
// every public arithmetic operation and the result is compared against a
// __float128 reference. Per-operation precision statistics are printed
// at the end. Built with g++ + libquadmath.

#include "multifloats.hh"

#include <quadmath.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>

namespace mf = multifloats;
using MF2 = mf::MultiFloat<double, 2>;

// =============================================================================
// __float128 reference helpers
// =============================================================================

using q_t = __float128;

static q_t to_q(MF2 const &x) {
  return (q_t)x._limbs[0] + (q_t)x._limbs[1];
}

static double q_rel_err(q_t got, q_t expected) {
  if (expected == 0) {
    return (got == 0) ? 0.0 : (double)fabsq(got);
  }
  q_t diff = got - expected;
  if (diff < 0) {
    diff = -diff;
  }
  q_t mag = expected < 0 ? -expected : expected;
  return (double)(diff / mag);
}

static char const *qstr(q_t v) {
  static char buf[64];
  quadmath_snprintf(buf, sizeof(buf), "%.30Qg", v);
  return buf;
}

// =============================================================================
// Stats
// =============================================================================

struct Stats {
  double max_rel = 0.0;
  double sum_rel = 0.0;
  long count = 0;
  void update(double rel) {
    if (!std::isfinite(rel)) {
      return;
    }
    if (rel > max_rel) {
      max_rel = rel;
    }
    sum_rel += rel;
    ++count;
  }
};

static long g_failures = 0;
static const long kFailureLimit = 50;

static void check_q(char const *op, MF2 const &got, q_t expected, double tol,
                    Stats &stats, MF2 const &a, MF2 const &b) {
  q_t got_q = to_q(got);
  double rel = q_rel_err(got_q, expected);
  stats.update(rel);
  if (!(rel <= tol)) {
    ++g_failures;
    std::fprintf(stderr, "FAIL [%s] rel_err=%g > tol=%g\n", op, rel, tol);
    std::fprintf(stderr, "  a        = (%a, %a)\n", a._limbs[0], a._limbs[1]);
    std::fprintf(stderr, "  b        = (%a, %a)\n", b._limbs[0], b._limbs[1]);
    std::fprintf(stderr, "  expected = %s\n", qstr(expected));
    std::fprintf(stderr, "  got      = %s\n", qstr(got_q));
    std::fprintf(stderr, "  got_limbs= (%a, %a)\n", got._limbs[0],
                 got._limbs[1]);
    if (g_failures >= kFailureLimit) {
      std::fprintf(stderr, "Too many failures, aborting.\n");
      std::exit(1);
    }
  }
}

// =============================================================================
// Tolerances
// =============================================================================

static constexpr double kAddTol = 1e-29;
static constexpr double kMulTol = 1e-29;
static constexpr double kDivTol = 1e-27;

// =============================================================================
// Random input generation
// =============================================================================

struct Rng {
  std::mt19937_64 engine;
  std::uniform_real_distribution<double> uniform{-1.0, 1.0};
  std::uniform_int_distribution<int> exponent{-200, 200};
  std::uniform_int_distribution<int> mode{0, 4};

  explicit Rng(uint64_t seed) : engine(seed) {}

  double random_double() {
    int e = exponent(engine);
    return std::ldexp(uniform(engine), e);
  }

  // Build a normalized double-double from an arbitrary __float128 value.
  MF2 random_mf() {
    int m = mode(engine);
    q_t v;
    if (m == 0) {
      v = (q_t)random_double();
    } else if (m == 1) {
      v = (q_t)random_double() + (q_t)random_double();
    } else if (m == 2) {
      double hi = random_double();
      double lo = std::ldexp(uniform(engine), exponent(engine) - 60);
      v = (q_t)hi + (q_t)lo;
    } else if (m == 3) {
      // Sum of three random doubles to exercise renormalization.
      v = (q_t)random_double() + (q_t)random_double() + (q_t)random_double();
    } else {
      // A scaled fraction with full precision.
      v = (q_t)random_double() / (q_t)random_double();
    }
    // Truncate to a normalized double-double exactly the way the kernel
    // would: high limb = round-to-nearest double of v, low limb = the
    // rounding error stored as a double.
    double hi = (double)v;
    if (!std::isfinite(hi)) {
      MF2 r;
      r._limbs[0] = hi;
      r._limbs[1] = hi;
      return r;
    }
    double lo = (double)(v - (q_t)hi);
    double s = hi + lo;
    double err = lo - (s - hi);
    MF2 r;
    r._limbs[0] = s;
    r._limbs[1] = err;
    return r;
  }
};

// =============================================================================
// Driver
// =============================================================================

static void print_stats(char const *name, Stats const &s) {
  if (s.count == 0) {
    std::printf("  %-12s n=%ld\n", name, s.count);
    return;
  }
  std::printf("  %-12s n=%-8ld max_rel=%.3e  mean_rel=%.3e\n", name, s.count,
              s.max_rel, s.sum_rel / s.count);
}

int main(int argc, char **argv) {
  long iterations = 100000;
  uint64_t seed = 0xC0FFEEULL;
  if (argc > 1) {
    iterations = std::atol(argv[1]);
  }
  if (argc > 2) {
    seed = std::strtoull(argv[2], nullptr, 0);
  }

  Rng rng(seed);
  std::printf("[multifloats_fuzz] iterations=%ld seed=0x%llx\n", iterations,
              (unsigned long long)seed);

  Stats neg, add, sub, mul, div, fmin_st, fmax_st, abs_st, cmp;

  for (long i = 0; i < iterations; ++i) {
    MF2 a = rng.random_mf();
    MF2 b = rng.random_mf();

    if (!mf::isfinite(a) || !mf::isfinite(b)) {
      continue;
    }

    q_t qa = to_q(a);
    q_t qb = to_q(b);

    check_q("neg", -a, -qa, 0.0, neg, a, b);

    q_t qabs = qa < 0 ? -qa : qa;
    check_q("abs", mf::abs(a), qabs, 0.0, abs_st, a, b);

    check_q("add", a + b, qa + qb, kAddTol, add, a, b);
    check_q("sub", a - b, qa - qb, kAddTol, sub, a, b);
    check_q("mul", a * b, qa * qb, kMulTol, mul, a, b);

    if (b._limbs[0] != 0.0) {
      check_q("div", a / b, qa / qb, kDivTol, div, a, b);
    }

    check_q("fmin", mf::fmin(a, b), qa < qb ? qa : qb, 0.0, fmin_st, a, b);
    check_q("fmax", mf::fmax(a, b), qa < qb ? qb : qa, 0.0, fmax_st, a, b);

    // Comparison consistency vs __float128.
    ++cmp.count;
    bool ok_eq = (a == b) == (qa == qb);
    bool ok_ne = (a != b) == (qa != qb);
    bool ok_lt = (a < b) == (qa < qb);
    bool ok_gt = (a > b) == (qa > qb);
    bool ok_le = (a <= b) == (qa <= qb);
    bool ok_ge = (a >= b) == (qa >= qb);
    if (!(ok_eq && ok_ne && ok_lt && ok_gt && ok_le && ok_ge)) {
      ++g_failures;
      std::fprintf(stderr, "FAIL [cmp] a=(%a,%a) b=(%a,%a)\n", a._limbs[0],
                   a._limbs[1], b._limbs[0], b._limbs[1]);
      if (g_failures >= kFailureLimit) {
        std::fprintf(stderr, "Too many failures, aborting.\n");
        std::exit(1);
      }
    }
  }

  std::printf("[multifloats_fuzz] failures=%ld\n", g_failures);
  std::printf("[multifloats_fuzz] precision report (relative error vs "
              "__float128):\n");
  print_stats("neg", neg);
  print_stats("abs", abs_st);
  print_stats("add", add);
  print_stats("sub", sub);
  print_stats("mul", mul);
  print_stats("div", div);
  print_stats("fmin", fmin_st);
  print_stats("fmax", fmax_st);
  std::printf("  cmp          n=%-8ld\n", cmp.count);

  return g_failures == 0 ? 0 : 1;
}
