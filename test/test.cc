// Unit tests for multifloats.hh.
//
// Each public operation is tested over a curated list of inputs and the
// result is compared against a __float128 reference. A per-operation
// precision report is printed at the end.
//
// Built with g++ + libquadmath (Apple Clang lacks __float128).

#include "multifloats.hh"

#include <quadmath.h>

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

namespace mf = multifloats;
using MF2 = mf::MultiFloat<double, 2>;

// =============================================================================
// __float128 reference helpers
// =============================================================================

using q_t = __float128;

static q_t to_q(MF2 const &x) {
  return (q_t)x._limbs[0] + (q_t)x._limbs[1];
}

static MF2 from_q(q_t v) {
  double hi = (double)v;
  if (!std::isfinite(hi)) {
    MF2 r;
    r._limbs[0] = hi;
    r._limbs[1] = hi;
    return r;
  }
  double lo = (double)(v - (q_t)hi);
  // fast_two_sum normalization (|hi| >= |lo| holds by construction).
  double s = hi + lo;
  double err = lo - (s - hi);
  MF2 r;
  r._limbs[0] = s;
  r._limbs[1] = err;
  return r;
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
// Stats / reporting
// =============================================================================

struct Stats {
  double max_rel = 0.0;
  double sum_rel = 0.0;
  int count = 0;
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

static int g_failures = 0;
static int g_checks = 0;

#define REQUIRE(cond)                                                          \
  do {                                                                         \
    ++g_checks;                                                                \
    if (!(cond)) {                                                             \
      ++g_failures;                                                            \
      std::fprintf(stderr, "FAIL %s:%d  REQUIRE(%s)\n", __FILE__, __LINE__,    \
                   #cond);                                                     \
    }                                                                          \
  } while (0)

// Compare a multifloat result against a __float128 reference, accumulate
// statistics, and fail loudly if the relative error exceeds `tol`.
static void check_q(char const *op, MF2 const &got, q_t expected, double tol,
                    Stats &stats, MF2 const *a = nullptr,
                    MF2 const *b = nullptr) {
  ++g_checks;
  q_t got_q = to_q(got);
  double rel = q_rel_err(got_q, expected);
  stats.update(rel);
  if (!(rel <= tol)) {
    ++g_failures;
    std::fprintf(stderr, "FAIL [%s] rel_err=%g > tol=%g\n", op, rel, tol);
    if (a) {
      std::fprintf(stderr, "  a        = (%a, %a)\n", a->_limbs[0],
                   a->_limbs[1]);
    }
    if (b) {
      std::fprintf(stderr, "  b        = (%a, %a)\n", b->_limbs[0],
                   b->_limbs[1]);
    }
    std::fprintf(stderr, "  expected = %s\n", qstr(expected));
    std::fprintf(stderr, "  got      = %s\n", qstr(got_q));
    std::fprintf(stderr, "  got_limbs= (%a, %a)\n", got._limbs[0],
                 got._limbs[1]);
  }
}

// =============================================================================
// Tolerances
// =============================================================================
//
// MultiFloat<double, 2> has ~104 bits of precision. The kernels lose a few
// bits to ordering/normalization, so we accept up to ~2^-100 ≈ 8e-31 for
// linear ops and ~2^-95 for division (one Newton step from 1/y[0]).

static constexpr double kAddTol = 1e-30;
static constexpr double kMulTol = 1e-30;
static constexpr double kDivTol = 1e-28;

// =============================================================================
// Inputs
// =============================================================================

// All inputs are constructed by rounding a __float128 value down to a
// normalized double-double. This guarantees that to_q() is the exact
// inverse of the construction, so the __float128 reference is exact for
// every operation.
static std::vector<MF2> sample_inputs() {
  std::vector<MF2> v;
  auto push_q = [&v](q_t value) { v.push_back(from_q(value)); };

  // Pure doubles.
  push_q((q_t)0);
  push_q((q_t)1);
  push_q((q_t)-1);
  push_q((q_t)2);
  push_q((q_t)0.5);
  push_q((q_t)3);
  push_q((q_t)0.1);
  push_q((q_t)-0.1);
  push_q((q_t)123456789.0);

  // Irrational doubles with a non-trivial lo limb (1/3, sqrt(2), pi, e).
  push_q((q_t)1 / (q_t)3);
  push_q(sqrtq((q_t)2));
  push_q(M_PIq);
  push_q(-M_Eq);

  // 1 + tiny perturbations exercising the lo limb.
  push_q((q_t)1 + scalbnq((q_t)1, -53));
  push_q((q_t)1 - scalbnq((q_t)1, -54));
  push_q((q_t)1 + scalbnq((q_t)1, -80));

  // Large and small magnitudes (still representable as a normalized
  // double-double whose total span fits within float128 precision).
  push_q(scalbnq((q_t)1, 100));
  push_q(scalbnq((q_t)1, -100));
  push_q(scalbnq(M_PIq, 50));
  push_q(scalbnq(-M_Eq, -50));

  return v;
}

// =============================================================================
// Construction / conversion / comparison (exact properties)
// =============================================================================

static void test_construction_and_conversion() {
  MF2 z;
  REQUIRE(z._limbs[0] == 0.0);
  REQUIRE(z._limbs[1] == 0.0);

  MF2 one(1.0);
  REQUIRE(one._limbs[0] == 1.0);
  REQUIRE(one._limbs[1] == 0.0);
  REQUIRE(static_cast<double>(one) == 1.0);
  REQUIRE(to_q(one) == (q_t)1);

  MF2 negz(-0.0);
  REQUIRE(static_cast<double>(negz) == 0.0);
  REQUIRE(std::signbit(negz._limbs[0]));
}

static void test_equality_and_ordering() {
  for (MF2 const &a : sample_inputs()) {
    for (MF2 const &b : sample_inputs()) {
      q_t qa = to_q(a);
      q_t qb = to_q(b);
      REQUIRE((a == b) == (qa == qb));
      REQUIRE((a != b) == (qa != qb));
      REQUIRE((a < b) == (qa < qb));
      REQUIRE((a > b) == (qa > qb));
      REQUIRE((a <= b) == (qa <= qb));
      REQUIRE((a >= b) == (qa >= qb));
    }
  }
}

// =============================================================================
// Arithmetic — compared against __float128
// =============================================================================

static void test_unary(Stats &stats) {
  for (MF2 const &x : sample_inputs()) {
    q_t qx = to_q(x);
    check_q("neg", -x, -qx, 0.0, stats);
    check_q("pos", +x, qx, 0.0, stats);
  }
}

static void test_addition(Stats &stats) {
  auto inputs = sample_inputs();
  for (MF2 const &a : inputs) {
    for (MF2 const &b : inputs) {
      q_t expected = to_q(a) + to_q(b);
      check_q("add", a + b, expected, kAddTol, stats, &a, &b);

      MF2 c = a;
      c += b;
      check_q("add+=", c, expected, kAddTol, stats, &a, &b);
    }
  }
}

static void test_subtraction(Stats &stats) {
  auto inputs = sample_inputs();
  for (MF2 const &a : inputs) {
    for (MF2 const &b : inputs) {
      q_t expected = to_q(a) - to_q(b);
      check_q("sub", a - b, expected, kAddTol, stats);

      MF2 c = a;
      c -= b;
      check_q("sub-=", c, expected, kAddTol, stats);
    }
  }
}

static void test_multiplication(Stats &stats) {
  auto inputs = sample_inputs();
  for (MF2 const &a : inputs) {
    for (MF2 const &b : inputs) {
      q_t expected = to_q(a) * to_q(b);
      check_q("mul", a * b, expected, kMulTol, stats);

      MF2 c = a;
      c *= b;
      check_q("mul*=", c, expected, kMulTol, stats);
    }
  }
}

static void test_division(Stats &stats) {
  auto inputs = sample_inputs();
  for (MF2 const &a : inputs) {
    for (MF2 const &b : inputs) {
      if (b._limbs[0] == 0.0) {
        continue;
      }
      q_t expected = to_q(a) / to_q(b);
      check_q("div", a / b, expected, kDivTol, stats);

      MF2 c = a;
      c /= b;
      check_q("div/=", c, expected, kDivTol, stats);
    }
  }
}

// =============================================================================
// cmath-style wrappers
// =============================================================================

static void test_abs_fmin_fmax(Stats &stats) {
  auto inputs = sample_inputs();
  for (MF2 const &x : inputs) {
    q_t qx = to_q(x);
    q_t qexp = qx < 0 ? -qx : qx;
    check_q("abs", mf::abs(x), qexp, 0.0, stats);
    check_q("fabs", mf::fabs(x), qexp, 0.0, stats);
  }
  for (MF2 const &a : inputs) {
    for (MF2 const &b : inputs) {
      q_t qa = to_q(a);
      q_t qb = to_q(b);
      check_q("fmin", mf::fmin(a, b), qa < qb ? qa : qb, 0.0, stats);
      check_q("fmax", mf::fmax(a, b), qa < qb ? qb : qa, 0.0, stats);
    }
  }
}

static void test_copysign(Stats &stats) {
  auto inputs = sample_inputs();
  for (MF2 const &x : inputs) {
    for (MF2 const &y : inputs) {
      q_t qx = to_q(x);
      bool xs = std::signbit(x._limbs[0]);
      bool ys = std::signbit(y._limbs[0]);
      q_t qexp = (xs == ys) ? qx : -qx;
      // Match the wrapper's exact bit-pattern semantics on the limb.
      check_q("copysign", mf::copysign(x, y), qexp, 0.0, stats);
    }
  }
}

static void test_classification() {
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  MF2 finite(1.5);
  MF2 zero(0.0);
  MF2 neg(-1.5);
  MF2 pinf(inf), ninf(-inf), mnan(nan);

  REQUIRE(!mf::signbit(finite));
  REQUIRE(mf::signbit(neg));
  REQUIRE(mf::isfinite(finite));
  REQUIRE(!mf::isfinite(pinf));
  REQUIRE(!mf::isfinite(mnan));
  REQUIRE(mf::isinf(pinf));
  REQUIRE(mf::isinf(ninf));
  REQUIRE(!mf::isinf(finite));
  REQUIRE(mf::isnan(mnan));
  REQUIRE(!mf::isnan(finite));
  REQUIRE(mf::fpclassify(finite) == FP_NORMAL);
  REQUIRE(mf::fpclassify(zero) == FP_ZERO);
  REQUIRE(mf::fpclassify(pinf) == FP_INFINITE);
  REQUIRE(mf::fpclassify(mnan) == FP_NAN);

  // Non-canonical DDs where hi is (signed) zero: the value lives in lo,
  // so signbit/abs must consult lo rather than blindly trusting hi.
  MF2 tiny_neg;
  tiny_neg._limbs[0] = 0.0;
  tiny_neg._limbs[1] = -1e-300;
  REQUIRE(mf::signbit(tiny_neg));
  REQUIRE(!mf::signbit(mf::abs(tiny_neg)));

  MF2 tiny_pos;
  tiny_pos._limbs[0] = -0.0;
  tiny_pos._limbs[1] = 1e-300;
  REQUIRE(!mf::signbit(tiny_pos));
  REQUIRE(!mf::signbit(mf::abs(tiny_pos)));
}

static void test_ldexp_scalbn_ilogb(Stats &stats) {
  auto inputs = sample_inputs();
  int const ns[] = {-50, -1, 0, 1, 25};
  for (MF2 const &x : inputs) {
    if (!mf::isfinite(x) || x == MF2(0.0)) {
      continue;
    }
    q_t qx = to_q(x);
    for (int n : ns) {
      q_t qexp = qx * scalbnq((q_t)1, n);
      check_q("ldexp", mf::ldexp(x, n), qexp, 0.0, stats);
      check_q("scalbn", mf::scalbn(x, n), qexp, 0.0, stats);
    }
    // mf::ilogb is documented as delegating to the leading limb (matches
    // Julia/Fortran exponent semantics).
    REQUIRE(mf::ilogb(x) == std::ilogb(x._limbs[0]));
  }
}

// =============================================================================
// Driver
// =============================================================================

static void print_stats(char const *name, Stats const &s) {
  if (s.count == 0) {
    std::printf("  %-12s n=%d\n", name, s.count);
    return;
  }
  std::printf("  %-12s n=%-6d max_rel=%.3e  mean_rel=%.3e\n", name, s.count,
              s.max_rel, s.sum_rel / s.count);
}

int main() {
  test_construction_and_conversion();
  test_equality_and_ordering();
  test_classification();

  Stats add, sub, mul, div, una, abs_fmm, csgn, ldx;
  test_unary(una);
  test_addition(add);
  test_subtraction(sub);
  test_multiplication(mul);
  test_division(div);
  test_abs_fmin_fmax(abs_fmm);
  test_copysign(csgn);
  test_ldexp_scalbn_ilogb(ldx);

  std::printf("[multifloats_test] %d checks, %d failures\n", g_checks,
              g_failures);
  std::printf("[multifloats_test] precision report (relative error vs "
              "__float128):\n");
  print_stats("unary", una);
  print_stats("add", add);
  print_stats("sub", sub);
  print_stats("mul", mul);
  print_stats("div", div);
  print_stats("abs/fmin/fmax", abs_fmm);
  print_stats("copysign", csgn);
  print_stats("ldexp/etc", ldx);

  return g_failures == 0 ? 0 : 1;
}
