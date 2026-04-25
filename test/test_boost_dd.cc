// Precision fuzz for boost::multiprecision::cpp_double_double, using the
// same __float128 (libquadmath) oracle as test/fuzz.cc. The point is to put
// boost's DD backend on the same scale as the multifloats kernels so the
// two `[*_fuzz] precision report` tables can be diffed side-by-side.
//
// Op coverage is restricted to what the boost backend (and Boost.Math)
// actually provides: arithmetic + standard libm transcendentals + erf /
// gamma. Bessel, complex, π-scaled trig, erfcx, array kernels — all out
// of scope (no native boost equivalent).

#include "test_common.hh"

#include <boost/multiprecision/cpp_double_fp.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>

#include <quadmath.h>

using boost::multiprecision::cpp_double_double;
using multifloats_test::q_t;
using multifloats_test::to_q;
using multifloats_test::from_q;
using multifloats_test::q_isfinite;
using multifloats_test::q_isnan;
using multifloats_test::qstr;

// =============================================================================
// boost <-> qp adapters. We renormalize through the existing fast_two_sum
// path in test_common.hh so the boost limbs satisfy the same |lo| <= 0.5 ULP
// of |hi| invariant the multifloats DD invariant maintains. Equivalent to
// constructing `cpp_double_double(hi) + cpp_double_double(lo)` (boost
// renormalizes on the constructor sum) but ~free.
// =============================================================================

static inline q_t bdd_to_q(cpp_double_double const &x) {
  auto const &r = x.backend().crep();
  return (q_t)r.first + (q_t)r.second;
}

static inline cpp_double_double bdd_from_q(q_t v) {
  multifloats::float64x2 mf = from_q(v);
  cpp_double_double out;
  auto &r = out.backend().rep();
  r.first  = mf.limbs[0];
  r.second = mf.limbs[1];
  return out;
}

static inline cpp_double_double bdd_from_d(double v) {
  cpp_double_double out;
  auto &r = out.backend().rep();
  r.first  = v;
  r.second = 0.0;
  return out;
}

// =============================================================================
// Per-op stats, classifier, failure printer — same shape as fuzz.cc.
// =============================================================================

struct StatEntry {
  char name[24] = {};
  double max_rel = 0.0;
  double sum_rel = 0.0;
  long count = 0;
};

static constexpr int kMaxStats = 256;
static StatEntry g_stats[kMaxStats];
static int g_nstats = 0;

static StatEntry *find_or_create_stat(char const *op) {
  for (int i = 0; i < g_nstats; ++i)
    if (std::strcmp(g_stats[i].name, op) == 0) return &g_stats[i];
  if (g_nstats >= kMaxStats) return nullptr;
  StatEntry *s = &g_stats[g_nstats++];
  std::snprintf(s->name, sizeof(s->name), "%s", op);
  return s;
}

static void update_stat(char const *op, double rel) {
  if (!std::isfinite(rel)) return;
  StatEntry *s = find_or_create_stat(op);
  if (!s) return;
  if (rel > s->max_rel) s->max_rel = rel;
  s->sum_rel += rel;
  ++s->count;
}

// Same tolerance tiers as fuzz.cc. Boost's generic transcendentals on top
// of the DD backend may not always meet the full-DD bar (~1e-26); we still
// classify them as full-DD here so the report shows whether they do — the
// pass/fail print is informational, not a CI gate.
static bool is_full_dd(char const *op) {
  static char const *kList[] = {
      "add", "sub", "mul", "div", "sqrt", "abs", "neg", "fma",
      "fmin", "fmax", "fdim", "copysign", "hypot", "fmod",
      "trunc", "round", "floor", "ceil",
      "exp", "exp2", "expm1", "log", "log2", "log10", "log1p", "pow",
      "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
      "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
      "erf", "erfc", "tgamma", "lgamma",
      nullptr};
  for (int i = 0; kList[i]; ++i)
    if (std::strcmp(kList[i], op) == 0) return true;
  return false;
}

static long g_failures = 0;
static long g_prints = 0;
static const long kPrintLimit = 50;
static constexpr double kSubnormalFloor = 1.0e-290;

static bool subnormal_range(q_t x) {
  q_t a = x < 0 ? -x : x;
  return a > (q_t)0 && a < (q_t)kSubnormalFloor;
}

static void report_fail(char const *op, char const *detail) {
  ++g_failures;
  if (g_prints < kPrintLimit) {
    std::fprintf(stderr, "FAIL [%s] %s\n", op, detail);
    ++g_prints;
  }
}

static void check(char const *op, cpp_double_double const &got, q_t expected,
                  q_t i1, q_t i2) {
  auto const &g = got.backend().crep();
  if (q_isnan(expected)) {
    if (!std::isnan(g.first)) {
      char buf[128];
      std::snprintf(buf, sizeof(buf), "expected NaN, got (%a, %a)",
                    g.first, g.second);
      report_fail(op, buf);
    }
    return;
  }
  if (!q_isfinite(expected)) {
    bool ok = std::isinf(g.first) && (std::signbit(g.first) == (expected < 0));
    if (!ok) {
      char buf[128];
      std::snprintf(buf, sizeof(buf), "expected inf, got (%a, %a)",
                    g.first, g.second);
      report_fail(op, buf);
    }
    return;
  }

  q_t got_q = bdd_to_q(got);
  q_t diff = got_q - expected;
  if (diff < 0) diff = -diff;

  q_t abs_i1 = i1 < 0 ? -i1 : i1;
  q_t abs_i2 = i2 < 0 ? -i2 : i2;
  q_t input_mag = abs_i1 > abs_i2 ? abs_i1 : abs_i2;
  if (input_mag < (q_t)1e-300q) input_mag = (q_t)1e-300q;

  q_t abs_exp = expected < 0 ? -expected : expected;
  double rel_err;
  double tol;
  if (abs_exp > input_mag * (q_t)1e-10q) {
    rel_err = (double)(diff / abs_exp);
    tol = is_full_dd(op) ? 1e-26 : 1e-15;
  } else {
    rel_err = (double)(diff / input_mag);
    tol = is_full_dd(op) ? 1e-28 : 1e-15;
  }

  if (subnormal_range(i1) || subnormal_range(i2) || subnormal_range(expected))
    return;

  update_stat(op, rel_err);

  if (rel_err > tol) {
    q_t huge_edge = (q_t)0x1.fffffffffffffp+1023q * (q_t)0.99q;
    if (abs_exp > huge_edge) return;
    ++g_failures;
    if (g_prints < kPrintLimit) {
      std::fprintf(stderr, "FAIL [%s] rel_err=%g > tol=%g\n", op, rel_err, tol);
      std::fprintf(stderr, "  i1       = %s\n", qstr(i1));
      std::fprintf(stderr, "  i2       = %s\n", qstr(i2));
      std::fprintf(stderr, "  expected = %s\n", qstr(expected));
      std::fprintf(stderr, "  got      = %s  (limbs %a, %a)\n", qstr(got_q),
                   g.first, g.second);
      ++g_prints;
    }
  }
}

// =============================================================================
// Random input generator — mirrors test/fuzz.cc Rng (verbatim).
// =============================================================================

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
  q_t wide(double r, double rexp) {
    int k = (int)(rexp * 60.0) - 30;
    q_t mag = powq((q_t)10.0q, (q_t)k);
    return (q_t)(r - 0.5) * mag;
  }
  void generate_pair(q_t &q1, q_t &q2) {
    double r1 = u(), r2 = u(), r3 = u(), r4 = u();
    int mode = (int)(r1 * 10.0);
    switch (mode) {
    case 0:
      q1 = pick_nonfinite(r2); q2 = pick_nonfinite(r3); break;
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
      int k = (int)(r3 * 10.0) + 20;
      q1 = (q_t)(r2 - 0.5) * 2.0q * powq((q_t)10.0q, (q_t)k);
      q2 = (q_t)(r4 - 0.5) * 2.0q * powq((q_t)10.0q, (q_t)k);
      break;
    }
    case 4: {
      int k = -((int)(r3 * 10.0) + 20);
      q1 = (q_t)(r2 - 0.5) * 2.0q * powq((q_t)10.0q, (q_t)k);
      q2 = (q_t)(r4 - 0.5) * 2.0q * powq((q_t)10.0q, (q_t)k);
      break;
    }
    default:
      q1 = wide(r2, r3); q2 = wide(r4, r1); break;
    }
  }
};

// =============================================================================
// Driver
// =============================================================================

static void print_all_stats() {
  std::printf("\nPer-operation precision report (relative error vs __float128):\n");
  std::printf("  %-16s %10s %14s %14s\n", "op", "n", "max_rel", "mean_rel");
  for (int i = 0; i < g_nstats; ++i) {
    StatEntry const &s = g_stats[i];
    if (s.count == 0) {
      std::printf("  %-16s %10ld     (no data)\n", s.name, s.count);
      continue;
    }
    std::printf("  %-16s %10ld  %14.3e %14.3e\n", s.name, s.count, s.max_rel,
                s.sum_rel / s.count);
  }
  std::printf("\n");
}

int main(int argc, char **argv) {
  long iterations = 1000000;
  uint64_t seed = 0x42ULL;
  if (argc > 1) iterations = std::atol(argv[1]);
  if (argc > 2) seed = std::strtoull(argv[2], nullptr, 0);

  std::printf("[boost_dd_fuzz] iterations=%ld seed=0x%llx (cpp_double_double)\n",
              iterations, (unsigned long long)seed);

  Rng rng(seed);
  using std::sqrt; using std::fabs; using std::fma;
  using std::exp; using std::exp2; using std::expm1;
  using std::log; using std::log2; using std::log10; using std::log1p;
  using std::pow; using std::hypot; using std::fmod;
  using std::sin; using std::cos; using std::tan;
  using std::asin; using std::acos; using std::atan; using std::atan2;
  using std::sinh; using std::cosh; using std::tanh;
  using std::asinh; using std::acosh; using std::atanh;
  using std::floor; using std::ceil; using std::trunc; using std::round;
  using std::fmin; using std::fmax; using std::copysign;
  using boost::math::erf; using boost::math::erfc;
  using boost::math::tgamma; using boost::math::lgamma;

  for (long i = 1; i <= iterations; ++i) {
    q_t q1, q2;
    rng.generate_pair(q1, q2);
    // Round-trip through DD so qp matches what DD can represent exactly.
    cpp_double_double f1 = bdd_from_q(q1);
    cpp_double_double f2 = bdd_from_q(q2);
    q1 = bdd_to_q(f1);
    q2 = bdd_to_q(f2);

    bool both_finite = q_isfinite(q1) && q_isfinite(q2);
    q_t aq1 = q1 < 0 ? -q1 : q1;
    q_t aq2 = q2 < 0 ? -q2 : q2;

    check("add", f1 + f2, q1 + q2, q1, q2);
    check("sub", f1 - f2, q1 - q2, q1, q2);
    if (both_finite)            check("mul", f1 * f2, q1 * q2, q1, q2);
    if (both_finite && q2 != 0) check("div", f1 / f2, q1 / q2, q1, q2);
    if (q_isfinite(q1) && q1 >= 0)
      check("sqrt", sqrt(f1), sqrtq(q1), q1, (q_t)0);
    check("abs", fabs(f1), q1 < 0 ? -q1 : q1, q1, (q_t)0);
    check("neg", -f1, -q1, q1, (q_t)0);

    if (q_isfinite(q1) && aq1 < (q_t)1e15q) {
      check("trunc", trunc(f1), truncq(q1), q1, (q_t)0);
      check("round", round(f1), roundq(q1), q1, (q_t)0);
      check("floor", floor(f1), floorq(q1), q1, (q_t)0);
      check("ceil",  ceil(f1),  ceilq(q1),  q1, (q_t)0);
    }

    // Periodic mixed-mode + binary (every 10th iter).
    if (i % 10 == 0 && both_finite) {
      check("fmin",     fmin(f1, f2),     fminq(q1, q2),     q1, q2);
      check("fmax",     fmax(f1, f2),     fmaxq(q1, q2),     q1, q2);
      check("copysign", copysign(f1, f2), copysignq(q1, q2), q1, q2);
      check("fdim",
            f1 > f2 ? f1 - f2 : cpp_double_double(0),
            q1 > q2 ? q1 - q2 : (q_t)0, q1, q2);
      check("hypot",    hypot(f1, f2),    hypotq(q1, q2),    q1, q2);
      if (q2 != 0 && aq1 < (q_t)1e20q && aq2 > (q_t)1e-20q)
        check("fmod",   fmod(f1, f2),     fmodq(q1, q2),     q1, q2);
      // fma: c = third input, use 1.0 for a clean reference.
      cpp_double_double f3 = bdd_from_d(1.0);
      check("fma",      fma(f1, f2, f3),  fmaq(q1, q2, (q_t)1.0q), q1, q2);
    }

    // Periodic transcendentals (every 100th iter).
    if (i % 100 == 0) {
      if (q_isfinite(q1) && q_isfinite(expq(q1)))   check("exp",   exp(f1),   expq(q1),   q1, (q_t)0);
      if (q_isfinite(q1) && q_isfinite(exp2q(q1)))  check("exp2",  exp2(f1),  exp2q(q1),  q1, (q_t)0);
      if (q_isfinite(q1) && q_isfinite(expm1q(q1))) check("expm1", expm1(f1), expm1q(q1), q1, (q_t)0);
      if (q_isfinite(q1) && q1 > 0)  check("log",   log(f1),   logq(q1),   q1, (q_t)0);
      if (q_isfinite(q1) && q1 > 0)  check("log2",  log2(f1),  log2q(q1),  q1, (q_t)0);
      if (q_isfinite(q1) && q1 > 0)  check("log10", log10(f1), log10q(q1), q1, (q_t)0);
      if (q_isfinite(q1) && q1 > -1) check("log1p", log1p(f1), log1pq(q1), q1, (q_t)0);

      if (q_isfinite(q1) && aq1 < (q_t)1e6q) {
        check("sin", sin(f1), sinq(q1), q1, (q_t)0);
        check("cos", cos(f1), cosq(q1), q1, (q_t)0);
        if (fabsq(cosq(q1)) > (q_t)1e-12q)
          check("tan", tan(f1), tanq(q1), q1, (q_t)0);
      }
      if (q_isfinite(q1) && aq1 <= 1) {
        check("asin", asin(f1), asinq(q1), q1, (q_t)0);
        check("acos", acos(f1), acosq(q1), q1, (q_t)0);
      }
      if (q_isfinite(q1)) check("atan", atan(f1), atanq(q1), q1, (q_t)0);
      if (both_finite)    check("atan2", atan2(f1, f2), atan2q(q1, q2), q1, q2);

      if (q_isfinite(q1) && aq1 < (q_t)700) {
        check("sinh", sinh(f1), sinhq(q1), q1, (q_t)0);
        check("cosh", cosh(f1), coshq(q1), q1, (q_t)0);
        check("tanh", tanh(f1), tanhq(q1), q1, (q_t)0);
      }
      if (q_isfinite(q1)) check("asinh", asinh(f1), asinhq(q1), q1, (q_t)0);
      if (q_isfinite(q1) && q1 >= 1)
        check("acosh", acosh(f1), acoshq(q1), q1, (q_t)0);
      if (q_isfinite(q1) && aq1 < 1)
        check("atanh", atanh(f1), atanhq(q1), q1, (q_t)0);

      // pow: positive base, bounded exponent.
      if (q_isfinite(q1) && q1 > 0 && q_isfinite(q2) && aq2 < 50) {
        cpp_double_double fp = bdd_from_q(q1);
        cpp_double_double fe = bdd_from_q(q2);
        check("pow", pow(fp, fe), powq(q1, q2), q1, q2);
      }

      // erf / erfc: bounded input.
      if (q_isfinite(q1) && aq1 < 6) {
        check("erf",  erf(f1),  erfq(q1),  q1, (q_t)0);
        check("erfc", erfc(f1), erfcq(q1), q1, (q_t)0);
      }
      // tgamma / lgamma: positive bounded argument, away from poles.
      if (q_isfinite(q1) && q1 > (q_t)0.1q && q1 < (q_t)170) {
        check("tgamma", tgamma(f1), tgammaq(q1), q1, (q_t)0);
        check("lgamma", lgamma(f1), lgammaq(q1), q1, (q_t)0);
      }
    }
  }

  print_all_stats();
  // Informational only — boost's generic transcendentals on cpp_double_double
  // can fall short of the multifloats hand-tuned bar, and that delta is the
  // point of the report. Always exit 0 so this isn't a CI gate.
  std::printf("[boost_dd_fuzz] non-tolerance failures (informational)=%ld\n",
              g_failures);
  return 0;
}
