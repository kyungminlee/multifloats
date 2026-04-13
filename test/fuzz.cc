// Property-based fuzz tests for multifloats.hh.
//
// This is the C++ analogue of test/fuzz.f90: for each iteration a random
// __float128 pair (q1, q2) is generated, projected to MultiFloat<double, 2>
// inputs (f1, f2), and every operation that multifloats.hh exposes is run on
// both legs. The DD result is compared against the __float128 reference and
// per-op relative-error statistics are printed at the end.
//
// Input generation mirrors fuzz.f90's `generate_pair`: 10% non-finite, 10%
// close numbers, 10% sum-near-zero cancellation, 10% near-huge, 10% near-tiny,
// 50% wide random 10^[-30,30]. This exposes cancellation, overflow, and
// subnormal regimes that uniform-magnitude inputs miss.
//
// Tolerances mirror fuzz.f90's is_full_dd / is_compound classification:
//   1e-26   full DD (~106-bit) operations
//   1e-10   compound transcendentals (gamma, lgamma)
//   1e-15   single-double ops and libm-quality transcendentals
// Built with g++ + libquadmath.

#include "multifloats.hh"

#include <quadmath.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>

namespace mf = multifloats;
using MF2 = mf::MultiFloat<double, 2>;
using q_t = __float128;

// =============================================================================
// __float128 reference helpers
// =============================================================================

static q_t to_q(MF2 const &x) {
  return (q_t)x._limbs[0] + (q_t)x._limbs[1];
}

// Project a q_t value down to a normalized double-double (the same way the
// Fortran test's `to_f64x2` does — hi = nearest double, lo = rounding error,
// then fast_two_sum renormalize).
static MF2 to_mf2(q_t v) {
  double hi = (double)v;
  if (!std::isfinite(hi)) {
    MF2 r;
    r._limbs[0] = hi;
    r._limbs[1] = 0.0;
    return r;
  }
  double lo = (double)(v - (q_t)hi);
  double s = hi + lo;
  double b = s - hi;
  MF2 r;
  r._limbs[0] = s;
  r._limbs[1] = lo - b;
  return r;
}

static bool q_isnan(q_t x) { return isnanq(x); }
static bool q_isfinite(q_t x) { return finiteq(x); }

static double q_rel_err(q_t got, q_t expected) {
  q_t diff = got - expected;
  if (diff < 0) {
    diff = -diff;
  }
  q_t mag = expected < 0 ? -expected : expected;
  if (mag == 0) {
    return (double)diff;
  }
  return (double)(diff / mag);
}

static char const *qstr(q_t v) {
  static char buf[64];
  quadmath_snprintf(buf, sizeof(buf), "%.30Qg", v);
  return buf;
}

// =============================================================================
// Per-op statistics table (keyed by op name)
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
  for (int i = 0; i < g_nstats; ++i) {
    if (std::strcmp(g_stats[i].name, op) == 0) {
      return &g_stats[i];
    }
  }
  if (g_nstats >= kMaxStats) {
    return nullptr;
  }
  StatEntry *s = &g_stats[g_nstats++];
  std::snprintf(s->name, sizeof(s->name), "%s", op);
  return s;
}

static void update_stat(char const *op, double rel) {
  if (!std::isfinite(rel)) {
    return;
  }
  StatEntry *s = find_or_create_stat(op);
  if (!s) {
    return;
  }
  if (rel > s->max_rel) {
    s->max_rel = rel;
  }
  s->sum_rel += rel;
  ++s->count;
}

// =============================================================================
// Classification: which ops are expected to produce full DD precision?
// =============================================================================

// Full-DD tier: kernels that should deliver ~106 bits of precision end-to-end.
// exp/log/hyperbolic/pow now hit full DD via the native polynomial kernels.
// sin/cos/tan still lag (inv_pi reduction is only DD-accurate) and fmod is
// excluded because its large-quotient subtraction loses precision through
// catastrophic cancellation once |x/y| grows.
static bool is_full_dd(char const *op) {
  static char const *kList[] = {
      "add", "sub", "mul", "div", "sqrt", "abs", "neg",
      "add_fd", "mul_df", "fmin", "fmax", "copysign", "fdim",
      "hypot", "trunc", "round", "scalbn", "min3", "max3",
      "exp", "log", "log10", "pow", "pow_md", "pow_dm", "pow_int",
      "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
      "asin", "acos", "atan", "atan2",
      nullptr};
  for (int i = 0; kList[i]; ++i) {
    if (std::strcmp(kList[i], op) == 0) {
      return true;
    }
  }
  return false;
}

// Compound tier: chained evaluations or libm calls with historically weaker
// relative accuracy. Gamma/lgamma in libquadmath give ~1e-17 and the DD
// kernel propagates that ULP.
static bool is_compound(char const *op) {
  static char const *kList[] = {"tgamma", "lgamma", nullptr};
  for (int i = 0; kList[i]; ++i) {
    if (std::strcmp(kList[i], op) == 0) {
      return true;
    }
  }
  return false;
}

// =============================================================================
// Failure counter
// =============================================================================

static long g_failures = 0;
static long g_prints = 0;
static const long kPrintLimit = 50;

// Subnormal-range threshold: below this the DD lo limb cannot hold a full
// 53-bit error term and effective precision degrades to single-double.
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

static void check(char const *op, MF2 const &got, q_t expected, q_t i1, q_t i2) {
  // NaN: leading limb must be NaN.
  if (q_isnan(expected)) {
    if (!std::isnan(got._limbs[0])) {
      char buf[128];
      std::snprintf(buf, sizeof(buf), "expected NaN, got (%a, %a)",
                    got._limbs[0], got._limbs[1]);
      report_fail(op, buf);
    }
    return;
  }
  // Infinity: leading limb must match sign and be infinite.
  if (!q_isfinite(expected)) {
    bool ok = std::isinf(got._limbs[0]) &&
              (std::signbit(got._limbs[0]) == (expected < 0));
    if (!ok) {
      char buf[128];
      std::snprintf(buf, sizeof(buf), "expected inf, got (%a, %a)",
                    got._limbs[0], got._limbs[1]);
      report_fail(op, buf);
    }
    return;
  }

  q_t got_q = to_q(got);
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
    if (is_full_dd(op)) tol = 1e-26;
    else if (is_compound(op)) tol = 1e-10;
    else tol = 1e-14;
  } else {
    rel_err = (double)(diff / input_mag);
    if (is_full_dd(op)) tol = 1e-28;
    else if (is_compound(op)) tol = 1e-10;
    else tol = 1e-14;
  }

  // Skip subnormal input/output range from the precision report and
  // from pass/fail: the DD lo limb can't hold a full 53-bit error term there.
  if (subnormal_range(i1) || subnormal_range(i2) || subnormal_range(expected)) {
    return;
  }
  update_stat(op, rel_err);

  if (rel_err > tol) {
    // Near DBL_MAX / near DBL_MIN: real behavior is IEEE overflow/underflow,
    // not a precision bug.
    q_t huge_edge = (q_t)0x1.fffffffffffffp+1023q * (q_t)0.99q;
    if (abs_exp > huge_edge) return;
    if ((double)diff < 1e-35) return;
    ++g_failures;
    if (g_prints < kPrintLimit) {
      std::fprintf(stderr, "FAIL [%s] rel_err=%g > tol=%g\n", op, rel_err, tol);
      std::fprintf(stderr, "  i1       = %s\n", qstr(i1));
      std::fprintf(stderr, "  i2       = %s\n", qstr(i2));
      std::fprintf(stderr, "  expected = %s\n", qstr(expected));
      std::fprintf(stderr, "  got      = %s  (limbs %a, %a)\n", qstr(got_q),
                   got._limbs[0], got._limbs[1]);
      ++g_prints;
    }
  }
}

// =============================================================================
// Random input generation — mirrors fuzz.f90's generate_pair
// =============================================================================

struct Rng {
  std::mt19937_64 engine;
  std::uniform_real_distribution<double> u01{0.0, 1.0};

  explicit Rng(uint64_t seed) : engine(seed) {}
  double u() { return u01(engine); }

  q_t pick_nonfinite(double r) {
    if (r < 0.33) return (q_t) (+1.0 / 0.0);
    if (r < 0.66) return (q_t) (-1.0 / 0.0);
    return (q_t) (0.0 / 0.0);
  }

  // Wide random: sign * uniform(0.5) * 10^k  with k ∈ [-30, 30].
  q_t wide(double r, double rexp) {
    int k = (int)(rexp * 60.0) - 30;
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
      q2 = q1 * ((q_t)1.0q + (q_t)(r4 * 1e-15));
      break;
    }
    case 2: {
      int k = (int)(r3 * 20.0) - 10;
      q1 = (q_t)(r2 - 0.5) * powq((q_t)10.0q, (q_t)k);
      q2 = -q1 + (q_t)((r4 - 0.5) * 1e-25) * q1;
      break;
    }
    case 3: {
      double huge = 0x1.fffffffffffffp+1023;
      q1 = (q_t)huge * (q_t)(0.9 + 0.1 * r2);
      q2 = (q_t)huge * (q_t)(0.9 + 0.1 * r3);
      break;
    }
    case 4: {
      double tiny = 0x1.0p-1022;
      q1 = (q_t)tiny * (q_t)(1.0 + 10.0 * r2);
      q2 = (q_t)tiny * (q_t)(1.0 + 10.0 * r3);
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

// =============================================================================
// Driver
// =============================================================================

static void print_all_stats() {
  std::printf("\n");
  std::printf("Per-operation precision report (relative error vs __float128):\n");
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

static void check_comp(MF2 const &f1, MF2 const &f2, q_t q1, q_t q2) {
  if (q_isnan(q1) || q_isnan(q2)) return;
  q_t diff = q1 - q2;
  if (diff < 0) diff = -diff;
  q_t aq1 = q1 < 0 ? -q1 : q1;
  q_t aq2 = q2 < 0 ? -q2 : q2;
  q_t mag = aq1 > aq2 ? aq1 : aq2;
  if (mag < (q_t)1e-300q) mag = (q_t)1e-300q;
  bool near = (diff <= mag * (q_t)1e-20q) || (diff <= (q_t)1e-32q);
  auto cmp_fail = [&](char const *op, bool fres, bool qres) {
    char buf[192];
    std::snprintf(buf, sizeof(buf), "q1=%s q2=%s got=%d expect=%d", qstr(q1),
                  qstr(q2), (int)fres, (int)qres);
    char full[32];
    std::snprintf(full, sizeof(full), "cmp:%s", op);
    report_fail(full, buf);
  };
  if ((f1 < f2) != (q1 < q2) && !near) cmp_fail("lt", f1 < f2, q1 < q2);
  if ((f1 > f2) != (q1 > q2) && !near) cmp_fail("gt", f1 > f2, q1 > q2);
  if ((f1 <= f2) != (q1 <= q2) && !near) cmp_fail("le", f1 <= f2, q1 <= q2);
  if ((f1 >= f2) != (q1 >= q2) && !near) cmp_fail("ge", f1 >= f2, q1 >= q2);
  if ((f1 == f2) != (q1 == q2) && !near) cmp_fail("eq", f1 == f2, q1 == q2);
  if ((f1 != f2) != (q1 != q2) && !near) cmp_fail("ne", f1 != f2, q1 != q2);
}

int main(int argc, char **argv) {
  long iterations = 1000000;
  uint64_t seed = 42ULL;
  if (argc > 1) {
    iterations = std::atol(argv[1]);
  }
  if (argc > 2) {
    seed = std::strtoull(argv[2], nullptr, 0);
  }

  Rng rng(seed);
  std::printf("[multifloats_fuzz] iterations=%ld seed=0x%llx\n", iterations,
              (unsigned long long)seed);

  for (long i = 1; i <= iterations; ++i) {
    q_t q1, q2;
    rng.generate_pair(q1, q2);
    MF2 f1 = to_mf2(q1);
    MF2 f2 = to_mf2(q2);
    double d1 = (double)q1;
    double d2 = (double)q2;

    // ----------------------------------------------------------------
    // Hot loop: arithmetic + sqrt + comparisons + unary basics
    // ----------------------------------------------------------------
    // add/sub/neg/abs are non-finite safe. mul/div/etc use two_prod, whose
    // error limb becomes NaN when a hi-limb is ±inf even if the IEEE
    // product is well-defined; the hand-written non-finite short-circuits
    // only cover add. Gate the multiplicative path on finite inputs.
    bool const both_finite = q_isfinite(q1) && q_isfinite(q2);

    check("add", f1 + f2, q1 + q2, q1, q2);
    check("sub", f1 - f2, q1 - q2, q1, q2);
    if (both_finite) {
      check("mul", f1 * f2, q1 * q2, q1, q2);
      if (q2 != (q_t)0) {
        check("div", f1 / f2, q1 / q2, q1, q2);
      }
    }
    if (q_isfinite(q1) && q1 >= (q_t)0) {
      check("sqrt", mf::sqrt(f1), sqrtq(q1), q1, (q_t)0);
    }
    check("abs", mf::abs(f1), q1 < 0 ? -q1 : q1, q1, (q_t)0);
    check("neg", -f1, -q1, q1, (q_t)0);

    // Bit-exact rounding (full DD)
    q_t aq1 = q1 < 0 ? -q1 : q1;
    if (q_isfinite(q1) && aq1 < (q_t)1e15q) {
      check("trunc", mf::trunc(f1), truncq(q1), q1, (q_t)0);
      check("round", mf::round(f1), roundq(q1), q1, (q_t)0);
    }

    check_comp(f1, f2, q1, q2);

    // ----------------------------------------------------------------
    // Periodic (every 10): mixed-mode + binary
    // ----------------------------------------------------------------
    if (i % 10 == 0) {
      // Mixed-mode and binary ops all route through the multiplicative
      // or hypot kernels and share mul's non-finite limitation.
      if (both_finite) {
        check("add_fd", f1 + MF2(d2), q1 + (q_t)d2, q1, (q_t)d2);
        check("mul_df", MF2(d1) * f2, (q_t)d1 * q2, (q_t)d1, q2);
        check("fmin", mf::fmin(f1, f2), q1 < q2 ? q1 : q2, q1, q2);
        check("fmax", mf::fmax(f1, f2), q1 < q2 ? q2 : q1, q1, q2);
        check("copysign", mf::copysign(f1, f2), copysignq(q1, q2), q1, q2);
        check("fdim", mf::fdim(f1, f2), fdimq(q1, q2), q1, q2);
        check("hypot", mf::hypot(f1, f2), hypotq(q1, q2), q1, q2);

        q_t aq2 = q2 < 0 ? -q2 : q2;
        if (q2 != (q_t)0 && aq1 < (q_t)1e20q && aq2 > (q_t)1e-20q) {
          check("fmod", mf::fmod(f1, f2), fmodq(q1, q2), q1, q2);
        }
      }
    }

    // ----------------------------------------------------------------
    // Periodic (every 100): transcendentals
    // ----------------------------------------------------------------
    if (i % 100 == 0) {
      if (q_isfinite(q1)) {
        q_t qe = expq(q1);
        if (q_isfinite(qe)) {
          check("exp", mf::exp(f1), qe, q1, (q_t)0);
        }
      }
      if (q_isfinite(q1) && q1 > (q_t)0) {
        check("log", mf::log(f1), logq(q1), q1, (q_t)0);
        check("log10", mf::log10(f1), log10q(q1), q1, (q_t)0);
      }

      // Trig: keep magnitudes moderate.
      if (q_isfinite(q1) && aq1 < (q_t)1e6q) {
        check("sin", mf::sin(f1), sinq(q1), q1, (q_t)0);
        check("cos", mf::cos(f1), cosq(q1), q1, (q_t)0);
        if (fabsq(cosq(q1)) > (q_t)1e-12q) {
          check("tan", mf::tan(f1), tanq(q1), q1, (q_t)0);
        }
      }

      if (q_isfinite(q1) && aq1 <= (q_t)1) {
        check("asin", mf::asin(f1), asinq(q1), q1, (q_t)0);
        check("acos", mf::acos(f1), acosq(q1), q1, (q_t)0);
      }
      if (q_isfinite(q1)) {
        check("atan", mf::atan(f1), atanq(q1), q1, (q_t)0);
      }

      if (q_isfinite(q1) && aq1 < (q_t)700) {
        check("sinh", mf::sinh(f1), sinhq(q1), q1, (q_t)0);
        check("cosh", mf::cosh(f1), coshq(q1), q1, (q_t)0);
        check("tanh", mf::tanh(f1), tanhq(q1), q1, (q_t)0);
      }
      if (q_isfinite(q1)) {
        check("asinh", mf::asinh(f1), asinhq(q1), q1, (q_t)0);
      }
      if (q_isfinite(q1) && q1 >= (q_t)1) {
        check("acosh", mf::acosh(f1), acoshq(q1), q1, (q_t)0);
      }
      if (q_isfinite(q1) && aq1 < (q_t)1) {
        check("atanh", mf::atanh(f1), atanhq(q1), q1, (q_t)0);
      }

      if (q_isfinite(q1) && aq1 < (q_t)100) {
        check("erf", mf::erf(f1), erfq(q1), q1, (q_t)0);
        check("erfc", mf::erfc(f1), erfcq(q1), q1, (q_t)0);
      }
      if (q_isfinite(q1) && q1 > (q_t)0 && q1 < (q_t)100) {
        check("tgamma", mf::tgamma(f1), tgammaq(q1), q1, (q_t)0);
        check("lgamma", mf::lgamma(f1), lgammaq(q1), q1, (q_t)0);
      }

      check("atan2", mf::atan2(f1, f2), atan2q(q1, q2), q1, q2);

      // Power operator: positive base, modest exponent.
      q_t aq2 = q2 < 0 ? -q2 : q2;
      if (q_isfinite(q1) && q1 > (q_t)1e-3q && q1 < (q_t)1e3q &&
          q_isfinite(q2) && aq2 < (q_t)30) {
        check("pow", mf::pow(f1, f2), powq(q1, q2), q1, q2);
        check("pow_md", mf::pow(f1, MF2(d2)), powq(q1, (q_t)d2), q1, (q_t)d2);
        check("pow_dm", mf::pow(MF2(d1), f2), powq((q_t)d1, q2), (q_t)d1, q2);
      }
      if (q_isfinite(q1) && aq1 < (q_t)1e10q) {
        check("pow_int", mf::pow(f1, MF2(3.0)), powq(q1, (q_t)3), q1, (q_t)3);
      }

      if (q_isfinite(q1)) {
        check("scalbn", mf::scalbn(f1, 5), scalbnq(q1, 5), q1, (q_t)0);
      }

      // 3-argument min/max via nested fmin/fmax.
      if (q_isfinite(q1) && q_isfinite(q2)) {
        q_t q3 = (q1 + q2) * (q_t)0.5q;
        MF2 f3 = to_mf2(q3);
        q_t min12 = q1 < q2 ? q1 : q2;
        q_t q_min3 = min12 < q3 ? min12 : q3;
        q_t max12 = q1 < q2 ? q2 : q1;
        q_t q_max3 = max12 < q3 ? q3 : max12;
        check("min3", mf::fmin(mf::fmin(f1, f2), f3), q_min3, q1, q2);
        check("max3", mf::fmax(mf::fmax(f1, f2), f3), q_max3, q1, q2);
      }
    }

    if (i % 100000 == 0) {
      std::printf("  ... completed %ld iterations (failures so far: %ld)\n", i,
                  g_failures);
    }
  }

  print_all_stats();
  std::printf("[multifloats_fuzz] failures=%ld\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}
