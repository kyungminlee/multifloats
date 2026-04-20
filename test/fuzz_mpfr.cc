// MPFR-based precision fuzz for multifloats.
//
// Per iteration: draw a random input pair, then for each op under test
// compute three values —
//   ref_mp :  mpreal at 200-bit precision (ground truth)
//   ref_q  :  libquadmath at 113-bit precision
//   got_dd :  multifloats DD (~106-bit)
// Report, per op, the max and mean of both rel_err(q vs mp) and
// rel_err(dd vs mp). This separates the DD kernel's true error from the
// float128 reference floor that scalar fuzz.cc cannot distinguish.
//
// Built only when -DBUILD_MPFR_TESTS=ON.

#include "multifloats.hh"
#include "test_common.hh"
#include "test_common_mpfr.hh"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>

namespace mf = multifloats;
using multifloats_test::q_t;
using multifloats_test::to_q;
using multifloats_test::from_q;
using multifloats_test::q_isfinite;
using multifloats_test::q_isnan;
using multifloats_test::mp_t;
using multifloats_test::to_mp;
using multifloats_test::mp_rel_err;
using multifloats_test::mp_isfinite;
using multifloats_test::mp_isnan;

// =============================================================================
// Per-op statistics — two columns, one for each reference comparison.
// =============================================================================

struct StatEntry {
  char name[24] = {};
  double max_q  = 0.0;  double sum_q  = 0.0;
  double max_dd = 0.0;  double sum_dd = 0.0;
  long count = 0;
};

static constexpr int kMaxStats = 128;
static StatEntry g_stats[kMaxStats];
static int g_nstats = 0;

static StatEntry *stat(char const *op) {
  for (int i = 0; i < g_nstats; ++i)
    if (std::strcmp(g_stats[i].name, op) == 0) return &g_stats[i];
  if (g_nstats >= kMaxStats) return nullptr;
  StatEntry *s = &g_stats[g_nstats++];
  std::snprintf(s->name, sizeof(s->name), "%s", op);
  return s;
}

static void record(char const *op, double rq, double rdd) {
  if (!std::isfinite(rq) || !std::isfinite(rdd)) return;
  StatEntry *s = stat(op);
  if (!s) return;
  if (rq  > s->max_q)  s->max_q  = rq;
  if (rdd > s->max_dd) s->max_dd = rdd;
  s->sum_q  += rq;
  s->sum_dd += rdd;
  ++s->count;
}

// Compare DD output and a float128 value against an mpreal reference.
// `op` labels the row; ref_mp must be finite for the sample to count.
//
// Skips three IEEE floors where reporting rel-err would be misleading:
//   1. ref_mp non-finite (NaN/Inf — no reference to compare).
//   2. |ref_mp| < 2^-1050 (deep subnormal; IEEE double rounds to zero
//      or the smallest subnormal, so rel-err up to 1.0 is expected and
//      is not a library bug).
//   3. Non-finite DD / quad outputs (overflow boundary).
static bool mp_below_subnormal_floor(mp_t const &v) {
  // DD precision is bounded by IEEE-double's exponent range: the lo limb
  // is itself a double, so it cannot go below 2^-1074. For the pair to
  // carry full 106-bit precision we need lo at ~hi·2^-53 to land in the
  // *normal* double range, i.e. hi ≥ 2^-969. Below that, DD loses bits
  // smoothly (at hi = 2^-1022 ≈ 53 bits of total precision survive).
  // Skipping samples with |ref| < 2^-969 keeps fuzz measuring kernel
  // quality rather than this format-level cliff. See doc/developer/AUDIT_TODO.md
  // (erfc deep-tail analysis, 2026-04-19) and compare erfcq, which has
  // no such cliff because __float128's exponent reaches ~2^-16494.
  static const mp_t floor_ = mpfr::pow(mp_t(2), mp_t(-969));
  return mpfr::abs(v) < floor_;
}

static void check(char const *op,
                  mf::float64x2 const &got_dd,
                  q_t ref_q,
                  mp_t const &ref_mp) {
  if (!mp_isfinite(ref_mp) || mp_isnan(ref_mp)) return;
  if (!q_isfinite(ref_q)) return;
  if (!std::isfinite(got_dd._limbs[0])) return;
  if (mp_below_subnormal_floor(ref_mp)) return;
  double rq  = mp_rel_err(to_mp(ref_q),          ref_mp);
  double rdd = mp_rel_err(to_mp(got_dd),          ref_mp);
  record(op, rq, rdd);
}

// =============================================================================
// Random input generation — finite inputs only.
// =============================================================================

struct Rng {
  std::mt19937_64 eng;
  std::uniform_real_distribution<double> u01{0.0, 1.0};
  explicit Rng(uint64_t s) : eng(s) {}
  double u() { return u01(eng); }

  q_t wide(double r, double rexp) {
    int k = (int)(rexp * 60.0) - 30;
    q_t mag = powq((q_t)10.0q, (q_t)k);
    return (q_t)(r - 0.5) * mag;
  }

  void generate_pair(q_t &q1, q_t &q2) {
    double r1 = u(), r2 = u(), r3 = u(), r4 = u();
    int mode = (int)(r1 * 8.0);
    switch (mode) {
    case 0: {  // near-equal (cancellation in sub)
      int k = (int)(r3 * 20.0) - 10;
      q1 = (q_t)(r2 - 0.5) * powq((q_t)10.0q, (q_t)k);
      q2 = q1 * ((q_t)1.0q + (q_t)(r4 * 1e-15));
      break;
    }
    case 1: {  // near-opposite (cancellation in add)
      int k = (int)(r3 * 20.0) - 10;
      q1 = (q_t)(r2 - 0.5) * powq((q_t)10.0q, (q_t)k);
      q2 = -q1 + (q_t)((r4 - 0.5) * 1e-25) * q1;
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
// Output
// =============================================================================

static void print_report() {
  std::printf("\n");
  std::printf("Per-operation 3-way precision report (reference = mpreal @ 200 bits):\n");
  std::printf("  %-10s %10s %14s %14s %14s %14s\n",
              "op", "n",
              "max_q",  "mean_q",
              "max_dd", "mean_dd");
  for (int i = 0; i < g_nstats; ++i) {
    StatEntry const &s = g_stats[i];
    if (s.count == 0) continue;
    double invn = 1.0 / (double)s.count;
    std::printf("  %-10s %10ld  %14.3e %14.3e %14.3e %14.3e\n",
                s.name, s.count,
                s.max_q,  s.sum_q  * invn,
                s.max_dd, s.sum_dd * invn);
  }
  std::printf("\n"
              "  Columns: q  = rel_err(libquadmath vs mpreal),\n"
              "           dd = rel_err(multifloats DD vs mpreal).\n"
              "  q ≈ 1e-33 is the float128 mantissa floor; dd above q means\n"
              "  the DD kernel — not the float128 reference — is the loss.\n");
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char **argv) {
  long iterations = 10000;
  uint64_t seed = 42ULL;
  if (argc > 1) iterations = std::atol(argv[1]);
  if (argc > 2) seed = std::strtoull(argv[2], nullptr, 0);

  mpfr::mpreal::set_default_prec(multifloats_test::kMpfrPrec);
  std::printf("[multifloats_fuzz_mpfr] iterations=%ld seed=0x%llx prec=%d bits\n",
              iterations, (unsigned long long)seed, (int)multifloats_test::kMpfrPrec);

  Rng rng(seed);

  for (long i = 1; i <= iterations; ++i) {
    q_t qseed1, qseed2;
    rng.generate_pair(qseed1, qseed2);
    if (!q_isfinite(qseed1) || !q_isfinite(qseed2)) continue;

    // Inputs live at DD precision: all three paths (DD / quad / mpreal)
    // must start from the SAME starting value or "error" is dominated by
    // input rounding, not by the op kernel under test.
    mf::float64x2 f1 = from_q(qseed1);
    mf::float64x2 f2 = from_q(qseed2);
    q_t q1 = to_q(f1);
    q_t q2 = to_q(f2);
    mp_t m1 = to_mp(f1);
    mp_t m2 = to_mp(f2);

    // ---------------- arithmetic ----------------
    check("add", f1 + f2, q1 + q2, m1 + m2);
    check("sub", f1 - f2, q1 - q2, m1 - m2);
    check("mul", f1 * f2, q1 * q2, m1 * m2);
    if (q2 != (q_t)0) {
      check("div", f1 / f2, q1 / q2, m1 / m2);
    }
    if (q1 >= (q_t)0) {
      check("sqrt", mf::sqrt(f1), sqrtq(q1), mpfr::sqrt(m1));
    }

    // ---------------- binary ops ----------------
    if (i % 5 == 0) {
      check("hypot", mf::hypot(f1, f2), hypotq(q1, q2), mpfr::hypot(m1, m2));
      check("fmin",  mf::fmin(f1, f2),  q1 < q2 ? q1 : q2, q1 < q2 ? m1 : m2);
      check("fmax",  mf::fmax(f1, f2),  q1 < q2 ? q2 : q1, q1 < q2 ? m2 : m1);
    }

    // ---------------- transcendentals (every 100) ----------------
    if (i % 100 != 0) continue;
    q_t aq1 = q1 < 0 ? -q1 : q1;

    // exp / log family
    {
      q_t qe = expq(q1);
      if (q_isfinite(qe)) {
        check("exp", mf::exp(f1), qe, mpfr::exp(m1));
      }
      q_t qem = expm1q(q1);
      if (q_isfinite(qem)) {
        check("expm1", mf::expm1(f1), qem, mpfr::expm1(m1));
      }
    }
    if (q1 > (q_t)0) {
      check("log",   mf::log(f1),   logq(q1),   mpfr::log(m1));
      check("log10", mf::log10(f1), log10q(q1), mpfr::log10(m1));
    }
    if (q1 > (q_t)-1) {
      check("log1p", mf::log1p(f1), log1pq(q1), mpfr::log1p(m1));
    }

    // trig
    if (aq1 < (q_t)1e6q) {
      check("sin", mf::sin(f1), sinq(q1), mpfr::sin(m1));
      check("cos", mf::cos(f1), cosq(q1), mpfr::cos(m1));
      if (fabsq(cosq(q1)) > (q_t)1e-12q) {
        check("tan", mf::tan(f1), tanq(q1), mpfr::tan(m1));
      }
    }
    if (aq1 <= (q_t)1) {
      check("asin", mf::asin(f1), asinq(q1), mpfr::asin(m1));
      check("acos", mf::acos(f1), acosq(q1), mpfr::acos(m1));
    }
    check("atan",  mf::atan(f1),      atanq(q1),      mpfr::atan(m1));
    check("atan2", mf::atan2(f1, f2), atan2q(q1, q2), mpfr::atan2(m1, m2));

    // hyperbolic
    if (aq1 < (q_t)700) {
      check("sinh", mf::sinh(f1), sinhq(q1), mpfr::sinh(m1));
      check("cosh", mf::cosh(f1), coshq(q1), mpfr::cosh(m1));
      check("tanh", mf::tanh(f1), tanhq(q1), mpfr::tanh(m1));
    }
    check("asinh", mf::asinh(f1), asinhq(q1), mpfr::asinh(m1));
    if (q1 >= (q_t)1) {
      check("acosh", mf::acosh(f1), acoshq(q1), mpfr::acosh(m1));
    }
    if (aq1 < (q_t)1) {
      check("atanh", mf::atanh(f1), atanhq(q1), mpfr::atanh(m1));
    }

    // special
    if (aq1 < (q_t)100) {
      check("erf",  mf::erf(f1),  erfq(q1),  mpfr::erf(m1));
      check("erfc", mf::erfc(f1), erfcq(q1), mpfr::erfc(m1));
    }
    if (q1 > (q_t)0 && q1 < (q_t)100) {
      check("tgamma", mf::tgamma(f1), tgammaq(q1), mpfr::gamma(m1));
      check("lgamma", mf::lgamma(f1), lgammaq(q1), mpfr::lngamma(m1));
    }

    // pow
    q_t aq2 = q2 < 0 ? -q2 : q2;
    if (q1 > (q_t)1e-3q && q1 < (q_t)1e3q && aq2 < (q_t)30) {
      check("pow", mf::pow(f1, f2), powq(q1, q2), mpfr::pow(m1, m2));
    }
  }

  print_report();
  return 0;
}
