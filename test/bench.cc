// Timing-only benchmark for every kernel category exposed by
// src/multifloats.hh, comparing MultiFloat<double, 2> against
// __float128 / libquadmath as a reference for the same operation.
// Precision measurement lives in test/fuzz.cc — this file is timing-only.
//
// Anti–dead-code-elimination strategy:
//   - Each rep computes the entire output array.
//   - At the end of each rep, a NOINLINE drain reads every element
//     into a plain `double` accumulator and feeds the result back into
//     one element of an input array. The NOINLINE barrier prevents the
//     optimizer from fusing the drain with the inner loop and eliding
//     reps after the first; the cross-rep dependency prevents loop
//     hoisting; the dp drain is much cheaper than the operation under
//     test, so the inner-loop work dominates the timing.
//   - Both legs use the same drain shape, so the qp/mf ratio is fair
//     even for one-cycle operations.

#include "multifloats.hh"

#include <quadmath.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>

namespace mf = multifloats;
using q_t = __float128;
using clk = std::chrono::steady_clock;

// =============================================================================
// Workspaces (file-static so the optimizer can't trivially track liveness
// through call boundaries) and rep counts.
// =============================================================================

static constexpr int N = 1024;
static constexpr int REPS_FAST = 400;       // +, -, *, /, abs, ...
static constexpr int REPS_TRIG = 40;        // exp, log, sin, ...
static constexpr int REPS_VERY_SLOW = 4;    // tgamma, erfc, ...

static q_t q1[N], q2[N], qpos[N], qsmall[N], qbnd[N], qres[N];
static mf::float64x2 f1[N], f2[N], fpos[N], fsmall[N], fbnd[N], fres[N];

static double q_sink = 0.0;
static double f_sink = 0.0;

// =============================================================================
// Plumbing
// =============================================================================

static mf::float64x2 to_mf2(q_t v) {
  double hi = (double)v;
  if (!std::isfinite(hi)) {
    mf::float64x2 r;
    r._limbs[0] = hi;
    r._limbs[1] = hi;
    return r;
  }
  double lo = (double)(v - (q_t)hi);
  double s = hi + lo;
  double err = lo - (s - hi);
  mf::float64x2 r;
  r._limbs[0] = s;
  r._limbs[1] = err;
  return r;
}

static void init_data() {
  std::mt19937_64 rng(42ULL);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  std::uniform_real_distribution<double> udd(-1.0, 1.0);

  // Build a genuine double-double from a leading double `hi` plus a random
  // perturbation at ~1 ULP of hi.  This gives a non-zero lo limb so that
  // precision measurements exercise the full DD path, matching the Fortran
  // fuzz test which also generates inputs with non-trivial lo limbs.
  auto make_dd = [&](double hi) -> std::pair<q_t, mf::float64x2> {
    double lo = hi * std::ldexp(udd(rng), -52);
    q_t v = (q_t)hi + (q_t)lo;
    return {v, to_mf2(v)};
  };

  // General arithmetic inputs in (-10, 10), nonzero (avoid div-by-0).
  auto away = [](double x) {
    return (x - 0.5) * 8.0 + (x >= 0.5 ? 0.25 : -0.25);
  };

  for (int i = 0; i < N; ++i) {
    double r1 = u(rng), r2 = u(rng), r3 = u(rng);
    auto [qv1, fv1] = make_dd(away(r1));    q1[i] = qv1;     f1[i] = fv1;
    auto [qv2, fv2] = make_dd(away(r2));    q2[i] = qv2;     f2[i] = fv2;
    // Positive in (0.1, 10) for log / sqrt / pow base / tgamma.
    auto [qvp, fvp] = make_dd(0.1 + 9.9 * r1);  qpos[i] = qvp; fpos[i] = fvp;
    // Small in (-3, 3) for trig / hyperbolic.
    auto [qvs, fvs] = make_dd((r2 - 0.5) * 6.0); qsmall[i] = qvs; fsmall[i] = fvs;
    // Bounded in (-0.9, 0.9) for asin / acos / atanh / pow exponent.
    auto [qvb, fvb] = make_dd((r3 - 0.5) * 1.8); qbnd[i] = qvb; fbnd[i] = fvb;
  }
}

// NOINLINE drain barriers. Reading every output element via a plain dp
// accumulator forces the inner loop to materialize qres / fres each rep.
// The feedback into one input element prevents loop hoisting.
__attribute__((noinline))
static void qfeed(q_t *in) {
  double s = 0.0;
  for (int i = 0; i < N; ++i) {
    s += (double)qres[i];
  }
  in[0] += (q_t)(s * 1e-30);
  q_sink += s;
}

__attribute__((noinline))
static void ffeed(mf::float64x2 *in) {
  double s = 0.0;
  for (int i = 0; i < N; ++i) {
    s += fres[i]._limbs[0];
  }
  in[0]._limbs[0] += s * 1e-30;
  f_sink += s;
}

static double seconds_since(clk::time_point t0) {
  return std::chrono::duration<double>(clk::now() - t0).count();
}

static void report(const char *name, long n_ops, double tq, double tf) {
  double speed = (tf > 0.0) ? (tq / tf) : 0.0;
  std::printf(" %-22s %12ld %10.4f %10.4f %9.2fx\n", name, n_ops, tq, tf,
              speed);
}

// =============================================================================
// Benchmark macro. Each invocation times one op for both legs.
//   NAME    : printed label
//   REPS    : repetition count
//   QEXPR   : per-element expression for the qp leg, may use `i`
//   MEXPR   : per-element expression for the mf leg, may use `i`
//   QFB     : input array (q_t*)  to feed back into for the qp leg
//   MFB     : input array (mf::float64x2*)  to feed back into for the mf leg
// =============================================================================

#define BENCH(NAME, REPS, QEXPR, MEXPR, QFB, MFB)                              \
  do {                                                                         \
    long n_ops = (long)N * (long)(REPS);                                       \
    /* Reinit before each leg so qp and mf start from identical clean          \
       inputs — drain feedback only affects within-leg reps. */                \
    init_data();                                                               \
    auto _t0 = clk::now();                                                     \
    for (int r = 0; r < (REPS); ++r) {                                         \
      for (int i = 0; i < N; ++i) {                                            \
        qres[i] = (QEXPR);                                                     \
      }                                                                        \
      qfeed(QFB);                                                              \
    }                                                                          \
    double tq = seconds_since(_t0);                                            \
    init_data();                                                               \
    _t0 = clk::now();                                                          \
    for (int r = 0; r < (REPS); ++r) {                                         \
      for (int i = 0; i < N; ++i) {                                            \
        fres[i] = (MEXPR);                                                     \
      }                                                                        \
      ffeed(MFB);                                                              \
    }                                                                          \
    double tf = seconds_since(_t0);                                            \
    report(NAME, n_ops, tq, tf);                                               \
  } while (0)

// =============================================================================
// Op categories
// =============================================================================

static void bench_arith() {
  BENCH("add",      REPS_FAST, q1[i] + q2[i],         f1[i] + f2[i],         q1, f1);
  BENCH("sub",      REPS_FAST, q1[i] - q2[i],         f1[i] - f2[i],         q1, f1);
  BENCH("mul",      REPS_FAST, q1[i] * q2[i],         f1[i] * f2[i],         q1, f1);
  BENCH("div",      REPS_FAST, q1[i] / q2[i],         f1[i] / f2[i],         q1, f1);
  BENCH("sqrt",     REPS_FAST, sqrtq(qpos[i]),        mf::sqrt(fpos[i]),     qpos, fpos);
  BENCH("cbrt",     REPS_FAST, cbrtq(q1[i]),          mf::cbrt(f1[i]),       q1, f1);
  BENCH("fma",      REPS_FAST, fmaq(q1[i], q2[i], qpos[i]),
                               mf::fma(f1[i], f2[i], fpos[i]),               q1, f1);
}

static void bench_unary_basic() {
  BENCH("abs",      REPS_FAST, fabsq(q1[i]),          mf::abs(f1[i]),        q1, f1);
  BENCH("neg",      REPS_FAST, -q1[i],                -f1[i],                q1, f1);
  BENCH("floor",    REPS_FAST, floorq(q1[i]),         mf::floor(f1[i]),      q1, f1);
  BENCH("ceil",     REPS_FAST, ceilq(q1[i]),          mf::ceil(f1[i]),       q1, f1);
  BENCH("trunc",    REPS_FAST, truncq(q1[i]),         mf::trunc(f1[i]),      q1, f1);
  BENCH("round",    REPS_FAST, roundq(q1[i]),         mf::round(f1[i]),      q1, f1);
  BENCH("rint",     REPS_FAST, rintq(q1[i]),          mf::rint(f1[i]),       q1, f1);
  BENCH("nearbyint", REPS_FAST, nearbyintq(q1[i]),    mf::nearbyint(f1[i]),  q1, f1);
}

static void bench_minmax_binary() {
  BENCH("fmin",     REPS_FAST, fminq(q1[i], q2[i]),   mf::fmin(f1[i], f2[i]), q1, f1);
  BENCH("fmax",     REPS_FAST, fmaxq(q1[i], q2[i]),   mf::fmax(f1[i], f2[i]), q1, f1);
  BENCH("fdim",     REPS_FAST, fdimq(q1[i], q2[i]),   mf::fdim(f1[i], f2[i]), q1, f1);
  BENCH("copysign", REPS_FAST, copysignq(q1[i], q2[i]),
                               mf::copysign(f1[i], f2[i]),                  q1, f1);
  BENCH("fmod",     REPS_FAST, fmodq(q1[i], q2[i]),   mf::fmod(f1[i], f2[i]), q1, f1);
  BENCH("hypot",    REPS_FAST, hypotq(q1[i], q2[i]),  mf::hypot(f1[i], f2[i]), q1, f1);
  BENCH("ldexp(.,5)", REPS_FAST, ldexpq(q1[i], 5),    mf::ldexp(f1[i], 5),   q1, f1);
}

static void bench_exp_log() {
  BENCH("exp",      REPS_TRIG, expq(qsmall[i]),       mf::exp(fsmall[i]),    qsmall, fsmall);
  BENCH("exp2",     REPS_TRIG, exp2q(qsmall[i]),      mf::exp2(fsmall[i]),   qsmall, fsmall);
  BENCH("expm1",    REPS_TRIG, expm1q(qsmall[i]),     mf::expm1(fsmall[i]),  qsmall, fsmall);
  BENCH("log",      REPS_TRIG, logq(qpos[i]),         mf::log(fpos[i]),      qpos, fpos);
  BENCH("log10",    REPS_TRIG, log10q(qpos[i]),       mf::log10(fpos[i]),    qpos, fpos);
  BENCH("log2",     REPS_TRIG, log2q(qpos[i]),        mf::log2(fpos[i]),     qpos, fpos);
  BENCH("log1p",    REPS_TRIG, log1pq(qpos[i]),       mf::log1p(fpos[i]),    qpos, fpos);
  BENCH("pow",      REPS_TRIG, powq(qpos[i], qbnd[i]),
                               mf::pow(fpos[i], fbnd[i]),                   qpos, fpos);
}

static void bench_trig() {
  BENCH("sin",      REPS_TRIG, sinq(qsmall[i]),       mf::sin(fsmall[i]),    qsmall, fsmall);
  BENCH("cos",      REPS_TRIG, cosq(qsmall[i]),       mf::cos(fsmall[i]),    qsmall, fsmall);
  BENCH("tan",      REPS_TRIG, tanq(qbnd[i]),         mf::tan(fbnd[i]),      qbnd, fbnd);
}

static void bench_inv_trig() {
  BENCH("asin",     REPS_TRIG, asinq(qbnd[i]),        mf::asin(fbnd[i]),     qbnd, fbnd);
  BENCH("acos",     REPS_TRIG, acosq(qbnd[i]),        mf::acos(fbnd[i]),     qbnd, fbnd);
  BENCH("atan",     REPS_TRIG, atanq(q1[i]),          mf::atan(f1[i]),       q1, f1);
  BENCH("atan2",    REPS_TRIG, atan2q(q1[i], q2[i]),  mf::atan2(f1[i], f2[i]), q1, f1);
}

static void bench_hyperbolic() {
  BENCH("sinh",     REPS_TRIG, sinhq(qsmall[i]),      mf::sinh(fsmall[i]),   qsmall, fsmall);
  BENCH("cosh",     REPS_TRIG, coshq(qsmall[i]),      mf::cosh(fsmall[i]),   qsmall, fsmall);
  BENCH("tanh",     REPS_TRIG, tanhq(qsmall[i]),      mf::tanh(fsmall[i]),   qsmall, fsmall);
}

static void bench_inv_hyperbolic() {
  BENCH("asinh",    REPS_TRIG, asinhq(q1[i]),         mf::asinh(f1[i]),      q1, f1);
  BENCH("acosh",    REPS_TRIG, acoshq(1.0q + qpos[i]),
                               mf::acosh(mf::float64x2(1.0) + fpos[i]),               qpos, fpos);
  BENCH("atanh",    REPS_TRIG, atanhq(qbnd[i]),       mf::atanh(fbnd[i]),    qbnd, fbnd);
}

static void bench_erf_gamma() {
  BENCH("erf",       REPS_VERY_SLOW, erfq(q1[i]),     mf::erf(f1[i]),        q1, f1);
  BENCH("erfc",      REPS_VERY_SLOW, erfcq(q1[i]),    mf::erfc(f1[i]),       q1, f1);
  BENCH("tgamma",    REPS_VERY_SLOW, tgammaq(qpos[i]),
                                     mf::tgamma(fpos[i]),                   qpos, fpos);
  BENCH("lgamma",    REPS_VERY_SLOW, lgammaq(qpos[i]),
                                     mf::lgamma(fpos[i]),                   qpos, fpos);
}

// =============================================================================
// main
// =============================================================================

int main() {
  init_data();

  std::printf("================================================================\n");
  std::printf(" multifloats C++ benchmark — N=%d elements per rep, batched\n", N);
  std::printf("================================================================\n\n");
  std::printf(" %-22s %12s %10s %10s %10s\n", "op", "n_ops", "qp [s]",
              "mf [s]", "speedup");
  std::printf(" ----------------------------------------------------------------\n");

  bench_arith();
  bench_unary_basic();
  bench_minmax_binary();
  bench_exp_log();
  bench_trig();
  bench_inv_trig();
  bench_hyperbolic();
  bench_inv_hyperbolic();
  bench_erf_gamma();

  std::printf(" ----------------------------------------------------------------\n");
  std::printf(" sinks (must remain live): qp_sink=%g  dd_sink=%g\n", q_sink,
              f_sink);
  return 0;
}
