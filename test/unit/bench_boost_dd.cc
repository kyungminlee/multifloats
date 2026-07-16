// Timing benchmark for boost::multiprecision::cpp_double_double, mirroring
// test/bench.cc structure: same workspace size, same rep counts, same anti-DCE
// drain pattern, same column layout. Both legs (qp via libquadmath, DD via
// boost) get identical inputs round-tripped through the boost backend so qp
// holds nothing past what DD can represent.
//
// Op coverage matches what boost actually provides natively or through
// Boost.Math; Bessel / complex / π-trig / array kernels are out of scope.

#include "test_common.hh"

#include <boost/multiprecision/cpp_double_fp.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/cbrt.hpp>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>

#include <quadmath.h>

using boost::multiprecision::cpp_double_double;
using multifloats_test::q_t;
using multifloats_test::from_q;
using clk = std::chrono::steady_clock;

// =============================================================================
// boost <-> qp adapters (same as test_boost_dd.cc).
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

// =============================================================================
// fx:: overload set — picks the boost cpp_double_double or qp leg from one
// op name. Unary / binary mirrors test/bench.cc's fx::WRAP_REAL_UNARY pattern.
// =============================================================================

namespace fx {

#define WRAP_UNARY_STD(NAME, QFN)                                              \
  inline cpp_double_double NAME(cpp_double_double const &x) {                  \
    using std::NAME; return NAME(x);                                           \
  }                                                                            \
  inline q_t NAME(q_t x) { return ::QFN(x); }

WRAP_UNARY_STD(sqrt,  sqrtq)
WRAP_UNARY_STD(trunc, truncq)
WRAP_UNARY_STD(round, roundq)
WRAP_UNARY_STD(floor, floorq)
WRAP_UNARY_STD(ceil,  ceilq)
WRAP_UNARY_STD(exp,   expq)
WRAP_UNARY_STD(exp2,  exp2q)
WRAP_UNARY_STD(expm1, expm1q)
WRAP_UNARY_STD(log,   logq)
WRAP_UNARY_STD(log2,  log2q)
WRAP_UNARY_STD(log10, log10q)
WRAP_UNARY_STD(log1p, log1pq)
WRAP_UNARY_STD(sin,   sinq)
WRAP_UNARY_STD(cos,   cosq)
WRAP_UNARY_STD(tan,   tanq)
WRAP_UNARY_STD(asin,  asinq)
WRAP_UNARY_STD(acos,  acosq)
WRAP_UNARY_STD(atan,  atanq)
WRAP_UNARY_STD(sinh,  sinhq)
WRAP_UNARY_STD(cosh,  coshq)
WRAP_UNARY_STD(tanh,  tanhq)
WRAP_UNARY_STD(asinh, asinhq)
WRAP_UNARY_STD(acosh, acoshq)
WRAP_UNARY_STD(atanh, atanhq)
#undef WRAP_UNARY_STD

inline cpp_double_double abs_op(cpp_double_double const &x) { using std::fabs; return fabs(x); }
inline q_t abs_op(q_t x) { return ::fabsq(x); }
inline cpp_double_double neg_op(cpp_double_double const &x) { return -x; }
inline q_t neg_op(q_t x) { return -x; }

#define WRAP_BINARY_STD(NAME, QFN)                                             \
  inline cpp_double_double NAME(cpp_double_double const &a,                    \
                                cpp_double_double const &b) {                  \
    using std::NAME; return NAME(a, b);                                        \
  }                                                                            \
  inline q_t NAME(q_t a, q_t b) { return ::QFN(a, b); }

WRAP_BINARY_STD(fmin,     fminq)
WRAP_BINARY_STD(fmax,     fmaxq)
WRAP_BINARY_STD(copysign, copysignq)
WRAP_BINARY_STD(fmod,     fmodq)
WRAP_BINARY_STD(hypot,    hypotq)
WRAP_BINARY_STD(pow,      powq)
WRAP_BINARY_STD(atan2,    atan2q)
#undef WRAP_BINARY_STD

inline cpp_double_double fdim_op(cpp_double_double const &a, cpp_double_double const &b) {
  return a > b ? cpp_double_double(a - b) : cpp_double_double(0);
}
inline q_t fdim_op(q_t a, q_t b) { return ::fdimq(a, b); }

inline cpp_double_double fma_op(cpp_double_double const &a,
                                cpp_double_double const &b,
                                cpp_double_double const &c) {
  using std::fma; return fma(a, b, c);
}
inline q_t fma_op(q_t a, q_t b, q_t c) { return ::fmaq(a, b, c); }

inline cpp_double_double cbrt_op(cpp_double_double const &x) { return boost::math::cbrt(x); }
inline q_t                cbrt_op(q_t x) { return ::cbrtq(x); }

// Bessel via Boost.Math (integer order).
inline cpp_double_double bj0(cpp_double_double const &x) { return boost::math::cyl_bessel_j(0, x); }
inline q_t                bj0(q_t x) { return ::j0q(x); }
inline cpp_double_double bj1(cpp_double_double const &x) { return boost::math::cyl_bessel_j(1, x); }
inline q_t                bj1(q_t x) { return ::j1q(x); }
inline cpp_double_double by0(cpp_double_double const &x) { return boost::math::cyl_neumann(0, x); }
inline q_t                by0(q_t x) { return ::y0q(x); }
inline cpp_double_double by1(cpp_double_double const &x) { return boost::math::cyl_neumann(1, x); }
inline q_t                by1(q_t x) { return ::y1q(x); }

// erf / erfc / tgamma / lgamma via Boost.Math.
inline cpp_double_double erf_op  (cpp_double_double const &x) { return boost::math::erf(x); }
inline q_t                erf_op  (q_t x) { return ::erfq(x); }
inline cpp_double_double erfc_op (cpp_double_double const &x) { return boost::math::erfc(x); }
inline q_t                erfc_op (q_t x) { return ::erfcq(x); }
inline cpp_double_double tgamma_op(cpp_double_double const &x) { return boost::math::tgamma(x); }
inline q_t                tgamma_op(q_t x) { return ::tgammaq(x); }
inline cpp_double_double lgamma_op(cpp_double_double const &x) { return boost::math::lgamma(x); }
inline q_t                lgamma_op(q_t x) { return ::lgammaq(x); }

}  // namespace fx

// =============================================================================
// Workspaces — same shape and population strategy as test/bench.cc.
// =============================================================================

static constexpr int N              = 1024;
static constexpr int REPS_FAST      = 400;
static constexpr int REPS_TRIG      = 40;
static constexpr int REPS_VERY_SLOW = 4;
static constexpr int    DD_LO_ULP_EXPONENT   = -52;
static constexpr double DRAIN_FEEDBACK_SCALE = 1e-30;

static q_t                q1[N], q2[N], qpos[N], qsmall[N], qbnd[N], qres[N];
static cpp_double_double  f1[N], f2[N], fpos[N], fsmall[N], fbnd[N], fres[N];

static double q_sink = 0.0, f_sink = 0.0;

static void init_data() {
  std::mt19937_64 rng(42ULL);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  std::uniform_real_distribution<double> udd(-1.0, 1.0);

  auto make_dd = [&](double hi) -> std::pair<q_t, cpp_double_double> {
    double lo = hi * std::ldexp(udd(rng), DD_LO_ULP_EXPONENT);
    q_t v = (q_t)hi + (q_t)lo;
    return {v, bdd_from_q(v)};
  };
  auto away = [](double x) {
    return (x - 0.5) * 8.0 + (x >= 0.5 ? 0.25 : -0.25);
  };

  for (int i = 0; i < N; ++i) {
    double r1 = u(rng), r2 = u(rng), r3 = u(rng);
    auto [qv1, fv1] = make_dd(away(r1));    q1[i] = qv1;    f1[i] = fv1;
    auto [qv2, fv2] = make_dd(away(r2));    q2[i] = qv2;    f2[i] = fv2;
    auto [qvp, fvp] = make_dd(0.1 + 9.9 * r1); qpos[i] = qvp; fpos[i] = fvp;
    auto [qvs, fvs] = make_dd((r2 - 0.5) * 6.0); qsmall[i] = qvs; fsmall[i] = fvs;
    auto [qvb, fvb] = make_dd((r3 - 0.5) * 1.8); qbnd[i] = qvb; fbnd[i] = fvb;
  }
}

__attribute__((noinline))
static void qfeed(q_t *in) {
  double s = 0.0;
  for (int i = 0; i < N; ++i) s += (double)qres[i];
  in[0] += (q_t)(s * DRAIN_FEEDBACK_SCALE);
  q_sink += s;
}
__attribute__((noinline))
static void ffeed(cpp_double_double *in) {
  double s = 0.0;
  for (int i = 0; i < N; ++i) {
    auto const &r = fres[i].backend().crep();
    s += r.first + r.second;
  }
  auto &r0 = in[0].backend().rep();
  r0.first += s * DRAIN_FEEDBACK_SCALE;
  f_sink += s;
}

static double seconds_since(clk::time_point t0) {
  return std::chrono::duration<double>(clk::now() - t0).count();
}
static void report(const char *name, long n_ops, double tq, double tf) {
  double speed = (tf > 0.0) ? (tq / tf) : 0.0;
  std::printf(" %-22s %12ld %10.4f %10.4f %9.2fx\n", name, n_ops, tq, tf, speed);
}

// =============================================================================
// BENCH macro — same shape as test/bench.cc.
// =============================================================================

#define BENCH(NAME, REPS, QEXPR, MEXPR, QFB, MFB)                              \
  do {                                                                         \
    long n_ops = (long)N * (long)(REPS);                                       \
    init_data();                                                               \
    auto _t0 = clk::now();                                                     \
    for (int r = 0; r < (REPS); ++r) {                                         \
      for (int i = 0; i < N; ++i) qres[i] = (QEXPR);                           \
      qfeed(QFB);                                                              \
    }                                                                          \
    double tq = seconds_since(_t0);                                            \
    init_data();                                                               \
    _t0 = clk::now();                                                          \
    for (int r = 0; r < (REPS); ++r) {                                         \
      for (int i = 0; i < N; ++i) fres[i] = (MEXPR);                           \
      ffeed(MFB);                                                              \
    }                                                                          \
    double tf = seconds_since(_t0);                                            \
    report(NAME, n_ops, tq, tf);                                               \
  } while (0)

#define BENCH1(NAME, REPS, QIN, FIN) \
  BENCH(#NAME, REPS, fx::NAME(QIN[i]), fx::NAME(FIN[i]), QIN, FIN)
#define BENCH2(NAME, REPS, QA, FA, QB, FB) \
  BENCH(#NAME, REPS, fx::NAME(QA[i], QB[i]), fx::NAME(FA[i], FB[i]), QA, FA)

// =============================================================================
// Op categories
// =============================================================================

static void bench_arith() {
  BENCH("add",  REPS_FAST, q1[i] + q2[i], f1[i] + f2[i], q1, f1);
  BENCH("sub",  REPS_FAST, q1[i] - q2[i], f1[i] - f2[i], q1, f1);
  BENCH("mul",  REPS_FAST, q1[i] * q2[i], f1[i] * f2[i], q1, f1);
  BENCH("div",  REPS_FAST, q1[i] / q2[i], f1[i] / f2[i], q1, f1);
  BENCH1(sqrt,  REPS_FAST, qpos, fpos);
  BENCH("fma",  REPS_FAST,
        fx::fma_op(q1[i], q2[i], qpos[i]),
        fx::fma_op(f1[i], f2[i], fpos[i]), q1, f1);
  BENCH("abs",  REPS_FAST, fx::abs_op(q1[i]), fx::abs_op(f1[i]), q1, f1);
  BENCH("neg",  REPS_FAST, fx::neg_op(q1[i]), fx::neg_op(f1[i]), q1, f1);
}

static void bench_rounding() {
  BENCH1(trunc, REPS_FAST, q1, f1);
  BENCH1(round, REPS_FAST, q1, f1);
  BENCH1(floor, REPS_FAST, q1, f1);
  BENCH1(ceil,  REPS_FAST, q1, f1);
}

static void bench_misc() {
  BENCH("cbrt", REPS_FAST, fx::cbrt_op(qpos[i]), fx::cbrt_op(fpos[i]), qpos, fpos);
  BENCH("ldexp", REPS_FAST,
        ::ldexpq(q1[i], 5),
        boost::multiprecision::ldexp(f1[i], 5), q1, f1);
  BENCH("scalbn", REPS_FAST,
        ::scalbnq(q1[i], 5),
        boost::multiprecision::scalbn(f1[i], 5), q1, f1);
}

static void bench_binary() {
  BENCH2(fmin,     REPS_FAST, q1, f1, q2, f2);
  BENCH2(fmax,     REPS_FAST, q1, f1, q2, f2);
  BENCH("fdim",    REPS_FAST,
        fx::fdim_op(q1[i], q2[i]), fx::fdim_op(f1[i], f2[i]), q1, f1);
  BENCH2(copysign, REPS_FAST, q1, f1, q2, f2);
  BENCH2(fmod,     REPS_FAST, q1, f1, q2, f2);
  BENCH2(hypot,    REPS_FAST, q1, f1, q2, f2);
}

static void bench_exp_log() {
  BENCH1(exp,    REPS_TRIG, qsmall, fsmall);
  BENCH1(exp2,   REPS_TRIG, qsmall, fsmall);
  BENCH1(expm1,  REPS_TRIG, qsmall, fsmall);
  BENCH1(log,    REPS_TRIG, qpos, fpos);
  BENCH1(log10,  REPS_TRIG, qpos, fpos);
  BENCH1(log2,   REPS_TRIG, qpos, fpos);
  BENCH1(log1p,  REPS_TRIG, qpos, fpos);
  BENCH2(pow,    REPS_TRIG, qpos, fpos, qbnd, fbnd);
}

static void bench_trig() {
  BENCH1(sin, REPS_TRIG, qsmall, fsmall);
  BENCH1(cos, REPS_TRIG, qsmall, fsmall);
  BENCH1(tan, REPS_TRIG, qbnd,   fbnd);
}

static void bench_inv_trig() {
  BENCH1(asin,  REPS_TRIG, qbnd, fbnd);
  BENCH1(acos,  REPS_TRIG, qbnd, fbnd);
  BENCH1(atan,  REPS_TRIG, q1,   f1);
  BENCH2(atan2, REPS_TRIG, q1, f1, q2, f2);
}

static void bench_hyperbolic() {
  BENCH1(sinh, REPS_TRIG, qsmall, fsmall);
  BENCH1(cosh, REPS_TRIG, qsmall, fsmall);
  BENCH1(tanh, REPS_TRIG, qsmall, fsmall);
}

static void bench_inv_hyperbolic() {
  BENCH1(asinh, REPS_TRIG, q1, f1);
  // acosh expects argument >= 1; shift fpos by +1.
  BENCH("acosh", REPS_TRIG,
        fx::acosh((q_t)1.0q + qpos[i]),
        fx::acosh(cpp_double_double(1) + fpos[i]), qpos, fpos);
  BENCH1(atanh, REPS_TRIG, qbnd, fbnd);
}

static void bench_bessel() {
  BENCH("bj0", REPS_VERY_SLOW, fx::bj0(q1[i]), fx::bj0(f1[i]), q1, f1);
  BENCH("bj1", REPS_VERY_SLOW, fx::bj1(q1[i]), fx::bj1(f1[i]), q1, f1);
  BENCH("bjn", REPS_VERY_SLOW,
        ::jnq(3, q1[i]), boost::math::cyl_bessel_j(3, f1[i]), q1, f1);
  BENCH("by0", REPS_VERY_SLOW, fx::by0(qpos[i]), fx::by0(fpos[i]), qpos, fpos);
  BENCH("by1", REPS_VERY_SLOW, fx::by1(qpos[i]), fx::by1(fpos[i]), qpos, fpos);
  BENCH("byn", REPS_VERY_SLOW,
        ::ynq(3, qpos[i]), boost::math::cyl_neumann(3, fpos[i]), qpos, fpos);
}

static void bench_erf_gamma() {
  BENCH("erf",    REPS_VERY_SLOW, fx::erf_op(q1[i]),    fx::erf_op(f1[i]),    q1, f1);
  BENCH("erfc",   REPS_VERY_SLOW, fx::erfc_op(q1[i]),   fx::erfc_op(f1[i]),   q1, f1);
  BENCH("tgamma", REPS_VERY_SLOW, fx::tgamma_op(qpos[i]), fx::tgamma_op(fpos[i]), qpos, fpos);
  BENCH("lgamma", REPS_VERY_SLOW, fx::lgamma_op(qpos[i]), fx::lgamma_op(fpos[i]), qpos, fpos);
}

int main() {
  init_data();

  std::printf("================================================================\n");
  std::printf(" boost::cpp_double_double benchmark — N=%d elements per rep\n", N);
  std::printf(" (DD leg: boost::multiprecision::cpp_double_double)\n");
  std::printf("================================================================\n\n");
  std::printf(" %-22s %12s %10s %10s %10s\n", "op", "n_ops", "qp [s]",
              "boost [s]", "speedup");
  std::printf(" ----------------------------------------------------------------\n");

  bench_arith();
  bench_rounding();
  bench_misc();
  bench_binary();
  bench_exp_log();
  bench_trig();
  bench_inv_trig();
  bench_hyperbolic();
  bench_inv_hyperbolic();
  bench_erf_gamma();
  bench_bessel();

  std::printf(" ----------------------------------------------------------------\n");
  std::printf(" sinks (must remain live): qp=%g  boost=%g\n", q_sink, f_sink);
  return 0;
}
