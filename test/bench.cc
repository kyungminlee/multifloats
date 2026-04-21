// Timing benchmark through the multifloats C ABI.
//
// This file times every kernel exposed in src/multifloats_c.h (sindd,
// cdd_muldd, j0dd, matmuldd_mm, …) against libquadmath as a reference
// for the SAME operation. The DD leg goes through the register-return
// C ABI (not the C++ header-only inline path), so the timings reflect
// what a real C or Fortran-bind-c caller sees.
//
// Precision measurement lives in test/fuzz_mpfr.cc; this file is
// timing-only.
//
// The three function families follow a consistent suffix scheme —
//   DD kernel      : NAMEdd     (C ABI)
//   libquadmath    : NAMEq      (113-bit __float128)
// — which we exploit by declaring overloads in `namespace fx` so a single
// BENCH1(NAME, ...) macro can dispatch both legs uniformly.

#include "multifloats.hh"
#include "multifloats_c.h"
#include "test_common.hh"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>

namespace mf = multifloats;
using multifloats_test::q_t;
using multifloats_test::from_q;
using clk = std::chrono::steady_clock;

// =============================================================================
// namespace fx — overloads of each op for (DD via C ABI) and (qp via libquadmath).
// =============================================================================

namespace fx {

// -- Real unary: NAMEdd and NAMEq -------------------------------------------
#define WRAP_REAL_UNARY(NAME)                                                  \
  inline mf::float64x2 NAME(mf::float64x2 const &x) {                          \
    return mf::detail::from_f64x2(::NAME##dd(mf::detail::to_f64x2(x)));        \
  }                                                                            \
  inline q_t NAME(q_t x) { return ::NAME##q(x); }

// neg: libquadmath uses the unary `-` operator on q_t, not a `negq` fn.
inline mf::float64x2 neg_op(mf::float64x2 const &x) {
  return mf::detail::from_f64x2(::negdd(mf::detail::to_f64x2(x)));
}
inline q_t neg_op(q_t x) { return -x; }

WRAP_REAL_UNARY(sqrt)
WRAP_REAL_UNARY(trunc)
WRAP_REAL_UNARY(round)
WRAP_REAL_UNARY(exp)
WRAP_REAL_UNARY(exp2)
WRAP_REAL_UNARY(expm1)
WRAP_REAL_UNARY(log)
WRAP_REAL_UNARY(log2)
WRAP_REAL_UNARY(log10)
WRAP_REAL_UNARY(log1p)
WRAP_REAL_UNARY(sin)
WRAP_REAL_UNARY(cos)
WRAP_REAL_UNARY(tan)
WRAP_REAL_UNARY(asin)
WRAP_REAL_UNARY(acos)
WRAP_REAL_UNARY(atan)
WRAP_REAL_UNARY(sinh)
WRAP_REAL_UNARY(cosh)
WRAP_REAL_UNARY(tanh)
WRAP_REAL_UNARY(asinh)
WRAP_REAL_UNARY(acosh)
WRAP_REAL_UNARY(atanh)
WRAP_REAL_UNARY(erf)
WRAP_REAL_UNARY(erfc)
// no ::erfcxq in libquadmath; omit erfcx from the bench.
WRAP_REAL_UNARY(tgamma)
WRAP_REAL_UNARY(lgamma)
// fabs: libquadmath is fabsq; DD is fabsdd. "abs" common name for readability.
inline mf::float64x2 abs_op(mf::float64x2 const &x) {
  return mf::detail::from_f64x2(::fabsdd(mf::detail::to_f64x2(x)));
}
inline q_t abs_op(q_t x) { return ::fabsq(x); }
#undef WRAP_REAL_UNARY

// -- π-scaled trig (DD: NAMEdd, qp: library fn composed with π·x or /π) ----
#define WRAP_PI_INPUT(NAME, FN)                                                \
  inline mf::float64x2 NAME(mf::float64x2 const &x) {                          \
    return mf::detail::from_f64x2(::NAME##dd(mf::detail::to_f64x2(x)));        \
  }                                                                            \
  inline q_t NAME(q_t x) { return ::FN##q(q_t(M_PIq) * x); }

WRAP_PI_INPUT(sinpi, sin)
WRAP_PI_INPUT(cospi, cos)
WRAP_PI_INPUT(tanpi, tan)
#undef WRAP_PI_INPUT

#define WRAP_PI_OUTPUT(NAME, FN)                                               \
  inline mf::float64x2 NAME(mf::float64x2 const &x) {                          \
    return mf::detail::from_f64x2(::NAME##dd(mf::detail::to_f64x2(x)));        \
  }                                                                            \
  inline q_t NAME(q_t x) { return ::FN##q(x) / q_t(M_PIq); }

WRAP_PI_OUTPUT(asinpi, asin)
WRAP_PI_OUTPUT(acospi, acos)
WRAP_PI_OUTPUT(atanpi, atan)
#undef WRAP_PI_OUTPUT

// -- Bessel (DD: NAMEdd, qp: NAMEq, both identical suffix scheme) ----------
#define WRAP_BESSEL01(NAME)                                                    \
  inline mf::float64x2 NAME(mf::float64x2 const &x) {                          \
    return mf::detail::from_f64x2(::NAME##dd(mf::detail::to_f64x2(x)));        \
  }                                                                            \
  inline q_t NAME(q_t x) { return ::NAME##q(x); }

WRAP_BESSEL01(j0)
WRAP_BESSEL01(j1)
WRAP_BESSEL01(y0)
WRAP_BESSEL01(y1)
#undef WRAP_BESSEL01

#define WRAP_BESSELN(NAME)                                                     \
  inline mf::float64x2 NAME(int n, mf::float64x2 const &x) {                   \
    return mf::detail::from_f64x2(::NAME##dd(n, mf::detail::to_f64x2(x)));     \
  }                                                                            \
  inline q_t NAME(int n, q_t x) { return ::NAME##q(n, x); }

WRAP_BESSELN(jn)
WRAP_BESSELN(yn)
#undef WRAP_BESSELN

// -- Real binary -----------------------------------------------------------
#define WRAP_REAL_BINARY(NAME)                                                 \
  inline mf::float64x2 NAME(mf::float64x2 const &a, mf::float64x2 const &b) {  \
    return mf::detail::from_f64x2(                                             \
        ::NAME##dd(mf::detail::to_f64x2(a), mf::detail::to_f64x2(b)));         \
  }                                                                            \
  inline q_t NAME(q_t a, q_t b) { return ::NAME##q(a, b); }

WRAP_REAL_BINARY(fmin)
WRAP_REAL_BINARY(fmax)
WRAP_REAL_BINARY(fdim)
WRAP_REAL_BINARY(copysign)
WRAP_REAL_BINARY(fmod)
WRAP_REAL_BINARY(hypot)
WRAP_REAL_BINARY(pow)
WRAP_REAL_BINARY(atan2)
#undef WRAP_REAL_BINARY

// fma: (a, b, c) → a·b + c
inline mf::float64x2 fma_op(mf::float64x2 const &a, mf::float64x2 const &b, mf::float64x2 const &c) {
  return mf::detail::from_f64x2(::fmadd(mf::detail::to_f64x2(a),
                                         mf::detail::to_f64x2(b),
                                         mf::detail::to_f64x2(c)));
}
inline q_t fma_op(q_t a, q_t b, q_t c) { return ::fmaq(a, b, c); }

// -- Complex unary: c<NAME>dd / c<NAME>q -----------------------------------
#define WRAP_COMPLEX_UNARY(NAME)                                               \
  inline complex64x2_t c##NAME(complex64x2_t z) { return ::c##NAME##dd(z); }   \
  inline __complex128  c##NAME(__complex128  z) { return ::c##NAME##q(z);  }

WRAP_COMPLEX_UNARY(sqrt)
WRAP_COMPLEX_UNARY(exp)
WRAP_COMPLEX_UNARY(log)
WRAP_COMPLEX_UNARY(sin)
WRAP_COMPLEX_UNARY(cos)
WRAP_COMPLEX_UNARY(tan)
WRAP_COMPLEX_UNARY(sinh)
WRAP_COMPLEX_UNARY(cosh)
WRAP_COMPLEX_UNARY(tanh)
WRAP_COMPLEX_UNARY(asin)
WRAP_COMPLEX_UNARY(acos)
WRAP_COMPLEX_UNARY(atan)
WRAP_COMPLEX_UNARY(asinh)
WRAP_COMPLEX_UNARY(acosh)
WRAP_COMPLEX_UNARY(atanh)
#undef WRAP_COMPLEX_UNARY

// cabs: complex → real
inline mf::float64x2 cabs(complex64x2_t z) { return mf::detail::from_f64x2(::cabsdd(z)); }
inline q_t           cabs(__complex128  z) { return ::cabsq(z); }

// conjg: complex → complex
inline complex64x2_t cconjg(complex64x2_t z) { return ::conjdd(z); }
inline __complex128  cconjg(__complex128  z) { return ::conjq(z); }

} // namespace fx

// =============================================================================
// Workspaces + drain (anti-DCE scaffolding). Kept similar to the original
// bench.cc so speedup numbers remain comparable across commits.
// =============================================================================

static constexpr int N                = 1024;
static constexpr int REPS_FAST        = 400;   // +, -, *, /, abs, ...
static constexpr int REPS_TRIG        = 40;    // exp, log, sin, ...
static constexpr int REPS_VERY_SLOW   = 4;     // tgamma, erfc, bessel, complex trans

// 2^-52: one double-ulp perturbation for DD lo-limb.
static constexpr int    DD_LO_ULP_EXPONENT = -52;
static constexpr double DRAIN_FEEDBACK_SCALE = 1e-30;

// Real workspaces.
static q_t           q1[N], q2[N], qpos[N], qsmall[N], qbnd[N], qres[N];
static mf::float64x2 f1[N], f2[N], fpos[N], fsmall[N], fbnd[N], fres[N];

// Complex workspaces.
static __complex128  zq1[N], zq2[N], zqres[N];
static complex64x2_t zf1[N], zf2[N], zfres[N];

static double q_sink = 0.0, f_sink = 0.0, zq_sink = 0.0, zf_sink = 0.0;

static void init_data() {
  std::mt19937_64 rng(42ULL);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  std::uniform_real_distribution<double> udd(-1.0, 1.0);

  // Genuine DD: hi + lo·2^-52 perturbation.
  auto make_dd = [&](double hi) -> std::pair<q_t, mf::float64x2> {
    double lo = hi * std::ldexp(udd(rng), DD_LO_ULP_EXPONENT);
    q_t v = (q_t)hi + (q_t)lo;
    return {v, from_q(v)};
  };
  auto away = [](double x) {
    return (x - 0.5) * 8.0 + (x >= 0.5 ? 0.25 : -0.25);
  };

  for (int i = 0; i < N; ++i) {
    double r1 = u(rng), r2 = u(rng), r3 = u(rng);
    auto [qv1, fv1] = make_dd(away(r1));    q1[i] = qv1;     f1[i] = fv1;
    auto [qv2, fv2] = make_dd(away(r2));    q2[i] = qv2;     f2[i] = fv2;
    auto [qvp, fvp] = make_dd(0.1 + 9.9 * r1);  qpos[i] = qvp; fpos[i] = fvp;
    auto [qvs, fvs] = make_dd((r2 - 0.5) * 6.0); qsmall[i] = qvs; fsmall[i] = fvs;
    auto [qvb, fvb] = make_dd((r3 - 0.5) * 1.8); qbnd[i] = qvb; fbnd[i] = fvb;
    // Complex: build from two fresh DD scalars per slot.
    auto [qre, fre] = make_dd(away(u(rng)));
    auto [qim, fim] = make_dd(away(u(rng)));
    __real__ zq1[i] = qre;  __imag__ zq1[i] = qim;
    zf1[i].re = mf::detail::to_f64x2(fre);
    zf1[i].im = mf::detail::to_f64x2(fim);
    auto [qre2, fre2] = make_dd(away(u(rng)));
    auto [qim2, fim2] = make_dd(away(u(rng)));
    __real__ zq2[i] = qre2; __imag__ zq2[i] = qim2;
    zf2[i].re = mf::detail::to_f64x2(fre2);
    zf2[i].im = mf::detail::to_f64x2(fim2);
  }
}

// Drain functions: sum one double per output into a sink; feed back into
// one input element with a tiny scale to keep the dep chain live.
__attribute__((noinline))
static void qfeed(q_t *in) {
  double s = 0.0;
  for (int i = 0; i < N; ++i) s += (double)qres[i];
  in[0] += (q_t)(s * DRAIN_FEEDBACK_SCALE);
  q_sink += s;
}
__attribute__((noinline))
static void ffeed(mf::float64x2 *in) {
  double s = 0.0;
  for (int i = 0; i < N; ++i) s += fres[i]._limbs[0];
  in[0]._limbs[0] += s * DRAIN_FEEDBACK_SCALE;
  f_sink += s;
}
__attribute__((noinline))
static void zqfeed(__complex128 *in) {
  double s = 0.0;
  for (int i = 0; i < N; ++i) s += (double)crealq(zqres[i]);
  __real__ in[0] += (q_t)(s * DRAIN_FEEDBACK_SCALE);
  zq_sink += s;
}
__attribute__((noinline))
static void zffeed(complex64x2_t *in) {
  double s = 0.0;
  for (int i = 0; i < N; ++i) s += zfres[i].re.hi;
  in[0].re.hi += s * DRAIN_FEEDBACK_SCALE;
  zf_sink += s;
}

static double seconds_since(clk::time_point t0) {
  return std::chrono::duration<double>(clk::now() - t0).count();
}
static void report(const char *name, long n_ops, double tq, double tf) {
  double speed = (tf > 0.0) ? (tq / tf) : 0.0;
  std::printf(" %-22s %12ld %10.4f %10.4f %9.2fx\n", name, n_ops, tq, tf, speed);
}

// =============================================================================
// BENCH macro — unchanged shape from the original bench.cc so downstream
// tooling (run_benchmarks.py parser) is unaffected. Per-element
// expressions reference `i` and the macro handles rep loops + drain
// + timing.
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

#define BENCHZ(NAME, REPS, QEXPR, MEXPR, QFB, MFB)                             \
  do {                                                                         \
    long n_ops = (long)N * (long)(REPS);                                       \
    init_data();                                                               \
    auto _t0 = clk::now();                                                     \
    for (int r = 0; r < (REPS); ++r) {                                         \
      for (int i = 0; i < N; ++i) zqres[i] = (QEXPR);                          \
      zqfeed(QFB);                                                             \
    }                                                                          \
    double tq = seconds_since(_t0);                                            \
    init_data();                                                               \
    _t0 = clk::now();                                                          \
    for (int r = 0; r < (REPS); ++r) {                                         \
      for (int i = 0; i < N; ++i) zfres[i] = (MEXPR);                          \
      zffeed(MFB);                                                             \
    }                                                                          \
    double tf = seconds_since(_t0);                                            \
    report(NAME, n_ops, tq, tf);                                               \
  } while (0)

// Unary convenience: BENCH1(sin, REPS_TRIG, qsmall, fsmall);
#define BENCH1(NAME, REPS, QIN, FIN) \
  BENCH(#NAME, REPS, fx::NAME(QIN[i]), fx::NAME(FIN[i]), QIN, FIN)
// Binary convenience: BENCH2(pow, REPS_TRIG, qpos, fpos, qbnd, fbnd);
#define BENCH2(NAME, REPS, QA, FA, QB, FB) \
  BENCH(#NAME, REPS, fx::NAME(QA[i], QB[i]), fx::NAME(FA[i], FB[i]), QA, FA)
// Bessel with fixed integer order n.
#define BENCH_BESN(LABEL, NAME, N_ORD, REPS, QIN, FIN) \
  BENCH(LABEL, REPS, fx::NAME(N_ORD, QIN[i]), fx::NAME(N_ORD, FIN[i]), QIN, FIN)
// Complex unary via C ABI: BENCHZ1(exp) → cdd_exp row.
#define BENCHZ1(STEM) \
  BENCHZ("cdd_" #STEM, REPS_VERY_SLOW, fx::c##STEM(zq1[i]), fx::c##STEM(zf1[i]), zq1, zf1)
#define BENCHZ1_FAST(STEM, REPS) \
  BENCHZ("cdd_" #STEM, REPS, fx::c##STEM(zq1[i]), fx::c##STEM(zf1[i]), zq1, zf1)
// Complex binary via operator (cadd, csub, cmul, cdiv).
#define BENCHZ2_OP(STEM, OP) \
  BENCHZ("cdd_" #STEM, REPS_FAST, zq1[i] OP zq2[i], c##STEM##dd(zf1[i], zf2[i]), zq1, zf1)

// =============================================================================
// Op categories
// =============================================================================

static void bench_arith() {
  BENCH("add",  REPS_FAST, q1[i] + q2[i], mf::detail::from_f64x2(adddd(mf::detail::to_f64x2(f1[i]), mf::detail::to_f64x2(f2[i]))), q1, f1);
  BENCH("sub",  REPS_FAST, q1[i] - q2[i], mf::detail::from_f64x2(subdd(mf::detail::to_f64x2(f1[i]), mf::detail::to_f64x2(f2[i]))), q1, f1);
  BENCH("mul",  REPS_FAST, q1[i] * q2[i], mf::detail::from_f64x2(muldd(mf::detail::to_f64x2(f1[i]), mf::detail::to_f64x2(f2[i]))), q1, f1);
  BENCH("div",  REPS_FAST, q1[i] / q2[i], mf::detail::from_f64x2(divdd(mf::detail::to_f64x2(f1[i]), mf::detail::to_f64x2(f2[i]))), q1, f1);
  BENCH1(sqrt,  REPS_FAST, qpos, fpos);
  BENCH("fma",  REPS_FAST, fx::fma_op(q1[i], q2[i], qpos[i]),
                           fx::fma_op(f1[i], f2[i], fpos[i]), q1, f1);
  BENCH("abs",  REPS_FAST, fx::abs_op(q1[i]), fx::abs_op(f1[i]), q1, f1);
  BENCH("neg",  REPS_FAST, fx::neg_op(q1[i]), fx::neg_op(f1[i]), q1, f1);
}

static void bench_rounding() {
  BENCH1(trunc, REPS_FAST, q1, f1);
  BENCH1(round, REPS_FAST, q1, f1);
}

static void bench_binary() {
  BENCH2(fmin,    REPS_FAST, q1, f1, q2, f2);
  BENCH2(fmax,    REPS_FAST, q1, f1, q2, f2);
  BENCH2(fdim,    REPS_FAST, q1, f1, q2, f2);
  BENCH2(copysign, REPS_FAST, q1, f1, q2, f2);
  BENCH2(fmod,    REPS_FAST, q1, f1, q2, f2);
  BENCH2(hypot,   REPS_FAST, q1, f1, q2, f2);
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
  BENCH1(sin,    REPS_TRIG, qsmall, fsmall);
  BENCH1(cos,    REPS_TRIG, qsmall, fsmall);
  BENCH1(tan,    REPS_TRIG, qbnd,   fbnd);
}

static void bench_pi_trig() {
  BENCH1(sinpi,  REPS_TRIG, qsmall, fsmall);
  BENCH1(cospi,  REPS_TRIG, qsmall, fsmall);
  BENCH1(tanpi,  REPS_TRIG, qbnd,   fbnd);
  BENCH1(asinpi, REPS_TRIG, qbnd,   fbnd);
  BENCH1(acospi, REPS_TRIG, qbnd,   fbnd);
  BENCH1(atanpi, REPS_TRIG, q1,     f1);
}

static void bench_inv_trig() {
  BENCH1(asin,   REPS_TRIG, qbnd, fbnd);
  BENCH1(acos,   REPS_TRIG, qbnd, fbnd);
  BENCH1(atan,   REPS_TRIG, q1,   f1);
  BENCH2(atan2,  REPS_TRIG, q1, f1, q2, f2);
}

static void bench_hyperbolic() {
  BENCH1(sinh,  REPS_TRIG, qsmall, fsmall);
  BENCH1(cosh,  REPS_TRIG, qsmall, fsmall);
  BENCH1(tanh,  REPS_TRIG, qsmall, fsmall);
}

static void bench_inv_hyperbolic() {
  BENCH1(asinh, REPS_TRIG, q1, f1);
  BENCH("acosh", REPS_TRIG,
        fx::acosh((q_t)1.0q + qpos[i]),
        fx::acosh(mf::float64x2(1.0) + fpos[i]), qpos, fpos);
  BENCH1(atanh, REPS_TRIG, qbnd, fbnd);
}

static void bench_erf_gamma() {
  BENCH1(erf,    REPS_VERY_SLOW, q1, f1);
  BENCH1(erfc,   REPS_VERY_SLOW, q1, f1);
  BENCH1(tgamma, REPS_VERY_SLOW, qpos, fpos);
  BENCH1(lgamma, REPS_VERY_SLOW, qpos, fpos);
}

static void bench_bessel() {
  BENCH("bj0", REPS_VERY_SLOW, fx::j0(q1[i]), fx::j0(f1[i]), q1, f1);
  BENCH("bj1", REPS_VERY_SLOW, fx::j1(q1[i]), fx::j1(f1[i]), q1, f1);
  BENCH_BESN("bjn", jn, 3, REPS_VERY_SLOW, q1, f1);
  BENCH("by0", REPS_VERY_SLOW, fx::y0(qpos[i]), fx::y0(fpos[i]), qpos, fpos);
  BENCH("by1", REPS_VERY_SLOW, fx::y1(qpos[i]), fx::y1(fpos[i]), qpos, fpos);
  BENCH_BESN("byn", yn, 3, REPS_VERY_SLOW, qpos, fpos);
}

static void bench_complex_arith() {
  BENCHZ2_OP(add, +);
  BENCHZ2_OP(sub, -);
  BENCHZ2_OP(mul, *);
  BENCHZ2_OP(div, /);
  // cabs returns real: use the real bench macro with complex inputs.
  BENCH("cdd_abs", REPS_FAST, fx::cabs(zq1[i]), fx::cabs(zf1[i]), q1, f1);
  BENCHZ("cdd_conjg", REPS_FAST, fx::cconjg(zq1[i]), fx::cconjg(zf1[i]), zq1, zf1);
}

static void bench_complex_trans() {
  BENCHZ1_FAST(sqrt, REPS_TRIG);
  BENCHZ1_FAST(exp,  REPS_TRIG);
  BENCHZ1_FAST(log,  REPS_TRIG);
  BENCHZ1_FAST(sin,  REPS_TRIG);
  BENCHZ1_FAST(cos,  REPS_TRIG);
  BENCHZ1_FAST(tan,  REPS_TRIG);
  BENCHZ1_FAST(sinh, REPS_TRIG);
  BENCHZ1_FAST(cosh, REPS_TRIG);
  BENCHZ1_FAST(tanh, REPS_TRIG);
  BENCHZ1_FAST(asin, REPS_VERY_SLOW);
  BENCHZ1_FAST(acos, REPS_VERY_SLOW);
  BENCHZ1_FAST(atan, REPS_VERY_SLOW);
  BENCHZ1_FAST(asinh, REPS_VERY_SLOW);
  BENCHZ1_FAST(acosh, REPS_VERY_SLOW);
  BENCHZ1_FAST(atanh, REPS_VERY_SLOW);
}

// =============================================================================
// main
// =============================================================================

int main() {
  init_data();

  std::printf("================================================================\n");
  std::printf(" multifloats C-ABI benchmark — N=%d elements per rep, batched\n", N);
  std::printf(" (timed through sindd / cdd_muldd / j0dd … register-return ABI)\n");
  std::printf("================================================================\n\n");
  std::printf(" %-22s %12s %10s %10s %10s\n", "op", "n_ops", "qp [s]", "mf [s]", "speedup");
  std::printf(" ----------------------------------------------------------------\n");

  bench_arith();
  bench_rounding();
  bench_binary();
  bench_exp_log();
  bench_trig();
  bench_pi_trig();
  bench_inv_trig();
  bench_hyperbolic();
  bench_inv_hyperbolic();
  bench_erf_gamma();
  bench_bessel();
  bench_complex_arith();
  bench_complex_trans();

  std::printf(" ----------------------------------------------------------------\n");
  std::printf(" sinks (must remain live): qp=%g  dd=%g  zq=%g  zf=%g\n",
              q_sink, f_sink, zq_sink, zf_sink);
  return 0;
}
