// Timing benchmark through the multifloats C ABI.
//
// This file times every kernel exposed in include/multifloats.h (sindd,
// cdd_muldd, j0dd, matmuldd_mm, …) against libquadmath as a reference
// for the SAME operation. The DD leg goes through the register-return
// C ABI (not the C++ header-only inline path), so the timings reflect
// what a real C or Fortran-bind-c caller sees.
//
// Precision measurement lives in test/fuzz.cc (compiled with -DUSE_MPFR);
// this file is timing-only.
//
// The three function families follow a consistent suffix scheme —
//   DD kernel      : NAMEdd     (C ABI)
//   libquadmath    : NAMEq      (113-bit __float128)
// — which we exploit by declaring overloads in `namespace fx` so a single
// BENCH1(NAME, ...) macro can dispatch both legs uniformly.

#include "multifloats.h"
#include "test_common.hh"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>

namespace mf = multifloats;
using multifloats::float64x2;
using multifloats::complex64x2;
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
    return mf::float64x2(NAME##dd(static_cast<float64x2>(x)));        \
  }                                                                            \
  inline q_t NAME(q_t x) { return ::NAME##q(x); }

// neg: libquadmath uses the unary `-` operator on q_t, not a `negq` fn.
inline mf::float64x2 neg_op(mf::float64x2 const &x) {
  return mf::float64x2(negdd(static_cast<float64x2>(x)));
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
// erfcx: libquadmath has no direct variant; compose from erfc and exp(x²).
// This is the same oracle fuzz.cc uses (in both its default and USE_MPFR modes).
inline mf::float64x2 erfcx(mf::float64x2 const &x) {
  return mf::float64x2(erfcxdd(static_cast<float64x2>(x)));
}
inline q_t erfcx(q_t x) { return ::expq(x * x) * ::erfcq(x); }
WRAP_REAL_UNARY(tgamma)
WRAP_REAL_UNARY(lgamma)
// fabs: libquadmath is fabsq; DD is fabsdd. "abs" common name for readability.
inline mf::float64x2 abs_op(mf::float64x2 const &x) {
  return mf::float64x2(fabsdd(static_cast<float64x2>(x)));
}
inline q_t abs_op(q_t x) { return ::fabsq(x); }
#undef WRAP_REAL_UNARY

// -- π-scaled trig (DD: NAMEdd, qp: library fn composed with π·x or /π) ----
#define WRAP_PI_INPUT(NAME, FN)                                                \
  inline mf::float64x2 NAME(mf::float64x2 const &x) {                          \
    return mf::float64x2(NAME##dd(static_cast<float64x2>(x)));        \
  }                                                                            \
  inline q_t NAME(q_t x) { return ::FN##q(q_t(M_PIq) * x); }

WRAP_PI_INPUT(sinpi, sin)
WRAP_PI_INPUT(cospi, cos)
WRAP_PI_INPUT(tanpi, tan)
#undef WRAP_PI_INPUT

#define WRAP_PI_OUTPUT(NAME, FN)                                               \
  inline mf::float64x2 NAME(mf::float64x2 const &x) {                          \
    return mf::float64x2(NAME##dd(static_cast<float64x2>(x)));        \
  }                                                                            \
  inline q_t NAME(q_t x) { return ::FN##q(x) / q_t(M_PIq); }

WRAP_PI_OUTPUT(asinpi, asin)
WRAP_PI_OUTPUT(acospi, acos)
WRAP_PI_OUTPUT(atanpi, atan)
#undef WRAP_PI_OUTPUT

// -- Bessel (DD: NAMEdd, qp: NAMEq, both identical suffix scheme) ----------
#define WRAP_BESSEL01(NAME)                                                    \
  inline mf::float64x2 NAME(mf::float64x2 const &x) {                          \
    return mf::float64x2(NAME##dd(static_cast<float64x2>(x)));        \
  }                                                                            \
  inline q_t NAME(q_t x) { return ::NAME##q(x); }

WRAP_BESSEL01(j0)
WRAP_BESSEL01(j1)
WRAP_BESSEL01(y0)
WRAP_BESSEL01(y1)
#undef WRAP_BESSEL01

#define WRAP_BESSELN(NAME)                                                     \
  inline mf::float64x2 NAME(int n, mf::float64x2 const &x) {                   \
    return mf::float64x2(NAME##dd(n, static_cast<float64x2>(x)));     \
  }                                                                            \
  inline q_t NAME(int n, q_t x) { return ::NAME##q(n, x); }

WRAP_BESSELN(jn)
WRAP_BESSELN(yn)
#undef WRAP_BESSELN

// -- Real binary -----------------------------------------------------------
#define WRAP_REAL_BINARY(NAME)                                                 \
  inline mf::float64x2 NAME(mf::float64x2 const &a, mf::float64x2 const &b) {  \
    return mf::float64x2(                                             \
        NAME##dd(static_cast<float64x2>(a), static_cast<float64x2>(b)));         \
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

// atan2pi: DD native, qp composed as atan2q(y, x) / M_PIq.
inline mf::float64x2 atan2pi(mf::float64x2 const &y, mf::float64x2 const &x) {
  return mf::float64x2(atan2pidd(
      static_cast<float64x2>(y), static_cast<float64x2>(x)));
}
inline q_t atan2pi(q_t y, q_t x) { return ::atan2q(y, x) / (q_t)M_PIq; }

// -- Fused sincos / sinhcosh: one DD call produces both outputs. Returning
//    a single value (hi+lo-limb sum) keeps the anti-DCE feedback pattern
//    usable while still forcing both outputs to be computed. The qp leg
//    does two separate sinq/cosq calls for a fair "single call with two
//    outputs vs two calls" comparison.
inline mf::float64x2 sincos_op(mf::float64x2 const &x) {
  float64x2 s, c;
  sincosdd(static_cast<float64x2>(x), &s, &c);
  return mf::float64x2(adddd(s, c));
}
inline q_t sincos_op(q_t x) { return ::sinq(x) + ::cosq(x); }

inline mf::float64x2 sinhcosh_op(mf::float64x2 const &x) {
  float64x2 s, c;
  sinhcoshdd(static_cast<float64x2>(x), &s, &c);
  return mf::float64x2(adddd(s, c));
}
inline q_t sinhcosh_op(q_t x) { return ::sinhq(x) + ::coshq(x); }

// yn_range: one yndd_range call produces 6 outputs; the qp leg runs 6
// separate ynq calls. Sum all 6 into one scalar so the drain consumes
// everything and the compiler can't elide any output.
inline mf::float64x2 yn_range_op(mf::float64x2 const &x) {
  float64x2 out[6];
  yndd_range(0, 5, static_cast<float64x2>(x), out);
  float64x2 acc = out[0];
  for (int k = 1; k < 6; ++k) acc = adddd(acc, out[k]);
  return mf::float64x2(acc);
}
inline q_t yn_range_op(q_t x) {
  q_t acc = ::ynq(0, x);
  for (int k = 1; k < 6; ++k) acc += ::ynq(k, x);
  return acc;
}

// fma: (a, b, c) → a·b + c
inline mf::float64x2 fma_op(mf::float64x2 const &a, mf::float64x2 const &b, mf::float64x2 const &c) {
  return mf::float64x2(fmadd(static_cast<float64x2>(a),
                                         static_cast<float64x2>(b),
                                         static_cast<float64x2>(c)));
}
inline q_t fma_op(q_t a, q_t b, q_t c) { return ::fmaq(a, b, c); }

// -- Complex unary: c<NAME>dd / c<NAME>q -----------------------------------
#define WRAP_COMPLEX_UNARY(NAME)                                               \
  inline complex64x2 c##NAME(complex64x2 z) { return c##NAME##dd(z); }   \
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
inline mf::float64x2 cabs(complex64x2 z) { return mf::float64x2(cabsdd(z)); }
inline q_t           cabs(__complex128  z) { return ::cabsq(z); }

// carg: complex → real (phase angle in (-π, π]).
inline mf::float64x2 carg(complex64x2 z) { return mf::float64x2(cargdd(z)); }
inline q_t           carg(__complex128  z) { return ::cargq(z); }

// conjg: complex → complex
inline complex64x2 cconjg(complex64x2 z) { return conjdd(z); }
inline __complex128  cconjg(__complex128  z) { return ::conjq(z); }

// cproj: Riemann sphere projection.
inline complex64x2 cproj(complex64x2 z) { return cprojdd(z); }
inline __complex128  cproj(__complex128  z) { return ::cprojq(z); }

// -- Complex composite-oracle ops: libquadmath has no direct variant
//    (cexpm1q, clog2q, clog10q, clog1pq, csinpiq, ccospiq are all
//    missing). Compose the qp reference from primitives.

inline complex64x2 cexpm1(complex64x2 z) { return cexpm1dd(z); }
inline __complex128  cexpm1(__complex128 z) {
  __complex128 one; __real__ one = 1; __imag__ one = 0;
  return ::cexpq(z) - one;
}

inline complex64x2 clog1p(complex64x2 z) { return clog1pdd(z); }
inline __complex128  clog1p(__complex128 z) {
  __complex128 one; __real__ one = 1; __imag__ one = 0;
  return ::clogq(z + one);
}

inline complex64x2 clog2(complex64x2 z) { return clog2dd(z); }
inline __complex128  clog2(__complex128 z) {
  __complex128 lg = ::clogq(z);
  q_t l2 = ::logq((q_t)2);
  __complex128 r;
  __real__ r = crealq(lg) / l2;
  __imag__ r = cimagq(lg) / l2;
  return r;
}

inline complex64x2 clog10(complex64x2 z) { return clog10dd(z); }
inline __complex128  clog10(__complex128 z) {
  __complex128 lg = ::clogq(z);
  q_t l10 = ::logq((q_t)10);
  __complex128 r;
  __real__ r = crealq(lg) / l10;
  __imag__ r = cimagq(lg) / l10;
  return r;
}

inline complex64x2 csinpi(complex64x2 z) { return csinpidd(z); }
inline __complex128  csinpi(__complex128 z) {
  __complex128 pi; __real__ pi = (q_t)M_PIq; __imag__ pi = 0;
  return ::csinq(pi * z);
}

inline complex64x2 ccospi(complex64x2 z) { return ccospidd(z); }
inline __complex128  ccospi(__complex128 z) {
  __complex128 pi; __real__ pi = (q_t)M_PIq; __imag__ pi = 0;
  return ::ccosq(pi * z);
}

// cpow: both libquadmath (cpowq) and DD (cpowdd) provide native binary
// entry points, so a generic two-input bench macro fits.
inline complex64x2 cpow(complex64x2 z, complex64x2 w) { return cpowdd(z, w); }
inline __complex128  cpow(__complex128 z, __complex128 w)   { return ::cpowq(z, w); }

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
static complex64x2 zf1[N], zf2[N], zfres[N];

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
    zf1[i] = complex64x2{fre, fim};
    auto [qre2, fre2] = make_dd(away(u(rng)));
    auto [qim2, fim2] = make_dd(away(u(rng)));
    __real__ zq2[i] = qre2; __imag__ zq2[i] = qim2;
    zf2[i] = complex64x2{fre2, fim2};
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
  for (int i = 0; i < N; ++i) s += fres[i].limbs[0] + fres[i].limbs[1];
  in[0].limbs[0] += s * DRAIN_FEEDBACK_SCALE;
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
static void zffeed(complex64x2 *in) {
  double s = 0.0;
  for (int i = 0; i < N; ++i)
    s += zfres[i].real().limbs[0] + zfres[i].real().limbs[1] + zfres[i].imag().limbs[0] + zfres[i].imag().limbs[1];
  mf::float64x2 r = in[0].real();
  r.limbs[0] += s * DRAIN_FEEDBACK_SCALE;
  in[0].real(r);
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
  BENCH("add",  REPS_FAST, q1[i] + q2[i], mf::float64x2(adddd(static_cast<float64x2>(f1[i]), static_cast<float64x2>(f2[i]))), q1, f1);
  BENCH("sub",  REPS_FAST, q1[i] - q2[i], mf::float64x2(subdd(static_cast<float64x2>(f1[i]), static_cast<float64x2>(f2[i]))), q1, f1);
  BENCH("mul",  REPS_FAST, q1[i] * q2[i], mf::float64x2(muldd(static_cast<float64x2>(f1[i]), static_cast<float64x2>(f2[i]))), q1, f1);
  BENCH("div",  REPS_FAST, q1[i] / q2[i], mf::float64x2(divdd(static_cast<float64x2>(f1[i]), static_cast<float64x2>(f2[i]))), q1, f1);
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
  BENCH("sincos", REPS_TRIG,
        fx::sincos_op(qsmall[i]), fx::sincos_op(fsmall[i]), qsmall, fsmall);
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
  BENCH2(atan2pi, REPS_TRIG, q1, f1, q2, f2);
}

static void bench_hyperbolic() {
  BENCH1(sinh,  REPS_TRIG, qsmall, fsmall);
  BENCH1(cosh,  REPS_TRIG, qsmall, fsmall);
  BENCH1(tanh,  REPS_TRIG, qsmall, fsmall);
  BENCH("sinhcosh", REPS_TRIG,
        fx::sinhcosh_op(qsmall[i]), fx::sinhcosh_op(fsmall[i]), qsmall, fsmall);
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
  // erfcx: use qbnd (|x| < 1) so the qp composite oracle (exp(x²)·erfc(x))
  // stays well inside qp range. The DD kernel's selling point — no tail
  // cancellation — kicks in harder at larger |x|, but those inputs need
  // the q-input range gated to keep exp(x²) representable.
  BENCH("erfcx", REPS_VERY_SLOW,
        fx::erfcx(qbnd[i]), fx::erfcx(fbnd[i]), qbnd, fbnd);
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

  // yn_range(0..5): one DD call fills six outputs vs six qp ynq calls.
  // Work-per-call is 6× a single yn, so n_ops reflects the speedup of
  // the fused sweep.
  BENCH("yn_range(0..5)", REPS_VERY_SLOW,
        fx::yn_range_op(qpos[i]), fx::yn_range_op(fpos[i]), qpos, fpos);
}

static void bench_complex_arith() {
  BENCHZ2_OP(add, +);
  BENCHZ2_OP(sub, -);
  BENCHZ2_OP(mul, *);
  BENCHZ2_OP(div, /);
  // cabs / carg return real: use the real bench macro with complex inputs.
  BENCH("cdd_abs", REPS_FAST, fx::cabs(zq1[i]), fx::cabs(zf1[i]), q1, f1);
  BENCH("cdd_arg", REPS_FAST, fx::carg(zq1[i]), fx::carg(zf1[i]), q1, f1);
  BENCHZ("cdd_conjg", REPS_FAST, fx::cconjg(zq1[i]), fx::cconjg(zf1[i]), zq1, zf1);
  BENCHZ("cdd_proj",  REPS_FAST, fx::cproj(zq1[i]),  fx::cproj(zf1[i]),  zq1, zf1);
}

static void bench_complex_trans() {
  BENCHZ1_FAST(sqrt, REPS_TRIG);
  BENCHZ1_FAST(exp,  REPS_TRIG);
  BENCHZ1_FAST(expm1, REPS_TRIG);
  BENCHZ1_FAST(log,  REPS_TRIG);
  BENCHZ1_FAST(log2, REPS_TRIG);
  BENCHZ1_FAST(log10, REPS_TRIG);
  BENCHZ1_FAST(log1p, REPS_TRIG);
  BENCHZ1_FAST(sin,  REPS_TRIG);
  BENCHZ1_FAST(cos,  REPS_TRIG);
  BENCHZ1_FAST(tan,  REPS_TRIG);
  BENCHZ1_FAST(sinpi, REPS_TRIG);
  BENCHZ1_FAST(cospi, REPS_TRIG);
  BENCHZ1_FAST(sinh, REPS_TRIG);
  BENCHZ1_FAST(cosh, REPS_TRIG);
  BENCHZ1_FAST(tanh, REPS_TRIG);
  BENCHZ1_FAST(asin, REPS_VERY_SLOW);
  BENCHZ1_FAST(acos, REPS_VERY_SLOW);
  BENCHZ1_FAST(atan, REPS_VERY_SLOW);
  BENCHZ1_FAST(asinh, REPS_VERY_SLOW);
  BENCHZ1_FAST(acosh, REPS_VERY_SLOW);
  BENCHZ1_FAST(atanh, REPS_VERY_SLOW);
  // cpow: binary complex. BENCHZ1_FAST assumes unary, so explicit form.
  BENCHZ("cdd_pow", REPS_VERY_SLOW,
         fx::cpow(zq1[i], zq2[i]), fx::cpow(zf1[i], zf2[i]), zq1, zf1);
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
