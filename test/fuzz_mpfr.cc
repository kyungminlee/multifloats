// MPFR-based precision fuzz for multifloats.
//
// Per iteration, draw random inputs; for each op compute three values:
//   ref_mp :  mpreal at 200-bit precision (ground truth)
//   ref_q  :  libquadmath at 113-bit precision
//   got_dd :  multifloats DD (~106-bit)
// and report, per op, the max and mean of both rel_err(q vs mp) and
// rel_err(dd vs mp). This separates the DD kernel's true error from the
// float128 reference floor that scalar fuzz.cc cannot distinguish.
//
// The three implementations follow a consistent naming pattern:
//   DD kernel      : NAMEdd    (multifloats C ABI)
//   libquadmath    : NAMEq     (113-bit __float128)
//   mpreal         : mpfr::NAME (200-bit reference)
// We exploit this by declaring overloads of NAME in namespace fx for all
// three argument types, so a CHECK(NAME) macro can call them uniformly.
//
// Built only when -DBUILD_MPFR_TESTS=ON.

#include "multifloats.hh"
#include "multifloats_c.h"
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
// Complex oracle at mpreal precision.
//
// MPFR proper has no complex functions. For each DD complex op we compose
// the reference by hand from real mpreal identities that match the C99
// Annex G branch cut convention (same as libquadmath). The DD complex
// kernels in src/multifloats_math.cc follow Annex G, so the oracle and
// the kernel share branch conventions.
// =============================================================================

struct CMp { mp_t re, im; };

static inline CMp c_add (CMp const &a, CMp const &b) { return {a.re + b.re, a.im + b.im}; }
static inline CMp c_sub (CMp const &a, CMp const &b) { return {a.re - b.re, a.im - b.im}; }
static inline CMp c_mul (CMp const &a, CMp const &b) {
  return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
static inline CMp c_div (CMp const &a, CMp const &b) {
  mp_t d = b.re * b.re + b.im * b.im;
  return {(a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d};
}
static inline CMp  c_conj(CMp const &z) { return {z.re, -z.im}; }
static inline mp_t c_abs (CMp const &z) { return mpfr::hypot(z.re, z.im); }
static inline mp_t c_arg (CMp const &z) { return mpfr::atan2(z.im, z.re); }

// C99 Annex G csqrt (Kahan form): branch cut on negative real axis,
// result in the right half-plane.
static CMp c_sqrt(CMp const &z) {
  mp_t az = c_abs(z);
  if (az == 0) return {mp_t(0), mp_t(0)};
  mp_t u, v;
  if (z.re >= 0) {
    u = mpfr::sqrt((az + z.re) / 2);
    v = z.im / (2 * u);
  } else {
    v = mpfr::sqrt((az - z.re) / 2);
    if (z.im < 0) v = -v;
    u = z.im / (2 * v);
  }
  return {u, v};
}

static CMp c_exp (CMp const &z) { mp_t e = mpfr::exp(z.re); return {e*mpfr::cos(z.im), e*mpfr::sin(z.im)}; }
static CMp c_log (CMp const &z) { return {mpfr::log(c_abs(z)), c_arg(z)}; }
static CMp c_sin (CMp const &z) { return {mpfr::sin(z.re)*mpfr::cosh(z.im),  mpfr::cos(z.re)*mpfr::sinh(z.im)}; }
static CMp c_cos (CMp const &z) { return {mpfr::cos(z.re)*mpfr::cosh(z.im), -mpfr::sin(z.re)*mpfr::sinh(z.im)}; }
static CMp c_sinh(CMp const &z) { return {mpfr::sinh(z.re)*mpfr::cos(z.im),  mpfr::cosh(z.re)*mpfr::sin(z.im)}; }
static CMp c_cosh(CMp const &z) { return {mpfr::cosh(z.re)*mpfr::cos(z.im),  mpfr::sinh(z.re)*mpfr::sin(z.im)}; }
static CMp c_tan (CMp const &z) { return c_div(c_sin(z),  c_cos(z)); }
static CMp c_tanh(CMp const &z) { return c_div(c_sinh(z), c_cosh(z)); }

// asin(z) = -i · log(iz + sqrt(1 - z²))
static CMp c_asin(CMp const &z) {
  CMp zz = c_mul(z, z);
  CMp arg = c_add({-z.im, z.re}, c_sqrt({mp_t(1) - zz.re, -zz.im}));
  CMp l = c_log(arg);
  return {l.im, -l.re};
}
// acos(z) = π/2 − asin(z)
static CMp c_acos(CMp const &z) {
  static const mp_t half_pi = mpfr::const_pi(multifloats_test::kMpfrPrec) / 2;
  CMp a = c_asin(z);
  return {half_pi - a.re, -a.im};
}
// atan(z) = (-i/2) · log((1 + iz)/(1 - iz))    [C99 Annex G principal]
static CMp c_atan(CMp const &z) {
  CMp num = {mp_t(1) - z.im,  z.re};
  CMp den = {mp_t(1) + z.im, -z.re};
  CMp l = c_log(c_div(num, den));
  return {l.im / 2, -l.re / 2};
}
// asinh(z) = log(z + sqrt(z² + 1))
static CMp c_asinh(CMp const &z) {
  CMp zz = c_mul(z, z);
  return c_log(c_add(z, c_sqrt({zz.re + 1, zz.im})));
}
// acosh(z) = log(z + sqrt(z-1)·sqrt(z+1))   [C99 Annex G principal]
static CMp c_acosh(CMp const &z) {
  CMp s = c_mul(c_sqrt({z.re - 1, z.im}), c_sqrt({z.re + 1, z.im}));
  return c_log(c_add(z, s));
}
// atanh(z) = ½ · log((1 + z)/(1 - z))
static CMp c_atanh(CMp const &z) {
  CMp l = c_log(c_div({mp_t(1) + z.re,  z.im}, {mp_t(1) - z.re, -z.im}));
  return {l.re / 2, l.im / 2};
}

// =============================================================================
// Input-type conversions (DD / qp / mpreal ↔ our carrier types).
// =============================================================================

// __complex128 helpers.
static inline __complex128 cq(q_t a, q_t b) {
  __complex128 z; __real__ z = a; __imag__ z = b; return z;
}
static inline complex64x2_t cd(mf::float64x2 const &a, mf::float64x2 const &b) {
  complex64x2_t z;
  z.re = mf::detail::to_f64x2(a);
  z.im = mf::detail::to_f64x2(b);
  return z;
}
static inline CMp to_cmp(complex64x2_t const &z) {
  return {to_mp(mf::detail::from_f64x2(z.re)),
          to_mp(mf::detail::from_f64x2(z.im))};
}
static inline CMp to_cmp(__complex128 z) { return {to_mp(crealq(z)), to_mp(cimagq(z))}; }

// =============================================================================
// namespace fx — overloaded wrappers over all three implementations.
//
// Each wrap macro below defines three inline overloads of the same short
// name, dispatching by argument type. Callers can then write fx::sin(x)
// regardless of whether x is DD / qp / mpreal. This is what drives the
// CHECK(NAME) macro pattern in the main loop.
// =============================================================================

namespace fx {

// -- Real unary: NAMEdd / NAMEq / mpfr::NAME ---------------------------------
#define WRAP_REAL_UNARY(NAME)                                                  \
  inline mf::float64x2 NAME(mf::float64x2 const &x) {                          \
    return mf::detail::from_f64x2(::NAME##dd(mf::detail::to_f64x2(x)));        \
  }                                                                            \
  inline q_t  NAME(q_t x)            { return ::NAME##q(x); }                  \
  inline mp_t NAME(mp_t const &x)    { return mpfr::NAME(x); }

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
WRAP_REAL_UNARY(exp)
WRAP_REAL_UNARY(expm1)
WRAP_REAL_UNARY(log)
WRAP_REAL_UNARY(log10)
WRAP_REAL_UNARY(log1p)
WRAP_REAL_UNARY(sqrt)
WRAP_REAL_UNARY(erf)
WRAP_REAL_UNARY(erfc)
#undef WRAP_REAL_UNARY

// -- Real unary with name exception (mpreal uses a different name) ----------
#define WRAP_REAL_UNARY_MPNAME(NAME, MPNAME)                                   \
  inline mf::float64x2 NAME(mf::float64x2 const &x) {                          \
    return mf::detail::from_f64x2(::NAME##dd(mf::detail::to_f64x2(x)));        \
  }                                                                            \
  inline q_t  NAME(q_t x)            { return ::NAME##q(x); }                  \
  inline mp_t NAME(mp_t const &x)    { return mpfr::MPNAME(x); }

WRAP_REAL_UNARY_MPNAME(tgamma, gamma)
WRAP_REAL_UNARY_MPNAME(lgamma, lngamma)
#undef WRAP_REAL_UNARY_MPNAME

// -- π-scaled trig: DD via C ABI, qp via mul-then-fn, mpreal via mul-then-fn.
//    {sin,cos,tan}pi(x) = fn(π·x);  {asin,acos,atan}pi(x) = fn(x)/π
#define WRAP_PI_INPUT(NAME, FN)                                                \
  inline mf::float64x2 NAME(mf::float64x2 const &x) {                          \
    return mf::detail::from_f64x2(::NAME##dd(mf::detail::to_f64x2(x)));        \
  }                                                                            \
  inline q_t  NAME(q_t x)            { return ::FN##q(q_t(M_PIq) * x); }       \
  inline mp_t NAME(mp_t const &x) {                                            \
    static const mp_t pi = mpfr::const_pi(multifloats_test::kMpfrPrec);        \
    return mpfr::FN(pi * x);                                                   \
  }

WRAP_PI_INPUT(sinpi, sin)
WRAP_PI_INPUT(cospi, cos)
WRAP_PI_INPUT(tanpi, tan)
#undef WRAP_PI_INPUT

#define WRAP_PI_OUTPUT(NAME, FN)                                               \
  inline mf::float64x2 NAME(mf::float64x2 const &x) {                          \
    return mf::detail::from_f64x2(::NAME##dd(mf::detail::to_f64x2(x)));        \
  }                                                                            \
  inline q_t  NAME(q_t x)            { return ::FN##q(x) / q_t(M_PIq); }       \
  inline mp_t NAME(mp_t const &x) {                                            \
    static const mp_t pi = mpfr::const_pi(multifloats_test::kMpfrPrec);        \
    return mpfr::FN(x) / pi;                                                   \
  }

WRAP_PI_OUTPUT(asinpi, asin)
WRAP_PI_OUTPUT(acospi, acos)
WRAP_PI_OUTPUT(atanpi, atan)
#undef WRAP_PI_OUTPUT

// -- Bessel (DD: {j0,j1,y0,y1}dd, qp: {j0,j1,y0,y1}q, mpreal: bessel*0/1) ---
#define WRAP_BESSEL01(NAME, MPNAME)                                            \
  inline mf::float64x2 NAME(mf::float64x2 const &x) {                          \
    return mf::detail::from_f64x2(::NAME##dd(mf::detail::to_f64x2(x)));        \
  }                                                                            \
  inline q_t  NAME(q_t x)            { return ::NAME##q(x); }                  \
  inline mp_t NAME(mp_t const &x)    { return mpfr::MPNAME(x); }

WRAP_BESSEL01(j0, besselj0)
WRAP_BESSEL01(j1, besselj1)
WRAP_BESSEL01(y0, bessely0)
WRAP_BESSEL01(y1, bessely1)
#undef WRAP_BESSEL01

// Bessel with integer order. DD takes int, qp's jnq takes int, mpreal
// expects long — accept int and widen.
#define WRAP_BESSELN(NAME, MPNAME)                                             \
  inline mf::float64x2 NAME(int n, mf::float64x2 const &x) {                   \
    return mf::detail::from_f64x2(::NAME##dd(n, mf::detail::to_f64x2(x)));     \
  }                                                                            \
  inline q_t  NAME(int n, q_t x)         { return ::NAME##q(n, x); }           \
  inline mp_t NAME(int n, mp_t const &x) { return mpfr::MPNAME((long)n, x); }

WRAP_BESSELN(jn, besseljn)
WRAP_BESSELN(yn, besselyn)
#undef WRAP_BESSELN

// -- Real binary ---------------------------------------------------------
#define WRAP_REAL_BINARY(NAME)                                                 \
  inline mf::float64x2 NAME(mf::float64x2 const &a, mf::float64x2 const &b) {  \
    return mf::detail::from_f64x2(                                             \
        ::NAME##dd(mf::detail::to_f64x2(a), mf::detail::to_f64x2(b)));         \
  }                                                                            \
  inline q_t  NAME(q_t a, q_t b)                 { return ::NAME##q(a, b); }   \
  inline mp_t NAME(mp_t const &a, mp_t const &b) { return mpfr::NAME(a, b); }

WRAP_REAL_BINARY(atan2)
WRAP_REAL_BINARY(hypot)
WRAP_REAL_BINARY(pow)
#undef WRAP_REAL_BINARY

// -- Complex unary: c<NAME>dd / c<NAME>q / c_<NAME> ---------------------
#define WRAP_COMPLEX_UNARY(NAME)                                               \
  inline complex64x2_t c##NAME(complex64x2_t z) { return ::c##NAME##dd(z); }   \
  inline __complex128  c##NAME(__complex128  z) { return ::c##NAME##q(z);  }   \
  inline CMp           c##NAME(CMp const  &z)   { return ::c_##NAME(z);    }

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

// -- Complex → real / real → complex ------------------------------------
// cabs: complex → real magnitude
inline mf::float64x2 cabs(complex64x2_t z) { return mf::detail::from_f64x2(::cabsdd(z)); }
inline q_t           cabs(__complex128  z) { return ::cabsq(z); }
inline mp_t          cabs(CMp const   &z)  { return c_abs(z); }

// cconjg: complex conjugate.
inline complex64x2_t cconjg(complex64x2_t z) { return ::conjdd(z); }
inline __complex128  cconjg(__complex128  z) { return ::conjq(z); }
inline CMp           cconjg(CMp const   &z)  { return c_conj(z); }

} // namespace fx

// =============================================================================
// Per-op statistics — two columns, one for each reference comparison.
// =============================================================================

struct StatEntry {
  char   name[24] = {};
  double max_q  = 0.0;  double sum_q  = 0.0;
  double max_dd = 0.0;  double sum_dd = 0.0;
  long   count = 0;
};

static constexpr int kMaxStats = 192;
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

// Subnormal floor. DD precision is bounded by IEEE-double's exponent range:
// the lo limb is itself a double, so it cannot go below 2^-1074. For the
// pair to carry full 106-bit precision we need lo at ~hi·2^-53 to land in
// the *normal* double range, i.e. hi ≥ 2^-969. Skipping samples with |ref|
// < 2^-969 keeps fuzz measuring kernel quality rather than this format-
// level cliff.
static bool mp_below_subnormal_floor(mp_t const &v) {
  static const mp_t floor_ = mpfr::pow(mp_t(2), mp_t(-969));
  return mpfr::abs(v) < floor_;
}

// Record a real scalar op.
static void check(char const *op, mf::float64x2 const &got_dd, q_t ref_q, mp_t const &ref_mp) {
  if (!mp_isfinite(ref_mp) || mp_isnan(ref_mp)) return;
  if (!q_isfinite(ref_q)) return;
  if (!std::isfinite(got_dd._limbs[0])) return;
  if (mp_below_subnormal_floor(ref_mp)) return;
  double rq  = mp_rel_err(to_mp(ref_q),   ref_mp);
  double rdd = mp_rel_err(to_mp(got_dd),  ref_mp);
  record(op, rq, rdd);
}

// Record a complex op. Emits "<op>_re" and "<op>_im" rows so ops.py's
// [("re", "<key>_re"), ("im", "<key>_im")] fuzz spec finds them directly.
static void check_cplx(char const *op, complex64x2_t const &got_dd,
                       __complex128 ref_q, CMp const &ref_mp) {
  if (!mp_isfinite(ref_mp.re) || !mp_isfinite(ref_mp.im)) return;
  q_t qre = crealq(ref_q), qim = cimagq(ref_q);
  if (!q_isfinite(qre) || !q_isfinite(qim)) return;
  if (!std::isfinite(got_dd.re.hi) || !std::isfinite(got_dd.im.hi)) return;

  char key_re[32], key_im[32];
  std::snprintf(key_re, sizeof(key_re), "%s_re", op);
  std::snprintf(key_im, sizeof(key_im), "%s_im", op);

  if (!mp_below_subnormal_floor(ref_mp.re)) {
    record(key_re,
           mp_rel_err(to_mp(qre),                                 ref_mp.re),
           mp_rel_err(to_mp(mf::detail::from_f64x2(got_dd.re)),   ref_mp.re));
  }
  if (!mp_below_subnormal_floor(ref_mp.im)) {
    record(key_im,
           mp_rel_err(to_mp(qim),                                 ref_mp.im),
           mp_rel_err(to_mp(mf::detail::from_f64x2(got_dd.im)),   ref_mp.im));
  }
}

// =============================================================================
// CHECK macros — read inputs from the enclosing scope (f1/f2/q1/q2/m1/m2 for
// scalars, zd1/zd2/zq1/zq2/zm1/zm2 for complex). Each macro calls all three
// typed overloads via fx::NAME and invokes check / check_cplx.
// =============================================================================

// Real unary (input: f1/q1/m1).
#define CHECK1(NAME) check(#NAME, fx::NAME(f1), fx::NAME(q1), fx::NAME(m1))
// Same, but override the printed label (used for ops whose bench_key
// doesn't match the generic name — e.g. Bessel: j0dd → "bj0").
#define CHECK1_LABEL(LABEL, NAME) check(LABEL, fx::NAME(f1), fx::NAME(q1), fx::NAME(m1))
// Real binary (inputs: f1,f2 / q1,q2 / m1,m2).
#define CHECK2(NAME) check(#NAME, fx::NAME(f1, f2), fx::NAME(q1, q2), fx::NAME(m1, m2))
// Bessel with integer order.
#define CHECK_BESN(LABEL, NAME, N) \
  check(LABEL, fx::NAME(N, f1), fx::NAME(N, q1), fx::NAME(N, m1))

// Complex unary (input: zd1/zq1/zm1).
#define CHECK_C1(STEM) check_cplx("cdd_" #STEM, \
  fx::c##STEM(zd1), fx::c##STEM(zq1), fx::c##STEM(zm1))
// Complex binary via operator (cadd, csub, cmul, cdiv).
#define CHECK_C2_OP(STEM, OP) check_cplx("cdd_" #STEM, \
  c##STEM##dd(zd1, zd2), zq1 OP zq2, c_##STEM(zm1, zm2))
// cabs: complex → real scalar.
#define CHECK_CABS() check("cdd_abs", fx::cabs(zd1), fx::cabs(zq1), fx::cabs(zm1))

// =============================================================================
// Random input generation. Real ops get wide magnitudes to stress the kernel;
// complex ops use a narrower band so mul/div don't hit catastrophic
// cancellation in the oracle's hand-rolled identities.
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
  // Narrow: |q| in [1e-3, 1e3]. Used for complex inputs.
  q_t narrow(double r, double rexp) {
    int k = (int)(rexp * 6.0) - 3;
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
  // Generate complex input pair (re, im) with narrow magnitudes.
  void generate_narrow_pair(q_t &re, q_t &im) {
    re = narrow(u(), u());
    im = narrow(u(), u());
  }
};

// =============================================================================
// Output
// =============================================================================

static void print_report() {
  std::printf("\n");
  std::printf("Per-operation 3-way precision report (reference = mpreal @ 200 bits):\n");
  std::printf("  %-14s %10s %14s %14s %14s %14s\n",
              "op", "n", "max_q", "mean_q", "max_dd", "mean_dd");
  for (int i = 0; i < g_nstats; ++i) {
    StatEntry const &s = g_stats[i];
    if (s.count == 0) continue;
    double invn = 1.0 / (double)s.count;
    std::printf("  %-14s %10ld  %14.3e %14.3e %14.3e %14.3e\n",
                s.name, s.count,
                s.max_q, s.sum_q * invn,
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
  uint64_t seed   = 42ULL;
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

    // Inputs live at DD precision: all three paths start from the SAME
    // value so "error" is the op kernel, not input rounding.
    mf::float64x2 f1 = from_q(qseed1), f2 = from_q(qseed2);
    q_t  q1 = to_q(f1),  q2 = to_q(f2);
    mp_t m1 = to_mp(f1), m2 = to_mp(f2);

    // ---------------- arithmetic ----------------
    check("add", f1 + f2, q1 + q2, m1 + m2);
    check("sub", f1 - f2, q1 - q2, m1 - m2);
    check("mul", f1 * f2, q1 * q2, m1 * m2);
    if (q2 != (q_t)0) check("div", f1 / f2, q1 / q2, m1 / m2);
    if (q1 >= (q_t)0) check("sqrt", mf::sqrt(f1), sqrtq(q1), mpfr::sqrt(m1));

    // ---------------- binary (every 5) ----------------
    if (i % 5 == 0) {
      CHECK2(hypot);
      check("fmin", mf::fmin(f1, f2), q1 < q2 ? q1 : q2, q1 < q2 ? m1 : m2);
      check("fmax", mf::fmax(f1, f2), q1 < q2 ? q2 : q1, q1 < q2 ? m2 : m1);
    }

    // ---------------- transcendentals (every 100) ----------------
    if (i % 100 != 0) continue;
    q_t aq1 = q1 < 0 ? -q1 : q1;

    // exp / log family
    {
      q_t qe = expq(q1);
      if (q_isfinite(qe)) CHECK1(exp);
      q_t qem = expm1q(q1);
      if (q_isfinite(qem)) CHECK1(expm1);
    }
    if (q1 > (q_t)0)  { CHECK1(log); CHECK1(log10); }
    if (q1 > (q_t)-1) { CHECK1(log1p); }

    // trig
    if (aq1 < (q_t)1e6q) {
      CHECK1(sin); CHECK1(cos);
      if (fabsq(cosq(q1)) > (q_t)1e-12q) CHECK1(tan);
    }
    if (aq1 <= (q_t)1) { CHECK1(asin); CHECK1(acos); }
    CHECK1(atan);
    CHECK2(atan2);

    // π-scaled trig (mpreal carries π cleanly; qp loses 1 ulp on the mul).
    if (aq1 < (q_t)1e6q) {
      CHECK1(sinpi); CHECK1(cospi);
      if (fabsq(cosq(q_t(M_PIq) * q1)) > (q_t)1e-12q) CHECK1(tanpi);
    }
    if (aq1 <= (q_t)1) { CHECK1(asinpi); CHECK1(acospi); }
    CHECK1(atanpi);

    // hyperbolic
    if (aq1 < (q_t)700) { CHECK1(sinh); CHECK1(cosh); CHECK1(tanh); }
    CHECK1(asinh);
    if (q1 >= (q_t)1) CHECK1(acosh);
    if (aq1 < (q_t)1) CHECK1(atanh);

    // special (erf / gamma)
    if (aq1 < (q_t)100)                    { CHECK1(erf); CHECK1(erfc); }
    if (q1 > (q_t)0 && q1 < (q_t)100)      { CHECK1(tgamma); CHECK1(lgamma); }

    // pow
    q_t aq2 = q2 < 0 ? -q2 : q2;
    if (q1 > (q_t)1e-3q && q1 < (q_t)1e3q && aq2 < (q_t)30) CHECK2(pow);

    // ---------------- Bessel (labels use ops.py's bench_key convention) ----------------
    if (aq1 < (q_t)200) {
      CHECK1_LABEL("bj0", j0);
      CHECK1_LABEL("bj1", j1);
      CHECK_BESN("bjn", jn, 3);
      if (q1 > (q_t)0) {
        CHECK1_LABEL("by0", y0);
        CHECK1_LABEL("by1", y1);
        CHECK_BESN("byn", yn, 3);
      }
    }

    // ---------------- Complex (every 100, narrower range) ----------------
    // Generate fresh narrow-magnitude inputs so cmul/cdiv don't cancel.
    // z_r1, z_i1, z_r2, z_i2 → complex pairs z1, z2.
    q_t qr1, qi1, qr2, qi2;
    rng.generate_narrow_pair(qr1, qi1);
    rng.generate_narrow_pair(qr2, qi2);
    if (!q_isfinite(qr1) || !q_isfinite(qi1) ||
        !q_isfinite(qr2) || !q_isfinite(qi2)) continue;

    mf::float64x2 fr1 = from_q(qr1), fi1 = from_q(qi1);
    mf::float64x2 fr2 = from_q(qr2), fi2 = from_q(qi2);
    q_t           qr1_ = to_q(fr1), qi1_ = to_q(fi1), qr2_ = to_q(fr2), qi2_ = to_q(fi2);
    mp_t          mr1 = to_mp(fr1), mi1 = to_mp(fi1), mr2 = to_mp(fr2), mi2 = to_mp(fi2);

    complex64x2_t zd1 = cd(fr1, fi1);
    complex64x2_t zd2 = cd(fr2, fi2);
    __complex128  zq1 = cq(qr1_, qi1_);
    __complex128  zq2 = cq(qr2_, qi2_);
    CMp           zm1 = {mr1, mi1};
    CMp           zm2 = {mr2, mi2};

    // Complex arithmetic (z1 op z2).
    CHECK_C2_OP(add, +);
    CHECK_C2_OP(sub, -);
    CHECK_C2_OP(mul, *);
    { mp_t dmag = mr2 * mr2 + mi2 * mi2;
      if (dmag > mp_t(0)) CHECK_C2_OP(div, /); }
    CHECK_CABS();
    CHECK_C1(conjg);

    // Complex transcendentals.
    CHECK_C1(sqrt);
    CHECK_C1(exp);
    CHECK_C1(sin); CHECK_C1(cos);
    CHECK_C1(sinh); CHECK_C1(cosh);
    {
      mp_t az = c_abs(zm1);
      if (az > mp_t(1e-200) && az < mp_t(1e100)) CHECK_C1(log);
    }
    if (c_abs(c_cos(zm1))  > mp_t(1e-10)) CHECK_C1(tan);
    if (c_abs(c_cosh(zm1)) > mp_t(1e-10)) CHECK_C1(tanh);
    CHECK_C1(asin);
    CHECK_C1(acos);
    CHECK_C1(atan);
    CHECK_C1(asinh);
    CHECK_C1(acosh);
    CHECK_C1(atanh);
  }

  print_report();
  return 0;
}
