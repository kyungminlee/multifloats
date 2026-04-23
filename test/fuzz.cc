// Property-based fuzz tests for multifloats.h.
//
// This is the C++ analogue of test/fuzz.f90: for each iteration a random
// __float128 pair (q1, q2) is generated, projected to float64x2 inputs
// (f1, f2), and every operation that multifloats.h exposes is run on
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
//   1e-15   single-double ops and libm-quality transcendentals (~5 ulp of dp)
// Built with g++ + libquadmath.
//
// ## USE_MPFR mode
//
// When compiled with `-DUSE_MPFR`, the binary additionally pulls in mpreal
// (MPFR at 200 bits) and reports a 3-way precision report: for each op, the
// per-iteration relative error of both the libquadmath reference and the
// multifloats DD result is measured against the 200-bit mpreal oracle, and
// max/mean columns of each are printed. Pass/fail gating is unchanged —
// the MPFR oracle is purely informational and does not participate in the
// regression-gate decision. The mpreal pieces — includes, the complex
// oracle, the extra stats columns, the extra check() parameter — are all
// reached through the `MP_PARAM` / `MP_ARG` / `MP_STMT` macros defined in
// one block near the top of this file, so the main loop body stays free
// of `#ifdef` sprinkles.

#include "multifloats.h"
#include "test_common.hh"

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
using multifloats_test::q_rel_err;
using multifloats_test::q_isfinite;
using multifloats_test::q_isnan;
using multifloats_test::qstr;

// Inline wrappers so complex arithmetic on __complex128 reads as
// caddq/csubq/cmulq/cdivq function calls, matching the cSTEMdd /
// cSTEMmp naming used by the DD kernels and mp oracle. This lets
// CHK_C2 treat the qp leg the same way as the other two and removes
// the need for a separate operator-taking macro.
static inline __complex128 caddq(__complex128 a, __complex128 b) { return a + b; }
static inline __complex128 csubq(__complex128 a, __complex128 b) { return a - b; }
static inline __complex128 cmulq(__complex128 a, __complex128 b) { return a * b; }
static inline __complex128 cdivq(__complex128 a, __complex128 b) { return a / b; }

// mpreal spells these gamma / lngamma; the C99 standard (and multifloats /
// libquadmath) use tgamma / lgamma. Since mpreal's names are also what
// the qp and DD sides can be made to match via one-line aliases, use
// the mpreal spelling everywhere so CHK1(gamma) / CHK1(lngamma) derive
// all three legs uniformly.
namespace multifloats {
inline float64x2 gamma  (float64x2 const &x) { return tgamma(x); }
inline float64x2 lngamma(float64x2 const &x) { return lgamma(x); }
}  // namespace multifloats
static inline q_t gammaq  (q_t x) { return tgammaq(x); }
static inline q_t lngammaq(q_t x) { return lgammaq(x); }

// =============================================================================
// USE_MPFR mode — all mpreal-specific types, helpers, the complex oracle,
// and the three call-site macros (MP_PARAM / MP_ARG / MP_STMT) live in this
// single block. Everything outside the block is mode-agnostic; the rest of
// the file compiles identically under both settings. The macros express
// "this argument only exists in one mode," which templates can't convey
// without stubbing every op in a parallel no-op class. __VA_ARGS__ lets the
// mp expression contain unparenthesized commas (e.g. `CMp{a, b}`).
// =============================================================================
#ifdef USE_MPFR

#include "test_common_mpfr.hh"

using multifloats_test::mp_t;
using multifloats_test::to_mp;
using multifloats_test::mp_rel_err;
using multifloats_test::mp_isfinite;
using multifloats_test::mp_isnan;

// mpreal covers fmin/fmax/fmod/copysign but is missing fdim. Fill the gap
// with a local wrapper so CHK2(fdim) derives mpfr::fdim(m1, m2) the same
// way as the other scalar binary ops.
namespace mpfr {
inline mp_t fdim(mp_t const &a, mp_t const &b) { return a > b ? a - b : mp_t(0); }
}  // namespace mpfr

#define MP_PARAM(...) , __VA_ARGS__
#define MP_ARG(...)   , (__VA_ARGS__)
#define MP_STMT(...)  __VA_ARGS__

// MPFR proper has no complex functions. For each DD complex op we compose
// the reference by hand from real mpreal identities that match the C99
// Annex G branch cut convention (same as libquadmath). The DD complex
// kernels in src/multifloats_math.cc follow Annex G, so the oracle and
// the kernel share branch conventions.
struct CMp { mp_t re, im; };

// Naming convention for the mp oracle: cNAMEmp(z), to mirror the DD
// kernel's cNAMEdd(z) and libquadmath's cNAMEq(z). This keeps CHK_C1 /
// CHK_COP below able to derive all three function names from a single
// STEM token via `## dd`, `## q`, `## mp`.
static inline CMp caddmp (CMp const &a, CMp const &b) { return {a.re + b.re, a.im + b.im}; }
static inline CMp csubmp (CMp const &a, CMp const &b) { return {a.re - b.re, a.im - b.im}; }
static inline CMp cmulmp (CMp const &a, CMp const &b) {
  return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
static inline CMp cdivmp (CMp const &a, CMp const &b) {
  mp_t d = b.re * b.re + b.im * b.im;
  return {(a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d};
}
static inline CMp  cconjgmp(CMp const &z) { return {z.re, -z.im}; }
static inline mp_t cabsmp (CMp const &z) { return mpfr::hypot(z.re, z.im); }
static inline mp_t cargmp (CMp const &z) { return mpfr::atan2(z.im, z.re); }

// C99 Annex G csqrt (Kahan form): branch cut on negative real axis,
// result in the right half-plane.
static CMp csqrtmp(CMp const &z) {
  mp_t az = cabsmp(z);
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

static CMp cexpmp (CMp const &z) { mp_t e = mpfr::exp(z.re); return {e*mpfr::cos(z.im), e*mpfr::sin(z.im)}; }
static CMp clogmp (CMp const &z) { return {mpfr::log(cabsmp(z)), cargmp(z)}; }
static CMp csinmp (CMp const &z) { return {mpfr::sin(z.re)*mpfr::cosh(z.im),  mpfr::cos(z.re)*mpfr::sinh(z.im)}; }
static CMp ccosmp (CMp const &z) { return {mpfr::cos(z.re)*mpfr::cosh(z.im), -mpfr::sin(z.re)*mpfr::sinh(z.im)}; }
static CMp csinhmp(CMp const &z) { return {mpfr::sinh(z.re)*mpfr::cos(z.im),  mpfr::cosh(z.re)*mpfr::sin(z.im)}; }
static CMp ccoshmp(CMp const &z) { return {mpfr::cosh(z.re)*mpfr::cos(z.im),  mpfr::sinh(z.re)*mpfr::sin(z.im)}; }
static CMp ctanmp (CMp const &z) { return cdivmp(csinmp(z),  ccosmp(z)); }
static CMp ctanhmp(CMp const &z) { return cdivmp(csinhmp(z), ccoshmp(z)); }

// casin(z) = -i · log(iz + sqrt(1 - z²))
static CMp casinmp(CMp const &z) {
  CMp zz = cmulmp(z, z);
  CMp arg = caddmp({-z.im, z.re}, csqrtmp({mp_t(1) - zz.re, -zz.im}));
  CMp l = clogmp(arg);
  return {l.im, -l.re};
}
// cacos(z) = π/2 − casin(z)
static CMp cacosmp(CMp const &z) {
  static const mp_t half_pi = mpfr::const_pi(multifloats_test::kMpfrPrec) / 2;
  CMp a = casinmp(z);
  return {half_pi - a.re, -a.im};
}
// catan(z) = (-i/2) · log((1 + iz)/(1 - iz))    [C99 Annex G principal]
static CMp catanmp(CMp const &z) {
  CMp num = {mp_t(1) - z.im,  z.re};
  CMp den = {mp_t(1) + z.im, -z.re};
  CMp l = clogmp(cdivmp(num, den));
  return {l.im / 2, -l.re / 2};
}
// casinh(z) = log(z + sqrt(z² + 1))
static CMp casinhmp(CMp const &z) {
  CMp zz = cmulmp(z, z);
  return clogmp(caddmp(z, csqrtmp({zz.re + 1, zz.im})));
}
// cacosh(z) = log(z + sqrt(z-1)·sqrt(z+1))   [C99 Annex G principal]
static CMp cacoshmp(CMp const &z) {
  CMp s = cmulmp(csqrtmp({z.re - 1, z.im}), csqrtmp({z.re + 1, z.im}));
  return clogmp(caddmp(z, s));
}
// catanh(z) = ½ · log((1 + z)/(1 - z))
static CMp catanhmp(CMp const &z) {
  CMp l = clogmp(cdivmp({mp_t(1) + z.re,  z.im}, {mp_t(1) - z.re, -z.im}));
  return {l.re / 2, l.im / 2};
}

// On a true 200-bit reference, the DD precision floor is |hi| ≥ 2^-969:
// below that, lo can't land in the normal-double range and the pair loses
// its full 106-bit error term.
static bool mp_below_subnormal_floor(mp_t const &v) {
  static const mp_t floor_ = mpfr::pow(mp_t(2), mp_t(-969));
  return mpfr::abs(v) < floor_;
}

#else  // !USE_MPFR

#define MP_PARAM(...)
#define MP_ARG(...)
#define MP_STMT(...)

#endif  // USE_MPFR

// =============================================================================
// Per-op statistics table (keyed by op name).
//
// Without USE_MPFR: one column — rel_err(DD vs qp).
// With USE_MPFR: two columns — rel_err(qp vs mp) and rel_err(DD vs mp), both
// against the 200-bit mpreal reference. Pass/fail still gates on DD vs qp
// regardless; the extra column is informational.
// =============================================================================

struct StatEntry {
  char name[24] = {};
#ifdef USE_MPFR
  double max_q  = 0.0; double sum_q  = 0.0;
  double max_dd = 0.0; double sum_dd = 0.0;
#else
  double max_rel = 0.0;
  double sum_rel = 0.0;
#endif
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

#ifdef USE_MPFR
static void update_stat(char const *op, double rq, double rdd) {
  if (!std::isfinite(rq) || !std::isfinite(rdd)) return;
  StatEntry *s = find_or_create_stat(op);
  if (!s) return;
  if (rq  > s->max_q)  s->max_q  = rq;
  if (rdd > s->max_dd) s->max_dd = rdd;
  s->sum_q  += rq;
  s->sum_dd += rdd;
  ++s->count;
}
#else
static void update_stat(char const *op, double rel) {
  if (!std::isfinite(rel)) return;
  StatEntry *s = find_or_create_stat(op);
  if (!s) return;
  if (rel > s->max_rel) s->max_rel = rel;
  s->sum_rel += rel;
  ++s->count;
}
#endif

// =============================================================================
// Classification: which ops are expected to produce full DD precision?
// =============================================================================

// Full-DD tier: kernels that should deliver ~106 bits of precision end-to-end.
// exp/log/trig/hyperbolic/pow/tgamma/lgamma hit full DD via the native
// polynomial kernels. fmod is excluded because its large-quotient subtraction
// loses precision through catastrophic cancellation once |x/y| grows.
//
// Complex ops emit "<op>_re" / "<op>_im" stat rows (one per component);
// matching both here keeps them on the full-DD tier.
// Bessel rows use "bj0"/"bj1"/"bjn"/"by0"/"by1"/"byn" labels (matches
// ops.py bench-key naming so the two surfaces line up).
// {asin,acos,atan}pi divide by π: DD and qp both divide by their own π
// (agree to ~1e-31), no amplification. Full-DD tolerance is fine.
// {sin,cos,tan}pi instead *multiply* by π inside the argument, so the DD
// vs qp π-representation mismatch gets amplified near zeros of the
// function — handled by is_pi_trig below.
static bool is_full_dd(char const *op) {
  static char const *kList[] = {
      "add", "sub", "mul", "div", "sqrt", "abs", "neg",
      "add_fd", "mul_df", "fmin", "fmax", "copysign", "fdim",
      "hypot", "trunc", "round", "scalbn", "min3", "max3", "fmadd",
      "floor", "ceil", "nearbyint", "rint",
      "lround", "llround", "lrint", "llrint",
      "logb", "ilogb", "ldexp", "scalbln",
      "frexp.frac", "frexp.exp", "modf.frac", "modf.int",
      "exp", "exp2", "expm1", "log", "log2", "log10", "log1p",
      "pow", "pow_md", "pow_dm", "pow_int",
      "sin", "cos", "tan", "sincos_s", "sincos_c",
      "sinh", "cosh", "tanh", "sinhcosh_s", "sinhcosh_c",
      "asinh", "acosh", "atanh",
      "asin", "acos", "atan", "atan2", "atan2pi",
      "asinpi", "acospi", "atanpi",
      "erf", "erfc", "erfcx",
      "gamma", "lngamma",
      "bj0", "bj1", "bjn", "by0", "by1", "byn", "yn_range",
      "cadd_re", "cadd_im", "csub_re", "csub_im",
      "cmul_re", "cmul_im", "cdiv_re", "cdiv_im",
      "cabs", "carg", "cconjg_re", "cconjg_im",
      "cproj_re", "cproj_im",
      "csqrt_re", "csqrt_im",
      "cexp_re", "cexp_im",
      // cexpm1 lives in is_reduced_dd — qp range-reduces im(z) against
      // its own π approximation and loses precision when im(z) is near nπ.
      "clog_re", "clog_im",
      "clog2_re", "clog2_im", "clog10_re", "clog10_im",
      "clog1p_re", "clog1p_im",
      "cpow_re", "cpow_im",
      "csin_re", "csin_im", "ccos_re", "ccos_im",
      "ctan_re", "ctan_im",
      "csinh_re", "csinh_im", "ccosh_re", "ccosh_im",
      "ctanh_re", "ctanh_im",
      "casinh_re", "casinh_im",
      // casin / cacos / catan / cacosh / catanh live in is_reduced_dd —
      // libquadmath's implementations lose precision on branch cuts.
      // csinpi / ccospi also live in is_reduced_dd — π-representation
      // mismatch amplified near zeros, same as scalar sinpi/cospi.
      nullptr};
  for (int i = 0; kList[i]; ++i) {
    if (std::strcmp(kList[i], op) == 0) {
      return true;
    }
  }
  return false;
}

// Reduced-DD tier: ops where the qp *oracle itself* cannot deliver
// full-DD (~1e-26) precision, so a stricter tolerance would measure the
// reference's precision floor rather than the DD kernel's. Two families:
//
//   1. Forward π-scaled trig (sinpi/cospi/tanpi). DD uses a 106-bit π, qp
//      uses 113-bit M_PIq — they differ by ~1e-33 relative. sin/cos/tan
//      amplify this by |cot(π·x)·π·x| near their zeros, pushing rel-err
//      past 1e-26 with no kernel bug.
//
//   2. Inverse complex transcendentals (casin/cacos/catan/cacosh/catanh).
//      libquadmath's c*q implementations evaluate the textbook
//      casin(z) = -i·log(iz + sqrt(1-z²)) family verbatim, which loses
//      precision on / near the branch cuts. The DD kernels use Kahan-
//      style formulations internally and are actually *more* accurate
//      than the qp reference at these edges. USE_MPFR mode is the right
//      tool to separate kernel error from reference floor; the default
//      mode just needs to not flag the reference floor as a regression.
static bool is_reduced_dd(char const *op) {
  static char const *kList[] = {
      "sinpi", "cospi", "tanpi",
      "casin_re", "casin_im",
      "cacos_re", "cacos_im",
      "catan_re", "catan_im",
      "cacosh_re", "cacosh_im",
      "catanh_re", "catanh_im",
      "csinpi_re", "csinpi_im",
      "ccospi_re", "ccospi_im",
      "cexpm1_re", "cexpm1_im",
      nullptr};
  for (int i = 0; kList[i]; ++i) {
    if (std::strcmp(kList[i], op) == 0) {
      return true;
    }
  }
  return false;
}

// Compound tier: no C++ operations currently fall into this tier.
// Kept as a placeholder for future chained-evaluation functions.
static bool is_compound(char const *) {
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

// CHK is the single check-site macro. Absorbs MP_ARG internally so the
// mpreal oracle expression never appears at a call site — keeping callers
// identical across USE_MPFR on/off.
//
// The mp slot is variadic (`...` / `__VA_ARGS__`), not a named positional
// arg: a call site like `CHK("cexpm1", ..., mag, (q_t)0, CMp{a, b})`
// has a comma inside the `CMp{...}` brace list that would otherwise get
// consumed by CHK itself (splitting the expression into two macro args).
// __VA_ARGS__ captures everything after the first five named slots as one
// unit, mirroring MP_ARG's own variadic contract.
//
// The expansion is a single `check(...)` expression-statement-ish call
// followed by the caller's `;`, so `if (cond) CHK(...);` is safe — no
// dangling-else hazard, no compound-statement wrapping needed. Keep it
// that way if you touch this.
#define CHK(op, mf_expr, q_expr, i1, i2, ...) \
    check(op, (mf_expr), (q_expr), (i1), (i2) MP_ARG(__VA_ARGS__))

// Shorthand wrappers over CHK. Each one token-pastes the op name into
// the DD / qp / mp function names following a convention, then forwards
// to CHK. The scope names they depend on (f1/q1/m1 for scalar ops,
// zd1/zq1/zm1/mag for complex ops) are set up once per main-loop
// iteration — relocating a call site out of that scope means reverting
// to raw CHK.
//
// Ops that don't follow the naming convention (abs, bjn, fmadd, sincos,
// π-scaled trig, erfcx, cconjg/cproj/carg/cabs/cpow/cexpm1/clog1p/clog2/
// clog10/csinpi/ccospi, mixed-mode add_fd/mul_df/pow_md/pow_dm/pow_int,
// min3/max3, plus arithmetic operators) keep using raw CHK.

// Scalar unary: mf::NAME(f1) / NAMEq(q1) / mpfr::NAME(m1).
#define CHK1(NAME) \
    CHK(#NAME, mf::NAME(f1), NAME##q(q1), q1, (q_t)0, mpfr::NAME(m1))

// Scalar binary: mf::NAME(f1, f2) / NAMEq(q1, q2) / mpfr::NAME(m1, m2).
#define CHK2(NAME) \
    CHK(#NAME, mf::NAME(f1, f2), NAME##q(q1, q2), q1, q2, mpfr::NAME(m1, m2))

// Complex unary: ::cSTEMdd(zd1) / cSTEMq(zq1) / cSTEMmp(zm1). Label is
// "cSTEM" (e.g. "csin", "csqrt") — matches libquadmath's cNAMEq naming
// with the q stripped.
#define CHK_C1(STEM) \
    CHK("c" #STEM, ::c##STEM##dd(zd1), c##STEM##q(zq1), \
        mag, (q_t)0, c##STEM##mp(zm1))

// Complex binary: ::cSTEMdd(zd1, zd2) / cSTEMq(zq1, zq2) / cSTEMmp(zm1, zm2).
// The qp leg reads as a function call thanks to the caddq/csubq/cmulq/cdivq
// inline wrappers above — this parallelizes cleanly with CHK_C1.
#define CHK_C2(STEM) \
    CHK("c" #STEM, ::c##STEM##dd(zd1, zd2), c##STEM##q(zq1, zq2), \
        mag, (q_t)0, c##STEM##mp(zm1, zm2))

// Gated-call variants: `CHK1_IF(sin, cond)` expands to the equivalent of
// `if (cond) CHK1(sin);` but wrapped in `do { } while (0)` so a stray
// dangling `else` at an outer call site can't bind to the inner `if`.
//
// CHK_IF takes COND as its *first* argument (not last) because CHK's tail
// slot is variadic for CMp{...} brace-comma safety; everything after COND
// is forwarded verbatim to CHK.
#define CHK_IF(COND, ...)      do { if (COND) CHK(__VA_ARGS__); }      while (0)
#define CHK1_IF(NAME, COND)    do { if (COND) CHK1(NAME); }    while (0)
#define CHK2_IF(NAME, COND)    do { if (COND) CHK2(NAME); }    while (0)
#define CHK_C1_IF(STEM, COND)  do { if (COND) CHK_C1(STEM); }  while (0)
#define CHK_C2_IF(STEM, COND)  do { if (COND) CHK_C2(STEM); }  while (0)

static void check(char const *op, mf::float64x2 const &got, q_t expected,
                  q_t i1, q_t i2 MP_PARAM(mp_t const &expected_mp)) {
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
    else if (is_reduced_dd(op)) tol = 1e-22;
    else if (is_compound(op)) tol = 1e-10;
    else tol = 1e-15;
  } else {
    rel_err = (double)(diff / input_mag);
    if (is_full_dd(op)) tol = 1e-28;
    else if (is_reduced_dd(op)) tol = 1e-22;
    else if (is_compound(op)) tol = 1e-10;
    else tol = 1e-15;
  }

  // Skip the subnormal range entirely — both stats and pass/fail. The DD
  // lo limb can't hold a full 53-bit error term below ~2^-969, so any
  // "error" reported there is measuring the format's cliff, not the
  // kernel. Matches the pre-merge fuzz.cc behavior; without this guard
  // every tiny-input division becomes a spurious pass/fail failure.
  if (subnormal_range(i1) || subnormal_range(i2) || subnormal_range(expected)) {
    return;
  }

#ifdef USE_MPFR
  // Extra USE_MPFR gate: record mp stats only when the mp oracle itself
  // is finite and above the same representational floor (expressed on the
  // 200-bit reference). Pass/fail below is unaffected and still runs on
  // rel_err(DD vs qp).
  if (mp_isfinite(expected_mp) && !mp_isnan(expected_mp) &&
      !mp_below_subnormal_floor(expected_mp)) {
    // Normalize absolute error by max(|expected|, input_mag) so ops whose
    // true result is intrinsically small relative to the inputs (fmod near
    // a zero, subtraction near cancellation) don't report astronomical
    // relative error from dividing by a tiny denominator. For amplifying
    // ops (mul, exp, …) |expected| ≥ input_mag so the max collapses to
    // |expected| and the reading is unchanged.
    mp_t denom_mp = mpfr::abs(expected_mp);
    mp_t input_mag_mp = to_mp(input_mag);
    if (denom_mp < input_mag_mp) denom_mp = input_mag_mp;
    double rq  = (mpfr::abs(to_mp(expected) - expected_mp) / denom_mp).toDouble();
    double rdd = (mpfr::abs(to_mp(got)      - expected_mp) / denom_mp).toDouble();
    update_stat(op, rq, rdd);
  }
#else
  update_stat(op, rel_err);
#endif

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

  // Narrow random: sign * uniform(0.5) * 10^k with k ∈ [-3, 3]. Used for
  // complex inputs — keeps re*re - im*im in c_mul / the cdivq division
  // well away from overflow and catastrophic-cancellation regimes so the
  // qp oracle stays a clean reference.
  q_t narrow(double r, double rexp) {
    int k = (int)(rexp * 6.0) - 3;
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
#ifdef USE_MPFR
  std::printf("Per-operation 3-way precision report (reference = mpreal @ 200 bits):\n");
  std::printf("  %-16s %10s %14s %14s %14s %14s\n",
              "op", "n", "max_q", "mean_q", "max_dd", "mean_dd");
  for (int i = 0; i < g_nstats; ++i) {
    StatEntry const &s = g_stats[i];
    if (s.count == 0) {
      std::printf("  %-16s %10ld     (no data)\n", s.name, s.count);
      continue;
    }
    double invn = 1.0 / (double)s.count;
    std::printf("  %-16s %10ld  %14.3e %14.3e %14.3e %14.3e\n",
                s.name, s.count,
                s.max_q,  s.sum_q  * invn,
                s.max_dd, s.sum_dd * invn);
  }
  std::printf("\n"
              "  Columns: q  = rel_err(libquadmath vs mpreal),\n"
              "           dd = rel_err(multifloats DD vs mpreal).\n"
              "  q ≈ 1e-33 is the float128 mantissa floor; dd above q means\n"
              "  the DD kernel — not the float128 reference — is the loss.\n");
#else
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
#endif
}

// =============================================================================
// Complex-op checker. Splits a DD complex result into its real and imag
// DD components and compares each against the qp __complex128 reference.
// Emits two stat rows per op ("<op>_re", "<op>_im") so the classifier in
// is_full_dd picks them up and so downstream tooling can key on component.
// =============================================================================

// Complex overload of check(). Splits the complex result into _re / _im
// components and dispatches each to the scalar overload. i1 carries the
// complex-input magnitude; i2 is unused in practice and set to (q_t)0 at
// call sites, but kept in the signature for uniformity with the scalar
// overload so a single CHK macro drives both.
static void check(char const *op, complex64x2_t const &got,
                  __complex128 expected, q_t i1, q_t i2
                  MP_PARAM(CMp const &expected_mp)) {
  mf::float64x2 got_re = mf::float64x2(got.re);
  mf::float64x2 got_im = mf::float64x2(got.im);
  q_t exp_re = crealq(expected);
  q_t exp_im = cimagq(expected);

  char key[32];
  std::snprintf(key, sizeof(key), "%s_re", op);
  CHK(key, got_re, exp_re, i1, i2, expected_mp.re);
  std::snprintf(key, sizeof(key), "%s_im", op);
  CHK(key, got_im, exp_im, i1, i2, expected_mp.im);
}

// Convert an mf::float64x2 pair to a C-ABI complex64x2_t for dispatching
// through the *dd kernels.
static inline complex64x2_t to_cdd(mf::float64x2 const &re, mf::float64x2 const &im) {
  complex64x2_t z;
  z.re = static_cast<float64x2_t>(re);
  z.im = static_cast<float64x2_t>(im);
  return z;
}

// Build a __complex128 from two q_t components.
static inline __complex128 to_cq(q_t re, q_t im) {
  __complex128 z;
  __real__ z = re;
  __imag__ z = im;
  return z;
}

static void check_comp(mf::float64x2 const &f1, mf::float64x2 const &f2, q_t q1, q_t q2) {
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
#ifdef USE_MPFR
  // 100× smaller iteration default — mpreal ops are ~100× slower than qp.
  long iterations = 10000;
#else
  long iterations = 1000000;
#endif
  uint64_t seed = 42ULL;
  if (argc > 1) {
    iterations = std::atol(argv[1]);
  }
  if (argc > 2) {
    seed = std::strtoull(argv[2], nullptr, 0);
  }

  Rng rng(seed);
#ifdef USE_MPFR
  mpfr::mpreal::set_default_prec(multifloats_test::kMpfrPrec);
  std::printf("[multifloats_fuzz] iterations=%ld seed=0x%llx prec=%d bits (USE_MPFR)\n",
              iterations, (unsigned long long)seed, (int)multifloats_test::kMpfrPrec);
#else
  std::printf("[multifloats_fuzz] iterations=%ld seed=0x%llx\n", iterations,
              (unsigned long long)seed);
#endif

  for (long i = 1; i <= iterations; ++i) {
    q_t q1, q2;
    rng.generate_pair(q1, q2);
    // Round-trip through DD so every qp input is already representable
    // in float64x2 exactly. Without this, qp carries 7 bits of precision
    // past what DD's 106-bit pair can hold, and those bits surface as
    // spurious DD-input error when DD and qp are measured against a
    // shared mpreal oracle (`sub(q1, q1·(1+1e-15))` would otherwise read
    // ~1e-17 qp rel-err instead of the actual 0.5 qp ulp). After the
    // clamp, to_mp(q1) == to_mp(from_q(q1)), so both systems see the
    // same math and the oracle can be sourced from either side.
    q1 = to_q(from_q(q1));
    q2 = to_q(from_q(q2));
    mf::float64x2 f1 = from_q(q1);
    mf::float64x2 f2 = from_q(q2);
    MP_STMT(mp_t m1 = to_mp(f1));
    MP_STMT(mp_t m2 = to_mp(f2));
    double d1 = (double)q1;
    double d2 = (double)q2;

    // Hot loop: arithmetic + sqrt + comparisons + unary basics.
    // add/sub/neg/abs are non-finite safe. mul/div/etc use two_prod, whose
    // error limb becomes NaN when a hi-limb is ±inf even if the IEEE product
    // is well-defined; gate the multiplicative path on finite inputs.
    bool const both_finite = q_isfinite(q1) && q_isfinite(q2);
    q_t aq1 = q1 < 0 ? -q1 : q1;

    CHK("add", f1 + f2, q1 + q2, q1, q2, m1 + m2);
    CHK("sub", f1 - f2, q1 - q2, q1, q2, m1 - m2);
    CHK_IF(both_finite, "mul", f1 * f2, q1 * q2, q1, q2, m1 * m2);
    CHK_IF(both_finite && q2 != (q_t)0, "div", f1 / f2, q1 / q2, q1, q2, m1 / m2);
    CHK1_IF(sqrt, q_isfinite(q1) && q1 >= (q_t)0);
    CHK("abs", mf::abs(f1), q1 < 0 ? -q1 : q1, q1, (q_t)0, mpfr::abs(m1));
    CHK("neg", -f1, -q1, q1, (q_t)0, -m1);
    CHK1_IF(trunc, q_isfinite(q1) && aq1 < (q_t)1e15q);
    CHK1_IF(round, q_isfinite(q1) && aq1 < (q_t)1e15q);
    // floor/ceil share trunc/round's magnitude gate: beyond 2^52 ≈ 4.5e15
    // every DD input is already integer-valued and the op is identity, so
    // the oracle comparison is uninformative. mpreal wraps both.
    CHK1_IF(floor, q_isfinite(q1) && aq1 < (q_t)1e15q);
    CHK1_IF(ceil,  q_isfinite(q1) && aq1 < (q_t)1e15q);
    // nearbyint / rint use the default rounding mode (round-half-to-even).
    // mpreal does not wrap nearbyint (no mpfr_nearbyint in MPFR), so pass
    // the qp result through to_mp as the MP oracle; pass/fail still runs
    // on DD-vs-qp which is what we care about. rint uses the same
    // passthrough for uniformity.
    CHK_IF(q_isfinite(q1) && aq1 < (q_t)1e15q, "nearbyint",
           mf::nearbyint(f1), nearbyintq(q1), q1, (q_t)0,
           to_mp(nearbyintq(q1)));
    CHK_IF(q_isfinite(q1) && aq1 < (q_t)1e15q, "rint",
           mf::rint(f1), rintq(q1), q1, (q_t)0, to_mp(rintq(q1)));

    // Integer-rounding ops (return long / long long). Gate at 2^50 ≈ 1e15
    // so the cast (long) → double in the adapter is exact (long has up to
    // 63 mantissa bits, double only 53). The same gate keeps the DD inputs
    // small enough that lround's half-integer correction never overflows.
    if (q_isfinite(q1) && aq1 < (q_t)1e15q) {
      auto dd_of_llong = [](long long v) {
        mf::float64x2 r;
        r._limbs[0] = (double)v;
        r._limbs[1] = 0.0;
        return r;
      };
      CHK("lround",  dd_of_llong(mf::lround(f1)),
          (q_t)lroundq(q1),  q1, (q_t)0, mp_t((long long)lroundq(q1)));
      CHK("llround", dd_of_llong(mf::llround(f1)),
          (q_t)llroundq(q1), q1, (q_t)0, mp_t((long long)llroundq(q1)));
      CHK("lrint",   dd_of_llong(mf::lrint(f1)),
          (q_t)lrintq(q1),   q1, (q_t)0, mp_t((long long)lrintq(q1)));
      CHK("llrint",  dd_of_llong(mf::llrint(f1)),
          (q_t)llrintq(q1),  q1, (q_t)0, mp_t((long long)llrintq(q1)));
    }

    // logb / ilogb: exponent extraction. Gate on finite nonzero — logb(0)
    // is -inf and ilogb(0) is FP_ILOGB0 (implementation-defined), both
    // oracle-matching on the qp side but not a useful precision check.
    // Use qp-passthrough for the MP oracle; mpreal's logb binding may or
    // may not exist depending on the mpreal version, and pass/fail runs
    // on rel_err(DD vs qp) regardless.
    CHK_IF(q_isfinite(q1) && q1 != (q_t)0, "logb",
           mf::logb(f1), logbq(q1), q1, (q_t)0, to_mp(logbq(q1)));
    if (q_isfinite(q1) && q1 != (q_t)0) {
      mf::float64x2 ilogb_dd;
      ilogb_dd._limbs[0] = (double)mf::ilogb(f1);
      ilogb_dd._limbs[1] = 0.0;
      CHK("ilogb", ilogb_dd, (q_t)ilogbq(q1), q1, (q_t)0,
          mp_t((long long)ilogbq(q1)));
    }

    check_comp(f1, f2, q1, q2);

    // Periodic (every 10): mixed-mode + binary. Mixed-mode and binary ops all
    // route through the multiplicative or hypot kernels — share mul's
    // non-finite limitation.
    if (i % 10 == 0 && both_finite) {
      CHK("add_fd", f1 + mf::float64x2(d2), q1 + (q_t)d2,
          q1, (q_t)d2, m1 + to_mp(d2));
      CHK("mul_df", mf::float64x2(d1) * f2, (q_t)d1 * q2,
          (q_t)d1, q2, to_mp(d1) * m2);
      CHK2(fmin);
      CHK2(fmax);
      CHK2(copysign);
      CHK2(fdim);
      CHK2(hypot);

      q_t aq2 = q2 < 0 ? -q2 : q2;
      CHK2_IF(fmod, q2 != (q_t)0 && aq1 < (q_t)1e20q && aq2 > (q_t)1e-20q);

      // fmadd(a, b, c) = a·b + c. Use d1 as the third input (DD-
      // representable) so the qp reference is exactly fmaq(q1, q2, d1).
      CHK("fmadd",
          mf::float64x2(::fmadd(
              static_cast<float64x2_t>(f1),
              static_cast<float64x2_t>(f2),
              static_cast<float64x2_t>(mf::float64x2(d1)))),
          fmaq(q1, q2, (q_t)d1), q1, q2, mpfr::fma(m1, m2, to_mp(d1)));
    }

    // Periodic (every 100): transcendentals.
    if (i % 100 == 0) {
      CHK1_IF(exp,   q_isfinite(q1) && q_isfinite(expq(q1)));
      CHK1_IF(log,   q_isfinite(q1) && q1 > (q_t)0);
      CHK1_IF(log10, q_isfinite(q1) && q1 > (q_t)0);
      CHK1_IF(expm1, q_isfinite(q1) && q_isfinite(expm1q(q1)));
      CHK1_IF(log1p, q_isfinite(q1) && q1 > (q_t)-1);
      CHK1_IF(exp2,  q_isfinite(q1) && q_isfinite(exp2q(q1)));
      CHK1_IF(log2,  q_isfinite(q1) && q1 > (q_t)0);

      // Trig: keep magnitudes moderate.
      if (q_isfinite(q1) && aq1 < (q_t)1e6q) {
        CHK1(sin);
        CHK1(cos);
        CHK1_IF(tan, fabsq(cosq(q1)) > (q_t)1e-12q);

        // Fused sincos: one range-reduction feeds both outputs.
        float64x2_t sc_s, sc_c;
        ::sincosdd(static_cast<float64x2_t>(f1), &sc_s, &sc_c);
        CHK("sincos_s", mf::float64x2(sc_s), sinq(q1),
            q1, (q_t)0, mpfr::sin(m1));
        CHK("sincos_c", mf::float64x2(sc_c), cosq(q1),
            q1, (q_t)0, mpfr::cos(m1));
      }
      CHK1_IF(asin, q_isfinite(q1) && aq1 <= (q_t)1);
      CHK1_IF(acos, q_isfinite(q1) && aq1 <= (q_t)1);
      CHK1_IF(atan, q_isfinite(q1));

      if (q_isfinite(q1) && aq1 < (q_t)700) {
        CHK1(sinh);
        CHK1(cosh);
        CHK1(tanh);

        // Fused sinhcosh.
        float64x2_t hc_s, hc_c;
        ::sinhcoshdd(static_cast<float64x2_t>(f1), &hc_s, &hc_c);
        CHK("sinhcosh_s", mf::float64x2(hc_s), sinhq(q1),
            q1, (q_t)0, mpfr::sinh(m1));
        CHK("sinhcosh_c", mf::float64x2(hc_c), coshq(q1),
            q1, (q_t)0, mpfr::cosh(m1));
      }
      CHK1_IF(asinh, q_isfinite(q1));
      CHK1_IF(acosh, q_isfinite(q1) && q1 >= (q_t)1);
      CHK1_IF(atanh, q_isfinite(q1) && aq1 <  (q_t)1);

      CHK1_IF(erf,  q_isfinite(q1) && aq1 < (q_t)100);
      CHK1_IF(erfc, q_isfinite(q1) && aq1 < (q_t)100);
      // erfcx(x) = exp(x²) * erfc(x). The DD kernel is precision-preserving
      // for the tail-tail cancellation; the qp oracle is the naive product.
      // Gate |x| < 26 so exp(x²) stays well inside qp's exponent range
      // and the product stays representable.
      CHK_IF(q_isfinite(q1) && aq1 < (q_t)26, "erfcx",
             mf::erfcx(f1), expq(q1 * q1) * erfcq(q1),
             q1, (q_t)0,
             mpfr::exp(m1 * m1) * mpfr::erfc(m1));
      if (q_isfinite(q1) && q1 > (q_t)0 && q1 < (q_t)100) {
        CHK1(gamma);
        CHK1(lngamma);
      }

      CHK2(atan2);

      // atan2pi(y, x) = atan2(y, x) / π. Division by π has low
      // amplification, so full-DD tolerance holds here (unlike forward
      // sin/cos/tanpi — see is_pi_trig).
      CHK("atan2pi",
          mf::float64x2(::atan2pidd(
              static_cast<float64x2_t>(f1), static_cast<float64x2_t>(f2))),
          atan2q(q1, q2) / (q_t)M_PIq,
          q1, q2,
          mpfr::atan2(m1, m2) / mpfr::const_pi(multifloats_test::kMpfrPrec));

      // Power: positive base, modest exponent.
      q_t aq2 = q2 < 0 ? -q2 : q2;
      if (q_isfinite(q1) && q1 > (q_t)1e-3q && q1 < (q_t)1e3q &&
          q_isfinite(q2) && aq2 < (q_t)30) {
        CHK2(pow);
        CHK("pow_md", mf::pow(f1, mf::float64x2(d2)), powq(q1, (q_t)d2),
            q1, (q_t)d2, mpfr::pow(m1, to_mp(d2)));
        CHK("pow_dm", mf::pow(mf::float64x2(d1), f2), powq((q_t)d1, q2),
            (q_t)d1, q2, mpfr::pow(to_mp(d1), m2));
      }
      CHK_IF(q_isfinite(q1) && aq1 < (q_t)1e10q, "pow_int",
             mf::pow(f1, mf::float64x2(3.0)), powq(q1, (q_t)3),
             q1, (q_t)3, mpfr::pow(m1, to_mp(3.0)));

      // scalbn(x, 5) = x · 2^5; the ·32 is exact at any precision.
      CHK_IF(q_isfinite(q1), "scalbn", mf::scalbn(f1, 5), scalbnq(q1, 5),
             q1, (q_t)0, m1 * mp_t(32));
      // ldexp / scalbln: POSIX aliases of scalbn under FLT_RADIX == 2
      // (IEEE 754 guarantees it). Same ·32 oracle. ldexp takes int,
      // scalbln takes long — distinct header-only wrappers in multifloats.h.
      CHK_IF(q_isfinite(q1), "ldexp",   mf::ldexp(f1, 5),    ldexpq(q1, 5),
             q1, (q_t)0, m1 * mp_t(32));
      CHK_IF(q_isfinite(q1), "scalbln", mf::scalbln(f1, 5L), scalblnq(q1, 5L),
             q1, (q_t)0, m1 * mp_t(32));

      // frexp(x, &e): DD returns fraction in [0.5, 1) for nonzero x and
      // exponent e such that x == fraction · 2^e. Emit two stat rows —
      // one for the fraction, one for the integer exponent — so per-
      // component regressions show up separately.
      if (q_isfinite(q1) && q1 != (q_t)0) {
        int e_mf = 0, e_q = 0;
        mf::float64x2 fr_mf = mf::frexp(f1, &e_mf);
        q_t fr_q = frexpq(q1, &e_q);
        CHK("frexp.frac", fr_mf, fr_q, q1, (q_t)0, to_mp(fr_q));
        mf::float64x2 eint_mf;
        eint_mf._limbs[0] = (double)e_mf;
        eint_mf._limbs[1] = 0.0;
        CHK("frexp.exp", eint_mf, (q_t)e_q, q1, (q_t)0,
            mp_t((long long)e_q));
      }

      // modf(x, &iptr): splits x into fractional and integer parts, both
      // preserving sign. Gate on |x| < 1e15 for the same reason as trunc:
      // beyond 2^52 the fractional part is zero and the test is identity.
      if (q_isfinite(q1) && aq1 < (q_t)1e15q) {
        mf::float64x2 ipart_mf;
        q_t ipart_q = 0;
        mf::float64x2 frac_mf = mf::modf(f1, &ipart_mf);
        q_t frac_q = modfq(q1, &ipart_q);
        CHK("modf.frac", frac_mf, frac_q, q1, (q_t)0, to_mp(frac_q));
        CHK("modf.int",  ipart_mf, ipart_q, q1, (q_t)0, to_mp(ipart_q));
      }

      // 3-argument min/max via nested fmin/fmax.
      if (both_finite) {
        q_t q3 = (q1 + q2) * (q_t)0.5q;
        // q1+q2 of two DD-precision values can spill into bit 107; clamp
        // back to DD so q3, f3, m3 stay fully round-trip consistent.
        q3 = to_q(from_q(q3));
        mf::float64x2 f3 = from_q(q3);
        q_t min12 = q1 < q2 ? q1 : q2;
        q_t q_min3 = min12 < q3 ? min12 : q3;
        q_t max12 = q1 < q2 ? q2 : q1;
        q_t q_max3 = max12 < q3 ? q3 : max12;
        MP_STMT(mp_t m3 = to_mp(f3));
        MP_STMT(mp_t mp_min12 = q1 < q2 ? m1 : m2);
        MP_STMT(mp_t mp_min3  = min12 < q3 ? mp_min12 : m3);
        MP_STMT(mp_t mp_max12 = q1 < q2 ? m2 : m1);
        MP_STMT(mp_t mp_max3  = max12 < q3 ? m3 : mp_max12);
        CHK("min3", mf::fmin(mf::fmin(f1, f2), f3), q_min3, q1, q2, mp_min3);
        CHK("max3", mf::fmax(mf::fmax(f1, f2), f3), q_max3, q1, q2, mp_max3);
      }

      // Bessel of the first (j) and second (y) kind. Labels match
      // ops.py's bench_key convention ("bj0" / "byn" / ...). y* is only
      // finite on positive arguments. Cap |x| at 200 so the DD recurrence
      // doesn't hit the qp reference's large-argument oscillation regime
      // where the two implementations can disagree by more than the kernel
      // error alone.
      if (q_isfinite(q1) && aq1 < (q_t)200) {
        CHK("bj0", mf::float64x2(::j0dd(static_cast<float64x2_t>(f1))),
            j0q(q1), q1, (q_t)0, mpfr::besselj0(m1));
        CHK("bj1", mf::float64x2(::j1dd(static_cast<float64x2_t>(f1))),
            j1q(q1), q1, (q_t)0, mpfr::besselj1(m1));

        static constexpr int kBesselOrders[] = {2, 3, 5, 8};
        for (int n : kBesselOrders) {
          CHK("bjn",
              mf::float64x2(::jndd(n, static_cast<float64x2_t>(f1))),
              jnq(n, q1), q1, (q_t)n, mpfr::besseljn((long)n, m1));
          CHK_IF(q1 > (q_t)0, "byn",
                 mf::float64x2(::yndd(n, static_cast<float64x2_t>(f1))),
                 ynq(n, q1), q1, (q_t)n, mpfr::besselyn((long)n, m1));
        }

        if (q1 > (q_t)0) {
          CHK("by0", mf::float64x2(::y0dd(static_cast<float64x2_t>(f1))),
              y0q(q1), q1, (q_t)0, mpfr::bessely0(m1));
          CHK("by1", mf::float64x2(::y1dd(static_cast<float64x2_t>(f1))),
              y1q(q1), q1, (q_t)0, mpfr::bessely1(m1));

          // yndd_range: single forward-recurrence sweep filling out[0..5].
          // Each output is compared individually against ynq(n, x). One
          // label "yn_range" covers all six so the stat row aggregates
          // error across the entire range sweep.
          float64x2_t yn_out[6];
          ::yndd_range(0, 5, static_cast<float64x2_t>(f1), yn_out);
          for (int n = 0; n <= 5; ++n)
            CHK("yn_range",
                mf::float64x2(yn_out[n]), ynq(n, q1),
                q1, (q_t)n, mpfr::besselyn((long)n, m1));
        }
      }

      // π-scaled trig. The qp oracle composes {sin,cos,tan}q(M_PIq·x) and
      // {asin,acos,atan}q(x)/M_PIq. The forward variants amplify the
      // DD-vs-qp π mismatch near zeros of the function — see is_pi_trig
      // commentary above — so we relax tolerance for those and still gate
      // aggressively on near-zero arguments.
      if (q_isfinite(q1) && aq1 < (q_t)1e6q) {
        q_t pix = (q_t)M_PIq * q1;
        q_t sp = sinq(pix);
        q_t cp = cosq(pix);
        // sinpi fails near integer x (sin(π·n)=0); cospi fails near
        // half-integer (cos(π·(n+1/2))=0). Gate each on its own magnitude.
        CHK_IF(fabsq(sp) > (q_t)1e-10q, "sinpi",
               mf::float64x2(::sinpidd(static_cast<float64x2_t>(f1))),
               sp, q1, (q_t)0,
               mpfr::sin(mpfr::const_pi(multifloats_test::kMpfrPrec) * m1));
        CHK_IF(fabsq(cp) > (q_t)1e-10q, "cospi",
               mf::float64x2(::cospidd(static_cast<float64x2_t>(f1))),
               cp, q1, (q_t)0,
               mpfr::cos(mpfr::const_pi(multifloats_test::kMpfrPrec) * m1));
        CHK_IF(fabsq(cp) > (q_t)1e-10q && fabsq(sp) > (q_t)1e-10q, "tanpi",
               mf::float64x2(::tanpidd(static_cast<float64x2_t>(f1))),
               sp / cp, q1, (q_t)0,
               mpfr::tan(mpfr::const_pi(multifloats_test::kMpfrPrec) * m1));
      }
      if (q_isfinite(q1) && aq1 <= (q_t)1) {
        CHK("asinpi",
            mf::float64x2(::asinpidd(static_cast<float64x2_t>(f1))),
            asinq(q1) / (q_t)M_PIq, q1, (q_t)0,
            mpfr::asin(m1) / mpfr::const_pi(multifloats_test::kMpfrPrec));
        CHK("acospi",
            mf::float64x2(::acospidd(static_cast<float64x2_t>(f1))),
            acosq(q1) / (q_t)M_PIq, q1, (q_t)0,
            mpfr::acos(m1) / mpfr::const_pi(multifloats_test::kMpfrPrec));
      }
      CHK_IF(q_isfinite(q1), "atanpi",
             mf::float64x2(::atanpidd(static_cast<float64x2_t>(f1))),
             atanq(q1) / (q_t)M_PIq, q1, (q_t)0,
             mpfr::atan(m1) / mpfr::const_pi(multifloats_test::kMpfrPrec));

      // Complex DD ops. Use narrow-magnitude inputs so the oracle's
      // re*re + im*im (in cdivq) and re*re - im*im (in cmulq) don't drift
      // into overflow or catastrophic cancellation.
      {
        q_t zr1, zi1, zr2, zi2;
        zr1 = rng.narrow(rng.u(), rng.u());
        zi1 = rng.narrow(rng.u(), rng.u());
        zr2 = rng.narrow(rng.u(), rng.u());
        zi2 = rng.narrow(rng.u(), rng.u());
        if (q_isfinite(zr1) && q_isfinite(zi1) &&
            q_isfinite(zr2) && q_isfinite(zi2)) {
          mf::float64x2 fr1 = from_q(zr1), fi1 = from_q(zi1);
          mf::float64x2 fr2 = from_q(zr2), fi2 = from_q(zi2);
          q_t qr1 = to_q(fr1), qi1 = to_q(fi1);
          q_t qr2 = to_q(fr2), qi2 = to_q(fi2);
          complex64x2_t zd1 = to_cdd(fr1, fi1);
          complex64x2_t zd2 = to_cdd(fr2, fi2);
          __complex128 zq1 = to_cq(qr1, qi1);
          __complex128 zq2 = to_cq(qr2, qi2);
          MP_STMT(CMp zm1 = {to_mp(fr1), to_mp(fi1)});
          MP_STMT(CMp zm2 = {to_mp(fr2), to_mp(fi2)});

          q_t mag1 = cabsq(zq1);
          q_t mag2 = cabsq(zq2);
          q_t mag  = mag1 > mag2 ? mag1 : mag2;
          if (mag < (q_t)1e-300q) mag = (q_t)1e-300q;

          CHK_C2(add);
          CHK_C2(sub);
          CHK_C2(mul);
          CHK_C2_IF(div, mag2 > (q_t)0);
          CHK("cabs",
              mf::float64x2(::cabsdd(zd1)), cabsq(zq1),
              mag, (q_t)0, cabsmp(zm1));
          CHK("cconjg", ::conjdd(zd1), conjq(zq1), mag, (q_t)0,
              cconjgmp(zm1));

          CHK_C1(sqrt);
          CHK_C1(exp);

          // cexpm1: the DD kernel preserves precision for small |z|. The
          // qp oracle (cexpq(z) - 1) loses it there, so gate |z| > 1e-5
          // to keep the oracle meaningful. The mp oracle comes as a
          // trailing CMp{...} — the brace list's comma is safely inside
          // CHK's variadic tail.
          if (mag1 > (q_t)1e-5q) {
            __complex128 one_c; __real__ one_c = (q_t)1; __imag__ one_c = (q_t)0;
            CHK("cexpm1", ::cexpm1dd(zd1), cexpq(zq1) - one_c, mag, (q_t)0,
                CMp{cexpmp(zm1).re - 1, cexpmp(zm1).im});
          }

          // clog: guard |z| away from 0 (and away from overflow) so the
          // oracle's arg computation stays accurate.
          if (mag1 > (q_t)1e-200q && mag1 < (q_t)1e100q) {
            CHK_C1(log);

            // clog2 / clog10 / clog1p: libquadmath has no direct variants,
            // compose from clogq. Both components are divided by the real
            // scalar log(2) / log(10) respectively.
            __complex128 lg = clogq(zq1);
            q_t log_2  = logq((q_t)2);
            q_t log_10 = logq((q_t)10);
            __complex128 lg2;  __real__ lg2  = crealq(lg) / log_2;
                               __imag__ lg2  = cimagq(lg) / log_2;
            __complex128 lg10; __real__ lg10 = crealq(lg) / log_10;
                               __imag__ lg10 = cimagq(lg) / log_10;
            MP_STMT(CMp clog2_oracle  = {clogmp(zm1).re / mpfr::log(mp_t(2)),
                                         clogmp(zm1).im / mpfr::log(mp_t(2))});
            MP_STMT(CMp clog10_oracle = {clogmp(zm1).re / mpfr::log(mp_t(10)),
                                         clogmp(zm1).im / mpfr::log(mp_t(10))});
            CHK("clog2",  ::clog2dd(zd1),  lg2,  mag, (q_t)0, clog2_oracle);
            CHK("clog10", ::clog10dd(zd1), lg10, mag, (q_t)0, clog10_oracle);
          }

          // clog1p: oracle = clogq(1 + z). Gate |1+z| to avoid oracle
          // blow-up on the branch cut at z=-1.
          {
            __complex128 one_c; __real__ one_c = (q_t)1; __imag__ one_c = (q_t)0;
            __complex128 one_plus_z = one_c + zq1;
            CHK_IF(cabsq(one_plus_z) > (q_t)1e-5q,
                   "clog1p", ::clog1pdd(zd1), clogq(one_plus_z), mag, (q_t)0,
                   clogmp({zm1.re + 1, zm1.im}));
          }

          // cpow: keep base and exponent modest — pow amplifies input
          // precision, so a wide exponent would swamp the tolerance.
          CHK_IF(mag1 > (q_t)1e-2q && mag1 < (q_t)1e2q && mag2 < (q_t)10,
                 "cpow", ::cpowdd(zd1, zd2), cpowq(zq1, zq2), mag, (q_t)0,
                 cexpmp(cmulmp(zm2, clogmp(zm1))));

          // csinpi / ccospi: forward π-scaled complex trig. Same π
          // representation mismatch as scalar {sin,cos}pi — reduced_dd
          // tolerance applies (see is_reduced_dd).
          {
            __complex128 pi_c; __real__ pi_c = (q_t)M_PIq; __imag__ pi_c = (q_t)0;
            __complex128 pi_z = pi_c * zq1;
            MP_STMT(mp_t mpi = mpfr::const_pi(multifloats_test::kMpfrPrec));
            CHK("csinpi", ::csinpidd(zd1), csinq(pi_z), mag, (q_t)0,
                csinmp({mpi * zm1.re, mpi * zm1.im}));
            CHK("ccospi", ::ccospidd(zd1), ccosq(pi_z), mag, (q_t)0,
                ccosmp({mpi * zm1.re, mpi * zm1.im}));
          }

          // cproj: finite inputs → identity; only meaningful distinction
          // from a plain copy happens on infinite inputs (which the
          // narrow generator doesn't produce). Still worth testing to
          // pin the finite path.
          CHK("cproj", ::cprojdd(zd1), cprojq(zq1), mag, (q_t)0, zm1);

          // carg: complex → real scalar phase. Near ±π on the branch cut
          // (negative real axis) DD and qp agree as long as imag-limbs
          // sign matches, which they do for DD-representable inputs.
          CHK("carg",
              mf::float64x2(::cargdd(zd1)), cargq(zq1), mag, (q_t)0,
              cargmp(zm1));

          CHK_C1(sin);
          CHK_C1(cos);
          CHK_C1(sinh);
          CHK_C1(cosh);

          // c{tan,tanh}: zeros of the corresponding {cos,cosh} cause the
          // division inside the oracle to blow up. Gate on oracle magnitude.
          __complex128 ccz  = ccosq(zq1);
          __complex128 ccsh = ccoshq(zq1);
          CHK_C1_IF(tan,  cabsq(ccz)  > (q_t)1e-10q);
          CHK_C1_IF(tanh, cabsq(ccsh) > (q_t)1e-10q);

          CHK_C1(asin);
          CHK_C1(acos);
          CHK_C1(atan);
          CHK_C1(asinh);
          CHK_C1(acosh);
          CHK_C1(atanh);
        }
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
