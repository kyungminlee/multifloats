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

#include "multifloats.hh"
#include "multifloats_c.h"
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

#define MP_PARAM(...) , __VA_ARGS__
#define MP_ARG(...)   , (__VA_ARGS__)
#define MP_STMT(...)  __VA_ARGS__

// MPFR proper has no complex functions. For each DD complex op we compose
// the reference by hand from real mpreal identities that match the C99
// Annex G branch cut convention (same as libquadmath). The DD complex
// kernels in src/multifloats_math.cc follow Annex G, so the oracle and
// the kernel share branch conventions.
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
      "exp", "exp2", "expm1", "log", "log2", "log10", "log1p",
      "pow", "pow_md", "pow_dm", "pow_int",
      "sin", "cos", "tan", "sincos_s", "sincos_c",
      "sinh", "cosh", "tanh", "sinhcosh_s", "sinhcosh_c",
      "asinh", "acosh", "atanh",
      "asin", "acos", "atan", "atan2", "atan2pi",
      "asinpi", "acospi", "atanpi",
      "erf", "erfc", "erfcx",
      "tgamma", "lgamma",
      "bj0", "bj1", "bjn", "by0", "by1", "byn", "yn_range",
      "cdd_add_re", "cdd_add_im", "cdd_sub_re", "cdd_sub_im",
      "cdd_mul_re", "cdd_mul_im", "cdd_div_re", "cdd_div_im",
      "cdd_abs", "cdd_arg", "cdd_conjg_re", "cdd_conjg_im",
      "cdd_proj_re", "cdd_proj_im",
      "cdd_sqrt_re", "cdd_sqrt_im",
      "cdd_exp_re", "cdd_exp_im",
      // cdd_expm1 lives in is_reduced_dd — qp range-reduces im(z) against
      // its own π approximation and loses precision when im(z) is near nπ.
      "cdd_log_re", "cdd_log_im",
      "cdd_log2_re", "cdd_log2_im", "cdd_log10_re", "cdd_log10_im",
      "cdd_log1p_re", "cdd_log1p_im",
      "cdd_pow_re", "cdd_pow_im",
      "cdd_sin_re", "cdd_sin_im", "cdd_cos_re", "cdd_cos_im",
      "cdd_tan_re", "cdd_tan_im",
      "cdd_sinh_re", "cdd_sinh_im", "cdd_cosh_re", "cdd_cosh_im",
      "cdd_tanh_re", "cdd_tanh_im",
      "cdd_asinh_re", "cdd_asinh_im",
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
      "cdd_asin_re", "cdd_asin_im",
      "cdd_acos_re", "cdd_acos_im",
      "cdd_atan_re", "cdd_atan_im",
      "cdd_acosh_re", "cdd_acosh_im",
      "cdd_atanh_re", "cdd_atanh_im",
      "cdd_sinpi_re", "cdd_sinpi_im",
      "cdd_cospi_re", "cdd_cospi_im",
      "cdd_expm1_re", "cdd_expm1_im",
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

// CHK1 / CHK2 / CHK1_IF / CHK2_IF collapse the common guard+check pattern.
// Expand inside the main loop where the locals `q1`, `q2`, `f1`, `f2`, and
// the 1-arg `(q_t)0` tail are all in scope. The explicit `i1, i2` trailing
// args on `check()` drive the input-magnitude floor inside check() itself.
// `mp_expr` is the mpreal oracle expression; it is dropped by MP_ARG when
// USE_MPFR is off, so call sites stay identical across modes.
#define CHK1(op, mf_expr, q_expr, mp_expr) \
    check(op, (mf_expr), (q_expr), q1, (q_t)0 MP_ARG(mp_expr))
#define CHK2(op, mf_expr, q_expr, mp_expr) \
    check(op, (mf_expr), (q_expr), q1, q2 MP_ARG(mp_expr))
#define CHK1_IF(cond, op, mf_expr, q_expr, mp_expr) \
    do { if (cond) CHK1(op, mf_expr, q_expr, mp_expr); } while (0)
#define CHK2_IF(cond, op, mf_expr, q_expr, mp_expr) \
    do { if (cond) CHK2(op, mf_expr, q_expr, mp_expr); } while (0)

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
    double rq  = mp_rel_err(to_mp(expected), expected_mp);
    double rdd = mp_rel_err(to_mp(got),      expected_mp);
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

static void check_cplx(char const *op, complex64x2_t const &got,
                       __complex128 expected, q_t input_mag
                       MP_PARAM(CMp const &expected_mp)) {
  mf::float64x2 got_re = mf::detail::from_f64x2(got.re);
  mf::float64x2 got_im = mf::detail::from_f64x2(got.im);
  q_t exp_re = crealq(expected);
  q_t exp_im = cimagq(expected);

  char key[32];
  std::snprintf(key, sizeof(key), "%s_re", op);
  check(key, got_re, exp_re, input_mag, (q_t)0 MP_ARG(expected_mp.re));
  std::snprintf(key, sizeof(key), "%s_im", op);
  check(key, got_im, exp_im, input_mag, (q_t)0 MP_ARG(expected_mp.im));
}

// Convert an mf::float64x2 pair to a C-ABI complex64x2_t for dispatching
// through the *dd kernels.
static inline complex64x2_t to_cdd(mf::float64x2 const &re, mf::float64x2 const &im) {
  complex64x2_t z;
  z.re = mf::detail::to_f64x2(re);
  z.im = mf::detail::to_f64x2(im);
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

    CHK2("add", f1 + f2, q1 + q2, m1 + m2);
    CHK2("sub", f1 - f2, q1 - q2, m1 - m2);
    CHK2_IF(both_finite, "mul", f1 * f2, q1 * q2, m1 * m2);
    CHK2_IF(both_finite && q2 != (q_t)0, "div", f1 / f2, q1 / q2, m1 / m2);
    CHK1_IF(q_isfinite(q1) && q1 >= (q_t)0, "sqrt",
            mf::sqrt(f1), sqrtq(q1), mpfr::sqrt(m1));
    CHK1("abs", mf::abs(f1), q1 < 0 ? -q1 : q1, mpfr::abs(m1));
    CHK1("neg", -f1, -q1, -m1);
    CHK1_IF(q_isfinite(q1) && aq1 < (q_t)1e15q, "trunc",
            mf::trunc(f1), truncq(q1), mpfr::trunc(m1));
    CHK1_IF(q_isfinite(q1) && aq1 < (q_t)1e15q, "round",
            mf::round(f1), roundq(q1), mpfr::round(m1));

    check_comp(f1, f2, q1, q2);

    // Periodic (every 10): mixed-mode + binary. Mixed-mode and binary ops all
    // route through the multiplicative or hypot kernels — share mul's
    // non-finite limitation.
    if (i % 10 == 0 && both_finite) {
      check("add_fd", f1 + mf::float64x2(d2), q1 + (q_t)d2, q1, (q_t)d2
            MP_ARG(m1 + to_mp(d2)));
      check("mul_df", mf::float64x2(d1) * f2, (q_t)d1 * q2, (q_t)d1, q2
            MP_ARG(to_mp(d1) * m2));
      CHK2("fmin", mf::fmin(f1, f2), q1 < q2 ? q1 : q2, q1 < q2 ? m1 : m2);
      CHK2("fmax", mf::fmax(f1, f2), q1 < q2 ? q2 : q1, q1 < q2 ? m2 : m1);
      CHK2("copysign", mf::copysign(f1, f2), copysignq(q1, q2),
           (m2 >= 0 ? mpfr::abs(m1) : -mpfr::abs(m1)));
      CHK2("fdim", mf::fdim(f1, f2), fdimq(q1, q2),
           (m1 > m2 ? m1 - m2 : mp_t(0)));
      CHK2("hypot", mf::hypot(f1, f2), hypotq(q1, q2), mpfr::hypot(m1, m2));

      q_t aq2 = q2 < 0 ? -q2 : q2;
      CHK2_IF(q2 != (q_t)0 && aq1 < (q_t)1e20q && aq2 > (q_t)1e-20q,
              "fmod", mf::fmod(f1, f2), fmodq(q1, q2), mpfr::fmod(m1, m2));

      // fmadd(a, b, c) = a·b + c. Use d1 as the third input (DD-
      // representable) so the qp reference is exactly fmaq(q1, q2, d1).
      check("fmadd",
            mf::detail::from_f64x2(::fmadd(
                mf::detail::to_f64x2(f1),
                mf::detail::to_f64x2(f2),
                mf::detail::to_f64x2(mf::float64x2(d1)))),
            fmaq(q1, q2, (q_t)d1), q1, q2
            MP_ARG(mpfr::fma(m1, m2, to_mp(d1))));
    }

    // Periodic (every 100): transcendentals.
    if (i % 100 == 0) {
      if (q_isfinite(q1)) {
        q_t qe = expq(q1);
        CHK1_IF(q_isfinite(qe), "exp", mf::exp(f1), qe, mpfr::exp(m1));
      }
      if (q_isfinite(q1) && q1 > (q_t)0) {
        CHK1("log",   mf::log(f1),   logq(q1),   mpfr::log(m1));
        CHK1("log10", mf::log10(f1), log10q(q1), mpfr::log10(m1));
      }
      if (q_isfinite(q1)) {
        q_t qem = expm1q(q1);
        CHK1_IF(q_isfinite(qem), "expm1", mf::expm1(f1), qem, mpfr::expm1(m1));
      }
      CHK1_IF(q_isfinite(q1) && q1 > (q_t)-1, "log1p",
              mf::log1p(f1), log1pq(q1), mpfr::log1p(m1));

      if (q_isfinite(q1)) {
        q_t qe2 = exp2q(q1);
        CHK1_IF(q_isfinite(qe2), "exp2",
                mf::detail::from_f64x2(::exp2dd(mf::detail::to_f64x2(f1))),
                qe2, mpfr::exp2(m1));
      }
      CHK1_IF(q_isfinite(q1) && q1 > (q_t)0, "log2",
              mf::detail::from_f64x2(::log2dd(mf::detail::to_f64x2(f1))),
              log2q(q1), mpfr::log2(m1));

      // Trig: keep magnitudes moderate.
      if (q_isfinite(q1) && aq1 < (q_t)1e6q) {
        CHK1("sin", mf::sin(f1), sinq(q1), mpfr::sin(m1));
        CHK1("cos", mf::cos(f1), cosq(q1), mpfr::cos(m1));
        CHK1_IF(fabsq(cosq(q1)) > (q_t)1e-12q, "tan",
                mf::tan(f1), tanq(q1), mpfr::tan(m1));

        // Fused sincos: one range-reduction feeds both outputs.
        float64x2_t sc_s, sc_c;
        ::sincosdd(mf::detail::to_f64x2(f1), &sc_s, &sc_c);
        check("sincos_s", mf::detail::from_f64x2(sc_s), sinq(q1),
              q1, (q_t)0 MP_ARG(mpfr::sin(m1)));
        check("sincos_c", mf::detail::from_f64x2(sc_c), cosq(q1),
              q1, (q_t)0 MP_ARG(mpfr::cos(m1)));
      }
      if (q_isfinite(q1) && aq1 <= (q_t)1) {
        CHK1("asin", mf::asin(f1), asinq(q1), mpfr::asin(m1));
        CHK1("acos", mf::acos(f1), acosq(q1), mpfr::acos(m1));
      }
      CHK1_IF(q_isfinite(q1), "atan", mf::atan(f1), atanq(q1), mpfr::atan(m1));

      if (q_isfinite(q1) && aq1 < (q_t)700) {
        CHK1("sinh", mf::sinh(f1), sinhq(q1), mpfr::sinh(m1));
        CHK1("cosh", mf::cosh(f1), coshq(q1), mpfr::cosh(m1));
        CHK1("tanh", mf::tanh(f1), tanhq(q1), mpfr::tanh(m1));

        // Fused sinhcosh.
        float64x2_t hc_s, hc_c;
        ::sinhcoshdd(mf::detail::to_f64x2(f1), &hc_s, &hc_c);
        check("sinhcosh_s", mf::detail::from_f64x2(hc_s), sinhq(q1),
              q1, (q_t)0 MP_ARG(mpfr::sinh(m1)));
        check("sinhcosh_c", mf::detail::from_f64x2(hc_c), coshq(q1),
              q1, (q_t)0 MP_ARG(mpfr::cosh(m1)));
      }
      CHK1_IF(q_isfinite(q1),                  "asinh",
              mf::asinh(f1), asinhq(q1), mpfr::asinh(m1));
      CHK1_IF(q_isfinite(q1) && q1 >= (q_t)1,  "acosh",
              mf::acosh(f1), acoshq(q1), mpfr::acosh(m1));
      CHK1_IF(q_isfinite(q1) && aq1 <  (q_t)1, "atanh",
              mf::atanh(f1), atanhq(q1), mpfr::atanh(m1));

      if (q_isfinite(q1) && aq1 < (q_t)100) {
        CHK1("erf",  mf::erf(f1),  erfq(q1),  mpfr::erf(m1));
        CHK1("erfc", mf::erfc(f1), erfcq(q1), mpfr::erfc(m1));
      }
      // erfcx(x) = exp(x²) * erfc(x). The DD kernel is precision-preserving
      // for the tail-tail cancellation; the qp oracle is the naive product.
      // Gate |x| < 26 so exp(x²) stays well inside qp's exponent range
      // and the product stays representable.
      if (q_isfinite(q1) && aq1 < (q_t)26) {
        CHK1("erfcx",
             mf::detail::from_f64x2(::erfcxdd(mf::detail::to_f64x2(f1))),
             expq(q1 * q1) * erfcq(q1),
             mpfr::exp(m1 * m1) * mpfr::erfc(m1));
      }
      if (q_isfinite(q1) && q1 > (q_t)0 && q1 < (q_t)100) {
        CHK1("tgamma", mf::tgamma(f1), tgammaq(q1), mpfr::gamma(m1));
        CHK1("lgamma", mf::lgamma(f1), lgammaq(q1), mpfr::lngamma(m1));
      }

      CHK2("atan2", mf::atan2(f1, f2), atan2q(q1, q2), mpfr::atan2(m1, m2));

      // atan2pi(y, x) = atan2(y, x) / π. Division by π has low
      // amplification, so full-DD tolerance holds here (unlike forward
      // sin/cos/tanpi — see is_pi_trig).
      CHK2("atan2pi",
           mf::detail::from_f64x2(::atan2pidd(
               mf::detail::to_f64x2(f1), mf::detail::to_f64x2(f2))),
           atan2q(q1, q2) / (q_t)M_PIq,
           mpfr::atan2(m1, m2) / mpfr::const_pi(multifloats_test::kMpfrPrec));

      // Power: positive base, modest exponent.
      q_t aq2 = q2 < 0 ? -q2 : q2;
      if (q_isfinite(q1) && q1 > (q_t)1e-3q && q1 < (q_t)1e3q &&
          q_isfinite(q2) && aq2 < (q_t)30) {
        CHK2("pow", mf::pow(f1, f2), powq(q1, q2), mpfr::pow(m1, m2));
        check("pow_md", mf::pow(f1, mf::float64x2(d2)), powq(q1, (q_t)d2),
              q1, (q_t)d2 MP_ARG(mpfr::pow(m1, to_mp(d2))));
        check("pow_dm", mf::pow(mf::float64x2(d1), f2), powq((q_t)d1, q2),
              (q_t)d1, q2 MP_ARG(mpfr::pow(to_mp(d1), m2)));
      }
      if (q_isfinite(q1) && aq1 < (q_t)1e10q) {
        check("pow_int", mf::pow(f1, mf::float64x2(3.0)), powq(q1, (q_t)3),
              q1, (q_t)3 MP_ARG(mpfr::pow(m1, to_mp(3.0))));
      }

      // scalbn(x, 5) = x · 2^5; the ·32 is exact at any precision.
      CHK1_IF(q_isfinite(q1), "scalbn",
              mf::scalbn(f1, 5), scalbnq(q1, 5), m1 * mp_t(32));

      // 3-argument min/max via nested fmin/fmax.
      if (both_finite) {
        q_t q3 = (q1 + q2) * (q_t)0.5q;
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
        CHK2("min3", mf::fmin(mf::fmin(f1, f2), f3), q_min3, mp_min3);
        CHK2("max3", mf::fmax(mf::fmax(f1, f2), f3), q_max3, mp_max3);
      }

      // Bessel of the first (j) and second (y) kind. Labels match
      // ops.py's bench_key convention ("bj0" / "byn" / ...). y* is only
      // finite on positive arguments. Cap |x| at 200 so the DD recurrence
      // doesn't hit the qp reference's large-argument oscillation regime
      // where the two implementations can disagree by more than the kernel
      // error alone.
      if (q_isfinite(q1) && aq1 < (q_t)200) {
        CHK1("bj0", mf::detail::from_f64x2(::j0dd(mf::detail::to_f64x2(f1))),
             j0q(q1), mpfr::besselj0(m1));
        CHK1("bj1", mf::detail::from_f64x2(::j1dd(mf::detail::to_f64x2(f1))),
             j1q(q1), mpfr::besselj1(m1));

        static constexpr int kBesselOrders[] = {2, 3, 5, 8};
        for (int n : kBesselOrders) {
          check("bjn",
                mf::detail::from_f64x2(::jndd(n, mf::detail::to_f64x2(f1))),
                jnq(n, q1), q1, (q_t)n MP_ARG(mpfr::besseljn((long)n, m1)));
          if (q1 > (q_t)0) {
            check("byn",
                  mf::detail::from_f64x2(::yndd(n, mf::detail::to_f64x2(f1))),
                  ynq(n, q1), q1, (q_t)n MP_ARG(mpfr::besselyn((long)n, m1)));
          }
        }

        if (q1 > (q_t)0) {
          CHK1("by0", mf::detail::from_f64x2(::y0dd(mf::detail::to_f64x2(f1))),
               y0q(q1), mpfr::bessely0(m1));
          CHK1("by1", mf::detail::from_f64x2(::y1dd(mf::detail::to_f64x2(f1))),
               y1q(q1), mpfr::bessely1(m1));

          // yn_rangedd: single forward-recurrence sweep filling out[0..5].
          // Each output is compared individually against ynq(n, x). One
          // label "yn_range" covers all six so the stat row aggregates
          // error across the entire range sweep.
          float64x2_t yn_out[6];
          ::yn_rangedd(0, 5, mf::detail::to_f64x2(f1), yn_out);
          for (int n = 0; n <= 5; ++n) {
            check("yn_range",
                  mf::detail::from_f64x2(yn_out[n]), ynq(n, q1), q1, (q_t)n
                  MP_ARG(mpfr::besselyn((long)n, m1)));
          }
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
        CHK1_IF(fabsq(sp) > (q_t)1e-10q, "sinpi",
                mf::detail::from_f64x2(::sinpidd(mf::detail::to_f64x2(f1))),
                sp, mpfr::sin(mpfr::const_pi(multifloats_test::kMpfrPrec) * m1));
        CHK1_IF(fabsq(cp) > (q_t)1e-10q, "cospi",
                mf::detail::from_f64x2(::cospidd(mf::detail::to_f64x2(f1))),
                cp, mpfr::cos(mpfr::const_pi(multifloats_test::kMpfrPrec) * m1));
        CHK1_IF(fabsq(cp) > (q_t)1e-10q && fabsq(sp) > (q_t)1e-10q, "tanpi",
                mf::detail::from_f64x2(::tanpidd(mf::detail::to_f64x2(f1))),
                sp / cp,
                mpfr::tan(mpfr::const_pi(multifloats_test::kMpfrPrec) * m1));
      }
      if (q_isfinite(q1) && aq1 <= (q_t)1) {
        CHK1("asinpi",
             mf::detail::from_f64x2(::asinpidd(mf::detail::to_f64x2(f1))),
             asinq(q1) / (q_t)M_PIq,
             mpfr::asin(m1) / mpfr::const_pi(multifloats_test::kMpfrPrec));
        CHK1("acospi",
             mf::detail::from_f64x2(::acospidd(mf::detail::to_f64x2(f1))),
             acosq(q1) / (q_t)M_PIq,
             mpfr::acos(m1) / mpfr::const_pi(multifloats_test::kMpfrPrec));
      }
      CHK1_IF(q_isfinite(q1), "atanpi",
              mf::detail::from_f64x2(::atanpidd(mf::detail::to_f64x2(f1))),
              atanq(q1) / (q_t)M_PIq,
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

          check_cplx("cdd_add", ::cadddd(zd1, zd2), zq1 + zq2, mag
                     MP_ARG(c_add(zm1, zm2)));
          check_cplx("cdd_sub", ::csubdd(zd1, zd2), zq1 - zq2, mag
                     MP_ARG(c_sub(zm1, zm2)));
          check_cplx("cdd_mul", ::cmuldd(zd1, zd2), zq1 * zq2, mag
                     MP_ARG(c_mul(zm1, zm2)));
          if (mag2 > (q_t)0) {
            check_cplx("cdd_div", ::cdivdd(zd1, zd2), zq1 / zq2, mag
                       MP_ARG(c_div(zm1, zm2)));
          }
          check("cdd_abs",
                mf::detail::from_f64x2(::cabsdd(zd1)), cabsq(zq1),
                mag, (q_t)0 MP_ARG(c_abs(zm1)));
          check_cplx("cdd_conjg", ::conjdd(zd1), conjq(zq1), mag
                     MP_ARG(c_conj(zm1)));

          check_cplx("cdd_sqrt", ::csqrtdd(zd1), csqrtq(zq1), mag
                     MP_ARG(c_sqrt(zm1)));
          check_cplx("cdd_exp",  ::cexpdd(zd1),  cexpq(zq1),  mag
                     MP_ARG(c_exp(zm1)));

          // cexpm1: the DD kernel preserves precision for small |z|. The
          // qp oracle (cexpq(z) - 1) loses it there, so gate |z| > 1e-5
          // to keep the oracle meaningful.
          if (mag1 > (q_t)1e-5q) {
            __complex128 one_c; __real__ one_c = (q_t)1; __imag__ one_c = (q_t)0;
            check_cplx("cdd_expm1", ::cexpm1dd(zd1), cexpq(zq1) - one_c, mag
                       MP_ARG(CMp{c_exp(zm1).re - 1, c_exp(zm1).im}));
          }

          // clog: guard |z| away from 0 (and away from overflow) so the
          // oracle's arg computation stays accurate.
          if (mag1 > (q_t)1e-200q && mag1 < (q_t)1e100q) {
            check_cplx("cdd_log", ::clogdd(zd1), clogq(zq1), mag
                       MP_ARG(c_log(zm1)));

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
            MP_STMT(CMp cdd_log2_mp  = {c_log(zm1).re / mpfr::log(mp_t(2)),
                                        c_log(zm1).im / mpfr::log(mp_t(2))});
            MP_STMT(CMp cdd_log10_mp = {c_log(zm1).re / mpfr::log(mp_t(10)),
                                        c_log(zm1).im / mpfr::log(mp_t(10))});
            check_cplx("cdd_log2",  ::clog2dd(zd1),  lg2,  mag MP_ARG(cdd_log2_mp));
            check_cplx("cdd_log10", ::clog10dd(zd1), lg10, mag MP_ARG(cdd_log10_mp));
          }

          // clog1p: oracle = clogq(1 + z). Gate |1+z| to avoid oracle
          // blow-up on the branch cut at z=-1.
          {
            __complex128 one_c; __real__ one_c = (q_t)1; __imag__ one_c = (q_t)0;
            __complex128 one_plus_z = one_c + zq1;
            if (cabsq(one_plus_z) > (q_t)1e-5q) {
              check_cplx("cdd_log1p", ::clog1pdd(zd1), clogq(one_plus_z), mag
                         MP_ARG(c_log({zm1.re + 1, zm1.im})));
            }
          }

          // cpow: keep base and exponent modest — pow amplifies input
          // precision, so a wide exponent would swamp the tolerance.
          if (mag1 > (q_t)1e-2q && mag1 < (q_t)1e2q && mag2 < (q_t)10) {
            check_cplx("cdd_pow", ::cpowdd(zd1, zd2), cpowq(zq1, zq2), mag
                       MP_ARG(c_exp(c_mul(zm2, c_log(zm1)))));
          }

          // csinpi / ccospi: forward π-scaled complex trig. Same π
          // representation mismatch as scalar {sin,cos}pi — reduced_dd
          // tolerance applies (see is_reduced_dd).
          {
            __complex128 pi_c; __real__ pi_c = (q_t)M_PIq; __imag__ pi_c = (q_t)0;
            __complex128 pi_z = pi_c * zq1;
            MP_STMT(mp_t mpi = mpfr::const_pi(multifloats_test::kMpfrPrec));
            check_cplx("cdd_sinpi", ::csinpidd(zd1), csinq(pi_z), mag
                       MP_ARG(c_sin({mpi * zm1.re, mpi * zm1.im})));
            check_cplx("cdd_cospi", ::ccospidd(zd1), ccosq(pi_z), mag
                       MP_ARG(c_cos({mpi * zm1.re, mpi * zm1.im})));
          }

          // cproj: finite inputs → identity; only meaningful distinction
          // from a plain copy happens on infinite inputs (which the
          // narrow generator doesn't produce). Still worth testing to
          // pin the finite path.
          check_cplx("cdd_proj", ::cprojdd(zd1), cprojq(zq1), mag MP_ARG(zm1));

          // carg: complex → real scalar phase. Near ±π on the branch cut
          // (negative real axis) DD and qp agree as long as imag-limbs
          // sign matches, which they do for DD-representable inputs.
          check("cdd_arg",
                mf::detail::from_f64x2(::cargdd(zd1)), cargq(zq1), mag, (q_t)0
                MP_ARG(c_arg(zm1)));

          check_cplx("cdd_sin",  ::csindd(zd1),  csinq(zq1),  mag MP_ARG(c_sin(zm1)));
          check_cplx("cdd_cos",  ::ccosdd(zd1),  ccosq(zq1),  mag MP_ARG(c_cos(zm1)));
          check_cplx("cdd_sinh", ::csinhdd(zd1), csinhq(zq1), mag MP_ARG(c_sinh(zm1)));
          check_cplx("cdd_cosh", ::ccoshdd(zd1), ccoshq(zq1), mag MP_ARG(c_cosh(zm1)));

          // c{tan,tanh}: zeros of the corresponding {cos,cosh} cause the
          // division inside the oracle to blow up. Gate on oracle magnitude.
          __complex128 ccz  = ccosq(zq1);
          __complex128 ccsh = ccoshq(zq1);
          if (cabsq(ccz)  > (q_t)1e-10q) {
            check_cplx("cdd_tan",  ::ctandd(zd1),  ctanq(zq1),  mag MP_ARG(c_tan(zm1)));
          }
          if (cabsq(ccsh) > (q_t)1e-10q) {
            check_cplx("cdd_tanh", ::ctanhdd(zd1), ctanhq(zq1), mag MP_ARG(c_tanh(zm1)));
          }

          check_cplx("cdd_asin",  ::casindd(zd1),  casinq(zq1),  mag MP_ARG(c_asin(zm1)));
          check_cplx("cdd_acos",  ::cacosdd(zd1),  cacosq(zq1),  mag MP_ARG(c_acos(zm1)));
          check_cplx("cdd_atan",  ::catandd(zd1),  catanq(zq1),  mag MP_ARG(c_atan(zm1)));
          check_cplx("cdd_asinh", ::casinhdd(zd1), casinhq(zq1), mag MP_ARG(c_asinh(zm1)));
          check_cplx("cdd_acosh", ::cacoshdd(zd1), cacoshq(zq1), mag MP_ARG(c_acosh(zm1)));
          check_cplx("cdd_atanh", ::catanhdd(zd1), catanhq(zq1), mag MP_ARG(c_atanh(zm1)));
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
