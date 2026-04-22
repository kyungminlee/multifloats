// Unit tests for multifloats.hh.
//
// Each public operation is tested over a curated list of inputs and the
// result is compared against a __float128 reference. A per-operation
// precision report is printed at the end.
//
// Built with g++ + libquadmath (Apple Clang lacks __float128).

#include "multifloats.hh"
#include "multifloats_td.hh"
#include "multifloats_c.h"
#include "test_common.hh"

#include <cmath>
#include <complex>
#include <cstdio>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace mf = multifloats;
using multifloats_test::q_t;
using multifloats_test::to_q;
using multifloats_test::from_q;
using multifloats_test::q_rel_err;
using multifloats_test::qstr;

// =============================================================================
// Stats / reporting
// =============================================================================

struct Stats {
  double max_rel = 0.0;
  double sum_rel = 0.0;
  int count = 0;
  void update(double rel) {
    if (!std::isfinite(rel)) {
      return;
    }
    if (rel > max_rel) {
      max_rel = rel;
    }
    sum_rel += rel;
    ++count;
  }
};

static int g_failures = 0;
static int g_checks = 0;

#define REQUIRE(cond)                                                          \
  do {                                                                         \
    ++g_checks;                                                                \
    if (!(cond)) {                                                             \
      ++g_failures;                                                            \
      std::fprintf(stderr, "FAIL %s:%d  REQUIRE(%s)\n", __FILE__, __LINE__,    \
                   #cond);                                                     \
    }                                                                          \
  } while (0)

// Compare a multifloat result against a __float128 reference, accumulate
// statistics, and fail loudly if the relative error exceeds `tol`.
static void check_q(char const *op, mf::float64x2 const &got, q_t expected, double tol,
                    Stats &stats, mf::float64x2 const *a = nullptr,
                    mf::float64x2 const *b = nullptr) {
  ++g_checks;
  q_t got_q = to_q(got);
  double rel = q_rel_err(got_q, expected);
  stats.update(rel);
  if (!(rel <= tol)) {
    ++g_failures;
    std::fprintf(stderr, "FAIL [%s] rel_err=%g > tol=%g\n", op, rel, tol);
    if (a) {
      std::fprintf(stderr, "  a        = (%a, %a)\n", a->_limbs[0],
                   a->_limbs[1]);
    }
    if (b) {
      std::fprintf(stderr, "  b        = (%a, %a)\n", b->_limbs[0],
                   b->_limbs[1]);
    }
    std::fprintf(stderr, "  expected = %s\n", qstr(expected));
    std::fprintf(stderr, "  got      = %s\n", qstr(got_q));
    std::fprintf(stderr, "  got_limbs= (%a, %a)\n", got._limbs[0],
                 got._limbs[1]);
  }
}

// =============================================================================
// Tolerances
// =============================================================================
//
// float64x2 has ~104 bits of precision. The kernels lose a few
// bits to ordering/normalization, so we accept up to ~2^-100 ≈ 8e-31 for
// linear ops and ~2^-95 for division (one Newton step from 1/y[0]).

static constexpr double kAddTol = 1e-30;
static constexpr double kMulTol = 1e-30;
static constexpr double kDivTol = 1e-28;

// =============================================================================
// Inputs
// =============================================================================

// All inputs are constructed by rounding a __float128 value down to a
// normalized double-double. This guarantees that to_q() is the exact
// inverse of the construction, so the __float128 reference is exact for
// every operation.
static std::vector<mf::float64x2> sample_inputs() {
  std::vector<mf::float64x2> v;
  auto push_q = [&v](q_t value) { v.push_back(from_q(value)); };

  // Pure doubles.
  push_q((q_t)0);
  push_q((q_t)1);
  push_q((q_t)-1);
  push_q((q_t)2);
  push_q((q_t)0.5);
  push_q((q_t)3);
  push_q((q_t)0.1);
  push_q((q_t)-0.1);
  push_q((q_t)123456789.0);

  // Irrational doubles with a non-trivial lo limb (1/3, sqrt(2), pi, e).
  push_q((q_t)1 / (q_t)3);
  push_q(sqrtq((q_t)2));
  push_q(M_PIq);
  push_q(-M_Eq);

  // 1 + tiny perturbations exercising the lo limb.
  push_q((q_t)1 + scalbnq((q_t)1, -53));
  push_q((q_t)1 - scalbnq((q_t)1, -54));
  push_q((q_t)1 + scalbnq((q_t)1, -80));

  // Large and small magnitudes (still representable as a normalized
  // double-double whose total span fits within float128 precision).
  push_q(scalbnq((q_t)1, 100));
  push_q(scalbnq((q_t)1, -100));
  push_q(scalbnq(M_PIq, 50));
  push_q(scalbnq(-M_Eq, -50));

  return v;
}

// =============================================================================
// Construction / conversion / comparison (exact properties)
// =============================================================================

static void test_construction_and_conversion() {
  mf::float64x2 z;
  REQUIRE(z._limbs[0] == 0.0);
  REQUIRE(z._limbs[1] == 0.0);

  mf::float64x2 one(1.0);
  REQUIRE(one._limbs[0] == 1.0);
  REQUIRE(one._limbs[1] == 0.0);
  REQUIRE(static_cast<double>(one) == 1.0);
  REQUIRE(to_q(one) == (q_t)1);

  mf::float64x2 negz(-0.0);
  REQUIRE(static_cast<double>(negz) == 0.0);
  REQUIRE(std::signbit(negz._limbs[0]));
}

static void test_equality_and_ordering() {
  for (mf::float64x2 const &a : sample_inputs()) {
    for (mf::float64x2 const &b : sample_inputs()) {
      q_t qa = to_q(a);
      q_t qb = to_q(b);
      REQUIRE((a == b) == (qa == qb));
      REQUIRE((a != b) == (qa != qb));
      REQUIRE((a < b) == (qa < qb));
      REQUIRE((a > b) == (qa > qb));
      REQUIRE((a <= b) == (qa <= qb));
      REQUIRE((a >= b) == (qa >= qb));
    }
  }
}

// =============================================================================
// Arithmetic — compared against __float128
// =============================================================================

static void test_unary(Stats &stats) {
  for (mf::float64x2 const &x : sample_inputs()) {
    q_t qx = to_q(x);
    check_q("neg", -x, -qx, 0.0, stats);
    check_q("pos", +x, qx, 0.0, stats);
  }
}

static void test_addition(Stats &stats) {
  auto inputs = sample_inputs();
  for (mf::float64x2 const &a : inputs) {
    for (mf::float64x2 const &b : inputs) {
      q_t expected = to_q(a) + to_q(b);
      check_q("add", a + b, expected, kAddTol, stats, &a, &b);

      mf::float64x2 c = a;
      c += b;
      check_q("add+=", c, expected, kAddTol, stats, &a, &b);
    }
  }
}

static void test_subtraction(Stats &stats) {
  auto inputs = sample_inputs();
  for (mf::float64x2 const &a : inputs) {
    for (mf::float64x2 const &b : inputs) {
      q_t expected = to_q(a) - to_q(b);
      check_q("sub", a - b, expected, kAddTol, stats);

      mf::float64x2 c = a;
      c -= b;
      check_q("sub-=", c, expected, kAddTol, stats);
    }
  }
}

static void test_multiplication(Stats &stats) {
  auto inputs = sample_inputs();
  for (mf::float64x2 const &a : inputs) {
    for (mf::float64x2 const &b : inputs) {
      q_t expected = to_q(a) * to_q(b);
      check_q("mul", a * b, expected, kMulTol, stats);

      mf::float64x2 c = a;
      c *= b;
      check_q("mul*=", c, expected, kMulTol, stats);
    }
  }
}

static void test_division(Stats &stats) {
  auto inputs = sample_inputs();
  for (mf::float64x2 const &a : inputs) {
    for (mf::float64x2 const &b : inputs) {
      if (b._limbs[0] == 0.0) {
        continue;
      }
      q_t expected = to_q(a) / to_q(b);
      check_q("div", a / b, expected, kDivTol, stats);

      mf::float64x2 c = a;
      c /= b;
      check_q("div/=", c, expected, kDivTol, stats);
    }
  }
}

// Bessel J_n(x) via Miller's backward recurrence is used when n > x. The
// CF start-depth is chosen so the truncated J_{n+k}/J_n is below the DD
// precision floor; before tightening the threshold (1e17 → 1e32) the DD
// answer only matched libquadmath to about double precision. Check that
// a Miller-path case now reaches full DD precision.
static void test_bessel_jn_miller_precision(Stats &stats) {
  struct Case { int n; double x; };
  Case const cases[] = {
      {10, 1.0}, {15, 2.0}, {20, 1.0}, {25, 5.0}, {40, 10.0},
  };
  for (Case c : cases) {
    float64x2_t xdd = {c.x, 0.0};
    float64x2_t got = jndd(c.n, xdd);
    mf::float64x2 g;
    g._limbs[0] = got.hi;
    g._limbs[1] = got.lo;
    q_t expected = jnq(c.n, (q_t)c.x);
    // 2e-31 ≈ 2^-102, a few ulps of DD.
    check_q("bessel_jn_miller", g, expected, 2e-31, stats);
  }
}

// C99 G.6.4.2: csqrt(±0 + ±0·i) returns +0 + ±0·i — the real part is always
// +0 regardless of the sign of the input real component, and the imaginary
// sign is preserved. Regression/guard test so a future "preserve sign of a"
// change (the audit's false-positive finding) doesn't silently violate C99.
static void test_csqrt_zero_branch() {
  auto make_z = [](double rh, double rl, double ih, double il) {
    complex64x2_t z;
    z.re.hi = rh; z.re.lo = rl;
    z.im.hi = ih; z.im.lo = il;
    return z;
  };
  double const pz = +0.0, nz = -0.0;

  complex64x2_t r;
  // csqrt(+0 + 0i) = +0 + 0i.
  r = csqrtdd(make_z(pz, 0.0, pz, 0.0));
  REQUIRE(!std::signbit(r.re.hi));
  REQUIRE(!std::signbit(r.im.hi));
  // csqrt(-0 + 0i) = +0 + 0i (NOT -0 + 0i).
  r = csqrtdd(make_z(nz, 0.0, pz, 0.0));
  REQUIRE(!std::signbit(r.re.hi));
  REQUIRE(!std::signbit(r.im.hi));
  // csqrt(-0 - 0i) = +0 - 0i (conjugate symmetry).
  r = csqrtdd(make_z(nz, 0.0, nz, 0.0));
  REQUIRE(!std::signbit(r.re.hi));
  REQUIRE(std::signbit(r.im.hi));
  // csqrt(+0 - 0i) = +0 - 0i.
  r = csqrtdd(make_z(pz, 0.0, nz, 0.0));
  REQUIRE(!std::signbit(r.re.hi));
  REQUIRE(std::signbit(r.im.hi));
}

// atan2 on zero-hi DDs must consult lo for the sign. If y = (+0, -subnormal)
// the effective value is slightly negative, so atan2 should return a negative
// quadrant, not +0. Regression for a bug where hi=0 caused libm atan2 to be
// called with raw +0 arguments, losing the lo-limb sign.
static void test_atan2_signed_zero() {
  auto make_dd = [](double hi, double lo) {
    mf::float64x2 r;
    r._limbs[0] = hi;
    r._limbs[1] = lo;
    return r;
  };
  double const pos_sub = 1e-320;
  double const neg_sub = -1e-320;
  mf::float64x2 y_neg = make_dd(0.0, neg_sub);
  mf::float64x2 y_pos = make_dd(0.0, pos_sub);
  mf::float64x2 x_neg = make_dd(0.0, neg_sub);
  mf::float64x2 x_pos = make_dd(0.0, pos_sub);

  // y<0, x>0 → answer in (−π/2, 0), specifically ≈ −π/4 for equal magnitudes.
  REQUIRE(mf::atan2(y_neg, x_pos)._limbs[0] < 0.0);
  // y>0, x<0 → answer in (π/2, π), specifically ≈ 3π/4.
  REQUIRE(mf::atan2(y_pos, x_neg)._limbs[0] > 1.5);
  // y<0, x<0 → answer in (−π, −π/2), specifically ≈ −3π/4.
  REQUIRE(mf::atan2(y_neg, x_neg)._limbs[0] < -1.5);
}

// to_string / operator<< round-trip and format sanity checks. The public
// API guarantees ~32 significant decimal digits; this test locks in the
// scientific-notation format, the signed-zero / nan / inf cases, and
// precision honoring via os.precision().
static void test_io_to_string_and_stream() {
  // Signed zero, nan, inf: exact textual match.
  REQUIRE(mf::to_string(mf::float64x2(0.0)) == "0e+00");
  mf::float64x2 nz;
  nz._limbs[0] = -0.0;
  REQUIRE(mf::to_string(nz) == "-0e+00");
  REQUIRE(mf::to_string(mf::float64x2(
              std::numeric_limits<double>::infinity())) == "inf");
  REQUIRE(mf::to_string(mf::float64x2(
              -std::numeric_limits<double>::infinity())) == "-inf");
  REQUIRE(mf::to_string(mf::float64x2(
              std::numeric_limits<double>::quiet_NaN())) == "nan");

  // Round-trip a DD value: compute π as DD, format to 32 digits, parse
  // back via __float128 strtoflt128, and require the relative error to be
  // below 1e-31 (one DD ulp floor).
  mf::float64x2 pi;
  pi._limbs[0] = 3.141592653589793;
  pi._limbs[1] = 1.2246467991473532e-16;  // DD π
  std::string s = mf::to_string(pi, 32);
  q_t parsed = strtoflt128(s.c_str(), nullptr);
  q_t pi_q = to_q(pi);
  q_t rel = fabsq((parsed - pi_q) / pi_q);
  REQUIRE(rel < 1e-31q);

  // Format shape: one leading digit, decimal point, 31 digits, 'e', sign,
  // two-digit exponent. Total length for 32 digits = 1 + 1 + 31 + 1 + 1 + 2 = 37.
  REQUIRE(s.size() == 37);
  REQUIRE(s[1] == '.');
  REQUIRE(s[33] == 'e');
  REQUIRE(s[34] == '+');

  // Precision clamping: request too few or too many, normalize to [1,34].
  REQUIRE(mf::to_string(mf::float64x2(1.5), 0).size() >= 5);
  REQUIRE(mf::to_string(mf::float64x2(1.5), 100).size() > 34);  // 34 digits fits

  // ostream: default precision of 6 should upgrade to 32 automatically.
  std::ostringstream os;
  os << pi;
  REQUIRE(os.str() == s);

  // Explicit precision override via os.precision() honored when > 17.
  std::ostringstream os2;
  os2.precision(20);
  os2 << pi;
  REQUIRE(os2.str().size() == 1 + 1 + 19 + 1 + 1 + 2);  // 25
  REQUIRE(os2.str()[21] == 'e');

  // Negative values: leading minus.
  REQUIRE(mf::to_string(mf::float64x2(-2.5), 3)[0] == '-');

  // Rounding at the boundary: 0.999999... with 3 digits → "1.00e+00".
  mf::float64x2 almost_one;
  almost_one._limbs[0] = 0.9999999999999999;
  almost_one._limbs[1] = 0.0;
  std::string r = mf::to_string(almost_one, 3);
  REQUIRE(r == "1.00e+00");
}

// Curated edge-input checks for log2 / log1p / expm1 / cbrt. These
// kernels were added in the Tier 2 API completion pass and are exercised
// by the fortran_fuzz loop across random inputs; the tests here lock in
// specific identities and precision-sensitive corner cases.
static void test_log_root_edges(Stats &stats) {
  auto mf_q = [&stats](char const *name, mf::float64x2 got, q_t expect,
                       q_t tol) {
    q_t err = fabsq(to_q(got) - expect);
    q_t mag = fabsq(expect) > 1 ? fabsq(expect) : (q_t)1;
    q_t rel = err / mag;
    if (rel > tol) {
      std::printf("FAIL %s rel=%g tol=%g\n", name, (double)rel, (double)tol);
      ++g_failures;
    }
    ++g_checks;
    stats.update(rel);
  };

  q_t const full_dd_tol = 1e-30q;   // ~1 DD ULP headroom.
  q_t const dp_tol      = 5e-15q;   // ~1 dp ULP (cbrt is single-Newton).

  // log2 identities at exact powers of 2.
  for (int k = -60; k <= 60; k += 5) {
    mf::float64x2 x = mf::ldexp(mf::float64x2(1.0), k);
    mf_q("log2(2^k)", mf::log2(x), (q_t)k, full_dd_tol);
  }
  mf_q("log2(1)", mf::log2(mf::float64x2(1.0)), (q_t)0, full_dd_tol);
  mf_q("log2(0.5)", mf::log2(mf::float64x2(0.5)), (q_t)-1, full_dd_tol);

  // log1p: small-|x| path (|x|<0.25 uses Taylor) must keep full DD.
  q_t tiny_vals[] = { (q_t)1e-20, (q_t)1e-10, (q_t)1e-6, (q_t)1e-3, (q_t)0.1 };
  for (q_t v : tiny_vals) {
    mf::float64x2 x = from_q(v);
    // log1p(x) ≈ x - x²/2 + x³/3 - ... Reference via __float128 log1pq.
    q_t ref = log1pq(v);
    mf_q("log1p(small)", mf::log1p(x), ref, full_dd_tol);
    mf_q("log1p(-small)", mf::log1p(from_q(-v)), log1pq(-v), full_dd_tol);
  }
  // log1p(0) must be exactly 0.
  mf::float64x2 lp0 = mf::log1p(mf::float64x2(0.0));
  REQUIRE(lp0._limbs[0] == 0.0 && lp0._limbs[1] == 0.0);

  // expm1: small-|x| path (|x|<0.5 uses Taylor) must keep full DD.
  for (q_t v : tiny_vals) {
    q_t ref = expm1q(v);
    mf_q("expm1(small)", mf::expm1(from_q(v)), ref, full_dd_tol);
    mf_q("expm1(-small)", mf::expm1(from_q(-v)), expm1q(-v), full_dd_tol);
  }
  // expm1(0) must be exactly 0.
  mf::float64x2 em0 = mf::expm1(mf::float64x2(0.0));
  REQUIRE(em0._limbs[0] == 0.0 && em0._limbs[1] == 0.0);
  // expm1(1): |x|>=0.5 routes through `exp(x) - 1`, inheriting the
  // exp kernel's precision (~1e-30 worst case, not Taylor-full-DD).
  mf_q("expm1(1)", mf::expm1(mf::float64x2(1.0)), M_Eq - 1, (q_t)1e-29);

  // cbrt: single-Newton result, ~1 dp ULP. Perfect cubes should come out
  // with relative error well under a dp ULP because the Newton correction
  // off a bit-exact `std::cbrt` seed is zero.
  double cubes[] = { 0.0, 1.0, 8.0, 27.0, 64.0, 125.0, -1.0, -8.0, -125.0 };
  double roots[] = { 0.0, 1.0, 2.0,  3.0,  4.0,   5.0, -1.0, -2.0,  -5.0 };
  for (int i = 0; i < 9; ++i) {
    mf::float64x2 y = mf::cbrt(mf::float64x2(cubes[i]));
    mf_q("cbrt(cube)", y, (q_t)roots[i], dp_tol);
  }
  // Non-cube inputs: compare against cbrtq.
  double non_cubes[] = { 2.0, 3.0, 7.0, 100.0, 1e-10, 1e10 };
  for (double v : non_cubes) {
    mf::float64x2 y = mf::cbrt(mf::float64x2(v));
    mf_q("cbrt", y, cbrtq((q_t)v), dp_tol);
  }
  // Sign preservation and signed-zero propagation.
  mf::float64x2 cb_nz = mf::cbrt(mf::float64x2(-0.0));
  REQUIRE(cb_nz._limbs[0] == 0.0);
}

// C99 Annex G branch-cut behavior for complex log, sqrt, asin, acos on
// the negative-real and imaginary axes. The sign of the 0 imaginary
// part selects which side of the cut to land on:
//   clog(-1 + 0i)   = (0, +π)
//   clog(-1 − 0i)   = (0, −π)
//   csqrt(-1 + 0i)  = (0, +1)
//   csqrt(-1 − 0i)  = (0, −1)
//   casin(2 + 0i)   = ( π/2, +log(2+√3))
//   casin(2 − 0i)   = ( π/2, −log(2+√3))
//   cacos(2 + 0i)   = (0, −log(2+√3))
//   cacos(2 − 0i)   = (0, +log(2+√3))
// Reference is libquadmath (clogq, csqrtq, casinq, cacosq).
static void test_complex_branch_cuts(Stats &stats) {
  using cq_t = __complex128;
  // NOTE: for branch-cut cases the sign of a zero imaginary part is the
  // whole point. The test-harness `from_q` normalizes via fast_two_sum,
  // which collapses (-0, +0) → (+0, +0) and loses signed zero. Bypass
  // it by laying the dp limbs down directly via to_cmf_dp.
  auto to_cmf_dp = [](double re_hi, double im_hi) {
    complex64x2_t z;
    z.re.hi = re_hi; z.re.lo = 0.0;
    z.im.hi = im_hi; z.im.lo = 0.0;
    return z;
  };
  auto check_c = [&stats](char const *name, complex64x2_t got, cq_t ref,
                          q_t tol) {
    q_t got_re = (q_t)got.re.hi + (q_t)got.re.lo;
    q_t got_im = (q_t)got.im.hi + (q_t)got.im.lo;
    q_t ref_re = crealq(ref);
    q_t ref_im = cimagq(ref);
    q_t mag = cabsq(ref); if (mag < 1) mag = 1;
    q_t err = fabsq(got_re - ref_re) + fabsq(got_im - ref_im);
    q_t rel = err / mag;
    if (rel > tol) {
      std::printf("FAIL %s  got=(%g,%g) ref=(%g,%g) rel=%g\n", name,
                  (double)got_re, (double)got_im,
                  (double)ref_re, (double)ref_im, (double)rel);
      ++g_failures;
    }
    ++g_checks;
    stats.update(rel);
  };
  auto check_sign = [](char const *name, double got, int want_sign) {
    bool got_neg = std::signbit(got);
    bool want_neg = (want_sign < 0);
    if (got_neg != want_neg) {
      std::printf("FAIL %s: got %g (signbit=%d), want sign=%d\n", name,
                  got, got_neg, want_sign);
      ++g_failures;
    }
    ++g_checks;
  };

  q_t const tol = 1e-30q;
  double const pz = +0.0;
  // Use copysign to defeat the compiler's constant folding of -0.0 at
  // -O3, which otherwise collapses -0.0 to +0 when it "knows" the
  // literal value (GCC treats -0.0 == 0.0 for constant propagation).
  double const nz = std::copysign(0.0, -1.0);

  // clog on the negative real axis: sign of imag input selects ±π.
  for (double x : {-1.0, -2.5, -1e-10, -1e10}) {
    complex64x2_t zp = to_cmf_dp(x, pz);
    complex64x2_t zn = to_cmf_dp(x, nz);
    cq_t ref_p = (cq_t)x; __imag__ ref_p = pz;
    cq_t ref_n = (cq_t)x; __imag__ ref_n = nz;
    check_c("clog(-x, +0)", clogdd(zp), clogq(ref_p), tol);
    check_c("clog(-x, -0)", clogdd(zn), clogq(ref_n), tol);
    check_sign("clog(-x, +0) imag>0", clogdd(zp).im.hi, +1);
    check_sign("clog(-x, -0) imag<0", clogdd(zn).im.hi, -1);
  }

  // csqrt on the negative real axis: imag takes the sign of imag input.
  for (double x : {-1.0, -4.0, -1e-10, -1e10}) {
    complex64x2_t zp = to_cmf_dp(x, pz);
    complex64x2_t zn = to_cmf_dp(x, nz);
    cq_t ref_p = (cq_t)x; __imag__ ref_p = pz;
    cq_t ref_n = (cq_t)x; __imag__ ref_n = nz;
    check_c("csqrt(-x, +0)", csqrtdd(zp), csqrtq(ref_p), tol);
    check_c("csqrt(-x, -0)", csqrtdd(zn), csqrtq(ref_n), tol);
    check_sign("csqrt(-x, +0) imag>0", csqrtdd(zp).im.hi, +1);
    check_sign("csqrt(-x, -0) imag<0", csqrtdd(zn).im.hi, -1);
  }

  // casin / cacos on the real axis with |Re z| > 1 have branch cuts at
  // the imag sign. C99 G.6.2.2:
  //   casin(x + 0i) for x > 1 = π/2 + i·log(x+√(x²−1))
  //   casin(x − 0i) for x > 1 = π/2 − i·log(x+√(x²−1))
  //   cacos(z) = π/2 − casin(z), so imag sign is the OPPOSITE of input.
  for (double x : {2.0, 5.0, -2.0, -5.0}) {
    complex64x2_t zp = to_cmf_dp(x, pz);
    complex64x2_t zn = to_cmf_dp(x, nz);
    cq_t ref_p = (cq_t)x; __imag__ ref_p = pz;
    cq_t ref_n = (cq_t)x; __imag__ ref_n = nz;
    check_c("casin(|x|>1,+0)", casindd(zp), casinq(ref_p), (q_t)1e-28);
    check_c("casin(|x|>1,-0)", casindd(zn), casinq(ref_n), (q_t)1e-28);
    check_c("cacos(|x|>1,+0)", cacosdd(zp), cacosq(ref_p), (q_t)1e-28);
    check_c("cacos(|x|>1,-0)", cacosdd(zn), cacosq(ref_n), (q_t)1e-28);
  }

  // catanh branch-point singularities at z = ±1 + 0i (C99 G.6.2.4).
  //   catanh(+1 + 0i) = +∞ + 0i
  //   catanh(−1 + 0i) = −∞ + 0i
  complex64x2_t r_p = catanhdd(to_cmf_dp(+1.0, pz));
  complex64x2_t r_n = catanhdd(to_cmf_dp(-1.0, pz));
  REQUIRE(std::isinf(r_p.re.hi) && r_p.re.hi > 0);
  REQUIRE(std::isinf(r_n.re.hi) && r_n.re.hi < 0);
  REQUIRE(r_p.im.hi == 0.0 && r_n.im.hi == 0.0);
}

// Huge-argument trig range-reduction sanity. sin(2π·k) is 0 exactly for
// any integer k; numerically it's bounded by whatever reduction precision
// the kernel uses. DD sin should keep the result well below 1 even out
// to 2^40 range — a failure here means the Cody-Waite / Payne-Hanek
// reduction has broken.
static void test_huge_argument_trig(Stats &stats) {
  // π in DD (standard split): hi = 3.141592..., lo ≈ 1.2246e-16.
  static const mf::float64x2 pi_dd(3.141592653589793,
                                   1.2246467991473532e-16);
  // 2π = 2 · π_dd, exact in DD since multiplication by 2 is a pure shift.
  static const mf::float64x2 two_pi_dd = mf::float64x2(2.0) * pi_dd;
  static const mf::float64x2 half_pi_dd = mf::float64x2(0.5) * pi_dd;

  // For each k ∈ {2^N}, compute y = 2π·k in DD (exact scaling) and check
  // that sin(y) / cos(y) are near 0 / 1 respectively. Tolerance grows
  // with k because DD(2π) is inexact at ~2^-106 relative, so absolute
  // reduction error is ~k·2^-106 — looser than 1 DD ULP but still tiny.
  for (int N = 0; N <= 40; N += 4) {
    double scale = std::ldexp(1.0, N);
    mf::float64x2 k(scale);
    mf::float64x2 y = two_pi_dd * k;       // 2π·k; scale exact, so this is
                                           // DD(2π) · 2^N, error ~2^-106 · 2π · 2^N.
    q_t tol = (q_t)scale * (q_t)1e-25;      // slack: ~100× the 1-DD-ulp floor.
    if (tol < (q_t)1e-25) tol = (q_t)1e-25;

    mf::float64x2 s = mf::sin(y);
    mf::float64x2 c = mf::cos(y);
    q_t s_mag = fabsq(to_q(s));
    q_t c_diff = fabsq(to_q(c) - 1);
    if (s_mag > tol) {
      std::printf("FAIL sin(2pi*2^%d) = %g, tol=%g\n", N,
                  (double)s_mag, (double)tol);
      ++g_failures;
    }
    if (c_diff > tol) {
      std::printf("FAIL cos(2pi*2^%d)-1 = %g, tol=%g\n", N,
                  (double)c_diff, (double)tol);
      ++g_failures;
    }
    g_checks += 2;
    stats.update(s_mag / (q_t)std::max(1.0, scale));
  }

  // sin(π/2 + 2π·k) must stay close to 1 for large k.
  for (int N = 10; N <= 30; N += 5) {
    double scale = std::ldexp(1.0, N);
    mf::float64x2 y = two_pi_dd * mf::float64x2(scale) + half_pi_dd;
    mf::float64x2 s = mf::sin(y);
    q_t err = fabsq(to_q(s) - 1);
    q_t tol = (q_t)scale * (q_t)1e-25;
    if (tol < (q_t)1e-25) tol = (q_t)1e-25;
    if (err > tol) {
      std::printf("FAIL sin(pi/2 + 2pi*2^%d) rel err %g, tol=%g\n", N,
                  (double)err, (double)tol);
      ++g_failures;
    }
    ++g_checks;
    stats.update(err);
  }

  // sinpi / cospi (C ABI only — no C++ template specialization). The
  // π-scaled form keeps the multiply by π implicit, so integer inputs
  // must yield exact zero (sin) / ±1 (cos) regardless of magnitude.
  for (int k = 1; k <= 20; ++k) {
    float64x2_t xk = {(double)k, 0.0};
    float64x2_t sp = sinpidd(xk);
    REQUIRE(sp.hi == 0.0 && sp.lo == 0.0);
  }
  for (int k = -5; k <= 5; ++k) {
    float64x2_t xk = {(double)k, 0.0};
    float64x2_t cp = cospidd(xk);
    REQUIRE(cp.hi == ((k & 1) ? -1.0 : 1.0));
    REQUIRE(cp.lo == 0.0);
  }
}

// Complex multiply must keep the 4-mul form (Re = ac − bd, Im = ad + bc).
// A naive Karatsuba rewrite (3 muls, 5 adds) would replace Im with
//   r − p − q  where r = (a+b)(c+d), p = ac, q = bd
// which catastrophically loses precision when the true Im has internal
// cancellation. For the classic witness a = (1, ε), b = (-1, ε) the
// true imaginary part is 0 exactly; 4-mul returns 0, Karatsuba returns
// −ε². See "Designs measured and rejected" in
// doc/developer/INTERNALS.md for the A/B study that rejected Karatsuba.
static void test_complex_mul_cancellation() {
  mf::float64x2 const one(1.0), neg_one(-1.0);
  mf::float64x2 eps;
  eps._limbs[0] = 1e-18;
  eps._limbs[1] = 0.0;

  std::complex<mf::float64x2> a(one, eps);
  std::complex<mf::float64x2> b(neg_one, eps);
  std::complex<mf::float64x2> prod = a * b;

  // Im(a·b) = 1·ε + ε·(−1) = 0 exactly. Both limbs must be +0.
  REQUIRE(prod.imag()._limbs[0] == 0.0);
  REQUIRE(prod.imag()._limbs[1] == 0.0);
  REQUIRE(!std::signbit(prod.imag()._limbs[0]));

  // Re(a·b) = (1)(−1) − ε·ε = −1 − ε² (full DD precision).
  q_t ref_re = (q_t)(-1.0) - (q_t)1e-18 * (q_t)1e-18;
  q_t got_re = to_q(prod.real());
  q_t err = fabsq(got_re - ref_re) / fabsq(ref_re);
  REQUIRE((double)err < 1e-31);
}

// Non-finite division must produce a DD where BOTH limbs are non-finite, so
// an `isnan(lo)` check on the result reports the same classification as the
// hi limb. Regression for a bug where NaN-producing division left lo == 0.
static void test_division_nonfinite() {
  double inf = std::numeric_limits<double>::infinity();
  mf::float64x2 one(1.0);
  mf::float64x2 zero(0.0);
  mf::float64x2 pinf(inf);

  mf::float64x2 nan_result = zero / zero;
  REQUIRE(std::isnan(nan_result._limbs[0]));
  REQUIRE(std::isnan(nan_result._limbs[1]));

  mf::float64x2 inf_result = one / zero;
  REQUIRE(std::isinf(inf_result._limbs[0]));
  REQUIRE(std::isinf(inf_result._limbs[1]));
  REQUIRE(inf_result._limbs[0] > 0 && inf_result._limbs[1] > 0);

  mf::float64x2 neg_inf = (mf::float64x2(-1.0)) / zero;
  REQUIRE(std::isinf(neg_inf._limbs[0]));
  REQUIRE(std::isinf(neg_inf._limbs[1]));
  REQUIRE(neg_inf._limbs[0] < 0 && neg_inf._limbs[1] < 0);

  // Finite / Inf is legitimately 0; default-init of lo=0 is correct.
  mf::float64x2 zero_result = one / pinf;
  REQUIRE(zero_result._limbs[0] == 0.0);
  REQUIRE(zero_result._limbs[1] == 0.0);
}

// =============================================================================
// cmath-style wrappers
// =============================================================================

static void test_abs_fmin_fmax(Stats &stats) {
  auto inputs = sample_inputs();
  for (mf::float64x2 const &x : inputs) {
    q_t qx = to_q(x);
    q_t qexp = qx < 0 ? -qx : qx;
    check_q("abs", mf::abs(x), qexp, 0.0, stats);
    check_q("fabs", mf::fabs(x), qexp, 0.0, stats);
  }
  for (mf::float64x2 const &a : inputs) {
    for (mf::float64x2 const &b : inputs) {
      q_t qa = to_q(a);
      q_t qb = to_q(b);
      check_q("fmin", mf::fmin(a, b), qa < qb ? qa : qb, 0.0, stats);
      check_q("fmax", mf::fmax(a, b), qa < qb ? qb : qa, 0.0, stats);
    }
  }
}

static void test_copysign(Stats &stats) {
  auto inputs = sample_inputs();
  for (mf::float64x2 const &x : inputs) {
    for (mf::float64x2 const &y : inputs) {
      q_t qx = to_q(x);
      bool xs = std::signbit(x._limbs[0]);
      bool ys = std::signbit(y._limbs[0]);
      q_t qexp = (xs == ys) ? qx : -qx;
      // Match the wrapper's exact bit-pattern semantics on the limb.
      check_q("copysign", mf::copysign(x, y), qexp, 0.0, stats);
    }
  }
}

static void test_classification() {
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  mf::float64x2 finite(1.5);
  mf::float64x2 zero(0.0);
  mf::float64x2 neg(-1.5);
  mf::float64x2 pinf(inf), ninf(-inf), mnan(nan);

  REQUIRE(!mf::signbit(finite));
  REQUIRE(mf::signbit(neg));
  REQUIRE(mf::isfinite(finite));
  REQUIRE(!mf::isfinite(pinf));
  REQUIRE(!mf::isfinite(mnan));
  REQUIRE(mf::isinf(pinf));
  REQUIRE(mf::isinf(ninf));
  REQUIRE(!mf::isinf(finite));
  REQUIRE(mf::isnan(mnan));
  REQUIRE(!mf::isnan(finite));
  REQUIRE(mf::fpclassify(finite) == FP_NORMAL);
  REQUIRE(mf::fpclassify(zero) == FP_ZERO);
  REQUIRE(mf::fpclassify(pinf) == FP_INFINITE);
  REQUIRE(mf::fpclassify(mnan) == FP_NAN);

  // Non-canonical DDs where hi is (signed) zero: the value lives in lo,
  // so signbit/abs must consult lo rather than blindly trusting hi.
  mf::float64x2 tiny_neg;
  tiny_neg._limbs[0] = 0.0;
  tiny_neg._limbs[1] = -1e-300;
  REQUIRE(mf::signbit(tiny_neg));
  REQUIRE(!mf::signbit(mf::abs(tiny_neg)));

  mf::float64x2 tiny_pos;
  tiny_pos._limbs[0] = -0.0;
  tiny_pos._limbs[1] = 1e-300;
  REQUIRE(!mf::signbit(tiny_pos));
  REQUIRE(!mf::signbit(mf::abs(tiny_pos)));
}

static void test_ldexp_scalbn_ilogb(Stats &stats) {
  auto inputs = sample_inputs();
  int const ns[] = {-50, -1, 0, 1, 25};
  for (mf::float64x2 const &x : inputs) {
    if (!mf::isfinite(x) || x == mf::float64x2(0.0)) {
      continue;
    }
    q_t qx = to_q(x);
    for (int n : ns) {
      q_t qexp = qx * scalbnq((q_t)1, n);
      check_q("ldexp", mf::ldexp(x, n), qexp, 0.0, stats);
      check_q("scalbn", mf::scalbn(x, n), qexp, 0.0, stats);
    }
    // mf::ilogb is documented as delegating to the leading limb (matches
    // Julia/Fortran exponent semantics).
    REQUIRE(mf::ilogb(x) == std::ilogb(x._limbs[0]));
  }
}

// =============================================================================
// lerp — C++20 specification properties
// =============================================================================
//
// `lerp(a, b, t)` must satisfy:
//   1. lerp(a, b, 0) == a  (exact)
//   2. lerp(a, b, 1) == b  (exact)
//   3. lerp(a, a, t) == a  (exact, for any finite t)
//   4. monotonicity: if t1 < t2 then lerp(a, b, t1) <= lerp(a, b, t2) when a<b
//   5. accuracy: for t ∈ [0, 1] and same-sign a,b the result has no cancellation
//   6. NaN propagation: if any argument is NaN the result is NaN
static void test_lerp(Stats &stats) {
  mf::float64x2 const cases[][2] = {
      {mf::float64x2(0.0), mf::float64x2(1.0)},
      {mf::float64x2(-1.0), mf::float64x2(1.0)},   // opposite signs
      {mf::float64x2(1.0), mf::float64x2(2.0)},    // same sign, ascending
      {mf::float64x2(-2.0), mf::float64x2(-1.0)},  // same sign, both negative
      {mf::float64x2(3.0), mf::float64x2(3.0)},    // a == b
      {from_q(M_PIq), from_q(M_Eq)},
      {from_q(scalbnq((q_t)1, 50)), from_q(scalbnq((q_t)1, -50))},
  };
  double const ts[] = {0.0, 0.25, 0.5, 0.75, 1.0, -0.5, 1.5, 2.0};

  for (auto const &p : cases) {
    mf::float64x2 const &a = p[0];
    mf::float64x2 const &b = p[1];
    q_t qa = to_q(a), qb = to_q(b);
    // Endpoint exactness (spec requirement, not tolerance-based).
    REQUIRE(mf::lerp(a, b, mf::float64x2(0.0)) == a);
    REQUIRE(mf::lerp(a, b, mf::float64x2(1.0)) == b);
    // a == b: every t returns a.
    if (a == b) {
      for (double t : ts) {
        REQUIRE(mf::lerp(a, b, mf::float64x2(t)) == a);
      }
      continue;
    }
    // Value check vs __float128 reference at intermediate t.
    for (double t : ts) {
      mf::float64x2 r = mf::lerp(a, b, mf::float64x2(t));
      q_t expected = qa + (q_t)t * (qb - qa);
      // For t in [0,1], lerp is exact to DD precision.
      // For extrapolation, allow the standard DD add/mul tolerance.
      check_q("lerp", r, expected, 1e-30, stats, &a, &b);
    }
  }

  // NaN propagation.
  mf::float64x2 nan_val;
  nan_val._limbs[0] = std::numeric_limits<double>::quiet_NaN();
  nan_val._limbs[1] = 0.0;
  REQUIRE(mf::isnan(mf::lerp(nan_val, mf::float64x2(1.0), mf::float64x2(0.5))));
  REQUIRE(mf::isnan(mf::lerp(mf::float64x2(0.0), nan_val, mf::float64x2(0.5))));
  REQUIRE(mf::isnan(mf::lerp(mf::float64x2(0.0), mf::float64x2(1.0), nan_val)));

  // Monotonicity sweep on a well-conditioned same-sign pair.
  mf::float64x2 a(1.0), b(2.0);
  mf::float64x2 prev = mf::lerp(a, b, mf::float64x2(0.0));
  for (int i = 1; i <= 100; ++i) {
    mf::float64x2 cur = mf::lerp(a, b, mf::float64x2((double)i / 100.0));
    REQUIRE(cur >= prev);
    prev = cur;
  }
}

// =============================================================================
// nextafter edge cases (power-of-2 boundaries)
// =============================================================================
//
// At |hi| = 2^k the downward ulp(hi) halves (2^(k-53) vs the upward 2^(k-52)),
// so the old formula eps = ldexp(ulp, -53) made the down-step 2× too small.
// The fix: always base eps on the UP-side ulp of |hi|, so the step is
// symmetric in both directions and nextafter(nextafter(x, +inf), -inf) == x.
static void test_nextafter_symmetry() {
  double const powers_of_two[] = {
      0x1p-50, 0x1p-10, 0x1p-1, 0x1p0, 0x1p1, 0x1p2, 0x1p3, 0x1p10, 0x1p50,
  };
  mf::float64x2 const pinf(std::numeric_limits<double>::infinity());
  mf::float64x2 const ninf(-std::numeric_limits<double>::infinity());

  for (double p2 : powers_of_two) {
    // Both signs — sign handling on the hi limb must not flip the ulp.
    for (double sign : {+1.0, -1.0}) {
      mf::float64x2 x(sign * p2);
      mf::float64x2 up = mf::nextafter(x, pinf);
      mf::float64x2 dn = mf::nextafter(x, ninf);

      // Both steps must be nonzero (no collapse to x).
      REQUIRE(up != x);
      REQUIRE(dn != x);
      // Step directions are correct.
      REQUIRE(up > x);
      REQUIRE(dn < x);
      // Round-trip identity: step and return lands exactly on x.
      REQUIRE(mf::nextafter(up, ninf) == x);
      REQUIRE(mf::nextafter(dn, pinf) == x);
      // Step magnitudes are symmetric in DD (not halved on one side).
      q_t qx = to_q(x);
      q_t qup = to_q(up);
      q_t qdn = to_q(dn);
      q_t up_step = qup - qx;
      q_t dn_step = qx - qdn;
      // Equality as q_t: we set up_step/dn_step to match bit-for-bit.
      REQUIRE(up_step == dn_step);
      // Step equals expected eps = 2^-105 * p2 (i.e., DD-rel ulp × |x|).
      q_t expected = scalbnq((q_t)p2, -105);
      REQUIRE(up_step == expected);
    }
  }

  // Non-power-of-2 values: step up and down should be equal (trivially —
  // ulp is unambiguous there) and round-trip still identity.
  double const non_p2[] = {
      1.5, 1.25, 1.125, 3.0, 5.0, 7.0, 0.75, 0.375, 0x1.23456789abcdep10,
  };
  for (double v : non_p2) {
    for (double sign : {+1.0, -1.0}) {
      mf::float64x2 x(sign * v);
      mf::float64x2 up = mf::nextafter(x, pinf);
      mf::float64x2 dn = mf::nextafter(x, ninf);
      REQUIRE(up != x);
      REQUIRE(dn != x);
      REQUIRE(mf::nextafter(up, ninf) == x);
      REQUIRE(mf::nextafter(dn, pinf) == x);
      q_t up_step = to_q(up) - to_q(x);
      q_t dn_step = to_q(x) - to_q(dn);
      REQUIRE(up_step == dn_step);
    }
  }

  // Identity: x == y implies nextafter returns y unchanged (exact).
  for (double p2 : powers_of_two) {
    mf::float64x2 x(p2);
    REQUIRE(mf::nextafter(x, x) == x);
  }

  // Non-trivial DD inputs (hi and lo both present): round-trip should hold.
  q_t const seeds[] = {
      (q_t)1 + scalbnq((q_t)1, -60),
      M_PIq,
      sqrtq((q_t)2) * scalbnq((q_t)1, 30),
      (q_t)1 / (q_t)7,
  };
  for (q_t seed : seeds) {
    mf::float64x2 x = from_q(seed);
    mf::float64x2 up = mf::nextafter(x, pinf);
    mf::float64x2 dn = mf::nextafter(x, ninf);
    REQUIRE(mf::nextafter(up, ninf) == x);
    REQUIRE(mf::nextafter(dn, pinf) == x);
  }
}

// =============================================================================
// Complex accessors and structural ops (cabs / carg / cproj / conj / real / imag)
// =============================================================================
//
// These are the C ABI symbols added to match libquadmath's surface.
// Reference is __complex128 from libquadmath. Each check runs the DD
// path (via the extern "C" symbol) and compares against cabsq / cargq
// / cprojq / conjq / crealq / cimagq.
static void test_complex_accessors(Stats &stats) {
  using cq_t = __complex128;
  auto to_cq = [](float64x2_t a) {
    return (q_t)a.hi + (q_t)a.lo;
  };

  // Diverse representative points: axes, all four quadrants, tiny and
  // huge magnitudes, signed zero / inf / nan in each component.
  struct Case { q_t re; q_t im; };
  std::vector<Case> cases = {
      {(q_t)0, (q_t)0},  {(q_t)1, (q_t)0},  {(q_t)-1, (q_t)0},
      {(q_t)0, (q_t)1},  {(q_t)0, (q_t)-1}, {(q_t)1, (q_t)1},
      {(q_t)-1, (q_t)1}, {(q_t)-1, (q_t)-1}, {(q_t)1, (q_t)-1},
      {(q_t)3, (q_t)4},  // classic 3-4-5.
      {scalbnq((q_t)1, 600), scalbnq((q_t)1, 600)},   // overflow-range hypot.
      {scalbnq((q_t)1, -600), scalbnq((q_t)1, -600)}, // underflow-range hypot.
      {M_PIq, M_Eq},
  };

  for (Case c : cases) {
    mf::float64x2 re = from_q(c.re);
    mf::float64x2 im = from_q(c.im);
    complex64x2_t z = {
        {re._limbs[0], re._limbs[1]},
        {im._limbs[0], im._limbs[1]},
    };

    cq_t cq;
    __real__ cq = c.re;
    __imag__ cq = c.im;
    // cabs: |z|
    {
      mf::float64x2 got;
      float64x2_t r = cabsdd(z);
      got._limbs[0] = r.hi;
      got._limbs[1] = r.lo;
      check_q("cabs", got, cabsq(cq), 1e-30, stats);
    }
    // carg: atan2(im, re); libquadmath's cargq takes __complex128.
    {
      mf::float64x2 got;
      float64x2_t r = cargdd(z);
      got._limbs[0] = r.hi;
      got._limbs[1] = r.lo;
      check_q("carg", got, cargq(cq), 1e-30, stats);
    }
    // creal / cimag: exact field accessors — must be bitwise equal.
    {
      float64x2_t rre = crealdd(z);
      float64x2_t rim = cimagdd(z);
      REQUIRE(rre.hi == z.re.hi && rre.lo == z.re.lo);
      REQUIRE(rim.hi == z.im.hi && rim.lo == z.im.lo);
    }
    // conj: (re, -im) bitwise (signs of both limbs flipped).
    {
      complex64x2_t c2 = conjdd(z);
      REQUIRE(c2.re.hi == z.re.hi && c2.re.lo == z.re.lo);
      // Signed zero: conj should flip +0 ↔ -0 on the imag part.
      REQUIRE(to_cq(c2.im) == -to_cq(z.im) ||
              (z.im.hi == 0.0 && z.im.lo == 0.0));
    }
  }

  // cproj / conj edge cases with non-finite inputs.
  double const pinf = std::numeric_limits<double>::infinity();
  double const nan = std::numeric_limits<double>::quiet_NaN();
  struct InfCase {
    double re;
    double im;
    double expect_re;
    double expect_im;  // cproj result (sign of im preserved).
  };
  InfCase const inf_cases[] = {
      {pinf, 1.0, pinf, +0.0},
      {1.0, pinf, pinf, +0.0},
      {-pinf, -2.0, pinf, -0.0},
      {1.0, -pinf, pinf, -0.0},
      {nan, pinf, pinf, +0.0},  // inf on any part collapses.
      {1.0, 2.0, 1.0, 2.0},     // finite: identity.
      {nan, nan, nan, nan},     // all-nan: identity (no inf).
  };
  for (InfCase ic : inf_cases) {
    complex64x2_t z = {{ic.re, 0.0}, {ic.im, 0.0}};
    complex64x2_t r = cprojdd(z);
    if (std::isnan(ic.expect_re)) {
      REQUIRE(std::isnan(r.re.hi));
    } else {
      REQUIRE(r.re.hi == ic.expect_re);
      // Signed-zero on imag is the whole point of cproj: verify bit-exact.
      if (ic.expect_im == 0.0) {
        REQUIRE(r.im.hi == 0.0);
        REQUIRE(std::signbit(r.im.hi) == std::signbit(ic.expect_im));
      } else if (std::isnan(ic.expect_im)) {
        REQUIRE(std::isnan(r.im.hi));
      } else {
        REQUIRE(r.im.hi == ic.expect_im);
      }
    }
  }
}

// =============================================================================
// Complex DD transcendentals added in the 1.0 API completion pass:
//   cexpm1, clog2, clog1p, csinpi, ccospi.
// (clog10 already existed.) We compare against libquadmath references:
//   cexpm1(z) ↔ cexpq(z) − 1
//   clog2(z)  ↔ clogq(z)/log(2)
//   clog1p(z) ↔ clogq(1 + z)
//   csinpi(z) ↔ csinq(π·z), with exact-integer z: exactly (0, 0).
//   ccospi(z) ↔ ccosq(π·z), with exact-integer z: exactly (±1, 0).
static void test_complex_new_transcendentals(Stats &stats) {
  using cq_t = __complex128;
  auto to_cmf = [](q_t re, q_t im) {
    mf::float64x2 r = from_q(re), i = from_q(im);
    complex64x2_t z;
    z.re.hi = r._limbs[0]; z.re.lo = r._limbs[1];
    z.im.hi = i._limbs[0]; z.im.lo = i._limbs[1];
    return z;
  };
  auto check_c = [&stats](char const *name, complex64x2_t got, cq_t ref,
                          q_t tol) {
    q_t re_got = (q_t)got.re.hi + (q_t)got.re.lo;
    q_t im_got = (q_t)got.im.hi + (q_t)got.im.lo;
    q_t ref_mag = cabsq(ref);
    q_t scale = ref_mag > 1 ? ref_mag : (q_t)1;
    q_t err = fabsq(re_got - crealq(ref)) + fabsq(im_got - cimagq(ref));
    q_t rel = err / scale;
    if (rel > tol) {
      std::printf("FAIL %s rel=%g tol=%g\n", name, (double)rel, (double)tol);
      ++g_failures;
    }
    ++g_checks;
    stats.update(rel);
  };

  q_t const pi_q = M_PIq;

  // Small-|z| checkpoints where cexpm1 / clog1p should give full DD
  // precision. Also a couple of larger magnitudes and the signed axes.
  struct Case { q_t re, im; };
  std::vector<Case> cases = {
      { (q_t)0,        (q_t)0        },
      { (q_t)1e-20,    (q_t)0        },
      { (q_t)0,        (q_t)1e-20    },
      { (q_t)1e-10,    (q_t)2e-10    },
      { (q_t)0.25,     (q_t)0.5      },
      { (q_t)-0.5,     (q_t)0.25     },
      { (q_t)1.5,      (q_t)-0.75    },
      { (q_t)-3.0,     (q_t)2.0      },
      { (q_t)0,        pi_q          },  // i·π: e^z = -1, so expm1 = -2.
  };

  for (Case c : cases) {
    complex64x2_t z = to_cmf(c.re, c.im);
    cq_t zq;
    __real__ zq = c.re;
    __imag__ zq = c.im;

    // cexpm1(z) vs cexpq(z) - 1. Tolerance 1e-28 — two fused sincos plus
    // half-angle, so ~2-3× the cexp worst case.
    complex64x2_t em1 = cexpm1dd(z);
    cq_t ref_em1 = cexpq(zq) - 1;
    check_c("cexpm1", em1, ref_em1, (q_t)1e-28);

    // clog1p(z) vs clogq(1 + z). Tolerance 1e-29 — log1p + atan2.
    complex64x2_t l1p = clog1pdd(z);
    cq_t one_plus_z = zq;
    __real__ one_plus_z = (q_t)1 + crealq(zq);
    cq_t ref_l1p = clogq(one_plus_z);
    // clog(0) is -inf; skip tolerance check there.
    q_t ref_l1p_re = crealq(ref_l1p);
    if (!isinfq(ref_l1p_re) && !isnanq(ref_l1p_re)) {
      check_c("clog1p", l1p, ref_l1p, (q_t)1e-29);
    }

    // clog2(z) vs clogq(z) / ln 2. Skip z=0 (log(0) = -inf).
    if (c.re != 0 || c.im != 0) {
      complex64x2_t l2 = clog2dd(z);
      cq_t ref_l2 = clogq(zq) / logq((q_t)2);
      check_c("clog2", l2, ref_l2, (q_t)1e-29);
    }
  }

  // Integer-argument sinpi/cospi: exact-zero imag for real integer z.
  // csinpi(n+0i) must be exactly (0, 0) (parity check only — real sinpi
  // already hits this; complex path must not add drift from π·b·sinhcosh).
  for (int n = -3; n <= 3; ++n) {
    complex64x2_t z = to_cmf((q_t)n, (q_t)0);
    complex64x2_t s = csinpidd(z);
    complex64x2_t c = ccospidd(z);
    // sin(n·π) = 0 exactly; cos(n·π) = (-1)^n.
    REQUIRE(s.re.hi == 0.0 && s.re.lo == 0.0);
    REQUIRE(s.im.hi == 0.0 && s.im.lo == 0.0);
    REQUIRE(c.re.hi == ((n & 1) ? -1.0 : 1.0));
    REQUIRE(c.re.lo == 0.0);
    REQUIRE(c.im.hi == 0.0 && c.im.lo == 0.0);
  }

  // Non-integer z: csinpi / ccospi vs csinq(π·z) / ccosq(π·z). The ref
  // multiplies in __float128 so it has ~30 bits of headroom.
  struct Case2 { q_t re, im; };
  // Note: huge real arguments (|Re z| >> 1) are deliberately excluded —
  // the __float128 reference computes csinq(M_PIq * z) which itself loses
  // bits when M_PIq * Re(z) overflows ~70 bits, drowning our precision
  // signal. csinpi's whole point is to avoid that on the DD side.
  std::vector<Case2> trig_cases = {
      { (q_t)0.25,  (q_t)0       },
      { (q_t)0.5,   (q_t)0.3     },
      { (q_t)10.75, (q_t)0       },
      { (q_t)-0.1,  (q_t)0.7     },
  };
  for (Case2 c : trig_cases) {
    complex64x2_t z = to_cmf(c.re, c.im);
    cq_t piz;
    __real__ piz = pi_q * c.re;
    __imag__ piz = pi_q * c.im;
    complex64x2_t s = csinpidd(z);
    complex64x2_t cc = ccospidd(z);
    check_c("csinpi", s, csinq(piz), (q_t)1e-28);
    check_c("ccospi", cc, ccosq(piz), (q_t)1e-28);
  }
}

// =============================================================================
// atan cutover edge cases
// =============================================================================
//
// The large-argument branch uses `t = -1/|x|` once `|x| >= 32/3`, which puts
// worst-case `|t| = 3/32 = 0.09375` at the fitted edge of the rational fit.
// The table branch spans k = 0..85 so the newly-reached cells (k = 83, 84, 85)
// and all crossings between them need to be walked.
static void test_atan_cutover(Stats &stats) {
  // Enumerate exact boundary values (k/8 centers where t = 0 exactly),
  // midpoints between them (worst table-branch |t|), and the new cutover.
  double const boundaries[] = {
      // Old cutover, now in table branch (k = 82 center).
      10.25,
      // New table entries: k/8 centers.
      10.375,                 // k = 83
      10.5,                   // k = 84
      10.625,                 // k = 85
      // Table-branch midpoints (maximum |t| within a cell).
      10.3125,                // between k=82 and k=83
      10.4375,                // between k=83 and k=84
      10.5625,                // between k=84 and k=85
      // New cutover: |x| = 32/3 -> large branch, |t| = 3/32 = 0.09375 exactly.
      32.0 / 3.0,
      // Large-branch values immediately beyond the cutover.
      10.7, 11.0, 12.0, 15.0, 100.0, 1e10, 1e100, 1e300,
      // Also pick up a value in the old 10.25..10.6667 gap that used to go
      // through the large branch at |t| > 0.09375.
      10.26, 10.3, 10.4, 10.5, 10.6, 10.65,
  };

  auto check_atan = [&](q_t qx) {
    mf::float64x2 x = from_q(qx);
    mf::float64x2 got = mf::atan(x);
    q_t expected = atanq(qx);
    // atan is always well-conditioned (|d/dx atan| <= 1), so tolerance is
    // the standard DD kernel tier.
    check_q("atan_cutover", got, expected, 1e-30, stats, &x);
  };

  for (double b : boundaries) {
    // The point itself, its bit-neighbors, and a ±16 ulp ring on the DD-side.
    q_t qb = (q_t)b;
    check_atan(qb);
    check_atan(-qb);
    q_t up = qb, dn = qb;
    for (int i = 0; i < 16; ++i) {
      up = nextafterq(up, (q_t)std::numeric_limits<double>::infinity());
      dn = nextafterq(dn, (q_t)-std::numeric_limits<double>::infinity());
      check_atan(up);
      check_atan(dn);
      check_atan(-up);
      check_atan(-dn);
    }
  }

  // Dense sweep across the newly-reached table range, with nontrivial lo
  // limbs so the kernel sees full-DD arguments (not just exact doubles).
  for (int i = 0; i < 2000; ++i) {
    // Uniform random in [10.25, 32/3).
    double frac = (double)i / 2000.0;
    q_t qx = (q_t)10.25 + ((q_t)(32.0 / 3.0) - (q_t)10.25) * (q_t)frac;
    // Add a sub-ulp perturbation so the lo limb is exercised.
    qx += scalbnq((q_t)1, -60) * (q_t)((i % 7) - 3);
    check_atan(qx);
    check_atan(-qx);
  }

  // Verify that the cutover value lands in the large branch: round-trip
  // identity arctan(x) + arctan(1/x) = pi/2 for x > 0 at |x| = 32/3.
  mf::float64x2 x = from_q((q_t)(32.0 / 3.0));
  mf::float64x2 one(1.0);
  mf::float64x2 y = one / x;
  mf::float64x2 sum = mf::atan(x) + mf::atan(y);
  q_t got = to_q(sum);
  q_t expected = M_PIq / (q_t)2;
  q_t diff = got - expected;
  if (diff < 0) diff = -diff;
  q_t mag = expected < 0 ? -expected : expected;
  double rel = (double)(diff / mag);
  ++g_checks;
  stats.update(rel);
  if (!(rel <= 1e-30)) {
    ++g_failures;
    std::fprintf(stderr,
                 "FAIL [atan+arccotan=pi/2 at 32/3] rel_err=%g\n", rel);
  }
}

// =============================================================================
// Triple-double scratch primitives (`detail::float64x3`, `td_mul_td`, …)
// =============================================================================
//
// Each primitive is verified against a __float128 reference. Inputs are
// chosen so the qp-exact operation fits inside qp's 113-bit significand
// — under that constraint, sum-preserving primitives (three_sum,
// td_add_double, td_add_dd, …) must match bit-exactly; td_mul_td is
// allowed qp-ulp rounding on the final-limb absorb. See
// doc/developer/TRIPLE_DOUBLE.md §3.

// Exact DD → qp: (q_t)(h) + (q_t)(m) + (q_t)(l) is exact when the limb
// span stays below 113 bits. Normalized TDs from our producers span at
// most ~159 bits, so for the carefully-chosen inputs below the sum fits.
static q_t to_q_td(mf::detail::float64x3 const &t) {
  return (q_t)t._limbs[0] + (q_t)t._limbs[1] + (q_t)t._limbs[2];
}

// Build a TD that represents `v` to ~159 bits by extracting residues.
static mf::detail::float64x3 td_from_q(q_t v) {
  double hi = (double)v;
  if (!std::isfinite(hi)) return {hi, 0.0, 0.0};
  double mid = (double)(v - (q_t)hi);
  double lo  = (double)(v - (q_t)hi - (q_t)mid);
  return {hi, mid, lo};
}

static void test_td_primitives(Stats &td_stats) {
  using mf::detail::float64x3;
  using mf::detail::three_sum;
  using mf::detail::renorm3;
  using mf::detail::td_add_double;
  using mf::detail::td_sub_double;
  using mf::detail::td_add_dd;
  using mf::detail::td_mul_td;
  using mf::detail::td_to_dd;
  using mf::detail::td_from_dd;

  auto expect_exact = [&](char const *name, q_t got, q_t expected) {
    ++g_checks;
    double rel = q_rel_err(got, expected);
    td_stats.update(rel);
    if (got != expected) {
      ++g_failures;
      std::fprintf(stderr,
                   "FAIL [%s] not exact: got=%s expected=%s (rel=%g)\n",
                   name, qstr(got), qstr(expected), rel);
    }
  };
  auto expect_near = [&](char const *name, q_t got, q_t expected, double tol) {
    ++g_checks;
    double rel = q_rel_err(got, expected);
    td_stats.update(rel);
    if (!(rel <= tol)) {
      ++g_failures;
      std::fprintf(stderr,
                   "FAIL [%s] rel=%g > tol=%g  got=%s expected=%s\n",
                   name, rel, tol, qstr(got), qstr(expected));
    }
  };

  // three_sum: widely separated magnitudes — each is exactly in qp, sum is exact.
  {
    double a = 1.0, b = 0x1p-60, c = 0x1p-120;
    double s, t, u;
    three_sum(a, b, c, s, t, u);
    expect_exact("three_sum wide",
                 (q_t)s + (q_t)t + (q_t)u,
                 (q_t)a + (q_t)b + (q_t)c);
  }
  // three_sum: cancellation — a + b = 0 exactly, c dominates.
  {
    double a = 1.0, b = -1.0, c = 0x1p-30;
    double s, t, u;
    three_sum(a, b, c, s, t, u);
    expect_exact("three_sum cancel",
                 (q_t)s + (q_t)t + (q_t)u,
                 (q_t)a + (q_t)b + (q_t)c);
  }
  // three_sum: reverse order (unordered inputs).
  {
    double a = 0x1p-120, b = 0x1p-60, c = 1.0;
    double s, t, u;
    three_sum(a, b, c, s, t, u);
    expect_exact("three_sum reversed",
                 (q_t)s + (q_t)t + (q_t)u,
                 (q_t)a + (q_t)b + (q_t)c);
  }

  // renorm3: sum-preserving and canonically ordered.
  {
    double h = 1.0, m = 0x1p-53, l = 0x1p-106;
    float64x3 r = renorm3(h, m, l);
    expect_exact("renorm3 sum",
                 (q_t)r._limbs[0] + (q_t)r._limbs[1] + (q_t)r._limbs[2],
                 (q_t)h + (q_t)m + (q_t)l);
    // Canonical ordering: |m| ≤ ulp(h), |l| ≤ ulp(m) (within 1 bit).
    REQUIRE(std::abs(r._limbs[1]) <=
            std::ldexp(std::abs(r._limbs[0]), -52));
    REQUIRE(std::abs(r._limbs[2]) <=
            std::ldexp(std::abs(r._limbs[1]), -52) ||
            r._limbs[1] == 0.0);
  }

  // td_add_double: sum preservation for a TD + exact double.
  {
    float64x3 td{1.0, 0x1p-53, 0x1p-106};
    double d = 0.5;
    float64x3 r = td_add_double(td, d);
    expect_exact("td_add_double",
                 (q_t)r._limbs[0] + (q_t)r._limbs[1] + (q_t)r._limbs[2],
                 (q_t)td._limbs[0] + (q_t)td._limbs[1] +
                     (q_t)td._limbs[2] + (q_t)d);
  }
  // td_sub_double: the cancellation case that motivates TD internals —
  //   cos(b)·e^a − 1 near the cancellation surface. Here we just verify
  //   the subtraction is exact (qp will not round a 107-bit operand).
  {
    // td representing 1 + 2^-80.
    float64x3 td{1.0, 0x1p-80, 0.0};
    float64x3 r = td_sub_double(td, 1.0);
    expect_exact("td_sub_double cancel",
                 (q_t)r._limbs[0] + (q_t)r._limbs[1] + (q_t)r._limbs[2],
                 (q_t)0x1p-80);
  }

  // td_add_dd: TD + DD with a non-trivial DD lo limb.
  {
    float64x3 td{1.0, 0x1p-53, 0x1p-106};
    mf::float64x2 dd(0x1p-30, 0x1p-83);
    float64x3 r = td_add_dd(td, dd);
    expect_exact("td_add_dd",
                 (q_t)r._limbs[0] + (q_t)r._limbs[1] + (q_t)r._limbs[2],
                 (q_t)td._limbs[0] + (q_t)td._limbs[1] +
                     (q_t)td._limbs[2] +
                     (q_t)dd._limbs[0] + (q_t)dd._limbs[1]);
  }

  // td_to_dd: round back to DD; the DD pair matches the TD sum to the
  //   DD resolution (~2^-105). A third limb below that threshold is
  //   absorbed into the leading pair with at most one DD-ulp of rounding.
  {
    float64x3 td{1.0, 0x1p-53, 0x1p-106};
    mf::float64x2 dd = td_to_dd(td);
    q_t got = to_q(dd);
    q_t expected = (q_t)1.0 + (q_t)0x1p-53 + (q_t)0x1p-106;
    expect_near("td_to_dd drop-3rd", got, expected, 1.3e-32);  // 1 DD ulp
  }
  // td_to_dd with a third limb that survives DD rounding: input
  //   (1, 2^-60, 2^-110) — the 2^-110 bit is captured exactly into
  //   the DD low limb through the final fold.
  {
    float64x3 td{1.0, 0x1p-60, 0x1p-110};
    mf::float64x2 dd = td_to_dd(td);
    q_t got = to_q(dd);
    q_t expected = (q_t)1.0 + (q_t)0x1p-60 + (q_t)0x1p-110;
    expect_exact("td_to_dd fold", got, expected);
  }

  // td_mul_td: two TDs representing irrational values (π and e-ish).
  //   qp can represent the 113-bit product; TD mul should match qp to
  //   ~2^-110 (qp rounding), well inside the qp precision band.
  {
    q_t va = M_PIq;
    q_t vb = M_Eq;
    float64x3 a = td_from_q(va);
    float64x3 b = td_from_q(vb);
    float64x3 p = td_mul_td(a, b);
    // qp product itself rounds at 2^-112; TD carries ~2^-159 of extra
    //   precision below that, so the match is limited by qp.
    expect_near("td_mul_td pi·e", to_q_td(p), va * vb, 2e-33);
  }
  // td_mul_td: near-cancellation product — a·b very close to 1. Inputs
  //   are chosen so the qp-exact product fits inside qp's 113-bit
  //   significand (residue at 2^-101, leading 1 bit and trailing sit
  //   within ~101 bits, well under the qp limit), which pins qp as the
  //   gold reference against the TD result.
  {
    q_t one = (q_t)1;
    q_t va = (q_t)2 + scalbnq(one, -50);
    q_t vb = (q_t)0.5 - scalbnq(one, -51);
    float64x3 a = td_from_q(va);
    float64x3 b = td_from_q(vb);
    float64x3 p = td_mul_td(a, b);
    expect_near("td_mul_td near-1", to_q_td(p), va * vb, 2e-33);

    // And the motivating cexpm1-style composition: (a*b) - 1 near the
    //   cancellation surface. TD subtract exact; qp reference rounds at
    //   2^-112, which is still well below TD's 2^-159.
    float64x3 diff = td_sub_double(p, 1.0);
    q_t diff_ref = va * vb - one;
    expect_near("td_mul·sub cancel", to_q_td(diff), diff_ref, 2e-33);
    // Key behavioural check: the DD projection of (TD·TD − 1) preserves
    //   the residue, whereas a naive DD mul → (DD·DD − 1) would have
    //   lost it entirely.
    mf::float64x2 diff_dd = td_to_dd(diff);
    REQUIRE(diff_dd._limbs[0] != 0.0 || diff_dd._limbs[1] != 0.0);
  }

  // td_from_dd round-trip: DD → TD → DD must be identity on any
  //   normalized DD.
  {
    mf::float64x2 dd = from_q(M_PIq);
    float64x3 td = td_from_dd(dd);
    mf::float64x2 back = td_to_dd(td);
    expect_exact("td_from_dd → td_to_dd",
                 to_q(back), to_q(dd));
  }
}

// =============================================================================
// Driver
// =============================================================================

static void print_stats(char const *name, Stats const &s) {
  if (s.count == 0) {
    std::printf("  %-12s n=%d\n", name, s.count);
    return;
  }
  std::printf("  %-12s n=%-6d max_rel=%.3e  mean_rel=%.3e\n", name, s.count,
              s.max_rel, s.sum_rel / s.count);
}

int main() {
  test_construction_and_conversion();
  test_equality_and_ordering();
  test_classification();
  test_nextafter_symmetry();

  Stats add, sub, mul, div, una, abs_fmm, csgn, ldx, atn, lrp, cxa, bjn, cxn, lre, htr, cbr, tdp;
  test_unary(una);
  test_addition(add);
  test_subtraction(sub);
  test_multiplication(mul);
  test_division(div);
  test_division_nonfinite();
  test_atan2_signed_zero();
  test_csqrt_zero_branch();
  test_complex_mul_cancellation();
  test_io_to_string_and_stream();
  test_abs_fmin_fmax(abs_fmm);
  test_copysign(csgn);
  test_ldexp_scalbn_ilogb(ldx);
  test_atan_cutover(atn);
  test_lerp(lrp);
  test_complex_accessors(cxa);
  test_complex_new_transcendentals(cxn);
  test_log_root_edges(lre);
  test_huge_argument_trig(htr);
  test_complex_branch_cuts(cbr);
  test_bessel_jn_miller_precision(bjn);
  test_td_primitives(tdp);

  std::printf("[multifloats_test] %d checks, %d failures\n", g_checks,
              g_failures);
  std::printf("[multifloats_test] precision report (relative error vs "
              "__float128):\n");
  print_stats("unary", una);
  print_stats("add", add);
  print_stats("sub", sub);
  print_stats("mul", mul);
  print_stats("div", div);
  print_stats("abs/fmin/fmax", abs_fmm);
  print_stats("copysign", csgn);
  print_stats("ldexp/etc", ldx);
  print_stats("atan cutover", atn);
  print_stats("lerp", lrp);
  print_stats("cx accessors", cxa);
  print_stats("cx new transc", cxn);
  print_stats("log/root edges", lre);
  print_stats("huge arg trig", htr);
  print_stats("cx branch cuts", cbr);
  print_stats("bessel jn (Miller)", bjn);
  print_stats("td primitives", tdp);

  // ----------------------------------------------------------------
  // Tolerance sensitivity ratchet (see §7 of doc/developer/INTERNALS.md).
  // ----------------------------------------------------------------
  // Pins observed max_rel for each category between [1/20×, 20×] of a
  // recorded expected value. Either direction fails:
  //   - observed > 20× expected → silent precision regression.
  //   - observed < 1/20× expected → the kernel silently got better; the
  //     pin should be tightened so the improvement is locked in.
  // The 20× band tolerates Estrin-style rounding variation (see Tier 3
  // sin/cos going 3.4e-32 → 5.0e-32) while still catching order-of-
  // magnitude drift in either direction.
  struct Pin { char const *name; Stats const *s; double expect_max; };
  Pin pins[] = {
      {"add",           &add,     1.1e-32},
      {"sub",           &sub,     1.3e-32},
      {"mul",           &mul,     2.1e-32},
      {"div",           &div,     3.4e-32},
      {"atan cutover",  &atn,     2.0e-32},
      {"cx accessors",  &cxa,     1.0e-32},
      {"cx new transc", &cxn,     1.2e-31},
      {"log/root edges",&lre,     5.6e-32},
      {"huge arg trig", &htr,     6.0e-33},
      {"cx branch cuts",&cbr,     1.5e-31},
      {"bessel Miller", &bjn,     1.0e-31},
  };
  std::printf("[multifloats_test] tolerance-ratchet (expected ± 20×):\n");
  int ratchet_fails = 0;
  for (Pin const &p : pins) {
    if (p.s->count == 0) continue;
    double ratio = p.s->max_rel / p.expect_max;
    bool too_high = ratio > 20.0;
    bool too_low  = ratio < 1.0 / 20.0 && p.s->max_rel > 0;
    char const *status = too_high ? "WORSE" : too_low ? "BETTER" : "ok";
    std::printf("  %-14s observed=%.3e expected=%.3e ratio=%.2f  [%s]\n",
                p.name, p.s->max_rel, p.expect_max, ratio, status);
    if (too_high || too_low) ++ratchet_fails;
  }
  if (ratchet_fails) {
    std::printf("[multifloats_test] %d tolerance-ratchet failures "
                "(update pins in test.cc if the kernel changed on purpose)\n",
                ratchet_fails);
    return 1;
  }

  return g_failures == 0 ? 0 : 1;
}
