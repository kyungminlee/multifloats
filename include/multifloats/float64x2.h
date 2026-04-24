/* multifloats — double-double arithmetic, unified C and C++ header.
 *
 * Structure:
 *   - A C-ABI section defining `float64x2` / `complex64x2` and every
 *     `extern "C"` `*dd` entry point. Valid in a plain-C translation unit
 *     (include from C, C++, or Fortran via iso_c_binding).
 *   - A C++-only section (guarded by `#ifdef __cplusplus`) providing the
 *     header-inline `multifloats::float64x2` class, <cmath>-style free
 *     functions, and the `std::complex<float64x2>` specializations.
 *
 * All DD math functions are defined as extern "C" in multifloats_math.cc.
 * They follow the libc / libquadmath convention: type suffix on the
 * function name (sindd, logdd, cexpdd) — mirroring `sinf`/`sinq`, `cexpf`/
 * `csinq`.
 */
#pragma once

/* ============================================================================
 * C ABI section — valid in a C or C++ translation unit.
 * ========================================================================= */

#include <stdint.h>
#include <stddef.h>  /* size_t */

/* ABI version. Bump on any breaking change to the `float64x2` layout,
 * the argument/return convention of any exported function, or removal of
 * an exported symbol. Additive changes (new *dd functions) keep the same
 * version. Callers can gate on this at compile time:
 *     #if !defined(MULTIFLOATS_ABI_VERSION) || MULTIFLOATS_ABI_VERSION < 2
 *     #error "multifloats 2.x required"
 *     #endif
 */
#define MULTIFLOATS_ABI_VERSION 2

/* Visibility attribute for every exported function. */
#if defined(__GNUC__) || defined(__clang__)
#  define MULTIFLOATS_API __attribute__((visibility("default")))
#else
#  define MULTIFLOATS_API
#endif

#ifdef __cplusplus

// Every error-free transformation here (two_prod, reduce_pi_half, erfc x-split,
// matmul mac_inl, ...) depends on std::fma being IEEE-compliant: a single
// multiply-add with one rounding. C99/C++ guarantees this for conforming
// implementations — including the software fallback when the target lacks
// a hardware FMA instruction. The one common way to break that guarantee is
// -ffast-math, which lets the compiler rewrite fma into a non-rounded sequence
// or fold it with surrounding expressions. Refuse to build in that mode: the
// DD arithmetic would silently collapse to plain double precision.
#if defined(__FAST_MATH__)
#  error "multifloats cannot be built with -ffast-math: it breaks the IEEE " \
         "semantics of std::fma that every error-free transformation relies " \
         "on. Drop -ffast-math (and -funsafe-math-optimizations) and rebuild."
#endif

#include <charconv>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iosfwd>
#include <limits>
#include <string>
#include <system_error>
#include <type_traits>

namespace multifloats {

namespace detail {

// =============================================================================
// Error-free transformations (inputs by value to avoid aliasing pitfalls).
// Kept T-generic: they're called with plain double from both float64x2 ops
// and the triple-double float64x3 internals, and the generic form documents
// them as EFTs rather than float64x2-specific helpers.
// =============================================================================

template <typename T>
constexpr void two_sum(T a, T b, T &sum, T &err) {
  T s = a + b;
  T a_prime = s - b;
  T b_prime = s - a_prime;
  err = (a - a_prime) + (b - b_prime);
  sum = s;
}

template <typename T>
constexpr void fast_two_sum(T a, T b, T &sum, T &err) {
  T s = a + b;
  T b_prime = s - a;
  err = b - b_prime;
  sum = s;
}

template <typename T> constexpr T one_prod(T a, T b) { return a * b; }

template <typename T>
constexpr void two_prod(T a, T b, T &prod, T &err) {
  T p = a * b;
  err = std::fma(a, b, -p);
  prod = p;
}

} // namespace detail
#endif /* __cplusplus */

/* =============================================================================
 * float64x2 — double-double (DD) POD and C ABI interchange type.
 *
 * Two IEEE-754 binary64 limbs. A canonical value satisfies
 * |limbs[1]| <= ulp(limbs[0])/2, giving ~104 bits of significand. Binary ops
 * are inlined as C++ constexpr members; transcendentals delegate to the
 * extern "C" `*dd` kernels defined in multifloats_math.cc.
 *
 * Layout is identical under C and C++; the Fortran `bind(c)` interop layer
 * assumes two back-to-back doubles (four for complex64x2).
 * ============================================================================= */
struct float64x2 {
  // The `= {}` member initializer (C++ only) is load-bearing: operator+,
  // operator-, operator*, and operator/ all have non-finite short-circuit
  // branches that write only `out.limbs[0]` and rely on `limbs[1]` being
  // zero-initialized. Dropping `= {}` to make the type trivially-default-
  // constructible would silently leak indeterminate bits through those
  // branches. If you ever want that, add explicit `out.limbs[1] = 0.0;`
  // lines to each short-circuit first.
  double limbs[2]
#ifdef __cplusplus
      = {}
#endif
      ;

#ifdef __cplusplus
  constexpr float64x2() = default;
  constexpr float64x2(float64x2 const &) = default;
  constexpr float64x2(float64x2 &&) = default;
  constexpr float64x2 &operator=(float64x2 const &) = default;
  constexpr float64x2 &operator=(float64x2 &&) = default;

  // One- and two-argument ctors in member-init form for symmetry. Non-explicit:
  // the scalar form is used implicitly throughout (e.g. `float64x2(1)` in
  // hypot/cbrt), and the two-arg form is the natural DD-literal spelling
  // `float64x2{hi, lo}` used inside class-method bodies.
  constexpr float64x2(double arg) : limbs{arg, 0.0} {}
  constexpr float64x2(double hi, double lo) : limbs{hi, lo} {}

  constexpr explicit operator double() const { return limbs[0]; }

  // Lexicographic comparison. NaN is unordered: IEEE `<` is false on either
  // side, so for NaN limbs both `<` checks fall through to the next limb,
  // yielding `operator<` = false — same as the previous _lex_compare loop.
  // For +0 vs -0 the leading comparisons are both false (IEEE +0 == -0), so
  // we correctly fall through to the next limb.
  constexpr bool operator==(float64x2 const &r) const {
    return limbs[0] == r.limbs[0] && limbs[1] == r.limbs[1];
  }
  constexpr bool operator!=(float64x2 const &r) const { return !(*this == r); }
  constexpr bool operator<(float64x2 const &r) const {
    if (limbs[0] < r.limbs[0]) return true;
    if (r.limbs[0] < limbs[0]) return false;
    return limbs[1] < r.limbs[1];
  }
  constexpr bool operator>(float64x2 const &r) const  { return  (r < *this); }
  constexpr bool operator<=(float64x2 const &r) const { return !(r < *this); }
  constexpr bool operator>=(float64x2 const &r) const { return !(*this < r); }

  constexpr float64x2 operator+() const { return *this; }
  constexpr float64x2 operator-() const { return {-limbs[0], -limbs[1]}; }

  // ---------------------------------------------------------------------------
  // Binary arithmetic — kernels inlined directly, translated from
  // MultiFloats.jl (mfadd / mfmul) and the Float64x2 division kernel.
  // ---------------------------------------------------------------------------

  constexpr float64x2 operator+(float64x2 const &rhs) const {
    float64x2 out;
    double s = limbs[0] + rhs.limbs[0];
    // Non-finite: the EFT below would propagate NaN into limbs[1];
    // short-circuit and let IEEE produce the correct leading limb.
    if (!std::isfinite(s)) {
      out.limbs[0] = s;
      return out;
    }
    // When both hi limbs are zero, two_sum loses the -0 sign
    // (IEEE 754: -0 + +0 = +0 in round-to-nearest).
    if (limbs[0] == 0.0 && rhs.limbs[0] == 0.0) {
      out.limbs[0] = s;
      out.limbs[1] = limbs[1] + rhs.limbs[1];
      return out;
    }
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
    detail::two_sum(limbs[0], rhs.limbs[0], a, b);
    detail::two_sum(limbs[1], rhs.limbs[1], c, d);
    detail::fast_two_sum(a, c, a, c);
    b += d;
    b += c;
    detail::fast_two_sum(a, b, out.limbs[0], out.limbs[1]);
    return out;
  }

  constexpr float64x2 operator-(float64x2 const &rhs) const {
    float64x2 out;
    double s = limbs[0] - rhs.limbs[0];
    if (!std::isfinite(s)) {
      out.limbs[0] = s;
      return out;
    }
    if (limbs[0] == 0.0 && rhs.limbs[0] == 0.0) {
      out.limbs[0] = s;
      out.limbs[1] = limbs[1] - rhs.limbs[1];
      return out;
    }
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
    detail::two_sum(limbs[0], -rhs.limbs[0], a, b);
    detail::two_sum(limbs[1], -rhs.limbs[1], c, d);
    detail::fast_two_sum(a, c, a, c);
    b += d;
    b += c;
    detail::fast_two_sum(a, b, out.limbs[0], out.limbs[1]);
    return out;
  }

  constexpr float64x2 operator*(float64x2 const &rhs) const {
    float64x2 out;
    double p00 = 0.0, e00 = 0.0;
    detail::two_prod(limbs[0], rhs.limbs[0], p00, e00);
    // Non-finite: two_prod's fma residual is NaN whenever p00 overflows
    // or either input is ±Inf/NaN (e.g. inf*2 → p00=inf, e00=fma(inf,
    // 2,-inf)=NaN). Short-circuit so the EFT residual doesn't contaminate
    // out.limbs[1]; IEEE supplies the correct leading limb in p00.
    if (!std::isfinite(p00)) {
      out.limbs[0] = p00;
      return out;
    }
    double p01 = detail::one_prod(limbs[0], rhs.limbs[1]);
    double p10 = detail::one_prod(limbs[1], rhs.limbs[0]);
    p01 += p10;
    e00 += p01;
    detail::fast_two_sum(p00, e00, out.limbs[0], out.limbs[1]);
    return out;
  }

  constexpr float64x2 operator/(float64x2 const &rhs) const {
    // Dekker-style: q1 = hi/rhs.hi, refine once.
    double q1 = limbs[0] / rhs.limbs[0];
    if (!std::isfinite(q1)) {
      // Mirror q1 into the lo limb so a non-finite result propagates
      // through both limbs. Otherwise isnan/isinf checks against the lo
      // limb would spuriously report "finite" on a NaN/Inf DD.
      float64x2 out;
      out.limbs[0] = q1;
      out.limbs[1] = q1;
      return out;
    }
    if (!std::isfinite(rhs.limbs[0])) {
      // Finite / ±Inf — q1 is ±0; the correct DD is {±0, 0}, which
      // default-initialization already gives us.
      float64x2 out;
      out.limbs[0] = q1;
      return out;
    }
    // r = this - q1 * rhs, computed as a full DD (q1 is a single-limb
    // scalar so q1*rhs is one two_prod + one one_prod = one DD).
    double p00 = 0.0, e00 = 0.0;
    detail::two_prod(q1, rhs.limbs[0], p00, e00);
    double p01 = detail::one_prod(q1, rhs.limbs[1]);
    double qhi = p00;
    double qlo = e00 + p01;
    // r = this - (qhi, qlo) via two_diff
    double r0 = 0.0, r0e = 0.0;
    detail::two_sum(limbs[0], -qhi, r0, r0e);
    double r1 = (limbs[1] - qlo) + r0e;
    double rh = 0.0, rl = 0.0;
    detail::fast_two_sum(r0, r1, rh, rl);
    // q2 = r.hi / rhs.hi
    double q2 = rh / rhs.limbs[0];
    float64x2 out;
    detail::fast_two_sum(q1, q2, out.limbs[0], out.limbs[1]);
    return out;
  }

  constexpr float64x2 &operator+=(float64x2 const &rhs) {
    return *this = *this + rhs;
  }
  constexpr float64x2 &operator-=(float64x2 const &rhs) {
    return *this = *this - rhs;
  }
  constexpr float64x2 &operator*=(float64x2 const &rhs) {
    return *this = *this * rhs;
  }
  constexpr float64x2 &operator/=(float64x2 const &rhs) {
    return *this = *this / rhs;
  }
#endif /* __cplusplus */
};

/* complex64x2 — double-double complex.
 *   • C:   a plain POD (two back-to-back float64x2s, four doubles).
 *   • C++: a type alias for std::complex<float64x2>. C++11 §26.4/4 only
 *          guarantees the `T[2]` layout for T ∈ {float, double, long double};
 *          for user-defined T the ABI match rests on implementation detail.
 *          libstdc++ and libc++ both store two `_M_value[2]` / `_M_real` +
 *          `_M_imag` members laid out as `T[2]`, and the static_asserts
 *          below pin sizeof to 4*sizeof(double). If a future standard
 *          library changes that, the asserts fire at compile time. */
#ifdef __cplusplus
using complex64x2 = std::complex<float64x2>;
#else
struct complex64x2 { struct float64x2 re, im; };
typedef struct float64x2 float64x2;
typedef struct complex64x2 complex64x2;
#endif

/* Guard against surprise padding — the entire C ABI and the Fortran
 * iso_c_binding layer assume float64x2 is exactly two back-to-back
 * doubles (and complex64x2 exactly four). C11 / C++11 required. */
#ifdef __cplusplus
static_assert(sizeof(float64x2) == 2 * sizeof(double),
              "float64x2 must be two back-to-back doubles with no padding");
static_assert(sizeof(complex64x2) == 4 * sizeof(double),
              "std::complex<float64x2> must be four back-to-back doubles with no padding");
#else
_Static_assert(sizeof(struct float64x2) == 2 * sizeof(double),
               "float64x2 must be two back-to-back doubles with no padding");
_Static_assert(sizeof(struct complex64x2) == 4 * sizeof(double),
               "complex64x2 must be four back-to-back doubles with no padding");
#endif

#ifdef __cplusplus
/* float64x2 is standard-layout + trivially-copyable, so its System V AMD64 /
 * ARM AAPCS calling convention matches a plain `struct { double, double }`.
 * Clang's -Wreturn-type-c-linkage conservatively warns because the struct has
 * user-provided constructors (non-POD in the C++98 sense) — the ABI is still
 * fine, silence the warning over the extern "C" block. */
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
#endif
extern "C" {
#endif

/* Arithmetic */
MULTIFLOATS_API float64x2 adddd(float64x2 a, float64x2 b);
MULTIFLOATS_API float64x2 subdd(float64x2 a, float64x2 b);
MULTIFLOATS_API float64x2 muldd(float64x2 a, float64x2 b);
MULTIFLOATS_API float64x2 divdd(float64x2 a, float64x2 b);

/* Unary */
MULTIFLOATS_API float64x2 negdd(float64x2 a);
MULTIFLOATS_API float64x2 fabsdd(float64x2 a);
MULTIFLOATS_API float64x2 sqrtdd(float64x2 a);

/* Rounding. truncdd matches C trunc / Fortran AINT (toward zero).
 * rounddd matches C round / Fortran ANINT (to nearest, halfway away). */
MULTIFLOATS_API float64x2 truncdd(float64x2 a);
MULTIFLOATS_API float64x2 rounddd(float64x2 a);

/* Binary */
MULTIFLOATS_API float64x2 fmindd(float64x2 a, float64x2 b);
MULTIFLOATS_API float64x2 fmaxdd(float64x2 a, float64x2 b);
MULTIFLOATS_API float64x2 hypotdd(float64x2 a, float64x2 b);
MULTIFLOATS_API float64x2 powdd(float64x2 a, float64x2 b);
/* Integer-exponent power via exponentiation-by-squaring: ⌈log2|n|⌉
 * squarings plus popcount(|n|) multiplies, exact up to DD precision
 * per step. One DD division at the end when n < 0. */
MULTIFLOATS_API float64x2 powidd(float64x2 a, int n);
MULTIFLOATS_API float64x2 fmoddd(float64x2 a, float64x2 b);
/* Floored modulo (Fortran `modulo` semantics): same sign as y.
 * modulodd(x, y) = fmod(x, y), plus y if the remainder and y have
 * opposite signs. Uses DD signbit (first-nonzero-limb), so
 * non-canonical DDs with hi==0 still get the right sign. */
MULTIFLOATS_API float64x2 modulodd(float64x2 a, float64x2 b);
MULTIFLOATS_API float64x2 fdimdd(float64x2 a, float64x2 b);
MULTIFLOATS_API float64x2 copysigndd(float64x2 a, float64x2 b);
MULTIFLOATS_API float64x2 fmadd(float64x2 a, float64x2 b, float64x2 c);

/* Exponential / logarithmic. `expm1dd` / `log1pdd` are the
 * cancellation-safe variants for arguments near zero — direct Taylor
 * / atanh-narrow kernels cover |x| below a threshold, larger |x|
 * falls through to the standard `exp` / `log` paths. */
MULTIFLOATS_API float64x2 expdd(float64x2 a);
MULTIFLOATS_API float64x2 exp2dd(float64x2 a);
MULTIFLOATS_API float64x2 expm1dd(float64x2 a);
MULTIFLOATS_API float64x2 logdd(float64x2 a);
MULTIFLOATS_API float64x2 log2dd(float64x2 a);
MULTIFLOATS_API float64x2 log10dd(float64x2 a);
MULTIFLOATS_API float64x2 log1pdd(float64x2 a);

/* Trigonometric */
MULTIFLOATS_API float64x2 sindd(float64x2 a);
MULTIFLOATS_API float64x2 cosdd(float64x2 a);
MULTIFLOATS_API float64x2 tandd(float64x2 a);
MULTIFLOATS_API float64x2 asindd(float64x2 a);
MULTIFLOATS_API float64x2 acosdd(float64x2 a);
MULTIFLOATS_API float64x2 atandd(float64x2 a);
MULTIFLOATS_API float64x2 atan2dd(float64x2 a, float64x2 b);

/* π-scaled trig: {sin,cos,tan}pidd(x) = fn(π·x),
 *                {asin,acos,atan}pidd(x) = fn(x)/π,
 *                atan2pidd(y,x) = atan2(y,x)/π. */
MULTIFLOATS_API float64x2 sinpidd(float64x2 a);
MULTIFLOATS_API float64x2 cospidd(float64x2 a);
MULTIFLOATS_API float64x2 tanpidd(float64x2 a);
MULTIFLOATS_API float64x2 asinpidd(float64x2 a);
MULTIFLOATS_API float64x2 acospidd(float64x2 a);
MULTIFLOATS_API float64x2 atanpidd(float64x2 a);
MULTIFLOATS_API float64x2 atan2pidd(float64x2 a, float64x2 b);

/* Hyperbolic */
MULTIFLOATS_API float64x2 sinhdd(float64x2 a);
MULTIFLOATS_API float64x2 coshdd(float64x2 a);
MULTIFLOATS_API float64x2 tanhdd(float64x2 a);
MULTIFLOATS_API float64x2 asinhdd(float64x2 a);
MULTIFLOATS_API float64x2 acoshdd(float64x2 a);
MULTIFLOATS_API float64x2 atanhdd(float64x2 a);

/* Error functions. erfcxdd is the scaled complementary error function
 *   erfcxdd(x) = exp(x^2) * erfcdd(x)
 * (standard name in Faddeeva/Julia/SciPy; Fortran calls it erfc_scaled). */
MULTIFLOATS_API float64x2 erfdd(float64x2 a);
MULTIFLOATS_API float64x2 erfcdd(float64x2 a);
MULTIFLOATS_API float64x2 erfcxdd(float64x2 a);

/* Gamma functions */
MULTIFLOATS_API float64x2 tgammadd(float64x2 a);
MULTIFLOATS_API float64x2 lgammadd(float64x2 a);

/* Bessel functions (POSIX naming). jndd / yndd: integer-order recurrence
 * seeded from j0dd / j1dd (Miller's backward recurrence for Jn when n > x).
 * Range variants fill `out[0..n2−n1]` with the order-n1 through order-n2
 * values; yndd_range uses a single stable forward-recurrence sweep (cheaper
 * than n2−n1+1 independent yndd calls), while jndd_range loops over jndd
 * because Jn's forward recurrence is unstable for n > x. */
MULTIFLOATS_API float64x2 j0dd(float64x2 a);
MULTIFLOATS_API float64x2 j1dd(float64x2 a);
MULTIFLOATS_API float64x2 y0dd(float64x2 a);
MULTIFLOATS_API float64x2 y1dd(float64x2 a);
MULTIFLOATS_API float64x2 jndd(int n, float64x2 a);
MULTIFLOATS_API float64x2 yndd(int n, float64x2 a);
MULTIFLOATS_API void jndd_range(int n1, int n2, float64x2 a, float64x2 *out);
MULTIFLOATS_API void yndd_range(int n1, int n2, float64x2 a, float64x2 *out);

/* Fused sincos / sinhcosh. One range-reduction + Taylor pair produces
 * both outputs, roughly halving the transcendental cost for call sites
 * that need both. Out-pointer style since C has no multi-value return. */
MULTIFLOATS_API void sincosdd(float64x2 a, float64x2 *s, float64x2 *c);
MULTIFLOATS_API void sinhcoshdd(float64x2 a, float64x2 *s, float64x2 *c);

/* Complex DD arithmetic. These are the canonical implementations; the
 * Fortran elemental `cdd + cdd`, `cdd * dp`, etc. routines wrap these via
 * bind(c). */
MULTIFLOATS_API complex64x2 cadddd(complex64x2 a, complex64x2 b);
MULTIFLOATS_API complex64x2 csubdd(complex64x2 a, complex64x2 b);
MULTIFLOATS_API complex64x2 cmuldd(complex64x2 a, complex64x2 b);
MULTIFLOATS_API complex64x2 cdivdd(complex64x2 a, complex64x2 b);

/* Complex DD transcendentals. Branch cuts match C99 Annex G (matching
 * libquadmath cexpq/clogq/csqrtq/...). Where the classic formula needs
 * both sin(y) and cos(y) or both sinh(y) and cosh(y), these use the fused
 * kernels internally so one range-reduction + Taylor pair covers both.
 * Precision and speed have been measured; the std::complex<float64x2>
 * specializations below delegate here where specialization wins (exp, sin,
 * cos, tan, sinh, cosh, tanh, atanh, acos). */
MULTIFLOATS_API complex64x2 cexpdd(complex64x2 z);
MULTIFLOATS_API complex64x2 cexpm1dd(complex64x2 z);
MULTIFLOATS_API complex64x2 clogdd(complex64x2 z);
MULTIFLOATS_API complex64x2 clog2dd(complex64x2 z);
MULTIFLOATS_API complex64x2 clog10dd(complex64x2 z);
MULTIFLOATS_API complex64x2 clog1pdd(complex64x2 z);
MULTIFLOATS_API complex64x2 cpowdd(complex64x2 z, complex64x2 w);
MULTIFLOATS_API complex64x2 csqrtdd(complex64x2 z);
MULTIFLOATS_API complex64x2 csindd(complex64x2 z);
MULTIFLOATS_API complex64x2 csinpidd(complex64x2 z);
MULTIFLOATS_API complex64x2 ccosdd(complex64x2 z);
MULTIFLOATS_API complex64x2 ccospidd(complex64x2 z);
MULTIFLOATS_API complex64x2 ctandd(complex64x2 z);
MULTIFLOATS_API complex64x2 casindd(complex64x2 z);
MULTIFLOATS_API complex64x2 cacosdd(complex64x2 z);
MULTIFLOATS_API complex64x2 catandd(complex64x2 z);
MULTIFLOATS_API complex64x2 csinhdd(complex64x2 z);
MULTIFLOATS_API complex64x2 ccoshdd(complex64x2 z);
MULTIFLOATS_API complex64x2 ctanhdd(complex64x2 z);
MULTIFLOATS_API complex64x2 casinhdd(complex64x2 z);
MULTIFLOATS_API complex64x2 cacoshdd(complex64x2 z);
MULTIFLOATS_API complex64x2 catanhdd(complex64x2 z);

/* Complex magnitude / argument / projection / conjugate / accessors.
 * cabsdd = |z|  (overflow-safe hypot).
 * cargdd = arg(z) in (-pi, pi].
 * cprojdd = projection onto the Riemann sphere (C99 7.3.9.4):
 *   infinities collapse to (+inf, copysign(0, imag)), else identity. */
MULTIFLOATS_API float64x2   cabsdd(complex64x2 z);
MULTIFLOATS_API float64x2   cargdd(complex64x2 z);
MULTIFLOATS_API complex64x2 cprojdd(complex64x2 z);
MULTIFLOATS_API complex64x2 conjdd(complex64x2 z);
MULTIFLOATS_API float64x2   crealdd(complex64x2 z);
MULTIFLOATS_API float64x2   cimagdd(complex64x2 z);

/* Matrix multiply (column-major, Fortran layout).
 *   matmuldd_mm: C(m,n) = A(m,k) * B(k,n)
 *   matmuldd_mv: y(m)   = A(m,k) * x(k)
 *   matmuldd_vm: y(n)   = x(k)   * B(k,n)
 * Leading dimensions equal the first extent (no strides).
 *
 * Scope vs BLAS GEMM. These are *not* a direct replacement for DGEMM:
 *
 *   - transa / transb: no. Always treats both operands in storage order;
 *     the caller must transpose the input arrays before the call if
 *     needed. (In Fortran, `matmul(transpose(a), b)` materialises the
 *     transpose, then the kernel walks the transposed buffer.)
 *   - alpha / beta:    no. The output is overwritten, not accumulated
 *     into. To compute `C := alpha*A*B + beta*C`, the caller must scale
 *     `A` or `B`, run matmul, then combine with C explicitly.
 *   - leading dimensions: always the first extent; no LDA/LDB/LDC.
 *
 * These constraints are deliberate: the compensated DD kernels use a
 * register-blocked panel design that assumes contiguous column-major
 * storage with the canonical shape. Adding transposed / strided / in-
 * place-accumulating variants is plausible but would require a new set
 * of panel dispatchers — tracked under "Deferred work" in
 * doc/developer/INTERNALS.md.
 *
 * renorm_interval: if > 0, renormalize accumulators every N reductions
 * (matches DD_FMA_RENORM_INTERVAL in the Fortran layer — keeps s_lo
 * bounded for large k). Pass 0 to renormalize only at the end. */
MULTIFLOATS_API void matmuldd_mm(const float64x2 *a, const float64x2 *b,
                         float64x2 *c,
                         int64_t m, int64_t k, int64_t n,
                         int64_t renorm_interval);
MULTIFLOATS_API void matmuldd_mv(const float64x2 *a, const float64x2 *x,
                         float64x2 *y,
                         int64_t m, int64_t k,
                         int64_t renorm_interval);
MULTIFLOATS_API void matmuldd_vm(const float64x2 *x, const float64x2 *b,
                         float64x2 *y,
                         int64_t k, int64_t n,
                         int64_t renorm_interval);

/* Formatted output. Scientific-notation decimal representation of `x`,
 * written into `[first, last)` with no NUL terminator — matching C++17
 * `std::to_chars` semantics. On success returns a pointer one past the
 * last byte written (so the written range is `[first, returned)`). On a
 * too-small buffer returns NULL and leaves the buffer unchanged. `precision`
 * is clamped to [1, 34] — DD carries ~32 significant digits, two extra for
 * round-half-to-even. Special values emit "nan", "inf", "-inf", "0e+00",
 * "-0e+00".
 *
 * A buffer of MULTIFLOATS_DD_CHARS_BUFSIZE bytes always fits any valid
 * output, so stack-allocating `char buf[MULTIFLOATS_DD_CHARS_BUFSIZE]` and
 * calling as `to_charsdd(x, prec, buf, buf + sizeof(buf))` is a safe use.
 * The C++ `multifloats::to_chars` / `to_string` / `operator<<` wrappers
 * below all layer on this so that `std::string` construction happens in
 * the consumer's TU — keeping the library free of `std::string` in its
 * ABI surface (libstdc++ `_GLIBCXX_USE_CXX11_ABI` dual-ABI safe). */
#define MULTIFLOATS_DD_CHARS_BUFSIZE 48
MULTIFLOATS_API char *to_charsdd(float64x2 x, int precision,
                                 char *first, char *last);

/* Comparison (return int: 1 = true, 0 = false) */
MULTIFLOATS_API int eqdd(float64x2 a, float64x2 b);
MULTIFLOATS_API int nedd(float64x2 a, float64x2 b);
MULTIFLOATS_API int ltdd(float64x2 a, float64x2 b);
MULTIFLOATS_API int ledd(float64x2 a, float64x2 b);
MULTIFLOATS_API int gtdd(float64x2 a, float64x2 b);
MULTIFLOATS_API int gedd(float64x2 a, float64x2 b);

/* libquadmath parity — every NAMEq in <quadmath.h> has a NAMEdd entry
 * here. Implementations are the canonical `multifloats::NAME` C++
 * helpers in the header section below; the C-ABI symbols in
 * src/float64x2_abi.inc are thin marshaling shims. */

/* Cube root, integer rounding (toward ±inf, nearest-current-mode), exponent
 * split / rescale, integer/fractional split, ulp-step, IEEE remainder. */
MULTIFLOATS_API float64x2 cbrtdd(float64x2 a);
MULTIFLOATS_API float64x2 ceildd(float64x2 a);
MULTIFLOATS_API float64x2 floordd(float64x2 a);
MULTIFLOATS_API float64x2 nearbyintdd(float64x2 a);
MULTIFLOATS_API float64x2 rintdd(float64x2 a);
MULTIFLOATS_API float64x2 logbdd(float64x2 a);
MULTIFLOATS_API float64x2 frexpdd(float64x2 a, int *exp);
MULTIFLOATS_API float64x2 modfdd(float64x2 a, float64x2 *iptr);
MULTIFLOATS_API float64x2 ldexpdd(float64x2 a, int n);
MULTIFLOATS_API float64x2 scalbndd(float64x2 a, int n);
MULTIFLOATS_API float64x2 scalblndd(float64x2 a, long n);
MULTIFLOATS_API float64x2 nextafterdd(float64x2 a, float64x2 b);
MULTIFLOATS_API float64x2 remainderdd(float64x2 a, float64x2 b);
MULTIFLOATS_API float64x2 remquodd(float64x2 a, float64x2 b, int *quo);
MULTIFLOATS_API int ilogbdd(float64x2 a);

/* Integer rounding — leading-limb's libm result with a half-integer
 * fixup that consults the lo limb (so e.g. lround((0.5, -ε)) = 0). */
MULTIFLOATS_API long lrounddd(float64x2 a);
MULTIFLOATS_API long long llrounddd(float64x2 a);
MULTIFLOATS_API long lrintdd(float64x2 a);
MULTIFLOATS_API long long llrintdd(float64x2 a);

/* Classification — C99 isnan/isinf/signbit/isfinite contract: non-zero
 * for true. `finitedd` mirrors libquadmath's legacy `finiteq` spelling
 * (alias for isfinitedd). */
MULTIFLOATS_API int isnandd(float64x2 a);
MULTIFLOATS_API int isinfdd(float64x2 a);
MULTIFLOATS_API int isfinitedd(float64x2 a);
MULTIFLOATS_API int finitedd(float64x2 a);
MULTIFLOATS_API int signbitdd(float64x2 a);
/* Signaling-NaN test on the leading limb. Uses the compiler builtin
 * where available (gcc >= 13, clang with __has_builtin); otherwise
 * returns 0 — multifloats never constructs sNaNs internally, so the
 * fallback only loses signal when the caller explicitly passed an
 * sNaN through the C ABI. */
MULTIFLOATS_API int issignalingdd(float64x2 a);

/* NaN with payload tag (libquadmath nanq parity). Calls std::nan(tagp)
 * for the leading limb; lo is set to 0. */
MULTIFLOATS_API float64x2 nandd(const char *tagp);

/* cis(x) = cos(x) + i sin(x). Mirrors libquadmath cexpiq. */
MULTIFLOATS_API complex64x2 cexpidd(float64x2 a);

#ifdef __cplusplus
}  /* extern "C" */
#if defined(__clang__)
#  pragma clang diagnostic pop
#endif
#endif

/* ============================================================================
 * C++-only section — free functions, transcendental decls, matmul, to_chars.
 * Still inside `namespace multifloats { }` opened at the top of the header.
 * ========================================================================= */
#ifdef __cplusplus

// =============================================================================
// <cmath>-style free functions (ADL on float64x2)
//
// Canonical-implementation convention: every header-inline definition below
// is the *single source of truth* for that operation. The extern "C"
// `NAMEdd` symbols in the C-ABI block above (declared, defined in
// src/float64x2_abi.inc) are thin marshaling shims that call the
// corresponding `multifloats::NAME` here — never separate implementations
// that could drift. Out-of-line transcendentals (exp, log, trig, …) flip
// the relationship: the header carries only the `multifloats::NAME`
// declaration, the body lives in src/multifloats_math.cc, and the C-ABI
// shim still calls through `multifloats::NAME`. In both directions the
// C-ABI surface is a forwarder, not a parallel implementation.
// =============================================================================

namespace detail {
// Index of the first nonzero limb, or 2 if every limb is a (possibly signed)
// zero. Used by abs / signbit to resolve the sign of non-canonical DDs like
// (+0, -eps), where signbit(hi) alone would misclassify.
constexpr std::size_t first_nonzero_limb_index(float64x2 const &x) {
  if (x.limbs[0] != 0.0) return 0;
  if (x.limbs[1] != 0.0) return 1;
  return 2;
}

// Triple-double scratch primitives (float64x3, td_add_*, td_mul_*, …)
// live in the internal multifloats_td.hh. Only kernels inside
// src/multifloats_math.cc and the direct-primitive tests pull them in;
// they are not part of the public ABI surface exposed here.

} // namespace detail

inline constexpr float64x2 fabs(float64x2 const &x) {
  std::size_t i = detail::first_nonzero_limb_index(x);
  // All limbs are a (possibly signed) zero: return canonical +0 so
  // fabs((-0, 0)) yields (+0, +0), matching IEEE fabs(-0.0) = +0.0.
  if (i == 2) return float64x2();
  return std::signbit(x.limbs[i]) ? -x : x;
}

inline constexpr float64x2 fmin(float64x2 const &a, float64x2 const &b) {
  return (a < b) ? a : b;
}

inline constexpr float64x2 fmax(float64x2 const &a, float64x2 const &b) {
  return (a < b) ? b : a;
}


inline constexpr float64x2 abs(float64x2 const &x) { return fabs(x); }
inline constexpr float64x2 min(float64x2 const &x, float64x2 const &y) { return fmin(x, y); }
inline constexpr float64x2 max(float64x2 const &x, float64x2 const &y) { return fmax(x, y); }

inline constexpr bool signbit(float64x2 const &x) {
  // For non-canonical zero-hi DDs (e.g. (+0, -eps)), the sign lives in
  // the first nonzero limb. Fall through to signbit(hi) when every limb
  // is a (possibly signed) zero, preserving IEEE -0 semantics.
  std::size_t i = detail::first_nonzero_limb_index(x);
  return std::signbit(x.limbs[i == 2 ? 0 : i]);
}

inline constexpr bool isfinite(float64x2 const &x) {
  // A DD with finite hi and non-finite lo is classified non-finite — this
  // matters for the operator+ short-circuit.
  for (std::size_t i = 0; i < 2; ++i) {
    if (!std::isfinite(x.limbs[i])) return false;
  }
  return true;
}

inline constexpr bool isinf(float64x2 const &x) {
  // Scan both limbs for symmetry with isnan / isfinite: a DD with a finite
  // hi and an inf lo (non-canonical but constructible via the C ABI) should
  // classify as inf, not "neither finite nor inf".
  for (std::size_t i = 0; i < 2; ++i) {
    if (std::isinf(x.limbs[i])) return true;
  }
  return false;
}

inline constexpr bool isnan(float64x2 const &x) {
  for (std::size_t i = 0; i < 2; ++i) {
    if (std::isnan(x.limbs[i])) return true;
  }
  return false;
}

inline constexpr int fpclassify(float64x2 const &x) {
  return std::fpclassify(x.limbs[0]);
}

inline constexpr float64x2 ldexp(float64x2 const &x, int n) {
  // Build the power-of-two scale once; multiplication by an exact power of
  // two is exact for every limb (no rounding, no renorm), avoiding the
  // two library calls of std::ldexp.
  double scale = std::ldexp(1.0, n);
  return {x.limbs[0] * scale, x.limbs[1] * scale};
}

inline constexpr float64x2 scalbn(float64x2 const &x, int n) {
  // POSIX alias of ldexp for FLT_RADIX == 2 (which is guaranteed by IEEE 754).
  return ldexp(x, n);
}

inline constexpr int ilogb(float64x2 const &x) {
  return std::ilogb(x.limbs[0]);
}

inline constexpr float64x2 copysign(float64x2 const &x, float64x2 const &y) {
  return (signbit(x) == signbit(y)) ? x : -x;
}

namespace detail {

// fast_two_sum (assumes |hi| >= |lo|), in-place renormalization helper for
// float64x2 rounding paths.
constexpr void renorm_fast(double &hi, double &lo) {
  double s = hi + lo;
  double b = s - hi;
  double e = lo - b;
  hi = s;
  lo = e;
}

} // namespace detail

// =============================================================================
// Rounding and integer-valued functions
// =============================================================================

inline constexpr float64x2 floor(float64x2 const &x) {
  float64x2 r;
  double fl_hi = std::floor(x.limbs[0]);
  if (fl_hi == x.limbs[0]) {
    // hi is already an integer; floor depends on the lo limb.
    r.limbs[0] = fl_hi;
    r.limbs[1] = std::floor(x.limbs[1]);
    detail::renorm_fast(r.limbs[0], r.limbs[1]);
  } else {
    r.limbs[0] = fl_hi;
    r.limbs[1] = 0.0;
  }
  return r;
}

inline constexpr float64x2 ceil(float64x2 const &x) {
  float64x2 r;
  double cl_hi = std::ceil(x.limbs[0]);
  if (cl_hi == x.limbs[0]) {
    r.limbs[0] = cl_hi;
    r.limbs[1] = std::ceil(x.limbs[1]);
    detail::renorm_fast(r.limbs[0], r.limbs[1]);
  } else {
    r.limbs[0] = cl_hi;
    r.limbs[1] = 0.0;
  }
  return r;
}

inline constexpr float64x2 trunc(float64x2 const &x) {
  return std::signbit(x.limbs[0]) ? -floor(-x) : floor(x);
}

inline constexpr float64x2 round(float64x2 const &x) {
  // Round half away from zero, matching std::round. Two half-integer
  // hazards handled here:
  //   * hi itself half-integer (e.g. 2.5): std::round jumps away from zero,
  //     undone when lo lies on the other side of the half-boundary.
  //   * hi exact integer with lo == ±0.5 (possible once ulp(hi) ≥ 1, i.e.
  //     |hi| ≥ 2^53): if sign(lo) opposes sign(hi), the true value is
  //     closer to zero, so the correct rounded value is hi itself rather
  //     than hi ± 1 that std::round(lo) would add.
  float64x2 r;
  double hi = std::round(x.limbs[0]);
  if (hi == x.limbs[0]) {
    double lo = x.limbs[1];
    double rlo = 0.0;
    if      (lo ==  0.5 && hi <  0.0) rlo = 0.0;
    else if (lo == -0.5 && hi >  0.0) rlo = 0.0;
    else                              rlo = std::round(lo);
    r.limbs[0] = hi;
    r.limbs[1] = rlo;
    detail::renorm_fast(r.limbs[0], r.limbs[1]);
  } else {
    double diff = x.limbs[0] - hi;
    if      (diff == -0.5 && x.limbs[1] < 0.0) hi -= 1.0;
    else if (diff ==  0.5 && x.limbs[1] > 0.0) hi += 1.0;
    r.limbs[0] = hi;
    r.limbs[1] = 0.0;
  }
  return r;
}

inline constexpr float64x2 nearbyint(float64x2 const &x) {
  float64x2 r;
  double hi = std::nearbyint(x.limbs[0]);
  if (hi == x.limbs[0]) {
    r.limbs[0] = hi;
    r.limbs[1] = std::nearbyint(x.limbs[1]);
    detail::renorm_fast(r.limbs[0], r.limbs[1]);
  } else {
    // Half-integer hi: std::nearbyint rounds to even and ignores lo, but
    // the true value = hi + lo can lie on the opposite side of the half
    // boundary. Mirror the adjustment in round() above.
    double diff = x.limbs[0] - hi;
    if      (diff ==  0.5 && x.limbs[1] > 0.0) hi += 1.0;
    else if (diff == -0.5 && x.limbs[1] < 0.0) hi -= 1.0;
    r.limbs[0] = hi;
    r.limbs[1] = 0.0;
  }
  return r;
}

inline constexpr float64x2 rint(float64x2 const &x) { return nearbyint(x); }

namespace detail {
// Shared half-integer correction for lround / llround. Given i = std::[l]lround(x_hi),
// adjust ±1 based on how lo crosses the half-integer boundary of the true value.
template <typename Int>
constexpr Int lround_adjust(float64x2 const &x, Int i) {
  double hi = x.limbs[0];
  double lo = x.limbs[1];
  double diff = hi - double(i);
  if (diff == 0.0) {
    // hi exact integer; lo (bounded by ulp(hi)/2) decides.
    if      (lo >   0.5)                 ++i;
    else if (lo <  -0.5)                 --i;
    else if (lo ==  0.5 && hi >=  0.0)   ++i;
    else if (lo == -0.5 && hi <=  0.0)   --i;
  } else if (diff == -0.5 && lo < 0.0) --i;
  else if   (diff ==  0.5 && lo > 0.0) ++i;
  return i;
}
} // namespace detail

inline constexpr long lround(float64x2 const &x) {
  return detail::lround_adjust<long>(x, std::lround(x.limbs[0]));
}

inline constexpr long long llround(float64x2 const &x) {
  return detail::lround_adjust<long long>(x, std::llround(x.limbs[0]));
}

inline constexpr long lrint(float64x2 const &x) {
  return std::lrint(rint(x).limbs[0]);
}

inline constexpr long long llrint(float64x2 const &x) {
  return std::llrint(rint(x).limbs[0]);
}

// =============================================================================
// Floating-point manipulation
// =============================================================================

inline constexpr float64x2 frexp(float64x2 const &x, int *exp) {
  float64x2 r;
  int e = 0;
  r.limbs[0] = std::frexp(x.limbs[0], &e);
  r.limbs[1] = std::ldexp(x.limbs[1], -e);
  *exp = e;
  return r;
}

inline constexpr float64x2 modf(float64x2 const &x, float64x2 *iptr) {
  *iptr = trunc(x);
  return x - *iptr;
}

inline constexpr float64x2 scalbln(float64x2 const &x, long n) {
  return {std::scalbln(x.limbs[0], n), std::scalbln(x.limbs[1], n)};
}

inline constexpr float64x2 logb(float64x2 const &x) {
  float64x2 r;
  r.limbs[0] = std::logb(x.limbs[0]);
  return r;
}

inline constexpr float64x2 nextafter(float64x2 const &x, float64x2 const &y) {
  if (x == y) return y;
  // One DD ulp ≈ ulp_up(|hi|) * 2^-53. Always use the upward ulp of |hi|;
  // the downward ulp halves at a power-of-2 boundary, so picking it there
  // would make the step 2× too small and break the round-trip identity
  // nextafter(nextafter(x, +inf), -inf) == x.
  double ax = std::abs(x.limbs[0]);
  double inf = std::numeric_limits<double>::infinity();
  // x.hi == 0 special: ulp_up(0) is min_subnormal ≈ 5e-324, and `ulp/2^53`
  // underflows to 0 — leaving us returning `x ± 0 == x`, which would loop
  // any caller that walks DD numbers via nextafter. Step by the smallest
  // representable double instead. (No DD lo limb can survive at this scale.)
  if (ax == 0.0) {
    double step = std::numeric_limits<double>::denorm_min();
    return (x < y) ? float64x2(step) : float64x2(-step);
  }
  double ulp = std::nextafter(ax, inf) - ax;
  double eps = std::ldexp(ulp, -53);
  return (x < y) ? x + float64x2(eps) : x - float64x2(eps);
}

inline constexpr float64x2 nexttoward(float64x2 const &x, float64x2 const &y) {
  return nextafter(x, y);
}

// libquadmath nanq parity. The tag string encodes the NaN payload via
// std::nan; lo is set to 0 so the result is a canonical DD NaN.
inline float64x2 nan(const char *tagp) {
  return float64x2(std::nan(tagp ? tagp : ""), 0.0);
}

// =============================================================================
// Basic arithmetic helpers
// =============================================================================

inline constexpr float64x2 fma(float64x2 const &x, float64x2 const &y,
                               float64x2 const &z) {
  // Not a hardware fma, but provides the cmath interface.
  return x * y + z;
}

inline constexpr float64x2 fmod(float64x2 const &x, float64x2 const &y) {
  // Reduction step picks q from the ilogb gap between r and ay:
  // gap ≤ 53 — scalar q fits in one double (the earlier all-scalar form
  // silently lost q's low integer bits past 2^53); gap > 53 — DD-level
  // trunc(r/ay) carries ~106 integer bits, so a single step drops gap by
  // ≥ 53 and iteration converges in O(gap/53) steps even past 2^106. DD
  // rounding of r − q·ay can leave a tiny negative residue; add-back
  // uses the same gap dispatch so recovery stays O(1).
  bool x_neg = x.limbs[0] < 0.0;
  float64x2 ax = x_neg ? -x : x;
  float64x2 ay = (y.limbs[0] < 0.0) ? -y : y;

  if (ax < ay) return x;

  float64x2 r = ax;
  while (true) {
    if (r.limbs[0] < 0.0) {
      double r_abs = -r.limbs[0];
      int gap = std::ilogb(r_abs) - std::ilogb(ay.limbs[0]);
      if (gap <= 0) {
        r = r + ay;
      } else if (gap <= 53) {
        double q = std::trunc(r_abs / ay.limbs[0]) + 1.0;
        r = r + ay * float64x2(q);
      } else {
        r = r + (trunc(-r / ay) + float64x2(1.0)) * ay;
      }
    } else if (r >= ay) {
      int gap = std::ilogb(r.limbs[0]) - std::ilogb(ay.limbs[0]);
      if (gap <= 53) {
        double q = std::trunc(r.limbs[0] / ay.limbs[0]);
        r = (q <= 1.0) ? (r - ay) : (r - ay * float64x2(q));
      } else {
        r = r - trunc(r / ay) * ay;
      }
    } else {
      break;
    }
    if (r.limbs[0] == 0.0 && r.limbs[1] == 0.0) break;
  }

  return x_neg ? -r : r;
}

// Floored modulo (matches Fortran `modulo` intrinsic): the result has
// the same sign as y, or is +0. Uses fmod for the truncated remainder
// then lifts to the floored form by adding y when the remainder and y
// have opposite signs. signbit() resolves the sign via the first
// nonzero limb, so non-canonical DDs (e.g. (+0, -eps)) get the correct
// sign adjustment — a subtlety the naive hi-only check misses.
inline constexpr float64x2 modulo(float64x2 const &x, float64x2 const &y) {
  float64x2 r = fmod(x, y);
  if (r.limbs[0] == 0.0 && r.limbs[1] == 0.0) return r;
  return (signbit(r) != signbit(y)) ? r + y : r;
}

// Integer-exponent power via exponentiation-by-squaring. Returns exact
// 1 for n == 0 (including for 0**0, matching C pow and Fortran). For
// n < 0 returns 1 / powi(base, -n), adding one DD division at the end.
inline constexpr float64x2 powi(float64x2 base, int n) {
  if (n == 0) return float64x2(1.0);
  bool neg = (n < 0);
  // Promote through long long so INT_MIN negates without overflowing int.
  unsigned long long u = neg
      ? static_cast<unsigned long long>(-static_cast<long long>(n))
      : static_cast<unsigned long long>(n);
  float64x2 result(1.0);
  float64x2 b = base;
  for (;;) {
    if (u & 1ull) result = result * b;
    u >>= 1;
    if (u == 0ull) break;
    b = b * b;
  }
  return neg ? (float64x2(1.0) / result) : result;
}

inline constexpr float64x2 remainder(float64x2 const &x, float64x2 const &y) {
  return x - round(x / y) * y;
}

inline constexpr float64x2 remquo(float64x2 const &x, float64x2 const &y, int *quo) {
  float64x2 q = round(x / y);
  *quo = static_cast<int>(q.limbs[0]);
  return x - q * y;
}

inline constexpr float64x2 fdim(float64x2 const &x, float64x2 const &y) {
  return (x > y) ? (x - y) : float64x2();
}

// C++20 std::lerp: exact at the endpoints, monotonic in t, and does not
// overshoot when a and b have the same sign and t is in [0, 1].
inline constexpr float64x2 lerp(float64x2 const &a, float64x2 const &b,
                      float64x2 const &t) {
  if ((a.limbs[0] <= 0.0 && b.limbs[0] >= 0.0) ||
      (a.limbs[0] >= 0.0 && b.limbs[0] <= 0.0)) {
    // Opposite signs (or one is zero): no cancellation risk.
    return t * b + (float64x2(1.0) - t) * a;
  }
  if (t == float64x2(1.0)) {
    return b;  // exact endpoint per C++20 spec
  }
  float64x2 x = a + t * (b - a);
  // Enforce monotonicity at the b end when t is past 1 or when rounding
  // nudges x beyond b — matches libstdc++/libc++ behavior.
  if ((t.limbs[0] > 1.0) == (b > a)) {
    return (b > x) ? b : x;
  }
  return (x > b) ? b : x;
}

// =============================================================================
// detail:: inline helpers for DD polynomial evaluation
//
// These are used both by the extern "C" implementations in multifloats_math.cc
// and by the header-only helpers (sqrt, trunc, etc.).
// =============================================================================

// Forward declaration so detail kernels below can use multifloats::sqrt via ADL.
inline constexpr float64x2 sqrt(float64x2 const &x);

namespace detail {

// Polynomial and conversion constants needed by the DD kernels are provided
// through dd_constants.hh which is included only in the .cc file. The
// erf/erfc rational-approximation constants also live in the .cc file.
//
// horner / neval / deval (DD polynomial evaluators) used to live here as
// inline definitions but have been moved to src/float64x2_poly.inc,
// which is pulled into namespace multifloats::detail inside
// multifloats_math.cc. Every call site lives in that TU, so inlining
// decisions are unchanged; the public header no longer carries ~250 lines
// of Estrin switch bodies.

} // namespace detail

// =============================================================================
// Roots
// =============================================================================

inline constexpr float64x2 sqrt(float64x2 const &x) {
  double s = std::sqrt(x.limbs[0]);
  // Bail on 0, -0, negative, NaN, +Inf — the Karp-Markstein refinement
  // would compute `inf - inf = NaN` in the residual step for +Inf, and
  // `0 / 0` for 0. Leading-limb sqrt handles every IEEE special case.
  if (!(x.limbs[0] > 0.0) || !std::isfinite(s)) {
    float64x2 r;
    r.limbs[0] = s;
    return r;
  }
  // Karp/Markstein: r = s + (x - s*s) / (2s), evaluated in DD. The
  // correction reduces the DD residual to a scalar via
  // `residual.limbs[0] * (0.5/s)`, so the residual's lo limb is
  // dropped on the floor. Two higher-fidelity variants were measured
  // (see doc/developer/INTERNALS.md anchor P1):
  //   (a) full DD divide `residual / (2*s_dd)` — sqrt worst case near
  //       perfect squares goes 0.76 → 0.39 ulp, but sqrt bench drops
  //       ~55% and hypot/acosh take a 10–25% hit.
  //   (b) `residual * float64x2(0.5/s)` (DD × scalar) — 0.76 → 0.58
  //       ulp, sqrt bench drops ~30%.
  // Baseline is already sub-1-ulp (0 ulp on exact k², ≤0.76 ulp with a
  // non-zero lo limb). The gain from (a)/(b) isn't worth the speed
  // regression for this library's usage pattern; keep baseline.
  const float64x2 s_dd(s);
  const float64x2 residual = x - s_dd * s_dd;
  const float64x2 correction(residual.limbs[0] * (0.5 / s));
  return s_dd + correction;
}

inline constexpr float64x2 cbrt(float64x2 const &x) {
  const double s = std::cbrt(x.limbs[0]);
  // Bail on 0/-0/±Inf/NaN — the residual would compute `inf - inf` for ±Inf
  // and `0/0` for ±0, poisoning both limbs. std::cbrt already returns the
  // correctly-signed special-value result for the leading limb.
  if (x.limbs[0] == 0.0 || !std::isfinite(s)) {
    float64x2 r;
    r.limbs[0] = s;
    return r;
  }
  const float64x2 s_dd(s);
  const float64x2 residual = x - s_dd * s_dd * s_dd;
  const float64x2 correction(residual.limbs[0] / (3.0 * s * s));
  return s_dd + correction;
}

inline constexpr float64x2 hypot(float64x2 const &x, float64x2 const &y) {
  // Defer to libm's hypot for non-finite so inf/NaN propagate correctly.
  if (!std::isfinite(x.limbs[0]) || !std::isfinite(y.limbs[0])) {
    float64x2 r;
    r.limbs[0] = std::hypot(x.limbs[0], y.limbs[0]);
    r.limbs[1] = 0.0;
    return r;
  }
  float64x2 ax = signbit(x) ? -x : x;
  float64x2 ay = signbit(y) ? -y : y;
  float64x2 big = (ax > ay) ? ax : ay;
  float64x2 small = (ax > ay) ? ay : ax;
  if (big.limbs[0] == 0.0) return float64x2();
  // Power-of-2 scale (exact) so big has exponent 0 before the square.
  // Replace per-limb ldexp calls with multiplies by 2^(±e): for a power-
  // of-2 multiplier and a non-subnormal result, `x * 2^k` and
  // `ldexp(x, k)` are bit-identical, but the multiply is one FP op while
  // ldexp is a libm call. `e` comes from `ilogb(big.hi)` on a finite
  // non-zero input, so `|e| ≤ 1023`; `2^(-e)` and `2^e` are both finite.
  int e = std::ilogb(big.limbs[0]);
  double down = std::ldexp(1.0, -e);
  big.limbs[0]   *= down; big.limbs[1]   *= down;
  small.limbs[0] *= down; small.limbs[1] *= down;
  float64x2 ratio = small / big;
  float64x2 root = big * sqrt(float64x2(1.0) + ratio * ratio);
  float64x2 r;
  double up = std::ldexp(1.0, e);
  r.limbs[0] = root.limbs[0] * up;
  r.limbs[1] = root.limbs[1] * up;
  // Overflow: true result exceeds double's range. Zero the trailing limb so
  // callers see a clean inf rather than (inf, NaN).
  if (!std::isfinite(r.limbs[0])) r.limbs[1] = 0.0;
  return r;
}

// =============================================================================
// Transcendentals — definitions live in multifloats_math.cc. The extern "C"
// `*dd` kernels are thin marshaling shims around these same functions, so
// C++ callers go straight to the C++ body with no ABI-crossing overhead.
// =============================================================================

// Power, exponential and logarithm
float64x2 exp   (float64x2 const &x);
float64x2 exp2  (float64x2 const &x);
float64x2 expm1 (float64x2 const &x);
float64x2 log   (float64x2 const &x);
float64x2 log10 (float64x2 const &x);
float64x2 log2  (float64x2 const &x);
float64x2 log1p (float64x2 const &x);
float64x2 pow   (float64x2 const &x, float64x2 const &y);

// Trigonometric
float64x2 sin   (float64x2 const &x);
float64x2 cos   (float64x2 const &x);
float64x2 tan   (float64x2 const &x);
float64x2 asin  (float64x2 const &x);
float64x2 acos  (float64x2 const &x);
float64x2 atan  (float64x2 const &x);
float64x2 atan2 (float64x2 const &y, float64x2 const &x);

// π-scaled trig: {sin,cos,tan}pi(x) = {sin,cos,tan}(π·x),
//                {asin,acos,atan}pi(x) = {asin,acos,atan}(x)/π,
//                atan2pi(y, x) = atan2(y, x)/π.
float64x2 sinpi   (float64x2 const &x);
float64x2 cospi   (float64x2 const &x);
float64x2 tanpi   (float64x2 const &x);
float64x2 asinpi  (float64x2 const &x);
float64x2 acospi  (float64x2 const &x);
float64x2 atanpi  (float64x2 const &x);
float64x2 atan2pi (float64x2 const &y, float64x2 const &x);

// Fused sincos / sinhcosh. One range-reduction + Taylor pair produces both
// outputs, roughly halving the transcendental cost when both are needed.
void sincos   (float64x2 const &x, float64x2 &s, float64x2 &c);
void sinhcosh (float64x2 const &x, float64x2 &s, float64x2 &c);

// Hyperbolic
float64x2 sinh  (float64x2 const &x);
float64x2 cosh  (float64x2 const &x);
float64x2 tanh  (float64x2 const &x);
float64x2 asinh (float64x2 const &x);
float64x2 acosh (float64x2 const &x);
float64x2 atanh (float64x2 const &x);

// Error and gamma. erfcx is the scaled complementary error function,
// erfcx(x) = exp(x^2) * erfc(x).
float64x2 erf    (float64x2 const &x);
float64x2 erfc   (float64x2 const &x);
float64x2 erfcx  (float64x2 const &x);
float64x2 tgamma (float64x2 const &x);
float64x2 lgamma (float64x2 const &x);

// Bessel functions of the first (j) and second (y) kind. Names follow
// POSIX <math.h>: j0/j1/jn for J_n, y0/y1/yn for Y_n. Note the signature
// differs from C++17 std::cyl_bessel_j / std::cyl_neumann — those take a
// double order, while the DD kernels take an integer order seeded off
// the fast-path j0/j1/y0/y1 rational fits and step up via Miller /
// forward recurrence for |n| ≥ 2. The 4-arg yn overload fills a single
// forward-recurrence sweep of n2-n1+1 outputs (cheaper than n2-n1+1
// independent 2-arg calls); the 4-arg jn overload loops the scalar
// Miller kernel because Jn forward-recurrence is unstable for n > x.
float64x2 j0 (float64x2 const &x);
float64x2 j1 (float64x2 const &x);
float64x2 y0 (float64x2 const &x);
float64x2 y1 (float64x2 const &x);
float64x2 jn (int n, float64x2 const &x);
float64x2 yn (int n, float64x2 const &x);
void      jn (int n1, int n2, float64x2 const &x, float64x2 *out);
void      yn (int n1, int n2, float64x2 const &x, float64x2 *out);

// =============================================================================
// Additional classification and ordered comparison
// =============================================================================

inline constexpr bool isnormal(float64x2 const &x) {
  return std::isnormal(x.limbs[0]);
}

inline constexpr bool isgreater(float64x2 const &x, float64x2 const &y) {
  return !isnan(x) && !isnan(y) && (x > y);
}

inline constexpr bool isgreaterequal(float64x2 const &x, float64x2 const &y) {
  return !isnan(x) && !isnan(y) && (x >= y);
}

inline constexpr bool isless(float64x2 const &x, float64x2 const &y) {
  return !isnan(x) && !isnan(y) && (x < y);
}

inline constexpr bool islessequal(float64x2 const &x, float64x2 const &y) {
  return !isnan(x) && !isnan(y) && (x <= y);
}

inline constexpr bool islessgreater(float64x2 const &x, float64x2 const &y) {
  return !isnan(x) && !isnan(y) && (x != y);
}

inline constexpr bool isunordered(float64x2 const &x, float64x2 const &y) {
  return isnan(x) || isnan(y);
}

// Layout sanity: float64x2 must remain standard-layout and trivially-copyable.
// The Fortran bind(c) interop, the C-ABI pass-by-value convention, and the
// matmul panel dispatchers all assume a two-double storage image with no
// hidden members, virtuals, or base classes. The user-provided constructors
// make it non-POD in the C++98 sense but preserve both properties.
static_assert(std::is_standard_layout<float64x2>::value,
              "float64x2 must be standard-layout");
static_assert(std::is_trivially_copyable<float64x2>::value,
              "float64x2 must be trivially copyable");

// ---- Compensated DD matrix multiply ----------------------------------------
//
// Column-major (Fortran-layout) GEMM/GEMV/VM panels with a compensated DD
// accumulator. Definitions live in src/multifloats_math.cc; the extern "C"
// `matmuldd_*` shims in the C section of this header forward directly to
// these (since `float64x2` is a single unified type across C and C++, no
// reinterpret_cast at the boundary). See the block comment on matmuldd_mm
// in the C section for the deliberate scope-down vs BLAS GEMM (no transpose
// flags, no alpha/beta, no leading-dimension arguments — always contiguous
// column-major).
//
// `renorm_interval` controls mid-accumulation renormalization of the (hi,
// lo) DD accumulator pair:
//   • 0  — renormalize only at the very end (cheapest; drifts on very large k).
//   • N>0 — two_sum every N reductions, keeping `lo` bounded. Matches
//           DD_FMA_RENORM_INTERVAL in the Fortran layer.
void matmul_mm(float64x2 const *a, float64x2 const *b, float64x2 *c,
               std::int64_t m, std::int64_t k, std::int64_t n,
               std::int64_t renorm_interval = 0);
void matmul_mv(float64x2 const *a, float64x2 const *x, float64x2 *y,
               std::int64_t m, std::int64_t k,
               std::int64_t renorm_interval = 0);
void matmul_vm(float64x2 const *x, float64x2 const *b, float64x2 *y,
               std::int64_t k, std::int64_t n,
               std::int64_t renorm_interval = 0);

// ---- Formatted I/O for float64x2 -------------------------------------------
//
// Scientific-notation decimal output at up to 34 significant digits (the
// natural DD significand limit is ~32; we allow two guard digits for
// round-half-to-even). Special values use the same textual forms as
// std::to_string for doubles ("nan", "inf", "-inf"). `operator<<` honors
// `os.precision()` when > 17; otherwise (including the C++ default of 6)
// the DD default of 32 digits is used.
//
// Layering. The core formatter is in src/multifloats_io.cc as a
// file-local helper, reached in the library TU by `operator<<`, by the
// out-of-line `multifloats::to_chars` (primary 5-arg form), and by the
// extern "C" `to_charsdd` shim in the C section above. The 3-arg and 4-arg
// `to_chars` overloads here are thin inline delegators to the 5-arg.
// `to_string` is header-inline — it layers on `to_chars` and constructs
// the returned `std::string` in the consumer's translation unit, so
// `libmultifloats.a` carries no `std::string` in its ABI and links cleanly
// against either libstdc++ string ABI (_GLIBCXX_USE_CXX11_ABI=0 or =1).
// `std::to_chars_result` in contrast is an ABI-stable pair `{char*, errc}`
// — no dual-ABI concern — so the 5-arg `to_chars` can safely be defined
// out-of-line in the library.

// `<charconv>`-style formatter for float64x2. Mirrors `std::to_chars` for
// the built-in floating-point types: writes into [first, last) without a
// NUL terminator and returns the one-past-the-last pointer together with a
// `std::errc`. On a too-small buffer, returns `{last, errc::value_too_large}`
// and leaves the buffer unchanged (nothing is written). On success, returns
// `{ptr, errc{}}` where `ptr == first + (bytes written)`.
//
// Scope: only `std::chars_format::scientific` is supported. `fixed`,
// `general`, and `hex` return `{first, errc::invalid_argument}` — the
// underlying formatter only implements scientific notation today. Precision
// defaults to 32 (DD's natural round-trip width); values outside [1, 34]
// are clamped by the core formatter.
std::to_chars_result
to_chars(char *first, char *last, float64x2 const &value,
         std::chars_format fmt, int precision) noexcept;
inline std::to_chars_result
to_chars(char *first, char *last, float64x2 const &value,
         std::chars_format fmt) noexcept {
  return to_chars(first, last, value, fmt, 32);
}
inline std::to_chars_result
to_chars(char *first, char *last, float64x2 const &value) noexcept {
  return to_chars(first, last, value, std::chars_format::scientific, 32);
}

inline std::string to_string(float64x2 const &x, int precision = 32) {
  char buf[MULTIFLOATS_DD_CHARS_BUFSIZE];
  auto res = to_chars(buf, buf + sizeof(buf), x,
                      std::chars_format::scientific, precision);
  // The caller buffer is sized for any legal precision, so value_too_large
  // cannot happen here — `res.ec` is always `errc{}` in practice.
  return std::string(buf, static_cast<std::size_t>(res.ptr - buf));
}
std::ostream &operator<<(std::ostream &os, float64x2 const &x);

} // namespace multifloats

// ---- std::complex<multifloats::float64x2> specializations ------------------
//
// Explicit specializations for the C++17 <complex> free-function templates
// on `std::complex<multifloats::float64x2>`. Definitions live in
// multifloats_math.cc (complex64x2_std.inc); the C-ABI
// `c*dd` symbols in the C section above are now thin wrappers that marshal
// complex64x2 ↔ std::complex<multifloats::float64x2> and call these.
//
// Covers:
//   • transcendentals: exp, log, sqrt, pow, sin, cos, tan, asin, acos,
//     atan, sinh, cosh, tanh, asinh, acosh, atanh (16)
//   • value accessors: abs, arg, proj (3)
//
// std::conj / std::real / std::imag ride the generic <complex> template
// unchanged — they compile to componentwise limb ops already. The
// Kind-D overloads that don't have std:: free functions (log1p, log2,
// expm1, sinpi, cospi) live in `namespace multifloats` below.

namespace std {

template <> complex<multifloats::float64x2> exp  (complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> log  (complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> sqrt (complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> pow  (complex<multifloats::float64x2> const &z,
                                                  complex<multifloats::float64x2> const &w);
template <> complex<multifloats::float64x2> sin  (complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> cos  (complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> tan  (complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> asin (complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> acos (complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> atan (complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> sinh (complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> cosh (complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> tanh (complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> asinh(complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> acosh(complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> atanh(complex<multifloats::float64x2> const &z);

template <> multifloats::float64x2          abs  (complex<multifloats::float64x2> const &z);
template <> multifloats::float64x2          arg  (complex<multifloats::float64x2> const &z);
template <> complex<multifloats::float64x2> proj (complex<multifloats::float64x2> const &z);

} // namespace std

// ---- Complex overloads in `namespace multifloats` --------------------------
//
// These are the C++ parity surface for c*dd kernels that either have no
// std::complex overload (expm1, log2, log1p) or whose names are outside
// the standard library (sinpi, cospi). They're `multifloats::sinpi(z)`,
// `multifloats::expm1(z)` etc. — not std:: specializations, because C++17
// doesn't have the corresponding free functions on std::complex to specialize.
// Definitions live in multifloats_math.cc alongside the scalar overloads.
namespace multifloats {
std::complex<float64x2> sinpi (std::complex<float64x2> const &z);
std::complex<float64x2> cospi (std::complex<float64x2> const &z);
std::complex<float64x2> expm1 (std::complex<float64x2> const &z);
std::complex<float64x2> log2  (std::complex<float64x2> const &z);
std::complex<float64x2> log10 (std::complex<float64x2> const &z);
std::complex<float64x2> log1p (std::complex<float64x2> const &z);

// cis(x) = cos(x) + i·sin(x). One fused sincos covers both components.
// Mirrors libquadmath cexpiq.
inline std::complex<float64x2> cexpi(float64x2 const &x) {
  float64x2 s, c;
  sincos(x, s, c);
  return std::complex<float64x2>(c, s);
}
} // namespace multifloats

#endif  /* __cplusplus */

/* Scope the visibility macro to this header so callers don't get a leaked
 * `MULTIFLOATS_API` identifier. Consumers wanting their own attribute
 * macros should define them in their own namespace. */
#undef MULTIFLOATS_API
