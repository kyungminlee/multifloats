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

#include <cfenv>
#include <charconv>
#include <cmath>
#include <complex>
#include <cstddef>
#include <climits>
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
/// @brief Double-double value: two IEEE-754 binary64 limbs (~104-bit
///        significand). A plain POD in C; carries the full class API in C++.
struct float64x2 {
  /// High and low limbs; `limbs[0]` is the leading binary64, `limbs[1]` the
  /// correction term with `|limbs[1]| <= ulp(limbs[0])/2` when canonical.
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
  /// @brief Construct from a double (low limb zero). Implicit on purpose.
  constexpr float64x2(double arg) : limbs{arg, 0.0} {}
  /// @brief Construct from explicit high and low limbs, `{hi, lo}`.
  constexpr float64x2(double hi, double lo) : limbs{hi, lo} {}

  /// @brief Truncate to the leading limb as a plain `double`.
  constexpr explicit operator double() const { return limbs[0]; }

  // Lexicographic by limb. The low-limb compare MUST be gated on an *equal*
  // high limb: a NaN high limb fails both `<` and `==`, so the operators
  // return false (NaN is unordered) — whereas a bare fall-through to the low
  // limb would let `(NaN, 0) < (x, lo>0)` wrongly return true. For +0 vs -0
  // the high limbs compare equal (IEEE +0 == -0), so we fall through to the
  // low limb. Same ordering for `<= / >=` below.
  /// @brief Lexicographic ordering by limb; NaN is unordered, +0 equals -0.
  constexpr bool operator==(float64x2 const &r) const {
    return limbs[0] == r.limbs[0] && limbs[1] == r.limbs[1];
  }
  constexpr bool operator!=(float64x2 const &r) const { return !(*this == r); }
  constexpr bool operator<(float64x2 const &r) const {
    if (limbs[0] < r.limbs[0]) return true;
    if (limbs[0] == r.limbs[0]) return limbs[1] < r.limbs[1];
    return false;
  }
  constexpr bool operator>(float64x2 const &r) const  { return  (r < *this); }
  // `<=` / `>=` are NOT `!(r < *this)` / `!(*this < r)`: negating `<` turns
  // NaN's unordered-false into a wrong `true`. Compare limbwise and require
  // an *equal* high limb (NaN fails `==`, so it falls through to false) — the
  // ordered IEEE result, at the same cost as `operator<`.
  constexpr bool operator<=(float64x2 const &r) const {
    if (limbs[0] < r.limbs[0]) return true;
    if (limbs[0] == r.limbs[0]) return limbs[1] <= r.limbs[1];
    return false;
  }
  constexpr bool operator>=(float64x2 const &r) const {
    if (r.limbs[0] < limbs[0]) return true;
    if (limbs[0] == r.limbs[0]) return limbs[1] >= r.limbs[1];
    return false;
  }

  /// @brief Unary plus (identity).
  constexpr float64x2 operator+() const { return *this; }
  /// @brief Unary negation, `-a`.
  constexpr float64x2 operator-() const { return {-limbs[0], -limbs[1]}; }

  // ---------------------------------------------------------------------------
  // Binary arithmetic — kernels inlined directly, translated from
  // MultiFloats.jl (mfadd / mfmul) and the Float64x2 division kernel.
  // ---------------------------------------------------------------------------

  /// @brief Double-double sum (compensated, error-free transformation).
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

  /// @brief Double-double difference (compensated, error-free transformation).
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

  /// @brief Double-double product (Dekker `two_prod` + FMA).
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

  /// @brief Double-double quotient (Dekker-style divide with one refinement).
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

  /// @brief Compound add-assign, `*this = *this + rhs`.
  constexpr float64x2 &operator+=(float64x2 const &rhs) {
    return *this = *this + rhs;
  }
  /// @brief Compound subtract-assign, `*this = *this - rhs`.
  constexpr float64x2 &operator-=(float64x2 const &rhs) {
    return *this = *this - rhs;
  }
  /// @brief Compound multiply-assign, `*this = *this * rhs`.
  constexpr float64x2 &operator*=(float64x2 const &rhs) {
    return *this = *this * rhs;
  }
  /// @brief Compound divide-assign, `*this = *this / rhs`.
  constexpr float64x2 &operator/=(float64x2 const &rhs) {
    return *this = *this / rhs;
  }
#endif /* __cplusplus */
};

/* complex64x2 — double-double complex, identical POD in C and C++.
 *
 * A distinct struct (not a typedef for `std::complex<float64x2>`) so that
 * LTO sees one and only one type identity for this struct at every C, C++,
 * and Fortran `bind(c)` ABI boundary. Keeping `std::complex<float64x2>` as
 * the authoritative kernel representation *inside* the library would make
 * IR-level signatures of `c*dd` entry points diverge between the
 * implementation TU and C/Fortran callers — LTO then rejects the link.
 *
 * For C++ callers who prefer `std::complex<float64x2>`, header-inline
 * overloads (declared after the extern "C" block) marshal between the two
 * representations. Layout compatibility is pinned by the static_asserts
 * below; the overloads forward by value and inline away cleanly. */
/// @brief Double-double complex value: a `float64x2` real and imaginary part.
///        A distinct POD (not `std::complex<float64x2>`) for a single ABI
///        type identity across C, C++, and Fortran `bind(c)`.
struct complex64x2 {
  struct float64x2 re;   ///< Real part.
  struct float64x2 im;   ///< Imaginary part.
};

#ifndef __cplusplus
typedef struct float64x2 float64x2;
typedef struct complex64x2 complex64x2;
#endif

/* Guard against surprise padding — the entire C ABI and the Fortran
 * iso_c_binding layer assume float64x2 is exactly two back-to-back
 * doubles (and complex64x2 exactly four). C11 / C++11 required. Also
 * pin that `std::complex<float64x2>` matches complex64x2's layout so
 * the C++ overload layer can rely on it for marshaling. */
#ifdef __cplusplus
static_assert(sizeof(float64x2) == 2 * sizeof(double),
              "float64x2 must be two back-to-back doubles with no padding");
static_assert(sizeof(complex64x2) == 4 * sizeof(double),
              "complex64x2 must be four back-to-back doubles with no padding");
static_assert(sizeof(std::complex<float64x2>) == 4 * sizeof(double),
              "std::complex<float64x2> must be four back-to-back doubles with no padding");
static_assert(std::is_standard_layout<complex64x2>::value,
              "complex64x2 must be standard-layout (ABI interchange)");
static_assert(std::is_trivially_copyable<complex64x2>::value,
              "complex64x2 must be trivially copyable (ABI interchange)");
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
/// @brief Sum of two double-doubles, `a + b`.
MULTIFLOATS_API float64x2 adddd(float64x2 a, float64x2 b);
/// @brief Difference of two double-doubles, `a - b`.
MULTIFLOATS_API float64x2 subdd(float64x2 a, float64x2 b);
/// @brief Product of two double-doubles, `a * b`.
MULTIFLOATS_API float64x2 muldd(float64x2 a, float64x2 b);
/// @brief Quotient of two double-doubles, `a / b`.
MULTIFLOATS_API float64x2 divdd(float64x2 a, float64x2 b);

/* Unary */
/// @brief Negation, `-a`.
MULTIFLOATS_API float64x2 negdd(float64x2 a);
/// @brief Absolute value, `|a|`.
MULTIFLOATS_API float64x2 fabsdd(float64x2 a);
/// @brief Square root.
MULTIFLOATS_API float64x2 sqrtdd(float64x2 a);

/* Rounding. truncdd matches C trunc / Fortran AINT (toward zero).
 * rounddd matches C round / Fortran ANINT (to nearest, halfway away). */
/// @brief Round toward zero (C `trunc` / Fortran `AINT`).
MULTIFLOATS_API float64x2 truncdd(float64x2 a);
/// @brief Round to nearest, halfway away from zero (C `round` / Fortran `ANINT`).
MULTIFLOATS_API float64x2 rounddd(float64x2 a);

/* Binary */
/// @brief Minimum of two double-doubles.
MULTIFLOATS_API float64x2 fmindd(float64x2 a, float64x2 b);
/// @brief Maximum of two double-doubles.
MULTIFLOATS_API float64x2 fmaxdd(float64x2 a, float64x2 b);
/// @brief Euclidean distance `sqrt(a^2 + b^2)`, with overflow-safe scaling.
MULTIFLOATS_API float64x2 hypotdd(float64x2 a, float64x2 b);
/// @brief Power `a**b`, via `exp(b * log(a))`.
MULTIFLOATS_API float64x2 powdd(float64x2 a, float64x2 b);
/* Integer-exponent power via exponentiation-by-squaring: ⌈log2|n|⌉
 * squarings plus popcount(|n|) multiplies, exact up to DD precision
 * per step. One DD division at the end when n < 0. */
/// @brief Integer power `a**n` by exponentiation-by-squaring.
MULTIFLOATS_API float64x2 powidd(float64x2 a, int n);
/// @brief Floating-point remainder of `a / b` (sign of `a`; C `fmod`).
MULTIFLOATS_API float64x2 fmoddd(float64x2 a, float64x2 b);
/* Floored modulo (Fortran `modulo` semantics): same sign as y.
 * modulodd(x, y) = fmod(x, y), plus y if the remainder and y have
 * opposite signs. Uses DD signbit (first-nonzero-limb), so
 * non-canonical DDs with hi==0 still get the right sign. */
/// @brief Floored remainder of `a / b` (sign of `b`; Fortran `modulo`).
MULTIFLOATS_API float64x2 modulodd(float64x2 a, float64x2 b);
/// @brief Positive difference, `max(a - b, 0)`.
MULTIFLOATS_API float64x2 fdimdd(float64x2 a, float64x2 b);
/// @brief Magnitude of `a` with the sign of `b`.
MULTIFLOATS_API float64x2 copysigndd(float64x2 a, float64x2 b);
/// @brief Fused multiply-add `a * b + c` with a single rounding.
MULTIFLOATS_API float64x2 fmadd(float64x2 a, float64x2 b, float64x2 c);

/* Exponential / logarithmic. `expm1dd` / `log1pdd` are the
 * cancellation-safe variants for arguments near zero — direct Taylor
 * / atanh-narrow kernels cover |x| below a threshold, larger |x|
 * falls through to the standard `exp` / `log` paths. */
/// @brief Base-e exponential, `e**a`.
MULTIFLOATS_API float64x2 expdd(float64x2 a);
/// @brief Base-2 exponential, `2**a`.
MULTIFLOATS_API float64x2 exp2dd(float64x2 a);
/// @brief Base-10 exponential, `10**a`.
MULTIFLOATS_API float64x2 exp10dd(float64x2 a);
/// @brief `e**a - 1`, accurate for small `a`.
MULTIFLOATS_API float64x2 expm1dd(float64x2 a);
/// @brief `2**a - 1`, accurate for small `a`.
MULTIFLOATS_API float64x2 exp2m1dd(float64x2 a);
/// @brief `10**a - 1`, accurate for small `a`.
MULTIFLOATS_API float64x2 exp10m1dd(float64x2 a);
/// @brief Natural (base-e) logarithm.
MULTIFLOATS_API float64x2 logdd(float64x2 a);
/// @brief Base-2 logarithm.
MULTIFLOATS_API float64x2 log2dd(float64x2 a);
/// @brief Base-10 logarithm.
MULTIFLOATS_API float64x2 log10dd(float64x2 a);
/// @brief `log(1 + a)`, accurate for small `a`.
MULTIFLOATS_API float64x2 log1pdd(float64x2 a);
/// @brief `log2(1 + a)`, accurate for small `a`.
MULTIFLOATS_API float64x2 log2p1dd(float64x2 a);
/// @brief `log10(1 + a)`, accurate for small `a`.
MULTIFLOATS_API float64x2 log10p1dd(float64x2 a);

/* Trigonometric */
/// @brief Sine.
MULTIFLOATS_API float64x2 sindd(float64x2 a);
/// @brief Cosine.
MULTIFLOATS_API float64x2 cosdd(float64x2 a);
/// @brief Tangent.
MULTIFLOATS_API float64x2 tandd(float64x2 a);
/// @brief Arc sine.
MULTIFLOATS_API float64x2 asindd(float64x2 a);
/// @brief Arc cosine.
MULTIFLOATS_API float64x2 acosdd(float64x2 a);
/// @brief Arc tangent.
MULTIFLOATS_API float64x2 atandd(float64x2 a);
/// @brief Arc tangent of `a/b`, using the signs of both to pick the quadrant.
MULTIFLOATS_API float64x2 atan2dd(float64x2 a, float64x2 b);

/* π-scaled trig: {sin,cos,tan}pidd(x) = fn(π·x),
 *                {asin,acos,atan}pidd(x) = fn(x)/π,
 *                atan2pidd(y,x) = atan2(y,x)/π. */
/// @brief `sin(pi * a)`.
MULTIFLOATS_API float64x2 sinpidd(float64x2 a);
/// @brief `cos(pi * a)`.
MULTIFLOATS_API float64x2 cospidd(float64x2 a);
/// @brief `tan(pi * a)`.
MULTIFLOATS_API float64x2 tanpidd(float64x2 a);
/// @brief `asin(a) / pi`.
MULTIFLOATS_API float64x2 asinpidd(float64x2 a);
/// @brief `acos(a) / pi`.
MULTIFLOATS_API float64x2 acospidd(float64x2 a);
/// @brief `atan(a) / pi`.
MULTIFLOATS_API float64x2 atanpidd(float64x2 a);
/// @brief `atan2(a, b) / pi`.
MULTIFLOATS_API float64x2 atan2pidd(float64x2 a, float64x2 b);

/* Hyperbolic */
/// @brief Hyperbolic sine.
MULTIFLOATS_API float64x2 sinhdd(float64x2 a);
/// @brief Hyperbolic cosine.
MULTIFLOATS_API float64x2 coshdd(float64x2 a);
/// @brief Hyperbolic tangent.
MULTIFLOATS_API float64x2 tanhdd(float64x2 a);
/// @brief Inverse hyperbolic sine.
MULTIFLOATS_API float64x2 asinhdd(float64x2 a);
/// @brief Inverse hyperbolic cosine.
MULTIFLOATS_API float64x2 acoshdd(float64x2 a);
/// @brief Inverse hyperbolic tangent.
MULTIFLOATS_API float64x2 atanhdd(float64x2 a);

/* Error functions. erfcxdd is the scaled complementary error function
 *   erfcxdd(x) = exp(x^2) * erfcdd(x)
 * (standard name in Faddeeva/Julia/SciPy; Fortran calls it erfc_scaled). */
/// @brief Error function.
MULTIFLOATS_API float64x2 erfdd(float64x2 a);
/// @brief Complementary error function, `1 - erf(a)`.
MULTIFLOATS_API float64x2 erfcdd(float64x2 a);
/// @brief Scaled complementary error function `exp(a^2) * erfc(a)` (Fortran `erfc_scaled`).
MULTIFLOATS_API float64x2 erfcxdd(float64x2 a);

/* Gamma functions */
/// @brief Gamma function, `Gamma(a)`.
MULTIFLOATS_API float64x2 tgammadd(float64x2 a);
/// @brief Natural log of `|Gamma(a)|`.
MULTIFLOATS_API float64x2 lgammadd(float64x2 a);

/* Bessel functions (POSIX naming). jndd / yndd: integer-order recurrence
 * seeded from j0dd / j1dd (Miller's backward recurrence for Jn when n > x).
 * Range variants fill `out[0..n2−n1]` with the order-n1 through order-n2
 * values; yndd_range uses a single stable forward-recurrence sweep (cheaper
 * than n2−n1+1 independent yndd calls), while jndd_range loops over jndd
 * because Jn's forward recurrence is unstable for n > x. */
/// @brief Bessel function of the first kind, order 0.
MULTIFLOATS_API float64x2 j0dd(float64x2 a);
/// @brief Bessel function of the first kind, order 1.
MULTIFLOATS_API float64x2 j1dd(float64x2 a);
/// @brief Bessel function of the second kind, order 0.
MULTIFLOATS_API float64x2 y0dd(float64x2 a);
/// @brief Bessel function of the second kind, order 1.
MULTIFLOATS_API float64x2 y1dd(float64x2 a);
/// @brief Bessel function of the first kind, integer order `n`.
MULTIFLOATS_API float64x2 jndd(int n, float64x2 a);
/// @brief Bessel function of the second kind, integer order `n`.
MULTIFLOATS_API float64x2 yndd(int n, float64x2 a);
/// @brief Fill `out` with first-kind Bessel values for orders `n1..n2`.
MULTIFLOATS_API void jndd_range(int n1, int n2, float64x2 a, float64x2 *out);
/// @brief Fill `out` with second-kind Bessel values for orders `n1..n2` (stable forward recurrence).
MULTIFLOATS_API void yndd_range(int n1, int n2, float64x2 a, float64x2 *out);

/* Fused sincos / sinhcosh. One range-reduction + Taylor pair produces
 * both outputs, roughly halving the transcendental cost for call sites
 * that need both. Out-pointer style since C has no multi-value return. */
/// @brief Compute `sin(a)` and `cos(a)` together into `*s` and `*c`.
MULTIFLOATS_API void sincosdd(float64x2 a, float64x2 *s, float64x2 *c);
/// @brief Compute `sinh(a)` and `cosh(a)` together into `*s` and `*c`.
MULTIFLOATS_API void sinhcoshdd(float64x2 a, float64x2 *s, float64x2 *c);

/* Complex DD arithmetic. These are the canonical implementations; the
 * Fortran elemental `cdd + cdd`, `cdd * dp`, etc. routines wrap these via
 * bind(c). */
/// @brief Complex sum of two double-double complex values.
MULTIFLOATS_API complex64x2 cadddd(complex64x2 a, complex64x2 b);
/// @brief Complex difference of two double-double complex values.
MULTIFLOATS_API complex64x2 csubdd(complex64x2 a, complex64x2 b);
/// @brief Complex product of two double-double complex values.
MULTIFLOATS_API complex64x2 cmuldd(complex64x2 a, complex64x2 b);
/// @brief Complex quotient of two double-double complex values.
MULTIFLOATS_API complex64x2 cdivdd(complex64x2 a, complex64x2 b);

/* Complex DD transcendentals. Branch cuts match C99 Annex G (matching
 * libquadmath cexpq/clogq/csqrtq/...). Where the classic formula needs
 * both sin(y) and cos(y) or both sinh(y) and cosh(y), these use the fused
 * kernels internally so one range-reduction + Taylor pair covers both.
 * Precision and speed have been measured; the std::complex<float64x2>
 * specializations below delegate here where specialization wins (exp, sin,
 * cos, tan, sinh, cosh, tanh, atanh, acos). */
/// @brief Complex base-e exponential.
MULTIFLOATS_API complex64x2 cexpdd(complex64x2 z);
/// @brief Complex `exp(z) - 1`, accurate for small `z`.
MULTIFLOATS_API complex64x2 cexpm1dd(complex64x2 z);
/// @brief Complex natural logarithm (principal branch).
MULTIFLOATS_API complex64x2 clogdd(complex64x2 z);
/// @brief Complex base-2 logarithm.
MULTIFLOATS_API complex64x2 clog2dd(complex64x2 z);
/// @brief Complex base-10 logarithm.
MULTIFLOATS_API complex64x2 clog10dd(complex64x2 z);
/// @brief Complex `log(1 + z)`, accurate for small `z`.
MULTIFLOATS_API complex64x2 clog1pdd(complex64x2 z);
/// @brief Complex power, `z**w`.
MULTIFLOATS_API complex64x2 cpowdd(complex64x2 z, complex64x2 w);
/// @brief Complex square root (principal branch).
MULTIFLOATS_API complex64x2 csqrtdd(complex64x2 z);
/// @brief Complex sine.
MULTIFLOATS_API complex64x2 csindd(complex64x2 z);
/// @brief Complex `sin(pi * z)`.
MULTIFLOATS_API complex64x2 csinpidd(complex64x2 z);
/// @brief Complex cosine.
MULTIFLOATS_API complex64x2 ccosdd(complex64x2 z);
/// @brief Complex `cos(pi * z)`.
MULTIFLOATS_API complex64x2 ccospidd(complex64x2 z);
/// @brief Complex tangent.
MULTIFLOATS_API complex64x2 ctandd(complex64x2 z);
/// @brief Complex arc sine.
MULTIFLOATS_API complex64x2 casindd(complex64x2 z);
/// @brief Complex arc cosine.
MULTIFLOATS_API complex64x2 cacosdd(complex64x2 z);
/// @brief Complex arc tangent.
MULTIFLOATS_API complex64x2 catandd(complex64x2 z);
/// @brief Complex hyperbolic sine.
MULTIFLOATS_API complex64x2 csinhdd(complex64x2 z);
/// @brief Complex hyperbolic cosine.
MULTIFLOATS_API complex64x2 ccoshdd(complex64x2 z);
/// @brief Complex hyperbolic tangent.
MULTIFLOATS_API complex64x2 ctanhdd(complex64x2 z);
/// @brief Complex inverse hyperbolic sine.
MULTIFLOATS_API complex64x2 casinhdd(complex64x2 z);
/// @brief Complex inverse hyperbolic cosine.
MULTIFLOATS_API complex64x2 cacoshdd(complex64x2 z);
/// @brief Complex inverse hyperbolic tangent.
MULTIFLOATS_API complex64x2 catanhdd(complex64x2 z);

/* Complex magnitude / argument / projection / conjugate / accessors.
 * cabsdd = |z|  (overflow-safe hypot).
 * cargdd = arg(z) in (-pi, pi].
 * cprojdd = projection onto the Riemann sphere (C99 7.3.9.4):
 *   infinities collapse to (+inf, copysign(0, imag)), else identity. */
/// @brief Complex magnitude `|z|` (overflow-safe).
MULTIFLOATS_API float64x2   cabsdd(complex64x2 z);
/// @brief Complex argument `arg(z)` in (-pi, pi].
MULTIFLOATS_API float64x2   cargdd(complex64x2 z);
/// @brief Projection of `z` onto the Riemann sphere (C99 `cproj`).
MULTIFLOATS_API complex64x2 cprojdd(complex64x2 z);
/// @brief Complex conjugate.
MULTIFLOATS_API complex64x2 conjdd(complex64x2 z);
/// @brief Real part.
MULTIFLOATS_API float64x2   crealdd(complex64x2 z);
/// @brief Imaginary part.
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
 * doc/dev/architecture.md.
 *
 * renorm_interval: if > 0, renormalize accumulators every N reductions
 * (matches DD_FMA_RENORM_INTERVAL in the Fortran layer — keeps s_lo
 * bounded for large k). Pass 0 to renormalize only at the end. */
/// @brief Column-major matrix product `C(m,n) = A(m,k) * B(k,n)`. `renorm_interval` bounds low-limb growth (0 renormalizes only at the end).
MULTIFLOATS_API void matmuldd_mm(const float64x2 *a, const float64x2 *b,
                         float64x2 *c,
                         int64_t m, int64_t k, int64_t n,
                         int64_t renorm_interval);
/// @brief Column-major matrix-vector product `y(m) = A(m,k) * x(k)`. See `renorm_interval`.
MULTIFLOATS_API void matmuldd_mv(const float64x2 *a, const float64x2 *x,
                         float64x2 *y,
                         int64_t m, int64_t k,
                         int64_t renorm_interval);
/// @brief Column-major vector-matrix product `y(n) = x(k) * B(k,n)`. See `renorm_interval`.
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
/// @brief Write the scientific-notation decimal form of `x` into `[first, last)` (no NUL), `std::to_chars`-style; returns one past the last byte, or NULL if the buffer is too small.
MULTIFLOATS_API char *to_charsdd(float64x2 x, int precision,
                                 char *first, char *last);

/* Comparison (return int: 1 = true, 0 = false) */
/// @brief Nonzero if `a == b`.
MULTIFLOATS_API int eqdd(float64x2 a, float64x2 b);
/// @brief Nonzero if `a != b`.
MULTIFLOATS_API int nedd(float64x2 a, float64x2 b);
/// @brief Nonzero if `a < b`.
MULTIFLOATS_API int ltdd(float64x2 a, float64x2 b);
/// @brief Nonzero if `a <= b`.
MULTIFLOATS_API int ledd(float64x2 a, float64x2 b);
/// @brief Nonzero if `a > b`.
MULTIFLOATS_API int gtdd(float64x2 a, float64x2 b);
/// @brief Nonzero if `a >= b`.
MULTIFLOATS_API int gedd(float64x2 a, float64x2 b);

/* libquadmath parity — every NAMEq in <quadmath.h> has a NAMEdd entry
 * here. Implementations are the canonical `multifloats::NAME` C++
 * helpers in the header section below; the C-ABI symbols in
 * src/float64x2/abi.inc are thin marshaling shims. */

/* Cube root, integer rounding (toward ±inf, nearest-current-mode), exponent
 * split / rescale, integer/fractional split, ulp-step, IEEE remainder. */
/// @brief Cube root.
MULTIFLOATS_API float64x2 cbrtdd(float64x2 a);
/// @brief Round toward +infinity.
MULTIFLOATS_API float64x2 ceildd(float64x2 a);
/// @brief Round toward -infinity.
MULTIFLOATS_API float64x2 floordd(float64x2 a);
/// @brief Round to integer in the current rounding mode, without raising inexact.
MULTIFLOATS_API float64x2 nearbyintdd(float64x2 a);
/// @brief Round to integer in the current rounding mode.
MULTIFLOATS_API float64x2 rintdd(float64x2 a);
/// @brief Unbiased radix-2 exponent of `a`, as a floating-point value.
MULTIFLOATS_API float64x2 logbdd(float64x2 a);
/// @brief Split `a` into a fraction in [0.5, 1) and an exponent stored in `*exp`.
MULTIFLOATS_API float64x2 frexpdd(float64x2 a, int *exp);
/// @brief Split `a` into integer part (stored in `*iptr`) and fractional part.
MULTIFLOATS_API float64x2 modfdd(float64x2 a, float64x2 *iptr);
/// @brief Scale by a power of two, `a * 2**n`.
MULTIFLOATS_API float64x2 ldexpdd(float64x2 a, int n);
/// @brief Scale by a power of two, `a * 2**n`.
MULTIFLOATS_API float64x2 scalbndd(float64x2 a, int n);
/// @brief Scale by a power of two, `a * 2**n` (`long` exponent).
MULTIFLOATS_API float64x2 scalblndd(float64x2 a, long n);
/// @brief Next representable double-double after `a` toward `b`.
MULTIFLOATS_API float64x2 nextafterdd(float64x2 a, float64x2 b);
/// @brief IEEE 754 remainder of `a / b` (round-to-nearest quotient).
MULTIFLOATS_API float64x2 remainderdd(float64x2 a, float64x2 b);
/// @brief IEEE 754 remainder of `a / b`; low quotient bits stored in `*quo`.
MULTIFLOATS_API float64x2 remquodd(float64x2 a, float64x2 b, int *quo);
/// @brief Unbiased radix-2 exponent of `a`, as an `int`.
MULTIFLOATS_API int ilogbdd(float64x2 a);

/* Integer rounding — leading-limb's libm result with a half-integer
 * fixup that consults the lo limb (so e.g. lround((0.5, -ε)) = 0). */
/// @brief Round to nearest, halfway away from zero, returning `long`.
MULTIFLOATS_API long lrounddd(float64x2 a);
/// @brief Round to nearest, halfway away from zero, returning `long long`.
MULTIFLOATS_API long long llrounddd(float64x2 a);
/// @brief Round to integer in the current rounding mode, returning `long`.
MULTIFLOATS_API long lrintdd(float64x2 a);
/// @brief Round to integer in the current rounding mode, returning `long long`.
MULTIFLOATS_API long long llrintdd(float64x2 a);

/* Classification — C99 isnan/isinf/signbit/isfinite/isnormal contract:
 * non-zero for true. `finitedd` mirrors libquadmath's legacy `finiteq`
 * spelling (alias for isfinitedd). `fpclassifydd` returns one of the
 * C99 FP_* values from <math.h> applied to the leading limb. */
/// @brief Nonzero if `a` is NaN.
MULTIFLOATS_API int isnandd(float64x2 a);
/// @brief Nonzero if `a` is infinite.
MULTIFLOATS_API int isinfdd(float64x2 a);
/// @brief Nonzero if `a` is finite.
MULTIFLOATS_API int isfinitedd(float64x2 a);
/// @brief Nonzero if `a` is finite (libquadmath `finiteq` spelling).
MULTIFLOATS_API int finitedd(float64x2 a);
/// @brief Nonzero if `a` is normal (finite, nonzero, not subnormal).
MULTIFLOATS_API int isnormaldd(float64x2 a);
/// @brief Nonzero if the sign bit of `a` is set.
MULTIFLOATS_API int signbitdd(float64x2 a);
/// @brief C99 floating-point classification of `a` (`FP_NAN`, `FP_INFINITE`, ...).
MULTIFLOATS_API int fpclassifydd(float64x2 a);

/* C99 ordered/unordered comparison predicates. `isunordereddd` returns
 * non-zero iff either argument is NaN; the others return non-zero iff
 * the named ordered relation holds (and silently return 0 for NaN
 * operands, matching the C99 macros). */
/// @brief Nonzero if `a` or `b` is NaN.
MULTIFLOATS_API int isunordereddd(float64x2 a, float64x2 b);
/// @brief Nonzero if `a > b` (quiet for NaN operands).
MULTIFLOATS_API int isgreaterdd(float64x2 a, float64x2 b);
/// @brief Nonzero if `a >= b` (quiet for NaN operands).
MULTIFLOATS_API int isgreaterequaldd(float64x2 a, float64x2 b);
/// @brief Nonzero if `a < b` (quiet for NaN operands).
MULTIFLOATS_API int islessdd(float64x2 a, float64x2 b);
/// @brief Nonzero if `a <= b` (quiet for NaN operands).
MULTIFLOATS_API int islessequaldd(float64x2 a, float64x2 b);
/// @brief Nonzero if `a < b` or `a > b` (quiet for NaN operands).
MULTIFLOATS_API int islessgreaterdd(float64x2 a, float64x2 b);
/* Signaling-NaN test on the leading limb. Uses the compiler builtin
 * where available (gcc >= 13, clang with __has_builtin); otherwise
 * returns 0 — multifloats never constructs sNaNs internally, so the
 * fallback only loses signal when the caller explicitly passed an
 * sNaN through the C ABI. */
/// @brief Nonzero if `a` is a signaling NaN.
MULTIFLOATS_API int issignalingdd(float64x2 a);

/* NaN with payload tag (libquadmath nanq parity). Calls std::nan(tagp)
 * for the leading limb; lo is set to 0. */
/// @brief Quiet NaN with the payload parsed from `tagp` (libquadmath `nanq`).
MULTIFLOATS_API float64x2 nandd(const char *tagp);

/* cis(x) = cos(x) + i sin(x). Mirrors libquadmath cexpiq. */
/// @brief `cos(a) + i*sin(a)` (libquadmath `cexpiq`).
MULTIFLOATS_API complex64x2 cexpidd(float64x2 a);

/* ----------------------------------------------------------------------------
 * C23-style `f64x2` aliases. Every libquadmath-style `*dd` entry point above
 * has a strong-symbol alias spelled with the `f64x2` type tag (e.g.
 * `sinf64x2`, `cabsf64x2`). The two names always refer to the same symbol;
 * pick whichever convention you prefer.
 *   - `dd`     — short, matches double-double literature and libquadmath's
 *                `q` style (e.g. `sindd`, `cabsdd`).
 *   - `f64x2`  — mirrors C23 / glibc's `_Float128` → `f128` suffix
 *                convention but truthfully says "two binary64 limbs"
 *                rather than (mis)labelling DD as quad precision.
 * -------------------------------------------------------------------------- */

/// @brief Sum of two double-doubles, `a + b`. (`f64x2`-tagged spelling of `adddd`)
MULTIFLOATS_API float64x2 addf64x2(float64x2 a, float64x2 b);
/// @brief Difference of two double-doubles, `a - b`. (`f64x2`-tagged spelling of `subdd`)
MULTIFLOATS_API float64x2 subf64x2(float64x2 a, float64x2 b);
/// @brief Product of two double-doubles, `a * b`. (`f64x2`-tagged spelling of `muldd`)
MULTIFLOATS_API float64x2 mulf64x2(float64x2 a, float64x2 b);
/// @brief Quotient of two double-doubles, `a / b`. (`f64x2`-tagged spelling of `divdd`)
MULTIFLOATS_API float64x2 divf64x2(float64x2 a, float64x2 b);
/// @brief Negation, `-a`. (`f64x2`-tagged spelling of `negdd`)
MULTIFLOATS_API float64x2 negf64x2(float64x2 a);
/// @brief Absolute value, `|a|`. (`f64x2`-tagged spelling of `fabsdd`)
MULTIFLOATS_API float64x2 fabsf64x2(float64x2 a);
/// @brief Square root. (`f64x2`-tagged spelling of `sqrtdd`)
MULTIFLOATS_API float64x2 sqrtf64x2(float64x2 a);
/// @brief Round toward zero (C `trunc` / Fortran `AINT`). (`f64x2`-tagged spelling of `truncdd`)
MULTIFLOATS_API float64x2 truncf64x2(float64x2 a);
/// @brief Round to nearest, halfway away from zero (C `round` / Fortran `ANINT`). (`f64x2`-tagged spelling of `rounddd`)
MULTIFLOATS_API float64x2 roundf64x2(float64x2 a);
/// @brief Minimum of two double-doubles. (`f64x2`-tagged spelling of `fmindd`)
MULTIFLOATS_API float64x2 fminf64x2(float64x2 a, float64x2 b);
/// @brief Maximum of two double-doubles. (`f64x2`-tagged spelling of `fmaxdd`)
MULTIFLOATS_API float64x2 fmaxf64x2(float64x2 a, float64x2 b);
/// @brief Euclidean distance `sqrt(a^2 + b^2)`, with overflow-safe scaling. (`f64x2`-tagged spelling of `hypotdd`)
MULTIFLOATS_API float64x2 hypotf64x2(float64x2 a, float64x2 b);
/// @brief Power `a**b`, via `exp(b * log(a))`. (`f64x2`-tagged spelling of `powdd`)
MULTIFLOATS_API float64x2 powf64x2(float64x2 a, float64x2 b);
/// @brief Integer power `a**n` by exponentiation-by-squaring. (`f64x2`-tagged spelling of `powidd`)
MULTIFLOATS_API float64x2 powif64x2(float64x2 a, int n);
/// @brief Floating-point remainder of `a / b` (sign of `a`; C `fmod`). (`f64x2`-tagged spelling of `fmoddd`)
MULTIFLOATS_API float64x2 fmodf64x2(float64x2 a, float64x2 b);
/// @brief Floored remainder of `a / b` (sign of `b`; Fortran `modulo`). (`f64x2`-tagged spelling of `modulodd`)
MULTIFLOATS_API float64x2 modulof64x2(float64x2 a, float64x2 b);
/// @brief Positive difference, `max(a - b, 0)`. (`f64x2`-tagged spelling of `fdimdd`)
MULTIFLOATS_API float64x2 fdimf64x2(float64x2 a, float64x2 b);
/// @brief Magnitude of `a` with the sign of `b`. (`f64x2`-tagged spelling of `copysigndd`)
MULTIFLOATS_API float64x2 copysignf64x2(float64x2 a, float64x2 b);
/// @brief Fused multiply-add `a * b + c` with a single rounding. (`f64x2`-tagged spelling of `fmadd`)
MULTIFLOATS_API float64x2 fmaf64x2(float64x2 a, float64x2 b, float64x2 c);
/// @brief Base-e exponential, `e**a`. (`f64x2`-tagged spelling of `expdd`)
MULTIFLOATS_API float64x2 expf64x2(float64x2 a);
/// @brief Base-2 exponential, `2**a`. (`f64x2`-tagged spelling of `exp2dd`)
MULTIFLOATS_API float64x2 exp2f64x2(float64x2 a);
/// @brief Base-10 exponential, `10**a`. (`f64x2`-tagged spelling of `exp10dd`)
MULTIFLOATS_API float64x2 exp10f64x2(float64x2 a);
/// @brief `e**a - 1`, accurate for small `a`. (`f64x2`-tagged spelling of `expm1dd`)
MULTIFLOATS_API float64x2 expm1f64x2(float64x2 a);
/// @brief `2**a - 1`, accurate for small `a`. (`f64x2`-tagged spelling of `exp2m1dd`)
MULTIFLOATS_API float64x2 exp2m1f64x2(float64x2 a);
/// @brief `10**a - 1`, accurate for small `a`. (`f64x2`-tagged spelling of `exp10m1dd`)
MULTIFLOATS_API float64x2 exp10m1f64x2(float64x2 a);
/// @brief Natural (base-e) logarithm. (`f64x2`-tagged spelling of `logdd`)
MULTIFLOATS_API float64x2 logf64x2(float64x2 a);
/// @brief Base-2 logarithm. (`f64x2`-tagged spelling of `log2dd`)
MULTIFLOATS_API float64x2 log2f64x2(float64x2 a);
/// @brief Base-10 logarithm. (`f64x2`-tagged spelling of `log10dd`)
MULTIFLOATS_API float64x2 log10f64x2(float64x2 a);
/// @brief `log(1 + a)`, accurate for small `a`. (`f64x2`-tagged spelling of `log1pdd`)
MULTIFLOATS_API float64x2 log1pf64x2(float64x2 a);
/// @brief `log2(1 + a)`, accurate for small `a`. (`f64x2`-tagged spelling of `log2p1dd`)
MULTIFLOATS_API float64x2 log2p1f64x2(float64x2 a);
/// @brief `log10(1 + a)`, accurate for small `a`. (`f64x2`-tagged spelling of `log10p1dd`)
MULTIFLOATS_API float64x2 log10p1f64x2(float64x2 a);
/// @brief Sine. (`f64x2`-tagged spelling of `sindd`)
MULTIFLOATS_API float64x2 sinf64x2(float64x2 a);
/// @brief Cosine. (`f64x2`-tagged spelling of `cosdd`)
MULTIFLOATS_API float64x2 cosf64x2(float64x2 a);
/// @brief Tangent. (`f64x2`-tagged spelling of `tandd`)
MULTIFLOATS_API float64x2 tanf64x2(float64x2 a);
/// @brief Arc sine. (`f64x2`-tagged spelling of `asindd`)
MULTIFLOATS_API float64x2 asinf64x2(float64x2 a);
/// @brief Arc cosine. (`f64x2`-tagged spelling of `acosdd`)
MULTIFLOATS_API float64x2 acosf64x2(float64x2 a);
/// @brief Arc tangent. (`f64x2`-tagged spelling of `atandd`)
MULTIFLOATS_API float64x2 atanf64x2(float64x2 a);
/// @brief Arc tangent of `a/b`, using the signs of both to pick the quadrant. (`f64x2`-tagged spelling of `atan2dd`)
MULTIFLOATS_API float64x2 atan2f64x2(float64x2 a, float64x2 b);
/// @brief `sin(pi * a)`. (`f64x2`-tagged spelling of `sinpidd`)
MULTIFLOATS_API float64x2 sinpif64x2(float64x2 a);
/// @brief `cos(pi * a)`. (`f64x2`-tagged spelling of `cospidd`)
MULTIFLOATS_API float64x2 cospif64x2(float64x2 a);
/// @brief `tan(pi * a)`. (`f64x2`-tagged spelling of `tanpidd`)
MULTIFLOATS_API float64x2 tanpif64x2(float64x2 a);
/// @brief `asin(a) / pi`. (`f64x2`-tagged spelling of `asinpidd`)
MULTIFLOATS_API float64x2 asinpif64x2(float64x2 a);
/// @brief `acos(a) / pi`. (`f64x2`-tagged spelling of `acospidd`)
MULTIFLOATS_API float64x2 acospif64x2(float64x2 a);
/// @brief `atan(a) / pi`. (`f64x2`-tagged spelling of `atanpidd`)
MULTIFLOATS_API float64x2 atanpif64x2(float64x2 a);
/// @brief `atan2(a, b) / pi`. (`f64x2`-tagged spelling of `atan2pidd`)
MULTIFLOATS_API float64x2 atan2pif64x2(float64x2 a, float64x2 b);
/// @brief Hyperbolic sine. (`f64x2`-tagged spelling of `sinhdd`)
MULTIFLOATS_API float64x2 sinhf64x2(float64x2 a);
/// @brief Hyperbolic cosine. (`f64x2`-tagged spelling of `coshdd`)
MULTIFLOATS_API float64x2 coshf64x2(float64x2 a);
/// @brief Hyperbolic tangent. (`f64x2`-tagged spelling of `tanhdd`)
MULTIFLOATS_API float64x2 tanhf64x2(float64x2 a);
/// @brief Inverse hyperbolic sine. (`f64x2`-tagged spelling of `asinhdd`)
MULTIFLOATS_API float64x2 asinhf64x2(float64x2 a);
/// @brief Inverse hyperbolic cosine. (`f64x2`-tagged spelling of `acoshdd`)
MULTIFLOATS_API float64x2 acoshf64x2(float64x2 a);
/// @brief Inverse hyperbolic tangent. (`f64x2`-tagged spelling of `atanhdd`)
MULTIFLOATS_API float64x2 atanhf64x2(float64x2 a);
/// @brief Error function. (`f64x2`-tagged spelling of `erfdd`)
MULTIFLOATS_API float64x2 erff64x2(float64x2 a);
/// @brief Complementary error function, `1 - erf(a)`. (`f64x2`-tagged spelling of `erfcdd`)
MULTIFLOATS_API float64x2 erfcf64x2(float64x2 a);
/// @brief Scaled complementary error function `exp(a^2) * erfc(a)` (Fortran `erfc_scaled`). (`f64x2`-tagged spelling of `erfcxdd`)
MULTIFLOATS_API float64x2 erfcxf64x2(float64x2 a);
/// @brief Gamma function, `Gamma(a)`. (`f64x2`-tagged spelling of `tgammadd`)
MULTIFLOATS_API float64x2 tgammaf64x2(float64x2 a);
/// @brief Natural log of `|Gamma(a)|`. (`f64x2`-tagged spelling of `lgammadd`)
MULTIFLOATS_API float64x2 lgammaf64x2(float64x2 a);
/// @brief Bessel function of the first kind, order 0. (`f64x2`-tagged spelling of `j0dd`)
MULTIFLOATS_API float64x2 j0f64x2(float64x2 a);
/// @brief Bessel function of the first kind, order 1. (`f64x2`-tagged spelling of `j1dd`)
MULTIFLOATS_API float64x2 j1f64x2(float64x2 a);
/// @brief Bessel function of the second kind, order 0. (`f64x2`-tagged spelling of `y0dd`)
MULTIFLOATS_API float64x2 y0f64x2(float64x2 a);
/// @brief Bessel function of the second kind, order 1. (`f64x2`-tagged spelling of `y1dd`)
MULTIFLOATS_API float64x2 y1f64x2(float64x2 a);
/// @brief Bessel function of the first kind, integer order `n`. (`f64x2`-tagged spelling of `jndd`)
MULTIFLOATS_API float64x2 jnf64x2(int n, float64x2 a);
/// @brief Bessel function of the second kind, integer order `n`. (`f64x2`-tagged spelling of `yndd`)
MULTIFLOATS_API float64x2 ynf64x2(int n, float64x2 a);
/// @brief Fill `out` with first-kind Bessel values for orders `n1..n2`. (`f64x2`-tagged spelling of `jndd_range`)
MULTIFLOATS_API void jnf64x2_range(int n1, int n2, float64x2 a, float64x2 *out);
/// @brief Fill `out` with second-kind Bessel values for orders `n1..n2` (stable forward recurrence). (`f64x2`-tagged spelling of `yndd_range`)
MULTIFLOATS_API void ynf64x2_range(int n1, int n2, float64x2 a, float64x2 *out);
/// @brief Compute `sin(a)` and `cos(a)` together into `*s` and `*c`. (`f64x2`-tagged spelling of `sincosdd`)
MULTIFLOATS_API void sincosf64x2(float64x2 a, float64x2 *s, float64x2 *c);
/// @brief Compute `sinh(a)` and `cosh(a)` together into `*s` and `*c`. (`f64x2`-tagged spelling of `sinhcoshdd`)
MULTIFLOATS_API void sinhcoshf64x2(float64x2 a, float64x2 *s, float64x2 *c);
/// @brief Complex sum of two double-double complex values. (`f64x2`-tagged spelling of `cadddd`)
MULTIFLOATS_API complex64x2 caddf64x2(complex64x2 a, complex64x2 b);
/// @brief Complex difference of two double-double complex values. (`f64x2`-tagged spelling of `csubdd`)
MULTIFLOATS_API complex64x2 csubf64x2(complex64x2 a, complex64x2 b);
/// @brief Complex product of two double-double complex values. (`f64x2`-tagged spelling of `cmuldd`)
MULTIFLOATS_API complex64x2 cmulf64x2(complex64x2 a, complex64x2 b);
/// @brief Complex quotient of two double-double complex values. (`f64x2`-tagged spelling of `cdivdd`)
MULTIFLOATS_API complex64x2 cdivf64x2(complex64x2 a, complex64x2 b);
/// @brief Complex base-e exponential. (`f64x2`-tagged spelling of `cexpdd`)
MULTIFLOATS_API complex64x2 cexpf64x2(complex64x2 z);
/// @brief Complex `exp(z) - 1`, accurate for small `z`. (`f64x2`-tagged spelling of `cexpm1dd`)
MULTIFLOATS_API complex64x2 cexpm1f64x2(complex64x2 z);
/// @brief Complex natural logarithm (principal branch). (`f64x2`-tagged spelling of `clogdd`)
MULTIFLOATS_API complex64x2 clogf64x2(complex64x2 z);
/// @brief Complex base-2 logarithm. (`f64x2`-tagged spelling of `clog2dd`)
MULTIFLOATS_API complex64x2 clog2f64x2(complex64x2 z);
/// @brief Complex base-10 logarithm. (`f64x2`-tagged spelling of `clog10dd`)
MULTIFLOATS_API complex64x2 clog10f64x2(complex64x2 z);
/// @brief Complex `log(1 + z)`, accurate for small `z`. (`f64x2`-tagged spelling of `clog1pdd`)
MULTIFLOATS_API complex64x2 clog1pf64x2(complex64x2 z);
/// @brief Complex power, `z**w`. (`f64x2`-tagged spelling of `cpowdd`)
MULTIFLOATS_API complex64x2 cpowf64x2(complex64x2 z, complex64x2 w);
/// @brief Complex square root (principal branch). (`f64x2`-tagged spelling of `csqrtdd`)
MULTIFLOATS_API complex64x2 csqrtf64x2(complex64x2 z);
/// @brief Complex sine. (`f64x2`-tagged spelling of `csindd`)
MULTIFLOATS_API complex64x2 csinf64x2(complex64x2 z);
/// @brief Complex `sin(pi * z)`. (`f64x2`-tagged spelling of `csinpidd`)
MULTIFLOATS_API complex64x2 csinpif64x2(complex64x2 z);
/// @brief Complex cosine. (`f64x2`-tagged spelling of `ccosdd`)
MULTIFLOATS_API complex64x2 ccosf64x2(complex64x2 z);
/// @brief Complex `cos(pi * z)`. (`f64x2`-tagged spelling of `ccospidd`)
MULTIFLOATS_API complex64x2 ccospif64x2(complex64x2 z);
/// @brief Complex tangent. (`f64x2`-tagged spelling of `ctandd`)
MULTIFLOATS_API complex64x2 ctanf64x2(complex64x2 z);
/// @brief Complex arc sine. (`f64x2`-tagged spelling of `casindd`)
MULTIFLOATS_API complex64x2 casinf64x2(complex64x2 z);
/// @brief Complex arc cosine. (`f64x2`-tagged spelling of `cacosdd`)
MULTIFLOATS_API complex64x2 cacosf64x2(complex64x2 z);
/// @brief Complex arc tangent. (`f64x2`-tagged spelling of `catandd`)
MULTIFLOATS_API complex64x2 catanf64x2(complex64x2 z);
/// @brief Complex hyperbolic sine. (`f64x2`-tagged spelling of `csinhdd`)
MULTIFLOATS_API complex64x2 csinhf64x2(complex64x2 z);
/// @brief Complex hyperbolic cosine. (`f64x2`-tagged spelling of `ccoshdd`)
MULTIFLOATS_API complex64x2 ccoshf64x2(complex64x2 z);
/// @brief Complex hyperbolic tangent. (`f64x2`-tagged spelling of `ctanhdd`)
MULTIFLOATS_API complex64x2 ctanhf64x2(complex64x2 z);
/// @brief Complex inverse hyperbolic sine. (`f64x2`-tagged spelling of `casinhdd`)
MULTIFLOATS_API complex64x2 casinhf64x2(complex64x2 z);
/// @brief Complex inverse hyperbolic cosine. (`f64x2`-tagged spelling of `cacoshdd`)
MULTIFLOATS_API complex64x2 cacoshf64x2(complex64x2 z);
/// @brief Complex inverse hyperbolic tangent. (`f64x2`-tagged spelling of `catanhdd`)
MULTIFLOATS_API complex64x2 catanhf64x2(complex64x2 z);
/// @brief Complex magnitude `|z|` (overflow-safe). (`f64x2`-tagged spelling of `cabsdd`)
MULTIFLOATS_API float64x2 cabsf64x2(complex64x2 z);
/// @brief Complex argument `arg(z)` in (-pi, pi]. (`f64x2`-tagged spelling of `cargdd`)
MULTIFLOATS_API float64x2 cargf64x2(complex64x2 z);
/// @brief Projection of `z` onto the Riemann sphere (C99 `cproj`). (`f64x2`-tagged spelling of `cprojdd`)
MULTIFLOATS_API complex64x2 cprojf64x2(complex64x2 z);
/// @brief Complex conjugate. (`f64x2`-tagged spelling of `conjdd`)
MULTIFLOATS_API complex64x2 conjf64x2(complex64x2 z);
/// @brief Real part. (`f64x2`-tagged spelling of `crealdd`)
MULTIFLOATS_API float64x2 crealf64x2(complex64x2 z);
/// @brief Imaginary part. (`f64x2`-tagged spelling of `cimagdd`)
MULTIFLOATS_API float64x2 cimagf64x2(complex64x2 z);
/// @brief Column-major matrix product `C(m,n) = A(m,k) * B(k,n)`. `renorm_interval` bounds low-limb growth (0 renormalizes only at the end). (`f64x2`-tagged spelling of `matmuldd_mm`)
MULTIFLOATS_API void matmulf64x2_mm(const float64x2 *a, const float64x2 *b, float64x2 *c, int64_t m, int64_t k, int64_t n, int64_t renorm_interval);
/// @brief Column-major matrix-vector product `y(m) = A(m,k) * x(k)`. See `renorm_interval`. (`f64x2`-tagged spelling of `matmuldd_mv`)
MULTIFLOATS_API void matmulf64x2_mv(const float64x2 *a, const float64x2 *x, float64x2 *y, int64_t m, int64_t k, int64_t renorm_interval);
/// @brief Column-major vector-matrix product `y(n) = x(k) * B(k,n)`. See `renorm_interval`. (`f64x2`-tagged spelling of `matmuldd_vm`)
MULTIFLOATS_API void matmulf64x2_vm(const float64x2 *x, const float64x2 *b, float64x2 *y, int64_t k, int64_t n, int64_t renorm_interval);
/// @brief Write the scientific-notation decimal form of `x` into `[first, last)` (no NUL), `std::to_chars`-style; returns one past the last byte, or NULL if the buffer is too small. (`f64x2`-tagged spelling of `to_charsdd`)
MULTIFLOATS_API char *to_charsf64x2(float64x2 x, int precision, char *first, char *last);
/// @brief Nonzero if `a == b`. (`f64x2`-tagged spelling of `eqdd`)
MULTIFLOATS_API int eqf64x2(float64x2 a, float64x2 b);
/// @brief Nonzero if `a != b`. (`f64x2`-tagged spelling of `nedd`)
MULTIFLOATS_API int nef64x2(float64x2 a, float64x2 b);
/// @brief Nonzero if `a < b`. (`f64x2`-tagged spelling of `ltdd`)
MULTIFLOATS_API int ltf64x2(float64x2 a, float64x2 b);
/// @brief Nonzero if `a <= b`. (`f64x2`-tagged spelling of `ledd`)
MULTIFLOATS_API int lef64x2(float64x2 a, float64x2 b);
/// @brief Nonzero if `a > b`. (`f64x2`-tagged spelling of `gtdd`)
MULTIFLOATS_API int gtf64x2(float64x2 a, float64x2 b);
/// @brief Nonzero if `a >= b`. (`f64x2`-tagged spelling of `gedd`)
MULTIFLOATS_API int gef64x2(float64x2 a, float64x2 b);
/// @brief Cube root. (`f64x2`-tagged spelling of `cbrtdd`)
MULTIFLOATS_API float64x2 cbrtf64x2(float64x2 a);
/// @brief Round toward +infinity. (`f64x2`-tagged spelling of `ceildd`)
MULTIFLOATS_API float64x2 ceilf64x2(float64x2 a);
/// @brief Round toward -infinity. (`f64x2`-tagged spelling of `floordd`)
MULTIFLOATS_API float64x2 floorf64x2(float64x2 a);
/// @brief Round to integer in the current rounding mode, without raising inexact. (`f64x2`-tagged spelling of `nearbyintdd`)
MULTIFLOATS_API float64x2 nearbyintf64x2(float64x2 a);
/// @brief Round to integer in the current rounding mode. (`f64x2`-tagged spelling of `rintdd`)
MULTIFLOATS_API float64x2 rintf64x2(float64x2 a);
/// @brief Unbiased radix-2 exponent of `a`, as a floating-point value. (`f64x2`-tagged spelling of `logbdd`)
MULTIFLOATS_API float64x2 logbf64x2(float64x2 a);
/// @brief Split `a` into a fraction in [0.5, 1) and an exponent stored in `*exp`. (`f64x2`-tagged spelling of `frexpdd`)
MULTIFLOATS_API float64x2 frexpf64x2(float64x2 a, int *exp);
/// @brief Split `a` into integer part (stored in `*iptr`) and fractional part. (`f64x2`-tagged spelling of `modfdd`)
MULTIFLOATS_API float64x2 modff64x2(float64x2 a, float64x2 *iptr);
/// @brief Scale by a power of two, `a * 2**n`. (`f64x2`-tagged spelling of `ldexpdd`)
MULTIFLOATS_API float64x2 ldexpf64x2(float64x2 a, int n);
/// @brief Scale by a power of two, `a * 2**n`. (`f64x2`-tagged spelling of `scalbndd`)
MULTIFLOATS_API float64x2 scalbnf64x2(float64x2 a, int n);
/// @brief Scale by a power of two, `a * 2**n` (`long` exponent). (`f64x2`-tagged spelling of `scalblndd`)
MULTIFLOATS_API float64x2 scalblnf64x2(float64x2 a, long n);
/// @brief Next representable double-double after `a` toward `b`. (`f64x2`-tagged spelling of `nextafterdd`)
MULTIFLOATS_API float64x2 nextafterf64x2(float64x2 a, float64x2 b);
/// @brief IEEE 754 remainder of `a / b` (round-to-nearest quotient). (`f64x2`-tagged spelling of `remainderdd`)
MULTIFLOATS_API float64x2 remainderf64x2(float64x2 a, float64x2 b);
/// @brief IEEE 754 remainder of `a / b`; low quotient bits stored in `*quo`. (`f64x2`-tagged spelling of `remquodd`)
MULTIFLOATS_API float64x2 remquof64x2(float64x2 a, float64x2 b, int *quo);
/// @brief Unbiased radix-2 exponent of `a`, as an `int`. (`f64x2`-tagged spelling of `ilogbdd`)
MULTIFLOATS_API int ilogbf64x2(float64x2 a);
/// @brief Round to nearest, halfway away from zero, returning `long`. (`f64x2`-tagged spelling of `lrounddd`)
MULTIFLOATS_API long lroundf64x2(float64x2 a);
/// @brief Round to nearest, halfway away from zero, returning `long long`. (`f64x2`-tagged spelling of `llrounddd`)
MULTIFLOATS_API long long llroundf64x2(float64x2 a);
/// @brief Round to integer in the current rounding mode, returning `long`. (`f64x2`-tagged spelling of `lrintdd`)
MULTIFLOATS_API long lrintf64x2(float64x2 a);
/// @brief Round to integer in the current rounding mode, returning `long long`. (`f64x2`-tagged spelling of `llrintdd`)
MULTIFLOATS_API long long llrintf64x2(float64x2 a);
/// @brief Nonzero if `a` is NaN. (`f64x2`-tagged spelling of `isnandd`)
MULTIFLOATS_API int isnanf64x2(float64x2 a);
/// @brief Nonzero if `a` is infinite. (`f64x2`-tagged spelling of `isinfdd`)
MULTIFLOATS_API int isinff64x2(float64x2 a);
/// @brief Nonzero if `a` is finite. (`f64x2`-tagged spelling of `isfinitedd`)
MULTIFLOATS_API int isfinitef64x2(float64x2 a);
/// @brief Nonzero if `a` is finite (libquadmath `finiteq` spelling). (`f64x2`-tagged spelling of `finitedd`)
MULTIFLOATS_API int finitef64x2(float64x2 a);
/// @brief Nonzero if `a` is normal (finite, nonzero, not subnormal). (`f64x2`-tagged spelling of `isnormaldd`)
MULTIFLOATS_API int isnormalf64x2(float64x2 a);
/// @brief Nonzero if the sign bit of `a` is set. (`f64x2`-tagged spelling of `signbitdd`)
MULTIFLOATS_API int signbitf64x2(float64x2 a);
/// @brief C99 floating-point classification of `a` (`FP_NAN`, `FP_INFINITE`, ...). (`f64x2`-tagged spelling of `fpclassifydd`)
MULTIFLOATS_API int fpclassifyf64x2(float64x2 a);
/// @brief Nonzero if `a` or `b` is NaN. (`f64x2`-tagged spelling of `isunordereddd`)
MULTIFLOATS_API int isunorderedf64x2(float64x2 a, float64x2 b);
/// @brief Nonzero if `a > b` (quiet for NaN operands). (`f64x2`-tagged spelling of `isgreaterdd`)
MULTIFLOATS_API int isgreaterf64x2(float64x2 a, float64x2 b);
/// @brief Nonzero if `a >= b` (quiet for NaN operands). (`f64x2`-tagged spelling of `isgreaterequaldd`)
MULTIFLOATS_API int isgreaterequalf64x2(float64x2 a, float64x2 b);
/// @brief Nonzero if `a < b` (quiet for NaN operands). (`f64x2`-tagged spelling of `islessdd`)
MULTIFLOATS_API int islessf64x2(float64x2 a, float64x2 b);
/// @brief Nonzero if `a <= b` (quiet for NaN operands). (`f64x2`-tagged spelling of `islessequaldd`)
MULTIFLOATS_API int islessequalf64x2(float64x2 a, float64x2 b);
/// @brief Nonzero if `a < b` or `a > b` (quiet for NaN operands). (`f64x2`-tagged spelling of `islessgreaterdd`)
MULTIFLOATS_API int islessgreaterf64x2(float64x2 a, float64x2 b);
/// @brief Nonzero if `a` is a signaling NaN. (`f64x2`-tagged spelling of `issignalingdd`)
MULTIFLOATS_API int issignalingf64x2(float64x2 a);
/// @brief Quiet NaN with the payload parsed from `tagp` (libquadmath `nanq`). (`f64x2`-tagged spelling of `nandd`)
MULTIFLOATS_API float64x2 nanf64x2(const char *tagp);
/// @brief `cos(a) + i*sin(a)` (libquadmath `cexpiq`). (`f64x2`-tagged spelling of `cexpidd`)
MULTIFLOATS_API complex64x2 cexpif64x2(float64x2 a);

#ifdef __cplusplus
}  /* extern "C" */
#if defined(__clang__)
#  pragma clang diagnostic pop
#endif

/* ----------------------------------------------------------------------------
 * C++ convenience overloads of the `c*dd` entry points taking
 * `std::complex<float64x2>`. These are header-only — they forward to the
 * POD-`complex64x2` extern "C" symbols above, so the ABI surface stays
 * LTO-compatible while C++ consumers keep an idiomatic std::complex API.
 *
 * Marshaling is one word-for-word copy into/out of `complex64x2` (two
 * float64x2 limbs). The `static_assert`s above pin layout compat; any
 * future divergence trips at compile time.
 * -------------------------------------------------------------------------- */
namespace detail {
inline complex64x2 cx_to_pod(std::complex<float64x2> const &z) {
  return { z.real(), z.imag() };
}
inline std::complex<float64x2> cx_from_pod(complex64x2 const &z) {
  return std::complex<float64x2>(z.re, z.im);
}
} // namespace detail

/// @brief Complex sum (`std::complex<float64x2>` overload).
inline std::complex<float64x2> cadddd(std::complex<float64x2> const &a,
                                      std::complex<float64x2> const &b) {
  return detail::cx_from_pod(cadddd(detail::cx_to_pod(a), detail::cx_to_pod(b)));
}
/// @brief Complex difference (`std::complex<float64x2>` overload).
inline std::complex<float64x2> csubdd(std::complex<float64x2> const &a,
                                      std::complex<float64x2> const &b) {
  return detail::cx_from_pod(csubdd(detail::cx_to_pod(a), detail::cx_to_pod(b)));
}
/// @brief Complex product (`std::complex<float64x2>` overload).
inline std::complex<float64x2> cmuldd(std::complex<float64x2> const &a,
                                      std::complex<float64x2> const &b) {
  return detail::cx_from_pod(cmuldd(detail::cx_to_pod(a), detail::cx_to_pod(b)));
}
/// @brief Complex quotient (`std::complex<float64x2>` overload).
inline std::complex<float64x2> cdivdd(std::complex<float64x2> const &a,
                                      std::complex<float64x2> const &b) {
  return detail::cx_from_pod(cdivdd(detail::cx_to_pod(a), detail::cx_to_pod(b)));
}

#define MULTIFLOATS_CDD_UNARY_OVERLOAD(name)                                 \
  inline std::complex<float64x2> name(std::complex<float64x2> const &z) {    \
    return detail::cx_from_pod(name(detail::cx_to_pod(z)));                  \
  }
/// @brief Complex base-e exponential (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(cexpdd)
/// @brief Complex `exp(z) - 1` (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(cexpm1dd)
/// @brief Complex natural logarithm (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(clogdd)
/// @brief Complex base-2 logarithm (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(clog2dd)
/// @brief Complex base-10 logarithm (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(clog10dd)
/// @brief Complex `log(1 + z)` (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(clog1pdd)
/// @brief Complex square root (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(csqrtdd)
/// @brief Complex sine (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(csindd)
/// @brief Complex `sin(pi*z)` (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(csinpidd)
/// @brief Complex cosine (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(ccosdd)
/// @brief Complex `cos(pi*z)` (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(ccospidd)
/// @brief Complex tangent (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(ctandd)
/// @brief Complex arc sine (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(casindd)
/// @brief Complex arc cosine (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(cacosdd)
/// @brief Complex arc tangent (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(catandd)
/// @brief Complex hyperbolic sine (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(csinhdd)
/// @brief Complex hyperbolic cosine (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(ccoshdd)
/// @brief Complex hyperbolic tangent (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(ctanhdd)
/// @brief Complex inverse hyperbolic sine (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(casinhdd)
/// @brief Complex inverse hyperbolic cosine (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(cacoshdd)
/// @brief Complex inverse hyperbolic tangent (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(catanhdd)
/// @brief Complex projection onto the Riemann sphere (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(cprojdd)
/// @brief Complex conjugate (`std::complex<float64x2>` overload).
MULTIFLOATS_CDD_UNARY_OVERLOAD(conjdd)
#undef MULTIFLOATS_CDD_UNARY_OVERLOAD

/// @brief Complex power `z**w` (`std::complex<float64x2>` overload).
inline std::complex<float64x2> cpowdd(std::complex<float64x2> const &z,
                                      std::complex<float64x2> const &w) {
  return detail::cx_from_pod(cpowdd(detail::cx_to_pod(z), detail::cx_to_pod(w)));
}

/// @brief Complex magnitude `|z|` (`std::complex<float64x2>` overload).
inline float64x2 cabsdd (std::complex<float64x2> const &z) { return cabsdd (detail::cx_to_pod(z)); }
/// @brief Complex argument `arg(z)` (`std::complex<float64x2>` overload).
inline float64x2 cargdd (std::complex<float64x2> const &z) { return cargdd (detail::cx_to_pod(z)); }
/// @brief Real part (`std::complex<float64x2>` overload).
inline float64x2 crealdd(std::complex<float64x2> const &z) { return z.real(); }
/// @brief Imaginary part (`std::complex<float64x2>` overload).
inline float64x2 cimagdd(std::complex<float64x2> const &z) { return z.imag(); }

// `cexpidd` takes a real `float64x2` so its overload set doesn't need a
// std::complex variant — `multifloats::cexpi(a)` already returns
// `std::complex<float64x2>` via the fused sincos path.

#endif  /* __cplusplus */

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
// src/float64x2/abi.inc) are thin marshaling shims that call the
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

/// @brief Absolute value, `|x|`.
inline constexpr float64x2 fabs(float64x2 const &x) {
  std::size_t i = detail::first_nonzero_limb_index(x);
  // All limbs are a (possibly signed) zero: return canonical +0 so
  // fabs((-0, 0)) yields (+0, +0), matching IEEE fabs(-0.0) = +0.0.
  if (i == 2) return float64x2();
  return std::signbit(x.limbs[i]) ? -x : x;
}

/// @brief Minimum of two double-doubles (NaN-as-missing-data, like C `fmin`).
inline constexpr float64x2 fmin(float64x2 const &a, float64x2 const &b) {
  // `b != b` (b is NaN) ⇒ return a. Combined with `a < b` being false for a
  // NaN, this returns the non-NaN operand if exactly one is NaN, with a single
  // extra comparison. `(value, NaN)` would otherwise wrongly return NaN.
  return (a < b || b.limbs[0] != b.limbs[0]) ? a : b;
}

/// @brief Maximum of two double-doubles (NaN-as-missing-data, like C `fmax`).
inline constexpr float64x2 fmax(float64x2 const &a, float64x2 const &b) {
  return (a < b || a.limbs[0] != a.limbs[0]) ? b : a;
}

// C23 min/max variants. fmax/fmin treat NaN as missing data ("if exactly
// one is NaN, return the other"); the C23 fmaximum/fminimum family
// instead propagates NaN, and also pins the signed-zero result so
// fmaximum(+0, -0) = +0 / fminimum(+0, -0) = -0 regardless of order.
// The *_num suffix variants restore the NaN-as-missing semantics on top
// of the signed-zero pinning, which is what most code wants when one
// operand may be a sentinel NaN. Predicates live on limbs[0] (the
// sign-bearing limb when the DD is a zero) so these inline before the
// dedicated isnan/signbit overloads later in the header.

inline constexpr float64x2 _nan_dd() {
  float64x2 r;
  r.limbs[0] = std::numeric_limits<double>::quiet_NaN();
  r.limbs[1] = 0.0;
  return r;
}

inline constexpr float64x2 _negzero_dd() {
  float64x2 r;
  r.limbs[0] = -0.0;
  r.limbs[1] = 0.0;
  return r;
}

/// @brief Maximum, propagating NaN and treating +0 as greater than -0 (C23).
inline constexpr float64x2 fmaximum(float64x2 const &a, float64x2 const &b) {
  if (std::isnan(a.limbs[0]) || std::isnan(b.limbs[0])) return _nan_dd();
  if (a.limbs[0] == 0.0 && b.limbs[0] == 0.0) {
    return (std::signbit(a.limbs[0]) && std::signbit(b.limbs[0]))
               ? _negzero_dd()
               : float64x2();
  }
  return (a < b) ? b : a;
}

/// @brief Minimum, propagating NaN and treating -0 as less than +0 (C23).
inline constexpr float64x2 fminimum(float64x2 const &a, float64x2 const &b) {
  if (std::isnan(a.limbs[0]) || std::isnan(b.limbs[0])) return _nan_dd();
  if (a.limbs[0] == 0.0 && b.limbs[0] == 0.0) {
    return (std::signbit(a.limbs[0]) || std::signbit(b.limbs[0]))
               ? _negzero_dd()
               : float64x2();
  }
  return (a < b) ? a : b;
}

/// @brief Maximum, ignoring NaN and treating +0 as greater than -0 (C23).
inline constexpr float64x2 fmaximum_num(float64x2 const &a,
                                         float64x2 const &b) {
  bool an = std::isnan(a.limbs[0]), bn = std::isnan(b.limbs[0]);
  if (an && bn) return _nan_dd();
  if (an) return b;
  if (bn) return a;
  if (a.limbs[0] == 0.0 && b.limbs[0] == 0.0) {
    return (std::signbit(a.limbs[0]) && std::signbit(b.limbs[0]))
               ? _negzero_dd()
               : float64x2();
  }
  return (a < b) ? b : a;
}

/// @brief Minimum, ignoring NaN and treating -0 as less than +0 (C23).
inline constexpr float64x2 fminimum_num(float64x2 const &a,
                                         float64x2 const &b) {
  bool an = std::isnan(a.limbs[0]), bn = std::isnan(b.limbs[0]);
  if (an && bn) return _nan_dd();
  if (an) return b;
  if (bn) return a;
  if (a.limbs[0] == 0.0 && b.limbs[0] == 0.0) {
    return (std::signbit(a.limbs[0]) || std::signbit(b.limbs[0]))
               ? _negzero_dd()
               : float64x2();
  }
  return (a < b) ? a : b;
}

// C23 by-magnitude variants: pick by |a| vs |b|; on equal magnitude
// (including ±0 / ±x ties) the tie-break matches fmaximum / fminimum.
// Compare full DD magnitudes via fabs — comparing only |a.hi| against
// |b.hi| misses cases where two DDs round to the same hi but the lo
// limbs flip the magnitude ordering (e.g. a = (-x, +eps), b = (+x, -eps)
// → |a| > |b| but |a.hi| == |b.hi|).

/// @brief Value with the larger magnitude (C23 `fmaximum_mag`).
inline constexpr float64x2 fmaximum_mag(float64x2 const &a,
                                         float64x2 const &b) {
  if (std::isnan(a.limbs[0]) || std::isnan(b.limbs[0])) return _nan_dd();
  float64x2 aa = fabs(a), bb = fabs(b);
  if (aa > bb) return a;
  if (aa < bb) return b;
  return fmaximum(a, b);
}

/// @brief Value with the smaller magnitude (C23 `fminimum_mag`).
inline constexpr float64x2 fminimum_mag(float64x2 const &a,
                                         float64x2 const &b) {
  if (std::isnan(a.limbs[0]) || std::isnan(b.limbs[0])) return _nan_dd();
  float64x2 aa = fabs(a), bb = fabs(b);
  if (aa < bb) return a;
  if (aa > bb) return b;
  return fminimum(a, b);
}

/// @brief Larger-magnitude value, ignoring NaN (C23).
inline constexpr float64x2 fmaximum_mag_num(float64x2 const &a,
                                             float64x2 const &b) {
  bool an = std::isnan(a.limbs[0]), bn = std::isnan(b.limbs[0]);
  if (an && bn) return _nan_dd();
  if (an) return b;
  if (bn) return a;
  float64x2 aa = fabs(a), bb = fabs(b);
  if (aa > bb) return a;
  if (aa < bb) return b;
  return fmaximum(a, b);
}

/// @brief Smaller-magnitude value, ignoring NaN (C23).
inline constexpr float64x2 fminimum_mag_num(float64x2 const &a,
                                             float64x2 const &b) {
  bool an = std::isnan(a.limbs[0]), bn = std::isnan(b.limbs[0]);
  if (an && bn) return _nan_dd();
  if (an) return b;
  if (bn) return a;
  float64x2 aa = fabs(a), bb = fabs(b);
  if (aa < bb) return a;
  if (aa > bb) return b;
  return fminimum(a, b);
}


/// @brief Absolute value, `|x|`.
inline constexpr float64x2 abs(float64x2 const &x) { return fabs(x); }
/// @brief Minimum of two double-doubles.
inline constexpr float64x2 min(float64x2 const &x, float64x2 const &y) { return fmin(x, y); }
/// @brief Maximum of two double-doubles.
inline constexpr float64x2 max(float64x2 const &x, float64x2 const &y) { return fmax(x, y); }

/// @brief True if the sign bit of `x` is set.
inline constexpr bool signbit(float64x2 const &x) {
  // For non-canonical zero-hi DDs (e.g. (+0, -eps)), the sign lives in
  // the first nonzero limb. Fall through to signbit(hi) when every limb
  // is a (possibly signed) zero, preserving IEEE -0 semantics.
  std::size_t i = detail::first_nonzero_limb_index(x);
  return std::signbit(x.limbs[i == 2 ? 0 : i]);
}

/// @brief True if `x` is finite.
inline constexpr bool isfinite(float64x2 const &x) {
  // A DD with finite hi and non-finite lo is classified non-finite — this
  // matters for the operator+ short-circuit.
  for (std::size_t i = 0; i < 2; ++i) {
    if (!std::isfinite(x.limbs[i])) return false;
  }
  return true;
}

/// @brief True if `x` is infinite.
inline constexpr bool isinf(float64x2 const &x) {
  // Scan both limbs for symmetry with isnan / isfinite: a DD with a finite
  // hi and an inf lo (non-canonical but constructible via the C ABI) should
  // classify as inf, not "neither finite nor inf".
  for (std::size_t i = 0; i < 2; ++i) {
    if (std::isinf(x.limbs[i])) return true;
  }
  return false;
}

/// @brief True if `x` is NaN.
inline constexpr bool isnan(float64x2 const &x) {
  for (std::size_t i = 0; i < 2; ++i) {
    if (std::isnan(x.limbs[i])) return true;
  }
  return false;
}

/// @brief C99 floating-point classification of `x` (`FP_NAN`, `FP_INFINITE`, ...).
inline constexpr int fpclassify(float64x2 const &x) {
  return std::fpclassify(x.limbs[0]);
}

/// @brief Scale by a power of two, `x * 2**n`.
inline constexpr float64x2 ldexp(float64x2 const &x, int n) {
  // Build the power-of-two scale once; multiplication by an exact power of
  // two is exact for every limb (no rounding, no renorm), avoiding the
  // two library calls of std::ldexp.
  double scale = std::ldexp(1.0, n);
  return {x.limbs[0] * scale, x.limbs[1] * scale};
}

/// @brief Scale by a power of two, `x * 2**n`.
inline constexpr float64x2 scalbn(float64x2 const &x, int n) {
  // POSIX alias of ldexp for FLT_RADIX == 2 (which is guaranteed by IEEE 754).
  return ldexp(x, n);
}

/// @brief Unbiased radix-2 exponent of `x`, as an `int`.
inline constexpr int ilogb(float64x2 const &x) {
  return std::ilogb(x.limbs[0]);
}

// C23 llogb: same as ilogb but returns long.
/// @brief Unbiased radix-2 exponent of `x`, as a `long`.
inline constexpr long llogb(float64x2 const &x) {
  return static_cast<long>(std::ilogb(x.limbs[0]));
}

/// @brief Magnitude of `x` with the sign of `y`.
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

/// @brief Round toward -infinity.
inline constexpr float64x2 floor(float64x2 const &x) {
  if (x.limbs[0] == 0.0 && x.limbs[1] == 0.0) return x;  // preserve ±0
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

/// @brief Round toward +infinity.
inline constexpr float64x2 ceil(float64x2 const &x) {
  if (x.limbs[0] == 0.0 && x.limbs[1] == 0.0) return x;  // preserve ±0
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

/// @brief Round toward zero.
inline constexpr float64x2 trunc(float64x2 const &x) {
  return std::signbit(x.limbs[0]) ? -floor(-x) : floor(x);
}

/// @brief Round to nearest, halfway away from zero.
inline constexpr float64x2 round(float64x2 const &x) {
  if (x.limbs[0] == 0.0 && x.limbs[1] == 0.0) return x;  // preserve ±0
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

/// @brief Round to integer in the current rounding mode, without raising inexact.
inline constexpr float64x2 nearbyint(float64x2 const &x) {
  if (x.limbs[0] == 0.0 && x.limbs[1] == 0.0) return x;  // preserve ±0 (rint/roundeven inherit)
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

/// @brief Round to integer in the current rounding mode.
inline constexpr float64x2 rint(float64x2 const &x) { return nearbyint(x); }

// C23 fromfp / ufromfp / fromfpx / ufromfpx: round x to integer per
// the requested rounding mode and return as intmax_t / uintmax_t IF
// it fits in `width` bits; otherwise raise FE_INVALID and return
// IM_MIN / IM_MAX (the *fp variants raise FE_INEXACT in addition,
// when the rounding actually moved the value).
//
// `rnd` is one of the C23 FP_INT_* rounding-direction macros (NOT the FE_*
// rounding-environment macros — they have different values):
//   FP_INT_UPWARD            → ceil  (toward +∞)
//   FP_INT_DOWNWARD          → floor (toward −∞)
//   FP_INT_TOWARDZERO        → trunc
//   FP_INT_TONEARESTFROMZERO → round (ties away from zero)
//   FP_INT_TONEAREST         → round to nearest, ties to even (default)
//
// For DD: each rounding direction has a native kernel; we dispatch on `rnd`,
// run the kernel, then narrow to integer.

// Fallback for a pre-C23 <math.h> that doesn't expose FP_INT_* (their values
// are fixed by the standard).
#ifndef FP_INT_UPWARD
#define FP_INT_UPWARD 0
#define FP_INT_DOWNWARD 1
#define FP_INT_TOWARDZERO 2
#define FP_INT_TONEARESTFROMZERO 3
#define FP_INT_TONEAREST 4
#endif

namespace detail {

inline float64x2 round_for_fromfp(float64x2 const &x, int rnd) {
  switch (rnd) {
  case FP_INT_UPWARD:            return ceil(x);
  case FP_INT_DOWNWARD:          return floor(x);
  case FP_INT_TOWARDZERO:        return trunc(x);
  case FP_INT_TONEARESTFROMZERO: return round(x);  // ties away from zero
  case FP_INT_TONEAREST:
  default: {  // round to nearest, ties to even
    int saved = std::fegetround();
    if (saved != FE_TONEAREST) std::fesetround(FE_TONEAREST);
    float64x2 r = nearbyint(x);
    if (saved != FE_TONEAREST) std::fesetround(saved);
    return r;
  }
  }
}

}  // namespace detail

/// @brief Round `x` to a signed integer of `width` bits using mode `rnd` (C23 `fromfp`).
inline intmax_t fromfp(float64x2 const &x, int rnd, unsigned int width) {
  if (std::isnan(x.limbs[0]) || std::isinf(x.limbs[0])) {
    std::feraiseexcept(FE_INVALID);
    return INTMAX_MIN;
  }
  if (width == 0 || width > 63) width = 63;
  float64x2 r = detail::round_for_fromfp(x, rnd);
  // Range check against [-2^(width-1), 2^(width-1) - 1]. Compare in DD.
  double upper_d = std::ldexp(1.0, static_cast<int>(width) - 1);  // 2^(w-1)
  float64x2 upper = float64x2(upper_d) - float64x2(1.0);
  float64x2 lower = -float64x2(upper_d);
  if (r > upper || r < lower) {
    std::feraiseexcept(FE_INVALID);
    return INTMAX_MIN;
  }
  return static_cast<intmax_t>(r.limbs[0]) +
         static_cast<intmax_t>(r.limbs[1]);
}

/// @brief Round `x` to an unsigned integer of `width` bits using mode `rnd` (C23 `ufromfp`).
inline uintmax_t ufromfp(float64x2 const &x, int rnd, unsigned int width) {
  if (std::isnan(x.limbs[0]) || std::isinf(x.limbs[0])) {
    std::feraiseexcept(FE_INVALID);
    return 0;
  }
  if (width == 0 || width > 64) width = 64;
  float64x2 r = detail::round_for_fromfp(x, rnd);
  if (r < float64x2()) {  // negative — out of unsigned range
    std::feraiseexcept(FE_INVALID);
    return 0;
  }
  double upper_d = std::ldexp(1.0, static_cast<int>(width));  // 2^w
  float64x2 upper = float64x2(upper_d) - float64x2(1.0);
  if (r > upper) {
    std::feraiseexcept(FE_INVALID);
    return 0;
  }
  // r.limbs[1] may be negative (canonical DD allows a low limb of either sign);
  // casting a negative double straight to uintmax_t is UB, so route the low
  // limb through intmax_t first. r.limbs[0] is a non-negative integer < 2^64.
  return static_cast<uintmax_t>(r.limbs[0]) +
         static_cast<uintmax_t>(static_cast<intmax_t>(r.limbs[1]));
}

/// @brief Like `fromfp`, also raising the inexact exception (C23 `fromfpx`).
inline intmax_t fromfpx(float64x2 const &x, int rnd, unsigned int width) {
  intmax_t out = fromfp(x, rnd, width);
  // Raise FE_INEXACT if the rounding moved the value (i.e. x wasn't already
  // an integer). Floating-point comparison against the rounded result.
  float64x2 r = detail::round_for_fromfp(x, rnd);
  if (!(r == x)) std::feraiseexcept(FE_INEXACT);
  return out;
}

/// @brief Like `ufromfp`, also raising the inexact exception (C23 `ufromfpx`).
inline uintmax_t ufromfpx(float64x2 const &x, int rnd, unsigned int width) {
  uintmax_t out = ufromfp(x, rnd, width);
  float64x2 r = detail::round_for_fromfp(x, rnd);
  if (!(r == x)) std::feraiseexcept(FE_INEXACT);
  return out;
}

// C23 roundeven: round to nearest, ties to even, INDEPENDENT of the
// current FE rounding mode (unlike nearbyint, which obeys FE_DOWNWARD
// etc.). For DD: nearbyint already implements round-to-even at the
// hardware level when FE_TONEAREST is set, so we use it as the kernel
// and explicitly clamp the rounding direction with a fenv guard so a
// caller running under FE_DOWNWARD still gets ties-to-even.
/// @brief Round to nearest integer, ties to even (C23 `roundeven`).
inline float64x2 roundeven(float64x2 const &x) {
#if defined(__cplusplus) && __cplusplus >= 201703L
  // Save / restore FE rounding mode around the kernel call. constexpr
  // can't reach fegetround/fesetround so this overload is non-constexpr.
  int saved = std::fegetround();
  if (saved != FE_TONEAREST) std::fesetround(FE_TONEAREST);
  float64x2 r = nearbyint(x);
  if (saved != FE_TONEAREST) std::fesetround(saved);
  return r;
#else
  return nearbyint(x);
#endif
}

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

/// @brief Round to nearest, halfway away from zero, returning `long`.
inline constexpr long lround(float64x2 const &x) {
  return detail::lround_adjust<long>(x, std::lround(x.limbs[0]));
}

/// @brief Round to nearest, halfway away from zero, returning `long long`.
inline constexpr long long llround(float64x2 const &x) {
  return detail::lround_adjust<long long>(x, std::llround(x.limbs[0]));
}

/// @brief Round to integer in the current rounding mode, returning `long`.
inline constexpr long lrint(float64x2 const &x) {
  return std::lrint(rint(x).limbs[0]);
}

/// @brief Round to integer in the current rounding mode, returning `long long`.
inline constexpr long long llrint(float64x2 const &x) {
  return std::llrint(rint(x).limbs[0]);
}

// =============================================================================
// Floating-point manipulation
// =============================================================================

/// @brief Split `x` into a fraction in [0.5, 1) and an exponent stored in `*exp`.
inline constexpr float64x2 frexp(float64x2 const &x, int *exp) {
  float64x2 r;
  int e = 0;
  r.limbs[0] = std::frexp(x.limbs[0], &e);
  r.limbs[1] = std::ldexp(x.limbs[1], -e);
  *exp = e;
  return r;
}

/// @brief Split `x` into integer part (stored in `*iptr`) and fractional part.
inline constexpr float64x2 modf(float64x2 const &x, float64x2 *iptr) {
  *iptr = trunc(x);
  // x - trunc(x) would flip a ±0 input to +0; the fractional part of ±0 is ±0.
  if (x.limbs[0] == 0.0 && x.limbs[1] == 0.0) return x;
  return x - *iptr;
}

/// @brief Scale by a power of two, `x * 2**n` (`long` exponent).
inline constexpr float64x2 scalbln(float64x2 const &x, long n) {
  return {std::scalbln(x.limbs[0], n), std::scalbln(x.limbs[1], n)};
}

/// @brief Unbiased radix-2 exponent of `x`, as a floating-point value.
inline constexpr float64x2 logb(float64x2 const &x) {
  float64x2 r;
  r.limbs[0] = std::logb(x.limbs[0]);
  return r;
}

/// @brief Next representable double-double after `x` toward `y`.
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

/// @brief Next representable double-double after `x` toward `y`.
inline constexpr float64x2 nexttoward(float64x2 const &x, float64x2 const &y) {
  return nextafter(x, y);
}

// C23 IEEE 754 metadata operations. For multifloats float64x2 these are
// thin shims over the lead-limb double semantics — the DD type has no
// payload of its own beyond what the leading limb's IEEE encoding
// carries, and the lo limb carries no NaN payload at all (any non-NaN
// canonical DD has a finite normalized lo).

// Always returns the input — DD has no non-canonical encoding for
// finite values (the lo-limb invariant |lo| <= 0.5 ulp(hi) is always
// maintained internally). NaN canonicalization is delegated to the
// underlying double: a sNaN hi becomes a qNaN.
/// @brief Canonicalize `x` to renormalized limbs (C23 `canonicalize`).
inline constexpr float64x2 canonicalize(float64x2 const &x) {
  if (std::isnan(x.limbs[0])) {
    float64x2 r;
    r.limbs[0] = x.limbs[0] + 0.0;  // sNaN → qNaN per IEEE 754
    r.limbs[1] = 0.0;
    return r;
  }
  return x;
}

// Signaling-equal predicate. multifloats does not distinguish sNaN /
// qNaN at the public API level (every NaN we emit is qNaN), so this
// is functionally `a == b`. Provided for C23 conformance.
/// @brief Signaling equality `x == y`, raising invalid for NaN (C23 `iseqsig`).
inline constexpr bool iseqsig(float64x2 const &a, float64x2 const &b) {
  return a == b;
}

// IEEE 754 totalOrder predicate generalized to DD: a lexicographic
// comparison over (hi, lo) treating ±0 / NaN sign-and-payload like the
// underlying double does. Captures the same monotonicity guarantees as
// std::totalOrder<double> on the leading limb, with the lo limb as a
// secondary key when hi limbs tie.
/// @brief IEEE 754 total-ordering predicate: true if `x` precedes or equals `y`.
inline bool totalorder(float64x2 const &a, float64x2 const &b) {
  // IEEE-754 totalOrder key transform:
  //   negative bit pattern → ~bits  (more-negative value → smaller key)
  //   positive bit pattern → bits ^ 0x80…0  (sets top bit so the key
  //     beats every negative key under unsigned compare)
  // Compare as unsigned so the sign-bit semantics line up.
  auto key = [](double d) -> unsigned long long {
    unsigned long long b;
    std::memcpy(&b, &d, sizeof b);
    long long sb;
    std::memcpy(&sb, &b, sizeof sb);
    return (sb < 0) ? ~b : b ^ 0x8000000000000000ULL;
  };
  unsigned long long ka = key(a.limbs[0]), kb = key(b.limbs[0]);
  if (ka != kb) return ka < kb;
  ka = key(a.limbs[1]); kb = key(b.limbs[1]);
  return ka <= kb;
}

/// @brief IEEE 754 total ordering on magnitudes, `totalorder(|x|, |y|)`.
inline bool totalordermag(float64x2 const &a, float64x2 const &b) {
  return totalorder(fabs(a), fabs(b));
}

// NaN payload extraction / injection. For DD we route through the
// leading limb's IEEE-754 NaN payload. setpayload returns 0 on
// success (matching libm convention), 1 on failure (payload too
// large to fit in the available NaN bits, or x not finite).
/// @brief Extract the payload of the NaN `x` (C23 `getpayload`).
inline float64x2 getpayload(float64x2 const &x) {
  if (!std::isnan(x.limbs[0])) return float64x2();
  unsigned long long bits;
  double hi = x.limbs[0];
  std::memcpy(&bits, &hi, sizeof bits);
  // Mask off sign + exponent + quiet bit; remaining 51 bits are the
  // payload. Convert through unsigned long long → double via the
  // float64x2 ctor (exact when bits fit in 53-bit mantissa).
  bits &= 0x0007ffffffffffffULL;
  return float64x2(static_cast<double>(bits));
}

/// @brief Make a quiet NaN with the given payload (C23 `setpayload`); nonzero on failure.
inline int setpayload(float64x2 &x, float64x2 const &payload) {
  // Payload must be a non-negative integer < 2^51 to fit in the qNaN
  // mantissa bits. Anything else: leave x unchanged, return 1.
  double ph = payload.limbs[0];
  if (!(ph >= 0.0) || ph != std::trunc(ph) || ph >= 0x1.0p51) {
    return 1;
  }
  unsigned long long pl = static_cast<unsigned long long>(ph);
  // qNaN encoding: exponent all-ones + quiet bit set + payload in low 51 bits.
  unsigned long long bits = 0x7ff8000000000000ULL | (pl & 0x0007ffffffffffffULL);
  double hi;
  std::memcpy(&hi, &bits, sizeof hi);
  x.limbs[0] = hi;
  x.limbs[1] = 0.0;
  return 0;
}

/// @brief Make a signaling NaN with the given payload (C23 `setpayloadsig`); nonzero on failure.
inline int setpayloadsig(float64x2 &x, float64x2 const &payload) {
  // Signaling NaN encoding: exponent all-ones + quiet bit clear +
  // payload in low 51 bits with at least one nonzero bit (else it's
  // an infinity). Caller must supply a nonzero payload.
  double ph = payload.limbs[0];
  if (!(ph > 0.0) || ph != std::trunc(ph) || ph >= 0x1.0p51) {
    return 1;
  }
  unsigned long long pl = static_cast<unsigned long long>(ph);
  unsigned long long bits = 0x7ff0000000000000ULL | (pl & 0x0007ffffffffffffULL);
  double hi;
  std::memcpy(&hi, &bits, sizeof hi);
  x.limbs[0] = hi;
  x.limbs[1] = 0.0;
  return 0;
}

// C23 nextup / nextdown: IEEE 754 next-toward-+inf / -inf, regardless
// of x. Forward to nextafter with a sentinel ±inf.
/// @brief Least representable double-double greater than `x` (C23 `nextup`).
inline constexpr float64x2 nextup(float64x2 const &x) {
  return nextafter(x, float64x2(std::numeric_limits<double>::infinity()));
}
/// @brief Greatest representable double-double less than `x` (C23 `nextdown`).
inline constexpr float64x2 nextdown(float64x2 const &x) {
  return nextafter(x, float64x2(-std::numeric_limits<double>::infinity()));
}

// libquadmath nanq parity. The tag string encodes the NaN payload via
// std::nan; lo is set to 0 so the result is a canonical DD NaN.
/// @brief Quiet NaN with the payload parsed from `tagp`.
inline float64x2 nan(const char *tagp) {
  return float64x2(std::nan(tagp ? tagp : ""), 0.0);
}

// =============================================================================
// Basic arithmetic helpers
// =============================================================================

/// @brief Fused multiply-add `x*y + z` with a single rounding.
inline constexpr float64x2 fma(float64x2 const &x, float64x2 const &y,
                               float64x2 const &z) {
  // Not a hardware fma, but provides the cmath interface.
  return x * y + z;
}

/// @brief Floating-point remainder of `x / y` (sign of `x`; C `fmod`).
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
/// @brief Floored remainder of `x / y` (sign of `y`; Fortran `modulo`).
inline constexpr float64x2 modulo(float64x2 const &x, float64x2 const &y) {
  float64x2 r = fmod(x, y);
  if (r.limbs[0] == 0.0 && r.limbs[1] == 0.0) return r;
  return (signbit(r) != signbit(y)) ? r + y : r;
}

// Integer-exponent power via exponentiation-by-squaring. Returns exact
// 1 for n == 0 (including for 0**0, matching C pow and Fortran). For
// n < 0 returns 1 / powi(base, -n), adding one DD division at the end.
// C23 pown: integer-exponent power. Identical kernel to powi (defined
// just below) — kept as a separate name so callers using the standard
// C23 spelling do not need to know the multifloats-internal alias.
/// @brief Integer power `x**n` (C23 `pown`).
inline constexpr float64x2 pown(float64x2 base, int n);

/// @brief Integer power `x**n` by exponentiation-by-squaring.
inline constexpr float64x2 powi(float64x2 base, int n) {
  if (n == 0) return float64x2(1.0);
  // ±0 base (IEEE pow): magnitude is 0 for n>0, +∞ for n<0; the sign is carried
  // only by an ODD exponent. Handle here because the compensated multiply below
  // flushes signed zero (so −0 would come out +0) and 1/(0,0) leaves a
  // non-canonical lo=∞.
  if (base.limbs[0] == 0.0 && base.limbs[1] == 0.0) {
    double s = (n & 1) ? std::copysign(1.0, base.limbs[0]) : 1.0;
    double mag = (n > 0) ? 0.0 : std::numeric_limits<double>::infinity();
    return float64x2(s * mag, 0.0);
  }
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

inline constexpr float64x2 pown(float64x2 base, int n) { return powi(base, n); }

/// @brief IEEE 754 remainder of `x / y` (round-to-nearest quotient).
// Not constexpr: IEEE 754 `remainder` rounds the quotient to nearest, ties to
// EVEN — `round` (ties away from zero) is wrong at half-integer quotients, e.g.
// remainder(1,2) must be +1, not -1. `roundeven` is the correct rounding but is
// non-constexpr (it pins the FE rounding mode), so these lose constexpr — as
// does the C library `remainder`/`remquo`, which were never constexpr anyway.
inline float64x2 remainder(float64x2 const &x, float64x2 const &y) {
  return x - roundeven(x / y) * y;
}

/// @brief IEEE 754 remainder of `x / y`; low quotient bits stored in `*quo`.
inline float64x2 remquo(float64x2 const &x, float64x2 const &y, int *quo) {
  float64x2 q = roundeven(x / y);
  // q can exceed the range of int for huge x/y; static_cast<int> of an
  // out-of-range double is UB, so clamp. (C only mandates the low-order bits;
  // any in-range value with the right sign satisfies the contract.)
  double qd = q.limbs[0];
  *quo = (qd >= 2147483648.0)    ? 2147483647
         : (qd <= -2147483648.0) ? (-2147483647 - 1)
                                 : static_cast<int>(qd);
  return x - q * y;
}

/// @brief Positive difference, `max(x - y, 0)`.
inline constexpr float64x2 fdim(float64x2 const &x, float64x2 const &y) {
  // C requires fdim to return NaN when either argument is NaN; the bare
  // `x > y` comparison is false for NaN and would wrongly yield 0.
  if (x.limbs[0] != x.limbs[0] || y.limbs[0] != y.limbs[0])
    return float64x2(std::numeric_limits<double>::quiet_NaN());
  return (x > y) ? (x - y) : float64x2();
}

// C++20 std::lerp: exact at the endpoints, monotonic in t, and does not
// overshoot when a and b have the same sign and t is in [0, 1].
/// @brief Linear interpolation `a + t*(b - a)` (`std::lerp`).
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
/// @brief Square root.
inline constexpr float64x2 sqrt(float64x2 const &x);

namespace detail {

// Polynomial and conversion constants needed by the DD kernels are provided
// through dd_constants.hh which is included only in the .cc file. The
// erf/erfc rational-approximation constants also live in the .cc file.
//
// horner / neval / deval (DD polynomial evaluators) used to live here as
// inline definitions but have been moved to src/float64x2/poly.inc,
// which is pulled into namespace multifloats::detail inside
// multifloats_math.cc. Every call site lives in that TU, so inlining
// decisions are unchanged; the public header no longer carries ~250 lines
// of Estrin switch bodies.

} // namespace detail

// =============================================================================
// Roots
// =============================================================================

inline constexpr float64x2 sqrt(float64x2 const &x) {
  double c = std::sqrt(x.limbs[0]);
  // Bail on 0, -0, negative, NaN, +Inf — the refinement step would compute
  // `inf - inf = NaN` for +Inf and `0 / 0` for 0. Leading-limb sqrt handles
  // every IEEE special case.
  if (!(x.limbs[0] > 0.0) || !std::isfinite(c)) {
    float64x2 r;
    r.limbs[0] = c;
    return r;
  }
  // Karp/Markstein: r = c + (x - c*c) / (2c), but evaluated entirely in
  // scalar after computing c*c exactly. two_prod(c, c) gives the full c*c
  // as (u, uu) in 2 FMAs; the residual then folds against both x limbs in
  // scalar and gets divided by 2c. Final fast_two_sum renormalizes.
  // Cheaper than building a DD `s_dd*s_dd` (avoids the DD multiply path
  // entirely) and matches Boost.Multiprecision's eval_sqrt formulation.
  double u = 0.0, uu = 0.0;
  detail::two_prod(c, c, u, uu);
  double cc = ((x.limbs[0] - u) - uu + x.limbs[1]) / (2.0 * c);
  float64x2 r;
  r.limbs[0] = c + cc;
  r.limbs[1] = (c - r.limbs[0]) + cc;
  return r;
}

// C23 rsqrt: 1/sqrt(x). Could be expressed as `1/sqrt(x)` but a direct
// formulation saves one DD-divide and dodges the spurious overflow at
// |x| just below DBL_MAX.
/// @brief Reciprocal square root, `1/sqrt(x)`.
inline constexpr float64x2 rsqrt(float64x2 const &x) {
  double c = 1.0 / std::sqrt(x.limbs[0]);
  // Bail on 0, -0, negative, NaN, +Inf for the same reasons sqrt does:
  // refinement would compute 0/0 or inf−inf otherwise.
  if (!(x.limbs[0] > 0.0) || !std::isfinite(c)) {
    float64x2 r;
    r.limbs[0] = c;
    return r;
  }
  // Newton step for 1/sqrt: y_new = y · (3 − x·y²) / 2 = y + y·(1 − x·y²)/2.
  // The residual (1 − x·c²) cancels to ~ulp_dd, so the dominant product
  // x.hi·c² must be captured exactly via two_prod — a plain mul would
  // round x.hi·c² to a double, leaving (1 − round(x.hi·c²)) at most one
  // ULP wide and capping the answer at double precision (1e-16).
  double u = 0.0, uu = 0.0;
  detail::two_prod(c, c, u, uu);                  // u + uu = c²  exact
  double p = 0.0, pe = 0.0;
  detail::two_prod(x.limbs[0], u, p, pe);         // p + pe = x.hi·c²  exact
  // 1 − p is exact whenever 1/2 ≤ p ≤ 2 (Sterbenz). Convergence guarantees
  // this for the rsqrt fixed point — c is correct to ~ulp(double) so
  // |x·c² − 1| ≲ 1e-15 ⇒ p ∈ [1−ulp, 1+ulp] and the subtraction is exact.
  double s = 1.0 - p;
  // Sub-ULP corrections: pe (~ulp_dd·|p|), x.hi·uu (~ulp_dd·|p|),
  // x.lo·u (~ulp_dd·|p|). Order doesn't matter at this scale; collapse
  // into a scalar.
  double resid = s - pe - x.limbs[0] * uu - x.limbs[1] * u;
  double dc = c * (0.5 * resid);
  float64x2 r;
  r.limbs[0] = c + dc;
  r.limbs[1] = (c - r.limbs[0]) + dc;
  return r;
}

/// @brief Cube root.
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

/// @brief Euclidean distance `sqrt(x^2 + y^2)`, with overflow-safe scaling.
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
/// @brief Base-e exponential, `e**x`.
float64x2 exp   (float64x2 const &x);
/// @brief Base-2 exponential, `2**x`.
float64x2 exp2  (float64x2 const &x);
/// @brief Base-10 exponential, `10**x`.
float64x2 exp10 (float64x2 const &x);
/// @brief `e**x - 1`, accurate for small `x`.
float64x2 expm1 (float64x2 const &x);
/// @brief `2**x - 1`, accurate for small `x`.
float64x2 exp2m1 (float64x2 const &x);
/// @brief `10**x - 1`, accurate for small `x`.
float64x2 exp10m1(float64x2 const &x);
/// @brief Natural (base-e) logarithm.
float64x2 log   (float64x2 const &x);
/// @brief Base-10 logarithm.
float64x2 log10 (float64x2 const &x);
/// @brief Base-2 logarithm.
float64x2 log2  (float64x2 const &x);
/// @brief `log(1 + x)`, accurate for small `x`.
float64x2 log1p (float64x2 const &x);
/// @brief `log2(1 + x)`, accurate for small `x`.
float64x2 log2p1 (float64x2 const &x);
/// @brief `log10(1 + x)`, accurate for small `x`.
float64x2 log10p1(float64x2 const &x);
/// @brief Power `x**y`, via `exp(y * log(x))`.
float64x2 pow   (float64x2 const &x, float64x2 const &y);
/// @brief Power `x**y` for `x >= 0`, via `exp(y * log(x))` (C23 `powr`).
float64x2 powr  (float64x2 const &x, float64x2 const &y);
/// @brief Real `n`-th root, `x**(1/n)` (C23 `rootn`).
float64x2 rootn (float64x2 const &x, int n);
/// @brief `(1 + x)**n` for integer `n`, accurate for small `x` (C23 `compoundn`).
float64x2 compoundn(float64x2 const &x, int n);

// Trigonometric
/// @brief Sine.
float64x2 sin   (float64x2 const &x);
/// @brief Cosine.
float64x2 cos   (float64x2 const &x);
/// @brief Tangent.
float64x2 tan   (float64x2 const &x);
/// @brief Arc sine.
float64x2 asin  (float64x2 const &x);
/// @brief Arc cosine.
float64x2 acos  (float64x2 const &x);
/// @brief Arc tangent.
float64x2 atan  (float64x2 const &x);
/// @brief Arc tangent of `y/x`, using the signs of both to pick the quadrant.
float64x2 atan2 (float64x2 const &y, float64x2 const &x);

// π-scaled trig: {sin,cos,tan}pi(x) = {sin,cos,tan}(π·x),
//                {asin,acos,atan}pi(x) = {asin,acos,atan}(x)/π,
//                atan2pi(y, x) = atan2(y, x)/π.
/// @brief `sin(pi * x)`.
float64x2 sinpi   (float64x2 const &x);
/// @brief `cos(pi * x)`.
float64x2 cospi   (float64x2 const &x);
/// @brief `tan(pi * x)`.
float64x2 tanpi   (float64x2 const &x);
/// @brief `asin(x) / pi`.
float64x2 asinpi  (float64x2 const &x);
/// @brief `acos(x) / pi`.
float64x2 acospi  (float64x2 const &x);
/// @brief `atan(x) / pi`.
float64x2 atanpi  (float64x2 const &x);
/// @brief `atan2(y, x) / pi`.
float64x2 atan2pi (float64x2 const &y, float64x2 const &x);

// Fused sincos / sinhcosh. One range-reduction + Taylor pair produces both
// outputs, roughly halving the transcendental cost when both are needed.
void sincos   (float64x2 const &x, float64x2 &s, float64x2 &c);
void sinhcosh (float64x2 const &x, float64x2 &s, float64x2 &c);

// Hyperbolic
/// @brief Hyperbolic sine.
float64x2 sinh  (float64x2 const &x);
/// @brief Hyperbolic cosine.
float64x2 cosh  (float64x2 const &x);
/// @brief Hyperbolic tangent.
float64x2 tanh  (float64x2 const &x);
/// @brief Inverse hyperbolic sine.
float64x2 asinh (float64x2 const &x);
/// @brief Inverse hyperbolic cosine.
float64x2 acosh (float64x2 const &x);
/// @brief Inverse hyperbolic tangent.
float64x2 atanh (float64x2 const &x);

// Error and gamma. erfcx is the scaled complementary error function,
// erfcx(x) = exp(x^2) * erfc(x).
/// @brief Error function.
float64x2 erf    (float64x2 const &x);
/// @brief Complementary error function, `1 - erf(x)`.
float64x2 erfc   (float64x2 const &x);
/// @brief Scaled complementary error function `exp(x^2) * erfc(x)`.
float64x2 erfcx  (float64x2 const &x);
/// @brief Gamma function, `Gamma(x)`.
float64x2 tgamma (float64x2 const &x);
/// @brief Natural log of `|Gamma(x)|`.
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
/// @brief Bessel function of the first kind, order 0.
float64x2 j0 (float64x2 const &x);
/// @brief Bessel function of the first kind, order 1.
float64x2 j1 (float64x2 const &x);
/// @brief Bessel function of the second kind, order 0.
float64x2 y0 (float64x2 const &x);
/// @brief Bessel function of the second kind, order 1.
float64x2 y1 (float64x2 const &x);
/// @brief Bessel function of the first kind, integer order `n`.
float64x2 jn (int n, float64x2 const &x);
/// @brief Bessel function of the second kind, integer order `n`.
float64x2 yn (int n, float64x2 const &x);
void      jn (int n1, int n2, float64x2 const &x, float64x2 *out);
void      yn (int n1, int n2, float64x2 const &x, float64x2 *out);

// =============================================================================
// Additional classification and ordered comparison
// =============================================================================

/// @brief True if `x` is normal (finite, nonzero, not subnormal).
inline constexpr bool isnormal(float64x2 const &x) {
  return std::isnormal(x.limbs[0]);
}

/// @brief True if `x > y` (quiet for NaN operands).
inline constexpr bool isgreater(float64x2 const &x, float64x2 const &y) {
  return !isnan(x) && !isnan(y) && (x > y);
}

/// @brief True if `x >= y` (quiet for NaN operands).
inline constexpr bool isgreaterequal(float64x2 const &x, float64x2 const &y) {
  return !isnan(x) && !isnan(y) && (x >= y);
}

/// @brief True if `x < y` (quiet for NaN operands).
inline constexpr bool isless(float64x2 const &x, float64x2 const &y) {
  return !isnan(x) && !isnan(y) && (x < y);
}

/// @brief True if `x <= y` (quiet for NaN operands).
inline constexpr bool islessequal(float64x2 const &x, float64x2 const &y) {
  return !isnan(x) && !isnan(y) && (x <= y);
}

/// @brief True if `x < y` or `x > y` (quiet for NaN operands).
inline constexpr bool islessgreater(float64x2 const &x, float64x2 const &y) {
  return !isnan(x) && !isnan(y) && (x != y);
}

/// @brief True if `x` or `y` is NaN.
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
/// @brief `std::to_chars`-style formatter for `float64x2` (scientific format
///        only). Writes into `[first, last)` without a NUL and returns
///        `{ptr, errc{}}` on success, or `{last, value_too_large}` if the
///        buffer is too small. `precision` is the number of significant digits
///        (defaults to 32; clamped to [1, 34]).
std::to_chars_result
to_chars(char *first, char *last, float64x2 const &value,
         std::chars_format fmt, int precision) noexcept;
/// @brief `to_chars` with the default precision of 32 digits.
inline std::to_chars_result
to_chars(char *first, char *last, float64x2 const &value,
         std::chars_format fmt) noexcept {
  return to_chars(first, last, value, fmt, 32);
}
/// @brief `to_chars` with scientific format and the default 32-digit precision.
inline std::to_chars_result
to_chars(char *first, char *last, float64x2 const &value) noexcept {
  return to_chars(first, last, value, std::chars_format::scientific, 32);
}

/// @brief Format `x` as a scientific-notation `std::string` (`precision`
///        significant digits, default 32).
inline std::string to_string(float64x2 const &x, int precision = 32) {
  char buf[MULTIFLOATS_DD_CHARS_BUFSIZE];
  auto res = to_chars(buf, buf + sizeof(buf), x,
                      std::chars_format::scientific, precision);
  // The caller buffer is sized for any legal precision, so value_too_large
  // cannot happen here — `res.ec` is always `errc{}` in practice.
  return std::string(buf, static_cast<std::size_t>(res.ptr - buf));
}
/// @brief Stream `x` in scientific notation (honors the stream's precision).
std::ostream &operator<<(std::ostream &os, float64x2 const &x);

} // namespace multifloats

// ---- std::complex<multifloats::float64x2> specializations ------------------
//
// Explicit specializations for the C++17 <complex> free-function templates
// on `std::complex<multifloats::float64x2>`. Definitions live in
// multifloats_math.cc (complex64x2/std.inc); the C-ABI
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

/// @brief Complex base-e exponential (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> exp  (complex<multifloats::float64x2> const &z);
/// @brief Complex natural logarithm (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> log  (complex<multifloats::float64x2> const &z);
/// @brief Complex square root (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> sqrt (complex<multifloats::float64x2> const &z);
/// @brief Complex power `z**w` (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> pow  (complex<multifloats::float64x2> const &z,
                                                  complex<multifloats::float64x2> const &w);
/// @brief Complex sine (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> sin  (complex<multifloats::float64x2> const &z);
/// @brief Complex cosine (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> cos  (complex<multifloats::float64x2> const &z);
/// @brief Complex tangent (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> tan  (complex<multifloats::float64x2> const &z);
/// @brief Complex arc sine (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> asin (complex<multifloats::float64x2> const &z);
/// @brief Complex arc cosine (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> acos (complex<multifloats::float64x2> const &z);
/// @brief Complex arc tangent (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> atan (complex<multifloats::float64x2> const &z);
/// @brief Complex hyperbolic sine (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> sinh (complex<multifloats::float64x2> const &z);
/// @brief Complex hyperbolic cosine (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> cosh (complex<multifloats::float64x2> const &z);
/// @brief Complex hyperbolic tangent (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> tanh (complex<multifloats::float64x2> const &z);
/// @brief Complex inverse hyperbolic sine (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> asinh(complex<multifloats::float64x2> const &z);
/// @brief Complex inverse hyperbolic cosine (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> acosh(complex<multifloats::float64x2> const &z);
/// @brief Complex inverse hyperbolic tangent (`std::complex<float64x2>` specialization).
template <> complex<multifloats::float64x2> atanh(complex<multifloats::float64x2> const &z);

/// @brief Complex magnitude `|z|` (`std::complex<float64x2>` specialization).
template <> multifloats::float64x2          abs  (complex<multifloats::float64x2> const &z);
/// @brief Complex argument `arg(z)` in (-pi, pi] (`std::complex<float64x2>`).
template <> multifloats::float64x2          arg  (complex<multifloats::float64x2> const &z);
/// @brief Projection onto the Riemann sphere (`std::complex<float64x2>`).
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
/// @brief Complex `sin(pi*z)` (`std::complex<float64x2>` overload).
std::complex<float64x2> sinpi (std::complex<float64x2> const &z);
/// @brief Complex `cos(pi*z)` (`std::complex<float64x2>` overload).
std::complex<float64x2> cospi (std::complex<float64x2> const &z);
/// @brief Complex `exp(z) - 1`, accurate for small `z` (`std::complex<float64x2>` overload).
std::complex<float64x2> expm1 (std::complex<float64x2> const &z);
/// @brief Complex base-2 logarithm (`std::complex<float64x2>` overload).
std::complex<float64x2> log2  (std::complex<float64x2> const &z);
/// @brief Complex base-10 logarithm (`std::complex<float64x2>` overload).
std::complex<float64x2> log10 (std::complex<float64x2> const &z);
/// @brief Complex `log(1 + z)`, accurate for small `z` (`std::complex<float64x2>` overload).
std::complex<float64x2> log1p (std::complex<float64x2> const &z);

// cis(x) = cos(x) + i·sin(x). One fused sincos covers both components.
// Mirrors libquadmath cexpiq.
/// @brief `cos(x) + i*sin(x)`, as a `std::complex<float64x2>`.
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
