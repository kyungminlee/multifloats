#pragma once

#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <string>
#include <type_traits>

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

// =============================================================================
// float64x2 — double-double (DD) arithmetic
//
// Two IEEE-754 binary64 limbs (hi, lo). A canonical value satisfies
// |lo| <= ulp(hi)/2, giving ~104 bits of significand. Binary ops are inlined
// here; transcendentals delegate to the extern "C" `*dd` kernels defined in
// multifloats_math.cc.
// =============================================================================

class float64x2 {
public:
  double _limbs[2] = {};

  constexpr float64x2() = default;
  constexpr float64x2(float64x2 const &) = default;
  constexpr float64x2(float64x2 &&) = default;
  constexpr float64x2 &operator=(float64x2 const &) = default;
  constexpr float64x2 &operator=(float64x2 &&) = default;

  // One- and two-argument ctors in member-init form for symmetry. Non-explicit:
  // the scalar form is used implicitly throughout (e.g. `float64x2(1)` in
  // hypot/cbrt), and the two-arg form is the natural DD-literal spelling
  // `float64x2{hi, lo}` used inside class-method bodies.
  constexpr float64x2(double arg) : _limbs{arg, 0.0} {}
  constexpr float64x2(double hi, double lo) : _limbs{hi, lo} {}

  constexpr explicit operator double() const { return _limbs[0]; }

  // Lexicographic comparison. NaN is unordered: IEEE `<` is false on either
  // side, so for NaN limbs both `<` checks fall through to the next limb,
  // yielding `operator<` = false — same as the previous _lex_compare loop.
  // For +0 vs -0 the leading comparisons are both false (IEEE +0 == -0), so
  // we correctly fall through to the next limb.
  constexpr bool operator==(float64x2 const &r) const {
    return _limbs[0] == r._limbs[0] && _limbs[1] == r._limbs[1];
  }
  constexpr bool operator!=(float64x2 const &r) const { return !(*this == r); }
  constexpr bool operator<(float64x2 const &r) const {
    if (_limbs[0] < r._limbs[0]) return true;
    if (r._limbs[0] < _limbs[0]) return false;
    return _limbs[1] < r._limbs[1];
  }
  constexpr bool operator>(float64x2 const &r) const  { return  (r < *this); }
  constexpr bool operator<=(float64x2 const &r) const { return !(r < *this); }
  constexpr bool operator>=(float64x2 const &r) const { return !(*this < r); }

  constexpr float64x2 operator+() const { return *this; }
  constexpr float64x2 operator-() const { return {-_limbs[0], -_limbs[1]}; }

  // ---------------------------------------------------------------------------
  // Binary arithmetic — kernels inlined directly, translated from
  // MultiFloats.jl (mfadd / mfmul) and the Float64x2 division kernel.
  // ---------------------------------------------------------------------------

  constexpr float64x2 operator+(float64x2 const &rhs) const {
    float64x2 out;
    double s = _limbs[0] + rhs._limbs[0];
    // Non-finite: the EFT below would propagate NaN into limbs[1];
    // short-circuit and let IEEE produce the correct leading limb.
    if (!std::isfinite(s)) {
      out._limbs[0] = s;
      return out;
    }
    // When both hi limbs are zero, two_sum loses the -0 sign
    // (IEEE 754: -0 + +0 = +0 in round-to-nearest).
    if (_limbs[0] == 0.0 && rhs._limbs[0] == 0.0) {
      out._limbs[0] = s;
      out._limbs[1] = _limbs[1] + rhs._limbs[1];
      return out;
    }
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
    detail::two_sum(_limbs[0], rhs._limbs[0], a, b);
    detail::two_sum(_limbs[1], rhs._limbs[1], c, d);
    detail::fast_two_sum(a, c, a, c);
    b += d;
    b += c;
    detail::fast_two_sum(a, b, out._limbs[0], out._limbs[1]);
    return out;
  }

  constexpr float64x2 operator-(float64x2 const &rhs) const {
    float64x2 out;
    double s = _limbs[0] - rhs._limbs[0];
    if (!std::isfinite(s)) {
      out._limbs[0] = s;
      return out;
    }
    if (_limbs[0] == 0.0 && rhs._limbs[0] == 0.0) {
      out._limbs[0] = s;
      out._limbs[1] = _limbs[1] - rhs._limbs[1];
      return out;
    }
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
    detail::two_sum(_limbs[0], -rhs._limbs[0], a, b);
    detail::two_sum(_limbs[1], -rhs._limbs[1], c, d);
    detail::fast_two_sum(a, c, a, c);
    b += d;
    b += c;
    detail::fast_two_sum(a, b, out._limbs[0], out._limbs[1]);
    return out;
  }

  constexpr float64x2 operator*(float64x2 const &rhs) const {
    float64x2 out;
    double p00 = 0.0, e00 = 0.0;
    detail::two_prod(_limbs[0], rhs._limbs[0], p00, e00);
    double p01 = detail::one_prod(_limbs[0], rhs._limbs[1]);
    double p10 = detail::one_prod(_limbs[1], rhs._limbs[0]);
    p01 += p10;
    e00 += p01;
    detail::fast_two_sum(p00, e00, out._limbs[0], out._limbs[1]);
    return out;
  }

  constexpr float64x2 operator/(float64x2 const &rhs) const {
    // Dekker-style: q1 = hi/rhs.hi, refine once.
    double q1 = _limbs[0] / rhs._limbs[0];
    if (!std::isfinite(q1)) {
      // Mirror q1 into the lo limb so a non-finite result propagates
      // through both limbs. Otherwise isnan/isinf checks against the lo
      // limb would spuriously report "finite" on a NaN/Inf DD.
      float64x2 out;
      out._limbs[0] = q1;
      out._limbs[1] = q1;
      return out;
    }
    if (!std::isfinite(rhs._limbs[0])) {
      // Finite / ±Inf — q1 is ±0; the correct DD is {±0, 0}, which
      // default-initialization already gives us.
      float64x2 out;
      out._limbs[0] = q1;
      return out;
    }
    // r = this - q1 * rhs, computed as a full DD (q1 is a single-limb
    // scalar so q1*rhs is one two_prod + one one_prod = one DD).
    double p00 = 0.0, e00 = 0.0;
    detail::two_prod(q1, rhs._limbs[0], p00, e00);
    double p01 = detail::one_prod(q1, rhs._limbs[1]);
    double qhi = p00;
    double qlo = e00 + p01;
    // r = this - (qhi, qlo) via two_diff
    double r0 = 0.0, r0e = 0.0;
    detail::two_sum(_limbs[0], -qhi, r0, r0e);
    double r1 = (_limbs[1] - qlo) + r0e;
    double rh = 0.0, rl = 0.0;
    detail::fast_two_sum(r0, r1, rh, rl);
    // q2 = r.hi / rhs.hi
    double q2 = rh / rhs._limbs[0];
    float64x2 out;
    detail::fast_two_sum(q1, q2, out._limbs[0], out._limbs[1]);
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
};

// =============================================================================
// <cmath>-style free functions (ADL on float64x2)
// =============================================================================

namespace detail {
// Index of the first nonzero limb, or 2 if every limb is a (possibly signed)
// zero. Used by abs / signbit to resolve the sign of non-canonical DDs like
// (+0, -eps), where signbit(hi) alone would misclassify.
constexpr std::size_t first_nonzero_limb_index(float64x2 const &x) {
  if (x._limbs[0] != 0.0) return 0;
  if (x._limbs[1] != 0.0) return 1;
  return 2;
}

// Triple-double scratch primitives (float64x3, td_add_*, td_mul_*, …)
// live in the internal multifloats_td.hh. Only kernels inside
// src/multifloats_math.cc and the direct-primitive tests pull them in;
// they are not part of the public ABI surface exposed here.

} // namespace detail

inline constexpr float64x2 fabs(float64x2 const &x) {
  std::size_t i = detail::first_nonzero_limb_index(x);
  if (i == 2) return x;
  return std::signbit(x._limbs[i]) ? -x : x;
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
  return std::signbit(x._limbs[i == 2 ? 0 : i]);
}

inline constexpr bool isfinite(float64x2 const &x) {
  // A DD with finite hi and non-finite lo is classified non-finite — this
  // matters for the operator+ short-circuit.
  for (std::size_t i = 0; i < 2; ++i) {
    if (!std::isfinite(x._limbs[i])) return false;
  }
  return true;
}

inline constexpr bool isinf(float64x2 const &x) {
  return std::isinf(x._limbs[0]);
}

inline constexpr bool isnan(float64x2 const &x) {
  for (std::size_t i = 0; i < 2; ++i) {
    if (std::isnan(x._limbs[i])) return true;
  }
  return false;
}

inline constexpr int fpclassify(float64x2 const &x) {
  return std::fpclassify(x._limbs[0]);
}

inline constexpr float64x2 ldexp(float64x2 const &x, int n) {
  // Build the power-of-two scale once; multiplication by an exact power of
  // two is exact for every limb (no rounding, no renorm), avoiding the
  // two library calls of std::ldexp.
  double scale = std::ldexp(1.0, n);
  return {x._limbs[0] * scale, x._limbs[1] * scale};
}

inline constexpr float64x2 scalbn(float64x2 const &x, int n) {
  // POSIX alias of ldexp for FLT_RADIX == 2 (which is guaranteed by IEEE 754).
  return ldexp(x, n);
}

inline constexpr int ilogb(float64x2 const &x) {
  return std::ilogb(x._limbs[0]);
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
  double fl_hi = std::floor(x._limbs[0]);
  if (fl_hi == x._limbs[0]) {
    // hi is already an integer; floor depends on the lo limb.
    r._limbs[0] = fl_hi;
    r._limbs[1] = std::floor(x._limbs[1]);
    detail::renorm_fast(r._limbs[0], r._limbs[1]);
  } else {
    r._limbs[0] = fl_hi;
    r._limbs[1] = 0.0;
  }
  return r;
}

inline constexpr float64x2 ceil(float64x2 const &x) {
  float64x2 r;
  double cl_hi = std::ceil(x._limbs[0]);
  if (cl_hi == x._limbs[0]) {
    r._limbs[0] = cl_hi;
    r._limbs[1] = std::ceil(x._limbs[1]);
    detail::renorm_fast(r._limbs[0], r._limbs[1]);
  } else {
    r._limbs[0] = cl_hi;
    r._limbs[1] = 0.0;
  }
  return r;
}

inline constexpr float64x2 trunc(float64x2 const &x) {
  return std::signbit(x._limbs[0]) ? -floor(-x) : floor(x);
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
  double hi = std::round(x._limbs[0]);
  if (hi == x._limbs[0]) {
    double lo = x._limbs[1];
    double rlo = 0.0;
    if      (lo ==  0.5 && hi <  0.0) rlo = 0.0;
    else if (lo == -0.5 && hi >  0.0) rlo = 0.0;
    else                              rlo = std::round(lo);
    r._limbs[0] = hi;
    r._limbs[1] = rlo;
    detail::renorm_fast(r._limbs[0], r._limbs[1]);
  } else {
    double diff = x._limbs[0] - hi;
    if      (diff == -0.5 && x._limbs[1] < 0.0) hi -= 1.0;
    else if (diff ==  0.5 && x._limbs[1] > 0.0) hi += 1.0;
    r._limbs[0] = hi;
    r._limbs[1] = 0.0;
  }
  return r;
}

inline constexpr float64x2 nearbyint(float64x2 const &x) {
  float64x2 r;
  double hi = std::nearbyint(x._limbs[0]);
  if (hi == x._limbs[0]) {
    r._limbs[0] = hi;
    r._limbs[1] = std::nearbyint(x._limbs[1]);
    detail::renorm_fast(r._limbs[0], r._limbs[1]);
  } else {
    r._limbs[0] = hi;
    r._limbs[1] = 0.0;
  }
  return r;
}

inline constexpr float64x2 rint(float64x2 const &x) { return nearbyint(x); }

namespace detail {
// Shared half-integer correction for lround / llround. Given i = std::[l]lround(x_hi),
// adjust ±1 based on how lo crosses the half-integer boundary of the true value.
template <typename Int>
constexpr Int lround_adjust(float64x2 const &x, Int i) {
  double hi = x._limbs[0];
  double lo = x._limbs[1];
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
  return detail::lround_adjust<long>(x, std::lround(x._limbs[0]));
}

inline constexpr long long llround(float64x2 const &x) {
  return detail::lround_adjust<long long>(x, std::llround(x._limbs[0]));
}

inline constexpr long lrint(float64x2 const &x) {
  return std::lrint(rint(x)._limbs[0]);
}

inline constexpr long long llrint(float64x2 const &x) {
  return std::llrint(rint(x)._limbs[0]);
}

// =============================================================================
// Floating-point manipulation
// =============================================================================

inline constexpr float64x2 frexp(float64x2 const &x, int *exp) {
  float64x2 r;
  int e = 0;
  r._limbs[0] = std::frexp(x._limbs[0], &e);
  r._limbs[1] = std::ldexp(x._limbs[1], -e);
  *exp = e;
  return r;
}

inline constexpr float64x2 modf(float64x2 const &x, float64x2 *iptr) {
  *iptr = trunc(x);
  return x - *iptr;
}

inline constexpr float64x2 scalbln(float64x2 const &x, long n) {
  return {std::scalbln(x._limbs[0], n), std::scalbln(x._limbs[1], n)};
}

inline constexpr float64x2 logb(float64x2 const &x) {
  float64x2 r;
  r._limbs[0] = std::logb(x._limbs[0]);
  return r;
}

inline constexpr float64x2 nextafter(float64x2 const &x, float64x2 const &y) {
  if (x == y) return y;
  // One DD ulp ≈ ulp_up(|hi|) * 2^-53. Always use the upward ulp of |hi|;
  // the downward ulp halves at a power-of-2 boundary, so picking it there
  // would make the step 2× too small and break the round-trip identity
  // nextafter(nextafter(x, +inf), -inf) == x.
  double ax = std::abs(x._limbs[0]);
  double inf = std::numeric_limits<double>::infinity();
  double ulp = std::nextafter(ax, inf) - ax;
  double eps = std::ldexp(ulp, -53);
  return (x < y) ? x + float64x2(eps) : x - float64x2(eps);
}

inline constexpr float64x2 nexttoward(float64x2 const &x, float64x2 const &y) {
  return nextafter(x, y);
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
  bool x_neg = x._limbs[0] < 0.0;
  float64x2 ax = x_neg ? -x : x;
  float64x2 ay = (y._limbs[0] < 0.0) ? -y : y;

  if (ax < ay) return x;

  float64x2 r = ax;
  while (true) {
    if (r._limbs[0] < 0.0) {
      double r_abs = -r._limbs[0];
      int gap = std::ilogb(r_abs) - std::ilogb(ay._limbs[0]);
      if (gap <= 0) {
        r = r + ay;
      } else if (gap <= 53) {
        double q = std::trunc(r_abs / ay._limbs[0]) + 1.0;
        r = r + ay * float64x2(q);
      } else {
        r = r + (trunc(-r / ay) + float64x2(1.0)) * ay;
      }
    } else if (r >= ay) {
      int gap = std::ilogb(r._limbs[0]) - std::ilogb(ay._limbs[0]);
      if (gap <= 53) {
        double q = std::trunc(r._limbs[0] / ay._limbs[0]);
        r = (q <= 1.0) ? (r - ay) : (r - ay * float64x2(q));
      } else {
        r = r - trunc(r / ay) * ay;
      }
    } else {
      break;
    }
    if (r._limbs[0] == 0.0 && r._limbs[1] == 0.0) break;
  }

  return x_neg ? -r : r;
}

inline constexpr float64x2 remainder(float64x2 const &x, float64x2 const &y) {
  return x - round(x / y) * y;
}

inline constexpr float64x2 remquo(float64x2 const &x, float64x2 const &y, int *quo) {
  float64x2 q = round(x / y);
  *quo = static_cast<int>(q._limbs[0]);
  return x - q * y;
}

inline constexpr float64x2 fdim(float64x2 const &x, float64x2 const &y) {
  return (x > y) ? (x - y) : float64x2();
}

// C++20 std::lerp: exact at the endpoints, monotonic in t, and does not
// overshoot when a and b have the same sign and t is in [0, 1].
inline constexpr float64x2 lerp(float64x2 const &a, float64x2 const &b,
                      float64x2 const &t) {
  if ((a._limbs[0] <= 0.0 && b._limbs[0] >= 0.0) ||
      (a._limbs[0] >= 0.0 && b._limbs[0] <= 0.0)) {
    // Opposite signs (or one is zero): no cancellation risk.
    return t * b + (float64x2(1.0) - t) * a;
  }
  if (t == float64x2(1.0)) {
    return b;  // exact endpoint per C++20 spec
  }
  float64x2 x = a + t * (b - a);
  // Enforce monotonicity at the b end when t is past 1 or when rounding
  // nudges x beyond b — matches libstdc++/libc++ behavior.
  if ((t._limbs[0] > 1.0) == (b > a)) {
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
// inline definitions but have been moved to src/multifloats_math_poly.inc,
// which is pulled into namespace multifloats::detail inside
// multifloats_math.cc. Every call site lives in that TU, so inlining
// decisions are unchanged; the public header no longer carries ~250 lines
// of Estrin switch bodies.

} // namespace detail

// =============================================================================
// C ABI handoff: include the extern "C" `*dd` function declarations at
// namespace global scope, then convert between float64x2 and float64x2_t.
// =============================================================================

} // namespace multifloats
#include "multifloats_c.h"
namespace multifloats {

namespace detail {
inline constexpr float64x2_t to_f64x2(float64x2 const &x) { return {x._limbs[0], x._limbs[1]}; }
inline constexpr float64x2 from_f64x2(float64x2_t x) { float64x2 r; r._limbs[0] = x.hi; r._limbs[1] = x.lo; return r; }
} // namespace detail

// =============================================================================
// Roots
// =============================================================================

inline constexpr float64x2 sqrt(float64x2 const &x) {
  double s = std::sqrt(x._limbs[0]);
  // Bail on 0, -0, negative, NaN, +Inf — the Karp-Markstein refinement
  // would compute `inf - inf = NaN` in the residual step for +Inf, and
  // `0 / 0` for 0. Leading-limb sqrt handles every IEEE special case.
  if (!(x._limbs[0] > 0.0) || !std::isfinite(s)) {
    float64x2 r;
    r._limbs[0] = s;
    return r;
  }
  // Karp/Markstein: r = s + (x - s*s) / (2s), evaluated in DD. The
  // correction reduces the DD residual to a scalar via
  // `residual._limbs[0] * (0.5/s)`, so the residual's lo limb is
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
  const float64x2 correction(residual._limbs[0] * (0.5 / s));
  return s_dd + correction;
}

inline constexpr float64x2 cbrt(float64x2 const &x) {
  if (x._limbs[0] == 0.0) return float64x2();
  const double s = std::cbrt(x._limbs[0]);
  const float64x2 s_dd(s);
  const float64x2 residual = x - s_dd * s_dd * s_dd;
  const float64x2 correction(residual._limbs[0] / (3.0 * s * s));
  return s_dd + correction;
}

inline constexpr float64x2 hypot(float64x2 const &x, float64x2 const &y) {
  // Defer to libm's hypot for non-finite so inf/NaN propagate correctly.
  if (!std::isfinite(x._limbs[0]) || !std::isfinite(y._limbs[0])) {
    float64x2 r;
    r._limbs[0] = std::hypot(x._limbs[0], y._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  float64x2 ax = signbit(x) ? -x : x;
  float64x2 ay = signbit(y) ? -y : y;
  float64x2 big = (ax > ay) ? ax : ay;
  float64x2 small = (ax > ay) ? ay : ax;
  if (big._limbs[0] == 0.0) return float64x2();
  // Power-of-2 scale (exact) so big has exponent 0 before the square.
  // Replace per-limb ldexp calls with multiplies by 2^(±e): for a power-
  // of-2 multiplier and a non-subnormal result, `x * 2^k` and
  // `ldexp(x, k)` are bit-identical, but the multiply is one FP op while
  // ldexp is a libm call. `e` comes from `ilogb(big.hi)` on a finite
  // non-zero input, so `|e| ≤ 1023`; `2^(-e)` and `2^e` are both finite.
  int e = std::ilogb(big._limbs[0]);
  double down = std::ldexp(1.0, -e);
  big._limbs[0]   *= down; big._limbs[1]   *= down;
  small._limbs[0] *= down; small._limbs[1] *= down;
  float64x2 ratio = small / big;
  float64x2 root = big * sqrt(float64x2(1.0) + ratio * ratio);
  float64x2 r;
  double up = std::ldexp(1.0, e);
  r._limbs[0] = root._limbs[0] * up;
  r._limbs[1] = root._limbs[1] * up;
  // Overflow: true result exceeds double's range. Zero the trailing limb so
  // callers see a clean inf rather than (inf, NaN).
  if (!std::isfinite(r._limbs[0])) r._limbs[1] = 0.0;
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

// Cylindrical Bessel. Names follow C++17 <cmath>: cyl_bessel_j for J_n,
// cyl_neumann for Y_n. Order-0 and order-1 have dedicated fast-path
// kernels; the integer-order dispatchers cyl_bessel_j/cyl_neumann call
// them internally and run Miller / forward recurrence for |n| ≥ 2. The
// 4-arg cyl_neumann overload fills a single forward-recurrence sweep of
// n2-n1+1 outputs (cheaper than n2-n1+1 independent 2-arg calls).
float64x2 cyl_bessel_j0 (float64x2 const &x);
float64x2 cyl_bessel_j1 (float64x2 const &x);
float64x2 cyl_neumann0  (float64x2 const &x);
float64x2 cyl_neumann1  (float64x2 const &x);
float64x2 cyl_bessel_j  (int n, float64x2 const &x);
float64x2 cyl_neumann   (int n, float64x2 const &x);
void      cyl_neumann   (int n1, int n2, float64x2 const &x, float64x2 *out);

// =============================================================================
// Additional classification and ordered comparison
// =============================================================================

inline constexpr bool isnormal(float64x2 const &x) {
  return std::isnormal(x._limbs[0]);
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

// ---- Formatted I/O for float64x2 -------------------------------------------
//
// Scientific-notation decimal output at up to 34 significant digits (the
// natural DD significand limit is ~32; we allow two guard digits for
// round-half-to-even). Special values use the same textual forms as
// std::to_string for doubles ("nan", "inf", "-inf"). `operator<<` honors
// `os.precision()` when > 17; otherwise (including the C++ default of 6)
// the DD default of 32 digits is used. Definitions live in
// src/multifloats_io.cc — inlining is not worth the header weight.
std::string to_string(float64x2 const &x, int precision = 32);
std::ostream &operator<<(std::ostream &os, float64x2 const &x);

} // namespace multifloats

// ---- std::complex<multifloats::float64x2> specializations ------------------
//
// For `exp, sin, cos, tan, sinh, cosh, tanh, atanh, acos` the generic
// <complex> template path is either slow (each call pair sin/cos or
// sinh/cosh goes through two separate extern-C kernels; the compiler
// can't fuse across the ABI boundary), or in the case of `acos`
// incorrect at DD precision on platforms where `long double == double`
// (libstdc++'s `(_Tp)1.5707963…L` truncates the π/2 low limb).
//
// These specializations delegate to the `c*dd` symbols in
// multifloats_math.cc, which use the fused sincos_full /
// sinhcosh_full kernels and hard-code the DD π/2 constant.
//
// The remaining eight functions (log, log10, pow, sqrt, asin, atan,
// asinh, acosh) are also exported via c*dd for Fortran / C callers,
// but we leave `std::log(complex<float64x2>)` etc. on the generic template:
// those paths offer no speedup and no correctness advantage here.
#include <complex>

namespace std {

#define MULTIFLOATS_CX_SPECIALIZE(fn)                                        \
  template <>                                                                \
  inline complex<multifloats::float64x2>                                     \
  fn(complex<multifloats::float64x2> const &z) {                             \
    ::complex64x2_t in = {multifloats::detail::to_f64x2(z.real()),           \
                          multifloats::detail::to_f64x2(z.imag())};          \
    ::complex64x2_t out = ::c##fn##dd(in);                                   \
    return complex<multifloats::float64x2>(                                  \
        multifloats::detail::from_f64x2(out.re),                             \
        multifloats::detail::from_f64x2(out.im));                            \
  }

MULTIFLOATS_CX_SPECIALIZE(exp)
MULTIFLOATS_CX_SPECIALIZE(sin)
MULTIFLOATS_CX_SPECIALIZE(cos)
MULTIFLOATS_CX_SPECIALIZE(tan)
MULTIFLOATS_CX_SPECIALIZE(sinh)
MULTIFLOATS_CX_SPECIALIZE(cosh)
MULTIFLOATS_CX_SPECIALIZE(tanh)
MULTIFLOATS_CX_SPECIALIZE(atanh)
MULTIFLOATS_CX_SPECIALIZE(acos)

#undef MULTIFLOATS_CX_SPECIALIZE

// Real-returning / identity specializations that delegate to the
// matching c*dd symbol. `abs` and `arg` match libquadmath's overflow-
// safe hypot / atan2 paths. `proj` handles the Riemann-sphere case.
template <>
inline multifloats::float64x2
abs(complex<multifloats::float64x2> const &z) {
  ::complex64x2_t in = {multifloats::detail::to_f64x2(z.real()),
                        multifloats::detail::to_f64x2(z.imag())};
  return multifloats::detail::from_f64x2(::cabsdd(in));
}

template <>
inline multifloats::float64x2
arg(complex<multifloats::float64x2> const &z) {
  ::complex64x2_t in = {multifloats::detail::to_f64x2(z.real()),
                        multifloats::detail::to_f64x2(z.imag())};
  return multifloats::detail::from_f64x2(::cargdd(in));
}

template <>
inline complex<multifloats::float64x2>
proj(complex<multifloats::float64x2> const &z) {
  ::complex64x2_t in = {multifloats::detail::to_f64x2(z.real()),
                        multifloats::detail::to_f64x2(z.imag())};
  ::complex64x2_t out = ::cprojdd(in);
  return complex<multifloats::float64x2>(
      multifloats::detail::from_f64x2(out.re),
      multifloats::detail::from_f64x2(out.im));
}

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
} // namespace multifloats
