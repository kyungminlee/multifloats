#pragma once

#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <ostream>
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
// Error-free transformations (inputs by value to avoid aliasing pitfalls)
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

template <typename T, std::size_t N> struct MultiFloat {
  static_assert(N == 1 || N == 2, "only N = 1 and N = 2 are implemented");
  T _limbs[N] = {};

  constexpr MultiFloat() = default;
  constexpr MultiFloat(MultiFloat const &) = default;
  constexpr MultiFloat(MultiFloat &&) = default;
  constexpr MultiFloat &operator=(MultiFloat const &) = default;
  constexpr MultiFloat &operator=(MultiFloat &&) = default;

  constexpr MultiFloat(T const &arg) { _limbs[0] = arg; }

  // Initialize every limb from a brace-init list: caller supplies exactly
  // N values, in order (hi, mid, ..., lo). No renormalization — the caller
  // is responsible for passing a canonical pair (|lo| <= ulp(hi)/2 for
  // N=2). Intended as a cheap factory for pre-split DD constants; use
  // arithmetic operators (or build from a scalar) if renormalization is
  // required.
  template <typename... Us,
            typename = std::enable_if_t<sizeof...(Us) == N &&
                                        (std::is_convertible_v<Us, T> && ...)>>
  constexpr MultiFloat(Us... args) : _limbs{T(args)...} {}

  constexpr explicit operator T() const { return _limbs[0]; }

  // Lexicographic three-way limb comparison. Returns -1 / 0 / +1 on the
  // first limb where `<` (or its swap) holds; NaN on either side yields 0
  // (fall-through — matching the original per-operator loops, which
  // preserves the existing DD-level "unordered" semantics for NaN).
  constexpr int _lex_compare(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (_limbs[i] < rhs._limbs[i]) return -1;
      if (rhs._limbs[i] < _limbs[i]) return +1;
    }
    return 0;
  }

  constexpr bool _equal_limbs(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (!(_limbs[i] == rhs._limbs[i])) return false;
    }
    return true;
  }

  constexpr bool operator==(MultiFloat const &rhs) const { return _equal_limbs(rhs); }
  constexpr bool operator!=(MultiFloat const &rhs) const { return !_equal_limbs(rhs); }
  constexpr bool operator<(MultiFloat const &rhs) const { return _lex_compare(rhs) < 0; }
  constexpr bool operator>(MultiFloat const &rhs) const { return _lex_compare(rhs) > 0; }
  constexpr bool operator<=(MultiFloat const &rhs) const { return _lex_compare(rhs) <= 0; }
  constexpr bool operator>=(MultiFloat const &rhs) const { return _lex_compare(rhs) >= 0; }

  constexpr MultiFloat operator+() const { return *this; }

  constexpr MultiFloat operator-() const {
    MultiFloat r;
    for (std::size_t i = 0; i < N; ++i) {
      r._limbs[i] = -_limbs[i];
    }
    return r;
  }

  // ---------------------------------------------------------------------------
  // Binary arithmetic — kernels inlined directly, translated from
  // MultiFloats.jl (mfadd / mfmul) and the Float64x2 division kernel in
  // fsrc/multifloats.fypp.
  // ---------------------------------------------------------------------------

  constexpr MultiFloat operator+(MultiFloat const &rhs) const {
    MultiFloat out;
    if constexpr (N == 1) {
      out._limbs[0] = _limbs[0] + rhs._limbs[0];
    } else { // N == 2
      T s = _limbs[0] + rhs._limbs[0];
      // Non-finite: the EFT below would propagate NaN into limbs[1];
      // short-circuit and let IEEE produce the correct leading limb.
      if (!std::isfinite(s)) {
        out._limbs[0] = s;
        return out;
      }
      // When both hi limbs are zero, two_sum loses the -0 sign
      // (IEEE 754: -0 + +0 = +0 in round-to-nearest).
      if (_limbs[0] == T(0) && rhs._limbs[0] == T(0)) {
        out._limbs[0] = s;
        out._limbs[1] = _limbs[1] + rhs._limbs[1];
        return out;
      }
      T a, b, c, d;
      detail::two_sum(_limbs[0], rhs._limbs[0], a, b);
      detail::two_sum(_limbs[1], rhs._limbs[1], c, d);
      detail::fast_two_sum(a, c, a, c);
      b += d;
      b += c;
      detail::fast_two_sum(a, b, out._limbs[0], out._limbs[1]);
    }
    return out;
  }

  constexpr MultiFloat operator-(MultiFloat const &rhs) const {
    MultiFloat out;
    if constexpr (N == 1) {
      out._limbs[0] = _limbs[0] - rhs._limbs[0];
    } else { // N == 2 — dedicated two_diff, mirrors operator+
      T s = _limbs[0] - rhs._limbs[0];
      if (!std::isfinite(s)) {
        out._limbs[0] = s;
        return out;
      }
      if (_limbs[0] == T(0) && rhs._limbs[0] == T(0)) {
        out._limbs[0] = s;
        out._limbs[1] = _limbs[1] - rhs._limbs[1];
        return out;
      }
      T a, b, c, d;
      detail::two_sum(_limbs[0], -rhs._limbs[0], a, b);
      detail::two_sum(_limbs[1], -rhs._limbs[1], c, d);
      detail::fast_two_sum(a, c, a, c);
      b += d;
      b += c;
      detail::fast_two_sum(a, b, out._limbs[0], out._limbs[1]);
    }
    return out;
  }

  constexpr MultiFloat operator*(MultiFloat const &rhs) const {
    MultiFloat out;
    if constexpr (N == 1) {
      out._limbs[0] = _limbs[0] * rhs._limbs[0];
    } else { // N == 2
      T p00, e00;
      detail::two_prod(_limbs[0], rhs._limbs[0], p00, e00);
      T p01 = detail::one_prod(_limbs[0], rhs._limbs[1]);
      T p10 = detail::one_prod(_limbs[1], rhs._limbs[0]);
      p01 += p10;
      e00 += p01;
      detail::fast_two_sum(p00, e00, out._limbs[0], out._limbs[1]);
    }
    return out;
  }

  constexpr MultiFloat operator/(MultiFloat const &rhs) const {
    if constexpr (N == 1) {
      MultiFloat out;
      out._limbs[0] = _limbs[0] / rhs._limbs[0];
      return out;
    } else { // N == 2 — Dekker-style: q1 = hi/rhs.hi, refine once.
      T q1 = _limbs[0] / rhs._limbs[0];
      if (!std::isfinite(q1)) {
        // Mirror q1 into the lo limb so a non-finite result propagates
        // through both limbs. Otherwise isnan/isinf checks against the lo
        // limb would spuriously report "finite" on a NaN/Inf DD.
        MultiFloat out;
        out._limbs[0] = q1;
        out._limbs[1] = q1;
        return out;
      }
      if (!std::isfinite(rhs._limbs[0])) {
        // Finite / ±Inf — q1 is ±0; the correct DD is {±0, 0}, which
        // default-initialization already gives us.
        MultiFloat out;
        out._limbs[0] = q1;
        return out;
      }
      // r = this - q1 * rhs, computed as a full DD (q1 is a single-limb
      // scalar so q1*rhs is one two_prod + one one_prod = one DD).
      T p00, e00;
      detail::two_prod(q1, rhs._limbs[0], p00, e00);
      T p01 = detail::one_prod(q1, rhs._limbs[1]);
      T qhi = p00;
      T qlo = e00 + p01;
      // r = this - (qhi, qlo) via two_diff
      T r0, r0e;
      detail::two_sum(_limbs[0], -qhi, r0, r0e);
      T r1 = (_limbs[1] - qlo) + r0e;
      T rh, rl;
      detail::fast_two_sum(r0, r1, rh, rl);
      // q2 = r.hi / rhs.hi
      T q2 = rh / rhs._limbs[0];
      MultiFloat out;
      detail::fast_two_sum(q1, q2, out._limbs[0], out._limbs[1]);
      return out;
    }
  }

  constexpr MultiFloat &operator+=(MultiFloat const &rhs) {
    return *this = *this + rhs;
  }
  constexpr MultiFloat &operator-=(MultiFloat const &rhs) {
    return *this = *this - rhs;
  }
  constexpr MultiFloat &operator*=(MultiFloat const &rhs) {
    return *this = *this * rhs;
  }
  constexpr MultiFloat &operator/=(MultiFloat const &rhs) {
    return *this = *this / rhs;
  }
};

using float64x2 = MultiFloat<double, 2>;

// =============================================================================
// <cmath>-style free functions (ADL on MultiFloat)
// =============================================================================

namespace detail {
// Index of the first nonzero limb, or N if every limb is a (possibly
// signed) zero. Used by abs / signbit to resolve the sign of non-canonical
// DDs like (+0, -eps), where signbit(hi) alone would misclassify.
template <typename T, std::size_t N>
constexpr std::size_t first_nonzero_limb_index(MultiFloat<T, N> const &x) {
  for (std::size_t i = 0; i < N; ++i) {
    if (x._limbs[i] != T(0)) return i;
  }
  return N;
}
} // namespace detail

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> abs(MultiFloat<T, N> const &x) {
  std::size_t i = detail::first_nonzero_limb_index(x);
  if (i == N) return x;
  return std::signbit(x._limbs[i]) ? -x : x;
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> fabs(MultiFloat<T, N> const &x) {
  return abs(x);
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> fmin(MultiFloat<T, N> const &a,
                                MultiFloat<T, N> const &b) {
  return (a < b) ? a : b;
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> fmax(MultiFloat<T, N> const &a,
                                MultiFloat<T, N> const &b) {
  return (a < b) ? b : a;
}

template <typename T, std::size_t N>
constexpr bool signbit(MultiFloat<T, N> const &x) {
  // For non-canonical zero-hi DDs (e.g. (+0, -eps)), the sign lives in
  // the first nonzero limb. Fall through to signbit(hi) when every limb
  // is a (possibly signed) zero, preserving IEEE -0 semantics.
  std::size_t i = detail::first_nonzero_limb_index(x);
  return std::signbit(x._limbs[i == N ? 0 : i]);
}

template <typename T, std::size_t N>
constexpr bool isfinite(MultiFloat<T, N> const &x) {
  for (std::size_t i = 0; i < N; ++i) {
    if (!std::isfinite(x._limbs[i])) {
      return false;
    }
  }
  return true;
}

template <typename T, std::size_t N>
constexpr bool isinf(MultiFloat<T, N> const &x) {
  return std::isinf(x._limbs[0]);
}

template <typename T, std::size_t N>
constexpr bool isnan(MultiFloat<T, N> const &x) {
  for (std::size_t i = 0; i < N; ++i) {
    if (std::isnan(x._limbs[i])) {
      return true;
    }
  }
  return false;
}

template <typename T, std::size_t N>
constexpr int fpclassify(MultiFloat<T, N> const &x) {
  return std::fpclassify(x._limbs[0]);
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> ldexp(MultiFloat<T, N> const &x, int n) {
  MultiFloat<T, N> r;
  // Build the power-of-two scale once; multiplication by an exact power of
  // two is exact for every limb (no rounding, no renorm), avoiding the
  // two library calls of std::ldexp.
  T scale = std::ldexp(T(1), n);
  for (std::size_t i = 0; i < N; ++i) {
    r._limbs[i] = x._limbs[i] * scale;
  }
  return r;
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> scalbn(MultiFloat<T, N> const &x, int n) {
  // POSIX alias of ldexp for FLT_RADIX == 2 (which is guaranteed by IEEE 754).
  return ldexp(x, n);
}

template <typename T, std::size_t N>
constexpr int ilogb(MultiFloat<T, N> const &x) {
  return std::ilogb(x._limbs[0]);
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> copysign(MultiFloat<T, N> const &x,
                                    MultiFloat<T, N> const &y) {
  return (signbit(x) == signbit(y)) ? x : -x;
}

namespace detail {

// fast_two_sum (assumes |hi| >= |lo|), in-place renormalization helper.
template <typename T>
constexpr void renorm_fast(T &hi, T &lo) {
  T s = hi + lo;
  T b = s - hi;
  T e = lo - b;
  hi = s;
  lo = e;
}

} // namespace detail

// =============================================================================
// Rounding and integer-valued functions
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> floor(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::floor(x._limbs[0]);
  } else { // N == 2
    T fl_hi = std::floor(x._limbs[0]);
    if (fl_hi == x._limbs[0]) {
      // hi is already an integer; floor depends on the lo limb.
      r._limbs[0] = fl_hi;
      r._limbs[1] = std::floor(x._limbs[1]);
      detail::renorm_fast(r._limbs[0], r._limbs[1]);
    } else {
      r._limbs[0] = fl_hi;
      r._limbs[1] = T(0);
    }
  }
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> ceil(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::ceil(x._limbs[0]);
  } else {
    T cl_hi = std::ceil(x._limbs[0]);
    if (cl_hi == x._limbs[0]) {
      r._limbs[0] = cl_hi;
      r._limbs[1] = std::ceil(x._limbs[1]);
      detail::renorm_fast(r._limbs[0], r._limbs[1]);
    } else {
      r._limbs[0] = cl_hi;
      r._limbs[1] = T(0);
    }
  }
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> trunc(MultiFloat<T, N> const &x) {
  return std::signbit(x._limbs[0]) ? -floor(-x) : floor(x);
}

template <typename T, std::size_t N>
MultiFloat<T, N> round(MultiFloat<T, N> const &x) {
  // Round half away from zero, matching std::round. Two half-integer
  // hazards handled here:
  //   * hi itself half-integer (e.g. 2.5): std::round jumps away from zero,
  //     undone when lo lies on the other side of the half-boundary.
  //   * hi exact integer with lo == ±0.5 (possible once ulp(hi) ≥ 1, i.e.
  //     |hi| ≥ 2^53): if sign(lo) opposes sign(hi), the true value is
  //     closer to zero, so the correct rounded value is hi itself rather
  //     than hi ± 1 that std::round(lo) would add.
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::round(x._limbs[0]);
  } else {
    T hi = std::round(x._limbs[0]);
    if (hi == x._limbs[0]) {
      T lo = x._limbs[1];
      T rlo;
      if      (lo == T( 0.5) && hi <  T(0)) rlo = T(0);
      else if (lo == T(-0.5) && hi >  T(0)) rlo = T(0);
      else                                  rlo = std::round(lo);
      r._limbs[0] = hi;
      r._limbs[1] = rlo;
      detail::renorm_fast(r._limbs[0], r._limbs[1]);
    } else {
      T diff = x._limbs[0] - hi;
      if (diff == T(-0.5) && x._limbs[1] < T(0)) hi -= T(1);
      else if (diff == T(0.5) && x._limbs[1] > T(0)) hi += T(1);
      r._limbs[0] = hi;
      r._limbs[1] = T(0);
    }
  }
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> nearbyint(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::nearbyint(x._limbs[0]);
  } else {
    T hi = std::nearbyint(x._limbs[0]);
    if (hi == x._limbs[0]) {
      r._limbs[0] = hi;
      r._limbs[1] = std::nearbyint(x._limbs[1]);
      detail::renorm_fast(r._limbs[0], r._limbs[1]);
    } else {
      r._limbs[0] = hi;
      r._limbs[1] = T(0);
    }
  }
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> rint(MultiFloat<T, N> const &x) {
  return nearbyint(x);
}

namespace detail {
// Shared half-integer correction for lround / llround. Given i = std::[l]lround(x_hi),
// adjust ±1 based on how lo crosses the half-integer boundary of the true value.
template <typename Int, typename T, std::size_t N>
Int lround_adjust(MultiFloat<T, N> const &x, Int i) {
  if constexpr (N == 2) {
    T hi = x._limbs[0];
    T lo = x._limbs[1];
    T diff = hi - T(i);
    if (diff == T(0)) {
      // hi exact integer; lo (bounded by ulp(hi)/2) decides.
      if      (lo >  T( 0.5))                 ++i;
      else if (lo <  T(-0.5))                 --i;
      else if (lo == T( 0.5) && hi >= T(0))   ++i;
      else if (lo == T(-0.5) && hi <= T(0))   --i;
    } else if (diff == T(-0.5) && lo < T(0)) --i;
    else if   (diff == T( 0.5) && lo > T(0)) ++i;
  }
  return i;
}
} // namespace detail

template <typename T, std::size_t N>
long lround(MultiFloat<T, N> const &x) {
  return detail::lround_adjust<long>(x, std::lround(x._limbs[0]));
}

template <typename T, std::size_t N>
long long llround(MultiFloat<T, N> const &x) {
  return detail::lround_adjust<long long>(x, std::llround(x._limbs[0]));
}

template <typename T, std::size_t N>
long lrint(MultiFloat<T, N> const &x) {
  return std::lrint(rint(x)._limbs[0]);
}

template <typename T, std::size_t N>
long long llrint(MultiFloat<T, N> const &x) {
  return std::llrint(rint(x)._limbs[0]);
}

// =============================================================================
// Floating-point manipulation
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> frexp(MultiFloat<T, N> const &x, int *exp) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::frexp(x._limbs[0], exp);
  } else {
    int e;
    r._limbs[0] = std::frexp(x._limbs[0], &e);
    r._limbs[1] = std::ldexp(x._limbs[1], -e);
    *exp = e;
  }
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> modf(MultiFloat<T, N> const &x, MultiFloat<T, N> *iptr) {
  *iptr = trunc(x);
  return x - *iptr;
}

template <typename T, std::size_t N>
MultiFloat<T, N> scalbln(MultiFloat<T, N> const &x, long n) {
  MultiFloat<T, N> r;
  for (std::size_t i = 0; i < N; ++i) {
    r._limbs[i] = std::scalbln(x._limbs[i], n);
  }
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> logb(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  r._limbs[0] = std::logb(x._limbs[0]);
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> nextafter(MultiFloat<T, N> const &x,
                           MultiFloat<T, N> const &y) {
  if (x == y) {
    return y;
  }
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::nextafter(x._limbs[0], y._limbs[0]);
  } else {
    // One DD ulp ≈ ulp_up(|hi|) * 2^-53.  Always use the upward ulp
    // of |hi|; the downward ulp halves at a power-of-2 boundary, so
    // picking it there would make the step 2× too small and break
    // the round-trip identity nextafter(nextafter(x, +inf), -inf) == x.
    T ax = std::abs(x._limbs[0]);
    T inf = std::numeric_limits<T>::infinity();
    T ulp = std::nextafter(ax, inf) - ax;
    T eps = std::ldexp(ulp, -53);
    return (x < y) ? x + MultiFloat<T, N>(eps) : x - MultiFloat<T, N>(eps);
  }
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> nexttoward(MultiFloat<T, N> const &x,
                            MultiFloat<T, N> const &y) {
  return nextafter(x, y);
}

// =============================================================================
// Basic arithmetic helpers
// =============================================================================

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> fma(MultiFloat<T, N> const &x,
                               MultiFloat<T, N> const &y,
                               MultiFloat<T, N> const &z) {
  // Not a hardware fma, but provides the cmath interface.
  return x * y + z;
}

template <typename T, std::size_t N>
MultiFloat<T, N> fmod(MultiFloat<T, N> const &x, MultiFloat<T, N> const &y) {
  // Reduction step picks q from the ilogb gap between r and ay:
  // gap ≤ 53 — scalar q fits in one double (the earlier all-scalar form
  // silently lost q's low integer bits past 2^53); gap > 53 — DD-level
  // trunc(r/ay) carries ~106 integer bits, so a single step drops gap by
  // ≥ 53 and iteration converges in O(gap/53) steps even past 2^106. DD
  // rounding of r − q·ay can leave a tiny negative residue; add-back
  // uses the same gap dispatch so recovery stays O(1).
  bool x_neg = x._limbs[0] < T(0);
  MultiFloat<T, N> ax = x_neg ? -x : x;
  MultiFloat<T, N> ay = (y._limbs[0] < T(0)) ? -y : y;

  if (ax < ay) return x;

  MultiFloat<T, N> r = ax;
  while (true) {
    if (r._limbs[0] < T(0)) {
      T r_abs = -r._limbs[0];
      int gap = std::ilogb(r_abs) - std::ilogb(ay._limbs[0]);
      if (gap <= 0) {
        r = r + ay;
      } else if (gap <= 53) {
        T q = std::trunc(r_abs / ay._limbs[0]) + T(1);
        r = r + ay * MultiFloat<T, N>(q);
      } else {
        r = r + (trunc(-r / ay) + MultiFloat<T, N>(T(1))) * ay;
      }
    } else if (r >= ay) {
      int gap = std::ilogb(r._limbs[0]) - std::ilogb(ay._limbs[0]);
      if (gap <= 53) {
        T q = std::trunc(r._limbs[0] / ay._limbs[0]);
        r = (q <= T(1)) ? (r - ay) : (r - ay * MultiFloat<T, N>(q));
      } else {
        r = r - trunc(r / ay) * ay;
      }
    } else {
      break;
    }
    if (r._limbs[0] == T(0) && r._limbs[1] == T(0)) break;
  }

  return x_neg ? -r : r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> remainder(MultiFloat<T, N> const &x,
                           MultiFloat<T, N> const &y) {
  return x - round(x / y) * y;
}

template <typename T, std::size_t N>
MultiFloat<T, N> remquo(MultiFloat<T, N> const &x, MultiFloat<T, N> const &y,
                        int *quo) {
  MultiFloat<T, N> q = round(x / y);
  *quo = static_cast<int>(q._limbs[0]);
  return x - q * y;
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> fdim(MultiFloat<T, N> const &x,
                                MultiFloat<T, N> const &y) {
  return (x > y) ? (x - y) : MultiFloat<T, N>();
}

// C++20 std::lerp: exact at the endpoints, monotonic in t, and does not
// overshoot when a and b have the same sign and t is in [0, 1].
template <typename T, std::size_t N>
MultiFloat<T, N> lerp(MultiFloat<T, N> const &a, MultiFloat<T, N> const &b,
                      MultiFloat<T, N> const &t) {
  if ((a._limbs[0] <= T(0) && b._limbs[0] >= T(0)) ||
      (a._limbs[0] >= T(0) && b._limbs[0] <= T(0))) {
    // Opposite signs (or one is zero): no cancellation risk.
    return t * b + (MultiFloat<T, N>(T(1)) - t) * a;
  }
  if (t == MultiFloat<T, N>(T(1))) {
    return b;  // exact endpoint per C++20 spec
  }
  MultiFloat<T, N> x = a + t * (b - a);
  // Enforce monotonicity at the b end when t is past 1 or when rounding
  // nudges x beyond b — matches libstdc++/libc++ behavior.
  if ((t._limbs[0] > T(1)) == (b > a)) {
    return (b > x) ? b : x;
  }
  return (x > b) ? b : x;
}

// =============================================================================
// detail:: inline helpers for DD polynomial evaluation
//
// These are used both by the extern "C" implementations in multifloats_math.cc
// and by the header-only template helpers (sqrt, trunc, etc.).
// =============================================================================

// Forward declaration so detail kernels below can use multifloats::sqrt via ADL.
template <typename T, std::size_t N>
MultiFloat<T, N> sqrt(MultiFloat<T, N> const &x);

namespace detail {

// Polynomial and conversion constants needed by inline helpers below are
// provided through dd_constants.hh which is included only in the .cc file.
// The erf/erfc rational-approximation constants also live in the .cc file.

// ---- Horner polynomial evaluation in DD ------------------------------------
inline float64x2 horner(float64x2 const &y, double const *hi, double const *lo,
                       int n) {
  float64x2 p(hi[n - 1], lo[n - 1]);
  for (int i = n - 2; i >= 0; --i) {
    p = p * y + float64x2(hi[i], lo[i]);
  }
  return p;
}

// ---- Estrin polynomial evaluation (neval/deval from libquadmath) ------------
// neval: evaluates P[n]*x^n + P[n-1]*x^(n-1) + ... + P[0]
// Uses Estrin's scheme for ILP when degree is known, Horner fallback otherwise.
//
// NOTE on convention: neval's `n` is the polynomial *degree*, so the
// coefficient array has n+1 entries (c[0..n]). This differs from
// `horner` above, which takes `n` as the coefficient COUNT. A common
// bug when migrating Horner → Estrin is forgetting to decrement n by 1.
//
// Tier 3 measurements (cpp_bench, relative to Horner baseline):
//   exp/exp2 (case 13), sin/cos (case 12), sinh/cosh (case 8),
//   asinh/atanh (case 14): 1.78×–1.88× speedup with no precision loss
//   outside 1 DD ULP. Extending with an x^16 case would benefit
//   expm1_taylor (degree 24) and log1p_taylor (degree 17); that
//   extension was skipped as the additional gain (~1.5× on small-|x|
//   branches only) did not justify the code size.
inline float64x2 neval(float64x2 const &x, double const *hi, double const *lo,
                      int n) {
  auto c = [&](int i) { return float64x2(hi[i], lo[i]); };
  float64x2 x2 = x * x;
  float64x2 x4 = x2 * x2;
  float64x2 x8 = x4 * x4;
  switch (n) {
  case 7: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4;
  }
  case 8: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4 + c(8) * x8;
  }
  case 9: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4 + p89 * x8;
  }
  case 10: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p810 = p89 + c(10) * x2;
    return p03 + p47 * x4 + p810 * x8;
  }
  case 11: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    return p03 + p47 * x4 + p811 * x8;
  }
  case 12: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    return p03 + p47 * x4 + (p811 + c(12) * x4) * x8;
  }
  case 13: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p1213 = c(12) + c(13) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    return p03 + p47 * x4 + (p811 + p1213 * x4) * x8;
  }
  case 14: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p1213 = c(12) + c(13) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    float64x2 p1214 = p1213 + c(14) * x2;
    return p03 + p47 * x4 + (p811 + p1214 * x4) * x8;
  }
  case 15: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p1213 = c(12) + c(13) * x;
    float64x2 p1415 = c(14) + c(15) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    float64x2 p1215 = p1213 + p1415 * x2;
    return p03 + p47 * x4 + (p811 + p1215 * x4) * x8;
  }
  case 17: {
    float64x2 x16 = x8 * x8;
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p1213 = c(12) + c(13) * x;
    float64x2 p1415 = c(14) + c(15) * x;
    float64x2 p1617 = c(16) + c(17) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    float64x2 p1215 = p1213 + p1415 * x2;
    float64x2 p0_15 = p03 + p47 * x4 + (p811 + p1215 * x4) * x8;
    return p0_15 + p1617 * x16;
  }
  case 24: {
    float64x2 x16 = x8 * x8;
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p1213 = c(12) + c(13) * x;
    float64x2 p1415 = c(14) + c(15) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    float64x2 p1215 = p1213 + p1415 * x2;
    float64x2 p0_15 = p03 + p47 * x4 + (p811 + p1215 * x4) * x8;
    float64x2 p1617 = c(16) + c(17) * x;
    float64x2 p1819 = c(18) + c(19) * x;
    float64x2 p2021 = c(20) + c(21) * x;
    float64x2 p2223 = c(22) + c(23) * x;
    float64x2 p1619 = p1617 + p1819 * x2;
    float64x2 p2023 = p2021 + p2223 * x2;
    float64x2 p16_24 = p1619 + p2023 * x4 + c(24) * x8;
    return p0_15 + p16_24 * x16;
  }
  default: {
    float64x2 y = c(n);
    for (int i = n - 1; i >= 0; --i)
      y = y * x + c(i);
    return y;
  }
  }
}

// deval: evaluates x^(n+1) + P[n]*x^n + ... + P[0] (monic leading term)
inline float64x2 deval(float64x2 const &x, double const *hi, double const *lo,
                      int n) {
  auto c = [&](int i) { return float64x2(hi[i], lo[i]); };
  float64x2 x2 = x * x;
  float64x2 x4 = x2 * x2;
  float64x2 x8 = x4 * x4;
  switch (n) {
  case 7: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4 + x8;
  }
  case 8: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4 + p89 * x8;
  }
  case 9: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p810 = p89 + x2;
    return p03 + p47 * x4 + p810 * x8;
  }
  case 10: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    return p03 + p47 * x4 + p811 * x8;
  }
  case 11: {
    float64x2 p01 = c(0) + c(1) * x;
    float64x2 p23 = c(2) + c(3) * x;
    float64x2 p45 = c(4) + c(5) * x;
    float64x2 p67 = c(6) + c(7) * x;
    float64x2 p89 = c(8) + c(9) * x;
    float64x2 p1011 = c(10) + c(11) * x;
    float64x2 p03 = p01 + p23 * x2;
    float64x2 p47 = p45 + p67 * x2;
    float64x2 p811 = p89 + p1011 * x2;
    return p03 + p47 * x4 + (p811 + x4) * x8;
  }
  default: {
    float64x2 y = x + c(n);
    for (int i = n - 1; i >= 0; --i)
      y = y * x + c(i);
    return y;
  }
  }
}

} // namespace detail

// =============================================================================
// C-ABI function declarations and DD conversion helpers
//
// All DD math functions are defined as extern "C" in multifloats_math.cc.
// The C++ templates below call these functions for MultiFloat<double, 2>.
// =============================================================================

} // namespace multifloats
#include "multifloats_c.h"
namespace multifloats {

namespace detail {
inline float64x2_t to_f64x2(float64x2 const &x) { return {x._limbs[0], x._limbs[1]}; }
inline float64x2 from_f64x2(float64x2_t x) { float64x2 r; r._limbs[0] = x.hi; r._limbs[1] = x.lo; return r; }
} // namespace detail

// =============================================================================
// Power, exponential and logarithm
//
// For T == double && N == 2 these call the extern "C" `*dd` functions
// which use Estrin polynomial evaluation and are compiled with full
// optimization.
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> sqrt(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::sqrt(x._limbs[0]);
    return r;
  } else {
    T s = std::sqrt(x._limbs[0]);
    // Bail on 0, -0, negative, NaN, +Inf — the Karp-Markstein refinement
    // would compute `inf - inf = NaN` in the residual step for +Inf, and
    // `0 / 0` for 0. Leading-limb sqrt handles every IEEE special case.
    if (!(x._limbs[0] > T(0)) || !std::isfinite(s)) {
      r._limbs[0] = s;
      return r;
    }
    // Karp/Markstein: r = s + (x - s*s) / (2s), evaluated in DD. The
    // correction reduces the DD residual to a scalar via
    // `residual._limbs[0] * (0.5/s)`, so the residual's lo limb is
    // dropped on the floor. Two higher-fidelity variants were measured
    // (see doc/developer/AUDIT_TODO.md P1):
    //   (a) full DD divide `residual / (2*s_dd)` — sqrt worst case near
    //       perfect squares goes 0.76 → 0.39 ulp, but sqrt bench drops
    //       ~55% and hypot/acosh take a 10–25% hit.
    //   (b) `residual * MultiFloat(0.5/s)` (DD × scalar) — 0.76 → 0.58
    //       ulp, sqrt bench drops ~30%.
    // Baseline is already sub-1-ulp (0 ulp on exact k², ≤0.76 ulp with a
    // non-zero lo limb). The gain from (a)/(b) isn't worth the speed
    // regression for this library's usage pattern; keep baseline.
    const MultiFloat<T, N> s_dd(s);
    const MultiFloat<T, N> residual = x - s_dd * s_dd;
    const MultiFloat<T, N> correction(residual._limbs[0] * (T(0.5) / s));
    return s_dd + correction;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> cbrt(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::cbrt(x._limbs[0]);
    return r;
  } else {
    if (x._limbs[0] == T(0)) {
      return MultiFloat<T, N>();
    }
    const T s = std::cbrt(x._limbs[0]);
    const MultiFloat<T, N> s_dd(s);
    const MultiFloat<T, N> residual = x - s_dd * s_dd * s_dd;
    const MultiFloat<T, N> correction(residual._limbs[0] / (T(3) * s * s));
    return s_dd + correction;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> hypot(MultiFloat<T, N> const &x, MultiFloat<T, N> const &y) {
  if constexpr (N == 1) {
    MultiFloat<T, N> r;
    r._limbs[0] = std::hypot(x._limbs[0], y._limbs[0]);
    return r;
  } else {
    // Defer to libm's hypot for non-finite so inf/NaN propagate correctly.
    if (!std::isfinite(x._limbs[0]) || !std::isfinite(y._limbs[0])) {
      MultiFloat<T, N> r;
      r._limbs[0] = std::hypot(x._limbs[0], y._limbs[0]);
      for (std::size_t i = 1; i < N; ++i) r._limbs[i] = T(0);
      return r;
    }
    MultiFloat<T, N> ax = signbit(x) ? -x : x;
    MultiFloat<T, N> ay = signbit(y) ? -y : y;
    MultiFloat<T, N> big = (ax > ay) ? ax : ay;
    MultiFloat<T, N> small = (ax > ay) ? ay : ax;
    if (big._limbs[0] == T(0)) return MultiFloat<T, N>();
    // Power-of-2 scale (exact) so big has exponent 0 before the square.
    // Replace per-limb ldexp calls with multiplies by 2^(±e): for a power-
    // of-2 multiplier and a non-subnormal result, `x * 2^k` and
    // `ldexp(x, k)` are bit-identical, but the multiply is one FP op while
    // ldexp is a libm call. `e` comes from `ilogb(big.hi)` on a finite
    // non-zero input, so `|e| ≤ 1023`; `2^(-e)` and `2^e` are both finite.
    int e = std::ilogb(big._limbs[0]);
    T down = std::ldexp(T(1), -e);
    for (std::size_t i = 0; i < N; ++i) {
      big._limbs[i] *= down;
      small._limbs[i] *= down;
    }
    MultiFloat<T, N> ratio = small / big;
    MultiFloat<T, N> root = big * sqrt(MultiFloat<T, N>(T(1)) + ratio * ratio);
    MultiFloat<T, N> r;
    T up = std::ldexp(T(1), e);
    for (std::size_t i = 0; i < N; ++i) {
      r._limbs[i] = root._limbs[i] * up;
    }
    // Overflow: true result exceeds T's range. Zero the trailing limbs so
    // callers see a clean inf rather than (inf, NaN).
    if (!std::isfinite(r._limbs[0])) {
      for (std::size_t i = 1; i < N; ++i) r._limbs[i] = T(0);
    }
    return r;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> exp(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::exp(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::expdd(detail::to_f64x2(x)));
  } else {
    T e = std::exp(x._limbs[0]);
    MultiFloat<T, N> e_dd(e);
    return e_dd + e_dd * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> exp2(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::exp2(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::exp2dd(detail::to_f64x2(x)));
  } else {
    T e = std::exp2(x._limbs[0]);
    const T ln2 = std::log(T(2));
    MultiFloat<T, N> e_dd(e);
    return e_dd + e_dd * MultiFloat<T, N>(x._limbs[1] * ln2);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> expm1(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::expm1(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::expm1dd(detail::to_f64x2(x)));
  } else {
    return exp(x) - MultiFloat<T, N>(T(1));
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> log(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::log(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::logdd(detail::to_f64x2(x)));
  } else {
    T l = std::log(x._limbs[0]);
    return MultiFloat<T, N>(l) +
           MultiFloat<T, N>(x._limbs[1] / x._limbs[0]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> log10(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::log10(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::log10dd(detail::to_f64x2(x)));
  } else {
    const T inv_ln10 = T(1) / std::log(T(10));
    return log(x) * MultiFloat<T, N>(inv_ln10);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> log2(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::log2(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::log2dd(detail::to_f64x2(x)));
  } else {
    const T inv_ln2 = T(1) / std::log(T(2));
    return log(x) * MultiFloat<T, N>(inv_ln2);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> log1p(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::log1p(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::log1pdd(detail::to_f64x2(x)));
  } else {
    return log(x + MultiFloat<T, N>(T(1)));
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> pow(MultiFloat<T, N> const &x, MultiFloat<T, N> const &y) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::pow(x._limbs[0], y._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::powdd(detail::to_f64x2(x), detail::to_f64x2(y)));
  } else {
    if (x._limbs[0] == T(0) && y._limbs[0] == T(0)) {
      return MultiFloat<T, N>(T(1));
    }
    return exp(y * log(x));
  }
}

// =============================================================================
// Trigonometric functions
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> sin(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::sin(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::sindd(detail::to_f64x2(x)));
  } else {
    T s = std::sin(x._limbs[0]);
    T c = std::cos(x._limbs[0]);
    return MultiFloat<T, N>(s) +
           MultiFloat<T, N>(c) * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> cos(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::cos(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::cosdd(detail::to_f64x2(x)));
  } else {
    T s = std::sin(x._limbs[0]);
    T c = std::cos(x._limbs[0]);
    return MultiFloat<T, N>(c) -
           MultiFloat<T, N>(s) * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> tan(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::tan(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::tandd(detail::to_f64x2(x)));
  } else {
    return sin(x) / cos(x);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> asin(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::asin(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::asindd(detail::to_f64x2(x)));
  } else {
    T a = std::asin(x._limbs[0]);
    MultiFloat<T, N> denom = sqrt(MultiFloat<T, N>(T(1)) - x * x);
    return MultiFloat<T, N>(a) + MultiFloat<T, N>(x._limbs[1]) / denom;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> acos(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::acos(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::acosdd(detail::to_f64x2(x)));
  } else {
    T a = std::acos(x._limbs[0]);
    MultiFloat<T, N> denom = sqrt(MultiFloat<T, N>(T(1)) - x * x);
    return MultiFloat<T, N>(a) - MultiFloat<T, N>(x._limbs[1]) / denom;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> atan(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::atan(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::atandd(detail::to_f64x2(x)));
  } else {
    T a = std::atan(x._limbs[0]);
    MultiFloat<T, N> denom = MultiFloat<T, N>(T(1)) + x * x;
    return MultiFloat<T, N>(a) + MultiFloat<T, N>(x._limbs[1]) / denom;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> atan2(MultiFloat<T, N> const &y, MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::atan2(y._limbs[0], x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::atan2dd(detail::to_f64x2(y), detail::to_f64x2(x)));
  } else {
    T a = std::atan2(y._limbs[0], x._limbs[0]);
    MultiFloat<T, N> num = x * MultiFloat<T, N>(y._limbs[1]) -
                           y * MultiFloat<T, N>(x._limbs[1]);
    MultiFloat<T, N> denom = x * x + y * y;
    return MultiFloat<T, N>(a) + num / denom;
  }
}

// =============================================================================
// Hyperbolic functions
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> sinh(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::sinh(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::sinhdd(detail::to_f64x2(x)));
  } else {
    T s = std::sinh(x._limbs[0]);
    T c = std::cosh(x._limbs[0]);
    return MultiFloat<T, N>(s) +
           MultiFloat<T, N>(c) * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> cosh(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::cosh(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::coshdd(detail::to_f64x2(x)));
  } else {
    T s = std::sinh(x._limbs[0]);
    T c = std::cosh(x._limbs[0]);
    return MultiFloat<T, N>(c) +
           MultiFloat<T, N>(s) * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> tanh(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::tanh(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::tanhdd(detail::to_f64x2(x)));
  } else {
    return sinh(x) / cosh(x);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> asinh(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::asinh(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::asinhdd(detail::to_f64x2(x)));
  } else {
    // d/dx asinh(x) = 1/sqrt(1 + x^2)
    T a = std::asinh(x._limbs[0]);
    MultiFloat<T, N> denom = sqrt(MultiFloat<T, N>(T(1)) + x * x);
    return MultiFloat<T, N>(a) + MultiFloat<T, N>(x._limbs[1]) / denom;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> acosh(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::acosh(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::acoshdd(detail::to_f64x2(x)));
  } else {
    // d/dx acosh(x) = 1/sqrt(x^2 - 1)
    T a = std::acosh(x._limbs[0]);
    MultiFloat<T, N> denom = sqrt(x * x - MultiFloat<T, N>(T(1)));
    return MultiFloat<T, N>(a) + MultiFloat<T, N>(x._limbs[1]) / denom;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> atanh(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::atanh(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::atanhdd(detail::to_f64x2(x)));
  } else {
    // d/dx atanh(x) = 1/(1 - x^2)
    T a = std::atanh(x._limbs[0]);
    MultiFloat<T, N> denom = MultiFloat<T, N>(T(1)) - x * x;
    return MultiFloat<T, N>(a) + MultiFloat<T, N>(x._limbs[1]) / denom;
  }
}

// =============================================================================
// Error and gamma functions
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> erf(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::erf(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::erfdd(detail::to_f64x2(x)));
  } else {
    // d/dx erf(x) = 2/sqrt(pi) * exp(-x^2)
    const T two_over_sqrt_pi = T(2) / std::sqrt(std::acos(T(-1)));
    T e = std::erf(x._limbs[0]);
    T deriv = two_over_sqrt_pi * std::exp(-x._limbs[0] * x._limbs[0]);
    return MultiFloat<T, N>(e) +
           MultiFloat<T, N>(deriv) * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> erfc(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::erfc(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::from_f64x2(::erfcdd(detail::to_f64x2(x)));
  } else {
    const T two_over_sqrt_pi = T(2) / std::sqrt(std::acos(T(-1)));
    T e = std::erfc(x._limbs[0]);
    T deriv = -two_over_sqrt_pi * std::exp(-x._limbs[0] * x._limbs[0]);
    return MultiFloat<T, N>(e) +
           MultiFloat<T, N>(deriv) * MultiFloat<T, N>(x._limbs[1]);
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> tgamma(MultiFloat<T, N> const &x) {
  if constexpr (N == 2 && std::is_same_v<T, double>) {
    return detail::from_f64x2(::tgammadd(detail::to_f64x2(x)));
  } else {
    MultiFloat<T, N> r;
    r._limbs[0] = std::tgamma(x._limbs[0]);
    return r;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> lgamma(MultiFloat<T, N> const &x) {
  if constexpr (N == 2 && std::is_same_v<T, double>) {
    return detail::from_f64x2(::lgammadd(detail::to_f64x2(x)));
  } else {
    MultiFloat<T, N> r;
    r._limbs[0] = std::lgamma(x._limbs[0]);
    return r;
  }
}

// =============================================================================
// Bessel functions
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> bessel_j0(MultiFloat<T, N> const &x) {
  if constexpr (N == 2 && std::is_same_v<T, double>) {
    return detail::from_f64x2(::j0dd(detail::to_f64x2(x)));
  } else {
    MultiFloat<T, N> r;
    r._limbs[0] = ::j0(x._limbs[0]);
    return r;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> bessel_j1(MultiFloat<T, N> const &x) {
  if constexpr (N == 2 && std::is_same_v<T, double>) {
    return detail::from_f64x2(::j1dd(detail::to_f64x2(x)));
  } else {
    MultiFloat<T, N> r;
    r._limbs[0] = ::j1(x._limbs[0]);
    return r;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> bessel_y0(MultiFloat<T, N> const &x) {
  if constexpr (N == 2 && std::is_same_v<T, double>) {
    return detail::from_f64x2(::y0dd(detail::to_f64x2(x)));
  } else {
    MultiFloat<T, N> r;
    r._limbs[0] = ::y0(x._limbs[0]);
    return r;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> bessel_y1(MultiFloat<T, N> const &x) {
  if constexpr (N == 2 && std::is_same_v<T, double>) {
    return detail::from_f64x2(::y1dd(detail::to_f64x2(x)));
  } else {
    MultiFloat<T, N> r;
    r._limbs[0] = ::y1(x._limbs[0]);
    return r;
  }
}

// =============================================================================
// Additional classification and ordered comparison
// =============================================================================

template <typename T, std::size_t N>
constexpr bool isnormal(MultiFloat<T, N> const &x) {
  return std::isnormal(x._limbs[0]);
}

template <typename T, std::size_t N>
constexpr bool isgreater(MultiFloat<T, N> const &x,
                         MultiFloat<T, N> const &y) {
  return !isnan(x) && !isnan(y) && (x > y);
}

template <typename T, std::size_t N>
constexpr bool isgreaterequal(MultiFloat<T, N> const &x,
                              MultiFloat<T, N> const &y) {
  return !isnan(x) && !isnan(y) && (x >= y);
}

template <typename T, std::size_t N>
constexpr bool isless(MultiFloat<T, N> const &x, MultiFloat<T, N> const &y) {
  return !isnan(x) && !isnan(y) && (x < y);
}

template <typename T, std::size_t N>
constexpr bool islessequal(MultiFloat<T, N> const &x,
                           MultiFloat<T, N> const &y) {
  return !isnan(x) && !isnan(y) && (x <= y);
}

template <typename T, std::size_t N>
constexpr bool islessgreater(MultiFloat<T, N> const &x,
                             MultiFloat<T, N> const &y) {
  return !isnan(x) && !isnan(y) && (x != y);
}

template <typename T, std::size_t N>
constexpr bool isunordered(MultiFloat<T, N> const &x,
                           MultiFloat<T, N> const &y) {
  return isnan(x) || isnan(y);
}

// ---- Formatted I/O for float64x2 -------------------------------------------
//
// Scientific-notation decimal output at up to 34 significant digits (the
// natural DD significand limit is ~32; we allow two guard digits for
// round-half-to-even). Special values use the same textual forms as
// std::to_string for doubles ("nan", "inf", "-inf"). `operator<<` honors
// `os.precision()` when > 17; otherwise (including the C++ default of 6)
// the DD default of 32 digits is used.
//
// These are inline — the library archive exports only extern "C" `*dd`
// symbols, so C++ helpers must live in the header to avoid link errors
// from the symbol-visibility strip in src/localize_symbols.sh.
namespace detail {

// Renormalize (hi, lo) into canonical DD form.
inline void io_renorm(double &hi, double &lo) {
  double s = hi + lo;
  double err = lo - (s - hi);
  hi = s;
  lo = err;
}

// (hi, lo) *= d, with d a double; result renormalized.
inline void io_dd_mul_d(double &hi, double &lo, double d) {
  double p = hi * d;
  double e = std::fma(hi, d, -p);
  lo = lo * d + e;
  hi = p;
  io_renorm(hi, lo);
}

} // namespace detail

inline std::string to_string(float64x2 const &x, int precision = 32) {
  double hi = x._limbs[0];
  double lo = x._limbs[1];

  if (std::isnan(hi) || std::isnan(lo)) return "nan";
  if (std::isinf(hi)) return hi > 0 ? "inf" : "-inf";
  if (hi == 0.0) return std::signbit(hi) ? "-0e+00" : "0e+00";

  if (precision < 1) precision = 1;
  if (precision > 34) precision = 34;

  const bool neg = hi < 0.0;
  if (neg) { hi = -hi; lo = -lo; }

  // Decimal exponent estimate. log10 may be off by 1 due to rounding; a
  // post-scale fixup corrects.
  int e10 = static_cast<int>(std::floor(std::log10(hi)));

  // Scale to [1, 10). Multiplying by 10 is exact in DD; multiplying by 0.1
  // is not, but each iteration retains ~106 bits of precision, so 308
  // iterations (max for dp range) lose at most a handful of ulps — well
  // below the visible 34 digits.
  int shift = -e10;
  while (shift > 0) { detail::io_dd_mul_d(hi, lo, 10.0); --shift; }
  while (shift < 0) { detail::io_dd_mul_d(hi, lo, 0.1);  ++shift; }
  // Drift correction: at most ±1 from the log10 estimate.
  if (hi >= 10.0)     { detail::io_dd_mul_d(hi, lo, 0.1);  ++e10; }
  else if (hi < 1.0)  { detail::io_dd_mul_d(hi, lo, 10.0); --e10; }

  // Extract precision+2 digits; the last two guard round-half-to-even.
  const int ndigits = precision + 2;
  const std::size_t prec = static_cast<std::size_t>(precision);
  char digits[36] = {};
  for (int i = 0; i < ndigits; ++i) {
    int d = static_cast<int>(hi);
    if (d < 0) d = 0;
    else if (d > 9) d = 9;
    digits[i] = static_cast<char>('0' + d);
    hi -= d;
    detail::io_renorm(hi, lo);
    detail::io_dd_mul_d(hi, lo, 10.0);
  }

  int guard = (digits[precision] - '0') * 10 + (digits[precision + 1] - '0');
  bool round_up = guard > 50 ||
                  (guard == 50 && ((digits[precision - 1] - '0') & 1));
  if (round_up) {
    for (int i = precision - 1; i >= 0; --i) {
      if (digits[i] < '9') { ++digits[i]; break; }
      digits[i] = '0';
      if (i == 0) {
        std::memmove(digits + 1, digits, prec - 1);
        digits[0] = '1';
        ++e10;
      }
    }
  }
  digits[precision] = '\0';

  std::string out;
  out.reserve(prec + 8);
  if (neg) out.push_back('-');
  out.push_back(digits[0]);
  if (precision > 1) {
    out.push_back('.');
    out.append(digits + 1, prec - 1);
  }
  out.push_back('e');
  int abs_e = e10 < 0 ? -e10 : e10;
  out.push_back(e10 < 0 ? '-' : '+');
  if (abs_e < 10) out.push_back('0');
  out += std::to_string(abs_e);
  return out;
}

inline std::ostream &operator<<(std::ostream &os, float64x2 const &x) {
  int p = static_cast<int>(os.precision());
  if (p <= 17) p = 32;
  os << to_string(x, p);
  return os;
}

} // namespace multifloats

// ---- std::complex<MultiFloat<double,2>> specializations -------------------
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
// but we leave `std::log(complex<MF>)` etc. on the generic template:
// those paths offer no speedup and no correctness advantage here.
#include <complex>

namespace std {

#define MULTIFLOATS_CX_SPECIALIZE(fn)                                        \
  template <>                                                                \
  inline complex<multifloats::MultiFloat<double, 2>>                         \
  fn(complex<multifloats::MultiFloat<double, 2>> const &z) {                 \
    ::complex64x2_t in = {multifloats::detail::to_f64x2(z.real()),           \
                          multifloats::detail::to_f64x2(z.imag())};          \
    ::complex64x2_t out = ::c##fn##dd(in);                                   \
    return complex<multifloats::MultiFloat<double, 2>>(                      \
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
inline multifloats::MultiFloat<double, 2>
abs(complex<multifloats::MultiFloat<double, 2>> const &z) {
  ::complex64x2_t in = {multifloats::detail::to_f64x2(z.real()),
                        multifloats::detail::to_f64x2(z.imag())};
  return multifloats::detail::from_f64x2(::cabsdd(in));
}

template <>
inline multifloats::MultiFloat<double, 2>
arg(complex<multifloats::MultiFloat<double, 2>> const &z) {
  ::complex64x2_t in = {multifloats::detail::to_f64x2(z.real()),
                        multifloats::detail::to_f64x2(z.imag())};
  return multifloats::detail::from_f64x2(::cargdd(in));
}

template <>
inline complex<multifloats::MultiFloat<double, 2>>
proj(complex<multifloats::MultiFloat<double, 2>> const &z) {
  ::complex64x2_t in = {multifloats::detail::to_f64x2(z.real()),
                        multifloats::detail::to_f64x2(z.imag())};
  ::complex64x2_t out = ::cprojdd(in);
  return complex<multifloats::MultiFloat<double, 2>>(
      multifloats::detail::from_f64x2(out.re),
      multifloats::detail::from_f64x2(out.im));
}

}  // namespace std
