#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

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

  constexpr explicit operator T() const { return _limbs[0]; }

  constexpr bool operator==(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (!(_limbs[i] == rhs._limbs[i])) {
        return false;
      }
    }
    return true;
  }

  constexpr bool operator!=(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (_limbs[i] != rhs._limbs[i]) {
        return true;
      }
    }
    return false;
  }

  constexpr bool operator<(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (_limbs[i] < rhs._limbs[i]) {
        return true;
      } else if (rhs._limbs[i] < _limbs[i]) {
        return false;
      }
    }
    return false;
  }

  constexpr bool operator>(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (_limbs[i] > rhs._limbs[i]) {
        return true;
      } else if (rhs._limbs[i] > _limbs[i]) {
        return false;
      }
    }
    return false;
  }

  constexpr bool operator<=(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (_limbs[i] < rhs._limbs[i]) {
        return true;
      } else if (rhs._limbs[i] < _limbs[i]) {
        return false;
      }
    }
    return true;
  }

  constexpr bool operator>=(MultiFloat const &rhs) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (_limbs[i] > rhs._limbs[i]) {
        return true;
      } else if (rhs._limbs[i] > _limbs[i]) {
        return false;
      }
    }
    return true;
  }

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
    return (*this) + (-rhs);
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
    } else { // N == 2 — single Newton refinement
      T s = _limbs[0] / rhs._limbs[0];
      // Non-finite quotient: y=0, 0/0, inf/finite, NaN.
      if (!std::isfinite(s)) {
        MultiFloat out;
        out._limbs[0] = s;
        return out;
      }
      // Infinite divisor: quotient is ±0 (already in s). The Newton
      // refinement would compute 0*inf = NaN; short-circuit.
      if (!std::isfinite(rhs._limbs[0])) {
        MultiFloat out;
        out._limbs[0] = s;
        return out;
      }
      MultiFloat u;
      u._limbs[0] = T(1) / rhs._limbs[0];
      MultiFloat quotient = (*this) * u;
      MultiFloat residual = quotient * rhs - (*this);
      MultiFloat correction = residual * u;
      return quotient - correction;
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

// =============================================================================
// <cmath>-style free functions (ADL on MultiFloat)
// =============================================================================

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> abs(MultiFloat<T, N> const &x) {
  return std::signbit(x._limbs[0]) ? -x : x;
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
  return std::signbit(x._limbs[0]);
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
  for (std::size_t i = 0; i < N; ++i) {
    r._limbs[i] = std::ldexp(x._limbs[i], n);
  }
  return r;
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> scalbn(MultiFloat<T, N> const &x, int n) {
  return ldexp(x, n);
}

template <typename T, std::size_t N>
constexpr int ilogb(MultiFloat<T, N> const &x) {
  return std::ilogb(x._limbs[0]);
}

template <typename T, std::size_t N>
constexpr MultiFloat<T, N> copysign(MultiFloat<T, N> const &x,
                                    MultiFloat<T, N> const &y) {
  bool xs = std::signbit(x._limbs[0]);
  bool ys = std::signbit(y._limbs[0]);
  return (xs == ys) ? x : -x;
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
  // Round half away from zero, matching std::round.
  MultiFloat<T, N> half(T(0.5));
  return std::signbit(x._limbs[0]) ? -floor(-x + half) : floor(x + half);
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

template <typename T, std::size_t N>
long lround(MultiFloat<T, N> const &x) {
  return std::lround(round(x)._limbs[0]);
}

template <typename T, std::size_t N>
long long llround(MultiFloat<T, N> const &x) {
  return std::llround(round(x)._limbs[0]);
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
    // Approximate one DD ulp ≈ ulp(hi) * 2^-53.
    T inf = std::numeric_limits<T>::infinity();
    T target = (x < y) ? inf : -inf;
    T next_hi = std::nextafter(x._limbs[0], target);
    T ulp = std::abs(next_hi - x._limbs[0]);
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
  // For small quotients (exponent gap <= 53): floor-multiple reduction
  // without DD divide. For large quotients: divide-based path.
  bool x_neg = x._limbs[0] < T(0);
  MultiFloat<T, N> ax = x_neg ? -x : x;
  MultiFloat<T, N> ay = (y._limbs[0] < T(0)) ? -y : y;

  if (ax < ay) return x;

  int diff = std::ilogb(ax._limbs[0]) - std::ilogb(ay._limbs[0]);
  MultiFloat<T, N> r;

  if (diff > 53) {
    r = ax - trunc(ax / ay) * ay;
  } else {
    T q = std::floor(ax._limbs[0] / ay._limbs[0]);
    if (q <= T(1)) {
      r = ax - ay;
    } else {
      r = ax - ay * MultiFloat<T, N>(q);
    }
    if (r._limbs[0] < T(0)) r = r + ay;
    if (!(r < ay)) r = r - ay;
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

// =============================================================================
// detail:: full-DD math kernels for MultiFloat<double, 2>
//
// These kernels mirror the Fortran implementation in fsrc/multifloats.fypp.
// The constant tables are the same hex-literal coefficients ported from
// MultiFloats.jl (see external/MultiFloats.jl/src/{exp,log,trig}.jl).
// =============================================================================

// Forward declaration so detail kernels below can use multifloats::sqrt via ADL.
template <typename T, std::size_t N>
MultiFloat<T, N> sqrt(MultiFloat<T, N> const &x);

namespace detail {

using MFD2 = MultiFloat<double, 2>;

// Polynomial and conversion constants needed by inline helpers below are
// provided through dd_constants.hh which is included only in the .cc file.
// The erf/erfc rational-approximation constants also live in the .cc file.
// ---- helper: build a DD pair from two doubles ------------------------------
inline MFD2 dd_pair(double hi, double lo) {
  MFD2 r;
  r._limbs[0] = hi;
  r._limbs[1] = lo;
  return r;
}

// ---- Horner polynomial evaluation in DD ------------------------------------
inline MFD2 dd_horner(MFD2 const &y, double const *hi, double const *lo,
                       int n) {
  MFD2 p = dd_pair(hi[n - 1], lo[n - 1]);
  for (int i = n - 2; i >= 0; --i) {
    p = p * y + dd_pair(hi[i], lo[i]);
  }
  return p;
}

// ---- Estrin polynomial evaluation (neval/deval from libquadmath) ------------
// neval: evaluates P[n]*x^n + P[n-1]*x^(n-1) + ... + P[0]
// Uses Estrin's scheme for ILP when degree is known, Horner fallback otherwise.
inline MFD2 dd_neval(MFD2 const &x, double const *hi, double const *lo,
                      int n) {
  auto c = [&](int i) { return dd_pair(hi[i], lo[i]); };
  MFD2 x2 = x * x;
  MFD2 x4 = x2 * x2;
  MFD2 x8 = x4 * x4;
  switch (n) {
  case 7: {
    MFD2 p01 = c(0) + c(1) * x;
    MFD2 p23 = c(2) + c(3) * x;
    MFD2 p45 = c(4) + c(5) * x;
    MFD2 p67 = c(6) + c(7) * x;
    MFD2 p03 = p01 + p23 * x2;
    MFD2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4;
  }
  case 8: {
    MFD2 p01 = c(0) + c(1) * x;
    MFD2 p23 = c(2) + c(3) * x;
    MFD2 p45 = c(4) + c(5) * x;
    MFD2 p67 = c(6) + c(7) * x;
    MFD2 p03 = p01 + p23 * x2;
    MFD2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4 + c(8) * x8;
  }
  case 9: {
    MFD2 p01 = c(0) + c(1) * x;
    MFD2 p23 = c(2) + c(3) * x;
    MFD2 p45 = c(4) + c(5) * x;
    MFD2 p67 = c(6) + c(7) * x;
    MFD2 p89 = c(8) + c(9) * x;
    MFD2 p03 = p01 + p23 * x2;
    MFD2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4 + p89 * x8;
  }
  case 10: {
    MFD2 p01 = c(0) + c(1) * x;
    MFD2 p23 = c(2) + c(3) * x;
    MFD2 p45 = c(4) + c(5) * x;
    MFD2 p67 = c(6) + c(7) * x;
    MFD2 p89 = c(8) + c(9) * x;
    MFD2 p03 = p01 + p23 * x2;
    MFD2 p47 = p45 + p67 * x2;
    MFD2 p810 = p89 + c(10) * x2;
    return p03 + p47 * x4 + p810 * x8;
  }
  case 11: {
    MFD2 p01 = c(0) + c(1) * x;
    MFD2 p23 = c(2) + c(3) * x;
    MFD2 p45 = c(4) + c(5) * x;
    MFD2 p67 = c(6) + c(7) * x;
    MFD2 p89 = c(8) + c(9) * x;
    MFD2 p1011 = c(10) + c(11) * x;
    MFD2 p03 = p01 + p23 * x2;
    MFD2 p47 = p45 + p67 * x2;
    MFD2 p811 = p89 + p1011 * x2;
    return p03 + p47 * x4 + p811 * x8;
  }
  default: {
    MFD2 y = c(n);
    for (int i = n - 1; i >= 0; --i)
      y = y * x + c(i);
    return y;
  }
  }
}

// deval: evaluates x^(n+1) + P[n]*x^n + ... + P[0] (monic leading term)
inline MFD2 dd_deval(MFD2 const &x, double const *hi, double const *lo,
                      int n) {
  auto c = [&](int i) { return dd_pair(hi[i], lo[i]); };
  MFD2 x2 = x * x;
  MFD2 x4 = x2 * x2;
  MFD2 x8 = x4 * x4;
  switch (n) {
  case 7: {
    MFD2 p01 = c(0) + c(1) * x;
    MFD2 p23 = c(2) + c(3) * x;
    MFD2 p45 = c(4) + c(5) * x;
    MFD2 p67 = c(6) + c(7) * x;
    MFD2 p03 = p01 + p23 * x2;
    MFD2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4 + x8;
  }
  case 8: {
    MFD2 p01 = c(0) + c(1) * x;
    MFD2 p23 = c(2) + c(3) * x;
    MFD2 p45 = c(4) + c(5) * x;
    MFD2 p67 = c(6) + c(7) * x;
    MFD2 p89 = c(8) + x;
    MFD2 p03 = p01 + p23 * x2;
    MFD2 p47 = p45 + p67 * x2;
    return p03 + p47 * x4 + p89 * x8;
  }
  case 9: {
    MFD2 p01 = c(0) + c(1) * x;
    MFD2 p23 = c(2) + c(3) * x;
    MFD2 p45 = c(4) + c(5) * x;
    MFD2 p67 = c(6) + c(7) * x;
    MFD2 p89 = c(8) + c(9) * x;
    MFD2 p03 = p01 + p23 * x2;
    MFD2 p47 = p45 + p67 * x2;
    MFD2 p810 = p89 + x2;
    return p03 + p47 * x4 + p810 * x8;
  }
  case 10: {
    MFD2 p01 = c(0) + c(1) * x;
    MFD2 p23 = c(2) + c(3) * x;
    MFD2 p45 = c(4) + c(5) * x;
    MFD2 p67 = c(6) + c(7) * x;
    MFD2 p89 = c(8) + c(9) * x;
    MFD2 p1011 = c(10) + x;
    MFD2 p03 = p01 + p23 * x2;
    MFD2 p47 = p45 + p67 * x2;
    MFD2 p811 = p89 + p1011 * x2;
    return p03 + p47 * x4 + p811 * x8;
  }
  case 11: {
    MFD2 p01 = c(0) + c(1) * x;
    MFD2 p23 = c(2) + c(3) * x;
    MFD2 p45 = c(4) + c(5) * x;
    MFD2 p67 = c(6) + c(7) * x;
    MFD2 p89 = c(8) + c(9) * x;
    MFD2 p1011 = c(10) + c(11) * x;
    MFD2 p03 = p01 + p23 * x2;
    MFD2 p47 = p45 + p67 * x2;
    MFD2 p811 = p89 + p1011 * x2;
    return p03 + p47 * x4 + (p811 + x4) * x8;
  }
  default: {
    MFD2 y = x + c(n);
    for (int i = n - 1; i >= 0; --i)
      y = y * x + c(i);
    return y;
  }
  }
}

// ---- exp2 / exp / log2 / log / log10 ---------------------------------------
MFD2 dd_exp2_kernel(MFD2 const &x);
MFD2 dd_exp2_full(MFD2 const &x);
MFD2 dd_exp_full(MFD2 const &x);
MFD2 dd_log2_kernel(MFD2 const &x);
MFD2 dd_log2_full(MFD2 const &x);
MFD2 dd_log_full(MFD2 const &x);
MFD2 dd_log10_full(MFD2 const &x);

// ---- sinpi / cospi / sin / cos / tan ---------------------------------------
MFD2 dd_sinpi_kernel(MFD2 const &x);
MFD2 dd_cospi_kernel(MFD2 const &x);
MFD2 dd_sinpi_full(MFD2 const &x);
MFD2 dd_cospi_full(MFD2 const &x);

// ---- sin / cos / tan -------------------------------------------------------
void dd_reduce_pi_half(MFD2 const &x, MFD2 &r, int &n_mod4);
MFD2 dd_sin_kernel(MFD2 const &x);
MFD2 dd_cos_kernel(MFD2 const &x);
MFD2 dd_sin_eval(MFD2 const &r);
MFD2 dd_cos_eval(MFD2 const &r);
MFD2 dd_sin_full(MFD2 const &x);
MFD2 dd_cos_full(MFD2 const &x);
MFD2 dd_tan_full(MFD2 const &x);

// ---- sinh / cosh / tanh ----------------------------------------------------
MFD2 dd_sinh_full(MFD2 const &x);
MFD2 dd_cosh_full(MFD2 const &x);
MFD2 dd_tanh_full(MFD2 const &x);

// ---- pow -------------------------------------------------------------------
MFD2 dd_pow_full(MFD2 const &x, MFD2 const &y);

// ---- asin / acos / atan / atan2 (Newton on full-DD forward) ---------------
MFD2 dd_asin_full(MFD2 const &x);
MFD2 dd_acos_full(MFD2 const &x);
MFD2 dd_atan_full(MFD2 const &x);
MFD2 dd_atan2_full(MFD2 const &y, MFD2 const &x);

// ---- asinh / acosh / atanh -------------------------------------------------
MFD2 dd_asinh_full(MFD2 const &x);
MFD2 dd_acosh_full(MFD2 const &x);
MFD2 dd_atanh_full(MFD2 const &x);

// ---- erf / erfc (piecewise rational, ported from libquadmath erfq.c) -------
MFD2 dd_erf_full(MFD2 const &x);
MFD2 dd_erfc_full(MFD2 const &x);

// ---- tgamma / lgamma -------------------------------------------------------
MFD2 dd_lgamma_stirling(MFD2 const &x);
MFD2 dd_lgamma_stirling_shift(MFD2 const &x);
MFD2 dd_lgamma_full(MFD2 const &x);
MFD2 dd_tgamma_full(MFD2 const &x);

} // namespace detail

// =============================================================================
// Power, exponential and logarithm
//
// For T == double && N == 2 these use the same polynomial / table-based
// kernels as the Fortran module (ported from MultiFloats.jl). Other (T, N)
// combinations fall back to the leading-limb std:: call (N == 1) or a
// derivative-corrected approach (N == 2 with non-double T).
// against the leading-limb std:: implementation, which yields ~2x the
// precision of double but does not necessarily reach the full ~106-bit
// double-double precision.
// =============================================================================

template <typename T, std::size_t N>
MultiFloat<T, N> sqrt(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::sqrt(x._limbs[0]);
    return r;
  } else {
    if (!(x._limbs[0] > T(0))) {
      r._limbs[0] = std::sqrt(x._limbs[0]); // 0, -0, NaN, or signaling
      return r;
    }
    // Karp/Markstein: r = s + (x - s*s) / (2s), evaluated in DD.
    T s = std::sqrt(x._limbs[0]);
    MultiFloat<T, N> s_dd(s);
    MultiFloat<T, N> residual = x - s_dd * s_dd;
    MultiFloat<T, N> correction(residual._limbs[0] * (T(0.5) / s));
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
    T s = std::cbrt(x._limbs[0]);
    MultiFloat<T, N> s_dd(s);
    MultiFloat<T, N> residual = x - s_dd * s_dd * s_dd;
    MultiFloat<T, N> correction(residual._limbs[0] / (T(3) * s * s));
    return s_dd + correction;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> hypot(MultiFloat<T, N> const &x, MultiFloat<T, N> const &y) {
  return sqrt(x * x + y * y);
}

template <typename T, std::size_t N>
MultiFloat<T, N> exp(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  if constexpr (N == 1) {
    r._limbs[0] = std::exp(x._limbs[0]);
    return r;
  } else if constexpr (std::is_same_v<T, double>) {
    return detail::dd_exp_full(x);
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
    return detail::dd_exp2_full(x);
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
    return detail::dd_log_full(x);
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
    return detail::dd_log10_full(x);
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
    return detail::dd_log2_full(x);
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
    return detail::dd_pow_full(x, y);
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
    return detail::dd_sin_full(x);
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
    return detail::dd_cos_full(x);
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
    return detail::dd_tan_full(x);
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
    return detail::dd_asin_full(x);
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
    return detail::dd_acos_full(x);
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
    return detail::dd_atan_full(x);
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
    return detail::dd_atan2_full(y, x);
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
    return detail::dd_sinh_full(x);
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
    return detail::dd_cosh_full(x);
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
    return detail::dd_tanh_full(x);
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
    return detail::dd_asinh_full(x);
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
    return detail::dd_acosh_full(x);
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
    return detail::dd_atanh_full(x);
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
    return detail::dd_erf_full(x);
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
    return detail::dd_erfc_full(x);
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
    return detail::dd_tgamma_full(x);
  } else {
    MultiFloat<T, N> r;
    r._limbs[0] = std::tgamma(x._limbs[0]);
    return r;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> lgamma(MultiFloat<T, N> const &x) {
  if constexpr (N == 2 && std::is_same_v<T, double>) {
    return detail::dd_lgamma_full(x);
  } else {
    MultiFloat<T, N> r;
    r._limbs[0] = std::lgamma(x._limbs[0]);
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

using float64x2 = MultiFloat<double, 2>;

} // namespace multifloats
