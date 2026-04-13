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

// ---- exp2 polynomial (14 coefs) -------------------------------------------
inline constexpr double exp2_coefs_hi[14] = {
    1.0,                       0.6931471805599453,
    0.2402265069591007,        5.5504108664821576e-02,
    9.618129107628477e-03,     1.3333558146428443e-03,
    1.5403530393381606e-04,    1.5252733804059840e-05,
    1.3215486790144305e-06,    1.0178086009239699e-07,
    7.0549116207971280e-09,    4.4455382718708807e-10,
    2.5678435993488194e-11,    1.3691488882450392e-12};
inline constexpr double exp2_coefs_lo[14] = {
    1.4286718760277223e-34,    2.319046813846301e-17,
    -9.396402954935143e-18,    -1.7580652375606938e-18,
    7.183164273073679e-20,     5.5605293410651576e-21,
    1.5044974304523385e-21,    -8.040605991629947e-22,
    -8.522826592655892e-23,    0.0,
    0.0,                       0.0,
    0.0,                       0.0};
inline constexpr double exp2_min_d = -1022.0;
inline constexpr double exp2_max_d = 1023.9999999999998;

// ---- conversion constants -------------------------------------------------
inline constexpr double log2_e_hi = 1.4426950408889634;
inline constexpr double log2_e_lo = 2.0355273740931033e-17;
inline constexpr double ln_2_hi = 0.6931471805599453;
inline constexpr double ln_2_lo = 2.319046813846301e-17;
inline constexpr double log10_2_hi = 0.3010299956639812;
inline constexpr double log10_2_lo = -2.8037281277851704e-18;

// ---- log2 polynomial: narrow (table path) ---------------------------------
inline constexpr double log2_narrow_hi[7] = {
    2.8853900817779268,    0.96179669392597557,
    0.57707801635558505,   0.41219858311113063,
    0.32192809488736207,   0.26261650020366453,
    0.22240777068641438};
inline constexpr double log2_narrow_lo[7] = {
    4.638031132243375e-17,    5.034626202368297e-17,
    5.288948275806925e-17,    1.0855475812016272e-18,
    0.0,                      0.0,
    0.0};

// ---- log2 polynomial: wide (direct path for x in [15/16, 17/16]) ----------
inline constexpr double log2_wide_hi[9] = {
    2.8853900817779268,     0.96179669392597557,
    0.57707801635558505,    0.41219858311113063,
    0.32192808510751845,    0.26206322127714287,
    0.22146410039507162,    0.19222806650807846,
    0.17031629527116502};
inline constexpr double log2_wide_lo[9] = {
    4.638031132243375e-17,    5.0346244164533024e-17,
    5.272175706110187e-17,    1.3727669404801418e-18,
    -1.0901881027318003e-18,  1.4691558006816344e-18,
    0.0,                      0.0,
    0.0};

// ---- log2 lookup table (32 entries) ---------------------------------------
inline constexpr double log2_centers[32] = {
    1.015625, 1.046875, 1.078125, 1.109375,
    1.140625, 1.171875, 1.203125, 1.234375,
    1.265625, 1.296875, 1.328125, 1.359375,
    1.390625, 1.421875, 1.453125, 1.484375,
    1.515625, 1.546875, 1.578125, 1.609375,
    1.640625, 1.671875, 1.703125, 1.734375,
    1.765625, 1.796875, 1.828125, 1.859375,
    1.890625, 1.921875, 1.953125, 1.984375};
inline constexpr double log2_values_hi[32] = {
    0.02236781302845451,  0.06608919045777244,
    0.10852445677816905,  0.14974711950468206,
    0.18982455888001723,  0.22881869049588088,
    0.2667865406949014,   0.30378074817710293,
    0.33985000288462475,  0.37503943134692475,
    0.4093909361377018,   0.4429434958487283,
    0.47573343096639775,  0.5077946401986962,
    0.5391588111080314,   0.5698556083309478,
    0.5999128421871277,   0.6293566200796096,
    0.6582114827517948,   0.6865005271832184,
    0.7142455176661227,   0.7414669864011469,
    0.7681843247769263,   0.794415866350106,
    0.8201789624151877,   0.8454900509443752,
    0.8703647195834046,   0.8948177633079435,
    0.9188632372745945,   0.9425145053392399,
    0.965784284662087,    0.9886846867721658};
inline constexpr double log2_values_lo[32] = {
    -1.593366605276194e-18,    -4.130247852756734e-18,
     5.4046572138033075e-18,    3.3957331682262494e-18,
    -2.362617117852667e-19,    -5.967894054218645e-18,
    -1.148454798555715e-17,    -8.333787019748188e-18,
    -2.0897960245560436e-17,   1.099000777384843e-17,
    -2.1361956385051908e-17,   2.7429379563921325e-17,
     2.6712179058256416e-18,   2.2368792763711565e-17,
    -4.246405680857825e-17,    3.494516357745965e-18,
    -2.4103897311490816e-17,   4.468163526988311e-17,
    -1.4783628552133162e-17,   -3.0880950164975563e-17,
    -1.670020420476703e-17,    4.3007535189465375e-18,
     2.7943000056050083e-17,   -6.972291600506703e-18,
    -2.0610765990304212e-17,   -1.7498130336849765e-17,
    -3.248732644336383e-17,    2.3416884695657537e-17,
    -7.610716771889941e-19,    1.3760330846314947e-17,
     5.439604524201502e-17,    4.274898271281587e-17};

// ---- 1/pi as DD (for sin/cos via sinpi) -----------------------------------
inline constexpr double inv_pi_hi = 0.3183098861837907;
inline constexpr double inv_pi_lo = -1.9678676675182486e-17;

// ---- sinpi polynomial (15 coefs, c[n] = -(-1)^n * pi^(2n-1)/(2n-1)!) ------
inline constexpr double sinpi_coefs_hi[15] = {
    3.141592653589793,        -5.16771278004997,
    2.5501640398773455,       -0.5992645293207921,
    0.08214588661112823,      -0.0073704309457143504,
    0.00046630280576761255,   -2.1915353447830217e-05,
    7.952054001475513e-07,    -2.2948428997269873e-08,
    5.392664662608129e-10,    -1.0518471716932065e-11,
    1.7302192458361107e-13,   -2.432561179993389e-15,
    2.9567015428549106e-17};
inline constexpr double sinpi_coefs_lo[15] = {
    1.2246467991473532e-16,   2.2665622825789447e-16,
    -7.931006345326556e-17,   2.845026112698218e-17,
    -3.847292805297656e-18,   -3.328281165603432e-19,
    1.0704561733683463e-20,   1.4648526682685598e-21,
    1.736540361519021e-23,    -7.376346207041088e-26,
    -4.6231664587063263e-26,  6.607471301444785e-28,
    4.02155341316903e-30,     1.1975701997015738e-31,
    -2.093244907518996e-34};

// ---- cospi polynomial -----------------------------------------------------
inline constexpr double cospi_coefs_hi[15] = {
    1.0,                      -4.934802200544679,
    4.0587121264167685,       -1.3352627688545895,
    0.2353306303588932,       -0.02580689139001406,
    0.0019295743094039231,    -0.0001046381049248457,
    4.303069587032947e-06,    -1.3878952462213771e-07,
    3.604730797462501e-09,    -7.700707130601354e-11,
    1.3768647280377414e-12,   -2.0906323353147685e-14,
    2.729327261598196e-16};
inline constexpr double cospi_coefs_lo[15] = {
    0.0,                      -3.1326477543698557e-16,
    -2.6602000824298645e-16,  3.1815237892149862e-18,
    -1.2583065576724427e-18,  1.170191067939226e-18,
    -9.669517939986956e-20,   -2.421206183964864e-21,
    -2.864010082936791e-22,   -7.479362090417238e-24,
    -1.833556774402799e-25,   4.7314468253686385e-27,
    -1.6034234137163717e-29,  -4.965817957054884e-32,
    -1.0546803731213643e-32};

// ---- sinh Taylor (9 coefs, n=0..8 → x, x³/6, ..., x¹⁷/17!) ----------------
inline constexpr double sinh_taylor_hi[9] = {
    1.0,                      0.16666666666666666,
    0.008333333333333333,     0.0001984126984126984,
    2.7557319223985893e-06,   2.505210838544172e-08,
    1.6059043836821613e-10,   7.647163731819816e-13,
    2.8114572543455206e-15};
inline constexpr double sinh_taylor_lo[9] = {
    0.0,                      9.25185853854297e-18,
    1.1564823173178714e-19,   1.7209558293420705e-22,
    -1.858393274046472e-22,   -1.448814070935912e-24,
    1.2585294588752098e-26,   7.03872877733453e-30,
    1.6508842730861433e-31};

// ---- asinh Taylor (15 coefs) ----------------------------------------------
inline constexpr double asinh_taylor_hi[15] = {
    1.0,                       -0.16666666666666666,
    0.075,                     -0.044642857142857144,
    0.030381944444444444,      -0.022372159090909092,
    0.017352764423076924,      -0.01396484375,
    0.011551800896139705,      -0.009761609529194078,
    0.008390335809616815,      -0.0073125258735988454,
    0.006447210311889649,      -0.005740037670841924,
    0.005153309682319905};
inline constexpr double asinh_taylor_lo[15] = {
    0.0,                       -9.25185853854297e-18,
    2.7755575615628915e-18,    9.912705577010326e-19,
    3.854941057726238e-19,     9.462128050782583e-19,
    -8.006416042969879e-19,    6.938893903907229e-19,
    8.163404592832033e-19,     -5.478074134663601e-19,
    4.130293990420969e-19,     3.394024192128536e-19,
    -3.1225022567582527e-19,   1.2849803525754126e-19,
    -3.888173308223878e-19};

// ---- atanh Taylor (15 coefs, c[n] = 1/(2n+1)) -----------------------------
inline constexpr double atanh_taylor_hi[15] = {
    1.0,                       0.3333333333333333,
    0.2,                       0.14285714285714285,
    0.1111111111111111,        0.09090909090909091,
    0.07692307692307693,       0.06666666666666667,
    0.058823529411764705,      0.05263157894736842,
    0.047619047619047616,      0.043478260869565216,
    0.04,                      0.037037037037037035,
    0.034482758620689655};
inline constexpr double atanh_taylor_lo[15] = {
    0.0,                       1.850371707708594e-17,
    -1.1102230246251566e-17,   7.93016446160826e-18,
    6.1679056923619804e-18,    -2.523234146875356e-18,
    -4.270088556250602e-18,    9.251858538542971e-19,
    8.163404592832033e-19,     2.921639538487254e-18,
    2.64338815386942e-18,      1.206764157201257e-18,
    -8.326672684688674e-19,    2.05596856412066e-18,
    4.785444071660157e-19};

// ---- erf/erfc rational approximation constants (from libquadmath erfq.c) ----
// Each constant is stored as a DD pair (hi, lo) in separate arrays.
// 2/sqrt(pi) - 1
inline constexpr double erf_efx_hi = 1.28379167095512586316e-01;
inline constexpr double erf_efx_lo = -1.24201160024630322884e-17;
// erf(1) split for near-1 approximation
inline constexpr double erf_const_hi = 8.45062911510467529297e-01;
inline constexpr double erf_const_lo = 0.0;
// erf(x) = x + x*P(x^2)/Q(x^2), |x| < 7/8
inline constexpr double erf_TN1_hi[9] = {
    -3.85825232425463714600e+10, 9.58031924859046478271e+10,
     1.30217051973488006592e+10, 2.92295695042639732361e+09,
     1.76431752078331947327e+08, 1.57343601460111867636e+07,
     4.02807738010572153144e+05, 1.64405680646728906140e+04,
     3.39086848005999144107e+01};
inline constexpr double erf_TN1_lo[9] = {
     2.14529208425578934281e-07,-1.00398476873047788400e-06,
    -8.83220675133069620841e-07, 9.41919230492270262881e-08,
    -7.53976783441145342803e-09,-4.62586250496415896819e-10,
    -1.42690214987251856082e-11, 5.45016441182450261786e-14,
     1.99168862057096300681e-15};
inline constexpr double erf_TD1_hi[9] = {
    -3.00535703069653320312e+11,-1.34260228312628280640e+11,
    -2.77715389335534095764e+10,-3.48382639103353214264e+09,
    -2.90632104707129955292e+08,-1.65334798572215419263e+07,
    -6.24552058156284852885e+05,-1.40212430417749874323e+04,
    -1.20936807247351069350e+02};
inline constexpr double erf_TD1_lo[9] = {
     2.75975114469310474536e-05,-2.10151742988399802589e-06,
    -3.64690946194366872661e-08, 1.45683540082177225790e-07,
    -3.27652420425108948098e-09, 3.01870884093293701209e-10,
    -2.49619401951340283333e-11,-8.53611593818100578213e-13,
     1.90054806458197120934e-15};
// erf(z+1) = erf_const + P(z)/Q(z), -0.125 <= z <= 0
inline constexpr double erf_TN2_hi[9] = {
    -4.08888969707748515248e+01, 7.15704643068180848786e+03,
    -2.19156191257440968911e+03, 2.18017491655531694050e+03,
     2.84857865804967048007e+02, 1.63036249095251292829e+02,
     6.31771235396186714439e+00, 2.45044103418349230594e+01,
     5.12766227770678817421e-01};
inline constexpr double erf_TN2_lo[9] = {
    -1.48533903324322527156e-15, 6.59808902688830162694e-14,
    -1.76436052090094233718e-13,-6.55139809198715350626e-14,
     1.88163317350736866056e-14,-9.15255322018434752951e-15,
    -1.70248078022225950350e-16, 1.28713277465725528564e-15,
    -5.09258099773952635238e-17};
inline constexpr double erf_TD2_hi[9] = {
     1.73102644592683391238e+04, 1.20968223900799039257e+04,
     1.16095029021799364273e+04, 5.39429464512712638680e+03,
     2.79123934053363245766e+03, 8.98936557133731866998e+02,
     2.97401649376634964028e+02, 6.14819275459037655196e+01,
     1.17850289249073849618e+01};
inline constexpr double erf_TD2_lo[9] = {
     9.58927979095155003479e-13,-2.17691489683133634456e-13,
    -1.40677647334490043822e-14, 1.91028917228274996716e-13,
     2.11778387732466363616e-13, 3.62967502970887454346e-14,
    -2.30556601665613793573e-14,-1.73219446212744015948e-15,
    -5.05229684170886269105e-16};
// erfc sub-interval constants: erfc(x0) = Ca + Cb
inline constexpr double erfc_Ca_hi[8] = {
    7.23663330078125000000e-01, 5.95870971679687500000e-01,
    4.79492187500000000000e-01, 3.76754760742187500000e-01,
    2.88833618164062500000e-01, 2.15911865234375000000e-01,
    1.57287597656250000000e-01, 1.11602783203125000000e-01};
inline constexpr double erfc_Ca_lo[8] = {0,0,0,0,0,0,0,0};
inline constexpr double erfc_Cb_hi[8] = {
    1.02797536380670146919e-05, 1.21188854902016758042e-05,
    7.93468695346231738495e-06, 4.35706939452755150554e-06,
    1.07481824223684007167e-05, 1.30737057653416849650e-05,
    1.16093940351306588923e-05, 8.98509516723592969827e-06};
inline constexpr double erfc_Cb_lo[8] = {
    2.39863816665404922098e-22, 3.70686184884307153505e-22,
   -1.31606354404668428198e-22,-1.46049487867006429167e-22,
    3.45472400088343275166e-22, 4.99328378791865841416e-22,
   -1.12906487312226461537e-22, 7.23284853368862606569e-22};
// erfc sub-interval center points
inline constexpr double erfc_x0[8] = {
    0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125};
// erfc region 2 numerator coefficients (8 sub-intervals, 9 coefficients each)
inline constexpr double erfc_RN_hi[8][9] = {
  {-2.35370709764128059760e+03, 3.87115965622874341534e+02,
   -3.88810513425826627554e+02,-2.12999853912006180678e+01,
   -8.12546226359403505057e+01, 8.15154909398350646654e+00,
   -5.03336203272920723606e+00,-4.25395662113513609026e-02,
   -8.09860287846385501487e-02},
  {-2.44616401640442609278e+03, 6.71875332449656411882e+02,
   -4.58163113804983595401e+02,-2.38284408898709223479e+01,
   -7.11923785240060027490e+01, 1.31360964610842021472e+01,
   -6.18860870208226465650e+00,-2.78711660110667822088e-02,
   -2.23039557057415384345e-02},
  {-2.62421241801118139847e+03, 8.47382890464782462914e+02,
   -5.28620745862838020912e+02,-3.89578123415531578644e+01,
   -6.20085790806516357065e+01, 1.46932461034692405377e+01,
   -6.96135652537065841017e+00, 5.14572438664116348778e-03,
    1.99025365594817961468e-02},
  {-2.34788794320068063826e+03, 8.00859066069210484784e+02,
   -5.25736331038412004091e+02,-4.47173771785780118648e+01,
   -4.84954038645257341500e+01, 1.14088526467713471391e+01,
   -6.73159108546026985209e+00, 1.37083165303304749250e-01,
    2.02295827998213890997e-02},
  {-1.76706873422027774723e+03, 6.69374664566524302245e+02,
   -4.74622424183727616764e+02,-2.27416063772878267457e+01,
   -3.54123226614093908893e+01, 6.98895051474705297778e+00,
   -5.80768721683654121080e+00, 3.63191598856734632061e-01,
   -1.48894548714963489283e-02},
  {-1.34204489908759342143e+03, 6.12722129422917305419e+02,
   -4.51982135652229146672e+02, 1.22327517782512877886e+01,
   -2.73078957138297120366e+01, 4.04518120492153876455e+00,
   -4.92514647787659232137e+00, 5.93387803661127977151e-01,
   -5.55764543585891626631e-02},
  {-1.13918093645415729043e+03, 6.13490312908689929827e+02,
   -4.62890902471532967866e+02, 4.16570238721073238253e+01,
   -2.28697991351522986747e+01, 1.87069525644987266766e+00,
   -4.17748660127310600387e+00, 7.53398037278964594066e-01,
   -8.62994543691775195526e-02},
  {-9.65270691645797342062e+02, 5.57706639605093300815e+02,
   -4.40633550884849682916e+02, 5.20289346649024295743e+01,
   -1.93131184766575785261e+01,-9.36431826874828815432e-02,
   -3.30639035128635283556e+00, 7.57380604528904433081e-01,
   -9.61174401148909285375e-02}};
inline constexpr double erfc_RN_lo[8][9] = {
  { 4.73170558809419951016e-14, 1.84653684177192327425e-14,
    8.33271221082055739908e-15, 1.38745372074413585866e-15,
    3.78104992321675881272e-15,-6.56424767749351984789e-16,
   -7.44032122400990630036e-17,-3.43964594862527500239e-22,
    2.25087681188738012879e-18},
  {-1.84800580322042305787e-13,-2.05425115189087899227e-14,
   -2.03418397743828986339e-14, 1.75496335496371180531e-17,
   -2.33022320291236327506e-15,-7.83838431159935867964e-16,
    2.67340728270014415322e-16,-6.63991075049570428079e-19,
   -1.19757229180367141697e-18},
  {-8.94553665141629588152e-14, 5.51936655467450532078e-14,
   -5.55979182060938498379e-14, 5.73498155990782167498e-16,
   -4.73867642051884384913e-16,-5.23788618568000253440e-16,
   -1.62632840435526800850e-16, 3.21818850554680407555e-19,
    9.87395168290367380940e-19},
  { 7.44713270246896720630e-14, 1.56940922085480358138e-14,
    3.12147692519952421399e-14,-4.39673879655837731306e-16,
    1.08290483286744236000e-15,-3.46305387930428683249e-16,
    4.04165349821486612008e-16,-5.21508605206189154568e-18,
   -1.54948432405526917843e-18},
  { 1.89932817942957236510e-14,-1.90022539120785261827e-14,
    2.09517661839342883154e-14,-5.74217768067705129051e-18,
    3.88331822499492179129e-16,-3.01387513182091940987e-16,
    3.79914633151257874730e-16, 1.18215679620870342081e-17,
    7.22943473328149096478e-19},
  { 2.40076774280190757211e-14,-5.66817987200589791875e-15,
    2.81096227299575779225e-14,-4.63659624252108973450e-16,
   -1.51963693402034040141e-15, 1.22328402309579182903e-16,
   -4.02035905754187183637e-16,-5.26851833405138632648e-17,
    2.40861448143013714199e-18},
  { 9.69313641057598085655e-14, 4.39242384025932109038e-14,
    1.16330861780438708458e-14,-2.99626347494855114681e-16,
    1.20267033308939329945e-15, 7.54078319160422241181e-17,
    2.50993312992860184735e-16, 1.99455869344006580758e-17,
   -4.77944534237082532037e-19},
  {-5.35746756780003163605e-14,-2.31471500641976833321e-14,
    1.15587334720672882398e-14,-2.23863787673675539397e-15,
   -6.07123826166575093736e-16, 4.90056602897308570878e-18,
    7.06667357252376399693e-17, 3.16912681475967571135e-17,
   -4.02769301156812748594e-19}};
// erfc region 2 denominator coefficients (8 sub-intervals, 8 coefficients each)
inline constexpr double erfc_RD_hi[8][8] = {
  { 2.22044879630669356629e+03, 1.89913325877957873900e+02,
    1.06190671228496103140e+03, 7.49708607230696770785e+01,
    2.14679611566267283251e+02, 1.12015600836257380735e+01,
    2.21101495207505251983e+01, 6.46965567532615026813e-01},
  { 2.49518743924186992444e+03, 2.50354944987292554970e+02,
    1.15903356098889548775e+03, 9.49375146654230519516e+01,
    2.27621492956235442762e+02, 1.36769752121906922326e+01,
    2.27698839599552833590e+01, 7.64774575364899678043e-01},
  { 2.98619076084797507065e+03, 5.28826275896107290464e+02,
    1.36364917807100687241e+03, 1.92170797564991602258e+02,
    2.58865110065102896897e+02, 2.62875292032145573273e+01,
    2.45564903588511427301e+01, 1.37882665359512857073e+00},
  { 3.07516617002483735632e+03, 8.73046894216079749640e+02,
    1.45847279916634056462e+03, 3.23042368756801977270e+02,
    2.80400987271989379224e+02, 4.46533422132322286302e+01,
    2.61272325968320586753e+01, 2.34152675118524422615e+00},
  { 2.74845752349815074922e+03, 1.02021339071347767913e+03,
    1.38885763593543265415e+03, 3.90336368114381798478e+02,
    2.78456834437813938621e+02, 5.55580083021676500721e+01,
    2.64621547095905036429e+01, 2.98490528210351735439e+00},
  { 2.55751800066170062564e+03, 1.07017143338288906307e+03,
    1.34484283442349305915e+03, 4.16114447844938183607e+02,
    2.76367025221985500139e+02, 5.99815348786894375621e+01,
    2.65769510843862875049e+01, 3.25214052439442191655e+00},
  { 2.74430344798113264915e+03, 1.26639635952618709780e+03,
    1.46673946142207341836e+03, 4.86871057075969417838e+02,
    2.99369430155975578600e+02, 6.86897681951025447233e+01,
    2.80150581624767731626e+01, 3.60443990919435020004e+00},
  { 3.03282962952014258917e+03, 1.65964847072196766931e+03,
    1.70354512865728452198e+03, 6.39346567773159904391e+02,
    3.48913139728103089965e+02, 8.84864173857078384344e+01,
    3.13226906255239292420e+01, 4.43013166329056318204e+00}};
inline constexpr double erfc_RD_lo[8][8] = {
  {-6.27441830403601723198e-14,-5.02040003320070854315e-15,
    7.87958289859515256886e-14, 2.57335362878924448224e-15,
   -3.66371719442709530328e-15,-7.06821692215733683415e-16,
    9.65789117453288412242e-16, 5.17558879533063371669e-17},
  {-1.91748249067394415999e-13, 3.03130186678051175078e-15,
   -6.05527361070248105269e-15,-7.03900555675233299669e-15,
   -9.93553881630375114594e-15, 5.70950974115243053483e-16,
    1.59153358834625447396e-15,-2.20589590972103180966e-17},
  {-1.27612476587595863438e-13, 1.61695532839709463566e-14,
    1.05940385331141582122e-13,-1.28333915527548085622e-14,
    5.41004105915573310737e-15,-1.26173887037390337903e-15,
    3.59637825960563935540e-16,-1.06345306333311082039e-16},
  {-1.40923354541151471183e-13, 5.35209025831859823487e-14,
   -8.48824065248633615298e-14,-6.32447496126085442220e-15,
   -1.80155579860175792787e-14, 8.03976668033019813376e-16,
    6.05719911577261535715e-16,-1.16427757814018246642e-16},
  {-7.25488371190404850714e-15, 7.64791900673302881914e-15,
   -3.21780435532953389353e-14,-2.33879524075479233945e-14,
    1.13009724882659258101e-14,-3.04428630600687946536e-15,
   -8.48553111096568905751e-16, 1.42695553857687075127e-16},
  {-3.68834648476162097763e-14,-6.81138302214246143968e-14,
    2.19061761356153048257e-14, 6.51371083551825690335e-15,
    1.96666938020146042168e-14,-4.79753188446535142992e-16,
    9.72417888775549545940e-16,-4.76308762493670926029e-17},
  { 5.22855317127000834380e-14,-3.25767377066117987484e-14,
   -6.68638610635560342950e-14,-2.22780189082456294172e-14,
    2.60476274521993943569e-14,-3.32592495196005085218e-15,
   -1.22775058540920737047e-15, 6.35137249807378684914e-17},
  {-2.50652335516190697564e-14, 5.06554262280507628754e-14,
    9.74245121093219960691e-14,-1.71409027222368703416e-14,
    4.77566742819381640423e-15,-4.36958420122609606540e-15,
    5.06296566118261205205e-16, 3.41888485601218206778e-16}};
// erfc asymptotic region: erfc(1/x) = (1/x)*exp(-1/x^2-0.5625+R(1/x^2))
// 8 sub-intervals in 1/x space, numerator coefficients
// Degrees: RNr1=9, RNr2=11, RNr3=11, RNr4=10, RNr5=10, RNr6=9, RNr7=9, RNr8=9
inline constexpr int erfc_asym_ndeg[8] = {9, 11, 11, 10, 10, 9, 9, 9};
inline constexpr int erfc_asym_ddeg[8] = {8, 10, 10, 10, 9, 9, 9, 8};
inline constexpr double erfc_asym_N_hi[8][12] = {
  {-4.25078088320236202508e-08,-5.37577705328861200733e-06,
   -2.57364594922089671026e-04,-6.19903292811354190289e-03,
   -8.26272119869340404552e-02,-6.24261522725732431738e-01,
   -2.60987473919959533930e+00,-5.58196756333667654104e+00,
   -5.12439892335602298346e+00,-1.29086524394429247309e+00,
    0,0},
  {-2.63891438342028712218e-08,-3.47919837026063403899e-06,
   -1.78398529533569765977e-04,-4.77787693312257608080e-03,
   -7.45063473898732453460e-02,-7.06831885487473332574e-01,
   -4.11391992193594457916e+00,-1.44044757322690628598e+01,
   -2.88348403153071828342e+01,-2.99088697432847645530e+01,
   -1.32528391491510486588e+01,-1.57243610622807028498e+00},
  {-1.95240112655120223790e-06,-2.13088174306637296923e-04,
   -8.37649395809019099712e-03,-1.65059264656098769741e-01,
   -1.83929081893331725084e+00,-1.21627871557088234056e+01,
   -4.81875934446236016129e+01,-1.12099466129747682430e+02,
   -1.45285076566231936113e+02,-9.48520785112895765678e+01,
   -2.56366385502579667843e+01,-1.78799594418756568892e+00},
  { 3.25853071202452781313e-03, 2.98705601687727806404e-01,
    8.73872908934020031779e+00, 1.20721116014864776389e+02,
    8.99755863248903324347e+02, 3.79802519769975742747e+03,
    9.11320366868308155972e+03, 1.20328589133993318683e+04,
    8.10064705791914002475e+03, 2.38388824990714510932e+03,
    2.12749357316645415494e+02, 0},
  {-3.33225892745528558259e-03,-2.69710075890028044832e-01,
   -6.08332855113962178706e+00,-6.11986352898330778771e+01,
   -3.17653528247559336251e+02,-8.93339517508056133011e+02,
   -1.36001950848847604902e+03,-1.07507557982818866549e+03,
   -4.01734656158601467268e+02,-5.85758136814526650937e+01,
   -2.07771592558783479987e+00, 0},
  { 1.64207687617683428805e+00, 1.20715000361111762572e+02,
    2.11926077931638974405e+03, 1.56294222773466335639e+04,
    5.65677918954971028143e+04, 1.05216624102148169186e+05,
    9.94979852478600078030e+04, 4.49179073408026524703e+04,
    8.37707409830153119401e+03, 4.50693480656798669770e+02,
    0, 0},
  { 1.68622219338598782201e+01, 1.17822454356760431438e+03,
    1.76455058429014934518e+04, 1.07375832189033477334e+05,
    3.13284074920594284777e+05, 4.60786493997409997974e+05,
    3.38978182010585209355e+05, 1.17404218711056513712e+05,
    1.66001360601116721227e+04, 6.70039395748066226588e+02,
    0, 0},
  { 3.58745148925535630724e+01, 5.40624974908734088785e+02,
    2.93130129062525111294e+03, 7.35925418524179531232e+03,
    9.20103184981063532177e+03, 5.74969709619319110061e+03,
    1.71041523441986078069e+03, 2.15075398254337869730e+02,
    8.74095358227214802582e+00, 4.87642297882871708636e-02,
    0, 0}};
inline constexpr double erfc_asym_N_lo[8][12] = {
  { 7.83791091856538706709e-25,-2.75159212843094247548e-22,
   -1.05949400481131210179e-20,-1.77370627751118470004e-19,
   -1.48590547647730991844e-19,-4.28993178067009580044e-17,
   -6.09244275626207188802e-17,-1.96104888489146111057e-16,
    3.73747713710190752080e-16, 1.02426809259807919963e-16,
    0,0},
  {-9.02204283809368586671e-25, 6.17337593656670016981e-23,
   -2.66074943137398301448e-21, 6.65371626226129166949e-20,
   -4.69472616791244541300e-18, 9.77184486105112586122e-19,
   -2.16603012006477360921e-16, 6.35636184521030092233e-16,
   -1.44996812168027541927e-15, 6.80203895586533049418e-16,
   -3.64602666517138111696e-18, 8.94744799591240228102e-17},
  { 2.92051220505938053031e-23, 1.67122052302002748923e-21,
    5.33858070720538207152e-20,-3.24882877426707872394e-19,
   -8.72729651085656826743e-17,-8.18528964875564981580e-16,
   -2.66324627288736927754e-15,-5.25023390952317401744e-15,
    9.69355786492914300685e-15, 5.48135275484624244360e-15,
    3.72167336972575267527e-16, 1.20830791212466141026e-17},
  { 2.19583494418401645087e-20,-1.34316484243471888421e-17,
   -5.67058141620254752945e-16, 1.85062981510161977197e-15,
   -3.41222862446749120244e-14,-2.01493512480007276643e-13,
   -5.84083406031399366527e-13, 5.17809905806477547606e-13,
    3.03785784513593951060e-13,-1.63480685534993422465e-13,
    9.42784892234639236005e-15, 0},
  { 1.24230613261092849426e-19, 4.56562249586093656816e-18,
    2.65641712967523494516e-16,-2.25258856306647001976e-15,
    1.89256795231336864350e-14, 4.04295418421603096804e-14,
    7.09631121572727401257e-14, 4.39462883421493879415e-14,
   -1.50139917998920314773e-14, 2.59863369536613786717e-15,
    1.93491159726021609125e-16, 0},
  { 1.02578017826786524228e-16, 6.32770277928578835669e-15,
    1.60689279048415465599e-13, 8.54148621818906978984e-13,
   -2.01438002940682022278e-12, 6.31692135445462718195e-15,
   -1.84673977835465696121e-12,-2.03623148452716109232e-12,
   -8.67740317399831899031e-13, 1.12391135245200087839e-14,
    0, 0},
  {-1.31225938630607544828e-15,-9.87733666081255936417e-14,
    1.21474811472372802790e-12, 4.86597929295448411595e-12,
    2.89853856963419782094e-11, 2.44879616115307300304e-11,
    2.09577688211568203204e-11, 6.57569102385174447691e-12,
   -6.82263893276342960730e-13,-3.28184226814005302671e-14,
    0, 0},
  {-5.64798739626397898355e-16,-4.55982846992354350432e-14,
   -2.26703528741816804503e-13, 2.71790836728392833168e-13,
    7.82338513289933378601e-13, 3.67141230258771804118e-13,
    4.50205662059899526887e-14,-1.16441248759464763781e-14,
   -6.90721851322408923462e-16, 1.33266161528873235095e-18,
    0, 0}};
// erfc asymptotic region denominator coefficients
inline constexpr double erfc_asym_D_hi[8][11] = {
  { 4.30897666174950863760e-06, 3.26539012643278011173e-04,
    9.81132883918704040704e-03, 1.51122251503602089695e-01,
    1.28926434191742989022e+00, 6.14764035618223037005e+00,
    1.57396687133773980349e+01, 1.95553412343509513960e+01,
    9.47261312136313549104e+00, 0, 0},
  { 2.67504272813673208157e-06, 2.17099786845181271866e-04,
    7.24996975268754010463e-03, 1.30204037585976872826e-01,
    1.38020248308291093586e+00, 8.92659411317416484621e+00,
    3.52108958478261655500e+01, 8.23354742753318191717e+01,
    1.07297157988580309507e+02, 6.94380311333796385043e+01,
    1.77569534103160791005e+01},
  { 1.97913068677034952557e-04, 1.15694171612848832609e-02,
    2.75265763430988630311e-01, 3.48224545724831857640e+00,
    2.56934706937269652371e+01, 1.14227900018045744446e+02,
    3.05650397719056456936e+02, 4.78084402092379491478e+02,
    4.10597272721255421857e+02, 1.72407218806374686437e+02,
    2.81593918346481828507e+01},
  {-3.30314198151454019303e-01,-1.35376862936360531364e+01,
   -2.20612763030362145855e+02,-1.86180033875806680044e+03,
   -8.88904877587260489236e+03,-2.46588810662794821837e+04,
   -3.93464221171077442705e+04,-3.45507725824225271936e+04,
   -1.52408397743969035218e+04,-2.81054188739798473762e+03,
   -1.34392955354115997579e+02},
  { 3.37787957041739916875e-01, 1.02196332274239072291e+01,
    1.20084764659294208400e+02, 7.11891552814292708717e+02,
    2.31815938006206624777e+03, 4.23872985353400963504e+03,
    4.27911490728482567647e+03, 2.25727718666326154562e+03,
    5.57047550128505463363e+02, 5.14218924385628923801e+01,
    0},
  {-1.66455764392826296216e+02,-3.80003590250765682867e+03,
   -3.27702819159173523076e+04,-1.38135947150288557168e+05,
   -3.08220428738258196972e+05,-3.69107148825673852116e+05,
   -2.30048244303834973834e+05,-6.87395530092763656285e+04,
   -8.26215881797833390010e+03,-2.51712225438443084613e+02,
    0},
  {-1.70930502471835893630e+03,-3.28003388748133293120e+04,
   -2.34528422802252200199e+05,-8.08675812309776432812e+05,
   -1.45690041451010876335e+06,-1.39165426488125510514e+06,
   -6.84236080186994047835e+05,-1.59743021444657351822e+05,
   -1.48887613060987659992e+04,-3.51176295093506041667e+02,
    0},
  { 6.35859313409690827257e+01, 9.90025381655245041657e+02,
    5.64292877785680138913e+03, 1.52419537519957084442e+04,
    2.11382964450000690704e+04, 1.52643856262646568211e+04,
    5.56137092214924177824e+03, 9.39403553017970466499e+02,
    6.14701959615039470464e+01, 0, 0}};
inline constexpr double erfc_asym_D_lo[8][11] = {
  { 3.97243123796133590860e-22, 7.23978614926339926734e-21,
    2.94858028245824543725e-19, 1.36459859037768887366e-17,
    6.86420395588609769671e-17, 3.99495641339572916525e-16,
   -1.89753586783063509343e-16,-7.23985869948104006900e-16,
   -1.87896812265910654812e-17, 0, 0},
  {-1.58020823740362359504e-22,-1.00778141680144524005e-21,
    1.84792891859379057074e-19,-5.36412456154939475181e-18,
   -4.69655422180910617019e-17, 5.06409949598045577615e-16,
   -8.26302782909874700806e-16,-5.41984454364874684207e-15,
   -6.19874066280096986385e-15, 6.19309418054138089925e-15,
   -1.71820327176194331122e-15},
  {-4.41083641755413185010e-21,-5.98493365038280273320e-19,
    3.33167989860232786036e-18, 2.10945533741279878963e-16,
   -1.65136294673291711692e-15,-2.47174009154103326126e-15,
   -2.75013624219138905579e-14,-9.31227787315541237817e-15,
    5.89285818227451989054e-15, 1.06492367178252564843e-14,
   -8.63689378460875922809e-16},
  {-8.11331474370500391680e-18, 1.29344886083329793567e-16,
   -6.33984055469487332284e-15, 1.03923038009302622689e-13,
   -8.15892818800572200478e-13, 7.89582092139533521049e-14,
   -6.78320427814793867282e-13,-2.55581733520039143535e-12,
    6.73641119219191933911e-13,-6.66198901608592701655e-14,
    4.19686602764693345181e-15},
  { 1.72803813491954528427e-17, 1.25218964946173042930e-16,
    1.11943430690240053381e-15, 1.69041528425541019063e-15,
    2.21612489324792338044e-13,-4.14015420175223832532e-13,
    2.09796430368067300454e-13,-1.45650043096089614785e-14,
   -3.40254805965876976888e-14,-2.56864717557594868212e-15,
    0},
  {-1.29716311879478228287e-14, 2.04084422308599885837e-13,
    3.02398735516117768713e-12, 1.25280373929831598611e-11,
    9.61905764413192768963e-12, 1.78150423772712648395e-11,
   -7.74061640863660987171e-12, 3.26153797429813149157e-12,
   -2.41981029997130034372e-13,-1.34991518130922011641e-15,
    0},
  { 6.15960150306902377950e-14,-2.68378559604889759029e-12,
    1.16899838736968820091e-11, 3.56195427247808446605e-11,
    4.49496124765592618778e-11, 3.67454359749143305480e-11,
    4.94680164990596574870e-11,-4.79630654403658214658e-12,
   -8.15021287216339079090e-13, 1.15264821214079067876e-14,
    0},
  { 7.83616802244381145894e-16,-3.42808814820722420962e-14,
   -3.68589243562754122898e-13, 2.37795774327328729815e-13,
   -1.57090126378563680620e-12, 2.41530036006392171817e-13,
   -3.21111445525845256551e-13, 3.86620976418577836784e-14,
   -1.26652535042228389060e-15, 0, 0}};

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
inline MFD2 dd_exp2_kernel(MFD2 const &x) {
  // n = nearest integer to x.hi, y = (x - n)/8
  double n_float = std::nearbyint(x._limbs[0]);
  MFD2 y = x + dd_pair(-n_float, 0.0);
  y._limbs[0] = std::ldexp(y._limbs[0], -3);
  y._limbs[1] = std::ldexp(y._limbs[1], -3);
  // poly(y) then cube via three squarings (to undo the /8)
  MFD2 p = dd_horner(y, exp2_coefs_hi, exp2_coefs_lo, 14);
  p = p * p;
  p = p * p;
  p = p * p;
  // multiply by 2^n via two ldexps (split to keep intermediates in range)
  int n = static_cast<int>(n_float);
  int half_n = n / 2;
  MFD2 r;
  r._limbs[0] = std::ldexp(p._limbs[0], half_n);
  r._limbs[1] = std::ldexp(p._limbs[1], half_n);
  r._limbs[0] = std::ldexp(r._limbs[0], n - half_n);
  r._limbs[1] = std::ldexp(r._limbs[1], n - half_n);
  return r;
}

inline MFD2 dd_exp2_full(MFD2 const &x) {
  MFD2 r;
  if (!std::isfinite(x._limbs[0])) {
    if (x._limbs[0] > 0.0) r._limbs[0] = x._limbs[0]; // +inf
    else if (x._limbs[0] < 0.0) r._limbs[0] = 0.0;    // -inf → 0
    else r._limbs[0] = x._limbs[0];                    // NaN
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] < exp2_min_d) return MFD2();
  if (x._limbs[0] > exp2_max_d) {
    r._limbs[0] = std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
    return r;
  }
  return dd_exp2_kernel(x);
}

inline MFD2 dd_exp_full(MFD2 const &x) {
  return dd_exp2_full(x * dd_pair(log2_e_hi, log2_e_lo));
}

inline MFD2 dd_log2_kernel(MFD2 const &x) {
  double hi = x._limbs[0];
  if (hi > 15.0 / 16.0 && hi < 17.0 / 16.0) {
    // direct path: t = (x - 1) / (x + 1)
    MFD2 num = x - MFD2(1.0);
    MFD2 den = x + MFD2(1.0);
    MFD2 t = num / den;
    MFD2 t_sq = t * t;
    MFD2 p = dd_horner(t_sq, log2_wide_hi, log2_wide_lo, 9);
    return t * p;
  }
  // table path
  int e = std::ilogb(hi); // matches Julia/IEEE convention (mantissa in [1,2))
  MFD2 m;
  m._limbs[0] = std::ldexp(x._limbs[0], -e);
  m._limbs[1] = std::ldexp(x._limbs[1], -e);
  int idx = static_cast<int>((m._limbs[0] - 1.0) * 32.0);
  if (idx < 0) idx = 0;
  if (idx > 31) idx = 31;
  MFD2 c = dd_pair(log2_centers[idx], 0.0);
  MFD2 v = dd_pair(log2_values_hi[idx], log2_values_lo[idx]);
  MFD2 num = m - c;
  MFD2 den = m + c;
  MFD2 t = num / den;
  MFD2 t_sq = t * t;
  MFD2 p = dd_horner(t_sq, log2_narrow_hi, log2_narrow_lo, 7);
  return dd_pair(static_cast<double>(e), 0.0) + v + t * p;
}

inline MFD2 dd_log2_full(MFD2 const &x) {
  MFD2 r;
  if (std::isnan(x._limbs[0]) || x._limbs[0] < 0.0) {
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] == 0.0) {
    r._limbs[0] = -std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
    return r;
  }
  if (!std::isfinite(x._limbs[0])) {
    r._limbs[0] = std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
    return r;
  }
  return dd_log2_kernel(x);
}

inline MFD2 dd_log_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::log(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  return dd_log2_full(x) * dd_pair(ln_2_hi, ln_2_lo);
}

inline MFD2 dd_log10_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::log10(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  return dd_log2_full(x) * dd_pair(log10_2_hi, log10_2_lo);
}

// ---- sinpi / cospi / sin / cos / tan ---------------------------------------
inline MFD2 dd_sinpi_kernel(MFD2 const &x) {
  return x * dd_horner(x * x, sinpi_coefs_hi, sinpi_coefs_lo, 15);
}
inline MFD2 dd_cospi_kernel(MFD2 const &x) {
  return dd_horner(x * x, cospi_coefs_hi, cospi_coefs_lo, 15);
}

inline MFD2 dd_sinpi_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  bool sign = x._limbs[0] < 0.0;
  MFD2 ax = sign ? -x : x;
  double n_float = std::nearbyint(2.0 * ax._limbs[0]);
  MFD2 rx = ax + dd_pair(-0.5 * n_float, 0.0);
  long long n_mod = static_cast<long long>(n_float) & 3LL;
  MFD2 res;
  switch (n_mod) {
  case 0: res = dd_sinpi_kernel(rx); break;
  case 1: res = dd_cospi_kernel(rx); break;
  case 2: res = -dd_sinpi_kernel(rx); break;
  default: res = -dd_cospi_kernel(rx); break;
  }
  return sign ? -res : res;
}

inline MFD2 dd_cospi_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 ax = (x._limbs[0] < 0.0) ? -x : x;
  double n_float = std::nearbyint(2.0 * ax._limbs[0]);
  MFD2 rx = ax + dd_pair(-0.5 * n_float, 0.0);
  long long n_mod = static_cast<long long>(n_float) & 3LL;
  switch (n_mod) {
  case 0: return dd_cospi_kernel(rx);
  case 1: return -dd_sinpi_kernel(rx);
  case 2: return -dd_cospi_kernel(rx);
  default: return dd_sinpi_kernel(rx);
  }
}

inline MFD2 dd_sin_full(MFD2 const &x) {
  return dd_sinpi_full(x * dd_pair(inv_pi_hi, inv_pi_lo));
}
inline MFD2 dd_cos_full(MFD2 const &x) {
  return dd_cospi_full(x * dd_pair(inv_pi_hi, inv_pi_lo));
}
inline MFD2 dd_tan_full(MFD2 const &x) {
  MFD2 y = x * dd_pair(inv_pi_hi, inv_pi_lo);
  return dd_sinpi_full(y) / dd_cospi_full(y);
}

// ---- sinh / cosh / tanh ----------------------------------------------------
inline MFD2 dd_sinh_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::sinh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) < 0.1) {
    // 9-term Taylor: x * (1 + x²/6 + x⁴/120 + ...)
    return x * dd_horner(x * x, sinh_taylor_hi, sinh_taylor_lo, 9);
  }
  MFD2 e = dd_exp_full(x);
  MFD2 ei = dd_exp_full(-x);
  return (e - ei) * dd_pair(0.5, 0.0);
}

inline MFD2 dd_cosh_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::cosh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 e = dd_exp_full(x);
  MFD2 ei = dd_exp_full(-x);
  return (e + ei) * dd_pair(0.5, 0.0);
}

inline MFD2 dd_tanh_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::tanh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] > 20.0) return MFD2(1.0);
  if (x._limbs[0] < -20.0) return MFD2(-1.0);
  if (std::abs(x._limbs[0]) < 0.5) {
    return dd_sinh_full(x) / dd_cosh_full(x);
  }
  // tanh(|x|) = 1 - 2/(e^(2|x|) + 1) — overflow-safe near ±1
  bool sign = x._limbs[0] < 0.0;
  MFD2 ax = sign ? -x : x;
  MFD2 two_ax = ax + ax;
  MFD2 e2 = dd_exp_full(two_ax);
  MFD2 res = MFD2(1.0) - dd_pair(2.0, 0.0) / (e2 + MFD2(1.0));
  return sign ? -res : res;
}

// ---- pow -------------------------------------------------------------------
inline MFD2 dd_pow_full(MFD2 const &x, MFD2 const &y) {
  if (x._limbs[0] == 0.0 && y._limbs[0] == 0.0) return MFD2(1.0);
  return dd_exp_full(y * dd_log_full(x));
}

// ---- asin / acos / atan / atan2 (Newton on full-DD forward) ---------------
inline MFD2 dd_asin_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::asin(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) > 1.0) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 y0 = dd_pair(std::asin(x._limbs[0]), 0.0);
  MFD2 sy = dd_sin_full(y0);
  MFD2 cy = dd_cos_full(y0);
  if (std::abs(cy._limbs[0]) < 1e-15) return y0;
  return y0 + (x - sy) / cy;
}

inline MFD2 dd_acos_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::acos(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) > 1.0) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 y0 = dd_pair(std::acos(x._limbs[0]), 0.0);
  MFD2 sy = dd_sin_full(y0);
  MFD2 cy = dd_cos_full(y0);
  if (std::abs(sy._limbs[0]) < 1e-15) return y0;
  return y0 - (x - cy) / sy;
}

inline MFD2 dd_atan_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::atan(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  // For |x| > 1 the Newton residual (x - tan(y0)) cancels to ~ulp(x),
  // capping precision at ~log2(x) bits lost. Apply the identity
  // atan(x) = sign(x)·π/2 - atan(1/x) to keep the Newton argument in
  // [-1, 1], where the residual stays well-scaled.
  bool use_recip = std::abs(x._limbs[0]) > 1.0;
  MFD2 arg = use_recip ? (dd_pair(1.0, 0.0) / x) : x;
  MFD2 y0 = dd_pair(std::atan(arg._limbs[0]), 0.0);
  MFD2 sy = dd_sin_full(y0);
  MFD2 cy = dd_cos_full(y0);
  MFD2 res;
  if (std::abs(cy._limbs[0]) < 1e-15) {
    res = y0;
  } else {
    MFD2 t = sy / cy;
    MFD2 c2 = cy * cy;
    res = y0 + c2 * (arg - t);
  }
  if (use_recip) {
    MFD2 pi_half_dd = dd_pair(1.5707963267948966, 6.123233995736766e-17);
    if (x._limbs[0] > 0.0) {
      res = pi_half_dd - res;
    } else {
      // atan(x) = -π/2 - atan(1/x), with 1/x < 0 giving res < 0.
      res._limbs[0] = -res._limbs[0];
      res._limbs[1] = -res._limbs[1];
      res = res - pi_half_dd;
    }
  }
  return res;
}

inline MFD2 dd_atan2_full(MFD2 const &y, MFD2 const &x) {
  if (!std::isfinite(x._limbs[0]) || !std::isfinite(y._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::atan2(y._limbs[0], x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] == 0.0 && y._limbs[0] == 0.0) {
    MFD2 r;
    r._limbs[0] = std::atan2(y._limbs[0], x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 pi_dd = dd_pair(3.141592653589793, 1.2246467991473532e-16);
  MFD2 pi_half_dd = dd_pair(1.5707963267948966, 6.123233995736766e-17);
  MFD2 res;
  if (std::abs(x._limbs[0]) >= std::abs(y._limbs[0])) {
    res = dd_atan_full(y / x);
    if (x._limbs[0] < 0.0) {
      if (y._limbs[0] >= 0.0) res = res + pi_dd;
      else res = res - pi_dd;
    }
  } else {
    res = dd_atan_full(x / y);
    if (y._limbs[0] > 0.0) res = pi_half_dd - res;
    else res = -pi_half_dd - res;
  }
  return res;
}

// ---- asinh / acosh / atanh -------------------------------------------------
inline MFD2 dd_asinh_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::asinh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) < 0.01) {
    return x * dd_horner(x * x, asinh_taylor_hi, asinh_taylor_lo, 15);
  }
  bool sign = x._limbs[0] < 0.0;
  MFD2 ax = sign ? -x : x;
  if (ax._limbs[0] > 1e150) {
    // log(2|x|) asymptotic to avoid x² overflow
    MFD2 r = dd_log_full(ax) + dd_pair(ln_2_hi, ln_2_lo);
    return sign ? -r : r;
  }
  MFD2 root = sqrt(ax * ax + MFD2(1.0));
  MFD2 res = dd_log_full(ax + root);
  return sign ? -res : res;
}

inline MFD2 dd_acosh_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::acosh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] < 1.0) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] > 1e150) {
    return dd_log_full(x) + dd_pair(ln_2_hi, ln_2_lo);
  }
  MFD2 root = sqrt(x * x - MFD2(1.0));
  return dd_log_full(x + root);
}

inline MFD2 dd_atanh_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::atanh(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) > 1.0) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  if (std::abs(x._limbs[0]) < 0.01) {
    return x * dd_horner(x * x, atanh_taylor_hi, atanh_taylor_lo, 15);
  }
  MFD2 num = MFD2(1.0) + x;
  MFD2 den = MFD2(1.0) - x;
  return dd_pair(0.5, 0.0) * dd_log_full(num / den);
}

// ---- erf / erfc (piecewise rational, ported from libquadmath erfq.c) -------
inline MFD2 dd_erfc_full(MFD2 const &x); // forward declaration

inline MFD2 dd_erf_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::erf(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  bool neg = x._limbs[0] < 0.0;
  MFD2 ax = neg ? -x : x;

  if (ax._limbs[0] >= 1.0) {
    if (ax._limbs[0] >= 16.0) return neg ? MFD2(-1.0) : MFD2(1.0);
    MFD2 res = MFD2(1.0) - dd_erfc_full(ax);
    return neg ? -res : res;
  }

  MFD2 y;
  if (ax._limbs[0] < 0.875) {
    if (ax._limbs[0] < 1e-18)
      return x + dd_pair(erf_efx_hi, erf_efx_lo) * x;
    MFD2 z = ax * ax;
    y = ax + ax * dd_neval(z, erf_TN1_hi, erf_TN1_lo, 8) /
                  dd_deval(z, erf_TD1_hi, erf_TD1_lo, 8);
  } else {
    MFD2 a = ax - MFD2(1.0);
    y = dd_pair(erf_const_hi, erf_const_lo) +
        dd_neval(a, erf_TN2_hi, erf_TN2_lo, 8) /
        dd_deval(a, erf_TD2_hi, erf_TD2_lo, 8);
  }
  return neg ? -y : y;
}

inline MFD2 dd_erfc_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::erfc(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  bool neg = x._limbs[0] < 0.0;
  MFD2 ax = neg ? -x : x;

  if (ax._limbs[0] < 0.25) {
    return MFD2(1.0) - dd_erf_full(x);
  }

  if (ax._limbs[0] > 107.0) {
    return neg ? MFD2(2.0) : MFD2();
  }

  MFD2 res;
  if (ax._limbs[0] < 1.25) {
    // 8 sub-intervals of width 0.125 in [0.25, 1.25)
    int i = static_cast<int>(8.0 * ax._limbs[0]) - 2;
    if (i < 0) i = 0;
    if (i > 7) i = 7;
    MFD2 z = ax - MFD2(erfc_x0[i]);
    MFD2 y = dd_neval(z, erfc_RN_hi[i], erfc_RN_lo[i], 8) /
             dd_deval(z, erfc_RD_hi[i], erfc_RD_lo[i], 7);
    res = y * z + dd_pair(erfc_Cb_hi[i], erfc_Cb_lo[i]) +
          dd_pair(erfc_Ca_hi[i], erfc_Ca_lo[i]);
  } else {
    // Asymptotic: erfc(x) = (1/x)*exp(-x^2 - 0.5625 + R(1/x^2))
    if (neg && ax._limbs[0] >= 9.0) return MFD2(2.0);

    MFD2 x2 = ax * ax;
    MFD2 z = MFD2(1.0) / x2;
    int i = static_cast<int>(8.0 / ax._limbs[0]);
    if (i > 7) i = 7;
    int nd = erfc_asym_ndeg[i];
    int dd = erfc_asym_ddeg[i];
    MFD2 p = dd_neval(z, erfc_asym_N_hi[i], erfc_asym_N_lo[i], nd) /
             dd_deval(z, erfc_asym_D_hi[i], erfc_asym_D_lo[i], dd);

    // Split x for accurate exp(-x^2): truncate hi limb to ~27 bits
    // so s^2 is exact in double precision
    union {
      double d;
      uint64_t u;
    } uu;
    uu.d = ax._limbs[0];
    uu.u &= 0xFFFFFFF800000000ULL; // keep top 27 bits of mantissa
    MFD2 s = dd_pair(uu.d, 0.0);

    MFD2 e1 = dd_exp_full(-(s * s) - MFD2(0.5625));
    MFD2 diff_sq = (s - ax) * (s + ax);
    MFD2 e2 = dd_exp_full(diff_sq + p);
    res = (e1 * e2) / ax;
  }

  return neg ? MFD2(2.0) - res : res;
}

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
  // No std::digamma in <cmath>, so the N == 2 case falls back to a
  // leading-limb evaluation (single-double precision).
  MultiFloat<T, N> r;
  r._limbs[0] = std::tgamma(x._limbs[0]);
  return r;
}

template <typename T, std::size_t N>
MultiFloat<T, N> lgamma(MultiFloat<T, N> const &x) {
  MultiFloat<T, N> r;
  r._limbs[0] = std::lgamma(x._limbs[0]);
  return r;
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
