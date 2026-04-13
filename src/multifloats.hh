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
    0.24022650695910072,       0.05550410866482158,
    0.009618129107628477,      0.0013333558146428443,
    0.0001540353039338161,     1.5252733804059841e-05,
    1.3215486790144443e-06,    1.0178086009239799e-07,
    7.054911615823195e-09,     4.445538268793799e-10,
    2.567936278800348e-11,     1.36919783052689e-12};
inline constexpr double exp2_coefs_lo[14] = {
    1.1483944385261558e-34,    2.3190468138462996e-17,
    -9.493931253185757e-18,    -3.165822290391908e-18,
    2.832460796447469e-19,     1.3928061121226897e-20,
    1.1765492220661074e-20,    -8.044711758039522e-22,
    -8.477254587913276e-23,    0.0,
    0.0,                       0.0,
    0.0,                       0.0};
inline constexpr double exp2_min_d = -1022.0;
inline constexpr double exp2_max_d = 1023.9999999999999;

// ---- conversion constants -------------------------------------------------
inline constexpr double log2_e_hi = 1.4426950408889634;
inline constexpr double log2_e_lo = 2.0355273740931033e-17;
inline constexpr double ln_2_hi = 0.6931471805599453;
inline constexpr double ln_2_lo = 2.3190468138462996e-17;
inline constexpr double log10_2_hi = 0.3010299956639812;
inline constexpr double log10_2_lo = -2.8037281277851704e-18;

// ---- log2 polynomial: narrow (table path) ---------------------------------
inline constexpr double log2_narrow_hi[7] = {
    2.8853900817779268,    0.9617966939259756,
    0.5770780163555853,    0.4121985831111324,
    0.3205988979754622,    0.26230818590911165,
    0.22199354561535342};
inline constexpr double log2_narrow_lo[7] = {
    4.0710547481862066e-17,   5.0577616648017525e-17,
    5.255105933880626e-17,    1.0864685980468599e-17,
    0.0,                      0.0,
    0.0};

// ---- log2 polynomial: wide (direct path for x in [15/16, 17/16]) ----------
inline constexpr double log2_wide_hi[9] = {
    2.8853900817779268,     0.9617966939259756,
    0.5770780163555853,     0.4121985831111324,
    0.32059889797532554,    0.26230818925164,
    0.22195308467900554,    0.1923579466882983,
    0.17044172319243897};
inline constexpr double log2_wide_lo[9] = {
    4.0710547481862066e-17,   5.0577616647866686e-17,
    5.255103712360494e-17,    1.3680430918811534e-17,
    -1.0910052064370687e-17,  1.46400932825286e-17,
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

// ---- pi/2 as a 3-part constant for Cody-Waite range reduction (~161 bits) -
inline constexpr double pi_half_cw1 = 1.5707963267948966;
inline constexpr double pi_half_cw2 = 6.123233995736766e-17;
inline constexpr double pi_half_cw3 = -1.4973849048591698e-33;

// ---- sin(x)/x Taylor coefficients (13 terms, c[k] = (-1)^k / (2k+1)!) -----
inline constexpr double sin_taylor_hi[13] = {
    1.0,                       -0.16666666666666666,
    0.008333333333333333,      -0.0001984126984126984,
    2.7557319223985893e-06,    -2.505210838544172e-08,
    1.6059043836821613e-10,    -7.647163731819816e-13,
    2.8114572543455206e-15,    -8.22063524662433e-18,
    1.9572941063391263e-20,    -3.868170170630684e-23,
    6.446950284384474e-26};
inline constexpr double sin_taylor_lo[13] = {
    0.0,                       -9.25185853854297e-18,
    1.1564823173178714e-19,    -1.7209558293420705e-22,
    -1.858393274046472e-22,    1.448814070935912e-24,
    1.2585294588752098e-26,    -7.03872877733453e-30,
    1.6508842730861433e-31,    -2.2141894119604265e-34,
    -1.3643503830087908e-36,   8.843177655482344e-40,
    -1.9330404233703465e-42};

// ---- cos(x) Taylor coefficients (13 terms, c[k] = (-1)^k / (2k)!) ---------
inline constexpr double cos_taylor_hi[13] = {
    1.0,                       -0.5,
    0.041666666666666664,      -0.001388888888888889,
    2.48015873015873e-05,      -2.755731922398589e-07,
    2.08767569878681e-09,      -1.1470745597729725e-11,
    4.779477332387385e-14,     -1.5619206968586225e-16,
    4.110317623312165e-19,     -8.896791392450574e-22,
    1.6117375710961184e-24};
inline constexpr double cos_taylor_lo[13] = {
    0.0,                       0.0,
    2.3129646346357427e-18,    5.300543954373577e-20,
    2.1511947866775882e-23,    -2.3767714622250297e-23,
    -1.20734505911326e-25,     -2.0655512752830745e-28,
    4.399205485834081e-31,     -1.1910679660273754e-32,
    1.4412973378659527e-36,    7.911402614872376e-38,
    -3.6846573564509766e-41};

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

// ---- erf Taylor (50 coefs, c[n] = (2/sqrt(pi))*(-1)^n/(n!*(2n+1))) -------
inline constexpr double erf_coefs_hi[50] = {
    1.1283791670955126,       -0.37612638903183754,
    0.11283791670955126,      -0.026866170645131252,
    0.005223977625442188,     -0.0008548327023450853,
    0.00012055332981789664,   -1.492565035840625e-05,
    1.6462114365889248e-06,   -1.6365844691234924e-07,
    1.4807192815879218e-08,   -1.2290555301717928e-09,
    9.422759064650411e-11,    -6.7113668551641105e-12,
    4.4632242632864775e-13,   -2.7835162072109215e-14,
    1.6342614095367152e-15,   -9.063970842808673e-17,
    4.763348040515068e-18,    -2.3784598852774293e-19,
    1.131218725924631e-20,    -5.136209054585811e-22,
    2.2308786802746453e-23,   -9.28672901131906e-25,
    3.71153285316323e-26,     -1.4263930180784176e-27,
    5.279103332510834e-29,    -1.8841244217042036e-30,
    6.492909974544561e-32,    -2.1630383901171243e-33,
    6.973730328792914e-35,    -2.1781748594796097e-36,
    6.597356545539203e-38,    -1.9395213725013487e-39,
    5.539127534424142e-41,    -1.538027363683162e-42,
    4.1552489658106737e-44,   -1.0930925207357809e-45,
    2.801843440026779e-47,    -7.002335114640117e-49,
    1.7073594878289175e-50,   -4.0639470618319806e-52,
    9.448392328628975e-54,    -2.1467878854142287e-55,
    4.769421502324767e-57,    -1.0365775670498273e-58,
    2.2049686442621386e-60,   -4.59265585479012e-62,
    9.370753999249602e-64,    -1.8737644566629794e-65};
inline constexpr double erf_coefs_lo[50] = {
    1.533545961316588e-17,    1.3391897206030649e-17,
    -4.017569161809194e-18,   4.6092880729453e-19,
    -8.962504586282528e-20,   5.0148896786169737e-20,
    6.480246840070509e-21,    -6.248427055364001e-22,
    -1.0547266132407653e-22,  1.5075323135139275e-24,
    -3.254656350443331e-25,   9.976105519856072e-26,
    -2.8231273253303265e-27,  1.0441838553137576e-28,
    -1.415652238799887e-29,   1.4189660114017355e-30,
    -1.159587670993415e-32,   3.004743561097066e-33,
    -1.856680962530252e-34,   -1.5875374697145176e-35,
    7.245880865990418e-38,    -1.4537189752115876e-38,
    -4.7283897669065125e-40,  -2.2114070452037343e-41,
    -2.864576959470938e-42,   -3.298130029636975e-44,
    2.3921205704042427e-45,   -1.106406475814513e-47,
    1.2592869837235584e-48,   -1.2253141605940307e-49,
    4.242571852198481e-51,    -1.3045967020940556e-52,
    4.177193607879865e-55,    5.869824242611732e-56,
    -2.1631332158779208e-57,  3.262876687823446e-60,
    -5.5037931014990846e-61,  5.540458322951681e-62,
    9.259808565197577e-64,    1.3200652030502501e-65,
    -9.490441382723132e-67,   -3.3631583324243266e-68,
    2.7994711326647915e-70,   1.7816791770026858e-71,
    1.5216055227045837e-74,   5.5040127417898383e-76,
    -1.3741888831536437e-76,  -1.6966388960870593e-78,
    -2.28824162542219e-80,    -7.09889581208805e-82};

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

// ---- sin / cos / tan -------------------------------------------------------
// Cody-Waite range reduction: n = nearest(x·2/π), r = x − n·(π/2) using a
// 3-part π/2 constant (~161 bits). Each n·pi_half_cwK is an exact DD
// product via FMA, then subtracted from r as full DD. This preserves ~106
// bits through the reduction — the critical difference from the previous
// sinpi(x·inv_pi) path, which lost the integer part of x/π for large |x|.
inline void dd_reduce_pi_half(MFD2 const &x, MFD2 &r, int &n_mod4) {
  double n_float = std::nearbyint(x._limbs[0] * 0.6366197723675814);
  r = x;
  MFD2 npi;
  npi._limbs[0] = n_float * pi_half_cw1;
  npi._limbs[1] = std::fma(n_float, pi_half_cw1, -npi._limbs[0]);
  r = r - npi;
  npi._limbs[0] = n_float * pi_half_cw2;
  npi._limbs[1] = std::fma(n_float, pi_half_cw2, -npi._limbs[0]);
  r = r - npi;
  npi._limbs[0] = n_float * pi_half_cw3;
  npi._limbs[1] = 0.0;
  r = r - npi;
  long long nn = static_cast<long long>(std::fabs(n_float)) & 3LL;
  if (n_float < 0.0) nn = (4 - nn) & 3;
  n_mod4 = static_cast<int>(nn);
}

// Taylor kernels for |x| ≤ π/8. sin evaluates as x·poly(x²) so the low
// limb of x is preserved losslessly.
inline MFD2 dd_sin_kernel(MFD2 const &x) {
  MFD2 x2 = x * x;
  return dd_horner(x2, sin_taylor_hi, sin_taylor_lo, 13) * x;
}
inline MFD2 dd_cos_kernel(MFD2 const &x) {
  MFD2 x2 = x * x;
  return dd_horner(x2, cos_taylor_hi, cos_taylor_lo, 13);
}

// Evaluate sin(r)/cos(r) for |r| ≤ π/4. For |r| > π/8 shift by π/4 and use
// the angle-addition identity — this halves the polynomial argument range
// (x² ≤ 0.154 vs 0.616) so the 13-term Taylor hits full DD at the boundary.
inline MFD2 dd_sin_eval(MFD2 const &r) {
  constexpr double pi8 = 0.392699081698724;
  if (std::fabs(r._limbs[0]) <= pi8) return dd_sin_kernel(r);
  MFD2 pi4_dd = dd_pair(0.7853981633974483, 3.061616997868383e-17);
  MFD2 inv_sqrt2 = dd_pair(0.7071067811865476, -4.833646656726457e-17);
  bool pos = r._limbs[0] > 0.0;
  MFD2 rp = pos ? (r - pi4_dd) : (r + pi4_dd);
  MFD2 sk = dd_sin_kernel(rp);
  MFD2 ck = dd_cos_kernel(rp);
  return pos ? (sk + ck) * inv_sqrt2 : (sk - ck) * inv_sqrt2;
}

inline MFD2 dd_cos_eval(MFD2 const &r) {
  constexpr double pi8 = 0.392699081698724;
  if (std::fabs(r._limbs[0]) <= pi8) return dd_cos_kernel(r);
  MFD2 pi4_dd = dd_pair(0.7853981633974483, 3.061616997868383e-17);
  MFD2 inv_sqrt2 = dd_pair(0.7071067811865476, -4.833646656726457e-17);
  bool pos = r._limbs[0] > 0.0;
  MFD2 rp = pos ? (r - pi4_dd) : (r + pi4_dd);
  MFD2 sk = dd_sin_kernel(rp);
  MFD2 ck = dd_cos_kernel(rp);
  return pos ? (ck - sk) * inv_sqrt2 : (ck + sk) * inv_sqrt2;
}

inline MFD2 dd_sin_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 r;
  int q;
  dd_reduce_pi_half(x, r, q);
  switch (q) {
  case 0: return dd_sin_eval(r);
  case 1: return dd_cos_eval(r);
  case 2: return -dd_sin_eval(r);
  default: return -dd_cos_eval(r);
  }
}

inline MFD2 dd_cos_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 r;
  int q;
  dd_reduce_pi_half(x, r, q);
  switch (q) {
  case 0: return dd_cos_eval(r);
  case 1: return -dd_sin_eval(r);
  case 2: return -dd_cos_eval(r);
  default: return dd_sin_eval(r);
  }
}

inline MFD2 dd_tan_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::numeric_limits<double>::quiet_NaN();
    r._limbs[1] = 0.0;
    return r;
  }
  MFD2 r;
  int q;
  dd_reduce_pi_half(x, r, q);
  switch (q) {
  case 0: return dd_sin_eval(r) / dd_cos_eval(r);
  case 1: return -dd_cos_eval(r) / dd_sin_eval(r);
  case 2: return dd_sin_eval(r) / dd_cos_eval(r);
  default: return -dd_cos_eval(r) / dd_sin_eval(r);
  }
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
  // For |x| > 37, exp(-2|x|) < 2^-106, so tanh(x) = ±1 to full DD.
  if (x._limbs[0] > 37.0) return MFD2(1.0);
  if (x._limbs[0] < -37.0) return MFD2(-1.0);
  if (std::abs(x._limbs[0]) < 0.5) {
    return dd_sinh_full(x) / dd_cosh_full(x);
  }
  // tanh(|x|) = (1 - em2) / (1 + em2) where em2 = exp(-2|x|) is small
  // for |x| > 0.5 — avoids the 1 − small cancellation of 1 − 2/(e^(2|x|)+1).
  bool sign = x._limbs[0] < 0.0;
  MFD2 neg_two_ax;
  neg_two_ax._limbs[0] = sign ? (2.0 * x._limbs[0]) : (-2.0 * x._limbs[0]);
  neg_two_ax._limbs[1] = sign ? (2.0 * x._limbs[1]) : (-2.0 * x._limbs[1]);
  MFD2 em2 = dd_exp_full(neg_two_ax);
  MFD2 res = (MFD2(1.0) - em2) / (MFD2(1.0) + em2);
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

// ---- erf / erfc ------------------------------------------------------------
inline MFD2 dd_erf_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::erf(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  bool sign = x._limbs[0] < 0.0;
  MFD2 ax = sign ? -x : x;
  MFD2 res;
  if (ax._limbs[0] < 2.0) {
    // 50-term Taylor
    res = ax * dd_horner(ax * ax, erf_coefs_hi, erf_coefs_lo, 50);
  } else {
    // 1 - libm_erfc(x) — DD subtraction preserves the lo limb
    res = MFD2(1.0) - dd_pair(std::erfc(ax._limbs[0]), 0.0);
  }
  return sign ? -res : res;
}

inline MFD2 dd_erfc_full(MFD2 const &x) {
  if (!std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = std::erfc(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  if (x._limbs[0] >= 6.0) {
    MFD2 r;
    r._limbs[0] = std::erfc(x._limbs[0]);
    r._limbs[1] = 0.0;
    return r;
  }
  return MFD2(1.0) - dd_erf_full(x);
}

// ---- tgamma / lgamma -------------------------------------------------------
// Route 1: derivative correction via a double-precision digamma.
//     tgamma(hi + lo) ≈ tgamma(hi) · (1 + ψ(hi) · lo)
//     lgamma(hi + lo) ≈ lgamma(hi) + ψ(hi) · lo
// The leading-limb libm call + DD correction gives ~106 bits of precision
// as long as ψ(hi) is good to ~1e-16 — which it is for the Stirling
// asymptotic + recurrence + reflection pipeline below.
//
// This lives in its own sub-namespace so a future full-DD Stirling
// implementation (route 2) can coexist and be benchmarked head-to-head.
namespace gamma_v1 {

// Double-precision digamma ψ(x) = d/dx ln Γ(x).
//   - x ≥ 12: 6-term Stirling asymptotic in 1/x² reaches ~6e-17 relative.
//   - 0 < x < 12: upward recurrence ψ(x) = ψ(x+1) − 1/x shifts into the
//     asymptotic range, accumulating the correction sum as we go.
//   - x < 0.5: reflection ψ(x) = ψ(1−x) − π·cot(π·x).
// Near non-positive integer poles, tgamma/lgamma themselves return
// inf/NaN, so the correction never runs there; we still guard `sin(πx)=0`.
inline double digamma_dp(double x) {
  if (x < 0.5) {
    double s = std::sin(M_PI * x);
    if (s == 0.0) return std::numeric_limits<double>::quiet_NaN();
    double c = std::cos(M_PI * x);
    return digamma_dp(1.0 - x) - M_PI * (c / s);
  }
  double result = 0.0;
  while (x < 12.0) {
    result -= 1.0 / x;
    x += 1.0;
  }
  double inv_x = 1.0 / x;
  double inv_x2 = inv_x * inv_x;
  // Coefficients B_{2k}/(2k) for k=1..6:
  //   1/12, -1/120, 1/252, -1/240, 1/132, -691/32760
  // Series: -Σ c_k · inv_x^{2k}  (Horner on inv_x2, outermost negation
  //         folds into the final subtraction).
  constexpr double c1 =  1.0 / 12.0;
  constexpr double c2 = -1.0 / 120.0;
  constexpr double c3 =  1.0 / 252.0;
  constexpr double c4 = -1.0 / 240.0;
  constexpr double c5 =  1.0 / 132.0;
  constexpr double c6 = -691.0 / 32760.0;
  double asymp = inv_x2 *
                 (c1 + inv_x2 *
                           (c2 + inv_x2 *
                                     (c3 + inv_x2 *
                                               (c4 + inv_x2 *
                                                         (c5 + inv_x2 * c6)))));
  result += std::log(x) - 0.5 * inv_x - asymp;
  return result;
}

inline MFD2 tgamma_full(MFD2 const &x) {
  double fhi = std::tgamma(x._limbs[0]);
  if (!std::isfinite(fhi) || !std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = fhi;
    r._limbs[1] = 0.0;
    return r;
  }
  double psi = digamma_dp(x._limbs[0]);
  double corr = fhi * psi * x._limbs[1];
  double s = fhi + corr;
  double bp = s - fhi;
  MFD2 r;
  r._limbs[0] = s;
  r._limbs[1] = corr - bp;
  return r;
}

inline MFD2 lgamma_full(MFD2 const &x) {
  double fhi = std::lgamma(x._limbs[0]);
  if (!std::isfinite(fhi) || !std::isfinite(x._limbs[0])) {
    MFD2 r;
    r._limbs[0] = fhi;
    r._limbs[1] = 0.0;
    return r;
  }
  double psi = digamma_dp(x._limbs[0]);
  double corr = psi * x._limbs[1];
  double s = fhi + corr;
  double bp = s - fhi;
  MFD2 r;
  r._limbs[0] = s;
  r._limbs[1] = corr - bp;
  return r;
}

} // namespace gamma_v1

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
    return detail::gamma_v1::tgamma_full(x);
  } else {
    MultiFloat<T, N> r;
    r._limbs[0] = std::tgamma(x._limbs[0]);
    return r;
  }
}

template <typename T, std::size_t N>
MultiFloat<T, N> lgamma(MultiFloat<T, N> const &x) {
  if constexpr (N == 2 && std::is_same_v<T, double>) {
    return detail::gamma_v1::lgamma_full(x);
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
