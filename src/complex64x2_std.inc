// std::complex<multifloats::float64x2> specializations for C++17 <complex>
// free-function templates. Pulled into multifloats_math.cc inside
// `namespace std { ... }` so the explicit specialization syntax lands in
// the right enclosing namespace.
//
// Branch cuts follow C99 Annex G (matching libquadmath cexpq/clogq/csqrtq).
// Where the real-axis formula needs both sin(y) and cos(y) (or both sinh
// and cosh), these call the fused multifloats::sincos / multifloats::sinhcosh
// kernels so one range reduction + Taylor pair covers both outputs.
// Functions built on top of std::sqrt / std::log (asin/acos/...) call the
// std:: specialization directly — within this TU the compiler sees the
// complete definition chain.
//
// Shared anon-namespace helpers `dd_cross_diff` and `dd_x2y2m1` are pulled
// in at file scope in multifloats_math.cc; unqualified name lookup from
// inside namespace std walks up to global scope to find them.

// ---- Exponential / logarithmic --------------------------------------------

// exp(a+bi) = e^a · (cos b + i·sin b). 2 transcendentals (exp + sincos).
// Annex G.6.3.1 special-value matrix: most edge cases would give NaN through
// the bare formula (e^∞·sin(0) = ∞·0 = NaN). Handle them up front.
template <>
complex<multifloats::float64x2>
exp(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 a = z.real(), b = z.imag();
  double a_hi = a._limbs[0], b_hi = b._limbs[0];
  bool b_zero = (b_hi == 0.0 && b._limbs[1] == 0.0);
  if (!std::isfinite(a_hi) || !std::isfinite(b_hi)) {
    double inf = std::numeric_limits<double>::infinity();
    double nan = std::numeric_limits<double>::quiet_NaN();
    if (std::isnan(a_hi)) {
      // NaN + i·0  → NaN + i·0; otherwise NaN + i·NaN.
      return b_zero ? complex<float64x2>(float64x2(nan, 0.0), b)
                    : complex<float64x2>(float64x2(nan, 0.0), float64x2(nan, 0.0));
    }
    if (std::isinf(b_hi)) {
      // finite + i·∞ or ±∞ + i·∞ → NaN + i·NaN (invalid).
      return complex<float64x2>(float64x2(nan, 0.0), float64x2(nan, 0.0));
    }
    if (std::isnan(b_hi)) {
      // ±∞ + i·NaN: +∞ → +∞ + i·NaN; -∞ → ±0 + i·±0.
      if (a_hi > 0.0)
        return complex<float64x2>(float64x2(inf, 0.0), float64x2(nan, 0.0));
      return complex<float64x2>(float64x2(0.0, 0.0), float64x2(0.0, 0.0));
    }
    // a is ±∞, b is finite.
    if (a_hi > 0.0) {
      // +∞ + i·0 → +∞ + i·0 (preserve b's sign).
      if (b_zero) return complex<float64x2>(float64x2(inf, 0.0), b);
      // +∞ · cis(b): each component is ±∞ with sign from cos(b)/sin(b).
      float64x2 s, c;
      multifloats::sincos(b, s, c);
      return complex<float64x2>(float64x2(std::copysign(inf, c._limbs[0]), 0.0),
                                float64x2(std::copysign(inf, s._limbs[0]), 0.0));
    }
    // -∞ + i·y, y finite: +0 · cis(y).
    if (b_zero) return complex<float64x2>(float64x2(0.0, 0.0), b);
    float64x2 s, c;
    multifloats::sincos(b, s, c);
    return complex<float64x2>(float64x2(std::copysign(0.0, c._limbs[0]), 0.0),
                              float64x2(std::copysign(0.0, s._limbs[0]), 0.0));
  }
  multifloats::float64x2 ea = multifloats::exp(a);
  multifloats::float64x2 s, c;
  multifloats::sincos(b, s, c);
  return complex<multifloats::float64x2>(ea * c, ea * s);
}

// log(a+bi) = log(|z|) + i·atan2(b, a). Two paths for log(|z|):
//   • |z| near 1 (big ∈ [0.5, 2]):  log|z| = 0.5·log1p(a²+b²−1) via
//     dd_x2y2m1 — the general formula loses bits when log(big) and
//     0.5·log(1+ratio²) cancel (both ≈ 0 ± ε).
//   • otherwise: log(big) + 0.5·log(1 + (small/big)²), which avoids
//     hypot overflow when |z| exceeds DBL_MAX yet log|z| is finite.
// atan2 handles the negative-real-axis cut and signed-zero propagation.
template <>
complex<multifloats::float64x2>
log(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 a = z.real(), b = z.imag();
  float64x2 ax = multifloats::abs(a);
  float64x2 ay = multifloats::abs(b);
  float64x2 big   = (ax > ay) ? ax : ay;
  float64x2 small = (ax > ay) ? ay : ax;
  float64x2 r;
  if (big._limbs[0] == 0.0) {
    r._limbs[0] = -std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
  } else if (big._limbs[0] >= 0.5 && big._limbs[0] <= 2.0) {
    float64x2 u = dd_x2y2m1(big, small);
    r = float64x2(0.5, 0.0) * multifloats::log1p(u);
  } else {
    float64x2 ratio = small / big;
    float64x2 one = float64x2(1.0, 0.0);
    float64x2 half = float64x2(0.5, 0.0);
    r = multifloats::log(big) + half * multifloats::log(one + ratio * ratio);
  }
  float64x2 phi = multifloats::atan2(b, a);
  return complex<multifloats::float64x2>(r, phi);
}

// sqrt(z). Principal branch; cut on negative real axis, continuous above.
// Uses Moshier's 2·Re·Im = Im identity to avoid cancellation in mod ± a.
template <>
complex<multifloats::float64x2>
sqrt(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 a = z.real(), b = z.imag();
  if (a._limbs[0] == 0.0 && b._limbs[0] == 0.0)
    return complex<multifloats::float64x2>(float64x2(0.0), b);  // preserves signed zero of imag
  float64x2 mod = multifloats::hypot(a, b);
  float64x2 half = float64x2(0.5, 0.0);
  if (a._limbs[0] >= 0.0) {
    float64x2 r = multifloats::sqrt((mod + a) * half);
    float64x2 i = b / (r + r);
    return complex<multifloats::float64x2>(r, i);
  } else {
    float64x2 s = multifloats::sqrt((mod - a) * half);
    float64x2 r = multifloats::abs(b) / (s + s);
    // Use signbit so b = -0 picks the −s branch (IEEE: -0 < 0 is false,
    // so `< 0.0` mishandles the negative-real-axis branch cut per C99 G.6.4.2).
    float64x2 i = std::signbit(b._limbs[0]) ? -s : s;
    return complex<multifloats::float64x2>(r, i);
  }
}

// pow(z, w) = exp(w · log(z)). C99 G.6.4.1; 0^w folds to 0.
template <>
complex<multifloats::float64x2>
pow(complex<multifloats::float64x2> const &z,
    complex<multifloats::float64x2> const &w) {
  using multifloats::float64x2;
  float64x2 a = z.real(), b = z.imag();
  // True zero requires *both* limbs of *both* components to be zero. A
  // hi-only test misclassifies subnormal DDs of the form (hi=+0, lo=ε)
  // as zero and returns 0 instead of exp(w · log(ε)).
  if (a._limbs[0] == 0.0 && a._limbs[1] == 0.0 &&
      b._limbs[0] == 0.0 && b._limbs[1] == 0.0) {
    return complex<multifloats::float64x2>(float64x2(0.0), float64x2(0.0));
  }
  // Renormalize each component if its hi limb is zero but lo is not —
  // otherwise std::log's log2 short-circuit would see hi == 0 and return
  // (−∞, …), propagating NaNs back through exp(w · log(z)).
  if (a._limbs[0] == 0.0 && a._limbs[1] != 0.0) { a = float64x2(a._limbs[1], 0.0); }
  if (b._limbs[0] == 0.0 && b._limbs[1] != 0.0) { b = float64x2(b._limbs[1], 0.0); }
  complex<multifloats::float64x2> zn(a, b);
  complex<multifloats::float64x2> l = std::log(zn);
  float64x2 lr = l.real(), li = l.imag();
  float64x2 wr = w.real(), wi = w.imag();
  // w * log(z) = (wr·lr − wi·li) + i·(wr·li + wi·lr).
  // The dominant error here is |Im(w·log z)|·ulp_dd amplified by exp;
  // this scales with |w|·|log z| ≈ a few hundred DD ulp at fuzz extremes.
  // Compensated cross-diff for the real part doesn't move the needle.
  complex<multifloats::float64x2> p(wr * lr - wi * li, wr * li + wi * lr);
  return std::exp(p);
}

// ---- Trigonometric ---------------------------------------------------------

// sin(a+bi) = sin(a)·cosh(b) + i·cos(a)·sinh(b). 2 fused transcendentals.
template <>
complex<multifloats::float64x2>
sin(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 a = z.real(), b = z.imag();
  float64x2 sa, ca, sb, cb;
  multifloats::sincos(a, sa, ca);
  multifloats::sinhcosh(b, sb, cb);
  return complex<multifloats::float64x2>(sa * cb, ca * sb);
}

// cos(a+bi) = cos(a)·cosh(b) − i·sin(a)·sinh(b). 2 fused transcendentals.
template <>
complex<multifloats::float64x2>
cos(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 a = z.real(), b = z.imag();
  float64x2 sa, ca, sb, cb;
  multifloats::sincos(a, sa, ca);
  multifloats::sinhcosh(b, sb, cb);
  return complex<multifloats::float64x2>(ca * cb, -(sa * sb));
}

// tan(a+bi) = (sin(a)·cos(a) + i·sinh(b)·cosh(b)) / (cos(a)² + sinh(b)²).
// libquadmath ctanq formula; 2 fused transcendentals + real div.
template <>
complex<multifloats::float64x2>
tan(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 a = z.real(), b = z.imag();
  float64x2 sa, ca, sb, cb;
  multifloats::sincos(a, sa, ca);
  multifloats::sinhcosh(b, sb, cb);
  float64x2 den = ca * ca + sb * sb;
  // Pole at a = π/2+kπ with b = 0 ⇒ den = 0. C99 Annex G.6.3.1
  // says ctan at a real-axis pole is (NaN, NaN).
  if (den._limbs[0] == 0.0) {
    double nan = std::numeric_limits<double>::quiet_NaN();
    return complex<multifloats::float64x2>(float64x2(nan, 0.0),
                                            float64x2(nan, 0.0));
  }
  // Overflow in sinh(b) for |b| ≳ 710 pushes den to +∞. The mathematical
  // limit is tan(a + i·b) → i·sign(b); bypass the division so the imag
  // part doesn't come out as ∞/∞ = NaN.
  if (!std::isfinite(den._limbs[0])) {
    float64x2 im_sign(std::copysign(1.0, b._limbs[0]), 0.0);
    return complex<multifloats::float64x2>(float64x2(0.0, 0.0), im_sign);
  }
  return complex<multifloats::float64x2>((sa * ca) / den, (sb * cb) / den);
}

// ---- Hyperbolic ------------------------------------------------------------

// sinh(a+bi) = sinh(a)·cos(b) + i·cosh(a)·sin(b). 2 fused transcendentals.
template <>
complex<multifloats::float64x2>
sinh(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 a = z.real(), b = z.imag();
  float64x2 sa, ca, sb, cb;
  multifloats::sinhcosh(a, sa, ca);
  multifloats::sincos(b, sb, cb);
  return complex<multifloats::float64x2>(sa * cb, ca * sb);
}

// cosh(a+bi) = cosh(a)·cos(b) + i·sinh(a)·sin(b). 2 fused transcendentals.
template <>
complex<multifloats::float64x2>
cosh(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 a = z.real(), b = z.imag();
  float64x2 sa, ca, sb, cb;
  multifloats::sinhcosh(a, sa, ca);
  multifloats::sincos(b, sb, cb);
  return complex<multifloats::float64x2>(ca * cb, sa * sb);
}

// tanh(a+bi) = (sinh(a)·cosh(a) + i·sin(b)·cos(b)) / (sinh(a)² + cos(b)²).
// libquadmath ctanhq formula.
template <>
complex<multifloats::float64x2>
tanh(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 a = z.real(), b = z.imag();
  float64x2 sa, ca, sb, cb;
  multifloats::sinhcosh(a, sa, ca);
  multifloats::sincos(b, sb, cb);
  float64x2 den = sa * sa + cb * cb;
  // Pole at a = 0 with b = π/2+kπ ⇒ den = 0. Annex G.6.2.6 says ctanh at
  // an imaginary-axis pole is (NaN, NaN).
  if (den._limbs[0] == 0.0) {
    double nan = std::numeric_limits<double>::quiet_NaN();
    return complex<multifloats::float64x2>(float64x2(nan, 0.0),
                                            float64x2(nan, 0.0));
  }
  // Overflow in sinh(a) for |a| ≳ 710 pushes den to +∞. Mathematical limit
  // is tanh(a + i·b) → sign(a); bypass division so it doesn't return NaN.
  if (!std::isfinite(den._limbs[0])) {
    float64x2 re_sign(std::copysign(1.0, a._limbs[0]), 0.0);
    return complex<multifloats::float64x2>(re_sign, float64x2(0.0, 0.0));
  }
  return complex<multifloats::float64x2>((sa * ca) / den, (sb * cb) / den);
}

// ---- Inverse hyperbolic ----------------------------------------------------

// asinh(z). Branch structure ported from libquadmath's
// __quadmath_kernel_casinhq: operates on |Re(z)|, |Im(z)|, restores
// signs at the end. The textbook log(z + sqrt(z² + 1)) form collapses
// to ~double precision once z sits near the unit circle, because
// arg = z + sqrt(1 + z²) and 1 carry the same leading digits — DD's
// 106 bits can't resolve log(|arg|) when |arg|² − 1 ≈ 0. Each branch
// rearranges the formula so the small result is computed from
// non-cancelling pieces, summed through log1p.
template <>
complex<multifloats::float64x2>
asinh(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 a = z.real(), b = z.imag();
  float64x2 rx = multifloats::abs(a);
  float64x2 ix = multifloats::abs(b);
  // DD eps = 2⁻¹⁰⁵ ≈ 2.46e-32; constants mirror the __float128 source.
  static constexpr double dd_inv_eps    = 0x1.0p105;
  static constexpr double dd_eps_over_8 = 0x1.0p-108;
  static constexpr double dd_eps_sq     = 0x1.0p-210;
  float64x2 one(1.0, 0.0), half(0.5, 0.0), two(2.0, 0.0);
  float64x2 re, im;

  if (rx._limbs[0] >= dd_inv_eps || ix._limbs[0] >= dd_inv_eps) {
    // |z| huge: z + sqrt(1+z²) ≈ 2z; avoid overflow from squaring.
    complex<multifloats::float64x2> l = std::log(complex<multifloats::float64x2>(rx, ix));
    re = l.real() + float64x2(ln_2_hi, ln_2_lo);
    im = l.imag();
  } else if (rx._limbs[0] >= 0.5 && ix._limbs[0] < dd_eps_over_8) {
    float64x2 s = multifloats::hypot(one, rx);
    re = multifloats::log(rx + s);
    im = multifloats::atan2(ix, s);
  } else if (rx._limbs[0] < dd_eps_over_8 && ix._limbs[0] >= 1.5) {
    float64x2 s = multifloats::sqrt((ix + one) * (ix - one));
    re = multifloats::log(ix + s);
    im = multifloats::atan2(s, rx);
  } else if (ix._limbs[0] > 1.0 && ix._limbs[0] < 1.5 && rx._limbs[0] < 0.5) {
    float64x2 ix2m1 = (ix + one) * (ix - one);
    if (rx._limbs[0] < dd_eps_sq) {
      float64x2 s = multifloats::sqrt(ix2m1);
      re = half * multifloats::log1p(two * (ix2m1 + ix * s));
      im = multifloats::atan2(s, rx);
    } else {
      float64x2 rx2 = rx * rx;
      float64x2 f = rx2 * (two + rx2 + two * ix * ix);
      float64x2 d = multifloats::sqrt(ix2m1 * ix2m1 + f);
      float64x2 dp = d + ix2m1;
      float64x2 dm = f / dp;
      float64x2 r1 = multifloats::sqrt((dm + rx2) * half);
      float64x2 r2 = rx * ix / r1;
      re = half * multifloats::log1p(rx2 + dp + two * (rx * r1 + ix * r2));
      im = multifloats::atan2(ix + r2, rx + r1);
    }
  } else if (ix._limbs[0] == 1.0 && ix._limbs[1] == 0.0 && rx._limbs[0] < 0.5) {
    if (rx._limbs[0] < dd_eps_over_8) {
      float64x2 sr = multifloats::sqrt(rx);
      re = half * multifloats::log1p(two * (rx + sr));
      im = multifloats::atan2(one, sr);
    } else {
      float64x2 rx2 = rx * rx;
      float64x2 d = rx * multifloats::sqrt(float64x2(4.0, 0.0) + rx2);
      float64x2 s1 = multifloats::sqrt((d + rx2) * half);
      float64x2 s2 = multifloats::sqrt((d - rx2) * half);
      re = half * multifloats::log1p(rx2 + d + two * (rx * s1 + s2));
      im = multifloats::atan2(one + s2, rx + s1);
    }
  } else if (ix._limbs[0] < 1.0 && rx._limbs[0] < 0.5) {
    if (ix._limbs[0] >= 0x1.0p-105) {
      if (rx._limbs[0] < dd_eps_sq) {
        float64x2 onemix2 = (one + ix) * (one - ix);
        float64x2 s = multifloats::sqrt(onemix2);
        re = half * multifloats::log1p(two * rx / s);
        im = multifloats::atan2(ix, s);
      } else {
        float64x2 onemix2 = (one + ix) * (one - ix);
        float64x2 rx2 = rx * rx;
        float64x2 f = rx2 * (two + rx2 + two * ix * ix);
        float64x2 d = multifloats::sqrt(onemix2 * onemix2 + f);
        float64x2 dp = d + onemix2;
        float64x2 dm = f / dp;
        float64x2 r1 = multifloats::sqrt((dp + rx2) * half);
        float64x2 r2 = rx * ix / r1;
        re = half * multifloats::log1p(rx2 + dm + two * (rx * r1 + ix * r2));
        im = multifloats::atan2(ix + r2, rx + r1);
      }
    } else {
      float64x2 s = multifloats::hypot(one, rx);
      re = half * multifloats::log1p(two * rx * (rx + s));
      im = multifloats::atan2(ix, s);
    }
  } else {
    // General case: rx ≥ 0.5 with ix moderate. Cancellation in 1+z² and
    // log|arg| is bounded here because |arg|² is well away from 1 — the
    // textbook form carries full DD.
    float64x2 yr = (rx - ix) * (rx + ix) + one;
    float64x2 yi = two * rx * ix;
    complex<multifloats::float64x2> root =
        std::sqrt(complex<multifloats::float64x2>(yr, yi));
    complex<multifloats::float64x2> arg(rx + root.real(), ix + root.imag());
    complex<multifloats::float64x2> l = std::log(arg);
    re = l.real();
    im = l.imag();
  }

  auto dd_signbit = [](double hi, double lo) {
    return (hi == 0.0 && lo != 0.0) ? std::signbit(lo) : std::signbit(hi);
  };
  if (dd_signbit(a._limbs[0], a._limbs[1])) {
    re._limbs[0] = -re._limbs[0];
    re._limbs[1] = -re._limbs[1];
  }
  if (dd_signbit(b._limbs[0], b._limbs[1])) {
    im._limbs[0] = -im._limbs[0];
    im._limbs[1] = -im._limbs[1];
  }
  return complex<multifloats::float64x2>(re, im);
}

// acosh(z). Domain-split:
//   |z| < 1:  delegate to std::acos via `acosh(z) = i · acos(z)` (Im(z) ≥ +0)
//             or the conj-reflection `acosh(z) = conj(i · acos(conj z))`
//             (Im(z) < 0, incl. −0). Component-wise:
//               Re(acosh z) = −Im(acos z')  (small-magnitude near z=0)
//               Im(acosh z) = ±Re(acos z')  (large-magnitude; +Re for
//                                            Im≥0, −Re for Im<0)
//             The swap puts the small output on the Im(acos) path, which
//             via std::asin is full-DD. std::acos's |z| ≤ 1 arm uses the
//             direct π/2 − asin form (no recursion back here).
//   |z| ≥ 1:  Kahan's 2·log(sqrt((z+1)/2) + sqrt((z−1)/2)) — stable near
//             z ≈ 1 (avoids z²−1 cancellation) and away from the small-
//             output regime that hurts near z ≈ 0.
//         C99 G.6.2.3 branch-cut convention on Re z < 1, Im z = ±0: the
//             delegated path routes through std::acos (which handles this
//             correctly) and applies the conj-reflection based on
//             signbit(z.imag().hi), so ±0 propagates symmetrically.
template <>
complex<multifloats::float64x2>
acosh(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 zr = z.real(), zi = z.imag();
  double az2 = zr._limbs[0] * zr._limbs[0] + zi._limbs[0] * zi._limbs[0];
  if (az2 < 1.0) {
    bool use_conj = std::signbit(zi._limbs[0]);
    complex<multifloats::float64x2> zc = use_conj
        ? complex<multifloats::float64x2>(zr, float64x2(-zi._limbs[0], -zi._limbs[1]))
        : z;
    complex<multifloats::float64x2> v = std::acos(zc);
    float64x2 vre = v.real(), vim = v.imag();
    float64x2 re(-vim._limbs[0], -vim._limbs[1]);
    float64x2 im = use_conj ? float64x2(-vre._limbs[0], -vre._limbs[1]) : vre;
    return complex<multifloats::float64x2>(re, im);
  }
  float64x2 a = zr, b = zi;
  float64x2 half = float64x2(0.5, 0.0);
  complex<multifloats::float64x2> zp((a + float64x2(1.0)) * half, b * half);
  complex<multifloats::float64x2> zm((a - float64x2(1.0)) * half, b * half);
  complex<multifloats::float64x2> s1 = std::sqrt(zp);
  complex<multifloats::float64x2> s2 = std::sqrt(zm);
  complex<multifloats::float64x2> sum(s1.real() + s2.real(),
                                      s1.imag() + s2.imag());
  complex<multifloats::float64x2> l = std::log(sum);
  return complex<multifloats::float64x2>(l.real() + l.real(), l.imag() + l.imag());
}

// atanh(z). Ported from libquadmath catanhq:
//   Re = 0.25 · log((1+a)² + b²)/((1−a)² + b²))
//   Im = 0.5 · atan2(2b, 1 − a² − b²)
// The Re-part uses log1p((4a)/den) when |num/den − 1| is small so the
// small answer survives cancellation in log(num/den). The Im-part's
// denominator picks between (1−|·|)(1+|·|)−(·)² and −dd_x2y2m1 based on
// which region is near the unit circle, avoiding the direct `1 − a² − b²`
// catastrophic subtraction near |z| = 1.
template <>
complex<multifloats::float64x2>
atanh(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 a = z.real(), b = z.imag();
  float64x2 one(1.0, 0.0), half(0.5, 0.0), quarter(0.25, 0.0);
  // Branch-point singularities at z = ±1 + 0i: C99 G.6.2.4 specifies
  //   catanh(+1 + 0i) = +∞ + 0i
  //   catanh(−1 + 0i) = −∞ + 0i
  // Without the short-circuit, num/den = 4/0 = inf in DD; the subsequent
  // 0.25 * inf DD multiply emits fma(0.25, inf, -inf) = NaN, producing
  // a wrong NaN-result at this exact branch point.
  if (b._limbs[0] == 0.0 && b._limbs[1] == 0.0 &&
      std::abs(a._limbs[0]) == 1.0 && a._limbs[1] == 0.0) {
    double inf = std::numeric_limits<double>::infinity();
    return complex<multifloats::float64x2>(
        float64x2(std::copysign(inf, a._limbs[0]), 0.0),
        float64x2(0.0, 0.0));
  }

  static constexpr double dd_eps_sq  = 0x1.0p-210;
  static constexpr double dd_eps_hlf = 0x1.0p-106;

  float64x2 absa = multifloats::abs(a);
  float64x2 absb = multifloats::abs(b);
  float64x2 re_res;

  if (absa._limbs[0] == 1.0 && absa._limbs[1] == 0.0 &&
      absb._limbs[0] < dd_eps_sq) {
    // Near-branch-point: catanh(±1 + ε·i) ≈ sign(a) · 0.5 · (ln2 − log|ε|).
    float64x2 ln2(ln_2_hi, ln_2_lo);
    float64x2 t = half * (ln2 - multifloats::log(absb));
    re_res = (a._limbs[0] < 0.0) ? -t : t;
  } else {
    // i2: drop b² when |b| below eps² so the DD round-off of b² doesn't
    // contaminate num/den. Matches the libquadmath guard.
    float64x2 i2 = (absb._limbs[0] >= dd_eps_sq) ? (b * b)
                                                 : float64x2(0.0, 0.0);
    float64x2 ap1 = one + a;
    float64x2 am1 = one - a;
    float64x2 num = i2 + ap1 * ap1;
    float64x2 den = i2 + am1 * am1;
    float64x2 f = num / den;
    if (f._limbs[0] < 0.5) {
      re_res = quarter * multifloats::log(f);
    } else {
      // 4a exact in DD (power-of-two scale); computing as (num−den)/den
      // would cancel when a is small — instead feed 4a/den into log1p.
      float64x2 four_a = (a + a) + (a + a);
      re_res = quarter * multifloats::log1p(four_a / den);
    }
  }

  // Imag-part denominator picks the non-cancelling form per region.
  float64x2 big   = (absa._limbs[0] > absb._limbs[0]) ? absa : absb;
  float64x2 small_= (absa._limbs[0] > absb._limbs[0]) ? absb : absa;
  float64x2 den_im;
  if (small_._limbs[0] < dd_eps_hlf) {
    den_im = (one - big) * (one + big);
  } else if (big._limbs[0] >= 1.0) {
    den_im = (one - big) * (one + big) - small_ * small_;
  } else if (big._limbs[0] >= 0.75 || small_._limbs[0] >= 0.5) {
    // Near unit circle: 1 − big² − small² = −dd_x2y2m1(big, small).
    float64x2 t = dd_x2y2m1(big, small_);
    den_im._limbs[0] = -t._limbs[0];
    den_im._limbs[1] = -t._limbs[1];
  } else {
    den_im = (one - big) * (one + big) - small_ * small_;
  }
  float64x2 im_res = half * multifloats::atan2(b + b, den_im);

  return complex<multifloats::float64x2>(re_res, im_res);
}

// ---- Inverse trigonometric -------------------------------------------------

// asin(z) = −i · asinh(i·z). Delegating to std::asinh inherits its
// libquadmath-ported 10-region branch schedule; the old textbook
// `−i · log(iz + sqrt(1 − z²))` loses precision near the real-axis
// branch cuts because 1 − z² cancels catastrophically when |Re z| ≈ 1.
// std::asinh's signed-zero trailer already reproduces C99 G.6.2.2 for
// |Re z| > 1 ± 0: i·z swaps limbs (a → Im, −b → Re) via exact DD
// negation, so sign(Im asin) = sign(Im z) carries through symmetrically.
// std::acos = π/2 − asin inherits the fix for free.
template <>
complex<multifloats::float64x2>
asin(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 zi = z.imag();
  complex<multifloats::float64x2> iz(
      float64x2(-zi._limbs[0], -zi._limbs[1]), z.real());
  complex<multifloats::float64x2> w = std::asinh(iz);
  float64x2 wre = w.real();
  return complex<multifloats::float64x2>(
      w.imag(), float64x2(-wre._limbs[0], -wre._limbs[1]));
}

// acos(z): domain-split between two full-DD formulations.
//   |z| ≤ 1:  π/2 − asin(z). std::asin (via std::asinh) is full-DD and
//             π/2 subtracts cleanly because asin(z) is bounded by π/2 with
//             room to spare — no cancellation when |z| ≤ 1.
//   |z| > 1:  −i · acosh(z) for Im(z) ≥ +0, else the conj-reflection
//             identity i · conj(acosh(conj z)) for Im(z) < 0 (including
//             −0). See std::acosh comment for the small/large output
//             swap that avoids the π/2 − asin cancellation capping the
//             old form at DD-absolute of π/2.
template <>
complex<multifloats::float64x2>
acos(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 zr = z.real(), zi = z.imag();
  double az2 = zr._limbs[0] * zr._limbs[0] + zi._limbs[0] * zi._limbs[0];
  if (az2 > 1.0) {
    bool use_conj = std::signbit(zi._limbs[0]);
    complex<multifloats::float64x2> zc = use_conj
        ? complex<multifloats::float64x2>(zr, float64x2(-zi._limbs[0], -zi._limbs[1]))
        : z;
    complex<multifloats::float64x2> w = std::acosh(zc);
    float64x2 wre = w.real(), wim = w.imag();
    float64x2 re = wim;
    float64x2 im = use_conj ? wre : float64x2(-wre._limbs[0], -wre._limbs[1]);
    return complex<multifloats::float64x2>(re, im);
  }
  static const float64x2 half_pi = float64x2(0x1.921fb54442d18p+0,
                                              0x1.1a62633145c07p-54);
  complex<multifloats::float64x2> s = std::asin(z);
  float64x2 sim = s.imag();
  return complex<multifloats::float64x2>(
      half_pi - s.real(), float64x2(-sim._limbs[0], -sim._limbs[1]));
}

// atan(z) = −i · atanh(i·z). Delegating to std::atanh reuses its Kahan-
// style Im-part (log1p(4·a/den) gated on f ≥ 0.5) and its unit-circle
// denominator picker (dd_x2y2m1 near |z| ≈ 1). The old direct form
// `0.25·log(num/den)` lost ~7 decimals of Im precision when num ≈ den
// (|b| small, |a| moderate); the identity form inherits full DD.
// i·z = (−b + a·i), and DD negation is exact, so the rotation is free.
template <>
complex<multifloats::float64x2>
atan(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  float64x2 zi = z.imag();
  complex<multifloats::float64x2> iz(
      float64x2(-zi._limbs[0], -zi._limbs[1]), z.real());
  complex<multifloats::float64x2> w = std::atanh(iz);
  float64x2 wre = w.real();
  return complex<multifloats::float64x2>(
      w.imag(), float64x2(-wre._limbs[0], -wre._limbs[1]));
}

// ---- Value accessors -------------------------------------------------------

// |z| via the scale-aware path shared with std::log: avoids overflow when
// |z| exceeds DBL_MAX even though hypot(a,b) itself would saturate.
template <>
multifloats::float64x2
abs(complex<multifloats::float64x2> const &z) {
  return multifloats::hypot(z.real(), z.imag());
}

// arg(z) = atan2(im, re), with the standard branch cut on the negative
// real axis; signed zero in im selects ±pi.
template <>
multifloats::float64x2
arg(complex<multifloats::float64x2> const &z) {
  return multifloats::atan2(z.imag(), z.real());
}

// Riemann-sphere projection (C99 7.3.9.4): if either component is
// infinite the result is (+inf, copysign(0, im)); otherwise z is
// returned unchanged. NaN on a finite part is left alone.
template <>
complex<multifloats::float64x2>
proj(complex<multifloats::float64x2> const &z) {
  using multifloats::float64x2;
  if (std::isinf(z.real()._limbs[0]) || std::isinf(z.imag()._limbs[0])) {
    return complex<multifloats::float64x2>(
        float64x2(std::numeric_limits<double>::infinity(), 0.0),
        float64x2(std::copysign(0.0, z.imag()._limbs[0]), 0.0));
  }
  return z;
}
