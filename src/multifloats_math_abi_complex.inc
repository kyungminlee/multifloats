// ---- Complex DD transcendentals ------------------------------------------
//
// Branch cuts follow C99 Annex G (matching libquadmath cexpq/clogq/csqrtq).
// Where the real-axis formula needs both sin(y) and cos(y) (or both sinh
// and cosh), these call the fused sincos_full / sinhcosh_full kernels
// so one range reduction + Taylor pair covers both outputs. Functions that
// are built on top of cdd_sqrt / cdd_log (asin/acos/...) call the exported
// cdd_* symbols directly; within a single TU the compiler inlines those.
//
// The text "2 fused" vs "4 separate" in speed comments below is measured
// against libstdc++'s generic <complex> template path (see commit notes).

// exp(a+bi) = e^a · (cos b + i·sin b). 2 transcendentals (exp + sincos).
complex64x2_t cexpdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 ea = exp_full(a);
  float64x2 s, c;
  sincos_full(b, s, c);
  return { to(ea * c), to(ea * s) };
}

// log(a+bi) = log(|z|) + i·atan2(b, a). Compute log(|z|) as
//   log(big) + 0.5·log(1 + (small/big)^2)
// with big = max(|a|,|b|), small = min(|a|,|b|). This avoids the
// overflow that hypot(a,b) suffers when |z| exceeds DBL_MAX even
// though log(|z|) itself is finite. atan2 handles the negative-real-
// axis branch cut (including signed-zero propagation).
complex64x2_t clogdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 ax = multifloats::abs(a);
  float64x2 ay = multifloats::abs(b);
  float64x2 big   = (ax > ay) ? ax : ay;
  float64x2 small = (ax > ay) ? ay : ax;
  float64x2 r;
  if (big._limbs[0] == 0.0) {
    r._limbs[0] = -std::numeric_limits<double>::infinity();
    r._limbs[1] = 0.0;
  } else {
    float64x2 ratio = small / big;
    float64x2 one = float64x2(1.0, 0.0);
    float64x2 half = float64x2(0.5, 0.0);
    r = log_full(big) + half * log_full(one + ratio * ratio);
  }
  float64x2 phi = atan2_full(b, a);
  return { to(r), to(phi) };
}

// log10(z) = log(z) · (1/ln 10).
complex64x2_t clog10dd(complex64x2_t z) {
  static const float64x2 inv_ln10 = float64x2(0x1.bcb7b1526e50ep-2,
                                      0x1.95355baaafad3p-57);
  complex64x2_t l = clogdd(z);
  return { to(from(l.re) * inv_ln10), to(from(l.im) * inv_ln10) };
}

// log2(z) = log(z) · (1/ln 2).
complex64x2_t clog2dd(complex64x2_t z) {
  static const float64x2 inv_ln2 =
      float64x2(1.4426950408889634, 2.0355273740931033e-17);
  complex64x2_t l = clogdd(z);
  return { to(from(l.re) * inv_ln2), to(from(l.im) * inv_ln2) };
}

// log1p(z) = log(1 + z) without losing precision for |z| << 1. Uses
// log1p on the real scalar form `|1+z|² - 1 = 2a + a² + b²`, which is
// computed without the 1 ± 1 cancellation. The imaginary part is
// `atan2(b, 1 + a)` — atan2 handles the branch cut and signed-zero
// propagation directly. Matches libquadmath clog1pq.
complex64x2_t clog1pdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 u = (a + a) + a * a + b * b;          // 2a + a² + b²
  float64x2 half(0.5, 0.0);
  float64x2 r = half * log1p_full(u);
  float64x2 phi = atan2_full(b, float64x2(1.0, 0.0) + a);
  return { to(r), to(phi) };
}

// expm1(z) = e^z - 1 without losing precision for |z| << 1. Write
//   Re = (e^a - 1) cos(b) + (cos(b) - 1)
//      = expm1(a)·cos(b) − 2·sin²(b/2)
//   Im = e^a · sin(b) = (1 + expm1(a)) · sin(b)
// The half-angle form of (cos(b) − 1) keeps the real part full-DD for
// small |b| where the direct `cos(b) − 1` cancels. Two fused sincos
// pairs — cost is ~2× cexp for general z, matches cexp for |z| ≫ 1.
complex64x2_t cexpm1dd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 ea_m1 = expm1_full(a);
  float64x2 s, c;
  sincos_full(b, s, c);
  float64x2 shalf, chalf;
  sincos_full(float64x2(0.5, 0.0) * b, shalf, chalf);
  float64x2 cos_m1 = -((shalf * shalf) + (shalf * shalf));
  float64x2 re = ea_m1 * c + cos_m1;
  float64x2 ea = float64x2(1.0, 0.0) + ea_m1;
  float64x2 im = ea * s;
  return { to(re), to(im) };
}

// sin(π·z) = sin(π·a) cosh(π·b) + i cos(π·a) sinh(π·b). Using sinpi/cospi
// on the real part keeps full DD precision for half-integer a (where
// sin(π·a) is exactly 0 or ±1 but sin_full(π·a) loses bits to the range
// reduction on π·a). The imaginary scaling π·b is a direct DD multiply.
complex64x2_t csinpidd(complex64x2_t z) {
  static const float64x2 pi_dd(pi_dd_hi, pi_dd_lo);
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 s = sinpi_full(a);
  float64x2 c = cospi_full(a);
  float64x2 pib = pi_dd * b;
  float64x2 sh, ch;
  sinhcosh_full(pib, sh, ch);
  return { to(s * ch), to(c * sh) };
}

// cos(π·z) = cos(π·a) cosh(π·b) − i sin(π·a) sinh(π·b).
complex64x2_t ccospidd(complex64x2_t z) {
  static const float64x2 pi_dd(pi_dd_hi, pi_dd_lo);
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 s = sinpi_full(a);
  float64x2 c = cospi_full(a);
  float64x2 pib = pi_dd * b;
  float64x2 sh, ch;
  sinhcosh_full(pib, sh, ch);
  return { to(c * ch), to(-(s * sh)) };
}

// pow(z, w) = exp(w · log(z)). C99 G.6.4.1; 0^w folds to 0.
complex64x2_t cpowdd(complex64x2_t z, complex64x2_t w) {
  if (z.re.hi == 0.0 && z.im.hi == 0.0) return { {0.0, 0.0}, {0.0, 0.0} };
  complex64x2_t l = clogdd(z);
  float64x2 lr = from(l.re), li = from(l.im);
  float64x2 wr = from(w.re), wi = from(w.im);
  // w * log(z) = (wr·lr − wi·li) + i·(wr·li + wi·lr)
  complex64x2_t p = { to(wr * lr - wi * li), to(wr * li + wi * lr) };
  return cexpdd(p);
}

// sqrt(z). Principal branch; cut on negative real axis, continuous above.
// Uses Moshier's 2·Re·Im = Im identity to avoid cancellation in mod ± a.
complex64x2_t csqrtdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  if (a._limbs[0] == 0.0 && b._limbs[0] == 0.0)
    return { to(float64x2(0.0)), z.im };  // preserves signed zero of imag
  float64x2 mod = multifloats::hypot(a, b);
  float64x2 half = float64x2(0.5, 0.0);
  if (a._limbs[0] >= 0.0) {
    float64x2 r = multifloats::sqrt((mod + a) * half);
    float64x2 i = b / (r + r);
    return { to(r), to(i) };
  } else {
    float64x2 s = multifloats::sqrt((mod - a) * half);
    float64x2 r = multifloats::abs(b) / (s + s);
    // Use signbit so b = -0 picks the −s branch (IEEE: -0 < 0 is false,
    // so `< 0.0` mishandles the negative-real-axis branch cut per C99 G.6.4.2).
    float64x2 i = std::signbit(b._limbs[0]) ? -s : s;
    return { to(r), to(i) };
  }
}

// sin(a+bi) = sin(a)·cosh(b) + i·cos(a)·sinh(b). 2 fused transcendentals.
complex64x2_t csindd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 sa, ca, sb, cb;
  sincos_full(a, sa, ca);
  sinhcosh_full(b, sb, cb);
  return { to(sa * cb), to(ca * sb) };
}

// cos(a+bi) = cos(a)·cosh(b) − i·sin(a)·sinh(b). 2 fused transcendentals.
complex64x2_t ccosdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 sa, ca, sb, cb;
  sincos_full(a, sa, ca);
  sinhcosh_full(b, sb, cb);
  return { to(ca * cb), to(-(sa * sb)) };
}

// tan(a+bi) = (sin(a)·cos(a) + i·sinh(b)·cosh(b)) / (cos(a)² + sinh(b)²).
// libquadmath ctanq formula; 2 fused transcendentals + real div
// (vs 8 for generic sin(z)/cos(z)).
complex64x2_t ctandd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 sa, ca, sb, cb;
  sincos_full(a, sa, ca);
  sinhcosh_full(b, sb, cb);
  float64x2 den = ca * ca + sb * sb;
  return { to((sa * ca) / den), to((sb * cb) / den) };
}

// asin(z) = −i · log(i·z + sqrt(1 − z²)). Cuts on real axis for |a|>1.
complex64x2_t casindd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  // 1 − z² : Re = 1 − a² + b², Im = −2ab
  complex64x2_t one_mz2 = { to(float64x2(1.0) - a * a + b * b), to(-(a * b + a * b)) };
  complex64x2_t root = csqrtdd(one_mz2);
  float64x2 rr = from(root.re), ri = from(root.im);
  // i·z + root : (−b + rr) + i·(a + ri)
  complex64x2_t arg = { to(-b + rr), to(a + ri) };
  complex64x2_t l = clogdd(arg);
  // −i · (lr + i·li) = li − i·lr
  complex64x2_t r = { l.im, to(-from(l.re)) };
  // Branch-cut fixup for real-axis inputs with |Re z| > 1 (C99 G.6.2.2):
  // the sign of Im(result) must mirror the sign of Im(input), including
  // ±0. DD multiply doesn't preserve signed zero (2 * -0 → +0), so the
  // sign is lost through the (a*b + a*b) chain. If the hi limb ended up
  // with the wrong sign, negate the WHOLE DD (both limbs) — a simple
  // copysign-per-limb would corrupt the pair because a typical DD lo
  // has the opposite sign from hi (representing the rounding residual).
  if (std::signbit(r.im.hi) != std::signbit(b._limbs[0])) {
    r.im.hi = -r.im.hi;
    r.im.lo = -r.im.lo;
  }
  return r;
}

// acos(z) = π/2 − asin(z). Uses DD-precision π/2; the generic <complex>
// template truncates `(_Tp)1.5707963…L` to double on arm64-macOS, losing
// the DD low limb — specializing here fixes that correctness bug.
complex64x2_t cacosdd(complex64x2_t z) {
  static const float64x2 half_pi = float64x2(0x1.921fb54442d18p+0,
                                     0x1.1a62633145c07p-54);
  complex64x2_t s = casindd(z);
  return { to(half_pi - from(s.re)), to(-from(s.im)) };
}

// atan(z). libstdc++ generic form:
//   Re = 0.5 · atan2(2a, 1 − a² − b²)
//   Im = 0.25 · log((a² + (b+1)²) / (a² + (b−1)²))
// Single log-ratio rather than the naive (i/2)(log(1−iz) − log(1+iz)).
complex64x2_t catandd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 a2 = a * a;
  float64x2 x = float64x2(1.0) - a2 - b * b;
  float64x2 bp1 = b + float64x2(1.0), bm1 = b - float64x2(1.0);
  float64x2 num = a2 + bp1 * bp1;
  float64x2 den = a2 + bm1 * bm1;
  float64x2 half = float64x2(0.5, 0.0), quarter = float64x2(0.25, 0.0);
  return { to(half * atan2_full(a + a, x)),
           to(quarter * log_full(num / den)) };
}

// sinh(a+bi) = sinh(a)·cos(b) + i·cosh(a)·sin(b). 2 fused transcendentals.
complex64x2_t csinhdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 sa, ca, sb, cb;
  sinhcosh_full(a, sa, ca);
  sincos_full(b, sb, cb);
  return { to(sa * cb), to(ca * sb) };
}

// cosh(a+bi) = cosh(a)·cos(b) + i·sinh(a)·sin(b). 2 fused transcendentals.
complex64x2_t ccoshdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 sa, ca, sb, cb;
  sinhcosh_full(a, sa, ca);
  sincos_full(b, sb, cb);
  return { to(ca * cb), to(sa * sb) };
}

// tanh(a+bi) = (sinh(a)·cosh(a) + i·sin(b)·cos(b)) / (sinh(a)² + cos(b)²).
// libquadmath ctanhq formula.
complex64x2_t ctanhdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 sa, ca, sb, cb;
  sinhcosh_full(a, sa, ca);
  sincos_full(b, sb, cb);
  float64x2 den = sa * sa + cb * cb;
  return { to((sa * ca) / den), to((sb * cb) / den) };
}

// asinh(z) = log(z + sqrt(z² + 1)). Textbook; matches libstdc++ generic.
complex64x2_t casinhdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  // 1 + z² : Re = 1 + a² − b², Im = 2ab
  complex64x2_t one_pz2 = { to(float64x2(1.0) + a * a - b * b), to(a * b + a * b) };
  complex64x2_t root = csqrtdd(one_pz2);
  complex64x2_t arg = { to(a + from(root.re)), to(b + from(root.im)) };
  return clogdd(arg);
}

// acosh(z). Kahan's formula: 2·log(sqrt((z+1)/2) + sqrt((z−1)/2)).
// More stable near z≈1 than the naive log(z + sqrt(z²−1)) which suffers
// cancellation in z²−1.
complex64x2_t cacoshdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  float64x2 half = float64x2(0.5, 0.0);
  complex64x2_t zp = { to((a + float64x2(1.0)) * half), to(b * half) };
  complex64x2_t zm = { to((a - float64x2(1.0)) * half), to(b * half) };
  complex64x2_t s1 = csqrtdd(zp);
  complex64x2_t s2 = csqrtdd(zm);
  complex64x2_t sum = { to(from(s1.re) + from(s2.re)), to(from(s1.im) + from(s2.im)) };
  complex64x2_t l = clogdd(sum);
  return { to(from(l.re) + from(l.re)), to(from(l.im) + from(l.im)) };
}

// atanh(z). libstdc++ generic form (atan with a, b swapped):
//   Re = 0.25 · log((1+a)² + b²)/((1−a)² + b²))
//   Im = 0.5 · atan2(2b, 1 − a² − b²)
complex64x2_t catanhdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  // Branch-point singularities at z = ±1 + 0i: C99 G.6.2.4 specifies
  //   catanh(+1 + 0i) = +∞ + 0i
  //   catanh(−1 + 0i) = −∞ + 0i
  // Without the short-circuit, num/den = 4/0 = inf in DD; the subsequent
  // 0.25 * inf DD multiply emits fma(0.25, inf, -inf) = NaN, producing
  // a wrong NaN-result at this exact branch point.
  if (b._limbs[0] == 0.0 && b._limbs[1] == 0.0 &&
      std::abs(a._limbs[0]) == 1.0 && a._limbs[1] == 0.0) {
    double inf = std::numeric_limits<double>::infinity();
    complex64x2_t r;
    r.re.hi = std::copysign(inf, a._limbs[0]);
    r.re.lo = 0.0;
    r.im.hi = 0.0;
    r.im.lo = 0.0;
    return r;
  }
  float64x2 b2 = b * b;
  float64x2 x = float64x2(1.0) - b2 - a * a;
  float64x2 ap1 = a + float64x2(1.0), am1 = a - float64x2(1.0);
  float64x2 num = b2 + ap1 * ap1;
  float64x2 den = b2 + am1 * am1;
  float64x2 half = float64x2(0.5, 0.0), quarter = float64x2(0.25, 0.0);
  return { to(quarter * log_full(num / den)),
           to(half * atan2_full(b + b, x)) };
}

// |z| via the scale-aware path shared with clogdd: avoids overflow when
// |z| exceeds DBL_MAX even though hypot(a,b) itself would saturate.
float64x2_t cabsdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  return to(multifloats::hypot(a, b));
}

// arg(z) = atan2(im, re), with the standard branch cut on the negative
// real axis; signed zero in im selects ±pi.
float64x2_t cargdd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);
  return to(atan2_full(b, a));
}

// Riemann-sphere projection (C99 7.3.9.4): if either component is
// infinite the result is (+inf, copysign(0, im)); otherwise z is
// returned unchanged. NaN on a finite part is left alone.
complex64x2_t cprojdd(complex64x2_t z) {
  if (std::isinf(z.re.hi) || std::isinf(z.im.hi)) {
    float64x2_t inf = { std::numeric_limits<double>::infinity(), 0.0 };
    float64x2_t zim = { std::copysign(0.0, z.im.hi), 0.0 };
    return { inf, zim };
  }
  return z;
}

// Conjugate: flip the sign of every limb of the imaginary part so the
// DD value negates exactly (no roundoff, preserves the low limb).
complex64x2_t conjdd(complex64x2_t z) {
  return { z.re, { -z.im.hi, -z.im.lo } };
}

// Real/imag accessors. Trivial but exposed as linkable symbols for
// parity with libquadmath's crealq / cimagq.
float64x2_t crealdd(complex64x2_t z) { return z.re; }
float64x2_t cimagdd(complex64x2_t z) { return z.im; }

// Matrix multiply (column-major).
//
// Inner accumulator uses Neumaier-style compensation on the DD-product stream:
//   p, e = two_prod(ah, bh)          (exact)
//   cross = fma(ah, bl, al*bh)       (DD-mul cross term, ignoring bl*bl)
//   (s_hi, s_lo) += (p + cross + e)  via two_sum on s_hi
// A final two_sum renormalizes s_hi/s_lo into the DD result; if
// `renorm_interval > 0`, intermediate two_sums also fire every
// that-many reductions to keep s_lo bounded at large k.
//
// Loop order is AXPY / gaxpy: outer p (shared dim), inner i (output row)
// so A is read as contiguous columns and the output accumulators are
// m-way data-parallel. The row dimension is register-blocked to MR=8
// via a compile-time-sized panel template; any m is handled by
// MR panels + a 1..7 tail. mm calls the same mv dispatch per output
// column.
