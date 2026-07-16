// =============================================================================
// Triple-double log and atan2.
//
// Each is one Newton step from the full-DD kernel, using the existing TD
// exp / sincos. The DD seed is accurate to ~1e-32 (relative), so a single
// quadratic-convergence step lands at the TD floor (~1e-47) — enough to make
// the *absolute* error of the result drop below the DD floor, which is what a
// downstream amplifier (e.g. cpow = exp(w·log z)) needs.
//
// These live in `multifloats::detail` (scratch primitives, not public ABI).
// Must be included after exp_log / trig / inv_trig (they call log, atan2,
// exp_full_td, sincos_td).
// =============================================================================
namespace detail {

// log(x) for a triple-double x > 0, returned as triple-double.
//   exp(−y_dd) = exp(−log x + r) = (1 + r)/x  ⇒  r ≈ x·exp(−y_dd) − 1,
// so log(x) = y_dd + (x·exp(−y_dd) − 1), with exp(−y_dd) carried in TD.
inline float64x3 log_td(float64x3 const &x) {
  float64x2 y_dd = log(td_to_dd(x));
  // For huge/tiny |x| (|log x| ≳ 600) exp(∓y_dd) would over/underflow; the DD
  // log is already full-DD relative there and the refinement gains nothing.
  if (!(std::fabs(y_dd.limbs[0]) < 600.0)) return td_from_dd(y_dd);
  float64x3 e = exp_full_td(float64x2(-y_dd.limbs[0], -y_dd.limbs[1]));
  float64x3 r = td_add_double(td_mul_td(e, x), -1.0);  // x·exp(−y_dd) − 1  (small)
  return td_add_dd(r, y_dd);
}

// atan2(y, x) to triple-double. One Newton step on the angle:
//   f(θ) = x·sinθ − y·cosθ = 0 at the true angle,  f'(θ) = x·cosθ + y·sinθ,
//   θ = θ_dd − f(θ_dd)/f'(θ_dd).
// The numerator cancels to ~|z|·(seed error), so it is formed in TD; the
// correction itself is tiny, so the final divide only needs DD precision.
inline float64x3 atan2_td(float64x3 const &y, float64x3 const &x) {
  float64x2 th_dd = atan2(td_to_dd(y), td_to_dd(x));
  if (!std::isfinite(th_dd.limbs[0])) return td_from_dd(th_dd);
  float64x3 s, c;
  sincos_td(th_dd, s, c);
  float64x3 num = td_add_td(td_mul_td(x, s), td_negate(td_mul_td(y, c)));  // x·sinθ − y·cosθ
  float64x3 den = td_add_td(td_mul_td(x, c), td_mul_td(y, s));             // x·cosθ + y·sinθ
  float64x2 den_dd = td_to_dd(den);
  if (den_dd.limbs[0] == 0.0) return td_from_dd(th_dd);
  float64x2 delta = float64x2(0.0) - td_to_dd(num) / den_dd;              // −f/f'  (small)
  return td_add_dd(td_from_dd(delta), th_dd);                             // θ_dd + delta
}

}  // namespace detail
