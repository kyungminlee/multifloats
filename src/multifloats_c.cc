// C-ABI wrappers around the C++ multifloats inline kernels.
// Each wrapper converts dd_t ↔ MultiFloat<double,2> and calls the
// inlined C++ operator / free function. The compiler inlines the C++
// body into the wrapper, so the only overhead vs a pure C++ call is
// the extern "C" ABI contract (which on ARM64 returns dd_t in d0/d1).

#include "multifloats.hh"
#include "multifloats_c.h"

namespace mf = multifloats;
using MF2 = mf::MultiFloat<double, 2>;

static inline MF2 from(dd_t x) {
  MF2 r;
  r._limbs[0] = x.hi;
  r._limbs[1] = x.lo;
  return r;
}

static inline dd_t to(MF2 const &x) {
  return {x._limbs[0], x._limbs[1]};
}

extern "C" {

dd_t dd_add(dd_t a, dd_t b) { return to(from(a) + from(b)); }
dd_t dd_sub(dd_t a, dd_t b) { return to(from(a) - from(b)); }
dd_t dd_mul(dd_t a, dd_t b) { return to(from(a) * from(b)); }
dd_t dd_div(dd_t a, dd_t b) { return to(from(a) / from(b)); }

dd_t dd_neg(dd_t a)  { return to(-from(a)); }
dd_t dd_abs(dd_t a)  { return to(mf::abs(from(a))); }
dd_t dd_sqrt(dd_t a) { return to(mf::sqrt(from(a))); }

dd_t dd_fmin(dd_t a, dd_t b)     { return to(mf::fmin(from(a), from(b))); }
dd_t dd_fmax(dd_t a, dd_t b)     { return to(mf::fmax(from(a), from(b))); }
dd_t dd_hypot(dd_t a, dd_t b)    { return to(mf::hypot(from(a), from(b))); }
dd_t dd_pow(dd_t a, dd_t b)      { return to(mf::pow(from(a), from(b))); }
dd_t dd_fmod(dd_t a, dd_t b)     { return to(mf::fmod(from(a), from(b))); }
dd_t dd_fdim(dd_t a, dd_t b)     { return to(mf::fdim(from(a), from(b))); }
dd_t dd_copysign(dd_t a, dd_t b) { return to(mf::copysign(from(a), from(b))); }
dd_t dd_fma(dd_t a, dd_t b, dd_t c) { return to(mf::fma(from(a), from(b), from(c))); }

dd_t dd_exp(dd_t a)   { return to(mf::exp(from(a))); }
dd_t dd_log(dd_t a)   { return to(mf::log(from(a))); }
dd_t dd_log10(dd_t a) { return to(mf::log10(from(a))); }
dd_t dd_sin(dd_t a)   { return to(mf::sin(from(a))); }
dd_t dd_cos(dd_t a)   { return to(mf::cos(from(a))); }
dd_t dd_tan(dd_t a)   { return to(mf::tan(from(a))); }
dd_t dd_asin(dd_t a)  { return to(mf::asin(from(a))); }
dd_t dd_acos(dd_t a)  { return to(mf::acos(from(a))); }
dd_t dd_atan(dd_t a)  { return to(mf::atan(from(a))); }
dd_t dd_atan2(dd_t a, dd_t b) { return to(mf::atan2(from(a), from(b))); }
dd_t dd_sinh(dd_t a)  { return to(mf::sinh(from(a))); }
dd_t dd_cosh(dd_t a)  { return to(mf::cosh(from(a))); }
dd_t dd_tanh(dd_t a)  { return to(mf::tanh(from(a))); }
dd_t dd_asinh(dd_t a) { return to(mf::asinh(from(a))); }
dd_t dd_acosh(dd_t a) { return to(mf::acosh(from(a))); }
dd_t dd_atanh(dd_t a) { return to(mf::atanh(from(a))); }
dd_t dd_erf(dd_t a)   { return to(mf::erf(from(a))); }
dd_t dd_erfc(dd_t a)  { return to(mf::erfc(from(a))); }
dd_t dd_tgamma(dd_t a) { return to(mf::tgamma(from(a))); }
dd_t dd_lgamma(dd_t a) { return to(mf::lgamma(from(a))); }
dd_t dd_bessel_j0(dd_t a) { return to(mf::bessel_j0(from(a))); }
dd_t dd_bessel_j1(dd_t a) { return to(mf::bessel_j1(from(a))); }
dd_t dd_bessel_y0(dd_t a) { return to(mf::bessel_y0(from(a))); }
dd_t dd_bessel_y1(dd_t a) { return to(mf::bessel_y1(from(a))); }

int dd_eq(dd_t a, dd_t b) { return from(a) == from(b) ? 1 : 0; }
int dd_ne(dd_t a, dd_t b) { return from(a) != from(b) ? 1 : 0; }
int dd_lt(dd_t a, dd_t b) { return from(a) <  from(b) ? 1 : 0; }
int dd_le(dd_t a, dd_t b) { return from(a) <= from(b) ? 1 : 0; }
int dd_gt(dd_t a, dd_t b) { return from(a) >  from(b) ? 1 : 0; }
int dd_ge(dd_t a, dd_t b) { return from(a) >= from(b) ? 1 : 0; }

}
