/* C-ABI wrappers around the C++ multifloats inline kernels.
 *
 * On ARM64, dd_t (16 bytes) is returned in d0/d1 — no hidden pointer.
 * This lets us isolate the cost of gfortran's derived-type ABI vs the
 * platform C calling convention for the same underlying DD arithmetic.
 */
#ifndef MULTIFLOATS_C_H
#define MULTIFLOATS_C_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct { double hi, lo; } dd_t;

/* Arithmetic */
dd_t dd_add(dd_t a, dd_t b);
dd_t dd_sub(dd_t a, dd_t b);
dd_t dd_mul(dd_t a, dd_t b);
dd_t dd_div(dd_t a, dd_t b);

/* Unary */
dd_t dd_neg(dd_t a);
dd_t dd_abs(dd_t a);
dd_t dd_sqrt(dd_t a);

/* Binary */
dd_t dd_fmin(dd_t a, dd_t b);
dd_t dd_fmax(dd_t a, dd_t b);
dd_t dd_hypot(dd_t a, dd_t b);
dd_t dd_pow(dd_t a, dd_t b);
dd_t dd_fmod(dd_t a, dd_t b);
dd_t dd_fdim(dd_t a, dd_t b);
dd_t dd_copysign(dd_t a, dd_t b);
dd_t dd_fma(dd_t a, dd_t b, dd_t c);

/* Transcendental */
dd_t dd_exp(dd_t a);
dd_t dd_log(dd_t a);
dd_t dd_log10(dd_t a);
dd_t dd_sin(dd_t a);
dd_t dd_cos(dd_t a);
dd_t dd_tan(dd_t a);
dd_t dd_asin(dd_t a);
dd_t dd_acos(dd_t a);
dd_t dd_atan(dd_t a);
dd_t dd_atan2(dd_t a, dd_t b);
dd_t dd_sinh(dd_t a);
dd_t dd_cosh(dd_t a);
dd_t dd_tanh(dd_t a);
dd_t dd_asinh(dd_t a);
dd_t dd_acosh(dd_t a);
dd_t dd_atanh(dd_t a);
dd_t dd_erf(dd_t a);
dd_t dd_erfc(dd_t a);
dd_t dd_tgamma(dd_t a);
dd_t dd_lgamma(dd_t a);
dd_t dd_bessel_j0(dd_t a);
dd_t dd_bessel_j1(dd_t a);
dd_t dd_bessel_y0(dd_t a);
dd_t dd_bessel_y1(dd_t a);

/* Comparison (return int: 1 = true, 0 = false) */
int dd_eq(dd_t a, dd_t b);
int dd_ne(dd_t a, dd_t b);
int dd_lt(dd_t a, dd_t b);
int dd_le(dd_t a, dd_t b);
int dd_gt(dd_t a, dd_t b);
int dd_ge(dd_t a, dd_t b);

#ifdef __cplusplus
}
#endif

#endif /* MULTIFLOATS_C_H */
