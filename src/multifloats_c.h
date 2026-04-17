/* Double-double arithmetic — C-ABI functions.
 *
 * All DD math functions are defined as extern "C" in multifloats_math.cc.
 * They follow the math.h naming convention: dd_sin, dd_j0, dd_tgamma, etc.
 * Include this header from C, C++, or Fortran (via iso_c_binding).
 */
#ifndef MULTIFLOATS_C_H
#define MULTIFLOATS_C_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Visibility: extern "C" functions are the only exported symbols. */
#if defined(__GNUC__) || defined(__clang__)
#  define DD_API __attribute__((visibility("default")))
#else
#  define DD_API
#endif

typedef struct { double hi, lo; } dd_t;

/* Arithmetic */
DD_API dd_t dd_add(dd_t a, dd_t b);
DD_API dd_t dd_sub(dd_t a, dd_t b);
DD_API dd_t dd_mul(dd_t a, dd_t b);
DD_API dd_t dd_div(dd_t a, dd_t b);

/* Unary */
DD_API dd_t dd_neg(dd_t a);
DD_API dd_t dd_abs(dd_t a);
DD_API dd_t dd_sqrt(dd_t a);

/* Binary */
DD_API dd_t dd_fmin(dd_t a, dd_t b);
DD_API dd_t dd_fmax(dd_t a, dd_t b);
DD_API dd_t dd_hypot(dd_t a, dd_t b);
DD_API dd_t dd_pow(dd_t a, dd_t b);
DD_API dd_t dd_fmod(dd_t a, dd_t b);
DD_API dd_t dd_fdim(dd_t a, dd_t b);
DD_API dd_t dd_copysign(dd_t a, dd_t b);
DD_API dd_t dd_fma(dd_t a, dd_t b, dd_t c);

/* Exponential / logarithmic */
DD_API dd_t dd_exp(dd_t a);
DD_API dd_t dd_exp2(dd_t a);
DD_API dd_t dd_log(dd_t a);
DD_API dd_t dd_log2(dd_t a);
DD_API dd_t dd_log10(dd_t a);

/* Trigonometric */
DD_API dd_t dd_sin(dd_t a);
DD_API dd_t dd_cos(dd_t a);
DD_API dd_t dd_tan(dd_t a);
DD_API dd_t dd_asin(dd_t a);
DD_API dd_t dd_acos(dd_t a);
DD_API dd_t dd_atan(dd_t a);
DD_API dd_t dd_atan2(dd_t a, dd_t b);

/* Hyperbolic */
DD_API dd_t dd_sinh(dd_t a);
DD_API dd_t dd_cosh(dd_t a);
DD_API dd_t dd_tanh(dd_t a);
DD_API dd_t dd_asinh(dd_t a);
DD_API dd_t dd_acosh(dd_t a);
DD_API dd_t dd_atanh(dd_t a);

/* Error functions */
DD_API dd_t dd_erf(dd_t a);
DD_API dd_t dd_erfc(dd_t a);

/* Gamma functions */
DD_API dd_t dd_tgamma(dd_t a);
DD_API dd_t dd_lgamma(dd_t a);

/* Bessel functions (POSIX naming) */
DD_API dd_t dd_j0(dd_t a);
DD_API dd_t dd_j1(dd_t a);
DD_API dd_t dd_y0(dd_t a);
DD_API dd_t dd_y1(dd_t a);

/* Matrix multiply (column-major, Fortran layout).
 *   dd_matmul_mm: C(m,n) = A(m,k) * B(k,n)
 *   dd_matmul_mv: y(m)   = A(m,k) * x(k)
 *   dd_matmul_vm: y(n)   = x(k)   * B(k,n)
 * Leading dimensions equal the first extent (no strides).
 *
 * renorm_interval: if > 0, renormalize accumulators every N reductions
 * (matches MF_FMA_RENORM_INTERVAL in the Fortran layer — keeps s_lo
 * bounded for large k). Pass 0 to renormalize only at the end. */
DD_API void dd_matmul_mm(const dd_t *a, const dd_t *b, dd_t *c,
                         int64_t m, int64_t k, int64_t n,
                         int64_t renorm_interval);
DD_API void dd_matmul_mv(const dd_t *a, const dd_t *x, dd_t *y,
                         int64_t m, int64_t k,
                         int64_t renorm_interval);
DD_API void dd_matmul_vm(const dd_t *x, const dd_t *b, dd_t *y,
                         int64_t k, int64_t n,
                         int64_t renorm_interval);

/* Comparison (return int: 1 = true, 0 = false) */
DD_API int dd_eq(dd_t a, dd_t b);
DD_API int dd_ne(dd_t a, dd_t b);
DD_API int dd_lt(dd_t a, dd_t b);
DD_API int dd_le(dd_t a, dd_t b);
DD_API int dd_gt(dd_t a, dd_t b);
DD_API int dd_ge(dd_t a, dd_t b);

#ifdef __cplusplus
}
#endif

#endif /* MULTIFLOATS_C_H */
