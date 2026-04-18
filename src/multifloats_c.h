/* Double-double arithmetic — C-ABI functions.
 *
 * All DD math functions are defined as extern "C" in multifloats_math.cc.
 * They follow the math.h naming convention: dd_sin, dd_j0, dd_tgamma, etc.
 * Include this header from C, C++, or Fortran (via iso_c_binding).
 */
#ifndef MULTIFLOATS_C_H
#define MULTIFLOATS_C_H

#include <stdint.h>

/* ABI version. Bump on any breaking change to the `dd_t` layout, the
 * argument/return convention of any dd_* function, or removal of an
 * exported symbol. Additive changes (new dd_* functions) keep the same
 * version. Callers can gate on this at compile time:
 *     #if !defined(MULTIFLOATS_ABI_VERSION) || MULTIFLOATS_ABI_VERSION < 1
 *     #error "multifloats 1.x required"
 *     #endif
 */
#define MULTIFLOATS_ABI_VERSION 1

#ifdef __cplusplus
extern "C" {
#endif

/* Visibility attribute for every exported function. */
#if defined(__GNUC__) || defined(__clang__)
#  define MULTIFLOATS_API __attribute__((visibility("default")))
#else
#  define MULTIFLOATS_API
#endif

typedef struct { double hi, lo; } dd_t;

/* Guard against surprise padding — the entire C ABI and the Fortran
 * iso_c_binding layer assume dd_t is exactly two back-to-back doubles.
 * C11 / C++11 are required for this assertion. */
#ifdef __cplusplus
static_assert(sizeof(dd_t) == 2 * sizeof(double),
              "dd_t must be two back-to-back doubles with no padding");
#else
_Static_assert(sizeof(dd_t) == 2 * sizeof(double),
               "dd_t must be two back-to-back doubles with no padding");
#endif

/* Arithmetic */
MULTIFLOATS_API dd_t dd_add(dd_t a, dd_t b);
MULTIFLOATS_API dd_t dd_sub(dd_t a, dd_t b);
MULTIFLOATS_API dd_t dd_mul(dd_t a, dd_t b);
MULTIFLOATS_API dd_t dd_div(dd_t a, dd_t b);

/* Unary */
MULTIFLOATS_API dd_t dd_neg(dd_t a);
MULTIFLOATS_API dd_t dd_abs(dd_t a);
MULTIFLOATS_API dd_t dd_sqrt(dd_t a);

/* Binary */
MULTIFLOATS_API dd_t dd_fmin(dd_t a, dd_t b);
MULTIFLOATS_API dd_t dd_fmax(dd_t a, dd_t b);
MULTIFLOATS_API dd_t dd_hypot(dd_t a, dd_t b);
MULTIFLOATS_API dd_t dd_pow(dd_t a, dd_t b);
MULTIFLOATS_API dd_t dd_fmod(dd_t a, dd_t b);
MULTIFLOATS_API dd_t dd_fdim(dd_t a, dd_t b);
MULTIFLOATS_API dd_t dd_copysign(dd_t a, dd_t b);
MULTIFLOATS_API dd_t dd_fma(dd_t a, dd_t b, dd_t c);

/* Exponential / logarithmic */
MULTIFLOATS_API dd_t dd_exp(dd_t a);
MULTIFLOATS_API dd_t dd_exp2(dd_t a);
MULTIFLOATS_API dd_t dd_log(dd_t a);
MULTIFLOATS_API dd_t dd_log2(dd_t a);
MULTIFLOATS_API dd_t dd_log10(dd_t a);

/* Trigonometric */
MULTIFLOATS_API dd_t dd_sin(dd_t a);
MULTIFLOATS_API dd_t dd_cos(dd_t a);
MULTIFLOATS_API dd_t dd_tan(dd_t a);
MULTIFLOATS_API dd_t dd_asin(dd_t a);
MULTIFLOATS_API dd_t dd_acos(dd_t a);
MULTIFLOATS_API dd_t dd_atan(dd_t a);
MULTIFLOATS_API dd_t dd_atan2(dd_t a, dd_t b);

/* π-scaled trig: {sin,cos,tan}pi(x) = fn(π·x),
 *                {asin,acos,atan}pi(x) = fn(x)/π,
 *                atan2pi(y,x) = atan2(y,x)/π. */
MULTIFLOATS_API dd_t dd_sinpi(dd_t a);
MULTIFLOATS_API dd_t dd_cospi(dd_t a);
MULTIFLOATS_API dd_t dd_tanpi(dd_t a);
MULTIFLOATS_API dd_t dd_asinpi(dd_t a);
MULTIFLOATS_API dd_t dd_acospi(dd_t a);
MULTIFLOATS_API dd_t dd_atanpi(dd_t a);
MULTIFLOATS_API dd_t dd_atan2pi(dd_t a, dd_t b);

/* Hyperbolic */
MULTIFLOATS_API dd_t dd_sinh(dd_t a);
MULTIFLOATS_API dd_t dd_cosh(dd_t a);
MULTIFLOATS_API dd_t dd_tanh(dd_t a);
MULTIFLOATS_API dd_t dd_asinh(dd_t a);
MULTIFLOATS_API dd_t dd_acosh(dd_t a);
MULTIFLOATS_API dd_t dd_atanh(dd_t a);

/* Error functions. dd_erfcx is the scaled complementary error function
 *   dd_erfcx(x) = exp(x^2) * dd_erfc(x)
 * (standard name in Faddeeva/Julia/SciPy; Fortran calls it erfc_scaled). */
MULTIFLOATS_API dd_t dd_erf(dd_t a);
MULTIFLOATS_API dd_t dd_erfc(dd_t a);
MULTIFLOATS_API dd_t dd_erfcx(dd_t a);

/* Gamma functions */
MULTIFLOATS_API dd_t dd_tgamma(dd_t a);
MULTIFLOATS_API dd_t dd_lgamma(dd_t a);

/* Bessel functions (POSIX naming) */
MULTIFLOATS_API dd_t dd_j0(dd_t a);
MULTIFLOATS_API dd_t dd_j1(dd_t a);
MULTIFLOATS_API dd_t dd_y0(dd_t a);
MULTIFLOATS_API dd_t dd_y1(dd_t a);

/* Matrix multiply (column-major, Fortran layout).
 *   dd_matmul_mm: C(m,n) = A(m,k) * B(k,n)
 *   dd_matmul_mv: y(m)   = A(m,k) * x(k)
 *   dd_matmul_vm: y(n)   = x(k)   * B(k,n)
 * Leading dimensions equal the first extent (no strides).
 *
 * renorm_interval: if > 0, renormalize accumulators every N reductions
 * (matches MF_FMA_RENORM_INTERVAL in the Fortran layer — keeps s_lo
 * bounded for large k). Pass 0 to renormalize only at the end. */
MULTIFLOATS_API void dd_matmul_mm(const dd_t *a, const dd_t *b, dd_t *c,
                         int64_t m, int64_t k, int64_t n,
                         int64_t renorm_interval);
MULTIFLOATS_API void dd_matmul_mv(const dd_t *a, const dd_t *x, dd_t *y,
                         int64_t m, int64_t k,
                         int64_t renorm_interval);
MULTIFLOATS_API void dd_matmul_vm(const dd_t *x, const dd_t *b, dd_t *y,
                         int64_t k, int64_t n,
                         int64_t renorm_interval);

/* Comparison (return int: 1 = true, 0 = false) */
MULTIFLOATS_API int dd_eq(dd_t a, dd_t b);
MULTIFLOATS_API int dd_ne(dd_t a, dd_t b);
MULTIFLOATS_API int dd_lt(dd_t a, dd_t b);
MULTIFLOATS_API int dd_le(dd_t a, dd_t b);
MULTIFLOATS_API int dd_gt(dd_t a, dd_t b);
MULTIFLOATS_API int dd_ge(dd_t a, dd_t b);

#ifdef __cplusplus
}
#endif

/* Scope the visibility macro to this header so callers don't get a leaked
 * `MULTIFLOATS_API` identifier. Consumers wanting their own attribute
 * macros should define them in their own namespace. */
#undef MULTIFLOATS_API

#endif /* MULTIFLOATS_C_H */
