/* Double-double arithmetic — C-ABI functions.
 *
 * All DD math functions are defined as extern "C" in multifloats_math.cc.
 * They follow the libc / libquadmath convention: type suffix on the
 * function name (sindd, logdd, cexpdd) — mirroring `sinf`/`sinq`, `cexpf`/
 * `csinq`. Include this header from C, C++, or Fortran (via iso_c_binding).
 */
#pragma once

#include <stdint.h>

/* ABI version. Bump on any breaking change to the `float64x2_t` layout,
 * the argument/return convention of any exported function, or removal of
 * an exported symbol. Additive changes (new *dd functions) keep the same
 * version. Callers can gate on this at compile time:
 *     #if !defined(MULTIFLOATS_ABI_VERSION) || MULTIFLOATS_ABI_VERSION < 2
 *     #error "multifloats 2.x required"
 *     #endif
 */
#define MULTIFLOATS_ABI_VERSION 2

#ifdef __cplusplus
extern "C" {
#endif

/* Visibility attribute for every exported function. */
#if defined(__GNUC__) || defined(__clang__)
#  define MULTIFLOATS_API __attribute__((visibility("default")))
#else
#  define MULTIFLOATS_API
#endif

typedef struct { double hi, lo; } float64x2_t;

/* Complex double-double: two back-to-back float64x2_t (four doubles, no
 * padding). Matches Fortran's `type(complex64x2) sequence` layout and is
 * used by the c*dd transcendentals below. */
typedef struct { float64x2_t re, im; } complex64x2_t;

/* Guard against surprise padding — the entire C ABI and the Fortran
 * iso_c_binding layer assume float64x2_t is exactly two back-to-back
 * doubles (and complex64x2_t exactly four). C11 / C++11 required. */
#ifdef __cplusplus
static_assert(sizeof(float64x2_t) == 2 * sizeof(double),
              "float64x2_t must be two back-to-back doubles with no padding");
static_assert(sizeof(complex64x2_t) == 4 * sizeof(double),
              "complex64x2_t must be four back-to-back doubles with no padding");
#else
_Static_assert(sizeof(float64x2_t) == 2 * sizeof(double),
               "float64x2_t must be two back-to-back doubles with no padding");
_Static_assert(sizeof(complex64x2_t) == 4 * sizeof(double),
               "complex64x2_t must be four back-to-back doubles with no padding");
#endif

/* Arithmetic */
MULTIFLOATS_API float64x2_t adddd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API float64x2_t subdd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API float64x2_t muldd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API float64x2_t divdd(float64x2_t a, float64x2_t b);

/* Unary */
MULTIFLOATS_API float64x2_t negdd(float64x2_t a);
MULTIFLOATS_API float64x2_t fabsdd(float64x2_t a);
MULTIFLOATS_API float64x2_t sqrtdd(float64x2_t a);

/* Rounding. truncdd matches C trunc / Fortran AINT (toward zero).
 * rounddd matches C round / Fortran ANINT (to nearest, halfway away). */
MULTIFLOATS_API float64x2_t truncdd(float64x2_t a);
MULTIFLOATS_API float64x2_t rounddd(float64x2_t a);

/* Binary */
MULTIFLOATS_API float64x2_t fmindd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API float64x2_t fmaxdd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API float64x2_t hypotdd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API float64x2_t powdd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API float64x2_t fmoddd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API float64x2_t fdimdd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API float64x2_t copysigndd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API float64x2_t fmadd(float64x2_t a, float64x2_t b, float64x2_t c);

/* Exponential / logarithmic. `expm1dd` / `log1pdd` are the
 * cancellation-safe variants for arguments near zero — direct Taylor
 * / atanh-narrow kernels cover |x| below a threshold, larger |x|
 * falls through to the standard `exp` / `log` paths. */
MULTIFLOATS_API float64x2_t expdd(float64x2_t a);
MULTIFLOATS_API float64x2_t exp2dd(float64x2_t a);
MULTIFLOATS_API float64x2_t expm1dd(float64x2_t a);
MULTIFLOATS_API float64x2_t logdd(float64x2_t a);
MULTIFLOATS_API float64x2_t log2dd(float64x2_t a);
MULTIFLOATS_API float64x2_t log10dd(float64x2_t a);
MULTIFLOATS_API float64x2_t log1pdd(float64x2_t a);

/* Trigonometric */
MULTIFLOATS_API float64x2_t sindd(float64x2_t a);
MULTIFLOATS_API float64x2_t cosdd(float64x2_t a);
MULTIFLOATS_API float64x2_t tandd(float64x2_t a);
MULTIFLOATS_API float64x2_t asindd(float64x2_t a);
MULTIFLOATS_API float64x2_t acosdd(float64x2_t a);
MULTIFLOATS_API float64x2_t atandd(float64x2_t a);
MULTIFLOATS_API float64x2_t atan2dd(float64x2_t a, float64x2_t b);

/* π-scaled trig: {sin,cos,tan}pidd(x) = fn(π·x),
 *                {asin,acos,atan}pidd(x) = fn(x)/π,
 *                atan2pidd(y,x) = atan2(y,x)/π. */
MULTIFLOATS_API float64x2_t sinpidd(float64x2_t a);
MULTIFLOATS_API float64x2_t cospidd(float64x2_t a);
MULTIFLOATS_API float64x2_t tanpidd(float64x2_t a);
MULTIFLOATS_API float64x2_t asinpidd(float64x2_t a);
MULTIFLOATS_API float64x2_t acospidd(float64x2_t a);
MULTIFLOATS_API float64x2_t atanpidd(float64x2_t a);
MULTIFLOATS_API float64x2_t atan2pidd(float64x2_t a, float64x2_t b);

/* Hyperbolic */
MULTIFLOATS_API float64x2_t sinhdd(float64x2_t a);
MULTIFLOATS_API float64x2_t coshdd(float64x2_t a);
MULTIFLOATS_API float64x2_t tanhdd(float64x2_t a);
MULTIFLOATS_API float64x2_t asinhdd(float64x2_t a);
MULTIFLOATS_API float64x2_t acoshdd(float64x2_t a);
MULTIFLOATS_API float64x2_t atanhdd(float64x2_t a);

/* Error functions. erfcxdd is the scaled complementary error function
 *   erfcxdd(x) = exp(x^2) * erfcdd(x)
 * (standard name in Faddeeva/Julia/SciPy; Fortran calls it erfc_scaled). */
MULTIFLOATS_API float64x2_t erfdd(float64x2_t a);
MULTIFLOATS_API float64x2_t erfcdd(float64x2_t a);
MULTIFLOATS_API float64x2_t erfcxdd(float64x2_t a);

/* Gamma functions */
MULTIFLOATS_API float64x2_t tgammadd(float64x2_t a);
MULTIFLOATS_API float64x2_t lgammadd(float64x2_t a);

/* Bessel functions (POSIX naming). jndd / yndd: integer-order recurrence
 * seeded from j0dd / j1dd (Miller's backward recurrence for Jn when n > x).
 * yn_rangedd: single forward-recurrence sweep filling n2 − n1 + 1 outputs
 * into `out`; cheaper than n2 − n1 + 1 independent yndd calls. */
MULTIFLOATS_API float64x2_t j0dd(float64x2_t a);
MULTIFLOATS_API float64x2_t j1dd(float64x2_t a);
MULTIFLOATS_API float64x2_t y0dd(float64x2_t a);
MULTIFLOATS_API float64x2_t y1dd(float64x2_t a);
MULTIFLOATS_API float64x2_t jndd(int n, float64x2_t a);
MULTIFLOATS_API float64x2_t yndd(int n, float64x2_t a);
MULTIFLOATS_API void yn_rangedd(int n1, int n2, float64x2_t a, float64x2_t *out);

/* Fused sincos / sinhcosh. One range-reduction + Taylor pair produces
 * both outputs, roughly halving the transcendental cost for call sites
 * that need both. Out-pointer style since C has no multi-value return. */
MULTIFLOATS_API void sincosdd(float64x2_t a, float64x2_t *s, float64x2_t *c);
MULTIFLOATS_API void sinhcoshdd(float64x2_t a, float64x2_t *s, float64x2_t *c);

/* Complex DD arithmetic. These are the canonical implementations; the
 * Fortran elemental `cdd + cdd`, `cdd * dp`, etc. routines wrap these via
 * bind(c). */
MULTIFLOATS_API complex64x2_t cadddd(complex64x2_t a, complex64x2_t b);
MULTIFLOATS_API complex64x2_t csubdd(complex64x2_t a, complex64x2_t b);
MULTIFLOATS_API complex64x2_t cmuldd(complex64x2_t a, complex64x2_t b);
MULTIFLOATS_API complex64x2_t cdivdd(complex64x2_t a, complex64x2_t b);

/* Complex DD transcendentals. Branch cuts match C99 Annex G (matching
 * libquadmath cexpq/clogq/csqrtq/...). Where the classic formula needs
 * both sin(y) and cos(y) or both sinh(y) and cosh(y), these use the fused
 * kernels internally so one range-reduction + Taylor pair covers both.
 * Precision and speed have been measured; the std::complex<MultiFloat<...>>
 * specializations in multifloats.hh delegate here where specialization
 * wins (exp, sin, cos, tan, sinh, cosh, tanh, atanh, acos). */
MULTIFLOATS_API complex64x2_t cexpdd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t cexpm1dd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t clogdd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t clog2dd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t clog10dd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t clog1pdd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t cpowdd(complex64x2_t z, complex64x2_t w);
MULTIFLOATS_API complex64x2_t csqrtdd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t csindd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t csinpidd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t ccosdd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t ccospidd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t ctandd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t casindd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t cacosdd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t catandd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t csinhdd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t ccoshdd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t ctanhdd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t casinhdd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t cacoshdd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t catanhdd(complex64x2_t z);

/* Complex magnitude / argument / projection / conjugate / accessors.
 * cabsdd = |z|  (overflow-safe hypot).
 * cargdd = arg(z) in (-pi, pi].
 * cprojdd = projection onto the Riemann sphere (C99 7.3.9.4):
 *   infinities collapse to (+inf, copysign(0, imag)), else identity. */
MULTIFLOATS_API float64x2_t   cabsdd(complex64x2_t z);
MULTIFLOATS_API float64x2_t   cargdd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t cprojdd(complex64x2_t z);
MULTIFLOATS_API complex64x2_t conjdd(complex64x2_t z);
MULTIFLOATS_API float64x2_t   crealdd(complex64x2_t z);
MULTIFLOATS_API float64x2_t   cimagdd(complex64x2_t z);

/* Matrix multiply (column-major, Fortran layout).
 *   matmuldd_mm: C(m,n) = A(m,k) * B(k,n)
 *   matmuldd_mv: y(m)   = A(m,k) * x(k)
 *   matmuldd_vm: y(n)   = x(k)   * B(k,n)
 * Leading dimensions equal the first extent (no strides).
 *
 * Scope vs BLAS GEMM. These are *not* a direct replacement for DGEMM:
 *
 *   - transa / transb: no. Always treats both operands in storage order;
 *     the caller must transpose the input arrays before the call if
 *     needed. (In Fortran, `matmul(transpose(a), b)` materialises the
 *     transpose, then the kernel walks the transposed buffer.)
 *   - alpha / beta:    no. The output is overwritten, not accumulated
 *     into. To compute `C := alpha*A*B + beta*C`, the caller must scale
 *     `A` or `B`, run matmul, then combine with C explicitly.
 *   - leading dimensions: always the first extent; no LDA/LDB/LDC.
 *
 * These constraints are deliberate: the compensated DD kernels use a
 * register-blocked panel design that assumes contiguous column-major
 * storage with the canonical shape. Adding transposed / strided / in-
 * place-accumulating variants is plausible but would require a new set
 * of panel dispatchers — tracked under "Deferred work" in
 * doc/developer/AUDIT_NOTES.md.
 *
 * renorm_interval: if > 0, renormalize accumulators every N reductions
 * (matches DD_FMA_RENORM_INTERVAL in the Fortran layer — keeps s_lo
 * bounded for large k). Pass 0 to renormalize only at the end. */
MULTIFLOATS_API void matmuldd_mm(const float64x2_t *a, const float64x2_t *b,
                         float64x2_t *c,
                         int64_t m, int64_t k, int64_t n,
                         int64_t renorm_interval);
MULTIFLOATS_API void matmuldd_mv(const float64x2_t *a, const float64x2_t *x,
                         float64x2_t *y,
                         int64_t m, int64_t k,
                         int64_t renorm_interval);
MULTIFLOATS_API void matmuldd_vm(const float64x2_t *x, const float64x2_t *b,
                         float64x2_t *y,
                         int64_t k, int64_t n,
                         int64_t renorm_interval);

/* Comparison (return int: 1 = true, 0 = false) */
MULTIFLOATS_API int eqdd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API int nedd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API int ltdd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API int ledd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API int gtdd(float64x2_t a, float64x2_t b);
MULTIFLOATS_API int gedd(float64x2_t a, float64x2_t b);

#ifdef __cplusplus
}
#endif

/* Scope the visibility macro to this header so callers don't get a leaked
 * `MULTIFLOATS_API` identifier. Consumers wanting their own attribute
 * macros should define them in their own namespace. */
#undef MULTIFLOATS_API
