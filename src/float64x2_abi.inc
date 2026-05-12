
// Every `NAMEdd` here is a thin marshaling shim around the matching
// `multifloats::NAME` free function (defined in multifloats_math.cc via
// the per-topic .inc files). `float64x2` / `complex64x2` are the unified
// C/C++ types, so the boundary is a plain call — no conversion needed.
// The MF_ABI_* macros collapse the boilerplate; shims that don't fit the
// pattern (multi-output fused kernels, range sweeps, matmul, integer-first
// Bessel) stay explicit below.

#define MF_ABI_UNARY(NAME) \
  float64x2 NAME##dd(float64x2 a) { return multifloats::NAME(a); }
#define MF_ABI_BINARY(NAME) \
  float64x2 NAME##dd(float64x2 a, float64x2 b) { \
    return multifloats::NAME(a, b); \
  }
#define MF_ABI_BINARY_OP(NAME, OP) \
  float64x2 NAME##dd(float64x2 a, float64x2 b) { return a OP b; }
#define MF_ABI_CMP(NAME, OP) \
  int NAME##dd(float64x2 a, float64x2 b) { return a OP b ? 1 : 0; }

// Arithmetic — C++ operators, not free functions.
MF_ABI_BINARY_OP(add, +)
MF_ABI_BINARY_OP(sub, -)
MF_ABI_BINARY_OP(mul, *)
MF_ABI_BINARY_OP(div, /)

// Unary
float64x2 negdd(float64x2 a) { return -a; }
MF_ABI_UNARY(fabs)
MF_ABI_UNARY(sqrt)

// Rounding (Fortran AINT / ANINT delegate here)
MF_ABI_UNARY(trunc)
MF_ABI_UNARY(round)

// Binary
MF_ABI_BINARY(fmin)
MF_ABI_BINARY(fmax)
MF_ABI_BINARY(hypot)
MF_ABI_BINARY(pow)
// Integer-exponent power: exponent is a plain int, not a DD, so it
// does not fit the MF_ABI_BINARY macro shape.
float64x2 powidd(float64x2 a, int n) { return multifloats::powi(a, n); }
MF_ABI_BINARY(fmod)
MF_ABI_BINARY(modulo)
MF_ABI_BINARY(fdim)
MF_ABI_BINARY(copysign)
float64x2 fmadd(float64x2 a, float64x2 b, float64x2 c) {
  return multifloats::fma(a, b, c);
}

// Exponential / logarithmic
MF_ABI_UNARY(exp)
MF_ABI_UNARY(exp2)
MF_ABI_UNARY(exp10)
MF_ABI_UNARY(expm1)
MF_ABI_UNARY(exp2m1)
MF_ABI_UNARY(exp10m1)
MF_ABI_UNARY(log)
MF_ABI_UNARY(log2)
MF_ABI_UNARY(log10)
MF_ABI_UNARY(log1p)
MF_ABI_UNARY(log2p1)
MF_ABI_UNARY(log10p1)

// Trigonometric
MF_ABI_UNARY(sin)
MF_ABI_UNARY(cos)
MF_ABI_UNARY(tan)
MF_ABI_UNARY(asin)
MF_ABI_UNARY(acos)
MF_ABI_UNARY(atan)
MF_ABI_BINARY(atan2)

// π-scaled trig: fn_pi(x) = fn(π·x) or fn(x)/π.
MF_ABI_UNARY(sinpi)
MF_ABI_UNARY(cospi)
MF_ABI_UNARY(tanpi)
MF_ABI_UNARY(asinpi)
MF_ABI_UNARY(acospi)
MF_ABI_UNARY(atanpi)
MF_ABI_BINARY(atan2pi)

// Hyperbolic
MF_ABI_UNARY(sinh)
MF_ABI_UNARY(cosh)
MF_ABI_UNARY(tanh)
MF_ABI_UNARY(asinh)
MF_ABI_UNARY(acosh)
MF_ABI_UNARY(atanh)

// Error functions
MF_ABI_UNARY(erf)
MF_ABI_UNARY(erfc)
MF_ABI_UNARY(erfcx)

// Gamma functions
MF_ABI_UNARY(tgamma)
MF_ABI_UNARY(lgamma)

// Bessel functions (math.h naming: j0, j1, y0, y1). Route to the
// order-specific C++ fast-path kernels in multifloats.h.
MF_ABI_UNARY(j0)
MF_ABI_UNARY(j1)
MF_ABI_UNARY(y0)
MF_ABI_UNARY(y1)

// Integer-order Bessel functions (recurrence from j0/j1, y0/y1).
float64x2 jndd(int n, float64x2 a) { return multifloats::jn(n, a); }
float64x2 yndd(int n, float64x2 a) { return multifloats::yn(n, a); }
// Range sweeps: Yn uses the stable forward-recurrence single-sweep kernel;
// Jn loops over Miller-backward per output (forward is unstable for n > x).
void jndd_range(int n1, int n2, float64x2 a, float64x2 *out) {
  if (n2 < n1) return;
  multifloats::jn(n1, n2, a, out);
}
void yndd_range(int n1, int n2, float64x2 a, float64x2 *out) {
  if (n2 < n1) return;
  multifloats::yn(n1, n2, a, out);
}

// Fused sincos / sinhcosh. Out-pointer style (C has no multi-value return).
void sincosdd(float64x2 a, float64x2 *s, float64x2 *c) {
  multifloats::sincos(a, *s, *c);
}

void sinhcoshdd(float64x2 a, float64x2 *s, float64x2 *c) {
  multifloats::sinhcosh(a, *s, *c);
}

// Matmul: thin wrappers around multifloats::matmul_* (implementations in
// multifloats_math.cc).
void matmuldd_mm(const float64x2 *a, const float64x2 *b,
                 float64x2 *c, int64_t m, int64_t k, int64_t n,
                 int64_t renorm_interval) {
  multifloats::matmul_mm(a, b, c, m, k, n, renorm_interval);
}

void matmuldd_mv(const float64x2 *a, const float64x2 *x,
                 float64x2 *y, int64_t m, int64_t k,
                 int64_t renorm_interval) {
  multifloats::matmul_mv(a, x, y, m, k, renorm_interval);
}

void matmuldd_vm(const float64x2 *x, const float64x2 *b,
                 float64x2 *y, int64_t k, int64_t n,
                 int64_t renorm_interval) {
  multifloats::matmul_vm(x, b, y, k, n, renorm_interval);
}

// Comparison (return int: 1 = true, 0 = false)
MF_ABI_CMP(eq, ==)
MF_ABI_CMP(ne, !=)
MF_ABI_CMP(lt, <)
MF_ABI_CMP(le, <=)
MF_ABI_CMP(gt, >)
MF_ABI_CMP(ge, >=)

// libquadmath parity — every NAMEq in <quadmath.h> has a matching NAMEdd
// shim below. Each delegates to the canonical multifloats:: helper.

// Floating-point manipulation
MF_ABI_UNARY(cbrt)
MF_ABI_UNARY(ceil)
MF_ABI_UNARY(floor)
MF_ABI_UNARY(nearbyint)
MF_ABI_UNARY(rint)
MF_ABI_UNARY(logb)
MF_ABI_BINARY(nextafter)
MF_ABI_BINARY(remainder)

// Out-pointer / int-arg variants — don't fit MF_ABI_UNARY/BINARY.
float64x2 frexpdd(float64x2 a, int *exp) { return multifloats::frexp(a, exp); }
float64x2 modfdd(float64x2 a, float64x2 *iptr) {
  return multifloats::modf(a, iptr);
}
float64x2 ldexpdd(float64x2 a, int n) { return multifloats::ldexp(a, n); }
float64x2 scalbndd(float64x2 a, int n) { return multifloats::scalbn(a, n); }
float64x2 scalblndd(float64x2 a, long n) { return multifloats::scalbln(a, n); }
float64x2 remquodd(float64x2 a, float64x2 b, int *quo) {
  return multifloats::remquo(a, b, quo);
}
int ilogbdd(float64x2 a) { return multifloats::ilogb(a); }

// Integer rounding
long      lrounddd (float64x2 a) { return multifloats::lround (a); }
long long llrounddd(float64x2 a) { return multifloats::llround(a); }
long      lrintdd  (float64x2 a) { return multifloats::lrint  (a); }
long long llrintdd (float64x2 a) { return multifloats::llrint (a); }

// Classification
int isnandd   (float64x2 a) { return multifloats::isnan   (a) ? 1 : 0; }
int isinfdd   (float64x2 a) { return multifloats::isinf   (a) ? 1 : 0; }
int isfinitedd(float64x2 a) { return multifloats::isfinite(a) ? 1 : 0; }
int finitedd  (float64x2 a) { return multifloats::isfinite(a) ? 1 : 0; }
int isnormaldd(float64x2 a) { return multifloats::isnormal(a) ? 1 : 0; }
int signbitdd (float64x2 a) { return multifloats::signbit (a) ? 1 : 0; }
int fpclassifydd(float64x2 a) { return multifloats::fpclassify(a); }

// C99 ordered/unordered comparison predicates
int isunordereddd   (float64x2 a, float64x2 b) { return multifloats::isunordered   (a, b) ? 1 : 0; }
int isgreaterdd     (float64x2 a, float64x2 b) { return multifloats::isgreater     (a, b) ? 1 : 0; }
int isgreaterequaldd(float64x2 a, float64x2 b) { return multifloats::isgreaterequal(a, b) ? 1 : 0; }
int islessdd        (float64x2 a, float64x2 b) { return multifloats::isless        (a, b) ? 1 : 0; }
int islessequaldd   (float64x2 a, float64x2 b) { return multifloats::islessequal   (a, b) ? 1 : 0; }
int islessgreaterdd (float64x2 a, float64x2 b) { return multifloats::islessgreater (a, b) ? 1 : 0; }

int issignalingdd(float64x2 a) {
#if defined(__has_builtin)
#  if __has_builtin(__builtin_issignaling)
  return __builtin_issignaling(a.limbs[0]) ? 1 : 0;
#  else
  (void)a; return 0;
#  endif
#elif defined(__GNUC__) && __GNUC__ >= 13 && !defined(__clang__)
  return __builtin_issignaling(a.limbs[0]) ? 1 : 0;
#else
  (void)a; return 0;
#endif
}

// NaN with payload tag
float64x2 nandd(const char *tagp) { return multifloats::nan(tagp); }

// cis(x) = cos(x) + i sin(x)
complex64x2 cexpidd(float64x2 a) {
  auto r = multifloats::cexpi(a);
  return { r.real(), r.imag() };
}

#undef MF_ABI_UNARY
#undef MF_ABI_BINARY
#undef MF_ABI_BINARY_OP
#undef MF_ABI_CMP
