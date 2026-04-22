
// Every `NAMEdd` here is a thin marshaling shim around the matching
// `multifloats::NAME` free function (defined in multifloats_math.cc via
// the per-topic .inc files). They cross `float64x2_t` ↔ `float64x2` at
// the boundary using the file-scope `from` / `to` helpers and otherwise
// don't add behavior. The MF_ABI_* macros collapse the boilerplate;
// shims that don't fit the pattern (multi-output fused kernels, range
// sweeps, matmul, integer-first Bessel) stay explicit below.

#define MF_ABI_UNARY(NAME) \
  float64x2_t NAME##dd(float64x2_t a) { return to(multifloats::NAME(from(a))); }
#define MF_ABI_UNARY_AS(CNAME, CXXNAME) \
  float64x2_t CNAME(float64x2_t a) { return to(multifloats::CXXNAME(from(a))); }
#define MF_ABI_BINARY(NAME) \
  float64x2_t NAME##dd(float64x2_t a, float64x2_t b) { \
    return to(multifloats::NAME(from(a), from(b))); \
  }
#define MF_ABI_BINARY_OP(NAME, OP) \
  float64x2_t NAME##dd(float64x2_t a, float64x2_t b) { \
    return to(from(a) OP from(b)); \
  }
#define MF_ABI_CMP(NAME, OP) \
  int NAME##dd(float64x2_t a, float64x2_t b) { \
    return from(a) OP from(b) ? 1 : 0; \
  }

// Arithmetic — C++ operators, not free functions.
MF_ABI_BINARY_OP(add, +)
MF_ABI_BINARY_OP(sub, -)
MF_ABI_BINARY_OP(mul, *)
MF_ABI_BINARY_OP(div, /)

// Unary
float64x2_t negdd(float64x2_t a) { return to(-from(a)); }
MF_ABI_UNARY_AS(fabsdd, abs)  // C's fabs is C++'s abs
MF_ABI_UNARY(sqrt)

// Rounding (Fortran AINT / ANINT delegate here)
MF_ABI_UNARY(trunc)
MF_ABI_UNARY(round)

// Binary
MF_ABI_BINARY(fmin)
MF_ABI_BINARY(fmax)
MF_ABI_BINARY(hypot)
MF_ABI_BINARY(pow)
MF_ABI_BINARY(fmod)
MF_ABI_BINARY(fdim)
MF_ABI_BINARY(copysign)
float64x2_t fmadd(float64x2_t a, float64x2_t b, float64x2_t c) {
  return to(multifloats::fma(from(a), from(b), from(c)));
}

// Exponential / logarithmic
MF_ABI_UNARY(exp)
MF_ABI_UNARY(exp2)
MF_ABI_UNARY(expm1)
MF_ABI_UNARY(log)
MF_ABI_UNARY(log2)
MF_ABI_UNARY(log10)
MF_ABI_UNARY(log1p)

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
// order-specific C++ fast-path kernels in multifloats.hh.
MF_ABI_UNARY(j0)
MF_ABI_UNARY(j1)
MF_ABI_UNARY(y0)
MF_ABI_UNARY(y1)

// Integer-order Bessel functions (recurrence from j0/j1, y0/y1).
float64x2_t jndd(int n, float64x2_t a) { return to(multifloats::jn(n, from(a))); }
float64x2_t yndd(int n, float64x2_t a) { return to(multifloats::yn(n, from(a))); }
// Range sweeps: Yn uses the stable forward-recurrence single-sweep kernel;
// Jn loops over Miller-backward per output (forward is unstable for n > x).
void jndd_range(int n1, int n2, float64x2_t a, float64x2_t *out) {
  int size = n2 - n1 + 1;
  if (size <= 0) return;
  std::vector<float64x2> buf(size);
  multifloats::jn(n1, n2, from(a), buf.data());
  for (int i = 0; i < size; ++i) out[i] = to(buf[i]);
}
void yndd_range(int n1, int n2, float64x2_t a, float64x2_t *out) {
  int size = n2 - n1 + 1;
  if (size <= 0) return;
  std::vector<float64x2> buf(size);
  multifloats::yn(n1, n2, from(a), buf.data());
  for (int i = 0; i < size; ++i) out[i] = to(buf[i]);
}

// Fused sincos / sinhcosh. Out-pointer style (C has no multi-value return).
void sincosdd(float64x2_t a, float64x2_t *s, float64x2_t *c) {
  float64x2 ss, cc;
  multifloats::sincos(from(a), ss, cc);
  *s = to(ss);
  *c = to(cc);
}

void sinhcoshdd(float64x2_t a, float64x2_t *s, float64x2_t *c) {
  float64x2 ss, cc;
  multifloats::sinhcosh(from(a), ss, cc);
  *s = to(ss);
  *c = to(cc);
}

// Matmul: thin wrappers around multifloats::matmul_* (implementations in
// multifloats_math.cc). The reinterpret_cast at the boundary is backed by a
// static_assert on layout compatibility in multifloats.hh.
void matmuldd_mm(const float64x2_t *a, const float64x2_t *b,
                 float64x2_t *c, int64_t m, int64_t k, int64_t n,
                 int64_t renorm_interval) {
  multifloats::matmul_mm(
      reinterpret_cast<const multifloats::float64x2 *>(a),
      reinterpret_cast<const multifloats::float64x2 *>(b),
      reinterpret_cast<multifloats::float64x2 *>(c),
      m, k, n, renorm_interval);
}

void matmuldd_mv(const float64x2_t *a, const float64x2_t *x,
                 float64x2_t *y, int64_t m, int64_t k,
                 int64_t renorm_interval) {
  multifloats::matmul_mv(
      reinterpret_cast<const multifloats::float64x2 *>(a),
      reinterpret_cast<const multifloats::float64x2 *>(x),
      reinterpret_cast<multifloats::float64x2 *>(y),
      m, k, renorm_interval);
}

void matmuldd_vm(const float64x2_t *x, const float64x2_t *b,
                 float64x2_t *y, int64_t k, int64_t n,
                 int64_t renorm_interval) {
  multifloats::matmul_vm(
      reinterpret_cast<const multifloats::float64x2 *>(x),
      reinterpret_cast<const multifloats::float64x2 *>(b),
      reinterpret_cast<multifloats::float64x2 *>(y),
      k, n, renorm_interval);
}

// Comparison (return int: 1 = true, 0 = false)
MF_ABI_CMP(eq, ==)
MF_ABI_CMP(ne, !=)
MF_ABI_CMP(lt, <)
MF_ABI_CMP(le, <=)
MF_ABI_CMP(gt, >)
MF_ABI_CMP(ge, >=)

#undef MF_ABI_UNARY
#undef MF_ABI_UNARY_AS
#undef MF_ABI_BINARY
#undef MF_ABI_BINARY_OP
#undef MF_ABI_CMP
