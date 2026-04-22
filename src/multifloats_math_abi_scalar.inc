
// Arithmetic
float64x2_t adddd(float64x2_t a, float64x2_t b) { return to(from(a) + from(b)); }
float64x2_t subdd(float64x2_t a, float64x2_t b) { return to(from(a) - from(b)); }
float64x2_t muldd(float64x2_t a, float64x2_t b) { return to(from(a) * from(b)); }
float64x2_t divdd(float64x2_t a, float64x2_t b) { return to(from(a) / from(b)); }
float64x2_t negdd(float64x2_t a)  { return to(-from(a)); }
float64x2_t fabsdd(float64x2_t a)  { return to(multifloats::abs(from(a))); }
float64x2_t sqrtdd(float64x2_t a) { return to(multifloats::sqrt(from(a))); }

// Rounding (Fortran AINT / ANINT delegate here)
float64x2_t truncdd(float64x2_t a) { return to(multifloats::trunc(from(a))); }
float64x2_t rounddd(float64x2_t a) { return to(multifloats::round(from(a))); }

// Binary
float64x2_t fmindd(float64x2_t a, float64x2_t b)     { return to(multifloats::fmin(from(a), from(b))); }
float64x2_t fmaxdd(float64x2_t a, float64x2_t b)     { return to(multifloats::fmax(from(a), from(b))); }
float64x2_t hypotdd(float64x2_t a, float64x2_t b)    { return to(multifloats::hypot(from(a), from(b))); }
float64x2_t powdd(float64x2_t a, float64x2_t b)      { return to(multifloats::pow(from(a), from(b))); }
float64x2_t fmoddd(float64x2_t a, float64x2_t b)     { return to(multifloats::fmod(from(a), from(b))); }
float64x2_t fdimdd(float64x2_t a, float64x2_t b)     { return to(multifloats::fdim(from(a), from(b))); }
float64x2_t copysigndd(float64x2_t a, float64x2_t b) { return to(multifloats::copysign(from(a), from(b))); }
float64x2_t fmadd(float64x2_t a, float64x2_t b, float64x2_t c) { return to(multifloats::fma(from(a), from(b), from(c))); }

// Exponential / logarithmic
float64x2_t expdd(float64x2_t a)   { return to(multifloats::exp(from(a))); }
float64x2_t exp2dd(float64x2_t a)  { return to(multifloats::exp2(from(a))); }
float64x2_t expm1dd(float64x2_t a) { return to(multifloats::expm1(from(a))); }
float64x2_t logdd(float64x2_t a)   { return to(multifloats::log(from(a))); }
float64x2_t log2dd(float64x2_t a)  { return to(multifloats::log2(from(a))); }
float64x2_t log10dd(float64x2_t a) { return to(multifloats::log10(from(a))); }
float64x2_t log1pdd(float64x2_t a) { return to(multifloats::log1p(from(a))); }

// Trigonometric
float64x2_t sindd(float64x2_t a)   { return to(multifloats::sin(from(a))); }
float64x2_t cosdd(float64x2_t a)   { return to(multifloats::cos(from(a))); }
float64x2_t tandd(float64x2_t a)   { return to(multifloats::tan(from(a))); }
float64x2_t asindd(float64x2_t a)  { return to(multifloats::asin(from(a))); }
float64x2_t acosdd(float64x2_t a)  { return to(multifloats::acos(from(a))); }
float64x2_t atandd(float64x2_t a)  { return to(multifloats::atan(from(a))); }
float64x2_t atan2dd(float64x2_t a, float64x2_t b) { return to(multifloats::atan2(from(a), from(b))); }

// π-scaled trig: fn_pi(x) = fn(π·x) or fn(x)/π.
float64x2_t sinpidd(float64x2_t a)  { return to(multifloats::sinpi(from(a))); }
float64x2_t cospidd(float64x2_t a)  { return to(multifloats::cospi(from(a))); }
float64x2_t tanpidd(float64x2_t a)  { return to(multifloats::tanpi(from(a))); }
float64x2_t asinpidd(float64x2_t a) { return to(multifloats::asinpi(from(a))); }
float64x2_t acospidd(float64x2_t a) { return to(multifloats::acospi(from(a))); }
float64x2_t atanpidd(float64x2_t a) { return to(multifloats::atanpi(from(a))); }
float64x2_t atan2pidd(float64x2_t a, float64x2_t b) { return to(multifloats::atan2pi(from(a), from(b))); }

// Hyperbolic
float64x2_t sinhdd(float64x2_t a)  { return to(multifloats::sinh(from(a))); }
float64x2_t coshdd(float64x2_t a)  { return to(multifloats::cosh(from(a))); }
float64x2_t tanhdd(float64x2_t a)  { return to(multifloats::tanh(from(a))); }
float64x2_t asinhdd(float64x2_t a) { return to(multifloats::asinh(from(a))); }
float64x2_t acoshdd(float64x2_t a) { return to(multifloats::acosh(from(a))); }
float64x2_t atanhdd(float64x2_t a) { return to(multifloats::atanh(from(a))); }

// Error functions
float64x2_t erfdd(float64x2_t a)   { return to(multifloats::erf(from(a))); }
float64x2_t erfcdd(float64x2_t a)  { return to(multifloats::erfc(from(a))); }
float64x2_t erfcxdd(float64x2_t a) { return to(multifloats::erfcx(from(a))); }

// Gamma functions
float64x2_t tgammadd(float64x2_t a) { return to(multifloats::tgamma(from(a))); }
float64x2_t lgammadd(float64x2_t a) { return to(multifloats::lgamma(from(a))); }

// Bessel functions (math.h naming: j0, j1, y0, y1). The fast-path kernels
// multifloats::bessel_j0/j1/y0/y1 are defined in multifloats_math_bessel.inc
// and live in namespace multifloats as internal helpers — the C ABI is the
// only public surface that reaches them directly.
float64x2_t j0dd(float64x2_t a) { return to(multifloats::bessel_j0(from(a))); }
float64x2_t j1dd(float64x2_t a) { return to(multifloats::bessel_j1(from(a))); }
float64x2_t y0dd(float64x2_t a) { return to(multifloats::bessel_y0(from(a))); }
float64x2_t y1dd(float64x2_t a) { return to(multifloats::bessel_y1(from(a))); }

// Integer-order Bessel functions (recurrence from j0/j1, y0/y1).
float64x2_t jndd(int n, float64x2_t a) { return to(multifloats::cyl_bessel_j(n, from(a))); }
float64x2_t yndd(int n, float64x2_t a) { return to(multifloats::cyl_neumann(n, from(a))); }
void yn_rangedd(int n1, int n2, float64x2_t a, float64x2_t *out) {
  int size = n2 - n1 + 1;
  if (size <= 0) return;
  std::vector<float64x2> buf(size);
  multifloats::cyl_neumann(n1, n2, from(a), buf.data());
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

void matmuldd_mm(const float64x2_t *__restrict__ a, const float64x2_t *__restrict__ b,
                  float64x2_t *__restrict__ c, int64_t m, int64_t k, int64_t n,
                  int64_t renorm_interval) {
  // NR-blocked mm: the MR-row × NR-col tile loads A[:,p] once per p and
  // reuses it across NR output columns, halving A-bandwidth vs a per-
  // column mv dispatch. Row-tail (1..MR-1 rows) and column-tail
  // (1..NR-1 cols) fall back to the single-column mv panel.
  constexpr int MR = 8;
  constexpr int NR = 2;
  int64_t j = 0;
  for (; j + NR <= n; j += NR) {
    int64_t i = 0;
    for (; i + MR <= m; i += MR) {
      gemm_panel<MR, NR>(a + i, b + j * k, c + j * m + i,
                            m, k, m, k, renorm_interval);
    }
    int tail_m = static_cast<int>(m - i);
    if (tail_m > 0) {
      for (int jj = 0; jj < NR; ++jj) {
        gaxpy_mv_tail(a + i, b + (j + jj) * k, c + (j + jj) * m + i,
                         tail_m, m, k, renorm_interval);
      }
    }
  }
  for (; j < n; ++j) {
    gaxpy_mv_dispatch(a, b + j * k, c + j * m, m, k, m, renorm_interval);
  }
}

void matmuldd_mv(const float64x2_t *__restrict__ a, const float64x2_t *__restrict__ x,
                  float64x2_t *__restrict__ y, int64_t m, int64_t k,
                  int64_t renorm_interval) {
  gaxpy_mv_dispatch(a, x, y, m, k, m, renorm_interval);
}

// vm: y[j] = sum_p x[p] * B[p, j]. Column-major B makes B[:, j]
// contiguous at fixed j, so one scalar accumulator per output is optimal.
void matmuldd_vm(const float64x2_t *__restrict__ x, const float64x2_t *__restrict__ b,
                  float64x2_t *__restrict__ y, int64_t k, int64_t n,
                  int64_t renorm_interval) {
  const bool simple = (renorm_interval <= 0) || (k <= renorm_interval);
  for (int64_t j = 0; j < n; ++j) {
    double s_hi = 0.0, s_lo = 0.0;
    const float64x2_t *__restrict__ bcol = b + j * k;
    if (simple) {
      for (int64_t p = 0; p < k; ++p) {
        mac_inl(x[p].hi, x[p].lo, bcol[p].hi, bcol[p].lo, s_hi, s_lo);
      }
    } else {
      const int64_t chunk = renorm_interval;
      int64_t p0 = 0;
      while (p0 < k) {
        int64_t pend = p0 + chunk;
        if (pend > k) pend = k;
        for (int64_t p = p0; p < pend; ++p) {
          mac_inl(x[p].hi, x[p].lo, bcol[p].hi, bcol[p].lo, s_hi, s_lo);
        }
        p0 = pend;
        if (p0 < k) renorm_inl(s_hi, s_lo);
      }
    }
    y[j] = finalize_inl(s_hi, s_lo);
  }
}

// Comparison
int eqdd(float64x2_t a, float64x2_t b) { return from(a) == from(b) ? 1 : 0; }
int nedd(float64x2_t a, float64x2_t b) { return from(a) != from(b) ? 1 : 0; }
int ltdd(float64x2_t a, float64x2_t b) { return from(a) <  from(b) ? 1 : 0; }
int ledd(float64x2_t a, float64x2_t b) { return from(a) <= from(b) ? 1 : 0; }
int gtdd(float64x2_t a, float64x2_t b) { return from(a) >  from(b) ? 1 : 0; }
int gedd(float64x2_t a, float64x2_t b) { return from(a) >= from(b) ? 1 : 0; }
