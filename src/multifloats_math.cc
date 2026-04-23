// Umbrella translation unit for the DD math kernels and their C-ABI
// wrappers. The implementation is split across a set of `.inc` files
// that are pulled in here as a single compiled TU; that way the
// compiler keeps every cross-kernel call inlineable (e.g.
// `cexpdd → multifloats::exp + multifloats::sincos`) without depending
// on -flto.
//
//   float64x2_exp_log.inc     exp / exp2 / expm1 / log /
//                                    log2 / log10 / log1p / pow
//   float64x2_trig.inc        sin / cos / tan / sincos /
//                                    sinpi / cospi / tanpi +
//                                    inverse π-scaled helpers
//   float64x2_hyp.inc         sinh / cosh / tanh / sinhcosh
//   float64x2_inv_trig.inc    atan / asin / acos / atan2 +
//                                    asinh / acosh / atanh
//   float64x2_special.inc     erf / erfc / erfc_scaled /
//                                    tgamma / lgamma
//   float64x2_bessel.inc      J0 / J1 / Y0 / Y1 / Jn / Yn
//   float64x2_matmul.inc      compensated GEMM panels
//   float64x2_abi.inc  extern "C" scalar wrappers,
//                                    matmul entry points, comparisons
//   complex64x2_abi.inc extern "C" complex DD kernels
//
// Public transcendental kernels live in `namespace multifloats` with
// matching declarations in multifloats.h. The `extern "C" *dd` shims
// are thin marshaling wrappers around those same C++ functions — no
// parallel `_full`-named body. Truly internal helpers (range reducers,
// polynomial evaluators, triple-double-output variants) stay in anon
// namespaces inside each .inc file.

#include "multifloats.h"
#include "multifloats_td.hh"
#include "dd_constants.hh"
#include <cstdint>
#include <cstring>
#include <vector>

// Triple-double primitive bodies (declarations in multifloats_td.hh) and
// DD polynomial evaluators (horner/neval/deval — previously inline in
// multifloats.h). Defined here so the kernel .inc files see the same
// same-TU inline bodies they did before.
namespace multifloats {
namespace detail {
#include "float64x2_td.inc"
#include "float64x2_poly.inc"
} // namespace detail
} // namespace multifloats

// Public kernels — definitions. Matching declarations live in multifloats.h.
// Internal helpers (range reducers, polynomial Estrin kernels, TD variants)
// are wrapped in `namespace detail { }` inside each .inc so they don't
// pollute the public `multifloats::` surface. The `using namespace
// multifloats::detail` below lets same-TU callers reach them without the
// qualifier.
namespace multifloats {
using namespace multifloats::detail;  // neval, deval, horner, float64x3, kernels
#include "float64x2_exp_log.inc"
#include "float64x2_trig.inc"
#include "float64x2_hyp.inc"
#include "float64x2_inv_trig.inc"
#include "float64x2_special.inc"
#include "float64x2_bessel.inc"
} // namespace multifloats

// Anon-namespace helpers `dd_cross_diff` and `dd_x2y2m1` are defined below
// (around line 140 originally); the std::complex specializations and C-ABI
// wrappers consume them. The kernels.inc complex log1p also wants
// dd_x2y2m1, so we forward-declare them at file scope here. The same-anon-ns
// definition below provides the body; the redeclaration is permitted within
// a single TU as both refer to the same internal-linkage entity.
namespace {
multifloats::float64x2 dd_cross_diff(multifloats::float64x2 a,
                                     multifloats::float64x2 b,
                                     multifloats::float64x2 c,
                                     multifloats::float64x2 d);
multifloats::float64x2 dd_x2y2m1(multifloats::float64x2 x,
                                 multifloats::float64x2 y);
}

// Complex overloads for Kind-D parity (sinpi/cospi/expm1/log2/log10/log1p
// on std::complex<float64x2>). Pulled in after <complex> is visible via
// the inclusion of multifloats.h at top.
namespace multifloats {
using namespace multifloats::detail;  // sincos_td + TD primitives
#include "complex64x2_abi_kernels.inc"
} // namespace multifloats

// =============================================================================
// C-ABI entry points — extern "C" functions following math.h naming convention.
// =============================================================================

// File-scope C ABI helpers and `using` aliases so both the matmul anon-ns
// block and the extern "C" shim blocks can call them. Kept non-inline / at
// global namespace so the extern "C" block can reach them by unqualified
// name (the ABI .inc files live inside extern "C", which is not a namespace
// but inherits enclosing lookup).
using multifloats::float64x2;
using multifloats::detail::float64x3;
using multifloats::detail::td_add_td;
using multifloats::detail::td_mul_td;
using multifloats::detail::td_sub_double;
using multifloats::detail::td_to_dd;
using multifloats::detail::td_from_dd;

static inline float64x2 from(float64x2_t x) { float64x2 r; r._limbs[0] = x.hi; r._limbs[1] = x.lo; return r; }
static inline float64x2_t to(float64x2 const &x) { return {x._limbs[0], x._limbs[1]}; }
static inline std::complex<float64x2> cc_from(complex64x2_t z) { return {from(z.re), from(z.im)}; }
static inline complex64x2_t cc_to(std::complex<float64x2> const &z) { return {to(z.real()), to(z.imag())}; }

// Anon-namespace helpers shared by the std:: complex specializations and
// the C-ABI cmul/cdiv wrappers (which keep compensated bodies because the
// std::complex operator* template is not specializable).
//   • dd_cross_diff(a,b,c,d) = a·b − c·d with compensated 14-term expansion
//     through a triple-double accumulator, DD-accurate even when |a·b|≈|c·d|.
//   • dd_x2y2m1(x,y)         = x² + y² − 1 without the 1 ± 1 cancellation;
//     mirrors libquadmath __quadmath_x2y2m1q on DD inputs.
namespace {

__attribute__((noinline, cold))
float64x2 dd_cross_diff(float64x2 a, float64x2 b, float64x2 c, float64x2 d) {
  double a0 = a._limbs[0], a1 = a._limbs[1];
  double b0 = b._limbs[0], b1 = b._limbs[1];
  double c0 = c._limbs[0], c1 = c._limbs[1];
  double d0 = d._limbs[0], d1 = d._limbs[1];
  double ab00_h = a0 * b0, ab00_l = std::fma(a0, b0, -ab00_h);
  double cd00_h = c0 * d0, cd00_l = std::fma(c0, d0, -cd00_h);
  double ab01_h = a0 * b1, ab01_l = std::fma(a0, b1, -ab01_h);
  double ab10_h = a1 * b0, ab10_l = std::fma(a1, b0, -ab10_h);
  double cd01_h = c0 * d1, cd01_l = std::fma(c0, d1, -cd01_h);
  double cd10_h = c1 * d0, cd10_l = std::fma(c1, d0, -cd10_h);
  double ab11 = a1 * b1;
  double cd11 = c1 * d1;
  double T0 = 0, T1 = 0, T2 = 0;
  #define TSUM(x) do { \
    double _x = (x), _s, _bb, _e; \
    _s = T0 + _x; _bb = _s - T0; _e = (T0 - (_s - _bb)) + (_x - _bb); \
    T0 = _s; _x = _e; \
    _s = T1 + _x; _bb = _s - T1; _e = (T1 - (_s - _bb)) + (_x - _bb); \
    T1 = _s; T2 += _e; \
  } while (0)
  TSUM(ab00_h);  TSUM(-cd00_h);
  TSUM(ab00_l);  TSUM(-cd00_l);
  TSUM(ab01_h);  TSUM(ab10_h);  TSUM(-cd01_h); TSUM(-cd10_h);
  TSUM(ab01_l);  TSUM(ab10_l);  TSUM(-cd01_l); TSUM(-cd10_l);
  TSUM(ab11);    TSUM(-cd11);
  #undef TSUM
  { double s = T1 + T2, bb = s - T1; T2 = (T1 - (s - bb)) + (T2 - bb); T1 = s; }
  { double s = T0 + T1, bb = s - T0; T1 = (T0 - (s - bb)) + (T1 - bb); T0 = s; }
  T1 += T2;
  double hi = T0 + T1;
  double lo = T1 - (hi - T0);
  return float64x2(hi, lo);
}

inline float64x2 dd_x2y2m1(float64x2 x, float64x2 y) {
  double x0 = x._limbs[0], x1 = x._limbs[1];
  double y0 = y._limbs[0], y1 = y._limbs[1];
  double xx00_h = x0 * x0, xx00_l = std::fma(x0, x0, -xx00_h);
  double yy00_h = y0 * y0, yy00_l = std::fma(y0, y0, -yy00_h);
  double twox0 = x0 + x0, twoy0 = y0 + y0;          // exact (×2)
  double xx01_h = twox0 * x1, xx01_l = std::fma(twox0, x1, -xx01_h);
  double yy01_h = twoy0 * y1, yy01_l = std::fma(twoy0, y1, -yy01_h);
  double xx11 = x1 * x1;
  double yy11 = y1 * y1;
  double T0 = 0, T1 = 0, T2 = 0;
  #define TSUM(v) do { \
    double _x = (v), _s, _bb, _e; \
    _s = T0 + _x; _bb = _s - T0; _e = (T0 - (_s - _bb)) + (_x - _bb); \
    T0 = _s; _x = _e; \
    _s = T1 + _x; _bb = _s - T1; _e = (T1 - (_s - _bb)) + (_x - _bb); \
    T1 = _s; T2 += _e; \
  } while (0)
  TSUM(xx00_h);  TSUM(yy00_h);  TSUM(-1.0);
  TSUM(xx01_h);  TSUM(yy01_h);
  TSUM(xx00_l);  TSUM(yy00_l);
  TSUM(xx01_l);  TSUM(yy01_l);
  TSUM(xx11);    TSUM(yy11);
  #undef TSUM
  { double s = T1 + T2, bb = s - T1; T2 = (T1 - (s - bb)) + (T2 - bb); T1 = s; }
  { double s = T0 + T1, bb = s - T0; T1 = (T0 - (s - bb)) + (T1 - bb); T0 = s; }
  T1 += T2;
  double hi = T0 + T1;
  double lo = T1 - (hi - T0);
  return float64x2(hi, lo);
}

#include "float64x2_matmul.inc"
} // anonymous namespace

// Public C++ matmul entry points. The panel dispatchers above operate on
// `float64x2_t` (C-ABI struct) since they were originally the bodies of the
// extern "C" `matmuldd_*` shims. `multifloats::float64x2` is layout-compatible
// (asserted in multifloats.h next to the declarations), so we reinterpret
// pointers once at the boundary and hand off to the same dispatchers. The
// `matmuldd_*` shims in float64x2_abi.inc are now one-line
// wrappers calling these in the opposite direction.
namespace multifloats {

void matmul_mm(float64x2 const *a, float64x2 const *b, float64x2 *c,
               std::int64_t m, std::int64_t k, std::int64_t n,
               std::int64_t renorm_interval) {
  auto *ac = reinterpret_cast<float64x2_t const *>(a);
  auto *bc = reinterpret_cast<float64x2_t const *>(b);
  auto *cc = reinterpret_cast<float64x2_t *>(c);
  // NR-blocked mm: the MR-row × NR-col tile loads A[:,p] once per p and
  // reuses it across NR output columns, halving A-bandwidth vs a per-
  // column mv dispatch. Row-tail (1..MR-1 rows) and column-tail
  // (1..NR-1 cols) fall back to the single-column mv panel.
  constexpr int MR = 8;
  constexpr int NR = 2;
  std::int64_t j = 0;
  for (; j + NR <= n; j += NR) {
    std::int64_t i = 0;
    for (; i + MR <= m; i += MR) {
      gemm_panel<MR, NR>(ac + i, bc + j * k, cc + j * m + i,
                         m, k, m, k, renorm_interval);
    }
    int tail_m = static_cast<int>(m - i);
    if (tail_m > 0) {
      for (int jj = 0; jj < NR; ++jj) {
        gaxpy_mv_tail(ac + i, bc + (j + jj) * k, cc + (j + jj) * m + i,
                      tail_m, m, k, renorm_interval);
      }
    }
  }
  for (; j < n; ++j) {
    gaxpy_mv_dispatch(ac, bc + j * k, cc + j * m, m, k, m, renorm_interval);
  }
}

void matmul_mv(float64x2 const *a, float64x2 const *x, float64x2 *y,
               std::int64_t m, std::int64_t k,
               std::int64_t renorm_interval) {
  gaxpy_mv_dispatch(reinterpret_cast<float64x2_t const *>(a),
                    reinterpret_cast<float64x2_t const *>(x),
                    reinterpret_cast<float64x2_t *>(y),
                    m, k, m, renorm_interval);
}

// vm: y[j] = sum_p x[p] * B[p, j]. Column-major B makes B[:, j]
// contiguous at fixed j, so one scalar accumulator per output is optimal.
void matmul_vm(float64x2 const *x, float64x2 const *b, float64x2 *y,
               std::int64_t k, std::int64_t n,
               std::int64_t renorm_interval) {
  auto *xc = reinterpret_cast<float64x2_t const *>(x);
  auto *bc = reinterpret_cast<float64x2_t const *>(b);
  auto *yc = reinterpret_cast<float64x2_t *>(y);
  const bool simple = (renorm_interval <= 0) || (k <= renorm_interval);
  for (std::int64_t j = 0; j < n; ++j) {
    double s_hi = 0.0, s_lo = 0.0;
    const float64x2_t *__restrict__ bcol = bc + j * k;
    if (simple) {
      for (std::int64_t p = 0; p < k; ++p) {
        mac_inl(xc[p].hi, xc[p].lo, bcol[p].hi, bcol[p].lo, s_hi, s_lo);
      }
    } else {
      const std::int64_t chunk = renorm_interval;
      std::int64_t p0 = 0;
      while (p0 < k) {
        std::int64_t pend = p0 + chunk;
        if (pend > k) pend = k;
        for (std::int64_t p = p0; p < pend; ++p) {
          mac_inl(xc[p].hi, xc[p].lo, bcol[p].hi, bcol[p].lo, s_hi, s_lo);
        }
        p0 = pend;
        if (p0 < k) renorm_inl(s_hi, s_lo);
      }
    }
    yc[j] = finalize_inl(s_hi, s_lo);
  }
}

} // namespace multifloats

// std:: complex template specializations — bodies for every C++ <complex>
// free function we override (exp/log/sqrt/pow, the trig/hyp triples and
// their inverses, abs/arg/proj). Declarations live in multifloats.h.
// Pulled into namespace std so the explicit specialization syntax is in
// the right enclosing namespace.
namespace std {
#include "complex64x2_std.inc"
} // namespace std

extern "C" {
#include "float64x2_abi.inc"
#include "complex64x2_abi.inc"
} // extern "C"
