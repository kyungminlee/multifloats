// Umbrella translation unit for the DD math kernels and their C-ABI
// wrappers. The implementation is split across a set of `.inc` files
// that are pulled in here as a single compiled TU; that way the
// compiler keeps every cross-kernel call inlineable (e.g.
// `cexpdd → multifloats::exp + multifloats::sincos`) without depending
// on -flto.
//
//   float64x2/exp_log.inc     exp / exp2 / expm1 / log /
//                                    log2 / log10 / log1p / pow
//   float64x2/trig.inc        sin / cos / tan / sincos /
//                                    sinpi / cospi / tanpi +
//                                    inverse π-scaled helpers
//   float64x2/hyp.inc         sinh / cosh / tanh / sinhcosh
//   float64x2/inv_trig.inc    atan / asin / acos / atan2 +
//                                    asinh / acosh / atanh
//   float64x2/special.inc     erf / erfc / erfc_scaled /
//                                    tgamma / lgamma
//   float64x2/bessel.inc      J0 / J1 / Y0 / Y1 / Jn / Yn
//   float64x2/matmul.inc      compensated GEMM panels
//   float64x2/abi.inc         extern "C" scalar wrappers,
//                                    matmul entry points, comparisons
//   complex64x2/abi.inc       extern "C" complex DD kernels
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

// Triple-double primitive bodies (declarations in multifloats_td.hh) and
// DD polynomial evaluators (horner/neval/deval — previously inline in
// multifloats.h). Defined here so the kernel .inc files see the same
// same-TU inline bodies they did before.
namespace multifloats {
namespace detail {
#include "float64x2/td.inc"
#include "float64x2/poly.inc"
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
#include "float64x2/exp_log.inc"
#include "float64x2/trig.inc"
#include "float64x2/hyp.inc"
#include "float64x2/inv_trig.inc"
#include "float64x2/td_log_atan2.inc"
#include "float64x2/special.inc"
#include "float64x2/bessel.inc"
} // namespace multifloats

// Complex overloads for Kind-D parity (sinpi/cospi/expm1/log2/log10/log1p
// on std::complex<float64x2>). Pulled in after <complex> is visible via
// the inclusion of multifloats.h at top.
namespace multifloats {
using namespace multifloats::detail;  // sincos_td + TD primitives
#include "complex64x2/abi_kernels.inc"
} // namespace multifloats

// =============================================================================
// C-ABI entry points — extern "C" functions following math.h naming convention.
// =============================================================================

// Pull the two public types into TU-scope so the anon namespace helpers,
// the matmul wrappers, and the extern "C" shim blocks below can all name
// them unqualified. `complex64x2` is a distinct POD (two back-to-back
// float64x2s); the C-ABI shims in complex64x2/abi.inc marshal it to/from
// `std::complex<float64x2>` at the kernel boundary so C/C++/Fortran all
// see the same struct at the ABI layer (LTO type-identity match).
using multifloats::float64x2;
using multifloats::complex64x2;

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
  double a0 = a.limbs[0], a1 = a.limbs[1];
  double b0 = b.limbs[0], b1 = b.limbs[1];
  double c0 = c.limbs[0], c1 = c.limbs[1];
  double d0 = d.limbs[0], d1 = d.limbs[1];
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
  double x0 = x.limbs[0], x1 = x.limbs[1];
  double y0 = y.limbs[0], y1 = y.limbs[1];
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

// Matmul is the only place DD arithmetic vectorizes to PACKED FMA (vfmadd…pd);
// the scalar transcendental EFTs are sequential chains that don't. So matmul is
// compiled twice: a baseline copy (mm_scalar) that runs on any supported CPU,
// and — on x86-64 GCC — an FMA copy (mm_vec) whose panels vectorize. The WHOLE
// implementation (panels + drivers) is recompiled for the FMA target via the
// block-level `#pragma GCC target`, so gemm_panel itself vectorizes; the public
// matmul_mm/mv/vm pick a copy by CPUID at runtime. Measured ~5x faster on an
// AVX2+FMA CPU than the baseline (and than a target_clones+flatten variant,
// whose clone did not actually vectorize the hot loop).
//
// GCC only: only GCC's vectorizer packs this compensated DD MAC loop — Clang/
// icx emit scalar FMA even at -march=haswell, so a vectorized copy buys them
// nothing and they use the baseline (still hardware *scalar* FMA at runtime via
// glibc's `fma` ifunc). Restricting to GCC also avoids the `#pragma clang
// attribute` path on Clang/icx. AArch64 has FMA in its base ISA (no dispatch);
// Apple/Mach-O excluded. Opt out with -DMULTIFLOATS_NO_MM_DISPATCH.
namespace mm_scalar {
#include "float64x2/matmul.inc"
}
#if defined(__x86_64__) && !defined(__APPLE__) && defined(__GNUC__) &&         \
    !defined(__clang__) && !defined(MULTIFLOATS_NO_MM_DISPATCH)
#define MULTIFLOATS_MM_DISPATCH_ACTIVE 1
#pragma GCC push_options
#pragma GCC target("arch=haswell")
namespace mm_vec {
#include "float64x2/matmul.inc"
}
#pragma GCC pop_options
#endif  // dispatch active
} // anonymous namespace

// Public C++ matmul entry points. Each dispatches to the baseline (mm_scalar)
// or the FMA-vectorized (mm_vec) copy of the driver by CPUID at runtime. The
// `matmuldd_*` shims in float64x2/abi.inc forward to these one-for-one.
namespace multifloats {

#ifdef MULTIFLOATS_MM_DISPATCH_ACTIVE
// CPUID is read once into a function-local static — evaluated on first call,
// well after libgcc's __cpu_model constructor runs — then cached.
static bool mm_use_vec() {
  static const bool v =
      __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
  return v;
}
#define MULTIFLOATS_MM_RUN(impl, ...)                                          \
  do {                                                                        \
    if (mm_use_vec()) ::mm_vec::impl(__VA_ARGS__);                            \
    else ::mm_scalar::impl(__VA_ARGS__);                                      \
  } while (0)
#else
#define MULTIFLOATS_MM_RUN(impl, ...) ::mm_scalar::impl(__VA_ARGS__)
#endif

void matmul_mm(float64x2 const *a, float64x2 const *b, float64x2 *c,
               std::int64_t m, std::int64_t k, std::int64_t n,
               std::int64_t renorm_interval) {
  MULTIFLOATS_MM_RUN(mm_impl, a, b, c, m, k, n, renorm_interval);
}

void matmul_mv(float64x2 const *a, float64x2 const *x, float64x2 *y,
               std::int64_t m, std::int64_t k, std::int64_t renorm_interval) {
  MULTIFLOATS_MM_RUN(mv_impl, a, x, y, m, k, renorm_interval);
}

void matmul_vm(float64x2 const *x, float64x2 const *b, float64x2 *y,
               std::int64_t k, std::int64_t n, std::int64_t renorm_interval) {
  MULTIFLOATS_MM_RUN(vm_impl, x, b, y, k, n, renorm_interval);
}

} // namespace multifloats

// std:: complex template specializations — bodies for every C++ <complex>
// free function we override (exp/log/sqrt/pow, the trig/hyp triples and
// their inverses, abs/arg/proj). Declarations live in multifloats.h.
// Pulled into namespace std so the explicit specialization syntax is in
// the right enclosing namespace.
namespace std {
#include "complex64x2/std.inc"
} // namespace std

// The C-ABI *dd shims live in `namespace multifloats` so their declarations
// (in multifloats.h) and definitions share the same C++ qualified name.
// `extern "C"` makes the linker symbol unmangled (`adddd`, not
// `_ZN11multifloats5adddd...`), matching the Fortran bind(c) and C clients;
// the namespace scope only affects C++ lookup/ADL.
namespace multifloats {
extern "C" {
#include "float64x2/abi.inc"
#include "complex64x2/abi.inc"
#include "float64x2/f64x2_aliases.inc"
} // extern "C"
} // namespace multifloats
