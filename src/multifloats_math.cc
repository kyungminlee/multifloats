// Umbrella translation unit for the DD math kernels and their C-ABI
// wrappers. The implementation is split across a set of `.inc` files
// that are pulled in here as a single compiled TU; that way the
// compiler keeps every cross-kernel call inlineable (e.g.
// `cexpdd → multifloats::exp + multifloats::sincos`) without depending
// on -flto.
//
//   multifloats_math_exp_log.inc     exp / exp2 / expm1 / log /
//                                    log2 / log10 / log1p / pow
//   multifloats_math_trig.inc        sin / cos / tan / sincos /
//                                    sinpi / cospi / tanpi +
//                                    inverse π-scaled helpers
//   multifloats_math_hyp.inc         sinh / cosh / tanh / sinhcosh
//   multifloats_math_inv_trig.inc    atan / asin / acos / atan2 +
//                                    asinh / acosh / atanh
//   multifloats_math_special.inc     erf / erfc / erfc_scaled /
//                                    tgamma / lgamma
//   multifloats_math_bessel.inc      J0 / J1 / Y0 / Y1 / Jn / Yn
//   multifloats_math_matmul.inc      compensated GEMM panels
//   multifloats_math_abi_scalar.inc  extern "C" scalar wrappers,
//                                    matmul entry points, comparisons
//   multifloats_math_abi_complex.inc extern "C" complex DD kernels
//
// Public transcendental kernels live in `namespace multifloats` with
// matching declarations in multifloats.hh. The `extern "C" *dd` shims
// are thin marshaling wrappers around those same C++ functions — no
// parallel `_full`-named body. Truly internal helpers (range reducers,
// polynomial evaluators, triple-double-output variants) stay in anon
// namespaces inside each .inc file.

#include "multifloats.hh"
#include "multifloats_td.hh"
#include "dd_constants.hh"
#include <cstdint>
#include <cstring>
#include <vector>

// Triple-double primitive bodies (declarations in multifloats_td.hh).
namespace multifloats {
namespace detail {
#include "multifloats_math_td.inc"
} // namespace detail
} // namespace multifloats

// Public kernels — definitions. Matching declarations live in multifloats.hh.
// Internal helpers that happen to share namespace scope (range reducers,
// polynomial Estrin kernels, TD variants) aren't header-declared; they leak
// as `multifloats::NAME` symbols but nobody outside this TU calls them.
// Once the symbol stripper is re-enabled in a follow-up commit, we'll
// tuck them behind `multifloats::detail::` or anon-namespace wrappers.
namespace multifloats {
using namespace multifloats::detail;  // neval, deval, horner, float64x3
#include "multifloats_math_exp_log.inc"
#include "multifloats_math_trig.inc"
#include "multifloats_math_hyp.inc"
#include "multifloats_math_inv_trig.inc"
#include "multifloats_math_special.inc"
#include "multifloats_math_bessel.inc"
} // namespace multifloats

// Complex overloads for Kind-D parity (sinpi/cospi/expm1/log2/log10/log1p
// on std::complex<float64x2>). Pulled in after <complex> is visible via
// the inclusion of multifloats.hh at top.
namespace multifloats {
#include "multifloats_math_abi_complex_kernels.inc"
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

namespace {
#include "multifloats_math_matmul.inc"
} // anonymous namespace

extern "C" {
#include "multifloats_math_abi_scalar.inc"
#include "multifloats_math_abi_complex.inc"
} // extern "C"
