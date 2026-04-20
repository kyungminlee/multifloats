// Umbrella translation unit for the DD math kernels and their C-ABI
// wrappers. The implementation is split across a set of `.inc` files
// that are pulled in here as a single compiled TU; that way the
// compiler keeps every cross-kernel call inlineable (e.g.
// `cexpdd → exp_full + sincos_full`) without depending on -flto.
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

#include "multifloats.hh"
#include "dd_constants.hh"
#include <cstdint>
#include <cstring>
#include <vector>

// All DD math kernels have internal linkage (anonymous namespace).
// Only the extern "C" `*dd` wrappers at the bottom are exported.
namespace {

using namespace multifloats::detail;  // neval, deval, horner
using multifloats::float64x2;

// Forward declarations for internal cross-references.
float64x2 exp_full(float64x2 const &x);
float64x2 exp2_full(float64x2 const &x);
float64x2 expm1_full(float64x2 const &x);
float64x2 log_full(float64x2 const &x);
float64x2 log1p_full(float64x2 const &x);
float64x2 sin_full(float64x2 const &x);
float64x2 cos_full(float64x2 const &x);
float64x2 sin_eval(float64x2 const &r);
float64x2 cos_eval(float64x2 const &r);
void sincos_eval(float64x2 const &r, float64x2 &s, float64x2 &c);
void sincos_full(float64x2 const &x, float64x2 &s, float64x2 &c);
void sinhcosh_full(float64x2 const &x, float64x2 &s, float64x2 &c);
float64x2 sinpi_full(float64x2 const &x);
float64x2 erfc_full(float64x2 const &x);
float64x2 lgamma_positive(float64x2 const &x);
float64x2 bessel_j0_full(float64x2 const &x);
float64x2 bessel_j1_full(float64x2 const &x);

#include "multifloats_math_exp_log.inc"
#include "multifloats_math_trig.inc"
#include "multifloats_math_hyp.inc"
#include "multifloats_math_inv_trig.inc"
#include "multifloats_math_special.inc"
#include "multifloats_math_bessel.inc"

} // anonymous namespace

// =============================================================================
// C-ABI entry points — extern "C" functions following math.h naming convention.
// These are the canonical DD implementations; the C++ template wrappers in
// multifloats.hh and the Fortran bind(C) interfaces both call these.
// =============================================================================

// multifloats_c.h is already pulled in via multifloats.hh.

namespace {
inline float64x2 from(float64x2_t x) { float64x2 r; r._limbs[0] = x.hi; r._limbs[1] = x.lo; return r; }
inline float64x2_t to(float64x2 const &x) { return {x._limbs[0], x._limbs[1]}; }

#include "multifloats_math_matmul.inc"
} // anonymous namespace

extern "C" {
#include "multifloats_math_abi_scalar.inc"
#include "multifloats_math_abi_complex.inc"
} // extern "C"
