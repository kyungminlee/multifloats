// AUTO-GENERATED — see comment in float64x2.h about the dual `dd` / `f64x2`
// naming surface. Each definition below is a thin wrapper that forwards
// to the existing `*dd` entry point under its `*f64x2` C23-style spelling.
//
// We previously used `__attribute__((alias))` for zero-overhead linker-level
// aliasing, but Mach-O (macOS) does not support strong aliases, and macOS
// GCC ICEs on the `weak, alias` form. Thin forwarders are portable across
// ELF/Mach-O/COFF; at -O2 (and certainly under LTO) the optimizer reduces
// them to a tail call, so the runtime cost is negligible.
//
// This file MUST be included once, inside the same `extern "C"` block as
// the originals (multifloats_math.cc handles that).
float64x2 addf64x2(float64x2 a, float64x2 b) { return adddd(a, b); }
float64x2 subf64x2(float64x2 a, float64x2 b) { return subdd(a, b); }
float64x2 mulf64x2(float64x2 a, float64x2 b) { return muldd(a, b); }
float64x2 divf64x2(float64x2 a, float64x2 b) { return divdd(a, b); }
float64x2 negf64x2(float64x2 a) { return negdd(a); }
float64x2 fabsf64x2(float64x2 a) { return fabsdd(a); }
float64x2 sqrtf64x2(float64x2 a) { return sqrtdd(a); }
float64x2 truncf64x2(float64x2 a) { return truncdd(a); }
float64x2 roundf64x2(float64x2 a) { return rounddd(a); }
float64x2 fminf64x2(float64x2 a, float64x2 b) { return fmindd(a, b); }
float64x2 fmaxf64x2(float64x2 a, float64x2 b) { return fmaxdd(a, b); }
float64x2 hypotf64x2(float64x2 a, float64x2 b) { return hypotdd(a, b); }
float64x2 powf64x2(float64x2 a, float64x2 b) { return powdd(a, b); }
float64x2 powif64x2(float64x2 a, int n) { return powidd(a, n); }
float64x2 fmodf64x2(float64x2 a, float64x2 b) { return fmoddd(a, b); }
float64x2 modulof64x2(float64x2 a, float64x2 b) { return modulodd(a, b); }
float64x2 fdimf64x2(float64x2 a, float64x2 b) { return fdimdd(a, b); }
float64x2 copysignf64x2(float64x2 a, float64x2 b) { return copysigndd(a, b); }
float64x2 fmaf64x2(float64x2 a, float64x2 b, float64x2 c) { return fmadd(a, b, c); }
float64x2 expf64x2(float64x2 a) { return expdd(a); }
float64x2 exp2f64x2(float64x2 a) { return exp2dd(a); }
float64x2 exp10f64x2(float64x2 a) { return exp10dd(a); }
float64x2 expm1f64x2(float64x2 a) { return expm1dd(a); }
float64x2 exp2m1f64x2(float64x2 a) { return exp2m1dd(a); }
float64x2 exp10m1f64x2(float64x2 a) { return exp10m1dd(a); }
float64x2 logf64x2(float64x2 a) { return logdd(a); }
float64x2 log2f64x2(float64x2 a) { return log2dd(a); }
float64x2 log10f64x2(float64x2 a) { return log10dd(a); }
float64x2 log1pf64x2(float64x2 a) { return log1pdd(a); }
float64x2 log2p1f64x2(float64x2 a) { return log2p1dd(a); }
float64x2 log10p1f64x2(float64x2 a) { return log10p1dd(a); }
float64x2 sinf64x2(float64x2 a) { return sindd(a); }
float64x2 cosf64x2(float64x2 a) { return cosdd(a); }
float64x2 tanf64x2(float64x2 a) { return tandd(a); }
float64x2 asinf64x2(float64x2 a) { return asindd(a); }
float64x2 acosf64x2(float64x2 a) { return acosdd(a); }
float64x2 atanf64x2(float64x2 a) { return atandd(a); }
float64x2 atan2f64x2(float64x2 a, float64x2 b) { return atan2dd(a, b); }
float64x2 sinpif64x2(float64x2 a) { return sinpidd(a); }
float64x2 cospif64x2(float64x2 a) { return cospidd(a); }
float64x2 tanpif64x2(float64x2 a) { return tanpidd(a); }
float64x2 asinpif64x2(float64x2 a) { return asinpidd(a); }
float64x2 acospif64x2(float64x2 a) { return acospidd(a); }
float64x2 atanpif64x2(float64x2 a) { return atanpidd(a); }
float64x2 atan2pif64x2(float64x2 a, float64x2 b) { return atan2pidd(a, b); }
float64x2 sinhf64x2(float64x2 a) { return sinhdd(a); }
float64x2 coshf64x2(float64x2 a) { return coshdd(a); }
float64x2 tanhf64x2(float64x2 a) { return tanhdd(a); }
float64x2 asinhf64x2(float64x2 a) { return asinhdd(a); }
float64x2 acoshf64x2(float64x2 a) { return acoshdd(a); }
float64x2 atanhf64x2(float64x2 a) { return atanhdd(a); }
float64x2 erff64x2(float64x2 a) { return erfdd(a); }
float64x2 erfcf64x2(float64x2 a) { return erfcdd(a); }
float64x2 erfcxf64x2(float64x2 a) { return erfcxdd(a); }
float64x2 tgammaf64x2(float64x2 a) { return tgammadd(a); }
float64x2 lgammaf64x2(float64x2 a) { return lgammadd(a); }
float64x2 j0f64x2(float64x2 a) { return j0dd(a); }
float64x2 j1f64x2(float64x2 a) { return j1dd(a); }
float64x2 y0f64x2(float64x2 a) { return y0dd(a); }
float64x2 y1f64x2(float64x2 a) { return y1dd(a); }
float64x2 jnf64x2(int n, float64x2 a) { return jndd(n, a); }
float64x2 ynf64x2(int n, float64x2 a) { return yndd(n, a); }
void jnf64x2_range(int n1, int n2, float64x2 a, float64x2 *out) { jndd_range(n1, n2, a, out); }
void ynf64x2_range(int n1, int n2, float64x2 a, float64x2 *out) { yndd_range(n1, n2, a, out); }
void sincosf64x2(float64x2 a, float64x2 *s, float64x2 *c) { sincosdd(a, s, c); }
void sinhcoshf64x2(float64x2 a, float64x2 *s, float64x2 *c) { sinhcoshdd(a, s, c); }
complex64x2 caddf64x2(complex64x2 a, complex64x2 b) { return cadddd(a, b); }
complex64x2 csubf64x2(complex64x2 a, complex64x2 b) { return csubdd(a, b); }
complex64x2 cmulf64x2(complex64x2 a, complex64x2 b) { return cmuldd(a, b); }
complex64x2 cdivf64x2(complex64x2 a, complex64x2 b) { return cdivdd(a, b); }
complex64x2 cexpf64x2(complex64x2 z) { return cexpdd(z); }
complex64x2 cexpm1f64x2(complex64x2 z) { return cexpm1dd(z); }
complex64x2 clogf64x2(complex64x2 z) { return clogdd(z); }
complex64x2 clog2f64x2(complex64x2 z) { return clog2dd(z); }
complex64x2 clog10f64x2(complex64x2 z) { return clog10dd(z); }
complex64x2 clog1pf64x2(complex64x2 z) { return clog1pdd(z); }
complex64x2 cpowf64x2(complex64x2 z, complex64x2 w) { return cpowdd(z, w); }
complex64x2 csqrtf64x2(complex64x2 z) { return csqrtdd(z); }
complex64x2 csinf64x2(complex64x2 z) { return csindd(z); }
complex64x2 csinpif64x2(complex64x2 z) { return csinpidd(z); }
complex64x2 ccosf64x2(complex64x2 z) { return ccosdd(z); }
complex64x2 ccospif64x2(complex64x2 z) { return ccospidd(z); }
complex64x2 ctanf64x2(complex64x2 z) { return ctandd(z); }
complex64x2 casinf64x2(complex64x2 z) { return casindd(z); }
complex64x2 cacosf64x2(complex64x2 z) { return cacosdd(z); }
complex64x2 catanf64x2(complex64x2 z) { return catandd(z); }
complex64x2 csinhf64x2(complex64x2 z) { return csinhdd(z); }
complex64x2 ccoshf64x2(complex64x2 z) { return ccoshdd(z); }
complex64x2 ctanhf64x2(complex64x2 z) { return ctanhdd(z); }
complex64x2 casinhf64x2(complex64x2 z) { return casinhdd(z); }
complex64x2 cacoshf64x2(complex64x2 z) { return cacoshdd(z); }
complex64x2 catanhf64x2(complex64x2 z) { return catanhdd(z); }
float64x2 cabsf64x2(complex64x2 z) { return cabsdd(z); }
float64x2 cargf64x2(complex64x2 z) { return cargdd(z); }
complex64x2 cprojf64x2(complex64x2 z) { return cprojdd(z); }
complex64x2 conjf64x2(complex64x2 z) { return conjdd(z); }
float64x2 crealf64x2(complex64x2 z) { return crealdd(z); }
float64x2 cimagf64x2(complex64x2 z) { return cimagdd(z); }
void matmulf64x2_mm(const float64x2 *a, const float64x2 *b, float64x2 *c, int64_t m, int64_t k, int64_t n, int64_t renorm_interval) { matmuldd_mm(a, b, c, m, k, n, renorm_interval); }
void matmulf64x2_mv(const float64x2 *a, const float64x2 *x, float64x2 *y, int64_t m, int64_t k, int64_t renorm_interval) { matmuldd_mv(a, x, y, m, k, renorm_interval); }
void matmulf64x2_vm(const float64x2 *x, const float64x2 *b, float64x2 *y, int64_t k, int64_t n, int64_t renorm_interval) { matmuldd_vm(x, b, y, k, n, renorm_interval); }
int eqf64x2(float64x2 a, float64x2 b) { return eqdd(a, b); }
int nef64x2(float64x2 a, float64x2 b) { return nedd(a, b); }
int ltf64x2(float64x2 a, float64x2 b) { return ltdd(a, b); }
int lef64x2(float64x2 a, float64x2 b) { return ledd(a, b); }
int gtf64x2(float64x2 a, float64x2 b) { return gtdd(a, b); }
int gef64x2(float64x2 a, float64x2 b) { return gedd(a, b); }
float64x2 cbrtf64x2(float64x2 a) { return cbrtdd(a); }
float64x2 ceilf64x2(float64x2 a) { return ceildd(a); }
float64x2 floorf64x2(float64x2 a) { return floordd(a); }
float64x2 nearbyintf64x2(float64x2 a) { return nearbyintdd(a); }
float64x2 rintf64x2(float64x2 a) { return rintdd(a); }
float64x2 logbf64x2(float64x2 a) { return logbdd(a); }
float64x2 frexpf64x2(float64x2 a, int *exp) { return frexpdd(a, exp); }
float64x2 modff64x2(float64x2 a, float64x2 *iptr) { return modfdd(a, iptr); }
float64x2 ldexpf64x2(float64x2 a, int n) { return ldexpdd(a, n); }
float64x2 scalbnf64x2(float64x2 a, int n) { return scalbndd(a, n); }
float64x2 scalblnf64x2(float64x2 a, long n) { return scalblndd(a, n); }
float64x2 nextafterf64x2(float64x2 a, float64x2 b) { return nextafterdd(a, b); }
float64x2 remainderf64x2(float64x2 a, float64x2 b) { return remainderdd(a, b); }
float64x2 remquof64x2(float64x2 a, float64x2 b, int *quo) { return remquodd(a, b, quo); }
int ilogbf64x2(float64x2 a) { return ilogbdd(a); }
long lroundf64x2(float64x2 a) { return lrounddd(a); }
long long llroundf64x2(float64x2 a) { return llrounddd(a); }
long lrintf64x2(float64x2 a) { return lrintdd(a); }
long long llrintf64x2(float64x2 a) { return llrintdd(a); }
int isnanf64x2(float64x2 a) { return isnandd(a); }
int isinff64x2(float64x2 a) { return isinfdd(a); }
int isfinitef64x2(float64x2 a) { return isfinitedd(a); }
int finitef64x2(float64x2 a) { return finitedd(a); }
int isnormalf64x2(float64x2 a) { return isnormaldd(a); }
int signbitf64x2(float64x2 a) { return signbitdd(a); }
int fpclassifyf64x2(float64x2 a) { return fpclassifydd(a); }
int isunorderedf64x2(float64x2 a, float64x2 b) { return isunordereddd(a, b); }
int isgreaterf64x2(float64x2 a, float64x2 b) { return isgreaterdd(a, b); }
int isgreaterequalf64x2(float64x2 a, float64x2 b) { return isgreaterequaldd(a, b); }
int islessf64x2(float64x2 a, float64x2 b) { return islessdd(a, b); }
int islessequalf64x2(float64x2 a, float64x2 b) { return islessequaldd(a, b); }
int islessgreaterf64x2(float64x2 a, float64x2 b) { return islessgreaterdd(a, b); }
int issignalingf64x2(float64x2 a) { return issignalingdd(a); }
float64x2 nanf64x2(const char *tagp) { return nandd(tagp); }
complex64x2 cexpif64x2(float64x2 a) { return cexpidd(a); }
