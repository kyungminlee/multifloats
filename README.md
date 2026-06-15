# multifloats

Double-double arithmetic for Fortran and C++.

`multifloats` provides two derived/template types,

| Fortran type | C / C++ type    | Storage             | Precision        |
| ------------ | --------------- | ------------------- | ---------------- |
| `real64x2`   | `float64x2`     | 2 × `real(dp)`      | ~106 bits (~32 decimal digits) |
| `cmplx64x2`  | `complex64x2`   | 2 × `float64x2`     | ~106 bits per component |

implemented natively from error-free transformations on regular IEEE
doubles. There are no quad-precision (`real(16)` / `__float128`) temporaries
in any arithmetic, transcendental, or array kernel — those are reserved
exclusively for the `to_qp` / `from_qp` conversion helpers and the defined
I/O routines.

The Fortran module ships in `fsrc/multifloats.fypp` (preprocessed via
[fypp]) and the unified C/C++ header in `include/multifloats/float64x2.h`
(with `include/multifloats.h` kept as a compatibility shim that simply
includes it). Both expose the same algorithmic surface so the same kernels
can be used from either language.

## Features

### Fortran (`fsrc/multifloats.fypp` → `multifloats` module)

- **Types** — `real64x2` and `cmplx64x2`. Both carry the `SEQUENCE`
  attribute so they can appear in `EQUIVALENCE` statements (required by
  some LAPACK routines such as `wlaln2` / `DLALN2`). The bind(c)
  interop types `float64x2` / `complex64x2` (matching the C struct
  tags) live in the same module for the C-ABI bridge.
- **Operators** — `+`, `-`, `*`, `/`, `**`, `==`, `/=`, `<`, `>`, `<=`,
  `>=` between every combination of {`real64x2`, `cmplx64x2`,
  `real(dp)`, `real(sp)`, `integer`, `complex(dp)`, `complex(sp)`}.
- **Constructors and assignment** — every supported numeric kind, in
  both directions, including the non-default integer kinds
  `integer(int8)`, `integer(int16)`, `integer(int64)`. Identity
  constructors `real64x2(real64x2)` and `cmplx64x2(cmplx64x2)`
  let generic code call the constructor regardless of input type.
- **`<cmath>`-style intrinsics** — `abs`, `sqrt`, `cbrt`, `exp`, `exp2`,
  `expm1`, `log`, `log2`, `log10`, `log1p`, `sin`, `cos`, `tan`, `asin`,
  `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`,
  `atanh`, `erf`, `erfc`, `erfc_scaled`, `gamma`, `log_gamma`, `pow`,
  `hypot`, `fmod`/`mod`, `modulo`, `dim`, `sign`, `min`, `max` (3..8
  arguments), `bessel_j0/j1/y0/y1/jn/yn`, `aint`, `anint`, `floor`,
  `ceiling`, `nint`, `int`, `dble`, `fraction`, `spacing`, `rrspacing`,
  `scale`, `set_exponent`, `exponent`, `nearest`, `epsilon`, `huge`,
  `tiny`, `precision`, `range`, `digits`, `radix`, `minexponent`,
  `maxexponent`, `storage_size`.
- **Reductions** — `sum`, `product`, `maxval`, `minval`, `maxloc`,
  `minloc`, `findloc`, `dot_product`, `norm2`, `matmul` for both
  `real64x2` and `cmplx64x2`, supporting all rank-1..7 forms with
  `dim=`, `mask=`, and `back=` arguments.
- **Other** — `random_number` (rank 0..7), defined formatted I/O for
  both types, and inquiry functions matching `real(dp)` semantics on the
  leading limb.

The `real64x2` interface is designed to mirror `REAL(KIND=16)` so that
existing quad-precision code can switch types with minimal changes.

### C / C++ (`include/multifloats/float64x2.h`)

A single C++17 header in the `multifloats` namespace, providing the
`float64x2` class with the same algorithmic kernels as the Fortran side.
`include/multifloats.h` is preserved as a thin compatibility shim that
forwards to it, so existing `#include "multifloats.h"` users keep working.

The full `<cmath>` double-name surface (`floor`, `ceil`, `trunc`, `round`,
`nearbyint`, `rint`, `lround`, `llround`, `lrint`, `llrint`, `frexp`,
`modf`, `scalbn`, `scalbln`, `ilogb`, `logb`, `nextafter`, `nexttoward`,
`copysign`, `fma`, `fmod`, `remainder`, `remquo`, `fdim`, `sqrt`, `cbrt`,
`hypot`, `pow`, `exp`, `exp2`, `expm1`, `log`, `log10`, `log2`, `log1p`,
`sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`,
`tanh`, `asinh`, `acosh`, `atanh`, `erf`, `erfc`, `tgamma`, `lgamma`,
`abs`, `fabs`, `fmin`, `fmax`, `fpclassify`, `isfinite`, `isinf`, `isnan`,
`isnormal`, `signbit`, `isgreater`, `isgreaterequal`, `isless`,
`islessequal`, `islessgreater`, `isunordered`) is implemented via ADL on
`float64x2`. The C++ side has zero external dependencies — `__float128`
is only used by the test harness as a high-precision reference.

## Precision

The `cpp_fuzz` / `fortran_fuzz` drivers run every operation through 1M random
inputs (adversarial strategies: subnormals, near-cancellation, near-overflow,
non-finite leading limbs) and report the per-op `(max_rel_err, mean_rel_err)`.
The C++ driver references `__float128`; the Fortran driver references
`real(16)`; and the optional `cpp_fuzz_mpfr` driver references a 200-bit
[MPFR](https://www.mpfr.org/) oracle, which separates the DD kernel's own error
from the ~113-bit float128 floor. A relative error of `0` means *exactly
bit-equal* to the reference for every input drawn.

One **DD ulp** is `2⁻¹⁰⁴ ≈ 5e-32`. "Full double-double" means a worst-case
relative error within a small constant of one DD ulp (mean error roughly
1/100 ulp). The numbers below are from a 1M-iteration run (seed 0); your build
lands in the same orders of magnitude.

### Full double-double (~1e-32, ≈ one DD ulp)

Nearly the whole surface — arithmetic, the polynomial `exp` / `log` families
ported from [MultiFloats.jl](https://github.com/dzhang314/MultiFloats.jl), and
every elementary transcendental.

| Op | max_rel | mean_rel |
| --- | --- | --- |
| `+`, `-` | 1.8e-32 | 1.2e-33 |
| `*`, `/` | 4.8e-32, 8.1e-32 | ~2.5e-33 |
| `sqrt`, `rsqrt`, `hypot` | 3.4e-32 – 9.2e-32 | ~4e-33 |
| `cbrt` | 4.3e-31 | 1.5e-32 |
| `exp`, `exp2`, `exp10`, `expm1`, … | 3.2e-32 – 5.5e-32 | ~3e-33 |
| `log`, `log2`, `log10`, `log1p`, … | 2.2e-32 – 6.7e-32 | ~4e-33 |
| `pow` (dd\*\*dd / dp), `powr`, `rootn` | ~5e-32 | ~5e-33 |
| `pow` (integer exponent) | 2.4e-31 | 8.4e-33 |
| `sin`, `cos`, `tan` | 3.7e-32 – 6.5e-32 | ~2e-33 |
| `sinpi`, `cospi`, `tanpi` (vs 200-bit MPFR) | 2.3e-32 – 4.7e-32 | ~3e-33 |
| `asin`, `acos`, `atan`, `atan2` | 1.4e-32 – 3.6e-32 | ~1.5e-33 |
| `sinh`, `cosh`, `tanh` | 3.2e-32 – 6.1e-32 | ~3e-33 |
| `asinh`, `acosh`, `atanh` | 3.1e-32 – 5.7e-32 | ~2e-33 |
| `erf`, `erfc`, `erfc_scaled` | 1.7e-32 – 5.5e-32 | ~2e-33 |
| `log_gamma` | 4.7e-32 | 6.4e-33 |
| `bessel_*` (`j0/j1/jn`, `y0/y1/yn`; vs 200-bit MPFR) | 5e-33 – 6e-32 | ~8e-33 |
| Complex `+ − × ÷`, and the `cdd_*` transcendentals (exp, log, sqrt, the trig / hyperbolic / inverse families, `sinpi`/`cospi`, `expm1`) | ~1e-32 – 7e-32 | ~1e-32 |
| `sum`, `dot_product`, `norm2`, `matmul` (n=8) | 2e-31 – 7e-30 | ~1e-32 |
| `product` (n=8) | 4.0e-50 | 4.0e-53 |

The three-argument `fma` is full DD on average (mean 2.9e-33) with a ~1.8e-30
worst case. Reduction error grows ~linearly with the element count `n`.

### Bit-exact (always 0)

`abs`, `neg`, `sign`, `aint`, `anint`, `trunc`, `round`, `floor`, `ceil`,
`nearbyint`, `rint`, `roundeven`, `lround`/`llround`/`lrint`/`llrint`,
`logb`/`ilogb`/`llogb`, `fraction`, `scale`, `set_exponent`, `fmin`/`fmax`
and the C23 `fmaximum*`/`fminimum*` family, `copysign`,
`scalbn`/`ldexp`/`scalbln`, `frexp`, `modf`, `min`/`max` (3..8 args), every
constructor and assignment, and the complex `conjg` / `proj` / `*`-real-part.

### Reduced precision

A handful of kernels fall short of full DD. None is near the single-double
floor — the worst case is ~1e-26.

Numbers here are against the **200-bit MPFR** oracle (the honest kernel
measure) — the float128 fuzz inflates several of these well past their true
error (see the caveat below).

| Op | max_rel | mean_rel | Why |
| --- | --- | --- | --- |
| `gamma` | 1.4e-31 | 1.0e-32 | `exp(lgamma)` amplifies lgamma's absolute error (`log_gamma` is full DD) |
| `mod`, `modulo`, `remainder` | ~3e-23 | ~1e-27 | cancellation when the remainder is near zero; full DD otherwise |
| `cdd_pow` | 1.8e-29 | ~5e-31 | `exp(w·log(z))` amplifies a large `w·log(z)` (gamma-like) |
| `cdd_div` (re), `cdd_log1p` (re) | 2.3e-31 / 4.6e-31 | ~1e-32 | cancellation in `(ac+bd)/(c²+d²)` and `½log((1+x)²+y²)` near `z→0` |

(Everything else — including the **whole Bessel family**, the **complex
inverse trig** `cdd_asin/acos/atan/asinh/acosh/atanh`, `cdd_sinpi`/`cospi`,
`cdd_expm1`, and the real `sinpi`/`cospi`/`tanpi` — is **full DD** against the
200-bit MPFR oracle. The float128 fuzz reports some of these at ~1e-26 – 1e-29:
that is the float128 *reference's* own floor (its 113-bit π differs from the
kernel's for π-scaled trig at large `x`; libquadmath loses precision on the
complex branch cuts), not kernel error.)

### What makes the full-DD kernels exact

- **Arithmetic** — error-free transformations: Knuth `two_sum`, Dekker
  `two_prod` + FMA, Karp/Markstein for `sqrt`.
- **`exp`** — 14-term DD polynomial in the 1/8th-reduced argument, cubed via
  three squarings (MultiFloats.jl).
- **`log` / `log10`** — 32-entry table indexed by the top 5 mantissa bits plus
  a narrow polynomial.
- **`pow`** — `exp(b·log(a))`; full DD because both halves are.
- **`atan2`** — full-DD `atan(y/x)` with a quadrant correction from
  high-precision DD `π` / `π/2`.
- **`erf` / `erfc`** — piecewise rational fits with DD coefficients (ported
  from libquadmath `erfq.c`), with an overflow-safe asymptotic split
  `exp(−x²) = exp(−s²−0.5625)·exp((s−x)(s+x)+R)` where `s` is `x` truncated to
  35 mantissa bits so `s²` is exact.
- **`erfc_scaled`** — avoids forming `exp(x²)` for `|x| ≥ 1.25`; the
  large-negative reflection `2·exp(x²) − erfc_scaled(|x|)` squares `x` to
  triple-double and folds the residual limb into the exponent, since a DD `x²`
  cannot hold the argument to the precision `exp` needs at large `|x|`.
- **`bessel_*`** — the Hankel asymptotic phase is computed via the exact
  identity `sin(x − π/4) = (sin x − cos x)/√2` from a triple-double `sincos x`,
  rather than rounding the angle `x − π/4` to a DD (which would lose `~x·2⁻¹⁰⁴`
  and amplify the relative error near the zeros of `y0`/`y1`). The integer-order
  `jn`/`yn` recurrence runs in triple-double with a corrected `2k/x` coefficient,
  so it stays full DD near the zeros of `Y_n` too.
- **`sinpi` / `cospi` / `tanpi`** — reduce `x` mod ½ *exactly* (the
  integer/half-integer part is removed with no rounding), then scale the small
  remainder by the DD π — so the result is full DD even near the zeros at large
  `x`. (The float128 fuzz can't see this: its 113-bit π differs from the
  kernel's, so it reports a spurious ~1e-26 there; MPFR confirms full DD.)

Only `gamma` (~2.6e-31 at large arguments) and the inherently
cancellation-bound cases above (`mod`/`remainder` near a zero remainder, the
complex `…−1` formulas) remain short of full DD against a clean reference.

## Matmul API and GEMM relationship

The `matmul` surface provides three shape-dispatched operations:

| Op                | Signature                                    |
| ----------------- | -------------------------------------------- |
| `matmul(A, B)`    | `A(m,k) * B(k,n) → C(m,n)`                   |
| `matmul(A, x)`    | `A(m,k) * x(k)   → y(m)`                     |
| `matmul(x, B)`    | `x(k)   * B(k,n) → y(n)`                     |

Both real (`real64x2`) and complex (`cmplx64x2`) operands are
supported through the same Fortran generic. The C ABI exposes the three
kernels directly as `matmuldd_mm` / `matmuldd_mv` / `matmuldd_vm`
(complex matmul is composed in Fortran from four real matmuls).

**Scope vs BLAS GEMM.** `multifloats::matmul` is *not* a drop-in
replacement for DGEMM / CGEMM:

- **No `transa` / `transb` flags.** Inputs are always interpreted in
  storage order. Callers who need a transposed operand must
  materialize the transpose first (`matmul(transpose(A), B)` in
  Fortran, or pre-swap indices when setting up the buffer).
- **No `alpha` / `beta` scaling.** The output is overwritten, not
  accumulated into. `C := α·A·B + β·C` must be built by scaling an
  operand and combining with the previous `C` externally.
- **No leading dimension (LDA/LDB/LDC).** Storage is contiguous
  column-major with leading dimension equal to the first extent.

These constraints are deliberate: the compensated DD kernels use a
register-blocked panel design that assumes contiguous column-major
storage with the canonical shape. Relaxing them (e.g. GEMM-style
trans/alpha/beta/LDA) is tracked under "Deferred / out-of-scope work"
in `doc/developer/INTERNALS.md` for a future release; it requires a
new set of panel dispatchers to cover the additional shapes.

**BLAS shims.** `blas/wgemm.f90` and `blas/wtrsm.f90` provide
`real(qp)`-mangled DGEMM/DTRSM-style wrappers that route to the DD
kernels for the non-transposed cases (for use with LAPACK routines that
need quad-precision substitutes).

**Renormalization interval.** Matmul and `dot_product` run a
compensated fused-multiply-accumulate loop. The low-limb accumulator
grows by ~1 ULP per iteration, so for large inner dimension `k` it must
be re-normalized periodically. The module constant
`DD_FMA_RENORM_INTERVAL` (default `8`) controls the frequency; set `0`
to disable periodic renorm (single renorm at the end only). A 4×k·k×4
matmul sweep shows:

| k     | ri=0      | ri=8      | ri=64     |
| ----- | --------- | --------- | --------- |
| 64    | 2.9e-31   | 2.9e-31   | 2.9e-31   |
| 1024  | 3.2e-29   | 4.9e-30   | 3.3e-30   |
| 65536 | 1.5e-28   | 3.2e-30   | 1.5e-29   |

`ri=8` is near-optimal at ~3-4% overhead. Call
`dd_set_fma_renorm_interval(n)` from Fortran (or pass `n` to the
`matmuldd_*` C kernels) to override per-call. For `k < 100`, `ri=0`
is fine; for `k > 10000`, keep `ri=8` for best accuracy; avoid
`ri > 32` at very large k since `s_lo` can alias off significant bits
before the next renormalization fires.

## Error handling

`multifloats` follows a strict *NaN-in-NaN-out* policy. Invalid inputs
propagate through the computation as IEEE 754 non-finite values rather
than interrupting control flow:

- **No `errno`.** `<cerrno>` is never read or written. `log(-1)` returns
  a DD NaN; `errno` is left untouched.
- **No `fenv` side effects.** The kernels do not call `feraiseexcept`,
  `feclearexcept`, or depend on the floating-point rounding mode.
- **No exceptions.** No C++ `throw`, no `std::terminate`, no
  `__builtin_trap`. No Fortran `error stop`.
- **No signaling NaN handling.** Any NaN (quiet or signalling) is
  treated as data, and arithmetic kernels propagate it.
- **No input validation.** Domain checks are limited to what IEEE 754
  naturally produces — e.g. `sqrt(-x)` on a negative DD returns NaN
  because the underlying `std::sqrt(hi)` returns NaN. There is no
  pre-check that raises an error.

This mirrors the behavior of `double` in C and `real(kind=16)` in
Fortran and makes the kernels safe to call from hot loops, vectorized
code, and parallel regions without synchronization on shared error
state.

## Minimal examples

### Fortran

```fortran
use multifloats                    ! generic sqrt/atan overloads live here
integer, parameter :: dp = 8
type(real64x2) :: x, y
x = real64x2(1.0_dp)               ! from a dp literal
y = sqrt(atan(x) * 4)              ! sqrt(pi) to full DD
print *, y                         ! defined I/O: ~32 digits
```

### C (via the C ABI)

```c
#include <multifloats/float64x2.h>  /* linked with -lmultifloats */
float64x2 x = {1.0, 0.0};
float64x2 pi4 = atandd(x);          /* pi/4 as a DD */
float64x2 pi = muldd((float64x2){4.0, 0.0}, pi4);
printf("pi.hi = %.17g  pi.lo = %.17g\n", pi.limbs[0], pi.limbs[1]);
```

### C++

```cpp
#include <multifloats/float64x2.h>  // same header; C++ section adds the class API
using namespace multifloats;
float64x2 x(1.0);
float64x2 pi = float64x2(4.0) * atan(x);
std::cout << to_string(pi, 32) << "\n";   // scientific, 32 digits
```

Link `libmultifloats.a` (C / C++) or `libmultifloatsf-<compiler>.a`
(Fortran); see the [Building](#building) section.

## Building

Requires:
- CMake ≥ 3.27
- A Fortran 2018 compiler with `REAL(KIND=16)` support and a C++17
  compiler with `__float128` / libquadmath. On macOS this means Homebrew
  GCC 13/14/15 (Apple Clang and Apple-shipped LLVM Flang are not
  sufficient). The build pins `g++-15` / `gfortran-15` automatically; pass
  `-DCMAKE_CXX_COMPILER=...` / `-DCMAKE_Fortran_COMPILER=...` to override.
- [`fypp`](https://fypp.readthedocs.io/) on `PATH` (the Fortran source is
  generated at configure time).

```sh
cmake -B build -S .
cmake --build build
```

A default build produces only the two installable libraries:

- `libmultifloats.a` (or `libmultifloats-<compiler>.a` when LTO is on —
  see `MULTIFLOATS_USE_LTO` below) — the C++ kernels (header-only API via
  `include/multifloats/float64x2.h`; this static archive holds the
  out-of-line math bodies and the extern "C" `*dd` entry points declared
  in the same header).
- `libmultifloatsf-<compiler>.a` — the Fortran module library (the
  compiler tag comes from `cmake/FortranCompiler.cmake`; the generated
  `.mod` files live under `build/fmod/`).

Both targets are exported as a CMake package; consumers should prefer
`find_package(multifloats)` / `find_package(multifloatsf)` over hard-coded
`-l` flags so the compiler-tagged archive name is resolved automatically.

Everything else — tests, benchmarks, the `blas-multifloat` smoke
target, MPFR / Boost comparison harnesses — is opt-in via CMake
options. None of them are pulled in by a default consumer build.

| Option                              | Default | What it adds |
|-------------------------------------|---------|--------------|
| `-DBUILD_TESTING=ON`                | OFF     | C++ + Fortran test/fuzz executables, ctest registrations, and the `libblas-multifloat.a` smoke target (`wgemm`, `wtrsm` BLAS shims for `real64x2` matrices). |
| `-DMULTIFLOATS_BUILD_BENCH=ON`      | OFF     | `cpp_bench`, `fortran_bench`, `fortran_bench_abi` micro-benchmarks. |
| `-DBUILD_MPFR_TESTS=ON`             | OFF     | `cpp_fuzz_mpfr` 3-way precision test (needs `libmpfr-dev`). Implies `BUILD_TESTING`. |
| `-DMULTIFLOATS_BUILD_BOOST_COMPARE=ON` | OFF  | `boost_dd_fuzz` / `boost_dd_bench` / `bjn_probe` against `boost::multiprecision::cpp_double_double` (fetches Boost ≥ 1.89 via FetchContent). |
| `-DMULTIFLOATS_USE_LTO=ON/OFF`      | ON      | LTO + fat-LTO objects on the installed archives. When ON, the C++ archive is named `libmultifloats-<compiler>.a` so per-compiler LTO builds can coexist; when OFF, it is the portable untagged `libmultifloats.a`. The generated `multifloatsConfig.cmake` selects the best match at `find_package` time. |
| `-DMULTIFLOATS_HIDDEN_VISIBILITY=ON/OFF` | ON | Apply `-fvisibility=hidden` + `-fvisibility-inlines-hidden` to the C++ kernels so only the `extern "C" dd_*` ABI is exported. Turn OFF if you need full C++ symbol visibility (debugging, profiling, or re-exporting through a downstream shared library). |
| `-DMULTIFLOATSF_INSTALL_PRECOMPILED_MOD=ON` | OFF | Additionally install a compiler-tagged precompiled `.mod` + tagged Fortran archive. The default ships the fypp-expanded `.f90` source under `<prefix>/share/multifloatsf/src/` and lets consumers compile the module with their own Fortran compiler at `find_package` time, sidestepping `.mod` format incompatibilities. |

To run the test suite:

```sh
cmake -B build -S . -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## Tests

```sh
ctest --test-dir build --output-on-failure
```

| Test                       | Language | What it covers |
| -------------------------- | -------- | -------------- |
| `precision_fortran`        | Fortran  | Targeted vs-quad precision checks for constructors, assignments, every arithmetic and reduction op, and edge-case sweeps (signed zero, infinities, NaN propagation, subnormal/huge boundary, ULP boundary, dim/mask/back reduction variants). |
| `fuzz_fortran`             | Fortran  | 1M random pairs through every public function for which random real input is meaningful. Adversarial input strategies cover subnormals, near-cancellation, overflow boundary, and non-finite limbs. Prints a per-op `(max_rel, mean_rel, count)` precision report at the end. |
| `precision_fortran_unit`   | Fortran  | Hand-written assertions for arithmetic, signed zero, NaN/Inf propagation, classification, math intrinsics, and rounding. |
| `precision_abi_equivalence`| Fortran  | Cross-checks the three DD entry paths (native Fortran ops / C wrapper via `multifloats.h` / hand-written `bind(c)` reimpl in `test/dd_bindc.f90`) produce bit-identical results. |
| `precision_cpp`            | C++      | Targeted vs-`__float128` precision checks for `include/multifloats.h`. |
| `fuzz_cpp`                 | C++      | C++ fuzz with the same precision-report machinery as the Fortran fuzz. |
| `fuzz_cpp_determinism` / `fuzz_fortran_determinism` | shell | Diff two runs of the fuzz binaries to catch non-deterministic state. |
| `dd_constants_up_to_date`  | Python   | Re-runs `scripts/gen_constants.py --check` to detect drift between `src/dd_constants.hh` and the generator. |
| `fortran_abi_sync`         | shell    | `scripts/check_fortran_abi_sync.sh`: every `bind(c, name=*dd*)` in the generated Fortran module must match a `MULTIFLOATS_API` entry in `multifloats.h`. |

### Optional: MPFR-based 3-way precision test

`-DBUILD_MPFR_TESTS=ON` enables `cpp_fuzz_mpfr` (registered as ctest
`precision_mpfr_cpp`). It compares each operation against an
arbitrary-precision [MPFR](https://www.mpfr.org/) reference at 200
bits, reporting — per op — both `rel_err(libquadmath vs mpreal)` and
`rel_err(multifloats DD vs mpreal)` side by side. This separates the
~113-bit float128 floor from the DD kernel's own error, which the
default `__float128`-based fuzz cannot distinguish.

Dependencies: `libmpfr-dev` (system package). The
[mpreal](https://github.com/advanpix/mpreal) C++ header is used if
installed system-wide; otherwise it is fetched via CMake FetchContent.

```sh
cmake -DBUILD_MPFR_TESTS=ON -S . -B build
cmake --build build --target cpp_fuzz_mpfr
build/cpp_fuzz_mpfr 10000 42   # iterations, seed
```

## Layout

```
fsrc/multifloats.fypp                  -- Fortran source (fypp template)
include/multifloats/float64x2.h        -- unified C / C++ public header
include/multifloats.h                  -- compatibility shim (forwards to the header above)
src/                                   -- C++ .cc sources and implementation-detail .inc fragments
blas/                                  -- BLAS shims for real64x2
bench/                                 -- microbenchmarks (opt-in via MULTIFLOATS_BUILD_BENCH)
test/                                  -- Fortran and C++ test suites
external/                              -- vendored references (MultiFloats.jl, LAPACK)
```

The C++ kernels follow the algorithms in
[Julia's MultiFloats.jl](https://github.com/dzhang314/MultiFloats.jl) (vendored
under `external/`); the Fortran kernels are direct translations of the same
double-double algorithms.

## License

See [LICENSE](LICENSE) (if present) or the source headers for licensing terms.

[fypp]: https://fypp.readthedocs.io/
