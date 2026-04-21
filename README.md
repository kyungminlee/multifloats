# multifloats

Double-double arithmetic for Fortran and C++.

`multifloats` provides two derived/template types,

| Type             | Storage             | Precision        |
| ---------------- | ------------------- | ---------------- |
| `float64x2`      | 2 × `real(dp)`      | ~106 bits (~32 decimal digits) |
| `complex64x2`   | 2 × `float64x2`     | ~106 bits per component |

implemented natively from error-free transformations on regular IEEE
doubles. There are no quad-precision (`real(16)` / `__float128`) temporaries
in any arithmetic, transcendental, or array kernel — those are reserved
exclusively for the `to_qp` / `from_qp` conversion helpers and the defined
I/O routines.

The Fortran module ships in `fsrc/multifloats.fypp` (preprocessed via
[fypp]) and the C++ header in `src/multifloats.hh`. Both expose the same
algorithmic surface so the same kernels can be used from either language.

## Features

### Fortran (`fsrc/multifloats.fypp` → `multifloats` module)

- **Types** — `float64x2` and `complex64x2`. Both carry the `SEQUENCE`
  attribute so they can appear in `EQUIVALENCE` statements (required by
  some LAPACK routines such as `wlaln2` / `DLALN2`).
- **Operators** — `+`, `-`, `*`, `/`, `**`, `==`, `/=`, `<`, `>`, `<=`,
  `>=` between every combination of {`float64x2`, `complex64x2`,
  `real(dp)`, `real(sp)`, `integer`, `complex(dp)`, `complex(sp)`}.
- **Constructors and assignment** — every supported numeric kind, in
  both directions, including the non-default integer kinds
  `integer(int8)`, `integer(int16)`, `integer(int64)`. Identity
  constructors `float64x2(float64x2)` and `complex64x2(complex64x2)`
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
  `float64x2` and `complex64x2`, supporting all rank-1..7 forms with
  `dim=`, `mask=`, and `back=` arguments.
- **Other** — `random_number` (rank 0..7), defined formatted I/O for
  both types, and inquiry functions matching `real(dp)` semantics on the
  leading limb.

The `float64x2` interface is designed to mirror `REAL(KIND=16)` so that
existing quad-precision code can switch types with minimal changes.

### C++ (`src/multifloats.hh`)

A single C++17 header in the `multifloats` namespace, providing
`MultiFloat<T, N>` (`N` ∈ {1, 2}) with the same algorithmic kernels and a
convenience alias `using float64x2 = MultiFloat<double, 2>;`.

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
`MultiFloat`. The C++ side has zero external dependencies — `__float128`
is only used by the test harness as a high-precision reference.

## Precision

The `cpp_fuzz` / `fortran_fuzz` drivers run every operation listed below
through 1M random input pairs (using adversarial input strategies that
include subnormals, near-cancellation pairs, near-overflow pairs, and
non-finite leading limbs) and report the per-op `(max_rel_err, mean_rel_err)`
against a quad-precision (`real(16)`) reference. The numbers in the tables
below are representative samples from a recent fuzz run (seed = 0, 1M
iterations); your build will land in the same orders of magnitude.

A relative error reported as `0` means *exactly bit-equal* to the qp
reference for every input the fuzz drew. Operations are organized by their
precision class.

### Full double-double (~1e-30 — 1e-32 max — close to one DD ulp)

The arithmetic kernels, the kernels whose result is a bit-exact rearrangement
of the input limbs, and the polynomial-based `exp` / `log` / `log10`
families ported from MultiFloats.jl. Worst-case error is a small constant
multiple of one DD ulp, mean error around 1/100 of a DD ulp.

| Op | max_rel | mean_rel |
| --- | --- | --- |
| `+` (dd±dd, dd±dp, etc.) | 1.5e-32 | 3.0e-34 |
| `-` | 6.2e-33 | 1.8e-34 |
| `*` | 3.5e-32 | 6.5e-34 |
| `/` | 5.3e-32 | 1.8e-33 |
| `sqrt` (Karp/Markstein iteration) | 5.0e-32 | 5.3e-33 |
| `mod`, `modulo` | 2.2e-32 | 9e-36 |
| `dim` | 6.2e-33 | 9.4e-35 |
| `min`, `max`, `min(a..h)`, `max(a..h)` | 6.2e-33 | 9.5e-35 |
| `hypot` (with overflow-safe scaling) | 7.8e-32 | 5.2e-33 |
| `exp` (14-term polynomial in 1/8th-reduced range, cubed) | 3.2e-30 | 2.2e-32 |
| `log`, `log10` (32-entry table + narrow polynomial) | 3.5e-32 | 3.3e-33 |
| `pow` (dd**dd, dd**dp, dp**dd via `exp(b·log(a))`) | 2.1e-30 | 5.3e-32 |
| `pow` (integer exponent, repeated multiplication) | 2.2e-32 | 1.3e-33 |
| `sinh` (Taylor for `|x|<0.1`, otherwise `(eˣ-e⁻ˣ)/2`) | 4.7e-30 | 4.1e-32 |
| `cosh` (`(eˣ+e⁻ˣ)/2`, well-conditioned) | 4.7e-30 | 4.8e-32 |
| `asin`, `acos` (one Newton step on full-DD `sin`/`cos`) | 3.7e-32 | 6.4e-33 |
| `acosh` (`log(x + sqrt(x²-1))` with large-`x` asymptotic) | 3.2e-32 | 2.0e-33 |
| `asinh` (Taylor for `|x|<0.01`, otherwise `log(x + sqrt(1+x²))`) | 4.3e-30 | 1.4e-32 |
| `atanh` (Taylor for `|x|<0.01`, otherwise `½ log((1+x)/(1-x))`) | 1.5e-30 | 1.1e-32 |
| `atan2` (full-DD `atan(y/x)` with quadrant correction) | 3.8e-32 | 1.5e-33 |
| Complex `+`, `-` (real and imag parts) | 1.4e-32 | 3.2e-34 |
| Complex `*` real part | 0 | 0 |
| Complex `*` imag part | 1.9e-32 | 1.5e-33 |
| Complex `/` real part | 4.5e-32 | 2.5e-33 |
| `cdd_sin`, `cdd_cos`, `cdd_sinh`, `cdd_cosh` (real and imag parts) | 1.0e-29 | 7e-32 |
| `cdd_tan`, `cdd_tanh` real / imag parts | 4e-30 | 3e-32 |
| `cdd_log` real and imag (overflow-safe formula) | 1.6e-31 | 2.7e-33 |
| `cdd_sqrt` real and imag (Kahan-style algorithm) | 6.6e-32 | 5.5e-33 |
| `cdd_atan`, `cdd_acos`, `cdd_acosh` (real and imag) | 4.3e-32 | 1.0e-32 |
| `cdd_asin_im`, `cdd_asinh_im` | 7.4e-32 | 1.8e-32 |
| `cdd_conjg`, `cdd_abs`, `cdd_aimag` | 6.7e-32 | 5.4e-33 |

### Near-DD precision (~1e-22 max, ~1e-25 mean)

| Op | max_rel | mean_rel | Why not full DD |
| --- | --- | --- | --- |
| `sin` | 9.2e-26 | 2.1e-28 | range reduction `x · 1/π` loses bits ∝ log₂|x| |
| `cos` | 2.3e-25 | 3.0e-28 | same |
| `tan` | 2.3e-25 | 5.2e-28 | same |
| `atan` | 1.1e-22 | 1.6e-25 | Newton step uses full-DD `sin`/`cos` but inherits their reduction precision |
| `tanh` | 4.0e-18 | 2.8e-21 | `1 - 2/(e²ˣ + 1)` near `|x| → ∞` rounds to 1 |

### Bit-exact (always 0 error)

Operations that are pure limb manipulation, sign flips, or trivial
promotions / truncations. The fuzz reports `max_rel = mean_rel = 0` over
1M iterations for every operation in this group.

- **Unary**: `abs`, `neg`, `sign`, `aint`, `anint`, `fraction`,
  `scale`, `set_exponent`
- **Mixed-mode arithmetic where one side is dp**: `dd + dp` (`add_fd`)
  reports 0 because the lo-limb error is in the dp ulp range
- **Every constructor**: `float64x2(...)` and `complex64x2(...)` for
  every supported numeric kind
- **Every assignment**: `dd ↔ {dp, sp, int, int8, int16, int64, cdp, csp}`
  and `cdd ↔ {dp, sp, int, int8, int16, int64, cdp, csp}`
- **Complex `*` real part**: bit-exact because the real part is computed
  as a single fma-style chain with no cancellation between terms

### Single-double, first-order derivative corrected (~1e-16 max)

`gamma` / `log_gamma` and the Bessel family are still computed as
`f(hi) + f'(hi) · lo` combined via `fast_two_sum`, which gives
roughly one double ulp of relative error — single-double precision.
There is no Julia polynomial port to crib from for these; getting
them to full DD would need dedicated polynomial / continued-fraction
approximations. `erfc_scaled` is in the same bucket.

`erf` and `erfc` have been upgraded to a hybrid kernel:

- **`|x| < 2`** — 50-term Taylor series with DD coefficients (full DD).
- **`|x| ≥ 2`** — `erf(x) = sign(x) · (1 − erfc_dp(|x|))` with the DD
  subtraction preserving the low limb down to the precision of libm's
  dp `erfc`. At `|x| = 2` that's ~5e-19; at `|x| ≥ 6` the low limb is
  below the DD ulp floor and the result is full DD.

Mean error across the fuzz drivers' random input distribution is
now ~2e-21 for `erf` (5 orders better than the previous
derivative-corrected version). Worst case is still ~2e-18 because
the `|x| ∈ [2, 6]` band inherits libm's dp precision for `erfc`.
Reaching full DD across the whole range would need a full-DD
asymptotic expansion (for `|x| > 6`, already trivial) and either a
many-term positive-series `e^(-x²) · Σ (2x²)^n / (2n+1)!!` for the
intermediate range or a much larger Taylor table.

| Op | max_rel | mean_rel |
| --- | --- | --- |
| `erf` | 2.2e-18 | 1.9e-21 |
| `erfc` | 5.0e-16 | 2.1e-18 |
| `erfc_scaled` | 4.8e-16 | 4.9e-17 |

### Compound — chained derivative corrections (~1e-12 to 1e-14 max)

Functions that internally chain two or more single-double-precision
transcendentals (`sin`, `cos`, `atan`, ...) or have cancellation in
their construction. Mean precision is still ~1e-15 but worst-case can
be a few orders looser.

| Op | max_rel | mean_rel | Notes |
| --- | --- | --- | --- |
| `gamma` | 1.2e-16 | 4.0e-17 | derivative correction unavailable; uses libm directly |
| `log_gamma` | 3.3e-16 | 5.0e-17 | likewise |
| `bessel_j0` | 3.9e-16 | 2.4e-17 | libm Bessel precision |
| `bessel_j1` | 2.5e-15 | 3.1e-17 | |
| `bessel_jn` (n=3 sample) | 5.0e-15 | 4.4e-17 | |
| `bessel_y0` | 2.3e-15 | 8.3e-17 | |
| `bessel_y1` | 7.9e-16 | 8.3e-17 | |
| `bessel_yn` (n=3 sample) | 1.4e-13 | 2.1e-16 | |
| `cdd_exp`, `cdd_sin`, `cdd_cos`, `cdd_sinh`, `cdd_cosh` (real and imag) | ~5e-29 max | ~1e-31 mean | full DD on average; near full DD worst-case (limited by `sin`/`cos` reduction) |
| `cdd_tan`, `cdd_tanh` (real and imag) | full DD | full DD | derived from `cdd_sin`/`cdd_cos` ratios |
| `cdd_div` imag part | 1.3e-16 | 5.0e-19 | fundamental cancellation in complex division (no fix) |
| `cdd_asin_re`, `cdd_asinh_re` | ~2e-24 max | ~3e-26 mean | bottlenecked by `atan2`'s precision floor (≈1e-22) |

### Array reductions (small-array fuzz, n = 8)

| Op | max_rel | mean_rel |
| --- | --- | --- |
| `sum`     | 4.0e-30 | 1.3e-32 |
| `product` | 1.1e-49 | 2.6e-52 |
| `maxval`  | 5.9e-33 | 1.5e-33 |
| `minval`  | 6.1e-33 | 1.4e-33 |
| `dot_product` | 1.2e-31 | 7.0e-33 |
| `norm2`   | 6.0e-32 | 1.0e-32 |
| `matmul` (mv, n=8) | 5.0e-31 | 6.6e-33 |

The reductions accumulate over `n` elements, so worst-case error grows
linearly with `n` while staying inside the full-DD regime.

### What's full DD vs not, and why

The full-DD kernels (~1e-32) are:

- **Arithmetic**: `+`, `-`, `*`, `/`, `sqrt`, `min`, `max`, `mod`, `modulo`,
  `dim`, `hypot`, `pow_int`. These use error-free transformations
  (Knuth's `two_sum`, Dekker's `two_prod` + FMA, Karp/Markstein for sqrt).
- **Bit-exact**: `abs`, `neg`, `sign`, `aint`, `anint`, `fraction`,
  `scale`, `set_exponent`, every constructor and assignment.
- **`exp`** — 14-term DD polynomial in the 1/8th-reduced argument,
  cubed via three squarings, ported from
  [Julia's MultiFloats.jl](https://github.com/dzhang314/MultiFloats.jl).
- **`log`, `log10`** — 32-entry lookup table indexed by the top 5
  mantissa bits, plus a polynomial in `t = (m - center)/(m + center)`.
  For `x ∈ [15/16, 17/16]` the kernel falls back to a direct polynomial
  in `t = (x - 1)/(x + 1)` with no table lookup. Ported from the same
  source.
- **`pow`** — `exp(b · log(a))`, full DD because both `exp` and `log`
  are full DD.
- **`sinh`, `cosh`** — derived from the new full-DD `exp`. `sinh` uses a
  9-term Taylor series with DD coefficients for `|x| < 0.1` (otherwise
  the `(eˣ - e⁻ˣ)/2` formula has cancellation), `cosh` uses
  `(eˣ + e⁻ˣ)/2` everywhere.
- **`asin`, `acos`** — one Newton step on the now-full-DD `sin`/`cos`,
  seeded by the libm leading-limb call. Quadratic convergence makes
  one step enough to reach DD precision.
- **`acosh`** — `log(x + sqrt(x² - 1))` directly, with a `log(2x)`
  asymptotic for huge `x` to avoid `x²` overflow.
- **`atan2`** — uses full-DD `atan(y/x)` with a quadrant correction
  using high-precision DD `π` and `π/2` constants. Picks the
  `atan(y/x)` or `±π/2 - atan(x/y)` form to keep `|argument| ≤ 1`.
- **Complex `+`, `-`, `*`, `/`, `conjg`, `abs`, `aimag`** — built from
  real DD ops directly, no transcendental dependencies.
- **`cdd_log`** real part — uses the overflow-safe formula
  `log(max(|a|,|b|)) + ½ log(1 + (min/max)²)` so that `|z| > huge`
  inputs don't overflow the intermediate `hypot`.
- **`cdd_sin`, `cdd_cos`, `cdd_sinh`, `cdd_cosh`, `cdd_tan`, `cdd_tanh`** —
  built from the new full-DD `sinh`/`cosh` and the near-DD
  `sin`/`cos`, so the leading-limb error of each component is
  ~1e-30.

The "near-DD" group (~1e-22 max, ~1e-25 mean): `sin`, `cos`, `tan`,
`atan`, `asinh`, `atanh`, `tanh`. These reach full DD on average but
have isolated worst-case inputs that drop toward single-double due to
range-reduction precision loss or formula cancellation. See the
near-DD precision table for the per-op explanations.

The single-double group (~1e-16): `erf`, `erfc`, `erfc_scaled`,
`gamma`, `log_gamma`, `bessel_*`. There is no Julia polynomial port
for these; getting them to full DD would need dedicated polynomial /
continued-fraction approximations.

The infrastructure (operators, reductions, complex / array support,
constructors, assignments) is already at full DD; only the per-function
kernels for the items above remain.

## Matmul API and GEMM relationship

The `matmul` surface provides three shape-dispatched operations:

| Op                | Signature                                    |
| ----------------- | -------------------------------------------- |
| `matmul(A, B)`    | `A(m,k) * B(k,n) → C(m,n)`                   |
| `matmul(A, x)`    | `A(m,k) * x(k)   → y(m)`                     |
| `matmul(x, B)`    | `x(k)   * B(k,n) → y(n)`                     |

Both real (`float64x2`) and complex (`complex64x2`) operands are
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
type(float64x2) :: x, y
x = float64x2(1.0_dp)              ! from a dp literal
y = sqrt(atan(x) * 4)              ! sqrt(pi) to full DD
print *, y                         ! defined I/O: ~32 digits
```

### C (via the C ABI)

```c
#include "multifloats_c.h"          /* linked with -lmultifloats */
float64x2_t x = {1.0, 0.0};
float64x2_t pi4 = atandd(x);        /* pi/4 as a DD */
float64x2_t pi = muldd((float64x2_t){4.0, 0.0}, pi4);
printf("pi.hi = %.17g  pi.lo = %.17g\n", pi.hi, pi.lo);
```

### C++

```cpp
#include "multifloats.hh"           // header-only public API
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
ctest --test-dir build --output-on-failure
```

The build produces:
- `libmultifloats.a` — the C++ kernels (header-only API via
  `src/multifloats.hh`; this static archive holds the out-of-line math
  bodies and the C-ABI entry points from `src/multifloats_c.h`).
- `libmultifloatsf-<compiler>.a` — the Fortran module library (the
  compiler tag comes from `cmake/FortranCompiler.cmake`; the generated
  `.mod` files live under `build/fmod/`).
- `libblas-multifloat.a` — `wgemm` / `wtrsm` BLAS shims that operate on
  `float64x2` matrices.
- Several test executables (see below).

## Tests

```sh
ctest --test-dir build --output-on-failure
```

| Test                       | Language | What it covers |
| -------------------------- | -------- | -------------- |
| `precision_fortran`        | Fortran  | Targeted vs-quad precision checks for constructors, assignments, every arithmetic and reduction op, and edge-case sweeps (signed zero, infinities, NaN propagation, subnormal/huge boundary, ULP boundary, dim/mask/back reduction variants). |
| `fuzz_fortran`             | Fortran  | 1M random pairs through every public function for which random real input is meaningful. Adversarial input strategies cover subnormals, near-cancellation, overflow boundary, and non-finite limbs. Prints a per-op `(max_rel, mean_rel, count)` precision report at the end. |
| `precision_fortran_unit`   | Fortran  | Hand-written assertions for arithmetic, signed zero, NaN/Inf propagation, classification, math intrinsics, and rounding. |
| `precision_abi_equivalence`| Fortran  | Cross-checks the three DD entry paths (native Fortran ops / C wrapper via `multifloats_c.h` / hand-written `bind(c)` reimpl in `test/dd_bindc.f90`) produce bit-identical results. |
| `precision_cpp`            | C++      | Targeted vs-`__float128` precision checks for `src/multifloats.hh`. |
| `fuzz_cpp`                 | C++      | C++ fuzz with the same precision-report machinery as the Fortran fuzz. |
| `fuzz_cpp_determinism` / `fuzz_fortran_determinism` | shell | Diff two runs of the fuzz binaries to catch non-deterministic state. |
| `dd_constants_up_to_date`  | Python   | Re-runs `scripts/gen_constants.py --check` to detect drift between `src/dd_constants.hh` and the generator. |
| `fortran_abi_sync`         | shell    | `scripts/check_fortran_abi_sync.sh`: every `bind(c, name=*dd*)` in the generated Fortran module must match a `MULTIFLOATS_API` entry in `multifloats_c.h`. |

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
fsrc/multifloats.fypp     -- Fortran source (fypp template)
src/multifloats.hh        -- C++17 header
blas/                     -- BLAS shims for float64x2
test/                     -- Fortran and C++ test suites
external/                 -- vendored references (MultiFloats.jl, LAPACK)
```

The C++ kernels follow the algorithms in
[Julia's MultiFloats.jl](https://github.com/dzhang314/MultiFloats.jl) (vendored
under `external/`); the Fortran kernels are direct translations of the same
double-double algorithms.

## License

See [LICENSE](LICENSE) (if present) or the source headers for licensing terms.

[fypp]: https://fypp.readthedocs.io/
