# Precision

The `cpp_fuzz` / `fortran_fuzz` drivers run every operation through 1M random
inputs — using adversarial strategies (subnormals, near-cancellation,
near-overflow, non-finite leading limbs) — and report the per-op
`(max_rel_err, mean_rel_err)`. The C++ driver references `__float128`
(libquadmath); the Fortran driver references `real(16)`; and the optional
`cpp_fuzz_mpfr` driver references a 200-bit [MPFR](https://www.mpfr.org/)
oracle, which separates the DD kernel's own error from the ~113-bit float128
floor. A relative error reported as `0` means *exactly bit-equal* to the
reference for every input drawn.

One **DD ulp** is `2⁻¹⁰⁴ ≈ 5e-32`, so "full double-double" means a worst-case
relative error that is a small constant multiple of `5e-32`.

The numbers below are from a representative 1M-iteration run (seed 0). Your
build will land in the same orders of magnitude.

## Precision classes at a glance

```{list-table}
:header-rows: 1

* - Class
  - Typical max rel err
  - Members
* - **Full double-double**
  - ~1e-32 (≈ 1 DD ulp)
  - `+ - * /`, `sqrt`, `cbrt`, `min`/`max`, `mod`, `dim`, `hypot`, `pow`,
    `exp`/`exp2`/`expm1`, `log`/`log2`/`log10`/`log1p`, the full trig set
    (`sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`), the full
    hyperbolic set (`sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`),
    `erf`, `erfc`, `erfc_scaled`, `lgamma`, complex `+ - * /`, and the
    `cdd_*` transcendentals.
* - **Bit-exact (always 0)**
  - 0
  - `abs`, `neg`, `sign`, `aint`, `anint`, `fraction`, `scale`,
    `set_exponent`, every constructor and assignment, complex `*` real part.
* - **Near full DD**
  - ~1e-31
  - `gamma` (a few DD ulp worst-case).
* - **Reduced — Bessel family**
  - ~1e-29 – 1e-27
  - `bessel_j0/j1/jn`, `bessel_y0/y1/yn`. libm-seeded; far better than one
    double ulp, but not full DD.
```

## Measured worst / mean relative error (1M, seed 0)

The transcendentals and special functions — the historically interesting
cases — all reach full DD except the Bessel family:

| Op | max_rel | mean_rel |
| --- | --- | --- |
| `sin` / `cos` / `tan` | 3.7e-32 / 4.1e-32 / 6.5e-32 | ~2e-33 |
| `atan` / `atan2` | 2.5e-32 | 1.3e-33 |
| `sinh` / `cosh` / `tanh` | ~6e-32 | ~3e-33 |
| `asinh` / `acosh` / `atanh` | ~5e-32 | ~2e-33 |
| `erf` | 1.7e-32 | 2.4e-33 |
| `erfc` | 2.5e-32 | 8.9e-34 |
| `erfc_scaled` (`erfcx`) | 5.5e-32 | 2.5e-33 |
| `lgamma` | 4.7e-32 | 6.4e-33 |
| `gamma` | 2.6e-31 | 1.0e-32 |
| `bessel_j0` / `j1` / `jn` | 1.1e-29 / 1.4e-29 / 1.4e-28 | ~1.5e-32 |
| `bessel_y0` / `y1` / `yn` | 7.9e-27 / 2.6e-29 / 3.6e-29 | ~3e-31 |

## Why the Bessel family is not full DD

`bessel_*` is seeded from libm's double-precision result and refined, with
`jn`/`yn` built on integer-order recurrences from `j0`/`j1`/`y0`/`y1`. The
worst case is `bessel_y0` near its zeros (~8e-27), where the function passes
through zero and the relative error is dominated by the libm seed. Reaching
full DD here would need dedicated DD polynomial / asymptotic kernels rather
than a libm-seeded refinement. `gamma` is a touch above the DD floor
(~2.6e-31) for the same reason at large arguments, while `lgamma` — computed
from a native DD Stirling series — is full DD.

## What makes the full-DD kernels exact

- **Arithmetic** uses error-free transformations: Knuth's `two_sum`, Dekker's
  `two_prod` + FMA, Karp/Markstein iteration for `sqrt`.
- **`exp`** — 14-term DD polynomial in the 1/8th-reduced argument, cubed via
  three squarings, ported from
  [MultiFloats.jl](https://github.com/dzhang314/MultiFloats.jl).
- **`log` / `log10`** — 32-entry lookup table indexed by the top 5 mantissa
  bits plus a narrow polynomial.
- **`pow`** — `exp(b · log(a))`, full DD because both halves are.
- **`atan2`** — full-DD `atan(y/x)` with a quadrant correction using
  high-precision DD `π` / `π/2` constants.
- **`erf` / `erfc`** — piecewise rational fits with DD coefficients (ported
  from libquadmath `erfq.c`), plus an overflow-safe asymptotic split
  `exp(-x²) = exp(-s²-0.5625)·exp((s-x)(s+x)+R)` where `s` is `x` truncated to
  35 mantissa bits, so `s²` is exact in a double.
- **`erfc_scaled`** — the asymptotic branch avoids forming `exp(x²)` for
  `|x| ≥ 1.25`. The large-negative reflection `2·exp(x²) − erfc_scaled(|x|)`
  squares `x` to **triple-double** and folds the residual limb back into the
  exponent (`exp(x²) = exp(x²_dd)·(1 + resid)`), since a DD `x²` cannot hold
  the argument to the absolute precision `exp` needs at large `|x|`.

```{important}
Every error-free transformation depends on `std::fma` being IEEE-compliant
(a single multiply-add with one rounding). The header refuses to compile under
`-ffast-math` / `-funsafe-math-optimizations`, which would silently collapse
the DD arithmetic to plain double precision.
```
