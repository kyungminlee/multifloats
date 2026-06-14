# Precision

The `cpp_fuzz` / `fortran_fuzz` drivers run every operation through 1M random
input pairs — using adversarial strategies (subnormals, near-cancellation,
near-overflow, non-finite leading limbs) — and report the per-op
`(max_rel_err, mean_rel_err)` against a quad-precision (`real(16)`) reference.
A relative error reported as `0` means *exactly bit-equal* to the qp reference
for every input drawn.

This page summarises which kernels reach **full double-double** precision and
which do not, and why. For the complete per-op tables see the
[benchmark report](../BENCHMARK.md) and the project README.

## Precision classes at a glance

```{list-table}
:header-rows: 1

* - Class
  - Typical max rel err
  - Members
* - **Full double-double**
  - ~1e-30 – 1e-32
  - `+ - * /`, `sqrt`, `min`/`max`, `mod`, `dim`, `hypot`, `pow_int`, `exp`,
    `log`, `log10`, `pow`, `sinh`, `cosh`, `asin`, `acos`, `acosh`, `atan2`,
    complex `+ - * /`, and the `cdd_*` transcendentals.
* - **Bit-exact (always 0)**
  - 0
  - `abs`, `neg`, `sign`, `aint`, `anint`, `fraction`, `scale`,
    `set_exponent`, every constructor and assignment, complex `*` real part.
* - **Near-DD**
  - ~1e-22 max, ~1e-25 mean
  - `sin`, `cos`, `tan`, `atan`, `asinh`, `atanh`, `tanh`.
* - **Single-double, derivative-corrected**
  - ~1e-16
  - `erf`, `erfc`, `erfc_scaled`, `gamma`, `log_gamma`, `bessel_*`.
```

## Why some kernels are not full DD

**Near-DD group** (`sin`, `cos`, `tan`, `atan`, `asinh`, `atanh`, `tanh`).
These reach full DD on average but have isolated worst-case inputs that drop
toward single-double — range reduction `x · 1/π` loses bits ∝ log₂|x|, and
formulas like `1 − 2/(e²ˣ+1)` round to 1 as `|x| → ∞`.

**Single-double group** (`erf`, `erfc`, `erfc_scaled`, `gamma`, `log_gamma`,
`bessel_*`). These are computed as `f(hi) + f'(hi)·lo` combined via
`fast_two_sum`, giving roughly one double ulp. There is no Julia polynomial
port to crib from for these; reaching full DD would need dedicated polynomial
or continued-fraction approximations.

`erf`/`erfc` use a hybrid kernel — a 50-term DD Taylor series for `|x| < 2`
(full DD) and `sign(x)·(1 − erfc_dp(|x|))` for `|x| ≥ 2`, which is full DD once
`|x| ≥ 6` and inherits libm's dp precision in the `[2, 6]` band.

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

```{important}
Every error-free transformation depends on `std::fma` being IEEE-compliant
(a single multiply-add with one rounding). The header refuses to compile under
`-ffast-math` / `-funsafe-math-optimizations`, which would silently collapse
the DD arithmetic to plain double precision.
```
