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
  - `+ - * /`, `fma`, `sqrt`, `cbrt`, `min`/`max`, `mod`, `dim`, `hypot`, `pow`,
    `exp`/`exp2`/`expm1`, `log`/`log2`/`log10`/`log1p`, the full trig set
    (`sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, and the π-scaled
    `sinpi`/`cospi`/`tanpi`), the full
    hyperbolic set (`sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`),
    `erf`, `erfc`, `erfc_scaled`, `lgamma`, the full Bessel family
    (`bessel_j0/j1/jn`, `bessel_y0/y1/yn`), complex `+ - * /`, and the
    `cdd_*` transcendentals.
* - **Bit-exact (always 0)**
  - 0
  - `abs`, `neg`, `sign`, `aint`, `anint`, `fraction`, `scale`,
    `set_exponent`, every constructor and assignment, complex `*` real part.
* - **Near full DD**
  - ~1e-31 – 3e-30
  - `gamma` (~1e-31) and the cancellation-bound `cdd_div`(re) (~3e-30) —
    a few to ~100 DD ulp.
* - **Reduced**
  - ~3e-29 – 1e-23
  - `cdd_log1p`(re) (~3e-29; `½log((1+x)²+y²)` cancellation near `z→0`), and
    `mod`/`modulo`/`remainder` (full DD normally — ~2e-32 — but ~1e-23 at very
    large arguments, ~2⁶⁵, where the integer-quotient reduction degrades).
```

## Measured worst / mean relative error (1M, seed 0)

The transcendentals and special functions — the historically interesting
cases — all reach full DD. The Bessel and π-scaled-trig columns list the
**200-bit MPFR** oracle (`max_dd`): the float128 reference's own precision
floor near the zeros of `y0` (and its 113-bit π differing from the kernel's
for `sinpi`/`cospi`/`tanpi` at large `x`, and its near-cancellation floor for
`fma`) limits *its* readings to ~1e-26 – 1e-30, so MPFR is the honest measure
of these kernels — `fma` itself is full DD (~4e-32).

| Op | max_rel | mean_rel |
| --- | --- | --- |
| `sin` / `cos` / `tan` | 3.7e-32 / 4.1e-32 / 6.5e-32 | ~2e-33 |
| `sinpi` / `cospi` / `tanpi` (MPFR) | 3.9e-32 / 2.3e-32 / 4.7e-32 | ~3e-33 |
| `atan` / `atan2` | 2.5e-32 | 1.3e-33 |
| `sinh` / `cosh` / `tanh` | ~6e-32 | ~3e-33 |
| `asinh` / `acosh` / `atanh` | ~5e-32 | ~2e-33 |
| `erf` | 1.7e-32 | 2.4e-33 |
| `erfc` | 2.5e-32 | 8.9e-34 |
| `erfc_scaled` (`erfcx`) | 5.5e-32 | 2.5e-33 |
| `lgamma` | 4.7e-32 | 6.4e-33 |
| `gamma` | 2.6e-31 | 1.0e-32 |
| `bessel_j0` / `j1` / `jn` (MPFR) | 5.9e-33 / 7.4e-33 / 1.1e-32 | ~3e-34 |
| `bessel_y0` / `y1` / `yn` (MPFR) | 4.9e-32 / 5.8e-32 / 5.9e-32 | ~8e-33 |

## The Bessel family

`bessel_*` evaluates the small-argument power series (`x ≤ 2`) or the Hankel
asymptotic (`x > 2`), with `jn`/`yn` built on integer-order recurrences from
`j0`/`j1`/`y0`/`y1`. The asymptotic phase `sin(x − π/4)` / `cos(x − π/4)` is
*not* formed by rounding the angle `x − π/4` to a DD — that would lose
`~x·2⁻¹⁰⁴` of the angle and blow up the relative error near `y0`'s zeros
(where `P·sin + Q·cos` cancels). Instead `sin x` / `cos x` are evaluated to
triple-double (proper `π/2` range reduction) and combined via the exact
identity `sin(x − π/4) = (sin x − cos x)/√2`. With this, all the kernels are
full DD (verified against the 200-bit MPFR oracle). `jn`/`yn` carry their
integer-order recurrence in triple-double — including a corrected `2k/x`
coefficient, whose plain-DD rounding would otherwise dominate near the zeros
of `Y_n` — so the recurrence does not accumulate above the DD floor. `gamma`
remains ~2.6e-31 at large arguments, while `lgamma` — a native DD Stirling
series — is full DD.

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
- **complex `pow`** — `exp(w · log z)` computes `log z` in **triple-double**
  (`log_td` / `atan2_td`, each one Newton step from the DD result using the TD
  `exp` / `sincos`), so `log z`'s absolute error no longer dominates after the
  complex exp's phase amplification.

```{important}
Every error-free transformation depends on `std::fma` being IEEE-compliant
(a single multiply-add with one rounding). The header refuses to compile under
`-ffast-math` / `-funsafe-math-optimizations`, which would silently collapse
the DD arithmetic to plain double precision.
```
