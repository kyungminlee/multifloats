# Benchmark Results

Comparison of multifloats double-double (DD) arithmetic against quad
precision (`real(16)` / `__float128` via libquadmath). The multifloats
implementations use explicit two-limb error-free transformations (EFTs)
on pairs of `double` values, achieving ~106 bits of significand — close
to quad precision's 113 bits — at a fraction of the cost.

## Systems

Short names are used as column labels throughout.

| Short name | CPU | OS | Compiler | Build |
|---|---|---|---|---|
| **M1 Max** | Apple M1 Max (arm64), 10 cores | macOS 26.3 (Darwin 25.3.0) | GNU Fortran 15.2.0 (Homebrew GCC 15.2.0_1) / g++ Apple clang version 17.0.0 (clang-1700.6.4.2) | CMake 4.3.1, `-O3 -flto`, STATIC library |

## Precision key

Precision is the maximum relative error vs an **mpreal @ 200-bit** oracle
over ~1M random inputs (fixed seed 42), reported in **DD ulps** (1 ulp ≈
2⁻¹⁰⁵ ≈ 2.46e-32 for the ~106-bit DD significand):

| Label | err (DD ulp) | Meaning |
|---|---|---|
| **full DD** | ~1 ulp | Full double-double EFT kernel (~106 bits) |
| **exact** | 0 | Bit-exact (no rounding involved) |
| **deriv-corrected** | ~1e7 to 1e14 ulp | `f(hi) + f'(hi)*lo` correction gives near-DD |
| **single-double** | ~1e15 to 1e17 ulp | Leading-limb libm call, no lo correction |

## Origin key

| Tag | Meaning |
|---|---|
| **Julia** | Ported from `external/MultiFloats.jl/src/` |
| **original** | Developed for this project |
| **sample** | Ad-hoc implementation in the benchmark harness, not exported from the library |

## Reference: quadmath vs MPFR

Baseline for reading the DD err numbers below. This table shows the
**quadmath leg's own** relative error against the same mpreal @ 200-bit
oracle, in **qp ulps** (1 qp ulp ≈ 2⁻¹¹² ≈ 1.93e-34 for the 113-bit
`__float128` significand). System-independent — qp precision is identical
on every machine.

Read the main tables as: *the DD kernel lands at N DD ulps vs the same
200-bit reference that puts quadmath at ~M qp ulps.* When the DD number
is close to its expected "full DD" band and the qp number is close to ½
qp ulp, both legs are behaving as specified.

| op | qp-vs-MPFR err |
|---|---|
| add | 0.5 qp ulp |
| sub | 0.5 qp ulp |
| mul | 0.5 qp ulp |
| div | 0.5 qp ulp |
| sqrt | 0.5 qp ulp |
| abs | exact |
| neg | exact |
| fmin | exact |
| fmax | exact |
| fdim | 0.5 qp ulp |
| copysign | exact |
| fmod | exact |
| hypot | 0.7 qp ulp |
| exp | 0.5 qp ulp |
| expm1 | 0.7 qp ulp |
| log | 0.5 qp ulp |
| log10 | 0.5 qp ulp |
| log1p | 0.6 qp ulp |
| pow | 0.5 qp ulp |
| sin | 0.5 qp ulp |
| cos | 0.3 qp ulp |
| sincos | 0.5 qp ulp (sin) / 0.3 qp ulp (cos) |
| sinpi | 1.9 qp ulp |
| cospi | 1.9 qp ulp |
| tan | 0.5 qp ulp |
| tanpi | 4330376 qp ulp |
| asin | 0.6 qp ulp |
| asinpi | 0.3 qp ulp |
| acos | 0.5 qp ulp |
| acospi | 0.9 qp ulp |
| atan | 0.7 qp ulp |
| atanpi | 0.3 qp ulp |
| atan2 | 0.6 qp ulp |
| atan2pi | 1.1 qp ulp |
| sinh | 1.2 qp ulp |
| cosh | 0.8 qp ulp |
| sinhcosh | 1.2 qp ulp (sinh) / 0.8 qp ulp (cosh) |
| tanh | 1.1 qp ulp |
| asinh | 1.2 qp ulp |
| acosh | 0.5 qp ulp |
| atanh | 1.1 qp ulp |
| erf | 0.7 qp ulp |
| erfc | 0.5 qp ulp |
| erfcx | 256 qp ulp |
| tgamma | 1.7 qp ulp |
| lgamma | 1.5 qp ulp |
| bessel\_j0 | 0.3 qp ulp |
| bessel\_j1 | 0.3 qp ulp |
| bessel\_jn(3,.) | 0.3 qp ulp |
| bessel\_y0 | 1.5 qp ulp |
| bessel\_y1 | 1.3 qp ulp |
| bessel\_yn(3,.) | 3.7 qp ulp |
| bessel\_yn\_range(0..5) | 3.1 qp ulp |
| cdd\_add | — |
| cdd\_sub | — |
| cdd\_mul | — |
| cdd\_div | — |
| cdd\_proj | — |
| cdd\_abs | — |
| cdd\_arg | — |
| cdd\_sqrt | — |
| cdd\_exp | — |
| cdd\_expm1 | — |
| cdd\_log | — |
| cdd\_log2 | — |
| cdd\_log10 | — |
| cdd\_log1p | — |
| cdd\_pow | — |
| cdd\_sin | — |
| cdd\_cos | — |
| cdd\_tan | — |
| cdd\_sinpi | — |
| cdd\_cospi | — |
| cdd\_sinh | — |
| cdd\_cosh | — |
| cdd\_tanh | — |
| cdd\_asin | — |
| cdd\_acos | — |
| cdd\_atan | — |
| cdd\_asinh | — |
| cdd\_acosh | — |
| cdd\_atanh | — |
| mul (dp\*dd) | 0.5 qp ulp |
| aint | — |
| anint | — |
| fraction | — |
| min | — |
| max | — |
| min3 | exact |
| max3 | exact |
| sign | — |
| dim | — |
| mod | — |
| modulo | — |
| pow\_int | 0.5 qp ulp |
| erfc\_scaled | — |
| cdd\_add | 0.5 qp ulp (re) / 0.5 qp ulp (im) |
| cdd\_sub | 0.5 qp ulp (re) / 0.5 qp ulp (im) |
| cdd\_mul | 0.9 qp ulp (re) / 0.9 qp ulp (im) |
| cdd\_div | 54 qp ulp (re) / 36 qp ulp (im) |
| cdd\_abs | 0.7 qp ulp |
| cdd\_sqrt | 1.0 qp ulp (re) / 0.9 qp ulp (im) |
| cdd\_exp | 1.2 qp ulp (re) / 1.1 qp ulp (im) |
| cdd\_expm1 | 1936 qp ulp (re) / 1.1 qp ulp (im) |
| cdd\_log | 1.1 qp ulp (re) / 0.9 qp ulp (im) |
| cdd\_log2 | 1.7 qp ulp (re) / 1.2 qp ulp (im) |
| cdd\_log10 | 1.2 qp ulp (re) / 1.2 qp ulp (im) |
| cdd\_log1p | 1129 qp ulp (re) / 1.1 qp ulp (im) |
| cdd\_pow | 112 qp ulp (re) / 4057 qp ulp (im) |
| cdd\_sin | 1.4 qp ulp (re) / 1.3 qp ulp (im) |
| cdd\_cos | 1.3 qp ulp (re) / 1.3 qp ulp (im) |
| cdd\_tan | 1.7 qp ulp (re) / 2.2 qp ulp (im) |
| cdd\_sinh | 1.5 qp ulp (re) / 1.4 qp ulp (im) |
| cdd\_cosh | 1.4 qp ulp (re) / 1.4 qp ulp (im) |
| cdd\_tanh | 2.2 qp ulp (re) / 1.9 qp ulp (im) |
| cdd\_asin | 1.3 qp ulp (re) / 1.5 qp ulp (im) |
| cdd\_acos | 1.1 qp ulp (re) / 1.5 qp ulp (im) |
| cdd\_atan | 1.2 qp ulp (re) / 2.1 qp ulp (im) |
| cdd\_asinh | 1.7 qp ulp (re) / 1.5 qp ulp (im) |
| cdd\_acosh | 1.5 qp ulp (re) / 1.1 qp ulp (im) |
| cdd\_atanh | 1.7 qp ulp (re) / 1.0 qp ulp (im) |
| arr\_sum (n=8) | — |
| arr\_product (n=8) | — |
| arr\_maxval (n=8) | — |
| arr\_minval (n=8) | — |
| arr\_dot (n=8) | — |
| arr\_norm2 (n=8) | — |
| arr\_matmul (8×8·8) | — |

## C: C ABI vs `__float128`

Ops timed through the **C ABI** entry points (`sindd`, `cdd_muldd`, `j0dd`,
…) — the register-return convention a C or Fortran-`bind(c)` caller sees,
not the C++ header-only inline path. Scope: every function in
`src/multifloats_c.h` (scalars, π-family, Bessel, complex arithmetic,
complex transcendentals). Shared ops also appear in the Fortran section
below under the Fortran elemental wrapper; the delta is the hidden-
pointer Fortran ABI overhead.

Each op is timed over 1024 elements × 400 repetitions (fast ops) or fewer
reps (transcendentals), with a NOINLINE drain after each rep to prevent
dead-code elimination. **×** = speedup (`qp_time / dd_time`, values > 1×
mean multifloats is faster); **err** = DD error in DD ulps from the 1M-
input MPFR fuzz run.

### Arithmetic

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| add | 0.7 ulp | **2.8×** | Julia: two\_sum EFT |
| sub | 0.7 ulp | **2.5×** | Julia: two\_sum EFT (negate + add) |
| mul | 2.0 ulp | **11×** | Julia: two\_prod EFT via FMA |
| div | 3.1 ulp | **7.5×** | original: Newton refinement (1/y seed, one step) |
| sqrt | 1.7 ulp | **44×** | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) |
| fma | — | **80×** | original: x\*y + z via DD ops |
| abs | exact | **2.1×** | original: sign-check + negate limbs |
| neg | exact | **2.6×** | original: negate both limbs |

### Rounding

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| trunc | exact | 1.7× | original: signbit ? −floor(−x) : floor(x) |
| round | exact | 1.9× | original: trunc(x + ½·sign(x)) |

### Binary

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| fmin | exact | **4.3×** | original: DD comparison + select |
| fmax | exact | **5.1×** | original: DD comparison + select |
| fdim | 0.7 ulp | **5.5×** | original: DD comparison, then subtract or zero |
| copysign | exact | 1.8× | original: sign-bit copy to hi, propagate to lo |
| fmod | 1.0 ulp | 0.96× | sample: floor-multiple reduction loop; fallback to div chain |
| hypot | 2.4 ulp | **10×** | original: scaled sqrt(x²+y²) |

### Exponential / logarithmic

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| exp | 1.4 ulp | **3.6×** | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction |
| exp2 | — | **4.9×** | Julia: exp2 polynomial (14-term Horner) |
| expm1 | 1.5 ulp | **5.5×** | original: exp(x) − 1 via DD sub |
| log | 1.4 ulp | **6.4×** | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) |
| log10 | 1.2 ulp | **8.7×** | Julia: log2 kernel × DD log10(2) |
| log2 | — | **8.3×** | Julia: log2 table lookup + polynomial |
| log1p | 3.7 ulp | **7.0×** | original: log(1 + x) via DD add |
| pow | 11 ulp | **5.7×** | Julia: exp(y × log(x)) |

### Trigonometric

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| sin | 1.4 ulp | **3.5×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| cos | 1.1 ulp | **3.4×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| sincos | 1.4 ulp (sin) / 1.1 ulp (cos) | **5.6×** | fused sin/cos: one range-reduction, two Taylor outputs |
| sinpi | 1.5 ulp | **4.5×** | Julia: sinpi Horner polynomial, direct |
| cospi | 1.4 ulp | **4.6×** | Julia: cospi Horner polynomial, direct |
| tan | 2.3 ulp | **2.5×** | original: sin/cos Taylor kernels + DD divide |
| tanpi | 1.9 ulp | **3.4×** | original: sinpi/cospi ratio |
| asin | 0.7 ulp | **6.0×** | original: piecewise rational P/Q (3 regions, from libquadmath asinq.c) |
| asinpi | 0.5 ulp | **6.0×** | original: asin(x)/π with exact-DD π division |
| acos | 0.6 ulp | **6.2×** | original: asin polynomial + half-angle identity |
| acospi | 1.2 ulp | **6.1×** | original: acos(x)/π with exact-DD π division |
| atan | 0.7 ulp | **4.1×** | original: 84-entry table lookup + rational P(t²)/Q(t²) (from libquadmath atanq.c) |
| atanpi | 0.6 ulp | **4.3×** | original: atan(x)/π with exact-DD π division |
| atan2 | 0.8 ulp | **3.5×** | original: table-based atan + quadrant correction |
| atan2pi | 1.4 ulp | **3.6×** | original: atan2(y,x)/π with exact-DD π division |

### Hyperbolic

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| sinh | 4.3 ulp | **3.1×** | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 |
| cosh | 1.0 ulp | **2.2×** | original: (exp+exp⁻¹)/2 |
| sinhcosh | 4.3 ulp (sinh) / 1.0 ulp (cosh) | **4.7×** | fused sinh/cosh: one range-reduction, two outputs |
| tanh | 4.1 ulp | **3.7×** | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) |
| asinh | 2.0 ulp | **7.2×** | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton |
| acosh | 0.8 ulp | **7.2×** | original: log(x+√(x²−1)) with Newton correction |
| atanh | 2.0 ulp | **6.4×** | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) |

### Error / special functions

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| erf | 0.7 ulp | **3.8×** | piecewise rational approx (ported from libquadmath erfq.c) |
| erfc | 0.5 ulp | **3.9×** | piecewise rational approx + split exp(-x^2) |
| erfcx | 328 ulp | **4.0×** | exp(x²)·erfc(x); scaled form avoiding tail cancellation |
| tgamma | 18 ulp | **8.5×** | piecewise rational approx + Stirling + reflection, exp(lgamma) |
| lgamma | 2.0 ulp | **4.8×** | piecewise rational approx + Stirling asymptotic |
| bessel\_j0 | 0.2 ulp | **6.5×** | piecewise rational + Hankel asymptotic (j0q.c) |
| bessel\_j1 | 0.4 ulp | **6.7×** | piecewise rational + Hankel asymptotic (j1q.c) |
| bessel\_jn(3,.) | 0.5 ulp | **4.4×** | forward/backward recurrence from j0/j1 |
| bessel\_y0 | 2.0 ulp | **7.0×** | piecewise rational + Hankel asymptotic (j0q.c) |
| bessel\_y1 | 1.5 ulp | **7.3×** | piecewise rational + Hankel asymptotic (j1q.c) |
| bessel\_yn(3,.) | 7.6 ulp | **6.9×** | forward recurrence from y0/y1 |
| bessel\_yn\_range(0..5) | 5.4 ulp | **29×** | single forward-recurrence sweep, 6 outputs / call |

### Complex arithmetic

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| cdd\_add | — | **2.4×** | original: component-wise DD add |
| cdd\_sub | — | **2.4×** | original: component-wise DD sub |
| cdd\_mul | — | **8.9×** | original: (ac−bd, ad+bc) via DD ops |
| cdd\_div | — | **5.4×** | original: (ac+bd, bc−ad)/(c²+d²) |
| cdd\_conjg | exact | 1.2× | original: negate im limbs |
| cdd\_proj | — | **2.7×** | C99 Annex G Riemann-sphere projection (identity for finite z) |
| cdd\_abs | — | **10×** | original: hypot(re, im) |
| cdd\_arg | — | **3.5×** | original: atan2(im, re) |

### Complex transcendentals

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| cdd\_sqrt | — | **7.8×** | original: Kahan-style (\|z\|+\|a\|)/2 with scaling |
| cdd\_exp | — | **3.9×** | original: exp(re)·(cos(im), sin(im)) |
| cdd\_expm1 | — | 1.7× | original: expm1(re)·cos(im) + (cos(im)−1) + i·exp(re)·sin(im) |
| cdd\_log | — | **4.8×** | original: (log(\|z\|), atan2(im,re)) |
| cdd\_log2 | — | **5.4×** | original: clog / log(2) component-wise |
| cdd\_log10 | — | **5.4×** | original: clog / log(10) component-wise |
| cdd\_log1p | — | **5.7×** | original: cancellation-safe log(1+z) near z=0 |
| cdd\_pow | — | **4.1×** | original: exp(w·log(z)) |
| cdd\_sin | — | **4.5×** | original: sin(re)cosh(im), cos(re)sinh(im) |
| cdd\_cos | — | **4.4×** | original: cos(re)cosh(im), −sin(re)sinh(im) |
| cdd\_tan | — | **4.5×** | original: complex sin/cos ratio |
| cdd\_sinpi | — | **4.1×** | original: csin(π·z) via π-scaled trig kernels |
| cdd\_cospi | — | **4.1×** | original: ccos(π·z) via π-scaled trig kernels |
| cdd\_sinh | — | **4.5×** | original: sinh(re)cos(im), cosh(re)sin(im) |
| cdd\_cosh | — | **4.6×** | original: cosh(re)cos(im), sinh(re)sin(im) |
| cdd\_tanh | — | **4.6×** | original: complex tanh via sinh/cosh |
| cdd\_asin | — | **5.0×** | original: −i·log(iz+√(1−z²)) |
| cdd\_acos | — | **4.1×** | original: π/2 − asin(z) |
| cdd\_atan | — | **4.9×** | original: (−i/2)·log((1+iz)/(1−iz)) |
| cdd\_asinh | — | **4.6×** | original: log(z+√(z²+1)) |
| cdd\_acosh | — | **4.3×** | original: log(z+√(z²−1)) |
| cdd\_atanh | — | **5.1×** | original: ½·log((1+z)/(1−z)) |

## Fortran: `float64x2` vs `real(16)`

Ops timed through the **Fortran elemental wrappers** on a `sequence`-
derived-type DD. gfortran passes/returns these via hidden pointers
(not FP registers), adding ~1.5× overhead vs the C ABI column even
after LTO. Scope: the full Fortran module surface — everything in the
C section above **plus** Fortran-only ops (array reductions, numeric
inquiry). **prec** column carries the precision label for the shared
surface; C-section duplicates omit it for table width.

### Arithmetic

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| add | full DD | 0.7 ulp | 2.0× | Julia: two\_sum EFT |
| sub | full DD | 0.7 ulp | **2.2×** | Julia: two\_sum EFT (negate + add) |
| mul | full DD | 2.0 ulp | **4.6×** | Julia: two\_prod EFT via FMA |
| div | full DD | 3.1 ulp | **2.8×** | original: Newton refinement (1/y seed, one step) |
| sqrt | full DD | 1.7 ulp | **42×** | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) |
| add (dd+dp) | exact | exact | 1.9× | Julia: two\_sum EFT |
| mul (dp\*dd) | full DD | 1.0 ulp | **5.0×** | Julia: two\_prod EFT via FMA |

### Unary

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| abs | exact | exact | 0.99× | original: sign-check + negate limbs |
| neg | exact | exact | 1.4× | original: negate both limbs |
| aint | exact | — | 1.0× | original: truncate hi, check DD fractional part |
| anint | exact | — | 1.5× | original: truncate hi, DD fractional part vs ±0.5 |
| fraction | exact | — | 1.2× | original: scale both limbs by −exponent |
| scale | exact | exact | 0.90× | original: ldexp on both limbs |
| set\_exponent | exact | exact | 1.5× | original: scale + set\_exponent on hi |

### Binary

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| min | full DD | — | 0.91× | original: DD comparison + select |
| max | full DD | — | 1.1× | original: DD comparison + select |
| min3 | full DD | exact | **2.7×** | original: chained min |
| max3 | full DD | exact | **2.9×** | original: chained max |
| sign | exact | — | 1.1× | original: sign-check + negate |
| dim | full DD | — | **2.8×** | original: DD comparison, then subtract or zero |
| hypot | full DD | 2.4 ulp | **4.2×** | original: scaled sqrt(x²+y²) |
| mod | full DD | — | 0.45× | sample: floor-multiple reduction loop; fallback to div chain |
| modulo | full DD | — | 0.94× | original: mod + sign adjustment |

### Exponential / logarithmic

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| exp | full DD | 1.4 ulp | **3.9×** | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction |
| log | full DD | 1.4 ulp | **6.5×** | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) |
| log10 | full DD | 1.2 ulp | **8.9×** | Julia: log2 kernel × DD log10(2) |
| pow | full DD | 11 ulp | **5.3×** | Julia: exp(y × log(x)) |
| pow\_int | full DD | 11 ulp | **7.6×** | original: repeated squaring via DD mul |

### Trigonometric

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| sin | full DD | 1.4 ulp | **3.5×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| cos | full DD | 1.1 ulp | **3.4×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| sincos | full DD | 1.4 ulp (sin) / 1.1 ulp (cos) | **5.4×** | fused sin/cos: one range-reduction, two Taylor outputs |
| sinpi | full DD | 1.5 ulp | **4.4×** | Julia: sinpi Horner polynomial, direct |
| cospi | full DD | 1.4 ulp | **4.5×** | Julia: cospi Horner polynomial, direct |
| tan | full DD | 2.3 ulp | **2.5×** | original: sin/cos Taylor kernels + DD divide |
| tanpi | full DD | 1.9 ulp | **3.3×** | original: sinpi/cospi ratio |
| asin | full DD | 0.7 ulp | **6.1×** | original: piecewise rational P/Q (3 regions, from libquadmath asinq.c) |
| asinpi | full DD | 0.5 ulp | **5.8×** | original: asin(x)/π with exact-DD π division |
| acos | full DD | 0.6 ulp | **6.1×** | original: asin polynomial + half-angle identity |
| acospi | full DD | 1.2 ulp | **6.1×** | original: acos(x)/π with exact-DD π division |
| atan | full DD | 0.7 ulp | **3.9×** | original: 84-entry table lookup + rational P(t²)/Q(t²) (from libquadmath atanq.c) |
| atanpi | full DD | 0.6 ulp | **4.5×** | original: atan(x)/π with exact-DD π division |
| atan2 | full DD | 0.8 ulp | **3.4×** | original: table-based atan + quadrant correction |
| atan2pi | full DD | 1.4 ulp | **3.6×** | original: atan2(y,x)/π with exact-DD π division |

### Hyperbolic

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| sinh | full DD | 4.3 ulp | **3.1×** | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 |
| cosh | full DD | 1.0 ulp | **2.2×** | original: (exp+exp⁻¹)/2 |
| sinhcosh | full DD | 4.3 ulp (sinh) / 1.0 ulp (cosh) | **4.7×** | fused sinh/cosh: one range-reduction, two outputs |
| tanh | full DD | 4.1 ulp | **3.5×** | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) |
| asinh | full DD | 2.0 ulp | **7.3×** | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton |
| acosh | full DD | 0.8 ulp | **6.8×** | original: log(x+√(x²−1)) with Newton correction |
| atanh | full DD | 2.0 ulp | **6.4×** | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) |

### Error / special functions

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| erf | full DD | 0.7 ulp | **3.8×** | piecewise rational approx (libquadmath erfq.c) |
| erfc | full DD | 0.5 ulp | **3.9×** | piecewise rational approx + split exp(-x^2) |
| erfc\_scaled | full DD | — | **5.9×** | exp(x^2)·erfc(x) with asymptotic cancellation |
| gamma | full DD | 18 ulp | **7.5×** | piecewise rational approx + Stirling + reflection |
| log\_gamma | full DD | 2.0 ulp | **5.2×** | piecewise rational approx + Stirling asymptotic |
| bessel\_j0 | full DD | 0.2 ulp | **6.8×** | piecewise rational + Hankel asymptotic (j0q.c) via C++ |
| bessel\_j1 | full DD | 0.4 ulp | **6.8×** | piecewise rational + Hankel asymptotic (j1q.c) via C++ |
| bessel\_jn(3,.) | full DD | 0.5 ulp | **4.2×** | forward/backward recurrence from j0/j1 |
| bessel\_y0 | full DD | 2.0 ulp | **7.0×** | piecewise rational + Hankel asymptotic (j0q.c) via C++ |
| bessel\_y1 | full DD | 1.5 ulp | **7.1×** | piecewise rational + Hankel asymptotic (j1q.c) via C++ |
| bessel\_yn(3,.) | full DD | 7.6 ulp | **6.9×** | forward recurrence from y0/y1 |

### Complex arithmetic

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| cdd\_add | full DD | 0.7 ulp (re) / 0.7 ulp (im) | **2.7×** | original: component-wise DD add |
| cdd\_sub | full DD | 0.7 ulp (re) / 0.6 ulp (im) | **2.7×** | original: component-wise DD sub |
| cdd\_mul | full DD | 1.7 ulp (re) / 1.2 ulp (im) | **2.1×** | original: (ac−bd, ad+bc) via DD ops |
| cdd\_div | full DD / deriv | 40 ulp (re) / 1.9 ulp (im) | **4.2×** | original: (ac+bd, bc−ad)/(c²+d²) |
| cdd\_conjg | exact | exact | **2.2×** | original: negate im limbs |
| cdd\_abs | full DD | 1.9 ulp | **5.1×** | original: hypot(re, im) |

### Complex transcendentals

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| cdd\_sqrt | full DD | 1.9 ulp (re) / 1.8 ulp (im) | **7.5×** | original: Kahan-style (\|z\|+\|a\|)/2 with scaling |
| cdd\_exp | full DD | 2.1 ulp (re) / 2.0 ulp (im) | **3.6×** | original: exp(re)·(cos(im), sin(im)) |
| cdd\_expm1 | reduced DD | 1.8 ulp (re) / 2.0 ulp (im) | 1.5× | original: expm1 + complex rotation with cancellation-safe Re |
| cdd\_log | full DD | 1.6 ulp (re) / 1.1 ulp (im) | **4.9×** | original: (log(\|z\|), atan2(im,re)) |
| cdd\_log2 | full DD | 1.7 ulp (re) / 1.2 ulp (im) | **5.7×** | original: clog / log(2) component-wise |
| cdd\_log10 | full DD | 1.5 ulp (re) / 1.4 ulp (im) | **5.7×** | original: clog / log(10) component-wise |
| cdd\_log1p | full DD | 20769 ulp (re) / 2.0 ulp (im) | **6.0×** | original: cancellation-safe log(1+z) near z=0 |
| cdd\_pow | full DD | 695 ulp (re) / 2498 ulp (im) | **4.0×** | original: exp(w·log(z)) |
| cdd\_sin | full DD | 2.8 ulp (re) / 4.2 ulp (im) | **4.4×** | original: sin(re)cosh(im), cos(re)sinh(im) |
| cdd\_cos | full DD | 1.8 ulp (re) / 3.0 ulp (im) | **4.2×** | original: cos(re)cosh(im), −sin(re)sinh(im) |
| cdd\_tan | full DD | 3.9 ulp (re) / 5.8 ulp (im) | **4.1×** | original: complex sin/cos ratio |
| cdd\_sinh | full DD | 3.3 ulp (re) / 2.2 ulp (im) | **4.4×** | original: sinh(re)cos(im), cosh(re)sin(im) |
| cdd\_cosh | full DD | 1.8 ulp (re) / 2.2 ulp (im) | **4.6×** | original: cosh(re)cos(im), sinh(re)sin(im) |
| cdd\_tanh | full DD | 5.2 ulp (re) / 4.7 ulp (im) | **4.7×** | original: complex tanh via sinh/cosh |
| cdd\_asin | deriv / full DD | 1.9 ulp (re) / 2.3 ulp (im) | **5.2×** | original: −i·log(iz+√(1−z²)) |
| cdd\_acos | full DD | 1.6 ulp (re) / 2.5 ulp (im) | **4.5×** | original: π/2 − asin(z) |
| cdd\_atan | full DD | 1.6 ulp (re) / 2.2 ulp (im) | **4.5×** | original: (i/2)·log((i+z)/(i−z)) |
| cdd\_asinh | deriv / full DD | 2.1 ulp (re) / 1.7 ulp (im) | **4.7×** | original: log(z+√(z²+1)) |
| cdd\_acosh | full DD | 2.5 ulp (re) / 1.6 ulp (im) | **4.7×** | original: log(z+√(z²−1)) |
| cdd\_atanh | deriv / full DD | 2.4 ulp (re) / 1.9 ulp (im) | **4.9×** | original: ½·log((1+z)/(1−z)) |

### Array reductions

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| arr\_sum (n=8) | full DD | — | **3.8×** | original: chained DD add |
| arr\_product (n=8) | full DD | — | **2.1×** | original: chained DD mul |
| arr\_maxval (n=8) | full DD | — | **2.8×** | original: chained DD compare |
| arr\_minval (n=8) | full DD | — | **2.8×** | original: chained DD compare |
| arr\_dot (n=8) | full DD | — | **4.1×** | original: fused multiply-accumulate with periodic renormalization |
| arr\_norm2 (n=8) | full DD | — | **6.1×** | original: sqrt(dot(x,x)) |
| arr\_matmul (8×8·8) | full DD | — | **2.5×** | original: AXPY-order C kernel, MR=8 register-blocked panels + 1..7 tail, periodic renorm |

## Notes

- **`mod` / `fmod`** is the only operation where quad precision is
  consistently faster. The DD `mod` uses a floor-multiple reduction loop
  for small quotients and a full DD divide chain for large quotients;
  libquadmath's `fmodq` uses a specialized bit-level remainder algorithm.

- **Trig range reduction** uses a 3-part π/2 constant (~161 bits) via
  Cody–Waite subtraction with DD arithmetic (FMA-captured product errors).
  Combined with the π/8 argument split, this gives full DD precision
  (~4e-32) for sin/cos/tan with the current 13-term Taylor kernels.
  For |x| > ~1e15, a Payne–Hanek reduction with a multi-word 2/π table
  would be needed.

- **tgamma / lgamma** (both C and Fortran) use a native double-double
  Stirling kernel shifting the argument up to x ≥ 25 via a product
  accumulator so a single `logdd(prod)` absorbs the recurrence. Small
  arguments (x < 0.5) use the reflection identity.

- **Inverse trigonometric** functions (asin, acos, atan, atan2) use
  piecewise rational polynomial approximations ported from libquadmath's
  `asinq.c` / `atanq.c`.

- **erf / erfc** use piecewise rational approximation coefficients from
  libquadmath's `erfq.c`, evaluated in full DD arithmetic via Estrin's
  scheme (`dd_neval` / `dd_deval`).

- The **Fortran** multifloats module uses `elemental` functions on a
  `sequence` derived type. gfortran's ABI passes/returns these via hidden
  pointers (not in FP registers), adding ~1.5× overhead vs the C ABI
  leg. See the performance note at the top of `fsrc/multifloats.fypp`
  for details and the `bind(c)` escape hatch.

- **Array reductions** (`dot_product`, `matmul`) use a fused multiply-
  accumulate kernel that computes the product's error-free representation
  and accumulates corrections into a scalar `s_lo`. `matmul` routes to a
  C kernel (`src/multifloats_math.cc`) in AXPY / gaxpy loop order; the
  kernel register-blocks any `m` via a strided panel template at `MR=8`
  plus a 1..7-row tail handler.

<!-- Auto-generated by bench/build_benchmark_md.py from per-system JSON
     results in bench/results/. Do not edit by hand. -->
