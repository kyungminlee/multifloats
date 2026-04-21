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
| abs | — |
| neg | — |
| fmin | exact |
| fmax | exact |
| fdim | — |
| copysign | — |
| fmod | — |
| hypot | 0.7 qp ulp |
| exp | 0.5 qp ulp |
| expm1 | 0.7 qp ulp |
| log | 0.5 qp ulp |
| log10 | 0.5 qp ulp |
| log1p | 0.7 qp ulp |
| pow | 0.5 qp ulp |
| sin | 0.8 qp ulp |
| cos | 0.7 qp ulp |
| sinpi | 33853776 qp ulp |
| cospi | 118903598 qp ulp |
| tan | 0.6 qp ulp |
| tanpi | 118903598 qp ulp |
| asin | 0.6 qp ulp |
| asinpi | 1.0 qp ulp |
| acos | 0.5 qp ulp |
| acospi | 0.9 qp ulp |
| atan | 0.8 qp ulp |
| atanpi | 1.1 qp ulp |
| atan2 | 0.9 qp ulp |
| sinh | 1.2 qp ulp |
| cosh | 0.8 qp ulp |
| tanh | 1.3 qp ulp |
| asinh | 1.9 qp ulp |
| acosh | 0.9 qp ulp |
| atanh | 1.1 qp ulp |
| erf | 0.7 qp ulp |
| erfc | 2.3 qp ulp |
| tgamma | 1.9 qp ulp |
| lgamma | 1.5 qp ulp |
| bessel\_j0 | 7.5 qp ulp |
| bessel\_j1 | 14 qp ulp |
| bessel\_jn(3,.) | 90 qp ulp |
| bessel\_y0 | 2.6 qp ulp |
| bessel\_y1 | 41 qp ulp |
| bessel\_yn(3,.) | 1071 qp ulp |
| cdd\_add | 0.5 qp ulp (re) / 0.5 qp ulp (im) |
| cdd\_sub | 0.5 qp ulp (re) / 0.5 qp ulp (im) |
| cdd\_mul | 536 qp ulp (re) / 34 qp ulp (im) |
| cdd\_div | 83 qp ulp (re) / 90 qp ulp (im) |
| cdd\_abs | 0.7 qp ulp |
| cdd\_sqrt | 1.3 qp ulp (re) / 1.1 qp ulp (im) |
| cdd\_exp | 1.5 qp ulp (re) / 1.3 qp ulp (im) |
| cdd\_log | 1.4 qp ulp (re) / 1.0 qp ulp (im) |
| cdd\_sin | 1.7 qp ulp (re) / 1.4 qp ulp (im) |
| cdd\_cos | 1.5 qp ulp (re) / 1.4 qp ulp (im) |
| cdd\_tan | 75392150385925696 qp ulp (re) / 2.5 qp ulp (im) |
| cdd\_sinh | 1.5 qp ulp (re) / 1.7 qp ulp (im) |
| cdd\_cosh | 1.4 qp ulp (re) / 1.5 qp ulp (im) |
| cdd\_tanh | 2.9 qp ulp (re) / 48179322550344664 qp ulp (im) |
| cdd\_asin | 2.6 qp ulp (re) / 1.9 qp ulp (im) |
| cdd\_acos | 1.8 qp ulp (re) / 1.9 qp ulp (im) |
| cdd\_atan | 1.6 qp ulp (re) / 2.4 qp ulp (im) |
| cdd\_asinh | 2.2 qp ulp (re) / 2.2 qp ulp (im) |
| cdd\_acosh | 1.9 qp ulp (re) / 1.8 qp ulp (im) |
| cdd\_atanh | 2.3 qp ulp (re) / 1.7 qp ulp (im) |
| mul (dp\*dd) | — |
| aint | — |
| anint | — |
| fraction | — |
| min | — |
| max | — |
| min3 | — |
| max3 | — |
| sign | — |
| dim | — |
| mod | — |
| modulo | — |
| pow\_int | — |
| erfc\_scaled | — |
| gamma | — |
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
| add | 0.7 ulp | **2.9×** | Julia: two\_sum EFT |
| sub | 0.7 ulp | **2.5×** | Julia: two\_sum EFT (negate + add) |
| mul | 2.0 ulp | **12×** | Julia: two\_prod EFT via FMA |
| div | 3.0 ulp | **7.0×** | original: Newton refinement (1/y seed, one step) |
| sqrt | 1.7 ulp | **31×** | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) |
| fma | — | **82×** | original: x\*y + z via DD ops |
| abs | — | **2.1×** | original: sign-check + negate limbs |
| neg | — | **2.7×** | original: negate both limbs |

### Rounding

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| trunc | exact | 1.6× | original: signbit ? −floor(−x) : floor(x) |
| round | exact | 1.9× | original: trunc(x + ½·sign(x)) |

### Binary

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| fmin | exact | **5.2×** | original: DD comparison + select |
| fmax | exact | **4.8×** | original: DD comparison + select |
| fdim | — | **5.2×** | original: DD comparison, then subtract or zero |
| copysign | — | 1.8× | original: sign-bit copy to hi, propagate to lo |
| fmod | — | 0.68× | sample: floor-multiple reduction loop; fallback to div chain |
| hypot | 2.4 ulp | **10×** | original: scaled sqrt(x²+y²) |

### Exponential / logarithmic

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| exp | 1.4 ulp | **3.8×** | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction |
| exp2 | — | **5.0×** | Julia: exp2 polynomial (14-term Horner) |
| expm1 | 1.5 ulp | **5.6×** | original: exp(x) − 1 via DD sub |
| log | 1.5 ulp | **6.7×** | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) |
| log10 | 1.2 ulp | **8.7×** | Julia: log2 kernel × DD log10(2) |
| log2 | — | **8.6×** | Julia: log2 table lookup + polynomial |
| log1p | 3.8 ulp | **6.6×** | original: log(1 + x) via DD add |
| pow | 16 ulp | **5.8×** | Julia: exp(y × log(x)) |

### Trigonometric

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| sin | 1.8 ulp | **3.6×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| cos | 1.4 ulp | **3.6×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| sinpi | 1.5 ulp | **4.5×** | Julia: sinpi Horner polynomial, direct |
| cospi | 1.5 ulp | **4.7×** | Julia: cospi Horner polynomial, direct |
| tan | 2.3 ulp | **2.6×** | original: sin/cos Taylor kernels + DD divide |
| tanpi | 2.1 ulp | **3.4×** | original: sinpi/cospi ratio |
| asin | 0.7 ulp | **6.0×** | original: piecewise rational P/Q (3 regions, from libquadmath asinq.c) |
| asinpi | 1.4 ulp | **6.2×** | original: asin(x)/π with exact-DD π division |
| acos | 0.6 ulp | **6.1×** | original: asin polynomial + half-angle identity |
| acospi | 1.4 ulp | **6.0×** | original: acos(x)/π with exact-DD π division |
| atan | 0.9 ulp | **4.0×** | original: 84-entry table lookup + rational P(t²)/Q(t²) (from libquadmath atanq.c) |
| atanpi | 1.8 ulp | **4.3×** | original: atan(x)/π with exact-DD π division |
| atan2 | 2.1 ulp | **3.6×** | original: table-based atan + quadrant correction |

### Hyperbolic

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| sinh | 4.3 ulp | **3.2×** | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 |
| cosh | 1.1 ulp | **2.2×** | original: (exp+exp⁻¹)/2 |
| tanh | 3.4 ulp | **3.6×** | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) |
| asinh | 2.0 ulp | **7.4×** | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton |
| acosh | 1.4 ulp | **7.0×** | original: log(x+√(x²−1)) with Newton correction |
| atanh | 2.4 ulp | **6.4×** | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) |

### Error / special functions

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| erf | 0.7 ulp | **4.1×** | piecewise rational approx (ported from libquadmath erfq.c) |
| erfc | 1.8 ulp | **4.3×** | piecewise rational approx + split exp(-x^2) |
| tgamma | 18 ulp | **8.7×** | piecewise rational approx + Stirling + reflection, exp(lgamma) |
| lgamma | 15 ulp | **5.0×** | piecewise rational approx + Stirling asymptotic |
| bessel\_j0 | 452 ulp | **6.9×** | piecewise rational + Hankel asymptotic (j0q.c) |
| bessel\_j1 | 1494 ulp | **7.0×** | piecewise rational + Hankel asymptotic (j1q.c) |
| bessel\_jn(3,.) | 598 ulp | **4.5×** | forward/backward recurrence from j0/j1 |
| bessel\_y0 | 162 ulp | **7.0×** | piecewise rational + Hankel asymptotic (j0q.c) |
| bessel\_y1 | 1884 ulp | **6.8×** | piecewise rational + Hankel asymptotic (j1q.c) |
| bessel\_yn(3,.) | 102142 ulp | **7.3×** | forward recurrence from y0/y1 |

### Complex arithmetic

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| cdd\_add | 0.7 ulp (re) / 0.7 ulp (im) | **2.4×** | original: component-wise DD add |
| cdd\_sub | 0.7 ulp (re) / 0.7 ulp (im) | **2.8×** | original: component-wise DD sub |
| cdd\_mul | 1231 ulp (re) / 102 ulp (im) | **12×** | original: (ac−bd, ad+bc) via DD ops |
| cdd\_div | 227 ulp (re) / 2.6 ulp (im) | **7.3×** | original: (ac+bd, bc−ad)/(c²+d²) |
| cdd\_conjg | exact | 1.2× | original: negate im limbs |
| cdd\_abs | 2.0 ulp | **10×** | original: hypot(re, im) |

### Complex transcendentals

| op | M1 Max err | M1 Max × | approach |
|---|---|---|---|
| cdd\_sqrt | 3.0 ulp (re) / 2.9 ulp (im) | **7.7×** | original: Kahan-style (\|z\|+\|a\|)/2 with scaling |
| cdd\_exp | 2.2 ulp (re) / 2.2 ulp (im) | **4.0×** | original: exp(re)·(cos(im), sin(im)) |
| cdd\_log | 1.8 ulp (re) / 1.6 ulp (im) | **4.9×** | original: (log(\|z\|), atan2(im,re)) |
| cdd\_sin | 2.8 ulp (re) / 8.0 ulp (im) | **4.6×** | original: sin(re)cosh(im), cos(re)sinh(im) |
| cdd\_cos | 2.2 ulp (re) / 7.3 ulp (im) | **4.5×** | original: cos(re)cosh(im), −sin(re)sinh(im) |
| cdd\_tan | 589001174890044 ulp (re) / 7.9 ulp (im) | **4.6×** | original: complex sin/cos ratio |
| cdd\_sinh | 5.5 ulp (re) / 2.8 ulp (im) | **4.5×** | original: sinh(re)cos(im), cosh(re)sin(im) |
| cdd\_cosh | 2.4 ulp (re) / 5.7 ulp (im) | **4.5×** | original: cosh(re)cos(im), sinh(re)sin(im) |
| cdd\_tanh | 5.6 ulp (re) / 376400957424568 ulp (im) | **4.5×** | original: complex tanh via sinh/cosh |
| cdd\_asin | 5290 ulp (re) / 2680523 ulp (im) | **4.9×** | original: −i·log(iz+√(1−z²)) |
| cdd\_acos | 33863511 ulp (re) / 2680523 ulp (im) | **5.0×** | original: π/2 − asin(z) |
| cdd\_atan | 6.1 ulp (re) / 237020239 ulp (im) | **6.1×** | original: (−i/2)·log((1+iz)/(1−iz)) |
| cdd\_asinh | 2.9 ulp (re) / 2.9 ulp (im) | **4.6×** | original: log(z+√(z²+1)) |
| cdd\_acosh | 2810331 ulp (re) / 2.6 ulp (im) | **4.2×** | original: log(z+√(z²−1)) |
| cdd\_atanh | 3.3 ulp (re) / 1.9 ulp (im) | **5.0×** | original: ½·log((1+z)/(1−z)) |

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
| sub | full DD | 0.7 ulp | **2.0×** | Julia: two\_sum EFT (negate + add) |
| mul | full DD | 2.0 ulp | **4.3×** | Julia: two\_prod EFT via FMA |
| div | full DD | 3.0 ulp | **2.8×** | original: Newton refinement (1/y seed, one step) |
| sqrt | full DD | 1.7 ulp | **30×** | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) |
| add (dd+dp) | exact | exact | 2.0× | Julia: two\_sum EFT |
| mul (dp\*dd) | full DD | — | **5.0×** | Julia: two\_prod EFT via FMA |

### Unary

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| abs | exact | — | 1.1× | original: sign-check + negate limbs |
| neg | exact | — | 1.4× | original: negate both limbs |
| aint | exact | — | 1.1× | original: truncate hi, check DD fractional part |
| anint | exact | — | 1.5× | original: truncate hi, DD fractional part vs ±0.5 |
| fraction | exact | — | 1.2× | original: scale both limbs by −exponent |
| scale | exact | exact | 0.89× | original: ldexp on both limbs |
| set\_exponent | exact | exact | 1.5× | original: scale + set\_exponent on hi |

### Binary

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| min | full DD | — | 1.8× | original: DD comparison + select |
| max | full DD | — | 0.85× | original: DD comparison + select |
| min3 | full DD | — | **2.5×** | original: chained min |
| max3 | full DD | — | **2.3×** | original: chained max |
| sign | exact | — | 1.1× | original: sign-check + negate |
| dim | full DD | — | 1.7× | original: DD comparison, then subtract or zero |
| hypot | full DD | 2.4 ulp | **9.7×** | original: scaled sqrt(x²+y²) |
| mod | full DD | — | 0.38× | sample: floor-multiple reduction loop; fallback to div chain |
| modulo | full DD | — | 1.1× | original: mod + sign adjustment |

### Exponential / logarithmic

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| exp | full DD | 1.4 ulp | **3.9×** | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction |
| log | full DD | 1.5 ulp | **6.8×** | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) |
| log10 | full DD | 1.2 ulp | **8.5×** | Julia: log2 kernel × DD log10(2) |
| pow | full DD | 16 ulp | **5.3×** | Julia: exp(y × log(x)) |
| pow\_int | full DD | — | **8.8×** | original: repeated squaring via DD mul |

### Trigonometric

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| sin | full DD | 1.8 ulp | **3.6×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| cos | full DD | 1.4 ulp | **3.5×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| sinpi | full DD | 1.5 ulp | **4.4×** | Julia: sinpi Horner polynomial, direct |
| cospi | full DD | 1.5 ulp | **4.6×** | Julia: cospi Horner polynomial, direct |
| tan | full DD | 2.3 ulp | **2.6×** | original: sin/cos Taylor kernels + DD divide |
| tanpi | full DD | 2.1 ulp | **3.5×** | original: sinpi/cospi ratio |
| asin | full DD | 0.7 ulp | **5.8×** | original: piecewise rational P/Q (3 regions, from libquadmath asinq.c) |
| asinpi | full DD | 1.4 ulp | **5.6×** | original: asin(x)/π with exact-DD π division |
| acos | full DD | 0.6 ulp | **6.1×** | original: asin polynomial + half-angle identity |
| acospi | full DD | 1.4 ulp | **6.0×** | original: acos(x)/π with exact-DD π division |
| atan | full DD | 0.9 ulp | **4.0×** | original: 84-entry table lookup + rational P(t²)/Q(t²) (from libquadmath atanq.c) |
| atanpi | full DD | 1.8 ulp | **4.2×** | original: atan(x)/π with exact-DD π division |
| atan2 | full DD | 2.1 ulp | **3.5×** | original: table-based atan + quadrant correction |

### Hyperbolic

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| sinh | full DD | 4.3 ulp | **3.2×** | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 |
| cosh | full DD | 1.1 ulp | **2.2×** | original: (exp+exp⁻¹)/2 |
| tanh | full DD | 3.4 ulp | **3.6×** | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) |
| asinh | full DD | 2.0 ulp | **7.5×** | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton |
| acosh | full DD | 1.4 ulp | **6.5×** | original: log(x+√(x²−1)) with Newton correction |
| atanh | full DD | 2.4 ulp | **5.7×** | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) |

### Error / special functions

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| erf | full DD | 0.7 ulp | **4.4×** | piecewise rational approx (libquadmath erfq.c) |
| erfc | full DD | 1.8 ulp | **4.2×** | piecewise rational approx + split exp(-x^2) |
| erfc\_scaled | full DD | — | **5.6×** | exp(x^2)·erfc(x) with asymptotic cancellation |
| gamma | full DD | — | **7.8×** | piecewise rational approx + Stirling + reflection |
| log\_gamma | full DD | 15 ulp | **4.8×** | piecewise rational approx + Stirling asymptotic |
| bessel\_j0 | full DD | 452 ulp | **6.8×** | piecewise rational + Hankel asymptotic (j0q.c) via C++ |
| bessel\_j1 | full DD | 1494 ulp | **6.6×** | piecewise rational + Hankel asymptotic (j1q.c) via C++ |
| bessel\_jn(3,.) | full DD | 598 ulp | **4.5×** | forward/backward recurrence from j0/j1 |
| bessel\_y0 | full DD | 162 ulp | **7.0×** | piecewise rational + Hankel asymptotic (j0q.c) via C++ |
| bessel\_y1 | full DD | 1884 ulp | **6.9×** | piecewise rational + Hankel asymptotic (j1q.c) via C++ |
| bessel\_yn(3,.) | full DD | 102142 ulp | **7.2×** | forward recurrence from y0/y1 |

### Complex arithmetic

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| cdd\_add | full DD | 0.7 ulp (re) / 0.7 ulp (im) | 1.1× | original: component-wise DD add |
| cdd\_sub | full DD | 0.7 ulp (re) / 0.7 ulp (im) | **2.5×** | original: component-wise DD sub |
| cdd\_mul | full DD | 1231 ulp (re) / 102 ulp (im) | **6.6×** | original: (ac−bd, ad+bc) via DD ops |
| cdd\_div | full DD / deriv | 227 ulp (re) / 2.6 ulp (im) | **5.4×** | original: (ac+bd, bc−ad)/(c²+d²) |
| cdd\_conjg | exact | exact | **2.0×** | original: negate im limbs |
| cdd\_abs | full DD | 2.0 ulp | **9.3×** | original: hypot(re, im) |

### Complex transcendentals

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| cdd\_sqrt | full DD | 3.0 ulp (re) / 2.9 ulp (im) | **7.6×** | original: Kahan-style (\|z\|+\|a\|)/2 with scaling |
| cdd\_exp | full DD | 2.2 ulp (re) / 2.2 ulp (im) | **3.7×** | original: exp(re)·(cos(im), sin(im)) |
| cdd\_log | full DD | 1.8 ulp (re) / 1.6 ulp (im) | **5.5×** | original: (log(\|z\|), atan2(im,re)) |
| cdd\_sin | full DD | 2.8 ulp (re) / 8.0 ulp (im) | **4.3×** | original: sin(re)cosh(im), cos(re)sinh(im) |
| cdd\_cos | full DD | 2.2 ulp (re) / 7.3 ulp (im) | **4.5×** | original: cos(re)cosh(im), −sin(re)sinh(im) |
| cdd\_tan | full DD | 589001174890044 ulp (re) / 7.9 ulp (im) | **4.4×** | original: complex sin/cos ratio |
| cdd\_sinh | full DD | 5.5 ulp (re) / 2.8 ulp (im) | **4.8×** | original: sinh(re)cos(im), cosh(re)sin(im) |
| cdd\_cosh | full DD | 2.4 ulp (re) / 5.7 ulp (im) | **4.7×** | original: cosh(re)cos(im), sinh(re)sin(im) |
| cdd\_tanh | full DD | 5.6 ulp (re) / 376400957424568 ulp (im) | **4.6×** | original: complex tanh via sinh/cosh |
| cdd\_asin | deriv / full DD | 5290 ulp (re) / 2680523 ulp (im) | **5.1×** | original: −i·log(iz+√(1−z²)) |
| cdd\_acos | full DD | 33863511 ulp (re) / 2680523 ulp (im) | **4.9×** | original: π/2 − asin(z) |
| cdd\_atan | full DD | 6.1 ulp (re) / 237020239 ulp (im) | **5.4×** | original: (i/2)·log((i+z)/(i−z)) |
| cdd\_asinh | deriv / full DD | 2.9 ulp (re) / 2.9 ulp (im) | **4.7×** | original: log(z+√(z²+1)) |
| cdd\_acosh | full DD | 2810331 ulp (re) / 2.6 ulp (im) | **4.4×** | original: log(z+√(z²−1)) |
| cdd\_atanh | deriv / full DD | 3.3 ulp (re) / 1.9 ulp (im) | **5.2×** | original: ½·log((1+z)/(1−z)) |

### Array reductions

| op | prec | M1 Max err | M1 Max × | approach |
|---|---|---|---|---|
| arr\_sum (n=8) | full DD | — | **4.3×** | original: chained DD add |
| arr\_product (n=8) | full DD | — | **2.0×** | original: chained DD mul |
| arr\_maxval (n=8) | full DD | — | **3.2×** | original: chained DD compare |
| arr\_minval (n=8) | full DD | — | **3.3×** | original: chained DD compare |
| arr\_dot (n=8) | full DD | — | **4.2×** | original: fused multiply-accumulate with periodic renormalization |
| arr\_norm2 (n=8) | full DD | — | **5.7×** | original: sqrt(dot(x,x)) |
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
