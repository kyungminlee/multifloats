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
| **pop-os** | 13th Gen Intel(R) Core(TM) i3-1315U (x86_64), 8 cores, 15 GB | Pop!_OS 24.04 LTS (Linux 6.17.9-76061709-generic) | GNU Fortran / g++ 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) | CMake 3.28.3, `-O3 -flto`, STATIC library |
| **sandbox** | Intel (family 6, model 207, AVX-512 + AMX), 16 vCPU, 21 GB (hypervisor-masked) | Ubuntu 24.04.4 LTS (Linux 4.4.0) | GNU Fortran / g++ 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) | CMake 3.28.3, `-O3 -flto`, STATIC library |

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

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| add | 0.7 ulp | 192 ulp | 192 ulp | **2.9×** | **8.9×** | **4.5×** | Julia: two\_sum EFT |
| sub | 0.7 ulp | 249 ulp | 249 ulp | **2.5×** | **4.6×** | **4.6×** | Julia: two\_sum EFT (negate + add) |
| mul | 2.0 ulp | 2.2 ulp | 2.2 ulp | **12×** | **7.3×** | **9.5×** | Julia: two\_prod EFT via FMA |
| div | 3.0 ulp | 3.1 ulp | 3.1 ulp | **7.0×** | **4.9×** | **5.8×** | original: Newton refinement (1/y seed, one step) |
| sqrt | 1.7 ulp | 1.8 ulp | 1.8 ulp | **31×** | **53×** | **55×** | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) |
| fma | — | — | — | **82×** | **168×** | **161×** | original: x\*y + z via DD ops |
| abs | — | 0.2 ulp | 0.2 ulp | **2.1×** | **7.4×** | **9.5×** | original: sign-check + negate limbs |
| neg | — | 0.2 ulp | 0.2 ulp | **2.7×** | **6.4×** | **8.8×** | original: negate both limbs |

### Rounding

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| trunc | exact | exact | exact | 1.6× | **4.7×** | **4.9×** | original: signbit ? −floor(−x) : floor(x) |
| round | exact | exact | exact | 1.9× | **2.3×** | **2.2×** | original: trunc(x + ½·sign(x)) |

### Binary

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| fmin | exact | 0.2 ulp | 0.2 ulp | **5.2×** | **8.7×** | **9.7×** | original: DD comparison + select |
| fmax | exact | 0.2 ulp | 0.2 ulp | **4.8×** | **10×** | **12×** | original: DD comparison + select |
| fdim | — | 159 ulp | 159 ulp | **5.2×** | **9.2×** | **8.6×** | original: DD comparison, then subtract or zero |
| copysign | — | 0.2 ulp | 0.2 ulp | 1.8× | **5.3×** | **5.4×** | original: sign-bit copy to hi, propagate to lo |
| fmod | — | 5115223702040951 ulp | 2485812121 ulp | 0.68× | 1.1× | 1.1× | sample: floor-multiple reduction loop; fallback to div chain |
| hypot | 2.4 ulp | 2.1 ulp | 2.1 ulp | **10×** | **13×** | **13×** | original: scaled sqrt(x²+y²) |

### Exponential / logarithmic

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| exp | 1.4 ulp | 1.3 ulp | 1.1 ulp | **3.8×** | **5.0×** | **5.1×** | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction |
| exp2 | — | — | — | **5.0×** | **6.0×** | **6.4×** | Julia: exp2 polynomial (14-term Horner) |
| expm1 | 1.5 ulp | 2.0 ulp | 2.0 ulp | **5.6×** | **5.0×** | **4.5×** | original: exp(x) − 1 via DD sub |
| log | 1.5 ulp | 1.2 ulp | 1.2 ulp | **6.7×** | **4.6×** | **4.8×** | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) |
| log10 | 1.2 ulp | 1.2 ulp | 1.2 ulp | **8.7×** | **6.1×** | **6.3×** | Julia: log2 kernel × DD log10(2) |
| log2 | — | — | — | **8.6×** | **5.9×** | **6.0×** | Julia: log2 table lookup + polynomial |
| log1p | 3.8 ulp | 2.4 ulp | 2.4 ulp | **6.6×** | **6.6×** | **6.7×** | original: log(1 + x) via DD add |
| pow | 16 ulp | 26 ulp | 26 ulp | **5.8×** | **6.2×** | **6.2×** | Julia: exp(y × log(x)) |

### Trigonometric

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| sin | 1.8 ulp | 1.7 ulp | 1.7 ulp | **3.6×** | **3.1×** | **3.5×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| cos | 1.4 ulp | 2.0 ulp | 2.0 ulp | **3.6×** | **3.2×** | **3.5×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| sinpi | 1.5 ulp | — | — | **4.5×** | — | — | Julia: sinpi Horner polynomial, direct |
| cospi | 1.5 ulp | — | — | **4.7×** | — | — | Julia: cospi Horner polynomial, direct |
| tan | 2.3 ulp | 2.3 ulp | 2.3 ulp | **2.6×** | **2.2×** | **2.5×** | original: sin/cos Taylor kernels + DD divide |
| tanpi | 2.1 ulp | — | — | **3.4×** | — | — | original: sinpi/cospi ratio |
| asin | 0.7 ulp | 0.8 ulp | 0.8 ulp | **6.0×** | **5.8×** | **6.3×** | original: piecewise rational P/Q (3 regions, from libquadmath asinq.c) |
| asinpi | 1.4 ulp | — | — | **6.2×** | — | — | original: asin(x)/π with exact-DD π division |
| acos | 0.6 ulp | 0.6 ulp | 0.6 ulp | **6.1×** | **6.3×** | **6.6×** | original: asin polynomial + half-angle identity |
| acospi | 1.4 ulp | — | — | **6.0×** | — | — | original: acos(x)/π with exact-DD π division |
| atan | 0.9 ulp | 0.8 ulp | 0.8 ulp | **4.0×** | **4.4×** | **4.4×** | original: 84-entry table lookup + rational P(t²)/Q(t²) (from libquadmath atanq.c) |
| atanpi | 1.8 ulp | — | — | **4.3×** | — | — | original: atan(x)/π with exact-DD π division |
| atan2 | 2.1 ulp | 1.1 ulp | 1.1 ulp | **3.6×** | **4.1×** | **4.2×** | original: table-based atan + quadrant correction |

### Hyperbolic

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| sinh | 4.3 ulp | 3.0 ulp | 3.0 ulp | **3.2×** | **2.7×** | **2.7×** | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 |
| cosh | 1.1 ulp | 1.3 ulp | 1.1 ulp | **2.2×** | **2.4×** | **2.6×** | original: (exp+exp⁻¹)/2 |
| tanh | 3.4 ulp | 2.7 ulp | 2.7 ulp | **3.6×** | **3.1×** | **3.2×** | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) |
| asinh | 2.0 ulp | 2.7 ulp | 2.7 ulp | **7.4×** | **8.8×** | **8.9×** | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton |
| acosh | 1.4 ulp | 1.3 ulp | 1.3 ulp | **7.0×** | **8.4×** | **8.8×** | original: log(x+√(x²−1)) with Newton correction |
| atanh | 2.4 ulp | 3.1 ulp | 3.1 ulp | **6.4×** | **6.6×** | **6.9×** | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) |

### Error / special functions

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| erf | 0.7 ulp | 0.8 ulp | 0.8 ulp | **4.1×** | **4.7×** | **5.0×** | piecewise rational approx (ported from libquadmath erfq.c) |
| erfc | 1.8 ulp | 1.6 ulp | 1.6 ulp | **4.3×** | **4.6×** | **5.1×** | piecewise rational approx + split exp(-x^2) |
| tgamma | 18 ulp | 9.9 ulp | 9.9 ulp | **8.7×** | **12×** | **13×** | piecewise rational approx + Stirling + reflection, exp(lgamma) |
| lgamma | 15 ulp | 7.3 ulp | 7.3 ulp | **5.0×** | **4.6×** | **5.2×** | piecewise rational approx + Stirling asymptotic |
| bessel\_j0 | 452 ulp | — | — | **6.9×** | — | — | piecewise rational + Hankel asymptotic (j0q.c) |
| bessel\_j1 | 1494 ulp | — | — | **7.0×** | — | — | piecewise rational + Hankel asymptotic (j1q.c) |
| bessel\_jn(3,.) | 598 ulp | — | — | **4.5×** | — | — | forward/backward recurrence from j0/j1 |
| bessel\_y0 | 162 ulp | — | — | **7.0×** | — | — | piecewise rational + Hankel asymptotic (j0q.c) |
| bessel\_y1 | 1884 ulp | — | — | **6.8×** | — | — | piecewise rational + Hankel asymptotic (j1q.c) |
| bessel\_yn(3,.) | 102142 ulp | — | — | **7.3×** | — | — | forward recurrence from y0/y1 |

### Complex arithmetic

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| cdd\_add | 0.7 ulp (re) / 0.7 ulp (im) | — | — | **2.4×** | — | — | original: component-wise DD add |
| cdd\_sub | 0.7 ulp (re) / 0.7 ulp (im) | — | — | **2.8×** | — | — | original: component-wise DD sub |
| cdd\_mul | 1231 ulp (re) / 102 ulp (im) | — | — | **12×** | — | — | original: (ac−bd, ad+bc) via DD ops |
| cdd\_div | 227 ulp (re) / 2.6 ulp (im) | — | — | **7.3×** | — | — | original: (ac+bd, bc−ad)/(c²+d²) |
| cdd\_conjg | exact | exact | exact | 1.2× | — | — | original: negate im limbs |
| cdd\_abs | 2.0 ulp | — | — | **10×** | — | — | original: hypot(re, im) |

### Complex transcendentals

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| cdd\_sqrt | 3.0 ulp (re) / 2.9 ulp (im) | — | — | **7.7×** | — | — | original: Kahan-style (\|z\|+\|a\|)/2 with scaling |
| cdd\_exp | 2.2 ulp (re) / 2.2 ulp (im) | — | — | **4.0×** | — | — | original: exp(re)·(cos(im), sin(im)) |
| cdd\_log | 1.8 ulp (re) / 1.6 ulp (im) | — | — | **4.9×** | — | — | original: (log(\|z\|), atan2(im,re)) |
| cdd\_sin | 2.8 ulp (re) / 8.0 ulp (im) | — | — | **4.6×** | — | — | original: sin(re)cosh(im), cos(re)sinh(im) |
| cdd\_cos | 2.2 ulp (re) / 7.3 ulp (im) | — | — | **4.5×** | — | — | original: cos(re)cosh(im), −sin(re)sinh(im) |
| cdd\_tan | 589001174890044 ulp (re) / 7.9 ulp (im) | — | — | **4.6×** | — | — | original: complex sin/cos ratio |
| cdd\_sinh | 5.5 ulp (re) / 2.8 ulp (im) | — | — | **4.5×** | — | — | original: sinh(re)cos(im), cosh(re)sin(im) |
| cdd\_cosh | 2.4 ulp (re) / 5.7 ulp (im) | — | — | **4.5×** | — | — | original: cosh(re)cos(im), sinh(re)sin(im) |
| cdd\_tanh | 5.6 ulp (re) / 376400957424568 ulp (im) | — | — | **4.5×** | — | — | original: complex tanh via sinh/cosh |
| cdd\_asin | 5290 ulp (re) / 2680523 ulp (im) | — | — | **4.9×** | — | — | original: −i·log(iz+√(1−z²)) |
| cdd\_acos | 33863511 ulp (re) / 2680523 ulp (im) | — | — | **5.0×** | — | — | original: π/2 − asin(z) |
| cdd\_atan | 6.1 ulp (re) / 237020239 ulp (im) | — | — | **6.1×** | — | — | original: (−i/2)·log((1+iz)/(1−iz)) |
| cdd\_asinh | 2.9 ulp (re) / 2.9 ulp (im) | — | — | **4.6×** | — | — | original: log(z+√(z²+1)) |
| cdd\_acosh | 2810331 ulp (re) / 2.6 ulp (im) | — | — | **4.2×** | — | — | original: log(z+√(z²−1)) |
| cdd\_atanh | 3.3 ulp (re) / 1.9 ulp (im) | — | — | **5.0×** | — | — | original: ½·log((1+z)/(1−z)) |

## Fortran: `float64x2` vs `real(16)`

Ops timed through the **Fortran elemental wrappers** on a `sequence`-
derived-type DD. gfortran passes/returns these via hidden pointers
(not FP registers), adding ~1.5× overhead vs the C ABI column even
after LTO. Scope: the full Fortran module surface — everything in the
C section above **plus** Fortran-only ops (array reductions, numeric
inquiry). **prec** column carries the precision label for the shared
surface; C-section duplicates omit it for table width.

### Arithmetic

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| add | full DD | 0.7 ulp | 0.6 ulp | 0.6 ulp | 2.0× | **5.6×** | **3.9×** | Julia: two\_sum EFT |
| sub | full DD | 0.7 ulp | 0.2 ulp | 0.2 ulp | **2.0×** | **3.4×** | **3.5×** | Julia: two\_sum EFT (negate + add) |
| mul | full DD | 2.0 ulp | 1.4 ulp | 1.4 ulp | **4.3×** | **3.9×** | **7.7×** | Julia: two\_prod EFT via FMA |
| div | full DD | 3.0 ulp | 2.2 ulp | 2.2 ulp | **2.8×** | 1.8× | 2.0× | original: Newton refinement (1/y seed, one step) |
| sqrt | full DD | 1.7 ulp | 0.7 ulp | 0.7 ulp | **30×** | **37×** | **37×** | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) |
| add (dd+dp) | exact | exact | exact | exact | 2.0× | **4.0×** | **4.2×** | Julia: two\_sum EFT |
| mul (dp\*dd) | full DD | — | 1.4 ulp | 1.4 ulp | **5.0×** | **4.3×** | **5.2×** | Julia: two\_prod EFT via FMA |

### Unary

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| abs | exact | — | exact | exact | 1.1× | 1.7× | 1.2× | original: sign-check + negate limbs |
| neg | exact | — | exact | exact | 1.4× | **2.3×** | **2.3×** | original: negate both limbs |
| aint | exact | — | exact | exact | 1.1× | 1.2× | 0.55× | original: truncate hi, check DD fractional part |
| anint | exact | — | exact | exact | 1.5× | 1.2× | 1.2× | original: truncate hi, DD fractional part vs ±0.5 |
| fraction | exact | — | exact | exact | 1.2× | 1.5× | 1.6× | original: scale both limbs by −exponent |
| scale | exact | exact | exact | exact | 0.89× | **5.2×** | **3.9×** | original: ldexp on both limbs |
| set\_exponent | exact | exact | exact | exact | 1.5× | **3.7×** | **3.9×** | original: scale + set\_exponent on hi |

### Binary

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| min | full DD | — | 0.2 ulp | 0.2 ulp | 1.8× | **3.3×** | **3.1×** | original: DD comparison + select |
| max | full DD | — | 0.2 ulp | 0.2 ulp | 0.85× | **3.1×** | **3.1×** | original: DD comparison + select |
| min3 | full DD | — | 0.2 ulp | 0.2 ulp | **2.5×** | **3.6×** | **3.6×** | original: chained min |
| max3 | full DD | — | 0.2 ulp | 0.2 ulp | **2.3×** | **3.7×** | **3.5×** | original: chained max |
| sign | exact | — | exact | exact | 1.1× | 1.6× | 1.5× | original: sign-check + negate |
| dim | full DD | — | 0.2 ulp | 0.2 ulp | 1.7× | **4.0×** | **4.3×** | original: DD comparison, then subtract or zero |
| hypot | full DD | 2.4 ulp | 2.5 ulp | 2.5 ulp | **9.7×** | **13×** | **14×** | original: scaled sqrt(x²+y²) |
| mod | full DD | — | 0.8 ulp | 0.8 ulp | 0.38× | 0.91× | 0.76× | sample: floor-multiple reduction loop; fallback to div chain |
| modulo | full DD | — | 0.8 ulp | 0.8 ulp | 1.1× | 1.3× | 1.1× | original: mod + sign adjustment |

### Exponential / logarithmic

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| exp | full DD | 1.4 ulp | 1.5 ulp | 1.1 ulp | **3.9×** | **5.2×** | **5.4×** | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction |
| log | full DD | 1.5 ulp | 1.8 ulp | 1.8 ulp | **6.8×** | **4.9×** | **5.1×** | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) |
| log10 | full DD | 1.2 ulp | 1.2 ulp | 1.2 ulp | **8.5×** | **6.5×** | **6.4×** | Julia: log2 kernel × DD log10(2) |
| pow | full DD | 16 ulp | 77 ulp | 77 ulp | **5.3×** | **6.1×** | **6.2×** | Julia: exp(y × log(x)) |
| pow\_int | full DD | — | 1.0 ulp | 1.0 ulp | **8.8×** | **4.1×** | **5.8×** | original: repeated squaring via DD mul |

### Trigonometric

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| sin | full DD | 1.8 ulp | 2.0 ulp | 2.0 ulp | **3.6×** | **3.2×** | **3.5×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| cos | full DD | 1.4 ulp | 1.9 ulp | 1.9 ulp | **3.5×** | **3.2×** | **3.6×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| sinpi | full DD | 1.5 ulp | — | — | **4.4×** | **3.9×** | **4.2×** | Julia: sinpi Horner polynomial, direct |
| cospi | full DD | 1.5 ulp | — | — | **4.6×** | **4.0×** | **4.4×** | Julia: cospi Horner polynomial, direct |
| tan | full DD | 2.3 ulp | 2.4 ulp | 2.4 ulp | **2.6×** | **2.3×** | **2.6×** | original: sin/cos Taylor kernels + DD divide |
| tanpi | full DD | 2.1 ulp | — | — | **3.5×** | — | — | original: sinpi/cospi ratio |
| asin | full DD | 0.7 ulp | 0.3 ulp | 0.3 ulp | **5.8×** | **5.7×** | **6.1×** | original: piecewise rational P/Q (3 regions, from libquadmath asinq.c) |
| asinpi | full DD | 1.4 ulp | — | — | **5.6×** | — | — | original: asin(x)/π with exact-DD π division |
| acos | full DD | 0.6 ulp | 0.5 ulp | 0.5 ulp | **6.1×** | **6.1×** | **6.4×** | original: asin polynomial + half-angle identity |
| acospi | full DD | 1.4 ulp | — | — | **6.0×** | — | — | original: acos(x)/π with exact-DD π division |
| atan | full DD | 0.9 ulp | 1.1 ulp | 1.1 ulp | **4.0×** | **4.0×** | **3.5×** | original: 84-entry table lookup + rational P(t²)/Q(t²) (from libquadmath atanq.c) |
| atanpi | full DD | 1.8 ulp | — | — | **4.2×** | — | — | original: atan(x)/π with exact-DD π division |
| atan2 | full DD | 2.1 ulp | 0.7 ulp | 0.7 ulp | **3.5×** | **3.8×** | **3.9×** | original: table-based atan + quadrant correction |

### Hyperbolic

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| sinh | full DD | 4.3 ulp | 4.5 ulp | 4.5 ulp | **3.2×** | **3.0×** | **3.2×** | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 |
| cosh | full DD | 1.1 ulp | 1.5 ulp | 1.1 ulp | **2.2×** | **2.6×** | **2.7×** | original: (exp+exp⁻¹)/2 |
| tanh | full DD | 3.4 ulp | 4.6 ulp | 4.6 ulp | **3.6×** | **3.1×** | **3.5×** | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) |
| asinh | full DD | 2.0 ulp | 2.7 ulp | 2.7 ulp | **7.5×** | **10×** | **10×** | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton |
| acosh | full DD | 1.4 ulp | 1.4 ulp | 1.4 ulp | **6.5×** | **6.9×** | **7.0×** | original: log(x+√(x²−1)) with Newton correction |
| atanh | full DD | 2.4 ulp | 2.5 ulp | 2.5 ulp | **5.7×** | **7.8×** | **7.5×** | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) |

### Error / special functions

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| erf | full DD | 0.7 ulp | 0.6 ulp | 0.6 ulp | **4.4×** | **4.4×** | **5.0×** | piecewise rational approx (libquadmath erfq.c) |
| erfc | full DD | 1.8 ulp | 1.1 ulp | 1.1 ulp | **4.2×** | **4.7×** | **5.0×** | piecewise rational approx + split exp(-x^2) |
| erfc\_scaled | full DD | — | 1.5 ulp | 1.5 ulp | **5.6×** | **7.0×** | **7.3×** | exp(x^2)·erfc(x) with asymptotic cancellation |
| gamma | full DD | — | 10 ulp | 11 ulp | **7.8×** | **7.7×** | **7.3×** | piecewise rational approx + Stirling + reflection |
| log\_gamma | full DD | 15 ulp | 2.1 ulp | 2.1 ulp | **4.8×** | **4.6×** | **5.1×** | piecewise rational approx + Stirling asymptotic |
| bessel\_j0 | full DD | 452 ulp | 851 ulp | 567 ulp | **6.8×** | **6.1×** | **6.9×** | piecewise rational + Hankel asymptotic (j0q.c) via C++ |
| bessel\_j1 | full DD | 1494 ulp | 1725 ulp | 945 ulp | **6.6×** | **5.8×** | **6.6×** | piecewise rational + Hankel asymptotic (j1q.c) via C++ |
| bessel\_jn(3,.) | full DD | 598 ulp | 2024 ulp | 2024 ulp | **4.5×** | **4.5×** | **4.9×** | forward/backward recurrence from j0/j1 |
| bessel\_y0 | full DD | 162 ulp | 584 ulp | 145 ulp | **7.0×** | **6.3×** | **7.2×** | piecewise rational + Hankel asymptotic (j0q.c) via C++ |
| bessel\_y1 | full DD | 1884 ulp | 4612 ulp | 1148 ulp | **6.9×** | **6.3×** | **7.2×** | piecewise rational + Hankel asymptotic (j1q.c) via C++ |
| bessel\_yn(3,.) | full DD | 102142 ulp | 10393 ulp | 5641 ulp | **7.2×** | **6.1×** | **6.9×** | forward recurrence from y0/y1 |

### Complex arithmetic

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| cdd\_add | full DD | 0.7 ulp (re) / 0.7 ulp (im) | 0.5 ulp (re) / 0.5 ulp (im) | 0.5 ulp (re) / 0.5 ulp (im) | 1.1× | **3.2×** | **4.4×** | original: component-wise DD add |
| cdd\_sub | full DD | 0.7 ulp (re) / 0.7 ulp (im) | 0.2 ulp (re) / 0.2 ulp (im) | 0.2 ulp (re) / 0.2 ulp (im) | **2.5×** | **5.2×** | **4.5×** | original: component-wise DD sub |
| cdd\_mul | full DD | 1231 ulp (re) / 102 ulp (im) | exact (re) / 0.8 ulp (im) | exact (re) / 0.8 ulp (im) | **6.6×** | **3.4×** | **5.1×** | original: (ac−bd, ad+bc) via DD ops |
| cdd\_div | full DD / deriv | 227 ulp (re) / 2.6 ulp (im) | 1.8 ulp (re) / 6031582967933934 ulp (im) | 1.7 ulp (re) / 2892068785384692 ulp (im) | **5.4×** | **2.6×** | **3.1×** | original: (ac+bd, bc−ad)/(c²+d²) |
| cdd\_conjg | exact | exact | exact | exact | **2.0×** | **2.7×** | **3.0×** | original: negate im limbs |
| cdd\_abs | full DD | 2.0 ulp | 1.5 ulp | 1.5 ulp | **9.3×** | **13×** | **14×** | original: hypot(re, im) |

### Complex transcendentals

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| cdd\_sqrt | full DD | 3.0 ulp (re) / 2.9 ulp (im) | 2.1 ulp (re) / 1.8 ulp (im) | 2.1 ulp (re) / 1.8 ulp (im) | **7.6×** | **12×** | **11×** | original: Kahan-style (\|z\|+\|a\|)/2 with scaling |
| cdd\_exp | full DD | 2.2 ulp (re) / 2.2 ulp (im) | 80 ulp (re) / 803 ulp (im) | 15 ulp (re) / 212 ulp (im) | **3.7×** | **4.2×** | **4.4×** | original: exp(re)·(cos(im), sin(im)) |
| cdd\_log | full DD | 1.8 ulp (re) / 1.6 ulp (im) | 2.2 ulp (re) / 0.7 ulp (im) | 2.2 ulp (re) / 0.7 ulp (im) | **5.5×** | **5.2×** | **6.1×** | original: (log(\|z\|), atan2(im,re)) |
| cdd\_sin | full DD | 2.8 ulp (re) / 8.0 ulp (im) | 8.5 ulp (re) / 8.2 ulp (im) | 8.5 ulp (re) / 8.2 ulp (im) | **4.3×** | **4.7×** | **5.0×** | original: sin(re)cosh(im), cos(re)sinh(im) |
| cdd\_cos | full DD | 2.2 ulp (re) / 7.3 ulp (im) | 8.2 ulp (re) / 9.0 ulp (im) | 8.2 ulp (re) / 9.0 ulp (im) | **4.5×** | **4.7×** | **5.0×** | original: cos(re)cosh(im), −sin(re)sinh(im) |
| cdd\_tan | full DD | 589001174890044 ulp (re) / 7.9 ulp (im) | 2.3 ulp (re) / 8.8 ulp (im) | 2.3 ulp (re) / 8.8 ulp (im) | **4.4×** | **4.7×** | **5.1×** | original: complex sin/cos ratio |
| cdd\_sinh | full DD | 5.5 ulp (re) / 2.8 ulp (im) | 88 ulp (re) / 803 ulp (im) | 15 ulp (re) / 212 ulp (im) | **4.8×** | **4.8×** | **5.0×** | original: sinh(re)cos(im), cosh(re)sin(im) |
| cdd\_cosh | full DD | 2.4 ulp (re) / 5.7 ulp (im) | 88 ulp (re) / 803 ulp (im) | 16 ulp (re) / 212 ulp (im) | **4.7×** | **4.5×** | **5.0×** | original: cosh(re)cos(im), sinh(re)sin(im) |
| cdd\_tanh | full DD | 5.6 ulp (re) / 376400957424568 ulp (im) | 4.2 ulp (re) / 3.6 ulp (im) | 4.2 ulp (re) / 3.6 ulp (im) | **4.6×** | **4.8×** | **5.2×** | original: complex tanh via sinh/cosh |
| cdd\_asin | deriv / full DD | 5290 ulp (re) / 2680523 ulp (im) | 1.3 ulp (re) / 8.9 ulp (im) | 1.3 ulp (re) / 7.0 ulp (im) | **5.1×** | **6.3×** | **6.6×** | original: −i·log(iz+√(1−z²)) |
| cdd\_acos | full DD | 33863511 ulp (re) / 2680523 ulp (im) | 0.6 ulp (re) / 8.9 ulp (im) | 0.6 ulp (re) / 7.0 ulp (im) | **4.9×** | **6.2×** | **6.5×** | original: π/2 − asin(z) |
| cdd\_atan | full DD | 6.1 ulp (re) / 237020239 ulp (im) | 1.4 ulp (re) / 23 ulp (im) | 1.4 ulp (re) / 23 ulp (im) | **5.4×** | **5.7×** | **5.7×** | original: (i/2)·log((i+z)/(i−z)) |
| cdd\_asinh | deriv / full DD | 2.9 ulp (re) / 2.9 ulp (im) | 10728989032 ulp (re) / 1.4 ulp (im) | 1.4 ulp (re) / 2.0 ulp (im) | **4.7×** | **6.1×** | **6.4×** | original: log(z+√(z²+1)) |
| cdd\_acosh | full DD | 2810331 ulp (re) / 2.6 ulp (im) | 12 ulp (re) / 1.3 ulp (im) | 9.4 ulp (re) / 1.3 ulp (im) | **4.4×** | **5.4×** | **6.1×** | original: log(z+√(z²−1)) |
| cdd\_atanh | deriv / full DD | 3.3 ulp (re) / 1.9 ulp (im) | 2436039088 ulp (re) / 1.1 ulp (im) | 1.3 ulp (re) / 1.1 ulp (im) | **5.2×** | **6.5×** | **5.7×** | original: ½·log((1+z)/(1−z)) |

### Array reductions

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| arr\_sum (n=8) | full DD | — | 269 ulp | 269 ulp | **4.3×** | **5.6×** | **5.8×** | original: chained DD add |
| arr\_product (n=8) | full DD | — | 0.0 ulp | 0.0 ulp | **2.0×** | 2.0× | **2.1×** | original: chained DD mul |
| arr\_maxval (n=8) | full DD | — | 0.2 ulp | 0.2 ulp | **3.2×** | **4.6×** | **5.4×** | original: chained DD compare |
| arr\_minval (n=8) | full DD | — | 0.2 ulp | 0.2 ulp | **3.3×** | **3.1×** | **6.8×** | original: chained DD compare |
| arr\_dot (n=8) | full DD | — | 8.7 ulp | 8.7 ulp | **4.2×** | **3.7×** | **4.5×** | original: fused multiply-accumulate with periodic renormalization |
| arr\_norm2 (n=8) | full DD | — | 1.3 ulp | 1.3 ulp | **5.7×** | **6.9×** | **6.9×** | original: sqrt(dot(x,x)) |
| arr\_matmul (8×8·8) | full DD | — | 288 ulp | 288 ulp | **2.5×** | 0.52× | 0.65× | original: AXPY-order C kernel, MR=8 register-blocked panels + 1..7 tail, periodic renorm |

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
