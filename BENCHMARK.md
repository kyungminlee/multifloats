# Benchmark Results

Comparison of multifloats double-double (DD) arithmetic against quad
precision (`real(16)` / `__float128` via libquadmath). The multifloats
implementations use explicit two-limb error-free transformations (EFTs)
on pairs of `double` values, achieving ~106 bits of significand — close
to quad precision's 113 bits — at a fraction of the cost.

## Systems

Three systems are benchmarked; short names are used as column labels
throughout.

| Short name | CPU | OS | Compiler | Build |
|---|---|---|---|---|
| **M1 Max** | Apple M1 Max (ARM64, 10 cores) | macOS 26.3 (Darwin 25.3.0) | GNU Fortran / g++ 15.2.0 (Homebrew GCC 15.2.0\_1) | CMake 4.3.1, `-O3 -flto`, OBJECT library |
| **Skylake** | Intel Xeon family 6 model 85 (Skylake-SP / Cascade Lake), 2.8 GHz, 16 cores, AVX-512 (KVM), 22 GB | Ubuntu 24.04.4 LTS | GNU Fortran / g++ 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) | CMake 3.28.3, `-O3 -flto`, OBJECT library |
| **Raptor Lake** | Intel Core i3-1315U (Raptor Lake), 4.5 GHz boost, 6 cores / 8 threads, 16 GB | Pop!\_OS 24.04 LTS (Linux 6.17.9) | GNU Fortran / g++ 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) | CMake 3.28.3, `-O3 -flto`, OBJECT library |

The x86-64 speedup numbers are generally higher than M1 for many
operations because gfortran's software `libquadmath` on x86-64 is 2–5×
slower per operation than Apple's ARM64 quad path, inflating the
multifloats speedup relative to M1 even though the DD kernels themselves
run at similar speed on all three machines.

## Precision key

Precision is measured as the maximum relative error vs the quad-precision
(`real(16)`) reference over ~1M random inputs:

| Label | max\_rel | Meaning |
|---|---|---|
| **full DD** | ~1e-32 | Full double-double EFT kernel (~106 bits) |
| **exact** | 0.0 | Bit-exact (no rounding involved) |
| **deriv-corrected** | ~1e-25 to 1e-18 | `f(hi) + f'(hi)*lo` correction gives near-DD |
| **single-double** | ~1e-16 to 1e-14 | Leading-limb libm call, no lo correction |

## Origin key

| Tag | Meaning |
|---|---|
| **Julia** | Ported from `external/MultiFloats.jl/src/` |
| **original** | Developed for this project |
| **sample** | Adapted from `external/float64x2-sample.cpp` |

## Fortran: `float64x2` vs `real(16)`

Each operation is timed over 1024 elements × 400 repetitions (fast ops)
or fewer reps (transcendentals), with a NOINLINE drain after each rep to
prevent dead-code elimination. **×** = speedup (`qp_time / mf_time`,
values > 1× mean multifloats is faster); **err** = max\_rel from the
1M-input fuzz run; **prec** = precision label.

### Arithmetic

| op | approach | prec | M1 Max × | M1 Max err | Skylake × | Skylake err | Raptor Lake × | Raptor Lake err |
|---|---|---|---|---|---|---|---|---|
| add | Julia: two\_sum EFT | full DD | 1.2× | 1.5e-32 | **3.6×** | 1.5e-32 | **4.4×** | 1.5e-32 |
| sub | Julia: two\_sum EFT (negate + add) | full DD | **2.0×** | 5.7e-33 | **3.1×** | 6.2e-33 | **5.6×** | 6.2e-33 |
| mul | Julia: two\_prod EFT via FMA | full DD | **11×** | 2.1e-32 | **5.9×** | 3.1e-32 | **6.8×** | 3.3e-32 |
| div | original: Newton refinement (1/y seed, one step) | full DD | **3.6×** | 3.8e-32 | **3.2×** | 6.1e-32 | 1.7× | 6.2e-32 |
| sqrt | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) | full DD | **36×** | 4.5e-32 | **17×** | 5.4e-32 | **39×** | 5.7e-32 |
| add (mf+dp) | Julia: two\_sum EFT | exact | **2.0×** | exact | **4.3×** | exact | **4.3×** | exact |
| mul (dp\*mf) | Julia: two\_prod EFT via FMA | full DD | **8.8×** | 2.1e-32 | **6.1×** | 3.1e-32 | **7.1×** | 2.8e-32 |

### Unary

| op | approach | prec | M1 Max × | M1 Max err | Skylake × | Skylake err | Raptor Lake × | Raptor Lake err |
|---|---|---|---|---|---|---|---|---|
| abs | original: sign-check + negate limbs | exact | 1.1× | exact | **2.4×** | exact | **2.8×** | exact |
| neg | original: negate both limbs | exact | **2.4×** | exact | **4.4×** | exact | **5.2×** | exact |
| aint | original: truncate hi, check DD fractional part | exact | 1.8× | exact | 1.4× | exact | 1.3× | exact |
| anint | original: truncate hi, DD fractional part vs ±0.5 | exact | **2.3×** | exact | **2.5×** | exact | **3.5×** | exact |
| fraction | original: scale both limbs by −exponent | exact | 1.5× | exact | 1.3× | exact | **2.4×** | exact |
| scale | original: ldexp on both limbs | exact | 1.0× | exact | **3.5×** | exact | **5.9×** | exact |
| set\_exponent | original: scale + set\_exponent on hi | exact | 1.7× | exact | **3.3×** | exact | **6.1×** | exact |

### Binary

| op | approach | prec | M1 Max × | M1 Max err | Skylake × | Skylake err | Raptor Lake × | Raptor Lake err |
|---|---|---|---|---|---|---|---|---|
| min | original: DD comparison + select | full DD | **2.2×** | 5.5e-33 | **4.2×** | 5.9e-33 | **4.9×** | 6.1e-33 |
| max | original: DD comparison + select | full DD | **2.3×** | 5.3e-33 | **4.3×** | 6.0e-33 | **5.1×** | 6.1e-33 |
| min3 | original: chained min | full DD | 1.8× | 5.7e-33 | **6.9×** | 5.8e-33 | **3.2×** | 5.8e-33 |
| max3 | original: chained max | full DD | **2.4×** | 6.1e-33 | **6.4×** | 5.7e-33 | **5.9×** | 5.9e-33 |
| sign | original: sign-check + negate | exact | 1.2× | exact | **2.3×** | exact | **2.3×** | exact |
| dim | original: DD comparison, then subtract or zero | full DD | 1.6× | 5.5e-33 | **4.3×** | 5.9e-33 | **3.3×** | 6.1e-33 |
| hypot | original: scaled sqrt(x²+y²) | full DD | **6.1×** | 6.4e-32 | **5.4×** | 7.2e-32 | **6.9×** | 7.7e-32 |
| mod | sample: floor-multiple reduction loop; fallback to div chain | full DD | 0.59× | 5.0e-33 | 1.1× | 3.2e-32 | 0.84× | 2.0e-32 |
| modulo | original: mod + sign adjustment | full DD | 1.2× | 5.0e-33 | 1.6× | 3.2e-32 | 1.5× | 2.0e-32 |

### Exponential / logarithmic

| op | approach | prec | M1 Max × | M1 Max err | Skylake × | Skylake err | Raptor Lake × | Raptor Lake err |
|---|---|---|---|---|---|---|---|---|
| exp | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction | full DD | **2.7×** | 4.3e-31 | **2.7×** | 2.8e-30 | **4.0×** | 3.4e-30 |
| log | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) | full DD | **4.5×** | 1.6e-32 | **3.8×** | 3.0e-32 | **5.7×** | 3.1e-32 |
| log10 | Julia: log2 kernel × DD log10(2) | full DD | **5.8×** | 1.5e-32 | **4.9×** | 3.4e-32 | **7.1×** | 2.9e-32 |
| pow | Julia: exp(y × log(x)) | full DD | **4.2×** | 5.2e-32 | **4.1×** | 2.2e-30 | **5.0×** | 1.9e-30 |
| pow\_int | original: repeated squaring via DD mul | full DD | **18×** | 1.7e-32 | **5.6×** | 2.2e-32 | **7.4×** | 2.1e-32 |

### Trigonometric

| op | approach | prec | M1 Max × | M1 Max err | Skylake × | Skylake err | Raptor Lake × | Raptor Lake err |
|---|---|---|---|---|---|---|---|---|
| sin | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split | full DD | 1.9× | 2.4e-32 | 1.5× | 3.8e-32 | **2.2×** | 4.2e-32 |
| cos | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split | full DD | 1.9× | 3.5e-32 | 1.5× | 4.2e-32 | **2.1×** | 3.6e-32 |
| sinpi | Julia: sinpi Horner polynomial, direct | full DD | **3.4×** | 4.9e-27 | **2.8×** | 4.9e-27 | **3.9×** | 4.9e-27 |
| cospi | Julia: cospi Horner polynomial, direct | full DD | **3.5×** | 8.2e-27 | **2.6×** | 8.2e-27 | **3.9×** | 8.2e-27 |
| tan | original: sin/cos Taylor kernels + DD divide | full DD | 0.89× | 4.7e-32 | 0.8× | 6.7e-32 | 1.1× | 6.2e-32 |
| asin | original: Newton step on sin, seeded by libm asin(hi) | full DD | 1.4× | 6.4e-33 | 1.2× | 4.3e-32 | 1.6× | 3.3e-32 |
| acos | original: Newton step on cos, seeded by libm acos(hi) | full DD | 1.5× | 3.6e-33 | 1.2× | 5.4e-32 | 1.6× | 1.1e-32 |
| atan | original: Newton on tan + atan(x)=π/2·sign(x)−atan(1/x) for \|x\|>1 | full DD | 0.77× | 3.6e-32 | 0.7× | 4.2e-32 | 0.91× | 5.5e-32 |
| atan2 | original: Newton step on atan + quadrant correction | full DD | 0.80× | 1.8e-32 | 0.8× | 3.0e-32 | 0.94× | 3.4e-32 |

All M1 Max values are post-fix (13-term Taylor + `atan(x) = π/2·sign(x) − atan(1/x)`).

### Hyperbolic

| op | approach | prec | M1 Max × | M1 Max err | Skylake × | Skylake err | Raptor Lake × | Raptor Lake err |
|---|---|---|---|---|---|---|---|---|
| sinh | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 | full DD | **2.3×** | 4.4e-31 | 1.9× | 3.7e-30 | **2.4×** | 3.4e-30 |
| cosh | original: (exp+exp⁻¹)/2 | full DD | 1.5× | 4.3e-31 | 1.6× | 3.7e-30 | **2.1×** | 3.4e-30 |
| tanh | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) | full DD | **2.4×** | 3.7e-31 | **2.3×** | 1.2e-30 | **2.7×** | 6.7e-31 |
| asinh | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton | full DD | **5.4×** | 1.4e-30 | **5.8×** | 2.7e-30 | **8.7×** | 3.9e-30 |
| acosh | original: log(x+√(x²−1)) with Newton correction | full DD | **5.4×** | 7.4e-33 | **4.9×** | 4.1e-32 | **7.1×** | 8.6e-32 |
| atanh | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) | full DD | **4.5×** | 4.0e-31 | **4.3×** | 1.2e-30 | **5.8×** | 8.6e-31 |

### Error / special functions

| op | approach | prec | M1 Max × | M1 Max err | Skylake × | Skylake err | Raptor Lake × | Raptor Lake err |
|---|---|---|---|---|---|---|---|---|
| erf | original: Taylor series (\|x\|<2) + libm erf(hi) deriv correction | deriv-corrected | **5.1×** | 3.2e-20 | **5.4×** | 7.8e-19 | **8.7×** | 3.1e-19 |
| erfc | original: 1−erf(x), or libm erfc(hi) for large x | deriv-corrected | **5.0×** | 2.4e-16 | **5.5×** | 2.5e-16 | **8.4×** | 1.7e-16 |
| erfc\_scaled | original: libm erfc\_scaled(hi), no lo correction | single-double | **185×** | 4.3e-16 | **148×** | 5.7e-16 | **193×** | 6.3e-16 |
| gamma | original: libm gamma(hi), no lo correction | single-double | **111×** | 8.4e-17 | **42×** | 4.1e-16 | **55×** | 3.2e-16 |
| log\_gamma | original: libm log\_gamma(hi), no lo correction | single-double | **69×** | 1.2e-16 | **40×** | 2.1e-16 | **53×** | 2.6e-16 |
| bessel\_j0 | original: libm bessel\_j0(hi), no lo correction | single-double | **99×** | 1.9e-16 | **68×** | 8.2e-16 | **94×** | 3.6e-16 |
| bessel\_j1 | original: libm bessel\_j1(hi), no lo correction | single-double | **102×** | 2.7e-16 | **69×** | 1.7e-13 | **96×** | 1.9e-15 |
| bessel\_jn(3,.) | original: libm bessel\_jn(3,hi), no lo correction | single-double | **104×** | 3.1e-16 | **67×** | 7.4e-15 | **84×** | 1.1e-14 |
| bessel\_y0 | original: libm bessel\_y0(hi), no lo correction | single-double | **136×** | 1.9e-16 | **75×** | 4.9e-16 | **113×** | 5.0e-16 |
| bessel\_y1 | original: libm bessel\_y1(hi), no lo correction | single-double | **157×** | 2.4e-16 | **71×** | 1.7e-15 | **47×** | 4.6e-15 |
| bessel\_yn(3,.) | original: libm bessel\_yn(3,hi), no lo correction | single-double | **144×** | 1.9e-15 | **75×** | 1.3e-14 | **111×** | 6.7e-15 |

### Complex arithmetic

| op | approach | prec | M1 Max × | M1 Max err | Skylake × | Skylake err | Raptor Lake × | Raptor Lake err |
|---|---|---|---|---|---|---|---|---|
| cx\_add | original: component-wise DD add | full DD | **2.4×** | 1.3e-32 | **3.6×** | 1.3e-32 | **2.9×** | 1.5e-32 |
| cx\_sub | original: component-wise DD sub | full DD | **2.2×** | 5.8e-33 | **3.4×** | 5.8e-33 | **2.8×** | 5.8e-33 |
| cx\_mul | original: (ac−bd, ad+bc) via DD ops | full DD | **6.8×** | 9.7e-33 | **4.0×** | 2.0e-32 | **4.6×** | 1.9e-32 |
| cx\_div | original: (ac+bd, bc−ad)/(c²+d²) | full DD / deriv | **7.2×** | 2.2e-32 (re) / 8.6e-17 (im) | **4.1×** | 4.6e-32 (re) / 1.1e-16 (im) | **4.3×** | 4.9e-32 (re) / 2.8e-16 (im) |
| cx\_conjg | original: negate im limbs | exact | 1.9× | exact | **3.5×** | exact | **3.9×** | exact |
| cx\_abs | original: hypot(re, im) | full DD | **3.9×** | 6.7e-32 | **4.9×** | 6.6e-32 | **5.9×** | 7.1e-32 |

### Complex transcendentals

| op | approach | prec | M1 Max × | M1 Max err | Skylake × | Skylake err | Raptor Lake × | Raptor Lake err |
|---|---|---|---|---|---|---|---|---|
| cx\_sqrt | original: Kahan-style (\|z\|+\|a\|)/2 with scaling | full DD | **5.7×** | 4.9e-32 | **4.6×** | 6.4e-32 | **7.0×** | 6.3e-32 |
| cx\_exp | original: exp(re)·(cos(im), sin(im)) | full DD | 1.5× | 1.0e-31 | 1.5× | 7.8e-31 | **2.4×** | 1.3e-30 |
| cx\_log | original: (log(\|z\|), atan2(im,re)) | full DD | 1.9× | 1.7e-32 | 1.8× | 2.0e-29 (re) / 2.8e-32 (im) | **2.5×** | 4.2e-32 |
| cx\_sin | original: sin(re)cosh(im), cos(re)sinh(im) | full DD | 1.5× | 9.5e-32 | 1.4× | 4.9e-31 | **2.0×** | 8.5e-31 |
| cx\_cos | original: cos(re)cosh(im), −sin(re)sinh(im) | full DD | 1.5× | 9.8e-32 | 1.4× | 4.9e-31 | **2.0×** | 8.5e-31 |
| cx\_tan | original: complex sin/cos ratio | full DD | 0.80× | 8.0e-31 | 0.7× | 3.1e-30 | 1.0× | 2.3e-30 |
| cx\_sinh | original: sinh(re)cos(im), cosh(re)sin(im) | full DD | 1.6× | 5.9e-32 | 1.5× | 1.2e-30 | **2.0×** | 1.4e-30 |
| cx\_cosh | original: cosh(re)cos(im), sinh(re)sin(im) | full DD | 1.6× | 7.3e-32 | 1.5× | 1.2e-30 | **2.0×** | 1.4e-30 |
| cx\_tanh | original: complex tanh via sinh/cosh | full DD | 0.84× | 8.2e-31 | 0.8× | 1.1e-30 | 1.1× | 4.2e-30 |
| cx\_asin | original: −i·log(iz+√(1−z²)) | deriv / full DD | **2.2×** | 2.3e-27 (re) / 1.4e-31 (im) | **2.4×** | 2.8e-23 (re) / 8.4e-31 (im) | **3.2×** | 4.0e-23 (re) / 1.0e-30 (im) |
| cx\_acos | original: π/2 − asin(z) | full DD | **2.3×** | 2.1e-33 (re) / 1.4e-31 (im) | **2.4×** | 2.9e-32 (re) / 8.4e-31 (im) | **3.0×** | 3.8e-32 (re) / 1.0e-30 (im) |
| cx\_atan | original: (i/2)·log((i+z)/(i−z)) | full DD | 1.6× | 1.0e-32 (re) / 1.1e-32 (im) | 1.4× | 3.8e-32 (re) / 4.8e-31 (im) | 1.8× | 2.1e-31 |
| cx\_asinh | original: log(z+√(z²+1)) | deriv / full DD | **2.2×** | 2.7e-26 (re) / 2.0e-32 (im) | **2.2×** | 6.7e-22 (re) / 6.6e-32 (im) | **2.9×** | 1.3e-22 (re) / 7.3e-32 (im) |
| cx\_acosh | original: log(z+√(z²−1)) | full DD | **2.1×** | 4.4e-32 (re) / 5.8e-33 (im) | **2.0×** | 7.2e-31 (re) / 2.9e-32 (im) | **2.6×** | 1.7e-30 (re) / 1.7e-32 (im) |
| cx\_atanh | original: ½·log((1+z)/(1−z)) | deriv / full DD | 1.7× | 1.5e-26 (re) / 2.5e-33 (im) | 1.5× | 7.2e-23 (re) / 6.7e-32 (im) | 1.9× | 1.5e-22 (re) / 4.5e-32 (im) |

### Array reductions

| op | approach | prec | M1 Max × | M1 Max err | Skylake × | Skylake err | Raptor Lake × | Raptor Lake err |
|---|---|---|---|---|---|---|---|---|
| arr\_sum (n=8) | original: chained DD add | full DD | 1.1× | 1.1e-32 | 1.4× | 2.5e-31 | **2.1×** | 5.5e-31 |
| arr\_product (n=8) | original: chained DD mul | full DD | **3.2×** | 1.1e-53 | **2.2×** | 4.3e-49 | **3.0×** | 5.8e-54 |
| arr\_maxval (n=8) | original: chained DD compare | full DD | **3.1×** | 2.9e-33 | **5.8×** | 6.1e-33 | **7.2×** | 6.1e-33 |
| arr\_minval (n=8) | original: chained DD compare | full DD | **3.4×** | 2.7e-33 | **5.2×** | 6.0e-33 | **5.5×** | 6.0e-33 |
| arr\_dot (n=8) | original: fused multiply-accumulate with periodic renormalization | full DD | **4.3×** | 2.5e-32 | **4.6×** | 1.0e-31 | **6.3×** | 1.1e-30 |
| arr\_norm2 (n=8) | original: sqrt(dot(x,x)) | full DD | **4.3×** | 1.4e-32 | **6.8×** | 6.4e-32 | **5.7×** | 5.8e-32 |
| arr\_matmul (8×8\*8) | original: fused multiply-accumulate with periodic renormalization | full DD | **2.1×** | 5.8e-32 | 0.8× | 8.8e-30 | 0.96× | 2.0e-29 |

## C++: `MultiFloat<double,2>` vs `__float128`

Header-only — all kernels inline into the call site. No LTO needed.
Precision characteristics are identical to the Fortran version (same
algorithms); see the Fortran tables for `max_rel` and `prec` of
corresponding operations. Only speedup is shown per system.

### Arithmetic

| op | approach | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|
| add | Julia: two\_sum EFT | **2.6×** | **4.1×** | **2.8×** |
| sub | Julia: two\_sum EFT (negate + add) | **2.8×** | **4.2×** | **7.2×** |
| mul | Julia: two\_prod EFT via FMA | **12×** | **7.7×** | **7.2×** |
| div | original: Newton refinement (1/y seed, one step) | 1.0× | **2.6×** | 1.3× |
| sqrt | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) | **52×** | **29×** | **56×** |
| cbrt | original: Newton correction on cbrt(hi) seed | **45×** | **13×** | **16×** |
| fma | original: x\*y + z via DD ops | **63×** | **69×** | **133×** |
| abs | original: sign-check + negate limbs | **2.0×** | **5.1×** | **5.8×** |
| neg | original: negate both limbs | **2.1×** | **3.6×** | **3.9×** |

### Rounding

| op | approach | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|
| floor | original: floor hi, adjust lo | **3.2×** | **4.2×** | **6.1×** |
| ceil | original: ceil hi, adjust lo | **3.2×** | **3.9×** | **5.5×** |
| trunc | original: signbit ? −floor(−x) : floor(x) | **2.5×** | **3.2×** | **4.5×** |
| round | original: trunc(x + ½·sign(x)) | 0.98× | 1.0× | 1.1× |
| rint | original: nearbyint on hi, adjust lo | **11×** | **11×** | **19×** |
| nearbyint | original: nearbyint on hi, adjust lo | **36×** | **65×** | **138×** |

### Binary

| op | approach | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|
| fmin | original: DD comparison + select | **5.5×** | **9.8×** | **9.9×** |
| fmax | original: DD comparison + select | **5.3×** | **13×** | **10×** |
| fdim | original: DD comparison, then subtract or zero | **6.1×** | **6.9×** | **8.4×** |
| copysign | original: sign-bit copy to hi, propagate to lo | 1.9× | **5.2×** | **7.2×** |
| fmod | sample: floor-multiple reduction loop; fallback to div chain | 0.85× | 0.76× | 0.82× |
| hypot | original: scaled sqrt(x²+y²) | **36×** | **20×** | **38×** |
| ldexp(.,5) | original: ldexp on both limbs | **2.0×** | **2.4×** | **4.1×** |

### Exponential / logarithmic

| op | approach | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|
| exp | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction | **3.1×** | **2.7×** | **3.9×** |
| exp2 | Julia: exp2 polynomial (14-term Horner) | **3.3×** | **3.1×** | **4.2×** |
| expm1 | original: exp(x) − 1 via DD sub | **4.2×** | **3.3×** | **3.9×** |
| log | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) | **4.4×** | **4.0×** | **5.1×** |
| log10 | Julia: log2 kernel × DD log10(2) | **5.9×** | **5.1×** | **6.9×** |
| log2 | Julia: log2 table lookup + polynomial | **5.6×** | **4.7×** | **6.5×** |
| log1p | original: log(1 + x) via DD add | **4.9×** | **4.3×** | **4.8×** |
| pow | Julia: exp(y × log(x)) | **4.3×** | **4.1×** | **4.9×** |

### Trigonometric

| op | approach | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|
| sin | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split | **2.8×** | **2.3×** | **3.0×** |
| cos | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split | **2.9×** | **2.3×** | **3.1×** |
| tan | original: sin/cos Taylor kernels + DD divide | 1.3× | 1.1× | 1.3× |
| asin | original: Newton step on sin, seeded by libm asin(hi) | 1.9× | 1.7× | **2.2×** |
| acos | original: Newton step on cos, seeded by libm acos(hi) | 1.9× | 1.7× | **2.3×** |
| atan | original: Newton on tan + atan(x)=π/2·sign(x)−atan(1/x) for \|x\|>1 | 1.1× | 1.0× | 1.2× |
| atan2 | original: Newton step on atan + quadrant correction | 1.0× | 1.1× | 1.2× |

### Hyperbolic

| op | approach | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|
| sinh | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 | **2.4×** | 1.9× | **2.3×** |
| cosh | original: (exp+exp⁻¹)/2 | 1.6× | 1.5× | **2.1×** |
| tanh | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) | **2.4×** | **2.2×** | **2.7×** |
| asinh | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton | **6.2×** | **6.0×** | **8.3×** |
| acosh | original: log(x+√(x²−1)) with Newton correction | **5.7×** | **5.5×** | **7.6×** |
| atanh | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) | **4.6×** | **4.4×** | **5.5×** |

### Error / special functions

| op | approach | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|
| erf | original: Taylor series (\|x\|<2) + libm erf(hi) deriv correction | **4.8×** | **5.7×** | **7.9×** |
| erfc | original: 1−erf(x), or libm erfc(hi) for large x | **5.1×** | **5.5×** | **7.1×** |
| tgamma | original: libm tgamma(hi), no lo correction | **103×** | **59×** | **84×** |
| lgamma | original: libm lgamma(hi), no lo correction | **72×** | **46×** | **54×** |

## Notes

- **`mod` / `fmod`** is the only operation where quad precision is
  consistently faster. The DD `mod` uses a floor-multiple reduction loop
  (adapted from `external/float64x2-sample.cpp`) for small quotients, or a
  full DD divide chain for large quotients; libquadmath's `fmodq` uses a
  specialized bit-level remainder algorithm. `modulo` (Fortran) now beats qp
  at 1.3–1.6× thanks to the iterative approach. Precision degrades as
  `~10^(log10(quotient) - 31)` for large quotients, which is the inherent DD
  precision limit.

- **Trig range reduction** uses a 3-part π/2 constant (~161 bits) via
  Cody–Waite subtraction with DD arithmetic (FMA-captured product errors).
  Combined with the π/8 argument split (which halves the polynomial
  evaluation range from x² ≤ 0.616 to x² ≤ 0.154), this gives full DD
  precision (~4e-32) for sin/cos/tan with the current 13-term Taylor
  kernels. The earlier 10-term kernels truncated at ~1e-23 at the boundary,
  which was below libm noise on Apple M1 but surfaced as ~1e-27 sin/cos
  max\_rel on x86-64; extending to 13 terms brings both targets to the
  DD floor. For |x| > ~1e15, a Payne–Hanek reduction with a multi-word
  2/π table would be needed.

- **x86-64 specifics.** On the Skylake-SP + glibc run, three issues that
  are masked on Apple M1 needed fixes to reach full DD precision:
  (1) `mod` used `floor` which returns default INTEGER and saturates for
  quotients > 2³¹, making the reduction loop run `|x/y|` iterations
  instead of one — one `mod(6.57e11, -44.65)` call was taking ~5 minutes.
  Fix: switch to `aint` (real-returning). (2) The 10-term sin/cos Taylor
  truncated at ~1e-23 at x² ≈ (π/8)², capping Newton-seeded asin/acos
  at ~1e-26. Fix: extend to 13 terms (Apple M1 was silently at ~5e-27 for
  the same reason but within the documented "near-DD" tolerance).
  (3) `atan(x)` for |x|>1 had ~log₂(x) bits of cancellation in the
  `x − tan(y₀)` Newton residual, giving ~1e-17 at |x|=10¹⁵. Fix: use
  `atan(x) = sign(x)·π/2 − atan(1/x)` so the Newton argument stays
  in [−1, 1]. All three fixes are in both the Fortran and C++ codebases
  and benefit both platforms.

- **Single-double precision** functions (gamma, bessel, erfc\_scaled, etc.)
  achieve 60–193× speedup by evaluating `f(hi)` via the leading-limb libm
  intrinsic without a lo-limb correction. This gives `double`-level accuracy
  (~15 digits) rather than full DD (~31 digits). For applications needing
  full DD precision on these functions, a polynomial or series expansion
  would be required (at significant implementation cost and reduced speedup).

- **Newton-corrected** functions (asin, acos, atan, atan2) compute a libm
  seed `f(hi)` and refine with one Newton step using the DD-accurate inverse
  function (sin, cos, tan). Quadratic convergence makes one step sufficient
  for full DD precision.

- **Deriv-corrected** functions (erf, erfc, atan for some ranges) use
  `f(hi) + f'(hi) * lo` to recover most of the DD precision. The correction
  is exact to first order, giving ~20–25 digits in typical cases.

- The **Fortran** multifloats module uses `elemental` functions on a
  `sequence` derived type. gfortran's ABI passes/returns these via hidden
  pointers (not in FP registers), adding ~1.5× overhead vs the C++ header-
  only version even with LTO inlining. See the performance note at the top
  of `fsrc/multifloats.fypp` for details and the `bind(c)` escape hatch.

- **Array reductions** (dot\_product, matmul) use a fused multiply-
  accumulate kernel that computes the product's error-free representation
  and accumulates corrections into a scalar `s_lo`, with periodic
  renormalization (configurable via `mf_set_fma_renorm_interval`).
