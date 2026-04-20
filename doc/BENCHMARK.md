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
| **M1 Max** | Apple M1 Max (arm64), 10 cores | macOS 26.3 (Darwin 25.3.0) | GNU Fortran 15.2.0 (Homebrew GCC 15.2.0_1) / g++ Apple clang version 17.0.0 (clang-1700.6.4.2) | CMake 4.3.1, `-O3` |
| **pop-os** | 13th Gen Intel(R) Core(TM) i3-1315U (x86_64), 8 cores, 15 GB | Pop!_OS 24.04 LTS (Linux 6.17.9-76061709-generic) | GNU Fortran / g++ 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) | CMake 3.28.3, `-O3 -flto`, STATIC library |
| **sandbox** | Intel (family 6, model 207, AVX-512 + AMX), 16 vCPU, 21 GB (hypervisor-masked) | Ubuntu 24.04.4 LTS (Linux 4.4.0) | GNU Fortran / g++ 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) | CMake 3.28.3, `-O3 -flto`, OBJECT library |

## Precision key

Precision is measured as the maximum relative error vs the quad-precision
(`real(16)`) reference over ~1M random inputs (fixed seed 42), reported in
**DD ulps** (1 ulp ≈ 2⁻¹⁰⁵ ≈ 2.46e-32 for the ~106-bit DD significand):

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

## Fortran: `float64x2` vs `real(16)`

Each operation is timed over 1024 elements × 400 repetitions (fast ops)
or fewer reps (transcendentals), with a NOINLINE drain after each rep to
prevent dead-code elimination. **×** = speedup (`qp_time / dd_time`,
values > 1× mean multifloats is faster); **err** = error in DD ulps
(1 ulp ≈ 2⁻¹⁰⁵ ≈ 2.46e-32) from the 1M-input fuzz run; **prec** =
precision label.

### Arithmetic

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| add | full DD | 0.6 ulp | 0.6 ulp | 0.6 ulp | **2.1×** | **5.6×** | **2.9×** | Julia: two\_sum EFT |
| sub | full DD | 0.2 ulp | 0.2 ulp | 0.2 ulp | **2.1×** | **3.4×** | **2.9×** | Julia: two\_sum EFT (negate + add) |
| mul | full DD | 1.4 ulp | 1.4 ulp | 1.4 ulp | **5.2×** | **3.9×** | **4.2×** | Julia: two\_prod EFT via FMA |
| div | full DD | 2.2 ulp | 2.2 ulp | 2.2 ulp | **3.1×** | 1.8× | 1.8× | original: Newton refinement (1/y seed, one step) |
| sqrt | full DD | 0.7 ulp | 0.7 ulp | 2.1 ulp | **30×** | **37×** | **39×** | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) |
| add (dd+dp) | exact | exact | exact | exact | **2.3×** | **4.0×** | **3.8×** | Julia: two\_sum EFT |
| mul (dp\*dd) | full DD | 1.4 ulp | 1.4 ulp | 1.4 ulp | **5.0×** | **4.3×** | **5.3×** | Julia: two\_prod EFT via FMA |

### Unary

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| abs | exact | exact | exact | exact | 1.1× | 1.7× | 1.3× | original: sign-check + negate limbs |
| neg | exact | exact | exact | exact | 1.4× | **2.3×** | **2.0×** | original: negate both limbs |
| aint | exact | exact | exact | exact | 0.83× | 1.2× | **2.1×** | original: truncate hi, check DD fractional part |
| anint | exact | exact | exact | exact | 1.5× | 1.2× | 1.8× | original: truncate hi, DD fractional part vs ±0.5 |
| fraction | exact | exact | exact | exact | 1.1× | 1.5× | 1.4× | original: scale both limbs by −exponent |
| scale | exact | exact | exact | exact | 0.88× | **5.2×** | **3.5×** | original: ldexp on both limbs |
| set\_exponent | exact | exact | exact | exact | 1.5× | **3.7×** | **3.6×** | original: scale + set\_exponent on hi |

### Binary

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| min | full DD | 0.2 ulp | 0.2 ulp | 0.2 ulp | 0.98× | **3.3×** | **2.9×** | original: DD comparison + select |
| max | full DD | 0.2 ulp | 0.2 ulp | 0.2 ulp | 0.95× | **3.1×** | **3.6×** | original: DD comparison + select |
| min3 | full DD | 0.2 ulp | 0.2 ulp | 0.2 ulp | **2.2×** | **3.6×** | **3.8×** | original: chained min |
| max3 | full DD | 0.2 ulp | 0.2 ulp | 0.2 ulp | **2.6×** | **3.7×** | **4.0×** | original: chained max |
| sign | exact | exact | exact | exact | 1.1× | 1.6× | 1.9× | original: sign-check + negate |
| dim | full DD | 0.2 ulp | 0.2 ulp | 0.2 ulp | 1.8× | **4.0×** | **4.5×** | original: DD comparison, then subtract or zero |
| hypot | full DD | 2.5 ulp | 2.5 ulp | 3.2 ulp | **5.8×** | **13×** | **7.2×** | original: scaled sqrt(x²+y²) |
| mod | full DD | 0.8 ulp | 0.8 ulp | 0.8 ulp | 0.51× | 0.91× | 1.4× | sample: floor-multiple reduction loop; fallback to div chain |
| modulo | full DD | 0.8 ulp | 0.8 ulp | 0.8 ulp | 1.1× | 1.3× | **2.1×** | original: mod + sign adjustment |

### Exponential / logarithmic

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| exp | full DD | 365 ulp | 1.5 ulp | 365 ulp | **5.4×** | **5.2×** | **3.7×** | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction |
| log | full DD | 1.8 ulp | 1.8 ulp | 1.8 ulp | **6.4×** | **4.9×** | **5.2×** | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) |
| log10 | full DD | 1.2 ulp | 1.2 ulp | 1.2 ulp | **8.8×** | **6.5×** | **6.6×** | Julia: log2 kernel × DD log10(2) |
| pow | full DD | 301 ulp | 77 ulp | 298 ulp | **6.2×** | **6.1×** | **4.9×** | Julia: exp(y × log(x)) |
| pow\_int | full DD | 1.0 ulp | 1.0 ulp | 1.0 ulp | **7.4×** | **4.1×** | **5.2×** | original: repeated squaring via DD mul |

### Trigonometric

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| sin | full DD | 2.0 ulp | 2.0 ulp | 1.5 ulp | **3.7×** | **3.2×** | **2.1×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| cos | full DD | 1.5 ulp | 1.9 ulp | 1.6 ulp | **3.7×** | **3.2×** | **2.1×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| sinpi | full DD | — | — | — | **4.6×** | **3.9×** | **3.9×** | Julia: sinpi Horner polynomial, direct |
| cospi | full DD | — | — | — | **4.8×** | **4.0×** | **4.1×** | Julia: cospi Horner polynomial, direct |
| tan | full DD | 2.4 ulp | 2.4 ulp | 2.5 ulp | **2.7×** | **2.3×** | 1.1× | original: sin/cos Taylor kernels + DD divide |
| asin | full DD | 0.3 ulp | 0.3 ulp | 0.7 ulp | **5.9×** | **5.7×** | **4.1×** | original: piecewise rational P/Q (3 regions, from libquadmath asinq.c) |
| acos | full DD | 0.5 ulp | 0.5 ulp | 0.8 ulp | **6.0×** | **6.1×** | **3.6×** | original: asin polynomial + half-angle identity |
| atan | full DD | 1.1 ulp | 1.1 ulp | 1.1 ulp | **4.0×** | **4.0×** | **3.5×** | original: 84-entry table lookup + rational P(t²)/Q(t²) (from libquadmath atanq.c) |
| atan2 | full DD | 0.8 ulp | 0.7 ulp | 1.2 ulp | **3.4×** | **3.8×** | **3.0×** | original: table-based atan + quadrant correction |

### Hyperbolic

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| sinh | full DD | 365 ulp | 4.5 ulp | 365 ulp | **4.5×** | **3.0×** | **2.3×** | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 |
| cosh | full DD | 365 ulp | 1.5 ulp | 365 ulp | **3.0×** | **2.6×** | 2.0× | original: (exp+exp⁻¹)/2 |
| tanh | full DD | 164 ulp | 4.6 ulp | 163 ulp | **5.1×** | **3.1×** | **2.9×** | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) |
| asinh | full DD | 868 ulp | 2.7 ulp | 868 ulp | **7.5×** | **10×** | **7.9×** | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton |
| acosh | full DD | 1.3 ulp | 1.4 ulp | 1.4 ulp | **6.6×** | **6.9×** | **6.5×** | original: log(x+√(x²−1)) with Newton correction |
| atanh | full DD | 2398 ulp | 2.5 ulp | 2400 ulp | **7.0×** | **7.8×** | **5.4×** | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) |

### Error / special functions

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| erf | full DD | 0.6 ulp | 0.6 ulp | 0.6 ulp | **5.0×** | **4.4×** | **3.3×** | piecewise rational approx (libquadmath erfq.c) |
| erfc | full DD | 267 ulp | 1.1 ulp | 269 ulp | **5.1×** | **4.7×** | **3.2×** | piecewise rational approx + split exp(-x^2) |
| erfc\_scaled | full DD | 312 ulp | 1.5 ulp | 312 ulp | **6.5×** | **7.0×** | **5.0×** | exp(x^2)·erfc(x) with asymptotic cancellation |
| gamma | full DD | 134209 ulp | 10 ulp | 134237 ulp | **9.0×** | **7.7×** | **5.7×** | piecewise rational approx + Stirling + reflection |
| log\_gamma | full DD | 5908 ulp | 2.1 ulp | 5909 ulp | **5.0×** | **4.6×** | **3.5×** | piecewise rational approx + Stirling asymptotic |
| bessel\_j0 | full DD | 851 ulp | 851 ulp | 567 ulp | **6.5×** | **6.1×** | **4.9×** | piecewise rational + Hankel asymptotic (j0q.c) via C++ |
| bessel\_j1 | full DD | 1725 ulp | 1725 ulp | 945 ulp | **6.4×** | **5.8×** | **4.7×** | piecewise rational + Hankel asymptotic (j1q.c) via C++ |
| bessel\_jn(3,.) | full DD | 2024 ulp | 2024 ulp | 2024 ulp | **4.1×** | **4.5×** | **3.8×** | forward/backward recurrence from j0/j1 |
| bessel\_y0 | full DD | 584 ulp | 584 ulp | 145 ulp | **7.0×** | **6.3×** | **4.9×** | piecewise rational + Hankel asymptotic (j0q.c) via C++ |
| bessel\_y1 | full DD | 4612 ulp | 4612 ulp | 1148 ulp | **7.3×** | **6.3×** | **5.0×** | piecewise rational + Hankel asymptotic (j1q.c) via C++ |
| bessel\_yn(3,.) | full DD | 10453 ulp | 10393 ulp | 5641 ulp | **7.3×** | **6.1×** | **5.0×** | forward recurrence from y0/y1 |

### Complex arithmetic

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| cdd\_add | full DD | 0.5 ulp (re) / 0.5 ulp (im) | 0.5 ulp (re) / 0.5 ulp (im) | 0.5 ulp (re) / 0.5 ulp (im) | **3.2×** | **3.2×** | **5.4×** | original: component-wise DD add |
| cdd\_sub | full DD | 0.2 ulp (re) / 0.2 ulp (im) | 0.2 ulp (re) / 0.2 ulp (im) | 0.2 ulp (re) / 0.2 ulp (im) | **3.2×** | **5.2×** | **5.2×** | original: component-wise DD sub |
| cdd\_mul | full DD | exact (re) / 0.8 ulp (im) | exact (re) / 0.8 ulp (im) | exact (re) / 0.8 ulp (im) | **5.3×** | **3.4×** | **5.3×** | original: (ac−bd, ad+bc) via DD ops |
| cdd\_div | full DD / deriv | 1.8 ulp (re) / 6031582967933934 ulp (im) | 1.8 ulp (re) / 6031582967933934 ulp (im) | 1.8 ulp (re) / 6031582967933934 ulp (im) | **5.0×** | **2.6×** | **3.1×** | original: (ac+bd, bc−ad)/(c²+d²) |
| cdd\_conjg | exact | exact | exact | exact | **2.0×** | **2.7×** | **3.1×** | original: negate im limbs |
| cdd\_abs | full DD | 1.5 ulp | 1.5 ulp | 3.2 ulp | **5.6×** | **13×** | **7.3×** | original: hypot(re, im) |

### Complex transcendentals

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| cdd\_sqrt | full DD | 2.1 ulp (re) / 1.8 ulp (im) | 2.1 ulp (re) / 1.8 ulp (im) | 2.6 ulp (re) / 2.4 ulp (im) | **7.2×** | **12×** | **5.8×** | original: Kahan-style (\|z\|+\|a\|)/2 with scaling |
| cdd\_exp | full DD | 283 ulp (re) / 811 ulp (im) | 80 ulp (re) / 803 ulp (im) | 286 ulp (re) / 278 ulp (im) | **4.4×** | **4.2×** | **2.1×** | original: exp(re)·(cos(im), sin(im)) |
| cdd\_log | full DD | 58 ulp (re) / 0.7 ulp (im) | 2.2 ulp (re) / 0.7 ulp (im) | 57 ulp (re) / 1.3 ulp (im) | **5.0×** | **5.2×** | **3.8×** | original: (log(\|z\|), atan2(im,re)) |
| cdd\_sin | full DD | 310 ulp (re) / 310 ulp (im) | 8.5 ulp (re) / 8.2 ulp (im) | 310 ulp (re) / 310 ulp (im) | **5.6×** | **4.7×** | 1.9× | original: sin(re)cosh(im), cos(re)sinh(im) |
| cdd\_cos | full DD | 310 ulp (re) / 310 ulp (im) | 8.2 ulp (re) / 9.0 ulp (im) | 310 ulp (re) / 310 ulp (im) | **5.9×** | **4.7×** | 1.7× | original: cos(re)cosh(im), −sin(re)sinh(im) |
| cdd\_tan | full DD | 266 ulp (re) / 521 ulp (im) | 2.3 ulp (re) / 8.8 ulp (im) | 63 ulp (re) / 27 ulp (im) | **5.5×** | **4.7×** | 0.97× | original: complex sin/cos ratio |
| cdd\_sinh | full DD | 298 ulp (re) / 811 ulp (im) | 88 ulp (re) / 803 ulp (im) | 298 ulp (re) / 278 ulp (im) | **6.0×** | **4.8×** | 2.0× | original: sinh(re)cos(im), cosh(re)sin(im) |
| cdd\_cosh | full DD | 283 ulp (re) / 811 ulp (im) | 88 ulp (re) / 803 ulp (im) | 286 ulp (re) / 298 ulp (im) | **6.3×** | **4.5×** | 1.9× | original: cosh(re)cos(im), sinh(re)sin(im) |
| cdd\_tanh | full DD | 494 ulp (re) / 259 ulp (im) | 4.2 ulp (re) / 3.6 ulp (im) | 27 ulp (re) / 35 ulp (im) | **6.0×** | **4.8×** | 0.96× | original: complex tanh via sinh/cosh |
| cdd\_asin | deriv / full DD | 1.4 ulp (re) / 8709 ulp (im) | 1.3 ulp (re) / 8.9 ulp (im) | 1031685047 ulp (re) / 8679 ulp (im) | **5.0×** | **6.3×** | **4.4×** | original: −i·log(iz+√(1−z²)) |
| cdd\_acos | full DD | 0.6 ulp (re) / 8709 ulp (im) | 0.6 ulp (re) / 8.9 ulp (im) | 0.5 ulp (re) / 8679 ulp (im) | **4.7×** | **6.2×** | **4.1×** | original: π/2 − asin(z) |
| cdd\_atan | full DD | 1.4 ulp (re) / 486 ulp (im) | 1.4 ulp (re) / 23 ulp (im) | 1.5 ulp (re) / 15 ulp (im) | **5.3×** | **5.7×** | **2.6×** | original: (i/2)·log((i+z)/(i−z)) |
| cdd\_asinh | deriv / full DD | 79778829935 ulp (re) / 1.4 ulp (im) | 10728989032 ulp (re) / 1.4 ulp (im) | 94159058344 ulp (re) / 2.9 ulp (im) | **4.7×** | **6.1×** | **3.7×** | original: log(z+√(z²+1)) |
| cdd\_acosh | full DD | 1078 ulp (re) / 1.3 ulp (im) | 12 ulp (re) / 1.3 ulp (im) | 251 ulp (re) / 0.8 ulp (im) | **3.9×** | **5.4×** | **3.1×** | original: log(z+√(z²−1)) |
| cdd\_atanh | deriv / full DD | 2436039088 ulp (re) / 1.1 ulp (im) | 2436039088 ulp (re) / 1.1 ulp (im) | 17395817069 ulp (re) / 2.3 ulp (im) | **5.8×** | **6.5×** | **3.1×** | original: ½·log((1+z)/(1−z)) |

### Array reductions

| op | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|---|
| arr\_sum (n=8) | full DD | 269 ulp | 269 ulp | 269 ulp | **4.0×** | **5.6×** | **5.6×** | original: chained DD add |
| arr\_product (n=8) | full DD | 0.0 ulp | 0.0 ulp | 0.0 ulp | **2.0×** | 2.0× | 2.0× | original: chained DD mul |
| arr\_maxval (n=8) | full DD | 0.2 ulp | 0.2 ulp | 0.2 ulp | **3.1×** | **4.6×** | **6.7×** | original: chained DD compare |
| arr\_minval (n=8) | full DD | 0.2 ulp | 0.2 ulp | 0.2 ulp | **3.4×** | **3.1×** | **6.7×** | original: chained DD compare |
| arr\_dot (n=8) | full DD | 8.7 ulp | 8.7 ulp | 8.7 ulp | **4.1×** | **3.7×** | **4.4×** | original: fused multiply-accumulate with periodic renormalization |
| arr\_norm2 (n=8) | full DD | 1.3 ulp | 1.3 ulp | 2.1 ulp | **5.4×** | **6.9×** | **6.1×** | original: sqrt(dot(x,x)) |
| arr\_matmul (8×8·8) | full DD | 288 ulp | 288 ulp | 288 ulp | **2.7×** | 0.52× | 0.56× | original: AXPY-order C kernel, MR=8 register-blocked panels + 1..7 tail, periodic renorm |

## C++: `MultiFloat<double,2>` vs `__float128`

Header-only — all kernels inline into the call site. No LTO needed.
See the Fortran tables for `prec` labels.

**Methodology.** Precision and timing are measured separately, mirroring
the Fortran split (`fortran_fuzz` / `fortran_bench`):

- **err** columns come from `cpp_fuzz` at 1M iterations, fixed seed 42.
- **×** columns come from `cpp_bench` at 1024 elements × reps.

### Arithmetic

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| add | 192 ulp | 192 ulp | 192 ulp | **3.1×** | **8.9×** | **4.9×** | Julia: two\_sum EFT |
| sub | 249 ulp | 249 ulp | 249 ulp | **3.2×** | **4.6×** | **4.9×** | Julia: two\_sum EFT (negate + add) |
| mul | 2.2 ulp | 2.2 ulp | 2.2 ulp | **14×** | **7.3×** | **8.1×** | Julia: two\_prod EFT via FMA |
| div | 3.1 ulp | 3.1 ulp | 3.1 ulp | **8.6×** | **4.9×** | **5.6×** | original: Newton refinement (1/y seed, one step) |
| sqrt | 1.8 ulp | 1.8 ulp | 1.8 ulp | **48×** | **53×** | **54×** | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) |
| cbrt | — | — | — | **39×** | **17×** | **17×** | original: Newton correction on cbrt(hi) seed |
| fma | — | — | — | **88×** | **168×** | **161×** | original: x\*y + z via DD ops |
| abs | 0.2 ulp | 0.2 ulp | 0.2 ulp | 1.6× | **7.4×** | **10×** | original: sign-check + negate limbs |
| neg | 0.2 ulp | 0.2 ulp | 0.2 ulp | **4.3×** | **6.4×** | **9.2×** | original: negate both limbs |

### Rounding

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| floor | exact | exact | exact | **3.0×** | **5.6×** | **5.5×** | original: floor hi, adjust lo |
| ceil | exact | exact | exact | **3.1×** | **5.4×** | **5.7×** | original: ceil hi, adjust lo |
| trunc | exact | exact | exact | **2.5×** | **4.7×** | **4.8×** | original: signbit ? −floor(−x) : floor(x) |
| round | exact | exact | exact | **2.1×** | **2.3×** | **2.4×** | original: trunc(x + ½·sign(x)) |
| rint | exact | exact | exact | **10×** | **18×** | **9.6×** | original: nearbyint on hi, adjust lo |
| nearbyint | exact | exact | exact | **33×** | **151×** | **92×** | original: nearbyint on hi, adjust lo |

### Binary

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| fmin | 0.2 ulp | 0.2 ulp | 0.2 ulp | **7.5×** | **8.7×** | **11×** | original: DD comparison + select |
| fmax | 0.2 ulp | 0.2 ulp | 0.2 ulp | **7.1×** | **10×** | **11×** | original: DD comparison + select |
| fdim | 159 ulp | 159 ulp | 159 ulp | **6.8×** | **9.2×** | **8.2×** | original: DD comparison, then subtract or zero |
| copysign | 0.2 ulp | 0.2 ulp | 0.2 ulp | **3.3×** | **5.3×** | **9.4×** | original: sign-bit copy to hi, propagate to lo |
| fmod | 5115223702040951 ulp | 5115223702040951 ulp | 5115223702040951 ulp | 1.1× | 1.1× | 1.2× | sample: floor-multiple reduction loop; fallback to div chain |
| hypot | 2.1 ulp | 2.1 ulp | 1.6 ulp | **6.4×** | **13×** | **38×** | original: scaled sqrt(x²+y²) |
| ldexp(.,5) | 0.2 ulp | 0.2 ulp | 0.2 ulp | **9.3×** | **24×** | **26×** | original: ldexp on both limbs |

### Exponential / logarithmic

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| exp | 315 ulp | 1.3 ulp | 302 ulp | **5.6×** | **5.0×** | **4.0×** | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction |
| exp2 | — | — | — | **6.0×** | **6.0×** | **4.3×** | Julia: exp2 polynomial (14-term Horner) |
| expm1 | — | — | — | **5.7×** | **5.0×** | **4.0×** | original: exp(x) − 1 via DD sub |
| log | 283 ulp | 1.2 ulp | 283 ulp | **6.7×** | **4.6×** | **4.7×** | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) |
| log10 | 283 ulp | 1.2 ulp | 283 ulp | **8.7×** | **6.1×** | **6.3×** | Julia: log2 kernel × DD log10(2) |
| log2 | — | — | — | **8.1×** | **5.9×** | **6.0×** | Julia: log2 table lookup + polynomial |
| log1p | — | — | — | **6.7×** | **6.6×** | **5.3×** | original: log(1 + x) via DD add |
| pow | 292 ulp | 26 ulp | 289 ulp | **6.6×** | **6.2×** | **5.5×** | Julia: exp(y × log(x)) |

### Trigonometric

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| sin | 1.6 ulp | 1.7 ulp | 1.5 ulp | **3.8×** | **3.1×** | **2.3×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| cos | 1.8 ulp | 2.0 ulp | 1.8 ulp | **3.9×** | **3.2×** | **2.3×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| tan | 2.3 ulp | 2.3 ulp | 2.0 ulp | **2.7×** | **2.2×** | 1.2× | original: sin/cos Taylor kernels + DD divide |
| asin | 0.8 ulp | 0.8 ulp | 0.8 ulp | **6.0×** | **5.8×** | **6.4×** | original: piecewise rational P/Q (3 regions, from libquadmath asinq.c) |
| acos | 0.6 ulp | 0.6 ulp | 0.6 ulp | **6.0×** | **6.3×** | **6.6×** | original: asin polynomial + half-angle identity |
| atan | 0.8 ulp | 0.8 ulp | 0.8 ulp | **4.1×** | **4.4×** | **4.5×** | original: 84-entry table lookup + rational P(t²)/Q(t²) (from libquadmath atanq.c) |
| atan2 | 1.1 ulp | 1.1 ulp | 1.1 ulp | **3.5×** | **4.1×** | **4.2×** | original: table-based atan + quadrant correction |

### Hyperbolic

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| sinh | 326 ulp | 3.0 ulp | 321 ulp | **4.5×** | **2.7×** | **2.3×** | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 |
| cosh | 315 ulp | 1.3 ulp | 302 ulp | **3.0×** | **2.4×** | 2.0× | original: (exp+exp⁻¹)/2 |
| tanh | 47 ulp | 2.7 ulp | 47 ulp | **5.0×** | **3.1×** | **2.6×** | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) |
| asinh | 714 ulp | 2.7 ulp | 714 ulp | **7.1×** | **8.8×** | **9.1×** | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton |
| acosh | 1.3 ulp | 1.3 ulp | 1.3 ulp | **7.2×** | **8.4×** | **8.0×** | original: log(x+√(x²−1)) with Newton correction |
| atanh | 2520 ulp | 3.1 ulp | 2520 ulp | **7.2×** | **6.6×** | **7.1×** | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) |

### Error / special functions

| op | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × | approach |
|---|---|---|---|---|---|---|---|
| erf | 13 ulp | 0.8 ulp | 14 ulp | **5.0×** | **4.7×** | **4.4×** | piecewise rational approx (ported from libquadmath erfq.c) |
| erfc | 208 ulp | 1.6 ulp | 208 ulp | **5.3×** | **4.6×** | **4.2×** | piecewise rational approx + split exp(-x^2) |
| tgamma | 4596 ulp | 9.9 ulp | 4564 ulp | **10×** | **12×** | **12×** | piecewise rational approx + Stirling + reflection, exp(lgamma) |
| lgamma | 165 ulp | 7.3 ulp | 164 ulp | **5.0×** | **4.6×** | **4.8×** | piecewise rational approx + Stirling asymptotic |

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

- **tgamma / lgamma** (both C++ and Fortran) use a native double-double
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
  pointers (not in FP registers), adding ~1.5× overhead vs the C++ header-
  only version even with LTO inlining. See the performance note at the top
  of `fsrc/multifloats.fypp` for details and the `bind(c)` escape hatch.

- **Array reductions** (`dot_product`, `matmul`) use a fused multiply-
  accumulate kernel that computes the product's error-free representation
  and accumulates corrections into a scalar `s_lo`. `matmul` routes to a
  C kernel (`src/multifloats_math.cc`) in AXPY / gaxpy loop order; the
  kernel register-blocks any `m` via a strided panel template at `MR=8`
  plus a 1..7-row tail handler.

<!-- Auto-generated by bench/build_benchmark_md.py from per-system JSON
     results in bench/results/. Do not edit by hand. -->
