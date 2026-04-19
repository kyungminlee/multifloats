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
(`real(16)`) reference over ~1M random inputs (fixed seed 42):

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
| **sample** | Ad-hoc implementation in the benchmark harness, not exported from the library |

## Fortran: `float64x2` vs `real(16)`

Each operation is timed over 1024 elements × 400 repetitions (fast ops)
or fewer reps (transcendentals), with a NOINLINE drain after each rep to
prevent dead-code elimination. **×** = speedup (`qp_time / dd_time`,
values > 1× mean multifloats is faster); **err** = max\_rel from the
1M-input fuzz run; **prec** = precision label.

### Arithmetic

| op | approach | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|---|
| add | Julia: two\_sum EFT | full DD | 1.5e-32 | 1.5e-32 | 1.5e-32 | **2.1×** | **3.9×** | **2.9×** |
| sub | Julia: two\_sum EFT (negate + add) | full DD | 6.2e-33 | 6.2e-33 | 6.2e-33 | **2.1×** | **3.4×** | **2.9×** |
| mul | Julia: two\_prod EFT via FMA | full DD | 3.3e-32 | 3.3e-32 | 3.3e-32 | **5.2×** | **4.0×** | **4.2×** |
| div | original: Newton refinement (1/y seed, one step) | full DD | 5.4e-32 | 5.4e-32 | 5.4e-32 | **3.1×** | 1.8× | 1.8× |
| sqrt | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) | full DD | 1.7e-32 | 1.7e-32 | 5.2e-32 | **30×** | **38×** | **39×** |
| add (dd+dp) | Julia: two\_sum EFT | exact | exact | exact | exact | **2.3×** | **4.1×** | **3.8×** |
| mul (dp\*dd) | Julia: two\_prod EFT via FMA | full DD | 3.3e-32 | 3.3e-32 | 3.3e-32 | **5.0×** | **4.0×** | **5.3×** |

### Unary

| op | approach | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|---|
| abs | original: sign-check + negate limbs | exact | exact | exact | exact | 1.1× | 1.3× | 1.3× |
| neg | original: negate both limbs | exact | exact | exact | exact | 1.4× | **2.3×** | **2.0×** |
| aint | original: truncate hi, check DD fractional part | exact | exact | exact | exact | 0.83× | 1.2× | **2.1×** |
| anint | original: truncate hi, DD fractional part vs ±0.5 | exact | exact | exact | exact | 1.5× | 1.2× | 1.8× |
| fraction | original: scale both limbs by −exponent | exact | exact | exact | exact | 1.1× | 1.4× | 1.4× |
| scale | original: ldexp on both limbs | exact | exact | exact | exact | 0.88× | **5.0×** | **3.5×** |
| set\_exponent | original: scale + set\_exponent on hi | exact | exact | exact | exact | 1.5× | **3.5×** | **3.6×** |

### Binary

| op | approach | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|---|
| min | original: DD comparison + select | full DD | 6.2e-33 | 6.2e-33 | 6.2e-33 | 0.98× | **2.4×** | **2.9×** |
| max | original: DD comparison + select | full DD | 6.1e-33 | 6.1e-33 | 6.1e-33 | 0.95× | **3.0×** | **3.6×** |
| min3 | original: chained min | full DD | 5.5e-33 | 5.5e-33 | 5.5e-33 | **2.2×** | **3.6×** | **3.8×** |
| max3 | original: chained max | full DD | 5.8e-33 | 5.8e-33 | 5.8e-33 | **2.6×** | **3.6×** | **4.0×** |
| sign | original: sign-check + negate | exact | exact | exact | exact | 1.1× | 1.4× | 1.9× |
| dim | original: DD comparison, then subtract or zero | full DD | 6.2e-33 | 6.2e-33 | 6.2e-33 | 1.8× | **4.0×** | **4.5×** |
| hypot | original: scaled sqrt(x²+y²) | full DD | 6.2e-32 | 6.2e-32 | 7.9e-32 | **5.8×** | **9.7×** | **7.2×** |
| mod | sample: floor-multiple reduction loop; fallback to div chain | full DD | 2.0e-32 | 2.0e-32 | 2.0e-32 | 0.51× | 1.5× | 1.4× |
| modulo | original: mod + sign adjustment | full DD | 2.0e-32 | 2.0e-32 | 2.0e-32 | 1.1× | **2.1×** | **2.1×** |

### Exponential / logarithmic

| op | approach | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|---|
| exp | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction | full DD | 9.0e-30 | 9.0e-30 | 9.0e-30 | **5.4×** | **5.7×** | **3.7×** |
| log | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) | full DD | 4.4e-32 | 4.4e-32 | 4.4e-32 | **6.4×** | **5.0×** | **5.2×** |
| log10 | Julia: log2 kernel × DD log10(2) | full DD | 2.9e-32 | 2.9e-32 | 2.9e-32 | **8.8×** | **6.3×** | **6.6×** |
| pow | Julia: exp(y × log(x)) | full DD | 7.4e-30 | 7.4e-30 | 7.4e-30 | **6.2×** | **6.5×** | **4.9×** |
| pow\_int | original: repeated squaring via DD mul | full DD | 2.4e-32 | 2.4e-32 | 2.4e-32 | **7.4×** | **4.5×** | **5.2×** |

### Trigonometric

| op | approach | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|---|
| sin | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split | full DD | 5.0e-32 | 5.0e-32 | 3.6e-32 | **3.7×** | **3.3×** | **2.1×** |
| cos | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split | full DD | 3.8e-32 | 4.7e-32 | 4.1e-32 | **3.7×** | **3.4×** | **2.1×** |
| sinpi | Julia: sinpi Horner polynomial, direct | full DD | — | — | — | **4.6×** | **4.1×** | **3.9×** |
| cospi | Julia: cospi Horner polynomial, direct | full DD | — | — | — | **4.8×** | **4.2×** | **4.1×** |
| tan | original: sin/cos Taylor kernels + DD divide | full DD | 5.8e-32 | 5.8e-32 | 6.1e-32 | **2.7×** | **2.4×** | 1.1× |
| asin | original: piecewise rational P/Q (3 regions, from libquadmath asinq.c) | full DD | 6.4e-33 | 6.4e-33 | 1.7e-32 | **5.9×** | **5.8×** | **4.1×** |
| acos | original: asin polynomial + half-angle identity | full DD | 1.3e-32 | 1.3e-32 | 1.9e-32 | **6.0×** | **6.1×** | **3.6×** |
| atan | original: 84-entry table lookup + rational P(t²)/Q(t²) (from libquadmath atanq.c) | full DD | 2.6e-32 | 2.6e-32 | 2.6e-32 | **4.0×** | **4.0×** | **3.5×** |
| atan2 | original: table-based atan + quadrant correction | full DD | 1.9e-32 | 1.8e-32 | 2.9e-32 | **3.4×** | **3.8×** | **3.0×** |

### Hyperbolic

| op | approach | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|---|
| sinh | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 | full DD | 9.0e-30 | 9.0e-30 | 9.0e-30 | **4.5×** | **3.4×** | **2.3×** |
| cosh | original: (exp+exp⁻¹)/2 | full DD | 9.0e-30 | 9.0e-30 | 9.0e-30 | **3.0×** | **3.0×** | 2.0× |
| tanh | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) | full DD | 4.0e-30 | 4.0e-30 | 4.0e-30 | **5.1×** | **4.0×** | **2.9×** |
| asinh | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton | full DD | 2.1e-29 | 2.1e-29 | 2.1e-29 | **7.5×** | **11×** | **7.9×** |
| acosh | original: log(x+√(x²−1)) with Newton correction | full DD | 3.2e-32 | 3.5e-32 | 3.5e-32 | **6.6×** | **6.8×** | **6.5×** |
| atanh | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) | full DD | 5.9e-29 | 5.9e-29 | 5.9e-29 | **7.0×** | **7.9×** | **5.4×** |

### Error / special functions

| op | approach | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|---|
| erf | piecewise rational approx (libquadmath erfq.c) | full DD | 1.5e-32 | 1.5e-32 | 1.5e-32 | **5.0×** | **4.9×** | **3.3×** |
| erfc | piecewise rational approx + split exp(-x^2) | full DD | 6.6e-30 | 6.6e-30 | 6.6e-30 | **5.1×** | **5.0×** | **3.2×** |
| erfc\_scaled | exp(x^2)·erfc(x) with asymptotic cancellation | full DD | 7.7e-30 | 7.7e-30 | 7.7e-30 | **6.5×** | **7.6×** | **5.0×** |
| gamma | piecewise rational approx + Stirling + reflection | full DD | 3.3e-27 | 3.3e-27 | 3.3e-27 | **9.0×** | **8.4×** | **5.7×** |
| log\_gamma | piecewise rational approx + Stirling asymptotic | full DD | 1.5e-28 | 1.5e-28 | 1.5e-28 | **5.0×** | **4.6×** | **3.5×** |
| bessel\_j0 | piecewise rational + Hankel asymptotic (j0q.c) via C++ | full DD | 2.1e-29 | 2.1e-29 | 1.4e-29 | **6.5×** | **6.1×** | **4.9×** |
| bessel\_j1 | piecewise rational + Hankel asymptotic (j1q.c) via C++ | full DD | 4.3e-29 | 4.3e-29 | 2.3e-29 | **6.4×** | **5.8×** | **4.7×** |
| bessel\_jn(3,.) | forward/backward recurrence from j0/j1 | full DD | 5.0e-29 | 5.0e-29 | 5.0e-29 | **4.1×** | **4.5×** | **3.8×** |
| bessel\_y0 | piecewise rational + Hankel asymptotic (j0q.c) via C++ | full DD | 1.4e-29 | 1.4e-29 | 3.6e-30 | **7.0×** | **6.6×** | **4.9×** |
| bessel\_y1 | piecewise rational + Hankel asymptotic (j1q.c) via C++ | full DD | 1.1e-28 | 1.1e-28 | 2.8e-29 | **7.3×** | **6.5×** | **5.0×** |
| bessel\_yn(3,.) | forward recurrence from y0/y1 | full DD | 2.6e-28 | 2.6e-28 | 1.4e-28 | **7.3×** | **6.3×** | **5.0×** |

### Complex arithmetic

| op | approach | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|---|
| cdd\_add | original: component-wise DD add | full DD | 1.2e-32 (re) / 1.2e-32 (im) | 1.2e-32 (re) / 1.2e-32 (im) | 1.2e-32 (re) / 1.2e-32 (im) | **3.2×** | **3.1×** | **5.4×** |
| cdd\_sub | original: component-wise DD sub | full DD | 5.6e-33 (re) / 5.6e-33 (im) | 5.6e-33 (re) / 5.6e-33 (im) | 5.6e-33 (re) / 5.6e-33 (im) | **3.2×** | **5.2×** | **5.2×** |
| cdd\_mul | original: (ac−bd, ad+bc) via DD ops | full DD | exact (re) / 1.9e-32 (im) | exact (re) / 1.9e-32 (im) | exact (re) / 1.9e-32 (im) | **5.3×** | **3.5×** | **5.3×** |
| cdd\_div | original: (ac+bd, bc−ad)/(c²+d²) | full DD / deriv | 4.5e-32 (re) / 1.5e-16 (im) | 4.5e-32 (re) / 1.5e-16 (im) | 4.5e-32 (re) / 1.5e-16 (im) | **5.0×** | **2.7×** | **3.1×** |
| cdd\_conjg | original: negate im limbs | exact | exact | exact | exact | **2.0×** | **2.9×** | **3.1×** |
| cdd\_abs | original: hypot(re, im) | full DD | 3.7e-32 | 3.7e-32 | 7.9e-32 | **5.6×** | **10.0×** | **7.3×** |

### Complex transcendentals

| op | approach | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|---|
| cdd\_sqrt | original: Kahan-style (\|z\|+\|a\|)/2 with scaling | full DD | 5.2e-32 (re) / 4.4e-32 (im) | 5.2e-32 (re) / 4.4e-32 (im) | 6.5e-32 (re) / 5.9e-32 (im) | **7.2×** | **11×** | **5.8×** |
| cdd\_exp | original: exp(re)·(cos(im), sin(im)) | full DD | 7.0e-30 (re) / 2.0e-29 (im) | 7.0e-30 (re) / 2.0e-29 (im) | 7.1e-30 (re) / 6.9e-30 (im) | **4.4×** | **4.4×** | **2.1×** |
| cdd\_log | original: (log(\|z\|), atan2(im,re)) | full DD | 1.4e-30 (re) / 1.8e-32 (im) | 1.4e-30 (re) / 1.8e-32 (im) | 1.4e-30 (re) / 3.1e-32 (im) | **5.0×** | **5.3×** | **3.8×** |
| cdd\_sin | original: sin(re)cosh(im), cos(re)sinh(im) | full DD | 7.7e-30 (re) / 7.7e-30 (im) | 7.7e-30 (re) / 7.7e-30 (im) | 7.7e-30 (re) / 7.7e-30 (im) | **5.6×** | **4.9×** | 1.9× |
| cdd\_cos | original: cos(re)cosh(im), −sin(re)sinh(im) | full DD | 7.7e-30 (re) / 7.7e-30 (im) | 7.7e-30 (re) / 7.7e-30 (im) | 7.7e-30 (re) / 7.7e-30 (im) | **5.9×** | **5.0×** | 1.7× |
| cdd\_tan | original: complex sin/cos ratio | full DD | 6.6e-30 (re) / 1.3e-29 (im) | 6.6e-30 (re) / 1.3e-29 (im) | 1.6e-30 (re) / 6.7e-31 (im) | **5.5×** | **5.0×** | 0.97× |
| cdd\_sinh | original: sinh(re)cos(im), cosh(re)sin(im) | full DD | 7.4e-30 (re) / 2.0e-29 (im) | 7.4e-30 (re) / 2.0e-29 (im) | 7.4e-30 (re) / 6.9e-30 (im) | **6.0×** | **5.0×** | 2.0× |
| cdd\_cosh | original: cosh(re)cos(im), sinh(re)sin(im) | full DD | 7.0e-30 (re) / 2.0e-29 (im) | 7.0e-30 (re) / 2.0e-29 (im) | 7.1e-30 (re) / 7.3e-30 (im) | **6.3×** | **5.0×** | 1.9× |
| cdd\_tanh | original: complex tanh via sinh/cosh | full DD | 1.2e-29 (re) / 6.4e-30 (im) | 1.2e-29 (re) / 6.4e-30 (im) | 6.7e-31 (re) / 8.6e-31 (im) | **6.0×** | **5.2×** | 0.96× |
| cdd\_asin | original: −i·log(iz+√(1−z²)) | deriv / full DD | 3.5e-32 (re) / 2.1e-28 (im) | 3.3e-32 (re) / 2.1e-28 (im) | 2.5e-23 (re) / 2.1e-28 (im) | **5.0×** | **6.2×** | **4.4×** |
| cdd\_acos | original: π/2 − asin(z) | full DD | 1.4e-32 (re) / 2.1e-28 (im) | 1.4e-32 (re) / 2.1e-28 (im) | 1.3e-32 (re) / 2.1e-28 (im) | **4.7×** | **6.1×** | **4.1×** |
| cdd\_atan | original: (i/2)·log((i+z)/(i−z)) | full DD | 3.5e-32 (re) / 1.2e-29 (im) | 3.5e-32 (re) / 1.2e-29 (im) | 3.6e-32 (re) / 3.8e-31 (im) | **5.3×** | **5.7×** | **2.6×** |
| cdd\_asinh | original: log(z+√(z²+1)) | deriv / full DD | 2.0e-21 (re) / 3.6e-32 (im) | 2.0e-21 (re) / 3.6e-32 (im) | 2.3e-21 (re) / 7.3e-32 (im) | **4.7×** | **6.0×** | **3.7×** |
| cdd\_acosh | original: log(z+√(z²−1)) | full DD | 2.7e-29 (re) / 3.1e-32 (im) | 2.7e-29 (re) / 3.1e-32 (im) | 6.2e-30 (re) / 1.9e-32 (im) | **3.9×** | **5.1×** | **3.1×** |
| cdd\_atanh | original: ½·log((1+z)/(1−z)) | deriv / full DD | 6.0e-23 (re) / 2.8e-32 (im) | 6.0e-23 (re) / 2.8e-32 (im) | 4.3e-22 (re) / 5.6e-32 (im) | **5.8×** | **6.4×** | **3.1×** |

### Array reductions

| op | approach | prec | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|---|
| arr\_sum (n=8) | original: chained DD add | full DD | 6.6e-30 | 6.6e-30 | 6.6e-30 | **4.0×** | **5.7×** | **5.6×** |
| arr\_product (n=8) | original: chained DD mul | full DD | 4.0e-50 | 4.0e-50 | 4.0e-50 | **2.0×** | 2.0× | 2.0× |
| arr\_maxval (n=8) | original: chained DD compare | full DD | 5.9e-33 | 5.9e-33 | 5.9e-33 | **3.1×** | **4.4×** | **6.7×** |
| arr\_minval (n=8) | original: chained DD compare | full DD | 6.0e-33 | 6.0e-33 | 6.0e-33 | **3.4×** | **5.7×** | **6.7×** |
| arr\_dot (n=8) | original: fused multiply-accumulate with periodic renormalization | full DD | 2.1e-31 | 2.1e-31 | 2.1e-31 | **4.1×** | **3.8×** | **4.4×** |
| arr\_norm2 (n=8) | original: sqrt(dot(x,x)) | full DD | 3.2e-32 | 3.2e-32 | 5.2e-32 | **5.4×** | **7.0×** | **6.1×** |
| arr\_matmul (8×8·8) | original: AXPY-order C kernel, MR=8 register-blocked panels + 1..7 tail, periodic renorm | full DD | 7.1e-30 | 7.1e-30 | 7.1e-30 | **2.7×** | 0.54× | 0.56× |

## C++: `MultiFloat<double,2>` vs `__float128`

Header-only — all kernels inline into the call site. No LTO needed.
See the Fortran tables for `prec` labels.

**Methodology.** Precision and timing are measured separately, mirroring
the Fortran split (`fortran_fuzz` / `fortran_bench`):

- **err** columns come from `cpp_fuzz` at 1M iterations, fixed seed 42.
- **×** columns come from `cpp_bench` at 1024 elements × reps.

### Arithmetic

| op | approach | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|
| add | Julia: two\_sum EFT | 4.7e-30 | 4.7e-30 | 4.7e-30 | **3.1×** | **4.0×** | **4.9×** |
| sub | Julia: two\_sum EFT (negate + add) | 6.1e-30 | 6.1e-30 | 6.1e-30 | **3.2×** | **4.7×** | **4.9×** |
| mul | Julia: two\_prod EFT via FMA | 5.5e-32 | 5.5e-32 | 5.5e-32 | **14×** | **7.2×** | **8.1×** |
| div | original: Newton refinement (1/y seed, one step) | 7.7e-32 | 7.7e-32 | 7.7e-32 | **8.6×** | **5.2×** | **5.6×** |
| sqrt | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) | 4.4e-32 | 4.4e-32 | 4.4e-32 | **48×** | **53×** | **54×** |
| cbrt | original: Newton correction on cbrt(hi) seed | — | — | — | **39×** | **17×** | **17×** |
| fma | original: x\*y + z via DD ops | — | — | — | **88×** | **167×** | **161×** |
| abs | original: sign-check + negate limbs | 6.2e-33 | 6.2e-33 | 6.2e-33 | 1.6× | **7.7×** | **10×** |
| neg | original: negate both limbs | 6.2e-33 | 6.2e-33 | 6.2e-33 | **4.3×** | **6.4×** | **9.2×** |

### Rounding

| op | approach | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|
| floor | original: floor hi, adjust lo | exact | exact | exact | **3.0×** | **5.6×** | **5.5×** |
| ceil | original: ceil hi, adjust lo | exact | exact | exact | **3.1×** | **5.7×** | **5.7×** |
| trunc | original: signbit ? −floor(−x) : floor(x) | exact | exact | exact | **2.5×** | **4.8×** | **4.8×** |
| round | original: trunc(x + ½·sign(x)) | exact | exact | exact | **2.1×** | **2.1×** | **2.4×** |
| rint | original: nearbyint on hi, adjust lo | exact | exact | exact | **10×** | **18×** | **9.6×** |
| nearbyint | original: nearbyint on hi, adjust lo | exact | exact | exact | **33×** | **139×** | **92×** |

### Binary

| op | approach | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|
| fmin | original: DD comparison + select | 6.2e-33 | 6.2e-33 | 6.2e-33 | **7.5×** | **9.3×** | **11×** |
| fmax | original: DD comparison + select | 6.2e-33 | 6.2e-33 | 6.2e-33 | **7.1×** | **11×** | **11×** |
| fdim | original: DD comparison, then subtract or zero | 3.9e-30 | 3.9e-30 | 3.9e-30 | **6.8×** | **8.7×** | **8.2×** |
| copysign | original: sign-bit copy to hi, propagate to lo | 6.2e-33 | 6.2e-33 | 6.2e-33 | **3.3×** | **3.9×** | **9.4×** |
| fmod | sample: floor-multiple reduction loop; fallback to div chain | 1.3e-16 | 1.3e-16 | 1.3e-16 | 1.1× | 0.97× | 1.2× |
| hypot | original: scaled sqrt(x²+y²) | 5.1e-32 | 5.1e-32 | 3.8e-32 | **6.4×** | **8.0×** | **38×** |
| ldexp(.,5) | original: ldexp on both limbs | 6.1e-33 | 6.1e-33 | 6.1e-33 | **9.3×** | **24×** | **26×** |

### Exponential / logarithmic

| op | approach | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|
| exp | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction | 7.8e-30 | 7.5e-30 | 7.5e-30 | **5.6×** | **6.1×** | **4.0×** |
| exp2 | Julia: exp2 polynomial (14-term Horner) | — | — | — | **6.0×** | **6.8×** | **4.3×** |
| expm1 | original: exp(x) − 1 via DD sub | — | — | — | **5.7×** | **4.6×** | **4.0×** |
| log | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) | 7.0e-30 | 7.0e-30 | 7.0e-30 | **6.7×** | **4.6×** | **4.7×** |
| log10 | Julia: log2 kernel × DD log10(2) | 7.0e-30 | 7.0e-30 | 7.0e-30 | **8.7×** | **6.1×** | **6.3×** |
| log2 | Julia: log2 table lookup + polynomial | — | — | — | **8.1×** | **5.9×** | **6.0×** |
| log1p | original: log(1 + x) via DD add | — | — | — | **6.7×** | **6.2×** | **5.3×** |
| pow | Julia: exp(y × log(x)) | 7.2e-30 | 7.2e-30 | 7.1e-30 | **6.6×** | **6.5×** | **5.5×** |

### Trigonometric

| op | approach | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|
| sin | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split | 4.0e-32 | 4.3e-32 | 3.6e-32 | **3.8×** | **3.4×** | **2.3×** |
| cos | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split | 4.4e-32 | 4.9e-32 | 4.4e-32 | **3.9×** | **3.5×** | **2.3×** |
| tan | original: sin/cos Taylor kernels + DD divide | 5.6e-32 | 5.6e-32 | 5.0e-32 | **2.7×** | **2.5×** | 1.2× |
| asin | original: piecewise rational P/Q (3 regions, from libquadmath asinq.c) | 1.9e-32 | 1.9e-32 | 1.9e-32 | **6.0×** | **6.2×** | **6.4×** |
| acos | original: asin polynomial + half-angle identity | 1.4e-32 | 1.4e-32 | 1.4e-32 | **6.0×** | **6.3×** | **6.6×** |
| atan | original: 84-entry table lookup + rational P(t²)/Q(t²) (from libquadmath atanq.c) | 1.9e-32 | 1.9e-32 | 1.9e-32 | **4.1×** | **4.3×** | **4.5×** |
| atan2 | original: table-based atan + quadrant correction | 2.8e-32 | 2.6e-32 | 2.6e-32 | **3.5×** | **4.0×** | **4.2×** |

### Hyperbolic

| op | approach | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|
| sinh | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 | 8.0e-30 | 8.0e-30 | 7.9e-30 | **4.5×** | **3.6×** | **2.3×** |
| cosh | original: (exp+exp⁻¹)/2 | 7.8e-30 | 7.5e-30 | 7.4e-30 | **3.0×** | **3.0×** | 2.0× |
| tanh | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) | 1.2e-30 | 1.2e-30 | 1.2e-30 | **5.0×** | **3.6×** | **2.6×** |
| asinh | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton | 1.8e-29 | 1.8e-29 | 1.8e-29 | **7.1×** | **9.3×** | **9.1×** |
| acosh | original: log(x+√(x²−1)) with Newton correction | 3.2e-32 | 3.2e-32 | 3.2e-32 | **7.2×** | **8.6×** | **8.0×** |
| atanh | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) | 6.2e-29 | 6.2e-29 | 6.2e-29 | **7.2×** | **7.1×** | **7.1×** |

### Error / special functions

| op | approach | M1 Max err | pop-os err | sandbox err | M1 Max × | pop-os × | sandbox × |
|---|---|---|---|---|---|---|---|
| erf | piecewise rational approx (ported from libquadmath erfq.c) | 3.3e-31 | 3.3e-31 | 3.4e-31 | **5.0×** | **5.4×** | **4.4×** |
| erfc | piecewise rational approx + split exp(-x^2) | 5.1e-30 | 5.1e-30 | 5.1e-30 | **5.3×** | **5.5×** | **4.2×** |
| tgamma | piecewise rational approx + Stirling + reflection, exp(lgamma) | 1.1e-28 | 1.1e-28 | 1.1e-28 | **10×** | **14×** | **12×** |
| lgamma | piecewise rational approx + Stirling asymptotic | 4.1e-30 | 4.1e-30 | 4.1e-30 | **5.0×** | **4.8×** | **4.8×** |

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
