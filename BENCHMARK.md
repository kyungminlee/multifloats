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
| **M1 Max** | Apple M1 Max (ARM64, 10 cores) | macOS 26.3 (Darwin 25.3.0) | GNU Fortran / g++ 15.2.0 (Homebrew GCC 15.2.0\_1) | CMake 4.3.1, `-O3 -flto`, STATIC library |
| **Skylake** | Intel Xeon family 6 model 85 (Skylake-SP / Cascade Lake), 2.8 GHz, 16 cores, AVX-512 (KVM), 22 GB | Ubuntu 24.04.4 LTS | GNU Fortran / g++ 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) | CMake 3.28.3, `-O3 -flto`, OBJECT library |
| **Raptor Lake** | Intel Core i3-1315U (Raptor Lake), 4.5 GHz boost, 6 cores / 8 threads, 16 GB | Pop!\_OS 24.04 LTS (Linux 6.17.9) | GNU Fortran / g++ 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) | CMake 3.28.3, `-O3 -flto`, OBJECT library |

The x86-64 speedup numbers are generally higher than M1 for many
operations because gfortran's software `libquadmath` on x86-64 is 2–5×
slower per operation than Apple's ARM64 quad path, inflating the
multifloats speedup relative to M1 even though the DD kernels themselves
run at similar speed on all three machines.

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
| **sample** | Adapted from `external/float64x2-sample.cpp` |

## Fortran: `float64x2` vs `real(16)`

Each operation is timed over 1024 elements × 400 repetitions (fast ops)
or fewer reps (transcendentals), with a NOINLINE drain after each rep to
prevent dead-code elimination. **×** = speedup (`qp_time / mf_time`,
values > 1× mean multifloats is faster); **err** = max\_rel from the
1M-input fuzz run; **prec** = precision label.

### Arithmetic

| op | approach | prec | M1 Max err | Skylake err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|---|---|
| add | Julia: two\_sum EFT | full DD | 1.5e-32 | 1.5e-32 | 1.5e-32 | 1.9× | **3.6×** | **5.6×** |
| sub | Julia: two\_sum EFT (negate + add) | full DD | 6.2e-33 | 6.2e-33 | 6.2e-33 | 1.3× | **3.1×** | **2.7×** |
| mul | Julia: two\_prod EFT via FMA | full DD | 3.3e-32 | 3.1e-32 | 3.3e-32 | **4.8×** | **5.9×** | **4.5×** |
| div | original: Newton refinement (1/y seed, one step) | full DD | 5.4e-32 | 6.1e-32 | 5.4e-32 | **3.3×** | **3.2×** | 1.7× |
| sqrt | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) | full DD | 5.2e-32 | 5.4e-32 | 5.2e-32 | **15×** | **17×** | **33×** |
| add (mf+dp) | Julia: two\_sum EFT | exact | exact | exact | exact | **2.0×** | **4.3×** | **4.0×** |
| mul (dp\*mf) | Julia: two\_prod EFT via FMA | full DD | 3.3e-32 | 3.1e-32 | 3.3e-32 | **4.9×** | **6.1×** | **3.9×** |

### Unary

| op | approach | prec | M1 Max err | Skylake err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|---|---|
| abs | original: sign-check + negate limbs | exact | exact | exact | exact | 1.1× | **2.4×** | 1.6× |
| neg | original: negate both limbs | exact | exact | exact | exact | 1.4× | **4.4×** | **2.4×** |
| aint | original: truncate hi, check DD fractional part | exact | exact | exact | exact | 1.7× | 1.4× | **2.5×** |
| anint | original: truncate hi, DD fractional part vs ±0.5 | exact | exact | exact | exact | 1.5× | **2.5×** | **2.4×** |
| fraction | original: scale both limbs by −exponent | exact | exact | exact | exact | 1.2× | 1.3× | 1.5× |
| scale | original: ldexp on both limbs | exact | exact | exact | exact | 0.91× | **3.5×** | **5.1×** |
| set\_exponent | original: scale + set\_exponent on hi | exact | exact | exact | exact | 1.5× | **3.3×** | **3.8×** |

### Binary

| op | approach | prec | M1 Max err | Skylake err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|---|---|
| min | original: DD comparison + select | full DD | 6.2e-33 | 5.9e-33 | 6.2e-33 | 0.90× | **4.2×** | **2.8×** |
| max | original: DD comparison + select | full DD | 6.1e-33 | 6.0e-33 | 6.1e-33 | 1.0× | **4.3×** | 1.8× |
| min3 | original: chained min | full DD | 5.5e-33 | 5.8e-33 | 5.5e-33 | **2.3×** | **6.9×** | **3.9×** |
| max3 | original: chained max | full DD | 5.8e-33 | 5.7e-33 | 5.8e-33 | **2.3×** | **6.4×** | **3.7×** |
| sign | original: sign-check + negate | exact | exact | exact | exact | 1.0× | **2.3×** | 1.5× |
| dim | original: DD comparison, then subtract or zero | full DD | 6.2e-33 | 5.9e-33 | 6.2e-33 | 1.9× | **4.3×** | **3.8×** |
| hypot | original: scaled sqrt(x²+y²) | full DD | 7.9e-32 | 7.2e-32 | 7.9e-32 | **3.9×** | **5.4×** | **7.3×** |
| mod | sample: floor-multiple reduction loop; fallback to div chain | full DD | 2.0e-32 | 3.2e-32 | 2.0e-32 | 0.57× | 1.1× | 1.5× |
| modulo | original: mod + sign adjustment | full DD | 2.0e-32 | 3.2e-32 | 2.0e-32 | 1.2× | 1.6× | **2.2×** |

### Exponential / logarithmic

| op | approach | prec | M1 Max err | Skylake err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|---|---|
| exp | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction | full DD | 9.0e-30 | 2.8e-30 | 9.0e-30 | **2.7×** | **2.7×** | **2.7×** |
| log | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) | full DD | 4.4e-32 | 3.0e-32 | 4.4e-32 | **4.1×** | **3.8×** | **5.1×** |
| log10 | Julia: log2 kernel × DD log10(2) | full DD | 2.9e-32 | 3.4e-32 | 2.9e-32 | **5.5×** | **4.9×** | **6.7×** |
| pow | Julia: exp(y × log(x)) | full DD | 7.4e-30 | 2.2e-30 | 7.4e-30 | **4.0×** | **4.1×** | **3.8×** |
| pow\_int | original: repeated squaring via DD mul | full DD | 2.4e-32 | 2.2e-32 | 2.4e-32 | **7.4×** | **5.6×** | **4.6×** |

### Trigonometric

| op | approach | prec | M1 Max err | Skylake err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|---|---|
| sin | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split | full DD | 3.4e-32 | 3.8e-32 | 3.6e-32 | 1.8× | 1.5× | **2.2×** |
| cos | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split | full DD | 3.0e-32 | 4.2e-32 | 4.1e-32 | 1.8× | 1.5× | **2.2×** |
| sinpi | Julia: sinpi Horner polynomial, direct | full DD | 4.9e-27 | 4.9e-27 | 4.9e-27 | **3.0×** | **2.8×** | **3.9×** |
| cospi | Julia: cospi Horner polynomial, direct | full DD | 8.2e-27 | 8.2e-27 | 8.2e-27 | **3.2×** | **2.6×** | **3.9×** |
| tan | original: sin/cos Taylor kernels + DD divide | full DD | 4.6e-32 | 6.7e-32 | 6.1e-32 | 0.85× | 0.8× | 1.1× |
| asin | original: piecewise rational P/Q (3 regions, from libquadmath asinq.c) | full DD | 1.7e-32 | 4.3e-32 | 1.7e-32 | **3.6×** | 1.2× | **4.2×** |
| acos | original: asin polynomial + half-angle identity | full DD | 1.9e-32 | 5.4e-32 | 1.9e-32 | **2.9×** | 1.2× | **3.7×** |
| atan | original: 84-entry table lookup + rational P(t²)/Q(t²) (from libquadmath atanq.c) | full DD | 2.6e-32 | 4.2e-32 | 2.6e-32 | **3.1×** | 0.7× | **3.5×** |
| atan2 | original: table-based atan + quadrant correction | full DD | 2.7e-32 | 3.0e-32 | 2.9e-32 | **2.6×** | 0.8× | **2.9×** |

### Hyperbolic

| op | approach | prec | M1 Max err | Skylake err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|---|---|
| sinh | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 | full DD | 9.0e-30 | 3.7e-30 | 9.0e-30 | **2.2×** | 1.9× | 1.5× |
| cosh | original: (exp+exp⁻¹)/2 | full DD | 9.0e-30 | 3.7e-30 | 9.0e-30 | 1.5× | 1.6× | 1.4× |
| tanh | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) | full DD | 4.0e-30 | 1.2e-30 | 4.0e-30 | **2.5×** | **2.3×** | **2.0×** |
| asinh | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton | full DD | 2.1e-29 | 2.7e-30 | 2.1e-29 | **5.3×** | **5.8×** | **7.9×** |
| acosh | original: log(x+√(x²−1)) with Newton correction | full DD | 3.2e-32 | 4.1e-32 | 3.5e-32 | **4.8×** | **4.9×** | **6.2×** |
| atanh | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) | full DD | 5.9e-29 | 1.2e-30 | 5.9e-29 | **4.3×** | **4.3×** | **5.3×** |

### Error / special functions

| op | approach | prec | M1 Max err | Skylake err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|---|---|
| erf | piecewise rational approx (libquadmath erfq.c) | full DD | 1.5e-32 | — | 1.5e-32 | **2.3×** | **5.4×** | **2.6×** |
| erfc | piecewise rational approx + split exp(-x^2) | full DD | 6.6e-30 | — | 6.6e-30 | **2.4×** | **5.5×** | **2.6×** |
| erfc\_scaled | exp(x^2)·erfc(x) with asymptotic cancellation | full DD | 7.7e-30 | — | 7.7e-30 | **3.3×** | **148×** | **3.8×** |
| gamma | piecewise rational approx + Stirling + reflection | full DD | 3.3e-27 | — | 3.3e-27 | **4.6×** | **42×** | **4.8×** |
| log\_gamma | piecewise rational approx + Stirling asymptotic | full DD | 1.5e-28 | — | 1.5e-28 | **2.9×** | **40×** | **3.4×** |
| bessel\_j0 | piecewise rational + Hankel asymptotic (j0q.c) via C++ | full DD | 1.4e-29 | — | 1.4e-29 | **3.8×** | **68×** | **4.6×** |
| bessel\_j1 | piecewise rational + Hankel asymptotic (j1q.c) via C++ | full DD | 2.3e-29 | — | 2.3e-29 | **4.2×** | **69×** | **4.0×** |
| bessel\_jn(3,.) | forward/backward recurrence from j0/j1 | full DD | 5.0e-29 | — | 5.0e-29 | **3.0×** | **67×** | **4.1×** |
| bessel\_y0 | piecewise rational + Hankel asymptotic (j0q.c) via C++ | full DD | 3.6e-30 | — | 3.6e-30 | **4.2×** | **75×** | **4.6×** |
| bessel\_y1 | piecewise rational + Hankel asymptotic (j1q.c) via C++ | full DD | 2.8e-29 | — | 2.8e-29 | **4.4×** | **71×** | **4.5×** |
| bessel\_yn(3,.) | forward recurrence from y0/y1 | full DD | 1.4e-28 | — | 1.4e-28 | **4.4×** | **75×** | **4.5×** |

### Complex arithmetic

| op | approach | prec | M1 Max err | Skylake err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|---|---|
| cx\_add | original: component-wise DD add | full DD | 1.2e-32 | 1.3e-32 | 1.2e-32 | **3.1×** | **3.6×** | **3.9×** |
| cx\_sub | original: component-wise DD sub | full DD | 5.6e-33 | 5.8e-33 | 5.6e-33 | **3.3×** | **3.4×** | **3.8×** |
| cx\_mul | original: (ac−bd, ad+bc) via DD ops | full DD | 1.9e-32 | 2.0e-32 | 1.9e-32 | **5.3×** | **4.0×** | **4.0×** |
| cx\_div | original: (ac+bd, bc−ad)/(c²+d²) | full DD / deriv | 4.5e-32 (re) / 1.5e-16 (im) | 4.6e-32 (re) / 1.1e-16 (im) | 4.5e-32 (re) / 1.5e-16 (im) | **5.3×** | **4.1×** | **2.7×** |
| cx\_conjg | original: negate im limbs | exact | exact | exact | exact | **2.0×** | **3.5×** | **2.4×** |
| cx\_abs | original: hypot(re, im) | full DD | 7.9e-32 | 6.6e-32 | 7.9e-32 | **3.5×** | **4.9×** | **7.3×** |

### Complex transcendentals

| op | approach | prec | M1 Max err | Skylake err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|---|---|
| cx\_sqrt | original: Kahan-style (\|z\|+\|a\|)/2 with scaling | full DD | 7.5e-32 | 6.4e-32 | 6.5e-32 | **4.3×** | **4.6×** | **5.7×** |
| cx\_exp | original: exp(re)·(cos(im), sin(im)) | full DD | 7.1e-30 | 7.8e-31 | 7.1e-30 | 1.5× | 1.5× | 1.9× |
| cx\_log | original: (log(\|z\|), atan2(im,re)) | full DD | 1.4e-30 (re) / 2.7e-32 (im) | 2.0e-29 (re) / 2.8e-32 (im) | 1.4e-30 (re) / 3.1e-32 (im) | 1.7× | 1.8× | **3.8×** |
| cx\_sin | original: sin(re)cosh(im), cos(re)sinh(im) | full DD | 7.7e-30 | 4.9e-31 | 7.7e-30 | 1.4× | 1.4× | 1.4× |
| cx\_cos | original: cos(re)cosh(im), −sin(re)sinh(im) | full DD | 7.7e-30 | 4.9e-31 | 7.7e-30 | 1.5× | 1.4× | 1.4× |
| cx\_tan | original: complex sin/cos ratio | full DD | 1.7e-30 | 3.1e-30 | 1.6e-30 | 0.77× | 0.7× | 0.72× |
| cx\_sinh | original: sinh(re)cos(im), cosh(re)sin(im) | full DD | 7.4e-30 | 1.2e-30 | 7.4e-30 | 1.6× | 1.5× | 1.5× |
| cx\_cosh | original: cosh(re)cos(im), sinh(re)sin(im) | full DD | 7.1e-30 | 1.2e-30 | 7.1e-30 | 1.6× | 1.5× | 1.6× |
| cx\_tanh | original: complex tanh via sinh/cosh | full DD | 1.2e-30 | 1.1e-30 | 8.6e-31 | 0.79× | 0.8× | 0.79× |
| cx\_asin | original: −i·log(iz+√(1−z²)) | deriv / full DD | 2.5e-23 (re) / 2.1e-28 (im) | 2.8e-23 (re) / 8.4e-31 (im) | 2.5e-23 (re) / 2.1e-28 (im) | **2.1×** | **2.4×** | **4.1×** |
| cx\_acos | original: π/2 − asin(z) | full DD | 1.7e-32 (re) / 2.1e-28 (im) | 2.9e-32 (re) / 8.4e-31 (im) | 1.3e-32 (re) / 2.1e-28 (im) | **2.1×** | **2.4×** | **4.1×** |
| cx\_atan | original: (i/2)·log((i+z)/(i−z)) | full DD | 5.0e-32 (re) / 3.1e-31 (im) | 3.8e-32 (re) / 4.8e-31 (im) | 3.6e-32 (re) / 3.8e-31 (im) | 1.4× | 1.4× | **2.6×** |
| cx\_asinh | original: log(z+√(z²+1)) | deriv / full DD | 2.3e-21 (re) / 6.8e-32 (im) | 6.7e-22 (re) / 6.6e-32 (im) | 2.3e-21 (re) / 7.3e-32 (im) | **2.0×** | **2.2×** | **4.0×** |
| cx\_acosh | original: log(z+√(z²−1)) | full DD | 6.2e-30 (re) / 2.2e-32 (im) | 7.2e-31 (re) / 2.9e-32 (im) | 6.2e-30 (re) / 1.9e-32 (im) | 1.8× | **2.0×** | **3.3×** |
| cx\_atanh | original: ½·log((1+z)/(1−z)) | deriv / full DD | 1.2e-22 (re) / 4.1e-32 (im) | 7.2e-23 (re) / 6.7e-32 (im) | 4.3e-22 (re) / 5.6e-32 (im) | 1.6× | 1.5× | **3.0×** |

### Array reductions

| op | approach | prec | M1 Max err | Skylake err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|---|---|
| arr\_sum (n=8) | original: chained DD add | full DD | 3.0e-30 | 2.5e-31 | 3.0e-30 | 1.2× | 1.4× | **2.6×** |
| arr\_product (n=8) | original: chained DD mul | full DD | 4.0e-50 | 4.3e-49 | 4.0e-50 | **2.0×** | **2.2×** | 2.0× |
| arr\_maxval (n=8) | original: chained DD compare | full DD | 5.9e-33 | 6.1e-33 | 5.9e-33 | **3.3×** | **5.8×** | **5.5×** |
| arr\_minval (n=8) | original: chained DD compare | full DD | 6.0e-33 | 6.0e-33 | 6.0e-33 | **3.0×** | **5.2×** | **5.5×** |
| arr\_dot (n=8) | original: fused multiply-accumulate with periodic renormalization | full DD | 2.1e-31 | 1.0e-31 | 2.1e-31 | **3.6×** | **4.6×** | **4.1×** |
| arr\_norm2 (n=8) | original: sqrt(dot(x,x)) | full DD | 5.2e-32 | 6.4e-32 | 5.2e-32 | **5.0×** | **6.8×** | **6.5×** |
| arr\_matmul (8×8\*8) | original: fused multiply-accumulate with periodic renormalization | full DD | 7.1e-30 | 8.8e-30 | 7.1e-30 | 1.0× | 0.8× | 0.58× |

## C++: `MultiFloat<double,2>` vs `__float128`

Header-only — all kernels inline into the call site. No LTO needed.
See the Fortran tables for `prec` labels.

**Methodology.** Precision and timing are now measured separately, to
mirror the Fortran split (`fortran_fuzz` / `fortran_bench`):

- **err** columns come from `cpp_fuzz` at 1M iterations, fixed seed 42.
  Inputs are sampled with the same corner-case biasing as `fuzz.f90`
  (10% non-finite, 10% close, 10% sum-near-zero cancellation,
  10% near-huge, 10% near-tiny, 50% wide random 10^[-30..30]). Each
  reported value is the per-op max relative error vs a `__float128`
  reference, aggregated across ~500k to ~1M finite samples for fast
  ops and ~5k for transcendentals.
- **×** columns come from `cpp_bench` at 1024 elements × reps, with an
  init-before-each-leg reset so drain feedback cannot drift the inputs
  between the qp and mf legs.

**Precision parity with Fortran.** The C++ header now hits full DD on
exp / log / pow, sin / cos / tan, and the whole hyperbolic family.
Earlier measurements showed ~1e-15 to ~1e-24 for these ops. Three
independent fixes landed together: stale `exp2_coefs` / `log2_narrow` /
`log2_wide` polynomial tables were resynced with the Fortran reference;
`dd_tanh_full` was rewritten to use the cancellation-free
`(1 − em2)/(1 + em2)` form; and `dd_sin_full` / `dd_cos_full` /
`dd_tan_full` now run 3-part Cody–Waite π/2 reduction with direct
sin/cos Taylor kernels rather than `sinpi(x · inv_pi)`. A subsequent
rewrite of tgamma / lgamma around a native DD Stirling kernel
(`detail::dd_lgamma_full` / `detail::dd_tgamma_full`) eliminated the
libm-seed floor for the gamma family. The Fortran API now also uses
the same Stirling approach (`dd_lgamma_stirling_shift` + reflection),
achieving full DD precision for both gamma and log\_gamma.
erf / erfc now use piecewise rational approximation coefficients from
libquadmath, achieving full DD precision.

All M1 Max err values are now measured with the proper fuzz/bench split
(`cpp_fuzz 1000000`, fixed seed 42). The Raptor Lake columns reflect
independent measurements on that system.

### Arithmetic

| op | approach | M1 Max err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|
| add | Julia: two\_sum EFT | 4.7e-30 | 4.7e-30 | **3.1×** | **4.1×** | **7.6×** |
| sub | Julia: two\_sum EFT (negate + add) | 6.1e-30 | 6.1e-30 | **3.2×** | **4.2×** | **3.9×** |
| mul | Julia: two\_prod EFT via FMA | 5.5e-32 | 5.5e-32 | **14×** | **7.7×** | **8.0×** |
| div | original: Newton refinement (1/y seed, one step) | 5.9e-32 | 5.9e-32 | **2.2×** | **2.6×** | 1.0× |
| sqrt | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) | 4.4e-32 | 4.4e-32 | **49×** | **29×** | **56×** |
| cbrt | original: Newton correction on cbrt(hi) seed | (not measured) | (not measured) | **40×** | **13×** | **18×** |
| fma | original: x\*y + z via DD ops | (not measured) | (not measured) | **88×** | **69×** | **141×** |
| abs | original: sign-check + negate limbs | 6.2e-33 | 6.2e-33 | **3.8×** | **5.1×** | **7.7×** |
| neg | original: negate both limbs | 6.2e-33 | 6.2e-33 | **4.3×** | **3.6×** | **6.1×** |

### Rounding

| op | approach | M1 Max err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|
| floor | original: floor hi, adjust lo | exact | exact | **3.1×** | **4.2×** | **5.5×** |
| ceil | original: ceil hi, adjust lo | exact | exact | **3.0×** | **3.9×** | **4.9×** |
| trunc | original: signbit ? −floor(−x) : floor(x) | exact | exact | **2.3×** | **3.2×** | **4.8×** |
| round | original: trunc(x + ½·sign(x)) | exact | exact | 1.2× | 1.0× | **2.1×** |
| rint | original: nearbyint on hi, adjust lo | exact | exact | **11×** | **11×** | **19×** |
| nearbyint | original: nearbyint on hi, adjust lo | exact | exact | **34×** | **65×** | **151×** |

### Binary

| op | approach | M1 Max err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|
| fmin | original: DD comparison + select | 6.2e-33 | 6.2e-33 | **7.6×** | **9.8×** | **8.1×** |
| fmax | original: DD comparison + select | 6.2e-33 | 6.2e-33 | **7.0×** | **13×** | **11×** |
| fdim | original: DD comparison, then subtract or zero | 3.9e-30 | 3.9e-30 | **6.7×** | **6.9×** | **7.9×** |
| copysign | original: sign-bit copy to hi, propagate to lo | 6.2e-33 | 6.2e-33 | **3.4×** | **5.2×** | **8.1×** |
| fmod | sample: floor-multiple reduction loop; fallback to div chain | 1.3e-16 | 1.3e-16 | 1.1× | 0.76× | 1.4× |
| hypot | original: scaled sqrt(x²+y²) | 3.9e-32 | 3.9e-32 | **38×** | **20×** | **37×** |
| ldexp(.,5) | original: ldexp on both limbs | 6.1e-33 | 6.1e-33 | **2.7×** | **2.4×** | **2.9×** |

### Exponential / logarithmic

| op | approach | M1 Max err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|
| exp | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction | 7.8e-30 | 7.5e-30 | **3.0×** | **2.7×** | **4.0×** |
| exp2 | Julia: exp2 polynomial (14-term Horner) | (not measured) | (not measured) | **3.1×** | **3.1×** | **4.5×** |
| expm1 | original: exp(x) − 1 via DD sub | (not measured) | (not measured) | **4.3×** | **3.3×** | **3.8×** |
| log | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) | 7.0e-30 | 7.0e-30 | **5.1×** | **4.0×** | **4.3×** |
| log10 | Julia: log2 kernel × DD log10(2) | 7.0e-30 | 7.0e-30 | **6.9×** | **5.1×** | **5.5×** |
| log2 | Julia: log2 table lookup + polynomial | (not measured) | (not measured) | **6.6×** | **4.7×** | **6.1×** |
| log1p | original: log(1 + x) via DD add | (not measured) | (not measured) | **5.6×** | **4.3×** | **3.3×** |
| pow | Julia: exp(y × log(x)) | 7.1e-30 | 7.1e-30 | **4.6×** | **4.1×** | **5.3×** |

### Trigonometric

| op | approach | M1 Max err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|
| sin | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split | 3.6e-32 | 3.6e-32 | **2.0×** | **2.3×** | **2.0×** |
| cos | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split | 4.4e-32 | 4.4e-32 | **2.1×** | **2.3×** | **2.1×** |
| tan | original: sin/cos Taylor kernels + DD divide | 5.2e-32 | 5.2e-32 | 0.97× | 1.1× | 1.2× |
| asin | original: piecewise rational P/Q (3 regions, from libquadmath asinq.c) | 1.9e-32 | 1.9e-32 | **5.6×** | 1.7× | **5.4×** |
| acos | original: asin polynomial + half-angle identity | 1.4e-32 | 1.4e-32 | **5.3×** | 1.7× | **5.5×** |
| atan | original: 84-entry table lookup + rational P(t²)/Q(t²) (from libquadmath atanq.c) | 2.5e-32 | 2.5e-32 | **3.5×** | 1.0× | **3.0×** |
| atan2 | original: table-based atan + quadrant correction | 3.2e-32 | 3.2e-32 | **2.9×** | 1.1× | **2.5×** |

### Hyperbolic

| op | approach | M1 Max err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|
| sinh | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 | 7.9e-30 | 7.9e-30 | **2.4×** | 1.9× | **2.2×** |
| cosh | original: (exp+exp⁻¹)/2 | 7.8e-30 | 7.4e-30 | 1.7× | 1.5× | 1.9× |
| tanh | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) | 1.2e-30 | 1.2e-30 | **2.5×** | **2.2×** | **2.4×** |
| asinh | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton | 1.8e-29 | 1.8e-29 | **6.9×** | **6.0×** | **8.1×** |
| acosh | original: log(x+√(x²−1)) with Newton correction | 3.2e-32 | 3.2e-32 | **6.2×** | **5.5×** | **8.4×** |
| atanh | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) | 6.2e-29 | 6.2e-29 | **5.2×** | **4.4×** | **4.0×** |

### Error / special functions

| op | approach | M1 Max err | Raptor Lake err | M1 Max × | Skylake × | Raptor Lake × |
|---|---|---|---|---|---|---|
| erf | piecewise rational approx (ported from libquadmath erfq.c) | 3.4e-31 | 3.4e-31 | **3.2×** | **5.7×** | **4.4×** |
| erfc | piecewise rational approx + split exp(-x^2) | 5.1e-30 | 5.1e-30 | **3.2×** | **5.5×** | **4.0×** |
| tgamma | piecewise rational approx + Stirling + reflection, exp(lgamma) | 1.1e-28 | 1.1e-28 | **7.2×** | **59×** | **11×** |
| lgamma | piecewise rational approx + Stirling asymptotic | 4.1e-30 | 4.1e-30 | **4.3×** | **46×** | **4.3×** |

## Notes

- **`mod` / `fmod`** is the only operation where quad precision is
  consistently faster. The DD `mod` uses a floor-multiple reduction loop
  (adapted from `external/float64x2-sample.cpp`) for small quotients, or a
  full DD divide chain for large quotients; libquadmath's `fmodq` uses a
  specialized bit-level remainder algorithm. `modulo` (Fortran) now beats qp
  at 1.4–1.6× thanks to the iterative approach. Precision degrades as
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
  truncated at ~1e-23 at x² ≈ (π/8)², capping precision. Fix: extend to
  13 terms. (3) Inverse trig functions (asin, acos, atan) were rewritten
  from Newton-on-sin/cos to piecewise rational polynomials from
  libquadmath (table-based atan, 3-region asin), eliminating the
  expensive sin/cos calls and improving both speed and precision.
  All three fixes are in both the Fortran and C++ codebases
  and benefit both platforms.

- **tgamma / lgamma** (both C++ and Fortran) use a native double-double
  Stirling kernel. `lgamma` evaluates the 13-term Stirling asymptotic
  `(x−½)·log x − x + ½·log(2π) + Σ B_{2k}/(2k(2k−1)·x^{2k−1})` in DD,
  after shifting the argument up to x ≥ 25 via a product accumulator so a
  single `dd_log(prod)` absorbs the recurrence. Small arguments (x < 0.5)
  use the reflection `log Γ(x) = log π − log|sin(πx)| − log Γ(1−x)`.
  `tgamma` / `gamma` derives from `exp(lgamma)`, with `π/(sin(πx)·Γ(1−x))`
  for negative x. Both deliver full DD precision (max\_rel ~9e-30 / ~5e-28
  on the 1M fuzz).

- **Inverse trigonometric** functions (asin, acos, atan, atan2) use
  piecewise rational polynomial approximations ported from libquadmath's
  `asinq.c` / `atanq.c`. `atan` uses an 84-entry lookup table for
  argument reduction followed by a small rational `P(t²)/Q(t²)`. `asin`
  uses three regions (|x|<0.5, 0.5–0.625 centered at 0.5625, 0.625–1 via
  half-angle identity). `acos` derives from the `asin` polynomial with
  domain-appropriate formulas to avoid cancellation near ±1.

- **erf / erfc** use piecewise rational approximation coefficients from
  libquadmath's `erfq.c`, evaluated in full DD arithmetic via Estrin's
  scheme (`dd_neval` / `dd_deval`). The asymptotic region splits `exp(-x^2)`
  into `exp(-s^2 - 0.5625) * exp((s-x)(s+x) + R)` with `s` truncated to
  ~27 mantissa bits so `s^2` is exact. Both deliver full DD precision
  (max\_rel ~1.5e-32 / ~6.6e-30 on the 1M fuzz).

- The **Fortran** multifloats module uses `elemental` functions on a
  `sequence` derived type. gfortran's ABI passes/returns these via hidden
  pointers (not in FP registers), adding ~1.5× overhead vs the C++ header-
  only version even with LTO inlining. See the performance note at the top
  of `fsrc/multifloats.fypp` for details and the `bind(c)` escape hatch.

- **Array reductions** (dot\_product, matmul) use a fused multiply-
  accumulate kernel that computes the product's error-free representation
  and accumulates corrections into a scalar `s_lo`, with periodic
  renormalization (configurable via `mf_set_fma_renorm_interval`).
