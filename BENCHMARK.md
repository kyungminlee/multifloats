# Benchmark Results

Comparison of multifloats double-double (DD) arithmetic against quad
precision (`real(16)` / `__float128` via libquadmath). The multifloats
implementations use explicit two-limb error-free transformations (EFTs)
on pairs of `double` values, achieving ~106 bits of significand — close
to quad precision's 113 bits — at a fraction of the cost.

## System

| | |
|---|---|
| **CPU** | Apple M1 Max (ARM64, 10 cores) |
| **RAM** | 64 GB |
| **OS** | macOS 26.3 (Darwin 25.3.0) |
| **Fortran** | GNU Fortran (Homebrew GCC 15.2.0\_1) 15.2.0 |
| **C++** | g++-15 (Homebrew GCC 15.2.0\_1) 15.2.0 |
| **Build** | CMake 4.3.1, `-O3 -flto`, OBJECT library (not STATIC — macOS `ar` strips GIMPLE from Mach-O, breaking LTO through `.a` archives) |

## Precision key

Precision is measured as the maximum relative error vs the quad-precision
(`real(16)`) reference over ~1M random inputs. The "precision" column
indicates the kernel implementation strategy:

| Label | max\_rel | Meaning |
|---|---|---|
| **full DD** | ~1e-30 | Full double-double EFT kernel (~106 bits) |
| **exact** | 0.0 | Bit-exact (no rounding involved) |
| **deriv-corrected** | ~1e-25 to 1e-18 | `f(hi) + f'(hi)*lo` correction gives near-DD |
| **single-double** | ~1e-16 to 1e-14 | Leading-limb libm call, no lo correction |

## Origin key

Each kernel's origin and algorithm are noted in the "approach" column:

| Tag | Meaning |
|---|---|
| **Julia** | Ported from `external/MultiFloats.jl/src/` |
| **original** | Developed for this project |
| **sample** | Adapted from `external/float64x2-sample.cpp` |

## Fortran: `float64x2` vs `real(16)`

Each operation is timed over 1024 elements × 400 repetitions (fast ops)
or fewer reps (transcendentals), with a NOINLINE drain after each rep to
prevent dead-code elimination. The "speedup" column is `qp_time / mf_time`:
values > 1× mean multifloats is faster.

### Arithmetic

| op | speedup | max\_rel | precision | approach |
|---|---|---|---|---|
| add | **1.9×** | 1.5e-32 | full DD | Julia: two\_sum EFT |
| sub | **2.1×** | 6.2e-33 | full DD | Julia: two\_sum EFT (negate + add) |
| mul | **8.6×** | 2.9e-32 | full DD | Julia: two\_prod EFT via FMA |
| div | **3.1×** | 5.1e-32 | full DD | original: Newton refinement (1/y seed, one step) |
| sqrt | **39×** | 5.5e-32 | full DD | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) |
| add (mf+dp) | **2.6×** | exact | full DD | Julia: two\_sum EFT |
| mul (dp\*mf) | **8.9×** | 2.9e-32 | full DD | Julia: two\_prod EFT via FMA |

### Unary

| op | speedup | max\_rel | precision | approach |
|---|---|---|---|---|
| abs | 1.4× | exact | exact | original: sign-check + negate limbs |
| neg | **2.2×** | exact | exact | original: negate both limbs |
| aint | 1.6× | exact | exact | original: truncate hi, check DD fractional part |
| anint | **2.4×** | exact | exact | original: truncate hi, DD fractional part vs ±0.5 |
| fraction | 1.4× | exact | exact | original: scale both limbs by -exponent |
| scale | 1.0× | exact | exact | original: ldexp on both limbs |
| set\_exponent | 1.8× | exact | exact | original: scale + set\_exponent on hi |

### Binary

| op | speedup | max\_rel | precision | approach |
|---|---|---|---|---|
| min | **2.2×** | 6.0e-33 | full DD | original: DD comparison + select |
| max | **2.1×** | 6.1e-33 | full DD | original: DD comparison + select |
| min3 | **2.2×** | 5.7e-33 | full DD | original: chained min |
| max3 | **4.3×** | 6.1e-33 | full DD | original: chained max |
| sign | 1.2× | exact | exact | original: sign-check + negate |
| dim | **2.8×** | 6.0e-33 | full DD | original: DD comparison, then subtract or zero |
| hypot | **4.7×** | 8.2e-32 | full DD | original: scaled sqrt(x²+y²) |
| mod | 0.53× | 1.9e-32 | full DD | sample: floor-multiple reduction loop; fallback to div chain for large quotients |
| modulo | **1.3×** | 1.9e-32 | full DD | original: mod + sign adjustment |

### Exponential / logarithmic

| op | speedup | max\_rel | precision | approach |
|---|---|---|---|---|
| exp | **2.9×** | 2.5e-30 | full DD | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction |
| log | **5.0×** | 4.0e-32 | full DD | Julia: log2 with table lookup (32 centers) + polynomial (7-term Horner) |
| log10 | **6.6×** | 2.9e-32 | full DD | Julia: log2 kernel × DD log10(2) |
| pow | **4.3×** | 2.9e-30 | full DD | Julia: exp(y × log(x)) |
| pow\_int | **16×** | 1.7e-32 | full DD | original: repeated squaring via DD mul |

### Trigonometric

| op | speedup | max\_rel | precision | approach |
|---|---|---|---|---|
| sin | **2.7×** | 4.7e-25 | near-DD | Julia: sinpi/cospi Horner polynomial (full DD) + 3-part 1/π Cody–Waite reduction (precision limited by reduction, not polynomial) |
| cos | **2.7×** | 3.4e-25 | near-DD | Julia: sinpi/cospi Horner polynomial (full DD) + 3-part 1/π Cody–Waite reduction (precision limited by reduction, not polynomial) |
| sinpi | **3.5×** | 4.9e-27 | full DD | Julia: sinpi Horner polynomial, direct (no 1/π reduction needed) |
| cospi | **3.5×** | 8.2e-27 | full DD | Julia: cospi Horner polynomial, direct (no 1/π reduction needed) |
| tan | 1.3× | 4.7e-25 | near-DD | Julia: sinpi/cospi kernels + DD divide; precision limited by 1/π reduction |
| asin | **2.0×** | 3.4e-32 | full DD | original: Newton step on sin, seeded by libm asin(hi) |
| acos | 1.9× | 2.6e-32 | full DD | original: Newton step on cos, seeded by libm acos(hi) |
| atan | 1.2× | 9.8e-23 | deriv-corrected | original: Newton step on tan, seeded by libm atan(hi) |
| atan2 | 1.2× | 4.3e-32 | full DD | original: Newton step on atan, with quadrant correction |

### Hyperbolic

| op | speedup | max\_rel | precision | approach |
|---|---|---|---|---|
| sinh | **2.3×** | 4.6e-30 | full DD | original: Taylor series (|x|<0.1) or (exp−exp⁻¹)/2 |
| cosh | 1.6× | 4.6e-30 | full DD | original: (exp+exp⁻¹)/2 |
| tanh | **2.7×** | 6.6e-31 | full DD | original: sinh/cosh (|x|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) |
| asinh | **6.1×** | 4.5e-30 | full DD | original: Taylor series (|x|<0.01) or log(x+√(x²+1)) with Newton |
| acosh | **5.6×** | 2.6e-32 | full DD | original: log(x+√(x²−1)) with Newton correction |
| atanh | **5.1×** | 9.8e-31 | full DD | original: Taylor series (|x|<0.01) or ½ log((1+x)/(1−x)) |

### Error / special functions

| op | speedup | max\_rel | precision | approach |
|---|---|---|---|---|
| erf | **4.8×** | 1.4e-18 | deriv-corrected | original: Taylor series (|x|<2) + libm erf(hi) derivative correction |
| erfc | **4.5×** | 6.0e-16 | deriv-corrected | original: 1 − erf(x), or libm erfc(hi) for large x |
| erfc\_scaled | **162×** | 4.8e-16 | single-double | original: libm erfc\_scaled(hi), no lo correction |
| gamma | **112×** | 1.1e-16 | single-double | original: libm gamma(hi), no lo correction |
| log\_gamma | **68×** | 1.8e-16 | single-double | original: libm log\_gamma(hi), no lo correction |
| bessel\_j0 | **96×** | 2.2e-15 | single-double | original: libm bessel\_j0(hi), no lo correction |
| bessel\_j1 | **101×** | 9.1e-16 | single-double | original: libm bessel\_j1(hi), no lo correction |
| bessel\_jn(3,.) | **104×** | 1.6e-14 | single-double | original: libm bessel\_jn(3,hi), no lo correction |
| bessel\_y0 | **144×** | 6.7e-16 | single-double | original: libm bessel\_y0(hi), no lo correction |
| bessel\_y1 | **145×** | 8.5e-16 | single-double | original: libm bessel\_y1(hi), no lo correction |
| bessel\_yn(3,.) | **165×** | 4.7e-15 | single-double | original: libm bessel\_yn(3,hi), no lo correction |

### Complex arithmetic

| op | speedup | max\_rel | precision | approach |
|---|---|---|---|---|
| cx\_add | **3.4×** | 1.3e-32 | full DD | original: component-wise DD add |
| cx\_sub | **3.5×** | 5.8e-33 | full DD | original: component-wise DD sub |
| cx\_mul | **7.4×** | 1.9e-32 | full DD | original: (ac−bd, ad+bc) via DD ops |
| cx\_div | **7.3×** | 4.2e-32 (re) / 6.4e-17 (im) | full DD / deriv | original: (ac+bd, bc−ad)/(c²+d²) |
| cx\_conjg | 1.8× | exact | exact | original: negate im limbs |
| cx\_abs | **4.0×** | 6.7e-32 | full DD | original: hypot(re, im) |

### Complex transcendentals

| op | speedup | max\_rel | precision | approach |
|---|---|---|---|---|
| cx\_sqrt | **5.8×** | 6.5e-32 | full DD | original: Kahan-style (|z|+|a|)/2 with scaling |
| cx\_exp | **2.1×** | 1.0e-29 | full DD | original: exp(re)·(cos(im), sin(im)) |
| cx\_log | **2.1×** | 3.6e-32 | full DD | original: (log(|z|), atan2(im,re)) |
| cx\_sin | 1.9× | 1.5e-30 (re) / 9.3e-30 (im) | full DD | original: sin(re)cosh(im), cos(re)sinh(im) |
| cx\_cos | 1.9× | 9.3e-30 (re) / 1.5e-30 (im) | full DD | original: cos(re)cosh(im), −sin(re)sinh(im) |
| cx\_tan | 0.96× | 2.4e-30 | full DD | original: sin/cos complex |
| cx\_sinh | **2.1×** | 1.0e-29 | full DD | original: sinh(re)cos(im), cosh(re)sin(im) |
| cx\_cosh | **2.0×** | 1.0e-29 | full DD | original: cosh(re)cos(im), sinh(re)sin(im) |
| cx\_tanh | 1.1× | 1.5e-30 | full DD | original: complex tanh via sinh/cosh |
| cx\_asin | **2.5×** | 1.2e-23 (re) / 5.3e-31 (im) | deriv / full DD | original: −i·log(iz+√(1−z²)) |
| cx\_acos | **2.5×** | 1.3e-32 (re) / 5.3e-31 (im) | full DD | original: π/2 − asin(z) |
| cx\_atan | 1.8× | 3.9e-32 (re) / 1.3e-31 (im) | full DD | original: (i/2)·log((i+z)/(i−z)) |
| cx\_asinh | **2.5×** | 1.4e-22 (re) / 9.3e-32 (im) | deriv / full DD | original: log(z+√(z²+1)) |
| cx\_acosh | **2.2×** | 7.6e-31 (re) / 3.2e-32 (im) | full DD | original: log(z+√(z²−1)) |
| cx\_atanh | **2.1×** | 9.4e-23 (re) / 3.7e-32 (im) | deriv / full DD | original: ½·log((1+z)/(1−z)) |

### Array reductions

| op | speedup | max\_rel | precision | approach |
|---|---|---|---|---|
| arr\_sum (n=8) | 1.2× | 2.3e-31 | full DD | original: chained DD add |
| arr\_product (n=8) | **3.4×** | 1.1e-53 | full DD | original: chained DD mul |
| arr\_maxval (n=8) | **4.6×** | 5.6e-33 | full DD | original: chained DD compare |
| arr\_minval (n=8) | **4.5×** | 5.6e-33 | full DD | original: chained DD compare |
| arr\_dot (n=8) | **5.0×** | 2.1e-31 | full DD | original: fused multiply-accumulate with periodic renormalization |
| arr\_norm2 (n=8) | **5.6×** | 5.4e-32 | full DD | original: sqrt(dot(x,x)) |
| arr\_matmul (8×8\*8) | **2.1×** | 1.9e-30 | full DD | original: fused multiply-accumulate with periodic renormalization |

## C++: `MultiFloat<double,2>` vs `__float128`

Header-only — all kernels inline into the call site. No LTO needed.
Precision characteristics are the same as the Fortran version (same
algorithms), so only speedup is shown.

### Arithmetic

| op | speedup |
|---|---|
| add | **3.0×** |
| sub | **4.4×** |
| mul | **11×** |
| div | **3.9×** |
| sqrt | **61×** |
| cbrt | **43×** |
| fma | **98×** |
| abs | **2.1×** |
| neg | **2.1×** |

### Rounding

| op | speedup |
|---|---|
| floor | **2.9×** |
| ceil | **3.2×** |
| trunc | **2.6×** |
| round | 1.0× |
| rint | **11×** |
| nearbyint | **34×** |

### Binary

| op | speedup |
|---|---|
| fmin | **5.2×** |
| fmax | **5.4×** |
| fdim | **6.5×** |
| copysign | 1.9× |
| fmod | 0.66× |
| hypot | **42×** |
| ldexp(.,5) | **2.1×** |

### Exponential / logarithmic

| op | speedup |
|---|---|
| exp | **3.1×** |
| exp2 | **3.6×** |
| expm1 | **4.6×** |
| log | **5.0×** |
| log10 | **6.7×** |
| log2 | **6.2×** |
| log1p | **5.5×** |
| pow | **4.5×** |

### Trigonometric

| op | speedup |
|---|---|
| sin | **3.2×** |
| cos | **3.1×** |
| tan | 1.4× |
| asin | **2.1×** |
| acos | **2.1×** |
| atan | 1.3× |
| atan2 | 1.1× |

### Hyperbolic

| op | speedup |
|---|---|
| sinh | **2.6×** |
| cosh | 1.8× |
| tanh | **2.8×** |
| asinh | **6.9×** |
| acosh | **6.3×** |
| atanh | **5.4×** |

### Error / special functions

| op | speedup |
|---|---|
| erf | **5.0×** |
| erfc | **5.2×** |
| tgamma | **105×** |
| lgamma | **71×** |

## Notes

- **`mod` / `fmod`** is the only operation where quad precision is
  consistently faster. The DD `mod` uses a floor-multiple reduction loop
  (adapted from `external/float64x2-sample.cpp`) for small quotients, or a
  full DD divide chain for large quotients; libquadmath's `fmodq` uses a
  specialized bit-level remainder algorithm. `modulo` (Fortran) now beats qp
  at 1.3× thanks to the iterative approach. Precision degrades as
  `~10^(log10(quotient) - 31)` for large quotients, which is the inherent DD
  precision limit.

- **Trig range reduction** uses a 3-part 1/π constant (~161 bits) via
  Cody–Waite reduction. This keeps sin/cos/tan accurate to ~1e-25 for
  |x| < 1e6 (the fuzz test range). For larger arguments, a Payne–Hanek
  reduction with a multi-word 2/π table would be needed.

- **Single-double precision** functions (gamma, bessel, erfc\_scaled, etc.)
  achieve 60–165× speedup by evaluating `f(hi)` via the leading-limb libm
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
