# Benchmark Results

Comparison of multifloats double-double (DD) arithmetic against quad
precision (`real(16)` / `__float128` via libquadmath). The multifloats
implementations use explicit two-limb error-free transformations (EFTs)
on pairs of `double` values, achieving ~106 bits of significand — close
to quad precision's 113 bits — at a fraction of the cost.

## System

| CPU | OS | Compiler | Build |
|---|---|---|---|
| Apple M1 Max (arm64), 10 cores | macOS 26.3 (Darwin 25.3.0) | GNU Fortran 15.2.0 (Homebrew GCC 15.2.0_1) / g++ Apple clang version 17.0.0 (clang-1700.6.4.2) | CMake 4.3.1, `-O3 -flto`, STATIC library |

## Precision key

DD precision is the maximum relative error vs an **mpreal @ 200-bit**
oracle over ~1M random inputs (fixed seed 42), reported in **DD ulps**
(1 ulp ≈ 2⁻¹⁰⁵ ≈ 2.46e-32 for the ~106-bit DD significand). The qp
column is the same metric for the libquadmath leg, in **qp ulps** (1 qp
ulp ≈ 2⁻¹¹² ≈ 1.93e-34).

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

## Operations

`scope` records which language surface the op is exposed through
(`C` = `include/multifloats.h`, `Fortran` = `fsrc/multifloats.fypp`,
`both` = both). The speedup column is `C-ABI × / Fortran-elemental ×`
(`qp_time / dd_time`); a `—` on one side means the op is not exposed
through that surface. The Fortran side is typically ~1.5× lower than
the C side for shared ops because gfortran returns DD via hidden
pointers rather than FP registers. **Bold** marks ≥ 2×.

| op | scope | prec | err_dd | err_qp | speedup (C / F) | approach |
|---|---|---|---|---|---|---|
| add | both | full DD | 0.7 ulp | 0.5 qp ulp | **2.8×** / 1.5× | Julia: two\_sum EFT |
| sub | both | full DD | 0.7 ulp | 0.5 qp ulp | **2.6×** / 1.4× | Julia: two\_sum EFT (negate + add) |
| mul | both | full DD | 2.0 ulp | 0.5 qp ulp | **10×** / **3.1×** | Julia: two\_prod EFT via FMA |
| div | both | full DD | 2.8 ulp | 0.5 qp ulp | **8.0×** / **3.8×** | original: Newton refinement (1/y seed, one step) |
| sqrt | both | full DD | 1.7 ulp | 0.5 qp ulp | **47×** / **35×** | Julia: Karp–Markstein (reciprocal sqrt seed + Newton) |
| fma | C |  | 2.2 ulp | 0.5 qp ulp | **78×** / — | original: x\*y + z via DD ops |
| abs | both | exact | exact | exact | **2.0×** / 0.68× | original: sign-check + negate limbs |
| neg | both | exact | exact | exact | **2.5×** / 0.94× | original: negate both limbs |
| trunc | C |  | exact | exact | 1.4× / — | original: signbit ? −floor(−x) : floor(x) |
| round | C |  | exact | exact | 1.8× / — | original: trunc(x + ½·sign(x)) |
| fmin | C |  | exact | exact | **4.3×** / — | original: DD comparison + select |
| fmax | C |  | exact | exact | **5.3×** / — | original: DD comparison + select |
| fdim | C |  | 0.7 ulp | 0.5 qp ulp | **3.3×** / — | original: DD comparison, then subtract or zero |
| copysign | C |  | exact | exact | 1.8× / — | original: sign-bit copy to hi, propagate to lo |
| fmod | C |  | 1.0 ulp | exact | 0.66× / — | sample: floor-multiple reduction loop; fallback to div chain |
| hypot | both | full DD | 2.2 ulp | 0.7 qp ulp | **9.7×** / **9.1×** | original: scaled sqrt(x²+y²) |
| exp | both | full DD | 1.6 ulp | 0.5 qp ulp | **3.7×** / **3.4×** | Julia: exp2 polynomial (14-term Horner) + ldexp reconstruction |
| exp2 | C |  | 1.2 ulp | 0.5 qp ulp | **4.8×** / — | Julia: exp2 polynomial (14-term Horner) |
| expm1 | C |  | 1.7 ulp | 0.7 qp ulp | **5.3×** / — | original: exp(x) − 1 via DD sub |
| log | both | full DD | 1.8 ulp | 0.6 qp ulp | **6.4×** / **6.2×** | Julia: log2 table lookup (32 centers) + polynomial (7-term Horner) |
| log10 | both | full DD | 1.2 ulp | 0.5 qp ulp | **8.4×** / **8.2×** | Julia: log2 kernel × DD log10(2) |
| log2 | C |  | 0.8 ulp | 0.5 qp ulp | **8.0×** / — | Julia: log2 table lookup + polynomial |
| log1p | C |  | 2.6 ulp | 0.6 qp ulp | **6.1×** / — | original: log(1 + x) via DD add |
| pow | both | full DD | 11 ulp | 0.5 qp ulp | **5.7×** / **5.4×** | Julia: exp(y × log(x)) |
| sin | both | full DD | 1.6 ulp | 0.5 qp ulp | **3.3×** / **3.2×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| cos | both | full DD | 0.7 ulp | 0.3 qp ulp | **3.2×** / **3.2×** | original: 13-term Taylor Horner + 3-part Cody–Waite π/2 + π/8 split |
| sincos | both | full DD | 1.6 ulp (sin) / 0.7 ulp (cos) | 0.5 qp ulp (sin) / 0.3 qp ulp (cos) | **5.3×** / **5.1×** | fused sin/cos: one range-reduction, two Taylor outputs |
| sinpi | both | full DD | 1.1 ulp | 1.9 qp ulp | **4.5×** / **4.4×** | Julia: sinpi Horner polynomial, direct |
| cospi | both | full DD | 1.4 ulp | 1.8 qp ulp | **4.2×** / **4.2×** | Julia: cospi Horner polynomial, direct |
| tan | both | full DD | 2.3 ulp | 0.5 qp ulp | **2.0×** / **2.2×** | original: sin/cos Taylor kernels + DD divide |
| tanpi | both | full DD | 2.0 ulp | 35287 qp ulp | **3.3×** / **3.1×** | original: sinpi/cospi ratio |
| asin | both | full DD | 0.7 ulp | 0.5 qp ulp | **5.7×** / **11×** | original: piecewise rational P/Q (3 regions, from libquadmath asinq.c) |
| asinpi | both | full DD | 0.5 ulp | 0.3 qp ulp | **5.7×** / **5.8×** | original: asin(x)/π with exact-DD π division |
| acos | both | full DD | 0.5 ulp | 0.6 qp ulp | **5.6×** / **6.2×** | original: asin polynomial + half-angle identity |
| acospi | both | full DD | 1.7 ulp | 1.0 qp ulp | **5.7×** / **5.8×** | original: acos(x)/π with exact-DD π division |
| atan | both | full DD | 0.9 ulp | 0.6 qp ulp | **3.7×** / **3.6×** | original: 84-entry table lookup + rational P(t²)/Q(t²) (from libquadmath atanq.c) |
| atanpi | both | full DD | 0.5 ulp | 0.3 qp ulp | **4.1×** / **4.3×** | original: atan(x)/π with exact-DD π division |
| atan2 | both | full DD | 1.6 ulp | 1.1 qp ulp | **3.3×** / **3.3×** | original: table-based atan + quadrant correction |
| atan2pi | both | full DD | 1.6 ulp | 1.2 qp ulp | **3.7×** / **3.5×** | original: atan2(y,x)/π with exact-DD π division |
| sinh | both | full DD | 4.0 ulp | 0.9 qp ulp | **3.2×** / **3.0×** | original: Taylor series (\|x\|<0.1) or (exp−exp⁻¹)/2 |
| cosh | both | full DD | 1.3 ulp | 0.9 qp ulp | **2.0×** / **2.0×** | original: (exp+exp⁻¹)/2 |
| sinhcosh | both | full DD | 4.0 ulp (sinh) / 1.3 ulp (cosh) | 0.9 qp ulp (sinh) / 0.9 qp ulp (cosh) | **4.7×** / **4.5×** | fused sinh/cosh: one range-reduction, two outputs |
| tanh | both | full DD | 3.4 ulp | 1.2 qp ulp | **3.6×** / **3.3×** | original: sinh/cosh (\|x\|<0.5) or (1−e⁻²ˣ)/(1+e⁻²ˣ) |
| asinh | both | full DD | 3.0 ulp | 1.6 qp ulp | **7.2×** / **7.2×** | original: Taylor series (\|x\|<0.01) or log(x+√(x²+1)) with Newton |
| acosh | both | full DD | 0.5 ulp | 0.5 qp ulp | **7.4×** / **6.5×** | original: log(x+√(x²−1)) with Newton correction |
| atanh | both | full DD | 2.0 ulp | 1.3 qp ulp | **6.4×** / **6.2×** | original: Taylor series (\|x\|<0.01) or ½·log((1+x)/(1−x)) |
| erf | both | full DD | 0.6 ulp | 0.7 qp ulp | **4.1×** / **4.0×** | piecewise rational approx (libquadmath erfq.c) |
| erfc | both | full DD | 0.4 ulp | 0.6 qp ulp | **3.2×** / **4.0×** | piecewise rational approx + split exp(-x^2) |
| erfcx | C |  | 290 ulp | 256 qp ulp | **3.3×** / — | exp(x²)·erfc(x); scaled form avoiding tail cancellation |
| tgamma | C |  | 13 ulp | 1.6 qp ulp | **8.8×** / — | piecewise rational approx + Stirling + reflection, exp(lgamma) |
| lgamma | C |  | 1.9 ulp | 1.3 qp ulp | **4.8×** / — | piecewise rational approx + Stirling asymptotic |
| bj0 | C |  | 0.2 ulp | 0.3 qp ulp | **6.6×** / — | piecewise rational + Hankel asymptotic (j0q.c) |
| bj1 | C |  | 0.4 ulp | 0.4 qp ulp | **7.1×** / — | piecewise rational + Hankel asymptotic (j1q.c) |
| bjn(3,.) | C |  | 0.5 ulp | 0.5 qp ulp | **4.3×** / — | forward/backward recurrence from j0/j1 |
| by0 | C |  | 1.8 ulp | 1.6 qp ulp | **5.2×** / — | piecewise rational + Hankel asymptotic (j0q.c) |
| by1 | C |  | 1.7 ulp | 1.5 qp ulp | **6.9×** / — | piecewise rational + Hankel asymptotic (j1q.c) |
| byn(3,.) | C |  | 7.4 ulp | 4.4 qp ulp | **6.8×** / — | forward recurrence from y0/y1 |
| yn\_range(0..5) | C |  | 6.0 ulp | 3.2 qp ulp | **26×** / — | single forward-recurrence sweep, 6 outputs / call |
| cdd\_add | both | full DD | 0.6 ulp (re) / 0.6 ulp (im) | 0.5 qp ulp (re) / 0.5 qp ulp (im) | **2.8×** / **2.5×** | original: component-wise DD add |
| cdd\_sub | both | full DD | 0.7 ulp (re) / 0.7 ulp (im) | 0.5 qp ulp (re) / 0.5 qp ulp (im) | **2.7×** / **2.5×** | original: component-wise DD sub |
| cdd\_mul | both | full DD | 1.8 ulp (re) / 2.1 ulp (im) | 0.8 qp ulp (re) / 0.9 qp ulp (im) | **8.0×** / 1.9× | original: (ac−bd, ad+bc) via DD ops |
| cdd\_div | both | full DD / deriv | 6.7 ulp (re) / 2.1 ulp (im) | 5.7 qp ulp (re) / 28 qp ulp (im) | **5.2×** / **4.0×** | original: (ac+bd, bc−ad)/(c²+d²) |
| cdd\_conjg | both | exact | exact | exact | 1.2× / 1.4× | original: negate im limbs |
| cdd\_proj | C |  | exact (re) / exact (im) | exact (re) / exact (im) | **2.7×** / — | C99 Annex G Riemann-sphere projection (identity for finite z) |
| cdd\_abs | both | full DD | 2.0 ulp | 0.7 qp ulp | **10×** / **9.3×** | original: hypot(re, im) |
| cdd\_arg | C |  | 1.7 ulp | 1.0 qp ulp | **3.1×** / — | original: atan2(im, re) |
| cdd\_sqrt | both | full DD | 1.7 ulp (re) / 2.1 ulp (im) | 0.9 qp ulp (re) / 1.0 qp ulp (im) | **6.9×** / **7.3×** | original: Kahan-style (\|z\|+\|a\|)/2 with scaling |
| cdd\_exp | both | full DD | 2.2 ulp (re) / 2.3 ulp (im) | 1.2 qp ulp (re) / 1.2 qp ulp (im) | **3.6×** / **3.8×** | original: exp(re)·(cos(im), sin(im)) |
| cdd\_expm1 | both | reduced DD | 1.7 ulp (re) / 2.3 ulp (im) | 1892 qp ulp (re) / 1.2 qp ulp (im) | 1.7× / 1.4× | original: expm1 + complex rotation with cancellation-safe Re |
| cdd\_log | both | full DD | 1.4 ulp (re) / 1.7 ulp (im) | 1.0 qp ulp (re) / 1.0 qp ulp (im) | **7.2×** / **4.8×** | original: (log(\|z\|), atan2(im,re)) |
| cdd\_log2 | both | full DD | 2.3 ulp (re) / 1.6 ulp (im) | 1.2 qp ulp (re) / 1.2 qp ulp (im) | **5.1×** / **5.0×** | original: clog / log(2) component-wise |
| cdd\_log10 | both | full DD | 1.6 ulp (re) / 1.8 ulp (im) | 0.9 qp ulp (re) / 1.3 qp ulp (im) | **4.7×** / **5.0×** | original: clog / log(10) component-wise |
| cdd\_log1p | both | full DD | 31 ulp (re) / 2.1 ulp (im) | 1188 qp ulp (re) / 1.1 qp ulp (im) | **5.5×** / **5.8×** | original: cancellation-safe log(1+z) near z=0 |
| cdd\_pow | both | full DD | 394 ulp (re) / 944 ulp (im) | 533 qp ulp (re) / 1435 qp ulp (im) | **3.6×** / **4.0×** | original: exp(w·log(z)) |
| cdd\_sin | both | full DD | 1.8 ulp (re) / 7.0 ulp (im) | 1.3 qp ulp (re) / 1.3 qp ulp (im) | **4.3×** / **4.0×** | original: sin(re)cosh(im), cos(re)sinh(im) |
| cdd\_cos | both | full DD | 1.7 ulp (re) / 1.8 ulp (im) | 1.3 qp ulp (re) / 1.5 qp ulp (im) | **4.2×** / **4.0×** | original: cos(re)cosh(im), −sin(re)sinh(im) |
| cdd\_tan | both | full DD | 2.7 ulp (re) / 8.1 ulp (im) | 1.4 qp ulp (re) / 2.1 qp ulp (im) | **4.4×** / **4.2×** | original: complex sin/cos ratio |
| cdd\_sinpi | C |  | 88 ulp (re) / 87 ulp (im) | 10722 qp ulp (re) / 2234 qp ulp (im) | **4.0×** / — | original: csin(π·z) via π-scaled trig kernels |
| cdd\_cospi | C |  | 87 ulp (re) / 88 ulp (im) | 2234 qp ulp (re) / 10722 qp ulp (im) | **3.9×** / — | original: ccos(π·z) via π-scaled trig kernels |
| cdd\_sinh | both | full DD | 3.8 ulp (re) / 2.3 ulp (im) | 1.7 qp ulp (re) / 1.5 qp ulp (im) | **4.2×** / **4.1×** | original: sinh(re)cos(im), cosh(re)sin(im) |
| cdd\_cosh | both | full DD | 2.3 ulp (re) / 2.3 ulp (im) | 1.5 qp ulp (re) / 1.3 qp ulp (im) | **4.2×** / **4.5×** | original: cosh(re)cos(im), sinh(re)sin(im) |
| cdd\_tanh | both | full DD | 4.7 ulp (re) / 1.8 ulp (im) | 1.8 qp ulp (re) / 1.7 qp ulp (im) | **4.5×** / **4.4×** | original: complex tanh via sinh/cosh |
| cdd\_asin | both | deriv / full DD | 1.9 ulp (re) / 1.8 ulp (im) | 1.5 qp ulp (re) / 1.7 qp ulp (im) | **4.8×** / **5.0×** | original: −i·log(iz+√(1−z²)) |
| cdd\_acos | both | full DD | 1.1 ulp (re) / 2.5 ulp (im) | 1.0 qp ulp (re) / 1.7 qp ulp (im) | **4.1×** / **4.4×** | original: π/2 − asin(z) |
| cdd\_atan | both | full DD | 1.4 ulp (re) / 2.3 ulp (im) | 1.1 qp ulp (re) / 1.4 qp ulp (im) | **4.0×** / **4.3×** | original: (i/2)·log((i+z)/(i−z)) |
| cdd\_asinh | both | deriv / full DD | 2.2 ulp (re) / 1.8 ulp (im) | 1.8 qp ulp (re) / 1.5 qp ulp (im) | **4.4×** / **4.9×** | original: log(z+√(z²+1)) |
| cdd\_acosh | both | full DD | 2.5 ulp (re) / 1.1 ulp (im) | 1.7 qp ulp (re) / 1.0 qp ulp (im) | **3.8×** / **4.5×** | original: log(z+√(z²−1)) |
| cdd\_atanh | both | deriv / full DD | 2.2 ulp (re) / 1.4 ulp (im) | 1.7 qp ulp (re) / 1.1 qp ulp (im) | **4.8×** / **5.0×** | original: ½·log((1+z)/(1−z)) |
| add (dd+dp) | Fortran | exact | exact | exact | — / 1.6× | Julia: two\_sum EFT |
| mul (dp\*dd) | Fortran | full DD | 1.1 ulp | 0.5 qp ulp | — / **2.6×** | Julia: two\_prod EFT via FMA |
| aint | Fortran | exact | exact | exact | — / 0.67× | original: truncate hi, check DD fractional part |
| anint | Fortran | exact | exact | exact | — / 1.1× | original: truncate hi, DD fractional part vs ±0.5 |
| fraction | Fortran | exact | exact | exact | — / 1.2× | original: scale both limbs by −exponent |
| scale | Fortran | exact | exact | exact | — / 0.74× | original: ldexp on both limbs |
| set\_exponent | Fortran | exact | exact | exact | — / 1.4× | original: scale + set\_exponent on hi |
| min | Fortran | full DD | exact | exact | — / 0.76× | original: DD comparison + select |
| max | Fortran | full DD | exact | exact | — / 1.4× | original: DD comparison + select |
| min3 | Fortran | full DD | exact | exact | — / 1.8× | original: chained min |
| max3 | Fortran | full DD | exact | exact | — / **2.4×** | original: chained max |
| sign | Fortran | exact | exact | exact | — / 0.61× | original: sign-check + negate |
| dim | Fortran | full DD | 0.7 ulp | 0.5 qp ulp | — / 1.9× | original: DD comparison, then subtract or zero |
| mod | Fortran | full DD | 1.0 ulp | exact | — / 0.41× | sample: floor-multiple reduction loop; fallback to div chain |
| modulo | Fortran | full DD | 1.0 ulp | exact | — / 0.80× | original: mod + sign adjustment |
| pow\_int | Fortran | full DD | 11 ulp | 0.5 qp ulp | — / **13×** | original: repeated squaring via DD mul |
| erfc\_scaled | Fortran | full DD | 290 ulp | 256 qp ulp | — / **5.5×** | exp(x^2)·erfc(x) with asymptotic cancellation |
| gamma | Fortran | full DD | 13 ulp | 1.6 qp ulp | — / **7.5×** | piecewise rational approx + Stirling + reflection |
| log\_gamma | Fortran | full DD | 1.9 ulp | 1.3 qp ulp | — / **4.9×** | piecewise rational approx + Stirling asymptotic |
| bessel\_j0 | Fortran | full DD | 0.2 ulp | 0.3 qp ulp | — / **6.3×** | piecewise rational + Hankel asymptotic (j0q.c) via C++ |
| bessel\_j1 | Fortran | full DD | 0.4 ulp | 0.4 qp ulp | — / **6.3×** | piecewise rational + Hankel asymptotic (j1q.c) via C++ |
| bessel\_jn(3,.) | Fortran | full DD | 0.5 ulp | 0.5 qp ulp | — / **4.0×** | forward/backward recurrence from j0/j1 |
| bessel\_y0 | Fortran | full DD | 1.8 ulp | 1.6 qp ulp | — / **6.5×** | piecewise rational + Hankel asymptotic (j0q.c) via C++ |
| bessel\_y1 | Fortran | full DD | 1.7 ulp | 1.5 qp ulp | — / **6.1×** | piecewise rational + Hankel asymptotic (j1q.c) via C++ |
| bessel\_yn(3,.) | Fortran | full DD | 7.4 ulp | 4.4 qp ulp | — / **6.0×** | forward recurrence from y0/y1 |
| arr\_sum (n=8) | Fortran | full DD | 1.3 ulp | 1.5 qp ulp | — / **3.2×** | original: chained DD add |
| arr\_product (n=8) | Fortran | full DD | 0.0 ulp | 0.0 qp ulp | — / 1.7× | original: chained DD mul |
| arr\_maxval (n=8) | Fortran | full DD | exact | exact | — / 1.9× | original: chained DD compare |
| arr\_minval (n=8) | Fortran | full DD | exact | exact | — / 1.6× | original: chained DD compare |
| arr\_dot (n=8) | Fortran | full DD | 0.2 ulp | 0.2 qp ulp | — / **3.6×** | original: fused multiply-accumulate with periodic renormalization |
| arr\_norm2 (n=8) | Fortran | full DD | 1.0 ulp | 0.9 qp ulp | — / **3.9×** | original: sqrt(dot(x,x)) |
| arr\_matmul (8×8·8) | Fortran | full DD | 0.4 ulp | 0.2 qp ulp | — / **2.6×** | original: AXPY-order C kernel, MR=8 register-blocked panels + 1..7 tail, periodic renorm |

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

<!-- Auto-generated by bench/build_benchmark_md.py from a per-system JSON
     result in bench/results/. Do not edit by hand. -->
