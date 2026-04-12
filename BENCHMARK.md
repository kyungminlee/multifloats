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

## Fortran: `float64x2` vs `real(16)`

Each operation is timed over 1024 elements × 400 repetitions (fast ops)
or fewer reps (transcendentals), with a NOINLINE drain after each rep to
prevent dead-code elimination. The "speedup" column is `qp_time / mf_time`:
values > 1× mean multifloats is faster.

### Arithmetic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| add | 0.0029 | 0.0016 | **1.9×** |
| sub | 0.0038 | 0.0018 | **2.1×** |
| mul | 0.0101 | 0.0012 | **8.6×** |
| div | 0.0138 | 0.0045 | **3.1×** |
| sqrt | 0.1532 | 0.0040 | **39×** |
| add (mf+dp) | 0.0028 | 0.0011 | **2.6×** |
| mul (dp\*mf) | 0.0067 | 0.0008 | **8.9×** |

### Unary

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| abs | 0.0011 | 0.0008 | 1.4× |
| neg | 0.0011 | 0.0005 | **2.2×** |
| aint | 0.0020 | 0.0012 | 1.6× |
| anint | 0.0022 | 0.0009 | **2.4×** |
| fraction | 0.0031 | 0.0022 | 1.4× |
| scale | 0.0018 | 0.0018 | 1.0× |
| set\_exponent | 0.0038 | 0.0022 | 1.8× |

### Binary

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| min | 0.0021 | 0.0010 | **2.2×** |
| max | 0.0021 | 0.0010 | **2.1×** |
| min3 | 0.0033 | 0.0015 | **2.2×** |
| max3 | 0.0057 | 0.0013 | **4.3×** |
| sign | 0.0015 | 0.0013 | 1.2× |
| dim | 0.0056 | 0.0020 | **2.8×** |
| hypot | 0.1997 | 0.0422 | **4.7×** |
| mod | 0.0027 | 0.0050 | 0.53× |
| modulo | 0.0085 | 0.0066 | **1.3×** |

### Exponential / logarithmic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| exp | 0.0216 | 0.0075 | **2.9×** |
| log | 0.0271 | 0.0054 | **5.0×** |
| log10 | 0.0358 | 0.0055 | **6.6×** |
| pow | 0.0682 | 0.0159 | **4.3×** |
| pow\_int | 0.0023 | 0.0001 | **16×** |

### Trigonometric

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| sin | 0.0215 | 0.0078 | **2.7×** |
| cos | 0.0217 | 0.0081 | **2.7×** |
| tan | 0.0212 | 0.0169 | 1.3× |
| asin | 0.0361 | 0.0183 | **2.0×** |
| acos | 0.0374 | 0.0201 | 1.9× |
| atan | 0.0226 | 0.0191 | 1.2× |
| atan2 | 0.0253 | 0.0214 | 1.2× |

### Hyperbolic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| sinh | 0.0364 | 0.0155 | **2.3×** |
| cosh | 0.0254 | 0.0160 | 1.6× |
| tanh | 0.0375 | 0.0139 | **2.7×** |
| asinh | 0.0538 | 0.0088 | **6.1×** |
| acosh | 0.0480 | 0.0085 | **5.6×** |
| atanh | 0.0417 | 0.0082 | **5.1×** |

### Error / special functions

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| erf | 0.0069 | 0.0014 | **4.8×** |
| erfc | 0.0067 | 0.0015 | **4.5×** |
| erfc\_scaled | 0.0091 | 0.0001 | **162×** |
| gamma | 0.0109 | 0.0001 | **112×** |
| log\_gamma | 0.0033 | 0.0000 | **68×** |
| bessel\_j0 | 0.0092 | 0.0001 | **96×** |
| bessel\_j1 | 0.0091 | 0.0001 | **101×** |
| bessel\_jn(3,.) | 0.0273 | 0.0003 | **104×** |
| bessel\_y0 | 0.0142 | 0.0001 | **144×** |
| bessel\_y1 | 0.0141 | 0.0001 | **145×** |
| bessel\_yn(3,.) | 0.0291 | 0.0002 | **165×** |

### Complex arithmetic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| cx\_add | 0.0088 | 0.0026 | **3.4×** |
| cx\_sub | 0.0097 | 0.0028 | **3.5×** |
| cx\_mul | 0.0288 | 0.0039 | **7.4×** |
| cx\_div | 0.0923 | 0.0126 | **7.3×** |
| cx\_conjg | 0.0019 | 0.0011 | 1.8× |
| cx\_abs | 0.1978 | 0.0490 | **4.0×** |

### Complex transcendentals

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| cx\_sqrt | 0.0040 | 0.0007 | **5.8×** |
| cx\_exp | 0.0049 | 0.0023 | **2.1×** |
| cx\_log | 0.0077 | 0.0037 | **2.1×** |
| cx\_sin | 0.0089 | 0.0048 | 1.9× |
| cx\_cos | 0.0090 | 0.0047 | 1.9× |
| cx\_tan | 0.0095 | 0.0098 | 0.96× |
| cx\_sinh | 0.0094 | 0.0046 | **2.1×** |
| cx\_cosh | 0.0095 | 0.0046 | **2.0×** |
| cx\_tanh | 0.0101 | 0.0096 | 1.1× |
| cx\_asin | 0.0118 | 0.0047 | **2.5×** |
| cx\_acos | 0.0118 | 0.0047 | **2.5×** |
| cx\_atan | 0.0072 | 0.0040 | 1.8× |
| cx\_asinh | 0.0120 | 0.0047 | **2.5×** |
| cx\_acosh | 0.0121 | 0.0056 | **2.2×** |
| cx\_atanh | 0.0084 | 0.0040 | **2.1×** |

### Array reductions

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| arr\_sum (n=8) | 0.0283 | 0.0238 | 1.2× |
| arr\_product (n=8) | 0.0395 | 0.0117 | **3.4×** |
| arr\_maxval (n=8) | 0.0086 | 0.0019 | **4.6×** |
| arr\_minval (n=8) | 0.0087 | 0.0020 | **4.5×** |
| arr\_dot (n=8) | 0.0304 | 0.0061 | **5.0×** |
| arr\_norm2 (n=8) | 0.1514 | 0.0270 | **5.6×** |
| arr\_matmul (8×8\*8) | 0.0308 | 0.0144 | **2.1×** |

## C++: `MultiFloat<double,2>` vs `__float128`

Header-only — all kernels inline into the call site. No LTO needed.

### Arithmetic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| add | 0.0031 | 0.0011 | **3.0×** |
| sub | 0.0048 | 0.0011 | **4.4×** |
| mul | 0.0075 | 0.0007 | **11×** |
| div | 0.0115 | 0.0029 | **3.9×** |
| sqrt | 0.1511 | 0.0025 | **61×** |
| cbrt | 0.2319 | 0.0054 | **43×** |
| fma | 0.1017 | 0.0010 | **98×** |
| abs | 0.0014 | 0.0007 | **2.1×** |
| neg | 0.0011 | 0.0005 | **2.1×** |

### Rounding

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| floor | 0.0020 | 0.0007 | **2.9×** |
| ceil | 0.0020 | 0.0006 | **3.2×** |
| trunc | 0.0023 | 0.0009 | **2.6×** |
| round | 0.0021 | 0.0021 | 1.0× |
| rint | 0.0067 | 0.0006 | **11×** |
| nearbyint | 0.0227 | 0.0007 | **34×** |

### Binary

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| fmin | 0.0048 | 0.0009 | **5.2×** |
| fmax | 0.0046 | 0.0008 | **5.4×** |
| fdim | 0.0059 | 0.0009 | **6.5×** |
| copysign | 0.0014 | 0.0008 | 1.9× |
| fmod | 0.0029 | 0.0044 | 0.66× |
| hypot | 0.2027 | 0.0048 | **42×** |
| ldexp(.,5) | 0.0037 | 0.0018 | **2.1×** |

### Exponential / logarithmic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| exp | 0.0220 | 0.0072 | **3.1×** |
| exp2 | 0.0236 | 0.0066 | **3.6×** |
| expm1 | 0.0334 | 0.0073 | **4.6×** |
| log | 0.0272 | 0.0054 | **5.0×** |
| log10 | 0.0362 | 0.0054 | **6.7×** |
| log2 | 0.0347 | 0.0056 | **6.2×** |
| log1p | 0.0317 | 0.0058 | **5.5×** |
| pow | 0.0687 | 0.0154 | **4.5×** |

### Trigonometric

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| sin | 0.0221 | 0.0070 | **3.2×** |
| cos | 0.0218 | 0.0071 | **3.1×** |
| tan | 0.0203 | 0.0147 | 1.4× |
| asin | 0.0363 | 0.0175 | **2.1×** |
| acos | 0.0369 | 0.0178 | **2.1×** |
| atan | 0.0226 | 0.0174 | 1.3× |
| atan2 | 0.0252 | 0.0232 | 1.1× |

### Hyperbolic

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| sinh | 0.0376 | 0.0143 | **2.6×** |
| cosh | 0.0252 | 0.0141 | 1.8× |
| tanh | 0.0366 | 0.0131 | **2.8×** |
| asinh | 0.0531 | 0.0077 | **6.9×** |
| acosh | 0.0478 | 0.0076 | **6.3×** |
| atanh | 0.0421 | 0.0078 | **5.4×** |

### Error / special functions

| op | qp \[s\] | mf \[s\] | speedup |
|---|---|---|---|
| erf | 0.0069 | 0.0014 | **5.0×** |
| erfc | 0.0070 | 0.0013 | **5.2×** |
| tgamma | 0.0105 | 0.0001 | **105×** |
| lgamma | 0.0032 | 0.0000 | **71×** |

## Notes

- **`mod` / `fmod`** is the only operation where quad precision is
  consistently faster. The DD `mod` uses a floor-multiple reduction loop
  (for small quotients) or a full DD divide chain (for large quotients);
  libquadmath's `fmodq` uses a specialized bit-level remainder algorithm.
  `modulo` (Fortran) now beats qp at 1.3× thanks to the iterative approach.
  Precision degrades as `~10^(log10(quotient) - 31)` for large quotients,
  which is the inherent DD precision limit.

- The **Fortran** multifloats module uses `elemental` functions on a
  `sequence` derived type. gfortran's ABI passes/returns these via hidden
  pointers (not in FP registers), adding ~1.5× overhead vs the C++ header-
  only version even with LTO inlining. See the performance note at the top
  of `fsrc/multifloats.fypp` for details and the `bind(c)` escape hatch.

- **Special functions** (gamma, bessel, erfc\_scaled) achieve 60–165×
  speedup because the DD kernels use a leading-limb libm call + derivative
  correction, which is inherently cheaper than libquadmath's full software
  quad implementation of these functions.

- **Array reductions** (dot\_product, matmul) use a fused multiply-
  accumulate kernel that computes the product's error-free representation
  and accumulates corrections into a scalar `s_lo`, with periodic
  renormalization (configurable via `mf_set_fma_renorm_interval`).
