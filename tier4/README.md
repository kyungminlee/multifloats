# tier4 — Branch-cut correctness scratch folder

Measurements and reproducers from the Tier 4 follow-up pass that fixed
two complex-kernel signed-zero bugs found by the Tier 4 branch-cut test
(`test_complex_branch_cuts` in `test/test.cc`).

## Files

| File                | What it is                                                                |
| ------------------- | ------------------------------------------------------------------------- |
| `baseline.txt`      | Pre-fix cpp_bench, cpp_test precision, fortran_fuzz for complex kernels.  |
| `fbench_before.txt` | Pre-fix fortran_bench complex-op timing snapshot.                         |
| `after.txt`         | Post-fix numbers on the same targets.                                     |
| `catanh_probe.cc`   | Tiny repro for `catanh(±1 + 0i) = NaN` (should be ±∞).                    |
| `casin_probe.cc`    | Signed-zero table: `casin({±2} × {+0, -0})`, principal case `(0.5, 0.5)`. |
| `casin_trace.cc`    | Step-by-step DD trace showing where the signed zero was lost.             |

## Fixes applied (committed together with this folder)

### `catanhdd(±1 + 0i)` returned NaN instead of ±∞

Root cause in `src/multifloats_math.cc:catanhdd`: at z = +1 + 0i the
formula `Re = 0.25 · log((b² + (a+1)²) / (b² + (a-1)²))` hits `num/den
= 4/0 = +inf` (DD division propagates correctly to {inf, inf}). But the
final `0.25 · inf` DD multiply emits `fma(0.25, inf, -inf) = NaN` in
the compensated-error term, turning the result into NaN. Added an
explicit short-circuit at the branch-point singularities `|a| == 1 &&
b == 0` that returns `(copysign(∞, a), +0)`.

### `casindd` / `cacosdd` lost the signed-zero imaginary input

Root cause in `src/multifloats_math.cc:casindd`: the formula
`Im(1 − z²) = −(a·b + a·b)` relies on DD multiply preserving signed
zero. But DD multiply of `(2, 0) * (-0, 0)` collapses to `(+0, 0)`
(`two_prod(2, -0)` loses the sign), so `csqrtdd` downstream picks the
+imag branch regardless of input sign. Fix: negate the entire DD pair
(both limbs) at the end of `casindd` when the computed hi sign
disagrees with `Im(z).hi`. A naive `copysign(lo, input_im)` would
corrupt the DD pair by ~1 dp ULP because the lo limb normally carries
the opposite sign from hi (the rounding residual). `cacosdd` delegates
to `casindd`, so the same fix covers it.

Both fixes are idempotent on inputs with `Im ≠ 0`, so
`multifloat_fuzz`'s random-input precision tables are byte-identical
before and after (see `after.txt`).

## Rebuilding the probes

```
cd tier4
g++-15 -std=c++17 -O2 catanh_probe.cc   -L../build/src -lmultifloats -o catanh_probe
g++-15 -std=c++17 -O2 casin_probe.cc    -L../build/src -lmultifloats -o casin_probe
g++-15 -std=c++17 -O2 casin_trace.cc    -L../build/src -lmultifloats -o casin_trace
```
