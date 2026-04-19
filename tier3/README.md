# tier3 — Performance-audit scratch folder

Raw measurements and reproducer sources from the Tier 3 performance
audit (commit range around `d0a8da4`..). Kept in-tree so the same
experiments are not re-run from scratch if someone revisits these
tradeoffs.

## Files

| File                   | What it is                                                    |
| ---------------------- | ------------------------------------------------------------- |
| `baseline.txt`         | Pre-Tier-3 cpp_bench + fortran_bench + fortran_fuzz dump.     |
| `step1.txt`            | After swapping `log2_wide` / `sinh_taylor` Horner → Estrin.   |
| `step2a.txt`           | After swapping `exp2` (14-coef) → Estrin case 13.             |
| `step2b.txt`           | After swapping `sin` / `cos` kernels → Estrin case 12.        |
| `step2c.txt`           | After swapping `asinh` / `atanh` Taylor → Estrin case 14.     |
| `step2_final.txt`      | Cumulative cpp_bench after all Estrin swaps.                  |
| `step2_summary.txt`    | Human summary table: per-op baseline → final speedup.         |
| `fuzz_final.txt`       | Post-Tier-3 fortran_fuzz precision report.                    |
| `karatsuba_probe.cc`   | A/B test: 4-mul vs Karatsuba 3-mul for DD complex multiply.   |
| `renorm_sweep.cc`      | Sweep over (k, ri) for matmul to characterize the tradeoff.   |
| `renorm_results.txt`   | Output of `renorm_sweep` — table used in README + fypp doc.   |

## Findings

### #10 Karatsuba complex multiply — REJECTED

`karatsuba_probe.cc` shows Karatsuba is **2.26× slower** for DD (extra
adds outweigh the saved mul) and **less precise** (catastrophic
cancellation: true Im = 0 for `(1,ε)·(-1,ε)`, Karatsuba returns
`-ε²`). Current 4-mul form is kept and now has a regression test
(`test_complex_mul_cancellation` in `test/test.cc`) plus inline warnings
in `fsrc/multifloats.fypp` (`cdd_arith_body` mul branch).

### #11 Estrin polynomial evaluation — ADOPTED

Swapped 7 hot Horner sites in `multifloats_math.cc` to Estrin via
`neval` (extended with cases 12–15). Measured speedups:

| op                      | baseline | after  | speedup |
| ----------------------- | -------- | ------ | ------- |
| `exp` / `exp2`          | 0.0077   | 0.0041 | 1.88×   |
| `sin` / `cos` / `tan`   | 0.0101   | 0.0056 | 1.78×   |
| `sinh` / `cosh` / `tanh`| 0.0156   | 0.0085 | 1.83×   |
| `pow`                   | 0.0151   | 0.0112 | 1.35×   |
| `expm1`                 | 0.0092   | 0.0062 | 1.48×   |

Precision: `sin` max_rel grew 3.4e-32 → 5.0e-32, `cos` 3.0e-32 →
3.8e-32; both still inside 1 DD ULP. All other ops identical.

`expm1_taylor` (degree 24) and `log1p_taylor` (degree 17) were not
swapped — they need an `x^16` Estrin level which is a substantial code
bulge for a ~1.5× gain on small-|x| branches only.

### #12 Renorm-interval tuning — DOCUMENTED

`renorm_sweep.cc` / `renorm_results.txt` sweep a 4×k·k×4 matmul over
`k ∈ {16..65536}` and `ri ∈ {0, 4, 8, 16, 32, 64}`. Default `ri = 8`
is near-optimal: ~50× precision win over `ri = 0` at `k = 65k` for
only ~3-4% runtime overhead. `ri = 4` costs 40–70% more runtime at
large k without a commensurate precision gain. The table is
reproduced in `README.md` (Matmul API section) and in the
`DD_FMA_RENORM_INTERVAL` module-level comment in
`fsrc/multifloats.fypp`.

## Rebuilding the probes

```
cd tier3
g++-15 -std=c++17 -O3 -march=native -fext-numeric-literals \
    karatsuba_probe.cc -lquadmath -o karatsuba_probe
g++-15 -std=c++17 -O3 -march=native -fext-numeric-literals \
    renorm_sweep.cc \
    -L../build/src -lmultifloats -lquadmath \
    -o renorm_sweep
```

(Needs `cmake --build ../build` to have produced `libmultifloats.a`.)
