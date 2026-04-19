# Audit TODO

Consolidated action plan from six parallel subagent audits (performance,
precision, API completeness, maintainability, testing, readability) run on
2026-04-18 against commit `44b3a64` (post Fortran prefix rename).

Every item below is sized to fit in a single PR. Effort tags: **S** = hours,
**M** = a day, **L** = multi-day.

---

## Tier 1 — Correctness (DONE 2026-04-18)

- [x] **1. Division preserves non-finite lo limb** — `src/multifloats.hh:227-234`
  Fixed: the `!isfinite(q1)` branch now mirrors `q1` into the lo limb so
  both limbs are non-finite together. Regression test
  `test_division_nonfinite` added in `test/test.cc`.
- [x] **2. `atan2dd(±0,±0)` sign from low limbs** — `src/multifloats_math.cc:705`
  Fixed: when both hi's are zero, pass the lo limbs to `std::atan2` so the
  effective sign (carried by lo when hi==0) determines the quadrant.
  Regression test `test_atan2_signed_zero` added.
- [x] **3. `csqrtdd` signed-zero at `a = -0`** — FALSE POSITIVE.
  C99 G.6.4.2 specifies `csqrt(±0 + 0i) → +0 + 0i`: the real part is
  always +0, not preserved. Current code is compliant; added a guard
  test `test_csqrt_zero_branch` so a future "fix" can't silently
  regress the behavior.
- [x] **4. Bessel Miller start-index threshold** — `src/multifloats_math.cc:1507`
  Fixed: threshold raised from `1e17` (~2^56, double-precision target) to
  `1e32` (~2^106, DD target) plus a `k < 1000` safety cap. Measured
  max_rel on the Miller path (n=10..40, x=1..10) drops from ~1e-17 to
  ~9e-32, near the DD floor. New test `test_bessel_jn_miller_precision`.

Perf sentinels before/after: `div` 0.0013s, `atan2` 0.0072→0.0073s,
`sqrt` 0.0032→0.0031s, `bessel_jn(3,·)` 0.0073s — all within noise.

## Tier 2 — API completeness (blocks 1.0)

- [x] **5. C++ I/O** — `to_string(float64x2, int precision=32)` + `operator<<`
  in `multifloats` namespace (inline in header; archive exports only
  extern "C" so C++ helpers must be header-resident). Scientific notation
  with up to 34 digits; round-half-to-even with 2-digit guard. Regression
  test `test_io_to_string_and_stream` covers nan/inf/signed-zero formats,
  round-trip vs `__float128`, precision clamping, carry rollover.
- [x] **6. Fortran `sincos` / `sinhcosh`** — exposed as pure elemental
  subroutines delegating to `sincosdd` / `sinhcoshdd`. Added bit-equal vs
  separate-call regression test (`test_sincos_sinhcosh` in
  `test/precision.f90`, 5 samples incl. 0, negative, large).
- [x] **7. Complex DD transcendentals** — added `clog2dd`, `clog1pdd`,
  `cexpm1dd`, `csinpidd`, `ccospidd` in C ABI (clog10dd already existed);
  all six exposed as Fortran generics (`log2`, `log10`, `log1p`, `expm1`,
  `sinpi`, `cospi`). Also filled pre-existing real gaps (log2, log1p,
  expm1 previously had `dd_*_full` kernels but no generic interface).
  New C++ test `test_complex_new_transcendentals` (max_rel 1.1e-31 — DD
  ulp) + Fortran `test_new_generic_intrinsics` covering all 9 identities.
- [x] **8. Matmul transA/transB/alpha/beta** — documented the current
  non-GEMM semantics explicitly in the C-ABI header and a new README
  "Matmul API and GEMM relationship" section (no trans flags, no
  alpha/beta, no LDA; contiguous column-major only). GEMM-style flag
  extension is left as a future item — requires new panel dispatchers
  for transposed shapes.
- [x] **9. Error-handling policy in README** — Added dedicated "Error
  handling" section documenting NaN-in-NaN-out, no errno, no fenv, no
  exceptions, no signalling NaN, no input validation.

## Tier 3 — Performance

- [x] **10. Complex Karatsuba 3-mul** — REJECTED via experiment.
  Karatsuba trades 1 DD mul (~11 dp ops) for 3 extra DD adds (~21 dp
  ops) — a net op-count loss at DD granularity. A/B probe on 200k
  random inputs: 4-mul 15.6 ns/op vs Karatsuba 35.3 ns/op (2.26×
  SLOWER). Precision also degrades: Im max_rel goes from 3.9e-32 to
  6.8e-32 (1.7× worse), with catastrophic cancellation on the classic
  witness a=(1,ε), b=(-1,ε) — true Im is 0 exactly, Karatsuba returns
  −ε². Added regression test `test_complex_mul_cancellation` that
  locks in the 4-mul semantics so a future contributor can't silently
  regress this.
- [x] **11. Estrin polynomial evaluation** — extended `neval` with
  Estrin cases 12–15 (`x⁸` recursion) and swapped 7 hot Horner sites:
  `exp2` (14 coefs), `sin`/`cos` kernels (13), `asinh`/`atanh` Taylor
  (15), `log2_wide` (9), `sinh` Taylor (9, ×2 callers).
  Measured speedups (cpp_bench): `exp` / `exp2` 1.88×, `sin` / `cos`
  / `tan` 1.78×, `sinh` / `cosh` / `tanh` 1.83×, `pow` 1.35×, `expm1`
  1.48× (indirect). Precision: `sin` max_rel 3.4e-32 → 5.0e-32, `cos`
  3.0e-32 → 3.8e-32 — both still < 1 DD ULP; all other ops identical.
  `expm1_taylor` (n=24) and `log1p_taylor` (n=17) left for a future
  pass (would require `x¹⁶` Estrin level).
- [x] **12. Renorm-interval tuning doc** — ran an empirical sweep over
  (k ∈ {16..65k}, ri ∈ {0, 4, 8, 16, 32, 64}) for 4×k·k×4 matmul.
  Default ri=8 confirmed as near-optimal (~50× precision gain over
  ri=0 at k=65k, only ~3-4% slower). Documented the tradeoff table in
  the Fortran fypp comment and the README matmul section with tuning
  guidance by k range.

## Tier 4 — Test coverage

- [ ] **13. Unit tests for log2/log1p/expm1/cbrt** with curated inputs. **S**
- [ ] **14. C-ABI correctness tests** — assert `bench_abi` / `dd_bindc`
  kernels bit-match C++. **S**
- [ ] **15. Matmul non-square / transposed shapes** — 3×5·5×2, outer product,
  LDA≠rows. **M**
- [ ] **16. Complex branch cuts** — `clog`, `csqrt`, `casin`, `cacos` vs
  libquadmath on negative-real / imaginary axes. **M**
- [ ] **17. Huge-argument trig** — `sin(2π·2^40)` style range-reduction
  sanity. **S**
- [ ] **18. Tolerance sensitivity sweep** — fail when observed error is far
  below hardcoded tolerance (detects drift in both directions). **M**
- [ ] **19. Fuzz seed determinism test** — same seed, same failure
  signature. **S**

## Tier 5 — Maintainability & readability

- [ ] **20. CHANGELOG.md + migration guide** for 44b3a64 prefix rename. **S**
- [ ] **21. Unify test helpers** — `to_q`, `from_mf2` into
  `test/test_common.hh`. **S**
- [ ] **22. `dd_constants.hh` TOC + per-block provenance** (Taylor / Remez /
  libquadmath citations) via `scripts/gen_constants.py`. **S**
- [ ] **23. Inline provenance comments** for polynomial evals in
  `multifloats_math.cc` (atan_P/Q, erf rationals). **S**
- [ ] **24. README minimal examples** — 5-line Fortran, C, C++ snippets. **S**
- [ ] **25. Delete or populate `dd_constants.f90.inc`** (7-line placeholder
  today). **S**
- [ ] **26. Split `multifloats_math.cc`** into `_exp_log.cc`, `_trig.cc`,
  `_special.cc`, or add a TOC header. **M**
- [ ] **27. Categorize ctest names** (`precision_*`, `fuzz_*`, `perf_*`). **S**
- [ ] **28. Clean up `work-gemini/`** + `external/` unused samples. **S**
- [ ] **29. PR CI workflow** (currently only tag builds). **S**
