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

- [ ] **5. C++ I/O** — `operator<<(ostream&, float64x2)` + `to_string`. `src/multifloats.hh`. **M**
- [ ] **6. Fortran `sincos` / `sinhcosh`** — expose fused pairs already in C ABI. **S**
- [ ] **7. Complex DD transcendentals** — `clog2dd`, `clog10dd`, `clog1pdd`,
  `cexpm1dd`, `csinpidd`, `ccospidd` in C ABI + Fortran. **M**
- [ ] **8. Matmul transA/transB/alpha/beta** — at minimum document current
  behavior; ideally add GEMM-style flags. **L**
- [ ] **9. Error-handling policy in README** — NaN in, NaN out; no errno; no
  signalling. **S**

## Tier 3 — Performance

- [ ] **10. Complex Karatsuba 3-mul** — replace 4-mul form in DD complex
  multiply. `fsrc/multifloats.fypp:2823-2824` + C++ mirror. **M**
- [ ] **11. Estrin polynomial evaluation** — alongside Horner for 6+ term
  polynomials (exp2/log2 hotspots). **M**
- [ ] **12. Renorm-interval tuning doc or auto-tune** — `k`-aware default
  for matmul. `src/multifloats_math.cc:1664` + README. **S**

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
