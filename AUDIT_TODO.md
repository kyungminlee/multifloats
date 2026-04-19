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

- [x] **13. Unit tests for log2/log1p/expm1/cbrt** — new
  `test_log_root_edges` in `test/test.cc` covers `log2(2^k)` identity at
  25 exponents, log1p / expm1 at tiny inputs down to 1e-20 (full DD
  precision), expm1(0)/log1p(0) exact-zero, expm1(1) (e-1 to kernel
  precision), cbrt on 9 perfect cubes + 6 non-cubes + signed zero.
- [x] **14. C-ABI correctness tests** — new `test/abi_equivalence.f90`
  asserts the three DD entry points (native operator / C wrapper
  `adddd`..`sqrtdd` / Fortran `bind(c)` reimplementation in
  `dd_bindc.f90`) agree on the same inputs: HI limb bit-exact, LO
  within 4 dp ULPs. Originally attempted as strict bit-equality but
  rejected: bindc's sqrt and the C divdd each pick slightly different
  compensated-error residuals (1 ULP in lo), all still full DD.
- [x] **15. Matmul non-square / transposed shapes** — added
  `test_matmul_shapes` in Fortran `precision.f90` covering rectangular
  (3×5·5×2), outer product (3×1·1×2), inner product (1×5·5×1),
  short-fat, tall-skinny, column-scalar, k-dominant. All against qp
  reference at full DD tolerance.
- [x] **16. Complex branch cuts** — `test_complex_branch_cuts` in
  `test/test.cc`. Fixed four real bugs found in the process:
  1. `atan2_full` used `y._limbs[0] >= 0.0` which is true for `-0.0`, so
     `clog(x + -0i)` for x<0 returned +π instead of -π; switched to
     `std::signbit(...)`.
  2. `csqrtdd` had the same `b._limbs[0] < 0.0` idiom for the imag-sign
     pick; same fix.
  3. `catanhdd(±1 + 0i)` returned NaN because the downstream DD multiply
     `0.25 · log(4/0) = 0.25 · ∞` emits `fma(0.25, ∞, −∞) = NaN` in
     the compensated-error term. Added explicit short-circuit at the
     branch-point singularity.
  4. `casindd` / `cacosdd` lost the signed-zero imag input because DD
     multiply (`2 · -0 → +0`) doesn't preserve the sign through the
     `−2ab` computation. Added a post-pass: if the computed Im(result)
     hi sign disagrees with Im(z) sign, negate the whole DD pair (both
     limbs — a limb-wise copysign would corrupt the pair since lo is
     typically the opposite sign of hi). Fuzz precision on random
     inputs is byte-identical before/after (the fix is idempotent for
     nonzero Im).
- [x] **17. Huge-argument trig** — `test_huge_argument_trig` in
  `test/test.cc` checks sin(2π·k) for k=2^N (N=0..40), sin(π/2+2π·k),
  exact-zero / ±1 at integer arguments of sinpi/cospi.
- [x] **18. Tolerance sensitivity sweep** — added a "tolerance-ratchet"
  post-pass in `main()` of `cpp_test` that pins observed max_rel for
  each test category into a [1/20×, 20×] band around an expected value.
  Fails in BOTH directions: too high = precision regression, too low =
  silent improvement (signal to tighten pin). Caught a real case during
  Tier 4 work: cx-branch-cuts pin was set at 2e-29 before running the
  test, the actual observation was 6.5e-33, ratchet flagged "BETTER"
  and the pin was tightened accordingly.
- [x] **19. Fuzz seed determinism test** — two new ctest cases
  (`cpp_fuzz_determinism`, `fortran_fuzz_determinism`) run the fuzz
  twice with the same seed (2000 iters for C++; the full 1M for
  Fortran) and diff the outputs. Any non-deterministic state (rogue
  static, uninitialized scratch, entropy-seeded RNG) would produce
  different output bits and fail CI.

## Tier 5 — Maintainability & readability

- [ ] **20. CHANGELOG.md + migration guide** for 44b3a64 prefix rename. **S**
- [x] **21. Unify test helpers** — consolidated `to_q` / `from_q` /
  `q_rel_err` / `qstr` / `q_isnan` / `q_isfinite` into
  `test/test_common.hh` (namespace `multifloats_test`). `test.cc`,
  `fuzz.cc`, `bench.cc` now include the shared header; `to_mf2`
  renamed to `from_q` for consistency. All 8 ctest targets pass.
- [x] **22. `dd_constants.hh` TOC + per-block provenance** — added
  `section(name, source, notes)` helper to
  `scripts/gen_constants.py`. Generated header now opens with a
  13-entry TOC and each block is framed by
  `=== SECTION: <NAME> ===` with explicit source citation
  (libquadmath erfq.c / atanq.c / asinq.c / lgammaq.c / j0q.c /
  j1q.c; in-house Taylor / Remez; mpmath for reference values).
  1736 constants verified at max conv err 5.8e-33.
- [x] **23. Inline provenance comments** — added formula + source
  one-liners at the polynomial sites that lacked them: erf_TN1/TD1,
  erf_TN2/TD2, erfc sub-interval, erfc asymptotic, bessel_j0/j1
  power-series and Hankel asymptotic branches. atan_P/Q and asin
  regions already had adequate inline formulas.
- [ ] **24. README minimal examples** — 5-line Fortran, C, C++ snippets. **S**
- [x] **25. Delete `dd_constants.f90.inc`** — was a 7-line
  empty-comment placeholder. The generator script (`F90_KEPT_NAMES =
  set()`) explicitly documented that the Fortran module no longer
  needs any generated constants (all `DD_*` named constants are
  materialized via the `DD_CONST` macro directly in the fypp, and
  every math routine delegates to C++). Removed the file, the
  `include 'dd_constants.f90.inc'` in `fsrc/multifloats.fypp`, the
  `DEPENDS` entry in `fsrc/CMakeLists.txt`, and the Fortran-output
  path (`F90_OUT`, `write_f90`) from `scripts/gen_constants.py`. Any
  future need for Fortran-side generated constants can be rehydrated
  in a follow-up commit.
- [x] **26. Split `multifloats_math.cc`** — carved 2358-line file into
  9 topical `.inc` files that the single compiled TU pulls in via
  `#include`. Preserves all cross-kernel inlining (critical for
  composed ops like `ctandd` / `casindd` / `ccoshdd` that chain
  through `sincos_full` / `sinhcosh_full` / `csqrtdd` / `clogdd`)
  even without `-flto`. Perf bit-identical to pre-split; all 8 ctest
  targets pass. Files: `multifloats_math_{exp_log, trig, hyp,
  inv_trig, special, bessel, matmul, abi_scalar, abi_complex}.inc`.
- [ ] **27. Categorize ctest names** (`precision_*`, `fuzz_*`, `perf_*`). **S**
- [ ] **28. Clean up `work-gemini/`** + `external/` unused samples. **S**
- [ ] **29. PR CI workflow** (currently only tag builds). **S**
