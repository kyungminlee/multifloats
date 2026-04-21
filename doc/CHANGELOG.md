# Changelog

All user-visible changes to `multifloats` are recorded here. Format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
the project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Dates are ISO-8601 UTC.

## [Unreleased]

### Breaking

- **Fortran-surface prefix rename** (commit `44b3a64`, 2026-04-17): the
  remaining `mf_` / `cx_` / `MF_` tokens in the Fortran module were
  normalized to `dd_` / `cdd_` / `DD_` to match the already-renamed C
  ABI (`sindd`, `matmuldd_mm`, …). The type `complex128x2` was also
  renamed to `complex64x2` so both types share the "number of dp
  limbs" naming convention (`float64x2` = 2 dp limbs,
  `complex64x2` = two `float64x2` components).

  This breaks any caller that imported these by name. The two derived
  type names `float64x2` / `complex64x2` and the generic-intrinsic
  names (`sin`, `matmul`, `sum`, `+`, `=`, …) are unaffected, so code
  that uses only generics — the vast majority of user code — needs no
  changes. Only callers that reached for the specific module procedure
  or the public constants need to adapt.

  **Migration table (search/replace, word-boundary):**

  | Old                       | New                        |
  | ------------------------- | -------------------------- |
  | `mf_*`                    | `dd_*`                     |
  | `cx_*`                    | `cdd_*`                    |
  | `MF_*`                    | `DD_*`                     |
  | `complex128x2`            | `complex64x2`              |
  | `mf_real`                 | `dd_real`                  |
  | `mf_to_double`            | `dd_to_double`             |
  | `mf_set_fma_renorm_interval` | `dd_set_fma_renorm_interval` |
  | `MF_ZERO` / `MF_ONE` / `MF_TWO` / `MF_HALF` / `MF_EIGHT` | `DD_ZERO` / … |
  | `MF_SAFMIN` / `MF_SAFMAX` | `DD_SAFMIN` / `DD_SAFMAX`  |
  | `MF_TSML` / `MF_TBIG` / `MF_SSML` / `MF_SBIG` | `DD_TSML` / … |
  | `MF_RTMIN` / `MF_RTMAX`   | `DD_RTMIN` / `DD_RTMAX`    |
  | `MF_FMA_RENORM_INTERVAL`  | `DD_FMA_RENORM_INTERVAL`   |

  A single `sed -i -E 's/\b(mf_|cx_|MF_)/\1/g'` with the obvious
  replacements is sufficient for most code bases. The rename is a
  textual change only; kernel behavior is bit-identical.

### Added

- C++ I/O: `to_string(float64x2, int precision=32)` and
  `operator<<(std::ostream&, float64x2 const&)` (audit #5).
- Fortran `sincos` / `sinhcosh` generic subroutines delegating to
  `sincosdd` / `sinhcoshdd` (audit #6).
- Complex DD transcendentals `clog2dd`, `clog1pdd`, `cexpm1dd`,
  `csinpidd`, `ccospidd` (audit #7) and matching Fortran generics
  `log2`, `log10`, `log1p`, `expm1`, `sinpi`, `cospi`.
- Ten new regression tests spanning log2/log1p/expm1/cbrt edges
  (#13), C-ABI equivalence (#14), non-square matmul shapes (#15),
  complex branch cuts (#16), huge-argument trig (#17), a tolerance
  ratchet (#18), and fuzz seed-determinism (#19).
- Minimal Fortran / C / C++ example snippets in `README.md` (#24).
- `CHANGELOG.md` with migration guide (this file, audit #20).
- TOC + per-section provenance in `src/dd_constants.hh`, emitted by
  `scripts/gen_constants.py` (#22). Coefficient groups now carry
  explicit citations to their libquadmath source (erfq / atanq /
  asinq / lgammaq / j0q / j1q), mpmath, or in-house Taylor/Remez.
- Inline formula comments for erf, erfc, and bessel polynomial
  evaluation sites (#23).
- `test/test_common.hh` unifying `to_q` / `from_q` / `q_rel_err` /
  `qstr` test helpers across `test.cc` / `fuzz.cc` / `bench.cc`
  (#21).
- `bench/` harness with `run_benchmarks.py` and auto-generated
  `BENCHMARK.md`.

### Changed

- **Precision**: division preserves non-finite lo limb (audit #1);
  `atan2dd(±0,±0)` picks the quadrant from low limbs when both hi's
  are zero (#2); Bessel Miller recurrence start index uses a DD-scale
  threshold (1e32 + `k < 1000` cap) instead of the old dp threshold
  (#4).
- **Performance**: Estrin polynomial evaluation replaces Horner at 7
  hot sites; measured speedups 1.78–1.88× on sin/cos/exp/tan and
  1.83× on hyperbolic (#11). Precision stays under 1 DD ulp.
- **Matmul**: `DD_FMA_RENORM_INTERVAL` default of 8 is now documented
  with an empirical tuning table (#12). Periodic renormalization
  gives ~50× precision gain at k=65k for only 3–4% overhead over
  ri=0.
- **README**: new "Error handling" section documenting the strict
  NaN-in-NaN-out policy (no errno, no fenv, no exceptions) (#9) and
  a "Matmul API and GEMM relationship" section (#8).
- **Test naming**: ctest `NAME` fields prefixed by category
  (`precision_*`, `fuzz_*`) so `ctest -R precision_` selects one
  group at a time (#27).
- **Source layout**: `multifloats_math.cc` was split into 9 topical
  `.inc` files pulled into a single compiled TU, preserving all
  cross-kernel inlining without requiring `-flto` (#26).

### Fixed

- `catandd` imaginary part lost ~7 decimals when `|b|` is small and
  `|a|` moderate (num/den ≈ 1, cancellation in `0.25·log(num/den)`):
  max_rel dropped from 3.3e-25 to 9.5e-32, mean_rel from 4.4e-28 to
  1.2e-32. Now delegates to `catanhdd` via `catan(z) = −i·catanh(iz)`,
  inheriting the Kahan `log1p(4a/den)` form and the unit-circle
  `dd_x2y2m1` denominator. Bench cost: ~6% slower (0.0016 → 0.0017 s
  / 4096 ops).
- `catanhdd(±1 + 0i)` returned NaN; now short-circuits at the
  branch-point singularity to `(copysign(inf, re), +0)` (#16.3).
- `casindd` / `cacosdd` lost the signed-zero imag-part sign because
  `2·(-0) = +0` in DD multiply; post-pass re-applies the input's sign
  to both limbs when they disagree (#16.4).
- `atan2_full` and `csqrtdd` used `x >= 0.0` which is `true` for
  `-0.0`; switched to `std::signbit` so `clog(x + -0i)` now returns
  the expected `-π` on the negative real axis (#16.1, #16.2).

### Removed

- Empty placeholder `src/dd_constants.f90.inc` (commit `c5fe91d`).
  The Fortran module materializes its named constants directly via
  the `DD_CONST` macro; the include file was a no-op left over from
  an earlier split.
- Scratch `work-gemini/` directory and unused `external/` samples
  (`erf.cpp`, `float64x2-sample.{cpp,hpp}`, `v3.12.1.tar.gz`) (#28).

## [0.1.0] — 2026-04-13

Initial tagged release. See the release workflow at
`.github/workflows/release.yml` for the build matrix (GCC 12–14,
LLVM Flang 19–21, Intel oneAPI 2024–2025). The tag predates the
Fortran-surface prefix rename: callers pinned to `v0.1.0` see the
old `mf_*` / `cx_*` names.
