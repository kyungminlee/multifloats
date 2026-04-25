# Changelog

All user-visible changes to `multifloats` are recorded here. Format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
the project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Dates are ISO-8601 UTC.

## [Unreleased]

### Changed

- **Default CMake build is now library-only.** A plain
  `cmake -S . -B build` produces just `multifloats` and
  `multifloatsf`. Tests, benchmarks, the `blas-multifloat` smoke
  target, and the MPFR / Boost comparison harnesses are all opt-in.
  Use `-DBUILD_TESTING=ON` for the test suite (which also brings in
  `blas-multifloat`), `-DMULTIFLOATS_BUILD_BENCH=ON` for the
  microbenchmarks, `-DBUILD_MPFR_TESTS=ON` for the 3-way MPFR
  precision check, and `-DMULTIFLOATS_BUILD_BOOST_COMPARE=ON` for
  the `boost::multiprecision::cpp_double_double` comparison
  (fetches Boost ≥ 1.89 via FetchContent). See README "Build" for
  the full table. Consumers that previously relied on the
  default-ON `BUILD_TESTING` (via `include(CTest)`) or on the
  default-on-top-level `MULTIFLOATS_BUILD_BENCH` need to opt in
  explicitly.

### Added

- `boost::multiprecision::cpp_double_double` precision and speed
  comparison harnesses (`boost_dd_fuzz`, `boost_dd_bench`) plus a
  `bjn_probe` regime sweep, gated behind
  `-DMULTIFLOATS_BUILD_BOOST_COMPARE=ON`. Documented in
  `doc/developer/BOOST_COMPARISON.md`.
- Worst-case input printer in both `cpp_fuzz` and `boost_dd_fuzz`
  precision reports — alongside `max_rel` / `mean_rel`, each op now
  prints the `(i1, i2, ref, got)` sample that produced its worst
  rel-err. Makes stochastic precision anomalies debuggable from a
  single run rather than requiring a custom probe rebuild.
- Fuzz coverage for `cbrt`, `fma_cxx`, `lerp`, `modulo`,
  `remainder`, `remquo`, `fpclassify`, `nextafter`, `nexttoward`,
  and `nan(tag)`.

### Fixed

- `multifloats::sqrt` Karp–Markstein refinement step now uses
  `two_prod(c, c)` (FMA-based exact scalar product) instead of a
  full DD multiply for the residual `x − c²`. ~2.2× faster, max_rel
  3.5e-32 vs 4.1e-32 prior.

## [0.3.4] — 2026-04-24

### Changed

- Bump CI action versions to the Node-24-native majors:
  `actions/checkout@v4` → `@v6`,
  `actions/upload-artifact@v4` → `@v7`,
  `actions/download-artifact@v4` → `@v8`.
  Supersedes the v0.3.3 `FORCE_JAVASCRIPT_ACTIONS_TO_NODE24` env opt-in,
  which is now redundant and removed. No workflow logic changes; the
  basic `name`/`path` inputs we use are stable across these majors.

## [0.3.3] — 2026-04-24

### Fixed

- Silence gfortran's `-Wlto-type-mismatch` on the Fortran ↔ C bind(c)
  boundary (≈1500 warnings per build). The Fortran bridge types
  `float64x2_t` / `complex64x2_t` are still layout-compatible with
  `multifloats::float64x2` / `std::complex<multifloats::float64x2>`
  (pinned by header static_asserts + verified by
  `precision_abi_equivalence` and `crosscheck_cpp_fortran`); only the
  type *names* differ, which gfortran's LTO type-checker flags.
  Scoped `-Wno-lto-type-mismatch` to GNU + LTO builds on the
  `multifloatsf` target.

### Changed

- CI / release workflows opt in to the Node.js 24 runtime
  (`FORCE_JAVASCRIPT_ACTIONS_TO_NODE24: true`) to silence the
  `actions/checkout@v4` / `actions/upload-artifact@v4` /
  `actions/download-artifact@v4` Node 20 deprecation warnings. Safe
  to drop after GitHub's 2026-06-02 force-default-to-Node-24
  transition.

## [0.3.2] — 2026-04-24

### Fixed

- Release-CI install failure on every non-GCC matrix entry (Intel
  oneAPI, LLVM Flang). The dual-variant `.f90` install in v0.3.0/0.3.1
  expected both `multifloats-quad.f90` and `multifloats-noquad.f90` to
  exist, but the non-GCC entries used `cmake --build --target
  multifloats --target multifloatsf`, which bypasses ALL-linked custom
  targets and only generated the variant the build host's REAL(16)
  probe had picked. `add_dependencies(multifloatsf
  multifloatsf_sources)` wires both fypp expansions as a side effect
  of building the library target.

## [0.3.1] — 2026-04-24

### Added

- **Dual-mode Fortran package distribution**: `find_package(multifloatsf)`
  now prefers a precompiled compiler-tagged archive + `.mod` when one
  matches the consumer's Fortran compiler (`multifloatsfTargets-<tag>.cmake`
  under `<prefix>/lib/cmake/multifloatsf/`), and falls back to the
  source-build path otherwise. Release tarballs now ship one tagged
  precompiled variant per build-matrix entry alongside the pre-expanded
  `.f90` source — consumers whose compiler exactly matches skip the
  module recompile entirely.
- `-DMULTIFLOATSF_INSTALL_PRECOMPILED_MOD=ON` wired into all release CI
  configure steps (`build-lto`, `build-lto-mac`) so the release artifact
  set carries per-compiler precompiled bundles.

### Changed

- `multifloatsfConfig.cmake` selector reworked from "always source" to
  "precompiled preferred, source fallback". The `@FC_DETECT_CODE@`
  placeholder in the template is now expanded at library-configure time
  from the `_FC_DETECT_CODE` snippet in `cmake/FortranCompiler.cmake`, so
  producer and consumer compute identical compiler tags.

## [0.3.0] — 2026-04-24

### Breaking

- **C++ type unification** (commit `6f1c472`): the parallel
  `float64x2_t` (C POD) and `multifloats::float64x2` (C++ class) types
  were collapsed into one `struct float64x2 { double limbs[2]; }` used
  by both languages. `_limbs` → `limbs`; `.hi` / `.lo` field accesses
  on the old C POD become `.limbs[0]` / `.limbs[1]`. Similarly
  `complex64x2_t` is gone — C sees `struct complex64x2 { float64x2 re, im; }`,
  C++ sees `using complex64x2 = std::complex<float64x2>;`. Fortran
  `bind(c)` interop types (`type dd_c`) switched from `hi, lo` fields
  to `limbs(2)`. Linker symbol set, calling convention, and on-disk
  layout are all unchanged; `MULTIFLOATS_ABI_VERSION` stays at 2.
- **C++ namespace move** (commit `6f1c472`): the public header no
  longer injects `float64x2` / `complex64x2` into the global namespace
  via `using multifloats::…`. External C++ consumers must spell the
  types as `multifloats::float64x2` / `multifloats::complex64x2`, or
  add their own `using multifloats::float64x2;` at TU scope. The
  extern "C" `*dd` prototypes now sit inside `namespace multifloats`
  — function calls resolve via ADL when the argument type is a
  `multifloats::float64x2` (common case, no user action) or via
  explicit qualification otherwise.
- **Fortran module distribution model** (commit `22ea804`): the
  default install no longer ships a compiler-tagged `.mod` + archive.
  Instead the pre-expanded `multifloats-quad.f90` and
  `multifloats-noquad.f90` sources are installed under
  `<prefix>/share/multifloatsf/src/` and a
  `multifloatsfConfig.cmake` lets consumers `find_package(multifloatsf)`
  + compile the module with their own Fortran compiler. Cross-compiler
  consumers are now viable (gfortran-13 library ↔ gfortran-15
  consumer, flang ↔ ifx, etc.). Precompiled `.mod` install is retained
  as an opt-in flag (`-DMULTIFLOATSF_INSTALL_PRECOMPILED_MOD=ON`) for
  distro packagers.

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

- **Fortran source-distribution package** (commit `22ea804`):
  `cmake/multifloatsfConfig.cmake.in`, both `multifloats-{quad,noquad}.f90`
  variants, and `test/consumer-fortran/` + `scripts/check_consumer_fortran.sh`
  as a new `consumer_fortran_smoke` ctest driving install → configure →
  build → run of a find_package-based consumer.
- **Layout static_asserts** on `multifloats::float64x2` (commit
  `6f1c472`): `std::is_standard_layout_v` + `std::is_trivially_copyable_v`,
  replacing the pre-refactor tautological `sizeof(a)==sizeof(a)`
  boundary check.
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

- **`isinf(float64x2)`** now scans both limbs (commit `6f1c472`),
  matching `isnan` / `isfinite`. A non-canonical DD with finite `hi`
  and infinite `lo` (constructible via the C ABI) previously
  classified as not-finite, not-inf, not-nan.
- `catandd` imaginary part lost ~7 decimals when `|b|` is small and
  `|a|` moderate (num/den ≈ 1, cancellation in `0.25·log(num/den)`):
  max_rel dropped from 3.3e-25 to 9.5e-32, mean_rel from 4.4e-28 to
  1.2e-32. Now delegates to `catanhdd` via `catan(z) = −i·catanh(iz)`,
  inheriting the Kahan `log1p(4a/den)` form and the unit-circle
  `dd_x2y2m1` denominator. Bench cost: ~6% slower (0.0016 → 0.0017 s
  / 4096 ops).
- `casindd` used the textbook `−i·log(iz + sqrt(1 − z²))` which lost
  ~3 decimals near the real-axis branch cuts because `1 − z²` cancels
  catastrophically when `|Re z| ≈ 1`. Now delegates to `casinhdd` via
  `asin(z) = −i·asinh(iz)`, inheriting the libquadmath-ported 10-region
  branch schedule. max_rel on cdd_asin_re dropped from 7.5e-29 to
  4.0e-32, cdd_asin_im from 7.6e-29 to 4.2e-32 (≈1800× better worst
  case). `cacosdd` = π/2 − casin inherits the Im-part fix, so
  cdd_acos_im also improved from 7.6e-29 to 4.2e-32; cdd_acos_re is
  unchanged (bottlenecked by the π/2 subtraction, not casindd).
  Bench cost: none (0.0034 → 0.0033 s / 4096 ops, within noise).
  casinhdd's own signed-zero trailer reproduces the C99 G.6.2.2
  branch-cut sign fixup for `|Re z| > 1 ± 0` (#16.4 still covered).
- `cacosdd` real part lost up to ~6 decimals of relative precision
  when |z| > 1 with tiny Im — the π/2 − asin.re subtraction was capped
  at DD-absolute of π/2 (~1e-32), which inflated to ~1e-27 relative when
  the true output was ~1e-6 (e.g. z ≈ 26.6 − 4e-5 i, acos(z).re ≈
  1.67e-6 vs the reference). Now domain-splits: |z| ≤ 1 keeps the direct
  π/2 − asin form (no cancellation there); |z| > 1 routes through
  cacoshdd via `acos(z) = −i·acosh(z)` so Re(acos) = Im(acosh), putting
  the small output on a non-cancelling path. max_rel on cdd_acos_re
  dropped from 1.04e-28 to 3.65e-32 (≈2800× better worst case, both
  vs `cacosq` and vs 200-bit MPFR). cdd_acos_im stays at 4.2e-32.
  Bench cost: ~10% slower (0.0034 → 0.0037 s / 4096 ops) from the
  cacoshdd Kahan-form overhead on the |z| > 1 half of inputs. Branch-
  cut signs for `|x| > 1 ± 0` still pass via the Im-signbit-directed
  conj-reflection (test/test.cc:562-565 still green).
- `cacoshdd` real part had the mirror problem: when |z| < 1, acosh(z)
  lies near iπ/2 (acosh(0) = iπ/2 exactly), so Re(acosh) is small and
  Kahan's `2·log(sqrt((z+1)/2) + sqrt((z−1)/2))` capped Re at
  DD-absolute of π/2 / 2 ≈ 1e-32. At the fuzz worst z ≈ −2.4e-4 +
  2.5e-7 i the output .re ≈ 2.5e-7 gave ~6e-26 relative. Symmetric fix:
  domain-split on |z|². |z| < 1 now delegates to cacosdd via the mirror
  identity `acosh(z) = i · acos(z)` for Im(z) ≥ +0, or
  `conj(i · acos(conj z))` for Im(z) < 0 (incl. −0). Component-wise
  (Re acosh, Im acosh) = (−Im acos, ±Re acos), routing the small
  component onto the full-DD Im(acos) path; |z| ≥ 1 keeps Kahan's
  original form (stable near z ≈ 1). The cacos↔cacosh mutual-
  delegation cycle terminates: cacos only routes for |z| > 1, cacosh
  only for |z| < 1, so neither recurses back. max_rel on cdd_acosh_re
  dropped from 2.95e-28 to 4.19e-32 (≈7000× better worst case, both
  vs `cacoshq` and vs 200-bit MPFR). All six complex inverse
  transcendentals (asin/acos/atan/asinh/acosh/atanh) × (re/im) are
  now at full DD (~4e-32 max_rel). Bench cost: within noise (0.0036 →
  0.0035 s / 4096 ops median over 3 runs).
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
