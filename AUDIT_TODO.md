# Audit follow-ups

Remaining findings from the 2026-04-17 five-lens audit. Applied safe fixes
already; these are the ones that needed explicit review.

## Correctness

- [x] **C1. Sign-of-zero in `abs` / `copysign` / `dd_bessel_j1`.**
      Fixed: `signbit`, `abs`, `copysign` in `src/multifloats.hh` now
      walk limbs until a non-zero one is found; `dd_bessel_j1_full` uses
      the same rule; lgamma reflection routes through `multifloats::abs`.
      Added a regression test in `test/test.cc::test_classification`.

- [x] **C2. Matmul periodic renorm uses unsafe `fast_two_sum`.**
      Fixed: introduced `dd_renorm_inl` (a full two_sum) in
      `src/multifloats_math.cc`; `dd_finalize_inl`, the periodic renorm
      in `dd_gaxpy_mv_panel`, and the periodic renorm in `matmul_vmdd`
      all route through it. Removed duplicated fast_two_sum open-code.
      Fuzz suite + benchmarks unchanged (tiny constant overhead on renorm
      boundary only).

## Parity (C++ ↔ Fortran)

- [x] **P1. Fortran overflow returns `huge()` instead of `+inf`.**
      Fixed: added `ieee_positive_inf` / `ieee_negative_inf` to the
      `ieee_arithmetic` use-list and replaced all overflow / pole
      returns in `dd_exp2_full`, `dd_log2_full`, `dd_log_gamma_full`,
      `dd_gamma_full`, `bessel_yn(n, 0)`, and complex `log` at `z=0`.
      Added `test_overflow_paths()` in `test/test.f90` covering each.
      Remaining `huge(0.0_dp)` sites are intentional (`huge()` inquiry
      and `nearest()` direction sentinel).

- [x] **P2. Fortran `sinpi` / `cospi` use legacy polynomial.**
      Fixed: Fortran `dd_sinpi_full` and `dd_cospi_full` now form
      `pr = π·rx` (DD-exact) and dispatch to the existing `dd_sin_eval`
      / `dd_cos_eval`, matching the C++ path (~4e-32 vs the legacy
      ~5e-27). The `dd_sinpi_kernel` / `dd_cospi_kernel` helpers are
      gone; `gen_sinpi_coefs` / `gen_cospi_coefs` were removed from
      `scripts/gen_constants.py`; regenerated `dd_constants.*` drops the
      four dead coefficient arrays. Added `test_sinpi_cospi_precision`
      asserting `sinpi(0.25) == sqrt(2)/2` at DD precision.

- [x] **P3. Uncompensated complex reductions.**
      Fixed via **option A** (4× real matmul): `cx_dot_product` delegates
      to 4× `mf_dot_product`, and `cx_matmul_mm/mv/vm` to 4× `mf_matmul_*`
      through split (re, im) parts. Inherits the compensated DD
      accumulator and renorm interval of the real path. Each function has
      a block comment noting option B (fused compensated kernel) for
      future revisit if cx reductions become a profile hotspot.
      `test_cx_matmul_dot` covers 8 identities in `test/test.f90`.

- [x] **P4. `mf_product` / `*_dim` reductions are naive.**
      Fixed: `mf_sum_${rank}d_dim` (rank>1), `cx_sum_${rank}d`, and
      `cx_sum_${rank}d_dim` (rank>1) now inline the Neumaier-compensated
      accumulator (two accumulators for complex) that was previously only
      in the rank-N no-dim `mf_sum`. Product paths stay naive — each DD×DD
      is already ~1 DD ulp, Neumaier compensation does not apply to
      multiplication (documented inline). Final renorm upgraded from
      fast_two_sum to full two_sum. `test_compensated_reductions` covers
      sum, sum(..., dim), cx sum, and cx sum(..., dim) with a
      cancellation-prone input.

## Speed

- [x] **S1. `matmul_mmdd` re-streams A per B column.**
      Fixed: introduced `dd_gemm_panel<MR=8, NR=2>` in
      `src/multifloats_math.cc`; each p-step loads A[:,p] once into
      registers and reuses across the 2-wide B tile. Row / column tails
      (1..MR-1 rows, trailing 1-col remainder) fall back to the existing
      `dd_gaxpy_mv_panel` path. Output is bit-identical to per-column mv
      (verified via `memcmp` on 9 shapes including chunked renorm paths).
      Measured +12–40% throughput on m,n,k ∈ {8..256} shapes.

- [x] **S2. Missing `__attribute__((always_inline))` on hot kernels.**
      Fixed: added `__attribute__((always_inline))` on the three leaf
      utilities `dd_mac_inl`, `dd_renorm_inl`, `dd_finalize_inl` in
      `src/multifloats_math.cc`. Panel templates (`dd_gaxpy_mv_panel`,
      `dd_gemm_panel`) kept as plain `static inline` — annotating them
      caused a measurable 10% regression on the 8×8×1024 shape under
      clang-arm64 (large body + template instantiation interacts badly
      with the compiler's register allocator when forced inline).
      Under clang-arm64 at `-O3` the leaf hint is a no-op (flat perf);
      defensive against `-O2` / other compilers as the note said.

- [x] **S3. Bessel asymptotic doubles the sin/cos work.**
      Fixed: added `dd_sincos_eval` (shares the π/4-shift Taylor kernels
      between sin and cos on |r| ≤ π/4) and `dd_sincos_full` (fused
      range reduction) in `src/multifloats_math.cc`. Wired the four
      Bessel asymptotic sites (j0/j1/y0/y1) and `dd_tan_full` (branches
      0/2 and 1/3 now evaluate once). Measured: j0/j1/y0/y1 asymptotic
      +28%, `tan` +40% throughput on the 1024-input sweep.
      Single-eval `sin`/`cos` paths unchanged.

- [x] **S4. `dd_reduce_pi_half` re-runs full-DD subtraction three times.**
      Fixed: inlined as three (two_sum + low-limb merge + fast_two_sum)
      steps instead of three full-DD `operator-` calls; folds the
      low-limb error into one `fast_two_sum` per step (~11 FLOPs vs
      ~20). Dropped the isfinite and zero-zero guards: the x-argument
      has already been checked for non-finite in the callers, and the
      accumulator cannot hit the zero-zero branch. Measured +5% on
      `sindd`/`cosdd` and +3% on `tandd`; Bessel asymptotic gains
      a further ~1-2% on top of S3. Full fuzz suite (C++ + Fortran)
      passes without tolerance adjustment.

## Maintainability

- [x] **M1. Naming drift: `exp2_min_d` (C++) vs `exp2_min` (Fortran).**
      Fixed: `scripts/gen_constants.py` now stores the numeric clamps
      (`min`, `max`) on the `exp2_clamp` group itself; both `write_cpp`
      and `write_f90` format from those fields. C++ symbols renamed
      `exp2_min_d`/`exp2_max_d` → `exp2_min`/`exp2_max` to match
      Fortran. Regenerated `dd_constants.{hh,f90.inc}`. Only callsite
      (`dd_exp2_full` in `src/multifloats_math.cc`) updated. All tests
      pass (C++ + Fortran + fuzz).

- [x] **M2. Fortran inlines Taylor coefficients already emitted.**
      Fixed: `dd_sinh_full`, `dd_asinh_full`, `dd_atanh_full` small-|x|
      branches now call `dd_horner(x*x, {sinh,asinh,atanh}_taylor_hi,
      *_taylor_lo)` using the generated arrays (removed ~93 lines of
      duplicated coefficient literals + manual Horner expansion).
      Benchmarks unchanged to within noise
      (sinh 0.0447→0.0443s, asinh 0.0639→0.0646s, atanh 0.0503→0.0498s
      on 40960 random inputs). All tests pass.

- [x] **M3. Arithmetic-dispatch mega-template in fypp.**
      Fixed: extracted five `#:def` macros at the top of the arithmetic
      block — `unpack_operand`, `dd_arith_addsub_body`,
      `dd_arith_mul_body`, `dd_arith_div_body`, `cx_arith_body`. The
      `n1_op_n2` function now just declares locals, unpacks operands,
      negates for `sub`, and dispatches to one macro. Generated
      Fortran is bit-identical (diff shows only relocated comments /
      blank lines, zero code changes across all 96 instantiations).
      All tests pass.

- [~] **M4. 8-way erfc dispatch duplicated.** **Obsolete:** erfc now
      lives in C++ only; the Fortran side picked up `erfcdd`/`erfcxdd`
      through `C_DELEGATE_UNARY_MAP` (see M3 refactor). No Fortran
      dispatch remains to deduplicate.

- [x] **M5. Section banners.**
      Fixed: `fsrc/multifloats.fypp` now uses a clean two-level scheme —
      `! ====` super-banners mark the three structural divisions (PUBLIC
      INTERFACE, MODULE DATA/CONSTANTS/C-ABI INTERFACES, PROCEDURE
      IMPLEMENTATIONS), and `! ----` sub-banners mark each topical
      section within them. Replaced 103 `! ====` lines with `! ----`,
      removed the stray duplicate banner before `dd_constants.f90.inc`,
      and inserted the three super-banner trios. No behavioral changes;
      full test suite passes.

## API / ABI

- [x] **A1. No ABI version macro.**
      Fixed: `#define MULTIFLOATS_ABI_VERSION 1` added near the top of
      `src/multifloats_c.h`; `VERSION 1.0.0` / `SOVERSION 1` set on the
      library target so future `BUILD_SHARED_LIBS=ON` builds record the
      soname automatically.

- [x] **A2. `VISIBILITY_INLINES_HIDDEN` not set.**
      Fixed: `VISIBILITY_INLINES_HIDDEN ON` is now on the library
      target's property list alongside `CXX_VISIBILITY_PRESET hidden`.

- [x] **A3. `localize_symbols.sh.in` Darwin/Linux asymmetry.**
      Fixed: both platforms now build an explicit keep-list by scanning
      `nm -g` for `T _dd_*` / `T dd_*` and feeding it to `nmedit -s`
      (Darwin) / `objcopy --keep-global-symbols` (Linux). Previously
      Linux used `objcopy --localize-hidden` which only touches
      already-hidden symbols — a regression vector when a newly
      introduced helper ever leaked as default-visibility.

- [x] **A4. `DD_API` macro leaks.**
      Fixed: renamed to `MULTIFLOATS_API` (project-scoped) AND
      `#undef`ed at the end of the header so consumers never see the
      identifier.

- [x] **A5. Add `_Static_assert(sizeof(float64x2_t) == 16)`.**
      Fixed: added `static_assert` (C++) / `_Static_assert` (C)
      branches immediately after the `float64x2_t` typedef. Asserts
      `sizeof(float64x2_t) == 2 * sizeof(double)` so the check is portable
      to any hypothetical non-IEEE double platform.

## Minor / nits (batch when convenient)

- [x] `src/multifloats_math.cc` — union type-pun in `dd_erfc_scaled_full` replaced
      with `std::memcpy` round-trip (C++17-safe; compiler folds to a reg move).
- [x] `src/multifloats_math.cc` — atan/asin tiny-|x| fast-path threshold
      widened from `2.4e-17` to `1.85e-16 ≈ sqrt(3)·2⁻⁵³`; cubic term still
      below 0.5 DD ulp. Inline comment records the bound.
- [x] `src/multifloats_math.cc` — Cody-Waite comment no longer cites "~161
      bits" (was confusing against the "~106 bits preserved" sentence right
      after); now describes the constant as "three back-to-back doubles".
- [x] `fsrc/multifloats.fypp` — `float64x2_t` bind(c) type uses `real(c_double)`
      instead of `real(dp)`.
- [~] `fsrc/multifloats.fypp:638-664` — dead `c_dd_tgamma`/`c_dd_lgamma`
      bind(C) declarations. **Obsolete:** gone after the M3 refactor
      (C_DELEGATE_UNARY_MAP covers them).
- [~] `fsrc/multifloats.fypp` product/dim — naive vs Neumaier. **Obsolete:**
      subsumed by fixed **P4**.

Deferred (flag for user review before applying):

- [ ] `src/multifloats.hh` `nextafter` — below a power of 2 the down-side
      `ulp` halves, so `eps = ldexp(ulp, -53)` undercounts by 2×.
      Fix changes semantics (may skip a DD-representable neighbor), so it
      wants a deliberate decision + a targeted test.
- [ ] `src/multifloats_math.cc` atan cutover `|x|≥10.25` — puts
      `|t|=1/|x|` past the rational's 0.09375 validity edge; shift to
      `|x|≥10+2/3`. Needs a fuzz re-sweep to confirm no new worst-case.
- [~] `fsrc/multifloats.fypp` `ax%limbs = -ax%limbs` sites. **Unfounded:**
      checked the generated assembly for both whole-array and per-limb
      forms on gfortran-15 arm64 at both `-O2` and `-O3`. Both emit
      identical code: a single 128-bit `ldr q31` + `fneg v31.2d, v31.2d`
      + `str q31`. The fixed-size 2-element array is recognized and
      vectorized with no temp allocation. Additionally, the hotter bessel
      sites at 1895 / 1901 fire once at function entry, not in an inner
      loop. Leave as-is — splitting to per-limb would just hurt
      readability without changing codegen.
- [x] `fsrc/multifloats.fypp` `mf_dot_product` — hoisted the `ri > 0`
      test out of the hot loop: the `ri <= 0` path is one tight FMA
      pass; the periodic-renorm path is a nested block loop with the
      fast_two_sum between blocks. Measured ~15% throughput gain across
      N=8/64/1024/16384 (3.02→2.57 ns per FMA at N=1024, gfortran-15
      arm64 M1 Max), no regression at any size. Tests pass.
