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
      in `dd_gaxpy_mv_panel`, and the periodic renorm in `dd_matmul_vm`
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

- [ ] **S1. `dd_matmul_mm` re-streams A per B column.**
      `src/multifloats_math.cc:1388-1394`. Add a second register-block
      dimension NR so each MR×NR µkernel loads A once per p and reuses
      it across NR b-column entries. Preserves the renorm-interval
      contract; only changes the outer dispatch.

- [ ] **S2. Missing `__attribute__((always_inline))` on hot kernels.**
      `src/multifloats_math.cc:1200-1212` (`dd_mac_inl`),
      `:1214-1218` (`dd_finalize_inl`), and `dd_gaxpy_mv_panel`. Some
      compilers refuse inlining at `-O2`. Benchmark before/after to
      confirm codegen win (it's a hint, not a correctness change).

- [ ] **S3. Bessel asymptotic doubles the sin/cos work.**
      `src/multifloats_math.cc:1054-1079, 1110-1144, 1146-1180`. Add a
      `dd_sincos_eval(angle)` helper that runs range reduction once and
      returns `(s, c)`; cut Bessel large-x cost ~40%. Also applies to
      `dd_tan_full` at `:308-310` (identical branches 0 and 2).

- [ ] **S4. `dd_reduce_pi_half` re-runs full-DD subtraction three times.**
      `src/multifloats_math.cc:204-220`. Restructure as three inline
      `two_sum` steps; save ~12 FLOPs per sin/cos/tan call.

## Maintainability

- [ ] **M1. Naming drift: `exp2_min_d` (C++) vs `exp2_min` (Fortran).**
      `scripts/gen_constants.py:1245, 1300`. Pick one; the `write_cpp` /
      `write_f90` split is the root cause — consolidate.

- [ ] **M2. Fortran inlines Taylor coefficients already emitted.**
      `fsrc/multifloats.fypp:1611-1620` (sinh), `:1858-1873` (asinh),
      `:2206-2220` (atanh). Replace with `dd_horner(x2, *_taylor_hi,
      *_taylor_lo, N)` using the generated arrays.

- [ ] **M3. Arithmetic-dispatch mega-template in fypp.**
      `fsrc/multifloats.fypp:4565-4744` — single ~180-line body expanded
      ~50×. Extract per-op bodies (`add` / `sub`, `mul`, `div`, complex)
      into named fypp `#:def` macros.

- [ ] **M4. 8-way erfc dispatch duplicated.**
      `fsrc/multifloats.fypp:2060-2085` vs `:2138-2163`. Extract
      `dd_erfc_asym_pq(idx, z, p)` helper.

- [ ] **M5. Section banners.**
      `fsrc/multifloats.fypp` is 5000 lines with three top-level
      banners. Promote the ~20 inline `! ====` sub-banners to a
      consistent style.

## API / ABI

- [ ] **A1. No ABI version macro.**
      Add `#define MULTIFLOATS_ABI_VERSION 1` in
      `src/multifloats_c.h`, or set `SOVERSION` via CMake.

- [ ] **A2. `VISIBILITY_INLINES_HIDDEN` not set.**
      `src/CMakeLists.txt:11`. Add `VISIBILITY_INLINES_HIDDEN ON` to
      reduce reliance on the post-build symbol scrub.

- [ ] **A3. `localize_symbols.sh.in` Darwin/Linux asymmetry.**
      Linux path uses `objcopy --localize-hidden` (only touches
      already-hidden symbols); Darwin uses `nmedit -s` with explicit
      preserve list. Mirror Darwin behavior on Linux (explicit keep list
      built from `nm -g | grep '^dd_'`).

- [ ] **A4. `DD_API` macro leaks.**
      `src/multifloats_c.h:17-21`. `#undef DD_API` at end of header, or
      rename to project-scoped `MULTIFLOATS_API`.

- [ ] **A5. Add `_Static_assert(sizeof(dd_t) == 16)`.**
      `src/multifloats_c.h`. Cheap defensive catch for surprise padding.

## Minor / nits (batch when convenient)

- [ ] `src/multifloats_math.cc:689-692` — union type-pun → `std::bit_cast<uint64_t>`.
- [ ] `src/multifloats.hh:530-544` — `nextafter` undercounts near powers of 2.
- [ ] `src/multifloats_math.cc:385, 430` — conservative `2.4e-17` threshold vs true `sqrt(3)·2^-53 ≈ 1.85e-16` in atan/asin (no accuracy loss, tiny perf).
- [ ] `src/multifloats_math.cc:392-400` — atan cutover `|x|≥10.25` puts `|t|=1/|x|` past rational's 0.09375 validity edge; shift to `|x|≥10+2/3`.
- [ ] `src/multifloats_math.cc:202` — comment overstates Cody-Waite reduction as "~161 bits"; actually ~106.
- [ ] `fsrc/multifloats.fypp:638-664` — dead `c_dd_tgamma` / `c_dd_lgamma` bind(C) declarations (Fortran reimplements; doesn't call).
- [ ] `fsrc/multifloats.fypp:641-643` — `real(dp)` in bind(c) type; idiomatically `real(c_double)`.
- [ ] `fsrc/multifloats.fypp:3521, 3527, 3605, 3643` — `ax%limbs = -ax%limbs` may allocate temps; split to per-limb.
- [ ] `fsrc/multifloats.fypp:4140-4173` — `if (ri > 0 .and. mod(k, ri) == 0)` inside hot loop; restructure as nested loops with `ri` as outer stride.
- [ ] `fsrc/multifloats.fypp:3966-4011` vs `fsrc/multifloats.fypp` product/dim — `sum` is Neumaier-compensated, other reductions naive (duplicate of P4).
