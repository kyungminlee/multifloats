# AUDIT_TODO

Roadmap of open items surfaced by the April-2026 source-tree audit of
`src/` and `fsrc/`. Each entry is referenced by `multifloats_c.h` / `README.md`
as the canonical tracking file.

Entries are grouped by category and ordered by severity. Checkboxes flip as
fixes land. File:line references are snapshots taken at the time of the audit.

---

## Correctness & edge cases

- [x] **C1 — `sinpi/cospi/tanpi` overflow is UB.** `src/multifloats_math_trig.inc:14,38,62`.
  `nearbyint(2.0 * ax._limbs[0])` overflows to ∞ for |x| ≳ 1e308, then the
  `long long` cast is undefined. Severity: **high**.
  _Resolved:_ added `pi_trig_arg_too_large(ax)` guard that returns NaN when
  `|x| ≥ 2^52`, the point where a double can no longer resolve the
  fractional part of `2x` (and well below the `∞` threshold that would
  cause the UB). All three entry points share the helper.
- [x] **C2 — `ctandd` divides by zero.** `src/multifloats_math_abi_complex.inc:185`.
  `den = ca*ca + sb*sb` has no zero guard; returns NaN/Inf via raw division
  instead of the C99 Annex G values. Severity: **high**.
  _Resolved:_ added two guards before the divide. `den.hi == 0` ⇒ return
  `(NaN, NaN)` (C99 Annex G.6.3.1 value at a real-axis pole). Non-finite
  `den` (sinh overflow for `|Im z| ≳ 710`) ⇒ return `(0, sign(b)·i)`, the
  mathematical limit of `tan(a + i·b)` as `|b| → ∞`.
- [x] **C3 — `cpowdd(0, w)` drops signed-zero info.** `src/multifloats_math_abi_complex.inc:128`.
  Zero check inspects only hi limbs; (−0, ε) is misclassified, sign of result
  violates C99 G.6.4.1. Severity: **medium**.
  _Resolved:_ two tightenings in `cpowdd`. (i) the "true zero" predicate
  now requires both limbs of both components to be zero, so a subnormal
  DD `(hi=0, lo=ε)` is no longer collapsed to zero. (ii) if either
  component has `hi == 0 && lo != 0` after that check, promote `lo` to
  `hi` (mirroring the log1p C5 fix) before calling `clogdd`, so the
  downstream `log2_full` zero short-circuit doesn't poison the result
  with NaNs. Verified via a targeted probe: `cpow((hi=0, lo=1e-300), 1)`
  now returns `(1e-300, 0)` instead of `NaN + NaN·i` (previous fix-step)
  or `0 + 0·i` (original bug).
- [x] **C4 — `casindd` branch-cut fixup ignores lo-limb sign.** `src/multifloats_math_abi_complex.inc:208`.
  `signbit(b._limbs[0])`-only; for `b=(+0,−ε)` imaginary sign flips wrong way
  on the real-axis cut. Severity: **medium**.
  _Resolved:_ introduced a DD-aware `dd_signbit(hi, lo)` helper that
  routes to `signbit(lo)` only when `hi == 0.0 && lo != 0.0` (so `±0`
  with a zero lo keeps `signbit(hi)`, preserving signed zero). Applied
  to both operands of the branch-cut comparison. Verified via a probe
  over `{+0, −0, (+0,−ε), (−0,+ε)}` for `Im(z)` with `Re(z)=2`:
  `Im(asin)` sign now matches the DD-level sign of the input in all
  four cases (previously `(+0,−ε)` and `(−0,+ε)` flipped the wrong
  way). cpp_fuzz: 0 failures, no tolerance regression.
- [x] **C5 — `log1p` subnormal fall-through.** `src/multifloats_math_exp_log.inc:171`.
  When `(1+x).hi == 0` but `lo ≠ 0`, control drops to `log_full(0) → −∞`
  instead of a large-negative finite value. Severity: **medium**.
  _Resolved:_ when the subnormal branch is taken, promote the lo limb to
  a normalized scalar via `float64x2(arg._limbs[1], 0.0)` before calling
  `log_full`, so the leading-limb zero check inside `log_full`/`log2_full`
  is not tripped.
- [x] **C6 — FMA assumed but not asserted.** `src/multifloats.hh:41`, trig.inc:97, special.inc:101.
  If `FP_FAST_FMA` is false, `two_prod` silently loses the error term and every
  DD operation degrades to double. Add `static_assert(FP_FAST_FMA)` or a runtime
  fallback. Severity: **medium** (silent correctness regression risk).
  _Resolved — narrower scope than originally framed:_ C99/C++ `std::fma`
  is always correctly rounded (even when emulated in software), so a
  hardware-FMA assert would reject valid conforming builds. The real
  silent-correctness risk is `-ffast-math`, which lets the compiler
  rewrite fma into a non-IEEE sequence. Added `#error` at the top of
  `src/multifloats.hh` when `__FAST_MATH__` is defined.

## Precision

- [x] **P1 — `sqrt` Karp–Markstein is single-limb.** `src/multifloats.hh:980-984`.
  `residual._limbs[0] * (0.5 / s)` only; full DD division gives ~1–2 ulp
  improvement near perfect squares, and feeds `asin/acos` near ±1. Severity:
  **medium**.
  _Closed — keep baseline, documented in code:_ the audit's "1–2 ulp loss
  near perfect squares" framing turned out to overstate the problem.
  A targeted probe (200k exact squares k² and 1.2M + 400k near-square
  samples, reference = libquadmath `sqrtq`, 1 ulp ≈ 2.47e-32) shows:
  ```
                            baseline   (b) DD×scalar   (a) full DD divide
    exact k²                0.00 ulp   0.00 ulp        0.00 ulp
    k² ± small lo limb      0.76 ulp   0.58 ulp        0.39 ulp
    k² ± ldexp(v, −54)      0.19 ulp   0.19 ulp        0.02 ulp
  ```
  Speed costs on pop-os (`cpp_bench` sqrt): baseline 52.6×, (b) ~36×
  (~30% slower), (a) ~24× (~55% slower); (a) also slows hypot ~13% and
  acosh ~24%. Baseline is already sub-1-ulp in every case, so the gain
  isn't worth the regression. A comment in `multifloats.hh` records the
  trade-off numbers so the optimization isn't re-attempted blindly.
  Cheap variant `(residual.hi + residual.lo) * (0.5/s)` was also tried
  and gives zero precision improvement (lo is below double rounding).
- [x] **P2 — `asin/acos` near ±1.** `src/multifloats_math_inv_trig.inc:87,121,128`.
  Region-3 half-angle path inherits the P1 error; `asin(±1)` only matches π/2
  to double. Severity: **medium** (gated by P1).
  _Closed — audit framing not reproducible:_ at exact `x = ±1` both
  `asin` and `acos` return within 0.04 ulp of the true value (π/2, 0,
  π). Global fuzz shows `asin` max 0.76 ulp, `acos` max 0.56 ulp —
  already sub-1-ulp. An initial denser probe over DD-extended inputs
  in [0.9, 1.0) appeared to show a 15.8-ulp spike for `acos` at
  `x = (0.99993902994472146, −2.13e-19)`, which was traced further
  and found to be a **float128 reference-precision artifact**, not a
  real implementation bug — see investigation note below.
  _Reference-noise note:_ when checking DD output against
  `1.0Q − ((__float128)xhi + (__float128)xlo)`, the inner
  `xhi + xlo` sum loses ~1–2 bits whenever `|xlo|` dips below
  ~2^-60 (below float128's 113-bit trailing-mantissa floor for
  magnitude-1 sums). Replacing the reference with the pairwise form
  `(1.0Q − (__float128)xhi) − (__float128)xlo` (each sub exact for
  double operands in float128) shows `got − ref = 0` bit-exact — the
  DD kernel was already optimal. The `acos` region-3 path and the
  underlying DD subtraction are correct for these inputs. Global
  fuzz uses `|xlo| ≈ 2^-53·|xhi|`, well above the float128 cliff, so
  its measurements were never polluted by this artifact.
- [x] **P3 — `pi_half_cw3` is a single double.** `src/dd_constants.hh:60-63`.
  Third Cody–Waite term limits reduction for |x| > 1e15. Severity: **low**.
  _Closed — reframed and capped:_ audit blamed cw3's single-double
  form, but empirically `reduce_pi_half` stays sub-1-ulp through
  `|x| ≤ 1.414·2^54 ≈ 2.5e16` (audit's 1e15 was conservative by ~25×).
  The true cliff is at `|x| ≥ 2^55`, where `nearbyint(x·2/π)` loses
  integer precision (ulp(x) ≥ 1 by then), so `n` is wrong by up to
  ±0.5 and no finite number of CW terms can recover the answer —
  adding a `cw4` wouldn't move the cliff. Worse, beyond the cliff
  `sin_full`/`cos_full`/`tan_full`/`sincos_full` returned
  progressively larger garbage (e.g. `sin(2^60)` ≈ -8.8e20) rather
  than bounded values, because `sin_kernel`'s Taylor series was
  evaluated on an unreduced argument. Added a `trig_arg_too_large`
  guard mirroring the existing `sinpi`/`cospi`/`tanpi` cap, returning
  NaN when `|x| ≥ 2^55`. Fuzz and bench on the supported range
  unchanged. Proper support for huge arguments would require
  Payne–Hanek reduction, which is out of scope for this library (see
  `std::sin` for double-precision Payne–Hanek).
- [x] **P4 — `erfc_full` asymptotic x-split asymmetry.** `src/multifloats_math_special.inc:100-104`.
  35-bit hi truncation makes `(s-ax)(s+ax)` asymmetric in precision. Severity:
  **low**.
  _Closed — audit cause incorrect; real cause is upstream:_ at the
  worst-case fuzz input (`x ≈ 7.29085`, erfc ≈ 328 ulp off), the
  `diff_sq = (s - ax)·(s + ax)` computed via the 35-bit-truncated `s`
  is **bit-exact** vs MPFR — the split trick is sound and its
  symmetry is not a real issue. The 300-ulp erfc tail is entirely
  inherited from `exp_full`: `exp(-53.72)` already has 279 ulp of
  relative error (verified against MPFR at 200-bit precision), and
  the asymptotic formula multiplies two `exp` calls, amplifying it.
  Global fuzz confirms the root cause sits at the exp/log level:
  `exp` max_rel 7.5e-30 (~304 ulp), `log` 6.97e-30 (~283 ulp),
  `log10` 6.98e-30 (~283 ulp). Tracked separately as **P6** below.
- [x] **P6 — `exp_full` / `log_full` precision tail (~300 ulp).**
  Discovered while investigating P4. The DD `exp`/`log` kernels emit
  worst-case relative errors around 300 ulp on clean single-double
  inputs (e.g. `exp(-53.72)`, `exp(1.0)`). Verified against MPFR,
  so it is not a float128 reference artifact. The erfc asymptotic
  branch, expm1, log10, pow, and every complex op that funnels
  through exp/log inherit this ceiling. Severity: **high**.
  _Closed — three-part fix:_
  (1) `exp2_min` lowered from −1022 to −1080 (`5eb1553`) so the
  subnormal-output branch goes through the kernel instead of
  short-circuiting to zero.
  (2) `exp2_kernel` Taylor degree bumped 13 → 15 and `log2_wide` degree
  8 → 10 (`0673187`) so the poly truncation drops well below the
  cube-amplification ceiling.
  (3) `exp_full` rewritten with a three-piece Cody–Waite split of ln2
  (`a997932`), removing the `y = x·log2e` cancellation that was the
  dominant source of the remaining ~220-ulp tail. The DD multiply
  injected ~|x|·ε_dd absolute error into `y`; subtracting
  `round(y.hi)` promoted that to ~1.8e-27 *relative*, which the cube
  then amplified. CW keeps r = x − n·ln2 at full DD precision.
  Result (cpp_fuzz_mpfr @ 200k iters, 200-bit mpreal reference):
    exp    7.5e-30 → 2.5e-31   (~304 ulp → ~10 ulp)
    expm1  5.2e-30 → 2.4e-31
    log    7.0e-30 → 3.0e-32   (~283 ulp → ~1 ulp)
    log10  7.0e-30 → 3.2e-32
    log1p  5.2e-30 → 6.3e-32
    sinh   5.2e-30 → 6.5e-31
    cosh   5.2e-30 → 2.3e-31
    tanh   5.2e-30 → 6.1e-31
  Speed cost: `exp` benchmark 6.10× → 5.45× vs libquadmath (~11%);
  log path unchanged. Residual tails in `pow`, `asinh`, `atanh` are
  now cancellation issues upstream of exp/log, not the exp/log
  kernels themselves — tracked as P7/P8/P9 below.
- [x] **P7 — `pow_full` ~100-ulp tail.** `src/multifloats_math_exp_log.inc:193`.
  Composition `exp(y · log(x))`: once P6 brought `log` to ~1 ulp_dd,
  the bottleneck became the DD multiply `y · log(x)` for large `|y|`.
  At the fuzz worst case `pow` sat at 2.4e-30 (~100 ulp_dd). Severity:
  **medium**.
  _Resolved:_ rewrote `pow_full` around `x^y = 2^(y·log2(x))` with a
  triple-double splitting trick. Factor `log2(x) = e_x + log2(m)` with
  `e_x = ilogb(x)` (integer) and `m ∈ [1, 2)`; then `y·e_x` is exact
  in triple-double form (`two_prod(y.hi, e_x)` + `two_prod(y.lo, e_x)`)
  and `y·log2(m)` is a DD product whose absolute error is bounded by
  `ulp_dd·|y|·|log2(m)| ≤ ulp_dd·|y|` regardless of |e_x|. Summing the
  pieces with cascaded `two_sum`s preserves every residual, so the
  fractional part fed into `exp2_kernel` has absolute error ~ulp_dd·|y|
  instead of ulp_dd·|y·log2(x)|. Fuzz @ 200k: `pow` 2.4e-30 → 3.9e-31
  (≈100 ulp_dd → ≈16 ulp_dd). The remaining tail is the exp2-kernel
  floor (cube-of-eighth + Taylor poly) at ~10 ulp_dd, shared with
  `exp_full`. Bonus: `pow` also got **faster**, 5.21× → 6.45× vs
  libquadmath (we now share a single range-reduction with exp2 instead
  of running log+exp back-to-back).
- [x] **P8 — `asinh_full` ~50-ulp tail.** *(Fixed 2026-04-19.)*
  The existing code already folded sign (`ax = |x|`), so the root cause
  wasn't `x + sqrt(x²+1)` cancelling — it was `log(ax + root)` where
  `ax + root ≈ 1` for small ax (log-of-near-1 cancellation). Fix: for
  `ax < 1`, compute `root - 1 = ax² / (root + 1)` accurately and return
  `log1p(ax + (root - 1))`. For `ax ≥ 1` the plain `log(ax + root)` path
  is preserved to avoid the log1p cost. Fuzz @ 200k: `asinh`
  1.2e-30 → 5.5e-32 (~48 → ~2.2 ulp_dd, 22×). Speed: 9.43× → ~8.7× vs
  libquadmath (~7% slower on the mixed-range bench).
- [x] **P9 — `atanh_full` ~40-ulp tail.** *(Fixed 2026-04-19.)*
  `0.5·log((1+x)/(1-x))` loses bits when the ratio ≈ 1 (|x| small but
  above the 0.01 Taylor threshold). Replaced with `0.5·log1p(2x/(1-x))`
  across the whole non-Taylor range — log1p routes to `log` internally
  once the argument grows, so the single formula handles |x| near 1 too.
  Fuzz @ 200k: `atanh` 9.4e-31 → 5.3e-32 (~38 → ~2.1 ulp_dd, 18×).
  Bench: 5.59× → ~6.15× vs libquadmath (10% **faster** — one fewer full
  `log_full` division round).

## Speed

- [ ] **S1 — `neval` Estrin ladder stops at degree 14.** `src/multifloats.hh:747-865`.
  `expm1_taylor` (degree 24) and `log1p_taylor` (degree 17) fall back to
  Horner, ~40–50% extra FLOPs on cold-tail paths. Severity: **low**.
- [ ] **S2 — Bessel `pq0/pq1` 7-deep branch trees.** `src/multifloats_math_bessel.inc:7-41`.
  Replace with small lookup keyed by `xinv_d/8`. Severity: **low** (cold path).
- [ ] **S3 — `hypot` per-limb `ldexp`.** `src/multifloats.hh:1027-1029`.
  Eight libm calls where one scale would do. Severity: **negligible**.
- [ ] **S4 — Division is single-refinement.** `src/multifloats.hh:227-263`.
  ~48–53 bits of precision, not full 106. Documented and intentional; tracked
  here so downstream callers aren't surprised. Severity: **informational**.

## Maintainability & hygiene

- [x] **M1 — `SOVERSION` hardcoded, not coupled to `MULTIFLOATS_ABI_VERSION`.**
  `src/CMakeLists.txt:14-15` vs `src/multifloats_c.h:21`. A manual bump of the
  header constant without a matching CMake edit produces a wrongly-tagged
  shared library. Severity: **medium**.
  _Resolved:_ `src/CMakeLists.txt` now uses `file(STRINGS …)` + a regex to
  pull `MULTIFLOATS_ABI_VERSION` out of `multifloats_c.h` at configure
  time and drives both `VERSION` and `SOVERSION` from it (fatal-errors if
  the header define can't be found). Verified to parse as `2` against
  the current header.
- [ ] **M2 — Fortran bindings hand-mirrored from C ABI.** `fsrc/multifloats.fypp:114-145`
  vs `src/multifloats_c.h:219-230`. No build-time sync check. Severity: **medium**.
- [ ] **M3 — `localize_symbols.sh` regex is fragile.** `src/CMakeLists.txt:25-27`,
  `src/localize_symbols.sh.in`. Pattern `dd(_[a-z]+)?$` silently mislabels any
  future non-`dd_*` export; Darwin vs Linux `nm` output differs. Severity:
  **medium**.
- [ ] **M4 — `dd_constants.hh` has no CMake regeneration rule.** `src/dd_constants.hh:1-4`.
  Stale constants can ship if `scripts/gen_constants.py` changes and the
  developer forgets to re-run it. Severity: **low**.
- [ ] **M5 — Inconsistent preamble comments across `.inc` files.**
  `multifloats_math_trig.inc` is documented; `bessel.inc`, `special.inc` are
  not. Severity: **low**.
- [ ] **M6 — Generated Fortran not in VCS.** Reviewers can't see the actual
  interfaces without a local build. Consider committing the generated file or
  diffing it in CI. Severity: **low**.

---

## Legend

- **Severity: high** — latent UB, silent correctness regression, or ABI break.
- **Severity: medium** — IEEE/C99 deviation, precision loss, or maintenance
  landmine that will eventually bite.
- **Severity: low / informational** — quality-of-life and doc items.

Fixes land one-per-commit so the benchmark/fuzz deltas stay bisectable.
