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
- [x] **P7b — `pow_full` residual ulp_dd floor.** *(Partially fixed 2026-04-19.)*
  After P7's TD range-reduction landed, `pow` sat at 3.86e-31 (≈15.7 ulp_dd)
  max at fuzz seed 0x2a. Root cause: `log2_kernel` is DD-output, and the
  downstream `y·Lf` DD multiply adds rounding at `~ulp(p01+p10)` which for
  `|y·Lf|~17` leaks ~2e-31. _Resolved (approach a, partial):_ offline-generated
  `log2_values_extra[32]` in `dd_constants.hh` (third-limb residual of each
  tabulated `log2(center)` via MPFR @300-bit); rewrote `pow_full`'s `y·log2(x)`
  accumulator to (i) fold `y·Lf_extra[idx]` in, (ii) replace the DD product
  `y·Lf` with explicit `two_prod` sub-products `y.hi·Lf.hi`, `y.hi·Lf.lo`,
  `y.lo·Lf.hi`, (iii) accumulate via a `two_sum` cascade with a carry limb `c`
  preserved through the integer-split via `two_sum(l, c, l, e)` (third-limb
  `e` survives into the fractional reduction). Fuzz @200k: `pow` 3.86e-31 →
  1.90e-31 (≈15.7 → 7.7 ulp_dd, **2.03× improvement**). Bench: 5.96× →
  6.11× vs libquadmath (unchanged within noise). The residual floor is now
  the DD-accuracy of `log2_kernel` itself (~ulp_dd·|y|·|Lf|): getting below
  ~7 ulp_dd requires a full TD-internal `log2_kernel` (approach c, large
  refactor), deferred. Severity: **low** (7.7 ulp_dd on a ~32-digit type
  is acceptable).
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
- [x] **P10 — `erfc` deep-tail "4.6e6 ulp_dd" is a DD-format floor, not a
  kernel bug.** *(Investigated 2026-04-19.)* Fuzz at seed 0x2a reported
  `erfc` max_dd ≈ 1.12e-25 (≈4.55e6 ulp_dd). A fine sweep of
  `erfc(x)` for x ∈ [25, 27] against 200-bit MPFR shows: ≤ 1 ulp_dd for
  `|result| ≥ 2^-969`, then smooth degradation to ~80 ulp_dd at 2^-979,
  ~1000 ulp_dd at 2^-983, 4.5e14 ulp_dd at 2^-1020, and 1.0 (full loss)
  past 2^-1074. **Root cause:** the DD pair is two `double`s, and a
  `double`'s exponent cannot go below 2^-1074. Once `hi < 2^-969`, the
  companion `lo` at ~hi·2^-53 is forced into the double-subnormal range
  and loses mantissa bits proportionally. At hi = 2^-1050, lo carries
  only ~24 bits. This is a *format* limit, not an algorithm one: a
  libquadmath cross-check shows `erfcq` holds ≤ 1 ulp_q across the same
  range because `__float128`'s 15-bit exponent reaches ~2^-16494, so
  `erfc(27) ≈ 5e-319` is still a fully-normal float128 with all 112
  bits intact.
  **Attempted fix that made things worse:** folding the asymptotic's
  `exp(-s²-0.5625) · exp(diff_sq+p)` into a single `exp` call on the
  combined DD argument regresses the normal range (0.2 → 12.7 ulp_dd at
  x=8) because a DD argument of magnitude ~700 has absolute error
  `ulp_dd·700`, which maps 1:1 to exp's relative output error. The
  two-exp split was already correct; it lets each exp run at full DD
  precision on a small argument.
  **Resolution:** tightened the fuzz subnormal-floor filter from
  `|ref| < 2^-1050` to `|ref| < 2^-969` so the measurement reflects
  kernel quality, not the format cliff. Post-filter `erfc` max_dd is
  ≤ 1 ulp_dd across the supported range. A genuine fix would require
  escaping DD (e.g. a separate integer-exponent path for deep-tail
  outputs, or a triple-double internal representation); deferred as
  out of scope until a concrete caller needs it.
- [x] **P11 — `lgamma`/`tgamma` kernel seam at x=13.5.** *(Fixed 2026-04-19.)*
  `lgamma_stirling` comment declared "x ≥ ~20" but `lgamma_positive` called it
  directly at x ≥ 13.5, where the 13-term asymptotic series is under-converged
  (`c_13/x^25` ≈ `10^-14` at x=13.5). Fuzz measured `lgamma` ≈ 1.5e3 ulp_dd and
  `tgamma` ≈ 3.8e4 ulp_dd, with the worst case at x ≈ 13.5 exactly. Kernel seam,
  not a DD-format cliff (results at ~21 and ~1.7e9 are fully normal).
  **Fix:** shift recurrence in `lgamma_positive` for 13.5 ≤ x < 25 —
  accumulate `prod = x·(x+1)·…·(x+N-1)` in DD (N ≤ 12), call
  `lgamma_stirling(x+N)`, subtract `log_full(prod)`. Result:
  `lgamma` 1.5e3 → ~2.3 ulp_dd, `tgamma` 3.8e4 → ~160 ulp_dd.
  Bench essentially unchanged (12.57× → 12.55× tgamma; 4.79× → 4.82× lgamma).
  **Follow-up (P11b, 2026-04-19):** direct-Stirling `tgamma` path replaces
  `exp(lgamma)` for x ≥ 13.5 (and for `Γ(1-x)` in the reflection branch),
  mirroring libquadmath's `tgammaq`. New `tgamma_stirling_direct`: shift
  to x_adj ≥ 24, factor `x = mant · 2^e` (mant ∈ [√½, √2)), compute
  `ret = pow(mant, x) · exp2(e·x_frac) · exp(-x) · sqrt(2π/x)`, then
  `ret · (1 + expm1(Σ B_2k/(2k(2k-1))/x^(2k-1)))`; the `2^(e·x_int)`
  factor is pulled out as `exp2_adj` and applied by `ldexp` at the end
  so no intermediate overflows. Biggest exp-argument drops from
  `log Γ(99) ≈ 359` to `x·log(mant) ≲ 63`, cutting the amplification
  floor ~10×. Fuzz @ 200k: `tgamma` 3.9e-30 → 3.5e-31 (~157 → ~14
  ulp_dd). Point probe at x=99: ~22 ulp_dd (was ~300); x=90: ~9;
  x ≤ 50: ≤ 4. Bench: 11.3× → 13.3× (direct path is cheaper than
  exp(lgamma_rational) at large x since it dodges the 14-piece
  rational schedule). All 9 ctest tests pass.
- [x] **P9 — `atanh_full` ~40-ulp tail.** *(Fixed 2026-04-19.)*
  `0.5·log((1+x)/(1-x))` loses bits when the ratio ≈ 1 (|x| small but
  above the 0.01 Taylor threshold). Replaced with `0.5·log1p(2x/(1-x))`
  across the whole non-Taylor range — log1p routes to `log` internally
  once the argument grows, so the single formula handles |x| near 1 too.
  Fuzz @ 200k: `atanh` 9.4e-31 → 5.3e-32 (~38 → ~2.1 ulp_dd, 18×).
  Bench: 5.59× → ~6.15× vs libquadmath (10% **faster** — one fewer full
  `log_full` division round).

## Speed

- [x] **S1 — `neval` Estrin ladder stops at degree 14.** `src/multifloats.hh:747-865`.
  `expm1_taylor` (degree 24) and `log1p_taylor` (degree 17) fell back to
  Horner. Severity: **low**.
  _Resolved (2026-04-19):_ added `neval` cases 17 and 24 (x16 = x8·x8 + a
  second-tier Estrin block), and switched `expm1_full` / `log1p_full` from
  `horner` to `neval`. Precision max_dd bit-identical on every op.
  Bench: `expm1` 3.92× → 4.33× (+10%), `log1p` 5.96× → 6.45× (+8%),
  `atanh` 6.16× → 6.67× (+8%, inherits log1p via P9),
  `asinh` 8.80× → 8.96× (+2%, hits log1p for |x|<1). All 9 ctests pass.
- [x] **S2 — Bessel `pq0/pq1` branch trees.** `src/multifloats_math_bessel.inc:7-41`.
  Audit said "7-deep"; the trees are actually **3-deep** balanced (8 leaves via
  `xinv_d ≤ 0.25 / 0.125 / 0.0625 …`). Severity: **low** (cold path).
  _Closed — won't-fix, measured (2026-04-19):_ replaced the 3-way tree with a
  computed-index switch (`idx = (int)(xinv_d·16.0)`, clamp to [0,7], then
  switch). Fuzz max_dd bit-identical across j0/j1/jn/y0/y1/yn (boundary shift
  is harmless because fits overlap to DD precision). Bench (mean of 5
  `fortran_bench` runs, vs libquadmath) before/after:
  ```
             tree      switch    Δ
    j0       6.23×     6.25×     +0.02×
    j1       5.90×     6.01×     +0.11×
    jn(3,·)  4.77×     4.48×     −0.29×
    y0       6.58×     6.19×     −0.39×
    y1       6.61×     6.40×     −0.21×
    yn(3,·)  6.27×     6.20×     −0.07×
  ```
  Run-to-run noise is ±0.15× on this bench; all deltas sit inside noise with
  a slight adverse drift. The 3-compare tree the compiler already emits is
  not the bottleneck — `neval`/`deval` over a 9–11-term DD rational plus
  `sincos_full` + DD `sqrt` dominate totally. Change rolled back.
- [x] **S3 — `hypot` per-limb `ldexp`.** `src/multifloats.hh:1093-1104`.
  Audit said "eight libm calls"; actual count at N=2 is **six**
  (downscale 2 limbs × {big, small} = 4, upscale 2 limbs of root = 2).
  Severity: **negligible** (audit) — turned out to be the hot path.
  _Resolved (2026-04-19):_ replaced per-limb `ldexp` with one
  `down = ldexp(1, -e)` + one `up = ldexp(1, e)` and multiplies against
  each limb. Power-of-2 multiplier + non-subnormal result ⇒ bit-identical
  to `ldexp`, but one FP op instead of a libm call. Bench (mean of 5
  `cpp_bench` runs vs libquadmath): **7.98× → 12.76× (+60%)**. Fuzz
  `max_rel` bit-identical at 3.990e-32 (1.62 ulp_dd). All 9 ctests pass.
- [x] **S4 — Division is single-refinement.** `src/multifloats.hh:236-277`.
  Audit claimed ~48–53 bits (i.e. barely above `double`). Severity:
  **informational**.
  _Closed — audit claim not reproducible (2026-04-19):_ the current operator
  is textbook Dekker: `q1 = hi/rhs.hi`, compute residual `r = this − q1·rhs`
  as a **full DD** (`two_prod(q1, rhs.hi)` + `one_prod(q1, rhs.lo)`, cascaded
  via `two_sum`/`fast_two_sum`), then `q2 = r.hi/rhs.hi`, final
  `fast_two_sum(q1, q2)`. Because the residual is DD (not just a single
  limb), the refinement recovers ~53 additional bits and the result is full
  DD precision. `cpp_fuzz_mpfr` (200-bit MPFR reference, 10k samples)
  confirms:
  ```
    op    max_rel          ulp_dd
    add   1.619e-32        0.66
    sub   1.740e-32        0.71
    mul   4.417e-32        1.79
    div   6.521e-32        2.65
  ```
  div sits at 2.65 ulp_dd — same order as mul, *four* orders of magnitude
  tighter than the audit-claimed 48-bit ceiling (which would be ~4e16 ulp_dd
  relative). No code change needed; AUDIT_TODO.md now reflects the verified
  numbers so future readers aren't misled.

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
- [x] **M2 — Fortran bindings hand-mirrored from C ABI.** `fsrc/multifloats.fypp:114-145`
  vs `src/multifloats_c.h:219-230`. No build-time sync check. Severity: **medium**.
  _Resolved (2026-04-19):_ added `scripts/check_fortran_abi_sync.sh` and
  wired it as ctest `fortran_abi_sync`. For every `bind(c, name='*dd*')`
  in the fypp-generated `build/fsrc/generated/multifloats.f90`, verify the
  target symbol is declared with `MULTIFLOATS_API` in `multifloats_c.h`.
  Verified catches drift: injecting `bind(c, name='sindd_BOGUS')` into
  the generated file makes the test fail with the expected diagnostic.
  The check is intentionally asymmetric — the C header remains a superset
  because Fortran reimplements basic arithmetic (add/sub/mul/div/fma,
  comparisons, complex helpers) natively.
- [x] **M3 — `localize_symbols.sh` regex is fragile.** `src/CMakeLists.txt:25-27`,
  `src/localize_symbols.sh.in`. Pattern `dd(_[a-z]+)?$` silently mislabels any
  future non-`dd_*` export; Darwin vs Linux `nm` output differs. Severity:
  **medium**.
  _Resolved (2026-04-19):_ rewrote the keep-list to come from
  `multifloats_c.h` directly — `grep -oE "^MULTIFLOATS_API[^(]+\("` gives
  the exact set of public symbols, bypassing any heuristic. The header
  path is baked in at configure time via `@CMAKE_CURRENT_SOURCE_DIR@`.
  Platform handling kept (Darwin/nmedit wants a leading-underscore list;
  Linux/objcopy wants plain names). Verified: `libmultifloats.a` exports
  exactly 95 T symbols, matching the 95 `MULTIFLOATS_API` declarations.
- [x] **M4 — `dd_constants.hh` has no CMake regeneration rule.** `src/dd_constants.hh:1-4`.
  Stale constants can ship if `scripts/gen_constants.py` changes and the
  developer forgets to re-run it. Severity: **low**.
  _Resolved (2026-04-19):_ added a `dd_constants_up_to_date` ctest that
  runs `scripts/gen_constants.py --check` via `uv run --with mpmath` when
  `uv` is available, falling back to the system `python3` otherwise. A
  drift between the script and the committed header now surfaces as a
  ctest failure, no manual remembered step needed.
- [x] **M5 — Inconsistent preamble comments across `.inc` files.**
  `multifloats_math_trig.inc` is documented; `bessel.inc`, `special.inc` are
  not. Severity: **low**.
  _Resolved (2026-04-19):_ `bessel.inc` already had a preamble; added a
  multi-paragraph block to `special.inc` covering the erfc Clenshaw–Burrus
  split, the tgamma direct-Stirling switchover at x ≥ 13.5, and the
  lgamma shift recurrence on [13.5, 25). Cross-references P4 / P11 /
  P11b for readers chasing the numerical rationale.
- [x] **M6 — Generated Fortran not in VCS.** Reviewers can't see the actual
  interfaces without a local build. Consider committing the generated file or
  diffing it in CI. Severity: **low**.
  _Closed — won't-fix with rationale (2026-04-19):_ committing generated
  output is a known anti-pattern — every non-trivial `.fypp` edit would
  produce a "regenerated file" diff alongside the real change, and a
  developer who forgets to regenerate pushes a silently-stale tree. The
  readability concern is real, though: the M2 fix above (ctest
  `fortran_abi_sync`) already catches the class of drift reviewers would
  most want to audit (C↔Fortran signature mismatch), and the fypp source
  itself is short enough to read directly. If a reviewer needs the
  expanded module for a one-off audit, `cmake --build build --target
  multifloatsf && cat build/fsrc/generated/multifloats.f90` produces it
  in seconds. No VCS pollution needed.

---

## Legend

- **Severity: high** — latent UB, silent correctness regression, or ABI break.
- **Severity: medium** — IEEE/C99 deviation, precision loss, or maintenance
  landmine that will eventually bite.
- **Severity: low / informational** — quality-of-life and doc items.

Fixes land one-per-commit so the benchmark/fuzz deltas stay bisectable.
