# Library internals

Cross-cutting reference for contributors working on the DD kernels,
build, and test surface. Distilled from the April-2026 source-tree
audit and the precision-improvement tiers that followed; the
item-by-item roadmap with resolutions, measurements, and commit
hashes lives in git history as the prior `doc/developer/AUDIT_TODO.md`
and `doc/developer/PROGRESS-PRECISION.md`.

Companion: `doc/developer/TRIPLE_DOUBLE.md` covers the narrow TD path
(`float64x3`, `exp_full_td`, `sincos_full_td`, `cexpm1dd` regime
split). This file points into it wherever the TD infrastructure is
the answer to a question raised here.

Source anchors `C1`–`C6`, `P1`–`P11b`, `S1`–`S4`, `M1`–`M6` in code
comments correspond to the historical roadmap — use
`git log --follow doc/developer/INTERNALS.md` to retrieve the
investigation behind any one of them.

---

## 1. Invariants the build relies on

The following are load-bearing; break one and the library silently
produces wrong answers or ships with the wrong ABI tag.

- **`-ffast-math` is forbidden.** `include/multifloats.h` emits an `#error`
  when `__FAST_MATH__` is defined. DD arithmetic relies on strict
  IEEE rounding of `fma` and on reassociation being blocked; fast-math
  silently violates both.
- **Hardware FMA is *not* required.** C99 `std::fma` is correctly
  rounded even in software. An early `static_assert(FP_FAST_FMA)` was
  rejected because it would refuse valid conforming builds; the real
  risk is `-ffast-math`, which is what we guard.
- **`SOVERSION` is derived from `MULTIFLOATS_ABI_VERSION`.**
  `src/CMakeLists.txt` greps the constant out of `include/multifloats.h`
  at configure time and drives both `VERSION` and `SOVERSION`. Bump
  the header; the CMake side follows. Manual `set(SOVERSION …)` edits
  are a bug.
- **Public symbols are the ones marked `MULTIFLOATS_API`.** The
  `localize_symbols.sh` helper extracts the export list by grepping
  `^MULTIFLOATS_API[^(]+\(` out of `multifloats.h` directly — it is
  *not* a regex over symbol names. New exports get the attribute;
  internal helpers do not.
- **Fortran bindings must mirror the C ABI.** ctest `fortran_abi_sync`
  (`scripts/check_fortran_abi_sync.sh`) walks `bind(c, name='*dd*')`
  declarations in the generated `multifloats.f90` and fails if the
  target symbol is missing from `multifloats.h`. The check is
  intentionally asymmetric — Fortran reimplements arithmetic natively,
  so the C header stays a superset.
- **`dd_constants.hh` regenerates cleanly.** ctest
  `dd_constants_up_to_date` runs `scripts/gen_constants.py --check`
  (via `uv` when available, falling back to system `python3`). Drift
  between script and header is a test failure.

## 2. Supported argument ranges

Where the kernels degrade and what they return outside the supported
range. Argument-range guards are narrower than the IEEE limits — they
exist to reject inputs the kernel mathematics cannot resolve, not just
to guard against overflow.

| Kernel                           | Guard                | Reason                                                                                  |
| -------------------------------- | -------------------- | --------------------------------------------------------------------------------------- |
| `sinpi` / `cospi` / `tanpi`      | `|x| < 2^52`, else NaN | `nearbyint(2x)` is exact below; above, the `(long long)` cast is UB.                  |
| `sin` / `cos` / `tan` / `sincos` | `|x| < 2^55`, else NaN | `nearbyint(x·2/π)` loses integer precision beyond 2^55; no finite number of Cody–Waite terms recovers the answer. |
| `erfc` tail                      | result `≥ 2^-969`    | DD's lo limb becomes subnormal below this threshold (format limit, not kernel).         |

Huge-argument trig via Payne–Hanek reduction is out of scope — `std::sin`
in libm is the workaround. The DD kernel intentionally refuses rather
than returning silently-wrong large values.

The fuzz harnesses filter by these same thresholds
(`cpp_fuzz_mpfr` drops `|ref| < 2^-969` for `erfc`), so reported max_ulp
numbers reflect kernel quality, not format cliffs.

## 3. Verified precision envelope

From `cpp_fuzz_mpfr` at 200-bit MPFR reference, post-audit:

| Op       | max_rel    | ≈ ulp_dd                                                      |
| -------- | ---------- | ------------------------------------------------------------- |
| `add`    | 1.6e-32    | 0.7                                                           |
| `sub`    | 1.7e-32    | 0.7                                                           |
| `mul`    | 4.4e-32    | 1.8                                                           |
| `div`    | 6.5e-32    | 2.7                                                           |
| `sqrt`   | —          | 0 ulp on exact `k²`; ≤ 0.76 ulp with a non-zero lo limb       |
| `exp`    | 2.5e-31    | ~10                                                           |
| `log`    | 3.0e-32    | ~1                                                            |
| `log1p`  | 6.3e-32    | ~3                                                            |
| `pow`    | 1.9e-31    | ~8                                                            |
| `asinh`  | 5.5e-32    | ~2                                                            |
| `atanh`  | 5.3e-32    | ~2                                                            |
| `tgamma` | 3.5e-31    | ~14                                                           |
| `lgamma` | —          | ~2 (after the `x ≥ 13.5` shift-recurrence fix)                |
| `cmul`   | 2.2e-31    | ~9 (Re and Im, after the 2⁻⁴ cancellation-gated compensation) |

1 ulp_dd ≈ 2⁻¹⁰⁵ ≈ 2.47e-32. These numbers are the bar: a change that
regresses any of them is rejected unless traded for a proportional win
elsewhere (documented in `doc/BENCHMARK.md`).

## 4. Designs measured and rejected

Each of these has been tried under audit-quality fuzz and bench; do not
re-attempt without new evidence that overrides the recorded trade-off.

- **Full DD divide in `sqrt`'s Karp–Markstein step** (`include/multifloats.h`
  near line 1200). Cuts worst-case residual 0.76 → 0.39 ulp near perfect
  squares, costs ~55% on the sqrt bench and 10–25% on hypot/acosh.
  Baseline is already sub-1-ulp. `residual * (0.5/s)` as DD×scalar was
  also tried: 0.76 → 0.58 ulp, 30% slower. A code comment at the sqrt
  site records the numbers.
- **Karatsuba complex multiply.** Saves one mul at the cost of a
  catastrophic cancellation in `Im = (a+b)(c+d) − ac − bd` when the
  true Im is near zero. Witness `a=(1,ε), b=(−1,ε)` — pinned in
  `test/test.cc::test_complex_mul_cancellation`. The 4-mul form stays.
- **Bessel `pq0/pq1` branch tree → computed switch.** A 3-deep balanced
  tree on `xinv_d` already compiles to near-optimal code; replacing
  with `(int)(xinv_d·16.0)` + switch runs inside bench noise
  (±0.15× on a 5-run mean), sometimes slower. The hot path is
  `neval`/`deval` plus `sincos_full` plus DD `sqrt`, not the dispatch.
- **`erfc` two-exp fold into one DD-arg `exp`.** Collapsing
  `exp(−s²−0.5625)·exp(diff_sq+p)` into a single `exp` on the DD sum
  regresses the normal range (0.2 → 12.7 ulp_dd at `x = 8`): a DD
  argument of magnitude ~700 carries `ulp_dd·700` absolute error,
  which maps 1:1 to the relative error of the output. Keep the split.
- **Bailey-style DD accumulation in `mac_inl`.** Zero precision
  improvement (the lo-limb tail is O(eps³·|s_hi|), below DD), ~9%
  slower on `arr_matmul`. The reported "288 ulp" on dot-product
  witnesses is a fuzz-metric artifact: `rel_err = abs / |result|`
  diverges when cancellation makes `|result|` tiny relative to the
  inputs.
- **Committing the fypp-generated `multifloats.f90` to VCS.** Every
  non-trivial `.fypp` edit produces a spurious "regenerated" diff, and
  forgetting to regenerate lets stale trees ship. The
  `fortran_abi_sync` ctest already catches the class of drift a
  reviewer cares about.

## 5. Pitfalls the audit uncovered

Traps that cost real time during the audit; they do not show up until
you hit exactly the right input.

- **float128 reference precision has a cliff.** When computing `(1 − x)`
  with `x = xhi + xlo` as a DD, `(__float128)xhi + (__float128)xlo`
  loses ~1–2 bits for `|xlo| < 2^-60`. Use the pairwise form
  `(1.0Q − (__float128)xhi) − (__float128)xlo` instead — each sub is
  exact for double operands inside float128. Ignore this and you'll
  chase ghost bugs in `acos` near `x = 1`.
- **Subnormal-lo signed zero.** `signbit(hi)` answers the wrong
  question when `hi == 0 && lo != 0`; use `dd_signbit(hi, lo)` for
  branch-cut code. Applied in `casindd`, used wherever the DD sign
  matters below the hi-limb's resolution.
- **`log_full(0) = −∞`.** If a subnormal DD `(hi=0, lo=ε)` drops into
  `log_full` via short-circuit, the answer is `−∞` instead of a finite
  large-negative. Promote lo → hi before calling (see `log1p_full`).
- **`cpow(0, w)` zero check must inspect both limbs.** `hi == 0`
  alone misclassifies `(hi=0, lo=ε)` as a true zero. See `cpowdd` —
  same pattern.
- **Large-argument trig UB.** `nearbyint(2.0·ax.hi)` overflows to `∞`
  for `|x| ≳ 1e308`; the subsequent `(long long)` cast is UB. Use
  `pi_trig_arg_too_large(ax)` / `trig_arg_too_large(ax)` — these
  guard at the precision cliff, well below the overflow point.
- **Fuzz ulp counts on cancellation witnesses are not kernel errors.**
  `fuzz.f90::check` normalizes by `|q|`; for dot-products where
  `|result| = 10⁻⁸·|inputs|`, any residual ε becomes `ε / |result|` ulp.
  Verify suspected blowups against `|inputs|` before rewriting kernels.
  The matmul dot-product "288 DD ulp" observed during precision tiering
  was the canonical example: at 100k trials a witness reports max
  12194 DD ulp, but normalizing by `|inputs|` instead of `|result|`
  shows max 0.447 DD ulp — the kernel was already at the `eps²` bound.
- **DD → DD complex transcendentals can collapse on the near-axis.**
  The textbook `casinh(z) = log(z + √(1 + z²))` loses all precision for
  `|Re z|` tiny relative to `|Im z|` because `|arg|² = 1` at the clog
  call. `catanh`'s `log((1+z)/(1−z))` has the same shape for `z` near
  the unit circle. Both were fixed by porting libquadmath's branching
  (`casinhdd`, `catanhdd` in `multifloats_math_abi_complex.inc`): a
  regime-picked `log1p` of a small positive sum replaces the
  cancelling `log`. Apply the same template before attempting TD
  internals on a complex transcendental.

## 6. Code landmarks

Entry points for contributors trying to find the subtle pieces.

| Landmark                                | What it solves                                                |
| --------------------------------------- | ------------------------------------------------------------- |
| `pi_trig_arg_too_large` (trig.inc)      | Unified range gate for `sin/cos/tan` (2^55) and `sinpi/cospi/tanpi` (2^52). |
| `dd_signbit(hi, lo)`                    | Signed-zero-aware sign test; use instead of `signbit(hi)` when branch cuts matter. |
| `dd_cross_diff(a,b,c,d)`                | `ab − cd` without catastrophic cancellation; expands to 14 scalar doubles through a TD accumulator. Marked `noinline, cold` so consumers (`cdivdd` always, `cmuldd` on cancellation) keep their hot path tight. |
| `dd_x2y2m1(x, y)`                       | `x² + y² − 1` without the `1 − a² − b²` cancellation; used by `catanh` and `cacosh` on the unit circle. |
| `reduce_pi_half` (trig.inc)             | DD Cody–Waite range reduction with `two_sum` on **both** hi and lo limbs per partial. The audit's "sloppy" add was the single worst precision leak. |
| `exp_full` (exp_log.inc)                | Three-piece CW split of `ln2`. Removes the `y = x·log2e` cancellation that caused the ~300-ulp `exp` tail. Do not replace with a single-limb split. |
| `pow_full` (exp_log.inc)                | Triple-double split: `log2(x) = e_x + log2(m)` with `e_x = ilogb(x)` exact; `y·log2(m)` bounded independent of `|e_x|`. Uses `log2_values_extra[32]` for third-limb residuals. |
| `tgamma_stirling_direct` (special.inc)  | Direct Stirling for `x ≥ 13.5`, avoiding the `exp(lgamma)` amplification. Mirrors libquadmath `tgammaq`. |
| `lgamma_positive` shift recurrence      | For `13.5 ≤ x < 25`, accumulates `Γ(x) = Γ(x+N) / (x·(x+1)·…·(x+N−1))` so Stirling runs at `x+N ≥ 25` where 13 terms converge. |
| Matmul `renorm_interval`                | Keeps compensated-FMA lo-limb bounded for large k. `ri = 8` is near-optimal; `README.md` table shows the trade-off. |
| `fmoddd` gap-aware reduction            | `gap = ilogb(r) − ilogb(y)` selects between scalar `trunc(r_hi/y_hi)` (gap ≤ 53) and DD-level `trunc(r/y)` (gap > 53). The DD path carries ~106 integer bits per step and converges even for huge gaps. Residual floor ~2.5e9 DD ulp is intrinsic to DD when gap > 106 — fixing it requires TD/QD, out of scope. |
| `casinhdd` / `catanhdd` branch schedules | Ported from libquadmath `__quadmath_kernel_casinhq` / `catanhq`. Each picks between ~8–10 regions on `(|Re|, |Im|)` and uses `log1p` of a positive sum (via `dd_x2y2m1` / `dd_cross_diff`) instead of `log(z + √(1+z²))` where the textbook form would collapse. `cacosh` / `cpow` on the unit circle share the `dd_x2y2m1` helper. |
| `cmuldd` cancellation-gated compensation | Computes the cheap 4-mul form first, then fires `dd_cross_diff` per component when the leading-limb ratio `|R|/max(|p|,|q|) < 2⁻⁴` (same shape as `cexpm1dd`'s `kCancelThresh`). 2⁻⁴ caps cheap-path rel-err at ~16 ulp_dd while firing on ~6% of uniform-quadrant inputs; tighter 2⁻¹ > halves the bench speedup, looser 2⁻⁸ lets 100+-ulp regimes through. Fuzz: 248 → 9 ulp_dd on Re, 183 → 9.6 ulp_dd on Im. ~40% bench cost — the floor is the compensated branch existing in the function body, not detector arithmetic. |
| `float64x3` + `td_mul_td` (multifloats.h) | Triple-double primitives for kernels whose DD intermediates cancel. Consumers today: `cexpm1dd` Re path via `exp_full_td` (`exp_log.inc`) + `sincos_full_td` (`trig.inc`). See `TRIPLE_DOUBLE.md` for the "constants must also be TD" rule and the `kCancelThresh` regime-split pattern. |

## 7. How to validate a precision or speed change

Every audit fix landed one-per-commit against this checklist so the
bisect history stays clean.

1. `ctest --test-dir build --output-on-failure` — all 9 must pass.
   `fortran_abi_sync`, `dd_constants_up_to_date`, and the determinism
   tests often catch drift the precision tests miss.
2. `cpp_fuzz_mpfr` — the truth for kernel-level precision; a 200-bit
   MPFR reference separates DD error from the float128-reference floor.
   Build with `-DBUILD_MPFR_TESTS=ON`.
3. `cpp_fuzz` / `fortran_fuzz` — fast, deterministic (`seed = 42`)
   fuzz against `__float128`. Good for quick iteration, weaker at the
   reference cliff.
4. `cpp_bench` / `fortran_bench` — mean of ≥ 5 runs vs libquadmath.
   Run-to-run noise is ≈ ±0.15× on most ops; smaller deltas are
   indistinguishable from noise.
5. For tolerance pins (`test/test.cc` "Tolerance sensitivity ratchet"),
   a ≥ 20× drop means the kernel got better — tighten the pin so the
   improvement is locked in. A ≥ 20× rise fails the build.

## 8. Deferred / out-of-scope work

Genuine limits of the DD representation or larger refactors that would
pay off but have not crossed the cost threshold yet.

- **Triple-double (TD) internal path for deep-tail `erfc`.** Outputs
  below `2^-969` lose bits because the DD lo limb goes subnormal. A
  separate integer-exponent path or TD internal would reach the
  float128-quality envelope. Deferred until a concrete caller asks.
- **TD-internal `log2_kernel` for sub-7-ulp `pow`.** The residual floor
  after the P7 fix is the DD accuracy of `log2_kernel` itself
  (~ulp_dd·|y|·|Lf|). Getting below ~7 ulp_dd wants a full TD rewrite
  of the kernel; large refactor, low payoff at the current spec.
- **TD input path for `csinh`/`ccosh`/`cexp` imaginary part.** Large
  `|Im z|` produces an amplification envelope of
  `2⁻¹⁰⁶·|z|/sin(reduced)` (~870 ulp at `z ≈ 47`); this is
  input-propagation, not a kernel defect. Fixing it requires DD→TD
  inputs, not a kernel change.
- **Compensated Horner (Graillat/Langlois) for Bessel's
  `neval`/`deval`.** Would drop polynomial-eval error to ~0.1 DD-ulp
  and the Bessel kernels with it. Blast radius is every
  fypp-generated caller; stopped after the T3a `reduce_pi_half` fix
  left Bessel already 4–7× better than baseline.
- **GEMM-style matmul (trans, alpha/beta, LDA).** Current kernels
  assume contiguous column-major with canonical shape. Extending
  requires a new set of panel dispatchers; tracked in `README.md` and
  `multifloats.h` as a future release item.
The `cexpm1dd` Re cancellation surface (`cos(b)·e^a = 1`, formerly
deferred) has been resolved by the TD path — see
`TRIPLE_DOUBLE.md` for the shipped regime split and the
precision / speed numbers.

---

Fixes land one-per-commit. Benchmark and fuzz deltas are the acceptance
criteria; precision is pinned via tolerance ratchets in `test/test.cc`.
