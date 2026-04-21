# Plan: TD internals for `cexpm1dd` Re path

Status: **draft** — written 2026-04-21 against commit on `main` (post
cdd_expm1 Im fix that brought Im from 4.3e15 → 2.2 DD ulp).

## Problem

`cdd_expm1_re` in `src/multifloats_math_abi_complex.inc:225` has
unbounded worst-case relative error, scaling roughly linearly with
sample count:

| fuzz_mpfr N | complex samples | max rel err | DD ulp |
|---|---|---|---|
| 10 000 | 100 | 2.5e-31 | 10 |
| 100 000 | 1 000 | 4.5e-30 | 184 |
| 10 000 000 | 100 000 | 2.65e-29 | 1 078 |

The "184 ulp" in `doc/BENCHMARK.md:286` is a sampling floor, not a
precision floor — the kernel will produce arbitrarily large relative
errors near the cancellation surface.

### Root cause

`Re(expm1(a+ib)) = 0` on the curve `cos(b)·e^a = 1` (for all `a > 0`
there's a `b = ±acos(e^-a)`). The current formula is

```
Re = expm1(a)·cos(b) + (cos(b) - 1)
   = expm1(a)·cos(b) - 2·sin²(b/2)        // half-angle for the 2nd term
```

Near the cancellation surface the two terms are ~0.9 each with opposite
signs; their sum is tiny. Each term is DD-accurate (abs err ~
|term|·2⁻¹⁰⁵ ≈ 2e-32). The sum inherits that absolute error, so

```
rel_err = abs_err / |Re| = O(2e-32 / |Re|)
```

which blows up as inputs approach the surface. This is fundamental to
DD-only internal arithmetic: you cannot get abs err below the
DD-ulp-of-the-intermediate floor.

## Approach: triple-double internals, localized

Compute `e^a · cos(b) - 1` using triple-double (TD, ~159-bit)
precision for the intermediates, then round back to DD at the end. No
precision loss at the subtraction, so abs err drops from ~2e-32 to
~2⁻¹⁵⁸ ≈ 3e-48, and `rel_err` returns to ~DD-ulp as long as
`|Re| ≫ 3e-48`.

Scope boundary: **do not** convert the DD library to TD wholesale. We
add a narrow, single-purpose TD path used only by `cexpm1dd`'s Re
component (Im is fine at 2.2 ulp and reuses the existing DD `exp_full`
and `sincos_full`).

### What the TD path needs

1. **TD exp** — compute `e^a` as a TD triple `(eh, em, el)` with
   `eh + em + el` close to true `e^a` at ~158 bits. We don't need a
   full TD-input TD-output exp — the input is DD, only the output is
   TD.

2. **TD cos / sin** — same: DD input, TD output. Needed because the
   cancellation comes from `cos(b)·e^a`; preserving precision in just
   `e^a` is not enough.

3. **TD × TD → TD multiply** — combine (eh, em, el) with (ch, cm, cl)
   into a TD product. 9-way scalar product expansion plus magnitude-
   descending compensated sum; pattern already demonstrated in
   `dd_cross_diff` / `dd_x2y2m1` at `multifloats_math_abi_complex.inc`.

4. **TD - 1, then round to DD.** Subtraction of an exact double from a
   TD is ~exact via a single `two_sum` on the leading limb, then
   propagating residuals. Final renormalization: TD → DD.

### Identity

```
Re(expm1(a+ib)) = TD(exp(a)) ⊗ TD(cos(b)) ⊖ 1   (as DD)
```

Im stays as-is:

```
Im(expm1(a+ib)) = exp_full(a) · sin_full(b)     (DD, unchanged)
```

Keeping the `(cos(b)-1) = -2sin²(b/2)` half-angle form as a fallback
for the `|z|` small regime — see "Regime split" below.

## Implementation steps

Ordered for independent verification at each stage.

### Step 1 — TD primitives (no domain kernels yet)

Add a minimal TD toolkit in a new `src/multifloats_td.inc` (or as
static-inline helpers in `multifloats_math_abi_complex.inc` if
genuinely one-use). Needed:

- `three_sum(a, b, c) → (s, t, u)` with `a+b+c = s+t+u` exact.
- `td_add_dd(td, dd) → td` — add a DD to a TD, renormalize.
- `td_sub_double(td, d) → td` — subtract an exact double.
- `td_mul_td(a, b) → td` — 9 scalar products + TD accumulator, final
  renormalize. Lift the TSUM pattern from `dd_cross_diff`.
- `td_to_dd(td) → dd` — renormalize and drop the third limb.

Unit-test each in isolation in `test/test.cc` against qp reference.
TD primitives should be bit-exact modulo final-limb rounding.

### Step 2 — TD exp (DD input → TD output)

Extend `exp2_from_reduced` in `src/multifloats_math_exp_log.inc` to
produce a third limb. Two places where residuals are currently dropped:

- The Horner polynomial `q = c1 + s·(c2 + …)` for `2^s - 1` runs at DD.
  Its dropped residual is ~ulp_DD ≈ 2^-105. Retaining it as a third
  limb costs a handful of FMAs per step.
- The table mul `T[k] · (1 + s·q)`: the `two_prod` of the two DD
  numbers has exactly-computable residual terms; keep them.

Deliverable: `exp_full_td(x) → (eh, em, el)`. Verified via
`exp_full_td` → round → DD must equal current `exp_full`, and TD error
vs MPFR should be ≲ 2^-158.

### Step 3 — TD sincos (DD input → TD outputs for both sin and cos)

Same shape as step 2 in `src/multifloats_math_trig.inc`. The payne-
hanek range reduction already produces a pair (rh, rl) with ~DD
precision; we need to extend it to (rh, rm, rl) for the near-π/2
cases where `cos(b)` swings through tiny values (that's where the
input-side precision matters).

Deliverable: `sincos_full_td(x) → (sh, sm, sl), (ch, cm, cl)`.

### Step 4 — `cexpm1dd` Re path

```cpp
complex64x2_t cexpm1dd(complex64x2_t z) {
  float64x2 a = from(z.re), b = from(z.im);

  // Im: unchanged.
  float64x2 ea = exp_full(a);
  float64x2 s, c;
  sincos_full(b, s, c);
  float64x2 im = ea * s;

  // Re: TD path near cancellation, DD fallback away from it.
  float64x2 re;
  if (/* regime-split test */) {
    // Current formula — fine when |Re| ≫ abs-err floor.
    float64x2 ea_m1 = expm1_full(a);
    float64x2 shalf, chalf;
    sincos_full(float64x2(0.5, 0.0) * b, shalf, chalf);
    float64x2 cos_m1 = -((shalf*shalf) + (shalf*shalf));
    re = ea_m1 * c + cos_m1;
  } else {
    // TD internal path.
    auto ea_td = exp_full_td(a);
    auto [s_td, c_td] = sincos_full_td(b);
    auto prod_td = td_mul_td(ea_td, c_td);
    auto diff_td = td_sub_double(prod_td, 1.0);
    re = td_to_dd(diff_td);
  }

  return { to(re), to(im) };
}
```

### Regime split

The TD path is only needed near the cancellation surface. Criterion
for falling through to the cheap DD formula:

```
near_cancel = (cos(b) > 0) &&
              |expm1(a)·cos(b) + (cos(b)-1)| < 2^-40 · max(|expm1(a)|, 1)
```

Threshold is tunable; `2^-40` keeps the DD path on the ~99.9% of
inputs that don't cancel, and the TD path catches the long tail where
current abs err of ~2e-32 relative to a tiny `|Re|` starts bleeding
through as visible DD-ulp. Final threshold chosen by fitting fuzz
output; placeholder for now.

## Verification

1. **Unit tests** for TD primitives in `test/test.cc` — every helper
   round-trips to DD bit-exactly.
2. **`fuzz_mpfr` at N=10M** — `cdd_expm1_re` worst case should drop
   from 1 078 DD ulp to ≲ 5 DD ulp (matching `cdd_exp_re`).
3. **`cdd_expm1_im` regression check** — must stay at 2.2 DD ulp
   (Im path unchanged).
4. **`cpp_fuzz` (qp reference)** — must continue to pass. The qp
   oracle itself has 387k qp-ulp on Re, so the fuzz tolerance tier
   (`is_reduced_dd` @ 1e-22) is insensitive; no threshold change
   expected.
5. **Bench regression** — `cexpm1` speedup vs qp will drop. Current
   1.9× → expect ~0.8–1.0× after TD path (2–3× runtime on the Re
   branch). If the regime split above is honored, the common case
   is unaffected.

## Effort estimate

- Step 1 (TD primitives + tests): ~1 day. Pattern well-established
  in `dd_cross_diff` / `dd_x2y2m1`.
- Step 2 (TD exp): ~1 day. `exp2_from_reduced` is already written to
  carry residuals in comments; promoting them to a real third limb
  is mostly plumbing.
- Step 3 (TD sincos): ~2 days. Harder because the range reduction
  for large |b| has to gain a third limb; pi_cw constants may need
  more entries. **This is the highest-risk step.**
- Step 4 (`cexpm1dd` rewrite + regime split + tuning): ~1 day.
- Regression sweep + tuning: ~1 day.

Total: **~1 week, single-engineer**, plus a review pass on the TD
sincos change since it interacts with the hot-path kernels.

## Risks / open questions

- **Scope creep on TD sincos.** Once TD sincos exists, there's
  pressure to extend it to other complex transcendentals (csin/ccos
  have similar cancellation structure on their own surfaces). Resist
  unless measured — `csin/ccos` benchmark at 2.8 / 7.3 ulp today,
  well inside DD band, so no motivation to convert.
- **Code size.** TD primitives are not free. Budget ~200 lines of
  `.inc` and ~30 KB of compiled code. Compare to current
  `multifloats_math_abi_complex.inc` at ~700 lines.
- **Threshold tuning is empirical.** The regime split criterion will
  need iteration. If the chosen threshold is too tight, we eat TD
  cost on benign inputs; too loose, we leak DD-precision worst cases
  back into the measurement.
- **MPFR oracle coverage.** `fuzz_mpfr` currently samples `cexpm1`
  at ~0.01 of its rate (every 100 iters, 1/1 of complex draws). To
  confidently measure ~5 DD ulp we need ~1M complex samples; consider
  raising the complex-draw rate in `fuzz_mpfr.cc` or adding a
  dedicated precision-cexpm1 driver similar to `find_worst.cc`.
- **Does TD genuinely fix abs err or just defer it?** Input `b` is
  DD, so the true value of `cos(b_DD)` has absolute error ≈
  |cos(b_DD)|·2⁻¹⁰⁵ inherent to the DD representation — independent of
  what precision we compute cos in. What TD buys us is that the
  *rounding error* of the cos/exp operation is pushed down to 2⁻¹⁵⁸
  instead of 2⁻¹⁰⁵. The input-precision floor is still 2⁻¹⁰⁵ relative
  to cos(b), which translates to abs err ~|cos·exp|·2⁻¹⁰⁵ ≈ 2⁻¹⁰⁵ at
  the cancellation surface. **This may cap the achievable
  improvement at ~1 DD ulp, not below it.** Confirm with a prototype
  before committing to the full TD sincos implementation.

## Decision gate

Before step 3 (the expensive part), prototype steps 1–2 plus a
*mock* TD cos (TD exp + fake TD cos = DD cos + zero third limb) and
measure cdd_expm1_re on N=1M samples. If the result is already within
1 DD ulp of full-DD, the open question above is answered negatively
and we skip TD sincos entirely — TD exp alone is enough. If not,
proceed to step 3.

## Baseline (2026-04-21, M1 Max, g++-15 -O3)

Locked in **after** the Im fix (ea from `1 + expm1(a)` → `exp_full(a)`
in `multifloats_math_abi_complex.inc:cexpm1dd`). Every step of the
plan below is measured against these numbers; Re is expected to
improve, Im must not regress.

### `fuzz_mpfr` — DD vs mpreal @ 200 bits

Reported as `max_rel`, converted to DD ulp via `ulp = max_rel / 2.46e-32`.

| op | N=10k (std) | N=100k | N=1M | N=10M |
|---|---|---|---|---|
| `cdd_expm1_re` | 2.5e-31 (10 ulp) | 4.5e-30 (184 ulp) | 4.5e-30 (184 ulp) | 2.65e-29 (1 078 ulp) |
| `cdd_expm1_im` | — | — | 5.465e-32 (2.2 ulp) | 6.434e-32 (2.6 ulp) |
| `expm1` (scalar, ref) | — | — | 3.77e-32 (1.5 ulp) | — |

Re worst case grows ~linearly with sample count; Im is flat at the
full-DD band (matches `cdd_exp_im` at 2.2 ulp) and is not a target of
this plan — recorded here to guard against regression only.

### `cpp_fuzz` — DD vs libquadmath (qp)

From `bench/results/benchmark-m1-max.json` (`c.fuzz` block, MPFR-run):

| op | DD-vs-MPFR (max_rel) | qp-vs-MPFR (max_rel) |
|---|---|---|
| `cdd_expm1_re` | 4.536e-30 | 7.455e-29 |
| `cdd_expm1_im` | 5.465e-32 | 2.5e-34 |

The qp Re is ~16× worse than DD Re — the libquadmath oracle itself
drifts on this op, which is why `fuzz.cc` (qp-referenced) places
cexpm1 in the `is_reduced_dd` tier @ 1e-22 tolerance. Any Re
improvement from this plan is measured against MPFR, not qp.

### `cpp_bench` — DD time vs qp time

From `bench/results/benchmark-m1-max.json` (`c.bench`):

| op | n_ops | dd_time (s) | qp_time (s) | speedup |
|---|---|---|---|---|
| `cdd_expm1` | 40 960 | 0.0281 | 0.0528 | 1.88× |
| `cdd_exp` (ref) | 40 960 | 0.0206 | 0.0728 | 3.54× |
| `expm1` (scalar) | 40 960 | 0.0061 | 0.0334 | 5.49× |

Runtime reference: `cdd_expm1` currently costs ~1.36× a `cdd_exp`
call (from the extra `exp_full` added in the Im fix). TD path will
push this higher on the fraction of inputs that take it. Regime-
split threshold will be chosen so the common DD path dominates.

### Success thresholds (restatement against these numbers)

- Re at N=10M: 1 078 ulp → **≲ 5 ulp** (target), **≲ 50 ulp**
  (acceptable — still an 20× win).
- Im at N=10M: 2.6 ulp → **≤ 2.6 ulp** (no regression).
- `cdd_expm1` bench speedup: 1.88× → **≥ 1.0×** (acceptable:
  DD-faster-than-qp preserved; TD branch amortized over regime
  split).

## Decision-gate outcome (2026-04-21)

**Result: negative.** Step 1 (TD primitives) and step 2 (TD exp) were
implemented and the prototype was measured. TD exp + mock TD cos did
not improve cdd_expm1_re — it *regressed* to 6.24e-27 (~250 000 DD
ulp) at N=1M, more than 1 000× worse than the half-angle-DD baseline.
Risk #5 materialized.

### What was built (and remains in-tree)

Re-usable, unit-tested TD infrastructure lives in
`src/multifloats.hh` (`multifloats::detail`):

- `float64x3` struct.
- `three_sum`, `renorm3`, `tsum3`, `tsum3_finalize` primitives.
- `td_add_double`, `td_sub_double`, `td_add_dd`, `td_mul_dd`,
  `td_mul_td`, `td_to_dd`, `td_from_dd`.

13 unit tests in `test/test.cc::test_td_primitives` exercise each
primitive against qp. Status: all pass, max_rel 1.23e-32.

`exp_full_td(x) → float64x3` in `src/multifloats_math_exp_log.inc`
was added but **is not currently wired into any kernel**. Left in
place as a working demonstration of the TD pipeline for future
reference; callers: none.

### Measurements (M1 Max, g++-15 -O3, N=1M)

| variant | cdd_expm1_re max_rel | cdd_expm1_im max_rel |
|---|---|---|
| baseline (half-angle DD) | **4.536e-30** (184 ulp) | 5.465e-32 (2.2 ulp) |
| TD exp · mock-TD cos − 1 (prototype) | 6.241e-27 (~254k ulp) | 5.465e-32 |
| pure DD with identity `e^a·cos(b) − 1` | 2.957e-26 (~1.2M ulp) | 5.465e-32 |
| `exp_full_td` |> `td_to_dd` fed into old identity | 4.536e-30 | 5.465e-32 |

The fourth row confirms `exp_full_td` round-trips to `exp_full` at DD
precision — TD exp itself is correct. The regression in rows 2–3
comes from the *identity* `e^a·cos(b) − 1`, which at equal arithmetic
precision is ~1 000× less precise than the half-angle form
`expm1(a)·cos(b) − 2sin²(b/2)` near the cancellation surface.

### Why the TD path did not help

The half-angle form owes its precision to *magnitude-scaled*
absolute error per term:

- `expm1(a)` is small when |a| is small; its DD abs err is
  ~|expm1(a)|·2⁻¹⁰⁵, scaling *with* the physical size of the
  quantity.
- `−2sin²(b/2)` is small when |b| is small; likewise abs err
  ~|sin²(b/2)|·2⁻¹⁰⁵.
- Their sum's abs err is bounded by max-term·2⁻¹⁰⁵ — which stays
  small when *both* terms are small.

The TD identity `e^a·cos(b) − 1` computes two **O(1)** quantities
first, then subtracts 1. The abs err of `e^a·cos(b)` before the
subtraction is ~|e^a·cos(b)|·ε_TD + |cos|·ε_cos, where ε_cos is the
precision of the cos representation. With *mock* TD cos
(DD-precision), ε_cos = 2⁻¹⁰⁵, so the subtraction inherits a fixed
~2⁻¹⁰⁵ abs err regardless of how small |Re| becomes — the same
ceiling as DD, with none of the per-term scaling.

### Would full TD sincos (step 3) fix this?

Partially. Replacing `mock_cos_td` with a TD cos that carries its
true 159-bit value would drop ε_cos from 2⁻¹⁰⁵ toward 2⁻¹⁵⁸ —
*provided* every constant inside the sincos kernel is also TD. But
the sincos kernel relies on DD constants (sincos coefficients, the
Payne-Hanek π-split constants `pi_cw1/cw2/cw3`, the DD `pi_dd`). Each
caps the output at ~2⁻¹⁰⁵. So a "TD sincos" built over DD constants
will still carry ε_cos ≈ 2⁻¹⁰⁵ and the cancellation-surface abs err
floor would not drop.

The analogous limitation applies to `exp_full_td` itself: `log2(e)`,
the `exp2_coefs`, the 257-entry `exp2_table`, and the Cody-Waite
`ln2_cw1/cw2/cw3` constants are all DD. The TD Horner and TD table-
mul eliminate *step-wise rounding* (each was ~2⁻¹⁰⁶), but the
constants impose a 2⁻¹⁰⁵ constant-precision floor that upper-bounds
the whole pipeline.

### Options for a real fix

1. **Upgrade every constant to TD.** (Chosen — see Phase A/B/C
   sections below.) Pre-compute 159-bit representations of `log2(e)`,
   `exp2_coefs[1..9]`, `exp2_table[257]`, `ln2_cw1..cw3` (or add
   `cw4`), sincos coefficients, `pi_dd`, and the Payne-Hanek π-split.
   Significant effort (~a week on top of step 3); payoff is genuine
   2⁻¹⁵⁸ precision throughout, which *should* drop cdd_expm1_re to
   ≲ 1 DD ulp. Result: dropped to ≲ 0.25 DD ulp.

2. **Rewrite the identity, not the arithmetic.** Keep DD throughout,
   but find an expression for Re whose per-term DD abs err scales
   with |Re| rather than max(|e^a|, |cos|). E.g.
   `Re = (cos(b) − 1) · e^a + (e^a − 1)` (rearranged),
   or a Taylor expansion around the cancellation surface `cos(b)·e^a
   = 1` parameterized by distance from that surface. Hard to derive;
   doesn't always exist.

3. **Accept the current baseline.** 184 DD ulp at N=10^5 is still
   within the DD-band for most users. The 1078-ulp figure at N=10^7
   is a sampling worst case on a multi-decade-wide input
   distribution, not typical.

4. **Back out the TD infrastructure.** `git reset` to
   `pre-expm1-optim`; drop the unused TD primitives, `exp_full_td`,
   and their tests. Nothing in main actually needs them today.

### Current state of the tree (post-gate — superseded by Phase A/B/C)

At the time the gate was recorded, `cexpm1dd` was reverted to the
half-angle DD baseline and the TD primitives / `exp_full_td` sat
inert. Phase A/B/C below then implemented option 1; see the final
"Tree state after Phase C" for the post-option-1 layout.

## Phase A / B / C — full TD path with TD constants (2026-04-21)

**Option 1 accepted.** After the decision-gate negative, every `exp` /
`sincos` constant was lifted to TD precision and the kernels rebuilt.
The `cexpm1dd` Re path now runs on `exp_full_td · sincos_full_td − 1`
round-tripped to DD.

### Phase A — TD exp with TD constants

- **A.1**: extended `scripts/gen_constants.py` with `to_td`,
  `scalar_td`, `array_td`, and `cw4`. Promoted `log2_e`, `pi_dd`,
  `inv_sqrt2`, `pi_quarter`, `exp2_coefs`, `exp2_table`, `sin_taylor`,
  `cos_taylor` to TD; extended `ln2_cw3` → `ln2_cw4` and `pi_half_cw3`
  → `pi_half_cw4`. Existing DD limbs (`_hi`, `_lo`) are byte-for-byte
  unchanged. Generator now verifies `max TD err 3.165e-49` (~2⁻¹⁶¹).
- **A.2**: rewrote `exp2_td_from_reduced` with TD Horner on 15 TD
  coefficients and TD table-mul; rewrote `exp_full_td` CW reduction
  with TSUM over 4 TD-split `ln2_cw` terms and a TD `log2(e)` multiply.
  Verified `td_to_dd(exp_full_td(a)) == exp_full(a)` at baseline
  precision via bisect in `cexpm1dd`.
- **A.3**: re-ran the decision gate (TD exp + mock TD cos). Still
  regressed — 6.24e-27 at N=1M, unchanged from the constants-
  untouched attempt. This pinned the remaining floor on *cos*
  precision, not exp. → proceed to Phase B.

### Phase B — TD sincos with TD constants

- Added `reduce_pi_half_td`, `sin_kernel_td`, `cos_kernel_td`,
  `sincos_eval_td`, `sincos_full_td` in `multifloats_math_trig.inc`.
  Added `td_add_td`, `td_negate` primitives.
- **First-cut bug**: `reduce_pi_half_td` computed
  `-(n_float * pi_half_cw1)` directly instead of via an FMA pair.
  `pi_half_cw1` has 33-bit significand; for |n| < 2²⁰ the product is
  exact in double, but the `trig_arg_too_large` cap allows |n| up to
  2⁵⁴, at which point the product silently rounds. Sampled regression:
  sincos_td(1000) → rel 3e-14; sincos_td(10¹⁰) → rel 1.6e-6. Fix:
  universal FMA capture for every `n·cw_k`, k = 1…4.
- After the fix, `sincos_full_td` measures sincos_s `5.929e-33`,
  sincos_c `5.909e-33` @ N=1M (~0.24 DD ulp each — actually *beats* the
  DD `sincos_full` kernel which sits at ~1.8 DD ulp).

### Phase C — wire full TD into cexpm1dd Re

```cpp
float64x3 ea_td = exp_full_td(a);
float64x3 s_td, c_td; sincos_full_td(b, s_td, c_td);
float64x3 prod_td = td_mul_td(ea_td, c_td);
float64x3 diff_td = td_sub_double(prod_td, 1.0);
float64x2 re = td_to_dd(diff_td);
float64x2 im = td_to_dd(ea_td) * td_to_dd(s_td);  // Im stays DD
```

### Results (2026-04-21, M1 Max, g++-15 -O3)

`fuzz_mpfr` — DD vs mpreal @ 200 bits (values in `max_rel`):

| op | baseline N=10M | TD path N=10M | improvement |
|---|---|---|---|
| `cdd_expm1_re` | 2.65e-29 (1078 ulp) | **6.160e-33 (~0.25 ulp)** | **≥ 4000×** |
| `cdd_expm1_im` | 6.434e-32 (2.6 ulp) | **4.573e-32 (1.9 ulp)** | 1.4× |

The TD Re worst case **stops growing with N** — the cancellation-
surface tail is fully closed. Im improves modestly as a side effect
of `sincos_full_td` being more precise than `sincos_full`.

`cpp_bench` cost (before regime split):

| op | baseline dd_time | TD path dd_time | vs qp |
|---|---|---|---|
| `cdd_expm1` | 0.0281 s | 0.0499 s (+1.78×) | 1.06× (was 1.88×) |

### Regime split — pay TD cost only on cancellation

Always compute the DD half-angle form; only fall through to the TD
pipeline when `|Re_DD| / max(|term1|, |cos_m1|) < kCancelThresh`. Each
bit of DD-path cancellation costs ~1 DD ulp, so the threshold is the
direct knob for worst-case DD rel-err. Tuning sweep (fuzz_mpfr N=10^7,
cpp_bench on M1 Max):

| kCancelThresh | Re DD ulp | dd_time (s) | vs qp | Re interpretation |
|---|---|---|---|---|
| *(no split — pure TD)* | **0.25** | 0.0499 | 1.06× | TD pipeline on every call |
| 2⁻⁴ | 28.3  | 0.0289 | 1.88× | 4-bit cancellation threshold |
| 2⁻² | 7.7   | 0.0298 | 1.78× | 2-bit cancellation threshold |
| **2⁻¹ (shipped)** | **4.3** | **0.0303** | **1.80×** | **1-bit cancellation threshold** |
| *(baseline, no TD at all)* | 1078  | 0.0281 | 1.88× | reference |

Shipped: **kCancelThresh = 2⁻¹**. Sits under the plan's ≤5 DD-ulp
ideal while adding only ~8% runtime vs the pure-DD baseline. The TD
path fires on ~1–2% of fuzz inputs (those close to the cancellation
surface), which is empirically enough to cap the worst case at the
desired band without reinflating the runtime.

### Success-threshold check (with regime split @ 2⁻¹)

| target | result | verdict |
|---|---|---|
| Re @ N=10M ≲ 5 DD ulp (ideal), ≲ 50 ulp (acceptable) | 4.3 ulp | ✓ under the ideal |
| Im @ N=10M ≤ 2.6 ulp (no regression) | 2.6 ulp | ✓ |
| `cdd_expm1` speedup vs qp ≥ 1.0× | 1.80× | ✓ (baseline was 1.88×) |
| ctest 11/11, cpp_fuzz failures=0, cpp_fuzz_mpfr failures=0 | all green | ✓ |

### Tree state after Phase C

- `src/multifloats.hh`: TD primitives (`float64x3`, `three_sum`,
  `renorm3`, `tsum3`, `tsum3_finalize`, `td_negate`, `td_add_double`,
  `td_sub_double`, `td_add_dd`, `td_add_td`, `td_mul_dd`, `td_mul_td`,
  `td_to_dd`, `td_from_dd`) — all in `multifloats::detail`, unit-tested.
- `src/multifloats_math_exp_log.inc`: `exp_full_td`,
  `exp2_td_from_reduced`.
- `src/multifloats_math_trig.inc`: `reduce_pi_half_td`,
  `sin_kernel_td`, `cos_kernel_td`, `sincos_eval_td`, `sincos_full_td`.
- `src/multifloats_math_abi_complex.inc::cexpm1dd`: TD Re path wired.
- `scripts/gen_constants.py`: `scalar_td`, `array_td`, `cw4` kinds.
- `src/dd_constants.hh`: auto-regenerated; `log2_e_lo2`,
  `exp2_coefs_lo2[16]`, `exp2_table_lo2[257]`, `sin_taylor_lo2[13]`,
  `cos_taylor_lo2[13]`, `pi_dd_lo2`, `pi_quarter_lo2`, `inv_sqrt2_lo2`,
  `ln2_cw4`, `pi_half_cw4`.

Plan closed.
