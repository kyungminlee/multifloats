# Triple-double internals

Reference for the narrow triple-double (TD) path that lives inside a
DD-first library. The TD infrastructure exists because DD-arithmetic
inherits a `~2⁻¹⁰⁵` absolute-error floor on every intermediate — fine
when the final answer scales with the input, catastrophic when the
answer is a small difference of two O(1) quantities. `cexpm1dd`'s Re
path is the one kernel that currently pays the TD cost; the primitives
are reusable for any future kernel with the same cancellation shape.

Source: `src/multifloats.hh` (primitives), `src/multifloats_math_exp_log.inc`
(TD exp), `src/multifloats_math_trig.inc` (TD sincos),
`src/multifloats_math_abi_complex.inc::cexpm1dd` (consumer),
`scripts/gen_constants.py` → `src/dd_constants.hh` (TD constants).

## 1. Why TD (and why only narrowly)

The DD representation's lo limb carries `|hi|·2⁻⁵³` worth of mantissa;
every arithmetic step injects up to one DD-ulp `≈ |result|·2⁻¹⁰⁵` of
rounding error. Two cases:

- **Magnitude-scaled answers** (`expm1`, `sin` near 0, `log1p`) —
  each term's abs err scales with the term's magnitude, which scales
  with `|answer|`, so relative error stays at ~1 DD ulp. DD is
  sufficient.
- **Cancellation between O(1) intermediates** (`Re(expm1(a+ib)) =
  e^a·cos(b) − 1` near the surface `e^a·cos(b) = 1`) — the two O(1)
  terms each carry `~2⁻¹⁰⁵` abs err; their difference is tiny, and the
  relative error of the difference blows up. No reordering of DD ops
  can escape this ceiling.

TD lifts the intermediate-arithmetic floor from `2⁻¹⁰⁵` to `~2⁻¹⁵⁸`.
That's enough to restore DD-quality relative error on the cancellation
surface when the *final* output is rounded back to DD.

**Scope boundary.** Do not TD-ify a kernel whose DD rel-err is already
bounded. `csin`/`ccos` benchmark at 2.8 / 7.3 DD ulp globally — well
inside the DD band — so lifting them to TD buys nothing. The TD path
exists for cancellation-surface kernels that measurably blow up at
high fuzz counts.

## 2. The constants-must-be-TD rule

The decisive lesson: **TD arithmetic over DD constants still clamps at
`2⁻¹⁰⁵`.** An initial `exp_full_td` built on DD `log2e` / `exp2_coefs` /
`exp2_table` / `ln2_cw1..cw3` did not help `cexpm1dd` at all — because
every constant inside the kernel imposed the DD-precision ceiling
regardless of how precise the arithmetic around them was.

Fix: promote every load-bearing constant to TD in `scripts/gen_constants.py`
and emit the extra limb as `*_lo2` alongside the existing `_hi` / `_lo`.
Current TD constants (listed here so future kernels know what is already
available before generating more):

- `log2_e_lo2`, `pi_dd_lo2`, `pi_quarter_lo2`, `inv_sqrt2_lo2`
- `exp2_coefs_lo2[16]`, `exp2_table_lo2[257]`
- `sin_taylor_lo2[13]`, `cos_taylor_lo2[13]`
- `ln2_cw4` (extends the Cody–Waite split by one limb)
- `pi_half_cw4` (same for the π/2 reduction)

The generator verifies each TD triple against MPFR at 300 bits; the
`dd_constants_up_to_date` ctest catches drift.

## 3. Primitive vocabulary

All in `multifloats::detail`. Named for what they do, not for any
particular caller. Unit-tested in
`test/test.cc::test_td_primitives` against `__float128`.

| Primitive                           | Contract                                                 |
| ----------------------------------- | -------------------------------------------------------- |
| `float64x3 {h, m, l}`               | Non-overlapping triple; `h + m + l` is the true value.   |
| `three_sum(a, b, c) → (s, t, u)`    | `a + b + c = s + t + u` exact; no normalization.         |
| `renorm3(h, m, l) → float64x3`      | Normalize a sloppy triple (magnitude-descending).        |
| `tsum3(…) / tsum3_finalize(…)`      | TSUM accumulator; pattern mirrors `dd_cross_diff`.       |
| `td_add_double`, `td_sub_double`    | TD + exact double; carries through all three limbs.      |
| `td_add_dd`, `td_add_td`            | TD + DD / TD + TD with full residue capture.             |
| `td_mul_dd`, `td_mul_td`            | 6- / 9-way scalar products through a TD accumulator.     |
| `td_to_dd`, `td_from_dd`            | Round to DD / zero-extend to TD.                         |
| `td_negate`                         | Negate all three limbs.                                  |

`td_mul_td` is the hot path: 9 scalar `two_prod`s + a magnitude-
descending compensated sum, then `renorm3`. Any derivation that can be
expressed as "TD ⊗ TD → TD ⊖ exact double → DD" will inherit
`~2⁻¹⁵⁸` residue through the arithmetic.

## 4. TD kernels

`exp_full_td(x) → float64x3` and `sincos_full_td(x) → (float64x3, float64x3)`
mirror their DD counterparts but carry a third limb through range
reduction, Horner / Taylor evaluation, and table multiplication.

Invariants both kernels preserve:

- `td_to_dd(exp_full_td(x)) == exp_full(x)` bit-exactly (fourth row of
  the Phase-C measurement table). The TD output is a strict superset
  of the DD output.
- Range reduction uses **FMA-captured** constant products. The first
  implementation wrote `-(n_float * pi_half_cw1)` directly; at
  `|n| ~ 2²⁰` the product silently rounded because `pi_half_cw1` has
  a 33-bit significand, and `sincos_full_td(10¹⁰)` degraded to
  `rel 1.6e-6`. Any `n·cw_k` product anywhere in a TD reduction must
  be captured through an FMA pair (`two_prod`), not a single multiply.

Adding a new TD kernel? The checklist:

1. Identify every constant the kernel loads; generate a TD companion
   limb in `scripts/gen_constants.py`.
2. Rewrite Horner / Taylor / table multiplies with `td_mul_td` +
   `tsum3` accumulators.
3. FMA-capture every `n·cw_k` range-reduction product.
4. Add a unit test that proves `td_to_dd(f_td(x)) == f(x)` across a
   targeted input set.
5. Add an `fuzz_mpfr` case and confirm the TD kernel sits at
   `≲ 0.3 DD ulp` — below the DD kernel's floor is the whole point.

## 5. The `cexpm1dd` Re path (current consumer)

Always compute the cheap DD half-angle form; fall through to TD only
when the two Re terms cancel enough for the DD floor to dominate.

```cpp
// DD path (fast, always run)
Re_dd = expm1(a)·cos(b) − 2·sin²(b/2)

// Cancellation detector: how many bits did the two terms agree on?
if (|Re_dd| ≥ kCancelThresh · max(|term1|, |cos_m1|))
    return Re_dd;                         // no cancellation, DD is enough

// TD path
Re_td = exp_full_td(a) ⊗ c_td ⊖ 1;        // all TD
return td_to_dd(Re_td);                   // fold back to DD
```

`kCancelThresh = 2⁻¹` is the shipped threshold. The tuning sweep
(`fuzz_mpfr` N = 10⁷, `cpp_bench` on M1 Max):

| kCancelThresh    | Re DD ulp | dd_time (s) | vs qp |
| ---------------- | --------- | ----------- | ----- |
| pure TD (no DD)  | 0.25      | 0.0499      | 1.06× |
| 2⁻⁴              | 28.3      | 0.0289      | 1.88× |
| 2⁻²              | 7.7       | 0.0298      | 1.78× |
| **2⁻¹ (ship)**   | **4.3**   | **0.0303**  | **1.80×** |
| no split (DD only, baseline) | 1078 | 0.0281 | 1.88× |

Each bit of cancellation tolerated costs ~1 DD ulp of DD-path rel-err.
The knob is the direct dial on worst-case precision; the threshold is
empirical per kernel.

Im stays pure DD — no cancellation structure there, and the existing
`exp_full(a)·sin(b)` sits at 2.6 DD ulp, bit-identical before and after
the Re rewrite. (`sincos_full_td` happens to be *more* precise than
`sincos_full` at ~0.24 vs ~1.8 DD ulp, which nudges `cdd_expm1_im` from
2.6 → 1.9 ulp as a side effect.)

## 6. When *not* to try TD

Options 2–4 from the cancellation-surface decision log remain valid
alternatives to reach for before building another TD kernel:

- **Rewrite the identity, not the arithmetic.** Some cancellation
  surfaces admit an equivalent DD-precise expression where per-term
  abs err scales with `|answer|`. The `expm1`/`log1p` pairings in
  `asinh_full` (P8) and `atanh_full` (P9) are the template — we beat
  40+ DD ulp tails with a reformulation, no TD needed.
- **Accept the band.** 184 DD ulp at `N=10⁵` is still inside the
  DD envelope for most users. A TD rewrite is warranted when a caller
  probes deep into the tail, not pre-emptively.
- **Narrow the input range.** The rest of the library refuses huge-
  argument trig (`|x| ≥ 2^55`) rather than silently returning wrong
  answers; the same pattern applies to any kernel whose DD floor is
  dictated by the input representation rather than the arithmetic.

The `cexpm1dd` Re path passed the threshold for TD because (a) the
worst case grew unboundedly with fuzz N, (b) no algebraic rewrite
recovered DD precision at the cancellation surface, and (c) the TD
infrastructure amortizes over the handful of future kernels that would
land in the same bucket.

## 7. Cost at a glance

Representative numbers for the shipped TD path (`cexpm1dd` Re),
measured on M1 Max / g++-15 -O3 at `fuzz_mpfr` N = 10⁷:

| metric                            | before TD | after TD (ship) |
| --------------------------------- | --------- | --------------- |
| `cdd_expm1_re` max_rel            | 2.65e-29  | 9.73e-32        |
| `cdd_expm1_re` DD ulp             | 1078      | 4.3             |
| `cdd_expm1_im` DD ulp             | 2.6       | 1.9             |
| `cdd_expm1` bench speedup vs qp   | 1.88×     | 1.80×           |
| `ctest` (11 tests)                | all green | all green       |

~4× improvement in the bench column is paid for a ≥ 250× improvement
in the worst-case precision column. Further TD kernels are expected
to ride similar ratios if they clear the scoping gate in section 6.
