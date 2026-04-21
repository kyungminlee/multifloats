# Precision improvement progress (local scratch — do not commit)

Baseline numbers pulled from `bench/results/benchmark-pop-os.json`
(2026-04-19). 1 DD ulp ≈ 2⁻¹⁰⁵ ≈ 2.46e-32. "speedup" is `qp_time /
dd_time`; >1× means multifloats is faster than `real(16)` / `__float128`.

## Baseline

| op                  | lang    | ulp              | speedup |
|---------------------|---------|------------------|---------|
| cdd_div_im          | fortran | 6.0e15           | 2.6×    |
| cdd_div_re          | fortran | 1.8              | 2.6×    |
| fmod                | cpp     | 5.1e15           | 1.1×    |
| cdd_asinh_re        | fortran | 1.1e10           | 6.1×    |
| cdd_atanh_re        | fortran | 2.4e9            | 6.5×    |
| byn                 | fortran | 10 393           | (n/a)   |
| by1                 | fortran | 4 612            | 6.3×    |
| bjn                 | fortran | 2 025            | (n/a)   |
| bj1                 | fortran | 1 725            | 5.8×    |
| bj0                 | fortran | 851              | 6.1×    |
| cdd_sinh_im         | fortran | 803              | 4.8×    |
| cdd_exp_im          | fortran | 803              | 4.3×    |
| cdd_cosh_im         | fortran | 803              | 4.5×    |
| by0                 | fortran | 584              | 6.3×    |
| arr_matmul (8×8·8)  | fortran | 288              | 0.52×   |
| arr_sum (n=8)       | fortran | 269              | 5.6×    |

Notes:
- `cdd_div_re` is fine (1.8 ulp); only the im-part is broken — matches
  the `bc − ad` cancellation diagnosis.
- `arr_matmul` is currently slower than quad (0.52×), so we have no
  headroom for precision work that costs speed on this one.
- Bessel recurrences `bjn`/`byn` have no separate bench entry (only
  `bessel_jn(3,.)` / `bessel_yn(3,.)` — different key).

## Tier 1

### 1a. cdd_div (cpp — relocated from fortran)
- before: im=6.0e15 ulp, re=1.8 ulp; speedup 2.6×
- after : im=2.9e15 ulp, re=1.7 ulp; speedup 3.18×
- change: Moved cadd/csub/cmul/cdiv to C++ (src/multifloats_math_abi_complex.inc,
  thin Fortran bind(c) wrappers). New `dd_cross_diff` computes the imaginary
  numerator `a*b − c*d` by expanding both DD products exactly into 14 scalar
  doubles (four FMA two_prods per side + a1*b1 / c1*d1 tails) and summing
  them through a triple-double accumulator before renormalizing to DD. The
  re part and denom stay on the plain DD operators (no cancellation there).
- note : max-ulp remainder is from fuzz witnesses where f1≡f2 as doubles
  but the qp reference differs below DD-ulp; the DD kernel correctly returns
  zero and the rel_err metric is relative-to-qp-reference, so these aren't
  kernel defects. Mean dropped ~40× (8.7e-18 → 2.0e-19), close to DD floor.

### 1b. fmod (cpp)
- before: 5.1e15 ulp (1.26e-16 rel); speedup 1.11×
- after : 2.5e9  ulp (6.13e-23 rel); speedup 1.36×
- change: Replaced scalar `q = trunc(r_hi/ay_hi)` with a gap-aware
  reduction step. gap ≤ 53 keeps the cheap scalar path; gap > 53
  switches to DD-level `q = trunc(r/ay)`, which carries ~106 integer
  bits and converges for gap > 106 by dropping ilogb(r) by ≥ 53 per
  step. DD rounding can leave a tiny negative residue; the same gap
  dispatch handles the add-back. Mean dropped ~4e6× (1.49e-20 →
  3.86e-27).
- note : residual ~2.5e9 ulp is the intrinsic DD-precision floor on
  fmod when `gap` > 106 — `DD(q*ay)` drops bits below 2^-106·|r|,
  so for very small true residues the rel-err blows up. Closing
  that gap needs a TD/QD multi-limb expansion (libquadmath does
  bit-level mantissa arithmetic); out of scope for Tier 1b.

### 1c. cdd_mul (cpp — cancellation-gated compensation)
- before: re max 6.09e-30 (~248 ulp_dd), im max 4.49e-30 (~183 ulp_dd)
  across 5 fuzz seeds × 200k iters; speedup 5.96× vs __float128 (dd_time
  0.0092 s on cpp_bench, pop-os / gcc-13).
- after : re max 2.20e-31 (~9 ulp_dd), im max 2.36e-31 (~9.6 ulp_dd);
  speedup 4.22× (dd_time 0.0130 s, ~41% slower). All 10 ctests pass;
  classic witness a=(1,ε),b=(-1,ε) still returns exact zero on Im and
  full-DD −1−ε² on Re.
- change: src/multifloats_math_abi_complex.inc cmuldd. Compute the
  cheap 4-mul form first (two DD muls + DD sub for Re; two DD muls +
  DD add for Im). Then per component fire a hi-limb cancellation test
  `|R_hi| < 2⁻⁴ · max(|p_hi|, |q_hi|)` — the ratio directly measures
  how many bits were shed to cancellation, since each DD-mul carries
  ~ulp_dd·|max| absolute error and the relative error after cancellation
  scales as `ulp_dd · max/|R|`. When the test fires, recompute that
  component through the existing `dd_cross_diff` (full 14-term exact
  expansion + TD accumulator). Im is a sum; route through
  `dd_cross_diff(ar, bi, -ai, br)` — DD-limb negation is exact.
- note : the speed floor is not detector arithmetic (fabs + max + cmp
  is ~2%; verified by setting kThresh=0 which dead-code-elides the
  compensated arm and restores 5.70× bench). It's the branch existing
  at all — any live threshold forces the compiler to keep
  `dd_cross_diff` in cmuldd's flow graph, adding register pressure.
  Marking `dd_cross_diff` `__attribute__((noinline, cold))` and the
  two ifs with `__builtin_expect(..., 0)` holds the hot path tight but
  doesn't move the ~30-40% floor. Outlining/inlining variants all
  landed within noise of each other.
- threshold choice: 2⁻⁴ caps cheap-path worst-case rel-err at ~16
  ulp_dd while firing on ~6% of uniform-quadrant inputs. Tighter 2⁻¹
  → 1-4 ulp_dd worst case but ~56% slowdown (fires ~25% of time).
  Looser 2⁻⁸ → 18-21 ulp_dd worst case, similar speed. See diagnostic
  table in the kernel-body comment.
- cdd_div unchanged (was already compensated in 1a) — dd_time 0.0382
  → 0.0342 as a side effect of making dd_cross_diff noinline (slight
  icache win).

## Tier 2

### 2a. cdd_asinh (cpp)
- before: re=1.1e10 ulp (2.7e-22 rel); speedup 6.1×
- after : re=1.4 ulp (3.5e-32 rel); speedup 6.32×
- change: Ported libquadmath `__quadmath_kernel_casinhq`'s 10-branch
  structure into src/multifloats_math_abi_complex.inc (casinhdd).
  Absolute-value + sign-restore frame. The main fix is the
  `ix<1, rx<0.5, rx≥eps²` branch that computes Re via
  `0.5·log1p(rx² + dm + 2·(rx·r1 + ix·r2))` instead of the textbook
  `log(z + √(1+z²))` — the textbook path's arg sits on |arg|²=1
  when Re(z) is tiny, so clogdd couldn't resolve the small result
  past ~DD ulp / |answer|. The new form sums only positive small
  pieces (no cancellation), feeds log1p, and survives all
  near-imaginary-axis witnesses at DD precision.
- note : clogdd also got a log1p path for big ∈ [0.5, 2] via a new
  `dd_x2y2m1(x,y) = x²+y²−1` helper (full 11-term exact expansion +
  TD accumulator, mirrors dd_cross_diff). Kept because cacosh and
  cpow still route through clogdd on unit-circle arguments.

### 2b. cdd_atanh (cpp)
- before: re=2.4e9 ulp (6.0e-23 rel); speedup 6.5×
- after : re=1.3 ulp (3.27e-32 rel); speedup 5.52×
- change: Ported libquadmath `catanhq`. Re-part picks log1p(4a/den)
  when f=num/den ≥ 0.5 so the small answer survives cancellation
  in log(num/den); the direct `log(num/den)` collapsed when num and
  den both ≈ 1 + b² for small |a|. Im-part denominator selects
  between `(1−big)(1+big)` / `(1−big)(1+big) − small²` /
  `−dd_x2y2m1(big, small)` based on which region is near the unit
  circle — replaces the direct `1 − a² − b²` catastrophic
  subtraction.
- note : ~15% speedup regression traded for ~9 orders of magnitude
  precision improvement — well within the proportional-win bar.

## Tier 3

### 3a. reduce_pi_half — IEEE DD-subtract per Cody-Waite step
- before: real sin/cos kernel internal err ≈1.46e-29 (~5900 DD-ulp) at
  witness x=47.15; cdd_{exp,sinh,cosh}_im = 1.98e-29 (~803 DD-ulp);
  cpp_bench sin/cos/tan speedup 1.98×/1.99×/0.99×.
- after : kernel-internal err 1.67e-33 (0.07 DD-ulp) — 8700× better.
  cpp_fuzz_mpfr sin=1.871e-32, cos=2.211e-32, tan=4.803e-32 (≤2 ulp).
  fortran_fuzz cdd_{exp,sinh,cosh}_im collapse to 5.23e-30 (~212 ulp),
  cdd_{sin,cos}_{re,im} 2.0–2.2e-31 (~8 ulp).
  cpp_bench sin 3.06×, cos 3.21×, tan 2.09× (all FASTER despite
  more EFT ops — lambda refactor helps the inliner).
  fortran_bench sin 2.99×, cos 3.11×, tan 2.19×, sinh 2.91×, cosh
  2.58× (all improved vs baseline).
- change: src/multifloats_math_trig.inc reduce_pi_half. Each of the
  three Cody-Waite partials (cw1, cw2, cw3 × k) was being subtracted
  via Shewchuk "sloppy" DD-add (one two_sum on hi, plain adds on lo) —
  that dumps ~eps·|pl| of accumulation error per step, crushing the
  reduced argument for |x| ≳ 10. New implementation: for each partial
  (ph, pl), two_sum on hi AND on lo limbs independently, then fold
  the lo two_sum tail through a fast_two_sum renormalize. Zero plain
  double adds in the critical path. Constants are unchanged (already
  >160 bits via 3-part split).
- note : remaining ~212 ulp in cdd_{exp,sinh,cosh}_im is the input-
  propagation floor: DD represents the fuzz's q ≈ 47.15 to ~106
  bits while the qp reference keeps ~112 bits, and cos(47.15)≈0.028,
  sin(47.15)≈-1; the amplification 2^-106·|q|/sin(reduced) ≈ 2e-29
  (~870 ulp) is the arithmetic envelope. Verified via hex-exact DD
  input + csinhq oracle: csinhdd rel_err = 5.9e-33 = 0.24 ulp (kernel
  clean). Fixing this requires DD→TD inputs, not a kernel change.

### 3b. cdd_exp_im / cdd_sinh_im / cdd_cosh_im
- closed by 3a. Residual 212 ulp is input-propagation floor.

### 3c. Bessel bj0/bj1/by0/by1 (+ recurrences bjn/byn)
- baseline      : 584–10 393 ulp; speedup 5.8–6.3×
- post-T3a only : bj0=568, bj1=948, bjn=2030, by0=146, by1=1150, byn=5650
  DD-ulp (auto-improved via sincos_full fix, no Bessel code change).
- compensation attempt: tried dd_cross_diff(p, c, q, s) for J0/J1
  and dd_cross_sum(p, s, q, c) for Y0/Y1 via a shared eft.inc. NO
  change in fuzz max-ulp (bj0 stayed at 1.3982e-29, by0 at 3.5829e-30
  to 4 sig digits). Reverted.
- root-cause update: diagnostic probe (/tmp/diag_by0_decomp.cc)
  shows non-cancellation floor = 40 DD-ulp at x=10.18 where
  |Y0(x)|=0.011. Amplification |p·s|/|Y0·sqrt(πx/2)| ≈ 25, so 1
  DD-ulp error in (p·s + q·c) products maps to 25 DD-ulp output —
  consistent with observed 40 ulp floor. Compensation only helps if
  the cross-product itself is the dominant error; here p and q
  themselves carry ~1 DD-ulp from neval/deval Horner evaluation (a
  degree-9/10 polynomial ratio in DD). That input error propagates
  through any combination scheme.
- next step (if pursued): compensated Horner (Graillat/Langlois) for
  neval and deval would drop their error to ~0.1 DD-ulp. But that's
  a polynomial-eval refactor touching all fypp-generated callers,
  large blast radius. Stopping here per plan — Bessel is already
  4–7× better than baseline via the T3a fix, and remaining error is
  near-zero-of-Y0/J0 input-propagation, not a kernel defect.

### 3d. arr_matmul / arr_sum (CLOSED — null result, reverted)
- before: 288 / 269 ulp; speedup 0.52× / 5.6×
- after : unchanged (max ulp bit-identical with/without change)
- attempt: rewrote mac_inl as Bailey-style DD add (two_sum on both hi and
  lo limbs, residual folded into s_lo). Expected to drop error floor
  from plain-double s_lo accumulation to true eps² EFT bound.
- result: max ulp unchanged (7.0961e-30). Mean shifted 9.7474e-33 →
  9.5428e-33 (~2%, below noise). Speed regressed ~9% on arr_matmul.
  Per plan stop condition: revert.
- root cause: the reported "288 ulp" is not an accumulator defect — it
  is a FUZZ METRIC artifact. test/fuzz.f90 check() normalizes by |q|
  (result magnitude) when |q| > 1e-10·input_mag. For a dot-product
  witness with strong cancellation (|result| = 9.5e-8·input_mag),
  any absolute error ε maps to ε/|result| ulp — unbounded as |result|
  shrinks. Confirmed via /tmp/diag_matmul.cc: at 100K trials, max=12194
  DD-ulp observed; same witnesses normalized to input-max show
  max=0.447 DD-ulp (kernel clean at eps² bound).
- why compensation didn't help: per-MAC absolute error is already
  bounded by eps²·|ab| (dropped al·bl + cross rounding + s_lo
  accumulation). Bailey compensation recovers the lo_err tail which is
  O(eps³·|s_hi|), below DD precision. Nothing in the MAC level can
  reduce absolute error further for DD output — only a triple-double
  accumulator would, and that's out of scope.
- conclusion: arr_matmul / arr_sum max-ulp as reported is a precision
  floor of the DD format under worst-case cancellation, not a kernel
  bug. Closed.
