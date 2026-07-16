# multifloats vs Boost.Multiprecision `cpp_double_double`

Side-by-side comparison of the multifloats `float64x2` kernels against
`boost::multiprecision::cpp_double_double` (Boost ≥ 1.86, header
`boost/multiprecision/cpp_double_fp.hpp`). Boost's backend provides
`+`, `-`, `*`, `/`, `sqrt`, `fabs`, `exp`, `log`, `pow`, `floor`, `ceil`,
`frexp`, `ldexp`, `ilogb` natively; everything else is supplied by
Boost.Math generics on top of the backend.

The two harnesses live in `test/test_boost_dd.cc` (precision fuzz) and
`test/bench_boost_dd.cc` (timing). Both gated by
`-DMULTIFLOATS_BUILD_BOOST_COMPARE=ON` (default OFF, fetches Boost via
`FetchContent`). Both use the same `__float128` (libquadmath) reference
oracle as `cpp_fuzz` and `cpp_bench`.

For the exhaustive op-by-op matrix (every operation either side
implements, with `max_rel` and bench timings side-by-side, including
multifloats-only surfaces such as the complex DD kernels and the
π-scaled trig family) see [`BOOST_OP_MATRIX.md`](BOOST_OP_MATRIX.md).
Regenerate it with `scripts/build_op_matrix.py` after re-running the
four harnesses.

## Methodology

- **Precision**: each op runs at 1,000,000 iterations against the same
  random input distribution as `test/fuzz.cc` (10% non-finite, 10% close
  pairs, 10% cancellation, 10% near-huge, 10% near-tiny, 50% wide
  random in 10⁶/10⁻⁶ × random sign). **Both harnesses use seed 42** —
  this is critical and was the source of an earlier misleading result
  (see "Sampling-variance correction" below).
- **Speed**: `bench` style with N=1024 inputs per rep, 400/40/4 reps for
  fast/trig/very-slow, anti-DCE drain after every rep. Both legs go
  through the same workspace; only the type changes.

Hardware used for the numbers below: i7-13700H (Raptor Lake, P-cores
4.9 GHz boost), GCC 13.3, `-O3`, single-threaded.

## Sampling-variance correction

An earlier comparison run had multifloats fuzzed at seed `42` and the
boost harness defaulting to seed `0x42 = 66`. The two RNGs draw from
the same distribution but land on different points; for ops where the
DD-precision floor is dominated by **how close x lands to a nearby
zero of the function** (Bessel, sin/cos near multiples of π/2, etc.),
that's the dominant source of variance in the reported max_rel.

Initial reading: boost wins `bj0`, `bj1`, `bjn` by 4–19×. **Wrong.**
With matched seeds, only `bjn` is a real boost win. `bj0` is a tie
(8.5e-29 vs 9.4e-29) and `bj1` flips to a multifloats win (3.7e-29
vs 1.4e-28).

Lesson: when comparing two stochastic precision sweeps that bottom
out at "how lucky did the RNG get near singularities," **always pin
the seed**. Sampling variance on max_rel from 1M iterations with rare
near-zero hits can be 4–10×.

## Headline (matched seed = 42, n = 1,000,000)

**Precision** (max_rel vs `__float128`):

| group | mf wins (≥1.4× better) | boost wins (≥1.4× better) | tied / exact |
|---|---|---|---|
| Arithmetic | add, sub, add_fd | — | mul, div, sqrt, abs, neg, mul_df |
| Rounding/scale | — | — | trunc, round, floor, ceil, ldexp, scalbn, frexp.frac, logb, ilogb |
| Comparison/binary | fmod (boost broken; see below) | — | fmin, fmax, fdim, copysign, hypot |
| exp/log | exp, exp2, expm1, log, log2, log10, log1p, pow | — | — |
| Trig | sin, cos, tan, asin, acos, atan, atan2 | — | — |
| Hyperbolic | sinh, cosh, tanh, asinh, acosh, atanh | — | — |
| Special | erf, erfc, tgamma, lgamma | — | — |
| Bessel | bj1, by0, by1, byn | — | bj0, **bjn (sampling variance — see §Bessel Jn revisit)** |

**Tally: 33 multifloats wins, 0 boost wins, 22 tied/close.**

**Speed** (boost time / multifloats time, lower = mf faster):

| group | mf wins (>1.4×) | boost wins (<0.7×) | tied |
|---|---|---|---|
| Arithmetic | mul 4.3×, div 3.3×, fma 3.9× | abs 0.13×, neg 0.26× | add, sub, sqrt |
| Rounding | round 1.6× | — | trunc |
| Binary | hypot 1.6×, fmod 2.4× | copysign 0.54× | fmin, fmax, fdim |
| exp/log | all 8 ops, 3.4× to 7.8× | — | — |
| Trig | all 7 ops, 5.1× to 45.5× | — | — |
| Hyperbolic | all 6 ops, 1.9× to 20.9× | — | — |
| Special | all 4 ops, 3.4× to 13.2× | — | — |
| Bessel | bj0 4.7×, bj1 4.2×, by0 43.7×, by1 46.6×, byn 24.4× | — | bjn |

**Tally: 36 multifloats wins, 3 boost wins, 8 tied.**

## Notable findings

### Boost's transcendentals run on Boost.Math generics
The backend only natively implements `exp`, `log`, `pow` plus the
arithmetic primitives. `sin`, `cos`, `tan`, `sinh`, `cosh`, etc. all
go through Boost.Math generic algorithms parameterized by precision.
Result:
- `sin` / `cos` / `tan`: **~10⁹× less precise** than multifloats
  (1e-25 vs 1e-32) and 5–8× slower. Boost's generic trig is using
  Taylor expansion with range reduction on the backend, which doesn't
  exploit the DD structure for either precision or speed.
- `asin` / `acos` / `atan` / `atan2`: 2–3× less precise, **23–46×
  slower** (atan is 0.21 s vs multifloats 4.6 ms over 40,960 ops).

### `boost::math::fmod(cpp_double_double)` is broken
Returns rel_err on order 10¹⁰. Failure mode (from fuzz failure prints):
when |dividend| < |divisor|, returns the divisor instead of the
dividend. Routine fmod identity `fmod(x, y) = x` for `|x| < |y|` is
violated. Worth filing upstream if anyone has time.

### `fma` is 4× slower on boost
Despite `fma` being one of the operations every DD library should
have specialized, boost has no native `eval_fma` and routes through
`boost::math::fma` which builds `a*b + c` as `(a*b) + c` — losing the
fused-rounding precision benefit and doing 2 DD ops where multifloats
does one combined two_prod-and-renorm.

### Bessel Jn revisit — the "bjn 19× win" was also sampling variance

Initial reading at seed 42 had mf bjn = 5.06e-28 vs boost = 4.14e-29
(after extending boost coverage to match multifloats' n ∈ {2,3,5,8}).
Looked like a clean ~12× boost win in the only Bessel slot left.

A regime-binned probe (`test/probe_bjn.cc`) checked all four jn
regimes — forward (n ≤ x/2), forward-near (x/2 < n ≤ x), Miller-near
(x < n ≤ 2x), Miller-far (n > 2x) — across both a deterministic grid
sweep and a 1M-sample fuzz-style random sweep. Multifloats wins or
ties boost in **every regime** including near-root (|J_n(x)| < 1e-3).

Sweep across 7 seeds (1M iters each):

| seed | mf | boost | winner |
|---|---|---|---|
| 1 | 2.6e-28 | 2.0e-28 | tied |
| 2 | 2.6e-29 | 1.5e-27 | mf 57× |
| 3 | 1.3e-29 | 9.0e-29 | mf 7× |
| 42 | 5.1e-28 | 4.1e-29 | boost 12× |
| 100 | 3.3e-28 | 2.6e-28 | tied |
| 1000 | 1.0e-28 | 9.2e-29 | tied |
| 12345 | 8.8e-28 | 1.3e-28 | boost 7× |

Geometric mean ~1.5e-28 for both — a wash. Same lesson as bj0/bj1:
max_rel near-zero is dominated by how close the RNG lands to a
J_n root, and seed-42 happened to land closer for multifloats.
Across the regime sweep multifloats hits the DD precision floor
(~1e-32 mean) just like boost.

No legitimate boost precision wins remain in the Bessel column.

### `abs` / `neg` / `copysign` speed: a measurement artifact, not a
### real gap

The bench numbers above show boost winning `abs` (0.13×), `neg`
(0.26×), and `copysign` (0.54×). **This is a methodology mismatch,
not a kernel-quality gap.** The multifloats kernels (`negdd`,
`fabsdd`, `copysigndd`) are exported from `libmultifloats` with
`MULTIFLOATS_API` = `__attribute__((visibility("default")))`; the
bench executable is built without `-flto`, so each call to a
1–3-instruction kernel pays a CALL + struct-spill/reload ABI dance
(~10 cycles of glue around 1 cycle of actual work).

Boost's `cpp_double_double::operator-()` and `fabs` are
`inline constexpr` in a header — fully inlined into the BENCH loop
body, no call, no spill.

Verified by rebuilding `cpp_bench` with `-flto`:

| op | cpp_bench default | cpp_bench `-flto` | boost_dd_bench |
|---|---|---|---|
| abs | 0.0030 s | **0.0006 s** | 0.0004 s |
| neg | 0.0027 s | **0.0003 s** | 0.0007 s |
| copysign | 0.0014 s | 0.0011 s | 0.0007 s |

With matched call-site inlining, multifloats *wins* `neg` by 2.3×
(boost's `negate()` checks isnan/isinf and renormalizes; multifloats
is one `xorpd`) and ties `abs` (boost's `fabs` skips the
canonical-zero branch multifloats keeps for `(±0, ±0)` correctness).

Why the bench keeps its current methodology: the file header at
`test/bench.cc:1` calls out that this is timing through the C ABI on
purpose, because that's what a Fortran or `bind(c)` consumer
*actually* sees. The win/loss verdict for those consumers is
faithfully reported. C++ users calling the inline `multifloats::abs(x)`
header path get the LTO-inlined cost, which matches or beats boost.

### Where boost is also slow
- `fmin`, `fmax`, `fdim` are tied or boost-faster on cpp_bench-side,
  but only by 1–2× — both are well-optimized.

## Things tried and rejected

### Taylor expansion of J0 around root_1
Hypothesis: boost's `bj0` precision win comes from explicitly
factoring `(x - root_n)` via the Sterbenz pair `((x - x11/256) - x12)`,
so the rational approximation only models the smooth part. Adopting
the same structure should close our `bj0` gap.

Implementation: degree-32 Taylor expansion of J0 around `root_1`
(generated by `mpmath.taylor` at 60-decimal precision), DD-encoded
as `(c_k_hi, c_k_lo)` pairs, gated for `x ∈ (root_1 ± 1.0)`.

Result: **null on the fuzz aggregate, regression in the targeted
window.** Probed at `x = root_1 ± 1e-6`:

| offset | Hankel rel_err (current) | Taylor rel_err (experiment) |
|---|---|---|
| −1e-3 | 4.2e-30 | 4.0e-30 |
| −1e-6 | 9.3e-28 | **1.0e-26 (worse)** |
| +1e-6 | 2.1e-27 | 3.9e-27 (worse) |
| −1e-9 | 1.4e-24 | **1.0e-23 (worse)** |

Reason: 32-step Horner accumulates ~32×ulp_DD of rounding error,
while the existing Hankel form `sqrt(2/πx) · (p·c − q·s)` is just
~3 DD ops. Both forms hit the same floor `DD_floor(x) / |x − root|`
near the root because that's the precision of the **input** itself,
not the kernel — but Hankel hits it with less added noise.

The earlier "boost wins by 4×" intuition was sampling variance, not
a real algorithmic gap. Reverted in commit 7c0785b ancestor.

## When boost.math precision *does* exceed DD precision
~~For `bjn`...~~ **Retracted.** This section originally claimed
boost's 4-regime dispatcher beat multifloats' 2-regime split by 19×.
After running the regime probe and a 7-seed sweep (see "Bessel Jn
revisit" above), the gap is sampling variance, not algorithmic.
Both implementations sit at the same ~1e-28 max_rel floor near
roots, dominated by input precision rather than recurrence error.
Multifloats' forward / Miller split is sufficient at DD precision —
no dispatcher rework warranted.

## How to reproduce

```bash
cmake -S . -B build -DMULTIFLOATS_BUILD_BOOST_COMPARE=ON
cmake --build build --target boost_dd_fuzz boost_dd_bench cpp_fuzz cpp_bench -j

./build/cpp_fuzz       1000000 42  > /tmp/cpp_fuzz.txt
./build/boost_dd_fuzz  1000000 42  > /tmp/boost_fuzz.txt
./build/cpp_bench                 > /tmp/cpp_bench.txt
./build/boost_dd_bench            > /tmp/boost_bench.txt
```

Both fuzz binaries take `iterations seed` as positional args. `42` is
the multifloats default; passing it explicitly to `boost_dd_fuzz`
overrides its `0x42` default and ensures matched sampling.
