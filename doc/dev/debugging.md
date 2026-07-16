# Debugging

A bug-fixing playbook: how to reproduce, isolate, and instrument a failure in
the DD kernels or the ABI surface. For the *catalogue* of numeric traps the
audit already uncovered — signed-zero limbs, cancellation on the near-axis,
range-guard UB — read [Architecture §6 (pitfalls)](architecture.md); this page
is about the mechanics of chasing a new one.

## Reproduce it

The fuzz drivers are **seeded** (`seed = 42`), so a failing input is
reproducible run to run. Re-run the exact driver and seed rather than hunting
for a fresh witness:

```sh
build/cpp_fuzz 1000000 42          # C++ vs __float128
build/cpp_fuzz_mpfr 10000 42       # vs 200-bit MPFR (needs -DBUILD_MPFR_TESTS=ON)
```

If a result looks non-deterministic, the `fuzz_*_determinism` ctests diff two
runs — start there before suspecting the kernel.

## Is it the kernel, or the reference?

A "regression" is often the **reference**, not the kernel. The `__float128`
oracle has a precision cliff near cancellation and near function zeros, and a
DD-normalized `rel_err = |abs| / |result|` blows up when `|result|` is tiny.
Before rewriting anything:

- Rebuild with `-DBUILD_MPFR_TESTS=ON` and check `cpp_fuzz_mpfr`. If the
  DD-vs-MPFR error is fine but DD-vs-float128 is not, the float128 leg is the
  floor — not a kernel defect.
- For dot-product / matmul witnesses, normalize by `|inputs|`, not `|result|`.
  The canonical "288 DD ulp" scare was 0.447 ulp once normalized correctly
  (see [Architecture §6](architecture.md)).

## Sanitizers

Build a `Debug` (or `RelWithDebInfo`) tree with LTO off and hidden-visibility
off so symbols and line info survive:

```sh
cmake -B build-asan -S . -DBUILD_TESTING=ON \
      -DCMAKE_BUILD_TYPE=Debug \
      -DMULTIFLOATS_USE_LTO=OFF -DMULTIFLOATS_HIDDEN_VISIBILITY=OFF \
      -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
      -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build-asan
ctest --test-dir build-asan --output-on-failure
```

UBSan is the one that pays off here: several range guards exist precisely
because an out-of-range `(long long)` cast on a huge trig argument is undefined
behavior (`trig_arg_too_large` / `pi_trig_arg_too_large` — see
[Architecture §2](architecture.md)). A UBSan hit on such a cast usually means a
guard is missing or fires too late.

Do **not** debug under `-ffast-math` — the header refuses to compile with
`__FAST_MATH__` set, because it silently collapses the DD arithmetic to plain
double.

## gdb on a DD value

A `float64x2` is two doubles; print both limbs, not the struct:

```
(gdb) p x.limbs[0]      # hi
(gdb) p x.limbs[1]      # lo
(gdb) p x.limbs[0] + x.limbs[1]   # value to double precision
```

The true value is `hi + lo` with `|lo| ≤ ½ ulp(hi)`; if `lo` is not far below
`hi` in magnitude the pair is un-normalized — usually the actual bug. For
signed-zero questions (`hi == 0 && lo != 0`) `signbit(hi)` lies; the kernels
use `dd_signbit(hi, lo)`, and so should your inspection.

## Symbol-leak and ABI-drift failures

These fail the build (and the release) rather than producing wrong numbers;
both name the offending symbol:

- **Symbol leak** — `check_exported_symbols.sh` (static archive) and
  `check_shared_exports.sh` (shared `.so`) fail if any strong global symbol
  outside the `MULTIFLOATS_API` export list plus the `multifloats::` C++ API
  leaks. A new export needs the `MULTIFLOATS_API` attribute; an internal
  helper must not have it.
- **Fortran ↔ C ABI drift** — `fortran_abi_sync` (`check_fortran_abi_sync.sh`)
  asserts every `bind(c, name=*dd*)` maps to a header entry. See
  [Codegen](codegen.md).

## Where the subtle pieces live

Architecture's ["Code landmarks"](architecture.md) table is the index of the
non-obvious kernels — `dd_cross_diff`, `reduce_pi_half`, `pow_full`, the
`casinhdd` / `catanhdd` branch schedules, the TD primitives. When a fix needs
one of them, start from that table rather than grepping cold.
