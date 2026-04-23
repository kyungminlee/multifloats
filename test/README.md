# Test suite

Tests, fuzz drivers, and benchmarks for `multifloats`. Everything here is
built (but not installed) by the top-level `CMakeLists.txt` and wired
into `ctest`.

The reference for precision is `__float128` / libquadmath by default, so
GCC + libquadmath is required to build (Apple Clang lacks `__float128`).
The optional MPFR driver raises the reference to 200 bits.

## Layout

### Curated unit tests

| File | Target | ctest name | Purpose |
|---|---|---|---|
| `test.cc` | `cpp_test` | `precision_cpp` | C++ per-op precision table over a curated input list, qp reference. |
| `test.f90` | `fortran_test` | `precision_fortran_unit` | Fortran operator / intrinsic coverage, qp reference. |
| `precision.f90` | `fortran_precision` | `precision_fortran` | Broader Fortran precision sweep (arithmetic, matmul shape coverage, signed zero, inf/NaN). |

### Randomized fuzz

| File | Target | ctest name | Purpose |
|---|---|---|---|
| `fuzz.cc` | `cpp_fuzz` | `fuzz_cpp` | 1 M-iteration C++ property fuzz over scalar, complex, Bessel, and π-scaled trig kernels; reports max / mean rel-err per op. |
| `fuzz.f90` | `fortran_fuzz` | `fuzz_fortran` | Fortran counterpart; seed hard-wired to 42 for determinism. |
| `fuzz.cc` + `-DUSE_MPFR` | `cpp_fuzz_mpfr` | `precision_mpfr_cpp` | Same source as `cpp_fuzz`, recompiled with `USE_MPFR` to pull in an MPFR 200-bit reference — separates DD kernel error from the `__float128` reference floor. Built only with `-DBUILD_MPFR_TESTS=ON`. |

Determinism is enforced by two extra ctest entries (`fuzz_cpp_determinism`,
`fuzz_fortran_determinism`) that diff two back-to-back runs.

### ABI equivalence

| File | Target | ctest name | Purpose |
|---|---|---|---|
| `abi_equivalence.f90` | `fortran_abi_equivalence` | `precision_abi_equivalence` | Pins native / C-ABI DD results to the same HI limb (bit-exact) and LO within 4 ulp on `add/sub/mul/div/sqrt`. |

### Benchmarks (not wired into ctest)

| File | Target | Purpose |
|---|---|---|
| `bench.cc` | `cpp_bench` | Times every C-ABI kernel (`sindd`, `cdd_muldd`, `j0dd`, ...); qp vs DD speedup. |
| `bench.f90` | `fortran_bench` | Fortran elemental counterpart — same kernels, Fortran hidden-pointer ABI. |
| `bench_abi.f90` | `fortran_bench_abi` | Compares Fortran-native operator path to C-ABI wrapper on the 5 core ops. |

All benchmarks use a cross-rep feedback drain inside the timed region so
`-O3` cannot hoist or elide the inner loop; see the header comment in
each file.

### Shared helpers

| File | Purpose |
|---|---|
| `test_common.hh` | `to_q` / `from_q` / `q_rel_err` / `qstr` bridge between `multifloats::float64x2` and `__float128`. Shared across C++ test / fuzz / bench. |
| `test_common_mpfr.hh` | Same bridge against `mpreal` at 200 bits. Pulled in by `fuzz.cc` when compiled with `-DUSE_MPFR`. |

## Running

```bash
cmake -S . -B build -DBUILD_MPFR_TESTS=ON   # MPFR flag optional
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Filter by category:

```bash
ctest -R precision_     # curated unit tests
ctest -R fuzz_          # randomized fuzz
```

The benchmark executables are under `build/test/` and are run ad-hoc (or
through `bench/run_benchmarks.py` — see `bench/README.md`).

## Fuzz input strategy

Both `fuzz.cc` and `fuzz.f90` share the same input recipe so results can
be compared across languages:

- 10% non-finite (inf, NaN, signed zero)
- 10% close numbers (probe cancellation)
- 10% sum-near-zero (a + ~−a)
- 10% near-huge, 10% near-tiny
- 50% uniform random magnitude in 10^[−30, 30]

Per-op tolerances come from a small `is_full_dd` / `is_compound` /
default classifier: full-DD ops must stay at ~1 ulp of the ~106-bit DD
format, compound transcendentals (gamma, lgamma) loosen to ~1 ulp of
`double`, and libm-quality scalars fall somewhere in between.

See `doc/developer/INTERNALS.md` for the per-op precision envelope
the tolerance thresholds pin against, and for the deferred items that
are intentionally out of scope.
