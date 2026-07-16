# Test

```sh
cmake -B build -S . -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

Keep testing on GCC — the C++ test executables need libquadmath and
`-fext-numeric-literals`, which LLVM and Intel toolchains lack.

## The ctest suite

| Test | Language | Covers |
| ---- | -------- | ------ |
| `precision_fortran` | Fortran | Targeted vs-quad checks for constructors, assignments, every arithmetic/reduction op, and edge-case sweeps. |
| `fuzz_fortran` | Fortran | 1M random pairs through every public function, with adversarial inputs. Prints a per-op precision report. |
| `precision_fortran_unit` | Fortran | Hand-written assertions for arithmetic, signed zero, NaN/Inf, classification, rounding. |
| `precision_abi_equivalence` | Fortran | Cross-checks native Fortran ops / C wrapper / `bind(c)` reimpl are bit-identical. |
| `precision_cpp` | C++ | Targeted vs-`__float128` checks. |
| `fuzz_cpp` | C++ | C++ fuzz with the same precision-report machinery. |
| `fuzz_*_determinism` | shell | Diff two fuzz runs to catch non-determinism. |
| `dd_constants_up_to_date` | Python | `codegen/gen_constants.py --check` drift detection. See [Codegen](codegen.md). |
| `fortran_abi_sync` | shell | Every `bind(c, name=*dd*)` matches a `MULTIFLOATS_API` entry. See [Codegen](codegen.md). |

The fuzz drivers are **seeded** (`seed = 42`), so runs are deterministic and
reproducible; the `fuzz_*_determinism` tests exist to catch any accidental
non-determinism.

## Optional: MPFR 3-way precision test

`-DBUILD_MPFR_TESTS=ON` enables `cpp_fuzz_mpfr`, which compares each op against
a 200-bit [MPFR](https://www.mpfr.org/) reference and reports both
`rel_err(libquadmath vs mpreal)` and `rel_err(multifloats DD vs mpreal)`,
separating the ~113-bit float128 floor from the DD kernel's own error. This is
the honest measure of kernel precision — the `__float128` reference has its own
cliff near cancellation and function zeros.

```sh
cmake -DBUILD_MPFR_TESTS=ON -S . -B build   # needs libmpfr-dev
cmake --build build --target cpp_fuzz_mpfr
build/cpp_fuzz_mpfr 10000 42   # iterations, seed
```

## Validating a precision or speed change

Every kernel change lands one-per-commit against this checklist so the bisect
history stays clean:

1. **`ctest … --output-on-failure`** — all must pass. `fortran_abi_sync`,
   `dd_constants_up_to_date`, and the determinism tests often catch drift the
   precision tests miss.
2. **`cpp_fuzz_mpfr`** — the truth for kernel-level precision (200-bit MPFR
   separates DD error from the reference floor). Build with
   `-DBUILD_MPFR_TESTS=ON`.
3. **`cpp_fuzz` / `fortran_fuzz`** — fast, deterministic fuzz against
   `__float128`. Good for quick iteration, weaker at the reference cliff.
4. **`cpp_bench` / `fortran_bench`** — mean of ≥ 5 runs vs libquadmath;
   run-to-run noise is ≈ ±0.15× on most ops. See [Benchmarks](benchmarks.md).
5. **Tolerance ratchets** (`test/test.cc` "Tolerance sensitivity ratchet") —
   a ≥ 20× drop means the kernel got better; tighten the pin to lock it in. A
   ≥ 20× rise fails the build.

The per-op precision envelope these tests defend is tabulated in
[Architecture](architecture.md); when in doubt about a regression, start
there.
