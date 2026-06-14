# Building

## Requirements

- CMake ≥ 3.27
- A Fortran 2018 compiler with `REAL(KIND=16)` support and a C++17 compiler
  with `__float128` / libquadmath. On macOS this means Homebrew GCC 13/14/15
  (Apple Clang and Apple-shipped LLVM Flang are not sufficient). The build
  pins `g++-15` / `gfortran-15` automatically; pass
  `-DCMAKE_CXX_COMPILER=...` / `-DCMAKE_Fortran_COMPILER=...` to override.
- [`fypp`](https://fypp.readthedocs.io/) on `PATH` — the Fortran source is
  generated from `fsrc/multifloats.fypp` at configure time.

```sh
cmake -B build -S .
cmake --build build
```

A default build produces only the two installable libraries (`libmultifloats.a`
and `libmultifloatsf-<compiler>.a`). Everything else — tests, benchmarks, the
BLAS smoke target, MPFR / Boost comparison harnesses — is opt-in and is **not**
pulled into a default consumer build.

## CMake options

| Option | Default | What it adds |
| ------ | ------- | ------------ |
| `-DBUILD_TESTING=ON` | OFF | C++ + Fortran test/fuzz executables, ctest registrations, and the `libblas-multifloat.a` smoke target. |
| `-DMULTIFLOATS_BUILD_BENCH=ON` | OFF | `cpp_bench`, `fortran_bench`, `fortran_bench_abi` micro-benchmarks. |
| `-DBUILD_MPFR_TESTS=ON` | OFF | `cpp_fuzz_mpfr` 3-way precision test (needs `libmpfr-dev`). Implies `BUILD_TESTING`. |
| `-DMULTIFLOATS_BUILD_BOOST_COMPARE=ON` | OFF | `boost_dd_fuzz` / `boost_dd_bench` / `bjn_probe` against `boost::multiprecision::cpp_double_double` (fetches Boost ≥ 1.89). |
| `-DMULTIFLOATS_USE_LTO=ON/OFF` | ON | LTO + fat-LTO objects on the installed archives. When ON the C++ archive is named `libmultifloats-<compiler>.a`; when OFF it is the portable untagged `libmultifloats.a`. |
| `-DMULTIFLOATS_HIDDEN_VISIBILITY=ON/OFF` | ON | `-fvisibility=hidden` so only the `extern "C" dd_*` ABI is exported. Turn OFF for full C++ symbol visibility (debugging, profiling, re-export). |
| `-DMULTIFLOATSF_INSTALL_PRECOMPILED_MOD=ON` | OFF | Additionally install a compiler-tagged precompiled `.mod` + tagged Fortran archive. By default the fypp-expanded `.f90` source is shipped and consumers compile the module with their own Fortran compiler, sidestepping `.mod` format incompatibilities. |

## Tests

```sh
cmake -B build -S . -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

| Test | Language | Covers |
| ---- | -------- | ------ |
| `precision_fortran` | Fortran | Targeted vs-quad checks for constructors, assignments, every arithmetic/reduction op, and edge-case sweeps. |
| `fuzz_fortran` | Fortran | 1M random pairs through every public function, with adversarial inputs. Prints a per-op precision report. |
| `precision_fortran_unit` | Fortran | Hand-written assertions for arithmetic, signed zero, NaN/Inf, classification, rounding. |
| `precision_abi_equivalence` | Fortran | Cross-checks native Fortran ops / C wrapper / `bind(c)` reimpl are bit-identical. |
| `precision_cpp` | C++ | Targeted vs-`__float128` checks. |
| `fuzz_cpp` | C++ | C++ fuzz with the same precision-report machinery. |
| `fuzz_*_determinism` | shell | Diff two fuzz runs to catch non-determinism. |
| `dd_constants_up_to_date` | Python | `scripts/gen_constants.py --check` drift detection. |
| `fortran_abi_sync` | shell | Every `bind(c, name=*dd*)` matches a `MULTIFLOATS_API` entry. |

### Optional: MPFR 3-way precision test

`-DBUILD_MPFR_TESTS=ON` enables `cpp_fuzz_mpfr`, which compares each op against
a 200-bit [MPFR](https://www.mpfr.org/) reference and reports both
`rel_err(libquadmath vs mpreal)` and `rel_err(multifloats DD vs mpreal)`,
separating the ~113-bit float128 floor from the DD kernel's own error.

```sh
cmake -DBUILD_MPFR_TESTS=ON -S . -B build
cmake --build build --target cpp_fuzz_mpfr
build/cpp_fuzz_mpfr 10000 42   # iterations, seed
```
