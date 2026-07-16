# Configure

Everything CMake decides before a single object is compiled: toolchain,
feature flags, and the generated-artifact guards. For the compile/link step
itself see [Build](build.md).

## Requirements

- **CMake ≥ 3.27.**
- **A Fortran 2018 compiler** with `REAL(KIND=16)` support **and a C++17
  compiler** with `__float128` / libquadmath. In practice this is GCC. On
  macOS it means Homebrew GCC 13/14/15 — Apple Clang and Apple-shipped LLVM
  Flang lack `__float128` / `REAL(KIND=16)` and are not sufficient. The build
  pins `g++-15` / `gfortran-15` automatically; override with
  `-DCMAKE_CXX_COMPILER=…` / `-DCMAKE_Fortran_COMPILER=…`.
- **[`fypp`](https://fypp.readthedocs.io/) on `PATH`** — required at configure
  time (`find_program(FYPP_EXECUTABLE fypp REQUIRED)`); it expands the Fortran
  source. See [Codegen](codegen.md).
- **`mpmath`** (Python) only if you edit the constant generator — see
  [Codegen](codegen.md). Not needed for a plain build.

Non-GCC toolchains (LLVM, Intel) can build the two libraries but **cannot**
build the C++ test executables, which need libquadmath and GCC's
`-fext-numeric-literals`. Keep testing on GCC.

## Configure a build

```sh
cmake -B build -S .
```

A default configure builds only the two installable libraries. Everything
else is opt-in through the options below.

## CMake options

| Option | Default | What it adds |
| ------ | ------- | ------------ |
| `-DBUILD_TESTING=ON` | OFF | C++ + Fortran test/fuzz executables and ctest registrations. |
| `-DMULTIFLOATS_BUILD_BENCH=ON` | OFF | `cpp_bench`, `fortran_bench`, `fortran_bench_abi` micro-benchmarks. See [Benchmarks](benchmarks.md). |
| `-DBUILD_MPFR_TESTS=ON` | OFF | `cpp_fuzz_mpfr` 3-way precision test (needs `libmpfr-dev`). Implies `BUILD_TESTING`. |
| `-DMULTIFLOATS_BUILD_BOOST_COMPARE=ON` | OFF | `boost_dd_fuzz` / `boost_dd_bench` / `bjn_probe` against `boost::multiprecision::cpp_double_double` (fetches Boost ≥ 1.89). See [BOOST_COMPARISON.md](BOOST_COMPARISON.md). |
| `-DMULTIFLOATS_USE_LTO=ON/OFF` | ON | LTO + fat-LTO objects on the installed archives. When ON the C++ archive is named `libmultifloats-lto-<compiler>.a`; when OFF it is the portable compiler-agnostic `libmultifloats-nolto.a`. |
| `-DMULTIFLOATS_HIDDEN_VISIBILITY=ON/OFF` | ON | `-fvisibility=hidden` so only the `extern "C" dd_*` ABI is exported. Turn OFF for full C++ symbol visibility (debugging, profiling, re-export). |
| `-DMULTIFLOATS_MM_DISPATCH=ON/OFF` | ON | Compile the matmul kernels twice (a portable baseline + an AVX2+FMA copy that vectorizes to packed `vfmadd…pd`) and pick one by CPUID at runtime — a single binary that uses packed FMA where the CPU supports it and runs anywhere otherwise. Only effective on x86-64 GCC (only GCC's vectorizer packs this loop); a no-op on Clang/icx, AArch64, and Apple. Turn OFF for a single baseline build. |
| `-DMULTIFLOATS_BUILD_SHARED=ON` | OFF | Also build a versioned shared library (`libmultifloats.so.N` / `libmultifloats.N.dylib`) from the same objects as the static archive, installed with the usual soname + dev symlinks. Its `SOVERSION` tracks `MULTIFLOATS_ABI_VERSION`; the ELF build links with full RELRO + `BIND_NOW` and `--no-undefined`. With `MULTIFLOATS_HIDDEN_VISIBILITY=ON` (the default) its dynamic symbol table exports only the `extern "C"` C ABI — for direct `-lmultifloats` / `dlopen` consumers (`find_package` still uses the static archive). The release ships one prebuilt shared object per platform. |
| `-DMULTIFLOATSF_INSTALL_PRECOMPILED_MOD=ON` | OFF | Additionally install a compiler-tagged precompiled `.mod` + tagged Fortran archive. By default the fypp-expanded `.f90` source is shipped and consumers compile the module with their own Fortran compiler, sidestepping `.mod` format incompatibilities. See the `.mod` note in [Build](build.md). |

## Build type

`CMAKE_BUILD_TYPE` defaults to `Release`; the cache exposes
`Debug` / `Release` / `RelWithDebInfo` / `MinSizeRel`. Use `Debug` (or
`RelWithDebInfo`) when running under a sanitizer or gdb — see
[Debugging](debugging.md).

## Configure-time guards

Two source-only ctests are registered at configure time (only when
`BUILD_TESTING=ON`) to catch committed generated artifacts drifting from
their generators: `dd_constants_up_to_date` and `fortran_abi_sync`. Both are
described in [Codegen](codegen.md) and [Test](test.md).
