# Installation

This page gets you from a fresh checkout to installed libraries you can link
against. For every build knob (tests, benchmarks, comparison harnesses, LTO,
shared library) see the contributor docs in
[`doc/dev/`](https://github.com/kyungminlee/multifloats/tree/main/doc/dev).

## Requirements

- CMake ≥ 3.27.
- A Fortran 2018 compiler with `REAL(KIND=16)` support and a C++17 compiler
  with `__float128` / libquadmath. On macOS this means Homebrew GCC 13/14/15
  (Apple Clang and Apple-shipped LLVM Flang are not sufficient).
- [`fypp`](https://fypp.readthedocs.io/) on `PATH` — the Fortran source is
  generated from `fsrc/multifloats.fypp` at configure time.

## Build the libraries

```sh
cmake -B build -S .
cmake --build build
```

A default build produces only the two installable libraries:

- **`libmultifloats-lto-<compiler>.a`** — the C / C++ kernels (LTO build,
  the default; the portable non-LTO build is `libmultifloats-nolto.a`). The
  API is header-only via `include/multifloats/float64x2.h`; this archive holds
  the out-of-line math bodies and the `extern "C"` `*dd` entry points.
- **`libmultifloatsf-<lto|nolto>-<compiler>.a`** — the Fortran module library
  (always compiler-tagged, since a Fortran `.mod` is version-locked).

Everything else — tests, benchmarks, MPFR / Boost comparison harnesses — is
opt-in and is **not** pulled into a default consumer build.

## Link against multifloats

Both libraries are exported as CMake packages. Prefer
`find_package(multifloats)` / `find_package(multifloatsf)` over hard-coded
`-l` flags so the compiler-tagged archive name is resolved automatically:

```cmake
find_package(multifloats REQUIRED)   # C / C++
find_package(multifloatsf REQUIRED)  # Fortran
target_link_libraries(app PRIVATE multifloats::multifloats)
```

For a direct compile, link the archive by name:

```sh
c++ -std=c++17 main.cc -lmultifloats
```

Next: [Usage](usage.md) walks through a first program in each language.
