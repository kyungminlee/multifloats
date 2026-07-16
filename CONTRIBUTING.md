# Contributing to multifloats

Thanks for your interest in improving multifloats. This file is the entry
point for contributors; the detailed guides live under
[`doc/dev/`](doc/dev/index.md).

## Quick start

```sh
cmake --preset debug      # configure (Ninja, tests on) into build/debug
cmake --build --preset debug
ctest --preset debug
```

The presets (`debug`, `release`, `asan`) are defined in
[`CMakePresets.json`](CMakePresets.json). The build pins Homebrew/GNU GCC by
default — Apple Clang lacks `__float128`/libquadmath. See
[`doc/dev/configure.md`](doc/dev/configure.md) for the full option matrix.

## Developer documentation

| Topic | Guide |
| ----- | ----- |
| Configure (CMake options, toolchains) | [`doc/dev/configure.md`](doc/dev/configure.md) |
| Build & install | [`doc/dev/build.md`](doc/dev/build.md) |
| Code generation (mpmath, fypp) | [`doc/dev/codegen.md`](doc/dev/codegen.md) |
| Test | [`doc/dev/test.md`](doc/dev/test.md) |
| Debugging (sanitizers, gdb, pitfalls) | [`doc/dev/debugging.md`](doc/dev/debugging.md) |
| Benchmarks | [`doc/dev/benchmark.md`](doc/dev/benchmark.md) |
| Architecture (internal design) | [`doc/dev/architecture.md`](doc/dev/architecture.md) |
| Release process | [`doc/dev/release.md`](doc/dev/release.md) |

## Conventions

- **C++ is the source of truth** for the numeric kernels; the Fortran side
  binds to them via the C ABI. Prefer moving a kernel to C++ over fighting the
  Fortran compiler.
- **One logical fix per commit**, with a message that explains *why* the old
  behavior was wrong (see the existing history for the style).
- **`-ffast-math` is forbidden** — the header refuses to compile under
  `__FAST_MATH__` because it collapses the compensated DD arithmetic.
- Run `clang-format` (config in [`.clang-format`](.clang-format)) on C/C++
  changes and keep the generated artifacts in sync (`ctest -R up_to_date`).

## Before you send a change

1. `ctest --preset debug` is green.
2. Generated sources are up to date (`dd_constants_up_to_date`,
   `fortran_abi_sync` pass).
3. For a precision- or performance-affecting change, follow the validation
   checklist in [`doc/dev/test.md`](doc/dev/test.md).
