# Developer docs

Internal documentation for people working *on* multifloats (not consumers of
it). These pages live in the source tree and are **excluded from the published
Sphinx site** (`doc/conf.py`); the user-facing manual is in `doc/user/`.

## The workflow, end to end

| Step | Where |
| ---- | ----- |
| **Configure** | [`configure.md`](configure.md) — requirements, toolchain, the CMake options table, LTO / shared-library / dispatch flags. |
| **Build** | [`build.md`](build.md) — the compile/install steps, the codegen prerequisite, and the generated `.mod` distribution model. |
| **Codegen** | [`codegen.md`](codegen.md) — the mpmath constant generator and the fypp Fortran template, plus their drift-detection ctests. |
| **Test** | [`test.md`](test.md) — every ctest target, the MPFR 3-way oracle, and the checklist for validating a precision/speed change. |
| **Debug** | [`debugging.md`](debugging.md) — reproduce, isolate kernel-vs-reference, sanitizers, gdb on DD values, symbol/ABI failures. |
| **Benchmark** | [`benchmarks.md`](benchmarks.md) — measured DD-vs-quad speedups and precision per op. |
| **Release** | [`release.md`](release.md) — versioning, cutting a tag, what the release workflow does, artifacts, dry-runs. |

## Internals & reference

- [`architecture.md`](architecture.md) — internal design and module map: the
  load-bearing invariants, supported argument ranges, the verified precision
  envelope, designs measured and rejected, pitfalls, and code landmarks.
- [`TRIPLE_DOUBLE.md`](TRIPLE_DOUBLE.md) — the narrow triple-double path
  (`float64x3`, `exp_full_td`, `sincos_full_td`, the `cexpm1dd` regime split).
- [`MATH_GAPS.md`](MATH_GAPS.md) — unimplemented surface vs C23 `<math.h>` and
  Fortran intrinsics: candidates to implement later.
- [`BOOST_COMPARISON.md`](BOOST_COMPARISON.md) +
  [`BOOST_OP_MATRIX.md`](BOOST_OP_MATRIX.md) — correctness/perf reference
  against `boost::multiprecision::cpp_double_double`.
- [`NAMING-SUGGESTION.md`](NAMING-SUGGESTION.md) — the intended naming scheme
  for `fsrc/multifloats.fypp`.

## Conventions

- **C++ is the source of truth** for the numeric kernels; port Fortran to C++
  when a compiler blocks optimization.
- **Fixes land one-per-commit**, gated by the fuzz/bench deltas and the
  tolerance ratchets in `test/test.cc` — see [`test.md`](test.md).
