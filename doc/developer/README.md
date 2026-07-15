# Developer docs

Internal notes for people working *on* multifloats (not consumers of it).
These pages are excluded from the published Sphinx site (`doc/conf.py`); user
docs live in `doc/guides/`.

## The workflow, end to end

| Step | Where |
| ---- | ----- |
| **Configure & build** | [`doc/guides/building.md`](../guides/building.md) — requirements, the CMake options table, LTO / shared-library / dispatch flags. |
| **Test** | [`doc/guides/building.md`](../guides/building.md#tests) — every ctest target and what it covers. Build with `-DBUILD_TESTING=ON`; add `-DBUILD_MPFR_TESTS=ON` for the 3-way precision oracle. |
| **Understand the internals** | [`INTERNALS.md`](INTERNALS.md) — DD kernels, build, and test surface. [`TRIPLE_DOUBLE.md`](TRIPLE_DOUBLE.md) — the narrow triple-double path. |
| **Fix / extend the math** | [`MATH_GAPS.md`](MATH_GAPS.md) — unimplemented surface vs C23 / Fortran intrinsics. [`BOOST_COMPARISON.md`](BOOST_COMPARISON.md) + [`BOOST_OP_MATRIX.md`](BOOST_OP_MATRIX.md) — correctness/perf reference against Boost. |
| **Deploy** | [`RELEASING.md`](RELEASING.md) — versioning, cutting a tag, what the release workflow does, artifacts, dry-runs. |

## Debugging a failure

- **Precision regression** — reproduce with the fuzz targets (seeded, so runs
  are deterministic); `fuzz_*_determinism` catches non-determinism. Use
  `-DBUILD_MPFR_TESTS=ON` → `cpp_fuzz_mpfr` to separate the DD kernel's error
  from the `__float128` oracle floor. See [`INTERNALS.md`](INTERNALS.md).
- **Symbol leak / ABI regression** — `check_exported_symbols.sh` (static
  archive) and `check_shared_exports.sh` (shared `.so`) run as ctests and in the
  release; both name the offending demangled symbol.
- **Fortran ↔ C ABI drift** — `check_fortran_abi_sync.sh` asserts every
  `bind(c, name=*dd*)` matches a `MULTIFLOATS_API` entry.

## Conventions

- C++ is the source of truth for the numeric kernels; port Fortran to C++ when a
  compiler blocks optimization.
- Naming: [`NAMING-SUGGESTION.md`](NAMING-SUGGESTION.md) documents the intended
  scheme for `fsrc/multifloats.fypp`.
