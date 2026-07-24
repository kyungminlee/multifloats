# Codegen

Two sources in this tree are generated rather than hand-written. Each has a
generator, a committed-or-build-time output, and a ctest that fails if the
output drifts from its generator.

| Generator | Output | Committed? | Drift ctest |
| --------- | ------ | ---------- | ----------- |
| `codegen/gen_constants.py` (mpmath) | `src/dd_constants.hh` | yes | `dd_constants_up_to_date` |
| `fsrc/multifloats.fypp` (fypp) | `multifloats-{quad,noquad}.f90` | no (build-time) | `fortran_abi_sync` |

The standalone generator (`gen_constants.py`) and its `requirements.txt` live
under [`codegen/`](https://github.com/kyungminlee/multifloats/tree/main/codegen). The `.fypp` sources stay next to the
artifacts they compile into — `fsrc/multifloats.fypp` is also installed as
source for downstream consumers, and the `test/*.fypp` harnesses belong with
the test suite — so only the constants generator was relocated.

## DD constants — `gen_constants.py`

`src/dd_constants.hh` holds the double-double polynomial coefficients and
conversion constants used by the C++ kernels (`include/multifloats.h`) — each
value split into a `(hi, lo)` pair. They are computed at 60 decimal digits
(`mp.dps = 60`, well above DD's ~32) with [mpmath](https://mpmath.org/), then
emitted as `inline constexpr` arrays.

```sh
python3 codegen/gen_constants.py          # regenerate src/dd_constants.hh
python3 codegen/gen_constants.py --check   # verify without writing (exit ≠ 0 on drift)
```

Requires `mpmath` (`pip install mpmath`). The header is **committed**: a plain
build uses the checked-in copy, so contributors don't need mpmath. You rerun
the generator only when adding or retuning a constant — then commit the
regenerated header in the same change.

The `dd_constants_up_to_date` ctest runs `--check` (via `uv run --with mpmath`
when `uv` is available, else system `python3`), so editing the generator
without regenerating the header is caught as a test failure.

The **Fortran** module does not pull these generated constants — its `DD_*`
named constants are defined directly via the `DD_CONST` macro in
`fsrc/multifloats.fypp`.

## Fortran source — `multifloats.fypp`

The Fortran module is a [fypp](https://fypp.readthedocs.io/) template. At
build time it is expanded into two variants (an `ALL`-attached custom target,
so a plain `cmake --build` runs it):

| Variant | fypp flags | When it's used |
| ------- | ---------- | -------------- |
| `multifloats-quad.f90`   | (none)                 | host has `REAL(KIND=16)` |
| `multifloats-noquad.f90` | `-DMULTIFLOATS_NO_QUAD` | no `REAL(16)` support |

`FORTRAN_HAS_REAL16` (detected at configure time) selects which variant the
in-tree `multifloatsf` library compiles; a stable
`generated/multifloats.f90` copy is also produced for the scripts that grep
the expanded source.

The expanded `.f90` is **deliberately not committed** — every non-trivial
template edit would otherwise produce a spurious "regenerated" diff, and a
stale checkout could ship. Edit the `.fypp`; the build re-expands.

### ABI sync guard

`fortran_abi_sync` (`scripts/check_fortran_abi_sync.sh`) walks every
`bind(c, name='*dd*')` declaration in the expanded module and fails if the
target symbol is missing from `include/multifloats.h`. The check is
intentionally asymmetric — Fortran reimplements the arithmetic natively, so
the C header stays a superset. When you add a `bind(c)` binding, add the
matching `MULTIFLOATS_API` entry to the header in the same change.

See [Architecture](architecture.md) for the invariants these guards protect.
