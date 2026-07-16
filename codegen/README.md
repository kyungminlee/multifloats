# codegen

Code generators for multifloats. Two sources in the tree are generated rather
than hand-written; each has a generator, an output, and a ctest that fails if
the output drifts from its generator. See [`doc/dev/codegen.md`](../doc/dev/codegen.md)
for the full description.

| Generator | Output | Committed? | Drift ctest |
| --------- | ------ | ---------- | ----------- |
| `codegen/gen_constants.py` (mpmath) | `src/dd_constants.hh` | yes | `dd_constants_up_to_date` |
| `fsrc/multifloats.fypp` (fypp) | `multifloats-{quad,noquad}.f90` | no (build-time) | `fortran_abi_sync` |

## Setup

```sh
python3 -m pip install -r codegen/requirements.txt
```

## DD constants

```sh
python3 codegen/gen_constants.py          # regenerate src/dd_constants.hh
python3 codegen/gen_constants.py --check   # verify without writing (exit != 0 on drift)
```

The header is committed, so a plain build needs neither mpmath nor this script;
rerun the generator only when adding or retuning a constant, and commit the
regenerated header in the same change.

## fypp sources

The `.fypp` templates are expanded at build time by CMake, so they are **not**
run from here — they live next to the artifacts they compile into:

- `fsrc/multifloats.fypp` — the Fortran module (also installed as source for
  downstream consumers, hence kept under `fsrc/`).
- `test/bench.fypp`, `test/fuzz.fypp` — test-suite harnesses.

`fypp` (listed in `requirements.txt`) is the preprocessor CMake invokes; you do
not normally run it by hand.
