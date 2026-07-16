# Benchmark

How multifloats measures double-double (DD) throughput against native quad
precision, and how to reproduce the numbers on a new machine.

The measured tables — one flat row per op, C-ABI and Fortran-elemental
speedups side by side — live in [`benchmark/results.md`](benchmark/results.md).
This page describes the harness that produces them.

## What is measured

Two executables per language, built only when `-DMULTIFLOATS_BUILD_BENCH=ON`
(see [Configure](configure.md)):

| Binary | Measures |
| ------ | -------- |
| `cpp_bench` / `fortran_bench` | wall-clock speedup of the DD kernel vs the `__float128` / `real(16)` reference, per op. |
| `cpp_fuzz` / `fortran_fuzz` | precision (`max_rel`, `mean_rel`) over 1M adversarial inputs per op. |

The bench and fuzz *sources* live with the test suite (`test/unit/bench.cc`,
`test/unit/bench.fypp`, `test/unit/bench_abi.f90`, `test/unit/fuzz.{cc,fypp}`); the harness
under `benchmark/` only drives them and renders the results.

## The harness

Everything under [`benchmark/`](../../benchmark/) is a Python analysis harness
(no compiled code):

```
benchmark/
├── run_benchmarks.py      # build + run the executables, emit one per-system JSON
├── build_benchmark_md.py  # render results.md from one JSON
├── ops.py                 # op metadata (labels, categories, precision notes)
├── BENCHMARK.md.j2        # Jinja2 template for the results page
└── data/                  # fresh run output (git-ignored, transient)
```

Curated, committed per-system JSONs are kept under
[`benchmark/baseline/`](benchmark/baseline/); `results.md` is rendered from
one of them.

## Reproducing

From the repo root, on the target machine:

```sh
python3 benchmark/run_benchmarks.py --name "M1 Max"   # -> benchmark/data/benchmark-m1-max.json
python3 benchmark/build_benchmark_md.py doc/dev/benchmark/baseline/benchmark-m1-max.json
```

`run_benchmarks.py` needs the build toolchain (CMake, GCC, fypp);
`build_benchmark_md.py` needs only Jinja2. See
[`benchmark/README.md`](../../benchmark/README.md) for the full option list,
the JSON schema, and how to adjust which ops appear in the tables.

To commit a refreshed run, promote the JSON from `benchmark/data/` into
`doc/dev/benchmark/baseline/`, regenerate `results.md`, and commit both.

See [Architecture](architecture.md) for the precision envelope these speedups
trade against.
