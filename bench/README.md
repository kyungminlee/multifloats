# Benchmark automation

Scripts that keep [`BENCHMARK.md`](../doc/BENCHMARK.md) in sync with
reproducible results collected on each target machine.

## Layout

```
bench/
├── ops.py                 # op metadata (categories, approaches, prec labels)
├── BENCHMARK.md.j2        # Jinja2 template for BENCHMARK.md
├── run_benchmarks.py      # builds + runs fuzz/bench, writes per-system JSON
├── build_benchmark_md.py  # renders BENCHMARK.md from per-system JSON files
└── results/
    ├── benchmark-m1-max.json
    ├── benchmark-skylake.json
    └── benchmark-raptor-lake.json
```

The JSON files under `results/` are the single source of truth for the
tables in `BENCHMARK.md`. Each file captures one system's bench speedups
(from `fortran_bench` / `cpp_bench`) plus its fuzz `max_rel` values (from
`fortran_fuzz` / `cpp_fuzz`), together with a short system description
(CPU, OS, compiler, build flags).

## Requirements

- Python 3.9+ with [Jinja2](https://jinja.palletsprojects.com/).
  `pip install jinja2` is enough.
- A working build toolchain for the repo: CMake ≥ 3.27, GCC 13+
  (Fortran + C++), and `fypp`. See the top-level `README.md`.

## Collecting results on a new machine

From the repo root:

```bash
python3 bench/run_benchmarks.py --name "M1 Max"
```

That will:

1. Configure and build `fortran_bench`, `fortran_fuzz`, `cpp_bench`,
   `cpp_fuzz` under `build/` (reusing an existing cache if present).
2. Run each binary, parse its stdout.
3. Auto-detect CPU / OS / compiler / build info.
4. Write `bench/results/benchmark-m1-max.json`.

The filename slug is derived from `--name` (lowercased, non-alphanumeric
→ `-`). Pass `--slug` to override.

### Overriding auto-detected fields

Auto-detection is best-effort; override any field that comes out wrong:

```bash
python3 bench/run_benchmarks.py \
    --name "Skylake" \
    --cpu  "Intel Xeon family 6 model 85 (Skylake-SP / Cascade Lake), 2.8 GHz, 16 cores, AVX-512 (KVM), 22 GB" \
    --os   "Ubuntu 24.04.4 LTS" \
    --compiler "GNU Fortran / g++ 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1)" \
    --build "CMake 3.28.3, \`-O3 -flto\`, OBJECT library"
```

### Skipping parts of the run

- `--no-build`                : assume binaries already exist under `--build-dir`.
- `--skip-fortran` / `--skip-cpp` : skip one language.
- `--skip-fuzz`               : skip the slow 1M-iteration precision runs (only
  refreshes speedups).
- `--build-dir <path>`        : point at an alternative build tree.

### Runtime

On a warm toolchain:

- `fortran_bench` + `cpp_bench`: ~1–2 min.
- `fortran_fuzz`: ~5–15 min (1M iterations × many ops).
- `cpp_fuzz`: ~2–5 min.

## Regenerating `BENCHMARK.md`

After collecting (or updating) a result JSON file, render the Markdown.
The renderer is **single-architecture** — pass exactly one JSON path:

```bash
python3 bench/build_benchmark_md.py bench/results/benchmark-m1-max.json
```

This writes `doc/BENCHMARK.md` from that one file as a single flat
table (one row per op, with C-ABI and Fortran-elemental measurements
merged side-by-side in a `speedup (C / F)` column).

Pass `--stdout` to preview without touching the file, or `-o <path>` to
write elsewhere:

```bash
python3 bench/build_benchmark_md.py bench/results/benchmark-m1-max.json --stdout | head -60
python3 bench/build_benchmark_md.py bench/results/benchmark-m1-max.json -o /tmp/B.md
```

## Typical end-to-end flow

```bash
# On the target machine
python3 bench/run_benchmarks.py --name "M1 Max"
python3 bench/build_benchmark_md.py bench/results/benchmark-m1-max.json
git add bench/results/benchmark-m1-max.json doc/BENCHMARK.md
git commit -m "bench: refresh m1-max results + regenerate BENCHMARK.md"
```

## JSON schema

One file per system. All numbers are plain JSON floats; missing data
surfaces in the rendered table as `—`.

```json
{
  "system": {
    "short_name": "M1 Max",
    "cpu": "Apple M1 Max (ARM64, 10 cores)",
    "os":  "macOS 26.3 (Darwin 25.3.0)",
    "compiler": "GNU Fortran / g++ 15.2.0",
    "build":    "CMake 4.3.1, `-O3 -flto`, STATIC library"
  },
  "fortran": {
    "bench": { "add": { "n_ops": 409600, "qp_time": 1.23, "dd_time": 0.65, "speedup": 1.90 }, "...": {} },
    "fuzz":  { "add": { "n": 1000000, "max_rel": 1.5e-32, "mean_rel": 4.5e-33 }, "...": {} }
  },
  "cpp": {
    "bench": { "...": {} },
    "fuzz":  { "...": {} }
  }
}
```

`bench` keys match the label printed by the bench programs (e.g.
`"add (dd+dp)"`, `"arr_matmul (8x8*8)"`); `fuzz` keys match the
`update_stats(op, ...)` labels used in `test/fuzz.{cc,f90}` (e.g.
`"cdd_div_re"` / `"cdd_div_im"` for the split complex ops).

## Adjusting what appears in the tables

Edit [`ops.py`](ops.py). Two `Section` lists feed the renderer:
`C_SECTIONS` (C ABI surface from `include/multifloats.h`) and
`FORTRAN_SECTIONS` (elemental wrappers from `fsrc/multifloats.fypp`).
`unified_rows()` merges entries that share a `bench_key` into one row
(scope `both`); ops with no counterpart on the other side stay separate
(scope `C` or `Fortran`). Fields on each `Op`:

| Field       | Used for                                                         |
|-------------|------------------------------------------------------------------|
| `bench_key` | Key into the `bench` dict (must match the bench program output). Also the merge key between C and Fortran sections. |
| `display`   | Markdown-escaped label shown in the `op` column.                 |
| `approach`  | Prose description printed in the `approach` column.              |
| `prec`      | Precision label (`full DD` / `exact` / …). Only set on the Fortran side. |
| `fuzz`      | How to assemble the `err_dd` / `err_qp` cells: `None`, `"exact"`, a key, or a list of `(label, key)` pairs for split re/im rows. |

The header prose, precision/origin keys, and the "Notes" section live
directly in [`BENCHMARK.md.j2`](BENCHMARK.md.j2); edit the template when
those need updating.
