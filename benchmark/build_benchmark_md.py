#!/usr/bin/env python3
"""Render ``BENCHMARK.md`` from one per-system JSON file produced by
``benchmark/run_benchmarks.py``.

    python3 benchmark/build_benchmark_md.py doc/dev/benchmark/baseline/benchmark-m1-max.json
    python3 benchmark/build_benchmark_md.py doc/dev/benchmark/baseline/benchmark-m1-max.json --stdout
    python3 benchmark/build_benchmark_md.py doc/dev/benchmark/baseline/benchmark-m1-max.json -o /tmp/BENCHMARK.md

The output is a single-architecture flat table — one row per op,
merging C-ABI and Fortran-elemental measurements side-by-side.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, StrictUndefined

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ops as ops_mod  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
TEMPLATE_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT = REPO_ROOT / "doc" / "dev" / "benchmark" / "results.md"

EMDASH = "—"  # —
TIMES = "×"   # ×


def _finite(v) -> bool:
    return isinstance(v, (int, float)) and math.isfinite(v)


_DD_ULP = 2.0 ** -105  # ~2.46e-32; 1 ulp of a ~106-bit DD significand
_QP_ULP = 2.0 ** -112  # ~1.93e-34; 1 ulp of a 113-bit __float128 significand


def _format_max_rel(v: float) -> str:
    if not _finite(v) or v < 0.0:
        return EMDASH
    if v == 0.0:
        return "exact"
    u = v / _DD_ULP
    if u < 10:
        return f"{u:.1f} ulp"
    return f"{u:.0f} ulp"


def _format_qp_err(v: float) -> str:
    if not _finite(v) or v < 0.0:
        return EMDASH
    if v == 0.0:
        return "exact"
    u = v / _QP_ULP
    if u < 10:
        return f"{u:.1f} qp ulp"
    return f"{u:.0f} qp ulp"


def _format_speedup_number(v: float) -> str:
    if v >= 10.0:
        return f"{int(round(v))}{TIMES}"
    if v >= 1.0:
        return f"{v:.1f}{TIMES}"
    return f"{v:.2f}{TIMES}"


def _format_speedup(v) -> str:
    if not _finite(v) or v <= 0.0:
        return EMDASH
    s = _format_speedup_number(float(v))
    return f"**{s}**" if v >= 2.0 else s


def _err_cell(spec, fuzz_data: dict, field: str, formatter) -> str:
    """Generic err-cell renderer parametrised on (re, im) split or single."""
    if spec is None:
        return EMDASH
    if spec == "exact":
        return "exact"
    if isinstance(spec, str):
        entry = fuzz_data.get(spec)
        if not entry or not _finite(entry.get(field)):
            return EMDASH
        return formatter(entry[field])
    parts: list[str] = []
    for label, key in spec:
        entry = fuzz_data.get(key)
        if entry and _finite(entry.get(field)):
            parts.append(f"{formatter(entry[field])} ({label})")
    if not parts:
        return EMDASH
    return " / ".join(parts)


def err_dd_filter(row: ops_mod.UnifiedRow, system: dict) -> str:
    """DD-vs-MPFR err cell, in DD ulps."""
    fuzz_data = (system.get("c") or {}).get("fuzz") or {}
    return _err_cell(row.fuzz, fuzz_data, "max_rel", _format_max_rel)


def err_qp_filter(row: ops_mod.UnifiedRow, system: dict) -> str:
    """quadmath-vs-MPFR err cell, in qp ulps. Same fuzz dict, different field."""
    fuzz_data = (system.get("c") or {}).get("fuzz") or {}
    return _err_cell(row.fuzz, fuzz_data, "max_qp", _format_qp_err)


def speedup_pair_filter(row: ops_mod.UnifiedRow, system: dict) -> str:
    """`<C-side> / <Fortran-side>` speedup pair (qp_time / dd_time)."""
    def lookup(lang: str, key: str | None) -> str:
        if key is None:
            return EMDASH
        bench_data = (system.get(lang) or {}).get("bench") or {}
        entry = bench_data.get(key)
        if not entry:
            return EMDASH
        return _format_speedup(entry.get("speedup"))

    return f"{lookup('c', row.c_bench_key)} / {lookup('fortran', row.f_bench_key)}"


def report_key_coverage(system: dict, rows: list[ops_mod.UnifiedRow]) -> None:
    """Warn (stderr, non-fatal) about ops.py <-> JSON key mismatches.

    A key expected by ops.py but absent from the JSON silently renders as
    EMDASH (the ``.get()`` fallbacks above), and a JSON key no row references
    is silently dropped — both usually mean a renamed bench label or a stale
    ops.py entry. An entirely absent dict (e.g. a ``--skip-fuzz`` run has no
    fuzz stats, a single-language JSON has no bench dict for the other side)
    is normal and not reported key-by-key. ``fuzz=None`` / ``"exact"`` ops
    and the one-language ``None`` bench keys expect nothing. Unreferenced
    fuzz keys are reported as a count only: fuzz.cc deliberately covers far
    more surface (ctors, assignments, …) than the table shows.
    """
    def fuzz_keys(spec: ops_mod.FuzzSpec) -> list[str]:
        if spec is None or spec == "exact":
            return []
        if isinstance(spec, str):
            return [spec]
        return [key for _label, key in spec]

    checks = [
        ("c bench",
         {r.c_bench_key for r in rows if r.c_bench_key is not None},
         (system.get("c") or {}).get("bench") or {}, True),
        ("fortran bench",
         {r.f_bench_key for r in rows if r.f_bench_key is not None},
         (system.get("fortran") or {}).get("bench") or {}, True),
        ("fuzz",
         {key for r in rows for key in fuzz_keys(r.fuzz)},
         (system.get("c") or {}).get("fuzz") or {}, False),
    ]
    warnings: list[str] = []
    for label, expected, data, list_unclaimed in checks:
        if not data:
            continue  # whole section absent from this JSON — legitimate
        missing = sorted(expected - data.keys())
        unclaimed = sorted(data.keys() - expected)
        if missing:
            warnings.append(f"{label}: {len(missing)} ops.py key(s) missing from"
                            f" the JSON (rendered as {EMDASH}): {', '.join(missing)}")
        if unclaimed and list_unclaimed:
            warnings.append(f"{label}: {len(unclaimed)} JSON key(s) not referenced"
                            f" by ops.py: {', '.join(unclaimed)}")
        elif unclaimed:
            warnings.append(f"{label}: {len(unclaimed)} JSON key(s) not shown in"
                            f" the table (expected: fuzz coverage is a superset)")
    if warnings:
        print(f"build_benchmark_md: {len(warnings)} key-coverage warning(s):",
              file=sys.stderr)
        for w in warnings:
            print(f"  - {w}", file=sys.stderr)


def load_system(path: Path) -> dict:
    blob = json.loads(path.read_text())
    return {
        **blob.get("system", {}),
        "fortran": blob.get("fortran", {}),
        "c": blob.get("c") or blob.get("cpp", {}),
    }


def render(system: dict, template_path: Path) -> str:
    env = Environment(
        loader=FileSystemLoader(str(template_path.parent)),
        undefined=StrictUndefined,
        trim_blocks=False,
        lstrip_blocks=False,
        keep_trailing_newline=True,
    )
    env.filters["err_dd"] = err_dd_filter
    env.filters["err_qp"] = err_qp_filter
    env.filters["speedup_pair"] = speedup_pair_filter
    tmpl = env.get_template(template_path.name)
    return tmpl.render(system=system, rows=ops_mod.unified_rows())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path,
                    help="per-system JSON file (e.g. doc/dev/benchmark/baseline/benchmark-m1-max.json)")
    ap.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT,
                    help=f"output path (default: {DEFAULT_OUTPUT.relative_to(REPO_ROOT)})")
    ap.add_argument("--template", type=Path,
                    default=TEMPLATE_DIR / "BENCHMARK.md.j2",
                    help="Jinja2 template (default: benchmark/BENCHMARK.md.j2)")
    ap.add_argument("--stdout", action="store_true", help="print to stdout instead of --output")
    args = ap.parse_args()

    if not args.input.is_file():
        print(f"error: {args.input} not found", file=sys.stderr)
        return 1

    system = load_system(args.input)
    out = render(system, args.template)

    if args.stdout:
        sys.stdout.write(out)
    else:
        args.output.write_text(out)
        print(f"wrote {args.output}", file=sys.stderr)

    report_key_coverage(system, ops_mod.unified_rows())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
