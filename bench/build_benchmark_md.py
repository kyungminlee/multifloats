#!/usr/bin/env python3
"""Render ``BENCHMARK.md`` from per-system JSON files produced by
``bench/run_benchmarks.py``.

    python3 bench/build_benchmark_md.py                     # uses bench/results/*.json
    python3 bench/build_benchmark_md.py bench/results/a.json bench/results/b.json
    python3 bench/build_benchmark_md.py -o /tmp/BENCHMARK.md

With no positional arguments, loads every ``benchmark-*.json`` under
``bench/results/``, sorted by filename. Columns in the output tables appear
in the same order.
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
RESULTS_DIR = REPO_ROOT / "bench" / "results"
DEFAULT_OUTPUT = REPO_ROOT / "doc" / "BENCHMARK.md"

EMDASH = "\u2014"  # —
TIMES = "\u00d7"   # ×


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
    """Format quadmath-vs-MPFR rel-err in qp ulps."""
    if not _finite(v) or v < 0.0:
        return EMDASH
    if v == 0.0:
        return "exact"
    u = v / _QP_ULP
    if u < 10:
        return f"{u:.1f} qp ulp"
    return f"{u:.0f} qp ulp"


def err_filter(op: ops_mod.Op, system: dict, lang: str) -> str:
    """Render the DD-vs-MPFR err cell for one op × one system."""
    spec = op.fuzz
    if spec is None:
        return EMDASH
    if spec == "exact":
        return "exact"
    fuzz_data = (system.get(lang) or {}).get("fuzz") or {}
    if isinstance(spec, str):
        entry = fuzz_data.get(spec)
        if not entry or not _finite(entry.get("max_rel")):
            return EMDASH
        return _format_max_rel(entry["max_rel"])
    parts: list[str] = []
    for label, key in spec:
        entry = fuzz_data.get(key)
        if entry and _finite(entry.get("max_rel")):
            parts.append(f"{_format_max_rel(entry['max_rel'])} ({label})")
    if not parts:
        return EMDASH
    return " / ".join(parts)


def qp_err_filter(op: ops_mod.Op, systems: list[dict]) -> str:
    """Render the quadmath-vs-MPFR reference err cell for one op.

    System-independent — qp precision is the same on every machine. Pulls
    `max_qp` from the first system that has an entry.
    """
    spec = op.fuzz
    if spec is None or spec == "exact":
        return EMDASH
    for sys in systems:
        fuzz_data = (sys.get("c") or sys.get("cpp") or {}).get("fuzz") or {}
        if isinstance(spec, str):
            entry = fuzz_data.get(spec)
            if entry and _finite(entry.get("max_qp")):
                return _format_qp_err(entry["max_qp"])
        else:
            parts: list[str] = []
            for label, key in spec:
                entry = fuzz_data.get(key)
                if entry and _finite(entry.get("max_qp")):
                    parts.append(f"{_format_qp_err(entry['max_qp'])} ({label})")
            if parts:
                return " / ".join(parts)
    return EMDASH


def _format_speedup_number(v: float) -> str:
    if v >= 10.0:
        return f"{int(round(v))}{TIMES}"
    if v >= 1.0:
        return f"{v:.1f}{TIMES}"
    return f"{v:.2f}{TIMES}"


def speedup_filter(op: ops_mod.Op, system: dict, lang: str) -> str:
    bench_data = (system.get(lang) or {}).get("bench") or {}
    entry = bench_data.get(op.bench_key)
    if not entry or not _finite(entry.get("speedup")):
        return EMDASH
    v = float(entry["speedup"])
    if v <= 0.0:
        return EMDASH
    s = _format_speedup_number(v)
    if v >= 2.0:
        return f"**{s}**"
    return s


def load_systems(paths: list[Path]) -> list[dict]:
    systems: list[dict] = []
    for p in paths:
        blob = json.loads(p.read_text())
        sysrec = {
            **blob.get("system", {}),
            "fortran": blob.get("fortran", {}),
            # Accept both the new "c" key and the legacy "cpp" key.
            "c": blob.get("c") or blob.get("cpp", {}),
        }
        # Ensure short_name is present for the template.
        if "short_name" not in sysrec:
            sysrec["short_name"] = p.stem.removeprefix("benchmark-")
        systems.append(sysrec)
    return systems


def render(systems: list[dict], template_path: Path) -> str:
    env = Environment(
        loader=FileSystemLoader(str(template_path.parent)),
        undefined=StrictUndefined,
        trim_blocks=False,
        lstrip_blocks=False,
        keep_trailing_newline=True,
    )
    env.filters["err"] = err_filter
    env.filters["speedup"] = speedup_filter
    # qp_err is a template-global helper (single-arg, system-independent).
    env.globals["qp_err"] = lambda op: qp_err_filter(op, systems)
    tmpl = env.get_template(template_path.name)
    return tmpl.render(
        systems=systems,
        fortran_sections=ops_mod.FORTRAN_SECTIONS,
        c_sections=ops_mod.C_SECTIONS,
        qp_ref_ops=ops_mod.qp_reference_ops(),
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="*", type=Path,
                    help="per-system JSON files (default: all bench/results/benchmark-*.json)")
    ap.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT,
                    help=f"output path (default: {DEFAULT_OUTPUT.relative_to(REPO_ROOT)})")
    ap.add_argument("--template", type=Path,
                    default=TEMPLATE_DIR / "BENCHMARK.md.j2",
                    help="Jinja2 template (default: bench/BENCHMARK.md.j2)")
    ap.add_argument("--stdout", action="store_true", help="print to stdout instead of --output")
    args = ap.parse_args()

    if args.inputs:
        paths = sorted(args.inputs)
    else:
        paths = sorted(RESULTS_DIR.glob("benchmark-*.json"))
    if not paths:
        print(f"error: no benchmark JSON files found (looked in {RESULTS_DIR})", file=sys.stderr)
        return 1

    systems = load_systems(paths)
    out = render(systems, args.template)

    if args.stdout:
        sys.stdout.write(out)
    else:
        args.output.write_text(out)
        print(f"wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
