#!/usr/bin/env python3
"""Run the fortran_fuzz / fortran_bench / cpp_fuzz / cpp_bench executables
and produce a per-system JSON blob suitable for ``bench/build_benchmark_md.py``.

Usage (simplest form):

    python3 bench/run_benchmarks.py --name "M1 Max"

The output lands at ``bench/results/benchmark-<slug>.json`` where ``<slug>``
is the ``--name`` slugified (lowercased, spaces → hyphens). All system
description fields auto-detect, but every one can be overridden with the
matching flag.

If you already have a build tree, pass ``--build-dir <path> --no-build``
to skip the cmake configure + build step.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_BUILD_DIR = REPO_ROOT / "build"
RESULTS_DIR = REPO_ROOT / "bench" / "results"


# ---------------------------------------------------------------------------
# System-info probing
# ---------------------------------------------------------------------------

def _run(cmd: list[str]) -> str:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return (out.stdout or out.stderr).strip()
    except (FileNotFoundError, OSError):
        return ""


def detect_cpu() -> str:
    if sys.platform == "darwin":
        brand = _run(["sysctl", "-n", "machdep.cpu.brand_string"])
        cores = _run(["sysctl", "-n", "hw.ncpu"])
        arch = platform.machine()
        core_part = f", {cores} cores" if cores else ""
        arch_part = f" ({arch})" if arch else ""
        return f"{brand}{arch_part}{core_part}".strip()
    if sys.platform.startswith("linux"):
        model = ""
        try:
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("model name"):
                        model = line.split(":", 1)[1].strip()
                        break
        except OSError:
            pass
        cores = _run(["nproc"])
        mem = ""
        try:
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemTotal"):
                        kb = int(line.split()[1])
                        mem = f", {kb // (1024 * 1024)} GB"
                        break
        except (OSError, ValueError):
            pass
        arch = platform.machine()
        core_part = f", {cores} cores" if cores else ""
        arch_part = f" ({arch})" if arch else ""
        return f"{model}{arch_part}{core_part}{mem}".strip().strip(",")
    return platform.processor() or platform.machine()


def detect_os() -> str:
    if sys.platform == "darwin":
        ver = _run(["sw_vers", "-productVersion"])
        name = _run(["sw_vers", "-productName"]) or "macOS"
        kernel = f"Darwin {platform.release()}" if platform.release() else ""
        base = f"{name} {ver}".strip()
        return f"{base} ({kernel})" if kernel else base
    if sys.platform.startswith("linux"):
        os_release: dict[str, str] = {}
        try:
            with open("/etc/os-release") as f:
                for line in f:
                    if "=" in line:
                        k, v = line.rstrip().split("=", 1)
                        os_release[k] = v.strip('"')
        except OSError:
            pass
        pretty = os_release.get("PRETTY_NAME", "")
        kernel = f"Linux {platform.release()}" if platform.release() else ""
        if pretty and kernel:
            return f"{pretty} ({kernel})"
        return pretty or kernel or "Linux"
    return f"{platform.system()} {platform.release()}".strip()


def detect_compiler(build_dir: Path) -> str:
    """Best-effort. Reads ``CMakeCache.txt`` if present, otherwise probes
    ``gfortran``/``g++`` on PATH."""
    fc, cxx = "", ""
    cache = build_dir / "CMakeCache.txt"
    if cache.exists():
        for line in cache.read_text().splitlines():
            if line.startswith("CMAKE_Fortran_COMPILER:"):
                fc = line.split("=", 1)[1]
            elif line.startswith("CMAKE_CXX_COMPILER:"):
                cxx = line.split("=", 1)[1]
    if not fc:
        fc = shutil.which("gfortran") or shutil.which("gfortran-15") or ""
    if not cxx:
        cxx = shutil.which("g++") or shutil.which("g++-15") or ""
    fc_ver = _version_line(fc) if fc else ""
    cxx_ver = _version_line(cxx) if cxx else ""
    if fc_ver and cxx_ver and fc_ver == cxx_ver:
        return f"GNU Fortran / g++ {fc_ver}"
    if fc_ver and cxx_ver:
        return f"GNU Fortran {fc_ver} / g++ {cxx_ver}"
    return fc_ver or cxx_ver or "unknown"


_VERSION_RE = re.compile(r"(\d+(?:\.\d+){1,2})\s*$")


def _version_line(exe: str) -> str:
    """Extract the trailing ``MAJOR.MINOR[.PATCH]`` from ``exe --version`` — and
    any parenthesized distribution tag on the same line — so callers get e.g.
    ``"13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1)"``."""
    out = _run([exe, "--version"])
    if not out:
        return ""
    first = out.splitlines()[0]
    m = _VERSION_RE.search(first)
    ver = m.group(1) if m else ""
    dist = ""
    lp, rp = first.find("("), first.rfind(")")
    if lp != -1 and rp > lp:
        dist = first[lp : rp + 1]
    if ver and dist:
        return f"{ver} {dist}"
    return ver or first.strip()


def detect_build(build_dir: Path) -> str:
    cmake_ver = _run(["cmake", "--version"]).splitlines()
    cmake = cmake_ver[0].replace("cmake version ", "CMake ") if cmake_ver else "CMake ?"

    # Extract optimization / LTO flags from the actual compile line so the
    # label reflects what was built, not a hardcoded guess.
    flags_file = build_dir / "src" / "CMakeFiles" / "multifloats.dir" / "flags.make"
    opt_tokens: list[str] = []
    try:
        for line in flags_file.read_text().splitlines():
            if line.startswith("CXX_FLAGS"):
                seen: set[str] = set()
                for tok in line.split("=", 1)[1].split():
                    if re.fullmatch(r"-O\d", tok) or tok.startswith("-flto"):
                        if tok not in seen:
                            seen.add(tok)
                            opt_tokens.append(tok)
                break
    except OSError:
        pass
    flags = " ".join(opt_tokens) if opt_tokens else "-O3"

    # STATIC vs OBJECT vs SHARED — read from CMake's target properties via
    # the generated archive/object artifacts.
    lib_type = "STATIC library"
    if (build_dir / "src" / "libmultifloats.a").exists():
        lib_type = "STATIC library"
    elif list(build_dir.glob("src/CMakeFiles/multifloats.dir/*.o")):
        lib_type = "OBJECT library"

    return f"{cmake}, `{flags}`, {lib_type}"


def slugify(name: str) -> str:
    s = name.strip().lower()
    s = re.sub(r"[^a-z0-9]+", "-", s).strip("-")
    return s


# ---------------------------------------------------------------------------
# Build + run
# ---------------------------------------------------------------------------

def cmake_build(build_dir: Path, targets: list[str]) -> None:
    build_dir.mkdir(parents=True, exist_ok=True)
    if not (build_dir / "CMakeCache.txt").exists():
        subprocess.run(
            ["cmake", "-S", str(REPO_ROOT), "-B", str(build_dir), "-DCMAKE_BUILD_TYPE=Release"],
            check=True,
        )
    subprocess.run(
        ["cmake", "--build", str(build_dir), "-j", "--target", *targets],
        check=True,
    )


def run_exe(exe: Path) -> str:
    print(f"  running {exe.name} ...", file=sys.stderr, flush=True)
    out = subprocess.run([str(exe)], capture_output=True, text=True, check=True)
    return out.stdout


# ---------------------------------------------------------------------------
# Stdout parsing
# ---------------------------------------------------------------------------

# A bench row looks like:
#   " op_name        n_ops  tq.xxxx  tf.xxxx  speed.xxx"
# with trailing ``x`` on the speedup. The op name can contain spaces /
# punctuation (e.g. "add (dd+dp)", "arr_matmul (8x8*8)", "ldexp(.,5)"), so we
# anchor on the four numeric trailing fields.
_BENCH_ROW_RE = re.compile(
    r"^\s+(?P<op>\S.+?)\s+"
    r"(?P<n_ops>\d+)\s+"
    r"(?P<tq>-?\d+\.\d+)\s+"
    r"(?P<tf>-?\d+\.\d+)\s+"
    r"(?P<speedup>-?\d+\.\d+)x\s*$"
)


def parse_bench(stdout: str) -> dict[str, dict]:
    results: dict[str, dict] = {}
    for line in stdout.splitlines():
        m = _BENCH_ROW_RE.match(line)
        if not m:
            continue
        op = m.group("op").strip()
        if op.lower() == "op":  # header line
            continue
        results[op] = {
            "n_ops": int(m.group("n_ops")),
            "qp_time": float(m.group("tq")),
            "dd_time": float(m.group("tf")),
            "speedup": float(m.group("speedup")),
        }
    return results


# A fuzz row looks like:
#    op_name       count   max_rel       mean_rel
# with the count sometimes followed by "   (no data)" when count=0.
_FUZZ_ROW_RE = re.compile(
    r"^\s+(?P<op>\S+)\s+"
    r"(?P<n>\d+)\s+"
    r"(?P<max_rel>-?\d+\.\d+[Ee][+-]?\d+)\s+"
    r"(?P<mean_rel>-?\d+\.\d+[Ee][+-]?\d+)\s*$"
)


def parse_fuzz(stdout: str) -> dict[str, dict]:
    results: dict[str, dict] = {}
    in_report = False
    for line in stdout.splitlines():
        if "Per-operation precision report" in line:
            in_report = True
            continue
        if not in_report:
            continue
        if not line.strip():
            # trailing blank line terminates the report
            if results:
                break
            continue
        m = _FUZZ_ROW_RE.match(line)
        if not m:
            continue
        op = m.group("op").strip()
        if op.lower() == "op":
            continue
        results[op] = {
            "n": int(m.group("n")),
            "max_rel": float(m.group("max_rel")),
            "mean_rel": float(m.group("mean_rel")),
        }
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True, help='short display name, e.g. "M1 Max"')
    ap.add_argument("--slug", default=None, help="filename slug override (default: slugified --name)")
    ap.add_argument("--cpu", default=None)
    ap.add_argument("--os", dest="os_", default=None)
    ap.add_argument("--compiler", default=None)
    ap.add_argument("--build", default=None, help="e.g. 'CMake 4.3.1, -O3 -flto, STATIC library'")
    ap.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIR)
    ap.add_argument("--no-build", action="store_true",
                    help="skip cmake configure + build (binaries must already exist)")
    ap.add_argument("--skip-fortran", action="store_true")
    ap.add_argument("--skip-cpp", action="store_true")
    ap.add_argument("--skip-fuzz", action="store_true",
                    help="skip the 1M-iteration fuzz runs (speedups only)")
    ap.add_argument("--output-dir", type=Path, default=RESULTS_DIR)
    args = ap.parse_args()

    targets: list[str] = []
    if not args.skip_fortran:
        targets += ["fortran_bench"]
        if not args.skip_fuzz:
            targets += ["fortran_fuzz"]
    if not args.skip_cpp:
        targets += ["cpp_bench"]
        if not args.skip_fuzz:
            targets += ["cpp_fuzz"]

    if not args.no_build:
        print("Building benchmark executables ...", file=sys.stderr)
        cmake_build(args.build_dir, targets)

    data: dict = {
        "system": {
            "short_name": args.name,
            "cpu": args.cpu or detect_cpu(),
            "os": args.os_ or detect_os(),
            "compiler": args.compiler or detect_compiler(args.build_dir),
            "build": args.build or detect_build(args.build_dir),
        },
        "fortran": {},
        "cpp": {},
    }

    if not args.skip_fortran:
        data["fortran"]["bench"] = parse_bench(run_exe(args.build_dir / "fortran_bench"))
        if not args.skip_fuzz:
            data["fortran"]["fuzz"] = parse_fuzz(run_exe(args.build_dir / "fortran_fuzz"))
    if not args.skip_cpp:
        data["cpp"]["bench"] = parse_bench(run_exe(args.build_dir / "cpp_bench"))
        if not args.skip_fuzz:
            data["cpp"]["fuzz"] = parse_fuzz(run_exe(args.build_dir / "cpp_fuzz"))

    slug = args.slug or slugify(args.name)
    out_path = args.output_dir / f"benchmark-{slug}.json"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    print(f"wrote {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
