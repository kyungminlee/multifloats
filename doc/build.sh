#!/usr/bin/env bash
# Build the multifloats documentation: Doxygen (XML) -> Sphinx/Breathe (HTML).
#
# One-time setup:
#     uv venv doc/.venv
#     uv pip install --python doc/.venv -r doc/requirements.txt
#     # plus the `doxygen` system package on PATH
#
# Usage:  ./doc/build.sh   (run from the repo root or from doc/)
set -euo pipefail

DOC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="$DOC_DIR/.venv"
SPHINX="$VENV/bin/sphinx-build"

if ! command -v doxygen >/dev/null 2>&1; then
    echo "error: doxygen not found on PATH (install the system package)." >&2
    exit 1
fi
if [[ ! -x "$SPHINX" ]]; then
    echo "error: $SPHINX not found. Run:" >&2
    echo "    uv venv doc/.venv && uv pip install --python doc/.venv -r doc/requirements.txt" >&2
    exit 1
fi

echo ">> doxygen (header -> XML)"
( cd "$DOC_DIR" && doxygen Doxyfile )

echo ">> sphinx-build (-> HTML)"
"$SPHINX" -b html "$DOC_DIR" "$DOC_DIR/_build/html" "$@"

echo ">> done: $DOC_DIR/_build/html/index.html"
