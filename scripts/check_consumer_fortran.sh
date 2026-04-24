#!/usr/bin/env bash
# Drive the test/consumer-fortran smoke test against a staged install of
# the in-tree build. Proves that the source-distributed multifloatsf
# package (shipped as pre-expanded .f90 + multifloatsfConfig.cmake) is
# find_package()-able by a fresh CMake project and that the consumer can
# exercise both the elemental Fortran path and the bind(c) C-ABI backend.
#
# Usage: scripts/check_consumer_fortran.sh [<main-build-dir>]
#   main-build-dir defaults to `build` relative to the repo root.
#
# Exit codes:
#   0 — success.
#   non-zero — anything that failed (install, consumer configure, build,
#              run, or output check).
set -euo pipefail

REPO="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
MAIN_BUILD="${1:-$REPO/build}"
STAGE="$(mktemp -d -t multifloatsf-consumer-install.XXXXXX)"
CONSUMER_BUILD="$(mktemp -d -t multifloatsf-consumer-build.XXXXXX)"
cleanup() { rm -rf "$STAGE" "$CONSUMER_BUILD"; }
trap cleanup EXIT

echo "==> Installing from $MAIN_BUILD into $STAGE"
cmake --install "$MAIN_BUILD" --prefix "$STAGE"

echo "==> Configuring consumer project"
# Reuse the same compilers the main build selected. The installed C
# library is tagged by compiler (LTO-tagged archive) so the consumer must
# use a matching toolchain or we'd need a non-LTO portable variant
# (-DMULTIFLOATS_USE_LTO=OFF at main build time). CMAKE_CXX_COMPILER is
# not exposed in CMakeCache.txt — pull it from the per-version
# CMakeFiles/<ver>/CMake<Lang>Compiler.cmake instead.
get_compiler() {
    local lang="$1"
    local files
    files=( "$MAIN_BUILD"/CMakeFiles/*/CMake"${lang}"Compiler.cmake )
    [[ -f "${files[0]}" ]] || { echo "!! cannot locate $lang compiler info"; exit 1; }
    # Lines look like: set(CMAKE_CXX_COMPILER "/path/to/compiler")
    sed -n "s|^set(CMAKE_${lang}_COMPILER \"\\(.*\\)\")$|\\1|p" "${files[0]}" | head -n1
}
MAIN_CXX="$(get_compiler CXX)"
MAIN_FC="$(get_compiler Fortran)"
echo "    using CXX=$MAIN_CXX"
echo "    using FC =$MAIN_FC"

cmake -S "$REPO/test/consumer-fortran" -B "$CONSUMER_BUILD" \
    -DCMAKE_PREFIX_PATH="$STAGE" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER="$MAIN_CXX" \
    -DCMAKE_Fortran_COMPILER="$MAIN_FC" \
    -Wno-dev

echo "==> Building consumer"
cmake --build "$CONSUMER_BUILD" -j

echo "==> Running consumer"
OUT="$("$CONSUMER_BUILD/consumer")"
echo "$OUT"

# The "ok" marker proves the program reached the final print. Exact bit
# patterns of the DD limbs depend on the compiler's IEEE conformance —
# don't lock them down here; the in-tree precision_cpp / precision_fortran
# ctests already pin those.
if ! grep -q '^ok$' <<<"$OUT"; then
    echo "!! consumer missing 'ok' marker — failure"
    exit 1
fi

# Sanity-check the hi-limb values. Match the first 6–8 significant digits
# only; the full ULP pattern depends on the compiler's print formatter.
if ! grep -Eq 'z *= *[[:space:]]*3\.0{5,}E\+00' <<<"$OUT"; then
    echo "!! z hi limb did not print 3.0 as expected"
    exit 1
fi
if ! grep -Eq 'sqrt2 *= *[[:space:]]*1\.41421356[0-9]*E\+00' <<<"$OUT"; then
    echo "!! sqrt(2) hi limb did not print the libm result"
    exit 1
fi

echo "==> consumer smoke PASSED"
