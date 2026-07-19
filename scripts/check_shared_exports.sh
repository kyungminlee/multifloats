#!/bin/bash
# Verify the shared library's DYNAMIC symbol table exports ONLY the declared
# extern "C" C ABI. This is the actual binary/DSO ABI a consumer links or
# dlopen()s, and it must be exactly the stable surface encoded by SOVERSION ==
# MULTIFLOATS_ABI_VERSION — nothing more.
#
# The shared object is built with -fvisibility=hidden + -fvisibility-inlines-
# hidden, so only the MULTIFLOATS_API-marked extern "C" functions have DEFAULT
# visibility and reach the dynsym; the multifloats:: C++ API is GLOBAL HIDDEN
# and the header-inline weak COMDAT copies are WEAK HIDDEN — neither is
# exported. A leak here means a visibility regression widened the DSO ABI: a
# C++ symbol (mangled _ZN11multifloats…) or a TU-internal helper escaped into
# the dynamic table. That is stricter than the static-archive audit
# (check_exported_symbols.sh), which also allows the multifloats:: C++ API
# because a static consumer links those out of the archive directly.
#
# ELF only (nm -D). The release builds the shared object in its Linux and macOS
# non-LTO jobs, but this audit runs on ELF; the Mach-O .dylib's export list is
# enforced by ld64's default treatment of hidden visibility.
#
# Usage:
#   check_shared_exports.sh <path/to/libmultifloats.so> <path/to/header>
set -e
. "$(dirname "$0")/lib/c_abi_symbols.sh"
LIB="$1"
HEADER="$2"
if [ -z "$LIB" ] || [ -z "$HEADER" ]; then
    echo "usage: $0 <libmultifloats.so> <header>" >&2
    exit 2
fi
if [ ! -f "$LIB" ] || [ ! -f "$HEADER" ]; then
    echo "check_shared_exports.sh: missing input file" >&2
    echo "  lib=$LIB" >&2
    echo "  header=$HEADER" >&2
    exit 2
fi

CABI=$(mktemp)
MANGLED=$(mktemp)
DEMANGLED=$(mktemp)
trap 'rm -f "$CABI" "$MANGLED" "$DEMANGLED"' EXIT

# extern "C" C-ABI keep-list — every `MULTIFLOATS_API <type> <name>(` in the
# header; the canonical scan shared with check_exported_symbols.sh /
# check_fortran_abi_sync.sh (lib/c_abi_symbols.sh).
list_c_abi_symbols "$HEADER" > "$CABI"

if [ ! -s "$CABI" ]; then
    echo "check_shared_exports.sh: no MULTIFLOATS_API symbols found in $HEADER" >&2
    exit 1
fi

# Exported, defined dynamic symbols. `nm -D --defined-only` lists the dynamic
# table; keep the exported binding types (T/t weak-or-strong text, plus any
# exported data D/B/R/V/W) and drop the address column.
nm -D --defined-only "$LIB" 2>/dev/null \
    | awk '$2 ~ /^[TWDBRV]$/ {print $3}' | sort -u > "$MANGLED"
c++filt < "$MANGLED" > "$DEMANGLED"

LEAKS=$(paste "$MANGLED" "$DEMANGLED" | awk -F'\t' -v cabi="$CABI" '
  BEGIN {
    while ((getline n < cabi) > 0) ok[n] = 1
    # Linker-synthesized markers the ELF linker may place in .dynsym. Not part
    # of the API; ignore if present (recent GCC hides most of these).
    split("_init _fini _end _edata __bss_start __end__ __data_start __dso_handle", g, " ")
    for (i in g) ok[g[i]] = 1
  }
  $1 in ok { next }   # declared C ABI or a known linker marker (by nm name)
  { print $2 }        # demangled form makes a C++ leak legible
')

if [ -n "$LEAKS" ]; then
    echo "check_shared_exports: undeclared dynamic export(s) in $LIB:" >&2
    echo "$LEAKS" | sed 's/^/  /' >&2
    echo "The shared library must export ONLY the extern \"C\" C ABI. A C++ or" >&2
    echo "internal symbol reaching the dynsym is a visibility regression —" >&2
    echo "check CXX_VISIBILITY_PRESET/VISIBILITY_INLINES_HIDDEN and MULTIFLOATS_API." >&2
    exit 1
fi

echo "check_shared_exports: OK ($(wc -l < "$MANGLED") dynamic exports, all declared C ABI)"
