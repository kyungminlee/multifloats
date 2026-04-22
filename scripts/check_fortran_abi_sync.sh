#!/bin/bash
# Verify that every `bind(c, name=<sym>)` binding in the generated Fortran
# module whose target ends in `dd` (or `dd_<suffix>`, e.g. matmuldd_mm)
# has a matching MULTIFLOATS_API declaration in include/multifloats.h.
# This catches a class of drift the audit called out under M2: if someone
# adds or renames an extern "C" entry point in the C header without
# updating the hand-mirrored Fortran interfaces, the binding would become
# a link-time UNDEF that surfaces only when Fortran code actually calls
# the renamed function. This check surfaces it as a ctest failure.
#
# The converse is intentionally NOT checked: Fortran reimplements basic
# arithmetic (adddd / subdd / muldd / divdd / fmadd, comparison ops,
# complex helpers) natively for speed, so the C header is a superset of
# the Fortran bindings.
#
# Usage:
#   check_fortran_abi_sync.sh <path/to/multifloats.h> <path/to/generated.f90>
set -e
HEADER="$1"
GENERATED_F90="$2"
if [ -z "$HEADER" ] || [ -z "$GENERATED_F90" ]; then
    echo "usage: $0 <header> <generated-f90>" >&2
    exit 2
fi
if [ ! -f "$HEADER" ] || [ ! -f "$GENERATED_F90" ]; then
    echo "check_fortran_abi_sync.sh: missing input file" >&2
    echo "  header=$HEADER" >&2
    echo "  fortran=$GENERATED_F90" >&2
    exit 2
fi

C_API=$(mktemp)
F_DD=$(mktemp)
trap 'rm -f "$C_API" "$F_DD"' EXIT

grep -oE "^MULTIFLOATS_API[^(]+\(" "$HEADER" \
    | awk '{print $NF}' | sed 's/($//' | sort -u > "$C_API"

grep -oE "bind\(c,\s*name\s*=\s*['\"][^'\"]+['\"]\)" "$GENERATED_F90" \
    | sed -E "s/.*name\s*=\s*['\"]([^'\"]+)['\"].*/\1/" \
    | grep -E 'dd($|_)' | sort -u > "$F_DD"

MISSING=$(comm -23 "$F_DD" "$C_API")
if [ -n "$MISSING" ]; then
    echo "check_fortran_abi_sync: Fortran bindings target symbols not declared in $HEADER:" >&2
    echo "$MISSING" >&2
    exit 1
fi
echo "check_fortran_abi_sync: OK ($(wc -l < "$F_DD") Fortran dd-bindings, all present in header)"
