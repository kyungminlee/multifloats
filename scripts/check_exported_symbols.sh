#!/bin/bash
# Verify the static archive exports ONLY the declared public API as strong-
# global (T) symbols. This is the non-destructive replacement for the old
# localize_symbols.sh post-build strip.
#
# The strip mutated the archive (objcopy --keep-global-symbols + strip
# --strip-unneeded) to localize stray globals. In practice it removed zero
# real leaks — every global in the archive is declared API — while its one
# concrete effect was to localize the weak COMDAT namespace copies
# (multifloats::hypot etc.), which GNU ld then reported as "relocation refers
# to discarded section" at the consumer's link. That forced STRIP_SYMBOLS=OFF
# on every non-LTO release job. This check keeps the "only declared API is
# global" guarantee without ever touching the archive: it asserts the surface
# and fails the build (naming the offender) if an undeclared global appears.
#
# Allowed strong-global (T) symbols:
#   (1) extern "C" C ABI — every `MULTIFLOATS_API <type> <name>(` in the header.
#   (2) multifloats:: / multifloats::detail:: free functions + operator<<
#       (the C++ API; detail:: primitives are exercised directly by tests).
#   (3) std::NAME<multifloats::float64x2>(...) specializations (std::exp, …).
# Weak (W) symbols are the header-inline COMDAT copies (operator/, hypot,
# the ctor, …) — ODR-merge artifacts, always harmless — and are NOT audited.
# The library defines no global data (D/R/B/G), so only T is checked.
#
# A leak is a TU-internal helper that escaped with external linkage. The fix
# is at the source: move it into an anonymous namespace / mark it `static` /
# put it in multifloats::detail::, or declare it MULTIFLOATS_API if it is
# genuinely public.
#
# Usage:
#   check_exported_symbols.sh <path/to/libmultifloats.a> <path/to/header>
set -e
LIB="$1"
HEADER="$2"
if [ -z "$LIB" ] || [ -z "$HEADER" ]; then
    echo "usage: $0 <archive.a> <header>" >&2
    exit 2
fi
if [ ! -f "$LIB" ] || [ ! -f "$HEADER" ]; then
    echo "check_exported_symbols.sh: missing input file" >&2
    echo "  archive=$LIB" >&2
    echo "  header=$HEADER" >&2
    exit 2
fi

CABI=$(mktemp)
MANGLED=$(mktemp)
DEMANGLED=$(mktemp)
trap 'rm -f "$CABI" "$MANGLED" "$DEMANGLED"' EXIT

# (1) extern "C" C-ABI keep-list. The trailing sed strips a `(` left by the
# grep and a leading `*`/`&` that hugs a pointer-returning name (e.g.
# `char *to_charsdd(`). Mirrors the scan in check_fortran_abi_sync.sh.
grep -oE "^MULTIFLOATS_API[^(]+\(" "$HEADER" \
    | awk '{print $NF}' | sed -e 's/($//' -e 's/^[*&]*//' | sort -u > "$CABI"

if [ ! -s "$CABI" ]; then
    echo "check_exported_symbols.sh: no MULTIFLOATS_API symbols found in $HEADER" >&2
    exit 1
fi

# Strong-global (T) defined symbols only. Per-object `file.o:` headers and
# blank lines in `nm` archive output have no $2 and are skipped by the filter.
nm --defined-only "$LIB" 2>/dev/null | awk '$2 == "T" {print $3}' > "$MANGLED"
c++filt < "$MANGLED" > "$DEMANGLED"

LEAKS=$(paste "$MANGLED" "$DEMANGLED" | awk -F'\t' -v cabi="$CABI" '
  BEGIN { while ((getline n < cabi) > 0) ok[n] = 1 }
  # (1) declared extern "C" ABI (nm name == demangled name for these).
  $2 in ok { next }
  # (2) C++ API: any function in namespace multifloats (incl. detail:: and
  #     operator<<). The demangled form starts with "multifloats::".
  $2 ~ /^multifloats::/ { next }
  # (3) std:: template specializations on our type (a return-type prefix may
  #     precede it in the demangled string — match the specialization anywhere).
  $2 ~ /std::[A-Za-z_][A-Za-z0-9_]*<multifloats::float64x2>/ { next }
  # Anything else is an undeclared global leak.
  { print $2 }
')

if [ -n "$LEAKS" ]; then
    echo "check_exported_symbols: undeclared global symbol(s) in $LIB:" >&2
    echo "$LEAKS" | sed 's/^/  /' >&2
    echo "Move each into an anonymous namespace / static / multifloats::detail::," >&2
    echo "or declare it MULTIFLOATS_API if it is genuinely public." >&2
    exit 1
fi

echo "check_exported_symbols: OK ($(wc -l < "$MANGLED") global symbols, all declared API)"
