#!/bin/bash
# Verify that the C23-style `*f64x2` alias surface is in sync: every
# `MULTIFLOATS_API <type> <name>f64x2(` declaration in the public header
# must have a hand-written forwarder definition in one of the given source
# files, and every forwarder must be declared in the header. The forwarders
# (src/float64x2/f64x2_aliases.inc, plus to_charsf64x2 in
# src/multifloats_io.cc) are maintained by hand in lockstep with the
# header — see the comment atop f64x2_aliases.inc — so a new `*f64x2`
# declaration without its forwarder would become a link-time UNDEF that
# surfaces only when client code actually calls the alias. This check
# surfaces the drift as a ctest failure instead, in both directions.
#
# Usage:
#   check_f64x2_aliases.sh <path/to/multifloats.h> <forwarder-src>...
set -e
. "$(dirname "$0")/lib/c_abi_symbols.sh"
HEADER="$1"
shift
if [ -z "$HEADER" ] || [ "$#" -eq 0 ]; then
    echo "usage: $0 <header> <forwarder-src>..." >&2
    exit 2
fi
for SRC in "$HEADER" "$@"; do
    if [ ! -f "$SRC" ]; then
        echo "check_f64x2_aliases.sh: missing input file" >&2
        echo "  file=$SRC" >&2
        exit 2
    fi
done

DECLARED=$(mktemp)
DEFINED=$(mktemp)
trap 'rm -f "$DECLARED" "$DEFINED"' EXIT

# `*f64x2` slice of the C-ABI keep-list (the canonical header scan shared
# with the other check_*.sh scripts via lib/c_abi_symbols.sh).
list_c_abi_symbols "$HEADER" | grep 'f64x2$' > "$DECLARED"

if [ ! -s "$DECLARED" ]; then
    echo "check_f64x2_aliases.sh: no *f64x2 MULTIFLOATS_API symbols found in $HEADER" >&2
    exit 1
fi

# Forwarder definitions: `<type> <name>f64x2(` at the start of a line. The
# same awk/sed tail as the header scan strips the `(` left by the grep and
# a leading `*`/`&` hugging a pointer-returning name (`char *to_charsf64x2(`).
grep -hoE '^[A-Za-z_][^(]*f64x2\(' "$@" \
    | awk '{print $NF}' | sed -e 's/($//' -e 's/^[*&]*//' | sort -u > "$DEFINED"

UNDEFINED=$(comm -23 "$DECLARED" "$DEFINED")
if [ -n "$UNDEFINED" ]; then
    echo "check_f64x2_aliases: declared in $HEADER but no forwarder defined in: $*" >&2
    echo "$UNDEFINED" >&2
    exit 1
fi

UNDECLARED=$(comm -13 "$DECLARED" "$DEFINED")
if [ -n "$UNDECLARED" ]; then
    echo "check_f64x2_aliases: forwarders defined but not declared in $HEADER:" >&2
    echo "$UNDECLARED" >&2
    exit 1
fi
echo "check_f64x2_aliases: OK ($(wc -l < "$DECLARED") *f64x2 aliases, header and forwarders in sync)"
