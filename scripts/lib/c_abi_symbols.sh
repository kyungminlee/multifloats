#!/bin/bash
# Canonical extern "C" C-ABI symbol scan, shared by the three audit scripts
# (check_fortran_abi_sync.sh / check_exported_symbols.sh /
# check_shared_exports.sh). This file is sourced, not executed — source it
# via "$(dirname "$0")/lib/c_abi_symbols.sh" so it resolves regardless of
# the caller's working directory (ctest and CI invoke the scripts from
# different CWDs).
#
# The pipeline used to be pasted into each script; the copies drifted (one
# kept the leading `*` of `char *to_charsdd(`, poisoning its keep-list).
# Keeping it here means every script audits the exact same symbol set.
#
# Usage:
#   list_c_abi_symbols <path/to/multifloats.h>
#
# Emits one bare symbol name per line, sorted and deduplicated: the name of
# every `MULTIFLOATS_API <type> <name>(` declaration in the header. The sed
# strips the `(` left by the grep and a leading `*`/`&` that hugs a
# pointer-returning name (e.g. `char *to_charsdd(` -> `to_charsdd`).
list_c_abi_symbols() {
    grep -oE "^MULTIFLOATS_API[^(]+\(" "$1" \
        | awk '{print $NF}' | sed -e 's/($//' -e 's/^[*&]*//' | sort -u
}
