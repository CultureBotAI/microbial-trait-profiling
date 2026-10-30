#!/usr/bin/env bash
# Purpose: refuse an output root under the system temporary directory, then run a command.
# Inputs: the output root and the command; outputs: whatever the command writes.
# Issue #1087; 2026-10-25. The recipes for modules that the arms also run use this, so
# those modules and the arms' recorded plans stay unchanged. --allow-temporary-output
# anywhere among the command's arguments admits a throwaway run and is not passed on.
# Usage: bash scripts/guard_output_root.sh OUTPUT COMMAND [ARG...]

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 OUTPUT COMMAND [ARG...]" >&2
    exit 2
fi

output="$1"
shift
allow=""
for argument; do
    shift
    if [[ "${argument}" == --allow-temporary-output ]]; then
        allow="${argument}"
    else
        set -- "$@" "${argument}"
    fi
done
uv run python -m mtp.workflow.output_roots "${output}" ${allow}
exec "$@"
