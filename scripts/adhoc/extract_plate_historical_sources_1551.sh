#!/usr/bin/env bash
# Preserve small committed source files for historical plate-number tests.
# Inputs: named Git blobs below. Outputs: exclusive test fixtures, no analysis rerun.
# Issues #1551/#1555; 2026-11-05. Run from this repository's worktree root.
set -euo pipefail
set -o noclobber

plate_source_dir=tests/fixtures/cellulose_tufte_historical_sources
mkdir -p "$plate_source_dir"

extract_source() {
    local revision_path=$1
    local output_name=$2
    local source_size
    source_size=$(git cat-file -s "$revision_path")
    if (( source_size > 262144 )); then
        printf 'Refusing historical source larger than 256 KiB: %s\n' "$revision_path" >&2
        return 1
    fi
    git show "$revision_path" > "$plate_source_dir/$output_name"
    test -s "$plate_source_dir/$output_name"
    test "$(git hash-object "$plate_source_dir/$output_name")" = \
        "$(git rev-parse "$revision_path")"
}

extract_source \
    7970cba2ddd40a636470ce767e22ed0b805bee33:data/analysis/cellulose_decr/cellulose_decr_summary.json \
    study_modal_7970cba.json
extract_source \
    a3c32cf9e32e46ff55b2077ed596b498e2c3e1cb:data/derived/integrated_evidence/integrated_evidence.md \
    integrated_evidence_a3c32cf.md
extract_source \
    a3c32cf9e32e46ff55b2077ed596b498e2c3e1cb:methods/integrated_evidence.md \
    methods_integrated_evidence_a3c32cf.md
