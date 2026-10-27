#!/usr/bin/env bash
# Purpose: qualify #1052 against one completed phenotype parent using merged code.
# Inputs: clean producer worktree, comparison manifest, label-support bundle,
#         fresh durable output root, expected producer commit.
# Outputs: parent table, scratch, command log and exit status; no replay.
# Issue: #1052. Date: 2026-10-25. Run through the coordination heavy_run wrapper.
set -euo pipefail
if [[ "$#" != 5 ]]; then
    printf '%s\n' 'usage: generate_delta_z_parent_1052.sh WORKTREE COMPARISON LABEL_SUPPORT OUTPUT COMMIT' >&2
    exit 2
fi
producer_worktree="$1"
comparison_manifest="$2"
label_support="$3"
output_root="$4"
expected_commit="$5"
cd "$producer_worktree"
test "$(git rev-parse HEAD)" = "$expected_commit"
git diff --quiet
git diff --cached --quiet
mkdir -p "$(dirname "$output_root")"
mkdir "$output_root"
set +e
just three-test-table --comparison-run "$comparison_manifest" \
    --label-support "$label_support" --output "$output_root/table" \
    --scratch-dir "$output_root/scratch" --code-commit "$expected_commit" \
    > "$output_root/generation.log" 2>&1
result=$?
set -e
printf '%s\n' "$result" > "$output_root/generation.exit"
exit "$result"
