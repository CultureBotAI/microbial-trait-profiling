#!/usr/bin/env bash
# Preserve the approved path-redacted ORA audit for historical Plate V source pins.
# Inputs: Git audit documents at fa5d2343 and 781e92e4. Output: exclusive test fixture.
# Issue #1993; 2026-11-10. No analysis execution or Git history mutation.
set -euo pipefail
set -o noclobber

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

git_read() {
    GIT_CONFIG_GLOBAL=/dev/null GIT_CONFIG_NOSYSTEM=1 git "$@"
}

source_path=docs/reviews/2026-10-19-ora-provenance-refresh.md
source_revision=781e92e4
earlier_revision=fa5d2343
output=tests/fixtures/cellulose_tufte_historical_sources/ora_provenance_refresh_781e92e4.md
source_size=$(git_read cat-file -s "$source_revision:$source_path")
if (( source_size > 262144 )); then
    printf 'Refusing historical source larger than 256 KiB\n' >&2
    exit 1
fi

git_read show "$source_revision:$source_path" > "$output"
test -s "$output"
test "$(git_read hash-object "$output")" = "$(git_read rev-parse "$source_revision:$source_path")"

pooled_rows() {
    rg '^\| Pooled (enrichment|depletion) \|'
}

cmp \
    <(git_read show "$earlier_revision:$source_path" | pooled_rows) \
    <(pooled_rows < "$output")
printf 'Source commit: %s\n' "$(git_read rev-parse "$source_revision")"
printf 'Fixture blob: %s\n' "$(git_read hash-object "$output")"
printf 'Fixture bytes: %s\n' "$source_size"
printf 'Pooled rows unchanged from %s:\n' "$earlier_revision"
pooled_rows < "$output"
