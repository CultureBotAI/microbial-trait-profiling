#!/usr/bin/env bash
# Purpose: run the real coverage CLI once against hash-bound pinned taxonomy (#1008).
# Inputs: core-run directory, KG edge path, expected edge SHA-256 and byte size.
# Outputs: this checkout's fresh .tmp/1008-pinned-validation only; 2026-10-23.
# Invoke through the shared heavy_run wrapper. This is validation, not published data.
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd -P)"
core_run="${1:?Pass the completed full-phenotype core-run directory}"
taxonomy_edges="${2:?Pass the pinned KG edge file}"
expected_sha="${3:?Pass its manifest SHA-256}"
expected_size="${4:?Pass its manifest byte size}"
if [ "$#" -ne 4 ]; then
    echo 'Expected exactly four arguments' >&2
    exit 2
fi
taxa="$core_run/steps/export_tables/derived/bugsigdb_taxa.tsv"
matrix="$core_run/steps/taxa_matrix/derived/taxa_trait_matrix.tsv"
profiles="$core_run/steps/count_profiles/derived/study_trait_profiles_wide.tsv"
experiments="$core_run/steps/export_tables/derived/experiments.tsv"
kg_provenance="$repo_root/data/kg-microbe/provenance.json"
inputs=("$taxonomy_edges" "$taxa" "$matrix" "$profiles" "$experiments"
    "$core_run/run.json" "$kg_provenance")
for input in "${inputs[@]}"; do
    test -f "$input"
done
# No fallback raw file is present in this isolated checkout. Avoid silently
# introducing a new input into this precisely scoped one-time validation.
if [ -e "$repo_root/data/raw/full_dump.csv" ]; then
    echo 'Unexpected fallback raw CSV; review the validation inputs first' >&2
    exit 2
fi
audit_root="$repo_root/.tmp/1008-pinned-validation"
mkdir -p "$repo_root/.tmp"
mkdir "$audit_root"  # Refuse an existing validation directory; never overwrite it.
cd "$repo_root"
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="$audit_root/matplotlib"
shasum -a 256 "${inputs[@]}" > "$audit_root/inputs-before.sha256"
read -r actual_sha _ < "$audit_root/inputs-before.sha256"
actual_size="$(wc -c < "$taxonomy_edges" | tr -d '[:space:]')"
printf '%s\n' "$actual_size" > "$audit_root/taxonomy-size-before.txt"
if [ "$actual_sha" != "$expected_sha" ] || [ "$actual_size" != "$expected_size" ]; then
    echo 'Pinned taxonomy identity mismatch; producer was not run' >&2
    exit 2
fi
shasum -a 256 scripts/create_taxa_coverage_histogram.py \
    scripts/adhoc/issue1008_validate_pinned_coverage.sh > "$audit_root/code.sha256"
set +e
/usr/bin/time -l just taxa-coverage \
    --strain-rollup phenotype --tables-only \
    --taxa "$taxa" --matrix "$matrix" --profiles "$profiles" \
    --experiments "$experiments" --taxonomy-edges "$taxonomy_edges" \
    --kg-provenance "$kg_provenance" --output-dir "$audit_root/results" \
    > "$audit_root/producer.log" 2>&1
producer_exit=$?
set -e
printf '%s\n' "$producer_exit" > "$audit_root/producer.exit"
shasum -a 256 "${inputs[@]}" > "$audit_root/inputs-after.sha256"
wc -c < "$taxonomy_edges" | tr -d '[:space:]' > "$audit_root/taxonomy-size-after.txt"
printf '\n' >> "$audit_root/taxonomy-size-after.txt"
cmp "$audit_root/inputs-before.sha256" "$audit_root/inputs-after.sha256"
cmp "$audit_root/taxonomy-size-before.txt" "$audit_root/taxonomy-size-after.txt"
shasum -a 256 -c "$audit_root/code.sha256"
exit "$producer_exit"
