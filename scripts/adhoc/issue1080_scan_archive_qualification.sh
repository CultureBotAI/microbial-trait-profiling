#!/usr/bin/env bash
# Purpose: real-recipe synthetic v2 archive/restore proof for #1080, 2026-10-28.
# Input: fresh absolute durable output root; output: logs, archives and parity.json.
# Invoke through the lead's heavy_run.sh. No real F data or historical scans are read.
set -euo pipefail
if [[ $# != 1 || "$1" != /* ]]; then
    echo 'Usage: issue1080_scan_archive_qualification.sh ABSOLUTE_FRESH_OUTPUT_ROOT' >&2
    exit 2
fi
qualification_root="$1"
repository_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$repository_root"
helper=scripts/adhoc/issue1080_scan_archive_qualification.py
uv run python "$helper" prepare "$qualification_root"
mkdir "$qualification_root/logs" "$qualification_root/archives"
for scope in phenotype off; do
    just core-run \
        --input "raw_dump=$qualification_root/original/inputs/full_dump.csv" \
        --input "kg_edges=$qualification_root/original/inputs/merged-kg_edges.tsv" \
        --input "kg_nodes=$qualification_root/original/inputs/merged-kg_nodes.tsv" \
        --parameter source_mode=fixture \
        --contrast-population disease_vs_healthy_tier1_2_subject_controls \
        --strain-rollup "$scope" \
        --output-root "$qualification_root/original/core-$scope" \
        > "$qualification_root/logs/core-$scope.log" 2>&1
    manifests=("$qualification_root/original/core-$scope"/run-*/run.json)
    [[ ${#manifests[@]} == 1 && -f "${manifests[0]}" ]]
    just cellulose-decr-selection-scan "${manifests[0]}" \
        "$qualification_root/original/scan-$scope" --permutations 200 --seed 3 \
        > "$qualification_root/logs/scan-$scope.log" 2>&1
    just evidence-archive "$qualification_root/archives/$scope" \
        --core-run "${manifests[0]}" \
        --consumer "$qualification_root/original/scan-$scope/core_consumption.json" \
        > "$qualification_root/logs/archive-$scope.log" 2>&1
done
just cellulose-decr-selection-scan-scopes \
    "$qualification_root/original/scan-phenotype" "$qualification_root/original/scan-off" \
    "$qualification_root/unmoved" > "$qualification_root/logs/unmoved.log" 2>&1
cp -R "$qualification_root/archives" "$qualification_root/restored"
mv "$qualification_root/archives" "$qualification_root/retired-archives"
mv "$qualification_root/original" "$qualification_root/retired-original"
for scope in phenotype off; do
    just evidence-resolve "$qualification_root/restored/$scope/archive.json" \
        "$qualification_root/original/scan-$scope/core_consumption.json" \
        > "$qualification_root/resolved-$scope.txt" \
        2> "$qualification_root/logs/resolve-$scope.log"
done
just cellulose-decr-selection-scan-scopes \
    "$qualification_root/original/scan-phenotype" "$qualification_root/original/scan-off" \
    "$qualification_root/relocated" \
    --archive "$qualification_root/original/scan-phenotype" \
        "$qualification_root/restored/phenotype/archive.json" \
    --archive "$qualification_root/original/scan-off" \
        "$qualification_root/restored/off/archive.json" \
    > "$qualification_root/logs/relocated.log" 2>&1
uv run python "$helper" compare "$qualification_root"
