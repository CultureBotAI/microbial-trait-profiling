#!/usr/bin/env bash
# Purpose: render one slice-2 development figure from pinned F evidence (#1046).
# Inputs: evidence base, fresh durable destination, figure name; 2026-10-29.
# Outputs: verified PDF/SVG/PNG/TSV/provenance; F archives remain read-only.
set -euo pipefail
if [[ $# -ne 3 ]]; then
    echo 'usage: issue1046_slice2_f.sh EVIDENCE_BASE NEW_OUTPUT FIGURE' >&2
    exit 2
fi
evidence_base="$1"
figure_output="$2"
figure_name="$3"
case "$figure_name" in
    rank-sensitivity|candidate-triples|trait-overlap) ;;
    *) echo 'Choose rank-sensitivity, candidate-triples, or trait-overlap' >&2; exit 2 ;;
esac
exec just publication-figures-slice2 \
    --three-test-run "$evidence_base/mtp-data/codex-1052-development/parents/0c9e0ea/tier12-phenotype/table" \
    --three-test-sha256 0547a89ccbdaddfed9e4a8c0fff3986141f42ea416d76c92b2845f52ca4d386d \
    --comparison-run "$evidence_base/mtp-runs/F-511c891/runs/tier12-phenotype/comparison" \
    --derived-run "$evidence_base/mtp-data/codex-1160-main-merge/9c9701fe21fd321b8d05586aefc3a1b5a3957e7a/phenotype" \
    --derived-sha256 b198dcbaed04a04caff2f0d820e2d00fbc258bf6ea2182a0a0aa53af766a2a5d \
    --scan-run "$evidence_base/mtp-runs/F-511c891/additive/tier12-phenotype/scan835" \
    --scan-sha256 b6307d6da48cfc552df5513ad9ff038574551773ab5d7c24cd273beeac5ac1ea \
    --matrix-run "$evidence_base/mtp-data/codex-1173-main-merge/0fd7e54d0c414447a951231e79b400b4eddc68f3/phenotype" \
    --matrix-sha256 551a73cb246856518b79aa0aea0941156b01ee78a5e9a1d3e094de45bf02687f \
    --trait-list methods/trait_lists/headline_ten_phenotype.txt \
    --code-commit "$(git rev-parse HEAD)" \
    --figure "$figure_name" --output "$figure_output"
