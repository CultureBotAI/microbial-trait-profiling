#!/usr/bin/env bash
# ==============================================================================
# Rerun Pipeline — Arm 2: Over-Representation Analysis (ORA) Enrichment & Depletion
# Branch Owner: Ilyas (IlyasMoutawwakil / ilyas)
# ==============================================================================
# Usage:
#   bash scripts/rerun_enrichment_ora_arm_ilyas.sh
#
# Prerequisites:
#   - Run from repository root after `just taxa-trait-matrix` and
#     `just study-trait-profiles-ora` (or with canonical
#     `data/derived/taxa_trait_matrix.tsv` and `data/derived/experiments.tsv`
#     already present on disk).
# ==============================================================================

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

echo "=== [Ilyas / ORA Enrichment Arm] Checking prerequisite derived tables ==="
for required_file in \
    "data/derived/taxa_trait_matrix.tsv" \
    "data/derived/experiments.tsv"; do
    if [[ ! -f "${required_file}" ]]; then
        echo "ERROR: Required input '${required_file}' not found." >&2
        echo "Run 'just taxa-trait-matrix' (after 'just download-kg', 'just pipeline', 'just study-tsvs') first." >&2
        exit 1
    fi
done

echo "=== [1/6] Building binary Homo sapiens ORA study-trait profiles ==="
just study-trait-profiles-ora

echo "=== [2/6] Running global-FDR ORA trait enrichment & exporting rankings (#119) ==="
just analyze

echo "=== [3/6] Running global-FDR ORA trait depletion & exporting rankings (#119) ==="
just analyze-depletion

echo "=== [4/6] Running ORA background stratification & annotation-density sensitivity (#118) ==="
just analyze-strata

echo "=== [5/6] Running ORA synthetic null & spike-in calibration harness (#263, #343, #348, #351) ==="
just calibrate

echo "=== [6/6] Rendering publication deck visualizations (Figures 1-5) ==="
just deck-figures

echo "=== [Ilyas / ORA Enrichment Arm] Rerun complete. Running verification checks ==="
just check
