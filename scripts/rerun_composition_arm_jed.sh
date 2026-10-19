#!/usr/bin/env bash
# ==============================================================================
# Rerun Pipeline — Arm 1: Paired Composition & Continuous Trait Profiling
# Branch Owner: Jed (jadhavjed-ux / jed)
# ==============================================================================
# Usage:
#   bash scripts/rerun_composition_arm_jed.sh
#
# Prerequisites:
#   - Run from repository root after `just study-trait-profiles` (or with canonical
#     `data/derived/study_trait_profiles_wide.tsv`,
#     `data/derived/taxa_trait_matrix.tsv`,
#     `data/derived/taxa_trait_matrix_columns.tsv`,
#     `data/derived/experiments.tsv`, and
#     `data/derived/bugsigdb_taxa.tsv` already present on disk).
#   - For Kaleido figure PNG export (`guild-maps`, `table1-trait-maps`), ensure
#     headless Chrome is installed (`just get-chrome`).
# ==============================================================================

set -euo pipefail
# Headless Agg, as every recipe uses, for the scripts this arm runs directly too (#719)
export MPLBACKEND=Agg

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

echo "=== [Jed / Paired Composition Arm] Checking prerequisite derived tables ==="
for required_file in \
    "data/raw/download_provenance.json" \
    "data/derived/study_trait_profiles_wide.tsv" \
    "data/derived/taxa_trait_matrix.tsv" \
    "data/derived/taxa_trait_matrix_columns.tsv" \
    "data/derived/experiments.tsv" \
    "data/derived/bugsigdb_taxa.tsv"; do
    if [[ ! -f "${required_file}" ]]; then
        echo "ERROR: Required input '${required_file}' not found." >&2
        echo "Run 'just study-trait-profiles' (after 'just download-kg', 'just pipeline', 'just study-tsvs', 'just taxa-trait-matrix') first." >&2
        exit 1
    fi
done
# The arms read the canonical human profile table committed as a .gz (#656).
just study-trait-profiles-check

echo "=== [1/9] Running paired within-study log-ratio analysis (15 Table 1 pairs) ==="
just table1-compositional-stats

echo "=== [2/9] Running predeclared log-ratio sensitivity scenarios (zero exclusions, fecal strata, LOO) ==="
just table1-log-ratio-sensitivity

echo "=== [3/9] Running paired composition empirical null calibration (100 replicates, seed 182) ==="
just calibrate-composition

echo "=== [4/9] Auditing annotation density confounding, hierarchy redundancy & PC1 adjustment ==="
just audit-trait-density

echo "=== [5/9] Auditing 6-trait biological coupling, partial correlations, VIF, LOTO & LCO ==="
just audit-trait-coupling

echo "=== [6/9] Regenerating clustered trait/condition/guild heatmaps & investigation outputs ==="
just trait-heatmaps
just trait-profile-investigation

echo "=== [7/9] Regenerating six-shared-trait Cytoscape heatmaps, oral view & Sankey diagrams ==="
just six-shared-trait-heatmaps
just oral-six-trait-heatmap
just six-shared-trait-combined-sankey
just six-shared-trait-high-support-sankeys

echo "=== [8/9] Regenerating Parkinson disease 6-trait taxon heatmaps, contributors & networks ==="
just pd-trait-taxa
just pd-trait-taxa-heatmap
just pd-trait-contributors
just pd-significant-ecological-network
just pd-ecological-reporter-combined

echo "=== [9/9] Rendering paired composition guild maps & Table 1 trait maps ==="
just guild-maps
just table1-trait-maps

echo "=== [Jed / Paired Composition Arm] Rerun complete. Running verification checks ==="
just check
