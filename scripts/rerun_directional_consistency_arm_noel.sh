#!/usr/bin/env bash
# ==============================================================================
# Rerun Pipeline — Arm 3: Directional Consistency (Sign Test, Wilcoxon & LOO)
# Branch Owner: Noel (noel)
# ==============================================================================
# Usage:
#   bash scripts/rerun_directional_consistency_arm_noel.sh
#
# Notes:
#   - Runs canonical directional screen, study-flip null calibration, contrast
#     orientation, taxon influence, validation evidence, and integrated evidence
#     recipes on `main`.
#   - Also executes the standalone directional-consistency CLI scripts
#     (`scripts/run_directional_trait_analysis.py`,
#     `scripts/run_validation_analyses.py`, `scripts/generate_manuscript_figures.py`)
#     when invoked on the `noel` branch.
# ==============================================================================

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

echo "=== [Noel / Directional Consistency Arm] Checking prerequisite derived tables ==="
for required_file in \
    "data/raw/download_provenance.json" \
    "data/derived/study_trait_profiles_wide.tsv" \
    "data/derived/taxa_trait_matrix.tsv" \
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

echo "=== [1/9] Running full-manifest study sign-test directional screen (#255, #265) ==="
just trait-direction-screen

echo "=== [2/9] Generating shared study-flip null assignment artefact (#263) ==="
just calibration-flip-assignments

echo "=== [3/9] Running widened directional screen under study-flip null (#263) ==="
just calibration-screen

echo "=== [4/9] Running empirical directional calibration report (null + power replicates, #182) ==="
just calibration-report

echo "=== [5/9] Classifying experiment contrast orientation & disease-vs-healthy subset (#188) ==="
just contrast-orientation

echo "=== [6/9] Computing source-count reconciliation & genus/family taxon influence (#190) ==="
just table1-taxon-influence-from "data/kg-microbe" "data/derived" "data/derived/table1_taxon_influence"

echo "=== [7/9] Inventorying validation evidence & building integrated cross-arm evidence table (#178, #191) ==="
just validation-evidence
just integrated-evidence

echo "=== [8/9] Regenerating unpooled Autism/ASD/Autoimmune follow-up (#128) ==="
uv run python scripts/investigate_facultative_aerobe_condition_followup.py

if [[ -f "scripts/run_directional_trait_analysis.py" ]]; then
    echo "=== [9/9] Running standalone noel-branch directional trait analysis & LOO validation scripts ==="
    uv run python scripts/run_directional_trait_analysis.py
    if [[ -f "scripts/run_validation_analyses.py" ]]; then
        uv run python scripts/run_validation_analyses.py
    fi
    if [[ -f "scripts/generate_manuscript_figures.py" ]]; then
        uv run python scripts/generate_manuscript_figures.py
    fi
else
    echo "=== [9/9] Standalone 'scripts/run_directional_trait_analysis.py' not on this branch; canonical main directional outputs updated ==="
fi

echo "=== [Noel / Directional Consistency Arm] Rerun complete. Running verification checks ==="
just check
