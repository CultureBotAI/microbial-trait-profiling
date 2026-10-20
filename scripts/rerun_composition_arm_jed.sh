#!/usr/bin/env bash
# Purpose: rerun paired composition, sensitivity, calibration, audits and standard figures from a verified core bundle.
# Inputs: core run.json and fresh output directory; outputs: results and arm_run.json.
# Issue #658; 2026-10-19. Historical branch-specific investigations are separate.
# Usage: bash scripts/rerun_composition_arm_jed.sh /path/run.json /path/new-output [--no-figures]

set -euo pipefail
# Keep rendering headless independently of the caller's session (#719).
export MPLBACKEND=Agg

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 CORE_RUN_JSON FRESH_OUTPUT_DIR [--replicates N] [--seed N] [--no-figures]" >&2
    echo "Prepare shared inputs with 'just rerun-core'; see docs/workflow.md." >&2
    exit 2
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"
exec just rerun-jed "$@"
