#!/usr/bin/env bash
# Purpose: render the 30 NCASM plates with the ORIGINAL, unported slide-wizard builders in
#   this checkout's Python environment, as the reference for the port in
#   src/mtp/ncasm_tufte_plates/ (compare with issue2264_compare_plate_renders.py).
# Inputs: SOURCE_DIR, the slide-wizard deck folder holding the six builder files, and
#   FIGURE_WIZARD_DIR, a figure-wizard checkout providing figure_wizard/tufte.py.
# Outputs: OUTPUT_DIR/tufte_plates/ (30 PNG and SVG, 20 PDF, manifest.json) and
#   OUTPUT_DIR/source_sha256.txt. Nothing is written into SOURCE_DIR.
# Issue #2264; 2026-11-11. Run from this repository's worktree root.
# ISSUE2264_BACKEND selects the matplotlib backend (default Agg, as the Justfile exports);
# the stored 2026-10-22 renders were drawn with the macOS backend, which resizes the
# 13.333 in figure to whole screen pixels (13.33 in), so only macosx reproduces them.
# Usage: bash scripts/adhoc/issue2264_render_original_builders.sh SOURCE_DIR FIGURE_WIZARD_DIR OUTPUT_DIR
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 SOURCE_DIR FIGURE_WIZARD_DIR OUTPUT_DIR" >&2
    exit 2
fi
repo_root="$(pwd)"
source_dir="$1"
figure_wizard_dir="$2"
output_dir="$3"

if [[ -e "${output_dir}" ]] && [[ -n "$(ls -A "${output_dir}")" ]]; then
    echo "Refusing to write into a non-empty directory: ${output_dir}" >&2
    exit 1
fi
mkdir -p "${output_dir}"
output_dir="$(cd "${output_dir}" && pwd)"

builders=(
    build_tufte_plates_20.py
    tufte_plates_common.py
    tufte_plates_part1.py
    tufte_plates_part2.py
    build_tufte_10_data_story_plots.py
    tufte_slides_common.py
)
for name in "${builders[@]}"; do
    cp "${source_dir}/${name}" "${output_dir}/${name}"
done
(cd "${output_dir}" && shasum -a 256 "${builders[@]}" > source_sha256.txt)

# The only edits: the three absolute local paths. The output base becomes this copy's
# own folder, and the figure-wizard checkout is the one named on the command line.
uv run python - "${output_dir}" "${figure_wizard_dir}" <<'PYTHON'
import re
import sys
from pathlib import Path

copy_dir = Path(sys.argv[1])
figure_wizard = str(Path(sys.argv[2]).resolve())
edits = {
    "build_tufte_plates_20.py": (
        r'Path\("/[^"]*/NCASM-disease-microbiome-traits"\)',
        "Path(__file__).resolve().parent",
    ),
    "tufte_plates_common.py": (r'Path\("/[^"]*/figure-wizard"\)', f"Path({figure_wizard!r})"),
    "tufte_slides_common.py": (r'Path\("/[^"]*/figure-wizard"\)', f"Path({figure_wizard!r})"),
}
for name, (pattern, replacement) in edits.items():
    path = copy_dir / name
    text, count = re.subn(pattern, lambda _match: replacement, path.read_text(encoding="utf-8"))
    if count != 1:
        raise SystemExit(f"{name}: expected one absolute path matching {pattern}, found {count}")
    path.write_text(text, encoding="utf-8")
PYTHON

(cd "${output_dir}" && MPLBACKEND="${ISSUE2264_BACKEND:-Agg}" MPLCONFIGDIR="${output_dir}/mplconfig" \
    uv run --project "${repo_root}" python build_tufte_plates_20.py)
test -s "${output_dir}/tufte_plates/manifest.json"
echo "Rendered $(find "${output_dir}/tufte_plates" -name '*.png' | wc -l | tr -d ' ') PNG files into ${output_dir}/tufte_plates"
