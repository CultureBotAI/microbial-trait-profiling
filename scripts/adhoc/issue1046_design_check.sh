#!/usr/bin/env bash
# Final-size checks for #1046, 2026-10-28. Inputs: figure bundle and installed CLI.
# Outputs: fresh durable audit directory. Source PDFs/PNGs stay read-only.
set -euo pipefail
if [ "$#" -ne 3 ]; then
    echo 'usage: issue1046_design_check.sh FIGURE_BUNDLE NEW_AUDIT_DIRECTORY DESIGN_CHECK' >&2
    exit 2
fi
bundle="$1"
audit="$2"
design_check="$3"
if [ ! -x "$design_check" ]; then
    echo 'Pass an already installed design-check executable' >&2
    exit 2
fi
uv run --no-sync python scripts/adhoc/issue1046_check_figure_bundle.py "$bundle" "$audit"
mkdir "$audit/design-check"
for pdf in "$bundle"/*/page-*.pdf; do
    [ -f "$pdf" ] || { echo 'No figure pages found' >&2; exit 2; }
    figure="$(basename "$(dirname "$pdf")")"
    page="$(basename "${pdf%.pdf}")"
    report="$audit/design-check/$figure-$page.yaml"
    if "$design_check" run --pdf "$pdf" --png "${pdf%.pdf}.png" \
        --target journal-double-column --out "$report" > "$audit/design-check/$figure-$page.log" 2>&1; then
        :
    else
        status=$?
        [ "$status" -eq 1 ] || exit "$status"
    fi
    "$design_check" validate "$report" >> "$audit/design-check/$figure-$page.log" 2>&1
done
uv run --no-sync python scripts/adhoc/issue1046_check_figure_bundle.py \
    "$bundle" "$audit/validated" --design-reports "$audit/design-check"
