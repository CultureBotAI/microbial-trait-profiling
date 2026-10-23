"""Render #947 narrow-segment/dense-rank fixtures for design review (2026-10-23).

Inputs: synthetic summary frames shared with the geometry regressions.
Outputs: native-size PNG/SVG/PDF figures beneath a new --output-dir.
No production data or canonical tracked outputs are read or regenerated.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

from matplotlib.figure import Figure

from mtp.paths import REPO_ROOT


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    root = parser.parse_args().output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    spec = importlib.util.spec_from_file_location(
        "coverage_geometry_fixture", REPO_ROOT / "tests/test_taxa_coverage_histogram.py"
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load the coverage geometry fixtures")
    fixture = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixture)
    coverage = fixture.coverage
    original_save = coverage.save_figure

    def save_with_pdf(figure: Figure, path: Path) -> None:
        figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
        original_save(figure, path)

    coverage.save_figure = save_with_pdf
    for count in (78, 2, 998):
        ranks, taxa, profiles = fixture.narrow_segment_frames(count)
        coverage.plot_occurrences_vs_unique(
            ranks, taxa, profiles, root / f"zero_{count}_of_1000.png"
        )
    inputs = root / "inputs"
    inputs.mkdir()
    paths = fixture.inputs.__wrapped__(inputs)
    _, original, _, profiles, _, _ = coverage.load_and_process_data(**paths)
    coverage.plot_rank_stratified_histogram(
        fixture.dense_rank_frame(original), root / "dense_rank.png", profiles
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
