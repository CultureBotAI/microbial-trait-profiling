"""Render coverage fixtures for the #912 design review, 2026-10-22.

Inputs: the repository's synthetic coverage fixture, with both taxa at species
rank to exercise mixed-profile labels. Outputs: native PNG/SVG, review PDFs and
the coverage tables/provenance in a new supplied directory. No production data.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import pandas as pd
from matplotlib.figure import Figure

from mtp.paths import REPO_ROOT


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.output_dir
    root.mkdir(parents=True, exist_ok=False)
    fixture_dir = root / "inputs"
    fixture_dir.mkdir()
    spec = importlib.util.spec_from_file_location(
        "coverage_review_fixture", REPO_ROOT / "tests/test_taxa_coverage_histogram.py"
    )
    if spec is None or spec.loader is None:
        raise ValueError("Coverage fixture module is unavailable")
    fixture = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixture)
    paths = fixture.inputs.__wrapped__(fixture_dir)
    matrix = pd.read_csv(paths["matrix"], sep="\t")
    matrix.loc[matrix.taxon_id.eq("t2"), "taxonomic_rank"] = "species"
    matrix.to_csv(paths["matrix"], sep="\t", index=False)
    coverage = fixture.coverage
    original = coverage.save_figure

    def with_pdf(fig: Figure, path: Path) -> None:
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
        original(fig, path)

    coverage.save_figure = with_pdf
    coverage.run_coverage(output_dir=root / "figures", **paths)


if __name__ == "__main__":
    main()
