"""Render #842 direction legends through the real CLI on fixtures (2026-10-23).

Inputs: synthetic three-study, three-trait profiles created here; no real data.
Outputs: isolated input TSV/CSV files and both PNG/PDF heatmaps for each basis.
Drawn-value and label validation belongs to the main()-level regression in
tests/test_study_trait_profile_heatmap.py, not to the basis-independent exports.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from matplotlib.figure import Figure


def write_inputs(root: Path) -> dict[str, Path]:
    """Create nonconstant profiles with all three required display traits."""
    root.mkdir(parents=True, exist_ok=True)
    paths = {
        "--study-trait-profiles-input": root / "profiles.tsv",
        "--trait-columns-input": root / "traits.tsv",
        "--studies-input": root / "studies.tsv",
        "--full-dump-input": root / "dump.csv",
    }
    traits = ["METPO:1000615", "METPO:1000698", "METPO:1000703"]
    values = [[8, 2, 1], [2, 5, 4], [4, 1, 8], [3, 7, 2], [1, 8, 3], [6, 2, 5]]
    profiles = []
    raw = []
    studies = []
    for index in range(3):
        study = f"bsdb:fixture{index + 1}"
        studies.append(
            {
                "study_id": study,
                "title": f"Synthetic cohort {index + 1}",
                "publication_year": "2026",
                "taxon_count": 3,
                "signature_count": 2,
            }
        )
        for signature, direction in enumerate(("inc", "dec"), 1):
            profiles.append(
                {
                    "study_id": study,
                    "experiment_id": f"{study}-1",
                    "inc_dec": direction,
                    **dict(zip(traits, values[index * 2 + signature - 1], strict=True)),
                }
            )
            raw.append(
                {
                    "BSDB ID": f"{study}/1/{signature}",
                    "Condition": f"Condition {index + 1}",
                    "Body site": "Feces",
                    "Location of subjects": "Synthetic",
                    "Host species": "Homo sapiens",
                    "Sequencing type": "16S",
                }
            )
    pd.DataFrame(profiles).to_csv(paths["--study-trait-profiles-input"], sep="\t", index=False)
    pd.DataFrame(studies).to_csv(paths["--studies-input"], sep="\t", index=False)
    pd.DataFrame(raw).to_csv(paths["--full-dump-input"], index=False)
    pd.DataFrame(
        {
            "column_id": traits,
            "object_name": ["mesophilic", "gram positive", "non motile"],
            "column_label": ["mesophilic", "gram positive", "non motile"],
            "object_category": ["biolink:PhenotypicQuality"] * 3,
            "subject_category": ["biolink:OrganismTaxon"] * 3,
            "object": traits,
        }
    ).to_csv(paths["--trait-columns-input"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    output = parser.parse_args().output_dir.resolve()
    source = Path(__file__).resolve().parents[1] / "create_study_trait_profile_heatmap.py"
    spec = importlib.util.spec_from_file_location("issue842_heatmap", source)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load profile heatmap renderer")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    inputs = write_inputs(output / "inputs")
    products = {
        "--output": "profiles.png",
        "--matrix-output": "matrix.tsv",
        "--study-metadata-output": "profiles.tsv",
        "--trait-metadata-output": "traits.tsv",
        "--paired-output": "paired.png",
        "--paired-matrix-output": "paired_matrix.tsv",
        "--paired-metadata-output": "paired_metadata.tsv",
        "--paired-required-trait-check-output": "required_traits.tsv",
    }
    original_save = Figure.savefig

    def save_both(figure: Figure, path: Path, *, dpi: int, bbox_inches: str) -> None:
        original_save(figure, path, dpi=dpi, bbox_inches=bbox_inches)
        original_save(figure, path.with_suffix(".pdf"), dpi=dpi, bbox_inches=bbox_inches)

    bases = ("group_1_vs_group_0", "disease_vs_comparison")
    for basis in bases:
        target = output / basis
        target.mkdir(parents=True, exist_ok=True)
        options = {**inputs, **{key: target / value for key, value in products.items()}}
        argv = [str(source), *(str(part) for item in options.items() for part in item)]
        argv += ["--direction-basis", basis, "--min-profile-support", "1"]
        with patch.object(sys, "argv", argv), patch.object(Figure, "savefig", save_both):
            module.main()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
