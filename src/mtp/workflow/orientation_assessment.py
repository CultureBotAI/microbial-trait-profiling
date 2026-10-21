"""Describe an already oriented core population without relabelling it again (#658)."""

from __future__ import annotations

import argparse
import csv
import json
from collections.abc import Sequence
from dataclasses import asdict
from pathlib import Path

from mtp import contrast_orientation
from mtp.contrast_orientation import (
    CLASSIFIER_VERSION,
    DIRECTION_CONVENTION,
    class_counts,
    classify,
    comparison_rows,
    control_tier_counts,
    read_experiment_records,
    read_review,
    write_manifest,
)
from mtp.paths import file_sha256
from mtp.trait_guilds import TABLE1_TRAITS, load_profiles, read_bugsigdb_vintage
from mtp.workflow.consumers import open_core_bundle, preflight_consumer_outputs
from mtp.workflow.population import ORIENTED_DIRECTION_BASIS, POPULATIONS

ASSESSMENT_VERSION = "selected-core-orientation-assessment-v1"


def _write_table(path: Path, rows: Sequence[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_assessment(
    *,
    core_run: Path,
    profile_input: Path,
    experiments_input: Path,
    full_dump_input: Path,
    studies_input: Path,
    source_studies_input: Path,
    selection_input: Path,
    population: str,
    output_dir: Path,
    review_input: Path,
    code_commit: str,
) -> Path:
    """Bind raw classification evidence to analyses of already oriented profiles."""
    if population == "all_paired" or population not in POPULATIONS:
        raise ValueError("The selected-core assessment requires a disease-oriented population")
    bundle = open_core_bundle(core_run)
    if (
        bundle.run.contrast_population != population
        or bundle.run.direction_basis != ORIENTED_DIRECTION_BASIS
    ):
        raise ValueError("Assessment population or direction basis differs from verified core")
    bindings = {
        "profile_input": (profile_input, "count_profiles"),
        "experiments_input": (experiments_input, "experiments"),
        "full_dump_input": (full_dump_input, "raw_dump"),
        "studies_input": (studies_input, "studies"),
        "source_studies_input": (source_studies_input, "source_studies"),
        "selection_input": (selection_input, "contrast_selection"),
    }
    for argument, (supplied, product) in bindings.items():
        if supplied.resolve() != bundle.path(product).resolve():
            raise ValueError(f"--core-run conflicts with explicit --{argument.replace('_', '-')}")
    selection = json.loads(selection_input.read_text(encoding="utf-8"))
    if (
        selection["population"] != bundle.run.contrast_population
        or selection["output_direction_basis"] != bundle.run.direction_basis
    ):
        raise ValueError(
            "Selection evidence population or direction basis differs from verified core"
        )
    selected = {
        (row["study_id"], row["experiment_id"]) for row in selection["selected_experiments"]
    }
    records = read_experiment_records(full_dump_input)
    classes = {key: classify(record) for key, record in records.items()}
    profiles, stats = load_profiles(
        profile_input, experiments_input, TABLE1_TRAITS, reject_missing_counts=True
    )
    if stats.unmatched_experiments:
        raise ValueError("Selected profiles have missing experiment metadata")
    profile_keys = [
        (profile.study_id, profile.experiment_id, profile.inc_dec) for profile in profiles
    ]
    if len(profile_keys) != len(set(profile_keys)):
        raise ValueError("Selected profiles contain duplicate study/experiment/direction keys")
    corpus = sorted({(profile.study_id, profile.experiment_id) for profile in profiles})
    if not corpus:
        raise ValueError("Selected profiles contain no human paired experiments")
    if not set(corpus).issubset(selected):
        raise ValueError("Analysed profiles contain experiments outside the selected population")
    if set(corpus) != selected:
        raise ValueError("Analysed profiles omit experiments from the selected population")
    if stats.profile_rows != stats.kept_profiles:
        raise ValueError("Selected profiles contain rows outside their human paired population")
    if not set(corpus).issubset(classes):
        raise ValueError("Analysed profiles have no original raw classification")

    validation = read_review(review_input, classes=classes)
    corpus_review = {key: validation[key] for key in corpus if key in validation}
    vintage = read_bugsigdb_vintage(full_dump_input)
    counts = class_counts(classes, corpus)
    tiers = control_tier_counts(classes, corpus)
    # Selection has already relabelled reversed signatures before profile creation.
    # Passing these profiles through oriented_profiles would undo that relabelling.
    comparisons = comparison_rows([(population, profiles)])
    output = output_dir.resolve()
    if output_dir.is_symlink() or (
        output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir()))
    ):
        raise FileExistsError(f"Orientation assessment output already exists: {output_dir}")
    preflight_consumer_outputs(
        bundle, output / "contrast_orientation_provenance.json", output_directory=output
    )
    output.mkdir(parents=True, exist_ok=True)
    write_manifest(output / "contrast_manifest.tsv", records, classes, corpus, vintage, validation)
    _write_table(output / "contrast_class_counts.tsv", counts)
    _write_table(output / "contrast_control_tier_counts.tsv", tiers)
    _write_table(output / "disease_control_comparison.tsv", comparisons)
    provenance = {
        "assessment_version": ASSESSMENT_VERSION,
        "classifier_version": CLASSIFIER_VERSION,
        "direction_convention": (
            "inc = higher abundance in the disease group than in the comparison group"
        ),
        "raw_direction_convention": DIRECTION_CONVENTION,
        "manifest_direction_basis": selection["input_direction_basis"],
        "profile_direction_basis": selection["output_direction_basis"],
        "profile_directions_relabelled_by_assessment": False,
        "contrast_population": population,
        "analysis_population": (
            f"{population}; human paired 16S and shotgun profiles already oriented by the "
            "verified core selection; inc means higher in the disease group"
        ),
        "bugsigdb_vintage": vintage,
        "code_commit": code_commit,
        "core_execution": {
            "run_id": bundle.run.run_id,
            "manifest_path": str(bundle.manifest),
            "manifest_sha256": file_sha256(bundle.manifest),
            "execution_fingerprint": bundle.run.execution_fingerprint,
            "contrast_population": bundle.run.contrast_population,
            "direction_basis": bundle.run.direction_basis,
            "bound_products": [product for _, product in bindings.values()],
        },
        "full_dump_sha256": file_sha256(full_dump_input),
        "profile_sha256": file_sha256(profile_input),
        "experiments_sha256": file_sha256(experiments_input),
        "analysis_source_sha256": file_sha256(Path(__file__)),
        "classifier_source_sha256": file_sha256(Path(contrast_orientation.__file__)),
        "selection_evidence": {
            "path": str(selection_input.resolve()),
            "sha256": file_sha256(selection_input),
            "schema_version": selection["schema_version"],
            "population": selection["population"],
            "input_direction_basis": selection["input_direction_basis"],
            "output_direction_basis": selection["output_direction_basis"],
            "reversed_experiments": selection["reversed_experiments"],
        },
        "analysed_experiments": len(corpus),
        "analysed_studies": len({key[0] for key in corpus}),
        "profile_filter_counts": asdict(stats),
        "class_counts": counts,
        "control_tier_counts": tiers,
        "review": {
            "path": str(review_input.resolve()),
            "sha256": file_sha256(review_input) if review_input.is_file() else None,
            "scope": "hand-reviewed experiments intersecting the analysed human paired population",
            "reviewed_experiments": len(corpus_review),
            "agreed": sum(status == "reviewed_agreed" for status in corpus_review.values()),
            "disagreed": sum(status == "reviewed_disagreed" for status in corpus_review.values()),
            "reviewed_experiments_outside_population": len(validation) - len(corpus_review),
            "not_reviewed_experiments": len(corpus) - len(corpus_review),
        },
        "limitations": [
            "lexical rules over group names; not a curated re-reading of each paper",
            "the export has no Group 0 definition",
            "original raw group labels and classifications are retained in the manifest; "
            "the comparison table uses the already oriented selected profiles",
            "classification and checksum agreement do not establish biological validity",
            "hand-review evidence covers only the reviewed experiments in this population",
        ],
    }
    path = output / "contrast_orientation_provenance.json"
    path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for option in (
        "core-run",
        "profile-input",
        "experiments-input",
        "full-dump-input",
        "studies-input",
        "source-studies-input",
        "selection-input",
        "output-dir",
        "review-input",
    ):
        parser.add_argument(f"--{option}", type=Path, required=True)
    parser.add_argument("--population", choices=POPULATIONS, required=True)
    parser.add_argument("--code-commit", required=True)
    args = parser.parse_args(argv)
    path = write_assessment(**vars(args))
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
