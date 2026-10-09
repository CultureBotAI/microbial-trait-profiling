"""Generate TSV tables from per-study BugSigDB YAML records."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from bsdbng.paths import DERIVED_DATA_DIR, STUDY_DATA_DIR


@dataclass(frozen=True)
class GeneratedTables:
    """Paths to the emitted TSV tables."""

    studies: Path
    experiments: Path
    experiment_ontology_terms: Path
    signatures: Path
    taxa: Path
    taxon_profiles_wide: Path
    taxon_profiles_long: Path


def _as_text(value: Any) -> str:
    if value is None:
        return ""
    return str(value)


def _ontology_terms(value: Any) -> list[str]:
    """Return the non-empty CURIEs held by a multivalued ontology-id slot."""
    if value is None:
        return []
    values = value if isinstance(value, list) else [value]
    return [text for item in values if (text := _as_text(item).strip())]


def _load_yaml(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        msg = f"Expected mapping at top level in {path}"
        raise ValueError(msg)
    return data


def _split_profile_direction(direction: str) -> str:
    if direction == "increased":
        return "inc"
    if direction == "decreased":
        return "dec"
    return direction


def _profile_feature(taxon: dict[str, Any], profile_key: str) -> str:
    preferred = _as_text(taxon.get(profile_key)).strip()
    if preferred:
        return preferred
    fallback = _as_text(taxon.get("taxon_name")).strip()
    if fallback:
        return fallback
    return _as_text(taxon.get("id")).strip()


def _write_dict_rows(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _as_text(row.get(field)) for field in fieldnames})


def _write_profile_wide(
    path: Path,
    profile_rows: dict[tuple[str, str, str], dict[str, int]],
    feature_names: list[str],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["study_id", "experiment_id", "inc_dec", *feature_names])
        for row_key in sorted(profile_rows):
            row = list(row_key) + [
                str(profile_rows[row_key].get(name, 0)) for name in feature_names
            ]
            writer.writerow(row)


def _write_profile_long(
    path: Path,
    profile_rows: dict[tuple[str, str, str], dict[str, int]],
    feature_names: list[str],
) -> None:
    row_keys = sorted(profile_rows)
    column_names = [
        f"{study_id}_{experiment_id}_{inc_dec}" for study_id, experiment_id, inc_dec in row_keys
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxon", *column_names])
        for feature_name in feature_names:
            row = [feature_name]
            for row_key in row_keys:
                row.append(str(profile_rows[row_key].get(feature_name, 0)))
            writer.writerow(row)


def generate_tables(
    study_dir: Path = STUDY_DATA_DIR,
    output_dir: Path = DERIVED_DATA_DIR,
    profile_key: str = "taxon_name",
) -> GeneratedTables:
    """Flatten YAML study records into multiple TSV tables."""
    output_dir.mkdir(parents=True, exist_ok=True)

    study_rows: list[dict[str, Any]] = []
    experiment_rows: list[dict[str, Any]] = []
    experiment_ontology_rows: list[dict[str, Any]] = []
    signature_rows: list[dict[str, Any]] = []
    taxon_rows: list[dict[str, Any]] = []
    profile_rows: dict[tuple[str, str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    feature_names: set[str] = set()

    for yaml_path in sorted(study_dir.glob("*.yaml")):
        study = _load_yaml(yaml_path)
        study_id = _as_text(study.get("id"))
        experiments_obj = study.get("experiments")
        experiments = experiments_obj if isinstance(experiments_obj, list) else []

        study_signature_count = 0
        study_taxon_count = 0

        for experiment in experiments:
            if not isinstance(experiment, dict):
                continue

            experiment_id = _as_text(experiment.get("id"))
            signatures_obj = experiment.get("signatures")
            signatures = signatures_obj if isinstance(signatures_obj, list) else []
            experiment_signature_count = 0
            experiment_taxon_count = 0

            for signature in signatures:
                if not isinstance(signature, dict):
                    continue

                signature_id = _as_text(signature.get("id"))
                direction = _as_text(signature.get("direction"))
                inc_dec = _split_profile_direction(direction)
                taxa_obj = signature.get("taxa")
                taxa = taxa_obj if isinstance(taxa_obj, list) else []

                experiment_signature_count += 1
                study_signature_count += 1
                experiment_taxon_count += len(taxa)
                study_taxon_count += len(taxa)

                signature_rows.append(
                    {
                        "study_id": study_id,
                        "experiment_id": experiment_id,
                        "signature_id": signature_id,
                        "direction": direction,
                        "taxon_count": len(taxa),
                    }
                )

                row_key = (study_id, experiment_id, inc_dec)
                for taxon in taxa:
                    if not isinstance(taxon, dict):
                        continue

                    taxon_id = _as_text(taxon.get("id"))
                    taxon_name = _as_text(taxon.get("taxon_name"))
                    taxonomic_rank = _as_text(taxon.get("taxonomic_rank"))
                    taxon_rows.append(
                        {
                            "study_id": study_id,
                            "experiment_id": experiment_id,
                            "signature_id": signature_id,
                            "direction": direction,
                            "taxon_id": taxon_id,
                            "taxon_name": taxon_name,
                            "taxonomic_rank": taxonomic_rank,
                        }
                    )

                    feature_name = _profile_feature(taxon, profile_key)
                    if feature_name:
                        profile_rows[row_key][feature_name] += 1
                        feature_names.add(feature_name)

            experiment_rows.append(
                {
                    "study_id": study_id,
                    "experiment_id": experiment_id,
                    "experiment_name": _as_text(experiment.get("experiment_name")),
                    "condition": _as_text(experiment.get("condition")),
                    "body_site": _as_text(experiment.get("body_site")),
                    "host_species": _as_text(experiment.get("host_species")),
                    "location_of_subjects": _as_text(experiment.get("location_of_subjects")),
                    "group_0_name": _as_text(experiment.get("group_0_name")),
                    "group_1_name": _as_text(experiment.get("group_1_name")),
                    "signature_count": experiment_signature_count,
                    "taxon_count": experiment_taxon_count,
                }
            )

            for slot in ("condition_ontology_id", "body_site_ontology_id"):
                for term_id in _ontology_terms(experiment.get(slot)):
                    experiment_ontology_rows.append(
                        {
                            "study_id": study_id,
                            "experiment_id": experiment_id,
                            "slot": slot,
                            "term_id": term_id,
                        }
                    )

        study_rows.append(
            {
                "study_id": study_id,
                "source_record_id": _as_text(study.get("source_record_id")),
                "pmid": _as_text(study.get("pmid")),
                "title": _as_text(study.get("title")),
                "publication_year": _as_text(study.get("publication_year")),
                "doi": _as_text(study.get("doi")),
                "url": _as_text(study.get("url")),
                "experiment_count": len(experiments),
                "signature_count": study_signature_count,
                "taxon_count": study_taxon_count,
                "input_file": yaml_path.name,
            }
        )

    feature_name_list = sorted(feature_names)
    tables = GeneratedTables(
        studies=output_dir / "studies.tsv",
        experiments=output_dir / "experiments.tsv",
        experiment_ontology_terms=output_dir / "experiment_ontology_terms.tsv",
        signatures=output_dir / "signatures.tsv",
        taxa=output_dir / "bugsigdb_taxa.tsv",
        taxon_profiles_wide=output_dir / "taxon_profiles_wide.tsv",
        taxon_profiles_long=output_dir / "taxon_profiles_long.tsv",
    )

    _write_dict_rows(
        tables.studies,
        [
            "study_id",
            "source_record_id",
            "pmid",
            "title",
            "publication_year",
            "doi",
            "url",
            "experiment_count",
            "signature_count",
            "taxon_count",
            "input_file",
        ],
        study_rows,
    )
    _write_dict_rows(
        tables.experiments,
        [
            "study_id",
            "experiment_id",
            "experiment_name",
            "condition",
            "body_site",
            "host_species",
            "location_of_subjects",
            "group_0_name",
            "group_1_name",
            "signature_count",
            "taxon_count",
        ],
        experiment_rows,
    )
    _write_dict_rows(
        tables.experiment_ontology_terms,
        ["study_id", "experiment_id", "slot", "term_id"],
        experiment_ontology_rows,
    )
    _write_dict_rows(
        tables.signatures,
        ["study_id", "experiment_id", "signature_id", "direction", "taxon_count"],
        signature_rows,
    )
    _write_dict_rows(
        tables.taxa,
        [
            "study_id",
            "experiment_id",
            "signature_id",
            "direction",
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
        ],
        taxon_rows,
    )
    _write_profile_wide(tables.taxon_profiles_wide, profile_rows, feature_name_list)
    _write_profile_long(tables.taxon_profiles_long, profile_rows, feature_name_list)
    return tables


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate TSV tables from BugSigDB per-study YAML files. "
            "Includes flattened tables plus wide/long taxon profile matrices."
        )
    )
    parser.add_argument(
        "--study-dir",
        type=Path,
        default=STUDY_DATA_DIR,
        help=f"Directory containing per-study YAML files (default: {STUDY_DATA_DIR})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DERIVED_DATA_DIR,
        help=f"Directory to write TSV outputs (default: {DERIVED_DATA_DIR})",
    )
    parser.add_argument(
        "--profile-key",
        choices=["taxon_name", "id"],
        default="taxon_name",
        help="Taxon field to use as columns/rows in the profile matrices.",
    )
    args = parser.parse_args()

    tables = generate_tables(args.study_dir, args.output_dir, args.profile_key)
    for path in (
        tables.studies,
        tables.experiments,
        tables.experiment_ontology_terms,
        tables.signatures,
        tables.taxa,
        tables.taxon_profiles_wide,
        tables.taxon_profiles_long,
    ):
        print(path)


if __name__ == "__main__":
    main()
