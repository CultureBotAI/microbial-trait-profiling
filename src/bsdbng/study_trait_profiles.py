"""Aggregate taxon trait counts by BugSigDB study profile."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import yaml

from bsdbng.paths import DERIVED_DATA_DIR, STUDY_DATA_DIR

DEFAULT_TAXA_TRAIT_MATRIX_PATH = DERIVED_DATA_DIR / "taxa_trait_matrix.tsv"
DEFAULT_OUTPUT_PATH = DERIVED_DATA_DIR / "study_trait_profiles_wide.tsv"
PROFILE_ID_SEPARATOR = "_"
PROFILE_FORMATS = ("wide", "long")

ProfileFormat = Literal["wide", "long"]
ProfileKey = tuple[str, str, str]


@dataclass(frozen=True)
class TaxaTraitMatrix:
    """Sparse taxon-to-trait count lookup loaded from a wide TSV matrix."""

    trait_columns: list[str]
    counts_by_taxon: dict[str, Counter[str]]


@dataclass(frozen=True)
class StudyTraitProfileStats:
    """Summary of a generated study trait profile table."""

    profiles: int
    trait_columns: int
    taxon_observations: int
    profile_taxa: int
    matched_profile_taxa: int
    unmatched_profile_taxa: int
    matrix_taxa: int
    output_path: Path
    study_dir: Path
    taxa_trait_matrix_path: Path
    output_format: ProfileFormat


def _as_text(value: Any) -> str:
    if value is None:
        return ""
    return str(value)


def _load_yaml(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Expected mapping at top level in {path}")
    return data


def direction_to_inc_dec(direction: str) -> str | None:
    """Return the compact output direction label for a BugSigDB signature direction."""
    normalized = direction.strip().lower()
    if normalized in {"increased", "inc"}:
        return "inc"
    if normalized in {"decreased", "dec"}:
        return "dec"
    return None


def profile_id(profile_key: ProfileKey) -> str:
    """Return a stable long-format column identifier for one profile row."""
    return PROFILE_ID_SEPARATOR.join(profile_key)


def _parse_count(value: str, *, path: Path, line_number: int, column_name: str) -> int:
    text = value.strip()
    if not text:
        return 0
    try:
        return int(text)
    except ValueError as exc:
        raise ValueError(
            f"Expected integer count in {path}:{line_number}, column {column_name!r}: {value!r}"
        ) from exc


def read_taxa_trait_matrix(path: Path) -> TaxaTraitMatrix:
    """Read a taxa-by-traits matrix keyed by NCBITaxon CURIE."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError(f"{path} is empty")
        if len(header) < 4:
            raise ValueError(f"{path} must have taxon_id, taxon_name, taxonomic_rank, and traits")
        if header[:3] != ["taxon_id", "taxon_name", "taxonomic_rank"]:
            raise ValueError(f"{path} first columns must be taxon_id, taxon_name, taxonomic_rank")

        trait_columns = header[3:]
        counts_by_taxon: dict[str, Counter[str]] = {}
        for line_number, row in enumerate(reader, start=2):
            if not row:
                continue

            taxon_id = row[0].strip()
            if not taxon_id:
                continue

            values = row[3:]
            if len(values) > len(trait_columns):
                raise ValueError(
                    f"{path}:{line_number} has {len(values)} trait values, "
                    f"but the header has {len(trait_columns)} trait columns"
                )

            counts: Counter[str] = Counter()
            for column_name, raw_value in zip(trait_columns, values, strict=False):
                count = _parse_count(
                    raw_value,
                    path=path,
                    line_number=line_number,
                    column_name=column_name,
                )
                if count:
                    counts[column_name] = count
            counts_by_taxon[taxon_id] = counts

    return TaxaTraitMatrix(trait_columns=trait_columns, counts_by_taxon=counts_by_taxon)


def collect_profile_taxa(
    study_dir: Path,
    *,
    deduplicate_taxa: bool = True,
) -> tuple[dict[ProfileKey, list[str]], int]:
    """Collect taxa from per-study YAML files by study, experiment, and direction."""
    profile_taxa: dict[ProfileKey, list[str]] = {}
    seen_taxa: defaultdict[ProfileKey, set[str]] = defaultdict(set)
    taxon_observations = 0

    for yaml_path in sorted(study_dir.glob("*.yaml")):
        study = _load_yaml(yaml_path)
        study_id = _as_text(study.get("id")).strip()
        experiments_obj = study.get("experiments")
        experiments = experiments_obj if isinstance(experiments_obj, list) else []

        for experiment in experiments:
            if not isinstance(experiment, dict):
                continue

            host = _as_text(experiment.get("host_species")).strip().lower()
            if host != "homo sapiens":
                continue

            experiment_id = _as_text(experiment.get("id")).strip()
            signatures_obj = experiment.get("signatures")
            signatures = signatures_obj if isinstance(signatures_obj, list) else []

            for signature in signatures:
                if not isinstance(signature, dict):
                    continue

                inc_dec = direction_to_inc_dec(_as_text(signature.get("direction")))
                if inc_dec is None:
                    continue

                key = (study_id, experiment_id, inc_dec)
                profile_taxa.setdefault(key, [])

                taxa_obj = signature.get("taxa")
                taxa = taxa_obj if isinstance(taxa_obj, list) else []
                for taxon in taxa:
                    if not isinstance(taxon, dict):
                        continue

                    taxon_id = _as_text(taxon.get("id")).strip()
                    if not taxon_id:
                        continue

                    taxon_observations += 1
                    if deduplicate_taxa and taxon_id in seen_taxa[key]:
                        continue

                    seen_taxa[key].add(taxon_id)
                    profile_taxa[key].append(taxon_id)

    return profile_taxa, taxon_observations


def aggregate_profile_counts(
    profile_taxa: dict[ProfileKey, list[str]],
    matrix: TaxaTraitMatrix,
) -> tuple[dict[ProfileKey, Counter[str]], int, int]:
    """Sum taxon trait vectors for each study profile."""
    profile_counts = {key: Counter[str]() for key in profile_taxa}
    matched_profile_taxa = 0
    unmatched_profile_taxa = 0

    for key, taxon_ids in profile_taxa.items():
        for taxon_id in taxon_ids:
            taxon_counts = matrix.counts_by_taxon.get(taxon_id)
            if taxon_counts is None:
                unmatched_profile_taxa += 1
                continue

            matched_profile_taxa += 1
            profile_counts[key].update(taxon_counts)

    return profile_counts, matched_profile_taxa, unmatched_profile_taxa


def write_wide_profiles(
    output_path: Path,
    profile_counts: dict[ProfileKey, Counter[str]],
    trait_columns: list[str],
) -> None:
    """Write one row per study, experiment, and inc/dec profile."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["study_id", "experiment_id", "inc_dec", *trait_columns])
        for key in sorted(profile_counts):
            counts = profile_counts[key]
            writer.writerow([*key, *[str(counts.get(column, 0)) for column in trait_columns]])


def write_long_profiles(
    output_path: Path,
    profile_counts: dict[ProfileKey, Counter[str]],
    trait_columns: list[str],
) -> None:
    """Write one row per trait with study profile columns."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    profile_keys = sorted(profile_counts)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait", *[profile_id(key) for key in profile_keys]])
        for trait_column in trait_columns:
            writer.writerow(
                [
                    trait_column,
                    *[str(profile_counts[key].get(trait_column, 0)) for key in profile_keys],
                ]
            )


def create_study_trait_profiles(
    study_dir: Path = STUDY_DATA_DIR,
    taxa_trait_matrix_path: Path = DEFAULT_TAXA_TRAIT_MATRIX_PATH,
    output_path: Path = DEFAULT_OUTPUT_PATH,
    output_format: ProfileFormat = "wide",
    deduplicate_taxa: bool = True,
) -> StudyTraitProfileStats:
    """Create a study profile TSV by summing taxa-trait vectors from YAML taxa."""
    if output_format not in PROFILE_FORMATS:
        raise ValueError(f"Unsupported output format: {output_format}")
    if not study_dir.exists():
        raise FileNotFoundError(f"Study YAML directory not found: {study_dir}")
    if not taxa_trait_matrix_path.exists():
        raise FileNotFoundError(f"Taxa trait matrix TSV not found: {taxa_trait_matrix_path}")

    matrix = read_taxa_trait_matrix(taxa_trait_matrix_path)
    profile_taxa, taxon_observations = collect_profile_taxa(
        study_dir,
        deduplicate_taxa=deduplicate_taxa,
    )
    profile_counts, matched_profile_taxa, unmatched_profile_taxa = aggregate_profile_counts(
        profile_taxa,
        matrix,
    )

    if output_format == "wide":
        write_wide_profiles(output_path, profile_counts, matrix.trait_columns)
    else:
        write_long_profiles(output_path, profile_counts, matrix.trait_columns)

    return StudyTraitProfileStats(
        profiles=len(profile_counts),
        trait_columns=len(matrix.trait_columns),
        taxon_observations=taxon_observations,
        profile_taxa=sum(len(taxon_ids) for taxon_ids in profile_taxa.values()),
        matched_profile_taxa=matched_profile_taxa,
        unmatched_profile_taxa=unmatched_profile_taxa,
        matrix_taxa=len(matrix.counts_by_taxon),
        output_path=output_path,
        study_dir=study_dir,
        taxa_trait_matrix_path=taxa_trait_matrix_path,
        output_format=output_format,
    )


def _default_output_path(output_format: ProfileFormat) -> Path:
    return DERIVED_DATA_DIR / f"study_trait_profiles_{output_format}.tsv"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate study trait profile TSVs by summing rows from "
            "data/derived/taxa_trait_matrix.tsv for taxa listed in per-study YAML."
        )
    )
    parser.add_argument(
        "--study-dir",
        type=Path,
        default=STUDY_DATA_DIR,
        help="Directory containing bsdbng per-study YAML files.",
    )
    parser.add_argument(
        "--taxa-trait-matrix",
        type=Path,
        default=DEFAULT_TAXA_TRAIT_MATRIX_PATH,
        help="Input taxa x traits TSV matrix.",
    )
    parser.add_argument(
        "--format",
        choices=PROFILE_FORMATS,
        default="wide",
        help="Output table shape.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output TSV path. Defaults to data/derived/study_trait_profiles_<format>.tsv.",
    )
    parser.add_argument(
        "--count-duplicates",
        action="store_true",
        help="Count repeated taxon IDs within the same study/experiment/direction profile.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    output_format = args.format
    output_path = args.output or _default_output_path(output_format)
    stats = create_study_trait_profiles(
        study_dir=args.study_dir,
        taxa_trait_matrix_path=args.taxa_trait_matrix,
        output_path=output_path,
        output_format=output_format,
        deduplicate_taxa=not args.count_duplicates,
    )

    print(f"Wrote {stats.output_format} study trait profiles: {stats.output_path}")
    print(f"Profiles: {stats.profiles}")
    print(f"Trait columns: {stats.trait_columns}")
    print(f"Taxon observations in YAML: {stats.taxon_observations}")
    print(f"Profile taxa summed: {stats.profile_taxa}")
    print(f"Matched profile taxa: {stats.matched_profile_taxa}")
    print(f"Unmatched profile taxa: {stats.unmatched_profile_taxa}")
    print(f"Matrix taxa: {stats.matrix_taxa}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
