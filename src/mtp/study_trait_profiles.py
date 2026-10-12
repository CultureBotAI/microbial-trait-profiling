"""Aggregate taxon trait counts by BugSigDB study profile."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import yaml

from mtp.paths import DERIVED_DATA_DIR, STUDY_DATA_DIR

DEFAULT_TAXA_TRAIT_MATRIX_PATH = DERIVED_DATA_DIR / "taxa_trait_matrix.tsv"
DEFAULT_OUTPUT_PATH = DERIVED_DATA_DIR / "study_trait_profiles_wide.tsv"
DEFAULT_ORA_OUTPUT_PATH = DERIVED_DATA_DIR / "study_trait_profiles_ora_wide.tsv"
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
    host_species: str | None = None,
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

            if host_species is not None:
                host = _as_text(experiment.get("host_species")).strip()
                if host != host_species.strip():
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
    *,
    binary: bool = False,
) -> tuple[dict[ProfileKey, Counter[str]], int, int]:
    """Aggregate taxon trait presence or counts for each study profile.

    When binary=True, counts the number of distinct taxa in the profile
    bearing each trait (0 or 1 per taxon). When binary=False (default), sums raw assertion counts.
    """
    profile_counts = {key: Counter[str]() for key in profile_taxa}
    matched_profile_taxa = 0
    unmatched_profile_taxa = 0

    for key, taxon_ids in profile_taxa.items():
        # In binary mode, evaluate distinct taxa so repeated taxa never increment binary presence
        taxa_to_process = list(dict.fromkeys(taxon_ids)) if binary else taxon_ids
        for taxon_id in taxa_to_process:
            taxon_counts = matrix.counts_by_taxon.get(taxon_id)
            if taxon_counts is None:
                unmatched_profile_taxa += 1
                continue

            matched_profile_taxa += 1
            if binary:
                # Binary indicator: each distinct taxon contributes 1 if it has the trait
                for trait, count in taxon_counts.items():
                    if count > 0:
                        profile_counts[key][trait] += 1
            else:
                profile_counts[key].update(taxon_counts)

    return profile_counts, matched_profile_taxa, unmatched_profile_taxa


def write_wide_profiles(
    output_path: Path,
    profile_counts: dict[ProfileKey, Counter[str]],
    trait_columns: list[str],
    profile_taxa_counts: dict[ProfileKey, int] | None = None,
) -> None:
    """Write one row per study, experiment, and inc/dec profile."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        extra_cols = ["taxa_count"] if profile_taxa_counts is not None else []
        writer.writerow(["study_id", "experiment_id", "inc_dec", *extra_cols, *trait_columns])
        for key in sorted(profile_counts):
            counts = profile_counts[key]
            taxa_val = [str(profile_taxa_counts[key])] if profile_taxa_counts is not None else []
            writer.writerow(
                [*key, *taxa_val, *[str(counts.get(column, 0)) for column in trait_columns]]
            )


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
    output_path: Path | None = None,
    output_format: ProfileFormat = "wide",
    deduplicate_taxa: bool = True,
    binary: bool = False,
    include_taxa_count: bool = False,
    host_species: str | None = None,
) -> StudyTraitProfileStats:
    """Create a study profile TSV by aggregating taxa-trait presence from YAML taxa."""
    if output_format not in PROFILE_FORMATS:
        raise ValueError(f"Unsupported output format: {output_format}")
    if include_taxa_count and output_format == "long":
        raise ValueError("--include-taxa-count is not supported with long format output")
    if not deduplicate_taxa and (binary or include_taxa_count):
        raise ValueError(
            "--count-duplicates cannot be combined with binary=True or include_taxa_count=True "
            "because the urn hypergeometric model requires distinct taxa drawn without replacement."
        )
    if not study_dir.exists():
        raise FileNotFoundError(f"Study YAML directory not found: {study_dir}")
    if not taxa_trait_matrix_path.exists():
        raise FileNotFoundError(f"Taxa trait matrix TSV not found: {taxa_trait_matrix_path}")

    is_ora = binary or include_taxa_count or (host_species is not None)
    if output_path is None:
        output_path = _default_output_path(output_format, is_ora=is_ora)
    elif output_path == DEFAULT_OUTPUT_PATH and is_ora:
        output_path = DEFAULT_ORA_OUTPUT_PATH

    matrix = read_taxa_trait_matrix(taxa_trait_matrix_path)
    profile_taxa, taxon_observations = collect_profile_taxa(
        study_dir,
        deduplicate_taxa=deduplicate_taxa,
        host_species=host_species,
    )
    profile_counts, matched_profile_taxa, unmatched_profile_taxa = aggregate_profile_counts(
        profile_taxa,
        matrix,
        binary=binary,
    )

    if output_format == "wide":
        profile_taxa_counts = None
        if include_taxa_count:
            # Issue #350 & #386: N must strictly count distinct matched signature
            # taxa in matrix universe M
            profile_taxa_counts = {
                key: sum(1 for taxon_id in set(taxa) if taxon_id in matrix.counts_by_taxon)
                for key, taxa in profile_taxa.items()
            }
        write_wide_profiles(
            output_path,
            profile_counts,
            matrix.trait_columns,
            profile_taxa_counts=profile_taxa_counts,
        )
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


def _default_output_path(output_format: ProfileFormat, *, is_ora: bool = False) -> Path:
    prefix = "study_trait_profiles_ora" if is_ora else "study_trait_profiles"
    return DERIVED_DATA_DIR / f"{prefix}_{output_format}.tsv"


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
    parser.add_argument(
        "--binary",
        action="store_true",
        default=False,
        help=(
            "Count binary taxon presence per trait (0 or 1 per taxon) "
            "rather than raw assertion counts."
        ),
    )
    parser.add_argument(
        "--include-taxa-count",
        action="store_true",
        default=False,
        help="Include matched signature taxa_count column in wide format output.",
    )
    parser.add_argument(
        "--host-species",
        type=str,
        default=None,
        help="Optional filter by host species (e.g. 'Homo sapiens').",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        stats = create_study_trait_profiles(
            study_dir=args.study_dir,
            taxa_trait_matrix_path=args.taxa_trait_matrix,
            output_path=args.output,
            output_format=args.format,
            deduplicate_taxa=not args.count_duplicates,
            binary=args.binary,
            include_taxa_count=args.include_taxa_count,
            host_species=args.host_species,
        )
    except ValueError as exc:
        parser.error(str(exc))

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
