"""Aggregate taxon trait counts by BugSigDB study profile."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import yaml

from bsdbng.count_inputs import (
    COMPLETE_COUNT_MODE,
    COUNT_INPUT_CONTRACT_VERSION,
    COUNT_INPUT_MODES,
    COVERAGE_COLUMNS,
    PROFILE_KEY_COLUMNS,
    CountInputMode,
    ProfileCoverage,
    annotation_coverage_status,
    taxa_trait_count_reader,
)
from bsdbng.paths import DERIVED_DATA_DIR, STUDY_DATA_DIR

DEFAULT_TAXA_TRAIT_MATRIX_PATH = DERIVED_DATA_DIR / "taxa_trait_matrix.tsv"
DEFAULT_OUTPUT_PATH = DERIVED_DATA_DIR / "study_trait_profiles_wide.tsv"
DEFAULT_COVERAGE_PATH = DERIVED_DATA_DIR / "study_trait_profiles_wide_coverage.tsv"
DEFAULT_PROVENANCE_PATH = DERIVED_DATA_DIR / "study_trait_profiles_wide_provenance.json"
PROFILE_ID_SEPARATOR = "_"
PROFILE_FORMATS = ("wide", "long")

ProfileFormat = Literal["wide", "long"]
ProfileKey = tuple[str, str, str]


@dataclass(frozen=True)
class TaxaTraitMatrix:
    """Sparse taxon-to-trait count lookup loaded from a wide TSV matrix."""

    trait_columns: list[str]
    counts_by_taxon: dict[str, Counter[str]]
    count_input_mode: CountInputMode
    blank_counts_interpreted_as_zero: int


@dataclass(frozen=True)
class ProfileAggregation:
    """Profile counts plus source-coverage information retained before analysis."""

    profile_counts: dict[ProfileKey, Counter[str]]
    coverage_by_profile: dict[ProfileKey, ProfileCoverage]
    matched_profile_taxa: int
    unmatched_profile_taxa: int


@dataclass(frozen=True)
class StudyTraitProfileStats:
    """Summary of a generated study trait profile table."""

    profiles: int
    trait_columns: int
    taxon_observations: int
    profile_taxa: int
    matched_profile_taxa: int
    unmatched_profile_taxa: int
    coverage_reported_distinct_taxa: int
    coverage_lookup_matched_taxa: int
    coverage_lookup_unmatched_taxa: int
    coverage_matched_taxa_with_annotations: int
    coverage_matched_taxa_without_annotations: int
    matrix_taxa: int
    output_path: Path
    coverage_output_path: Path
    provenance_output_path: Path
    study_dir: Path
    taxa_trait_matrix_path: Path
    output_format: ProfileFormat
    count_input_mode: CountInputMode
    blank_counts_interpreted_as_zero: int


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


def read_taxa_trait_matrix(
    path: Path, *, count_input_mode: CountInputMode = COMPLETE_COUNT_MODE
) -> TaxaTraitMatrix:
    """Read a complete taxa-by-traits count matrix keyed by unique taxon ID."""
    reader = taxa_trait_count_reader(path, mode=count_input_mode)
    trait_columns = list(reader.count_columns)
    counts_by_taxon: dict[str, Counter[str]] = {}
    for row in reader:
        counts_by_taxon[row.key[0]] = Counter(
            {
                column: count
                for column, count in zip(trait_columns, row.counts, strict=True)
                if count
            }
        )
    return TaxaTraitMatrix(
        trait_columns=trait_columns,
        counts_by_taxon=counts_by_taxon,
        count_input_mode=count_input_mode,
        blank_counts_interpreted_as_zero=reader.summary.blank_counts_interpreted_as_zero,
    )


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
    """Sum taxon trait vectors, retaining the historical three-value API."""
    result = aggregate_profile_counts_with_coverage(profile_taxa, matrix)
    return result.profile_counts, result.matched_profile_taxa, result.unmatched_profile_taxa


def _trait_filter_id(trait_columns: list[str]) -> str:
    content = "\n".join(trait_columns).encode()
    digest = hashlib.sha256(content).hexdigest()[:16]
    return f"columns:{len(trait_columns)}:sha256:{digest}"


def aggregate_profile_counts_with_coverage(
    profile_taxa: dict[ProfileKey, list[str]],
    matrix: TaxaTraitMatrix,
    *,
    coverage_trait_columns: list[str] | None = None,
) -> ProfileAggregation:
    """Sum taxon vectors and retain lookup/annotation coverage per profile."""
    trait_columns = (
        matrix.trait_columns if coverage_trait_columns is None else coverage_trait_columns
    )
    if not trait_columns:
        raise ValueError("Coverage requires at least one declared trait column")
    if len(set(trait_columns)) != len(trait_columns):
        raise ValueError("Coverage trait columns must be unique")
    unknown = sorted(set(trait_columns).difference(matrix.trait_columns))
    if unknown:
        raise ValueError(f"Coverage trait columns are absent from the matrix: {unknown}")
    trait_filter_id = _trait_filter_id(trait_columns)

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

    coverage_by_profile: dict[ProfileKey, ProfileCoverage] = {}
    for key, taxon_ids in profile_taxa.items():
        distinct_ids = sorted(set(taxon_ids))
        unmatched_ids: list[str] = []
        without_annotation_ids: list[str] = []
        with_annotations = 0
        for taxon_id in distinct_ids:
            taxon_counts = matrix.counts_by_taxon.get(taxon_id)
            if taxon_counts is None:
                unmatched_ids.append(taxon_id)
            elif any(taxon_counts.get(column, 0) > 0 for column in trait_columns):
                with_annotations += 1
            else:
                without_annotation_ids.append(taxon_id)
        matched = len(distinct_ids) - len(unmatched_ids)
        without_annotations = len(without_annotation_ids)
        status = annotation_coverage_status(
            reported=len(distinct_ids),
            matched=matched,
            unmatched=len(unmatched_ids),
            with_annotations=with_annotations,
            without_annotations=without_annotations,
        )
        coverage_by_profile[key] = ProfileCoverage(
            study_id=key[0],
            experiment_id=key[1],
            inc_dec=key[2],
            reported_distinct_taxa=len(distinct_ids),
            lookup_matched_taxa=matched,
            lookup_unmatched_taxa=len(unmatched_ids),
            matched_taxa_with_annotations=with_annotations,
            matched_taxa_without_annotations=without_annotations,
            coverage_status=status,
            unmatched_taxon_ids=tuple(unmatched_ids),
            matched_without_annotation_taxon_ids=tuple(without_annotation_ids),
            trait_filter_id=trait_filter_id,
        )

    return ProfileAggregation(
        profile_counts=profile_counts,
        coverage_by_profile=coverage_by_profile,
        matched_profile_taxa=matched_profile_taxa,
        unmatched_profile_taxa=unmatched_profile_taxa,
    )


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


def write_profile_coverage(
    output_path: Path, coverage_by_profile: dict[ProfileKey, ProfileCoverage]
) -> None:
    """Write source coverage separately from numeric profile counts."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COVERAGE_COLUMNS)
        for key in sorted(coverage_by_profile):
            coverage = coverage_by_profile[key]
            writer.writerow(
                [
                    coverage.study_id,
                    coverage.experiment_id,
                    coverage.inc_dec,
                    coverage.reported_distinct_taxa,
                    coverage.lookup_matched_taxa,
                    coverage.lookup_unmatched_taxa,
                    coverage.matched_taxa_with_annotations,
                    coverage.matched_taxa_without_annotations,
                    coverage.coverage_status,
                    json.dumps(list(coverage.unmatched_taxon_ids), separators=(",", ":")),
                    json.dumps(
                        list(coverage.matched_without_annotation_taxon_ids),
                        separators=(",", ":"),
                    ),
                    coverage.trait_filter_id,
                ]
            )


def _sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_profile_provenance(
    output_path: Path,
    *,
    matrix: TaxaTraitMatrix,
    profile_output_path: Path,
    coverage_output_path: Path,
    coverage_by_profile: dict[ProfileKey, ProfileCoverage],
    taxa_trait_matrix_path: Path,
) -> None:
    """Record the count mode, exact inputs, output hashes, and coverage totals."""
    coverage_rows = list(coverage_by_profile.values())
    trait_filter_ids = sorted({row.trait_filter_id for row in coverage_rows}) or [
        _trait_filter_id(matrix.trait_columns)
    ]
    if len(trait_filter_ids) != 1:
        raise ValueError(f"Expected one declared coverage trait filter, got {trait_filter_ids}")
    payload = {
        "count_input_contract": COUNT_INPUT_CONTRACT_VERSION,
        "taxa_trait_matrix": {
            "path": str(taxa_trait_matrix_path),
            "sha256": _sha256(taxa_trait_matrix_path),
            "rows": len(matrix.counts_by_taxon),
            "trait_columns": len(matrix.trait_columns),
            "count_input_mode": matrix.count_input_mode,
            "blank_counts_interpreted_as_zero": matrix.blank_counts_interpreted_as_zero,
        },
        "profile_output": {
            "path": str(profile_output_path),
            "sha256": _sha256(profile_output_path),
        },
        "coverage_output": {
            "path": str(coverage_output_path),
            "sha256": _sha256(coverage_output_path),
            "profile_key": list(PROFILE_KEY_COLUMNS),
            "trait_filter_id": trait_filter_ids[0],
            "status_is_separate_from": [
                "observed_numeric_zero",
                "analytical_exclusion",
                "test_eligibility",
                "display_masking",
                "statistical_significance",
            ],
            "totals": {
                column: sum(getattr(row, column) for row in coverage_rows)
                for column in (
                    "reported_distinct_taxa",
                    "lookup_matched_taxa",
                    "lookup_unmatched_taxa",
                    "matched_taxa_with_annotations",
                    "matched_taxa_without_annotations",
                )
            },
        },
        "zero_semantics": (
            "A numeric zero is an observed count under the declared trait matrix. "
            "Lookup and annotation coverage are recorded only in the coverage sidecar."
        ),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def create_study_trait_profiles(
    study_dir: Path = STUDY_DATA_DIR,
    taxa_trait_matrix_path: Path = DEFAULT_TAXA_TRAIT_MATRIX_PATH,
    output_path: Path = DEFAULT_OUTPUT_PATH,
    output_format: ProfileFormat = "wide",
    deduplicate_taxa: bool = True,
    count_input_mode: CountInputMode = COMPLETE_COUNT_MODE,
    coverage_output_path: Path | None = None,
    provenance_output_path: Path | None = None,
) -> StudyTraitProfileStats:
    """Create a study profile TSV by summing taxa-trait vectors from YAML taxa."""
    if output_format not in PROFILE_FORMATS:
        raise ValueError(f"Unsupported output format: {output_format}")
    if not study_dir.exists():
        raise FileNotFoundError(f"Study YAML directory not found: {study_dir}")
    if not taxa_trait_matrix_path.exists():
        raise FileNotFoundError(f"Taxa trait matrix TSV not found: {taxa_trait_matrix_path}")

    resolved_coverage_path = coverage_output_path or output_path.with_name(
        f"{output_path.stem}_coverage.tsv"
    )
    resolved_provenance_path = provenance_output_path or output_path.with_name(
        f"{output_path.stem}_provenance.json"
    )
    matrix = read_taxa_trait_matrix(taxa_trait_matrix_path, count_input_mode=count_input_mode)
    profile_taxa, taxon_observations = collect_profile_taxa(
        study_dir,
        deduplicate_taxa=deduplicate_taxa,
    )
    aggregation = aggregate_profile_counts_with_coverage(profile_taxa, matrix)

    if output_format == "wide":
        write_wide_profiles(output_path, aggregation.profile_counts, matrix.trait_columns)
    else:
        write_long_profiles(output_path, aggregation.profile_counts, matrix.trait_columns)
    write_profile_coverage(resolved_coverage_path, aggregation.coverage_by_profile)
    write_profile_provenance(
        resolved_provenance_path,
        matrix=matrix,
        profile_output_path=output_path,
        coverage_output_path=resolved_coverage_path,
        coverage_by_profile=aggregation.coverage_by_profile,
        taxa_trait_matrix_path=taxa_trait_matrix_path,
    )

    coverage_rows = list(aggregation.coverage_by_profile.values())

    return StudyTraitProfileStats(
        profiles=len(aggregation.profile_counts),
        trait_columns=len(matrix.trait_columns),
        taxon_observations=taxon_observations,
        profile_taxa=sum(len(taxon_ids) for taxon_ids in profile_taxa.values()),
        matched_profile_taxa=aggregation.matched_profile_taxa,
        unmatched_profile_taxa=aggregation.unmatched_profile_taxa,
        coverage_reported_distinct_taxa=sum(row.reported_distinct_taxa for row in coverage_rows),
        coverage_lookup_matched_taxa=sum(row.lookup_matched_taxa for row in coverage_rows),
        coverage_lookup_unmatched_taxa=sum(row.lookup_unmatched_taxa for row in coverage_rows),
        coverage_matched_taxa_with_annotations=sum(
            row.matched_taxa_with_annotations for row in coverage_rows
        ),
        coverage_matched_taxa_without_annotations=sum(
            row.matched_taxa_without_annotations for row in coverage_rows
        ),
        matrix_taxa=len(matrix.counts_by_taxon),
        output_path=output_path,
        coverage_output_path=resolved_coverage_path,
        provenance_output_path=resolved_provenance_path,
        study_dir=study_dir,
        taxa_trait_matrix_path=taxa_trait_matrix_path,
        output_format=output_format,
        count_input_mode=count_input_mode,
        blank_counts_interpreted_as_zero=matrix.blank_counts_interpreted_as_zero,
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
    parser.add_argument(
        "--count-input-mode",
        choices=COUNT_INPUT_MODES,
        default=COMPLETE_COUNT_MODE,
        help=(
            "Count table format. 'complete' rejects blanks; 'legacy_sparse_zero' "
            "explicitly interprets in-width blanks as zeros and records their count."
        ),
    )
    parser.add_argument(
        "--coverage-output",
        type=Path,
        help="Profile coverage sidecar path. Defaults beside --output.",
    )
    parser.add_argument(
        "--provenance-output",
        type=Path,
        help="Count-contract provenance JSON path. Defaults beside --output.",
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
        count_input_mode=args.count_input_mode,
        coverage_output_path=args.coverage_output,
        provenance_output_path=args.provenance_output,
    )

    print(f"Wrote {stats.output_format} study trait profiles: {stats.output_path}")
    print(f"Wrote profile coverage: {stats.coverage_output_path}")
    print(f"Wrote count provenance: {stats.provenance_output_path}")
    print(f"Profiles: {stats.profiles}")
    print(f"Trait columns: {stats.trait_columns}")
    print(f"Taxon observations in YAML: {stats.taxon_observations}")
    print(f"Profile taxa summed: {stats.profile_taxa}")
    print(f"Matched profile taxa: {stats.matched_profile_taxa}")
    print(f"Unmatched profile taxa: {stats.unmatched_profile_taxa}")
    print(f"Distinct profile taxa in coverage: {stats.coverage_reported_distinct_taxa}")
    print(f"Lookup-matched distinct taxa: {stats.coverage_lookup_matched_taxa}")
    print(f"Lookup-unmatched distinct taxa: {stats.coverage_lookup_unmatched_taxa}")
    print(f"Matched distinct taxa with annotations: {stats.coverage_matched_taxa_with_annotations}")
    print(
        "Matched distinct taxa without annotations: "
        f"{stats.coverage_matched_taxa_without_annotations}"
    )
    print(f"Count input mode: {stats.count_input_mode}")
    print(f"Blank counts interpreted as zero: {stats.blank_counts_interpreted_as_zero}")
    print(f"Matrix taxa: {stats.matrix_taxa}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
