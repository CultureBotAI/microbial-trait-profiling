"""Strict, shared readers for taxon and profile count tables.

Count tables are data contracts, not forgiving import formats.  A missing,
fractional, negative, or malformed count must not become an observed zero before
the statistical methods see it.  The reader streams rows so the wide profile
matrix does not need to be materialized twice in memory.

Adapted for issue #629 from the unmerged #184 / PR #224 implementation
(commit 94039ed), preserving its count-input-v1 compatibility policy.
"""

from __future__ import annotations

import csv
import gzip
import json
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Literal

CountInputMode = Literal["complete", "legacy_sparse_zero"]

COMPLETE_COUNT_MODE: CountInputMode = "complete"
LEGACY_SPARSE_ZERO_MODE: CountInputMode = "legacy_sparse_zero"
COUNT_INPUT_MODES = (COMPLETE_COUNT_MODE, LEGACY_SPARSE_ZERO_MODE)
COUNT_INPUT_CONTRACT_VERSION = "count-input-v1"

PROFILE_KEY_COLUMNS = ("study_id", "experiment_id", "inc_dec")
TAXON_MATRIX_PREFIX = ("taxon_id", "taxon_name", "taxonomic_rank")

COVERAGE_COUNT_COLUMNS = (
    "reported_distinct_taxa",
    "lookup_matched_taxa",
    "lookup_unmatched_taxa",
    "matched_taxa_with_annotations",
    "matched_taxa_without_annotations",
)
COVERAGE_COLUMNS = (
    *PROFILE_KEY_COLUMNS,
    *COVERAGE_COUNT_COLUMNS,
    "coverage_status",
    "unmatched_taxon_ids",
    "matched_without_annotation_taxon_ids",
    "trait_filter_id",
    "lineage_resolution",
)

COVERAGE_STATUSES = (
    "no_reported_taxa",
    "no_lookup_matches",
    "partial_lookup_coverage",
    "matched_without_annotations",
    "partial_annotation_coverage",
    "complete_annotation_coverage",
)


def open_input_text(path: Path) -> IO[str]:
    """Open an uncompressed or gzip-compressed TSV as UTF-8 text."""
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def parse_nonnegative_integer(
    value: str,
    *,
    path: Path,
    line_number: int,
    column_name: str,
    mode: CountInputMode = COMPLETE_COUNT_MODE,
) -> int:
    """Parse one count without coercion, truncation, or implicit missingness."""
    if mode not in COUNT_INPUT_MODES:
        raise ValueError(f"Unsupported count input mode: {mode!r}")
    text = value.strip()
    if not text:
        if mode == LEGACY_SPARSE_ZERO_MODE:
            return 0
        raise ValueError(
            f"Missing count in {path}:{line_number}, column {column_name!r}; "
            "missing is not an observed zero"
        )
    if not text.isascii() or not text.isdecimal():
        raise ValueError(
            f"Expected a finite nonnegative integer count in {path}:{line_number}, "
            f"column {column_name!r}; got {value!r}"
        )
    return int(text)


@dataclass(frozen=True)
class CountTableRow:
    """One validated count-table row."""

    line_number: int
    key: tuple[str, ...]
    values: tuple[str, ...]
    counts: tuple[int, ...]


@dataclass(frozen=True)
class CountReadSummary:
    """Validation facts captured while consuming a count table."""

    rows: int
    count_columns: int
    blank_counts_interpreted_as_zero: int
    mode: CountInputMode


@dataclass(frozen=True)
class ProfileCoverage:
    """Source annotation coverage for one profile, separate from test eligibility."""

    study_id: str
    experiment_id: str
    inc_dec: str
    reported_distinct_taxa: int
    lookup_matched_taxa: int
    lookup_unmatched_taxa: int
    matched_taxa_with_annotations: int
    matched_taxa_without_annotations: int
    coverage_status: str
    unmatched_taxon_ids: tuple[str, ...]
    matched_without_annotation_taxon_ids: tuple[str, ...]
    trait_filter_id: str
    lineage_resolution: str = "unresolved"

    @property
    def key(self) -> tuple[str, str, str]:
        return self.study_id, self.experiment_id, self.inc_dec

    @property
    def is_complete(self) -> bool:
        return self.coverage_status == "complete_annotation_coverage"


def annotation_coverage_status(
    *,
    reported: int,
    matched: int,
    unmatched: int,
    with_annotations: int,
    without_annotations: int,
) -> str:
    """Classify source coverage without treating an analytical zero as missing."""
    if min(reported, matched, unmatched, with_annotations, without_annotations) < 0:
        raise ValueError("Coverage counts must be nonnegative")
    if matched + unmatched != reported:
        raise ValueError("Lookup-matched plus unmatched taxa must equal reported taxa")
    if with_annotations + without_annotations != matched:
        raise ValueError("Annotated plus unannotated matched taxa must equal matched taxa")
    if reported == 0:
        return "no_reported_taxa"
    if matched == 0:
        return "no_lookup_matches"
    if unmatched:
        return "partial_lookup_coverage"
    if with_annotations == 0:
        return "matched_without_annotations"
    if without_annotations:
        return "partial_annotation_coverage"
    return "complete_annotation_coverage"


class CountTableReader:
    """Validate and stream a TSV whose keys and count columns are declared.

    ``leading_columns`` must match the beginning of the header exactly.
    ``key_columns`` define row uniqueness and must occur in the header.  When
    ``count_columns`` is omitted, every column after ``leading_columns`` is a
    count.  Supplying it explicitly is useful for sidecars that also carry
    nonnumeric status/provenance fields.
    """

    def __init__(
        self,
        path: Path,
        *,
        leading_columns: Sequence[str],
        key_columns: Sequence[str],
        count_columns: Sequence[str] | None = None,
        exact_columns: Sequence[str] | None = None,
        mode: CountInputMode = COMPLETE_COUNT_MODE,
    ) -> None:
        if mode not in COUNT_INPUT_MODES:
            raise ValueError(f"Unsupported count input mode: {mode!r}")
        self.path = path
        self.leading_columns = tuple(leading_columns)
        self.key_columns = tuple(key_columns)
        self.mode = mode
        self._rows = 0
        self._blank_counts = 0

        with open_input_text(path) as handle:
            header = next(csv.reader(handle, delimiter="\t"), None)
        if header is None:
            raise ValueError(f"{path} is empty")
        self.header = tuple(header)
        self._validate_header(exact_columns)

        positions = {column: index for index, column in enumerate(self.header)}
        missing_keys = [column for column in self.key_columns if column not in positions]
        if missing_keys:
            raise ValueError(f"{path} is missing key columns: {missing_keys}")
        self._key_positions = tuple(positions[column] for column in self.key_columns)

        if count_columns is None:
            inferred = self.header[len(self.leading_columns) :]
            self.count_columns = tuple(inferred)
        else:
            self.count_columns = tuple(count_columns)
            missing_counts = [column for column in self.count_columns if column not in positions]
            if missing_counts:
                raise ValueError(f"{path} is missing count columns: {missing_counts}")
        if not self.count_columns:
            raise ValueError(f"{path} has no count columns")
        self._count_positions = tuple(positions[column] for column in self.count_columns)

    def _validate_header(self, exact_columns: Sequence[str] | None) -> None:
        if not self.header:
            raise ValueError(f"{self.path} has an empty header")
        for index, column in enumerate(self.header, start=1):
            if not column.strip():
                raise ValueError(f"{self.path} has a blank header at column {index}")
        seen: dict[str, int] = {}
        for index, column in enumerate(self.header, start=1):
            if column in seen:
                raise ValueError(
                    f"{self.path} has duplicate column {column!r} at columns "
                    f"{seen[column]} and {index}"
                )
            seen[column] = index
        actual_prefix = self.header[: len(self.leading_columns)]
        if actual_prefix != self.leading_columns:
            raise ValueError(
                f"{self.path} must start with {self.leading_columns}; got {actual_prefix}"
            )
        if exact_columns is not None and self.header != tuple(exact_columns):
            raise ValueError(
                f"{self.path} must have exactly these columns: {tuple(exact_columns)}; "
                f"got {self.header}"
            )

    @property
    def summary(self) -> CountReadSummary:
        return CountReadSummary(
            rows=self._rows,
            count_columns=len(self.count_columns),
            blank_counts_interpreted_as_zero=self._blank_counts,
            mode=self.mode,
        )

    def __iter__(self) -> Iterator[CountTableRow]:
        self._rows = 0
        self._blank_counts = 0
        seen_keys: dict[tuple[str, ...], int] = {}
        with open_input_text(self.path) as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = tuple(next(reader, ()))
            if header != self.header:
                raise ValueError(f"{self.path} changed after its header was validated")
            for raw_row in reader:
                line_number = reader.line_num
                if len(raw_row) != len(self.header):
                    raise ValueError(
                        f"Wrong row width in {self.path}:{line_number}; expected "
                        f"{len(self.header)} fields, got {len(raw_row)}"
                    )
                values = tuple(raw_row)
                key_parts: list[str] = []
                for column, position in zip(self.key_columns, self._key_positions, strict=True):
                    value = values[position].strip()
                    if not value:
                        raise ValueError(
                            f"Missing identifier in {self.path}:{line_number}, column {column!r}"
                        )
                    key_parts.append(value)
                key = tuple(key_parts)
                if key in seen_keys:
                    raise ValueError(
                        f"Duplicate key {key!r} in {self.path}:{line_number}; "
                        f"first seen at row {seen_keys[key]}"
                    )
                seen_keys[key] = line_number

                parsed: list[int] = []
                for column, position in zip(self.count_columns, self._count_positions, strict=True):
                    value = values[position]
                    if not value.strip() and self.mode == LEGACY_SPARSE_ZERO_MODE:
                        self._blank_counts += 1
                    parsed.append(
                        parse_nonnegative_integer(
                            value,
                            path=self.path,
                            line_number=line_number,
                            column_name=column,
                            mode=self.mode,
                        )
                    )
                self._rows += 1
                yield CountTableRow(line_number, key, values, tuple(parsed))


def profile_count_reader(
    path: Path, *, mode: CountInputMode = COMPLETE_COUNT_MODE
) -> CountTableReader:
    """Return the shared strict reader for study profile count matrices."""
    return CountTableReader(
        path,
        leading_columns=PROFILE_KEY_COLUMNS,
        key_columns=PROFILE_KEY_COLUMNS,
        mode=mode,
    )


def taxa_trait_count_reader(
    path: Path, *, mode: CountInputMode = COMPLETE_COUNT_MODE
) -> CountTableReader:
    """Return the shared strict reader for taxon-by-trait count matrices."""
    return CountTableReader(
        path,
        leading_columns=TAXON_MATRIX_PREFIX,
        key_columns=("taxon_id",),
        mode=mode,
    )


def _parse_identifier_list(
    value: str, *, path: Path, line_number: int, column: str
) -> tuple[str, ...]:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Expected a JSON identifier list in {path}:{line_number}, column {column!r}"
        ) from exc
    if not isinstance(parsed, list) or any(
        not isinstance(item, str) or not item for item in parsed
    ):
        raise ValueError(f"Expected a JSON string list in {path}:{line_number}, column {column!r}")
    if parsed != sorted(set(parsed)):
        raise ValueError(
            f"Identifiers must be unique and sorted in {path}:{line_number}, column {column!r}"
        )
    return tuple(parsed)


def read_profile_coverage(path: Path) -> dict[tuple[str, str, str], ProfileCoverage]:
    """Read and reconcile a generated profile coverage sidecar."""
    reader = CountTableReader(
        path,
        leading_columns=PROFILE_KEY_COLUMNS,
        key_columns=PROFILE_KEY_COLUMNS,
        count_columns=COVERAGE_COUNT_COLUMNS,
        exact_columns=COVERAGE_COLUMNS,
    )
    positions = {column: index for index, column in enumerate(reader.header)}
    coverage: dict[tuple[str, str, str], ProfileCoverage] = {}
    for row in reader:
        counts = dict(zip(reader.count_columns, row.counts, strict=True))
        status = row.values[positions["coverage_status"]].strip()
        if status not in COVERAGE_STATUSES:
            raise ValueError(f"Unknown coverage status in {path}:{row.line_number}: {status!r}")
        expected = annotation_coverage_status(
            reported=counts["reported_distinct_taxa"],
            matched=counts["lookup_matched_taxa"],
            unmatched=counts["lookup_unmatched_taxa"],
            with_annotations=counts["matched_taxa_with_annotations"],
            without_annotations=counts["matched_taxa_without_annotations"],
        )
        if status != expected:
            raise ValueError(
                f"Coverage status mismatch in {path}:{row.line_number}: "
                f"recorded {status!r}, expected {expected!r}"
            )
        unmatched_ids = _parse_identifier_list(
            row.values[positions["unmatched_taxon_ids"]],
            path=path,
            line_number=row.line_number,
            column="unmatched_taxon_ids",
        )
        no_annotation_ids = _parse_identifier_list(
            row.values[positions["matched_without_annotation_taxon_ids"]],
            path=path,
            line_number=row.line_number,
            column="matched_without_annotation_taxon_ids",
        )
        if len(unmatched_ids) != counts["lookup_unmatched_taxa"]:
            raise ValueError(f"Unmatched identifier count mismatch in {path}:{row.line_number}")
        if len(no_annotation_ids) != counts["matched_taxa_without_annotations"]:
            raise ValueError(f"Unannotated identifier count mismatch in {path}:{row.line_number}")
        if set(unmatched_ids).intersection(no_annotation_ids):
            raise ValueError(
                f"Matched and unmatched identifiers overlap in {path}:{row.line_number}"
            )
        if row.key[2] not in {"inc", "dec"}:
            raise ValueError(
                f"Invalid profile direction in {path}:{row.line_number}: {row.key[2]!r}"
            )
        trait_filter_id = row.values[positions["trait_filter_id"]].strip()
        if not trait_filter_id:
            raise ValueError(f"Missing trait_filter_id in {path}:{row.line_number}")
        lineage_resolution = row.values[positions["lineage_resolution"]]
        if lineage_resolution != "unresolved":
            raise ValueError(
                f"Aggregate matrix coverage cannot resolve direct/descendant lineage in {path}"
            )
        record = ProfileCoverage(
            study_id=row.key[0],
            experiment_id=row.key[1],
            inc_dec=row.key[2],
            reported_distinct_taxa=counts["reported_distinct_taxa"],
            lookup_matched_taxa=counts["lookup_matched_taxa"],
            lookup_unmatched_taxa=counts["lookup_unmatched_taxa"],
            matched_taxa_with_annotations=counts["matched_taxa_with_annotations"],
            matched_taxa_without_annotations=counts["matched_taxa_without_annotations"],
            coverage_status=status,
            unmatched_taxon_ids=unmatched_ids,
            matched_without_annotation_taxon_ids=no_annotation_ids,
            trait_filter_id=trait_filter_id,
            lineage_resolution=lineage_resolution,
        )
        coverage[record.key] = record
    return coverage
