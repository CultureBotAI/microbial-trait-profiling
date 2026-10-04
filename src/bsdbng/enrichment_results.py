"""Shared contract for trait enrichment (ORA) result files.

The analyzer (``scripts/trait_enrichment.py``) writes a results TSV of the
selected (significant) tests, preceded by ``# Key: Value`` metadata lines. The
ranking exporter and the enrichment summary read that file. This module is the
single place that defines the required columns, validates the file, and
classifies the result so that an empty table stays interpretable (#296):

* ``discoveries``         -- at least one selected row.
* ``zero_discoveries``    -- tests were performed, none were selected.
* ``zero_eligible_tests`` -- no tests could be performed.

Malformed input (missing file, missing columns, non-numeric statistics,
row counts that disagree with the metadata, or an empty table that does not
say how many tests were performed) raises :class:`MalformedResultsError`
instead of yielding an empty successful report.
"""

from __future__ import annotations

import io
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

STATUS_DISCOVERIES = "discoveries"
STATUS_ZERO_DISCOVERIES = "zero_discoveries"
STATUS_ZERO_ELIGIBLE_TESTS = "zero_eligible_tests"
RESULT_STATUSES: tuple[str, ...] = (
    STATUS_DISCOVERIES,
    STATUS_ZERO_DISCOVERIES,
    STATUS_ZERO_ELIGIBLE_TESTS,
)

# Metadata keys written by the analyzer and read downstream.
META_TESTS_PERFORMED = "Total-Tests-Performed"
META_SIGNIFICANT_HITS = "Significant-Hits"
META_RESULT_STATUS = "Result-Status"

ID_COLUMNS: tuple[str, ...] = ("study_id", "experiment_id", "inc_dec", "trait")
PROBABILITY_COLUMNS: tuple[str, ...] = ("p_value", "fdr_q")
FLOAT_COLUMNS: tuple[str, ...] = ("fold", "expected")
COUNT_COLUMNS: tuple[str, ...] = ("count", "testable_profiles")
RESULT_COLUMNS: tuple[str, ...] = (
    "study_id",
    "experiment_id",
    "inc_dec",
    "trait",
    "p_value",
    "fdr_q",
    "fold",
    "count",
    "expected",
    "testable_profiles",
)

_REGENERATE_HINT = (
    "Regenerate it with `scripts/trait_enrichment.py --tsv-output <path>` (or `just analyze`)."
)


class MalformedResultsError(ValueError):
    """Raised when an enrichment results file cannot be interpreted."""


def classify_result_status(tests_performed: int, selected_rows: int) -> str:
    """Classify a result set from its test count and its number of selected rows."""
    if tests_performed < 0 or selected_rows < 0:
        raise MalformedResultsError(
            f"Counts must be non-negative (tests_performed={tests_performed}, "
            f"selected_rows={selected_rows})."
        )
    if selected_rows > tests_performed:
        raise MalformedResultsError(
            f"{selected_rows} selected rows exceed the {tests_performed} tests performed."
        )
    if tests_performed == 0:
        return STATUS_ZERO_ELIGIBLE_TESTS
    if selected_rows == 0:
        return STATUS_ZERO_DISCOVERIES
    return STATUS_DISCOVERIES


@dataclass(frozen=True)
class ResultsTable:
    """A validated enrichment results TSV."""

    path: Path
    frame: Any  # pandas.DataFrame with RESULT_COLUMNS (ids as str, statistics numeric)
    metadata: dict[str, str]
    tests_performed: int | None
    significant_hits: int | None
    status: str


def parse_metadata_lines(lines: list[str]) -> dict[str, str]:
    """Parse ``# Key: Value`` lines into a dict (other comment lines are ignored)."""
    metadata: dict[str, str] = {}
    for line in lines:
        content = line.lstrip("#").strip()
        key, sep, value = content.partition(": ")
        if sep and key:
            metadata[key.strip()] = value.strip()
    return metadata


def _parse_metadata_count(metadata: dict[str, str], key: str, path: Path) -> int | None:
    raw = metadata.get(key)
    if raw is None:
        return None
    try:
        value = int(raw)
    except ValueError:
        raise MalformedResultsError(
            f"{path}: metadata '# {key}: {raw}' is not an integer. {_REGENERATE_HINT}"
        ) from None
    if value < 0:
        raise MalformedResultsError(f"{path}: metadata '# {key}: {raw}' is negative.")
    return value


def _validate_numeric(frame: Any, path: Path, first_data_line: int) -> Any:
    for column in PROBABILITY_COLUMNS + FLOAT_COLUMNS + COUNT_COLUMNS:
        converted = pd.to_numeric(frame[column], errors="coerce")
        bad = converted.isna()
        if column in PROBABILITY_COLUMNS:
            bad = bad | (converted < 0) | (converted > 1)
        if column in COUNT_COLUMNS:
            bad = bad | (converted < 0) | (converted.map(lambda v: not float(v).is_integer()))
        if bool(bad.any()):
            position = int(bad.to_numpy().nonzero()[0][0])
            raise MalformedResultsError(
                f"{path}, line {first_data_line + position}: invalid value "
                f"{frame[column].iloc[position]!r} in column '{column}'. {_REGENERATE_HINT}"
            )
        frame[column] = converted.astype("int64") if column in COUNT_COLUMNS else converted
    for column in ID_COLUMNS:
        empty = frame[column].map(lambda v: not str(v).strip())
        if bool(empty.any()):
            position = int(empty.to_numpy().nonzero()[0][0])
            raise MalformedResultsError(
                f"{path}, line {first_data_line + position}: empty value in column '{column}'."
            )
    return frame


def load_results_tsv(path: Path) -> ResultsTable:
    """Read and validate an analyzer results TSV.

    Raises ``FileNotFoundError`` if the file is missing and
    :class:`MalformedResultsError` if it cannot be interpreted.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Enrichment results TSV not found: {path}. {_REGENERATE_HINT}")

    lines = path.read_text(encoding="utf-8").splitlines()
    header_index = 0
    while header_index < len(lines) and lines[header_index].startswith("#"):
        header_index += 1
    metadata = parse_metadata_lines(lines[:header_index])

    body = lines[header_index:]
    if not body or not body[0].strip():
        raise MalformedResultsError(f"{path}: no column header line found. {_REGENERATE_HINT}")
    header = body[0].split("\t")
    missing = [column for column in RESULT_COLUMNS if column not in header]
    if missing:
        raise MalformedResultsError(
            f"{path}: missing required columns {missing} (found {header}). {_REGENERATE_HINT}"
        )

    frame = pd.read_csv(
        io.StringIO("\n".join(body)),
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    # Line numbers are 1-based: metadata lines, then the header, then data.
    frame = _validate_numeric(frame, path, first_data_line=header_index + 2)
    for column in ID_COLUMNS:
        frame[column] = frame[column].astype(str)

    n_rows = len(frame)
    tests_performed = _parse_metadata_count(metadata, META_TESTS_PERFORMED, path)
    significant_hits = _parse_metadata_count(metadata, META_SIGNIFICANT_HITS, path)

    if significant_hits is not None and significant_hits != n_rows:
        raise MalformedResultsError(
            f"{path}: metadata declares {significant_hits} significant hits but the file "
            f"contains {n_rows} rows; the file may be truncated. {_REGENERATE_HINT}"
        )

    if tests_performed is None:
        if n_rows == 0:
            raise MalformedResultsError(
                f"{path}: contains zero result rows and no '# {META_TESTS_PERFORMED}' metadata, "
                "so zero discoveries cannot be distinguished from zero eligible tests or a "
                f"failed run. {_REGENERATE_HINT}"
            )
        status = STATUS_DISCOVERIES
    else:
        try:
            status = classify_result_status(tests_performed, n_rows)
        except MalformedResultsError as exc:
            raise MalformedResultsError(f"{path}: {exc} {_REGENERATE_HINT}") from None

    declared = metadata.get(META_RESULT_STATUS)
    if declared is not None and declared != status:
        raise MalformedResultsError(
            f"{path}: metadata declares status '{declared}' but the rows and test count "
            f"imply '{status}'. {_REGENERATE_HINT}"
        )

    return ResultsTable(
        path=path,
        frame=frame,
        metadata=metadata,
        tests_performed=tests_performed,
        significant_hits=significant_hits,
        status=status,
    )
