"""Compatibility shim forwarding bsdbng.enrichment_results to mtp.enrichment_results."""

import sys
from typing import Any

from mtp.enrichment_results import (
    _REGENERATE_HINT,
    COUNT_COLUMNS,
    FLOAT_COLUMNS,
    ID_COLUMNS,
    META_RESULT_STATUS,
    META_SIGNIFICANT_HITS,
    META_TESTS_PERFORMED,
    PROBABILITY_COLUMNS,
    RESULT_COLUMNS,
    RESULT_STATUSES,
    STATUS_DISCOVERIES,
    STATUS_ZERO_DISCOVERIES,
    STATUS_ZERO_ELIGIBLE_TESTS,
    MalformedResultsError,
    ResultsTable,
    _parse_metadata_count,
    _validate_numeric,
    classify_result_status,
    load_results_tsv,
    parse_metadata_lines,
)

__all__ = [
    "COUNT_COLUMNS",
    "FLOAT_COLUMNS",
    "ID_COLUMNS",
    "META_RESULT_STATUS",
    "META_SIGNIFICANT_HITS",
    "META_TESTS_PERFORMED",
    "PROBABILITY_COLUMNS",
    "RESULT_COLUMNS",
    "RESULT_STATUSES",
    "STATUS_DISCOVERIES",
    "STATUS_ZERO_DISCOVERIES",
    "STATUS_ZERO_ELIGIBLE_TESTS",
    "_REGENERATE_HINT",
    "MalformedResultsError",
    "ResultsTable",
    "_parse_metadata_count",
    "_validate_numeric",
    "classify_result_status",
    "load_results_tsv",
    "parse_metadata_lines",
]

import mtp.enrichment_results as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.enrichment_results", run_name="__main__")
