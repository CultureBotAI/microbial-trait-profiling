"""Compatibility shim forwarding bsdbng.ingest to mtp.ingest."""

import sys
from typing import Any

from mtp.ingest import (
    METAPHLAN_RANK_MAP,
    REQUIRED_EXPORTS,
    _build_rows_from_separate_csvs,
    _build_study_record,
    _clean,
    _clean_float,
    _clean_int,
    _clean_list,
    _parse_experiment_id,
    _parse_metaphlan_name,
    _parse_signature_id,
    _parse_study_id,
    _parse_taxa,
    _read_csv,
    _taxon_name_to_rank,
    assert_required_exports,
    ingest,
)

__all__ = [
    "METAPHLAN_RANK_MAP",
    "REQUIRED_EXPORTS",
    "_build_rows_from_separate_csvs",
    "_build_study_record",
    "_clean",
    "_clean_float",
    "_clean_int",
    "_clean_list",
    "_parse_experiment_id",
    "_parse_metaphlan_name",
    "_parse_signature_id",
    "_parse_study_id",
    "_parse_taxa",
    "_read_csv",
    "_taxon_name_to_rank",
    "assert_required_exports",
    "ingest",
]

import mtp.ingest as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.ingest", run_name="__main__")
