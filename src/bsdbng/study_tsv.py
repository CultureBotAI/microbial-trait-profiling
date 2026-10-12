"""Compatibility shim forwarding bsdbng.study_tsv to mtp.study_tsv."""

import sys
from typing import Any

from mtp.study_tsv import (
    GeneratedTables,
    _as_text,
    _load_yaml,
    _ontology_terms,
    _profile_feature,
    _split_profile_direction,
    _write_dict_rows,
    _write_profile_long,
    _write_profile_wide,
    generate_tables,
    main,
)

__all__ = [
    "GeneratedTables",
    "_as_text",
    "_load_yaml",
    "_ontology_terms",
    "_profile_feature",
    "_split_profile_direction",
    "_write_dict_rows",
    "_write_profile_long",
    "_write_profile_wide",
    "generate_tables",
    "main",
]

import mtp.study_tsv as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.study_tsv", run_name="__main__")
