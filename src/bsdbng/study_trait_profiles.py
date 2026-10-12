"""Compatibility shim forwarding bsdbng.study_trait_profiles to mtp.study_trait_profiles."""

import sys
from typing import Any

from mtp.study_trait_profiles import (
    DEFAULT_ORA_OUTPUT_PATH,
    DEFAULT_OUTPUT_PATH,
    DEFAULT_TAXA_TRAIT_MATRIX_PATH,
    PROFILE_FORMATS,
    PROFILE_ID_SEPARATOR,
    ProfileFormat,
    ProfileKey,
    StudyTraitProfileStats,
    TaxaTraitMatrix,
    _as_text,
    _default_output_path,
    _load_yaml,
    _parse_count,
    aggregate_profile_counts,
    build_parser,
    collect_profile_taxa,
    create_study_trait_profiles,
    direction_to_inc_dec,
    main,
    profile_id,
    read_taxa_trait_matrix,
    write_long_profiles,
    write_wide_profiles,
)

__all__ = [
    "DEFAULT_ORA_OUTPUT_PATH",
    "DEFAULT_OUTPUT_PATH",
    "DEFAULT_TAXA_TRAIT_MATRIX_PATH",
    "PROFILE_FORMATS",
    "PROFILE_ID_SEPARATOR",
    "ProfileFormat",
    "ProfileKey",
    "StudyTraitProfileStats",
    "TaxaTraitMatrix",
    "_as_text",
    "_default_output_path",
    "_load_yaml",
    "_parse_count",
    "aggregate_profile_counts",
    "build_parser",
    "collect_profile_taxa",
    "create_study_trait_profiles",
    "direction_to_inc_dec",
    "main",
    "profile_id",
    "read_taxa_trait_matrix",
    "write_long_profiles",
    "write_wide_profiles",
]

import mtp.study_trait_profiles as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.study_trait_profiles", run_name="__main__")
