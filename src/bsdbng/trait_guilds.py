"""Compatibility shim forwarding bsdbng.trait_guilds to mtp.trait_guilds."""

import sys
from typing import Any

from mtp.trait_guilds import (
    DEFAULT_EXPERIMENTS_PATH,
    DEFAULT_FULL_DUMP_PATH,
    DEFAULT_PROFILE_PATH,
    DIRECTIONS,
    GRAM_GUILDS,
    HUMAN_HOST_SPECIES,
    MEDIADIVE_MARKER,
    OXYGEN_GUILDS,
    PROFILE_ID_COLUMNS,
    PROFILE_METADATA_COLUMNS,
    TABLE1_TRAIT_COLORS,
    TABLE1_TRAIT_LABELS,
    TABLE1_TRAITS,
    CorpusStats,
    ExperimentMetadata,
    PairExclusion,
    Profile,
    TraitColumnSelection,
    _apply_pair_filter,
    _read_profile_header,
    _select_trait_columns,
    column_term_id,
    drop_zero_total_pairs,
    format_pair_exclusion,
    format_stats_report,
    load_profiles,
    read_bugsigdb_vintage,
    resolve_trait_columns,
    retained_trait_columns,
    write_stats_report,
)

__all__ = [
    "DEFAULT_EXPERIMENTS_PATH",
    "DEFAULT_FULL_DUMP_PATH",
    "DEFAULT_PROFILE_PATH",
    "DIRECTIONS",
    "GRAM_GUILDS",
    "HUMAN_HOST_SPECIES",
    "MEDIADIVE_MARKER",
    "OXYGEN_GUILDS",
    "PROFILE_ID_COLUMNS",
    "PROFILE_METADATA_COLUMNS",
    "TABLE1_TRAITS",
    "TABLE1_TRAIT_COLORS",
    "TABLE1_TRAIT_LABELS",
    "CorpusStats",
    "ExperimentMetadata",
    "PairExclusion",
    "Profile",
    "TraitColumnSelection",
    "_apply_pair_filter",
    "_read_profile_header",
    "_select_trait_columns",
    "column_term_id",
    "drop_zero_total_pairs",
    "format_pair_exclusion",
    "format_stats_report",
    "load_profiles",
    "read_bugsigdb_vintage",
    "resolve_trait_columns",
    "retained_trait_columns",
    "write_stats_report",
]

import mtp.trait_guilds as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.trait_guilds", run_name="__main__")
