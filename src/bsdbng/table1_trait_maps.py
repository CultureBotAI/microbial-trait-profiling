"""Compatibility shim forwarding bsdbng.table1_trait_maps to mtp.table1_trait_maps."""

import sys
from typing import Any

from mtp.table1_trait_maps import (
    CONFIDENCE_LEVEL,
    FOCUS_PANELS,
    OUTPUT_PREFIX,
    TRAIT_ORDER,
    FocusPanel,
    MissingConditions,
    ShiftResult,
    TraitShift,
    _absent_site_reason,
    build_parser,
    compute_direction_shifts,
    main,
    missing_condition_note,
    missing_conditions,
    study_trait_contrasts,
    tally_direction_signs,
    trait_colorscale,
    write_shift_figure,
    write_shift_table,
)

__all__ = [
    "CONFIDENCE_LEVEL",
    "FOCUS_PANELS",
    "OUTPUT_PREFIX",
    "TRAIT_ORDER",
    "FocusPanel",
    "MissingConditions",
    "ShiftResult",
    "TraitShift",
    "_absent_site_reason",
    "build_parser",
    "compute_direction_shifts",
    "main",
    "missing_condition_note",
    "missing_conditions",
    "study_trait_contrasts",
    "tally_direction_signs",
    "trait_colorscale",
    "write_shift_figure",
    "write_shift_table",
]

import mtp.table1_trait_maps as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.table1_trait_maps", run_name="__main__")
