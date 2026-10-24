"""Compatibility shim forwarding bsdbng.trait_screen to mtp.trait_screen."""

import sys
from typing import Any

from mtp.trait_screen import (
    ALL_TIED,
    CONFIDENCE_LEVEL,
    CONTRAST_BASIS,
    DEFAULT_COLUMN_LABELS_PATH,
    DEFAULT_OUTPUT_DIR,
    EFFECTIVE_ESTIMATOR,
    MANIFEST_PATH,
    MIN_DEPENDENCE_TRAITS,
    MIN_DEPENDENCE_UNITS,
    NO_ELIGIBLE_STUDIES,
    PROFILE_BASIS,
    Q_THRESHOLDS,
    TESTED,
    VERSION,
    ColumnLabel,
    Dependence,
    ScreenResult,
    TraitScreenRow,
    _cell,
    _dependence_paragraphs,
    _family_rows,
    _report,
    _sha256,
    candidate_slots,
    cheverud_note,
    dependence,
    effective_family_size,
    harmonic,
    load_column_labels,
    main,
    read_header,
    run_screen,
    screen_candidates,
    significant,
    step_up,
    step_up_array,
    table1_slug,
    write_screen_outputs,
)

__all__ = [
    "ALL_TIED",
    "CONFIDENCE_LEVEL",
    "CONTRAST_BASIS",
    "DEFAULT_COLUMN_LABELS_PATH",
    "DEFAULT_OUTPUT_DIR",
    "EFFECTIVE_ESTIMATOR",
    "MANIFEST_PATH",
    "MIN_DEPENDENCE_TRAITS",
    "MIN_DEPENDENCE_UNITS",
    "NO_ELIGIBLE_STUDIES",
    "PROFILE_BASIS",
    "Q_THRESHOLDS",
    "TESTED",
    "VERSION",
    "ColumnLabel",
    "Dependence",
    "ScreenResult",
    "TraitScreenRow",
    "_cell",
    "_dependence_paragraphs",
    "_family_rows",
    "_report",
    "_sha256",
    "candidate_slots",
    "cheverud_note",
    "dependence",
    "effective_family_size",
    "harmonic",
    "load_column_labels",
    "main",
    "read_header",
    "run_screen",
    "screen_candidates",
    "significant",
    "step_up",
    "step_up_array",
    "table1_slug",
    "write_screen_outputs",
]

import mtp.trait_screen as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.trait_screen", run_name="__main__")
