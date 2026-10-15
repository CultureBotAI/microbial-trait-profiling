"""Compatibility shim.

Forwards bsdbng.composition_calibration.composition
to mtp.composition_calibration.composition.
"""

from __future__ import annotations

import sys
from typing import Any

from mtp.composition_calibration.composition import (
    CORRECTIONS,
    DEFAULT_COLUMNS_PATH,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_PROFILES_PATH,
    DEFAULT_REPLICATES,
    SEED,
    STREAM,
    THRESHOLDS,
    UNIT,
    ReplicateCallCount,
    SummaryRow,
    Table1CalibrationRow,
    build_parser,
    load_composition_trait_columns,
    main,
    run_composition_calibration,
    save_composition_calibration_artifacts,
)

__all__ = [
    "CORRECTIONS",
    "DEFAULT_COLUMNS_PATH",
    "DEFAULT_OUTPUT_DIR",
    "DEFAULT_PROFILES_PATH",
    "DEFAULT_REPLICATES",
    "SEED",
    "STREAM",
    "THRESHOLDS",
    "UNIT",
    "ReplicateCallCount",
    "SummaryRow",
    "Table1CalibrationRow",
    "build_parser",
    "load_composition_trait_columns",
    "main",
    "run_composition_calibration",
    "save_composition_calibration_artifacts",
]

import mtp.composition_calibration.composition as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.composition_calibration.composition", run_name="__main__")
