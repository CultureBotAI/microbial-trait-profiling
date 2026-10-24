"""Compatibility shim forwarding bsdbng.composition_calibration to mtp.composition_calibration."""

from __future__ import annotations

import sys
from typing import Any

from mtp.composition_calibration import (
    CORRECTIONS,
    DEFAULT_COLUMNS_PATH,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_PROFILES_PATH,
    DEFAULT_PSEUDOCOUNT_SCALE,
    DEFAULT_REPLICATES,
    SEED,
    STREAM,
    THRESHOLDS,
    UNIT,
    CompositionContrastData,
    JedCompositionAdapter,
    MethodAdapter,
    MethodInferenceResult,
    ReplicateCallCount,
    SummaryRow,
    Table1CalibrationRow,
    build_parser,
    compute_composition_contrasts,
    compute_profile_zscores,
    filter_paired_profiles,
    generate_study_flip_null_profiles,
    load_composition_trait_columns,
    main,
    replicate_rng,
    run_composition_calibration,
    save_composition_calibration_artifacts,
    stream_key,
    study_flip_draw_order,
)

__all__ = [
    "CORRECTIONS",
    "DEFAULT_COLUMNS_PATH",
    "DEFAULT_OUTPUT_DIR",
    "DEFAULT_PROFILES_PATH",
    "DEFAULT_PSEUDOCOUNT_SCALE",
    "DEFAULT_REPLICATES",
    "SEED",
    "STREAM",
    "THRESHOLDS",
    "UNIT",
    "CompositionContrastData",
    "JedCompositionAdapter",
    "MethodAdapter",
    "MethodInferenceResult",
    "ReplicateCallCount",
    "SummaryRow",
    "Table1CalibrationRow",
    "build_parser",
    "compute_composition_contrasts",
    "compute_profile_zscores",
    "filter_paired_profiles",
    "generate_study_flip_null_profiles",
    "load_composition_trait_columns",
    "main",
    "replicate_rng",
    "run_composition_calibration",
    "save_composition_calibration_artifacts",
    "stream_key",
    "study_flip_draw_order",
]

import mtp.composition_calibration as _target
import mtp.composition_calibration.adapters as _adapters
import mtp.composition_calibration.adapters.jed_composition as _jed_composition
import mtp.composition_calibration.composition as _composition
import mtp.composition_calibration.protocol as _protocol
import mtp.composition_calibration.scenarios as _scenarios

sys.modules[__name__] = _target
sys.modules[f"{__name__}.adapters"] = _adapters
sys.modules[f"{__name__}.adapters.jed_composition"] = _jed_composition
sys.modules[f"{__name__}.composition"] = _composition
sys.modules[f"{__name__}.protocol"] = _protocol
sys.modules[f"{__name__}.scenarios"] = _scenarios


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.composition_calibration", run_name="__main__")
