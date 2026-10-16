"""Paired composition arm calibration package (#263)."""

from __future__ import annotations

from mtp.composition_calibration.adapters.jed_composition import (
    DEFAULT_PSEUDOCOUNT_SCALE,
    CompositionContrastData,
    JedCompositionAdapter,
    compute_composition_contrasts,
    compute_profile_zscores,
    filter_paired_profiles,
)
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
from mtp.composition_calibration.protocol import MethodAdapter, MethodInferenceResult
from mtp.composition_calibration.scenarios import (
    generate_study_flip_null_profiles,
    replicate_rng,
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
