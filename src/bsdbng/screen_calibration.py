"""Compatibility shim forwarding bsdbng.screen_calibration to mtp.screen_calibration."""

import sys
from typing import Any

from mtp.screen_calibration import (
    CORRECTIONS,
    DEFAULT_OUTPUT_DIR,
    STREAM,
    THRESHOLDS,
    UNIT,
    VERSION,
    ReplicateCount,
    count_calls,
    main,
    replicate_flipped_profiles,
    run_null_screens,
    summarize,
)

__all__ = [
    "CORRECTIONS",
    "DEFAULT_OUTPUT_DIR",
    "STREAM",
    "THRESHOLDS",
    "UNIT",
    "VERSION",
    "ReplicateCount",
    "count_calls",
    "main",
    "replicate_flipped_profiles",
    "run_null_screens",
    "summarize",
]

import mtp.screen_calibration as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.screen_calibration", run_name="__main__")
