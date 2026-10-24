"""Compatibility shim.

Forwards bsdbng.composition_calibration.adapters.jed_composition
to mtp.composition_calibration.adapters.jed_composition.
"""

from __future__ import annotations

import sys
from typing import Any

from mtp.composition_calibration.adapters.jed_composition import (
    DEFAULT_PSEUDOCOUNT_SCALE,
    RESULT_COLUMNS,
    CompositionContrastData,
    JedCompositionAdapter,
    compute_composition_contrasts,
    compute_profile_zscores,
    filter_paired_profiles,
)

__all__ = [
    "DEFAULT_PSEUDOCOUNT_SCALE",
    "RESULT_COLUMNS",
    "CompositionContrastData",
    "JedCompositionAdapter",
    "compute_composition_contrasts",
    "compute_profile_zscores",
    "filter_paired_profiles",
]

import mtp.composition_calibration.adapters.jed_composition as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.composition_calibration.adapters.jed_composition", run_name="__main__")
