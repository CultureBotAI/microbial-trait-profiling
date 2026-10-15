"""Compatibility shim.

Forwards bsdbng.composition_calibration.adapters
to mtp.composition_calibration.adapters.
"""

from __future__ import annotations

import sys
from typing import Any

from mtp.composition_calibration.adapters import (
    DEFAULT_PSEUDOCOUNT_SCALE,
    CompositionContrastData,
    JedCompositionAdapter,
    compute_composition_contrasts,
    compute_profile_zscores,
    filter_paired_profiles,
)

__all__ = [
    "DEFAULT_PSEUDOCOUNT_SCALE",
    "CompositionContrastData",
    "JedCompositionAdapter",
    "compute_composition_contrasts",
    "compute_profile_zscores",
    "filter_paired_profiles",
]

import mtp.composition_calibration.adapters as _target
import mtp.composition_calibration.adapters.jed_composition as _jed_composition

sys.modules[__name__] = _target
sys.modules[f"{__name__}.jed_composition"] = _jed_composition


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.composition_calibration.adapters", run_name="__main__")
