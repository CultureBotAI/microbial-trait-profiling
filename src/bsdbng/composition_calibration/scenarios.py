"""Compatibility shim.

Forwards bsdbng.composition_calibration.scenarios
to mtp.composition_calibration.scenarios.
"""

from __future__ import annotations

import sys
from typing import Any

from mtp.composition_calibration.scenarios import (
    generate_study_flip_null_profiles,
    replicate_rng,
    stream_key,
)

__all__ = [
    "generate_study_flip_null_profiles",
    "replicate_rng",
    "stream_key",
]

import mtp.composition_calibration.scenarios as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.composition_calibration.scenarios", run_name="__main__")
