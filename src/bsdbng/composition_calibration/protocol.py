"""Compatibility shim.

Forwards bsdbng.composition_calibration.protocol
to mtp.composition_calibration.protocol.
"""

from __future__ import annotations

import sys
from typing import Any

from mtp.composition_calibration.protocol import (
    MethodAdapter,
    MethodInferenceResult,
)

__all__ = [
    "MethodAdapter",
    "MethodInferenceResult",
]

import mtp.composition_calibration.protocol as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.composition_calibration.protocol", run_name="__main__")
