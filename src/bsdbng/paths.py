"""Compatibility shim forwarding bsdbng.paths to mtp.paths."""

import sys
from typing import Any

from mtp.paths import (
    DERIVED_DATA_DIR,
    KG_MICROBE_DIR,
    RAW_DATA_DIR,
    REPO_ROOT,
    SCHEMA_PATH,
    STUDY_DATA_DIR,
)

__all__ = [
    "DERIVED_DATA_DIR",
    "KG_MICROBE_DIR",
    "RAW_DATA_DIR",
    "REPO_ROOT",
    "SCHEMA_PATH",
    "STUDY_DATA_DIR",
]

import mtp.paths as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.paths", run_name="__main__")
