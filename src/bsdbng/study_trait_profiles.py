"""Compatibility shim forwarding bsdbng.study_trait_profiles to mtp.study_trait_profiles."""

import sys
from typing import Any

import mtp.study_trait_profiles as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.study_trait_profiles", run_name="__main__")
