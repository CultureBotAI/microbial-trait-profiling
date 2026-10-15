"""Compatibility shim forwarding bsdbng.study_tsv to mtp.study_tsv."""

import sys
from typing import Any

import mtp.study_tsv as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.study_tsv", run_name="__main__")
