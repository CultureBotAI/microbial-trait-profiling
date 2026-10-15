"""Compatibility shim forwarding bsdbng.taxa_trait_matrix to mtp.taxa_trait_matrix."""

import sys
from typing import Any

import mtp.taxa_trait_matrix as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.taxa_trait_matrix", run_name="__main__")
