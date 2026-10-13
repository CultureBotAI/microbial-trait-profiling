"""Compatibility shim forwarding bsdbng.contribution_inputs to mtp.contribution_inputs."""

import sys
from typing import Any

from mtp.contribution_inputs import (
    RANKS,
    Assertion,
    Base,
    Graph,
    Key,
    Pair,
    count,
    digest,
    encoded,
    load_base,
    load_graph,
    rows,
    verify_kg,
)

__all__ = [
    "RANKS",
    "Assertion",
    "Base",
    "Graph",
    "Key",
    "Pair",
    "count",
    "digest",
    "encoded",
    "load_base",
    "load_graph",
    "rows",
    "verify_kg",
]

import mtp.contribution_inputs as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.contribution_inputs", run_name="__main__")
