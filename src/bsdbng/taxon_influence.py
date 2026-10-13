"""Compatibility shim forwarding bsdbng.taxon_influence to mtp.taxon_influence."""

import sys
from typing import Any

from mtp.taxon_influence import (
    TIE,
    adjust,
    comparison,
    concentration,
    critical_t,
    describe,
    directional_counts,
    directional_test,
    effect_rows,
    overlap,
    pair_manifest,
    pair_values,
    removal,
    run_influence,
    scenarios,
    scope_keys,
    source_occurrences,
    study_means,
    taxonomic_counts,
    trait_leading_genera,
)

__all__ = [
    "TIE",
    "adjust",
    "comparison",
    "concentration",
    "critical_t",
    "describe",
    "directional_counts",
    "directional_test",
    "effect_rows",
    "overlap",
    "pair_manifest",
    "pair_values",
    "removal",
    "run_influence",
    "scenarios",
    "scope_keys",
    "source_occurrences",
    "study_means",
    "taxonomic_counts",
    "trait_leading_genera",
]

import mtp.taxon_influence as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.taxon_influence", run_name="__main__")
