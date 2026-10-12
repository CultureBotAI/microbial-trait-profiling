"""Compatibility shim forwarding bsdbng.datamodel to mtp.datamodel."""

import sys
from typing import Any

from mtp.datamodel import (
    BiolinkCategoryEnum,
    BiolinkPredicateEnum,
    DirectionEnum,
    Edge,
    Experiment,
    KgTraitAssertion,
    Node,
    Phenotype,
    Signature,
    Study,
    Taxon,
    TaxonomicRankEnum,
    TaxonTraitMatrixRow,
)

__all__ = [
    "BiolinkCategoryEnum",
    "BiolinkPredicateEnum",
    "DirectionEnum",
    "Edge",
    "Experiment",
    "KgTraitAssertion",
    "Node",
    "Phenotype",
    "Signature",
    "Study",
    "Taxon",
    "TaxonTraitMatrixRow",
    "TaxonomicRankEnum",
]

import mtp.datamodel as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.datamodel", run_name="__main__")
