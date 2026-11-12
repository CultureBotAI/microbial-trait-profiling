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
from mtp.datamodel import core_products_pydantic as _core_products
from mtp.datamodel import mtp_pydantic as _models
from mtp.datamodel import multiple_testing_pydantic as _multiple_testing
from mtp.datamodel import workflow_pydantic as _workflow

# Keep this package's search path and alias children before Python can execute
# generated models under a second name. The legacy generated file remains a
# separate schema artifact; public imports use the canonical model classes.
for _name, _module in (
    ("mtp_pydantic", _models),
    ("bsdbng_pydantic", _models),
    ("core_products_pydantic", _core_products),
    ("multiple_testing_pydantic", _multiple_testing),
    ("workflow_pydantic", _workflow),
):
    sys.modules[f"{__name__}.{_name}"] = _module
    globals()[_name] = _module


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.datamodel", run_name="__main__")
