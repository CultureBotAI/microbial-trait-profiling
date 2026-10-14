"""Compatibility shim forwarding bsdbng.ontology_audit to mtp.ontology_audit."""

import sys
from typing import Any

from mtp.ontology_audit import (
    CANONICAL_REDUNDANCY_CLUSTERS,
    DEFAULT_COMPARISON_OUT,
    DEFAULT_PAIRS_OUT,
    DEFAULT_REDUNDANCY_SUMMARY_OUT,
    RedundancyCluster,
    RedundantPair,
    build_redundant_pairs_table,
    collapse_trait_matrix,
    get_canonical_redundancy_clusters,
    main,
    run_ontology_redundancy_audit,
)

__all__ = [
    "CANONICAL_REDUNDANCY_CLUSTERS",
    "DEFAULT_COMPARISON_OUT",
    "DEFAULT_PAIRS_OUT",
    "DEFAULT_REDUNDANCY_SUMMARY_OUT",
    "RedundancyCluster",
    "RedundantPair",
    "build_redundant_pairs_table",
    "collapse_trait_matrix",
    "get_canonical_redundancy_clusters",
    "main",
    "run_ontology_redundancy_audit",
]

import mtp.ontology_audit as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.ontology_audit", run_name="__main__")
