"""Compatibility shim forwarding bsdbng.label_support_audit to mtp.label_support_audit."""

import sys
from typing import Any

from mtp.label_support_audit import (
    AUDIT_TRAITS,
    COASSERTION_PAIRS,
    FOCAL_TAXA,
    FOCAL_TRAITS,
    MADIN_RULES,
    MADIN_SOURCE,
    METADATA_ONLY_MADIN_SOURCES,
    PROFILE_KEYS,
    Row,
    Scope,
    bind_occurrences,
    build_coassertion_census,
    build_label_support_tables,
    build_support_threshold_sensitivity,
    load_kg_provenance_graph,
    load_madin_subsources,
    main,
    metadata_source_status,
    run_label_support_audit,
)

__all__ = [
    "AUDIT_TRAITS",
    "COASSERTION_PAIRS",
    "FOCAL_TAXA",
    "FOCAL_TRAITS",
    "MADIN_RULES",
    "MADIN_SOURCE",
    "METADATA_ONLY_MADIN_SOURCES",
    "PROFILE_KEYS",
    "Row",
    "Scope",
    "bind_occurrences",
    "build_coassertion_census",
    "build_label_support_tables",
    "build_support_threshold_sensitivity",
    "load_kg_provenance_graph",
    "load_madin_subsources",
    "main",
    "metadata_source_status",
    "run_label_support_audit",
]

import mtp.label_support_audit as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.label_support_audit", run_name="__main__")
