"""Compatibility shim forwarding bsdbng.trait_terms to mtp.trait_terms."""

import sys
from typing import Any

from mtp.trait_terms import (
    _HASH_BLOCK_BYTES,
    DEFAULT_NODES_PATH,
    DEFAULT_PROVENANCE_PATH,
    EXPECTED_LABELS,
    NEGATION_PREFIXES,
    NodeTableIdentity,
    TermCheck,
    check_term_labels,
    declared_node_md5,
    file_md5,
    grouped_identifiers,
    identify_node_table,
    main,
    unlisted_negations,
)

__all__ = [
    "DEFAULT_NODES_PATH",
    "DEFAULT_PROVENANCE_PATH",
    "EXPECTED_LABELS",
    "NEGATION_PREFIXES",
    "_HASH_BLOCK_BYTES",
    "NodeTableIdentity",
    "TermCheck",
    "check_term_labels",
    "declared_node_md5",
    "file_md5",
    "grouped_identifiers",
    "identify_node_table",
    "main",
    "unlisted_negations",
]

import mtp.trait_terms as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.trait_terms", run_name="__main__")
