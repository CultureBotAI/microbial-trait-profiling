"""Compatibility shim forwarding bsdbng.taxa_trait_matrix to mtp.taxa_trait_matrix."""

import sys
from typing import Any

from mtp.taxa_trait_matrix import (
    DEFAULT_KG_DIR,
    DEFAULT_OUTPUT_PATH,
    DEFAULT_TAXA_PATH,
    EDGES_FILENAME,
    NODES_FILENAME,
    SUBCLASS_OF,
    TAXON_PREFIX,
    TRAIT_KEY_SEPARATOR,
    TRAIT_PREFIX,
    MatrixStats,
    NodeMetadata,
    ReportedTaxon,
    _clean_label,
    _is_taxon,
    _is_trait,
    _node_label,
    _target_ancestors,
    build_parser,
    column_metadata_rows,
    compute_trait_counts,
    create_taxa_trait_matrix,
    main,
    make_trait_key,
    matrix_rows,
    node_ids_for_trait_keys,
    read_kg_edges,
    read_node_metadata,
    read_reported_taxa,
    resolve_kg_file,
    split_trait_key,
    trait_column_id,
    trait_column_label,
    write_matrix,
)

__all__ = [
    "DEFAULT_KG_DIR",
    "DEFAULT_OUTPUT_PATH",
    "DEFAULT_TAXA_PATH",
    "EDGES_FILENAME",
    "NODES_FILENAME",
    "SUBCLASS_OF",
    "TAXON_PREFIX",
    "TRAIT_KEY_SEPARATOR",
    "TRAIT_PREFIX",
    "MatrixStats",
    "NodeMetadata",
    "ReportedTaxon",
    "_clean_label",
    "_is_taxon",
    "_is_trait",
    "_node_label",
    "_target_ancestors",
    "build_parser",
    "column_metadata_rows",
    "compute_trait_counts",
    "create_taxa_trait_matrix",
    "main",
    "make_trait_key",
    "matrix_rows",
    "node_ids_for_trait_keys",
    "read_kg_edges",
    "read_node_metadata",
    "read_reported_taxa",
    "resolve_kg_file",
    "split_trait_key",
    "trait_column_id",
    "trait_column_label",
    "write_matrix",
]

import mtp.taxa_trait_matrix as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.taxa_trait_matrix", run_name="__main__")
