"""Reconcile published taxon matrices to consumed KG evidence for issue #652.

Validation traverses the canonical edge records independently of the matrix
writer. Each annotated source taxon contributes once per feature to itself and
its reported ancestors; duplicate edges and multiple taxonomy paths add nothing.
Under a strain rollup scope, the BacDive strain records it admits count at their
linked NCBITaxon nodes by the matrix builder's own rule, scoped_taxon_traits (#788).
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

from mtp.count_input import taxa_trait_count_reader
from mtp.taxa_trait_matrix import (
    STRAIN_PREFIX,
    SUBCLASS_OF,
    ReportedTaxon,
    StrainRollup,
    assertion_trait_key,
    column_metadata_rows,
    node_ids_for_trait_keys,
    read_assertion_term_metadata,
    read_node_metadata,
    read_reported_taxa,
    scoped_taxon_traits,
    strain_trait_rolls_up,
)
from mtp.workflow.kg_inputs import iter_kg_edges


def _source_counts(
    edges_path: Path, nodes_path: Path, reported: list[ReportedTaxon], strain_rollup: StrainRollup
) -> dict[str, Counter[str]]:
    parents: defaultdict[str, set[str]] = defaultdict(set)
    annotations: defaultdict[str, set[str]] = defaultdict(set)
    strain_taxa: defaultdict[str, set[str]] = defaultdict(set)
    strain_traits: defaultdict[str, set[str]] = defaultdict(set)
    terms = read_assertion_term_metadata(nodes_path, edges_path)
    for record in iter_kg_edges(edges_path):
        subject, predicate, obj = (
            record[column].strip() for column in ("subject", "predicate", "object")
        )
        if subject.startswith("NCBITaxon:"):
            if predicate == SUBCLASS_OF and obj.startswith("NCBITaxon:"):
                parents[subject].add(obj)
            elif obj.startswith("METPO:") or predicate.startswith("METPO:"):
                annotations[subject].add(assertion_trait_key(record, terms))
        elif subject.startswith(STRAIN_PREFIX):
            if predicate == SUBCLASS_OF and obj.startswith("NCBITaxon:"):
                strain_taxa[subject].add(obj)
            elif obj.startswith("METPO:") or predicate.startswith("METPO:"):
                key = assertion_trait_key(record, terms)
                if strain_trait_rolls_up(key, strain_rollup):
                    strain_traits[subject].add(key)
    scoped = scoped_taxon_traits(
        dict(annotations), dict(strain_taxa), dict(strain_traits), strain_rollup
    )
    counts = {taxon.taxon_id: Counter[str]() for taxon in reported}
    for source, features in scoped.items():
        pending = [source]
        visited: set[str] = set()
        while pending:
            current = pending.pop()
            if current in visited:
                continue
            visited.add(current)
            if current in counts:
                counts[current].update(features)
            pending.extend(parents.get(current, ()))
    return counts


def _expected_metadata(
    counts: dict[str, Counter[str]], nodes_path: Path
) -> tuple[list[str], list[list[str]]]:
    features = {feature for values in counts.values() for feature, count in values.items() if count}
    nodes = read_node_metadata(nodes_path, node_ids_for_trait_keys(features))
    return column_metadata_rows(counts, nodes)


def validate_matrix_sources(
    path: Path,
    reported_path: Path,
    edges_path: Path,
    nodes_path: Path,
    metadata_path: Path | None,
    *,
    strain_rollup: StrainRollup = "off",
) -> dict[str, object]:
    """Check all matrix cells, biological row identities and KG feature metadata.

    The expected cells apply the strain rollup scope the producing run declared (#788).
    """
    reported = read_reported_taxa(reported_path)
    reported_by_id = {taxon.taxon_id: taxon for taxon in reported}
    counts = _source_counts(edges_path, nodes_path, reported, strain_rollup)
    metadata_header, metadata_rows = _expected_metadata(counts, nodes_path)
    features = tuple(row[0] for row in metadata_rows)
    reader = taxa_trait_count_reader(path)
    if reader.count_columns != features:
        raise ValueError("Matrix feature vocabulary/order differs from consumed KG evidence")
    seen: set[str] = set()
    for row in reader:
        taxon_id = row.key[0]
        if taxon_id not in reported_by_id:
            raise ValueError(f"Matrix contains an unknown reported taxon: {taxon_id}")
        taxon = reported_by_id[taxon_id]
        if row.values[:3] != (taxon.taxon_id, taxon.taxon_name, taxon.taxonomic_rank):
            raise ValueError(f"Matrix label/rank differs from reported taxon: {taxon_id}")
        expected = tuple(counts[taxon_id][feature] for feature in features)
        if row.counts != expected:
            raise ValueError(f"Matrix count differs from consumed KG evidence for {taxon_id}")
        seen.add(taxon_id)
    if seen != set(reported_by_id):
        raise ValueError("Matrix membership differs from reported taxon source")
    if metadata_path is not None:
        with metadata_path.open(encoding="utf-8", newline="") as handle:
            observed = list(csv.reader(handle, delimiter="\t"))
        if observed != [metadata_header, *metadata_rows]:
            raise ValueError("Feature metadata differs from consumed KG evidence and nodes")
    return {
        "kg_source_aggregation_checked": True,
        "kg_source_cells_checked": reader.summary.rows * len(features),
        "reported_taxon_identity_checked": True,
        "kg_feature_metadata_checked": metadata_path is not None,
    }


def validate_feature_node_metadata(path: Path, nodes_path: Path) -> None:
    """Check feature descriptions against the supplied node snapshot, without a graph rebuild."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    counts = {"features": Counter({row["column_id"]: 1 for row in rows})}
    header, expected = _expected_metadata(counts, nodes_path)
    observed = [[row[column] for column in header] for row in rows]
    if observed != expected:
        raise ValueError("Feature metadata differs from consumed KG nodes")
