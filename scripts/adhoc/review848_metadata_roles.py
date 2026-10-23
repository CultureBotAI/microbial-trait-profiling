"""Review #848 predicate-role guards; synthetic TSV inputs, JSON evidence output (2026-10-22).

Run with the reviewed worktree's locked environment/PYTHONPATH. This helper only
creates synthetic input files and exercises public readers; no production KG is read.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from mtp.taxa_trait_matrix import read_kg_edges
from mtp.trait_terms import unlisted_negations


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.output_dir
    root.mkdir(parents=True, exist_ok=True)
    nodes = root / "nodes.tsv"
    with nodes.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "category", "name", "synonym"])
        writer.writerow(
            ["METPO:1999999", "biolink:PhenotypicQuality", "lacks flagella", "does not swim"]
        )
        writer.writerow(
            ["example:negative_relation", "biolink:OntologyClass", "does not exhibit", ""]
        )
        writer.writerow(["METPO:1000602", "biolink:PhenotypicQuality", "aerobic", ""])
    edges = root / "edges.tsv"
    edges.write_text(
        "id\tsubject\tpredicate\tobject\n"
        "phenotype\tNCBITaxon:1\tbiolink:has_phenotype\tMETPO:1999999\n"
        "negative\tNCBITaxon:1\texample:negative_relation\tMETPO:1000602\n",
        encoding="utf-8",
    )
    flagged = unlisted_negations(nodes)
    try:
        _, full = read_kg_edges(edges, strain_rollup="full", nodes_path=nodes)
        outcome = {node: sorted(keys) for node, keys in full.items()}
        error = None
    except ValueError as exception:
        outcome = {}
        error = str(exception)
    record = {
        "phenotype_only_label_flagged_as_predicate": "METPO:1999999" in flagged,
        "standalone_unlisted_negations": flagged,
        "full_scope_traits": outcome,
        "full_scope_error": error,
    }
    (root / "evidence.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
