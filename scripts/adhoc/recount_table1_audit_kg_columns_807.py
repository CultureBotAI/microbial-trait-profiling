"""Recount the Table 1 audit's KG columns and scope counts for #807 (2026-10-23).

Purpose: recompute the numbers in docs/table1_trait_semantics_audit.md that come from the
pinned KG-Microbe edges (#807, #1004):
- the edges on each of the six Table 1 trait objects, by subject kind, predicate and
  primary_knowledge_source;
- the NCBITaxon nodes that carry each trait directly (off) and under the phenotype
  strain rollup, and the BacDive strains that carry it;
- the direct coccus/sphere overlap;
- the reported BugSigDB taxa that bear each trait under off and under phenotype. This
  is the matrix's own rollup, so it reproduces the matrix counts of each scope.
The node and taxon counts use mtp.taxa_trait_matrix's own reader and scope rule.

Inputs:
- --edges: default data/kg-microbe/merged-kg_edges.tsv. The expected md5 555ecf5b is the
  one data/kg-microbe/provenance.json records.
- --reported-taxa: default data/derived/bugsigdb_taxa.tsv. The expected sha256 prefix
  e7fa982b is the pin in data/analysis/cellulose_decr/cellulose_decr_provenance.json.

Output: plain-text counts on stdout. The edge file is read twice (about 1 GB each time),
so this is a heavy run.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from mtp.taxa_trait_matrix import (
    compute_trait_counts,
    read_kg_trait_edges,
    read_reported_taxa,
    scoped_taxon_traits,
)

SIX = {
    "METPO:1000608": "facultatively aerobic",
    "METPO:2000103|GO:0030245": "cellulose catabolic process",
    "METPO:1000668": "coccus shaped",
    "METPO:1000683": "sphere shaped",
    "METPO:2000103|pathways:arsenate_reduction": "arsenate reduction",
    "METPO:1000623": "moderately halophilic",
}
COCCUS, SPHERE = "METPO:1000668", "METPO:1000683"
# Positions in merged-kg_edges.tsv. Its header line holds a stray carriage return, so
# lines are split on tabs by position rather than read by header name.
SUBJECT, PREDICATE, OBJECT, PRIMARY_SOURCE = 1, 2, 3, 10


def edge_counts(edges_path: Path) -> Counter[tuple[str, str, str, str]]:
    objects = {key.split("|")[-1]: key for key in SIX}
    counts: Counter[tuple[str, str, str, str]] = Counter()
    with edges_path.open(encoding="utf-8", newline="\n") as handle:
        next(handle)
        for line in handle:
            row = line.rstrip("\n").split("\t")
            key = objects.get(row[OBJECT])
            if key is None:
                continue
            subject = row[SUBJECT]
            if subject.startswith("NCBITaxon:"):
                kind = "NCBITaxon"
            elif subject.startswith("kgmicrobe.strain:"):
                kind = "strain"
            else:
                kind = subject.split(":")[0]
            source = row[PRIMARY_SOURCE] if len(row) > PRIMARY_SOURCE else ""
            counts[(key, kind, row[PREDICATE], source)] += 1
    return counts


def nodes_with(sets: dict[str, set[str]], key: str) -> set[str]:
    return {node for node, keys in sets.items() if key in keys}


def recount(edges_path: Path, reported_taxa_path: Path) -> None:
    for (key, kind, predicate, source), n in sorted(edge_counts(edges_path).items()):
        print(f"  {SIX[key]:28s} {kind:10s} {predicate:28s} {source:24s} {n}")

    edges = read_kg_trait_edges(edges_path, strains=True)
    direct = edges.taxon_traits
    phenotype = scoped_taxon_traits(
        edges.taxon_traits, edges.strain_taxa, edges.strain_traits, "phenotype"
    )
    for key, label in SIX.items():
        own, rolled = nodes_with(direct, key), nodes_with(phenotype, key)
        strains = [s for s, keys in edges.strain_traits.items() if key in keys]
        linked = [s for s in strains if edges.strain_taxa.get(s)]
        print(
            f"{label:28s} direct NCBITaxon nodes={len(own)} phenotype nodes={len(rolled)} "
            f"added={len(rolled - own)} strains={len(strains)} linked={len(linked)}"
        )
    coccus, sphere = nodes_with(direct, COCCUS), nodes_with(direct, SPHERE)
    shared = len(coccus & sphere)
    print(
        f"direct coccus={len(coccus)} sphere={len(sphere)} intersection={shared} "
        f"sphere-in-coccus={100 * shared / len(sphere):.1f}% "
        f"coccus-in-sphere={100 * shared / len(coccus):.1f}% "
        f"Jaccard={shared / len(coccus | sphere):.4f}"
    )

    reported = read_reported_taxa(reported_taxa_path)
    print(f"reported taxa: {len(reported)}")
    for scope, sets in (("off", direct), ("phenotype", phenotype)):
        counts = compute_trait_counts(reported, edges.child_to_parents, sets)
        bearing = [sum(1 for c in counts.values() if c[key] > 0) for key in SIX]
        print(f"reported taxa bearing each trait, {scope}: {bearing}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edges", type=Path, default=Path("data/kg-microbe/merged-kg_edges.tsv"))
    parser.add_argument(
        "--reported-taxa", type=Path, default=Path("data/derived/bugsigdb_taxa.tsv")
    )
    args = parser.parse_args()
    recount(args.edges, args.reported_taxa)


if __name__ == "__main__":
    main()
