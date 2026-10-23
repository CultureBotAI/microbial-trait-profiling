"""Compare trait coverage across the BacDive strain rollup scopes (#774, #788).

Inputs: the reported taxa the matrix is built from, the KG edges and nodes, and
optionally a contrast-selected study directory. Outputs, in --output-dir:
trait_coverage.tsv, one row per trait key that any NCBITaxon node or BacDive strain
carries, and summary.json with the strain links (orphan nodes included), and per scope
the vocabulary, the focal traits and, with --study-dir, the signature taxa that gain a
trait label.

The scopes are taxa_trait_matrix's --strain-rollup values: "off" (no strain records),
"phenotype" (bare METPO phenotypes without handling and clinical classes) and "full"
(every strain edge that does not assert an absence, with every negated key dropped).
The evaluable vocabulary is what retained_trait_columns keeps of the matrix columns:
every trait some reported taxon carries, without MediaDive growth media. All counts are
annotation coverage.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from collections.abc import Iterable
from pathlib import Path

from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES
from mtp.paths import file_sha256
from mtp.study_trait_profiles import collect_profile_taxa
from mtp.taxa_trait_matrix import (
    PHENOTYPE_EXCLUDED_TRAITS,
    STRAIN_ROLLUPS,
    KgTraitEdges,
    NodeMetadata,
    StrainRollup,
    compute_trait_counts,
    is_negated,
    node_ids_for_trait_keys,
    read_kg_trait_edges,
    read_node_metadata,
    read_reported_taxa,
    scoped_taxon_traits,
    split_trait_key,
    strain_trait_rolls_up,
    trait_column_label,
)
from mtp.trait_guilds import HUMAN_HOST_SPECIES, MEDIADIVE_MARKER, retained_trait_columns

# The six Table 1 traits and facultatively anaerobic, the oxygen class named in #774.
FOCAL_TRAITS = {
    "METPO:2000103|GO:0030245": "cellulose catabolic process",
    "METPO:1000608": "facultatively aerobic",
    "METPO:1000623": "moderately halophilic",
    "METPO:1000668": "coccus shaped",
    "METPO:1000683": "sphere shaped",
    "METPO:2000103|pathways:arsenate_reduction": "arsenate reduction",
    "METPO:1000605": "facultatively anaerobic",
}
COVERAGE_COLUMNS = (
    "trait_key",
    "label",
    "predicate_id",
    "predicate_name",
    "negated_predicate",
    "culture_medium",
    "excluded_phenotype",
    "focal_trait",
    "linked_strains",
    "unlinked_strains",
    *(f"ncbitaxon_nodes_{scope}" for scope in STRAIN_ROLLUPS),
    *(f"reported_taxa_{scope}" for scope in STRAIN_ROLLUPS),
    *(f"vocabulary_{scope}" for scope in STRAIN_ROLLUPS),
)

Counts = dict[str, Counter[str]]


def carried(trait_sets: Iterable[set[str]]) -> Counter[str]:
    """How many nodes, strains or taxa carry each trait key."""
    return Counter(key for traits in trait_sets for key in traits)


def present(counts: Counts) -> dict[str, set[str]]:
    """The traits each reported taxon carries at or below itself."""
    return {taxon: {key for key, n in row.items() if n > 0} for taxon, row in counts.items()}


def retained(columns: Iterable[str]) -> set[str]:
    """The evaluable vocabulary, by the rule every analysis arm applies."""
    _, kept, _ = retained_trait_columns(sorted(columns))
    return set(kept)


def predicate_name(trait_key: str, metadata: dict[str, NodeMetadata]) -> str:
    predicate, _ = split_trait_key(trait_key)
    if predicate is None:
        return ""
    return metadata.get(predicate, NodeMetadata("", "")).name


def orphan_links(edges: KgTraitEdges, strains: Iterable[str], targets: set[str]) -> dict[str, int]:
    """Linked NCBITaxon nodes without an NCBITaxon parent, and the strains reaching them.

    Such a node's strains reach a matrix row only when the node is itself reported.
    """
    linked = {strain: edges.strain_taxa.get(strain, set()) for strain in strains}
    nodes = {node for links in linked.values() for node in links}
    orphans = {node for node in nodes if node not in edges.child_to_parents}
    return {
        "linked_nodes": len(nodes),
        "orphan_nodes": len(orphans),
        "orphan_strains": sum(bool(links & orphans) for links in linked.values()),
        "orphan_nodes_reported": len(orphans & targets),
    }


def signature_taxa(
    study_dir: Path,
    host_species: str | None,
    labels_off: dict[str, set[str]],
    labels: dict[str, set[str]],
) -> dict[str, object]:
    """Distinct taxa of the 16S and shotgun profiles, and the trait labels they gain."""
    profiles, _ = collect_profile_taxa(
        study_dir, host_species=host_species, sequencing_types=CORPUS_SEQUENCING_TYPES
    )
    taxa = {taxon for observed in profiles.values() for taxon in observed}
    matched = sorted(taxa.intersection(labels_off))
    gaining = [taxon for taxon in matched if labels[taxon] - labels_off[taxon]]
    return {
        "host_species": host_species or "all",
        "profiles": len(profiles),
        "signature_taxa": len(taxa),
        "matched_to_matrix_rows": len(matched),
        "with_a_retained_trait": sum(bool(labels[taxon]) for taxon in matched),
        "gaining_a_trait_label": len(gaining),
        "gaining_a_first_trait_label": sum(not labels_off[taxon] for taxon in gaining),
        "retained_traits_carried": len({k for taxon in matched for k in labels[taxon]}),
        "focal_traits": {
            key: sum(key in labels[taxon] for taxon in matched) for key in FOCAL_TRAITS
        },
    }


def audit(
    reported_taxa_path: Path,
    edges_path: Path,
    nodes_path: Path,
    study_dir: Path | None,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Return the per-trait coverage rows and the summary, checking every scope's rule."""
    reported = read_reported_taxa(reported_taxa_path)
    targets = {taxon.taxon_id for taxon in reported}
    edges = read_kg_trait_edges(edges_path, strains=True)
    nodes: dict[StrainRollup, dict[str, set[str]]] = {}
    labels: dict[StrainRollup, dict[str, set[str]]] = {}
    vocabularies: dict[StrainRollup, set[str]] = {}
    columns: dict[StrainRollup, set[str]] = {}
    for scope in STRAIN_ROLLUPS:
        nodes[scope] = scoped_taxon_traits(
            edges.taxon_traits, edges.strain_taxa, edges.strain_traits, scope
        )
        labels[scope] = present(
            compute_trait_counts(reported, edges.child_to_parents, nodes[scope])
        )
        columns[scope] = set(carried(labels[scope].values()))
        vocabularies[scope] = retained(columns[scope])
    for scope in ("phenotype", "full"):
        dropped = {
            taxon: {k for k in labels["off"][taxon] - labels[scope][taxon] if not is_negated(k)}
            for taxon in labels["off"]
        }
        if any(dropped.values()):
            raise ValueError(f"The {scope} rollup removed a non-negated label from a taxon")
        if not vocabularies["off"] <= vocabularies[scope]:
            raise ValueError(f"A trait left the vocabulary under the {scope} rollup")
    if any(is_negated(key) for keys in nodes["full"].values() for key in keys):
        raise ValueError("A node carries a negated key under the full rollup")
    if set(PHENOTYPE_EXCLUDED_TRAITS) & columns["phenotype"]:
        raise ValueError("An excluded handling or clinical class entered the phenotype matrix")

    node_counts = {scope: carried(nodes[scope].values()) for scope in STRAIN_ROLLUPS}
    taxa_counts = {scope: carried(labels[scope].values()) for scope in STRAIN_ROLLUPS}
    linked = [s for s in edges.strain_traits if edges.strain_taxa.get(s)]
    unlinked = [s for s in edges.strain_traits if not edges.strain_taxa.get(s)]
    strains_linked = carried(edges.strain_traits[s] for s in linked)
    strains_unlinked = carried(edges.strain_traits[s] for s in unlinked)
    keys = set(node_counts["full"]) | set(node_counts["off"]) | set(strains_linked)
    keys |= set(strains_unlinked)
    metadata = read_node_metadata(nodes_path, node_ids_for_trait_keys(keys))
    rows: list[dict[str, object]] = []
    for key in sorted(keys, key=lambda k: (trait_column_label(k, metadata).casefold(), k)):
        predicate, _ = split_trait_key(key)
        rows.append(
            {
                "trait_key": key,
                "label": trait_column_label(key, metadata),
                "predicate_id": predicate or "",
                "predicate_name": predicate_name(key, metadata),
                "negated_predicate": is_negated(key),
                "culture_medium": MEDIADIVE_MARKER in key.lower(),
                "excluded_phenotype": key in PHENOTYPE_EXCLUDED_TRAITS,
                "focal_trait": FOCAL_TRAITS.get(key, ""),
                "linked_strains": strains_linked[key],
                "unlinked_strains": strains_unlinked[key],
                **{f"ncbitaxon_nodes_{s}": node_counts[s][key] for s in STRAIN_ROLLUPS},
                **{f"reported_taxa_{s}": taxa_counts[s][key] for s in STRAIN_ROLLUPS},
                **{f"vocabulary_{s}": key in vocabularies[s] for s in STRAIN_ROLLUPS},
            }
        )

    link_counts = Counter(len(edges.strain_taxa.get(s, ())) for s in edges.strain_traits)
    scopes: dict[str, object] = {}
    for scope in STRAIN_ROLLUPS:
        entering = sorted(vocabularies[scope] - vocabularies["off"])
        names = Counter(
            predicate_name(key, metadata) or "phenotype (METPO object)" for key in entering
        )
        admitted = [
            s
            for s in edges.strain_traits
            if any(strain_trait_rolls_up(key, scope) for key in edges.strain_traits[s])
        ]
        scopes[scope] = {
            "ncbitaxon_nodes_with_a_trait": len(nodes[scope]),
            "matrix_columns": len(columns[scope]),
            "retained_vocabulary": len(vocabularies[scope]),
            "entering": len(entering),
            "leaving": len(vocabularies["off"] - vocabularies[scope]),
            "entering_by_predicate": dict(names.most_common()),
            "strains_with_admitted_edges": len(admitted),
            "orphan_links": orphan_links(edges, admitted, targets),
            "focal_traits": {
                key: {
                    "label": label,
                    "ncbitaxon_nodes": node_counts[scope][key],
                    "reported_taxa": taxa_counts[scope][key],
                    "in_vocabulary": key in vocabularies[scope],
                }
                for key, label in FOCAL_TRAITS.items()
            },
        }
    summary: dict[str, object] = {
        "inputs": {
            str(path): file_sha256(path) for path in (reported_taxa_path, edges_path, nodes_path)
        },
        "reported_taxa": len(reported),
        "strains": {
            "with_trait_edges": len(edges.strain_traits),
            "by_ncbitaxon_link_count": {str(n): c for n, c in sorted(link_counts.items())},
            "linked_without_trait_edges": len(set(edges.strain_taxa) - set(edges.strain_traits)),
            "strain_trait_pairs": sum(len(t) for t in edges.strain_traits.values()),
            "strain_trait_pairs_without_ncbitaxon_link": sum(
                len(edges.strain_traits[s]) for s in unlinked
            ),
            "orphan_links": orphan_links(edges, linked, targets),
        },
        "excluded_phenotype_traits": PHENOTYPE_EXCLUDED_TRAITS,
        "scopes": scopes,
    }
    if study_dir is not None:
        retained_labels = {
            scope: {t: found & vocabularies[scope] for t, found in labels[scope].items()}
            for scope in STRAIN_ROLLUPS
        }
        summary["signature_taxa"] = {
            "study_dir": str(study_dir),
            "cohorts": [
                {
                    scope: signature_taxa(
                        study_dir, host, retained_labels["off"], retained_labels[scope]
                    )
                    for scope in STRAIN_ROLLUPS
                }
                for host in (None, HUMAN_HOST_SPECIES)
            ],
        }
    return rows, summary


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=COVERAGE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reported-taxa", type=Path, required=True)
    parser.add_argument("--edges", type=Path, required=True)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--study-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    rows, summary = audit(args.reported_taxa, args.edges, args.nodes, args.study_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_rows(args.output_dir / "trait_coverage.tsv", rows)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: summary[key] for key in ("strains", "scopes")}, indent=2))


if __name__ == "__main__":
    main()
