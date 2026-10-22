"""Compare trait coverage with and without the BacDive strain rollup (#774, 2026-10-21).

Inputs: the reported taxa the matrix is built from, the KG edges and nodes, and
optionally a contrast-selected study directory. Outputs, in --output-dir:
trait_coverage.tsv, one row per trait key that any NCBITaxon node or BacDive strain
carries, and summary.json with the strain links, the vocabulary sizes, the focal
traits and, with --study-dir, the signature taxa that gain a trait label.

"Off" is the default matrix; "on" counts strain trait edges at the strains'
NCBITaxon nodes (taxa_trait_matrix --strain-rollup). The evaluable vocabulary is
what retained_trait_columns keeps of the matrix columns: every trait some reported
taxon carries, without MediaDive growth media. All counts are annotation coverage.
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
    NodeMetadata,
    compute_trait_counts,
    node_ids_for_trait_keys,
    read_kg_trait_edges,
    read_node_metadata,
    read_reported_taxa,
    roll_up_strain_traits,
    split_trait_key,
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
NEGATED_PREDICATE_PREFIXES = ("does not ", "is not ")
COVERAGE_COLUMNS = (
    "trait_key",
    "label",
    "predicate_id",
    "predicate_name",
    "negated_predicate",
    "culture_medium",
    "focal_trait",
    "ncbitaxon_nodes_off",
    "ncbitaxon_nodes_on",
    "linked_strains",
    "unlinked_strains",
    "reported_taxa_off",
    "reported_taxa_on",
    "vocabulary_off",
    "vocabulary_on",
    "vocabulary_change",
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


def vocabulary_change(off: bool, on: bool) -> str:
    if off and on:
        return "stays"
    if on:
        return "enters"
    return "leaves" if off else "outside"


def predicate_name(trait_key: str, metadata: dict[str, NodeMetadata]) -> str:
    predicate, _ = split_trait_key(trait_key)
    if predicate is None:
        return ""
    return metadata.get(predicate, NodeMetadata("", "")).name


def signature_taxa(
    study_dir: Path,
    host_species: str | None,
    labels_off: dict[str, set[str]],
    labels_on: dict[str, set[str]],
) -> dict[str, object]:
    """Distinct taxa of the 16S and shotgun profiles, and the trait labels they gain."""
    profiles, _ = collect_profile_taxa(
        study_dir, host_species=host_species, sequencing_types=CORPUS_SEQUENCING_TYPES
    )
    taxa = {taxon for observed in profiles.values() for taxon in observed}
    matched = sorted(taxa.intersection(labels_off))
    gaining = [taxon for taxon in matched if labels_on[taxon] - labels_off[taxon]]
    return {
        "host_species": host_species or "all",
        "profiles": len(profiles),
        "signature_taxa": len(taxa),
        "matched_to_matrix_rows": len(matched),
        "with_a_retained_trait_off": sum(bool(labels_off[taxon]) for taxon in matched),
        "with_a_retained_trait_on": sum(bool(labels_on[taxon]) for taxon in matched),
        "gaining_a_trait_label": len(gaining),
        "gaining_a_first_trait_label": sum(not labels_off[taxon] for taxon in gaining),
        "retained_traits_carried_off": len({k for taxon in matched for k in labels_off[taxon]}),
        "retained_traits_carried_on": len({k for taxon in matched for k in labels_on[taxon]}),
        "focal_traits": {
            key: {
                "label": label,
                "taxa_off": sum(key in labels_off[taxon] for taxon in matched),
                "taxa_on": sum(key in labels_on[taxon] for taxon in matched),
            }
            for key, label in FOCAL_TRAITS.items()
        },
    }


def audit(
    reported_taxa_path: Path,
    edges_path: Path,
    nodes_path: Path,
    study_dir: Path | None,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Return the per-trait coverage rows and the summary, checking that on only adds."""
    reported = read_reported_taxa(reported_taxa_path)
    edges = read_kg_trait_edges(edges_path, strains=True)
    taxon_traits_on = roll_up_strain_traits(
        edges.taxon_traits, edges.strain_taxa, edges.strain_traits
    )
    counts_off = compute_trait_counts(reported, edges.child_to_parents, edges.taxon_traits)
    counts_on = compute_trait_counts(reported, edges.child_to_parents, taxon_traits_on)
    for taxon, row in counts_off.items():
        if any(counts_on[taxon][key] < n for key, n in row.items()):
            raise ValueError(f"The rollup lowered a count of {taxon}")

    nodes_off = carried(edges.taxon_traits.values())
    nodes_on = carried(taxon_traits_on.values())
    linked = [s for s in edges.strain_traits if edges.strain_taxa.get(s)]
    unlinked = [s for s in edges.strain_traits if not edges.strain_taxa.get(s)]
    strains_linked = carried(edges.strain_traits[s] for s in linked)
    strains_unlinked = carried(edges.strain_traits[s] for s in unlinked)
    labels_off_all = present(counts_off)
    labels_on_all = present(counts_on)
    taxa_off = carried(labels_off_all.values())
    taxa_on = carried(labels_on_all.values())
    columns_off, columns_on = set(taxa_off), set(taxa_on)
    vocabulary_off, vocabulary_on = retained(columns_off), retained(columns_on)
    if not vocabulary_off <= vocabulary_on:
        raise ValueError("A trait left the vocabulary under the rollup")

    keys = set(nodes_on) | set(strains_linked) | set(strains_unlinked)
    metadata = read_node_metadata(nodes_path, node_ids_for_trait_keys(keys))
    rows: list[dict[str, object]] = []
    for key in sorted(keys, key=lambda k: (trait_column_label(k, metadata).casefold(), k)):
        predicate, _ = split_trait_key(key)
        name = predicate_name(key, metadata)
        rows.append(
            {
                "trait_key": key,
                "label": trait_column_label(key, metadata),
                "predicate_id": predicate or "",
                "predicate_name": name,
                "negated_predicate": name.startswith(NEGATED_PREDICATE_PREFIXES),
                "culture_medium": MEDIADIVE_MARKER in key.lower(),
                "focal_trait": FOCAL_TRAITS.get(key, ""),
                "ncbitaxon_nodes_off": nodes_off[key],
                "ncbitaxon_nodes_on": nodes_on[key],
                "linked_strains": strains_linked[key],
                "unlinked_strains": strains_unlinked[key],
                "reported_taxa_off": taxa_off[key],
                "reported_taxa_on": taxa_on[key],
                "vocabulary_off": key in vocabulary_off,
                "vocabulary_on": key in vocabulary_on,
                "vocabulary_change": vocabulary_change(key in vocabulary_off, key in vocabulary_on),
            }
        )

    entering = [row for row in rows if row["vocabulary_change"] == "enters"]
    by_predicate = Counter(
        str(row["predicate_name"] or "phenotype (METPO object)") for row in entering
    )
    link_counts = Counter(len(edges.strain_taxa.get(s, ())) for s in edges.strain_traits)
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
            "ncbitaxon_nodes_linked_from_trait_bearing_strains": len(
                {taxon for s in linked for taxon in edges.strain_taxa[s]}
            ),
        },
        "ncbitaxon_nodes_with_a_trait": {
            "off": len(edges.taxon_traits),
            "on": len(taxon_traits_on),
        },
        "vocabulary": {
            "matrix_columns_off": len(columns_off),
            "matrix_columns_on": len(columns_on),
            "retained_off": len(vocabulary_off),
            "retained_on": len(vocabulary_on),
            "entering": len(entering),
            "leaving": len(vocabulary_off - vocabulary_on),
            "entering_with_negated_predicate": sum(bool(r["negated_predicate"]) for r in entering),
            "entering_by_predicate": dict(by_predicate.most_common()),
        },
        "focal_traits": {
            key: {
                "label": label,
                "ncbitaxon_nodes_off": nodes_off[key],
                "ncbitaxon_nodes_on": nodes_on[key],
                "linked_strains": strains_linked[key],
                "unlinked_strains": strains_unlinked[key],
                "reported_taxa_off": taxa_off[key],
                "reported_taxa_on": taxa_on[key],
                "vocabulary_change": vocabulary_change(key in vocabulary_off, key in vocabulary_on),
            }
            for key, label in FOCAL_TRAITS.items()
        },
    }
    if study_dir is not None:
        labels_off = {t: labels & vocabulary_off for t, labels in labels_off_all.items()}
        labels_on = {t: labels & vocabulary_on for t, labels in labels_on_all.items()}
        summary["signature_taxa"] = {
            "study_dir": str(study_dir),
            "cohorts": [
                signature_taxa(study_dir, host, labels_off, labels_on)
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
    print(json.dumps({k: summary[k] for k in ("strains", "vocabulary")}, indent=2))


if __name__ == "__main__":
    main()
