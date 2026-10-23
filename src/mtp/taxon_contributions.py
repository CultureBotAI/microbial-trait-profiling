"""Reconcile source counts and export an auditable Table 1 contribution ledger."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from mtp.contribution_inputs import (
    RANKS,
    Base,
    Graph,
    Key,
    digest,
    encoded,
    load_base,
    load_graph,
    verify_kg,
)
from mtp.paths import DERIVED_DATA_DIR, REPO_ROOT
from mtp.taxa_trait_matrix import STRAIN_ROLLUPS

Record = dict[str, Any]


def write_table(path: Path, fields: list[str], records: Iterable[Record]) -> int:
    """Write explicit headers even for empty results; gzip large relational tables."""
    opener = gzip.open if path.suffix == ".gz" else open
    number = 0
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in records:
            writer.writerow(row)
            number += 1
    return number


@dataclass
class Ledger:
    sources: dict[str, set[str]]
    reconstructed: dict[str, Counter[str]]
    mismatches: list[Record]


def reconcile(base: Base, graph: Graph) -> Ledger:
    """Mirror the original deduplicated source-trait rollup, then profile summation."""
    targets = set(base.taxa) | {t for taxa in base.reported.values() for t in taxa}
    sources: dict[str, set[str]] = {t: set() for t in targets}
    reconstructed: dict[str, Counter[str]] = {t: Counter() for t in targets}
    for source, traits in graph.traits.items():
        for target in graph.ancestors(source) & targets:
            sources[target].add(source)
            reconstructed[target].update(traits)
    mismatches = []
    for target in sorted(targets):
        expected = base.taxa.get(target, {})
        for trait in base.columns:
            actual = reconstructed[target][trait]
            if actual != expected.get(trait, 0):
                mismatches.append(
                    {
                        "level": "taxon_matrix",
                        "key": target,
                        "trait_key": trait,
                        "expected": expected.get(trait, 0),
                        "reconstructed": actual,
                    }
                )
    for key, expected in sorted(base.profiles.items()):
        actual_profile: Counter[str] = Counter()
        for target in base.reported[key]:
            actual_profile.update(reconstructed[target])
        for trait in base.columns:
            if actual_profile[trait] != expected[trait]:
                mismatches.append(
                    {
                        "level": "profile",
                        "key": encoded(key),
                        "trait_key": trait,
                        "expected": expected[trait],
                        "reconstructed": actual_profile[trait],
                    }
                )
    return Ledger(sources, reconstructed, mismatches)


def profile_record(base: Base, key: Key) -> Record:
    meta = base.metadata[key[:2]]
    return {
        "study_id": key[0],
        "experiment_id": key[1],
        "direction": key[2],
        "condition": meta["condition"],
        "body_site": meta["body_site"],
    }


PROFILE_FIELDS = ["study_id", "experiment_id", "direction", "condition", "body_site"]


def export_ledger(base: Base, graph: Graph, ledger: Ledger, output: Path) -> dict[str, int]:
    """Normalized joins preserve all contributors without replicating edge metadata."""
    six = set(base.traits.values())

    def links() -> Iterator[Record]:
        for key, taxa in sorted(base.reported.items()):
            for target, signatures in sorted(taxa.items()):
                yield {
                    **profile_record(base, key),
                    "reported_taxon": target,
                    "reported_name": base.names.get(target, ""),
                    "reported_rank": base.ranks.get(target, ""),
                    "signature_ids": encoded(sorted(signatures)),
                }

    def units() -> Iterator[Record]:
        for target, sources in sorted(ledger.sources.items()):
            for source in sorted(sources):
                path, multiple = graph.path(source, target)
                for trait in sorted(graph.traits[source] & six):
                    evidence = graph.assertions.get((source, trait), [])
                    yield dict(
                        reported_taxon=target,
                        source_taxon=source,
                        source_name=graph.names.get(source, ""),
                        trait_key=trait,
                        count_contribution=1,
                        kg_edge_count=len(evidence),
                        primary_source_ids=encoded(
                            sorted({a.primary_source for a in evidence if a.primary_source})
                        ),
                        original_assertion_identity="unresolved_KG_edges_are_not_original_assertions",
                        rollup_path=encoded(path),
                        alternative_parent_paths=multiple,
                        **{rank: graph.ranked(source, rank, base.ranks)[0] for rank in RANKS},
                    )

    def assertions() -> Iterator[Record]:
        for (source, trait), evidence in sorted(graph.assertions.items()):
            for assertion in evidence:
                yield dict(source_taxon=source, trait_key=trait, **asdict(assertion))

    def profile_counts() -> Iterator[Record]:
        for key, taxa in sorted(base.reported.items()):
            for trait in base.traits.values():
                occurrences: Counter[str] = Counter()
                reported_count = 0
                for target in taxa:
                    contributing = {s for s in ledger.sources[target] if trait in graph.traits[s]}
                    occurrences.update(contributing)
                    reported_count += bool(contributing)
                yield {
                    **profile_record(base, key),
                    "trait_key": trait,
                    "profile_count": base.profiles[key][trait],
                    "distinct_reported_taxa": reported_count,
                    "distinct_source_taxa": len(occurrences),
                    "descendant_source_taxa": sum(s not in taxa for s in occurrences),
                    "unique_KG_edges": len(
                        {
                            (a.edge_id, a.source_row)
                            for s in occurrences
                            for a in graph.assertions.get((s, trait), [])
                        }
                    ),
                    "repeated_propagated_counts": sum(occurrences.values()) - len(occurrences),
                    "original_assertion_count": "",
                    "full_non_mediadive_total": sum(base.profiles[key].values()),
                }

    return {
        "profile_reported_taxa": write_table(
            output / "profile_reported_taxa.tsv.gz",
            [*PROFILE_FIELDS, "reported_taxon", "reported_name", "reported_rank", "signature_ids"],
            links(),
        ),
        "source_contributions": write_table(
            output / "source_contributions.tsv.gz",
            [
                "reported_taxon",
                "source_taxon",
                "source_name",
                "trait_key",
                "count_contribution",
                "kg_edge_count",
                "primary_source_ids",
                "original_assertion_identity",
                "rollup_path",
                "alternative_parent_paths",
                *RANKS,
            ],
            units(),
        ),
        "source_assertions": write_table(
            output / "source_assertions.tsv.gz",
            ["source_taxon", "trait_key", *assertion_fields()],
            assertions(),
        ),
        "profile_counts": write_table(
            output / "profile_counts.tsv.gz",
            [
                *PROFILE_FIELDS,
                "trait_key",
                "profile_count",
                "distinct_reported_taxa",
                "distinct_source_taxa",
                "descendant_source_taxa",
                "unique_KG_edges",
                "repeated_propagated_counts",
                "original_assertion_count",
                "full_non_mediadive_total",
            ],
            profile_counts(),
        ),
    }


def assertion_fields() -> list[str]:
    return [
        "edge_id",
        "predicate",
        "object_id",
        "primary_source",
        "knowledge_source",
        "knowledge_level",
        "source_row",
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kg-dir", type=Path, default=REPO_ROOT / "data/kg-microbe")
    parser.add_argument("--inputs", type=Path, default=DERIVED_DATA_DIR)
    parser.add_argument("--output", type=Path, default=DERIVED_DATA_DIR / "table1_taxon_influence")
    parser.add_argument(
        "--protocol", type=Path, default=REPO_ROOT / "methods/table1_taxon_influence_protocol.json"
    )
    parser.add_argument(
        "--strain-rollup",
        choices=STRAIN_ROLLUPS,
        default="off",
        help=(
            "The BacDive strain rollup scope the input matrix was built with, so that the "
            "ledger reconstructs its counts by the same rule (#788)."
        ),
    )
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    output: Path = args.output
    output.mkdir(parents=True, exist_ok=True)
    print("Verifying pinned KG and input counts", flush=True)
    hashes = verify_kg(args.kg_dir, protocol["kg_md5"])
    names = [
        "study_trait_profiles_wide.tsv",
        "experiments.tsv",
        "bugsigdb_taxa.tsv",
        "taxa_trait_matrix.tsv",
    ]
    base = load_base(*(args.inputs / name for name in names))
    hashes.update({name: digest(args.inputs / name) for name in names})
    hashes[args.protocol.name] = digest(args.protocol)
    print(
        f"Loading KG for {len(base.profiles)} profiles, {len(base.columns)} denominator traits",
        flush=True,
    )
    graph = load_graph(
        args.kg_dir / "merged-kg_edges.tsv",
        args.kg_dir / "merged-kg_nodes.tsv",
        set(base.columns),
        set(base.traits.values()),
        strain_rollup=args.strain_rollup,
    )
    print("Reconciling all matrix and profile counts", flush=True)
    ledger = reconcile(base, graph)
    write_table(
        output / "reconciliation.tsv",
        ["level", "key", "trait_key", "expected", "reconstructed"],
        ledger.mismatches,
    )
    write_table(
        output / "input_eligibility.tsv",
        ["study_id", "experiment_id", "direction", "status", "orientation_status"],
        base.eligibility,
    )
    manifest: Record = {
        "completed": False,
        "input_sha256": hashes,
        "protocol": protocol,
        "traits": base.traits,
        "denominator_columns": base.columns,
        "strain_rollup": args.strain_rollup,
        "profiles": len(base.profiles),
        "studies": len({k[0] for k in base.profiles}),
        "header_continuation_repaired": graph.header_repaired,
        "agent_field_continued_records": graph.continued_records,
        "mismatches": len(ledger.mismatches),
        "code_sha256": {
            name: digest(Path(__file__).parent / name)
            for name in ("contribution_inputs.py", "taxon_contributions.py", "taxon_influence.py")
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if ledger.mismatches:
        raise ValueError(
            f"{len(ledger.mismatches)} count mismatches; influence calculations blocked. "
            f"Inspect {output / 'reconciliation.tsv'}"
        )
    print("Exporting reconciled source ledger", flush=True)
    manifest["rows"] = export_ledger(base, graph, ledger, output)
    # Import here to keep count reconstruction independently usable and testable.
    from mtp.taxon_influence import run_influence

    manifest["influence"] = run_influence(base, graph, ledger, protocol, output)
    manifest["completed"] = True
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote audit to {output}", flush=True)


if __name__ == "__main__":
    main()


__all__ = [
    "PROFILE_FIELDS",
    "Ledger",
    "Record",
    "assertion_fields",
    "export_ledger",
    "main",
    "profile_record",
    "reconcile",
    "write_table",
]
