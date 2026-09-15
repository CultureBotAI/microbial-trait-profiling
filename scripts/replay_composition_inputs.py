"""Reconcile the frozen #174 matrix against pinned raw KG and BugSigDB inputs."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from bsdbng.contribution_inputs import digest, load_graph, verify_kg
from bsdbng.ingest import _parse_experiment_id, _parse_study_id, _parse_taxa, _read_csv
from bsdbng.study_trait_profiles import direction_to_inc_dec
from bsdbng.taxa_trait_matrix import ReportedTaxon, compute_trait_counts
from bsdbng.trait_composition import DENOMINATOR_DEFINITION, retained_trait_counts

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / ".tmp/issue174"
REPORT = ROOT / "trait-visualization/denominator-harmonization"
KG_MD5 = {
    "merged-kg_edges.tsv": "555ecf5ba5d89288cdc44e368375374f",
    "merged-kg_nodes.tsv": "dd76b4656f9341db25a7c4afa96cd58e",
}
FULL_DUMP_SHA = "af84bfca7bfe76913a20c1fedbaadbfad4e4bf78605a2cfd7d00eb0c0fdd5788"


def replay(data: Path, kg: Path) -> None:
    dump = data / "raw/full_dump.csv"
    if digest(dump) != FULL_DUMP_SHA:
        raise ValueError("BugSigDBExports v1.3.1 full_dump.csv hash mismatch")
    hashes = verify_kg(kg, KG_MD5)
    source = WORK / "study_trait_profiles_wide.tsv"
    if digest(source) != "c07d589f6cc88724e43c683ebd132bcd5c7f316f36b95df3ff9acad24ebbf62f":
        raise ValueError("Issue #174 profile matrix hash mismatch")
    counts = retained_trait_counts(
        pd.read_csv(source, sep="\t").set_index(["study_id", "experiment_id", "inc_dec"])
    )
    reported: defaultdict[tuple[str, str, str], set[str]] = defaultdict(set)
    taxa: dict[str, ReportedTaxon] = {}
    parse_log: list[dict[str, str]] = []
    for row in _read_csv(dump):
        direction = direction_to_inc_dec(row["Abundance in Group 1"])
        if direction is None:
            continue
        signature = row["BSDB ID"]
        key = (_parse_study_id(signature), _parse_experiment_id(signature), direction)
        if key not in counts.index:
            continue
        reported[key]  # Preserve explicitly empty signatures.
        for taxon in _parse_taxa(
            row["NCBI Taxonomy IDs"], row["MetaPhlAn taxon names"], signature, parse_log
        ):
            reported[key].add(taxon["id"])
            taxa[taxon["id"]] = ReportedTaxon(
                taxon["id"], taxon["taxon_name"], taxon.get("taxonomic_rank", "")
            )
    skipped = [entry for entry in parse_log if entry["level"] == "skip"]
    if skipped:
        raise ValueError(f"Unparsed reported taxa: {skipped[:3]}")
    print(
        f"Replaying {len(taxa)} reported taxa, {len(counts)} profiles, "
        f"{len(counts.columns)} retained traits",
        flush=True,
    )
    graph = load_graph(
        kg / "merged-kg_edges.tsv", kg / "merged-kg_nodes.tsv", set(counts.columns), set()
    )
    taxon_counts = compute_trait_counts(list(taxa.values()), graph.parents, graph.traits)
    paired = pd.read_csv(
        ROOT / "trait-visualization/study_trait_profile_paired_enrichment_metadata.tsv",
        sep="\t",
        index_col=0,
    )
    paired_ids = set(paired["increased_profile_ids"]) | set(paired["decreased_profile_ids"])
    records = []
    discrepancies = []
    for key, row in counts.iterrows():
        total: Counter[str] = Counter()
        for taxon in reported[key]:
            total.update(taxon_counts[taxon])
        changed = [trait for trait in counts.columns if total[trait] != row[trait]]
        records.append(
            {
                "study_id": key[0],
                "experiment_id": key[1],
                "inc_dec": key[2],
                "paired": "|".join(key) in paired_ids,
                "reported_taxa": len(reported[key]),
                "frozen_retained_total": int(row.sum()),
                "replayed_retained_total": sum(total.values()),
                "mismatching_traits": len(changed),
            }
        )
        for trait in changed:
            discrepancies.append(
                {
                    "profile_id": "|".join(key),
                    "trait_id": trait,
                    "frozen_count": int(row[trait]),
                    "replayed_count": total[trait],
                }
            )
    REPORT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(REPORT / "raw_input_replay_profiles.tsv", sep="\t", index=False)
    pd.DataFrame(
        discrepancies, columns=["profile_id", "trait_id", "frozen_count", "replayed_count"]
    ).to_csv(REPORT / "raw_input_replay_discrepancies.tsv", sep="\t", index=False)
    manifest = {
        "code_commit": (WORK / "code_commit.txt").read_text().strip(),
        "source_sha256": {
            "scripts/replay_composition_inputs.py": digest(Path(__file__)),
            "src/bsdbng/contribution_inputs.py": digest(ROOT / "src/bsdbng/contribution_inputs.py"),
        },
        "kg_sha256": hashes,
        "full_dump_sha256": digest(dump),
        "profile_matrix_sha256": digest(source),
        "denominator": DENOMINATOR_DEFINITION,
        "profile_count": len(counts),
        "paired_profile_count": len(paired_ids),
        "retained_trait_count": len(counts.columns),
        "mismatching_cells": len(discrepancies),
        "profiles_with_mismatches": sum(record["mismatching_traits"] > 0 for record in records),
        "paired_profiles_with_mismatches": sum(
            record["paired"] and record["mismatching_traits"] > 0 for record in records
        ),
        "kg_header_continuation_repaired": graph.header_repaired,
        "kg_agent_field_continuations_repaired": graph.continued_records,
        "contract": "Deduplicate reported taxa within each study/experiment/direction. "
        "Each distinct KG trait-bearing taxon contributes once to each reported ancestor. "
        "Preserve nested reported taxa and original multiplicity; no influence model change.",
    }
    (REPORT / "raw_input_replay_provenance.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    if discrepancies:
        raise ValueError("Raw-input replay differs from the frozen matrix; see discrepancy tables")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--kg-dir", type=Path, required=True)
    args = parser.parse_args()
    replay(args.data_root.resolve(), args.kg_dir.resolve())


if __name__ == "__main__":
    main()
