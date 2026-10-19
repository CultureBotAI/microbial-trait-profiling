#!/usr/bin/env python3
"""Audit #629 production ingest without modifying scientific inputs (2026-10-19).

Inputs: raw BugSigDB CSV, new study collection/log, historical collection/log.
Output: JSON reason counts, raw-to-YAML leaf-taxon reconciliation, and baseline
membership differences. All directory reads include gitignored YAML files.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

import yaml


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def log_summary(path: Path) -> dict[str, object]:
    events = json.loads(path.read_text(encoding="utf-8"))
    counts = Counter((event["level"], event["entity"], event["reason"]) for event in events)
    groups = [
        {"level": level, "entity": entity, "reason": reason, "events": count}
        for (level, entity, reason), count in counts.most_common()
    ]
    return {
        "path": str(path),
        "sha256": digest(path),
        "events": len(events),
        "reasons": groups,
        "non_lineage_examples": [e for e in events if not e["reason"].startswith("lineage with ")][
            :12
        ],
    }


def signature_id(raw: str) -> str:
    return raw.replace("/", "-")


def collection(path: Path) -> dict[str, Any]:
    signatures: dict[str, list[str]] = {}
    studies: set[str] = set()
    experiments: set[str] = set()
    ranks: Counter[str] = Counter()
    hosts: Counter[str] = Counter()
    placeholders: list[dict[str, str]] = []
    files = sorted(path.glob("*.yaml"))
    for source in files:
        study = yaml.safe_load(source.read_text(encoding="utf-8"))
        studies.add(study["id"])
        for experiment in study.get("experiments", []):
            experiments.add(experiment["id"])
            hosts[str(experiment.get("host_species"))] += 1
            for signature in experiment.get("signatures", []):
                signatures[signature["id"]] = [taxon["id"] for taxon in signature.get("taxa", [])]
                for taxon in signature.get("taxa", []):
                    ranks[str(taxon.get("taxonomic_rank"))] += 1
                    if not taxon.get("taxon_name") or taxon["taxon_name"].startswith("taxon_"):
                        placeholders.append({"signature": signature["id"], **taxon})
    return {
        "files": len(files),
        "study_ids": studies,
        "experiment_ids": experiments,
        "signatures": signatures,
        "ranks": dict(ranks),
        "hosts": dict(hosts),
        "placeholder_or_blank_names": placeholders,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--studies", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--baseline-studies", type=Path, required=True)
    parser.add_argument("--baseline-log", type=Path, required=True)
    args = parser.parse_args()
    with args.raw.open(encoding="utf-8", newline="") as handle:
        raw = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    source_signatures: dict[str, list[str]] = {}
    source_studies: set[str] = set()
    source_experiments: set[str] = set()
    blank_taxa: list[dict[str, str]] = []
    chain_mismatches: list[dict[str, object]] = []
    lineage_levels: Counter[int] = Counter()
    for row in raw:
        key = row["BSDB ID"]
        source_studies.add(key.split("/")[0])
        source_experiments.add("-".join(key.split("/")[:2]))
        text = row["NCBI Taxonomy IDs"].strip()
        ids = (
            []
            if text in {"", "NA"}
            else [value.strip() for value in text.split(";") if value.strip()]
        )
        names = row.get("MetaPhlAn taxon names", "") or row.get("Taxon", "") or ""
        name_chains = [name.strip() for name in names.split(",") if name.strip()]
        leaves: list[str] = []
        for chain in ids:
            parts = [part.strip() for part in chain.split("|") if part.strip()]
            lineage_levels[len(parts)] += 1
            if parts:
                leaf = parts[-1]
                leaves.append(leaf if leaf.startswith("NCBITaxon:") else f"NCBITaxon:{leaf}")
        source_signatures[signature_id(key)] = leaves
        if not ids:
            blank_taxa.append(
                {
                    field: row.get(field, "")
                    for field in (
                        "BSDB ID",
                        "NCBI Taxonomy IDs",
                        "MetaPhlAn taxon names",
                        "Taxon",
                        "Abundance in Group 1",
                        "Title",
                        "Source",
                        "Description",
                    )
                }
            )
        elif len(ids) != len(name_chains):
            chain_mismatches.append(
                {"id": key, "id_chains": len(ids), "name_chains": len(name_chains)}
            )
    current = collection(args.studies)
    baseline = collection(args.baseline_studies)
    mismatches = [
        key for key, leaves in current["signatures"].items() if leaves != source_signatures.get(key)
    ]
    skipped = sorted(set(source_signatures) - set(current["signatures"]))
    source_log = log_summary(args.log)
    baseline_log = log_summary(args.baseline_log)
    current_names = {path.name for path in args.studies.glob("*.yaml")}
    baseline_names = {path.name for path in args.baseline_studies.glob("*.yaml")}
    shared = sorted(current_names.intersection(baseline_names))
    changed = [
        name
        for name in shared
        if digest(args.studies / name) != digest(args.baseline_studies / name)
    ]
    result = {
        "raw": {
            "path": str(args.raw),
            "sha256": digest(args.raw),
            "rows": len(raw),
            "studies": len(source_studies),
            "experiments": len(source_experiments),
            "distinct_signatures": len(source_signatures),
            "taxon_lineage_chains": sum(lineage_levels.values()),
            "lineage_levels": dict(sorted(lineage_levels.items())),
            "id_name_chain_count_mismatches": chain_mismatches,
        },
        "current_log": source_log,
        "baseline_log": baseline_log,
        "current_collection": {
            "files": current["files"],
            "studies": len(current["study_ids"]),
            "experiments": len(current["experiment_ids"]),
            "signatures": len(current["signatures"]),
            "taxon_observations": sum(len(value) for value in current["signatures"].values()),
            "ranks": current["ranks"],
            "hosts": current["hosts"],
            "placeholder_or_blank_names": current["placeholder_or_blank_names"],
            "leaf_id_mismatch_signatures": mismatches,
            "missing_source_studies": sorted(source_studies - current["study_ids"]),
            "missing_source_experiments": sorted(source_experiments - current["experiment_ids"]),
            "missing_source_signatures": skipped,
            "missing_signatures_with_nonempty_source_taxa": [
                key for key in skipped if source_signatures[key]
            ],
        },
        "blank_taxa_source_rows": blank_taxa,
        "baseline_collection": {
            "files": baseline["files"],
            "studies": len(baseline["study_ids"]),
            "experiments": len(baseline["experiment_ids"]),
            "signatures": len(baseline["signatures"]),
            "shared_files": len(shared),
            "changed_shared_files": changed,
            "only_current_files": sorted(current_names - baseline_names),
            "only_baseline_files": sorted(baseline_names - current_names),
        },
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
