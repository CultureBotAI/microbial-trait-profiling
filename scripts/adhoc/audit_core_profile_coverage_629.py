#!/usr/bin/env python3
"""Inspect completed #629 core annotation products without rerunning analyses (2026-10-19).

Input: completed production run.json and its published matrix, catalog, profiles,
and coverage TSVs. Output: JSON dimensions, annotation coverage, and selected
profile cells. Counts describe annotations/lookup memberships, not organisms.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def audit(manifest: Path) -> dict[str, object]:
    run = json.loads(manifest.read_text(encoding="utf-8"))
    if run["status"] != "completed":
        raise ValueError("Inspect only a completed run's published artifacts")
    artifacts = {record["product_id"]: record for record in run["artifacts"]}

    def product_path(product: str) -> Path:
        path = Path(artifacts[product]["path"])
        return path if path.is_absolute() else manifest.parent / path

    with product_path("feature_catalog").open(encoding="utf-8", newline="") as handle:
        catalog = list(csv.DictReader(handle, delimiter="\t"))
    retained = [row["column_id"] for row in catalog if row["feature_class"] != "culture_medium"]
    selected = {"METPO:1000608", "METPO:1000668", "METPO:1000683"}
    selected_columns = [row["column_id"] for row in catalog if row["column_id"] in selected]
    matrix_status: Counter[str] = Counter()
    matrix_status_by_rank: Counter[tuple[str, str]] = Counter()
    with product_path("taxa_matrix").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            nonmedia = any(int(row[column]) for column in retained)
            all_features = nonmedia or any(int(row[column["column_id"]]) for column in catalog)
            status = (
                "retained_annotation" if nonmedia else "media_only" if all_features else "all_zero"
            )
            matrix_status[status] += 1
            matrix_status_by_rank[(row["taxonomic_rank"], status)] += 1
    profiles: dict[str, object] = {}
    for prefix in ("count", "ora"):
        path = product_path(f"{prefix}_profiles")
        zeros: Counter[str] = Counter()
        selected_rows: list[dict[str, object]] = []
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                zeros["profiles"] += 1
                nonmedia = any(int(row[column]) for column in retained)
                all_features = nonmedia or any(int(row[column["column_id"]]) for column in catalog)
                zeros["zero_all_features"] += not all_features
                zeros["zero_nonmedia_features"] += not nonmedia
                zeros["media_only_profiles"] += all_features and not nonmedia
                if row["study_id"] == "bsdb:27115405":
                    selected_rows.append(
                        {
                            "study_id": row["study_id"],
                            "experiment_id": row["experiment_id"],
                            "direction": row["inc_dec"],
                            "taxa_count": int(row["taxa_count"]) if prefix == "ora" else None,
                            "selected_feature_values": {
                                column: int(row[column]) for column in selected_columns
                            },
                        }
                    )
        coverage_path = product_path(f"{prefix}_coverage")
        with coverage_path.open(encoding="utf-8", newline="") as handle:
            coverage = list(csv.DictReader(handle, delimiter="\t"))
        count_columns = [
            "reported_distinct_taxa",
            "lookup_matched_taxa",
            "lookup_unmatched_taxa",
            "matched_taxa_with_annotations",
            "matched_taxa_without_annotations",
        ]
        profiles[prefix] = {
            "profile_path": str(path),
            "profile_sha256": digest(path),
            "coverage_path": str(coverage_path),
            "coverage_sha256": digest(coverage_path),
            "cohort": artifacts[f"{prefix}_profiles"]["cohort"],
            "profile_zero_status": dict(zeros),
            "coverage_status_profiles": dict(Counter(row["coverage_status"] for row in coverage)),
            "profile_taxon_memberships": {
                column: sum(int(row[column]) for row in coverage) for column in count_columns
            },
            "selected_study_cells": selected_rows,
            "feature_filter_ids": sorted({row["trait_filter_id"] for row in coverage}),
            "lineage_resolution": sorted({row["lineage_resolution"] for row in coverage}),
        }
    return {
        "run_id": run["run_id"],
        "manifest": str(manifest.resolve()),
        "manifest_sha256": digest(manifest),
        "execution_fingerprint": run["execution_fingerprint"],
        "features": len(catalog),
        "retained_nonmedia_features": len(retained),
        "culture_medium_features": len(catalog) - len(retained),
        "matrix_annotation_status_rows": dict(matrix_status),
        "matrix_annotation_status_by_rank": [
            {"rank": rank, "annotation_status": status, "rows": count}
            for (rank, status), count in sorted(matrix_status_by_rank.items())
        ],
        "selected_feature_labels": {
            row["column_id"]: row["column_label"] for row in catalog if row["column_id"] in selected
        },
        "profiles": profiles,
        "scope": (
            "Core input cohort; no pairing or zero-total exclusions applied. "
            "Membership sums reuse taxa across profiles."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.manifest), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
