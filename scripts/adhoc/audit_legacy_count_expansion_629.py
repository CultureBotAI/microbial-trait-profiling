#!/usr/bin/env python3
"""Explain production count-profile expansion for issue #629 (2026-10-19).

Inputs: saved core comparison JSON, its completed run manifest, actual baseline
and current count profiles, published experiment table and feature catalog.
Output: JSON cohort reconciliation and feature-level counts of nonzero cells.
This read-only audit distinguishes scope expansion from shared-cell changes;
it does not infer the historical producer's provenance from matching values.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

KEY_COLUMNS = ("study_id", "experiment_id", "inc_dec")
HUMAN = "Homo sapiens"


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_keys(path: Path) -> tuple[list[str], set[tuple[str, ...]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        keys = set()
        for row in reader:
            key = tuple(row[column] for column in KEY_COLUMNS)
            if key in keys:
                raise ValueError(f"Duplicate profile key in {path}: {key}")
            keys.add(key)
    return [column for column in columns if column not in KEY_COLUMNS], keys


def audit(comparison_path: Path) -> dict[str, object]:
    comparison = json.loads(comparison_path.read_text(encoding="utf-8"))
    manifest = Path(comparison["run_manifest"])
    if digest(manifest) != comparison["run_manifest_sha256"]:
        raise ValueError("Run manifest changed since the saved comparison")
    run = json.loads(manifest.read_text(encoding="utf-8"))
    if run["status"] != "completed":
        raise ValueError("Run must be completed")
    artifacts = {record["product_id"]: record for record in run["artifacts"]}

    def product_path(product: str) -> Path:
        path = Path(artifacts[product]["path"])
        return path if path.is_absolute() else manifest.parent / path

    count_comparison = comparison["comparisons"]["count_profiles"]
    baseline_path = Path(count_comparison["baseline_path"])
    current_path = Path(count_comparison["current_path"])
    for prefix, path in (("baseline", baseline_path), ("current", current_path)):
        if digest(path) != count_comparison[f"{prefix}_sha256"]:
            raise ValueError(f"{prefix} count profile changed since the saved comparison")
    experiment_path = product_path("experiments")
    with experiment_path.open(encoding="utf-8", newline="") as handle:
        experiments = {
            (row["study_id"], row["experiment_id"]): row["host_species"].strip()
            for row in csv.DictReader(handle, delimiter="\t")
        }
    catalog_path = product_path("feature_catalog")
    with catalog_path.open(encoding="utf-8", newline="") as handle:
        catalog = {row["column_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    baseline_columns, baseline_keys = read_keys(baseline_path)
    current_columns, current_keys = read_keys(current_path)
    added_keys = current_keys - baseline_keys
    removed_keys = baseline_keys - current_keys
    added_columns = sorted(set(current_columns) - set(baseline_columns))
    if added_keys != {tuple(key) for key in count_comparison["added_rows"]}:
        raise ValueError("Added profile keys differ from saved comparison")
    if removed_keys != {tuple(key) for key in count_comparison["removed_rows"]}:
        raise ValueError("Removed profile keys differ from saved comparison")
    if added_columns != sorted(count_comparison["added_features"]):
        raise ValueError("Added features differ from saved comparison")

    def host(key: tuple[str, ...]) -> str:
        return experiments[key[:2]]

    def host_counts(keys: set[tuple[str, ...]]) -> dict[str, int]:
        return dict(sorted(Counter(host(key) for key in keys).items()))

    human_keys = {key for key in current_keys if host(key) == HUMAN}
    feature_counts: dict[str, Counter[str]] = {column: Counter() for column in added_columns}
    examples: list[dict[str, object]] = []
    with current_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = tuple(row[column] for column in KEY_COLUMNS)
            cohort = "human" if host(key) == HUMAN else "nonhuman"
            for column in added_columns:
                value = int(row[column])
                if value:
                    feature_counts[column][f"{cohort}_nonzero_profiles"] += 1
                    feature_counts[column][f"{cohort}_count_sum"] += value
                    feature_counts[column][f"{cohort}_max_count"] = max(
                        feature_counts[column][f"{cohort}_max_count"], value
                    )
                    if len(examples) < 8:
                        examples.append(
                            {
                                "profile_key": key,
                                "host": host(key),
                                "feature": column,
                                "count": value,
                            }
                        )
    features = [
        {
            "column_id": column,
            "label": catalog[column]["column_label"],
            "feature_class": catalog[column]["feature_class"],
            **{
                f"{cohort}_{metric}": feature_counts[column][f"{cohort}_{metric}"]
                for cohort in ("human", "nonhuman")
                for metric in ("nonzero_profiles", "count_sum", "max_count")
            },
        }
        for column in added_columns
    ]
    return {
        "run_id": run["run_id"],
        "comparison_path": str(comparison_path.resolve()),
        "comparison_sha256": digest(comparison_path),
        "manifest_path": str(manifest.resolve()),
        "manifest_sha256": digest(manifest),
        "experiment_table_sha256": digest(experiment_path),
        "catalog_sha256": digest(catalog_path),
        "baseline_path": str(baseline_path),
        "baseline_sha256": digest(baseline_path),
        "current_path": str(current_path),
        "current_sha256": digest(current_path),
        "baseline_rows": len(baseline_keys),
        "current_rows": len(current_keys),
        "added_rows": len(added_keys),
        "removed_rows": len(removed_keys),
        "baseline_profiles_by_host": host_counts(baseline_keys),
        "added_profiles_by_host": host_counts(added_keys),
        "shared_profiles_by_host": host_counts(baseline_keys & current_keys),
        "baseline_keys_equal_current_human_keys": baseline_keys == human_keys,
        "added_rows_all_nonhuman": all(host(key) != HUMAN for key in added_keys),
        "added_rows_missing_host": sum(not host(key) for key in added_keys),
        "baseline_features": len(baseline_columns),
        "current_features": len(current_columns),
        "added_features": len(added_columns),
        "added_feature_classes": dict(Counter(row["feature_class"] for row in features)),
        "added_features_all_zero_in_humans": all(
            counts["human_nonzero_profiles"] == 0 for counts in feature_counts.values()
        ),
        "added_features_all_zero_in_all_profiles": sum(
            counts["human_nonzero_profiles"] == 0 and counts["nonhuman_nonzero_profiles"] == 0
            for counts in feature_counts.values()
        ),
        "added_feature_nonzero_cell_examples": examples,
        "added_feature_details": features,
        "shared_compared_cells": count_comparison["compared_cells"],
        "shared_changed_cells": count_comparison["changed_cells"],
        "interpretation": (
            "Cohort membership is joined to the published experiment table. "
            "A zero means no annotation count in this product, not biological absence. "
            "Historical baseline provenance remains imported_unverified."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("comparison", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.comparison), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
