"""Extract the values quoted in the final rerun record from the run directories.

Purpose: reproduce every run ID, count, scope, verdict and additive figure that
docs/reviews/2026-10-23-final-rerun-record.md quotes, for #1004.
Inputs: --runs-root, the directory holding one subdirectory per unit
(core/run-*/run.json, arms/{ilyas,jed,noel}/arm_run.json and
comparison/comparison_run.json); --additive-root, the Tier 1+2 x off additive
directory (scan835/, influence-extension/, label-support/).
Output: JSON on stdout.
Date: 2026-10-23. Read-only.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
from typing import Any

ARMS = ("ilyas", "jed", "noel")
CELLULOSE = "METPO:2000103|GO:0030245"
LINEAGES = (
    "Lachnospiraceae",
    "Oscillospiraceae",
    "Lachnospiraceae + Oscillospiraceae",
    "Ruminococcus",
    "Acetivibrio",
    "Ruminiclostridium",
    "Clostridiaceae",
)


def read_json(path: Path) -> dict[str, Any]:
    """Read one JSON object."""
    with path.open(encoding="utf-8") as handle:
        data: dict[str, Any] = json.load(handle)
    return data


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a tab-separated table, gzipped or not."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def unit_values(unit_dir: Path) -> dict[str, Any]:
    """Summarise one verified unit: identity, scope, counts and verdicts."""
    core_path = next(unit_dir.glob("core/run-*/run.json"))
    core = read_json(core_path)
    parameters = {item["name"]: item["value"] for item in core["parameters"]}
    profiles = read_tsv(
        core_path.parent / "steps/count_profiles/derived/study_trait_profiles_wide.tsv"
    )
    arms = {}
    for arm in ARMS:
        record = read_json(unit_dir / "arms" / arm / "arm_run.json")
        arms[arm] = {
            "id": record["id"],
            "status": record["status"],
            "core_run_id": record["core_run_id"],
            "git_commit": record["git_commit"],
            "git_dirty": record["git_dirty"],
        }
    comparison_path = unit_dir / "comparison/comparison_run.json"
    comparison = read_json(comparison_path) if comparison_path.exists() else None
    values: dict[str, Any] = {
        "core_run_id": core["run_id"],
        "core_status": core["status"],
        "git_commit": core["git_commit"],
        "git_dirty": core["git_dirty"],
        "contrast_population": core["contrast_population"],
        "strain_rollup": parameters.get("strain_rollup"),
        "count_host_species": parameters.get("count_host_species"),
        "ora_host_species": parameters.get("ora_host_species"),
        "profiles": len(profiles),
        "experiments": len({(row["study_id"], row["experiment_id"]) for row in profiles}),
        "studies": len({row["study_id"] for row in profiles}),
        "arms": arms,
        "comparison": None,
    }
    if comparison is not None:
        values["comparison"] = {
            "id": comparison["id"],
            "status": comparison["status"],
            "core_run_id": comparison["core_run_id"],
            "git_commit": comparison["git_commit"],
            "git_dirty": comparison["git_dirty"],
            "bases": {
                basis["id"]: [basis["retained_feature_count"], basis["rollup_policy"]]
                for basis in comparison["analysis_bases"]
            },
            "verdicts": {
                verdict["id"]: verdict["outcome"]
                for verdict in comparison.get("calibration_verdicts", [])
            },
        }
    return values


def scan_values(scan_dir: Path) -> dict[str, Any]:
    """Summarise the condition-triple scan for the cellulose trait."""
    summary = read_json(scan_dir / "selection_scan_summary.json")
    conditions = read_tsv(scan_dir / "selection_scan_conditions.tsv")
    eligible = [row for row in conditions if row["eligible"] == "True"]
    by_studies = sorted(eligible, key=lambda row: -int(row["stool_studies"]))
    rule_triple = {row["condition"] for row in by_studies[:3]}
    triples = read_tsv(scan_dir / "selection_scan_triples.tsv.gz")
    tested = [row for row in triples if row["u6_p_mean_distance"] not in ("", "NA")]
    lowest = min(tested, key=lambda row: float(row["u6_p_mean_distance"]))

    def triple_row(names: set[str]) -> dict[str, str]:
        for row in triples:
            found = {row["condition_a"], row["condition_b"], row["condition_c"]}
            if row["trait_id"] == CELLULOSE and found == names:
                return row
        raise LookupError(f"no cellulose row for {sorted(names)}")

    def ranks(row: dict[str, str]) -> dict[str, str]:
        keys = (
            "mean_pairwise_jaccard_distance",
            "u6_expected_mean_distance",
            "u6_p_mean_distance",
            "rank_mean_distance",
            "u6_p_core",
            "rank_core",
        )
        return {key: row[key] for key in keys}

    lead = summary["lead_triple"]["conditions"]
    return {
        "counts": summary["counts"],
        "rule": {key: summary["rule"][key] for key in ("permutations", "seed", "min_studies")},
        "lead_triple": {"conditions": lead, **ranks(triple_row(set(lead)))},
        "rule_triple": {
            "conditions": sorted(rule_triple),
            "stool_studies": {row["condition"]: row["stool_studies"] for row in by_studies[:3]},
            **ranks(triple_row(rule_triple)),
        },
        "lowest_p_trait_triple": {
            "trait_id": lowest["trait_id"],
            "conditions": [lowest["condition_a"], lowest["condition_b"], lowest["condition_c"]],
            "u6_p_mean_distance": lowest["u6_p_mean_distance"],
            "tested_trait_triples": len(tested),
        },
    }


def lineage_values(extension_dir: Path) -> dict[str, list[str]]:
    """Return pooled and matched cellulose effect reductions for the record's lineages."""
    rows = read_tsv(extension_dir / "lineage_fraction_effects.tsv")
    values = {}
    for row in rows:
        if row["scope"] == "pooled" and row["trait_key"] == CELLULOSE and row["name"] in LINEAGES:
            values[row["name"]] = [
                row["effect_reduction_percent"],
                row["matched_effect_reduction_percent"],
            ]
    provenance = read_json(extension_dir / "extension_provenance.json")
    values["reconciliation_mismatches"] = [str(provenance["reconciliation_mismatches"])]
    return values


def label_support_values(label_dir: Path) -> dict[str, Any]:
    """Return the label-support audit's verification fields."""
    provenance = read_json(label_dir / "provenance.json")
    keys = ("verification_status", "strain_rollup", "madin_status", "scope")
    return {key: provenance[key] for key in keys}


def main() -> None:
    """Print the record's values as JSON."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-root", type=Path, required=True)
    parser.add_argument("--additive-root", type=Path, required=True)
    parser.add_argument(
        "--units",
        nargs="+",
        default=["tier12-off", "tier12-phenotype", "full-phenotype", "full-off"],
    )
    args = parser.parse_args()
    output = {
        "units": {unit: unit_values(args.runs_root / unit) for unit in args.units},
        "scan": scan_values(args.additive_root / "scan835"),
        "lineage": lineage_values(args.additive_root / "influence-extension"),
        "label_support": label_support_values(args.additive_root / "label-support"),
    }
    print(json.dumps(output, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
