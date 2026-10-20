"""Audit retained assessment refresh outputs for #729/#730/#734 (2026-10-19).

Inputs: before/after calibration and integrated-output directories. Output: JSON
identities, row comparisons and explicit allowed assessment-column changes.
Scientific integrated fields and non-guild calibration results must stay exact;
this comparison does not independently validate any scientific procedure.
Assessment-source compatibility is enforced by integration; consumer-record
verification and provenance/assessment sidecar inspection remain separate.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from mtp.paths import file_sha256

ASSESSMENT_COLUMNS = frozenset(
    {
        "calibration_status",
        "sensitivity_status",
        "coverage_status",
        "orientation_status",
        "taxon_influence_status",
    }
)
# This refresh changes calibration evidence only. Other assessments must retain
# their recorded status; binding metadata is recorded in a separate sidecar.
ALLOWED_CHANGED_COLUMNS = frozenset({"calibration_status"})
CALIBRATION_TABLES = (
    "calibration_results.tsv",
    "calibration_replicates.tsv",
    "calibration_criteria.tsv",
)
INTEGRATION_TABLES = ("integrated_evidence.tsv", "integrated_evidence_method_results.tsv")
GUILD_IDENTITY = {
    "calibration_results.tsv": ("scenario", "procedure", "label", "role"),
    "calibration_replicates.tsv": ("scenario", "procedure", "replicate"),
    "calibration_criteria.tsv": ("criterion", "procedure"),
}
GUILD_FIXED_POLICY = {
    "calibration_results.tsv": ("truth", "replicates"),
    "calibration_replicates.tsv": (),
    "calibration_criteria.tsv": (
        "requirement",
        "tolerance",
        "replicates_required",
        "role",
        "expected_verdict",
    ),
}


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        if not header or len(set(header)) != len(header):
            raise ValueError(f"Missing or repeated columns: {path}")
        rows = [dict(zip(header, row, strict=True)) for row in reader]
    if not rows:
        raise ValueError(f"Empty table: {path}")
    return header, rows


def compare_table(
    before: Path,
    after: Path,
    *,
    assessment_columns: bool = False,
    exclude_guild: bool = False,
) -> dict[str, object]:
    header, left = read_table(before)
    other_header, right = read_table(after)
    if header != other_header:
        raise ValueError(f"Column identity/order differs: {before} vs {after}")
    ignored = ALLOWED_CHANGED_COLUMNS if assessment_columns else frozenset()
    if assessment_columns and not ASSESSMENT_COLUMNS.issubset(header):
        raise ValueError("Integrated table lacks required assessment columns")
    total_rows = {"total_rows_before": len(left), "total_rows_after": len(right)}
    guild_identity: dict[str, object] = {}
    if exclude_guild:
        identity = GUILD_IDENTITY[before.name]
        policy = GUILD_FIXED_POLICY[before.name]
        if not {*identity, *policy}.issubset(header):
            raise ValueError("Calibration table lacks guild identity or fixed-policy columns")
        old_guild = [row for row in left if row["procedure"] == "guild_mcnemar"]
        new_guild = [row for row in right if row["procedure"] == "guild_mcnemar"]
        old_keys = [tuple(row[name] for name in identity) for row in old_guild]
        new_keys = [tuple(row[name] for name in identity) for row in new_guild]
        guild_identity = {
            "guild_rows_before": len(old_keys),
            "guild_rows_after": len(new_keys),
            "guild_identity_columns": list(identity),
            "guild_fixed_policy_columns": list(policy),
            "guild_fixed_policy_equal": (
                [tuple(row[name] for name in policy) for row in old_guild]
                == [tuple(row[name] for name in policy) for row in new_guild]
            ),
            "guild_row_identity_equal": (
                bool(old_keys) and old_keys == new_keys and len(set(old_keys)) == len(old_keys)
            ),
        }
        left = [row for row in left if row["procedure"] != "guild_mcnemar"]
        right = [row for row in right if row["procedure"] != "guild_mcnemar"]
        if not left or not right:
            raise ValueError("Calibration comparison has no non-guild rows")
    changed = Counter(
        name
        for old, new in zip(left, right, strict=False)
        for name in header
        if old[name] != new[name]
    )
    scientific = set(header) - ignored
    return {
        "before": str(before.resolve()),
        "after": str(after.resolve()),
        "before_sha256": file_sha256(before),
        "after_sha256": file_sha256(after),
        "rows_before": len(left),
        "rows_after": len(right),
        **total_rows,
        **guild_identity,
        "rows_equal_outside_assessments": (
            len(left) == len(right) and not scientific.intersection(changed)
        ),
        "changed_columns": dict(sorted(changed.items())),
        "allowed_assessment_columns": sorted(ignored),
        "guild_rows_excluded": exclude_guild,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "before-calibration",
        "after-calibration",
        "before-integration",
        "after-integration",
        "output",
    ):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        parser.error("Audit output must be a new file")
    pairs = [
        (args.before_calibration / name, args.after_calibration / name)
        for name in CALIBRATION_TABLES
    ] + [
        (args.before_integration / name, args.after_integration / name)
        for name in INTEGRATION_TABLES
    ]
    unchanged = [
        (args.before_calibration / name, args.after_calibration / name)
        for name in ("calibration_flip_assignments.tsv", "calibration_injection.tsv")
    ] + [
        (args.before_integration / name, args.after_integration / name)
        for name in ("integrated_evidence_features.tsv", "integrated_evidence_families.tsv")
    ]
    for pair in (*pairs, *unchanged):
        for path in pair:
            if args.output.resolve() == path.resolve() or (
                args.output.exists() and path.samefile(args.output)
            ):
                parser.error("Audit output must not overwrite an input")
    comparisons = [
        compare_table(
            before,
            after,
            exclude_guild=before.name in CALIBRATION_TABLES,
            assessment_columns=before.name in INTEGRATION_TABLES,
        )
        for before, after in pairs
    ]
    identical = [
        {
            "before": str(before.resolve()),
            "after": str(after.resolve()),
            "before_sha256": file_sha256(before),
            "after_sha256": file_sha256(after),
            "byte_identical": file_sha256(before) == file_sha256(after),
        }
        for before, after in unchanged
    ]
    passed = all(
        row["rows_equal_outside_assessments"]
        and row.get("guild_row_identity_equal", True)
        and row.get("guild_fixed_policy_equal", True)
        for row in comparisons
    ) and all(row["byte_identical"] for row in identical)
    result = {
        "status": "ok" if passed else "failed",
        "scope": "Exact scientific table comparison, allowing only listed assessment-column "
        "updates and allowing changed guild statistical outcomes while preserving "
        "guild row identity, declared replicate counts and fixed policy. "
        "Not independent calibration or biological validation. Assessment-source binding, "
        "consumer verification and provenance/assessment sidecar inspection are separate checks.",
        "audit_source_sha256": file_sha256(Path(__file__)),
        "table_comparisons": comparisons,
        "unchanged_artifacts": identical,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
