"""Compare core-bound strata with the preserved direct-path outputs (#727, 2026-10-19).

Inputs: historical and regenerated stratification directories. Output: JSON to
stdout with exact file comparisons and complete/disjoint partition accounting.
Only the output-directory references in the text sensitivity report may differ.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from mtp.paths import file_sha256
from mtp.stratification import PRIMARY_STRATA


def compare(before: Path, after: Path) -> dict[str, object]:
    """Fail on changed scientific products or incomplete/overlapping partitions."""
    reserved = {"core_consumption.json", "core_consumption.coverage.tsv"}
    old_files = {path.name for path in before.iterdir() if path.is_file()}
    new_files = {path.name for path in after.iterdir() if path.is_file()} - reserved
    if len(old_files) != 24 or old_files != new_files:
        raise ValueError("Expected the same 24 stratification products in both directories")
    comparisons = []
    for name in sorted(old_files):
        old, new = before / name, after / name
        old_hash, new_hash = file_sha256(old), file_sha256(new)
        equal = old_hash == new_hash
        if not equal and name == "stratification_sensitivity_summary.txt":
            equal = old.read_text().replace(str(before), "<OUTPUT>") == new.read_text().replace(
                str(after), "<OUTPUT>"
            )
        if not equal:
            raise ValueError(f"Scientific output differs: {name}")
        comparisons.append(
            {
                "file": name,
                "historical_sha256": old_hash,
                "current_sha256": new_hash,
                "byte_identical": old_hash == new_hash,
                "scientifically_identical": equal,
            }
        )
    memberships: dict[str, set[tuple[str, str, str]]] = {}
    for stratum in PRIMARY_STRATA:
        with (after / f"study_trait_profiles_{stratum}.tsv").open() as handle:
            rows = csv.DictReader(handle, delimiter="\t")
            keys = [(row["study_id"], row["experiment_id"], row["inc_dec"]) for row in rows]
        if len(keys) != len(set(keys)):
            raise ValueError(f"Duplicate profile keys in {stratum}")
        memberships[stratum] = set(keys)
    feces, other = memberships["feces"], memberships["feces_excluded"]
    if feces & other or feces | other != memberships["pooled"]:
        raise ValueError("Feces partitions must be disjoint and cover the pooled profiles")
    return {
        "issue": 727,
        "historical": str(before),
        "current": str(after),
        "passed": True,
        "partition_profiles": {key: len(value) for key, value in memberships.items()},
        "complete_disjoint_partitions": True,
        "products": comparisons,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("historical", type=Path)
    parser.add_argument("current", type=Path)
    args = parser.parse_args()
    print(json.dumps(compare(args.historical.resolve(), args.current.resolve()), indent=2))


if __name__ == "__main__":
    main()
