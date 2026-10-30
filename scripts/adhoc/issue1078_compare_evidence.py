"""Compare pinned table bytes without changing either input (#1078; 2026-10-26).

Inputs: accepted and candidate per-trait TSVs. Output: JSON with SHA256 pins,
changed columns, counts and exact example values for the qualification review.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from collections import Counter
from pathlib import Path


def read_table(path: Path) -> tuple[list[str], dict[str, dict[str, str]], str]:
    raw = path.read_bytes()
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8")), delimiter="\t")
    fields = list(reader.fieldnames or [])
    rows = list(reader)
    keyed = {row["trait"]: row for row in rows}
    if not fields or len(keyed) != len(rows):
        raise ValueError("Missing columns or duplicate traits")
    return fields, keyed, hashlib.sha256(raw).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("accepted", type=Path)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()
    fields, before, accepted_sha = read_table(args.accepted)
    columns, after, candidate_sha = read_table(args.candidate)
    if fields != columns or before.keys() != after.keys():
        raise ValueError("Table columns or trait identities differ")
    changes: Counter[str] = Counter()
    examples: dict[str, list[dict[str, str | float | bool]]] = {}
    for trait, original in before.items():
        for column in fields:
            old, new = original[column], after[trait][column]
            if old == new:
                continue
            changes[column] += 1
            example: dict[str, str | float | bool] = {
                "trait": trait,
                "accepted": old,
                "candidate": new,
            }
            try:
                left, right = float(old), float(new)
            except ValueError:
                example["numeric"] = False
            else:
                if not math.isfinite(left) or not math.isfinite(right):
                    raise ValueError("Changed numerical value is nonfinite")
                example["numeric"] = True
                example["absolute_difference"] = abs(right - left)
                example["within_rel_1e_12_abs_1e_15"] = math.isclose(
                    left, right, rel_tol=1e-12, abs_tol=1e-15
                )
            if len(examples.setdefault(column, [])) < 8:
                examples[column].append(example)
    print(
        json.dumps(
            {
                "accepted_sha256": accepted_sha,
                "candidate_sha256": candidate_sha,
                "traits": len(before),
                "changed_columns": dict(changes),
                "examples": examples,
            },
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
