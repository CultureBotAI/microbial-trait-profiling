"""Prepare independent historical policy inventory fixtures (#1129, 2026-10-27).

Inputs: explicit raw-policy and implicit-union Git references in this repository.
Output: deterministic JSON containing full source commits and actual matrix blob
SHA-256/size. Tests read this committed fixture without requiring Git history.
No historical code is executed and no network access is used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from git import Repo

MATRIX_SOURCE = "src/mtp/taxa_trait_matrix.py"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-commit", required=True)
    parser.add_argument("--union-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = []
    with Repo(Path(__file__).resolve().parents[2]) as repo:
        for reference, policy in (
            (args.raw_commit, "raw_object_v1"),
            (args.union_commit, "facultative_union_v1"),
        ):
            commit = repo.commit(reference)
            payload = (commit.tree / MATRIX_SOURCE).data_stream.read()
            records.append(
                {
                    "source_reference": reference,
                    "source_commit": commit.hexsha,
                    "policy": policy,
                    "code_inventory": {
                        "path": MATRIX_SOURCE,
                        "sha256": hashlib.sha256(payload).hexdigest(),
                        "size_bytes": len(payload),
                    },
                }
            )
    result = {"format": "historical-trait-policy-inventories-v1", "records": records}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, sort_keys=True, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
