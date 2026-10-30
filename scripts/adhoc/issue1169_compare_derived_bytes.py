"""Compare eight completed derived TSVs without replay (#1169, 2026-10-29).

Inputs: baseline/current manifest paths and SHA256 pins, expected current commit.
Output: JSON on stdout only; no input or report file is written. This checks
declared contracts and exact output bytes, not native upstream verification.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

TABLES = frozenset(
    {
        "composition_pair_delta_z.tsv",
        "composition_study_delta_z.tsv",
        "trait_summary.tsv",
        "quadrant_pairs.tsv",
        "quadrant_summary.tsv",
        "stool_disease_shifts.tsv",
        "stool_study_sensitivity.tsv",
        "site_disease_shifts.tsv",
    }
)
CONTRACT_FIELDS = (
    "version",
    "summary_only",
    "parent_run_id",
    "parent_manifest_sha256",
    "parent_producer",
    "core_run_id",
    "comparison_run_id",
    "source_run_ids",
    "source_commits",
    "strain_rollup",
    "scope_role",
    "contrast_population",
    "direction_basis",
    "analysis_bases",
    "family_size",
    "retained_traits",
    "counts",
    "trait_list",
    "policies",
    "parent_input_inventory",
    "inputs",
)


def identity(raw: bytes) -> dict[str, str | int]:
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate manifest field: {key}")
        result[key] = value
    return result


def read_bound(path: Path, snapshots: dict[Path, tuple[bytes, Path]]) -> bytes:
    target = path.resolve(strict=True)
    if not path.is_file() or path.is_symlink():
        raise ValueError("Bundle inputs must be regular nonsymlink files")
    raw = path.read_bytes()
    if path.resolve(strict=True) != target:
        raise ValueError("Input target changed while reading")
    if path in snapshots and snapshots[path] != (raw, target):
        raise ValueError("Input bytes or target changed during parity check")
    snapshots[path] = raw, target
    return raw


def load_bundle(
    path: Path,
    pin: str,
    snapshots: dict[Path, tuple[bytes, Path]],
) -> tuple[dict[str, Any], dict[str, bytes]]:
    raw = read_bound(path, snapshots)
    if re.fullmatch(r"[0-9a-f]{64}", pin) is None or identity(raw)["sha256"] != pin:
        raise ValueError(f"Explicit manifest pin differs: {path}")
    record = json.loads(raw, object_pairs_hook=unique_object)
    if not isinstance(record, dict) or record.get("version") != "delta-z-derived-v1":
        raise ValueError("Expected a completed derived v1 manifest")
    if record.get("summary_only") is not False:
        raise ValueError("Expected a complete non-summary derived bundle")
    if not isinstance(record.get("outputs"), dict) or set(record["outputs"]) != TABLES:
        raise ValueError("Declared inventory must contain exactly the eight derived TSVs")
    if {entry.name for entry in path.parent.iterdir()} != {path.name, *TABLES}:
        raise ValueError("Physical bundle inventory differs from the declared complete bundle")
    tables: dict[str, bytes] = {}
    for name in sorted(TABLES):
        content = read_bound(path.parent / name, snapshots)
        expected = record["outputs"][name]
        if (
            not isinstance(expected, dict)
            or type(expected.get("size_bytes")) is not int
            or expected != identity(content)
        ):
            raise ValueError(f"Bound table identity differs: {name}")
        tables[name] = content
    return record, tables


def compare_bundles(
    baseline: Path,
    current: Path,
    *,
    baseline_sha256: str,
    current_sha256: str,
    expected_current_commit: str,
) -> dict[str, Any]:
    snapshots: dict[Path, tuple[bytes, Path]] = {}
    old, old_tables = load_bundle(baseline, baseline_sha256, snapshots)
    new, new_tables = load_bundle(current, current_sha256, snapshots)
    if (
        re.fullmatch(r"[0-9a-f]{40}", expected_current_commit) is None
        or new.get("code_commit") != expected_current_commit
        or new.get("code_dirty") is not False
    ):
        raise ValueError("Current producer must declare the expected clean commit")
    for name in CONTRACT_FIELDS:
        if name not in old or name not in new or old[name] != new[name]:
            raise ValueError(f"Analysis contract differs or is absent: {name}")
    for name in sorted(TABLES):
        if old_tables[name] != new_tables[name]:
            raise ValueError(f"Derived table bytes differ: {name}")
    for path in snapshots:
        read_bound(path, snapshots)
    for path in (baseline, current):
        if {entry.name for entry in path.parent.iterdir()} != {path.name, *TABLES}:
            raise ValueError("Bundle inventory changed during parity check")
    return {
        "version": "derived-byte-parity-v1",
        "status": "PASS",
        "native_verification_performed": False,
        "verification_scope": "Pinned output bytes and declared contracts; no upstream replay",
        "baseline": {"path": str(baseline), **identity(snapshots[baseline][0]), "record": old},
        "current": {"path": str(current), **identity(snapshots[current][0]), "record": new},
        "expected_current_commit": expected_current_commit,
        "contract_fields": list(CONTRACT_FIELDS),
        "tables": {name: identity(old_tables[name]) for name in sorted(TABLES)},
        "helper": {"path": str(Path(__file__).resolve()), **identity(Path(__file__).read_bytes())},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for role in ("baseline", "current"):
        parser.add_argument(f"--{role}-manifest", type=Path, required=True)
        parser.add_argument(f"--{role}-sha256", required=True)
    parser.add_argument("--expected-current-commit", required=True)
    args = parser.parse_args()
    result = compare_bundles(
        args.baseline_manifest,
        args.current_manifest,
        baseline_sha256=args.baseline_sha256,
        current_sha256=args.current_sha256,
        expected_current_commit=args.expected_current_commit,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
