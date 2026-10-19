"""Audit regenerated complete ORA exports for #720 (2026-10-19).

Inputs: a JSON manifest of historical/current TSV pairs and their explicit profile
and taxa-matrix inputs. Output: a JSON audit; nonzero exit on unexpected changes.
Only new input fingerprint metadata is allowed. Historical exports are read-only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections.abc import Iterator
from dataclasses import asdict, dataclass
from itertools import chain
from pathlib import Path

from mtp.enrichment_results import (
    META_RESULT_STATUS,
    META_SIGNIFICANT_HITS,
    META_TESTS_PERFORMED,
    RESULT_COLUMNS,
    classify_result_status,
)
from mtp.paths import file_sha256

FINGERPRINTS = {"Profiles-SHA256": "profiles", "Taxa-Matrix-SHA256": "taxa_matrix"}
PAIR_FIELDS = {"label", "historical", "current", "profiles", "taxa_matrix"}


def same_path(left: Path, right: Path) -> bool:
    """Include existing hard links when protecting input files from overwrite."""
    return left.resolve() == right.resolve() or (
        left.exists() and right.exists() and left.samefile(right)
    )


@dataclass(frozen=True)
class TableAudit:
    path: str
    sha256: str
    body_sha256: str
    body_bytes: int
    rows: int
    significant_rows: int
    metadata: dict[str, str]
    errors: list[str]


def read_table(path: Path) -> TableAudit:
    """Stream a complete export and reconcile its declared and observed counts."""
    metadata: dict[str, str] = {}
    body_digest = hashlib.sha256()
    body_bytes = 0
    rows = 0
    significant = 0
    errors: list[str] = []
    with path.open("rb") as handle:
        first_line = handle.readline()
        while first_line.startswith(b"#"):
            key, separator, value = first_line.decode("utf-8").strip()[1:].strip().partition(": ")
            if not separator or not key or key in metadata:
                raise ValueError(f"{path}: malformed or duplicate metadata: {first_line!r}")
            metadata[key] = value
            first_line = handle.readline()

        def lines() -> Iterator[str]:
            nonlocal body_bytes
            for line in chain((first_line,), handle):
                body_digest.update(line)
                body_bytes += len(line)
                yield line.decode("utf-8")

        reader = csv.reader(lines(), delimiter="\t", strict=True)
        header = next(reader, [])
        required = {*RESULT_COLUMNS, "is_significant"}
        if len(set(header)) != len(header) or not required.issubset(header):
            raise ValueError(f"{path}: missing or duplicate complete-result columns")
        selected_index = header.index("is_significant")
        for row in reader:
            rows += 1
            if len(row) != len(header):
                raise ValueError(f"{path}: malformed table row {rows}")
            selected = row[selected_index]
            if selected not in {"True", "False"}:
                raise ValueError(f"{path}: invalid is_significant value in row {rows}")
            significant += selected == "True"
    for key, observed in ((META_TESTS_PERFORMED, rows), (META_SIGNIFICANT_HITS, significant)):
        try:
            declared = int(metadata[key])
        except (KeyError, ValueError) as exc:
            raise ValueError(f"{path}: missing or invalid {key}") from exc
        if declared != observed:
            errors.append(f"{key}: declared {declared}, observed {observed}")
    status = classify_result_status(rows, significant)
    if metadata.get(META_RESULT_STATUS) != status:
        errors.append(
            f"{META_RESULT_STATUS}: declared {metadata.get(META_RESULT_STATUS)!r}, "
            f"observed {status!r}"
        )
    return TableAudit(
        path=str(path),
        sha256=file_sha256(path),
        body_sha256=body_digest.hexdigest(),
        body_bytes=body_bytes,
        rows=rows,
        significant_rows=significant,
        metadata=metadata,
        errors=errors,
    )


def audit_pair(pair: dict[str, str], base: Path) -> dict[str, object]:
    paths = {key: (base / value).resolve() for key, value in pair.items() if key != "label"}
    if same_path(paths["historical"], paths["current"]):
        raise ValueError(f"{pair['label']}: historical and current paths must differ")
    historical = read_table(paths["historical"])
    current = read_table(paths["current"])
    old_meta = historical.metadata
    new_meta = current.metadata
    inputs = {
        key: {"path": str(paths[key]), "sha256": file_sha256(paths[key])}
        for key in FINGERPRINTS.values()
    }
    changes = {
        key: {"historical": old_meta.get(key), "current": new_meta.get(key)}
        for key in sorted(old_meta.keys() | new_meta.keys())
        if old_meta.get(key) != new_meta.get(key)
    }
    errors = [
        f"{name}: {error}"
        for name, table in (("historical", historical), ("current", current))
        for error in table.errors
    ]
    body_equal = (
        historical.body_sha256 == current.body_sha256
        and historical.body_bytes == current.body_bytes
    )
    if not body_equal:
        errors.append("Table body bytes differ")
    for key in changes:
        if key not in FINGERPRINTS or key in old_meta:
            errors.append(f"Unexpected metadata change: {key}")
    for key, input_key in FINGERPRINTS.items():
        if new_meta.get(key) != inputs[input_key]["sha256"]:
            errors.append(f"{key}: current export does not match supplied input")
    return {
        "label": pair["label"],
        "historical": asdict(historical),
        "current": asdict(current),
        "inputs": inputs,
        "metadata_changes": changes,
        "table_body_equal": body_equal,
        "passed": not errors,
        "errors": errors,
    }


def audit_manifest(manifest: Path, output: Path) -> dict[str, object]:
    raw = json.loads(manifest.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or set(raw) != {"pairs"}:
        raise ValueError("Manifest must contain exactly the 'pairs' key")
    pairs = raw["pairs"]
    if not isinstance(pairs, list) or not pairs:
        raise ValueError("Manifest 'pairs' must be a nonempty list")
    labels: set[str] = set()
    for pair in pairs:
        if (
            not isinstance(pair, dict)
            or set(pair) != PAIR_FIELDS
            or any(not isinstance(value, str) or not value.strip() for value in pair.values())
        ):
            raise ValueError(f"Each pair must supply nonempty strings for {sorted(PAIR_FIELDS)}")
        if pair["label"] in labels:
            raise ValueError(f"Duplicate pair label: {pair['label']}")
        labels.add(pair["label"])
        for key in PAIR_FIELDS - {"label"}:
            if same_path(manifest.parent / pair[key], output):
                raise ValueError("Audit output must not overwrite any input")
    if same_path(manifest, output):
        raise ValueError("Audit output must not overwrite manifest")
    results = []
    for pair in pairs:
        try:
            results.append(audit_pair(pair, manifest.parent))
        except (OSError, ValueError, csv.Error) as exc:
            results.append({"label": pair["label"], "passed": False, "errors": [str(exc)]})
    return {
        "issue": 720,
        "manifest": str(manifest.resolve()),
        "manifest_sha256": file_sha256(manifest),
        "passed": all(result["passed"] for result in results),
        "pairs": results,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit_manifest(args.manifest, args.output)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"ORA refresh audit: {'PASS' if result['passed'] else 'FAIL'} ({args.output})")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
