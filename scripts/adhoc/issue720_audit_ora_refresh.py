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
import math
from collections.abc import Iterator
from dataclasses import asdict, dataclass
from itertools import chain
from pathlib import Path

from mtp.enrichment_results import (
    COUNT_COLUMNS,
    FLOAT_COLUMNS,
    ID_COLUMNS,
    META_RESULT_STATUS,
    META_SIGNIFICANT_HITS,
    META_TESTS_PERFORMED,
    PROBABILITY_COLUMNS,
    RESULT_COLUMNS,
    classify_result_status,
)
from mtp.paths import file_sha256

FINGERPRINTS = {"Profiles-SHA256": "profiles", "Taxa-Matrix-SHA256": "taxa_matrix"}
PAIR_FIELDS = {"label", "historical", "current", "profiles", "taxa_matrix"}
BACKGROUND_COLUMNS = ("background_m", "background_k")
COMPLETE_COLUMNS = {
    *RESULT_COLUMNS,
    *BACKGROUND_COLUMNS,
    "is_significant",
    "correction_family",
    "eligibility_status",
}


def number(
    value: str,
    field: str,
    *,
    integer: bool = False,
    probability: bool = False,
    positive_infinity: bool = False,
) -> float:
    """Apply the shared numeric domains without materializing an ORA table (#724)."""
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"invalid {field}: {value!r}") from exc
    if (
        (not math.isfinite(parsed) and not (positive_infinity and parsed == math.inf))
        or parsed < 0
        or (integer and not parsed.is_integer())
        or (probability and parsed > 1)
    ):
        raise ValueError(f"invalid {field}: {value!r}")
    return parsed


@dataclass(frozen=True)
class RowContract:
    """Validate each complete row against the producer's declared analysis settings."""

    indices: dict[str, int]
    upper_tail: bool
    fdr: bool
    family: str
    threshold: float
    fold_cutoff: float
    background_m: float
    testable_profiles: float

    def selected(self, row: list[str]) -> bool:
        indices = self.indices
        for column in ID_COLUMNS:
            if not row[indices[column]].strip():
                raise ValueError(f"empty {column}")
        if row[indices["inc_dec"]] not in {"inc", "dec"}:
            raise ValueError("invalid inc_dec")
        values = {}
        for column in PROBABILITY_COLUMNS + FLOAT_COLUMNS + COUNT_COLUMNS + BACKGROUND_COLUMNS:
            values[column] = number(
                row[indices[column]],
                column,
                integer=column in COUNT_COLUMNS + BACKGROUND_COLUMNS,
                probability=column in PROBABILITY_COLUMNS,
            )
        if (
            values["background_m"] != self.background_m
            or values["testable_profiles"] != self.testable_profiles
        ):
            raise ValueError("row background_m or testable_profiles differs from metadata")
        if not 0 < values["background_k"] <= values["background_m"]:
            raise ValueError("background_k must be positive and not exceed background_m")
        if values["count"] > values["background_k"]:
            raise ValueError("count exceeds background_k")
        if values["expected"] <= 0 or values["testable_profiles"] <= 0:
            raise ValueError("testable rows require positive expected and testable_profiles")
        if row[indices["eligibility_status"]] != "testable":
            raise ValueError("invalid eligibility_status")
        if row[indices["correction_family"]] != self.family:
            raise ValueError("correction_family differs from FDR-Scope")
        if not self.fdr and values["fdr_q"] != values["p_value"]:
            raise ValueError("fdr_q differs from p_value with FDR disabled")
        selected = row[indices["is_significant"]]
        if selected not in {"True", "False"}:
            raise ValueError("invalid is_significant")
        statistic = values["fdr_q" if self.fdr else "p_value"]
        fold_passes = (
            values["fold"] >= self.fold_cutoff
            if self.upper_tail
            else values["fold"] <= self.fold_cutoff
        )
        if (selected == "True") != (statistic <= self.threshold and fold_passes):
            raise ValueError("is_significant disagrees with declared thresholds")
        return selected == "True"


def row_contract(header: list[str], metadata: dict[str, str]) -> RowContract:
    """Read settings once, accepting both tails, correction scopes and FDR modes."""
    for key, allowed in (
        ("Tail", {"upper", "lower"}),
        ("FDR", {"True", "False"}),
        ("FDR-Scope", {"global", "profile"}),
        ("Testable-Set", {"K>0 & N>0"}),
    ):
        if metadata.get(key) not in allowed:
            raise ValueError(f"missing or invalid {key}")
    fold_key = "Min-Fold" if metadata["Tail"] == "upper" else "Max-Fold"
    try:
        total_profiles = number(
            metadata["Total-Profiles-Input"], "Total-Profiles-Input", integer=True
        )
        testable_profiles = number(metadata["Testable-Profiles"], "Testable-Profiles", integer=True)
        if testable_profiles > total_profiles:
            raise ValueError("Testable-Profiles exceeds Total-Profiles-Input")
        return RowContract(
            indices={column: index for index, column in enumerate(header)},
            upper_tail=metadata["Tail"] == "upper",
            fdr=metadata["FDR"] == "True",
            family=metadata["FDR-Scope"],
            threshold=number(metadata["P-Threshold"], "P-Threshold", probability=True),
            fold_cutoff=number(metadata[fold_key], fold_key, positive_infinity=True),
            background_m=number(
                metadata["Background-Taxa-Universe-M"], "Background-Taxa-Universe-M", integer=True
            ),
            testable_profiles=testable_profiles,
        )
    except KeyError as exc:
        raise ValueError(f"missing metadata: {exc.args[0]}") from exc


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
        if len(set(header)) != len(header) or not COMPLETE_COLUMNS.issubset(header):
            raise ValueError(f"{path}: missing or duplicate complete-result columns")
        try:
            contract = row_contract(header, metadata)
        except ValueError as exc:
            raise ValueError(f"{path}: {exc}") from exc
        for row in reader:
            rows += 1
            if len(row) != len(header):
                raise ValueError(f"{path}: malformed table row {rows}")
            try:
                significant += contract.selected(row)
            except ValueError as exc:
                raise ValueError(f"{path}: row {rows}: {exc}") from exc
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
