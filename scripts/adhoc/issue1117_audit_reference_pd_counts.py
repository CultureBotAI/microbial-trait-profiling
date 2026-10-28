"""Audit archived PD pair support counts for #1117/#1165 (2026-10-28).

Inputs: an explicitly pinned reference consumer manifest and its colocated PD
summary/pair TSVs. Outputs: a count TSV and provenance JSON in a fresh durable
directory. Hash and parse the same bytes; do not reconstruct taxon sets, rerun
scientific tests, or claim producer/core re-verification.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from typing import Any

from mtp.workflow.output_roots import require_durable_output_path

SUMMARY = "parkinsons_stool_study_overlap_summary.tsv"
PAIRS = "parkinsons_stool_study_pairs.tsv"
COUNTS = ("n_studies", "n_study_pairs", "nonzero_pairs", "zero_pairs")


def identity(raw: bytes) -> dict[str, str | int]:
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def rows(raw: bytes) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8")), delimiter="\t")
    names = reader.fieldnames
    if not names or len(names) != len(set(names)):
        raise ValueError("Missing or duplicate TSV columns")
    result = list(reader)
    if any(
        set(row) != set(names) or any(value is None for value in row.values()) for row in result
    ):
        raise ValueError("Malformed TSV row")
    return result


def integer(row: dict[str, str], name: str) -> int:
    value = row[name]
    if not value or any(character not in "0123456789" for character in value):
        raise ValueError(f"Expected a nonnegative integer for {name}: {value!r}")
    return int(value)


def support_counts(
    summary: list[dict[str, str]], pairs: list[dict[str, str]]
) -> list[dict[str, str | int]]:
    scopes = [row["scope"] for row in summary]
    if not scopes or len(scopes) != len(set(scopes)) or any(not scope for scope in scopes):
        raise ValueError("Summary scope identities must be nonempty and unique")
    by_scope: dict[str, list[dict[str, str]]] = {scope: [] for scope in scopes}
    seen: set[tuple[str, str, str]] = set()
    for row in pairs:
        scope = row["scope"]
        if scope not in by_scope:
            raise ValueError("Pair ledger contains a scope absent from the summary")
        a, b = row["study_a_id"], row["study_b_id"]
        if not a or not b or a.strip() != a or b.strip() != b or a == b:
            raise ValueError("Pair ledger contains a blank, malformed or self-paired study key")
        pair = (scope, min(a, b), max(a, b))
        if pair in seen:
            raise ValueError("Pair ledger contains a duplicate unordered study pair")
        seen.add(pair)
        na, nb = integer(row, "study_a_taxa"), integer(row, "study_b_taxa")
        intersection, union = integer(row, "intersection_taxa"), integer(row, "union_taxa")
        if union <= 0 or intersection > min(na, nb) or union != na + nb - intersection:
            raise ValueError("Pair ledger contains inconsistent integer intersection/union counts")
        by_scope[scope].append(row)
    result: list[dict[str, str | int]] = []
    for row in summary:
        scope = row["scope"]
        selected = by_scope[scope]
        studies = {row[key] for row in selected for key in ("study_a_id", "study_b_id")}
        nonzero = sum(integer(pair, "intersection_taxa") > 0 for pair in selected)
        actual: dict[str, str | int] = {
            "scope": scope,
            "n_studies": len(studies),
            "n_study_pairs": len(selected),
            "nonzero_pairs": nonzero,
            "zero_pairs": len(selected) - nonzero,
        }
        if len(selected) != len(studies) * (len(studies) - 1) // 2:
            raise ValueError(f"Incomplete unordered study-pair family: {scope}")
        if any(actual[name] != integer(row, name) for name in COUNTS):
            raise ValueError(f"Integer pair ledger disagrees with archived summary: {scope}")
        result.append(actual)
    return result


def audit(reference: Path, expected_sha256: str, output_dir: Path) -> dict[str, Any]:
    reference = reference.resolve(strict=True)
    output_dir = require_durable_output_path(output_dir, allowed=False)
    source_root = reference.parent
    if output_dir.is_relative_to(source_root) or source_root.is_relative_to(output_dir):
        raise ValueError("Audit output must be disjoint from the reference directory")
    if output_dir.exists():
        raise ValueError("Audit output directory must be fresh")
    reference_bytes = reference.read_bytes()
    if identity(reference_bytes)["sha256"] != expected_sha256:
        raise ValueError("Reference consumer manifest differs from its explicit pin")
    record: dict[str, Any] = json.loads(reference_bytes)
    if record["consumer"] != "matrix-overlap-coverage":
        raise ValueError("Reference has the wrong consumer identity")
    snapshots = {reference: reference_bytes}
    inputs: dict[str, Any] = {
        "reference_consumer": {"path": str(reference), **identity(reference_bytes)},
    }
    parsed = {}
    for name in (SUMMARY, PAIRS):
        matches = [entry for entry in record["outputs"] if Path(entry["path"]).name == name]
        if len(matches) != 1:
            raise ValueError(f"Reference output basename is missing or ambiguous: {name}")
        entry = matches[0]
        path = source_root / name
        raw = path.read_bytes()
        if identity(raw) != {key: entry[key] for key in ("sha256", "size_bytes")}:
            raise ValueError(f"Reference output bytes differ from consumer inventory: {name}")
        snapshots[path] = raw
        parsed[name] = rows(raw)
        inputs[name] = {
            "recorded_path": entry["path"],
            "read_path": str(path),
            "path_policy": "explicit_basename_relocation_to_reference_manifest_directory",
            **identity(raw),
        }
    results = support_counts(parsed[SUMMARY], parsed[PAIRS])
    for name, path in {
        "audit_script": Path(__file__).resolve(),
        "durable_output_policy": Path(__file__).resolve().parents[2]
        / "src/mtp/workflow/output_roots.py",
    }.items():
        raw = path.read_bytes()
        snapshots[path] = raw
        inputs[name] = {"path": str(path), **identity(raw)}
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer, fieldnames=["scope", *COUNTS], delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(results)
    table_bytes = buffer.getvalue().encode("utf-8")
    output_identities = {"pd_support_counts.tsv": identity(table_bytes)}
    payload = {
        "issue": 1165,
        "parent_issue": 1117,
        "audit_kind": "archived_reference_integer_pair_support_only",
        "reference_consumer_expected_sha256": expected_sha256,
        "inputs": inputs,
        "checks": {
            "same_bytes_hashed_and_parsed": True,
            "consumer_bound_summary_and_pair_ledger": True,
            "unique_unordered_nonself_pairs": True,
            "positive_consistent_unions": True,
            "complete_pair_families": True,
            "summary_support_counts_match_pair_ledger": True,
        },
        "scopes": results,
        "producer_core_reverification_performed": False,
        "taxon_sets_reconstructed": False,
        "scientific_tests_run": False,
        "limitation": "Reference output identity and internal integer-count consistency only; "
        "the archived consumer predates separate sidecar producer binding. "
        "No independent clean-producer or current core-verification claim is made.",
        "outputs": output_identities,
    }
    for path, raw in snapshots.items():
        if path.read_bytes() != raw:
            raise ValueError(f"Input changed during the audit: {path}")
    output_dir.mkdir(parents=True, exist_ok=False)
    with (output_dir / "pd_support_counts.tsv").open("xb") as handle:
        handle.write(table_bytes)
    report_bytes = (json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    with (output_dir / "audit.json").open("xb") as handle:
        handle.write(report_bytes)
    return {
        "output_directory": str(output_dir),
        "outputs": {**output_identities, "audit.json": identity(report_bytes)},
        "scopes": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-consumer", type=Path, required=True)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.reference_consumer, args.reference_sha256, args.output_dir)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
