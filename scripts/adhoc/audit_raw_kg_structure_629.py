#!/usr/bin/env python3
"""Audit physical KG TSV records before the #629 core workflow preflight.

Inputs: KG edge/node TSV paths. Output: JSON row-width counts and bounded malformed
examples on stdout. Serves issue #629; 2026-10-19. Does not modify input files.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

from mtp.paths import file_sha256
from mtp.workflow.kg_inputs import audit_kg_edges, iter_kg_edges


def native_core_stream_parity(path: Path) -> dict[str, object]:
    """Compare every native-consumed taxonomy/trait tuple, retaining order and duplicates."""

    def relevant(subject: str, predicate: str, obj: str) -> bool:
        return subject.startswith("NCBITaxon:") and (
            (predicate == "biolink:subclass_of" and obj.startswith("NCBITaxon:"))
            or predicate.startswith("METPO:")
            or obj.startswith("METPO:")
        )

    native_hash, canonical_hash = hashlib.sha256(), hashlib.sha256()
    native_count = canonical_count = 0
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) < 4:
                continue
            triple = tuple(field.strip() for field in row[1:4])
            if relevant(*triple):
                native_hash.update((json.dumps(triple, separators=(",", ":")) + "\n").encode())
                native_count += 1
    for row in iter_kg_edges(path):
        triple = tuple(row[name].strip() for name in ("subject", "predicate", "object"))
        if relevant(*triple):
            canonical_hash.update((json.dumps(triple, separators=(",", ":")) + "\n").encode())
            canonical_count += 1
    equal = (
        native_count == canonical_count and native_hash.hexdigest() == canonical_hash.hexdigest()
    )
    if not equal:
        raise ValueError("Canonical KG parser changed native-consumed taxonomy/trait tuples")
    return {
        "native_core_edge_records": native_count,
        "canonical_core_edge_records": canonical_count,
        "native_core_tuple_sha256": native_hash.hexdigest(),
        "canonical_core_tuple_sha256": canonical_hash.hexdigest(),
        "identical_native_numeric_input_stream": equal,
    }


def audit(path: Path) -> dict[str, object]:
    widths: Counter[int] = Counter()
    malformed_ids: Counter[str] = Counter()
    examples: list[dict[str, object]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        for row in reader:
            widths[len(row)] += 1
            if len(row) != len(header):
                head = row[0] if row else ""
                category = (
                    "identifier_prefix" if re.match(r"^[A-Za-z][A-Za-z0-9_.-]*:", head) else "other"
                )
                malformed_ids[category] += 1
                if len(examples) < 12:
                    examples.append({"line": reader.line_num, "width": len(row), "head": row[:4]})
    result: dict[str, object] = {
        "path": str(path.resolve()),
        "sha256": file_sha256(path),
        "header": header,
        "width_counts": dict(sorted(widths.items())),
        "malformed_first_field_categories": dict(malformed_ids),
        "examples": examples,
    }
    if header[:4] == ["id", "subject", "predicate", "object"]:
        result["canonical_contract"] = audit_kg_edges(path)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", type=Path, nargs="+")
    parser.add_argument("--verify-native-core-stream", action="store_true")
    args = parser.parse_args()
    results = []
    for path in args.paths:
        result = audit(path)
        if args.verify_native_core_stream and "canonical_contract" in result:
            result["numeric_input_parity"] = native_core_stream_parity(path)
        results.append(result)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
