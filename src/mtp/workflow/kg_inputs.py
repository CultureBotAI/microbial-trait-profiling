"""Versioned KG edge serialization contract, preserving native count semantics (#629).

The KG export can contain an unquoted line break in its header and in the
second agent_type field. These exact producer encodings are reconstructed,
audited, and tested against the native first-four-column counting parser.
Unknown malformed records fail; no taxon or trait row is silently skipped.
"""

from __future__ import annotations

import csv
from collections.abc import Iterator
from pathlib import Path

KG_EDGE_CONTRACT = "kg-edge-export-v1"
_SPLIT_HEADER = ["id", "subject", "predicate", "object", "relation", "agent_type", "agent_type"]
_HEADER_CONTINUATION = [
    "",
    "knowledge_level",
    "knowledge_source",
    "meta",
    "primary_knowledge_source",
]


def iter_kg_edges(path: Path, audit: dict[str, object] | None = None) -> Iterator[dict[str, str]]:
    """Read complete evidence records and record every recognized serialization repair.

    Both duplicate-named agent fields are preserved as agent_type and
    agent_type_secondary. The four columns used by native numerical aggregation
    are unchanged. Physical source line numbers remain available in each record.
    """
    facts = audit if audit is not None else {}
    records = continuations = header_repairs = duplicate_fields = 0
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if header[:4] != ["id", "subject", "predicate", "object"]:
            raise ValueError(
                f"{path}: expected leading KG edge columns id/subject/predicate/object"
            )
        if header == _SPLIT_HEADER:
            if next(reader, []) != _HEADER_CONTINUATION:
                raise ValueError(f"{path}: invalid known KG header continuation")
            header = [*header[:6], "agent_type_secondary", *_HEADER_CONTINUATION[1:]]
            header_repairs = 1
            duplicate_fields = 1
        elif header == [*_SPLIT_HEADER, *_HEADER_CONTINUATION[1:]]:
            header[6] = "agent_type_secondary"
            duplicate_fields = 1
        elif len(set(header)) != len(header) or any(not name for name in header):
            raise ValueError(f"{path}: duplicate or missing KG header names")
        previous_end = reader.line_num
        for row in reader:
            start_line = previous_end + 1
            if (
                len(row) == 7
                and len(header) == 11
                and row[6] == "manual_validation_of_automated_agent"
            ):
                continuation = next(reader, [])
                if len(continuation) != 5 or continuation[0] != "":
                    raise ValueError(f"{path}:{start_line}: invalid known agent-field continuation")
                row += continuation[1:]
                continuations += 1
            if len(row) != len(header):
                raise ValueError(f"{path}:{start_line}: unknown malformed KG edge record")
            if any(not field.strip() for field in row[:4]):
                raise ValueError(
                    f"{path}:{start_line}: empty KG edge identity/subject/predicate/object"
                )
            record = dict(zip(header, row, strict=True))
            record["source_physical_row_start"] = str(start_line)
            record["source_physical_row_end"] = str(reader.line_num)
            previous_end = reader.line_num
            records += 1
            yield record
    facts.update(
        {
            "parser_contract": KG_EDGE_CONTRACT,
            "records": records,
            "header_continuations_reconstructed": header_repairs,
            "duplicate_agent_fields_preserved": duplicate_fields,
            "agent_field_continuations_reconstructed": continuations,
            "excluded_records": 0,
            "native_count_columns": ["subject", "predicate", "object"],
        }
    )


def audit_kg_edges(path: Path) -> dict[str, object]:
    """Validate a whole raw edge export without retaining its graph in memory."""
    facts: dict[str, object] = {}
    for _ in iter_kg_edges(path, facts):
        pass
    if not facts["records"]:
        raise ValueError(f"{path}: empty KG edge export")
    return facts
