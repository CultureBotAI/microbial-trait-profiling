"""Disk-backed membership checks and one full-family numeric correction (#853)."""

from __future__ import annotations

import hashlib
import json
import sqlite3
import tempfile
from collections.abc import Iterable, Iterator
from contextlib import closing
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

from mtp.datamodel.multiple_testing_pydantic import (
    AdjustedHypothesis,
    ExternalFamilyCorrection,
    HypothesisDefinition,
    HypothesisFamilyMetadata,
    RawHypothesis,
    SourceArtifact,
)
from mtp.multiple_testing import (
    adjusted_hypothesis,
    adjusted_probabilities,
    canonical_sha256,
    load_policy,
    scientific_hypothesis_sha256,
    validate_family_metadata,
    validate_hypothesis_definition,
    validate_raw_hypothesis,
)
from mtp.paths import file_sha256


def json_text(value: Any) -> str:
    """Canonical record encoding shared with the in-memory correction engine."""
    return json.dumps(
        value, sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(",", ":")
    )


def _array_hash(records: Iterable[str]) -> str:
    digest = hashlib.sha256()
    digest.update(b"[")
    for index, record in enumerate(records):
        if index:
            digest.update(b",")
        digest.update(record.encode("utf-8"))
    digest.update(b"]")
    return digest.hexdigest()


def _family_hash(metadata: HypothesisFamilyMetadata, db: sqlite3.Connection) -> str:
    values = metadata.model_dump(mode="json")
    digest = hashlib.sha256()
    digest.update(b"{")
    for index, key in enumerate(sorted([*values, "members"])):
        if index:
            digest.update(b",")
        digest.update((json_text(key) + ":").encode("utf-8"))
        if key == "members":
            digest.update(b"[")
            for member_index, (declaration,) in enumerate(
                db.execute("SELECT declaration FROM tests ORDER BY id")
            ):
                if member_index:
                    digest.update(b",")
                digest.update(declaration.encode("utf-8"))
            digest.update(b"]")
        else:
            digest.update(json_text(values[key]).encode("utf-8"))
    digest.update(b"}")
    return digest.hexdigest()


def _write_lines(path: Path, records: Iterable[str]) -> SourceArtifact:
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(record + "\n")
    return SourceArtifact(identity=path.name, sha256=file_sha256(path))


def read_corrected_inventory(path: Path) -> Iterator[AdjustedHypothesis]:
    """Read one validated corrected record at a time, without materializing a table."""
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            row = AdjustedHypothesis.model_validate_json(line)
            validate_raw_hypothesis(
                RawHypothesis.model_validate(
                    {name: getattr(row, name) for name in RawHypothesis.model_fields}
                )
            )
            yield row


def correct_inventory(
    metadata: HypothesisFamilyMetadata,
    declarations: Iterable[HypothesisDefinition],
    observations: Iterable[RawHypothesis],
    output: Path,
    *,
    recorded_at: datetime | None = None,
) -> ExternalFamilyCorrection:
    """Validate membership on disk and correct one global family, never chunks.

    SQLite enforces unique coordinates and exact expected/observed membership.
    Only the three numeric p/BY/BH vectors are retained in memory. JSON Lines
    inventories and their canonical identities match the small-family contract.
    """
    metadata = validate_family_metadata(metadata)
    policy = load_policy()
    instant = recorded_at if recorded_at is not None else datetime.now(UTC)
    offset = instant.utcoffset()
    if offset is None or offset.total_seconds() != 0:
        raise ValueError("recorded_at must be timezone-aware UTC")
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise FileExistsError("correction output directory must be empty")
    with (
        tempfile.TemporaryDirectory(prefix="mtp-multiplicity-") as scratch,
        closing(sqlite3.connect(Path(scratch) / "inventory.sqlite")) as db,
        db,
    ):
        db.execute(
            "CREATE TABLE tests (id TEXT PRIMARY KEY, scientific_id TEXT NOT NULL UNIQUE, "
            "declaration TEXT NOT NULL, raw TEXT, p REAL)"
        )
        for declared in declarations:
            member = validate_hypothesis_definition(declared)
            try:
                db.execute(
                    "INSERT INTO tests(id, scientific_id, declaration) VALUES (?, ?, ?)",
                    (
                        member.hypothesis_id,
                        scientific_hypothesis_sha256(member),
                        json_text(member.model_dump(mode="json")),
                    ),
                )
            except sqlite3.IntegrityError as error:
                raise ValueError("duplicate scientific hypothesis or declared identity") from error
        count = db.execute("SELECT COUNT(*) FROM tests").fetchone()[0]
        if count != metadata.expected_family_size:
            raise ValueError("expected family size differs from declared inventory")
        for observed in observations:
            row = validate_raw_hypothesis(observed)
            cursor = db.execute(
                "UPDATE tests SET raw=?, p=? WHERE id=? AND raw IS NULL",
                (
                    json_text(row.model_dump(mode="json")),
                    row.p_value if row.p_value is not None else 1.0,
                    row.hypothesis_id,
                ),
            )
            if cursor.rowcount != 1:
                raise ValueError("duplicate or unexpected raw hypothesis membership")
        if db.execute("SELECT COUNT(*) FROM tests WHERE raw IS NULL").fetchone()[0]:
            raise ValueError("raw hypothesis membership differs from declared inventory")
        probabilities = np.fromiter(
            (r[0] for r in db.execute("SELECT p FROM tests ORDER BY id")),
            dtype=np.float64,
            count=count,
        )
        by, bh = adjusted_probabilities(probabilities, policy=policy)
        del probabilities
        family_hash = _family_hash(metadata, db)
        raw_hash = _array_hash(row[0] for row in db.execute("SELECT raw FROM tests ORDER BY id"))
        declared_artifact = _write_lines(
            output / "declared.jsonl",
            (row[0] for row in db.execute("SELECT declaration FROM tests ORDER BY id")),
        )
        raw_artifact = _write_lines(
            output / "raw.jsonl",
            (row[0] for row in db.execute("SELECT raw FROM tests ORDER BY id")),
        )
        counts = {"tested": 0, "untestable": 0}

        def corrected_records() -> Iterator[str]:
            for (raw_record,), q_by, q_bh in zip(
                db.execute("SELECT raw FROM tests ORDER BY id"), by, bh, strict=True
            ):
                raw_row = RawHypothesis.model_validate_json(raw_record)
                counts[raw_row.status] += 1
                result = adjusted_hypothesis(raw_row, float(q_by), float(q_bh), policy=policy)
                yield json_text(result.model_dump(mode="json"))

        corrected_artifact = _write_lines(output / "adjusted.jsonl", corrected_records())
        with (output / "adjusted.jsonl").open(encoding="utf-8") as handle:
            output_hash = _array_hash(line.rstrip("\n") for line in handle)
    return ExternalFamilyCorrection(
        metadata=metadata,
        declared=declared_artifact,
        raw=raw_artifact,
        corrected=corrected_artifact,
        recorded_at=instant,
        policy_sha256=canonical_sha256(policy.model_dump(mode="json")),
        family_sha256=family_hash,
        input_sha256=raw_hash,
        output_sha256=output_hash,
        tested=counts["tested"],
        untestable=counts["untestable"],
    )
