"""Reconcile source provenance claims against consumed bytes (issue #655)."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from mtp.datamodel.core_products_pydantic import SourceIdentity, VerifiedRelease
from mtp.download import (
    BUGSIGDB_EXPORTS_DOI,
    BUGSIGDB_EXPORTS_TAG_SHA,
    BUGSIGDB_EXPORTS_VERSION,
    PINNED_FULL_DUMP_BYTES,
    PINNED_FULL_DUMP_URL,
    IntegrityError,
    verify_pinned_dump,
)
from mtp.paths import file_sha256


def validate_source_identity(path: Path, context: Mapping[str, Path]) -> dict[str, object]:
    """Check source identifiers, pin claims and snapshot header independently of labels."""
    source = SourceIdentity.model_validate_json(path.read_text(encoding="utf-8"))
    if set(source.consumed_files) != {"raw_dump", "kg_edges", "kg_nodes"}:
        raise ValueError("Source identity must identify exactly the consumed raw and KG files")
    for name, record in source.consumed_files.items():
        actual = context[name]
        if (
            record.product_id != name
            or record.sha256 != file_sha256(actual)
            or record.size_bytes != actual.stat().st_size
            or record.path != str(actual.resolve())
        ):
            raise ValueError(f"Source identity does not match consumed file: {name}")

    raw = context["raw_dump"]
    pinned = False
    if raw.stat().st_size == PINNED_FULL_DUMP_BYTES:
        try:
            verify_pinned_dump(raw.read_bytes())
        except IntegrityError:
            pass
        else:
            pinned = True
    expected_status = "pinned_bytes_verified" if pinned else "origin_unverified"
    if source.pinned_verified != pinned or source.verification_status != expected_status:
        raise ValueError("Source identity pin/status claims do not match consumed bytes")
    if source.declared_source_mode == "pinned" and not pinned:
        raise ValueError("Source identity declares pinned mode for unverified bytes")
    expected_release = (
        VerifiedRelease(
            version=BUGSIGDB_EXPORTS_VERSION,
            doi=BUGSIGDB_EXPORTS_DOI,
            commit=BUGSIGDB_EXPORTS_TAG_SHA,
            url=PINNED_FULL_DUMP_URL,
        )
        if pinned
        else None
    )
    if source.verified_release != expected_release:
        raise ValueError("Source identity release claim does not match consumed bytes")
    with raw.open(encoding="utf-8") as handle:
        first_line = handle.readline().rstrip("\n")
    expected_header = first_line if first_line.startswith("#") else None
    if source.snapshot_header != expected_header:
        raise ValueError("Source identity snapshot header does not match consumed raw input")
    return {"pinned_verified": pinned, "source_claims_reconciled": True}
