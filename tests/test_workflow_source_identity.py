"""Reject contradictory source provenance claims from issue #655 (2026-10-19)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from mtp.download import (
    BUGSIGDB_EXPORTS_DOI,
    BUGSIGDB_EXPORTS_TAG_SHA,
    BUGSIGDB_EXPORTS_VERSION,
    PINNED_FULL_DUMP_URL,
)
from mtp.workflow import source_identity
from mtp.workflow.definition import read_definition
from mtp.workflow.operations import identify_sources, validate_core_product

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "workflow"
INPUTS = {
    "raw_dump": FIXTURE / "full_dump.csv",
    "kg_edges": FIXTURE / "merged-kg_edges.tsv",
    "kg_nodes": FIXTURE / "merged-kg_nodes.tsv",
}
PRODUCT = next(
    item
    for item in read_definition(ROOT / "workflows" / "bugsigdb_core.yaml").products
    if item.id == "source_identity"
)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("pinned_verified", True),
        ("verification_status", "pinned_bytes_verified"),
        ("declared_source_mode", "pinned"),
        ("snapshot_header", "# fabricated snapshot"),
        (
            "verified_release",
            {"version": "fabricated", "doi": "wrong", "commit": "wrong", "url": "wrong"},
        ),
    ],
)
def test_contradictory_source_claims_fail(tmp_path: Path, field: str, value: object) -> None:
    output = tmp_path / "source_identity.json"
    identify_sources(INPUTS, {"source_identity": output}, {"source_mode": "fixture"})
    record = json.loads(output.read_text(encoding="utf-8"))
    record[field] = value
    output.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(ValueError, match="Source identity"):
        validate_core_product(
            PRODUCT,
            output,
            INPUTS,
            trait_key_policy="raw_object_v1",
            include_policy_metadata=False,
        )


@pytest.mark.parametrize("mode", ["unknown", "fixture", "live"])
def test_valid_unverified_declarations_remain_explicit(tmp_path: Path, mode: str) -> None:
    output = tmp_path / "source_identity.json"
    identify_sources(INPUTS, {"source_identity": output}, {"source_mode": mode})
    result = validate_core_product(
        PRODUCT,
        output,
        INPUTS,
        trait_key_policy="raw_object_v1",
        include_policy_metadata=False,
    )
    assert result[0].status == "passed"
    record = json.loads(output.read_text(encoding="utf-8"))
    assert record["declared_source_mode"] == mode
    assert record["verification_status"] == "origin_unverified"
    assert not record["pinned_verified"]
    assert "verified_release" not in record


@pytest.mark.parametrize("change", [None, "downgrade", "version", "doi", "commit", "url"])
def test_verified_pin_requires_exact_release_claims(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, change: str | None
) -> None:
    output = tmp_path / "source_identity.json"
    identify_sources(INPUTS, {"source_identity": output}, {"source_mode": "fixture"})
    raw_bytes = INPUTS["raw_dump"].read_bytes()

    # Exercise release consistency without downloading the large pinned source;
    # the byte-integrity checker itself is covered by the download tests.
    def accept_fixture(data: bytes) -> None:
        assert data == raw_bytes

    monkeypatch.setattr(source_identity, "PINNED_FULL_DUMP_BYTES", len(raw_bytes))
    monkeypatch.setattr(source_identity, "verify_pinned_dump", accept_fixture)
    record = json.loads(output.read_text(encoding="utf-8"))
    record.update(
        pinned_verified=True,
        verification_status="pinned_bytes_verified",
        verified_release={
            "version": BUGSIGDB_EXPORTS_VERSION,
            "doi": BUGSIGDB_EXPORTS_DOI,
            "commit": BUGSIGDB_EXPORTS_TAG_SHA,
            "url": PINNED_FULL_DUMP_URL,
        },
    )
    if change == "downgrade":
        record.update(pinned_verified=False, verification_status="origin_unverified")
    elif change is not None:
        record["verified_release"][change] = "fabricated"
    output.write_text(json.dumps(record), encoding="utf-8")
    if change is None:
        assert (
            validate_core_product(
                PRODUCT,
                output,
                INPUTS,
                trait_key_policy="raw_object_v1",
                include_policy_metadata=False,
            )[0].status
            == "passed"
        )
    else:
        with pytest.raises(ValueError, match="Source identity"):
            validate_core_product(
                PRODUCT,
                output,
                INPUTS,
                trait_key_policy="raw_object_v1",
                include_policy_metadata=False,
            )
