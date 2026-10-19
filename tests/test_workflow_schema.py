"""Scientific distinctions enforced by the generated issue-629 LinkML models."""

from __future__ import annotations

from typing import Any

import pytest
from pydantic import ValidationError

from mtp.datamodel.core_products_pydantic import CountProfileRow, ORAProfileRow, SourceIdentity
from mtp.workflow.models import LegacyImportRecord


def profile() -> dict[str, Any]:
    return {
        "study_id": "bsdb:1",
        "experiment_id": "bsdb:1/1",
        "inc_dec": "inc",
        "trait_filter_id": "sha256:" + "a" * 64,
        "feature_counts": [{"feature_id": "METPO:0001", "count": 7}],
    }


def test_profile_unit_and_denominator_contracts_are_distinct() -> None:
    assert CountProfileRow.model_validate(profile()).unit == "rolled_up_source_taxon_count"
    with pytest.raises(ValidationError, match="sample_size"):
        ORAProfileRow.model_validate(profile())
    value = {**profile(), "sample_size": 9}
    assert ORAProfileRow.model_validate(value).unit == "distinct_matched_reported_taxon_count"
    with pytest.raises(ValidationError, match="sample_size"):
        CountProfileRow.model_validate(value)


@pytest.mark.parametrize("invalid", [-1, 1.5, "unknown"])
def test_profile_counts_cannot_be_negative_fractional_or_missing(invalid: object) -> None:
    value = profile()
    value["feature_counts"][0]["count"] = invalid
    with pytest.raises(ValidationError, match="count"):
        CountProfileRow.model_validate(value)


def test_source_identity_requires_explicit_verified_byte_bindings() -> None:
    value = {
        "recorded_at_utc": "2026-10-19T00:00:00Z",
        "declared_source_mode": "fixture",
        "pinned_verified": False,
        "verification_status": "origin_unverified",
        "interpretation": "Annotation evidence",
        "consumed_files": {
            "raw_dump": {
                "product_id": "raw_dump",
                "path": "fixture.csv",
                "sha256": "a" * 64,
                "size_bytes": 1,
            },
        },
    }
    assert not SourceIdentity.model_validate(value).pinned_verified
    invalid = {**value, "verification_status": "trusted_without_checks"}
    with pytest.raises(ValidationError, match="verification_status"):
        SourceIdentity.model_validate(invalid)


def test_legacy_import_cannot_assert_a_successful_execution() -> None:
    value = {
        "id": "legacy:1",
        "source_path": "old.yaml",
        "source_sha256": "a" * 64,
        "imported_at": "2026-10-19T00:00:00Z",
        "verification_status": "completed",
        "claims_json": "{}",
    }
    with pytest.raises(ValidationError, match="verification_status"):
        LegacyImportRecord.model_validate(value)
