"""Issue #976: versioned reasons retain legacy callers and reject semantic tampering."""

import json
from pathlib import Path

import pytest

from mtp.datamodel import multiple_testing_pydantic as models
from mtp.datamodel.multiple_testing_pydantic import (
    CorrectionRun,
    HypothesisFamilyMetadata,
    RawHypothesis,
)
from mtp.multiple_testing import (
    canonical_sha256,
    correct_family,
    validate_raw_hypothesis,
    verify_correction,
)
from mtp.multiplicity.inventory import correct_inventory, read_corrected_inventory
from test_multiple_testing import STAMP, family

LEGACY = "legacy_free_text_v1"
MODERN = "analysis_reason_v1"
ANALYSIS_REASONS = (
    "all_tied",
    "no_eligible_studies",
    "no_positive_pairs",
    "insufficient_support",
    "insufficient_variance",
    "condition_below_minimum_studies",
    "shared_study_between_conditions",
    "undefined_distance_in_null",
    "all_sets_empty",
    "pool_too_small",
    "no_units",
)


def observation(
    *,
    reason: str | None = "all_tied",
    vocabulary: str | None = MODERN,
    identifier: str = "h0",
    status: str = "untestable",
    p_value: float | None = None,
) -> RawHypothesis:
    return RawHypothesis.model_validate(
        {
            "hypothesis_id": identifier,
            "status": status,
            "p_value": p_value,
            "reason": reason,
            "reason_vocabulary": vocabulary,
        }
    )


def test_generated_reason_enums_declare_exact_supported_contract() -> None:
    assert {member.value for member in models.ReasonVocabulary} == {LEGACY, MODERN}
    assert {member.value for member in models.AnalysisReasonV1} == set(ANALYSIS_REASONS)


@pytest.mark.parametrize("reason", ["no non-tied eligible studies", "  historical explanation  "])
def test_omitted_vocabulary_preserves_legacy_text_and_round_trip(reason: str) -> None:
    row = RawHypothesis(hypothesis_id="h0", status="untestable", reason=reason)
    checked = validate_raw_hypothesis(row)
    assert checked.model_dump(mode="json")["reason_vocabulary"] == LEGACY
    assert checked.reason == reason
    run = correct_family(family(1), [row], recorded_at=STAMP)
    restored = CorrectionRun.model_validate_json(run.model_dump_json())
    assert restored == run
    assert restored.results[0].reason == reason
    assert restored.results[0].model_dump(mode="json")["reason_vocabulary"] == LEGACY
    verify_correction(restored, family=family(1), hypotheses=[row])


@pytest.mark.parametrize("reason", ANALYSIS_REASONS)
def test_every_declared_modern_reason_is_accepted(reason: str) -> None:
    row = validate_raw_hypothesis(observation(reason=reason))
    assert row.reason == reason
    assert row.model_dump(mode="json")["reason_vocabulary"] == MODERN
    run = correct_family(family(1), [row], recorded_at=STAMP)
    result = run.results[0]
    assert result.reason == reason
    assert result.p_value is result.q_by is result.q_bh is None
    assert result.model_dump(mode="json")["reason_vocabulary"] == MODERN
    verify_correction(run, family=family(1), hypotheses=[row])


@pytest.mark.parametrize(
    "reason", ["no_eligible_study", "All_tied", " all_tied", "all_tied ", "all_tied\n", "", None]
)
def test_modern_reasons_reject_unknown_whitespace_and_missing_values(reason: str | None) -> None:
    with pytest.raises(ValueError, match=r"reason|Reason"):
        validate_raw_hypothesis(observation(reason=reason))


@pytest.mark.parametrize("vocabulary", ["analysis_reason_v2", " analysis_reason_v1", "", None])
def test_vocabulary_is_versioned_and_nonnull(vocabulary: str | None) -> None:
    with pytest.raises(ValueError, match="reason_vocabulary"):
        observation(vocabulary=vocabulary)


@pytest.mark.parametrize("method", ["model_copy", "model_construct"])
@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("reason", "unsupported_reason"),
        ("reason", " all_tied"),
        ("reason_vocabulary", "analysis_reason_v2"),
        ("reason_vocabulary", None),
    ],
)
def test_unvalidated_model_paths_cannot_bypass_reason_contract(
    method: str, field: str, value: str | None
) -> None:
    original = observation()
    if method == "model_copy":
        altered = original.model_copy(update={field: value})
    else:
        fields = original.model_dump(mode="json")
        fields[field] = value
        altered = RawHypothesis.model_construct(**fields)
    with pytest.raises(ValueError, match=r"reason|Reason"):
        validate_raw_hypothesis(altered)


@pytest.mark.parametrize("vocabulary", [LEGACY, MODERN])
def test_tested_rows_forbid_reasons_under_both_vocabularies(vocabulary: str) -> None:
    row = observation(vocabulary=vocabulary, status="tested", p_value=0.02, reason=None)
    assert validate_raw_hypothesis(row).p_value == 0.02
    with pytest.raises(ValueError, match="tested hypothesis"):
        validate_raw_hypothesis(row.model_copy(update={"reason": "all_tied"}))


def test_mixed_vocabularies_have_identical_memory_and_disk_results(tmp_path: Path) -> None:
    declared = family(4)
    values = [
        observation(identifier="h0", vocabulary=LEGACY, status="tested", p_value=0.01, reason=None),
        observation(identifier="h1", vocabulary=LEGACY, reason="historical undefined explanation"),
        observation(identifier="h2", reason="all_tied"),
        observation(identifier="h3", status="tested", p_value=0.1, reason=None),
    ]
    expected = correct_family(declared, values, recorded_at=STAMP)
    metadata = HypothesisFamilyMetadata.model_validate(declared.model_dump(exclude={"members"}))
    actual = correct_inventory(
        metadata, iter(declared.members), iter(values), tmp_path, recorded_at=STAMP
    )
    observed = list(read_corrected_inventory(tmp_path / "adjusted.jsonl"))
    assert observed == expected.results
    assert actual.family_sha256 == expected.family_sha256
    assert actual.input_sha256 == expected.input_sha256
    assert actual.output_sha256 == expected.output_sha256
    assert [row.model_dump(mode="json")["reason_vocabulary"] for row in observed] == [
        LEGACY,
        LEGACY,
        MODERN,
        MODERN,
    ]
    raw_records = [json.loads(line) for line in (tmp_path / "raw.jsonl").read_text().splitlines()]
    assert [row["reason_vocabulary"] for row in raw_records] == [LEGACY, LEGACY, MODERN, MODERN]


@pytest.mark.parametrize(
    "change",
    [
        {"reason": "unsupported_reason"},
        {"reason": " all_tied"},
        {"reason": None},
        {"p_value": 0.1},
        {"status": "tested", "p_value": 0.1},
    ],
)
def test_corrected_reader_validates_raw_semantics(
    tmp_path: Path, change: dict[str, str | float | None]
) -> None:
    run = correct_family(family(1), [observation()], recorded_at=STAMP)
    record = run.results[0].model_dump(mode="json")
    record.update(change)
    path = tmp_path / "adjusted.jsonl"
    path.write_text(json.dumps(record) + "\n")
    with pytest.raises(ValueError, match=r"reason|Reason|hypothesis"):
        list(read_corrected_inventory(path))


@pytest.mark.parametrize(("original", "replacement"), [(MODERN, LEGACY), (LEGACY, MODERN)])
def test_rehashed_vocabulary_tampering_fails_against_trusted_inputs(
    original: str, replacement: str
) -> None:
    values = [observation(vocabulary=original)]
    run = correct_family(family(1), values, recorded_at=STAMP)
    payload = run.model_dump(mode="json")
    payload["hypotheses"][0]["reason_vocabulary"] = replacement
    payload["results"][0]["reason_vocabulary"] = replacement
    payload["input_sha256"] = canonical_sha256(payload["hypotheses"])
    payload["output_sha256"] = canonical_sha256(payload["results"])
    assert payload["input_sha256"] != run.input_sha256
    assert payload["output_sha256"] != run.output_sha256
    envelope = {
        key: value
        for key, value in payload.items()
        if key not in {"policy", "family", "hypotheses", "results", "run_sha256"}
    }
    envelope["recorded_at"] = STAMP.isoformat()
    payload["run_sha256"] = canonical_sha256(envelope)
    altered = CorrectionRun.model_validate(payload)
    with pytest.raises(ValueError, match="verification failed"):
        verify_correction(altered, family=family(1), hypotheses=values)
