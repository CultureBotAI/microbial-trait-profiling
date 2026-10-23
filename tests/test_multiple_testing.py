"""Issue #852: independent arithmetic and adversarial family verification."""

from collections.abc import Sequence
from datetime import UTC, datetime
from math import inf, nan
from pathlib import Path

import pytest
import yaml

from mtp.datamodel.multiple_testing_pydantic import (
    CorrectionRun,
    HypothesisCoordinates,
    HypothesisDefinition,
    HypothesisFamily,
    RawHypothesis,
    SourceArtifact,
)
from mtp.multiple_testing import canonical_sha256, correct_family, load_policy, verify_correction

STAMP = datetime(2026, 10, 22, tzinfo=UTC)
DIGEST = "1" * 64


@pytest.mark.parametrize("value", [nan, inf, -inf])
@pytest.mark.parametrize("nested", [False, True])
def test_canonical_hash_rejects_nonfinite_numbers(value: float, nested: bool) -> None:
    payload = {"membership": [value]} if nested else value
    with pytest.raises(ValueError, match="Out of range float"):
        canonical_sha256(payload)


def family(size: int) -> HypothesisFamily:
    return HypothesisFamily(
        family_id="fixture:study-sign",
        analysis_id="fixture:screen",
        analysis_basis_sha256=DIGEST,
        source_membership_sha256=DIGEST,
        retained_feature_count=size,
        strain_rollup="off",
        feature_filter_id="sha256:" + DIGEST,
        expected_family_size=size,
        population="fixture_paired_human",
        direction_basis="curated_group1_vs_group0",
        profile_unit="count_profiles",
        inference_unit="study",
        hypothesis_unit="trait",
        estimand="study-level sign of increased-minus-decreased annotation fraction",
        selection_rule="all retained fixture traits",
        selection_history=["fixture declared before outcomes"],
        sources=[SourceArtifact(identity="fixture:raw-table", sha256=DIGEST)],
        members=[
            HypothesisDefinition(
                hypothesis_id=f"h{i}",
                feature_ids=[f"fixture:trait{i}"],
                coordinates=HypothesisCoordinates(),
                test="binomial_sign",
                tail="two_sided",
                eligibility_sha256=DIGEST,
            )
            for i in range(size)
        ],
    )


def raw(values: Sequence[float | None]) -> list[RawHypothesis]:
    return [
        RawHypothesis(
            hypothesis_id=f"h{i}",
            status="tested" if value is not None else "untestable",
            p_value=value,
            reason=None if value is not None else "no non-tied eligible studies",
        )
        for i, value in enumerate(values)
    ]


def test_independent_reference_arithmetic_and_comparison_calls() -> None:
    # BH: min_{j>=i}(m p[j]/j); BY adds H4=25/12, clipped to one.
    values = raw([0.01, 0.04, 0.04, 0.9])
    run = correct_family(family(4), values, recorded_at=STAMP)
    assert [r.q_bh for r in run.results] == pytest.approx([0.04, 4 / 75, 4 / 75, 0.9])
    assert [r.q_by for r in run.results] == pytest.approx([1 / 12, 1 / 9, 1 / 9, 1])
    assert run.results[0].reject_bh_005 is True
    assert run.results[0].reject_by_005 is False
    verify_correction(run, family=family(4), hypotheses=values)


@pytest.mark.parametrize("values", [[], [0.0], [1.0], [0.05], [0.001], [0.0, 0.0, 1.0]])
def test_empty_singleton_ties_and_endpoints(values: list[float]) -> None:
    run = correct_family(family(len(values)), raw(values), recorded_at=STAMP)
    assert len(run.results) == len(values)
    verify_correction(run, family=family(len(values)), hypotheses=raw(values))
    if len(values) == 1:
        assert run.results[0].q_by == values[0]
        assert run.results[0].q_bh == values[0]
        assert run.results[0].reject_by_005 == (values[0] <= 0.05)
        assert run.results[0].reject_by_0001 == (values[0] <= 0.001)


def test_untestable_is_reserved_but_exported_undefined() -> None:
    run = correct_family(family(2), raw([0.01, None]), recorded_at=STAMP)
    assert run.results[0].q_bh == pytest.approx(0.02)
    assert run.results[0].q_by == pytest.approx(0.03)
    missing = run.results[1]
    assert missing.p_value is None
    assert missing.q_by is None
    assert missing.q_bh is None
    assert missing.reject_by_005 is None
    assert missing.reject_bh_0001 is None
    assert missing.reason == "no non-tied eligible studies"


def test_row_order_is_not_identity_but_ratio_order_is() -> None:
    declared = family(3)
    values = raw([0.01, 0.2, None])
    run = correct_family(declared, values, recorded_at=STAMP)
    declared.members.reverse()
    assert correct_family(declared, list(reversed(values)), recorded_at=STAMP) == run
    declared.members[0].feature_ids = ["numerator", "denominator"]
    ratio = correct_family(declared, values, recorded_at=STAMP)
    declared.members[0].feature_ids.reverse()
    assert correct_family(declared, values, recorded_at=STAMP).family_sha256 != ratio.family_sha256


@pytest.mark.parametrize("value", [-0.001, 1.001, inf, -inf, nan])
def test_invalid_p_values_fail(value: float) -> None:
    with pytest.raises(ValueError, match=r"p_value|finite probability"):
        correct_family(family(1), raw([value]))


@pytest.mark.parametrize(
    ("status", "value", "reason"),
    [
        ("tested", None, None),
        ("tested", 0.1, "undefined"),
        ("untestable", 1.0, "none"),
        ("untestable", None, None),
        ("untestable", None, " "),
    ],
)
def test_test_status_and_missingness_must_agree(
    status: str, value: float | None, reason: str | None
) -> None:
    with pytest.raises(ValueError, match="hypothesis needs"):
        correct_family(
            family(1),
            [RawHypothesis(hypothesis_id="h0", status=status, p_value=value, reason=reason)],
        )


def test_missing_extra_and_duplicate_members_are_errors() -> None:
    values = raw([0.01, 0.2])
    with pytest.raises(ValueError, match="membership"):
        correct_family(family(2), values[:1])
    with pytest.raises(ValueError, match="membership"):
        correct_family(family(1), values)
    with pytest.raises(ValueError, match="duplicate"):
        correct_family(family(2), [values[0], values[0]])
    declared = family(2)
    declared.members = [declared.members[0], declared.members[0]]
    with pytest.raises(ValueError, match="duplicate"):
        correct_family(declared, values)
    declared = family(2)
    declared.expected_family_size = 1
    with pytest.raises(ValueError, match="family size"):
        correct_family(declared, values)


def test_feature_count_is_distinct_from_hypothesis_count() -> None:
    declared = family(2)
    declared.retained_feature_count = 17
    run = correct_family(declared, raw([0.01, 0.2]))
    assert run.family.retained_feature_count == 17
    assert run.family.expected_family_size == 2
    assert run.results[0].q_by == pytest.approx(0.03)


def test_policy_cannot_be_changed_under_supported_version() -> None:
    policy = load_policy()
    policy.primary_method = "bh"
    with pytest.raises(ValueError, match="policy"):
        correct_family(family(1), raw([0.2]), policy=policy)


@pytest.mark.parametrize("field", ["q_by", "q_bh", "reject_by_005", "p_value", "reason"])
def test_tampered_results_rejected_even_with_recomputed_envelope(field: str) -> None:
    run = correct_family(family(1), raw([0.01]), recorded_at=STAMP)
    payload = run.model_dump(mode="json")
    payload["results"][0][field] = {"reject_by_005": False, "reason": "wrong"}.get(field, 0.9)
    payload["output_sha256"] = canonical_sha256(payload["results"])
    envelope = {
        key: value
        for key, value in payload.items()
        if key not in {"policy", "family", "hypotheses", "results", "run_sha256"}
    }
    envelope["recorded_at"] = STAMP.isoformat()
    payload["run_sha256"] = canonical_sha256(envelope)
    tampered = CorrectionRun.model_validate(payload)
    with pytest.raises(ValueError, match="verification failed"):
        verify_correction(tampered, family=family(1), hypotheses=raw([0.01]))


@pytest.mark.parametrize("field", ["population", "estimand", "selection_rule", "direction_basis"])
def test_self_consistent_metadata_tampering_rejected_against_trusted_family(field: str) -> None:
    altered = family(1)
    setattr(altered, field, "changed metadata")
    run = correct_family(altered, raw([0.01]), recorded_at=STAMP)
    with pytest.raises(ValueError, match="family"):
        verify_correction(run, family=family(1), hypotheses=raw([0.01]))


def test_round_trip_and_raw_input_tampering() -> None:
    run = correct_family(family(2), raw([0.01, None]), recorded_at=STAMP)
    round_trip = CorrectionRun.model_validate_json(run.model_dump_json())
    assert round_trip == run
    verify_correction(round_trip, family=family(2), hypotheses=raw([0.01, None]))
    with pytest.raises(ValueError, match=r"input|hypotheses"):
        verify_correction(round_trip, family=family(2), hypotheses=raw([0.02, None]))


def test_non_utc_or_naive_timestamp_rejected() -> None:
    with pytest.raises(ValueError, match="UTC"):
        correct_family(family(1), raw([0.2]), recorded_at=datetime(2026, 10, 22))


def test_checked_in_family_fixture_and_all_untestable() -> None:
    path = Path(__file__).parent / "fixtures" / "multiple_testing_family.yaml"
    declared = HypothesisFamily.model_validate(yaml.safe_load(path.read_text(encoding="utf-8")))
    run = correct_family(declared, raw([None, None]), recorded_at=STAMP)
    assert len(run.results) == 2
    assert all(row.q_by is None for row in run.results)
    verify_correction(run, family=declared, hypotheses=raw([None, None]))


@pytest.mark.parametrize(
    "field",
    [
        "code_sha256",
        "input_sha256",
        "output_sha256",
        "run_sha256",
        "numpy_version",
        "scipy_version",
    ],
)
def test_execution_hash_tampering_rejected(field: str) -> None:
    run = correct_family(family(1), raw([0.01]), recorded_at=STAMP)
    setattr(run, field, "0" * 64)
    with pytest.raises(ValueError, match="verification failed"):
        verify_correction(run, family=family(1), hypotheses=raw([0.01]))


def test_source_metadata_and_dropped_results_are_verified() -> None:
    run = correct_family(family(2), raw([0.01, 0.2]), recorded_at=STAMP)
    run.results.pop()
    with pytest.raises(ValueError, match="results"):
        verify_correction(run, family=family(2), hypotheses=raw([0.01, 0.2]))
    declared = family(2)
    declared.sources[0].sha256 = "0" * 64
    run = correct_family(declared, raw([0.01, 0.2]), recorded_at=STAMP)
    with pytest.raises(ValueError, match="family"):
        verify_correction(run, family=family(2), hypotheses=raw([0.01, 0.2]))


def test_regex_suffixes_cannot_bypass_status_or_digest_validation() -> None:
    malformed = RawHypothesis(hypothesis_id="h0", status="untestable\n", reason="no observations")
    with pytest.raises(ValueError, match="status"):
        correct_family(family(1), [malformed])
    declared = family(1)
    declared.analysis_basis_sha256 = DIGEST + "\n"
    with pytest.raises(ValueError, match="SHA256"):
        correct_family(declared, raw([0.2]))


def test_renamed_scientific_duplicates_are_rejected() -> None:
    declared = family(2)
    declared.members[1] = declared.members[0].model_copy(
        update={"hypothesis_id": "h1", "eligibility_sha256": "2" * 64}
    )
    with pytest.raises(ValueError, match="duplicate scientific hypothesis"):
        correct_family(declared, raw([0.01, 0.01]))


def test_negative_zero_probability_is_rejected() -> None:
    with pytest.raises(ValueError, match="finite probability"):
        correct_family(family(1), raw([-0.0]))


@pytest.mark.parametrize("field", ["q_by", "q_bh", "p_value"])
def test_result_signed_zero_cannot_disagree_with_recorded_output_hash(field: str) -> None:
    run = correct_family(family(1), raw([0.0]), recorded_at=STAMP)
    setattr(run.results[0], field, -0.0)
    with pytest.raises(ValueError, match="output_sha256"):
        verify_correction(run, family=family(1), hypotheses=raw([0.0]))


def test_same_count_different_membership_is_rejected() -> None:
    values = raw([0.01, 0.2])
    values[1].hypothesis_id = "other"
    with pytest.raises(ValueError, match="membership"):
        correct_family(family(2), values)


@pytest.mark.parametrize(
    "field",
    [
        "family_id",
        "analysis_id",
        "population",
        "direction_basis",
        "profile_unit",
        "inference_unit",
        "hypothesis_unit",
        "estimand",
        "selection_rule",
    ],
)
def test_blank_family_metadata_is_rejected(field: str) -> None:
    declared = family(1)
    setattr(declared, field, " ")
    with pytest.raises(ValueError, match="blank family metadata"):
        correct_family(declared, raw([0.2]))


@pytest.mark.parametrize("history", [[], [" "]])
def test_selection_history_is_required(history: list[str]) -> None:
    declared = family(1)
    declared.selection_history = history
    with pytest.raises(ValueError, match="selection history"):
        correct_family(declared, raw([0.2]))


def test_source_identities_must_be_unique() -> None:
    declared = family(1)
    declared.sources.append(declared.sources[0].model_copy())
    with pytest.raises(ValueError, match="duplicate source identity"):
        correct_family(declared, raw([0.2]))


def test_eligibility_digest_requires_full_match() -> None:
    declared = family(1)
    declared.members[0].eligibility_sha256 = DIGEST + "\n"
    with pytest.raises(ValueError, match="SHA256"):
        correct_family(declared, raw([0.2]))


@pytest.mark.parametrize("features", [[], [" "], ["trait:x", ""]])
def test_feature_identities_must_be_nonempty(features: list[str]) -> None:
    declared = family(1)
    declared.members[0].feature_ids = features
    with pytest.raises(ValueError, match="feature identities"):
        correct_family(declared, raw([0.2]))


@pytest.mark.parametrize("coordinate", ["profile_id", "stratum_id", "split_id", "condition_ids"])
def test_explicit_coordinates_distinguish_scientific_tests(coordinate: str) -> None:
    declared = family(2)
    declared.members[1].feature_ids = declared.members[0].feature_ids.copy()
    for i, member in enumerate(declared.members):
        value = [f"source:condition{i}"] if coordinate == "condition_ids" else f"source:{i}"
        member.coordinates = HypothesisCoordinates.model_validate({coordinate: value})
    run = correct_family(declared, raw([0.01, 0.02]), recorded_at=STAMP)
    verify_correction(run, family=declared, hypotheses=raw([0.01, 0.02]))
    declared.members[1].coordinates = declared.members[0].coordinates.model_copy(deep=True)
    with pytest.raises(ValueError, match="duplicate scientific hypothesis"):
        correct_family(declared, raw([0.01, 0.02]))


def test_condition_set_order_is_canonical_but_membership_is_identity() -> None:
    declared = family(1)
    declared.members[0].coordinates.condition_ids = ["CRC", "HIV", "PD"]
    run = correct_family(declared, raw([0.01]), recorded_at=STAMP)
    declared.members[0].coordinates.condition_ids.reverse()
    assert correct_family(declared, raw([0.01]), recorded_at=STAMP) == run
    declared.members[0].coordinates.condition_ids[0] = "other"
    with pytest.raises(ValueError, match="family"):
        verify_correction(run, family=declared, hypotheses=raw([0.01]))


@pytest.mark.parametrize("coordinate", ["profile_id", "stratum_id", "split_id", "condition_ids"])
def test_blank_coordinates_are_rejected(coordinate: str) -> None:
    declared = family(1)
    value = [" "] if coordinate == "condition_ids" else " "
    declared.members[0].coordinates = HypothesisCoordinates.model_validate({coordinate: value})
    with pytest.raises(ValueError, match="blank"):
        correct_family(declared, raw([0.01]))


def test_duplicate_conditions_are_rejected() -> None:
    declared = family(1)
    declared.members[0].coordinates.condition_ids = ["PD", "PD"]
    with pytest.raises(ValueError, match="duplicate condition"):
        correct_family(declared, raw([0.01]))


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("strain_rollup", "full"),
        ("strain_rollup", "off\n"),
        ("feature_filter_id", "sha256:" + DIGEST + "\n"),
    ],
)
def test_scope_and_filter_identity_are_strict(field: str, value: str) -> None:
    declared = family(1).model_copy(update={field: value})
    with pytest.raises(ValueError, match=r"strain_rollup|strain rollup|feature filter"):
        correct_family(declared, raw([0.01]))


@pytest.mark.parametrize("field", ["strain_rollup", "feature_filter_id"])
def test_scope_and_filter_tampering_is_rejected(field: str) -> None:
    declared = family(1)
    changed = "phenotype" if field == "strain_rollup" else "sha256:" + "2" * 64
    setattr(declared, field, changed)
    run = correct_family(declared, raw([0.01]), recorded_at=STAMP)
    with pytest.raises(ValueError, match="family"):
        verify_correction(run, family=family(1), hypotheses=raw([0.01]))


def test_modified_policy_file_is_rejected(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import mtp.multiple_testing as engine

    path = tmp_path / "policy.yaml"
    content = engine.POLICY_PATH.read_text(encoding="utf-8").replace("- 0.05", "- 0.1")
    path.write_text(content, encoding="utf-8")
    monkeypatch.setattr(engine, "POLICY_PATH", path)
    with pytest.raises(ValueError, match="policy file"):
        load_policy(path)
