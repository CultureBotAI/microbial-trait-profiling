"""Whole-study split selection and held-out inference regressions (#854)."""

from __future__ import annotations

from dataclasses import replace

import pandas as pd
import pytest

from mtp.cellulose_decr.study_split import (
    Candidate,
    Protocol,
    evaluate_candidates,
    make_splits,
    partition_frames,
    select_candidates,
)
from mtp.trait_guilds import Profile
from mtp.workflow.coverage import feature_filter_id


def corpus() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[Profile]]:
    profile_rows = []
    taxa_rows = []
    counts = []
    for index in range(18):
        study = f"s{index:02}"
        condition = ("A", "B", "C")[index % 3]
        for experiment in ("e1", "e2"):
            for direction in ("inc", "dec"):
                key = f"{study}|{experiment}|{direction}"
                row = {
                    "study_id": study,
                    "experiment_id": experiment,
                    "condition": condition,
                    "body_site": "Feces",
                    "compartment": "Feces",
                    "dir": direction,
                }
                profile_rows.append(dict(profile_id=key, **row))
                taxa_rows.append(dict(profile_id=key, taxon_id="t1", **row))
                counts.append(
                    Profile(
                        study,
                        experiment,
                        direction,
                        condition,
                        "Feces",
                        {"trait": 1 if direction == "inc" else 3, "other": 3},
                        4 if direction == "inc" else 6,
                    )
                )
    return (
        pd.DataFrame(profile_rows).set_index("profile_id"),
        pd.DataFrame(taxa_rows),
        pd.DataFrame({"trait": [1, 0], "other": [0, 1]}, index=["t1", "t2"]),
        counts,
    )


def test_splits_are_complete_disjoint_reproducible_and_order_independent() -> None:
    studies = [f"s{i}" for i in range(11)]
    forward = make_splits(studies, splits=5, seed=7)
    assert forward == make_splits(list(reversed(studies)), splits=5, seed=7)
    for split in forward:
        assert len(split.training) == 5
        assert not set(split.training) & set(split.testing)
        assert set(split.training) | set(split.testing) == set(studies)
    assert forward != make_splits(studies, splits=5, seed=8)


def test_partition_keeps_every_experiment_and_taxon_of_each_study() -> None:
    profiles, taxa, _, _ = corpus()
    split = make_splits(profiles.study_id.tolist(), splits=1, seed=9)[0]
    training_profiles, training_taxa = partition_frames(profiles, taxa, split.training)
    assert set(training_profiles.study_id) == set(split.training)
    assert set(training_taxa.study_id) == set(split.training)
    assert training_profiles.groupby("study_id").size().eq(4).all()
    assert training_taxa.groupby("study_id").size().eq(4).all()


def test_training_selection_does_not_force_historical_lead_or_use_heldout() -> None:
    profiles, taxa, matrix, counts = corpus()
    split = make_splits(profiles.study_id.tolist(), splits=1, seed=3)[0]
    selected = select_candidates(
        profiles,
        taxa,
        matrix,
        counts,
        ("trait", "other"),
        split,
        Protocol(min_studies=1, permutations=19),
    )
    assert len(selected) == 1
    assert selected[0].trait == "trait"
    assert selected[0].conditions == ("A", "B", "C")
    heldout = set(split.testing)
    changed = [p for p in counts if p.study_id not in heldout]
    profiles = profiles.copy()
    profiles.loc[profiles.study_id.isin(heldout), "condition"] = "HIDDEN"
    taxa = taxa[~taxa.study_id.isin(heldout)]
    assert selected == select_candidates(
        profiles,
        taxa,
        matrix,
        changed,
        ("trait", "other"),
        split,
        Protocol(min_studies=1, permutations=19),
    )


def test_missing_heldout_condition_reserves_both_tests_without_replacement() -> None:
    profiles, taxa, matrix, _ = corpus()
    candidate = Candidate("trait", ("A", "B", "missing"), 0.0)
    tests, details = evaluate_candidates(
        [candidate],
        profiles,
        taxa,
        matrix,
        Protocol(min_studies=1, permutations=19),
        seed=4,
    )
    assert len(tests) == 2
    assert all(row.status == "untestable" and row.p_value is None for row in tests)
    assert {row.reason for row in tests} == {"condition_below_minimum_studies"}
    assert len(details) == 2


def test_actual_scanner_runs_on_heldout_pool_and_keeps_both_statistics() -> None:
    profiles, taxa, matrix, _ = corpus()
    tests, details = evaluate_candidates(
        [Candidate("trait", ("A", "B", "C"), 0.0)],
        profiles,
        taxa,
        matrix,
        Protocol(min_studies=1, permutations=19),
        seed=4,
    )
    assert len(tests) == 2
    assert all(row.status == "tested" and row.p_value == 1.0 for row in tests)
    assert {row["pool_studies"] for row in details} == {18}


@pytest.mark.parametrize("studies", [[], ["only"], ["", "s2"]])
def test_invalid_study_inventory_fails(studies: list[str]) -> None:
    with pytest.raises(ValueError, match="two nonblank study"):
        make_splits(studies, splits=1, seed=2)


def test_shared_heldout_study_is_untestable_instead_of_using_disjoint_null() -> None:
    profiles, taxa, matrix, _ = corpus()
    mixed = profiles.study_id.eq("s00") & profiles.experiment_id.eq("e2")
    profiles.loc[mixed, "condition"] = "B"
    taxa.loc[taxa.study_id.eq("s00") & taxa.experiment_id.eq("e2"), "condition"] = "B"
    tests, _ = evaluate_candidates(
        [Candidate("trait", ("A", "B", "C"), 0.0)],
        profiles,
        taxa,
        matrix,
        Protocol(min_studies=1, permutations=19),
        seed=4,
    )
    assert len(tests) == 2
    assert all(row.reason == "shared_study_between_conditions" for row in tests)


def test_undefined_null_distance_is_not_counted_as_non_extreme() -> None:
    profiles, taxa, matrix, _ = corpus()
    extra = [f"s{i:02}" for i in range(9, 18)]
    profiles.loc[profiles.study_id.isin(extra), "condition"] = "D"
    taxa.loc[taxa.study_id.isin(extra), "condition"] = "D"
    taxa.loc[~taxa.study_id.eq("s00"), "taxon_id"] = "t2"
    tests, _ = evaluate_candidates(
        [Candidate("trait", ("A", "B", "C"), 1.0)],
        profiles,
        taxa,
        matrix,
        Protocol(min_studies=1, permutations=99),
        seed=4,
    )
    assert tests[0].status == "untestable"
    assert tests[0].reason == "undefined_distance_in_null"
    assert tests[1].status == "tested"
    assert tests[1].p_value == 1.0


def test_training_study_leakage_is_rejected_before_selection() -> None:
    from mtp.cellulose_decr.study_split import StudySplit

    with pytest.raises(ValueError, match="leakage"):
        StudySplit(0, ("s1",), ("s1",))


def test_family_keeps_missing_tests_and_rejects_tampering() -> None:
    from mtp.cellulose_decr.study_split import StudySplit, declare_family
    from mtp.datamodel.multiple_testing_pydantic import RawHypothesis, SourceArtifact
    from mtp.multiple_testing import correct_family, verify_correction
    from mtp.workflow.downstream_contracts import build_analysis_basis

    basis = build_analysis_basis(
        basis_id="fixture",
        feature_catalog_sha256="a" * 64,
        feature_filter_id=feature_filter_id(("trait", "other")),
        retained_feature_count=2,
        profile_unit="count_profiles",
        aggregation_policy="sum",
        rollup_policy="ncbi_taxonomy_closure_rollup_off",
        denominator_policy="retained",
        background_inventory_id="none",
        weighting_policy="equal_study",
    )
    candidates = [Candidate("trait", ("A", "B", "C"), 0.0)]
    family = declare_family(
        candidates,
        StudySplit(0, ("s1",), ("s2",)),
        basis,
        [SourceArtifact(identity="fixture", sha256="c" * 64)],
        analysis_id="fixture",
        strain_rollup="off",
    )
    assert family.expected_family_size == 2
    assert family.retained_feature_count == 2
    assert family.strain_rollup == "off"
    assert family.feature_filter_id == basis.feature_filter_id
    assert family.profile_unit == basis.profile_unit
    assert all(member.coordinates.condition_ids == ["A", "B", "C"] for member in family.members)
    assert all(member.coordinates.split_id == "split-000" for member in family.members)
    raw = [
        RawHypothesis(
            hypothesis_id=member.hypothesis_id,
            status="untestable",
            reason="condition_below_minimum_studies",
        )
        for member in family.members
    ]
    with pytest.raises(ValueError, match="membership"):
        correct_family(family, raw[:1])
    run = correct_family(family, raw)
    assert all(row.p_value is None and row.q_by is None for row in run.results)
    run.results[0].q_by = 0.01
    with pytest.raises(ValueError, match="verification failed"):
        verify_correction(run, family=family, hypotheses=raw)


@pytest.mark.parametrize("selection", ["selected", "high_minimum", "sparse_conditions"])
def test_repeated_split_report_has_no_aggregate_p_value(selection: str) -> None:
    from datetime import UTC, datetime

    from mtp.cellulose_decr.study_split import run_splits, summarize_stability
    from mtp.datamodel.multiple_testing_pydantic import SourceArtifact
    from mtp.workflow.downstream_contracts import build_analysis_basis

    profiles, taxa, matrix, counts = corpus()
    if selection == "sparse_conditions":
        profiles["condition"] = profiles["study_id"]
        taxa["condition"] = taxa["study_id"]
        counts = [replace(profile, condition=profile.study_id) for profile in counts]
    minimum = {"selected": 1, "high_minimum": 100, "sparse_conditions": 3}[selection]
    basis = build_analysis_basis(
        basis_id="fixture",
        feature_catalog_sha256="a" * 64,
        feature_filter_id=feature_filter_id(("trait", "other")),
        retained_feature_count=2,
        profile_unit="count_profiles",
        aggregation_policy="sum",
        rollup_policy="ncbi_taxonomy_closure_rollup_off",
        denominator_policy="retained",
        background_inventory_id="none",
        weighting_policy="equal_study",
    )
    result = run_splits(
        profiles,
        taxa,
        matrix,
        counts,
        ("trait", "other"),
        basis,
        [SourceArtifact(identity="fixture", sha256="c" * 64)],
        Protocol(splits=4, min_studies=minimum, permutations=19),
        analysis_id="fixture",
        strain_rollup="off",
        recorded_at=datetime(2026, 10, 22, tzinfo=UTC),
    )
    assert len(result) == 4
    for record in result:
        training = set(record["split"]["training"])
        testing = set(record["split"]["testing"])
        assert not training & testing
        assert len(training | testing) == 18
        assert record["correction"]["family"]["expected_family_size"] == 2 * len(
            record["selected_candidates"]
        )
    summary = summarize_stability(result)
    if selection != "selected":
        assert summary == []
        assert all(record["selected_candidates"] == [] for record in result)
        assert all(record["correction"]["results"] == [] for record in result)
        assert all(
            record["selection_status"] == "fewer_than_three_eligible_conditions"
            for record in result
        )
        return
    assert summary
    assert all(record["selection_status"] == "selected" for record in result)
    assert all(row["total_splits"] == 4 for row in summary)
    assert all("p_value" not in row and "q_by" not in row for row in summary)
