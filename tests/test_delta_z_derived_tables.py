"""Independent synthetic contracts for derived Delta-z tables (#1052, 2026-10-25)."""

from __future__ import annotations

import math
from typing import Literal

import numpy as np
import pandas as pd
import pytest
from numpy.typing import NDArray
from scipy import __version__ as SCIPY_VERSION
from scipy import stats

from mtp.delta_z_derived import tables as tables_module
from mtp.delta_z_derived.tables import build_tables, trait_summary

FA = "METPO:1000608"
OA = "METPO:1000607"
ANAEROBIC = "METPO:1000603"
SIGN = "TEST:threshold"
SHIFT = "TEST:shift"
ZERO = "TEST:zero"
TRAITS = [FA, OA, ANAEROBIC, SIGN, SHIFT, ZERO]
Inputs = tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]
REVIEWED_SCIPY_CONTRACTS = {
    "1.17.1": "scipy-1.17.1-wilcoxon-auto-wilcox-v1",
    "1.18.1": "scipy-1.18.1-wilcoxon-auto-wilcox-v1",
}


@pytest.fixture
def inputs() -> Inputs:
    """Asymmetric pairs, repeated experiment IDs across studies and exact site labels."""
    index = pd.MultiIndex.from_tuples(
        [("s1", "e1"), ("s1", "e2"), ("s1", "e3")] + [(f"s{i}", "e1") for i in range(2, 7)],
        names=["study_id", "experiment_id"],
    )
    pairs = pd.DataFrame(
        {
            FA: [6, 2, 2, -3, 0, 1e-13, -1e-13, 0],
            OA: [-3, 1, -4, 2, -2, -1e-13, -1e-13, 0],
            ANAEROBIC: [0.0] * 8,
            SIGN: [2e-12, 2e-12, 2e-12, -2e-12, 1e-12, -1e-12, 0, 0],
            SHIFT: [1, 2, 4, -8, 16, -32, 64, 0],
            ZERO: [0.0] * 8,
        },
        index=index,
        dtype=float,
    )
    studies = pairs.groupby(level="study_id").mean()
    evidence = pd.DataFrame(
        {
            "trait": TRAITS,
            "trait_label": [
                "facultatively aerobic",
                "obligately anaerobic",
                "anaerobic",
                "threshold",
                "shift",
                "zero",
            ],
            "experiment_mean_delta_z": pairs.mean().reindex(TRAITS).to_numpy(),
            "study_mean_delta_z": studies.mean().reindex(TRAITS).to_numpy(),
            "noel_studies_up": [3] * len(TRAITS),
            "noel_studies_down": [2] * len(TRAITS),
            "noel_studies_tied": [1] * len(TRAITS),
        }
    )
    for name in ("experiment_mean_delta_z", "study_mean_delta_z"):
        evidence[f"{name}_rank_descending"] = evidence[name].rank(method="min", ascending=False)
    metadata = index.to_frame(index=False).assign(
        body_site=["Feces", "Feces", "Feces", "Feces", "Oral", "feces", "Feces", ""],
        condition=[
            "Disease B",
            "Disease B",
            "Disease B",
            "Disease A",
            "Disease B",
            "Disease A",
            "Disease C",
            "",
        ],
        body_site_category=["stool"] * 8,
    )
    return pairs, studies, evidence, metadata


def test_study_signs_use_supplied_study_ledger_and_strict_threshold(inputs: Inputs) -> None:
    pairs, studies, evidence, metadata = inputs
    result = build_tables(pairs, studies, evidence, metadata, [SHIFT, ZERO])
    summary = result.trait_summary.set_index("trait")
    assert list(summary.index) == TRAITS
    assert summary.loc[FA, "sdz_up"] == 1
    assert summary.loc[FA, "sdz_down"] == 1
    assert summary.loc[FA, "sdz_tied"] == 4
    assert summary.loc[FA, "sdz_active_studies"] == 2
    assert summary.loc[SIGN, "sdz_up"] == 1
    assert summary.loc[SIGN, "sdz_down"] == 1
    assert summary.loc[SIGN, "sdz_tied"] == 4
    assert summary.loc[ZERO, "sdz_active_studies"] == 0
    assert set(summary.sdz_studies_total) == {6}
    assert set(summary.sdz_active_threshold) == {1e-12}
    assert summary.loc[FA, "noel_studies_up"] == 3
    assert summary.loc[FA, "noel_studies_down"] == 2
    assert summary.loc[SHIFT, "experiment_mean_delta_z"] == pytest.approx(47 / 8)
    assert summary.loc[SHIFT, "study_mean_delta_z"] == pytest.approx((7 / 3 - 8 + 16 - 32 + 64) / 6)
    pd.testing.assert_frame_equal(result.trait_summary[evidence.columns], evidence)


def test_summary_only_has_no_oxygen_or_headline_list_requirement(inputs: Inputs) -> None:
    _, studies, evidence, _ = inputs
    selected = [SHIFT, ZERO]
    summary = trait_summary(studies[selected], evidence.loc[evidence.trait.isin(selected)])
    assert list(summary.trait) == selected
    assert list(summary.sdz_active_studies) == [5, 0]


def test_active_study_support_boundary_is_exactly_fifty(inputs: Inputs) -> None:
    _, _, evidence, _ = inputs
    studies = pd.DataFrame(
        {
            trait: [1.0] * n + [0.0] * (51 - n)
            for trait, n in zip(TRAITS, [49, 50, 51, 0, 0, 0], strict=True)
        },
        index=pd.Index([f"study{i}" for i in range(51)], name="study_id"),
    )
    evidence = evidence.copy()
    evidence["study_mean_delta_z"] = studies.mean().reindex(TRAITS).to_numpy()
    summary = trait_summary(studies, evidence)
    assert list(summary.sdz_active_studies) == [49, 50, 51, 0, 0, 0]
    assert summary.sdz_active_at_least_50.tolist() == [False, True, True, False, False, False]


def test_quadrants_use_obligate_label_exact_signs_and_all_values_for_pearson(
    inputs: Inputs,
) -> None:
    pairs, studies, evidence, metadata = inputs
    result = build_tables(pairs, studies, evidence, metadata, [SHIFT, ZERO])
    summary = result.quadrant_summary.set_index("population")
    all_rows = summary.loc["all"]
    assert all_rows["fa_up_oa_down"] == 3
    assert all_rows["fa_down_oa_up"] == 1
    assert all_rows["both_up"] == 1
    assert all_rows["both_down"] == 1
    assert all_rows["axis"] == 2
    assert all_rows["n_comparisons"] == 8
    assert all_rows["n_studies"] == 6
    assert all_rows["sign_test_n"] == 4
    assert all_rows["sign_test_p"] == pytest.approx(0.625)
    assert all_rows["pearson_r"] == pytest.approx(stats.pearsonr(pairs[FA], pairs[OA]).statistic)
    assert set(summary.fa_trait) == {FA}
    assert set(summary.oa_trait) == {OA}
    assert set(summary.sign_threshold) == {0.0}
    assert set(summary.p_adjustment) == {"none_nominal"}
    stool = summary.loc["stool"]
    assert stool["n_comparisons"] == 5
    assert stool["n_studies"] == 3
    assert stool["fa_up_oa_down"] == 2
    assert stool["fa_down_oa_up"] == 1
    assert stool["sign_test_p"] == 1.0
    per_pair = result.quadrant_pairs.set_index(["study_id", "experiment_id"])
    assert per_pair.loc[("s4", "e1"), "quadrant"] == "fa_up_oa_down"
    assert not per_pair.loc[("s4", "e1"), "is_stool"]
    assert per_pair.loc[("s5", "e1"), "quadrant"] == "both_down"
    assert per_pair.loc[("s6", "e1"), "quadrant"] == "axis"


def test_stool_uses_exact_site_top_diseases_and_nominal_comparison_wilcoxon(inputs: Inputs) -> None:
    pairs, studies, evidence, metadata = inputs
    result = build_tables(pairs, studies, evidence, metadata, [SHIFT, ZERO], top_diseases=2)
    table = result.stool_disease_shifts
    assert list(zip(table.condition, table.trait, strict=True)) == [
        ("Disease B", SHIFT),
        ("Disease B", ZERO),
        ("Disease A", SHIFT),
        ("Disease A", ZERO),
    ]
    b = table.loc[table.condition.eq("Disease B") & table.trait.eq(SHIFT)].iloc[0]
    assert b["n_comparisons"] == 3
    assert b["n_studies"] == 1
    assert b["mean_delta_z"] == pytest.approx(7 / 3)
    assert b["wilcoxon_statistic"] == 0.0
    assert b["wilcoxon_p"] == pytest.approx(0.25)
    assert b["wilcoxon_status"] == "tested"
    assert set(table.test_unit) == {"comparison"}
    assert set(table.p_adjustment) == {"none_nominal"}
    zero = table.loc[table.trait == ZERO]
    assert set(zero.wilcoxon_status) == {"all_tied"}
    assert zero.wilcoxon_p.isna().all()
    assert zero.wilcoxon_statistic.isna().all()


def test_row_two_quadrants_exclude_multisite_nonaxis_pairs(inputs: Inputs) -> None:
    pairs, studies, evidence, metadata = inputs
    original = build_tables(pairs, studies, evidence, metadata, [SHIFT, ZERO])
    metadata = metadata.copy()
    multisite = metadata.study_id.eq("s4") & metadata.experiment_id.eq("e1")
    metadata.loc[multisite, "body_site"] = "Feces,Blood"
    result = build_tables(pairs, studies, evidence, metadata, [SHIFT, ZERO])
    summary = result.quadrant_summary.set_index("population")
    counts = ["fa_up_oa_down", "fa_down_oa_up", "both_up", "both_down", "axis"]
    all_rows = summary.loc[["all"]].iloc[0]
    stool = summary.loc[["stool"]].iloc[0]
    assert all_rows["n_comparisons"] == 8
    assert all_rows["n_studies"] == 6
    assert all_rows[counts].tolist() == [3, 1, 1, 1, 2]
    assert stool["n_comparisons"] == 5
    assert stool["n_studies"] == 3
    assert stool[counts].tolist() == [2, 1, 1, 1, 0]
    per_pair = result.quadrant_pairs.set_index(["study_id", "experiment_id"])
    assert set(per_pair.loc[per_pair.is_stool].index) == {
        ("s1", "e1"),
        ("s1", "e2"),
        ("s1", "e3"),
        ("s2", "e1"),
        ("s5", "e1"),
    }
    excluded = per_pair.loc[[("s4", "e1")]].iloc[0]
    assert excluded["body_site"] == "Feces,Blood"
    assert excluded["fa_delta_z"] > 0
    assert excluded["oa_delta_z"] < 0
    assert excluded["quadrant"] == "fa_up_oa_down"
    assert not excluded["is_stool"]
    pd.testing.assert_series_equal(
        all_rows, original.quadrant_summary.set_index("population").loc[["all"]].iloc[0]
    )
    pd.testing.assert_frame_equal(result.stool_disease_shifts, original.stool_disease_shifts)
    pd.testing.assert_frame_equal(result.stool_study_sensitivity, original.stool_study_sensitivity)


def test_row_five_preserves_every_exact_cell_and_reports_two_counts(inputs: Inputs) -> None:
    pairs, studies, evidence, metadata = inputs
    result = build_tables(pairs, studies, evidence, metadata, [SHIFT])
    table = result.site_disease_shifts
    assert len(table) == 6 * len(TRAITS)
    assert set(table.trait) == set(TRAITS)
    cells = table.set_index(["body_site", "condition", "trait"])
    feces = cells.loc[[("Feces", "Disease B", SHIFT)]].iloc[0]
    assert feces["n_comparisons"] == 3
    assert feces["n_studies"] == 1
    assert feces["mean_delta_z"] == pytest.approx(7 / 3)
    assert cells.loc[("Oral", "Disease B", SHIFT), "mean_delta_z"] == 16
    assert cells.loc[("feces", "Disease A", SHIFT), "mean_delta_z"] == -32
    assert cells.loc[("", "", SHIFT), "n_comparisons"] == 1
    assert cells.loc[("", "", SHIFT), "mean_delta_z"] == 0.0


def test_group_means_and_wilcoxon_count_comparisons_not_study_means(inputs: Inputs) -> None:
    pairs, _, evidence, metadata = inputs
    pairs = pairs.copy()
    pairs.loc[pairs.index[:4], SHIFT] = [1.0, 2.0, 3.0, -10.0]
    studies = pairs.groupby(level="study_id").mean()
    evidence = evidence.copy()
    evidence["experiment_mean_delta_z"] = pairs.mean().reindex(TRAITS).to_numpy()
    evidence["study_mean_delta_z"] = studies.mean().reindex(TRAITS).to_numpy()
    metadata = metadata.copy()
    metadata.loc[:3, "condition"] = "Combined condition"
    result = build_tables(pairs, studies, evidence, metadata, [SHIFT], top_diseases=1)
    stool = result.stool_disease_shifts.iloc[0]
    assert stool["mean_delta_z"] == -1.0
    assert stool["n_comparisons"] == 4
    assert stool["n_studies"] == 2
    assert stool["wilcoxon_statistic"] == 4.0
    assert stool["wilcoxon_p"] == pytest.approx(0.875)
    cell = result.site_disease_shifts.loc[
        result.site_disease_shifts.body_site.eq("Feces")
        & result.site_disease_shifts.condition.eq("Combined condition")
        & result.site_disease_shifts.trait.eq(SHIFT)
    ].iloc[0]
    assert cell["mean_delta_z"] == -1.0
    assert cell["n_comparisons"] == 4
    assert cell["n_studies"] == 2


def test_absent_stool_population_is_explicit_and_keeps_output_schema(inputs: Inputs) -> None:
    pairs, studies, evidence, metadata = inputs
    metadata = metadata.assign(body_site="Saliva")
    result = build_tables(pairs, studies, evidence, metadata, [SHIFT])
    stool = result.quadrant_summary.loc[result.quadrant_summary.population.eq("stool")].iloc[0]
    assert stool["n_comparisons"] == 0
    assert stool["n_studies"] == 0
    assert stool["sign_test_status"] == "no_discordant_pairs"
    assert stool["pearson_status"] == "insufficient_support"
    assert result.stool_disease_shifts.empty
    assert {"trait", "n_comparisons", "wilcoxon_p"}.issubset(result.stool_disease_shifts.columns)


def test_inputs_are_not_changed_and_metadata_order_is_irrelevant(inputs: Inputs) -> None:
    pairs, studies, evidence, metadata = inputs
    originals = [frame.copy(deep=True) for frame in inputs]
    expected = build_tables(*inputs, [SHIFT, ZERO])
    extra = metadata.iloc[[0]].assign(study_id="unpaired", experiment_id="unused")
    observed = build_tables(
        pairs.iloc[::-1],
        studies.iloc[::-1],
        evidence,
        pd.concat([extra, metadata.iloc[::-1]]),
        [SHIFT, ZERO],
    )
    for field in (
        "trait_summary",
        "quadrant_summary",
        "stool_disease_shifts",
        "site_disease_shifts",
    ):
        pd.testing.assert_frame_equal(getattr(expected, field), getattr(observed, field))
    for actual, original in zip(inputs, originals, strict=True):
        pd.testing.assert_frame_equal(actual, original)


def test_small_round_trip_noise_is_allowed_but_changed_study_values_are_not(inputs: Inputs) -> None:
    pairs, studies, evidence, metadata = inputs
    rounded = studies.copy()
    rounded.loc["s1", SHIFT] = 7 / 3 + 1e-15
    build_tables(pairs, rounded, evidence, metadata, [SHIFT])
    altered = studies.copy()
    altered.loc["s1", SHIFT] = 7 / 3 + 0.01
    with pytest.raises(ValueError, match=r"study|Study"):
        build_tables(pairs, altered, evidence, metadata, [SHIFT])


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf, "not numeric"])
def test_nonfinite_or_non_numeric_delta_z_is_rejected(inputs: Inputs, bad: float | str) -> None:
    pairs, studies, evidence, metadata = inputs
    damaged = pairs.astype(object)
    damaged.iloc[0, 0] = bad
    with pytest.raises(ValueError, match=r"finite|numeric"):
        build_tables(damaged, studies, evidence, metadata, [SHIFT])


@pytest.mark.parametrize("which", ["pairs", "studies", "evidence", "metadata"])
def test_duplicate_identity_is_rejected(inputs: Inputs, which: str) -> None:
    frames = list(inputs)
    position = {"pairs": 0, "studies": 1, "evidence": 2, "metadata": 3}[which]
    frames[position] = pd.concat([frames[position], frames[position].iloc[[0]]])
    with pytest.raises(ValueError, match=r"[Dd]uplicate|unique"):
        build_tables(frames[0], frames[1], frames[2], frames[3], [SHIFT])


@pytest.mark.parametrize("damage", ["missing_pair", "blank_key", "duplicate_unpaired"])
def test_metadata_identity_guard_rejects_loss_or_ambiguity(inputs: Inputs, damage: str) -> None:
    pairs, studies, evidence, metadata = inputs
    metadata = metadata.copy()
    if damage == "missing_pair":
        metadata = metadata.iloc[1:]
    elif damage == "blank_key":
        metadata.loc[0, "study_id"] = " "
    else:
        extra = metadata.iloc[[0]].assign(study_id="unpaired")
        metadata = pd.concat([metadata, extra, extra])
    with pytest.raises(ValueError, match=r"metadata|identity|key|[Dd]uplicate|unique"):
        build_tables(pairs, studies, evidence, metadata, [SHIFT])


@pytest.mark.parametrize("traits", [[], [SHIFT, SHIFT], ["TEST:absent"]])
def test_explicit_trait_list_must_be_unique_nonempty_and_present(
    inputs: Inputs, traits: list[str]
) -> None:
    with pytest.raises(ValueError, match=r"trait|Trait"):
        build_tables(*inputs, traits)


def test_oxygen_label_cannot_be_substituted_and_catalogs_must_match(inputs: Inputs) -> None:
    pairs, studies, evidence, metadata = inputs
    with pytest.raises(ValueError, match=r"oxygen|obligate|METPO:1000607|catalog"):
        build_tables(
            pairs.drop(columns=OA),
            studies.drop(columns=OA),
            evidence.loc[evidence.trait != OA],
            metadata,
            [SHIFT],
        )
    with pytest.raises(ValueError, match=r"catalog|trait|Trait"):
        build_tables(pairs, studies.drop(columns=ZERO), evidence, metadata, [SHIFT])
    bad_evidence = evidence.copy()
    bad_evidence.loc[bad_evidence.trait == SHIFT, "experiment_mean_delta_z"] = 999
    with pytest.raises(ValueError, match=r"mean|evidence"):
        build_tables(pairs, studies, bad_evidence, metadata, [SHIFT])


def test_zero_only_oxygen_vectors_are_explicitly_untestable(inputs: Inputs) -> None:
    pairs, studies, evidence, metadata = inputs
    pairs[[FA, OA]] = 0.0
    studies[[FA, OA]] = 0.0
    evidence.loc[
        evidence.trait.isin([FA, OA]), ["experiment_mean_delta_z", "study_mean_delta_z"]
    ] = 0.0
    result = build_tables(pairs, studies, evidence, metadata, [SHIFT])
    assert set(result.quadrant_summary.sign_test_status) == {"no_discordant_pairs"}
    assert result.quadrant_summary.sign_test_p.isna().all()
    assert set(result.quadrant_summary.pearson_status) == {"constant_input"}
    assert result.quadrant_summary.pearson_r.isna().all()


def _cell_inputs(
    values: list[float],
    *,
    study_ids: list[str] | None = None,
    sites: list[str] | None = None,
) -> Inputs:
    """One condition with arbitrary comparison/study multiplicity and exact sites."""
    count = len(values)
    ids = study_ids if study_ids is not None else [f"s{index}" for index in range(count)]
    index = pd.MultiIndex.from_arrays(
        [ids, [f"e{number}" for number in range(count)]], names=["study_id", "experiment_id"]
    )
    pairs = pd.DataFrame({FA: np.ones(count), OA: -np.ones(count), SHIFT: values}, index=index)
    studies = pairs.groupby(level="study_id").mean()
    evidence = pd.DataFrame(
        {
            "trait": [FA, OA, SHIFT],
            "trait_label": ["verified facultative", "verified obligate", "shift"],
            "experiment_mean_delta_z": pairs.mean().to_numpy(),
            "study_mean_delta_z": studies.mean().to_numpy(),
            "noel_studies_up": [0, 0, 0],
            "noel_studies_down": [0, 0, 0],
            "noel_studies_tied": [len(studies)] * 3,
        }
    )
    metadata = index.to_frame(index=False).assign(
        body_site=sites if sites is not None else ["Feces"] * count,
        condition="One condition",
    )
    return pairs, studies, evidence, metadata


def test_mixed_zero_wilcoxon_pins_wilcox_and_no_continuity_correction() -> None:
    # Eleven nonzero ranks 1..11; negative ranks 3+6+9=18. With three
    # exact zeros, original n=14 forces auto's asymptotic path. No tie
    # correction: mean rank sum=33, variance=11*12*23/24=126.5.
    values = [0.0] * 3 + [1.0, 2.0, -3.0, 4.0, 5.0, -6.0, 7.0, 8.0, -9.0, 10.0, 11.0]
    row = build_tables(*_cell_inputs(values), [SHIFT]).stool_disease_shifts.iloc[0]
    assert row["wilcoxon_statistic"] == 18.0
    assert row["wilcoxon_p"] == pytest.approx(math.erfc(15 / math.sqrt(253)), rel=1e-13)
    assert row["wilcoxon_nonzero_comparisons"] == 11
    assert row["wilcoxon_zero_comparisons"] == 3
    assert row["wilcoxon_zero_method"] == "wilcox"
    assert row["wilcoxon_method"] == "asymptotic"
    assert row["wilcoxon_method_requested"] == "auto"
    assert not row["wilcoxon_correction"]
    assert row["wilcoxon_method_contract"] == REVIEWED_SCIPY_CONTRACTS[SCIPY_VERSION]
    assert row["wilcoxon_scipy_version"] == SCIPY_VERSION


def test_locked_environment_records_its_reviewed_wilcoxon_contract() -> None:
    """Exercise the installed SciPy; changing only a version constant proves nothing."""
    assert SCIPY_VERSION in REVIEWED_SCIPY_CONTRACTS
    result = build_tables(*_cell_inputs([1.0, -2.0, 3.0]), [SHIFT])
    for table in (result.stool_disease_shifts, result.stool_study_sensitivity):
        row = table.iloc[0]
        assert row["wilcoxon_statistic"] == 2.0
        assert row["wilcoxon_p"] == 0.75
        assert row["wilcoxon_method"] == "exact"
        assert row["wilcoxon_scipy_version"] == SCIPY_VERSION
        assert row["wilcoxon_method_contract"] == REVIEWED_SCIPY_CONTRACTS[SCIPY_VERSION]


def _all_positive_asymptotic_p(nonzero: int) -> float:
    """Independent normal-tail oracle for untied ranks 1..n without correction."""
    mean = nonzero * (nonzero + 1) / 4
    variance = nonzero * (nonzero + 1) * (2 * nonzero + 1) / 24
    return math.erfc(mean / math.sqrt(2 * variance))


@pytest.mark.parametrize(
    ("values", "expected_method", "expected_statistic", "expected_p"),
    [
        ([1.0, -2.0, 3.0], "exact", 2.0, 0.75),
        ([0.0, 1.0, -2.0], "permutation", 1.0, 1.0),
        ([1.0, -1.0, 2.0], "permutation", 1.5, 0.75),
        ([0.0, 1e-13, -1e-13], "permutation", 1.5, 1.0),
        ([0.0, *map(float, range(1, 13))], "permutation", 0.0, 2 / 2**12),
        (
            [0.0, *map(float, range(1, 14))],
            "asymptotic",
            0.0,
            _all_positive_asymptotic_p(13),
        ),
        (list(map(float, range(1, 51))), "exact", 0.0, 2 / 2**50),
        (list(map(float, range(1, 52))), "asymptotic", 0.0, _all_positive_asymptotic_p(51)),
        ([0.0, *map(float, range(1, 51))], "asymptotic", 0.0, _all_positive_asymptotic_p(50)),
        ([0.0, 0.0], "not_tested", math.nan, math.nan),
    ],
)
def test_resolved_auto_method_uses_original_size_zeros_and_absolute_ties(
    values: list[float], expected_method: str, expected_statistic: float, expected_p: float
) -> None:
    result = build_tables(*_cell_inputs(values), [SHIFT])
    row = result.stool_disease_shifts.iloc[0]
    assert row["wilcoxon_method"] == expected_method
    assert row["wilcoxon_method_requested"] == "auto"
    assert row["wilcoxon_zero_comparisons"] == values.count(0.0)
    assert row["wilcoxon_nonzero_comparisons"] == len(values) - values.count(0.0)
    if expected_method == "not_tested":
        assert row["wilcoxon_status"] == "all_tied"
        assert pd.isna(row["wilcoxon_statistic"])
        assert pd.isna(row["wilcoxon_p"])
    else:
        assert row["wilcoxon_statistic"] == expected_statistic
        assert row["wilcoxon_p"] == pytest.approx(expected_p, rel=1e-13, abs=0.0)
        # Reporting the resolved method must not alter the historical auto result.
        expected = stats.wilcoxon(
            values, zero_method="wilcox", correction=False, alternative="two-sided", method="auto"
        )
        assert row["wilcoxon_statistic"] == expected.statistic
        assert row["wilcoxon_p"] == expected.pvalue


def test_quadrant_labels_come_from_verified_evidence_even_with_empty_stool() -> None:
    inputs = _cell_inputs([1.0, -2.0], sites=["Saliva", "Saliva"])
    result = build_tables(*inputs, [SHIFT])
    assert result.quadrant_summary.fa_label.tolist() == ["verified facultative"] * 2
    assert result.quadrant_summary.oa_label.tolist() == ["verified obligate"] * 2


def test_row_three_counts_every_trait_cell_against_exact_raw_stool_metadata(inputs: Inputs) -> None:
    pairs, studies, evidence, metadata = inputs
    metadata = metadata.copy()
    metadata.loc[7, ["body_site", "condition"]] = ["Feces,Blood", "Disease A"]
    result = build_tables(pairs, studies, evidence, metadata, [SHIFT, ZERO])
    exact = metadata.loc[metadata["body_site"].eq("Feces")]
    expected_counts = exact.groupby("condition").size().to_dict()
    assert expected_counts == {"Disease A": 1, "Disease B": 3, "Disease C": 1}
    for row in result.stool_disease_shifts.itertuples(index=False):
        assert row.n_comparisons == expected_counts[row.condition]
        rows = exact.loc[exact["condition"].eq(row.condition)]
        assert row.n_studies == rows["study_id"].nunique()
    assert result.stool_disease_shifts.groupby("trait")["n_comparisons"].sum().to_dict() == {
        SHIFT: 5,
        ZERO: 5,
    }


def test_study_sensitivity_averages_within_cell_then_weights_studies_equally() -> None:
    inputs = _cell_inputs(
        [1.0, 2.0, 3.0, -10.0, 100.0],
        study_ids=["s1", "s1", "s1", "s2", "s1"],
        sites=["Feces", "Feces", "Feces", "Feces", "Oral"],
    )
    result = build_tables(*inputs, [SHIFT])
    comparison = result.stool_disease_shifts.iloc[0]
    study = result.stool_study_sensitivity.iloc[0]
    assert comparison["mean_delta_z"] == -1.0
    assert comparison["wilcoxon_statistic"] == 4.0
    assert comparison["wilcoxon_p"] == 0.875
    assert study["mean_delta_z"] == -4.0  # equally weighted within-cell means [2, -10]
    assert study["n_comparisons"] == 4
    assert study["n_studies"] == 2
    assert study["wilcoxon_statistic"] == 1.0
    assert study["wilcoxon_p"] == 1.0
    assert study["wilcoxon_nonzero_studies"] == 2
    assert study["wilcoxon_zero_studies"] == 0
    assert study["wilcoxon_method"] == "exact"
    assert study["test_unit"] == "study"
    assert study["mean_weighting"] == "equal_weight_within_cell_study_means"
    assert study["p_adjustment"] == "none_nominal"
    assert list(result.stool_study_sensitivity.condition) == list(
        result.stool_disease_shifts.condition
    )
    assert list(result.stool_study_sensitivity.trait) == list(result.stool_disease_shifts.trait)
    assert list(result.stool_study_sensitivity.disease_rank) == list(
        result.stool_disease_shifts.disease_rank
    )


def test_study_sensitivity_retains_within_study_cancellation_as_an_exact_zero() -> None:
    result = build_tables(*_cell_inputs([1.0, -1.0, 4.0], study_ids=["s1", "s1", "s2"]), [SHIFT])
    row = result.stool_study_sensitivity.iloc[0]
    assert row["mean_delta_z"] == 2.0
    assert row["n_comparisons"] == 3
    assert row["n_studies"] == 2
    assert row["wilcoxon_nonzero_studies"] == 1
    assert row["wilcoxon_zero_studies"] == 1
    assert row["wilcoxon_method"] == "permutation"
    assert row["wilcoxon_p"] == 1.0


def test_study_sensitivity_exact_zero_is_distinct_from_activity_threshold() -> None:
    inputs = _cell_inputs([0.0, 0.5e-12, 1e-12, 2e-12])
    before = tuple(frame.copy(deep=True) for frame in inputs)
    result = build_tables(*inputs, [SHIFT])
    activity = result.trait_summary.set_index("trait").loc[SHIFT]
    sensitivity = result.stool_study_sensitivity.iloc[0]
    assert activity["sdz_active_studies"] == 1
    assert activity["sdz_up"] == 1
    assert activity["sdz_down"] == 0
    assert activity["sdz_tied"] == 3
    assert sensitivity["n_studies"] == sensitivity["n_comparisons"] == 4
    assert sensitivity["wilcoxon_nonzero_studies"] == 3, "Exact-zero study support must be retained"
    assert sensitivity["wilcoxon_zero_studies"] == 1
    assert sensitivity["wilcoxon_statistic"] == 0.0
    # All three nonzero ranks have the same sign: two extreme assignments / 2**3.
    assert sensitivity["wilcoxon_p"] == 0.25
    assert sensitivity["wilcoxon_method"] == "permutation"
    assert sensitivity["wilcoxon_method_contract"] == REVIEWED_SCIPY_CONTRACTS[SCIPY_VERSION]
    assert sensitivity["mean_delta_z"] == pytest.approx(0.875e-12, rel=1e-15, abs=0)
    assert result.stool_disease_shifts.iloc[0]["wilcoxon_nonzero_comparisons"] == 3
    for original, after in zip(before, inputs, strict=True):
        pd.testing.assert_frame_equal(original, after)


@pytest.mark.parametrize("sign", [-1.0, 1.0])
def test_study_sensitivity_preserves_signed_nextafter_activity_boundaries(sign: float) -> None:
    values = [
        0.0,
        sign * np.nextafter(1e-12, 0.0),
        sign * 1e-12,
        sign * np.nextafter(1e-12, np.inf),
    ]
    result = build_tables(*_cell_inputs(values), [SHIFT])
    activity = result.trait_summary.set_index("trait").loc[SHIFT]
    sensitivity = result.stool_study_sensitivity.iloc[0]
    assert activity["sdz_up"] == int(sign > 0)
    assert activity["sdz_down"] == int(sign < 0)
    assert activity["sdz_active_studies"] == 1
    assert activity["sdz_tied"] == 3
    assert sensitivity["wilcoxon_nonzero_studies"] == 3
    assert sensitivity["wilcoxon_zero_studies"] == 1
    assert sensitivity["wilcoxon_statistic"] == 0.0
    assert sensitivity["wilcoxon_p"] == 0.25
    assert sensitivity["wilcoxon_method"] == "permutation"


def test_exact_zero_contract_rejects_in_memory_tolerance_unification(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Exercise the plausible study-only dead-band change without editing source files."""
    original = tables_module._wilcoxon
    changed: list[int] = []

    def unified(
        values: NDArray[np.float64], *, unit: Literal["comparison", "study"] = "comparison"
    ) -> dict[str, object]:
        if unit == "study":
            values = np.where(np.abs(values) <= 1e-12, 0.0, values)
            changed.append(int(np.count_nonzero(values)))
        return original(values, unit=unit)

    monkeypatch.setattr(tables_module, "_wilcoxon", unified)
    with pytest.raises(AssertionError, match="Exact-zero study support must be retained"):
        test_study_sensitivity_exact_zero_is_distinct_from_activity_threshold()
    assert changed == [1]


def test_study_sensitivity_empty_population_keeps_its_schema() -> None:
    result = build_tables(*_cell_inputs([1.0, -2.0], sites=["Saliva", "Saliva"]), [SHIFT])
    assert result.stool_study_sensitivity.empty
    assert {
        "trait",
        "n_studies",
        "n_comparisons",
        "wilcoxon_nonzero_studies",
        "test_unit",
    }.issubset(result.stool_study_sensitivity.columns)


@pytest.mark.parametrize("version", ["1.18.0", "1.19.0", "1.18.1.dev0"])
def test_auto_method_reporting_rejects_an_unreviewed_scipy_version(
    monkeypatch: pytest.MonkeyPatch,
    version: str,
) -> None:
    monkeypatch.setattr(tables_module, "SCIPY_VERSION", version)
    with pytest.raises(ValueError, match=r"Wilcoxon auto-method contract requires .*SciPy"):
        build_tables(*_cell_inputs([1.0, -2.0]), [SHIFT])
