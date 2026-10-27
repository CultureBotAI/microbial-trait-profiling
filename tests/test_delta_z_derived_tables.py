"""Independent synthetic contracts for derived Delta-z tables (#1052, 2026-10-25)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from mtp.delta_z_derived.tables import build_tables, trait_summary

FA = "METPO:1000608"
OA = "METPO:1000607"
ANAEROBIC = "METPO:1000603"
SIGN = "TEST:threshold"
SHIFT = "TEST:shift"
ZERO = "TEST:zero"
TRAITS = [FA, OA, ANAEROBIC, SIGN, SHIFT, ZERO]
Inputs = tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]


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
    assert int(summary.sdz_active_studies.ge(50).sum()) == 2


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
