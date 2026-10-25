"""Full-catalog composition and secondary signed-rank contracts for issue #1044."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from mtp.composition_calibration.adapters.jed_composition import JedCompositionAdapter
from mtp.three_test_table.composition import composition_table


@pytest.fixture
def profiles() -> pd.DataFrame:
    """Unequal experiments per study, an absent trait, and a MediaDive-only pair."""
    records: list[dict[str, object]] = []
    for study, experiment, increased, decreased in (
        ("s1", "e1", 9, 1),
        ("s1", "e2", 8, 2),
        ("s1", "e3", 7, 3),
        ("s2", "e4", 2, 8),
        ("s3", "e5", 6, 4),
        ("s4", "e6", 0, 0),
    ):
        for direction, value in (("inc", increased), ("dec", decreased)):
            records.append(
                {
                    "study_id": study,
                    "experiment_id": experiment,
                    "inc_dec": direction,
                    "trait_A": value,
                    "trait_B": 10 - value if experiment != "e6" else 0,
                    "trait_zero": 0,
                    "mediadive.medium:1": 1000,
                }
            )
    return pd.DataFrame(records)


def test_native_observed_t_is_reused_with_full_family_and_study_weights(
    profiles: pd.DataFrame,
) -> None:
    traits = ["trait_A", "trait_B", "trait_zero"]
    result = composition_table(profiles, traits)
    table = result.table.set_index("trait")
    adapter = JedCompositionAdapter()
    for correction in ("bh", "by"):
        native = adapter.run_inference(
            profiles, pd.DataFrame(), {"traits": traits, "correction": correction}
        ).results_df.set_index("trait")
        for trait in traits[:2]:
            assert table.loc[trait, "study_t_p"] == native.loc[trait, "p_value"]
            assert table.loc[trait, "study_t_statistic"] == native.loc[trait, "t_statistic"]
            assert table.loc[trait, f"study_t_q_{correction}"] == native.loc[trait, "q_value"]
            assert (
                table.loc[trait, "study_mean_delta_z_ci95_low"] == native.loc[trait, "ci95_lower"]
            )
            assert (
                table.loc[trait, "study_mean_delta_z_ci95_high"] == native.loc[trait, "ci95_upper"]
            )
    assert set(table["composition_family_size"]) == {3}
    assert set(table["composition_studies"]) == {4}
    assert set(table["composition_pairs"]) == {6}
    assert set(table["composition_zero_total_profiles"]) == {2}
    assert table.loc["trait_A", "experiment_mean_delta_z"] != pytest.approx(
        table.loc["trait_A", "study_mean_delta_z"]
    )
    pd.testing.assert_frame_equal(
        result.study_contrasts,
        result.pair_contrasts.groupby(level="study_id").mean(),
    )


def test_transform_denominator_matches_declared_deck_formula(profiles: pd.DataFrame) -> None:
    traits = ["trait_A", "trait_B", "trait_zero"]
    result = composition_table(profiles, traits)
    counts = profiles[traits]
    fractions = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    logs = np.log1p(1000 * fractions)
    zscores = (logs - logs.mean()) / logs.std(ddof=0).replace(0, np.nan)
    expected = (
        zscores.fillna(0)
        .assign(
            study_id=profiles.study_id,
            experiment_id=profiles.experiment_id,
            inc_dec=profiles.inc_dec,
        )
        .set_index(["study_id", "experiment_id", "inc_dec"])
    )
    expected = expected.xs("inc", level="inc_dec") - expected.xs("dec", level="inc_dec")
    pd.testing.assert_frame_equal(result.pair_contrasts, expected, atol=1e-14, rtol=1e-14)


def test_secondary_wilcoxon_counts_studies_and_uses_full_catalog(profiles: pd.DataFrame) -> None:
    result = composition_table(profiles, ["trait_A", "trait_B", "trait_zero"])
    table = result.table.set_index("trait")
    raw = []
    for trait in ("trait_A", "trait_B"):
        values = result.study_contrasts[trait]
        expected = stats.wilcoxon(
            values, zero_method="wilcox", correction=False, alternative="two-sided", method="approx"
        )
        raw.append(expected.pvalue)
        assert table.loc[trait, "study_wilcoxon_statistic"] == expected.statistic
        assert table.loc[trait, "study_wilcoxon_p"] == expected.pvalue
        assert table.loc[trait, "study_wilcoxon_zero_studies"] == 1
        assert table.loc[trait, "study_wilcoxon_nonzero_studies"] == 3
    expected_by = stats.false_discovery_control([*raw, 1.0], method="by")
    assert table.loc["trait_A", "study_wilcoxon_q_by"] == expected_by[0]
    assert table.loc["trait_zero", "study_wilcoxon_reason"] == "all_tied"
    assert pd.isna(table.loc["trait_zero", "study_wilcoxon_q_by"])
    assert pd.isna(table.loc["trait_zero", "study_t_p"])
    assert table.loc["trait_zero", "study_t_reason"] == "all_tied"
    assert table.loc["trait_zero", "study_mean_delta_z"] == 0


def test_rank_ties_and_reordering_are_explicit(profiles: pd.DataFrame) -> None:
    profiles["trait_twin"] = profiles.trait_A
    traits = ["trait_A", "trait_twin", "trait_B", "trait_zero"]
    expected = composition_table(profiles, traits).table.set_index("trait")
    observed = composition_table(profiles.iloc[::-1], traits[::-1]).table.set_index("trait")
    for prefix in ("experiment_mean_delta_z", "study_mean_delta_z"):
        for order in ("ascending", "descending", "absolute_descending"):
            column = f"{prefix}_rank_{order}"
            assert expected.loc["trait_A", column] == expected.loc["trait_twin", column]
            assert expected.loc["trait_zero", f"{prefix}_rank_absolute_descending"] == 4
            pd.testing.assert_series_equal(
                expected[column].sort_index(), observed[column].sort_index()
            )


def test_unpaired_and_single_study_preserve_rows_and_undefined_statistics(
    profiles: pd.DataFrame,
) -> None:
    selected = profiles.loc[profiles.study_id.isin(["s1", "s2"])].copy()
    selected = selected.loc[~((selected.study_id == "s2") & (selected.inc_dec == "dec"))]
    result = composition_table(selected, ["trait_A", "trait_B", "trait_zero"])
    assert len(result.table) == 3
    assert set(result.table.composition_studies) == {1}
    assert set(result.table.composition_unpaired_profiles_excluded) == {1}
    assert set(result.table.study_t_reason) == {"insufficient_support"}
    assert result.table.study_t_p.isna().all()
    assert result.table.study_mean_delta_z_ci95_low.isna().all()
    assert result.table.experiment_mean_delta_z.notna().all()


def test_no_pairs_retains_entire_untestable_catalog(profiles: pd.DataFrame) -> None:
    result = composition_table(
        profiles.loc[profiles.inc_dec == "inc"], ["trait_A", "trait_B", "trait_zero"]
    )
    assert len(result.table) == 3
    assert set(result.table.composition_studies) == {0}
    assert result.table.study_mean_delta_z.isna().all()
    assert result.table.study_t_p.isna().all()
    assert result.table.study_wilcoxon_p.isna().all()
    assert result.table.study_mean_delta_z_rank_ascending.isna().all()


@pytest.mark.parametrize(
    "problem", ["duplicate", "blank_id", "missing_id", "bad_direction", "alias_duplicate"]
)
def test_invalid_profile_identity_fails(profiles: pd.DataFrame, problem: str) -> None:
    if problem == "duplicate":
        profiles = pd.concat([profiles, profiles.iloc[[0]]], ignore_index=True)
    elif problem == "blank_id":
        profiles.loc[0, "study_id"] = " "
    elif problem == "missing_id":
        profiles.loc[0, "experiment_id"] = None
    elif problem == "bad_direction":
        profiles.loc[0, "inc_dec"] = "unknown"
    else:
        duplicate = profiles.iloc[[0]].copy()
        duplicate["inc_dec"] = "increased"
        profiles = pd.concat([profiles, duplicate], ignore_index=True)
    with pytest.raises(ValueError, match=r"identity|direction"):
        composition_table(profiles, ["trait_A", "trait_B", "trait_zero"])


@pytest.mark.parametrize("bad_count", [float("nan"), float("inf"), -1, 0.5, "bad", 2**53])
def test_invalid_counts_fail_before_adapter_coercion(
    profiles: pd.DataFrame, bad_count: float | str
) -> None:
    profiles["trait_A"] = profiles.trait_A.astype(object)
    profiles.loc[0, "trait_A"] = bad_count
    with pytest.raises(ValueError, match="count"):
        composition_table(profiles, ["trait_A", "trait_B", "trait_zero"])


@pytest.mark.parametrize("traits", [[], ["trait_A", "trait_A"], ["trait_A"], ["missing"]])
def test_full_retained_catalog_is_required(profiles: pd.DataFrame, traits: list[str]) -> None:
    with pytest.raises(ValueError, match=r"trait|catalog"):
        composition_table(profiles, traits)


def test_constant_nonzero_study_shifts_are_not_a_fabricated_t_discovery(
    profiles: pd.DataFrame,
) -> None:
    pair = profiles.loc[profiles.experiment_id == "e1"]
    repeated = pd.concat(
        [pair.assign(study_id=f"study{i}", experiment_id=f"pair{i}") for i in range(3)],
        ignore_index=True,
    )
    result = composition_table(repeated, ["trait_A", "trait_B", "trait_zero"])
    row = result.table.set_index("trait").loc["trait_A"]
    assert row["study_mean_delta_z"] > 0
    assert row["study_t_status"] == "untestable"
    assert row["study_t_reason"] == "insufficient_variance"
    assert pd.isna(row["study_t_p"])
    assert pd.isna(row["study_t_q_by"])
    assert pd.isna(row["study_mean_delta_z_ci95_low"])
    assert row["study_wilcoxon_status"] == "tested"


def test_duplicate_count_columns_fail(profiles: pd.DataFrame) -> None:
    duplicate = pd.concat([profiles, profiles[["trait_A"]]], axis=1)
    with pytest.raises(ValueError, match=r"[Dd]uplicate"):
        composition_table(duplicate, ["trait_A", "trait_B", "trait_zero"])


def test_unsafe_count_total_is_rejected(profiles: pd.DataFrame) -> None:
    profiles.loc[0, ["trait_A", "trait_B"]] = 2**52
    with pytest.raises(ValueError, match="count"):
        composition_table(profiles, ["trait_A", "trait_B", "trait_zero"])


def test_all_trait_paired_profile_support_and_primary_correction_are_explicit(
    profiles: pd.DataFrame,
) -> None:
    profiles["trait_ubiquitous"] = 1
    profiles["trait_rare"] = 0
    profiles.loc[0, "trait_rare"] = 1
    unpaired = profiles.iloc[[0]].assign(study_id="unpaired", experiment_id="unpaired")
    traits = ["trait_A", "trait_B", "trait_zero", "trait_ubiquitous", "trait_rare"]
    result = composition_table(pd.concat([profiles, unpaired], ignore_index=True), traits)
    table = result.table.set_index("trait")
    assert set(table.study_t_primary_correction) == {"BY"}
    assert table.composition_profile_support.to_dict() == {
        "trait_A": 10,
        "trait_B": 10,
        "trait_zero": 0,
        "trait_ubiquitous": 12,
        "trait_rare": 1,
    }
    np.testing.assert_array_equal(
        table.composition_profile_support_fraction, [10 / 12, 10 / 12, 0, 1, 1 / 12]
    )
    assert set(table.composition_profiles) == {12}
    assert set(table.composition_family_size) == {5}


def test_no_pair_support_fraction_is_undefined(profiles: pd.DataFrame) -> None:
    result = composition_table(
        profiles.loc[profiles.inc_dec == "inc"], ["trait_A", "trait_B", "trait_zero"]
    )
    assert result.table.composition_profile_support.eq(0).all()
    assert result.table.composition_profile_support_fraction.isna().all()


def test_full_family_by_is_distinct_from_bh_and_keeps_untestable_member() -> None:
    """Strong synthetic contrasts distinguish both corrections and the family size (#1077)."""
    records: list[dict[str, object]] = []
    for index in range(12):
        for direction, first, second in (
            ("inc", 55 + index % 5, 25 + (3 * index) % 7),
            ("dec", 10 + (2 * index) % 5, 10 + index % 4),
        ):
            records.append(
                {
                    "study_id": f"study{index}",
                    "experiment_id": "experiment",
                    "inc_dec": direction,
                    "trait_A": first,
                    "trait_B": second,
                    "trait_C": 100 - first - second,
                    "trait_zero": 0,
                }
            )
    tested = ["trait_A", "trait_B", "trait_C"]
    result = composition_table(pd.DataFrame.from_records(records), [*tested, "trait_zero"])
    table = result.table.set_index("trait")
    raw = [float(stats.ttest_1samp(result.study_contrasts[trait], 0.0).pvalue) for trait in tested]
    assert all(0 < value < 0.01 for value in raw)
    probabilities = np.array([*raw, 1.0])
    expected_by = stats.false_discovery_control(probabilities, method="by")
    expected_bh = stats.false_discovery_control(probabilities, method="bh")
    np.testing.assert_allclose(table.loc[tested, "study_t_p"], raw, rtol=1e-12, atol=0)
    np.testing.assert_allclose(
        table.loc[tested, "study_t_q_by"], expected_by[:-1], rtol=1e-12, atol=0
    )
    np.testing.assert_allclose(
        table.loc[tested, "study_t_q_bh"], expected_bh[:-1], rtol=1e-12, atol=0
    )
    assert np.all(expected_by[:-1] > expected_bh[:-1])
    assert np.all(expected_by[:-1] > stats.false_discovery_control(raw, method="by"))
    assert set(table.composition_family_size) == {4}
    assert table.loc["trait_zero", "study_t_status"] == "untestable"
    assert table.loc["trait_zero", "study_t_reason"] == "all_tied"
    assert table.loc["trait_zero", ["study_t_p", "study_t_q_by", "study_t_q_bh"]].isna().all()
