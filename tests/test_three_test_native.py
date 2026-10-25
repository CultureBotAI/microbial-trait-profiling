"""Issue #1044: import complete native families without changing their estimands."""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control

from mtp.multiplicity.families import rows
from mtp.multiplicity.sources import (
    VerifiedSources,
    bind_recorded_hashes,
    recorded_arm_output,
    recorded_comparison_output,
)
from mtp.paths import file_sha256
from mtp.three_test_table.composition import composition_table
from mtp.three_test_table.native import native_tables, validate_native_overlap
from multiplicity_fixture import fixture_sources, write_table


def _native_path(sources: VerifiedSources) -> Path:
    return recorded_comparison_output(
        sources.comparison, "trait_profile_paired_contrast_trait_stats_pooled.tsv"
    )


def _noel_path(sources: VerifiedSources) -> Path:
    return recorded_arm_output(
        sources.arms["noel"], "directional_screen", "trait_direction_screen.tsv"
    )


@pytest.mark.parametrize("subset", [False, True])
def test_native_import_preserves_observed_values_and_full_catalog_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, subset: bool
) -> None:
    sources = fixture_sources(tmp_path, monkeypatch, subset=subset)
    noel = {row["trait"]: row for row in rows(_noel_path(sources))}
    native = {row["trait"]: row for row in rows(_native_path(sources))}
    frame, consumed = native_tables(sources)
    assert frame["trait"].tolist() == list(sources.features["count_profiles"])
    assert frame["trait"].is_unique
    assert {_noel_path(sources), _native_path(sources)}.issubset(consumed)
    assert tuple(sorted(set(consumed))) == consumed
    assert all(path.is_file() for path in consumed)
    result = frame.set_index("trait")
    for trait, row in noel.items():
        for field in ("studies_up", "studies_down", "studies_tied", "studies_scored"):
            assert result.loc[trait, f"noel_{field}"] == int(row[field])
        assert result.loc[trait, "noel_status"] == row["test_status"]
        assert result.loc[trait, "noel_family_size"] == len(noel)
        for target, source in (
            ("share_up", "share_up"),
            ("ci_low", "ci_low"),
            ("ci_high", "ci_high"),
            ("p_value", "binomial_p"),
            ("q_by_nominal", "by_q_nominal"),
            ("q_bh_nominal", "bh_q_nominal"),
        ):
            value = result.loc[trait, f"noel_{target}"]
            assert value == float(row[source]) if row[source] else pd.isna(value)
    for trait, row in native.items():
        assert result.loc[trait, "native_selected_status"] == row["study_test_status"]
        assert result.loc[trait, "native_selected_family_size"] == len(native)
        for target, source in (
            ("study_count", "study_count"),
            ("study_mean_dz", "study_mean_increased_minus_decreased"),
            ("study_ci95_lower", "study_ci95_lower"),
            ("study_ci95_upper", "study_ci95_upper"),
            ("study_t_statistic", "study_t_statistic"),
            ("study_p_value", "study_p_value"),
            ("study_q_bh", "study_q_value"),
            ("experiment_count", "pair_count"),
            ("experiment_mean_dz", "mean_increased_minus_decreased"),
        ):
            assert result.loc[trait, f"native_selected_{target}"] == float(row[source])
    if not subset:
        extra = result.loc["qualified|extra"]
        assert extra["native_selected_status"] == "not_selected"
        assert pd.isna(extra["native_selected_study_p_value"])
        assert pd.isna(extra["native_selected_study_q_by"])
        assert pd.isna(extra["native_selected_study_mean_dz"])
        assert extra["noel_status"] in {"tested", "all_tied"}


def test_native_by_uses_selected_family_not_the_retained_catalog(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = fixture_sources(tmp_path, monkeypatch)
    records = list(rows(_native_path(sources)))
    # Isolate correction-family arithmetic with controlled imported probabilities.
    # The preservation tests above use the unchanged actual producer outputs.
    p_values = [0.001 * (index + 1) for index in range(len(records))]
    native_bh = false_discovery_control(p_values, method="bh")
    for row, p_value, q_value in zip(records, p_values, native_bh, strict=True):
        row["study_p_value"] = str(p_value)
        row["study_q_value"] = str(q_value)
    write_table(_native_path(sources), records)
    expected = false_discovery_control(p_values, method="by")
    wrong_family = false_discovery_control([*p_values, 1.0], method="by")[:-1]
    assert np.any(expected != wrong_family), "Fixture must distinguish the two families"
    table, _ = native_tables(sources)
    result = table.set_index("trait")
    actual = [result.loc[row["trait"], "native_selected_study_q_by"] for row in records]
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("producer", ["noel", "native"])
@pytest.mark.parametrize("mutation", ["truncate", "duplicate", "replace_member"])
def test_native_import_rejects_incomplete_or_duplicate_observed_membership(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, producer: str, mutation: str
) -> None:
    sources = fixture_sources(tmp_path, monkeypatch)
    path = _noel_path(sources) if producer == "noel" else _native_path(sources)
    records = list(rows(path))
    if mutation == "truncate":
        records.pop()
    elif mutation == "duplicate":
        records.append(records[0].copy())
    else:
        records[-1] = records[0].copy()
    write_table(path, records)
    with pytest.raises(ValueError, match=r"(membership|duplicate|inventory)"):
        native_tables(sources)


def test_native_import_matches_traits_by_identity_not_export_row_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = fixture_sources(tmp_path, monkeypatch)
    expected, _ = native_tables(sources)
    for path in (_noel_path(sources), _native_path(sources)):
        write_table(path, list(reversed(list(rows(path)))))
    actual, _ = native_tables(sources)
    pd.testing.assert_frame_equal(actual, expected)


def test_native_import_rejects_profile_order_different_from_bound_feature_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = fixture_sources(tmp_path, monkeypatch)
    changed = replace(
        sources,
        features={
            **sources.features,
            "count_profiles": tuple(reversed(sources.features["count_profiles"])),
        },
    )
    with pytest.raises(ValueError, match="retained feature inventory"):
        native_tables(changed)


def test_native_import_checks_upstream_recorded_table_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = fixture_sources(tmp_path, monkeypatch)
    bound = bind_recorded_hashes(sources, file_sha256(Path(sources.comparison.manifest_path)))
    table = _noel_path(sources)
    records = list(rows(table))
    records[0]["share_up"] = "0.123456"
    write_table(table, records)
    with pytest.raises(ValueError, match=r"Source.*(record|hash|bytes)"):
        native_tables(bound)


def _overlap_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, subset: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sources = fixture_sources(tmp_path, monkeypatch, subset=subset)
    profiles = pd.read_csv(sources.product("count_profiles"), sep="\t")
    composition = composition_table(profiles, sources.features["count_profiles"]).table
    native, _ = native_tables(sources)
    return composition, native


@pytest.mark.parametrize("subset", [False, True])
def test_native_overlap_matches_actual_producer_exports(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, subset: bool
) -> None:
    composition, native = _overlap_inputs(tmp_path, monkeypatch, subset=subset)
    checked = validate_native_overlap(composition, native)
    selected = checked.native_selected_status != "not_selected"
    assert selected.any()
    assert checked.loc[selected, "native_selected_overlap_status"].eq("tested_match").all()
    assert checked.loc[~selected, "native_selected_overlap_status"].eq("not_selected").all()
    assert "native_selected_overlap_status" not in native


@pytest.mark.parametrize(
    "column",
    [
        "native_selected_study_t_statistic",
        "native_selected_study_p_value",
        "native_selected_study_mean_dz",
        "native_selected_study_ci95_lower",
        "native_selected_study_ci95_upper",
        "native_selected_experiment_mean_dz",
        "native_selected_study_sd",
        "native_selected_study_count",
        "native_selected_experiment_count",
        "native_selected_profile_support",
    ],
)
def test_native_overlap_rejects_changed_raw_statistics_or_population(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, column: str
) -> None:
    composition, native = _overlap_inputs(tmp_path, monkeypatch)
    index = native.index[native.native_selected_status == "tested"][0]
    native.loc[index, column] += 1
    with pytest.raises(ValueError, match="Native overlap"):
        validate_native_overlap(composition, native)


def test_native_overlap_never_compares_q_values_from_different_families(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    composition, native = _overlap_inputs(tmp_path, monkeypatch)
    native["native_selected_study_q_bh"] = 0.123
    native["native_selected_study_q_by"] = 0.987
    checked = validate_native_overlap(composition.iloc[::-1], native)
    assert checked.trait.tolist() == native.trait.tolist()
    assert checked.native_selected_study_q_bh.eq(0.123).all()
    assert checked.native_selected_study_q_by.eq(0.987).all()


def test_native_overlap_preserves_relative_precision_for_tiny_p_values(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    composition, native = _overlap_inputs(tmp_path, monkeypatch)
    index = native.index[native.native_selected_status == "tested"][0]
    trait = native.loc[index, "trait"]
    composition.loc[composition.trait == trait, "study_t_p"] = 1e-100
    native.loc[index, "native_selected_study_p_value"] = 1e-100
    validate_native_overlap(composition, native)
    native.loc[index, "native_selected_study_p_value"] = 1e-101
    with pytest.raises(ValueError, match="Native overlap"):
        validate_native_overlap(composition, native)


def test_native_overlap_exposes_legitimate_dispersion_threshold_difference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    composition, native = _overlap_inputs(tmp_path, monkeypatch)
    index = native.index[native.native_selected_status == "tested"][0]
    trait = native.loc[index, "trait"]
    composition.loc[composition.trait == trait, "study_delta_z_sd"] = 5e-9
    native.loc[index, "native_selected_study_sd"] = 5e-9
    native.loc[index, "native_selected_status"] = "insufficient_variance"
    for column in ("study_t_statistic", "study_p_value", "study_ci95_lower", "study_ci95_upper"):
        native.loc[index, f"native_selected_{column}"] = np.nan
    checked = validate_native_overlap(composition, native)
    assert checked.loc[index, "native_selected_overlap_status"] == "native_sd_threshold_difference"
    assert "1e-12" in checked.loc[index, "native_selected_overlap_reason"]
    assert "1e-8" in checked.loc[index, "native_selected_overlap_reason"]
    assert composition.loc[composition.trait == trait, "study_t_status"].item() == "tested"


def test_native_overlap_rejects_unexplained_status_difference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    composition, native = _overlap_inputs(tmp_path, monkeypatch)
    index = native.index[native.native_selected_status == "tested"][0]
    native.loc[index, "native_selected_status"] = "insufficient_variance"
    with pytest.raises(ValueError, match="Native overlap"):
        validate_native_overlap(composition, native)


@pytest.mark.parametrize("reason", ["all_tied", "insufficient_variance", "insufficient_support"])
def test_native_overlap_preserves_matching_untestable_descriptions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, reason: str
) -> None:
    composition, native = _overlap_inputs(tmp_path, monkeypatch)
    index = native.index[native.native_selected_status == "tested"][0]
    selected = composition.trait == native.loc[index, "trait"]
    composition.loc[selected, "study_t_status"] = "untestable"
    composition.loc[selected, "study_t_reason"] = reason
    composition.loc[selected, "study_delta_z_sd"] = 0.0
    native.loc[index, "native_selected_status"] = reason
    native.loc[index, "native_selected_study_sd"] = 0.0
    if reason == "all_tied":
        composition.loc[selected, ["study_mean_delta_z", "experiment_mean_delta_z"]] = 0.0
        native.loc[
            index, ["native_selected_study_mean_dz", "native_selected_experiment_mean_dz"]
        ] = 0.0
    elif reason == "insufficient_support":
        composition.loc[selected, "composition_studies"] = 1
        composition.loc[selected, "study_delta_z_sd"] = np.nan
        native.loc[index, "native_selected_study_count"] = 1
    for full_column, native_column in (
        ("study_t_statistic", "study_t_statistic"),
        ("study_t_p", "study_p_value"),
        ("study_mean_delta_z_ci95_low", "study_ci95_lower"),
        ("study_mean_delta_z_ci95_high", "study_ci95_upper"),
    ):
        composition.loc[selected, full_column] = np.nan
        native.loc[index, f"native_selected_{native_column}"] = np.nan
    checked = validate_native_overlap(composition, native)
    assert checked.loc[index, "native_selected_overlap_status"] == "untestable_match"
    assert checked.loc[index, "native_selected_overlap_reason"] == reason
    assert pd.isna(checked.loc[index, "native_selected_study_p_value"])


@pytest.mark.parametrize("mutation", ["duplicate", "missing", "unexpected"])
def test_native_overlap_rejects_ambiguous_trait_membership(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mutation: str
) -> None:
    composition, native = _overlap_inputs(tmp_path, monkeypatch)
    if mutation == "duplicate":
        native = pd.concat([native, native.iloc[[0]]], ignore_index=True)
    elif mutation == "missing":
        native = native.iloc[:-1]
    else:
        native.loc[0, "trait"] = "unexpected|trait"
    with pytest.raises(ValueError, match="Native overlap"):
        validate_native_overlap(composition, native)
