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
from mtp.three_test_table.native import native_tables
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
