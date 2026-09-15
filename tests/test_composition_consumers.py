"""Downstream composition must use the same universe and keyed caches."""

import importlib
from pathlib import Path
from types import ModuleType

import numpy as np
import pandas as pd
import pytest

from bsdbng.trait_composition import (
    DENOMINATOR_DEFINITION,
    load_retained_profiles,
    trait_fractions,
    validate_cached_totals,
)


def consumer(name: str, monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    return importlib.import_module(name)


def fixture_counts() -> pd.DataFrame:
    return pd.DataFrame(
        {"A": [10, 20, 0, 0], "B": [90, 20, 0, 0], "C": [1, 0, 0, 0]},
        index=[f"s{i}|e{i}|inc" for i in range(4)],
    )


def write_profiles(path: Path, counts: pd.DataFrame) -> None:
    keys = pd.DataFrame(
        [key.split("|") for key in counts.index],
        columns=["study_id", "experiment_id", "inc_dec"],
        index=counts.index,
    )
    keys.join(counts).to_csv(path, sep="\t", index=False)


def test_cached_totals_reject_stale_and_incompatible_inputs(tmp_path: Path) -> None:
    counts = fixture_counts()
    counts["mediadive.medium:1"] = [900, 10, 700, 0]
    metadata = pd.DataFrame({"total_trait_signal": [101, 40, 0, 0]}, index=counts.index)
    validate_cached_totals(counts, metadata)
    path = tmp_path / "profiles.tsv"
    write_profiles(path, counts.iloc[::-1])
    pd.testing.assert_frame_equal(
        load_retained_profiles(path, metadata), counts.drop(columns="mediadive.medium:1")
    )
    with pytest.raises(ValueError, match="Stale or incompatible"):
        validate_cached_totals(counts, metadata.assign(total_trait_signal=counts.sum(axis=1)))
    with pytest.raises(ValueError, match="Stale or incompatible"):
        validate_cached_totals(counts[["A"]], metadata)
    with pytest.raises(ValueError, match="missing"):
        validate_cached_totals(counts.iloc[1:], metadata)
    with pytest.raises(ValueError, match="unique"):
        validate_cached_totals(pd.concat([counts, counts]), metadata)


def test_sankey_and_pd_use_all_retained_traits(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sankey = consumer("create_six_shared_trait_sankeys", monkeypatch)
    pd_analysis = consumer("create_pd_cellulose_analysis", monkeypatch)
    counts = fixture_counts()
    six = [trait.trait_id for trait in sankey.SHARED_TRAITS]
    counts = counts.rename(columns={"A": six[0], "B": six[1]})
    for trait in six[2:]:
        counts[trait] = 0
    metadata = pd.DataFrame({"total_trait_signal": counts.sum(axis=1)}, index=counts.index)
    path = tmp_path / "profiles.tsv"
    monkeypatch.setattr(sankey, "PROFILE_MATRIX_PATH", path)
    monkeypatch.setattr(pd_analysis, "PROFILE_MATRIX_PATH", path)
    scaled = np.log1p(trait_fractions(counts)[six] * 1000)
    expected = ((scaled - scaled.mean()) / scaled.std(ddof=0).replace(0, np.nan)).fillna(0)
    for media in ([900, 2, 700, 0], [0, 100000, 0, 999]):
        with_media = counts.assign(**{"mediadive.medium:1": media})
        write_profiles(path, with_media)
        loaded = sankey.load_raw_trait_counts(sankey.SHARED_TRAITS, metadata)
        pd.testing.assert_frame_equal(sankey.build_zscore_matrix(loaded, metadata), expected)
        pd.testing.assert_frame_equal(pd_analysis.build_zscore_matrix(metadata, six), expected)
    assert expected.iloc[2].equals(expected.iloc[3])


def test_module_selection_and_aggregation_are_media_invariant(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = consumer("create_condition_trait_module_heatmap", monkeypatch)
    counts = fixture_counts()
    labels = pd.DataFrame(
        dict.fromkeys(
            ["display_label", "column_label", "object_name", "object_category", "display_category"],
            counts.columns,
        ),
        index=counts.columns,
    )
    baseline = module.select_traits(counts, labels, 3, 1, 1.0, False)
    path = tmp_path / "profiles.tsv"
    write_profiles(path, counts)
    study, metadata = module.load_aggregated_study_profiles(path)
    for amount in [1, 1000]:
        media = counts.assign(**{"mediadive.medium:1": np.array([900, 2, 700, 0]) * amount})
        write_profiles(path, media)
        actual_study, actual_metadata = module.load_aggregated_study_profiles(path)
        pd.testing.assert_frame_equal(study, actual_study)
        pd.testing.assert_frame_equal(metadata, actual_metadata)
        selected, stats = module.select_traits(media, labels, 3, 1, 1.0, False)
        pd.testing.assert_index_equal(selected, baseline[0])
        pd.testing.assert_frame_equal(stats, baseline[1])
        pd.testing.assert_frame_equal(
            module.build_plot_matrix(media, selected), module.build_plot_matrix(counts, selected)
        )
        groups = metadata.assign(
            condition=["x", "y", "x", "y"],
            body_site="site",
            location="place",
            host_species="Homo sapiens",
        )
        expected_groups = module.build_condition_modules(
            study, groups, ["condition"], 1, 1, 0, 0, False, 10
        )
        actual_groups = module.build_condition_modules(
            study.assign(**{"mediadive.medium:1": [10000, 0, 90000, 2]}),
            groups,
            ["condition"],
            1,
            1,
            0,
            0,
            False,
            10,
        )
        for expected_group, actual_group in zip(expected_groups, actual_groups, strict=True):
            pd.testing.assert_frame_equal(expected_group, actual_group)


def test_followup_metadata_and_guild_fraction_are_media_invariant(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    followup = consumer("analyze_enrichment_infographic_followup", monkeypatch)
    monkeypatch.setattr(followup, "TARGET_CONDITIONS", ["condition"])
    monkeypatch.setattr(followup, "GUILD_TRAITS", {"guild": {"A": "A"}})
    counts = fixture_counts()
    metadata = pd.DataFrame(
        {
            "condition": "condition",
            "study_id": ["s0", "s1", "s2", "s3"],
            "experiment_id": ["e0", "e1", "e2", "e3"],
            "inc_dec": "inc",
        },
        index=counts.index,
    )
    expected = followup.summarize_guild_connections(counts, metadata)
    assert expected.iloc[0]["mean_relative_guild_signal"] == pytest.approx((10 / 101 + 0.5) / 4)
    with_media = counts.assign(**{"mediadive.medium:1": [900, 10, 700, 0]})
    pd.testing.assert_frame_equal(
        expected, followup.summarize_guild_connections(with_media, metadata)
    )
    experiments = metadata[["study_id", "experiment_id", "condition"]]
    studies = metadata[["study_id"]]
    path = tmp_path / "profiles.tsv"
    write_profiles(path, counts)
    before = followup.load_profiles(path, experiments, studies)
    write_profiles(path, with_media)
    after = followup.load_profiles(path, experiments, studies)
    pd.testing.assert_frame_equal(before[0], after[0])
    pd.testing.assert_frame_equal(before[1], after[1])
    assert before[2] == after[2]
    neighbor_metadata = after[1].copy()
    neighbor_metadata.loc[counts.index[1], "condition"] = "other"
    expected_neighbors = followup.compute_profile_neighbors(counts, neighbor_metadata, 2, 1)
    actual_neighbors = followup.compute_profile_neighbors(with_media, neighbor_metadata, 2, 1)
    for expected_table, actual_table in zip(expected_neighbors, actual_neighbors, strict=True):
        pd.testing.assert_frame_equal(expected_table, actual_table)


def test_guild_cache_aligns_keys_and_rejects_old_denominator(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    guild = consumer("create_condition_body_site_guild_maps", monkeypatch)
    metadata = pd.DataFrame(
        {
            "study_id": ["s1", "s2"],
            "experiment_id": ["e1", "e2"],
            "direction_label": ["increased", "decreased"],
            "denominator_definition": DENOMINATOR_DEFINITION,
            "total_trait_signal": [10, 10],
        },
        index=["s1|e1|inc", "s2|e2|dec"],
    )
    matrix = pd.DataFrame(
        {"Trait": [1.0, -1.0]}, index=["s2 | e2 | decreased", "s1 | e1 | increased"]
    )
    for setting, frame in [
        ("MATRIX_PATH", matrix),
        ("METADATA_PATH", metadata),
        (
            "TRAIT_METADATA_PATH",
            pd.DataFrame({"display_label": ["Trait"], "object_name": ["A"]}, index=["A"]),
        ),
    ]:
        path = tmp_path / f"{setting}.tsv"
        frame.to_csv(path, sep="\t")
        monkeypatch.setattr(guild, setting, path)
    raw_path = tmp_path / "profiles.tsv"
    write_profiles(raw_path, pd.DataFrame({"A": [1, 3], "B": [9, 7]}, index=metadata.index))
    monkeypatch.setattr(guild, "PROFILE_MATRIX_PATH", raw_path)
    actual, _ = guild.load_study_heatmap()
    assert actual.index.tolist() == metadata.index.tolist()
    assert actual["A"].tolist() == [-1.0, 1.0]
    (matrix * 2).to_csv(guild.MATRIX_PATH, sep="\t")
    with pytest.raises(ValueError, match="Stale transformed"):
        guild.load_study_heatmap()
    metadata.drop(columns="denominator_definition").to_csv(guild.METADATA_PATH, sep="\t")
    with pytest.raises(ValueError, match="Regenerate"):
        guild.load_study_heatmap()


def test_heatmap_comparison_preserves_native_q_values_and_missing_cells(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    audit = consumer("audit_composition_downstream", monkeypatch)
    report = tmp_path / "report"
    report.mkdir()
    monkeypatch.setattr(audit, "REPORT", report)
    for stage, q_value in [("before", 0.04), ("after", 0.08)]:
        root = tmp_path / stage
        directory = root / "six-shared-traits"
        directory.mkdir(parents=True)
        monkeypatch.setattr(audit, stage.upper(), root)
        manifest = []
        for suffix in ["trait_pattern_matrix", "direction_contrasts"]:
            traits = ["A", "B"] if stage == "before" else ["A"]
            cells = pd.DataFrame(
                {
                    "view_key": suffix,
                    "condition_body_site_pair": "condition | site",
                    "trait_id": traits,
                    "fisher_q_value": q_value,
                }
            )
            cells.to_csv(
                directory / f"six_shared_trait_condition_body_site_{suffix}.tsv",
                sep="\t",
                index=False,
            )
            manifest.append(
                {
                    "view_key": suffix,
                    "tested_cell_count": len(traits),
                    "displayed_cell_count": 1,
                    "bh_selected_cell_count": int(q_value < 0.05),
                }
            )
        pd.DataFrame(manifest).to_csv(
            directory / "six_shared_trait_cytoscape_heatmap_manifest.tsv", sep="\t", index=False
        )
    audit.compare_heatmap_families()
    cells = pd.read_csv(report / "heatmap_cells_before_after.tsv", sep="\t")
    matched = cells[cells["_merge"] == "both"]
    assert matched.fisher_q_value_before.eq(0.04).all()
    assert matched.fisher_q_value_after.eq(0.08).all()
    assert cells["_merge"].eq("left_only").sum() == 2
    families = pd.read_csv(report / "heatmap_family_comparison.tsv", sep="\t")
    assert families.tested_cell_count_before.eq(2).all()
    assert families.tested_cell_count_after.eq(1).all()
