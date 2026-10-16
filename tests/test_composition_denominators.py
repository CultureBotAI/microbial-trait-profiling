"""Unit tests for composition_denominators module (Refs #174, #202).

Verifies non-MediaDive denominator normalization, MediaDive column exclusion,
zero-signal handling, cache validation, and key alignment.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from mtp.composition_denominators import (
    DENOMINATOR_DEFINITION,
    load_retained_profiles,
    retained_trait_counts,
    trait_fractions,
    validate_cached_totals,
)

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"


def _load_script_module(module_name: str) -> ModuleType:
    if str(SCRIPTS_DIR) not in sys.path:
        sys.path.insert(0, str(SCRIPTS_DIR))
    module_path = SCRIPTS_DIR / f"{module_name}.py"
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def test_denominator_definition_constant() -> None:
    assert "retained non-MediaDive" in DENOMINATOR_DEFINITION


def test_retained_trait_counts_filters_mediadive() -> None:
    counts = pd.DataFrame(
        {
            "METPO:1000615": [1, 2],
            "METPO:2000517|mediadive.medium:1": [10, 20],
            "METPO:2000517|MediaDive.Medium:2": [30, 40],
            "METPO:1000616": [5, 6],
        },
        index=["p1", "p2"],
    )
    retained = retained_trait_counts(counts)
    assert list(retained.columns) == ["METPO:1000615", "METPO:1000616"]
    pd.testing.assert_frame_equal(
        retained,
        pd.DataFrame(
            {"METPO:1000615": [1, 2], "METPO:1000616": [5, 6]},
            index=["p1", "p2"],
        ),
    )


def test_trait_fractions_normalization() -> None:
    counts = pd.DataFrame(
        {
            "A": [10, 0, 0, 20],
            "B": [90, 0, 0, 80],
            "mediadive.medium:1": [900, 0, 700, 5],
        },
        index=["p1", "p2", "p3", "p4"],
    )
    fractions = trait_fractions(counts)
    assert list(fractions.columns) == ["A", "B"]
    assert fractions["A"].tolist() == pytest.approx([0.1, 0.0, 0.0, 0.2])
    assert fractions["B"].tolist() == pytest.approx([0.9, 0.0, 0.0, 0.8])
    assert fractions.iloc[1].eq(0.0).all()
    assert fractions.iloc[2].eq(0.0).all()
    assert not fractions.isna().any().any()


def test_trait_fractions_all_mediadive_matrix() -> None:
    only_media = pd.DataFrame(
        {"mediadive.medium:1": [10, 0], "mediadive.medium:2": [5, 20]},
        index=["p1", "p2"],
    )
    assert retained_trait_counts(only_media).shape == (2, 0)
    fractions = trait_fractions(only_media)
    assert fractions.shape == (2, 0)


def test_validate_cached_totals_success() -> None:
    counts = pd.DataFrame(
        {"A": [5, 10], "B": [15, 20], "mediadive.medium:1": [100, 200]},
        index=["s1|e1|inc", "s2|e2|dec"],
    )
    metadata = pd.DataFrame(
        {"total_trait_signal": [20, 30]},
        index=["s1|e1|inc", "s2|e2|dec"],
    )
    # Should pass without raising
    validate_cached_totals(counts, metadata)


def test_validate_cached_totals_detects_stale_cache() -> None:
    counts = pd.DataFrame(
        {"A": [5, 10], "B": [15, 20], "mediadive.medium:1": [100, 200]},
        index=["s1|e1|inc", "s2|e2|dec"],
    )
    stale_metadata = pd.DataFrame(
        {"total_trait_signal": [120, 230]},  # Includes MediaDive, stale!
        index=["s1|e1|inc", "s2|e2|dec"],
    )
    with pytest.raises(ValueError, match="Stale or incompatible non-MediaDive denominator cache"):
        validate_cached_totals(counts, stale_metadata)


def test_validate_cached_totals_rejects_duplicate_index() -> None:
    dup_counts = pd.DataFrame({"A": [1, 2]}, index=["dup", "dup"])
    metadata = pd.DataFrame({"total_trait_signal": [1, 2]}, index=["dup", "dup"])
    with pytest.raises(ValueError, match="Profile keys must be unique"):
        validate_cached_totals(dup_counts, metadata)


def test_validate_cached_totals_rejects_missing_keys() -> None:
    counts = pd.DataFrame({"A": [1, 2]}, index=["k1", "k2"])
    metadata = pd.DataFrame({"total_trait_signal": [1, 2]}, index=["k1", "missing_k"])
    with pytest.raises(ValueError, match="Cached metadata contains profile keys missing"):
        validate_cached_totals(counts, metadata)


def test_load_retained_profiles(tmp_path: Path) -> None:
    data = {
        "study_id": ["s1", "s2", "s3"],
        "experiment_id": ["e1", "e2", "e3"],
        "inc_dec": ["inc", "dec", "inc"],
        "trait_x": [10, 20, 0],
        "trait_y": [30, 40, 0],
        "mediadive.medium:1": [999, 888, 777],
    }
    df = pd.DataFrame(data)
    tsv_path = tmp_path / "profiles.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)

    metadata = pd.DataFrame(
        {"total_trait_signal": [40, 60]},
        index=["s1|e1|inc", "s2|e2|dec"],
    )

    loaded = load_retained_profiles(tsv_path, metadata)
    assert list(loaded.columns) == ["trait_x", "trait_y"]
    assert loaded.index.tolist() == ["s1|e1|inc", "s2|e2|dec"]
    assert loaded.loc["s1|e1|inc", "trait_x"] == 10
    assert loaded.loc["s2|e2|dec", "trait_y"] == 40


def test_load_retained_profiles_rejects_duplicate_keys(tmp_path: Path) -> None:
    data = {
        "study_id": ["s1", "s1", "s2"],
        "experiment_id": ["e1", "e1", "e2"],
        "inc_dec": ["inc", "inc", "dec"],
        "trait_x": [5, 5, 20],
        "trait_y": [15, 15, 40],
    }
    tsv_path = tmp_path / "dup_profiles.tsv"
    pd.DataFrame(data).to_csv(tsv_path, sep="\t", index=False)
    metadata = pd.DataFrame(
        {"total_trait_signal": [40, 60]},
        index=["s1|e1|inc", "s2|e2|dec"],
    )
    with pytest.raises(ValueError, match="Profile keys must be unique"):
        load_retained_profiles(tsv_path, metadata)


def test_study_heatmap_caller_excludes_mediadive_from_denominator() -> None:
    study_heatmap = _load_script_module("create_study_trait_profile_heatmap")

    selected = pd.Index(["METPO:1000608", "METPO:1000668"])
    base = pd.DataFrame(
        {
            "METPO:1000608": [10.0, 20.0, 30.0],
            "METPO:1000668": [30.0, 20.0, 10.0],
        },
        index=["p1", "p2", "p3"],
    )
    with_media = base.assign(**{"METPO:2000517|mediadive.medium:1": [900.0, 5.0, 500.0]})

    pd.testing.assert_frame_equal(
        study_heatmap.build_plot_matrix(base, selected),
        study_heatmap.build_plot_matrix(with_media, selected),
    )


def test_infographic_caller_excludes_mediadive_from_denominator() -> None:
    infographic = _load_script_module("create_trait_profile_infographic")

    selected = pd.Index(["METPO:1000608", "METPO:1000668"])
    base = pd.DataFrame(
        {
            "METPO:1000608": [10.0, 20.0, 30.0],
            "METPO:1000668": [30.0, 20.0, 10.0],
        },
        index=["p1", "p2", "p3"],
    )
    with_media = base.assign(**{"METPO:2000517|mediadive.medium:1": [900.0, 5.0, 500.0]})

    pd.testing.assert_frame_equal(
        infographic.build_zscore_matrix(base, selected),
        infographic.build_zscore_matrix(with_media, selected),
    )


def test_condition_module_heatmap_caller_excludes_mediadive() -> None:
    module_heatmap = _load_script_module("create_condition_trait_module_heatmap")

    selected = pd.Index(["METPO:1000608", "METPO:1000668"])
    base = pd.DataFrame(
        {
            "METPO:1000608": [10.0, 20.0, 30.0, 40.0],
            "METPO:1000668": [30.0, 20.0, 10.0, 0.0],
        },
        index=["s1|inc", "s2|inc", "s3|dec", "s4|dec"],
    )
    with_media = base.assign(**{"METPO:2000517|mediadive.medium:1": [900.0, 5.0, 500.0, 100.0]})

    pd.testing.assert_frame_equal(
        module_heatmap.build_plot_matrix(base, selected),
        module_heatmap.build_plot_matrix(with_media, selected),
    )


def test_enrichment_followup_caller_excludes_mediadive() -> None:
    followup = _load_script_module("analyze_enrichment_infographic_followup")

    target_cond = followup.TARGET_CONDITIONS[0]
    other_cond = followup.TARGET_CONDITIONS[1]
    guild_trait = next(iter(next(iter(followup.GUILD_TRAITS.values()))))
    metadata = pd.DataFrame(
        {
            "profile_id": ["p1", "p2", "p3"],
            "study_id": ["s1", "s2", "s3"],
            "experiment_id": ["e1", "e2", "e3"],
            "inc_dec": ["inc", "dec", "inc"],
            "direction_label": ["increased", "decreased", "increased"],
            "condition": [target_cond, target_cond, other_cond],
            "title": ["t1", "t2", "t3"],
            "experiment_label": ["l1", "l2", "l3"],
            "body_site": ["Feces", "Feces", "Feces"],
            "location": ["USA", "USA", "USA"],
            "profile_trait_count": [2, 2, 2],
            "total_trait_signal": [40.0, 40.0, 40.0],
        },
        index=["p1", "p2", "p3"],
    )
    base = pd.DataFrame(
        {
            guild_trait: [10.0, 20.0, 15.0],
            "METPO:1000668": [30.0, 20.0, 25.0],
        },
        index=["p1", "p2", "p3"],
    )
    with_media = base.assign(**{"METPO:2000517|mediadive.medium:1": [900.0, 5.0, 500.0]})

    neighbors_base, centroids_base, cross_base = followup.compute_profile_neighbors(
        base, metadata, top_n=2, min_condition_profiles=1
    )
    neighbors_media, centroids_media, cross_media = followup.compute_profile_neighbors(
        with_media, metadata, top_n=2, min_condition_profiles=1
    )
    pd.testing.assert_frame_equal(neighbors_base, neighbors_media)
    pd.testing.assert_frame_equal(centroids_base, centroids_media)
    pd.testing.assert_frame_equal(cross_base, cross_media)
    pd.testing.assert_frame_equal(
        followup.summarize_guild_connections(base, metadata),
        followup.summarize_guild_connections(with_media, metadata),
    )


def test_six_shared_sankeys_and_pd_callers_exclude_mediadive(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pd_cellulose = _load_script_module("create_pd_cellulose_analysis")
    pd_taxon = _load_script_module("create_pd_six_trait_taxon_heatmap")
    sankeys = _load_script_module("create_six_shared_trait_sankeys")

    shared_ids = [t.trait_id for t in sankeys.SHARED_TRAITS]
    rows = {
        "study_id": ["s1", "s1", "s2"],
        "experiment_id": ["e1", "e1", "e2"],
        "inc_dec": ["inc", "dec", "inc"],
        shared_ids[0]: [10.0, 20.0, 30.0],
        shared_ids[1]: [5.0, 15.0, 25.0],
        shared_ids[2]: [5.0, 5.0, 5.0],
        shared_ids[3]: [10.0, 10.0, 10.0],
        shared_ids[4]: [5.0, 10.0, 15.0],
        shared_ids[5]: [5.0, 0.0, 5.0],
        "OTHER_TRAIT": [60.0, 40.0, 10.0],
    }
    base_tsv = tmp_path / "base_profiles.tsv"
    media_tsv = tmp_path / "media_profiles.tsv"
    pd.DataFrame(rows).to_csv(base_tsv, sep="\t", index=False)
    pd.DataFrame({**rows, "mediadive.medium:1": [900.0, 10.0, 400.0]}).to_csv(
        media_tsv, sep="\t", index=False
    )

    metadata = pd.DataFrame(
        {"total_trait_signal": [100.0, 100.0, 100.0]},
        index=["s1|e1|inc", "s1|e1|dec", "s2|e2|inc"],
    )

    monkeypatch.setattr(sankeys, "PROFILE_MATRIX_PATH", base_tsv)
    raw_base = sankeys.load_raw_trait_counts(sankeys.SHARED_TRAITS, metadata)
    z_base = sankeys.build_zscore_matrix(raw_base, metadata)

    monkeypatch.setattr(sankeys, "PROFILE_MATRIX_PATH", media_tsv)
    raw_media = sankeys.load_raw_trait_counts(sankeys.SHARED_TRAITS, metadata)
    z_media = sankeys.build_zscore_matrix(raw_media, metadata)
    assert "OTHER_TRAIT" in raw_media.columns
    assert "mediadive.medium:1" not in raw_media.columns
    pd.testing.assert_frame_equal(z_base, z_media)

    monkeypatch.setattr(pd_cellulose, "PROFILE_MATRIX_PATH", base_tsv)
    pd_z_base = pd_cellulose.build_zscore_matrix(metadata, shared_ids[:2])
    monkeypatch.setattr(pd_cellulose, "PROFILE_MATRIX_PATH", media_tsv)
    pd_z_media = pd_cellulose.build_zscore_matrix(metadata, shared_ids[:2])
    pd.testing.assert_frame_equal(pd_z_base, pd_z_media)

    counts_base, cols_base = pd_taxon.read_profile_counts(base_tsv)
    counts_media, cols_media = pd_taxon.read_profile_counts(media_tsv)
    assert cols_base == cols_media
    pd.testing.assert_frame_equal(counts_base, counts_media)

    ref_base = pd_taxon.fixed_reference(counts_base, shared_ids[:2], metadata.index)
    with_injected_media = counts_base.assign(**{"mediadive.medium:1": [900.0, 10.0, 400.0]})
    ref_media = pd_taxon.fixed_reference(with_injected_media, shared_ids[:2], metadata.index)
    pd.testing.assert_frame_equal(ref_base[0], ref_media[0])
    pd.testing.assert_series_equal(ref_base[1], ref_media[1])


def test_jed_composition_excludes_mediadive_from_denominator() -> None:
    from mtp.composition_calibration.adapters.jed_composition import (
        JedCompositionAdapter,
        compute_profile_zscores,
    )

    base = pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s2", "s2"],
            "experiment_id": ["e1", "e1", "e2", "e2"],
            "inc_dec": ["inc", "dec", "inc", "dec"],
            "trait_a": [10.0, 20.0, 15.0, 25.0],
            "trait_b": [30.0, 20.0, 25.0, 15.0],
        },
        index=["p1", "p2", "p3", "p4"],
    )
    with_media = base.assign(**{"mediadive.medium:1": [900.0, 5.0, 500.0, 10.0]})

    z_base = compute_profile_zscores(base, ["trait_a", "trait_b"])
    z_media = compute_profile_zscores(with_media, ["trait_a", "trait_b", "mediadive.medium:1"])
    pd.testing.assert_frame_equal(z_base, z_media)

    adapter = JedCompositionAdapter()
    res_base = adapter.run_inference(base, pd.DataFrame(), {"traits": ["trait_a", "trait_b"]})
    res_media = adapter.run_inference(
        with_media,
        pd.DataFrame(),
        {"traits": ["trait_a", "trait_b", "mediadive.medium:1"]},
    )
    pd.testing.assert_frame_equal(res_base.results_df, res_media.results_df)


def test_pd_network_and_calibration_and_cache_loaders_exclude_mediadive(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from mtp.composition_calibration.composition import (
        load_composition_trait_columns,
        run_composition_calibration,
    )
    from mtp.density_adjustment import (
        get_non_media_trait_columns,
        normalize_traits_by_profile_composition,
    )

    pd_network = _load_script_module("create_pd_significant_condition_taxon_trait_network")
    followup = _load_script_module("analyze_enrichment_infographic_followup")
    aerobe = _load_script_module("investigate_facultative_aerobe_expansion")
    anaerobe = _load_script_module("investigate_gut_anaerobe_depletion")

    # 1. PD significant condition taxon-trait network
    shared_ids = ["METPO:1000608", "METPO:1000668"]
    rows = {
        "study_id": ["s1", "s1", "s2"],
        "experiment_id": ["e1", "e1", "e2"],
        "inc_dec": ["inc", "dec", "inc"],
        shared_ids[0]: [10.0, 20.0, 30.0],
        shared_ids[1]: [30.0, 20.0, 10.0],
        "OTHER_TRAIT": [60.0, 60.0, 60.0],
    }
    base_tsv = tmp_path / "pd_net_base.tsv"
    media_tsv = tmp_path / "pd_net_media.tsv"
    pd.DataFrame(rows).to_csv(base_tsv, sep="\t", index=False)
    pd.DataFrame({**rows, "mediadive.medium:1": [900.0, 5.0, 500.0]}).to_csv(
        media_tsv, sep="\t", index=False
    )
    metadata = pd.DataFrame(
        {"total_trait_signal": [100.0, 100.0, 100.0]},
        index=["s1|e1|inc", "s1|e1|dec", "s2|e2|inc"],
    )
    traits_df = pd.DataFrame({"trait_id": shared_ids})
    monkeypatch.setattr(pd_network, "PROFILE_MATRIX_PATH", base_tsv)
    z_base = pd_network.build_zscore_matrix(metadata, traits_df)
    monkeypatch.setattr(pd_network, "PROFILE_MATRIX_PATH", media_tsv)
    z_media = pd_network.build_zscore_matrix(metadata, traits_df)
    pd.testing.assert_frame_equal(z_base, z_media)

    # 2. Enrichment followup recomputes trait count and total signal from retained counts
    target_cond = followup.TARGET_CONDITIONS[0]
    other_cond = followup.TARGET_CONDITIONS[1]
    guild_trait = next(iter(next(iter(followup.GUILD_TRAITS.values()))))
    stale_meta = pd.DataFrame(
        {
            "profile_id": ["p1", "p2", "p3"],
            "study_id": ["s1", "s2", "s3"],
            "experiment_id": ["e1", "e2", "e3"],
            "inc_dec": ["inc", "dec", "inc"],
            "direction_label": ["increased", "decreased", "increased"],
            "condition": [target_cond, target_cond, other_cond],
            "title": ["t1", "t2", "t3"],
            "experiment_label": ["l1", "l2", "l3"],
            "body_site": ["Feces", "Feces", "Feces"],
            "location": ["USA", "USA", "USA"],
            "profile_trait_count": [99, 99, 99],
            "total_trait_signal": [999.0, 999.0, 999.0],
        },
        index=["p1", "p2", "p3"],
    )
    vals = pd.DataFrame(
        {
            guild_trait: [10.0, 20.0, 15.0],
            "METPO:1000668": [30.0, 0.0, 25.0],
            "METPO:2000517|mediadive.medium:1": [500.0, 600.0, 700.0],
        },
        index=["p1", "p2", "p3"],
    )
    neighbors, _, _ = followup.compute_profile_neighbors(
        vals, stale_meta, top_n=2, min_condition_profiles=1
    )
    p1_rows = neighbors[neighbors["neighbor_profile_id"] == "p1"]
    assert (p1_rows["profile_trait_count"] == 2).all()
    assert (p1_rows["total_trait_signal"] == 40.0).all()

    # 3. Composition calibration & density adjustment filter mediadive.medium: columns
    cat_tsv = tmp_path / "cats.tsv"
    pd.DataFrame(
        {
            "column_id": ["trait_a", "trait_b", "METPO:2000517|mediadive.medium:1"],
            "object_category": ["METPO:1000608", "METPO:1000668", "METPO:9999999"],
        }
    ).to_csv(cat_tsv, sep="\t", index=False)
    assert load_composition_trait_columns(cat_tsv) == ["trait_a", "trait_b"]
    assert get_non_media_trait_columns(cat_tsv) == ["trait_a", "trait_b"]

    calib_profiles = pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s2", "s2"],
            "experiment_id": ["e1", "e1", "e2", "e2"],
            "inc_dec": ["inc", "dec", "inc", "dec"],
            "condition": ["cond1", "cond1", "cond1", "cond1"],
            "trait_a": [10.0, 20.0, 15.0, 25.0],
            "trait_b": [30.0, 20.0, 25.0, 15.0],
            "mediadive.medium:1": [900.0, 5.0, 500.0, 10.0],
        }
    )
    calib_meta, calib_summary, _, _, _ = run_composition_calibration(
        calib_profiles,
        ["trait_a", "trait_b", "mediadive.medium:1"],
        replicates=4,
        seed=42,
    )
    assert calib_meta["n_traits"] == 2
    assert (calib_summary["n_tests"] == 2).all()

    norm_base = normalize_traits_by_profile_composition(
        calib_profiles.drop(columns=["mediadive.medium:1"]),
        ["trait_a", "trait_b"],
    )
    norm_media = normalize_traits_by_profile_composition(
        calib_profiles,
        ["trait_a", "trait_b", "mediadive.medium:1"],
    )
    assert norm_base.shape == norm_media.shape == (4, 2)
    assert norm_base == pytest.approx(norm_media)

    # 4. Cached heatmap loaders in aerobe & anaerobe scripts reject missing denominator_definition
    vis_dir = tmp_path / "vis"
    vis_dir.mkdir()
    pd.DataFrame({"display_a": [0.5, -0.5]}, index=["p1", "p2"]).to_csv(
        vis_dir / "study_trait_profile_clustered_heatmap_matrix.tsv", sep="\t"
    )
    pd.DataFrame(
        {
            "study_id": ["s1", "s2"],
            "experiment_id": ["e1", "e2"],
            "inc_dec": ["inc", "dec"],
            "total_trait_signal": [40.0, 40.0],
        },
        index=["p1", "p2"],
    ).to_csv(vis_dir / "study_trait_profile_study_metadata.tsv", sep="\t")
    pd.DataFrame(
        {
            "trait_id": ["trait_a"],
            "display_label": ["display_a"],
            "object_name": ["Trait A"],
        },
        index=["trait_a"],
    ).to_csv(vis_dir / "study_trait_profile_trait_metadata.tsv", sep="\t")

    monkeypatch.setattr(aerobe, "VIS_DIR", vis_dir)
    monkeypatch.setattr(anaerobe, "VIS_DIR", vis_dir)
    with pytest.raises(ValueError, match="missing the harmonized non-MediaDive denominator marker"):
        aerobe.load_study_heatmap()
    with pytest.raises(ValueError, match="missing the harmonized non-MediaDive denominator marker"):
        anaerobe.load_study_heatmap()
