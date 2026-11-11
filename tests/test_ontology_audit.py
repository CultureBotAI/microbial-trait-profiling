"""Tests for METPO ontology redundancy audit and branch collapsing (Issue #264 Task 2)."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from mtp.density_adjustment import (
    DEFAULT_COLUMNS_PATH,
    DEFAULT_PROFILES_PATH,
)
from mtp.ontology_audit import (
    RedundancyCluster,
    build_redundant_pairs_table,
    collapse_trait_matrix,
    get_canonical_redundancy_clusters,
    run_ontology_redundancy_audit,
)


def test_canonical_redundancy_clusters_structure() -> None:
    """Verify that canonical redundancy clusters are properly defined and consistent."""
    clusters = get_canonical_redundancy_clusters()
    assert len(clusters) == 10

    cluster_ids = [c.cluster_id for c in clusters]
    assert len(cluster_ids) == len(set(cluster_ids)), "Cluster IDs must be unique"

    for cluster in clusters:
        assert isinstance(cluster, RedundancyCluster)
        assert cluster.surviving_trait
        assert len(cluster.collapsed_traits) >= 1
        assert cluster.surviving_trait not in cluster.collapsed_traits
        assert cluster.tier in {
            "ancestor_descendant",
            "lexical_synonym",
            "synonym",
            "ancestor_descendant_and_synonym",
            "synonym_and_ancestor_descendant",
            "sibling_near_twin",
        }
        assert cluster.selection_rule
        assert cluster.selection_rule.startswith("max_support")

    by_id = {c.cluster_id: c for c in clusters}
    co_ox = by_id["co_oxidation_pathway"]
    assert co_ox.surviving_trait == "METPO:2000103|pathways:carbon_monoxide_oxidation"
    assert co_ox.collapsed_traits == ["METPO:2000103|pathways:carbonmonoxide_oxidation"]
    assert (
        "13,744 profile counts across 2,752 nonzero profiles"
        in by_id["bacillus_rod_shape"].rationale
    )


def test_collapse_trait_matrix_synthetic() -> None:
    """Verify matrix collapsing under both 'native' and 'max' count semantics."""
    df = pd.DataFrame(
        {
            "study_id": ["s1", "s2"],
            "experiment_id": ["e1", "e2"],
            "inc_dec": ["inc", "dec"],
            "t_parent": [10.0, 5.0],
            "t_child1": [8.0, 12.0],
            "t_unrelated": [1.0, 2.0],
        }
    )
    cols = ["t_parent", "t_child1", "t_unrelated"]
    cluster = RedundancyCluster(
        cluster_id="test_cluster",
        cluster_label="Test Cluster",
        tier="ancestor_descendant",
        surviving_trait="t_parent",
        surviving_label="parent",
        collapsed_traits=["t_child1"],
        selection_rule="max_support",
        rationale="test",
    )

    # Native semantics
    mat_nat, cols_nat = collapse_trait_matrix(df, cols, [cluster], count_semantics="native")
    assert cols_nat == ["t_parent", "t_unrelated"]
    assert np.array_equal(mat_nat[:, 0], np.array([10.0, 5.0]))
    assert np.array_equal(mat_nat[:, 1], np.array([1.0, 2.0]))

    # Max semantics: for s2, max(5.0, 12.0) = 12.0
    mat_max, cols_max = collapse_trait_matrix(df, cols, [cluster], count_semantics="max")
    assert cols_max == ["t_parent", "t_unrelated"]
    assert np.array_equal(mat_max[:, 0], np.array([10.0, 12.0]))
    assert np.array_equal(mat_max[:, 1], np.array([1.0, 2.0]))


def test_collapse_trait_matrix_invalid_semantics() -> None:
    """Verify error raised on invalid count semantics."""
    df = pd.DataFrame({"t1": [1.0], "t2": [2.0]})
    cluster = RedundancyCluster(
        cluster_id="c",
        cluster_label="c",
        tier="synonym",
        surviving_trait="t1",
        surviving_label="t1",
        collapsed_traits=["t2"],
        selection_rule="rule",
        rationale="rat",
    )
    invalid_semantics: Any = "invalid"
    with pytest.raises(ValueError, match="Unknown count_semantics"):
        collapse_trait_matrix(df, ["t1", "t2"], [cluster], count_semantics=invalid_semantics)


def test_build_redundant_pairs_table_synthetic(tmp_path: Path) -> None:
    """Verify pairwise table construction and metric calculation."""
    df = pd.DataFrame(
        {
            "t1": [1.0, 2.0, 3.0],
            "t2": [2.0, 4.0, 6.0],  # perfect correlation r=1.0
        }
    )
    cluster = RedundancyCluster(
        cluster_id="test",
        cluster_label="Test",
        tier="synonym",
        surviving_trait="t1",
        surviving_label="trait 1",
        collapsed_traits=["t2"],
        selection_rule="rule",
        rationale="rat",
    )

    # Minimal columns DataFrame
    tmp_cols = pd.DataFrame(
        {
            "column_id": ["t1", "t2"],
            "column_label": ["Trait One", "Trait Two"],
            "object_category": ["cat", "cat"],
        }
    )

    cols_path = tmp_path / "columns.tsv"
    tmp_cols.to_csv(cols_path, sep="\t", index=False)

    pairs_df = build_redundant_pairs_table(df, columns_path=cols_path, clusters=[cluster])
    assert len(pairs_df) == 1
    row = pairs_df.iloc[0]
    assert row["trait_1"] == "t1"
    assert row["trait_2"] == "t2"
    assert pytest.approx(row["correlation"]) == 1.0
    assert row["support_trait_1"] == 6
    assert row["support_trait_2"] == 12
    assert row["surviving_trait"] == "t1"


def _has_real_derived_data() -> bool:
    profiles_exist = (
        DEFAULT_PROFILES_PATH.exists()
        or Path(f"{DEFAULT_PROFILES_PATH}.gz").exists()
        or DEFAULT_PROFILES_PATH.with_name(f"{DEFAULT_PROFILES_PATH.name}.gz").exists()
    )
    columns_exist = (
        DEFAULT_COLUMNS_PATH.exists()
        or Path(f"{DEFAULT_COLUMNS_PATH}.gz").exists()
        or DEFAULT_COLUMNS_PATH.with_name(f"{DEFAULT_COLUMNS_PATH.name}.gz").exists()
    )
    return profiles_exist and columns_exist


@pytest.mark.skipif(
    not _has_real_derived_data(),
    reason="Requires derived study trait profiles and columns tables (.tsv or .tsv.gz)",
)
def test_ontology_redundancy_audit_integration(tmp_path: Path) -> None:
    """Verify the full ontology redundancy audit pipeline on real derived data."""
    prof_path = (
        DEFAULT_PROFILES_PATH
        if DEFAULT_PROFILES_PATH.exists()
        else DEFAULT_PROFILES_PATH.with_name(f"{DEFAULT_PROFILES_PATH.name}.gz")
    )
    cols_path = (
        DEFAULT_COLUMNS_PATH
        if DEFAULT_COLUMNS_PATH.exists()
        else DEFAULT_COLUMNS_PATH.with_name(f"{DEFAULT_COLUMNS_PATH.name}.gz")
    )
    results = run_ontology_redundancy_audit(
        profiles_path=prof_path,
        columns_path=cols_path,
        output_dir=tmp_path,
    )

    # Check generated output files
    assert results["pairs_tsv"].exists()
    assert results["summary_tsv"].exists()
    assert results["comparison_tsv"].exists()

    # Check standalone counts bounding trait space
    assert results["n_clusters"] == 10
    assert results["n_redundant_pairs"] == 24
    assert results["n_collapsed_features"] == 16
    assert results["bounded_trait_space"] == 301

    # Verify that raw traits were 317
    assert results["raw_metrics"].n_traits == 317
    assert results["native_metrics"].n_traits == 301
    assert results["max_metrics"].n_traits == 301
    assert results["adj_collapsed_metrics"].n_traits == 301

    # Key scientific assertions:
    # 1. PC1 variance in collapsed matrix remains > 60% (curation density, not redundancy!)
    assert results["max_metrics"].pc1_variance_pct > 60.0

    # 2. PC1 regression residualization drops PC1 variance substantially (from ~66.9% to ~38.7%)
    assert results["adj_collapsed_metrics"].pc1_variance_pct < 40.0

    # 3. Post-PC1 residualization median pairwise r is centered near 0
    assert pytest.approx(results["adj_collapsed_metrics"].median_pairwise_r, abs=0.05) == 0.0

    # 4. In comparison table, verify all rows are populated
    comp_df = pd.read_csv(results["comparison_tsv"], sep="\t")
    assert len(comp_df) == 4
    assert set(comp_df["configuration"]) == {
        "raw_uncollapsed_317",
        "collapsed_native_survivor_301",
        "collapsed_deduplicated_max_301",
        "collapsed_max_plus_pc1_adjusted_301",
    }

    # 5. near_twin_r is finite in raw and NaN in collapsed configurations (#475)
    raw_row = comp_df[comp_df["configuration"] == "raw_uncollapsed_317"].iloc[0]
    assert not np.isnan(raw_row["near_twin_r"])
    for cfg in [
        "collapsed_native_survivor_301",
        "collapsed_deduplicated_max_301",
        "collapsed_max_plus_pc1_adjusted_301",
    ]:
        collapsed_row = comp_df[comp_df["configuration"] == cfg].iloc[0]
        assert np.isnan(collapsed_row["near_twin_r"])

    # 6. Surviving trait has at least as much support as each collapsed member (#473)
    pairs_df = results["df_pairs"]
    survivor_pairs = pairs_df[pairs_df["trait_1"] == pairs_df["surviving_trait"]]
    assert not survivor_pairs.empty
    assert (survivor_pairs["support_trait_1"] >= survivor_pairs["support_trait_2"]).all()


def test_build_redundant_pairs_table_gz_columns(tmp_path: Path) -> None:
    """Verify build_redundant_pairs_table handles .gz columns file and fallback (#475)."""
    df = pd.DataFrame(
        {
            "study_id": ["s1", "s2", "s3"],
            "experiment_id": ["e1", "e2", "e3"],
            "t1": [1.0, 2.0, 3.0],
            "t2": [2.0, 4.0, 6.0],
        }
    )
    cluster = RedundancyCluster(
        cluster_id="c1",
        cluster_label="Test",
        tier="synonym",
        surviving_trait="t1",
        surviving_label="Trait One",
        collapsed_traits=["t2"],
        selection_rule="rule",
        rationale="rat",
    )
    tmp_cols = pd.DataFrame(
        {
            "column_id": ["t1", "t2"],
            "column_label": ["Trait One", "Trait Two"],
            "object_category": ["cat", "cat"],
        }
    )
    gz_cols_path = tmp_path / "columns.tsv.gz"
    tmp_cols.to_csv(gz_cols_path, sep="\t", index=False, compression="gzip")

    # 1. Direct .gz path
    pairs_direct = build_redundant_pairs_table(df, columns_path=gz_cols_path, clusters=[cluster])
    assert len(pairs_direct) == 1
    assert pairs_direct.iloc[0]["trait_1"] == "t1"

    # 2. Fallback when uncompressed path given but only .gz exists
    uncompressed_path = tmp_path / "columns.tsv"
    pairs_fallback = build_redundant_pairs_table(
        df, columns_path=uncompressed_path, clusters=[cluster]
    )
    assert len(pairs_fallback) == 1
    assert pairs_fallback.iloc[0]["trait_1"] == "t1"

    # 3. Missing file raises FileNotFoundError
    non_existent = tmp_path / "does_not_exist.tsv"
    with pytest.raises(FileNotFoundError, match="Columns file not found"):
        build_redundant_pairs_table(df, columns_path=non_existent, clusters=[cluster])


def test_compute_correlation_metrics_near_twin_nan() -> None:
    """Verify near_twin_r evaluates to NaN when near-twin pair is missing (#475)."""
    from mtp.density_adjustment import (
        TWIN_A,
        TWIN_B,
        compute_correlation_metrics,
    )

    # 1. When both twins are present
    rng = np.random.default_rng(42)
    mat = rng.standard_normal((50, 3))
    metrics_present = compute_correlation_metrics(mat, [TWIN_A, TWIN_B, "other_trait"])
    assert not np.isnan(metrics_present.near_twin_r)

    # 2. When TWIN_B is absent (collapsed out)
    metrics_absent = compute_correlation_metrics(mat[:, [0, 2]], [TWIN_A, "other_trait"])
    assert np.isnan(metrics_absent.near_twin_r)


def test_write_consolidated_report_dynamic_numbers(tmp_path: Path) -> None:
    """Verify write_consolidated_report formats numbers dynamically in Section 6 (#475)."""
    import importlib.util

    from mtp.density_adjustment import (
        ConditionalTwinSummary,
        CorrelationMetrics,
        Table1TraitEffect,
    )

    script_path = (
        Path(__file__).resolve().parents[1] / "scripts" / "audit_trait_density_correlation.py"
    )
    spec = importlib.util.spec_from_file_location("audit_trait_density_correlation", script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    dummy_metrics = CorrelationMetrics(
        n_samples=100,
        n_traits=10,
        constant_traits_excluded=0,
        median_pairwise_r=0.1,
        mean_pairwise_r=0.1,
        pct_pairs_gt_05=5.0,
        pct_pairs_gt_09=1.0,
        pc1_variance_pct=55.55,
        pc1_5_variance_pct=75.0,
        m_eff_li_ji=8.0,
        m_eff_galwey=6.0,
        m_eff_nyholt=9.0,
        test_pair_r=0.3,
        near_twin_r=float("nan"),
    )

    density_res = {
        "raw_pearson": dummy_metrics,
        "raw_spearman": dummy_metrics,
        "adj_pc1_pearson": dummy_metrics,
        "table1_effects": [
            Table1TraitEffect(
                trait_id="METPO:1000608",
                trait_label="facultatively aerobic",
                raw_mean_diff=0.05,
                raw_wilcoxon_p=0.000345,
                adj_mean_diff=0.01,
                adj_wilcoxon_p=0.456,
                relative_change_pct=-80.0,
                status="collapsed (attenuated to nonsignificant)",
                magnitude_change_pct=-80.0,
                sign_flipped=False,
                raw_wilcoxon_bh_q=0.00069,
                adj_wilcoxon_bh_q=0.456,
                bh_family_size=2,
            ),
            Table1TraitEffect(
                trait_id="METPO:2000103|GO:0030245",
                trait_label="cellulose catabolic process",
                raw_mean_diff=-0.08,
                raw_wilcoxon_p=1e-8,
                adj_mean_diff=-0.09,
                adj_wilcoxon_p=3.4e-12,
                relative_change_pct=-12.5,
                status="enhanced (significant, larger magnitude)",
                magnitude_change_pct=12.5,
                sign_flipped=False,
                raw_wilcoxon_bh_q=1e-8,
                adj_wilcoxon_bh_q=6.8e-12,
                bh_family_size=2,
            ),
        ],
        "twin_summary": ConditionalTwinSummary(
            trait_a="METPO:1000668",
            trait_b="METPO:1000683",
            trait_a_label="coccus shaped",
            trait_b_label="sphere shaped",
            profile_correlation=0.92,
            profile_partial_r_depth=0.85,
            vif_profile=6.12,
            profile_slope_a_on_b=0.95,
            profile_intercept_a_on_b=0.01,
            n_paired_experiments=1234,
            paired_correlation=0.88,
            vif_paired=4.5,
            paired_slope_a_on_b=0.91,
            paired_slope_pvalue=1e-15,
        ),
    }

    ontology_res = {
        "n_redundant_pairs": 18,
        "n_clusters": 8,
        "n_collapsed_features": 12,
        "bounded_trait_space": 250,
        "native_metrics": dummy_metrics,
        "max_metrics": dummy_metrics,
        "adj_collapsed_metrics": dummy_metrics,
    }

    report_file = module.write_consolidated_report(
        output_dir=tmp_path,
        density_res=density_res,
        ontology_res=ontology_res,
    )
    content = report_file.read_text(encoding="utf-8")

    # Assert dynamic formatting
    assert "N=1,234" in content
    assert "VIF = 6.12" in content
    # Section 6 reads the BH q-values and labels of the reported effects (#1321).
    assert (
        "After PC1 regression the facultatively aerobic paired contrast is labelled "
        "'collapsed (attenuated to nonsignificant)'"
    ) in content
    assert "(raw BH q=6.90e-04, adjusted BH q=4.56e-01)" in content
    assert (
        "Significant before and after PC1 regression in the same direction "
        "(BH q <= 0.05; survived or enhanced): cellulose catabolic process."
    ) in content
    assert "(nominal p=3.45e-04, q=6.90e-04)" in content
    assert "(nominal p=3.40e-12, q=6.80e-12)" in content
    assert "over all 2 reported traits" in content
    assert "18 redundant pairs (12 features collapsed into 8 clusters)" in content
    assert "bounding distinct trait space to 250 features" in content

    # The report never states a family it does not report: no rows, or a family size
    # other than the number of rows, is an error rather than "over all 0 traits" (#1324).
    effects = density_res["table1_effects"]
    assert isinstance(effects, list)
    for bad_effects in ([], [replace(effect, bh_family_size=6) for effect in effects]):
        with pytest.raises(ValueError, match="family"):
            module.write_consolidated_report(
                output_dir=tmp_path / "bad",
                density_res={**density_res, "table1_effects": bad_effects},
                ontology_res=ontology_res,
            )


def test_ontology_redundancy_audit_counts_observed_clusters_only(tmp_path: Path) -> None:
    """#1402: n_clusters counts clusters with at least one redundant pair present in the inputs."""
    traits = ["METPO:1000668", "METPO:1000683", "METPO:1000608"]
    profiles = tmp_path / "profiles.tsv"
    columns = tmp_path / "columns.tsv"
    pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s2", "s2"],
            "experiment_id": ["e1", "e1", "e2", "e2"],
            "inc_dec": ["inc", "dec", "inc", "dec"],
            "METPO:1000668": [2.0, 4.0, 6.0, 8.0],
            "METPO:1000683": [1.0, 2.0, 3.0, 4.0],
            "METPO:1000608": [5.0, 1.0, 2.0, 3.0],
        }
    ).to_csv(profiles, sep="\t", index=False)
    pd.DataFrame(
        {
            "column_id": traits,
            "column_label": ["coccus shaped", "sphere shaped", "facultatively aerobic"],
            "object_name": ["coccus shaped", "sphere shaped", "facultatively aerobic"],
            "object_category": ["biolink:PhenotypicQuality"] * len(traits),
        }
    ).to_csv(columns, sep="\t", index=False)

    res = run_ontology_redundancy_audit(
        profiles_path=profiles,
        columns_path=columns,
        output_dir=tmp_path / "out",
    )
    assert res["n_redundant_pairs"] == 1
    assert res["n_clusters"] == 1
    assert res["n_collapsed_features"] == 1
    summary_df = pd.read_csv(res["summary_tsv"], sep="\t").set_index("metric")["value"].to_dict()
    assert summary_df["redundant_clusters_count"] == 1.0
    assert summary_df["redundant_pairs_count"] == 1.0

    # Zero observed pairs -> 0 observed clusters
    no_pair_traits = ["METPO:1000668", "METPO:1000608"]
    profiles_zero = tmp_path / "profiles_zero.tsv"
    columns_zero = tmp_path / "columns_zero.tsv"
    pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s2", "s2"],
            "experiment_id": ["e1", "e1", "e2", "e2"],
            "inc_dec": ["inc", "dec", "inc", "dec"],
            "METPO:1000668": [2.0, 4.0, 6.0, 8.0],
            "METPO:1000608": [5.0, 1.0, 2.0, 3.0],
        }
    ).to_csv(profiles_zero, sep="\t", index=False)
    pd.DataFrame(
        {
            "column_id": no_pair_traits,
            "column_label": ["coccus shaped", "facultatively aerobic"],
            "object_name": ["coccus shaped", "facultatively aerobic"],
            "object_category": ["biolink:PhenotypicQuality"] * len(no_pair_traits),
        }
    ).to_csv(columns_zero, sep="\t", index=False)

    res_zero = run_ontology_redundancy_audit(
        profiles_path=profiles_zero,
        columns_path=columns_zero,
        output_dir=tmp_path / "out_zero",
    )
    assert res_zero["n_redundant_pairs"] == 0
    assert res_zero["n_clusters"] == 0

    # Multiple redundant pairs in one cluster (#1825 K02): motility has survivor METPO:1000702
    # plus three collapsed traits (3 pairs in 1 cluster when 3 members are present).
    multi_pair_traits = ["METPO:1000702", "METPO:1000705", "METPO:1000704", "METPO:1000608"]
    profiles_multi = tmp_path / "profiles_multi.tsv"
    columns_multi = tmp_path / "columns_multi.tsv"
    pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s2", "s2"],
            "experiment_id": ["e1", "e1", "e2", "e2"],
            "inc_dec": ["inc", "dec", "inc", "dec"],
            "METPO:1000702": [3.0, 6.0, 9.0, 12.0],
            "METPO:1000705": [2.0, 4.0, 6.0, 8.0],
            "METPO:1000704": [1.0, 2.0, 3.0, 4.0],
            "METPO:1000608": [5.0, 1.0, 2.0, 3.0],
        }
    ).to_csv(profiles_multi, sep="\t", index=False)
    pd.DataFrame(
        {
            "column_id": multi_pair_traits,
            "column_label": [
                "motile",
                "axially filamented",
                "flagellated",
                "facultatively aerobic",
            ],
            "object_name": ["motile", "axially filamented", "flagellated", "facultatively aerobic"],
            "object_category": ["biolink:PhenotypicQuality"] * len(multi_pair_traits),
        }
    ).to_csv(columns_multi, sep="\t", index=False)

    res_multi = run_ontology_redundancy_audit(
        profiles_path=profiles_multi,
        columns_path=columns_multi,
        output_dir=tmp_path / "out_multi",
    )
    assert res_multi["n_redundant_pairs"] == 3
    assert res_multi["n_clusters"] == 1
    assert res_multi["n_collapsed_features"] == 2
    assert res_multi["n_orphan_features"] == 0

    # Paired cluster (coccus + sphere) plus orphan-removed member (flagellated without motile,
    # #1857)
    mixed_traits = ["METPO:1000668", "METPO:1000683", "METPO:1000704", "METPO:1000608"]
    profiles_mixed = tmp_path / "profiles_mixed.tsv"
    columns_mixed = tmp_path / "columns_mixed.tsv"
    pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s2", "s2"],
            "experiment_id": ["e1", "e1", "e2", "e2"],
            "inc_dec": ["inc", "dec", "inc", "dec"],
            "METPO:1000668": [2.0, 4.0, 6.0, 8.0],
            "METPO:1000683": [1.0, 2.0, 3.0, 4.0],
            "METPO:1000704": [1.0, 2.0, 3.0, 4.0],
            "METPO:1000608": [5.0, 1.0, 2.0, 3.0],
        }
    ).to_csv(profiles_mixed, sep="\t", index=False)
    pd.DataFrame(
        {
            "column_id": mixed_traits,
            "column_label": [
                "coccus shaped",
                "sphere shaped",
                "flagellated",
                "facultatively aerobic",
            ],
            "object_name": [
                "coccus shaped",
                "sphere shaped",
                "flagellated",
                "facultatively aerobic",
            ],
            "object_category": ["biolink:PhenotypicQuality"] * len(mixed_traits),
        }
    ).to_csv(columns_mixed, sep="\t", index=False)

    res_mixed = run_ontology_redundancy_audit(
        profiles_path=profiles_mixed,
        columns_path=columns_mixed,
        output_dir=tmp_path / "out_mixed",
    )
    assert res_mixed["n_redundant_pairs"] == 1
    assert res_mixed["n_clusters"] == 1
    assert res_mixed["n_collapsed_features"] == 2
    assert res_mixed["n_orphan_features"] == 1

    # Multi-orphan case (#1857, PR #1901 review): motility survivor (METPO:1000702) is absent
    # while two of its collapsed members (METPO:1000705, METPO:1000704) are present alongside
    # coccus + sphere. The member-member motility pair is in df_pairs, but n_clusters must count
    # only clusters whose survivor is present (1 cluster: morphology_sphere, not 2).
    multi_orphan_traits = [
        "METPO:1000668",
        "METPO:1000683",
        "METPO:1000705",
        "METPO:1000704",
        "METPO:1000608",
    ]
    profiles_multi_orphan = tmp_path / "profiles_multi_orphan.tsv"
    columns_multi_orphan = tmp_path / "columns_multi_orphan.tsv"
    pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s2", "s2"],
            "experiment_id": ["e1", "e1", "e2", "e2"],
            "inc_dec": ["inc", "dec", "inc", "dec"],
            "METPO:1000668": [2.0, 4.0, 6.0, 8.0],
            "METPO:1000683": [1.0, 2.0, 3.0, 4.0],
            "METPO:1000705": [2.0, 4.0, 6.0, 8.0],
            "METPO:1000704": [1.0, 2.0, 3.0, 4.0],
            "METPO:1000608": [5.0, 1.0, 2.0, 3.0],
        }
    ).to_csv(profiles_multi_orphan, sep="\t", index=False)
    pd.DataFrame(
        {
            "column_id": multi_orphan_traits,
            "column_label": [
                "coccus shaped",
                "sphere shaped",
                "axially filamented",
                "flagellated",
                "facultatively aerobic",
            ],
            "object_name": [
                "coccus shaped",
                "sphere shaped",
                "axially filamented",
                "flagellated",
                "facultatively aerobic",
            ],
            "object_category": ["biolink:PhenotypicQuality"] * len(multi_orphan_traits),
        }
    ).to_csv(columns_multi_orphan, sep="\t", index=False)

    res_multi_orphan = run_ontology_redundancy_audit(
        profiles_path=profiles_multi_orphan,
        columns_path=columns_multi_orphan,
        output_dir=tmp_path / "out_multi_orphan",
    )
    assert res_multi_orphan["n_redundant_pairs"] == 2
    assert res_multi_orphan["n_clusters"] == 1
    assert res_multi_orphan["n_collapsed_features"] == 3
    assert res_multi_orphan["n_orphan_features"] == 2


def test_build_redundant_pairs_table_applies_display_trait_label(tmp_path: Path) -> None:
    """Verify build_redundant_pairs_table applies display_trait_label (#1250)."""
    df = pd.DataFrame(
        {
            "METPO:1000608": [1.0, 2.0, 3.0],
            "METPO:1000607": [2.0, 4.0, 6.0],
        }
    )
    cluster_survivor_608 = RedundancyCluster(
        cluster_id="facultative_test_1",
        cluster_label="Facultative Test 1",
        tier="synonym",
        surviving_trait="METPO:1000608",
        surviving_label="facultatively aerobic",
        collapsed_traits=["METPO:1000607"],
        selection_rule="max_support_canonical",
        rationale="test",
    )
    cluster_collapsed_608 = RedundancyCluster(
        cluster_id="facultative_test_2",
        cluster_label="Facultative Test 2",
        tier="synonym",
        surviving_trait="METPO:1000607",
        surviving_label="obligately anaerobic",
        collapsed_traits=["METPO:1000608"],
        selection_rule="max_support_canonical",
        rationale="test",
    )

    union_cols_path = tmp_path / "union_columns.tsv"
    pd.DataFrame(
        {
            "column_id": ["METPO:1000608", "METPO:1000607"],
            "column_label": ["facultatively aerobic", "obligately anaerobic"],
            "object_category": ["biolink:PhenotypicQuality", "biolink:PhenotypicQuality"],
            "trait_key_policy": ["facultative_union_v1", "facultative_union_v1"],
            "trait_member_ids_json": [
                '["METPO:1000605","METPO:1000608","METPO:1000612"]',
                '["METPO:1000607"]',
            ],
        }
    ).to_csv(union_cols_path, sep="\t", index=False)

    raw_cols_path = tmp_path / "raw_columns.tsv"
    pd.DataFrame(
        {
            "column_id": ["METPO:1000608", "METPO:1000607"],
            "column_label": ["facultatively aerobic", "obligately anaerobic"],
            "object_category": ["biolink:PhenotypicQuality", "biolink:PhenotypicQuality"],
            "trait_key_policy": ["raw_object_v1", "raw_object_v1"],
            "trait_member_ids_json": [
                '["METPO:1000608"]',
                '["METPO:1000607"]',
            ],
        }
    ).to_csv(raw_cols_path, sep="\t", index=False)

    union_pairs = build_redundant_pairs_table(
        df,
        columns_path=union_cols_path,
        clusters=[cluster_survivor_608, cluster_collapsed_608],
    )
    assert list(union_pairs["label_1"]) == [
        "facultative (aerobe or anaerobe)",
        "obligately anaerobic",
    ]
    assert list(union_pairs["label_2"]) == [
        "obligately anaerobic",
        "facultative (aerobe or anaerobe)",
    ]

    raw_pairs = build_redundant_pairs_table(
        df,
        columns_path=raw_cols_path,
        clusters=[cluster_survivor_608, cluster_collapsed_608],
    )
    assert list(raw_pairs["label_1"]) == [
        "facultatively aerobic",
        "obligately anaerobic",
    ]
    assert list(raw_pairs["label_2"]) == [
        "obligately anaerobic",
        "facultatively aerobic",
    ]
