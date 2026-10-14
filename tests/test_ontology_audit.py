"""Tests for METPO ontology redundancy audit and branch collapsing (Issue #264 Task 2)."""

from __future__ import annotations

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


@pytest.mark.skipif(
    not DEFAULT_PROFILES_PATH.exists() or not DEFAULT_COLUMNS_PATH.exists(),
    reason="Requires derived study trait profiles and columns tables",
)
def test_ontology_redundancy_audit_integration(tmp_path: Path) -> None:
    """Verify the full ontology redundancy audit pipeline on real derived data."""
    results = run_ontology_redundancy_audit(
        profiles_path=DEFAULT_PROFILES_PATH,
        columns_path=DEFAULT_COLUMNS_PATH,
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

    # 2. PC1 regression residualization successfully drops PC1 variance down to ~24%
    assert results["adj_collapsed_metrics"].pc1_variance_pct < 30.0

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
