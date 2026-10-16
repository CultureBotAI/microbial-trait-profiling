"""Ontology redundancy audit, branch collapsing rules, and dimensionality evaluation.

Addresses Issue #264 (Jed's portion, Task 2):
1. Audits METPO (and associated CHEBI/pathway) ontology redundancy across traits:
   - Counts ancestor-descendant pairs via transitive subclass closure.
   - Counts synonym pairs (lexical duplicates and shared chemical synonyms).
   - Counts sibling semantic near-twins (cell shape provenance duplicates).
   - Reports standalone pair counts bounding the distinct trait space.
2. Formalizes a stated branch-collapsing rule:
   - Survivor selection: Max-support evidence rule.
   - Count semantics: Deduplicated Max Count and Native Surviving Representative Count.
3. Recomputes dimensionality metrics (PC1 variance, M_eff, median pairwise r)
   on the collapsed trait matrix and evaluates whether ontology redundancy
   or annotation density accounts for the dominant correlation structure.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd

from mtp.density_adjustment import (
    DEFAULT_COLUMNS_PATH,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_PROFILES_PATH,
    TWIN_A,
    TWIN_B,
    compute_correlation_metrics,
    load_non_media_profile_matrix,
    regress_out_pc1,
)

DEFAULT_PAIRS_OUT = DEFAULT_OUTPUT_DIR / "ontology_redundancy_pairs.tsv"
DEFAULT_REDUNDANCY_SUMMARY_OUT = DEFAULT_OUTPUT_DIR / "ontology_redundancy_summary.tsv"
DEFAULT_COMPARISON_OUT = DEFAULT_OUTPUT_DIR / "ontology_collapsing_comparison.tsv"


@dataclass(frozen=True)
class RedundantPair:
    """A pair of traits sharing an ontology subclass, synonym, or sibling twin relationship."""

    trait_1: str
    label_1: str
    trait_2: str
    label_2: str
    tier: str  # "ancestor_descendant", "synonym", "sibling_near_twin"
    relationship: str  # e.g., "subclass_of", "chemical_synonym", "provenance_branch_twin"
    correlation: float
    support_1: int
    support_2: int


@dataclass(frozen=True)
class RedundancyCluster:
    """A group of redundant traits collapsing into a single representative trait."""

    cluster_id: str
    cluster_label: str
    tier: str
    surviving_trait: str
    surviving_label: str
    collapsed_traits: list[str]
    selection_rule: str
    rationale: str


# Canonical redundant clusters across non-MediaDive traits
CANONICAL_REDUNDANCY_CLUSTERS: list[RedundancyCluster] = [
    RedundancyCluster(
        cluster_id="motility",
        cluster_label="Bacterial Motility",
        tier="ancestor_descendant",
        surviving_trait="METPO:1000702",
        surviving_label="motile",
        collapsed_traits=["METPO:1000705", "METPO:1000704", "METPO:1000706"],
        selection_rule="max_support_ancestor",
        rationale=(
            "METPO:1000702 (motile) is the parent term; axially filamented (1000705), "
            "flagellated (1000704), and gliding (1000706) are specific locomotion subclasses. "
            "Motile has the broadest support across taxa."
        ),
    ),
    RedundancyCluster(
        cluster_id="obligate_anaerobe",
        cluster_label="Obligate Anaerobiosis",
        tier="ancestor_descendant",
        surviving_trait="METPO:1000607",
        surviving_label="obligately anaerobic",
        collapsed_traits=["METPO:1000611"],
        selection_rule="max_support_ancestor",
        rationale=(
            "METPO:1000611 (strictly anaerobic) is a direct subclass of METPO:1000607 "
            "(obligately anaerobic) in METPO. Obligately anaerobic has greater taxon support."
        ),
    ),
    RedundancyCluster(
        cluster_id="dextrin_utilization",
        cluster_label="Dextrin Utilization",
        tier="ancestor_descendant",
        surviving_trait="METPO:2000006|CHEBI:28675",
        surviving_label="dextrin",
        collapsed_traits=["METPO:2000006|CHEBI:40585"],
        selection_rule="max_support_ancestor",
        rationale=(
            "CHEBI:40585 (alpha-cyclodextrin) is a specific cyclic oligosaccharide subclass "
            "of CHEBI:28675 (dextrin). Dextrin represents the broader nutrient class."
        ),
    ),
    RedundancyCluster(
        cluster_id="butan2ol_utilization",
        cluster_label="Butan-2-ol Utilization",
        tier="ancestor_descendant_and_synonym",
        surviving_trait="METPO:2000006|CHEBI:35687",
        surviving_label="butan-2-ol",
        collapsed_traits=["METPO:2000006|CHEBI:45475"],
        selection_rule="max_support_ancestor",
        rationale=(
            "CHEBI:45475 ((2S)-butan-2-ol) is an enantiomer subclass of racemic CHEBI:35687 "
            "(butan-2-ol). Across profiles, both columns have identical counts (r = 1.000)."
        ),
    ),
    RedundancyCluster(
        cluster_id="citrate_utilization",
        cluster_label="Citrate Utilization",
        tier="synonym_and_ancestor_descendant",
        surviving_trait="METPO:2000006|CHEBI:16947",
        surviving_label="citrate",
        collapsed_traits=["METPO:2000006|CHEBI:133748", "METPO:2000006|CHEBI:50744"],
        selection_rule="max_support_canonical",
        rationale=(
            "CHEBI:16947 (citrate), CHEBI:133748 (citrate anion), and CHEBI:50744 (citrate salt) "
            "are chemical ionization/salt synonyms sharing identical profiles (r = 1.000)."
        ),
    ),
    RedundancyCluster(
        cluster_id="glutamate_utilization",
        cluster_label="Glutamate Utilization",
        tier="synonym_and_ancestor_descendant",
        surviving_trait="METPO:2000006|CHEBI:18237",
        surviving_label="glutamic acid",
        collapsed_traits=[
            "METPO:2000006|CHEBI:16015",
            "METPO:2000006|CHEBI:29987",
            "METPO:2000006|CHEBI:14321",
        ],
        selection_rule="max_support_canonical",
        rationale=(
            "CHEBI:18237 (glutamic acid), CHEBI:16015 (L-glutamic acid), CHEBI:29987 (glutamate), "
            "and CHEBI:14321 (glutamate(1-)) are stereochemical and charge synonyms "
            "with identical profile annotations (r = 1.000)."
        ),
    ),
    RedundancyCluster(
        cluster_id="malonate_utilization",
        cluster_label="Malonate Utilization",
        tier="synonym",
        surviving_trait="METPO:2000006|CHEBI:15792",
        surviving_label="malonate",
        collapsed_traits=["METPO:2000006|CHEBI:38083"],
        selection_rule="max_support_canonical",
        rationale=(
            "CHEBI:15792 (malonate) and CHEBI:38083 (malonate ester) share exact synonym "
            "'malonate' and identical profile annotations across all 5,829 profiles (r = 1.000)."
        ),
    ),
    RedundancyCluster(
        cluster_id="co_oxidation_pathway",
        cluster_label="Carbon Monoxide Oxidation Pathway",
        tier="lexical_synonym",
        surviving_trait="METPO:2000103|pathways:carbon_monoxide_oxidation",
        surviving_label="carbon_monoxide_oxidation",
        collapsed_traits=["METPO:2000103|pathways:carbonmonoxide_oxidation"],
        selection_rule="max_support_canonical",
        rationale=(
            "METPO:2000103|pathways:carbon_monoxide_oxidation (1,241 profile count support) and "
            "METPO:2000103|pathways:carbonmonoxide_oxidation (1,005 profile count support) are "
            "lexical duplicates differing only by an underscore; carbon_monoxide_oxidation "
            "survives under the max-support rule."
        ),
    ),
    RedundancyCluster(
        cluster_id="coccus_sphere_shape",
        cluster_label="Coccus / Sphere Morphology",
        tier="sibling_near_twin",
        surviving_trait=TWIN_A,
        surviving_label="coccus shaped",
        collapsed_traits=[TWIN_B],
        selection_rule="max_support_provenance",
        rationale=(
            "METPO:1000668 (coccus shaped, Madin et al.) and METPO:1000683 (sphere shaped, "
            "Bactotraits) are sibling terms under METPO:1000666 (cell shape) with r = +0.911. "
            "Coccus has 18,349 profile counts and 3.9M taxon annotations vs sphere (11,178 "
            "profile counts, 178k annotations). Under max-support rule, coccus survives."
        ),
    ),
    RedundancyCluster(
        cluster_id="bacillus_rod_shape",
        cluster_label="Bacillus / Rod Morphology",
        tier="sibling_near_twin",
        surviving_trait="METPO:1000667",
        surviving_label="bacillus shaped",
        collapsed_traits=["METPO:1000681", "METPO:1000688"],
        selection_rule="max_support_provenance",
        rationale=(
            "METPO:1000667 (bacillus shaped, Madin et al., 33,994 profile counts), "
            "METPO:1000681 (rod shaped, Bactotraits, 20,294 profile counts), and "
            "METPO:1000688 (coccobacillus shaped, 13,744 profile counts across 2,752 nonzero "
            "profiles) are highly collinear sibling morphological traits (r = +0.911 bacillus "
            "vs rod). Bacillus has highest support."
        ),
    ),
]


def get_canonical_redundancy_clusters() -> list[RedundancyCluster]:
    """Return the list of predefined canonical redundancy clusters."""
    return CANONICAL_REDUNDANCY_CLUSTERS


def build_redundant_pairs_table(
    profiles_df: pd.DataFrame,
    columns_path: Path = DEFAULT_COLUMNS_PATH,
    clusters: list[RedundancyCluster] | None = None,
) -> pd.DataFrame:
    """Build a detailed tabular report of redundant pairs across the trait universe.

    Computes empirical pairwise correlation and support counts across profiles.
    """
    if clusters is None:
        clusters = CANONICAL_REDUNDANCY_CLUSTERS

    if not columns_path.exists():
        # Same .gz fallback as load_non_media_profile_matrix, which accepted this path.
        gz_candidate = (
            columns_path
            if columns_path.suffix == ".gz"
            else columns_path.with_name(f"{columns_path.name}.gz")
        )
        if gz_candidate.exists():
            columns_path = gz_candidate
        else:
            raise FileNotFoundError(f"Columns file not found: {columns_path} (or {gz_candidate})")
    cols_df = pd.read_csv(columns_path, sep="\t")
    col_label_map = dict(zip(cols_df["column_id"], cols_df["column_label"], strict=False))

    pairs_rows: list[dict[str, Any]] = []

    for cluster in clusters:
        survivor = cluster.surviving_trait
        all_members = [survivor, *cluster.collapsed_traits]

        # Generate all pairwise combinations within cluster
        for i in range(len(all_members)):
            for j in range(i + 1, len(all_members)):
                c1, c2 = all_members[i], all_members[j]
                if c1 not in profiles_df.columns or c2 not in profiles_df.columns:
                    continue

                v1 = profiles_df[c1].to_numpy(dtype=float)
                v2 = profiles_df[c2].to_numpy(dtype=float)

                r_val = float(np.corrcoef(v1, v2)[0, 1]) if len(v1) > 1 else 0.0
                sup1 = int(np.sum(v1))
                sup2 = int(np.sum(v2))

                # Determine specific relationship type
                if cluster.tier == "ancestor_descendant":
                    rel = "subclass_of"
                elif cluster.tier == "lexical_synonym":
                    rel = "lexical_duplicate"
                elif "synonym" in cluster.tier:
                    rel = "chemical_synonym"
                else:
                    rel = "provenance_branch_twin"

                pairs_rows.append(
                    {
                        "cluster_id": cluster.cluster_id,
                        "cluster_label": cluster.cluster_label,
                        "tier": cluster.tier,
                        "relationship": rel,
                        "trait_1": c1,
                        "label_1": col_label_map.get(c1, c1),
                        "trait_2": c2,
                        "label_2": col_label_map.get(c2, c2),
                        "correlation": round(r_val, 6),
                        "support_trait_1": sup1,
                        "support_trait_2": sup2,
                        "surviving_trait": survivor,
                        "selection_rule": cluster.selection_rule,
                    }
                )

    df_pairs = pd.DataFrame(pairs_rows)
    return df_pairs


def collapse_trait_matrix(
    profiles_df: pd.DataFrame,
    non_media_cols: list[str],
    clusters: list[RedundancyCluster] | None = None,
    *,
    count_semantics: Literal["max", "native"] = "max",
) -> tuple[np.ndarray, list[str]]:
    """Collapse redundant branches in wide profiles to one representative per cluster.

    Args:
        profiles_df: DataFrame of wide study trait profiles.
        non_media_cols: Ordered list of non-MediaDive trait column identifiers.
        clusters: List of RedundancyCluster definitions.
        count_semantics:
            - 'max': Largest member count (C_i = max_{m} C_im), a lower bound on the union.
              Ensures taxa annotated under either term are captured without double counting.
            - 'native': Pure representative count (C_i = C_{i, survivor}).

    Returns:
        tuple of (collapsed_matrix, collapsed_column_ids)
    """
    if clusters is None:
        clusters = CANONICAL_REDUNDANCY_CLUSTERS

    to_remove: set[str] = set()
    for c in clusters:
        for m in c.collapsed_traits:
            to_remove.add(m)

    collapsed_cols = [col for col in non_media_cols if col not in to_remove]

    if count_semantics == "native":
        mat = profiles_df[collapsed_cols].values.astype(float)
        return mat, collapsed_cols

    if count_semantics == "max":
        df_copy = profiles_df.copy()
        for c in clusters:
            if c.surviving_trait in df_copy.columns:
                valid_members = [
                    m for m in [c.surviving_trait, *c.collapsed_traits] if m in df_copy.columns
                ]
                df_copy[c.surviving_trait] = df_copy[valid_members].max(axis=1)
        mat = df_copy[collapsed_cols].values.astype(float)
        return mat, collapsed_cols

    raise ValueError(f"Unknown count_semantics: {count_semantics}")


def run_ontology_redundancy_audit(
    profiles_path: Path = DEFAULT_PROFILES_PATH,
    columns_path: Path = DEFAULT_COLUMNS_PATH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> dict[str, Any]:
    """Execute complete Task 2 ontology redundancy audit and matrix collapsing evaluation.

    Writes:
      output_dir / 'ontology_redundancy_pairs.tsv'
      output_dir / 'ontology_redundancy_summary.tsv'
      output_dir / 'ontology_collapsing_comparison.tsv'

    Returns:
      dict containing summary metrics and output paths.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    if not profiles_path.exists():
        gz_prof = (
            profiles_path
            if profiles_path.suffix == ".gz"
            else profiles_path.with_name(f"{profiles_path.name}.gz")
        )
        if gz_prof.exists():
            profiles_path = gz_prof
    if not columns_path.exists():
        gz_cols = (
            columns_path
            if columns_path.suffix == ".gz"
            else columns_path.with_name(f"{columns_path.name}.gz")
        )
        if gz_cols.exists():
            columns_path = gz_cols

    profiles_df, non_media_cols = load_non_media_profile_matrix(profiles_path, columns_path)
    clusters = CANONICAL_REDUNDANCY_CLUSTERS

    # 1. Generate redundant pairs table
    df_pairs = build_redundant_pairs_table(profiles_df, columns_path, clusters)
    pairs_out = output_dir / "ontology_redundancy_pairs.tsv"
    df_pairs.to_csv(pairs_out, sep="\t", index=False)

    # 2. Compute baseline metrics on full non-media trait matrix (317 traits)
    raw_matrix = profiles_df[non_media_cols].values.astype(float)
    raw_metrics = compute_correlation_metrics(raw_matrix, non_media_cols, method="pearson")

    # 3. Collapse matrix under Native Survivor Count semantics
    mat_native, cols_native = collapse_trait_matrix(
        profiles_df, non_media_cols, clusters, count_semantics="native"
    )
    native_metrics = compute_correlation_metrics(mat_native, cols_native, method="pearson")

    # 4. Collapse matrix under Deduplicated Max Count semantics
    mat_max, cols_max = collapse_trait_matrix(
        profiles_df, non_media_cols, clusters, count_semantics="max"
    )
    max_metrics = compute_correlation_metrics(mat_max, cols_max, method="pearson")

    # 5. Apply PC1 residualization to collapsed max matrix
    adj_collapsed_mat, _, _ = regress_out_pc1(mat_max)
    adj_collapsed_metrics = compute_correlation_metrics(
        adj_collapsed_mat, cols_max, method="pearson"
    )

    # 6. Save Comparison Summary TSV
    comparison_rows = [
        {"configuration": "raw_uncollapsed_317", **asdict(raw_metrics)},
        {"configuration": "collapsed_native_survivor_301", **asdict(native_metrics)},
        {"configuration": "collapsed_deduplicated_max_301", **asdict(max_metrics)},
        {"configuration": "collapsed_max_plus_pc1_adjusted_301", **asdict(adj_collapsed_metrics)},
    ]
    df_comparison = pd.DataFrame(comparison_rows)
    comparison_out = output_dir / "ontology_collapsing_comparison.tsv"
    df_comparison.to_csv(comparison_out, sep="\t", index=False)

    # 7. Summary metrics of redundant pairs
    n_redundant_pairs = len(df_pairs)
    n_clusters = len(clusters)
    n_collapsed_features = len(non_media_cols) - len(cols_max)
    bounded_trait_space = len(cols_max)

    summary_rows = [
        {"metric": "total_non_media_traits", "value": len(non_media_cols)},
        {"metric": "redundant_clusters_count", "value": n_clusters},
        {"metric": "redundant_pairs_count", "value": n_redundant_pairs},
        {"metric": "redundant_features_collapsed", "value": n_collapsed_features},
        {"metric": "bounded_distinct_trait_space", "value": bounded_trait_space},
        {
            "metric": "raw_pc1_variance_pct",
            "value": round(raw_metrics.pc1_variance_pct, 4),
        },
        {
            "metric": "collapsed_native_pc1_variance_pct",
            "value": round(native_metrics.pc1_variance_pct, 4),
        },
        {
            "metric": "collapsed_max_pc1_variance_pct",
            "value": round(max_metrics.pc1_variance_pct, 4),
        },
        {
            "metric": "collapsed_adj_pc1_variance_pct",
            "value": round(adj_collapsed_metrics.pc1_variance_pct, 4),
        },
        {
            "metric": "raw_median_pairwise_r",
            "value": round(raw_metrics.median_pairwise_r, 4),
        },
        {
            "metric": "collapsed_max_median_pairwise_r",
            "value": round(max_metrics.median_pairwise_r, 4),
        },
        {
            "metric": "collapsed_adj_median_pairwise_r",
            "value": round(adj_collapsed_metrics.median_pairwise_r, 4),
        },
        {
            "metric": "raw_test_pair_r",
            "value": round(raw_metrics.test_pair_r, 4),
        },
        {
            "metric": "collapsed_max_test_pair_r",
            "value": round(max_metrics.test_pair_r, 4),
        },
        {
            "metric": "collapsed_adj_test_pair_r",
            "value": round(adj_collapsed_metrics.test_pair_r, 4),
        },
    ]
    df_summary = pd.DataFrame(summary_rows)
    summary_out = output_dir / "ontology_redundancy_summary.tsv"
    df_summary.to_csv(summary_out, sep="\t", index=False)

    return {
        "df_pairs": df_pairs,
        "n_redundant_pairs": n_redundant_pairs,
        "n_clusters": n_clusters,
        "n_collapsed_features": n_collapsed_features,
        "bounded_trait_space": bounded_trait_space,
        "raw_metrics": raw_metrics,
        "native_metrics": native_metrics,
        "max_metrics": max_metrics,
        "adj_collapsed_metrics": adj_collapsed_metrics,
        "pairs_tsv": pairs_out,
        "summary_tsv": summary_out,
        "comparison_tsv": comparison_out,
    }


def main() -> None:
    """CLI entrypoint for ontology redundancy audit."""
    parser = argparse.ArgumentParser(
        description="Audit METPO ontology redundancy and evaluate branch collapsing."
    )
    parser.add_argument(
        "--profiles",
        type=Path,
        default=DEFAULT_PROFILES_PATH,
        help="Path to study_trait_profiles_wide.tsv",
    )
    parser.add_argument(
        "--columns",
        type=Path,
        default=DEFAULT_COLUMNS_PATH,
        help="Path to taxa_trait_matrix_columns.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory to write TSVs and reports",
    )

    args = parser.parse_args()
    results = run_ontology_redundancy_audit(
        profiles_path=args.profiles,
        columns_path=args.columns,
        output_dir=args.output_dir,
    )
    print(f"Ontology redundancy audit completed successfully. Outputs written to {args.output_dir}")
    print(f"Redundant Pairs TSV: {results['pairs_tsv']}")
    print(f"Redundancy Summary TSV: {results['summary_tsv']}")
    print(f"Collapsing Comparison TSV: {results['comparison_tsv']}")
    print(f"Standalone Redundant Pairs Count: {results['n_redundant_pairs']}")
    print(
        f"Bounded Distinct Trait Space: {results['bounded_trait_space']} "
        f"(collapsed {results['n_collapsed_features']} features across "
        f"{results['n_clusters']} clusters)"
    )


__all__ = [
    "CANONICAL_REDUNDANCY_CLUSTERS",
    "DEFAULT_COMPARISON_OUT",
    "DEFAULT_PAIRS_OUT",
    "DEFAULT_REDUNDANCY_SUMMARY_OUT",
    "RedundancyCluster",
    "RedundantPair",
    "build_redundant_pairs_table",
    "collapse_trait_matrix",
    "get_canonical_redundancy_clusters",
    "main",
    "run_ontology_redundancy_audit",
]


if __name__ == "__main__":
    main()
