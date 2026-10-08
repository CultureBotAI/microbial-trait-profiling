#!/usr/bin/env python3
"""Create stacked bar plots for trait enrichment across X (traits) and Y (profiles) axes."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform


def parse_args():
    parser = argparse.ArgumentParser(description="Plot enrichment stacked bar graphs.")
    parser.add_argument(
        "--results",
        type=Path,
        default=Path("data/derived/trait_enrichment_results.tsv"),
        help="Path to trait enrichment results TSV.",
    )
    parser.add_argument(
        "--universe",
        type=Path,
        default=Path("data/derived/study_trait_profiles_wide.tsv"),
        help="Path to the full study trait profiles wide TSV.",
    )
    parser.add_argument(
        "--experiments",
        type=Path,
        default=Path("data/derived/experiments.tsv"),
        help="Path to the experiments metadata TSV.",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/derived"), help="Directory to save the plots."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Load enrichment results
    print(f"Loading enrichment results from {args.results}...")
    if not args.results.exists():
        print(f"Error: Results file {args.results} not found.")
        return
    results_df = pd.read_csv(args.results, sep="\t", comment="#")
    # Exclude MediaDive growth-medium traits for consistency across analysis arms
    results_df = results_df[~results_df["trait"].str.contains("mediadive.medium:", na=False)]

    # Load experiments metadata for stacking by condition/anatomy
    print(f"Loading experiment metadata from {args.experiments}...")
    if args.experiments.exists():
        exp_df = pd.read_csv(args.experiments, sep="\t")
        # Ensure we only have required columns to avoid duplicates
        exp_df = exp_df[["study_id", "experiment_id", "condition", "body_site"]].drop_duplicates()

        # Merge with results
        results_with_meta = results_df.merge(exp_df, on=["study_id", "experiment_id"], how="left")

        # Fill missing with "Unknown"
        results_with_meta["condition"] = results_with_meta["condition"].fillna("Unknown")
        results_with_meta["body_site"] = results_with_meta["body_site"].fillna("Unknown")

        # Create combined category
        results_with_meta["category"] = (
            results_with_meta["condition"] + " (" + results_with_meta["body_site"] + ")"
        )
    else:
        print(
            f"Warning: Experiments file {args.experiments} not found. "
            "Skipping metadata-based plots."
        )
        results_with_meta = None

    # 1. Y-Axis: Profile Depth (How many traits are enriched in each profile?)
    # Group by profile and direction
    # We define a unique profile as the combination of study, experiment, and direction.
    results_df["profile_id"] = (
        results_df["study_id"].astype(str)
        + "|"
        + results_df["experiment_id"].astype(str)
        + "|"
        + results_df["inc_dec"].astype(str)
    )
    profile_trait_counts = (
        results_df.groupby("profile_id")
        .agg(
            {"study_id": "first", "experiment_id": "first", "inc_dec": "first", "trait": "nunique"}
        )
        .reset_index()
    )
    profile_trait_counts.columns = [
        "profile_id",
        "study_id",
        "experiment_id",
        "inc_dec",
        "trait_count",
    ]

    # Filter out profiles with 0 hits (remove tails)
    active_profiles = profile_trait_counts[profile_trait_counts["trait_count"] > 0].copy()

    # Define bins for the profile trait counts
    max_traits = active_profiles["trait_count"].max()
    num_bins = 30
    bins_y = np.linspace(1, max_traits, num_bins + 1)

    # Assign each profile to a bin
    active_profiles["bin_idx"] = pd.cut(
        active_profiles["trait_count"], bins=bins_y, labels=False, include_lowest=True
    )

    # Group by bin and inc_dec to get counts for the stack
    profile_binned = active_profiles.groupby(["bin_idx", "inc_dec"]).size().unstack(fill_value=0)

    # Ensure both directions exist in columns
    for d in ["inc", "dec"]:
        if d not in profile_binned.columns:
            profile_binned[d] = 0

    # Sort columns for consistent stacking
    profile_binned = profile_binned[["inc", "dec"]]

    # 2. X-Axis: Trait Frequency (How many profiles is each trait enriched in?)
    # Count unique profiles for each trait separated by direction.
    # This ensures traits enriched in both directions are counted twice.
    trait_inc = results_df[results_df["inc_dec"] == "inc"].groupby("trait")["profile_id"].nunique()
    trait_dec = results_df[results_df["inc_dec"] == "dec"].groupby("trait")["profile_id"].nunique()
    trait_total = results_df.groupby("trait")["profile_id"].nunique()

    trait_data = pd.DataFrame({"inc": trait_inc, "dec": trait_dec, "total": trait_total}).fillna(0)

    # Filter out traits with 0 hits (remove tails)
    active_traits = trait_data[trait_data["total"] > 0].copy()

    # Define bins for trait frequency
    max_freq = active_traits["total"].max()
    bins_x = np.linspace(1, max_freq, num_bins + 1)

    # Assign each trait to a bin based on total frequency
    active_traits["bin_idx"] = pd.cut(
        active_traits["total"], bins=bins_x, labels=False, include_lowest=True
    )

    # For each bin, sum the 'inc' and 'dec' enrichment hits
    trait_binned = (
        active_traits.groupby("bin_idx")[["inc", "dec"]]
        .sum()
        .rename(columns={"inc": "Increased", "dec": "Decreased"})[["Increased", "Decreased"]]
    )
    profile_binned = profile_binned.rename(columns={"inc": "Increased", "dec": "Decreased"})[
        ["Increased", "Decreased"]
    ]

    # Plotting
    sns.set_theme(style="whitegrid")
    # Harmonized palette: Red/Coral for Increased, Blue for Decreased
    colors = {"Increased": "#c0392b", "Decreased": "#2874a6"}
    tick_indices = range(0, num_bins, 2)

    def save_histogram(
        df,
        title,
        xlabel,
        ylabel,
        path,
        use_log=True,
        is_stacked=True,
        color_map=None,
        bin_labels=None,
        legend_title="Signature Direction",
        subtitle=(
            "Discrete urn model (M=5,184). Annotation representation "
            "in curated signatures, not patient-level abundance."
        ),
    ):
        plt.figure(figsize=(12, 7))
        df.plot(kind="bar", stacked=is_stacked, width=0.8, color=color_map, ax=plt.gca())
        full_title = f"{title}\n{subtitle}" if subtitle else title
        plt.title(full_title, pad=15, fontsize=12)
        plt.xlabel(xlabel, fontsize=11)
        plt.ylabel(ylabel, fontsize=11)
        if use_log:
            plt.yscale("log")

        if bin_labels:
            plt.xticks(tick_indices, [bin_labels[i] for i in tick_indices], rotation=45)

        plt.legend(
            title=legend_title,
            bbox_to_anchor=(1.05, 1) if len(df.columns) > 5 else None,
            loc="upper left" if len(df.columns) > 5 else None,
        )
        plt.grid(axis="y", which="both", linestyle="-", alpha=0.3)
        plt.savefig(path, dpi=300, bbox_inches="tight")
        print(f"Wrote plot to {path}")
        plt.close()

    # 1. Profile Depth Plots (Y-Axis)
    bin_labels_y = [f"{int(bins_y[i])}-{int(bins_y[i + 1])}" for i in range(len(bins_y) - 1)]
    profile_binned = profile_binned.reindex(range(num_bins), fill_value=0)
    save_histogram(
        profile_binned,
        f"Profile Enrichment Depth (Y-Axis)\n{len(active_profiles)} active profiles",
        "Number of Enriched Traits (Binned)",
        "Number of Profiles",
        args.output_dir / "profile_enrichment_y_histogram.png",
        color_map=[colors["Increased"], colors["Decreased"]],
        bin_labels=bin_labels_y,
    )

    # 2. Trait Frequency Plots (X-Axis)
    bin_labels_x = [f"{int(bins_x[i])}-{int(bins_x[i + 1])}" for i in range(len(bins_x) - 1)]
    trait_binned = trait_binned.reindex(range(num_bins), fill_value=0)
    save_histogram(
        trait_binned,
        f"Trait Enrichment Frequency (X-Axis)\n{len(active_traits)} active traits",
        "Number of Profiles Enriched In (Binned)",
        "Number of Profile Hits",
        args.output_dir / "trait_enrichment_x_histogram.png",
        color_map=[colors["Increased"], colors["Decreased"]],
        bin_labels=bin_labels_x,
    )

    # 3. X-Axis: Trait Frequency Stacked by Condition/Anatomy
    if results_with_meta is not None:
        print("Creating Condition/Anatomy split plots...")
        top_n = 12
        cat_counts = results_with_meta["category"].value_counts()
        top_cats = cat_counts.head(top_n).index.tolist()

        results_with_meta["plot_category"] = results_with_meta["category"].apply(
            lambda x: x if x in top_cats else "Other"
        )

        trait_to_total = active_traits["total"].to_dict()
        results_with_meta["trait_freq"] = results_with_meta["trait"].map(trait_to_total)
        results_with_meta["bin_idx"] = pd.cut(
            results_with_meta["trait_freq"], bins=bins_x, labels=False, include_lowest=True
        )

        # Helper to plot both count and 100% stacked for a subset
        def plot_category_histograms(df_subset, suffix):
            # Group by bin and plot_category
            binned = df_subset.groupby(["bin_idx", "plot_category"]).size().unstack(fill_value=0)
            binned = binned.reindex(range(num_bins), fill_value=0)

            # Sort columns
            other_col = [c for c in binned.columns if c == "Other"]
            main_cols = sorted([c for c in binned.columns if c != "Other"])
            binned = binned[main_cols + other_col]

            cat_colors = sns.color_palette("husl", len(binned.columns))
            disp_suffix = {
                "all": "All Signatures",
                "inc": "Increased Signatures",
                "dec": "Decreased Signatures",
            }.get(suffix, suffix)

            # 1. Normal count plot (log scale)
            save_histogram(
                binned,
                f"Trait Enrichment Frequency by Condition ({disp_suffix})\nStacked Counts",
                "Number of Profiles Enriched In (Binned)",
                "Number of Profile Hits",
                args.output_dir / f"trait_enrichment_x_condition_{suffix}_counts.png",
                color_map=cat_colors,
                bin_labels=bin_labels_x,
                legend_title="Condition (Body Site)",
            )

            # 2. 100% Stacked (stretched)
            binned_stretched = binned.div(binned.sum(axis=1), axis=0).fillna(0)
            title_str = (
                f"Trait Enrichment Frequency by Condition ({disp_suffix})\n"
                "100% Stacked (Relative Proportion)"
            )
            save_histogram(
                binned_stretched,
                title_str,
                "Number of Profiles Enriched In (Binned)",
                "Relative Proportion (0-1)",
                args.output_dir / f"trait_enrichment_x_condition_{suffix}_stretched.png",
                use_log=False,
                color_map=cat_colors,
                bin_labels=bin_labels_x,
                legend_title="Condition (Body Site)",
            )
            return binned

        # Plot for All, Inc, Dec
        print("Generating 'All' condition plots...")
        plot_category_histograms(results_with_meta, "all")

        print("Generating 'Increased' condition plots...")
        plot_category_histograms(results_with_meta[results_with_meta["inc_dec"] == "inc"], "inc")

        print("Generating 'Decreased' condition plots...")
        plot_category_histograms(results_with_meta[results_with_meta["inc_dec"] == "dec"], "dec")

        # 4. Correlation Matrix (Euclidean Distance Heatmap)
        print("Creating Euclidean distance heatmap for categories (including Direction)...")
        # We want to compare (Category + Direction) profiles with readable direction names
        dir_labels = results_with_meta["inc_dec"].map({"inc": "Increased", "dec": "Decreased"})
        results_with_meta["cat_dir"] = results_with_meta["plot_category"] + " [" + dir_labels + "]"

        # Group by bin and cat_dir
        cat_dir_binned = (
            results_with_meta.groupby(["bin_idx", "cat_dir"]).size().unstack(fill_value=0)
        )
        cat_dir_binned = cat_dir_binned.reindex(range(num_bins), fill_value=0)

        # Normalize vectors to compare 'shape' rather than absolute magnitude
        # This prevents high-count categories from dominating the distance calculation
        cat_dir_norm = cat_dir_binned.div(cat_dir_binned.sum(axis=0), axis=1).fillna(0)

        # Transpose so rows are observations (categories)
        obs_matrix = cat_dir_norm.T

        if not obs_matrix.empty and len(obs_matrix) > 1:
            # Calculate Euclidean distance on normalized vectors
            dist_vec = pdist(obs_matrix, metric="euclidean")
            dist_matrix = squareform(dist_vec)
            dist_df = pd.DataFrame(dist_matrix, index=obs_matrix.index, columns=obs_matrix.index)

            # Perform linkage
            linkage_matrix = linkage(dist_vec, method="average")

            plt.figure(figsize=(16, 14))
            g = sns.clustermap(
                dist_df,
                row_linkage=linkage_matrix,
                col_linkage=linkage_matrix,
                cmap="rocket_r",
                linewidths=0.5,
                annot=len(dist_df) < 25,
                fmt=".2f",
                cbar_kws={"label": "Euclidean Distance", "orientation": "horizontal"},
                cbar_pos=(0.03, 0.91, 0.12, 0.018),
            )
            g.ax_heatmap.set_xlabel("")
            g.ax_heatmap.set_ylabel("")
            g.fig.suptitle(
                "Similarity of (Condition + Direction) Enrichment Profiles\n"
                "Euclidean Distance on Normalized Frequency Vectors\n"
                "Discrete urn model (M=5,184). Annotation representation in curated signatures.",
                y=0.99,
                fontsize=13,
            )
            g.fig.subplots_adjust(top=0.89, right=0.90)

            heatmap_path = (
                args.output_dir / "trait_enrichment_category_direction_distance_heatmap.png"
            )
            plt.savefig(heatmap_path, dpi=300, bbox_inches="tight")
            print(f"Wrote detailed distance heatmap to {heatmap_path}")
            plt.close()

    # Save summary counts
    active_traits.to_csv(args.output_dir / "trait_enrichment_counts_stacked.tsv", sep="\t")
    active_profiles.to_csv(args.output_dir / "profile_enrichment_counts_stacked.tsv", sep="\t")


if __name__ == "__main__":
    main()
