#!/usr/bin/env python3
"""Create stacked bar plots for trait enrichment across X (traits) and Y (profiles) axes."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(description="Plot enrichment stacked bar graphs.")
    parser.add_argument(
        "--results",
        type=Path,
        default=Path("data/derived/trait_enrichment_results.tsv"),
        help="Path to trait enrichment results TSV."
    )
    parser.add_argument(
        "--universe",
        type=Path,
        default=Path("data/derived/study_trait_profiles_wide.tsv"),
        help="Path to the full study trait profiles wide TSV."
    )
    parser.add_argument(
        "--experiments",
        type=Path,
        default=Path("data/derived/experiments.tsv"),
        help="Path to the experiments metadata TSV."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/derived"),
        help="Directory to save the plots."
    )
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Load enrichment results
    print(f"Loading enrichment results from {args.results}...")
    if not args.results.exists():
        print(f"Error: Results file {args.results} not found.")
        return
    results_df = pd.read_csv(args.results, sep="\t")

    # Load experiments metadata for stacking by condition/anatomy
    print(f"Loading experiment metadata from {args.experiments}...")
    if args.experiments.exists():
        exp_df = pd.read_csv(args.experiments, sep="\t")
        # Ensure we only have the columns we need to avoid duplicates if study_tsv was run multiple times
        exp_df = exp_df[["study_id", "experiment_id", "condition", "body_site"]].drop_duplicates()
        
        # Merge with results
        results_with_meta = results_df.merge(exp_df, on=["study_id", "experiment_id"], how="left")
        
        # Fill missing with "Unknown"
        results_with_meta["condition"] = results_with_meta["condition"].fillna("Unknown")
        results_with_meta["body_site"] = results_with_meta["body_site"].fillna("Unknown")
        
        # Create combined category
        results_with_meta["category"] = results_with_meta["condition"] + " (" + results_with_meta["body_site"] + ")"
    else:
        print(f"Warning: Experiments file {args.experiments} not found. Skipping metadata-based plots.")
        results_with_meta = None

    # 1. Y-Axis: Profile Depth (How many traits are enriched in each profile?)
    # Group by profile and direction
    profile_trait_counts = results_df.groupby(["study_id", "experiment_id", "inc_dec"])["trait"].nunique().reset_index()
    profile_trait_counts.columns = ["study_id", "experiment_id", "inc_dec", "trait_count"]
    
    # Filter out profiles with 0 hits (remove tails)
    active_profiles = profile_trait_counts[profile_trait_counts["trait_count"] > 0].copy()
    
    # Define bins for the profile trait counts
    max_traits = active_profiles["trait_count"].max()
    num_bins = 30
    bins_y = np.linspace(1, max_traits, num_bins + 1)
    
    # Assign each profile to a bin
    active_profiles["bin_idx"] = pd.cut(active_profiles["trait_count"], bins=bins_y, labels=False, include_lowest=True)
    
    # Group by bin and inc_dec to get counts for the stack
    profile_binned = active_profiles.groupby(["bin_idx", "inc_dec"]).size().unstack(fill_value=0)
    
    # Ensure both directions exist in columns
    for d in ["inc", "dec"]:
        if d not in profile_binned.columns:
            profile_binned[d] = 0
    
    # Sort columns for consistent stacking
    profile_binned = profile_binned[["inc", "dec"]]
    
    # 2. X-Axis: Trait Frequency (How many profiles is each trait enriched in?)
    # Count occurrences for each trait separated by profile direction
    trait_inc = results_df[results_df["inc_dec"] == "inc"].groupby("trait")["experiment_id"].nunique()
    trait_dec = results_df[results_df["inc_dec"] == "dec"].groupby("trait")["experiment_id"].nunique()
    trait_total = results_df.groupby("trait")["experiment_id"].nunique()
    
    trait_data = pd.DataFrame({
        "inc": trait_inc,
        "dec": trait_dec,
        "total": trait_total
    }).fillna(0)
    
    # Filter out traits with 0 hits (remove tails)
    active_traits = trait_data[trait_data["total"] > 0].copy()
    
    # Define bins for trait frequency
    max_freq = active_traits["total"].max()
    bins_x = np.linspace(1, max_freq, num_bins + 1)
    
    # Assign each trait to a bin based on total frequency
    active_traits["bin_idx"] = pd.cut(active_traits["total"], bins=bins_x, labels=False, include_lowest=True)
    
    # For each bin, sum the 'inc' and 'dec' enrichment hits
    trait_binned = active_traits.groupby("bin_idx")[["inc", "dec"]].sum()
    
    # Plotting
    sns.set_theme(style="whitegrid")
    colors = {"inc": "#3498db", "dec": "#e74c3c"} # Blue for Increase, Red for Decrease
    
    # Profile Depth Plot (Y-Axis)
    plt.figure(figsize=(12, 7))
    bin_labels_y = [f"{int(bins_y[i])}-{int(bins_y[i+1])}" for i in range(len(bins_y)-1)]
    
    # Reindex to ensure all bins are represented even if empty (within range)
    profile_binned = profile_binned.reindex(range(num_bins), fill_value=0)
    
    ax_y = profile_binned.plot(
        kind="bar", 
        stacked=True, 
        width=0.8, 
        color=[colors["inc"], colors["dec"]],
        ax=plt.gca()
    )
    
    plt.title(f"Profile Enrichment Depth (Y-Axis)\n{len(active_profiles)} active profiles (zero-hit tail removed)", pad=15)
    plt.xlabel("Number of Enriched Traits (Binned)")
    plt.ylabel("Number of Profiles (log10 scale)")
    plt.yscale("log")
    
    # Set x-labels every 5 bins to avoid overcrowding
    tick_indices = range(0, num_bins, 2)
    plt.xticks(tick_indices, [bin_labels_y[i] for i in tick_indices], rotation=45)
    
    plt.legend(title="Profile Direction")
    plt.grid(axis='y', which='both', linestyle='-', alpha=0.3)
    
    profile_plot_path = args.output_dir / "profile_enrichment_y_histogram.png"
    plt.savefig(profile_plot_path, dpi=300, bbox_inches="tight")
    print(f"Wrote stacked profile plot to {profile_plot_path}")
    plt.close()
    
    # Trait Frequency Plot (X-Axis)
    plt.figure(figsize=(12, 7))
    bin_labels_x = [f"{int(bins_x[i])}-{int(bins_x[i+1])}" for i in range(len(bins_x)-1)]
    
    trait_binned = trait_binned.reindex(range(num_bins), fill_value=0)
    
    ax_x = trait_binned.plot(
        kind="bar", 
        stacked=True, 
        width=0.8, 
        color=[colors["inc"], colors["dec"]],
        ax=plt.gca()
    )
    
    plt.title(f"Trait Enrichment Frequency (X-Axis)\n{len(active_traits)} active traits (zero-hit tail removed)", pad=15)
    plt.xlabel("Number of Profiles Enriched In (Binned)")
    plt.ylabel("Total Enrichment Hits (log10 scale)")
    plt.yscale("log")
    
    plt.xticks(tick_indices, [bin_labels_x[i] for i in tick_indices], rotation=45)
    
    plt.legend(title="Enrichment Direction")
    plt.grid(axis='y', which='both', linestyle='-', alpha=0.3)
    
    trait_plot_path = args.output_dir / "trait_enrichment_x_histogram.png"
    plt.savefig(trait_plot_path, dpi=300, bbox_inches="tight")
    print(f"Wrote stacked trait plot to {trait_plot_path}")
    plt.close()
    
    # 3. X-Axis: Trait Frequency Stacked by Condition/Anatomy (Separate Plot)
    if results_with_meta is not None:
        print("Creating Trait Frequency plot stacked by Condition and Anatomy...")
        plt.figure(figsize=(14, 8))
        
        # Identify top categories for the legend, group the rest as "Other"
        top_n = 12
        cat_counts = results_with_meta["category"].value_counts()
        top_cats = cat_counts.head(top_n).index.tolist()
        
        results_with_meta["plot_category"] = results_with_meta["category"].apply(
            lambda x: x if x in top_cats else "Other"
        )
        
        # We need the same bins as the previous X plot
        # active_traits["total"] contains the frequency per trait
        trait_to_total = active_traits["total"].to_dict()
        results_with_meta["trait_freq"] = results_with_meta["trait"].map(trait_to_total)
        
        # Assign to bins
        results_with_meta["bin_idx"] = pd.cut(
            results_with_meta["trait_freq"],
            bins=bins_x,
            labels=False,
            include_lowest=True
        )
        
        # Group by bin and plot_category
        trait_cat_binned = results_with_meta.groupby(["bin_idx", "plot_category"]).size().unstack(fill_value=0)
        
        # Reindex to ensure all bins are represented
        trait_cat_binned = trait_cat_binned.reindex(range(num_bins), fill_value=0)
        
        # Sort columns so "Other" is last and others are consistent (e.g. by total frequency)
        other_col = [c for c in trait_cat_binned.columns if c == "Other"]
        main_cols = sorted([c for c in trait_cat_binned.columns if c != "Other"])
        trait_cat_binned = trait_cat_binned[main_cols + other_col]
        
        # Use a more diverse color palette for many categories
        cat_colors = sns.color_palette("husl", len(trait_cat_binned.columns))
        
        ax_x_cat = trait_cat_binned.plot(
            kind="bar",
            stacked=True,
            width=0.8,
            color=cat_colors,
            ax=plt.gca()
        )
        
        plt.title(f"Trait Enrichment Frequency by Condition & Anatomy (X-Axis)\nStacked by top {top_n} categories", pad=15)
        plt.xlabel("Number of Profiles Enriched In (Binned)")
        plt.ylabel("Total Enrichment Hits (log10 scale)")
        plt.yscale("log")
        
        plt.xticks(tick_indices, [bin_labels_x[i] for i in tick_indices], rotation=45)
        
        # Place legend outside to the right
        plt.legend(title="Condition (Anatomy)", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(axis='y', which='both', linestyle='-', alpha=0.3)
        
        trait_cat_plot_path = args.output_dir / "trait_enrichment_x_condition_site_histogram.png"
        plt.savefig(trait_cat_plot_path, dpi=300, bbox_inches="tight")
        print(f"Wrote stacked trait condition plot to {trait_cat_plot_path}")
        plt.close()

    # Save summary counts
    active_traits.to_csv(args.output_dir / "trait_enrichment_counts_stacked.tsv", sep="\t")
    active_profiles.to_csv(args.output_dir / "profile_enrichment_counts_stacked.tsv", sep="\t")

if __name__ == "__main__":
    main()

