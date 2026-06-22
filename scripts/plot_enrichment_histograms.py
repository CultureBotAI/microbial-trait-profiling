#!/usr/bin/env python3
"""Create histogram plots for trait enrichment across X (traits) and Y (profiles) axes."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(description="Plot enrichment histograms.")
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
        "--output-dir",
        type=Path,
        default=Path("data/derived"),
        help="Directory to save the plots."
    )
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Load universe to get all possible profiles and traits
    print(f"Loading universe from {args.universe}...")
    # Read only the header first to get trait columns
    universe_header = pd.read_csv(args.universe, sep="\t", nrows=0)
    all_traits = [col for col in universe_header.columns if col not in {"study_id", "experiment_id", "inc_dec"}]
    
    # Read universe to get all profiles
    universe_df = pd.read_csv(args.universe, sep="\t", usecols=["study_id", "experiment_id", "inc_dec"])
    universe_df["profile_id"] = universe_df["study_id"] + "|" + universe_df["experiment_id"] + "|" + universe_df["inc_dec"]
    all_profiles = universe_df["profile_id"].unique()
    
    print(f"Universe contains {len(all_traits)} traits and {len(all_profiles)} profiles.")
    
    # Load enrichment results
    print(f"Loading enrichment results from {args.results}...")
    results_df = pd.read_csv(args.results, sep="\t")
    results_df["profile_id"] = results_df["study_id"] + "|" + results_df["experiment_id"] + "|" + results_df["inc_dec"]
    
    # 1. X-Axis: Enriched Traits (how many profiles is each trait enriched in?)
    trait_counts = results_df.groupby("trait")["profile_id"].nunique()
    # Align with all traits in universe
    trait_counts = trait_counts.reindex(all_traits, fill_value=0)
    
    # 2. Y-Axis: Enriched Profiles (how many traits are enriched in each profile?)
    profile_counts = results_df.groupby("profile_id")["trait"].nunique()
    # Align with all profiles in universe
    profile_counts = profile_counts.reindex(all_profiles, fill_value=0)
    
    # Plotting
    sns.set_theme(style="whitegrid")

    # Trait Histogram
    plt.figure(figsize=(12, 7))
    trait_zeros = (trait_counts == 0).sum()
    
    # Use plt.hist with rwidth to force gaps between bars and avoid "step plot" look
    plt.hist(
        trait_counts,
        bins=30,
        color="#3498db",
        edgecolor="black",
        linewidth=1.0,
        rwidth=0.8,
        log=True
    )
    
    plt.title(f"Trait Enrichment Frequency (X-Axis)\n{len(all_traits)} traits total ({trait_zeros} have 0 hits)", pad=15)
    plt.xlabel("Number of Profiles Enriched In")
    plt.ylabel("Number of Traits (log10 scale)")
    plt.grid(axis='y', which='major', linestyle='-', alpha=0.3)
    plt.grid(axis='y', which='minor', linestyle=':', alpha=0.2)
    
    trait_plot_path = args.output_dir / "trait_enrichment_x_histogram.png"
    plt.savefig(trait_plot_path, dpi=300, bbox_inches="tight")
    print(f"Wrote trait histogram to {trait_plot_path}")
    plt.close()
    
    # Profile Histogram
    plt.figure(figsize=(12, 7))
    profile_zeros = (profile_counts == 0).sum()
    
    plt.hist(
        profile_counts,
        bins=30,
        color="#e74c3c",
        edgecolor="black",
        linewidth=1.0,
        rwidth=0.8,
        log=True
    )
    
    plt.title(f"Profile Enrichment Depth (Y-Axis)\n{len(all_profiles)} profiles total ({profile_zeros} have 0 hits)", pad=15)
    plt.xlabel("Number of Enriched Traits")
    plt.ylabel("Number of Profiles (log10 scale)")
    plt.grid(axis='y', which='major', linestyle='-', alpha=0.3)
    plt.grid(axis='y', which='minor', linestyle=':', alpha=0.2)
    
    profile_plot_path = args.output_dir / "profile_enrichment_y_histogram.png"
    plt.savefig(profile_plot_path, dpi=300, bbox_inches="tight")
    print(f"Wrote profile histogram to {profile_plot_path}")
    plt.close()
    
    # Save raw counts for verification
    trait_counts_df = trait_counts.reset_index()
    trait_counts_df.columns = ["trait", "enrichment_count"]
    trait_counts_path = args.output_dir / "trait_enrichment_counts.tsv"
    trait_counts_df.to_csv(trait_counts_path, sep="\t", index=False)
    print(f"Wrote trait counts to {trait_counts_path}")
    
    profile_counts_df = profile_counts.reset_index()
    profile_counts_df.columns = ["profile_id", "enrichment_count"]
    profile_counts_path = args.output_dir / "profile_enrichment_counts.tsv"
    profile_counts_df.to_csv(profile_counts_path, sep="\t", index=False)
    print(f"Wrote profile counts to {profile_counts_path}")

if __name__ == "__main__":
    main()
