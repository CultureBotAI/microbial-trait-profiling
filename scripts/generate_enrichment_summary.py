#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
import re

def main():
    # Paths
    trait_counts_path = Path("bsdbng/data/derived/trait_enrichment_counts.tsv")
    profile_counts_path = Path("bsdbng/data/derived/profile_enrichment_counts.tsv")
    mapping_path = Path("metpo_phenotype_mapping.tsv")
    matrix_path = Path("bsdbng/data/derived/trait_enrichment_clustered_heatmap_matrix.tsv")
    experiments_path = Path("bsdbng/data/raw/experiments.csv")
    output_text_path = Path("bsdbng/data/derived/enrichment_summary_top10.txt")
    output_img_path = Path("bsdbng/data/derived/enrichment_infographic.png")

    # 1. Load Data
    print("Loading data...")
    trait_counts = pd.read_csv(trait_counts_path, sep="\t")
    profile_counts = pd.read_csv(profile_counts_path, sep="\t")
    mapping = pd.read_csv(mapping_path, sep="\t").set_index("METPO_ID")["Description"].to_dict()
    
    # Load basic metadata from matrix
    matrix_meta = pd.read_csv(matrix_path, sep="\t", usecols=["profile_id", "study_id", "experiment_id", "inc_dec"])
    
    # Load raw experiments for Body Site and Condition
    if experiments_path.exists():
        print("Enriching metadata from experiments.csv...")
        exdf = pd.read_csv(experiments_path, dtype=str)
        # Normalize columns
        exdf.columns = [c.strip().replace('"', '') for c in exdf.columns]
        
        # Create a key that matches matrix IDs
        exdf['study_key'] = "bsdb:" + exdf['Study'].astype(str).str.strip()
        # Extract experiment number, e.g., "Experiment 1" -> "1"
        exdf['exp_num'] = exdf['Experiment page name'].str.extract(r'(\d+)', expand=False)
        exdf['exp_key'] = exdf['study_key'] + "-" + exdf['exp_num'].fillna("1")
        
        # Keep only one row per experiment key to avoid non-unique index error
        exdf_unique = exdf.drop_duplicates(subset='exp_key')
        lookup = exdf_unique.set_index('exp_key')[['Body site', 'Condition']].to_dict('index')
        
        matrix_meta['body_site'] = matrix_meta['experiment_id'].map(lambda x: lookup.get(x, {}).get('Body site', 'Unknown'))
        matrix_meta['condition'] = matrix_meta['experiment_id'].map(lambda x: lookup.get(x, {}).get('Condition', 'Unknown'))
    else:
        print("experiments.csv not found, using 'Unknown' for metadata.")
        matrix_meta['body_site'] = "Unknown"
        matrix_meta['condition'] = "Unknown"

    matrix_meta = matrix_meta.set_index("profile_id")

    # 2. Top 10 Traits
    top_traits = trait_counts.sort_values("enrichment_count", ascending=False).head(10).copy()
    top_traits["description"] = top_traits["trait"].map(lambda x: mapping.get(x.split("|")[0], x))
    
    # 3. Top 10 Profiles
    top_profiles = profile_counts.sort_values("enrichment_count", ascending=False).head(10).copy()
    top_profiles = top_profiles.join(matrix_meta, on="profile_id")

    # 4. Additional Finding: Body Site Distribution
    all_profile_activity = profile_counts[profile_counts["enrichment_count"] > 0].join(matrix_meta, on="profile_id")
    # Clean body site names for the pie chart
    all_profile_activity['body_site'] = all_profile_activity['body_site'].fillna("Unknown")
    body_site_activity = all_profile_activity.groupby("body_site")["enrichment_count"].sum().sort_values(ascending=False)
    
    # 5. Additional Finding: Directionality Bias
    direction_activity = all_profile_activity.groupby("inc_dec")["enrichment_count"].sum()

    # 6. Write Text File
    print(f"Writing summary to {output_text_path}...")
    with open(output_text_path, "w", encoding="utf-8") as f:
        f.write("=== TOP 10 MOST ENRICHED TRAITS ===\n")
        f.write(f"{'Trait ID':<40} {'Enrichment Count':<20} {'Description'}\n")
        f.write("-" * 100 + "\n")
        for _, row in top_traits.iterrows():
            f.write(f"{row['trait']:<40} {row['enrichment_count']:<20} {row['description']}\n")
        
        f.write("\n\n=== TOP 10 PROFILES WITH MOST ENRICHED TRAITS ===\n")
        f.write(f"{'Profile ID':<80} {'Count':<10} {'Condition':<30} {'Body Site'}\n")
        f.write("-" * 140 + "\n")
        for _, row in top_profiles.iterrows():
            cond = str(row['condition'])[:28]
            site = str(row['body_site'])
            f.write(f"{row['profile_id']:<80} {row['enrichment_count']:<10} {cond:<30} {site}\n")
            
        f.write("\n\n=== SCIENTIFIC FINDINGS SUMMARY ===\n")
        f.write(f"Total Enriched Traits across all profiles: {trait_counts['enrichment_count'].sum():,}\n")
        if not all_profile_activity.empty:
            f.write(f"Average enrichment count per profile (active only): {all_profile_activity['enrichment_count'].mean():.2f}\n")
        
        if not body_site_activity.empty:
            f.write(f"Top Body Site by Total Enrichment: {body_site_activity.index[0]} ({body_site_activity.iloc[0]} hits)\n")
        
        bias_inc = direction_activity.get("inc", 0)
        bias_dec = direction_activity.get("dec", 0)
        bias_str = "Increased" if bias_inc > bias_dec else "Decreased"
        f.write(f"Directionality Bias: More traits found in '{bias_str}' profiles ({bias_inc:,} inc vs {bias_dec:,} dec).\n")

    # 7. Generate Infographic
    print(f"Generating infographic to {output_img_path}...")
    plt.style.use('seaborn-v0_8-whitegrid')
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.4, wspace=0.3)

    # Panel A: Top 10 Traits
    ax1 = fig.add_subplot(gs[0, 0])
    sns.barplot(data=top_traits, x="enrichment_count", y="description", hue="description", palette="viridis", ax=ax1, legend=False)
    ax1.set_title("A. Top 10 Most Enriched Traits", fontsize=14, fontweight='bold')
    ax1.set_xlabel("Number of Profiles Enriched In")
    ax1.set_ylabel("")

    # Panel B: Top 10 Profiles
    ax2 = fig.add_subplot(gs[0, 1])
    top_profiles["short_label"] = top_profiles.apply(
        lambda r: f"{str(r['condition'])[:15]}... ({str(r['study_id']).split(':')[-1]})", axis=1
    )
    sns.barplot(data=top_profiles, x="enrichment_count", y="short_label", hue="short_label", palette="magma", ax=ax2, legend=False)
    ax2.set_title("B. Profiles with Highest Enrichment Depth", fontsize=14, fontweight='bold')
    ax2.set_xlabel("Number of Enriched Traits")
    ax2.set_ylabel("")

    # Panel C: Body Site Distribution
    ax3 = fig.add_subplot(gs[1, 0])
    if not body_site_activity.empty:
        plot_sites = body_site_activity.head(8)
        ax3.pie(plot_sites, labels=plot_sites.index, autopct='%1.1f%%', colors=sns.color_palette("pastel"), startangle=140)
        ax3.set_title("C. Enrichment Activity by Body Site", fontsize=14, fontweight='bold')
    else:
        ax3.text(0.5, 0.5, "No Body Site Data Available", ha='center', va='center')

    # Panel D: Scientific Summary Text
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')
    
    top_trait_name = top_traits.iloc[0]['description'] if not top_traits.empty else "N/A"
    top_trait_count = top_traits.iloc[0]['enrichment_count'] if not top_traits.empty else 0
    top_site_name = body_site_activity.index[0] if not body_site_activity.empty else "N/A"

    summary_text = (
        "KEY SCIENTIFIC FINDINGS\n\n"
        f"• Global Scale: {len(trait_counts):,} traits analyzed across {len(profile_counts):,} study profiles.\n\n"
        f"• Prevalence: The most enriched trait, '{top_trait_name}',\n"
        f"  appears in {top_trait_count:,} distinct biological contexts.\n\n"
        f"• Hotspots: Microbial signatures in the '{top_site_name}'\n"
        "  show high density of trait-level changes.\n\n"
        f"• Consistency: {bias_inc:,} traits were enriched in 'Increased'\n"
        f"  taxa vs {bias_dec:,} in 'Decreased' taxa, suggesting a\n"
        "  potential gain-of-function dominance in dysbiosis."
    )
    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes, fontsize=12, va='top', ha='left',
             bbox=dict(boxstyle='round,pad=1', facecolor='aliceblue', alpha=0.5))
    ax4.set_title("D. Summary & Implications", fontsize=14, fontweight='bold')

    plt.suptitle("CultureBotAI: Trait Enrichment Landscape Analysis", fontsize=20, fontweight='bold', y=0.98)
    
    plt.savefig(output_img_path, dpi=300, bbox_inches="tight")
    print(f"Infographic saved to {output_img_path}")
    print("Done.")

if __name__ == "__main__":
    main()
