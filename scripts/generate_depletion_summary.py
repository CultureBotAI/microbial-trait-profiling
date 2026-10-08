#!/usr/bin/env python3
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def shorten_grows_in(desc):
    """Shortens long 'grows in' medium labels for presentation."""
    if "grows in:" in str(desc):
        match = re.search(r"grows in:\s*(.*?)(?:\s+mediadive\.medium:|\s+DSMZ Medium|$)", str(desc))
        if match:
            clean_name = match.group(1).strip()
            words = clean_name.split()
            if len(words) > 3:
                clean_name = " ".join(words[:3])
            return f"Grows in: {clean_name}"
    return desc


def main():
    # Paths
    results_path = Path("data/derived/trait_depletion_results.tsv")
    mapping_path = Path("metpo_phenotype_mapping.tsv")
    experiments_path = Path("data/raw/full_dump.csv")
    trait_meta_path = Path("data/derived/taxa_trait_matrix_columns.tsv")
    output_text_path = Path("data/derived/depletion_summary_top10.txt")
    output_img_path = Path("data/derived/depletion_infographic.png")

    # 1. Load Data
    print("Loading data...")
    if not results_path.exists():
        print(f"Error: {results_path} not found.", file=sys.stderr)
        return 1

    results_df = pd.read_csv(results_path, sep="\t", comment="#")
    # Exclude growth mediums
    results_df = results_df[~results_df["trait"].str.contains("mediadive.medium:", na=False)]

    # Calculate trait counts (how many profiles is it depleted in?)
    trait_counts = results_df.groupby("trait").size().reset_index(name="depletion_count")

    # Calculate profile counts
    results_df["profile_id"] = (
        results_df["study_id"].astype(str)
        + "|"
        + results_df["experiment_id"].astype(str)
        + "|"
        + results_df["inc_dec"].astype(str)
    )
    profile_counts = results_df.groupby("profile_id").size().reset_index(name="depletion_count")

    # Load phenotype mapping
    mapping = {}
    if mapping_path.exists():
        mapping = pd.read_csv(mapping_path, sep="\t").set_index("METPO_ID")["Description"].to_dict()

    # Load more detailed trait mapping from matrix columns
    if trait_meta_path.exists():
        print("Enriching trait labels from matrix columns...")
        trait_meta = pd.read_csv(trait_meta_path, sep="\t")
        for _, row in trait_meta.iterrows():
            mapping[row["column_id"]] = row["column_label"]

    # Load and filter experiments metadata (strictly human only)
    print(f"Loading and filtering metadata from {experiments_path}...")
    if experiments_path.exists():
        exdf = pd.read_csv(experiments_path, comment="#", dtype=str)
        exdf.columns = [c.strip().replace('"', "") for c in exdf.columns]
        exdf = exdf[exdf["Host species"].str.strip().str.lower() == "homo sapiens"].copy()
        exdf["study_key"] = "bsdb:" + exdf["Study"].astype(str).str.strip()
        exdf["exp_num"] = exdf["Experiment"].str.extract(r"(\d+)", expand=False)
        exdf["experiment_id"] = exdf["study_key"] + "-" + exdf["exp_num"].fillna("1")
        exdf["sample_size_0"] = pd.to_numeric(exdf["Group 0 sample size"], errors="coerce").fillna(
            0
        )
        exdf["sample_size_1"] = pd.to_numeric(exdf["Group 1 sample size"], errors="coerce").fillna(
            0
        )
        exdf["total_samples"] = exdf["sample_size_0"] + exdf["sample_size_1"]
        exdf_unique = exdf.drop_duplicates(subset="experiment_id")
        all_exp_meta = exdf_unique[
            ["study_key", "experiment_id", "Body site", "Condition", "total_samples"]
        ].copy()
        all_exp_meta.columns = [
            "study_id",
            "experiment_id",
            "body_site",
            "condition",
            "total_samples",
        ]
    else:
        print("experiments.csv not found, metadata-based analysis will be skipped.")
        all_exp_meta = pd.DataFrame()

    # 2. Top 10 Traits (Filtered)
    trait_counts["description"] = trait_counts["trait"].map(
        lambda x: mapping.get(x, mapping.get(x.split("|")[0], x))
    )
    filtered_traits = trait_counts[
        ~trait_counts["description"].str.contains("grows in:", case=False, na=False)
    ].copy()
    top_traits = filtered_traits.sort_values("depletion_count", ascending=False).head(10).copy()
    top_traits["description"] = top_traits["description"].apply(
        lambda x: str(x).split(" [METPO")[0] if " [METPO" in str(x) else x
    )
    top_traits["description"] = top_traits["description"].apply(
        lambda x: str(x).split(" METPO:")[0] if " METPO:" in str(x) else x
    )
    top_traits["description"] = top_traits["description"].apply(shorten_grows_in)

    # 3. Robust Condition-Level Analysis (Panel B)
    if not all_exp_meta.empty:
        meta_inc = all_exp_meta.copy()
        meta_inc["inc_dec"] = "inc"
        meta_dec = all_exp_meta.copy()
        meta_dec["inc_dec"] = "dec"
        full_meta = pd.concat([meta_inc, meta_dec])
        full_meta["profile_id"] = (
            full_meta["study_id"].astype(str)
            + "|"
            + full_meta["experiment_id"].astype(str)
            + "|"
            + full_meta["inc_dec"].astype(str)
        )

        full_profile_data = full_meta.merge(profile_counts, on="profile_id", how="left").fillna(
            {"depletion_count": 0}
        )
        full_profile_data["condition"] = full_profile_data["condition"].fillna("Unknown")
        full_profile_data["body_site"] = full_profile_data["body_site"].fillna("Unknown")
        full_profile_data["total_samples"] = full_profile_data["total_samples"].fillna(0)

        condition_stats = (
            full_profile_data.groupby(["condition", "inc_dec"])["depletion_count"]
            .agg(["sum", "count"])
            .reset_index()
        )
        condition_stats["avg_depletion"] = condition_stats["sum"] / condition_stats["count"]

        cond_pivot = condition_stats.pivot(
            index="condition", columns="inc_dec", values="avg_depletion"
        ).fillna(0)
        cond_pivot["total_avg"] = cond_pivot.sum(axis=1)
        top_conditions_avg = cond_pivot.sort_values("total_avg", ascending=False).head(10)

        # 4. Normalized Body Site Analysis (Panel C)
        valid_site_data = full_profile_data[
            (full_profile_data["total_samples"] > 0)
            & (~full_profile_data["body_site"].str.contains(",", na=False))
        ].copy()

        if not valid_site_data.empty:
            site_stats = (
                valid_site_data.groupby("body_site")
                .agg({"depletion_count": "sum", "total_samples": "sum", "profile_id": "count"})
                .reset_index()
            )
            site_stats["density"] = (
                site_stats["depletion_count"] / site_stats["total_samples"]
            ) * 100
            body_site_density = site_stats.sort_values("density", ascending=False)
        else:
            body_site_density = pd.DataFrame(
                columns=["body_site", "depletion_count", "total_samples", "density"]
            )
    else:
        top_conditions_avg = pd.DataFrame()
        body_site_density = pd.DataFrame()
        full_profile_data = pd.DataFrame()

    # 5. Directionality Bias
    all_profile_activity = full_profile_data[full_profile_data["depletion_count"] > 0]
    direction_activity = all_profile_activity.groupby("inc_dec")["depletion_count"].sum()

    # 6. Write Text File
    print(f"Writing summary to {output_text_path}...")
    with open(output_text_path, "w", encoding="utf-8") as f:
        f.write("=== TOP 10 MOST DEPLETED TRAITS ===\n")
        f.write(f"{'Trait ID':<40} {'Depletion Count':<20} {'Description'}\n")
        f.write("-" * 100 + "\n")
        for _, row in top_traits.iterrows():
            f.write(f"{row['trait']:<40} {row['depletion_count']:<20} {row['description']}\n")

        f.write("\n\n=== TOP 10 CONDITIONS BY AVERAGE DEPLETION DEPTH ===\n")
        f.write(f"{'Condition':<40} {'Avg Total':<15} {'Avg Inc':<15} {'Avg Dec'}\n")
        f.write("-" * 100 + "\n")
        for cond, row in top_conditions_avg.iterrows():
            inc = row.get("inc", 0)
            dec = row.get("dec", 0)
            f.write(f"{cond!s:<40} {row['total_avg']:<15.2f} {inc:<15.2f} {dec:.2f}\n")

        f.write("\n\n=== TOP 10 BODY SITES BY DEPLETION DENSITY (Normalized by Sample Size) ===\n")
        f.write(
            "Note: Sites with missing sample size data are excluded from density normalization.\n"
        )
        f.write(
            f"{'Body Site':<40} {'Calls/100 Reported Samples':<30} "
            f"{'Total Calls':<15} {'Total Samples'}\n"
        )
        f.write("-" * 100 + "\n")
        for _, row in body_site_density.head(10).iterrows():
            bs_str = f"{row['body_site']!s:<40}"
            den_str = f"{row['density']:<30.2f}"
            cnt_str = f"{row['depletion_count']:<15.0f}"
            tot_str = f"{row['total_samples']:.0f}"
            f.write(f"{bs_str} {den_str} {cnt_str} {tot_str}\n")

        f.write("\n\n=== SCIENTIFIC FINDINGS SUMMARY ===\n")
        total_dep = trait_counts["depletion_count"].sum()
        f.write(f"Total Depleted Traits across all profiles: {total_dep:,}\n")
        f.write("Functional scope: 317 metabolic/physiological traits analyzed.\n")
        if not full_profile_data.empty:
            avg_dep = full_profile_data["depletion_count"].mean()
            f.write(f"Average depletion count per profile (all human): {avg_dep:.2f}\n")

        if not body_site_density.empty:
            top_site = body_site_density.iloc[0]["body_site"]
            top_density = body_site_density.iloc[0]["density"]
            f.write(
                f"Top Body Site by Depletion Density: {top_site} "
                f"({top_density:.2f} calls/100 reported samples; "
                "aggregate call ratio, not patient loss)\n"
            )

        bias_inc = int(results_df[results_df["inc_dec"] == "inc"].shape[0])
        bias_dec = int(results_df[results_df["inc_dec"] == "dec"].shape[0])
        if bias_inc == bias_dec:
            f.write(
                f"Directionality Bias: none ({bias_inc:,} inc vs {bias_dec:,} dec selected rows).\n"
            )
        else:
            bias_str = "Increased" if bias_inc > bias_dec else "Decreased"
            f.write(
                f"Directionality Bias: More traits depleted in '{bias_str}' "
                f"profiles ({bias_inc:,} inc vs {bias_dec:,} dec).\n"
            )
            f.write(
                "Directional Distribution Note: Raw inc/dec hit totals reflect compositional "
                "profile counts rather than evidence of taxa dominance or "
                "biological directionality bias.\n"
            )

    # 7. Generate Infographic
    print(f"Generating infographic to {output_img_path}...")
    plt.style.use("seaborn-v0_8-whitegrid")
    fig = plt.figure(figsize=(20, 15))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.35)

    # Panel A: Top 10 Traits (Shortened Labels) - Cool palette
    ax1 = fig.add_subplot(gs[0, 0])
    sns.barplot(
        data=top_traits,
        x="depletion_count",
        y="description",
        hue="description",
        palette="mako",
        ax=ax1,
        legend=False,
    )
    ax1.set_title("A. Top Depleted Traits", fontsize=16, fontweight="bold")
    ax1.set_xlabel("Number of Profiles Depleted In")
    ax1.set_ylabel("")

    # Panel B: Condition Average Depletion (Stacked) - Blues/Teals
    ax2 = fig.add_subplot(gs[0, 1])
    plot_data_b = top_conditions_avg.drop(columns="total_avg").copy()
    plot_data_b.index = [str(i)[:25] for i in plot_data_b.index]
    plot_data_b = plot_data_b.iloc[::-1]
    plot_data_b.plot(kind="barh", stacked=True, color=["#2c3e50", "#2980b9"], ax=ax2, alpha=0.8)

    for i, (_idx, row) in enumerate(plot_data_b.iterrows()):
        total = row.sum()
        if total > 0:
            ax2.text(total + 0.5, i, f"{round(total)}", va="center", fontsize=10, fontweight="bold")

    ax2.set_title("B. Avg Depleted Traits per Condition", fontsize=16, fontweight="bold")
    ax2.set_xlabel("Average Number of Depleted Traits")
    ax2.set_ylabel("")
    ax2.legend(
        ["Decreased", "Increased"],
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        title="Taxonomic Signature",
    )

    # Panel C: Anatomical Representation
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.axis("off")

    if not body_site_density.empty:
        site_coords = {
            "Nasopharynx,Oropharynx": (0.55, 0.96),
            "Mucosa of oral region,Mucosa of oropharynx": (0.57, 0.93),
            "Surface of tongue": (0.56, 0.90),
            "Saliva,Dental plaque": (0.54, 0.87),
            "Subgingival dental plaque,Saliva": (0.56, 0.84),
            "Saliva,Supragingival dental plaque,Buccal mucosa,Tongue": (0.54, 0.81),
            "Dental plaque": (0.54, 0.87),
            "Meconium,Feces": (0.48, 0.35),
            "Meconium": (0.48, 0.35),
            "Mesenteric lymph node": (0.38, 0.42),
            "Large intestine": (0.42, 0.47),
            "Jejunum": (0.44, 0.50),
            "Ileum,Jejunum": (0.46, 0.53),
            "Mucosa of ascending colon": (0.42, 0.44),
            "Caecum": (0.46, 0.42),
            "Rectum": (0.48, 0.30),
            "Skin of abdomen": (0.22, 0.62),
            "Axilla skin": (0.20, 0.65),
            "Skin of body": (0.22, 0.58),
            "Anatomical surface": (0.2, 0.58),
            "Body external integument structure": (0.2, 0.54),
            "Bulbar conjunctiva": (0.22, 0.72),
            "Pinna": (0.22, 0.76),
            "Middle nasal meatus": (0.45, 0.88),
            "Venous blood": (0.32, 0.55),
            "Ovary": (0.52, 0.32),
            "Rumen": (0.42, 0.38),
            "Mucus,Exoskeleton": (0.78, 0.65),
            "Tadpole": (0.22, 0.32),
        }

        t = np.linspace(0, 2 * np.pi, 100)
        ax3.plot(0.5 + 0.07 * np.cos(t), 0.9 + 0.07 * np.sin(t), color="gray", lw=2.0, alpha=0.8)
        ax3.plot(
            [0.38, 0.62, 0.65, 0.6, 0.4, 0.35, 0.38],
            [0.8, 0.8, 0.3, 0.25, 0.25, 0.3, 0.8],
            color="gray",
            lw=2.0,
            alpha=0.8,
        )
        ax3.plot([0.4, 0.37, 0.45, 0.48], [0.25, 0.05, 0.05, 0.25], color="gray", lw=2.0, alpha=0.8)
        ax3.plot([0.6, 0.63, 0.55, 0.52], [0.25, 0.05, 0.05, 0.25], color="gray", lw=2.0, alpha=0.8)
        ax3.plot([0.38, 0.22, 0.25, 0.38], [0.8, 0.6, 0.58, 0.78], color="gray", lw=2.0, alpha=0.8)
        ax3.plot([0.62, 0.78, 0.75, 0.62], [0.8, 0.6, 0.58, 0.78], color="gray", lw=2.0, alpha=0.8)

        plot_sites = body_site_density.head(8).copy()
        max_d = plot_sites["density"].max()
        cmap = plt.get_cmap("GnBu")

        ready_to_plot = []
        for _, row in plot_sites.iterrows():
            site = str(row["body_site"])
            site_lower = site.lower()
            coord = None
            if site in site_coords:
                coord = site_coords[site]
            else:
                for key in sorted(site_coords, key=len, reverse=True):
                    if key.lower() in site_lower or site_lower in key.lower():
                        coord = site_coords[key]
                        break
            if coord:
                ready_to_plot.append(
                    {"name": site, "x": coord[0], "y": coord[1], "density": row["density"]}
                )

        for item in ready_to_plot:
            radius = (item["density"] / max_d) * 0.05 + 0.015
            color = cmap(item["density"] / max_d)
            circle = plt.Circle(
                (item["x"], item["y"]), radius, color=color, alpha=0.6, ec="black", lw=1.5, zorder=3
            )
            ax3.add_patch(circle)
            item["radius"] = radius

        left_side = [i for i in ready_to_plot if i["x"] < 0.5]
        right_side = [i for i in ready_to_plot if i["x"] >= 0.5]
        left_side.sort(key=lambda item: item["y"], reverse=True)
        right_side.sort(key=lambda item: item["y"], reverse=True)

        def place_side_labels(side_items, x_target, ha):
            last_y = 1.02
            min_gap = 0.12
            for item in side_items:
                target_y = min(last_y - min_gap, item["y"])
                target_y = max(0.02, min(0.98, target_y))
                label_text = item["name"].split(",")[0].strip()
                if len(label_text) > 20:
                    label_text = label_text[:18] + ".."
                ax3.annotate(
                    f"{label_text}\n({item['density']:.1f})",
                    xy=(item["x"], item["y"]),
                    xytext=(x_target, target_y),
                    fontsize=9,
                    fontweight="bold",
                    va="center",
                    ha=ha,
                    arrowprops={
                        "arrowstyle": "-",
                        "color": "gray",
                        "lw": 0.8,
                        "alpha": 0.5,
                        "shrinkA": 0,
                        "shrinkB": 5,
                    },
                    bbox={
                        "facecolor": "white",
                        "alpha": 0.9,
                        "edgecolor": "lightgray",
                        "boxstyle": "round,pad=0.3",
                    },
                    zorder=5,
                )
                last_y = target_y

        place_side_labels(left_side, x_target=0.03, ha="left")
        place_side_labels(right_side, x_target=0.97, ha="right")

        ax3.set_xlim(0, 1)
        ax3.set_ylim(0, 1)
        ax3.set_title("C. Anatomical Under-Representation", fontsize=16, fontweight="bold")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=max_d))
        fig.colorbar(sm, ax=ax3, shrink=0.5, label="Calls / 100 Reported Samples")
    else:
        ax3.text(0.5, 0.5, "No Body Site Data Available", ha="center", va="center")

    # Panel D: Scientific Summary
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis("off")

    top_3_sites = body_site_density.head(3)
    top_3_str = ", ".join(
        [f"'{row['body_site']}' ({row['density']:.1f})" for _, row in top_3_sites.iterrows()]
    )

    summary_text = (
        "UNDER-REPRESENTATION SUMMARY\n\n"
        f"• Under-Represented Trait Scope: {len(trait_counts):,} metabolic traits analyzed.\n"
        f"  Total of {trait_counts['depletion_count'].sum():,} "
        "significant under-representations.\n\n"
        f"• Condition Profile Rates: Panel B shows mean under-representation calls per profile,\n"
        "  summarizing trait counts with significantly fewer annotations than "
        "background expectation.\n\n"
        f"• Under-Representation Representation: Top 3 sites by call density\n"
        f"  (calls/100 reported samples): {top_3_str}.\n\n"
        f"• Signature Direction Contrast: {direction_activity.get('inc', 0):,} 'Increased' vs "
        f"{direction_activity.get('dec', 0):,} 'Decreased'\n"
        "  signature calls indicate relative annotation distribution across signature directions."
    )
    ax4.text(
        0.05,
        0.95,
        summary_text,
        transform=ax4.transAxes,
        fontsize=13,
        va="top",
        ha="left",
        bbox={"boxstyle": "round,pad=1", "facecolor": "ghostwhite", "alpha": 0.7},
    )
    ax4.set_title("D. Scientific Implications", fontsize=16, fontweight="bold")

    plt.suptitle(
        "CultureBotAI: Meta-Analysis of Microbial Trait Depletion",
        fontsize=24,
        fontweight="bold",
        y=0.98,
    )
    plt.savefig(output_img_path, dpi=300, bbox_inches="tight")
    print(f"Infographic saved to {output_img_path}")
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
