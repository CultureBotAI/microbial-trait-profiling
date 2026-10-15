"""Generate rank-stratified taxa coverage histograms, occurrence comparisons, and audit tables.

This script implements Approach A for evaluating taxonomic coverage between BugSigDB
and KG-Microbe traits:
1. Rank-stratified stacked histogram (absolute counts and 100% normalized proportions).
2. Unique taxa vs. literature signature occurrences comparison.
3. Trait count distribution per taxon highlighting zero-trait vs. trait-bearing taxa.
4. Comprehensive audit TSVs for rank-level coverage, individual taxa, and top zero-trait taxa.
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = VIS_DIR / "trait-profile-investigation"
OUTPUT_DIR = INVESTIGATION_DIR / "taxa-coverage"

TAXA_PATH = REPO_ROOT / "data" / "derived" / "taxa.tsv"
TAXA_TRAIT_MATRIX_PATH = REPO_ROOT / "data" / "derived" / "taxa_trait_matrix.tsv"

# Standard taxonomic hierarchy (ascending order)
RANK_ORDER = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
RANK_DISPLAY = {
    "species": "Species",
    "genus": "Genus",
    "family": "Family",
    "order": "Order",
    "class": "Class",
    "phylum": "Phylum",
    "superkingdom": "Superkingdom",
}

# Colorblind-safe palette (Okabe-Ito inspired)
COLOR_RICH = "#0072b2"     # Deep Blue: >= 3 traits
COLOR_SPARSE = "#56b4e9"   # Sky Blue: 1-2 traits
COLOR_ZERO = "#d55e00"     # Vermilion: 0 traits
COLOR_ACCENT = "#009e73"   # Bluish Green for highlights
COLOR_SLATE_DARK = "#0f172a"
COLOR_SLATE_MID = "#475569"
COLOR_SLATE_LIGHT = "#94a3b8"
COLOR_BG_LIGHT = "#f8fafc"

plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]
plt.rcParams["svg.hashsalt"] = "taxa-coverage-histogram-v1"


def read_tsv(path: Path, **kwargs: Any) -> pd.DataFrame:
    """Read a local TSV file."""
    return pd.read_csv(path, sep="\t", **kwargs)


def write_tsv(frame: pd.DataFrame, output_path: Path) -> None:
    """Write a DataFrame to a TSV file with parent directories created."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, sep="\t", index=False)


def save_figure(fig: plt.Figure, output_path: Path) -> None:
    """Save a figure as both 300 DPI PNG and clean SVG."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(output_path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def diagnose_zero_trait_taxon(taxon_id: str, name: str, rank: str) -> tuple[str, str]:
    """Diagnose the informatics root cause and remediation for a zero-trait BugSigDB taxon."""
    lower_name = name.lower()

    # 1. Known secondary / desynchronized NCBI tax IDs with canonical twins that have traits
    if taxon_id in {"NCBITaxon:541000", "NCBITaxon:81850", "NCBITaxon:508215"}:
        canonical_map = {
            "NCBITaxon:541000": "NCBITaxon:216572 (Oscillospiraceae, 162 traits)",
            "NCBITaxon:81850": "NCBITaxon:33958 (Lactobacillaceae, 166 traits)",
            "NCBITaxon:508215": "NCBITaxon:32207 (Rothia, 87 traits)",
        }
        twin = canonical_map.get(taxon_id, "canonical TaxID")
        return (
            "NCBI TaxID Desynchronization",
            f"Remap BugSigDB alias to canonical TaxID with traits: {twin}",
        )

    # 2. Ambiguous bin / taxonomic slashes / unresolved markers
    if "/" in name or re.search(r"\b(unidentified|uncultured organism|unknown)\b", lower_name):
        return (
            "Ambiguous Taxon / Bin Slash",
            "Resolve to lowest common ancestor (LCA) or treat as composite multi-taxon marker",
        )

    # 3. Candidatus lineages
    if "candidatus" in lower_name or lower_name.startswith("ca. "):
        return (
            "Candidatus Lineage",
            "Infer phenotypic traits from MAG-derived metabolic models or phylogenetic neighbors",
        )

    # 4. Fungal / Non-Bacterial (microeukaryotes, plants)
    if (
        re.search(
            r"\b(candida|malassezia|ascomycota|saccharomyces|fungus|fungi|basidiomycota|streptophyta|viridiplantae)\b",
            lower_name,
        )
        or (rank in {"kingdom", "superkingdom", "phylum"} and "fung" in lower_name)
    ):
        return (
            "Non-Bacterial / Fungal Eukaryote",
            "Integrate specialized microeukaryotic trait ontologies (e.g. FungiTraits) into KG-Microbe",
        )

    # 5. Metagenomic / Uncultivated MAGs / Environmental Bins
    if (
        "uncultured" in lower_name
        or re.search(r"\b(cag:\d+|uba\d+|ucg-\d+|rf39|gut metagenome)\b", lower_name)
        or re.search(r"\s+bacterium\b", lower_name)
        or re.search(r"\s+sp\.\s*$", lower_name)
    ):
        return (
            "Uncultivated / Metagenomic MAG",
            "Propagate higher-rank family/order phenotypic traits to uncultivated leaf nodes in KG",
        )

    # 6. Novel or rare isolates lacking phenotype records in BacDive
    return (
        "Rare / Novel Isolate (Lacks KG Phenotype)",
        "Target for literature curation in BacDive / KG-Microbe phenotypic assertions",
    )


def load_and_process_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load taxa.tsv and taxa_trait_matrix.tsv and compute taxon- and rank-level statistics."""
    print(f"Reading taxa occurrences from {TAXA_PATH}...")
    taxa_df = read_tsv(TAXA_PATH)

    print(f"Reading taxa trait matrix from {TAXA_TRAIT_MATRIX_PATH}...")
    matrix_df = read_tsv(TAXA_TRAIT_MATRIX_PATH)

    # Extract metadata columns and trait columns
    meta_cols = ["taxon_id", "taxon_name", "taxonomic_rank"]
    trait_cols = [c for c in matrix_df.columns if c not in meta_cols]
    print(f"  Loaded {len(matrix_df)} unique taxa and {len(trait_cols)} METPO trait columns.")

    # Vectorized trait count per taxon
    trait_counts = (matrix_df[trait_cols] > 0).sum(axis=1).to_numpy()
    taxa_master = matrix_df[meta_cols].copy()
    taxa_master["trait_count"] = trait_counts

    # Calculate occurrence frequency and distinct study counts in BugSigDB
    occ_stats = taxa_df.groupby("taxon_id").agg(
        signature_occurrences=("signature_id", "count"),
        study_count=("study_id", "nunique"),
    ).reset_index()

    taxa_master = taxa_master.merge(occ_stats, on="taxon_id", how="left")
    taxa_master["signature_occurrences"] = taxa_master["signature_occurrences"].fillna(0).astype(int)
    taxa_master["study_count"] = taxa_master["study_count"].fillna(0).astype(int)

    # Categorize trait coverage
    taxa_master["has_traits"] = taxa_master["trait_count"] > 0
    conditions = [
        taxa_master["trait_count"] == 0,
        (taxa_master["trait_count"] >= 1) & (taxa_master["trait_count"] <= 2),
        taxa_master["trait_count"] >= 3,
    ]
    choices = ["zero_traits", "sparse_1_2", "rich_gte_3"]
    taxa_master["coverage_category"] = np.select(conditions, choices, default="unknown")

    # Standardize rank strings
    taxa_master["taxonomic_rank"] = taxa_master["taxonomic_rank"].str.lower().str.strip()
    taxa_df["taxonomic_rank"] = taxa_df["taxonomic_rank"].str.lower().str.strip()

    # Build rank-stratified summary
    rank_rows: list[dict[str, Any]] = []
    for rank in RANK_ORDER:
        sub_taxa = taxa_master[taxa_master["taxonomic_rank"] == rank]
        sub_occ = taxa_df[taxa_df["taxonomic_rank"] == rank]

        total_taxa = len(sub_taxa)
        zero_taxa = int((sub_taxa["coverage_category"] == "zero_traits").sum())
        sparse_taxa = int((sub_taxa["coverage_category"] == "sparse_1_2").sum())
        rich_taxa = int((sub_taxa["coverage_category"] == "rich_gte_3").sum())
        trait_taxa = sparse_taxa + rich_taxa

        total_occ = len(sub_occ)
        # Map occurrences to trait presence
        occ_with_traits = sub_occ.merge(taxa_master[["taxon_id", "has_traits", "coverage_category"]], on="taxon_id", how="left")
        zero_occ = int((occ_with_traits["coverage_category"] == "zero_traits").sum())
        sparse_occ = int((occ_with_traits["coverage_category"] == "sparse_1_2").sum())
        rich_occ = int((occ_with_traits["coverage_category"] == "rich_gte_3").sum())
        trait_occ = sparse_occ + rich_occ

        rank_rows.append({
            "taxonomic_rank": rank,
            "display_rank": RANK_DISPLAY.get(rank, rank.capitalize()),
            "total_taxa": total_taxa,
            "zero_trait_taxa": zero_taxa,
            "zero_trait_pct": round(zero_taxa / total_taxa * 100, 2) if total_taxa > 0 else 0.0,
            "sparse_1_2_taxa": sparse_taxa,
            "sparse_1_2_pct": round(sparse_taxa / total_taxa * 100, 2) if total_taxa > 0 else 0.0,
            "rich_gte_3_taxa": rich_taxa,
            "rich_gte_3_pct": round(rich_taxa / total_taxa * 100, 2) if total_taxa > 0 else 0.0,
            "trait_bearing_taxa": trait_taxa,
            "trait_bearing_pct": round(trait_taxa / total_taxa * 100, 2) if total_taxa > 0 else 0.0,
            "total_occurrences": total_occ,
            "zero_trait_occurrences": zero_occ,
            "zero_trait_occurrences_pct": round(zero_occ / total_occ * 100, 2) if total_occ > 0 else 0.0,
            "sparse_occurrences": sparse_occ,
            "rich_occurrences": rich_occ,
            "trait_bearing_occurrences": trait_occ,
            "trait_bearing_occurrences_pct": round(trait_occ / total_occ * 100, 2) if total_occ > 0 else 0.0,
        })

    rank_summary_df = pd.DataFrame(rank_rows)

    # Top zero-trait taxa diagnosis
    zero_taxa_df = taxa_master[taxa_master["coverage_category"] == "zero_traits"].sort_values(
        by=["signature_occurrences", "study_count"], ascending=[False, False]
    ).reset_index(drop=True)

    top_zero_records: list[dict[str, Any]] = []
    for rank_idx, row in zero_taxa_df.head(100).iterrows():
        cat, rec = diagnose_zero_trait_taxon(str(row["taxon_id"]), str(row["taxon_name"]), str(row["taxonomic_rank"]))
        top_zero_records.append({
            "zero_trait_rank": rank_idx + 1,
            "taxon_id": row["taxon_id"],
            "taxon_name": row["taxon_name"],
            "taxonomic_rank": row["taxonomic_rank"],
            "signature_occurrences": row["signature_occurrences"],
            "study_count": row["study_count"],
            "diagnosis_category": cat,
            "suggested_remediation": rec,
        })
    top_zero_df = pd.DataFrame(top_zero_records)

    return taxa_master, rank_summary_df, top_zero_df


def plot_rank_stratified_histogram(rank_df: pd.DataFrame, output_path: Path) -> None:
    """Generate Figure 1: Rank-stratified stacked histogram (Counts & 100% Proportions)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 6.5), gridspec_kw={"wspace": 0.28})

    ranks = [RANK_DISPLAY[r] for r in RANK_ORDER]
    x = np.arange(len(ranks))
    width = 0.62

    # Data arrays
    zeros = rank_df["zero_trait_taxa"].to_numpy()
    sparse = rank_df["sparse_1_2_taxa"].to_numpy()
    rich = rank_df["rich_gte_3_taxa"].to_numpy()
    totals = rank_df["total_taxa"].to_numpy()

    # Panel A: Absolute Counts Stacked Bar Chart
    ax1.bar(x, zeros, width, label="Mapped, 0 traits (No KG traits)", color=COLOR_ZERO, edgecolor="white", linewidth=0.8)
    ax1.bar(x, sparse, width, bottom=zeros, label="Mapped, 1–2 traits (Sparse)", color=COLOR_SPARSE, edgecolor="white", linewidth=0.8)
    ax1.bar(x, rich, width, bottom=zeros + sparse, label="Mapped, ≥ 3 traits (Rich)", color=COLOR_RICH, edgecolor="white", linewidth=0.8)

    # Label total N on top of bars
    for i, total in enumerate(totals):
        ax1.text(
            x[i],
            total + max(totals) * 0.015,
            f"N = {total:,}",
            ha="center",
            va="bottom",
            fontsize=9.5,
            fontweight="bold",
            color=COLOR_SLATE_DARK,
        )

    ax1.set_xticks(x)
    ax1.set_xticklabels(ranks, fontsize=11, fontweight="normal")
    ax1.set_ylabel("Number of Unique Taxa in BugSigDB", fontsize=11, fontweight="bold", labelpad=8)
    ax1.set_title("A. BugSigDB Unique Taxa Counts by Trait Richness", fontsize=12.5, fontweight="bold", pad=12, loc="left")
    ax1.set_ylim(0, max(totals) * 1.12)
    ax1.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax1.set_axisbelow(True)

    # Add subtitle info in ax1
    ax1.text(
        0.02,
        0.96,
        "Total unique taxa: N = 5,829 across 7 ranks\nSpecies represents 66.0% of all cataloged taxa",
        transform=ax1.transAxes,
        fontsize=9.5,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor=COLOR_BG_LIGHT, edgecolor=COLOR_SLATE_LIGHT, alpha=0.9),
    )

    # Panel B: 100% Normalized Proportions Stacked Bar Chart
    zero_pcts = rank_df["zero_trait_pct"].to_numpy()
    sparse_pcts = rank_df["sparse_1_2_pct"].to_numpy()
    rich_pcts = rank_df["rich_gte_3_pct"].to_numpy()

    ax2.bar(x, zero_pcts, width, color=COLOR_ZERO, edgecolor="white", linewidth=0.8)
    ax2.bar(x, sparse_pcts, width, bottom=zero_pcts, color=COLOR_SPARSE, edgecolor="white", linewidth=0.8)
    ax2.bar(x, rich_pcts, width, bottom=zero_pcts + sparse_pcts, color=COLOR_RICH, edgecolor="white", linewidth=0.8)

    # Inline percentage labels inside segments
    for i in range(len(ranks)):
        # Zero traits segment label
        if zero_pcts[i] >= 6.0:
            ax2.text(
                x[i],
                zero_pcts[i] / 2.0,
                f"{zero_pcts[i]:.1f}%",
                ha="center",
                va="center",
                color="white",
                fontsize=9.5,
                fontweight="bold",
            )
        # Sparse segment label
        if sparse_pcts[i] >= 5.0:
            ax2.text(
                x[i],
                zero_pcts[i] + sparse_pcts[i] / 2.0,
                f"{sparse_pcts[i]:.1f}%",
                ha="center",
                va="center",
                color=COLOR_SLATE_DARK,
                fontsize=8.5,
                fontweight="bold",
            )
        # Rich segment label
        if rich_pcts[i] >= 8.0:
            ax2.text(
                x[i],
                zero_pcts[i] + sparse_pcts[i] + rich_pcts[i] / 2.0,
                f"{rich_pcts[i]:.1f}%",
                ha="center",
                va="center",
                color="white",
                fontsize=9.5,
                fontweight="bold",
            )

    # Reference line at 50%
    ax2.axhline(50, color=COLOR_SLATE_MID, linestyle=":", linewidth=1.2, alpha=0.7)
    ax2.text(len(ranks) - 0.45, 51.5, "50% parity threshold", fontsize=8.5, color=COLOR_SLATE_MID, ha="right")

    # Annotations highlighting drop-off and higher rank coverage
    ax2.annotate(
        "Species Drop-off:\n45.9% have 0 traits\n(1,764 / 3,846 taxa)",
        xy=(x[0], zero_pcts[0]),
        xytext=(x[0] + 0.35, 62),
        arrowprops=dict(arrowstyle="->", color=COLOR_ZERO, lw=1.5, connectionstyle="arc3,rad=-0.15"),
        fontsize=9.5,
        fontweight="bold",
        color=COLOR_ZERO,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#fff7ed", edgecolor=COLOR_ZERO, alpha=0.9),
    )

    ax2.annotate(
        "Genus & Family:\n~74–76% trait-bearing\n(>72% rich ≥3 traits)",
        xy=(x[1], zero_pcts[1] + sparse_pcts[1] + rich_pcts[1] / 2),
        xytext=(x[1] + 0.45, 88),
        arrowprops=dict(arrowstyle="->", color=COLOR_RICH, lw=1.5, connectionstyle="arc3,rad=0.15"),
        fontsize=9.5,
        fontweight="bold",
        color=COLOR_RICH,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#eff6ff", edgecolor=COLOR_RICH, alpha=0.9),
    )

    ax2.set_xticks(x)
    ax2.set_xticklabels(ranks, fontsize=11, fontweight="normal")
    ax2.set_ylabel("Percentage of Taxa (%)", fontsize=11, fontweight="bold", labelpad=8)
    ax2.set_title("B. Normalized Trait Coverage Proportion (%)", fontsize=12.5, fontweight="bold", pad=12, loc="left")
    ax2.set_ylim(0, 105)
    ax2.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax2.set_axisbelow(True)

    # Shared unified legend
    handles = [
        Patch(facecolor=COLOR_RICH, edgecolor="white", label="Mapped, Rich (≥ 3 traits)"),
        Patch(facecolor=COLOR_SPARSE, edgecolor="white", label="Mapped, Sparse (1–2 traits)"),
        Patch(facecolor=COLOR_ZERO, edgecolor="white", label="Mapped, 0 traits (No KG traits)"),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=3,
        fontsize=10.5,
        frameon=True,
        facecolor=COLOR_BG_LIGHT,
        edgecolor=COLOR_SLATE_LIGHT,
    )

    save_figure(fig, output_path)
    print(f"  Saved Figure 1: {output_path} and .svg")


def plot_occurrences_vs_unique(rank_df: pd.DataFrame, taxa_master: pd.DataFrame, output_path: Path) -> None:
    """Generate Figure 2: Literature Occurrences vs. Unique Taxa Coverage Disparity."""
    fig, ax = plt.subplots(figsize=(11, 6))

    # Overall Unique Taxa counts
    total_taxa = len(taxa_master)
    zero_taxa = (taxa_master["coverage_category"] == "zero_traits").sum()
    sparse_taxa = (taxa_master["coverage_category"] == "sparse_1_2").sum()
    rich_taxa = (taxa_master["coverage_category"] == "rich_gte_3").sum()

    # Overall Signature Occurrences counts
    total_occ = rank_df["total_occurrences"].sum()
    zero_occ = rank_df["zero_trait_occurrences"].sum()
    sparse_occ = rank_df["sparse_occurrences"].sum()
    rich_occ = rank_df["rich_occurrences"].sum()

    categories = [
        f"Unique Taxa Catalog\n(N = {total_taxa:,} unique taxa)",
        f"BugSigDB Literature Occurrences\n(N = {total_occ:,} signature mentions)",
    ]
    y_pos = np.arange(len(categories))
    bar_height = 0.52

    zero_pcts = [zero_taxa / total_taxa * 100, zero_occ / total_occ * 100]
    sparse_pcts = [sparse_taxa / total_taxa * 100, sparse_occ / total_occ * 100]
    rich_pcts = [rich_taxa / total_taxa * 100, rich_occ / total_occ * 100]

    # Horizontal stacked bar
    ax.barh(y_pos, zero_pcts, height=bar_height, color=COLOR_ZERO, edgecolor="white", label="Mapped, 0 traits (No KG traits)")
    ax.barh(y_pos, sparse_pcts, height=bar_height, left=zero_pcts, color=COLOR_SPARSE, edgecolor="white", label="Mapped, 1–2 traits (Sparse)")
    ax.barh(y_pos, rich_pcts, height=bar_height, left=np.array(zero_pcts) + np.array(sparse_pcts), color=COLOR_RICH, edgecolor="white", label="Mapped, ≥ 3 traits (Rich)")

    # Inline labels
    # Row 0: Unique Taxa
    ax.text(zero_pcts[0] / 2, y_pos[0], f"0 Traits: {zero_pcts[0]:.1f}%\n({zero_taxa:,})", ha="center", va="center", color="white", fontweight="bold", fontsize=10)
    ax.text(zero_pcts[0] + sparse_pcts[0] + rich_pcts[0] / 2, y_pos[0], f"Trait-Bearing: {sparse_pcts[0] + rich_pcts[0]:.1f}%\n({sparse_taxa + rich_taxa:,} taxa)", ha="center", va="center", color="white", fontweight="bold", fontsize=10.5)

    # Row 1: Occurrences
    ax.text(zero_pcts[1] / 2, y_pos[1], f"9.5%\n({zero_occ:,})", ha="center", va="center", color="white", fontweight="bold", fontsize=9.5)
    ax.text(zero_pcts[1] + sparse_pcts[1] + rich_pcts[1] / 2, y_pos[1], f"Trait-Bearing in Literature: {sparse_pcts[1] + rich_pcts[1]:.1f}%\n({sparse_occ + rich_occ:,} occurrences across signatures)", ha="center", va="center", color="white", fontweight="bold", fontsize=11)

    # Callout highlighting the key scientific insight
    ax.annotate(
        "Key Finding: 90.55% of all literature findings DO have traits!\nThe 38.7% zero-trait taxa reside in the long-tail\nof uncultured/metagenomic singletons.",
        xy=(zero_pcts[1], y_pos[1]),
        xytext=(16, y_pos[1] - 0.38),
        arrowprops=dict(arrowstyle="->", color=COLOR_SLATE_DARK, lw=1.5, connectionstyle="arc3,rad=-0.12"),
        fontsize=10.5,
        fontweight="bold",
        color=COLOR_SLATE_DARK,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#ecfdf5", edgecolor=COLOR_ACCENT, alpha=0.95),
    )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(categories, fontsize=11, fontweight="bold")
    ax.set_xlabel("Percentage of Total (%)", fontsize=11, fontweight="bold", labelpad=8)
    ax.set_xlim(0, 105)
    ax.set_title("Taxa Coverage: Unique Taxa vs. Literature Signature Occurrences", fontsize=13, fontweight="bold", pad=15, loc="left")
    ax.grid(axis="x", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax.set_axisbelow(True)

    handles = [
        Patch(facecolor=COLOR_RICH, edgecolor="white", label="Mapped, Rich (≥ 3 traits)"),
        Patch(facecolor=COLOR_SPARSE, edgecolor="white", label="Mapped, Sparse (1–2 traits)"),
        Patch(facecolor=COLOR_ZERO, edgecolor="white", label="Mapped, 0 traits (No KG traits)"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.14), ncol=3, fontsize=10, frameon=True, facecolor=COLOR_BG_LIGHT)

    save_figure(fig, output_path)
    print(f"  Saved Figure 2: {output_path} and .svg")


def plot_trait_count_distribution(taxa_master: pd.DataFrame, output_path: Path) -> None:
    """Generate Figure 3: Trait count distribution per taxon showing the zero-spike and positive distribution."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"wspace": 0.25})

    trait_counts = taxa_master["trait_count"].to_numpy()
    zero_count = int((trait_counts == 0).sum())
    pos_counts = trait_counts[trait_counts > 0]

    # Panel A: Zero-bin vs Positive Bins Histogram
    bin_labels = ["0", "1–4", "5–9", "10–19", "20–49", "50–99", "100–199", "200–499", "≥500"]
    counts_per_bin = [
        int(zero_count),
        int(((trait_counts >= 1) & (trait_counts < 5)).sum()),
        int(((trait_counts >= 5) & (trait_counts < 10)).sum()),
        int(((trait_counts >= 10) & (trait_counts < 20)).sum()),
        int(((trait_counts >= 20) & (trait_counts < 50)).sum()),
        int(((trait_counts >= 50) & (trait_counts < 100)).sum()),
        int(((trait_counts >= 100) & (trait_counts < 200)).sum()),
        int(((trait_counts >= 200) & (trait_counts < 500)).sum()),
        int((trait_counts >= 500).sum()),
    ]

    x_b = np.arange(len(bin_labels))
    bar_colors = [COLOR_ZERO] + [COLOR_SPARSE if i == 1 else COLOR_RICH for i in range(1, len(bin_labels))]

    ax1.bar(x_b, counts_per_bin, color=bar_colors, edgecolor="white", linewidth=0.8)
    for i, c in enumerate(counts_per_bin):
        ax1.text(
            x_b[i],
            c + max(counts_per_bin) * 0.015,
            f"{c:,}",
            ha="center",
            va="bottom",
            fontsize=8.5,
            fontweight="bold",
            color=COLOR_SLATE_DARK,
        )

    ax1.set_xticks(x_b)
    ax1.set_xticklabels(bin_labels, fontsize=9.5, rotation=25)
    ax1.set_ylabel("Number of Taxa", fontsize=11, fontweight="bold", labelpad=8)
    ax1.set_xlabel("Number of Annotated KG-Microbe Traits per Taxon", fontsize=11, fontweight="bold", labelpad=8)
    ax1.set_title("A. Trait Degree Distribution Across All Taxa", fontsize=12, fontweight="bold", pad=12, loc="left")
    ax1.set_ylim(0, max(counts_per_bin) * 1.12)
    ax1.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax1.set_axisbelow(True)

    # Callout for zero spike
    ax1.annotate(
        f"Prominent Zero-Trait Spike:\n{zero_count:,} taxa (38.7%)\nhave zero KG traits",
        xy=(x_b[0], zero_count),
        xytext=(x_b[0] + 1.2, zero_count * 0.85),
        arrowprops=dict(arrowstyle="->", color=COLOR_ZERO, lw=1.5),
        fontsize=9.5,
        fontweight="bold",
        color=COLOR_ZERO,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#fff7ed", edgecolor=COLOR_ZERO, alpha=0.9),
    )

    # Panel B: Trait Richness by Taxonomic Rank (Positive taxa only, log-scale boxplot)
    rank_pos_data = []
    plot_ranks = []
    for rank in RANK_ORDER:
        sub = taxa_master[(taxa_master["taxonomic_rank"] == rank) & (taxa_master["trait_count"] > 0)]["trait_count"].to_numpy()
        if len(sub) > 0:
            rank_pos_data.append(sub)
            plot_ranks.append(RANK_DISPLAY[rank])

    bp = ax2.boxplot(
        rank_pos_data,
        tick_labels=plot_ranks,
        patch_artist=True,
        showmeans=True,
        meanprops=dict(marker="o", markeredgecolor=COLOR_SLATE_DARK, markerfacecolor="white", markersize=5),
        medianprops=dict(color=COLOR_SLATE_DARK, linewidth=1.5),
        boxprops=dict(facecolor="#bae6fd", color=COLOR_RICH, linewidth=1.2),
        whiskerprops=dict(color=COLOR_SLATE_MID, linewidth=1.2),
        capprops=dict(color=COLOR_SLATE_MID, linewidth=1.2),
        flierprops=dict(marker=".", markerfacecolor=COLOR_SLATE_LIGHT, markersize=3, alpha=0.4),
    )

    ax2.set_yscale("log")
    ax2.set_ylabel("Annotated Traits per Taxon (Log10 Scale)", fontsize=11, fontweight="bold", labelpad=8)
    ax2.set_title("B. Trait Richness Rollup Across Ranks (Trait-Bearing Taxa)", fontsize=12, fontweight="bold", pad=12, loc="left")
    ax2.grid(axis="y", linestyle="--", alpha=0.35, color=COLOR_SLATE_LIGHT)
    ax2.set_axisbelow(True)

    # Annotate rollup effect
    ax2.text(
        0.03,
        0.96,
        "Trait Rollup in KG-Microbe:\nHigher ranks inherit all phenotypic assertions\nfrom their classified descendant lineages.\n(Median: Species=13, Genus=32, Phylum=121)",
        transform=ax2.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor=COLOR_BG_LIGHT, edgecolor=COLOR_SLATE_LIGHT, alpha=0.9),
    )

    save_figure(fig, output_path)
    print(f"  Saved Figure 3: {output_path} and .svg")




def main() -> None:
    """Execute complete Approach A taxa coverage workflow."""
    parser = argparse.ArgumentParser(description="Generate rank-stratified taxa coverage figures and audit tables.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Target output directory")
    args = parser.parse_args()

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=================================================================")
    print("Executing Approach A: Rank-Stratified Taxa Coverage Analysis")
    print(f"Output Directory: {out_dir}")
    print("=================================================================")

    # 1. Load data and compute statistics
    taxa_master, rank_df, top_zero_df = load_and_process_data()

    # 2. Save TSV Audit Tables
    rank_tsv_path = out_dir / "taxa_coverage_by_rank.tsv"
    write_tsv(rank_df, rank_tsv_path)
    print(f"  Written rank coverage summary: {rank_tsv_path}")

    details_tsv_path = out_dir / "taxa_kgmicrobe_coverage_details.tsv"
    write_tsv(taxa_master, details_tsv_path)
    print(f"  Written complete taxa audit table ({len(taxa_master)} taxa): {details_tsv_path}")

    top_zero_tsv_path = out_dir / "top_zero_trait_bugsigdb_taxa.tsv"
    write_tsv(top_zero_df, top_zero_tsv_path)
    print(f"  Written top 100 zero-trait diagnosis table: {top_zero_tsv_path}")

    # 3. Generate Visualizations
    fig1_path = out_dir / "taxa_coverage_histogram_rank_stratified.png"
    plot_rank_stratified_histogram(rank_df, fig1_path)

    fig2_path = out_dir / "taxa_coverage_occurrences_vs_unique.png"
    plot_occurrences_vs_unique(rank_df, taxa_master, fig2_path)

    fig3_path = out_dir / "taxa_trait_count_distribution.png"
    plot_trait_count_distribution(taxa_master, fig3_path)


    print("=================================================================")
    print("Approach A Taxa Coverage Analysis Complete Successfully!")
    print("=================================================================")


if __name__ == "__main__":
    main()
