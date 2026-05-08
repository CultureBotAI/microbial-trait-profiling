"""Create a clustered heatmap from the BugSigDB YAML-derived taxa TSV."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "trait-visualization/bsdb-taxa-only"


def build_binary_matrix(
    taxa_path: Path,
    top_taxa: int,
    max_profiles: int,
) -> tuple[pd.DataFrame, pd.Series]:
    """Build a binary profile x genus matrix from taxa.tsv."""
    taxa = pd.read_csv(
        taxa_path,
        sep="\t",
        usecols=["study_id", "experiment_id", "direction", "taxon_name", "taxonomic_rank"],
    )
    taxa = taxa[taxa["taxonomic_rank"] == "genus"].copy()
    taxa["inc_dec"] = (
        taxa["direction"].map({"increased": "inc", "decreased": "dec"}).fillna(taxa["direction"])
    )
    taxa["profile_id"] = taxa["study_id"] + "|" + taxa["experiment_id"] + "|" + taxa["inc_dec"]

    # Binary presence/absence: one genus counts once per profile.
    taxa = taxa[["profile_id", "inc_dec", "taxon_name"]].drop_duplicates()

    top_taxon_names = (
        taxa.groupby("taxon_name")["profile_id"]
        .nunique()
        .sort_values(ascending=False)
        .head(top_taxa)
        .index
    )
    taxa = taxa[taxa["taxon_name"].isin(top_taxon_names)].copy()

    matrix = (
        taxa.assign(value=1)
        .pivot(index="profile_id", columns="taxon_name", values="value")
        .fillna(0)
        .astype(int)
    )
    matrix = matrix.loc[:, matrix.sum(axis=0).sort_values(ascending=False).index]

    if len(matrix) > max_profiles:
        profile_richness = matrix.sum(axis=1)
        selected_profiles = (
            profile_richness.sort_values(ascending=False).head(max_profiles).sort_index().index
        )
        matrix = matrix.loc[selected_profiles]

    inc_dec = pd.Series(
        [profile_id.rsplit("|", 1)[-1] for profile_id in matrix.index],
        index=matrix.index,
        name="inc_dec",
    )
    return matrix, inc_dec


def save_clustered_heatmap(
    matrix: pd.DataFrame,
    inc_dec: pd.Series,
    output_path: Path,
) -> None:
    """Render and save a clustered heatmap."""
    sns.set_theme(style="white")
    row_colors = inc_dec.map({"inc": "#c0392b", "dec": "#2874a6"}).fillna("#7f8c8d")

    cluster = sns.clustermap(
        matrix,
        cmap=sns.color_palette(["#f7f7f7", "#0f766e"], as_cmap=True),
        method="average",
        metric="jaccard",
        row_colors=row_colors,
        linewidths=0,
        xticklabels=True,
        yticklabels=False,
        figsize=(24, 18),
        cbar_pos=(0.02, 0.82, 0.02, 0.12),
        dendrogram_ratio=(0.14, 0.12),
        colors_ratio=(0.03, 0.03),
    )

    cluster.ax_heatmap.set_xlabel("Genus")
    cluster.ax_heatmap.set_ylabel("Profiles")
    cluster.ax_heatmap.set_title(
        "BugSigDB Genus Presence/Absence Clustered Heatmap\n"
        "Top prevalent genera; rows are study|experiment|direction profiles",
        pad=18,
    )

    for label in cluster.ax_heatmap.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(8)

    cluster.fig.text(0.02, 0.96, "Row annotation", fontsize=10, fontweight="bold")
    cluster.fig.text(0.02, 0.94, "red = increased", color="#c0392b", fontsize=9)
    cluster.fig.text(0.02, 0.92, "blue = decreased", color="#2874a6", fontsize=9)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cluster.fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(cluster.fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Create a clustered heatmap from taxa.tsv using a binary profile x genus matrix. "
            "Profiles are filtered to the most information-dense rows for readability."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=REPO_ROOT / "data/derived/taxa.tsv",
        help="Input taxa TSV path.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_DIR / "clustered_heatmap_genus_top100.png",
        help="Output heatmap image path.",
    )
    parser.add_argument(
        "--matrix-output",
        type=Path,
        default=OUTPUT_DIR / "clustered_heatmap_genus_top100_matrix.tsv",
        help="Output path for the filtered binary matrix used to plot the heatmap.",
    )
    parser.add_argument(
        "--top-taxa",
        type=int,
        default=100,
        help="Number of most prevalent genera to keep.",
    )
    parser.add_argument(
        "--max-profiles",
        type=int,
        default=300,
        help="Maximum number of profiles to include in the rendered heatmap.",
    )
    args = parser.parse_args()

    matrix, inc_dec = build_binary_matrix(args.input, args.top_taxa, args.max_profiles)

    matrix_to_save = matrix.copy()
    matrix_to_save.insert(0, "inc_dec", inc_dec)
    matrix_to_save.to_csv(args.matrix_output, sep="\t")

    save_clustered_heatmap(matrix, inc_dec, args.output)

    print(f"Wrote matrix: {args.matrix_output}")
    print(f"Wrote heatmap: {args.output}")
    print(f"Profiles plotted: {len(matrix)}")
    print(f"Genera plotted: {len(matrix.columns)}")


if __name__ == "__main__":
    main()
