"""Create a clustered heatmap from the BugSigDB YAML-derived taxa TSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "trait-visualization/bsdb-taxa-only"
HUMAN_HOST_SPECIES = "Homo sapiens"


def clean_value(value: object) -> str:
    """Normalize empty BugSigDB metadata values."""
    text = str(value or "").strip()
    if not text or text.upper() == "NA" or text.lower() == "nan":
        return ""
    return text


def study_id_from_bsdb_id(bsdb_id: str) -> str:
    """Convert a full BugSigDB signature ID into the derived study ID format."""
    compact_id = bsdb_id.strip().removeprefix("bsdb:").split("/", 1)[0]
    return f"bsdb:{compact_id}"


def experiment_id_from_bsdb_id(bsdb_id: str) -> str:
    """Convert a full BugSigDB signature ID into the derived experiment ID format."""
    parts = bsdb_id.strip().removeprefix("bsdb:").split("/")
    if len(parts) >= 2:
        return f"bsdb:{parts[0]}-{parts[1]}"
    return bsdb_id


def load_human_experiment_ids(full_dump_path: Path) -> set[tuple[str, str]]:
    """Return study/experiment IDs for Homo sapiens 16S and WMS rows in full_dump.csv.

    The raw export is required. The sequencing type keeps the analysis corpus's
    rule; a blank one is left out (#742, #787).
    """
    if full_dump_path is None:
        raise ValueError("full_dump.csv is required to select the human 16S/WMS corpus")
    if not full_dump_path.is_file():
        raise FileNotFoundError(f"{full_dump_path} is missing; full_dump.csv is required")

    experiment_ids: set[tuple[str, str]] = set()
    with full_dump_path.open(encoding="utf-8", newline="") as handle:
        rows = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(rows)
        for row in reader:
            bsdb_id = clean_value(row.get("BSDB ID"))
            host_species = clean_value(row.get("Host species"))
            sequencing_type = clean_value(row.get("Sequencing type"))
            if (
                bsdb_id
                and host_species == HUMAN_HOST_SPECIES
                and sequencing_type in CORPUS_SEQUENCING_TYPES
            ):
                experiment_ids.add(
                    (
                        study_id_from_bsdb_id(bsdb_id),
                        experiment_id_from_bsdb_id(bsdb_id),
                    )
                )
    return experiment_ids


def build_binary_matrix(
    taxa_path: Path,
    top_taxa: int,
    max_profiles: int,
    full_dump_path: Path,
) -> tuple[pd.DataFrame, pd.Series]:
    """Build a binary profile x genus matrix from taxa.tsv."""
    human_experiment_ids = load_human_experiment_ids(full_dump_path)
    # A missing or empty corpus must never become an all-host fallback (#753, #787).
    if not human_experiment_ids:
        msg = (
            f"{full_dump_path} has no {HUMAN_HOST_SPECIES} 16S or WMS rows; check its "
            "Host species and Sequencing type columns"
        )
        raise ValueError(msg)
    taxa = pd.read_csv(
        taxa_path,
        sep="\t",
        usecols=["study_id", "experiment_id", "direction", "taxon_name", "taxonomic_rank"],
    )
    taxa = taxa[
        [
            (study_id, experiment_id) in human_experiment_ids
            for study_id, experiment_id in zip(
                taxa["study_id"],
                taxa["experiment_id"],
                strict=True,
            )
        ]
    ].copy()
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
        "BugSigDB Human-Host Genus Presence/Absence Clustered Heatmap\n"
        "Top prevalent genera; rows are human-host 16S and shotgun "
        "study|experiment|direction profiles",
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
        default=REPO_ROOT / "data/derived/bugsigdb_taxa.tsv",
        help="Input taxa TSV path.",
    )
    parser.add_argument(
        "--full-dump-input",
        type=Path,
        default=REPO_ROOT / "data/raw/full_dump.csv",
        help="Required BugSigDB full_dump.csv used to select Homo sapiens 16S/WMS profiles.",
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

    matrix, inc_dec = build_binary_matrix(
        args.input,
        args.top_taxa,
        args.max_profiles,
        args.full_dump_input,
    )

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
