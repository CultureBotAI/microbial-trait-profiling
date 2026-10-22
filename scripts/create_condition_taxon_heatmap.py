"""Create a clustered heatmap of condition x genus prevalence."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
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


def load_study_condition_pairs(full_dump_path: Path) -> pd.DataFrame:
    """Extract unique human-host study/experiment condition pairs from full_dump.csv.

    Only 16S and WMS experiments are kept, the analysis corpus's rule; a blank
    sequencing type is left out (#742).
    """
    experiment_conditions: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    with full_dump_path.open(encoding="utf-8", newline="") as handle:
        rows = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(rows)
        for row in reader:
            bsdb_id = clean_value(row.get("BSDB ID"))
            condition = clean_value(row.get("Condition"))
            host_species = clean_value(row.get("Host species"))
            sequencing_type = clean_value(row.get("Sequencing type"))
            if (
                not bsdb_id
                or not condition
                or host_species != HUMAN_HOST_SPECIES
                or sequencing_type not in CORPUS_SEQUENCING_TYPES
            ):
                continue
            study_id = study_id_from_bsdb_id(bsdb_id)
            experiment_id = experiment_id_from_bsdb_id(bsdb_id)
            experiment_conditions[(study_id, experiment_id)].add(condition)

    pairs = [
        {
            "study_id": study_id,
            "experiment_id": experiment_id,
            "condition": condition,
        }
        for (study_id, experiment_id), conditions in experiment_conditions.items()
        for condition in sorted(conditions)
    ]
    return pd.DataFrame(pairs)


def rank_by_study_count(frame: pd.DataFrame, column: str) -> pd.Series:
    """Count each value's unique studies, most first, ties in ascending name order.

    The stable sort makes a top-N cut that falls inside a tie keep the same values
    on every run (#794).
    """
    counts = frame.groupby(column)["study_id"].nunique().sort_index()
    return counts.sort_values(ascending=False, kind="stable")


def build_condition_genus_matrix(
    taxa_path: Path,
    full_dump_path: Path,
    top_conditions: int,
    top_taxa: int,
) -> tuple[pd.DataFrame, pd.Series]:
    """Build a condition x genus prevalence matrix.

    The taxa table's experiments are the population: dump condition labels are kept
    only for them, so a contrast-population subset ranks conditions and counts each
    condition's studies from its own experiments, not the whole corpus (#768).
    """
    study_conditions = load_study_condition_pairs(full_dump_path)
    if study_conditions.empty:
        msg = "No condition labels found in full_dump.csv"
        raise ValueError(msg)

    taxa = pd.read_csv(
        taxa_path,
        sep="\t",
        usecols=["study_id", "experiment_id", "taxon_name", "taxonomic_rank"],
    )
    study_conditions = study_conditions.merge(
        taxa[["study_id", "experiment_id"]].drop_duplicates(), on=["study_id", "experiment_id"]
    )
    if study_conditions.empty:
        msg = f"No full_dump.csv condition labels for the experiments in {taxa_path}"
        raise ValueError(msg)
    taxa = taxa[taxa["taxonomic_rank"] == "genus"].copy()
    taxa = taxa[["study_id", "experiment_id", "taxon_name"]].drop_duplicates()

    condition_sizes = rank_by_study_count(study_conditions, "condition").head(top_conditions)
    study_conditions = study_conditions[
        study_conditions["condition"].isin(condition_sizes.index)
    ].copy()

    merged = study_conditions.merge(taxa, on=["study_id", "experiment_id"], how="inner")
    if merged.empty:
        msg = "No overlap between condition labels and genus taxa"
        raise ValueError(msg)

    top_taxon_names = rank_by_study_count(merged, "taxon_name").head(top_taxa).index
    merged = merged[merged["taxon_name"].isin(top_taxon_names)].copy()

    unique_condition_taxa = merged[["condition", "study_id", "taxon_name"]].drop_duplicates()

    counts = (
        unique_condition_taxa.groupby(["condition", "taxon_name"])["study_id"]
        .nunique()
        .rename("study_count")
        .reset_index()
    )
    counts["condition_study_count"] = counts["condition"].map(condition_sizes)
    counts["prevalence"] = counts["study_count"] / counts["condition_study_count"]

    matrix = (
        counts.pivot(index="condition", columns="taxon_name", values="prevalence")
        .fillna(0.0)
        .astype(float)
    )
    matrix = matrix.loc[condition_sizes.index.intersection(matrix.index)]
    matrix = matrix.loc[:, matrix.mean(axis=0).sort_values(ascending=False).index]
    row_annotation = condition_sizes.loc[matrix.index]
    return matrix, row_annotation


def save_condition_heatmap(
    matrix: pd.DataFrame, study_counts: pd.Series, output_path: Path
) -> None:
    """Render and save the clustered heatmap."""
    sns.set_theme(style="white")
    display_index = [
        f"{condition} (n={int(study_counts.loc[condition])})" for condition in matrix.index
    ]
    plot_matrix = matrix.copy()
    plot_matrix.index = display_index

    cluster = sns.clustermap(
        plot_matrix,
        cmap="YlGnBu",
        method="average",
        metric="euclidean",
        linewidths=0,
        xticklabels=True,
        yticklabels=True,
        figsize=(22, 16),
        cbar_pos=(0.02, 0.82, 0.02, 0.12),
        dendrogram_ratio=(0.14, 0.12),
    )

    cluster.ax_heatmap.set_xlabel("Genus")
    cluster.ax_heatmap.set_ylabel("Condition")
    cluster.ax_heatmap.set_title(
        "BugSigDB Human-Host Condition x Genus Clustered Heatmap\n"
        "Cell values are the fraction of human-host 16S and shotgun studies in each condition "
        "containing the genus",
        pad=18,
    )

    for label in cluster.ax_heatmap.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(8)

    for label in cluster.ax_heatmap.get_yticklabels():
        label.set_fontsize(9)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cluster.fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(cluster.fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Create a clustered heatmap of condition x genus prevalence using taxa.tsv "
            "and full_dump.csv."
        )
    )
    parser.add_argument(
        "--taxa-input",
        type=Path,
        default=REPO_ROOT / "data/derived/bugsigdb_taxa.tsv",
        help=(
            "Input taxa TSV path. Its experiments are the population whose conditions are "
            "ranked and whose studies form each condition's denominator."
        ),
    )
    parser.add_argument(
        "--full-dump-input",
        type=Path,
        default=REPO_ROOT / "data/raw/full_dump.csv",
        help=(
            "Input full_dump.csv path. It supplies the condition labels and keeps Homo "
            "sapiens 16S and WMS experiments."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_DIR / "condition_taxon_clustered_heatmap_top40x100.png",
        help="Output heatmap image path.",
    )
    parser.add_argument(
        "--matrix-output",
        type=Path,
        default=OUTPUT_DIR / "condition_taxon_clustered_heatmap_top40x100_matrix.tsv",
        help="Output matrix TSV path.",
    )
    parser.add_argument(
        "--top-conditions",
        type=int,
        default=40,
        help="Number of most common conditions to keep.",
    )
    parser.add_argument(
        "--top-taxa",
        type=int,
        default=100,
        help="Number of most prevalent genera to keep.",
    )
    args = parser.parse_args()

    matrix, study_counts = build_condition_genus_matrix(
        args.taxa_input,
        args.full_dump_input,
        args.top_conditions,
        args.top_taxa,
    )

    matrix_to_save = matrix.copy()
    matrix_to_save.insert(0, "study_count", study_counts.astype(int))
    matrix_to_save.to_csv(args.matrix_output, sep="\t")

    save_condition_heatmap(matrix, study_counts, args.output)

    print(f"Wrote matrix: {args.matrix_output}")
    print(f"Wrote heatmap: {args.output}")
    print(f"Conditions plotted: {len(matrix)}")
    print(f"Genera plotted: {len(matrix.columns)}")


if __name__ == "__main__":
    main()
