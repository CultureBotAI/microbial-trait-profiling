"""Create a clustered heatmap of BugSigDB studies by aggregate trait profiles."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "trait-visualization"
PROFILE_METADATA_COLUMNS = ["study_id", "experiment_id", "inc_dec"]


def clean_value(value: object) -> str:
    """Normalize empty BugSigDB metadata values."""
    text = str(value or "").strip()
    if not text or text.upper() == "NA" or text.lower() == "nan":
        return "unknown"
    return text


def compact_label(value: str, max_length: int = 42) -> str:
    """Shorten long labels for heatmap axes and legends."""
    value = clean_value(value)
    if len(value) <= max_length:
        return value
    return f"{value[: max_length - 3]}..."


def study_id_from_bsdb_id(bsdb_id: str) -> str:
    """Convert a full BugSigDB signature ID into the derived study ID format."""
    compact_id = bsdb_id.strip().removeprefix("bsdb:").split("/", 1)[0]
    return f"bsdb:{compact_id}"


def most_common_label(values: list[str]) -> str:
    """Return the most common non-empty label, with deterministic ties."""
    labels = [clean_value(value) for value in values if clean_value(value) != "unknown"]
    if not labels:
        return "unknown"
    counts = Counter(labels)
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]


def load_raw_study_metadata(full_dump_path: Path) -> pd.DataFrame:
    """Extract condition and body-site labels from BugSigDB full_dump.csv."""
    values_by_study: dict[str, dict[str, list[str]]] = defaultdict(
        lambda: {
            "condition": [],
            "body_site": [],
            "location": [],
            "host_species": [],
        }
    )
    with full_dump_path.open(encoding="utf-8", newline="") as handle:
        rows = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(rows)
        for row in reader:
            bsdb_id = clean_value(row.get("BSDB ID"))
            if bsdb_id == "unknown":
                continue
            study_id = study_id_from_bsdb_id(bsdb_id)
            values_by_study[study_id]["condition"].append(clean_value(row.get("Condition")))
            values_by_study[study_id]["body_site"].append(clean_value(row.get("Body site")))
            values_by_study[study_id]["location"].append(
                clean_value(row.get("Location of subjects"))
            )
            values_by_study[study_id]["host_species"].append(clean_value(row.get("Host species")))

    records = []
    for study_id, metadata in values_by_study.items():
        records.append(
            {
                "study_id": study_id,
                "condition": most_common_label(metadata["condition"]),
                "body_site": most_common_label(metadata["body_site"]),
                "location": most_common_label(metadata["location"]),
                "host_species": most_common_label(metadata["host_species"]),
                "condition_count": len(
                    {value for value in metadata["condition"] if clean_value(value) != "unknown"}
                ),
                "body_site_count": len(
                    {value for value in metadata["body_site"] if clean_value(value) != "unknown"}
                ),
            }
        )
    return pd.DataFrame.from_records(records).set_index("study_id")


def load_study_titles(studies_path: Path) -> pd.DataFrame:
    """Load compact title metadata from derived studies.tsv."""
    studies = pd.read_csv(
        studies_path,
        sep="\t",
        usecols=["study_id", "title", "publication_year", "taxon_count", "signature_count"],
        dtype=str,
    )
    for column in ["title", "publication_year"]:
        studies[column] = studies[column].map(clean_value)
    for column in ["taxon_count", "signature_count"]:
        studies[column] = pd.to_numeric(studies[column], errors="coerce").fillna(0).astype(int)
    return studies.set_index("study_id")


def load_trait_labels(columns_path: Path) -> pd.DataFrame:
    """Load display labels and categories for trait columns."""
    columns = pd.read_csv(columns_path, sep="\t", dtype=str).fillna("")
    columns = columns.set_index("column_id", drop=False)
    columns["display_label"] = [
        format_trait_label(column_id, row) for column_id, row in columns.iterrows()
    ]
    columns["display_category"] = columns["object_category"].str.removeprefix("biolink:")
    columns.loc[columns["display_category"] == "", "display_category"] = "unknown"
    return columns


def format_trait_label(column_id: str, row: pd.Series) -> str:
    """Build a compact two-line trait label."""
    name = row.get("object_name", "") or row.get("column_label", "") or column_id
    name = compact_label(str(name).replace("_", " "), 34)
    short_id = column_id.split("|")[-1].replace("pathways:", "").replace("carbon_substrates:", "")
    short_id = compact_label(short_id, 22)
    return f"{name}\n{short_id}"


def group_top_categories(series: pd.Series, max_categories: int) -> pd.Series:
    """Keep the most common categories and group the rest as Other."""
    cleaned = series.map(clean_value)
    counts = cleaned[cleaned != "unknown"].value_counts()
    keep = set(counts.head(max_categories).index)
    grouped = cleaned.where(cleaned.isin(keep), "Other")
    grouped = grouped.where(cleaned != "unknown", "unknown")
    return grouped


def load_aggregated_profiles(profile_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Aggregate inc/dec experiment-level trait profiles to study-level profiles."""
    profiles = pd.read_csv(profile_path, sep="\t")
    trait_columns = [
        column for column in profiles.columns if column not in PROFILE_METADATA_COLUMNS
    ]
    values = profiles[trait_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    study_matrix = values.groupby(profiles["study_id"], sort=False).sum()

    profile_counts = profiles.groupby("study_id", sort=False).agg(
        profile_count=("experiment_id", "size"),
        experiment_count=("experiment_id", "nunique"),
    )
    direction_counts = pd.crosstab(profiles["study_id"], profiles["inc_dec"])
    for direction in ["inc", "dec"]:
        if direction not in direction_counts.columns:
            direction_counts[direction] = 0
    metadata = profile_counts.join(direction_counts[["inc", "dec"]], how="left").fillna(0)
    metadata = metadata.rename(
        columns={"inc": "increased_profile_count", "dec": "decreased_profile_count"}
    )
    metadata["increased_profile_count"] = metadata["increased_profile_count"].astype(int)
    metadata["decreased_profile_count"] = metadata["decreased_profile_count"].astype(int)
    metadata["total_trait_signal"] = study_matrix.sum(axis=1)
    metadata["trait_count"] = study_matrix.gt(0).sum(axis=1)
    return study_matrix, metadata, trait_columns


def select_studies(
    study_matrix: pd.DataFrame,
    study_metadata: pd.DataFrame,
    max_studies: int,
    min_profiles: int,
    min_traits: int,
    min_total_signal: int,
) -> pd.Index:
    """Select studies with enough profile and trait signal for a readable heatmap."""
    eligible = study_metadata[
        (study_metadata["profile_count"] >= min_profiles)
        & (study_metadata["trait_count"] >= min_traits)
        & (study_metadata["total_trait_signal"] >= min_total_signal)
    ].copy()
    if len(eligible) < 2:
        msg = "Fewer than two studies passed the selection filters."
        raise ValueError(msg)
    eligible = eligible.loc[study_matrix.index.intersection(eligible.index)]
    eligible = eligible.sort_values(
        ["profile_count", "total_trait_signal", "trait_count"],
        ascending=False,
    )
    return eligible.head(max_studies).index


def select_traits(
    selected_study_matrix: pd.DataFrame,
    trait_labels: pd.DataFrame,
    max_traits: int,
    min_study_support: int,
    max_study_support_fraction: float,
    include_media_traits: bool,
) -> tuple[pd.Index, pd.DataFrame]:
    """Select traits that vary across the selected study profiles."""
    row_sums = selected_study_matrix.sum(axis=1).replace(0, pd.NA)
    composition = selected_study_matrix.div(row_sums, axis=0).fillna(0)
    support = selected_study_matrix.gt(0).sum(axis=0).astype(int)
    total_signal = selected_study_matrix.sum(axis=0).astype(float)
    variance = composition.var(axis=0)

    trait_metadata = pd.DataFrame(
        {
            "trait_id": composition.columns,
            "study_support": support,
            "total_signal": total_signal,
            "composition_variance": variance,
        },
        index=composition.columns,
    ).join(
        trait_labels[
            [
                "display_label",
                "column_label",
                "object_name",
                "object_category",
                "display_category",
            ]
        ],
        how="left",
    )
    trait_metadata[["display_label", "display_category"]] = trait_metadata[
        ["display_label", "display_category"]
    ].fillna("")
    trait_metadata["selection_score"] = [
        variance_value * math.log1p(signal)
        for variance_value, signal in zip(
            trait_metadata["composition_variance"],
            trait_metadata["total_signal"],
            strict=True,
        )
    ]

    max_support = max(2, int(len(selected_study_matrix) * max_study_support_fraction))
    eligible = trait_metadata[
        (trait_metadata["study_support"] >= min_study_support)
        & (trait_metadata["study_support"] <= max_support)
        & (trait_metadata["composition_variance"] > 0)
    ].copy()
    if not include_media_traits:
        eligible = eligible[
            ~eligible.index.to_series().str.contains("mediadive.medium", regex=False)
            & ~eligible["display_label"].str.contains("mediadive.medium", regex=False)
        ].copy()
    eligible = eligible.sort_values(
        ["selection_score", "composition_variance", "study_support"],
        ascending=False,
    ).head(max_traits)
    if len(eligible) < 2:
        msg = "Fewer than two traits passed the selection filters."
        raise ValueError(msg)
    return eligible.index, trait_metadata.loc[eligible.index]


def build_plot_matrix(
    selected_study_matrix: pd.DataFrame,
    selected_traits: pd.Index,
) -> pd.DataFrame:
    """Create column-scaled log composition values for plotting."""
    row_sums = selected_study_matrix.sum(axis=1).replace(0, pd.NA)
    composition = selected_study_matrix[selected_traits].div(row_sums, axis=0).fillna(0)
    scaled = (composition * 1000).apply(lambda column: column.map(math.log1p))
    zscores = (scaled - scaled.mean(axis=0)) / scaled.std(axis=0).replace(0, pd.NA)
    return zscores.fillna(0).clip(-2.5, 2.5)


def build_display_matrix(
    plot_matrix: pd.DataFrame,
    study_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Apply readable labels for the plotted matrix."""
    display = plot_matrix.copy()
    display.index = [
        " | ".join(
            [
                study_id,
                compact_label(study_metadata.loc[study_id, "body_site_group"], 24),
                compact_label(study_metadata.loc[study_id, "condition_group"], 36),
            ]
        )
        for study_id in display.index
    ]
    display.columns = [
        trait_metadata.loc[trait_id, "display_label"] or trait_id for trait_id in display.columns
    ]
    return display


def palette_for_categories(categories: pd.Series, palette_name: str) -> dict[str, object]:
    """Build a stable palette for row or column annotation categories."""
    ordered = [
        category for category in sorted(categories.unique()) if category not in {"Other", "unknown"}
    ]
    colors: dict[str, object] = {}
    if ordered:
        colors = dict(
            zip(
                ordered,
                sns.color_palette(palette_name, n_colors=len(ordered)),
                strict=True,
            )
        )
    colors["Other"] = "#b0b0b0"
    colors["unknown"] = "#e5e7eb"
    return colors


def legend_handles(palette: dict[str, object], present_categories: pd.Series) -> list[Patch]:
    """Build legend handles only for categories present in the plot."""
    present = set(present_categories)
    return [Patch(color=color, label=label) for label, color in palette.items() if label in present]


def save_heatmap(
    display_matrix: pd.DataFrame,
    study_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    output_path: Path,
) -> None:
    """Render and save the clustered study profile heatmap."""
    sns.set_theme(style="white")
    display_to_study_id = {
        display_label: display_label.split(" | ", 1)[0] for display_label in display_matrix.index
    }
    display_to_trait_id = dict(zip(display_matrix.columns, trait_metadata.index, strict=True))

    body_palette = palette_for_categories(study_metadata["body_site_group"], "tab20")
    condition_palette = palette_for_categories(study_metadata["condition_group"], "Set3")
    trait_palette = palette_for_categories(trait_metadata["display_category"], "Set2")

    row_study_ids = pd.Index([display_to_study_id[label] for label in display_matrix.index])
    row_colors = pd.DataFrame(
        {
            "body_site": study_metadata.loc[row_study_ids, "body_site_group"]
            .map(body_palette)
            .to_numpy(),
            "condition": study_metadata.loc[row_study_ids, "condition_group"]
            .map(condition_palette)
            .to_numpy(),
        },
        index=display_matrix.index,
    )
    col_trait_ids = pd.Index([display_to_trait_id[label] for label in display_matrix.columns])
    col_colors = pd.DataFrame(
        {
            "trait_category": trait_metadata.loc[col_trait_ids, "display_category"]
            .map(trait_palette)
            .to_numpy(),
        },
        index=display_matrix.columns,
    )

    cluster = sns.clustermap(
        display_matrix,
        cmap="vlag",
        center=0,
        method="average",
        metric="correlation",
        row_colors=row_colors,
        col_colors=col_colors,
        linewidths=0.12,
        linecolor="#e5e7eb",
        xticklabels=True,
        yticklabels=True,
        figsize=(24, 18),
        cbar_pos=(0.02, 0.80, 0.02, 0.14),
        dendrogram_ratio=(0.12, 0.10),
        colors_ratio=(0.03, 0.03),
    )
    cluster.fig.suptitle(
        "Study Trait Profile Clustering\nBugSigDB Studies by Aggregate Trait Composition",
        y=1.08,
        fontsize=16,
        fontweight="bold",
    )
    cluster.ax_heatmap.set_xlabel("Trait")
    cluster.ax_heatmap.set_ylabel("Study | body site | condition")
    for label in cluster.ax_heatmap.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(7)
    for label in cluster.ax_heatmap.get_yticklabels():
        label.set_fontsize(7)

    cluster.ax_row_dendrogram.legend(
        handles=legend_handles(body_palette, study_metadata["body_site_group"]),
        title="Body site",
        loc="center",
        bbox_to_anchor=(-0.25, 0.65),
        frameon=False,
    )
    cluster.fig.legend(
        handles=legend_handles(trait_palette, trait_metadata["display_category"]),
        title="Trait category",
        loc="upper center",
        bbox_to_anchor=(0.5, 1.16),
        ncol=4,
        frameon=False,
        fontsize=9,
        title_fontsize=10,
    )
    cluster.fig.legend(
        handles=legend_handles(condition_palette, study_metadata["condition_group"]),
        title="Condition",
        loc="lower left",
        bbox_to_anchor=(0.02, 0.02),
        ncol=1,
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cluster.fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(cluster.fig)


def write_outputs(
    plot_matrix: pd.DataFrame,
    display_matrix: pd.DataFrame,
    study_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    output_path: Path,
    matrix_output_path: Path,
    study_metadata_output_path: Path,
    trait_metadata_output_path: Path,
) -> None:
    """Write the plotted matrix and source metadata."""
    matrix_output_path.parent.mkdir(parents=True, exist_ok=True)
    display_matrix.to_csv(matrix_output_path, sep="\t")
    study_ids = plot_matrix.index
    trait_ids = plot_matrix.columns
    study_metadata.loc[study_ids].to_csv(study_metadata_output_path, sep="\t")
    trait_metadata.loc[trait_ids].to_csv(trait_metadata_output_path, sep="\t")
    print(f"Wrote heatmap: {output_path}")
    print(f"Wrote matrix: {matrix_output_path}")
    print(f"Wrote study metadata: {study_metadata_output_path}")
    print(f"Wrote trait metadata: {trait_metadata_output_path}")
    print(f"Studies plotted: {len(study_ids)}")
    print(f"Traits plotted: {len(trait_ids)}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Cluster studies by aggregate trait profiles to inspect whether similar "
            "conditions or body sites converge on similar microbial traits."
        )
    )
    parser.add_argument(
        "--study-trait-profiles-input",
        type=Path,
        default=REPO_ROOT / "data/derived/study_trait_profiles_wide.tsv",
        help="Input wide study trait profiles TSV path.",
    )
    parser.add_argument(
        "--trait-columns-input",
        type=Path,
        default=REPO_ROOT / "data/derived/taxa_trait_matrix_columns.tsv",
        help="Input trait column metadata TSV path.",
    )
    parser.add_argument(
        "--studies-input",
        type=Path,
        default=REPO_ROOT / "data/derived/studies.tsv",
        help="Input derived studies TSV path.",
    )
    parser.add_argument(
        "--full-dump-input",
        type=Path,
        default=REPO_ROOT / "data/raw/full_dump.csv",
        help="Input BugSigDB full_dump.csv path.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_DIR / "study_trait_profile_clustered_heatmap.png",
        help="Output heatmap image path.",
    )
    parser.add_argument(
        "--matrix-output",
        type=Path,
        default=OUTPUT_DIR / "study_trait_profile_clustered_heatmap_matrix.tsv",
        help="Output plotted z-score matrix TSV path.",
    )
    parser.add_argument(
        "--study-metadata-output",
        type=Path,
        default=OUTPUT_DIR / "study_trait_profile_study_metadata.tsv",
        help="Output selected study metadata TSV path.",
    )
    parser.add_argument(
        "--trait-metadata-output",
        type=Path,
        default=OUTPUT_DIR / "study_trait_profile_trait_metadata.tsv",
        help="Output selected trait metadata TSV path.",
    )
    parser.add_argument("--max-studies", type=int, default=80)
    parser.add_argument("--max-traits", type=int, default=60)
    parser.add_argument("--min-profiles", type=int, default=2)
    parser.add_argument("--min-traits", type=int, default=8)
    parser.add_argument("--min-total-signal", type=int, default=20)
    parser.add_argument("--min-study-support", type=int, default=6)
    parser.add_argument("--max-study-support-fraction", type=float, default=0.95)
    parser.add_argument("--max-body-site-categories", type=int, default=8)
    parser.add_argument("--max-condition-categories", type=int, default=8)
    parser.add_argument(
        "--include-media-traits",
        action="store_true",
        help="Include MediaDive growth-medium traits in the selected traits.",
    )
    args = parser.parse_args()

    study_matrix, study_profile_metadata, _trait_columns = load_aggregated_profiles(
        args.study_trait_profiles_input
    )
    study_metadata = study_profile_metadata.join(load_study_titles(args.studies_input), how="left")
    study_metadata = study_metadata.join(load_raw_study_metadata(args.full_dump_input), how="left")
    for column in [
        "title",
        "publication_year",
        "condition",
        "body_site",
        "location",
        "host_species",
    ]:
        study_metadata[column] = study_metadata[column].map(clean_value)

    selected_studies = select_studies(
        study_matrix,
        study_metadata,
        args.max_studies,
        args.min_profiles,
        args.min_traits,
        args.min_total_signal,
    )
    selected_study_matrix = study_matrix.loc[selected_studies]
    study_metadata = study_metadata.loc[selected_studies].copy()
    study_metadata["body_site_group"] = group_top_categories(
        study_metadata["body_site"],
        args.max_body_site_categories,
    )
    study_metadata["condition_group"] = group_top_categories(
        study_metadata["condition"],
        args.max_condition_categories,
    )

    trait_labels = load_trait_labels(args.trait_columns_input)
    selected_traits, trait_metadata = select_traits(
        selected_study_matrix,
        trait_labels,
        args.max_traits,
        args.min_study_support,
        args.max_study_support_fraction,
        args.include_media_traits,
    )
    plot_matrix = build_plot_matrix(selected_study_matrix, selected_traits)
    display_matrix = build_display_matrix(plot_matrix, study_metadata, trait_metadata)
    save_heatmap(display_matrix, study_metadata, trait_metadata, args.output)
    write_outputs(
        plot_matrix,
        display_matrix,
        study_metadata,
        trait_metadata,
        args.output,
        args.matrix_output,
        args.study_metadata_output,
        args.trait_metadata_output,
    )


if __name__ == "__main__":
    main()
