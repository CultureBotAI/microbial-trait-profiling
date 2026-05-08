"""Create clustered heatmaps for condition/body-site trait modules."""

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


def compact_label(value: str, max_length: int = 46) -> str:
    """Shorten long labels for axes and legends."""
    value = clean_value(value)
    if len(value) <= max_length:
        return value
    return f"{value[: max_length - 3]}..."


def study_id_from_bsdb_id(bsdb_id: str) -> str:
    """Convert a full BugSigDB signature ID into the derived study ID format."""
    compact_id = bsdb_id.strip().removeprefix("bsdb:").split("/", 1)[0]
    return f"bsdb:{compact_id}"


def most_common_label(values: list[str]) -> str:
    """Return the most common non-empty label with deterministic ties."""
    labels = [clean_value(value) for value in values if clean_value(value) != "unknown"]
    if not labels:
        return "unknown"
    counts = Counter(labels)
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]


def load_raw_study_metadata(full_dump_path: Path) -> pd.DataFrame:
    """Extract per-study condition and body-site labels from BugSigDB full_dump.csv."""
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
    """Load compact study metadata from derived studies.tsv."""
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
    """Load readable labels and categories for trait columns."""
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


def load_aggregated_study_profiles(profile_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate inc/dec experiment profiles to study-level trait profiles."""
    profiles = pd.read_csv(profile_path, sep="\t")
    trait_columns = [
        column for column in profiles.columns if column not in PROFILE_METADATA_COLUMNS
    ]
    values = profiles[trait_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    study_matrix = values.groupby(profiles["study_id"], sort=False).sum()

    profile_metadata = profiles.groupby("study_id", sort=False).agg(
        profile_count=("experiment_id", "size"),
        experiment_count=("experiment_id", "nunique"),
    )
    direction_counts = pd.crosstab(profiles["study_id"], profiles["inc_dec"])
    for direction in ["inc", "dec"]:
        if direction not in direction_counts.columns:
            direction_counts[direction] = 0
    profile_metadata = profile_metadata.join(direction_counts[["inc", "dec"]], how="left").fillna(0)
    profile_metadata = profile_metadata.rename(
        columns={"inc": "increased_profile_count", "dec": "decreased_profile_count"}
    )
    profile_metadata["total_trait_signal"] = study_matrix.sum(axis=1)
    profile_metadata["trait_count"] = study_matrix.gt(0).sum(axis=1)
    return study_matrix, profile_metadata


def normalize_group_fields(group_by: list[str]) -> list[str]:
    """Validate supported metadata fields for grouping."""
    supported = {"condition", "body_site", "location", "host_species"}
    unknown = [field for field in group_by if field not in supported]
    if unknown:
        msg = f"Unsupported group field(s): {', '.join(unknown)}"
        raise ValueError(msg)
    return group_by


def build_condition_modules(
    study_matrix: pd.DataFrame,
    study_metadata: pd.DataFrame,
    group_by: list[str],
    min_studies: int,
    min_profiles: int,
    min_traits: int,
    min_total_signal: int,
    include_unknown_groups: bool,
    max_groups: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate study profiles into condition/body-site modules."""
    metadata = study_metadata.copy()
    for field in group_by:
        metadata[field] = metadata[field].map(clean_value)
    if not include_unknown_groups:
        for field in group_by:
            metadata = metadata[metadata[field] != "unknown"].copy()

    metadata = metadata.loc[metadata.index.intersection(study_matrix.index)].copy()
    metadata["module_id"] = metadata[group_by].agg(" | ".join, axis=1)
    group_matrix = study_matrix.loc[metadata.index].groupby(metadata["module_id"], sort=False).sum()

    group_metadata = metadata.groupby("module_id", sort=False).agg(
        study_count=("profile_count", "size"),
        profile_count=("profile_count", "sum"),
        experiment_count=("experiment_count", "sum"),
        condition=("condition", most_common_label),
        body_site=("body_site", most_common_label),
        location=("location", most_common_label),
        host_species=("host_species", most_common_label),
    )
    group_metadata["total_trait_signal"] = group_matrix.sum(axis=1)
    group_metadata["trait_count"] = group_matrix.gt(0).sum(axis=1)

    eligible = group_metadata[
        (group_metadata["study_count"] >= min_studies)
        & (group_metadata["profile_count"] >= min_profiles)
        & (group_metadata["trait_count"] >= min_traits)
        & (group_metadata["total_trait_signal"] >= min_total_signal)
    ].copy()
    if len(eligible) < 2:
        msg = "Fewer than two condition modules passed the selection filters."
        raise ValueError(msg)

    eligible = eligible.sort_values(
        ["study_count", "profile_count", "total_trait_signal", "trait_count"],
        ascending=False,
    ).head(max_groups)
    return group_matrix.loc[eligible.index], eligible


def group_top_categories(series: pd.Series, max_categories: int) -> pd.Series:
    """Keep the most common categories and group the rest as Other."""
    cleaned = series.map(clean_value)
    counts = cleaned[cleaned != "unknown"].value_counts()
    keep = set(counts.head(max_categories).index)
    grouped = cleaned.where(cleaned.isin(keep), "Other")
    grouped = grouped.where(cleaned != "unknown", "unknown")
    return grouped


def select_traits(
    group_matrix: pd.DataFrame,
    trait_labels: pd.DataFrame,
    max_traits: int,
    min_group_support: int,
    max_group_support_fraction: float,
    include_media_traits: bool,
) -> tuple[pd.Index, pd.DataFrame]:
    """Select condition-discriminating trait columns."""
    row_sums = group_matrix.sum(axis=1).replace(0, pd.NA)
    composition = group_matrix.div(row_sums, axis=0).fillna(0)
    support = group_matrix.gt(0).sum(axis=0).astype(int)
    total_signal = group_matrix.sum(axis=0).astype(float)
    variance = composition.var(axis=0)

    trait_metadata = pd.DataFrame(
        {
            "trait_id": group_matrix.columns,
            "module_support": support,
            "total_signal": total_signal,
            "composition_variance": variance,
        },
        index=group_matrix.columns,
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

    max_support = max(2, int(len(group_matrix) * max_group_support_fraction))
    eligible = trait_metadata[
        (trait_metadata["module_support"] >= min_group_support)
        & (trait_metadata["module_support"] <= max_support)
        & (trait_metadata["composition_variance"] > 0)
    ].copy()
    if not include_media_traits:
        eligible = eligible[
            ~eligible.index.to_series().str.contains("mediadive.medium", regex=False)
            & ~eligible["display_label"].str.contains("mediadive.medium", regex=False)
        ].copy()
    eligible = eligible.sort_values(
        ["selection_score", "composition_variance", "module_support"],
        ascending=False,
    ).head(max_traits)
    if len(eligible) < 2:
        msg = "Fewer than two traits passed the selection filters."
        raise ValueError(msg)
    return eligible.index, trait_metadata.loc[eligible.index]


def build_plot_matrix(group_matrix: pd.DataFrame, selected_traits: pd.Index) -> pd.DataFrame:
    """Create column-scaled log composition values for plotting."""
    row_sums = group_matrix.sum(axis=1).replace(0, pd.NA)
    composition = group_matrix[selected_traits].div(row_sums, axis=0).fillna(0)
    scaled = (composition * 1000).apply(lambda column: column.map(math.log1p))
    zscores = (scaled - scaled.mean(axis=0)) / scaled.std(axis=0).replace(0, pd.NA)
    return zscores.fillna(0).clip(-2.5, 2.5)


def build_display_matrix(
    plot_matrix: pd.DataFrame,
    group_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Apply readable labels to rows and columns."""
    display = plot_matrix.copy()
    display.index = [
        " | ".join(
            [
                compact_label(str(group_metadata.loc[module_id, "condition"]), 38),
                compact_label(str(group_metadata.loc[module_id, "body_site"]), 24),
                f"n={int(group_metadata.loc[module_id, 'study_count'])}",
            ]
        )
        for module_id in display.index
    ]
    display.columns = [
        trait_metadata.loc[trait_id, "display_label"] or trait_id for trait_id in display.columns
    ]
    return display


def palette_for_categories(categories: pd.Series, palette_name: str) -> dict[str, object]:
    """Build a stable palette for annotation categories."""
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
    group_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    output_path: Path,
) -> None:
    """Render and save the condition trait module heatmap."""
    sns.set_theme(style="white")
    display_to_module_id = dict(zip(display_matrix.index, group_metadata.index, strict=True))
    display_to_trait_id = dict(zip(display_matrix.columns, trait_metadata.index, strict=True))

    body_palette = palette_for_categories(group_metadata["body_site_group"], "tab20")
    size_palette = palette_for_categories(group_metadata["module_size_group"], "Set2")
    trait_palette = palette_for_categories(trait_metadata["display_category"], "Set2")

    row_module_ids = pd.Index([display_to_module_id[label] for label in display_matrix.index])
    row_colors = pd.DataFrame(
        {
            "body_site": group_metadata.loc[row_module_ids, "body_site_group"]
            .map(body_palette)
            .to_numpy(),
            "module_size": group_metadata.loc[row_module_ids, "module_size_group"]
            .map(size_palette)
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
        "Condition/Disease Trait Modules\nBugSigDB Condition Groups by Aggregate Trait Composition",
        y=1.08,
        fontsize=16,
        fontweight="bold",
    )
    cluster.ax_heatmap.set_xlabel("Trait")
    cluster.ax_heatmap.set_ylabel("Condition | body site | study count")
    for label in cluster.ax_heatmap.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(7)
    for label in cluster.ax_heatmap.get_yticklabels():
        label.set_fontsize(7)

    cluster.ax_row_dendrogram.legend(
        handles=legend_handles(body_palette, group_metadata["body_site_group"]),
        title="Body site",
        loc="center",
        bbox_to_anchor=(-0.28, 0.65),
        frameon=False,
    )
    cluster.fig.legend(
        handles=legend_handles(size_palette, group_metadata["module_size_group"]),
        title="Study count",
        loc="lower left",
        bbox_to_anchor=(0.02, 0.02),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
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

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cluster.fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(cluster.fig)


def write_outputs(
    plot_matrix: pd.DataFrame,
    display_matrix: pd.DataFrame,
    group_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    output_path: Path,
    matrix_output_path: Path,
    group_metadata_output_path: Path,
    trait_metadata_output_path: Path,
) -> None:
    """Write the plotted matrix and source metadata."""
    matrix_output_path.parent.mkdir(parents=True, exist_ok=True)
    display_matrix.to_csv(matrix_output_path, sep="\t")
    group_metadata.loc[plot_matrix.index].to_csv(group_metadata_output_path, sep="\t")
    trait_metadata.loc[plot_matrix.columns].to_csv(trait_metadata_output_path, sep="\t")
    print(f"Wrote heatmap: {output_path}")
    print(f"Wrote matrix: {matrix_output_path}")
    print(f"Wrote group metadata: {group_metadata_output_path}")
    print(f"Wrote trait metadata: {trait_metadata_output_path}")
    print(f"Condition modules plotted: {len(plot_matrix)}")
    print(f"Traits plotted: {len(plot_matrix.columns)}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Group BugSigDB studies by condition, body site, or other metadata, "
            "then cluster aggregate trait profiles to reveal condition-specific signatures."
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
        default=OUTPUT_DIR / "condition_trait_module_clustered_heatmap.png",
        help="Output heatmap image path.",
    )
    parser.add_argument(
        "--matrix-output",
        type=Path,
        default=OUTPUT_DIR / "condition_trait_module_clustered_heatmap_matrix.tsv",
        help="Output plotted z-score matrix TSV path.",
    )
    parser.add_argument(
        "--group-metadata-output",
        type=Path,
        default=OUTPUT_DIR / "condition_trait_module_group_metadata.tsv",
        help="Output selected condition module metadata TSV path.",
    )
    parser.add_argument(
        "--trait-metadata-output",
        type=Path,
        default=OUTPUT_DIR / "condition_trait_module_trait_metadata.tsv",
        help="Output selected trait metadata TSV path.",
    )
    parser.add_argument(
        "--group-by",
        nargs="+",
        default=["condition", "body_site"],
        help="Metadata fields for grouping: condition, body_site, location, host_species.",
    )
    parser.add_argument("--max-groups", type=int, default=80)
    parser.add_argument("--max-traits", type=int, default=60)
    parser.add_argument("--min-studies", type=int, default=2)
    parser.add_argument("--min-profiles", type=int, default=4)
    parser.add_argument("--min-traits", type=int, default=8)
    parser.add_argument("--min-total-signal", type=int, default=40)
    parser.add_argument("--min-group-support", type=int, default=4)
    parser.add_argument("--max-group-support-fraction", type=float, default=0.95)
    parser.add_argument("--max-body-site-categories", type=int, default=10)
    parser.add_argument(
        "--include-unknown-groups",
        action="store_true",
        help="Include groups where a grouping field is unknown.",
    )
    parser.add_argument(
        "--include-media-traits",
        action="store_true",
        help="Include MediaDive growth-medium traits in the selected traits.",
    )
    args = parser.parse_args()

    group_by = normalize_group_fields(args.group_by)
    study_matrix, profile_metadata = load_aggregated_study_profiles(args.study_trait_profiles_input)
    study_metadata = profile_metadata.join(load_study_titles(args.studies_input), how="left")
    study_metadata = study_metadata.join(load_raw_study_metadata(args.full_dump_input), how="left")
    for column in ["condition", "body_site", "location", "host_species"]:
        study_metadata[column] = study_metadata[column].map(clean_value)

    group_matrix, group_metadata = build_condition_modules(
        study_matrix,
        study_metadata,
        group_by,
        args.min_studies,
        args.min_profiles,
        args.min_traits,
        args.min_total_signal,
        args.include_unknown_groups,
        args.max_groups,
    )
    group_metadata["body_site_group"] = group_top_categories(
        group_metadata["body_site"],
        args.max_body_site_categories,
    )
    group_metadata["module_size_group"] = pd.cut(
        group_metadata["study_count"],
        bins=[0, 2, 5, 10, float("inf")],
        labels=["2", "3-5", "6-10", "11+"],
        include_lowest=True,
    ).astype(str)

    trait_labels = load_trait_labels(args.trait_columns_input)
    selected_traits, trait_metadata = select_traits(
        group_matrix,
        trait_labels,
        args.max_traits,
        args.min_group_support,
        args.max_group_support_fraction,
        args.include_media_traits,
    )
    plot_matrix = build_plot_matrix(group_matrix, selected_traits)
    display_matrix = build_display_matrix(plot_matrix, group_metadata, trait_metadata)
    save_heatmap(display_matrix, group_metadata, trait_metadata, args.output)
    write_outputs(
        plot_matrix,
        display_matrix,
        group_metadata,
        trait_metadata,
        args.output,
        args.matrix_output,
        args.group_metadata_output,
        args.trait_metadata_output,
    )


if __name__ == "__main__":
    main()
