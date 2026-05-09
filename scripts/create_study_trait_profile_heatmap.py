"""Create a clustered heatmap of BugSigDB study/experiment/direction trait profiles."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import leaves_list, linkage

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "trait-visualization"
PROFILE_METADATA_COLUMNS = ["study_id", "experiment_id", "inc_dec"]
DIRECTION_LABELS = {"inc": "increased", "dec": "decreased"}


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


def profile_id_from_parts(study_id: str, experiment_id: str, inc_dec: str) -> str:
    """Build a stable row ID for one study/experiment/direction profile."""
    return f"{study_id}|{experiment_id}|{inc_dec}"


def load_profile_matrix(profile_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Load one row per study, experiment, and inc/dec trait profile."""
    profiles = pd.read_csv(profile_path, sep="\t")
    trait_columns = [
        column for column in profiles.columns if column not in PROFILE_METADATA_COLUMNS
    ]
    values = profiles[trait_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    profile_ids = [
        profile_id_from_parts(str(study_id), str(experiment_id), str(inc_dec))
        for study_id, experiment_id, inc_dec in profiles[PROFILE_METADATA_COLUMNS].itertuples(
            index=False,
            name=None,
        )
    ]
    values.index = profile_ids
    profile_matrix = values.groupby(level=0, sort=False).sum()

    metadata = profiles[PROFILE_METADATA_COLUMNS].copy()
    metadata["profile_id"] = profile_ids
    metadata = metadata.drop_duplicates(subset=["profile_id"], keep="first").set_index("profile_id")
    metadata["direction_label"] = (
        metadata["inc_dec"].map(DIRECTION_LABELS).fillna(metadata["inc_dec"])
    )
    metadata["total_trait_signal"] = profile_matrix.sum(axis=1)
    metadata["trait_count"] = profile_matrix.gt(0).sum(axis=1)
    return profile_matrix, metadata, trait_columns


def select_profiles(
    profile_matrix: pd.DataFrame,
    profile_metadata: pd.DataFrame,
    max_profiles: int,
    min_traits: int,
    min_total_signal: int,
) -> pd.Index:
    """Select profile rows for the heatmap; by default this keeps every profile."""
    eligible = profile_metadata[
        (profile_metadata["trait_count"] >= min_traits)
        & (profile_metadata["total_trait_signal"] >= min_total_signal)
    ].copy()
    if len(eligible) < 2:
        msg = "Fewer than two profiles passed the selection filters."
        raise ValueError(msg)
    eligible = eligible.loc[profile_matrix.index.intersection(eligible.index)]
    if max_profiles <= 0:
        return eligible.index

    eligible = eligible.sort_values(
        ["total_trait_signal", "trait_count"],
        ascending=False,
    )
    return eligible.head(max_profiles).index


def select_traits(
    selected_profile_matrix: pd.DataFrame,
    trait_labels: pd.DataFrame,
    max_traits: int,
    min_profile_support: int,
    max_profile_support_fraction: float,
    include_media_traits: bool,
) -> tuple[pd.Index, pd.DataFrame]:
    """Select traits that vary across the selected study profiles."""
    row_sums = selected_profile_matrix.sum(axis=1).replace(0, pd.NA)
    composition = selected_profile_matrix.div(row_sums, axis=0).fillna(0)
    support = selected_profile_matrix.gt(0).sum(axis=0).astype(int)
    total_signal = selected_profile_matrix.sum(axis=0).astype(float)
    variance = composition.var(axis=0)

    trait_metadata = pd.DataFrame(
        {
            "trait_id": composition.columns,
            "profile_support": support,
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

    max_support = max(2, int(len(selected_profile_matrix) * max_profile_support_fraction))
    eligible = trait_metadata[
        (trait_metadata["profile_support"] >= min_profile_support)
        & (trait_metadata["profile_support"] <= max_support)
        & (trait_metadata["composition_variance"] > 0)
    ].copy()
    if not include_media_traits:
        eligible = eligible[
            ~eligible.index.to_series().str.contains("mediadive.medium", regex=False)
            & ~eligible["display_label"].str.contains("mediadive.medium", regex=False)
        ].copy()
    eligible = eligible.sort_values(
        ["selection_score", "composition_variance", "profile_support"],
        ascending=False,
    ).head(max_traits)
    if len(eligible) < 2:
        msg = "Fewer than two traits passed the selection filters."
        raise ValueError(msg)
    return eligible.index, trait_metadata.loc[eligible.index]


def build_plot_matrix(
    selected_profile_matrix: pd.DataFrame,
    selected_traits: pd.Index,
) -> pd.DataFrame:
    """Create column z-scores from log-scaled trait fractions."""
    row_sums = selected_profile_matrix.sum(axis=1).replace(0, pd.NA)
    composition = selected_profile_matrix[selected_traits].div(row_sums, axis=0).fillna(0)
    scaled = (composition * 1000).apply(lambda column: column.map(math.log1p))
    zscores = (scaled - scaled.mean(axis=0)) / scaled.std(axis=0, ddof=0).replace(0, pd.NA)
    return zscores.fillna(0)


def cluster_index(indexed_matrix: pd.DataFrame) -> pd.Index:
    """Return a clustered row order for a block, falling back to current order."""
    if len(indexed_matrix) < 2:
        return indexed_matrix.index
    values = indexed_matrix.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        return indexed_matrix.index
    try:
        row_linkage = linkage(values, method="average", metric="euclidean")
    except ValueError:
        return indexed_matrix.index
    return indexed_matrix.index[leaves_list(row_linkage)]


def order_profiles_by_direction(
    plot_matrix: pd.DataFrame,
    profile_metadata: pd.DataFrame,
) -> pd.Index:
    """Cluster rows within increased/decreased sections and keep sections separate."""
    ordered: list[str] = []
    for direction in ["inc", "dec"]:
        direction_index = profile_metadata.index[profile_metadata["inc_dec"] == direction]
        direction_index = direction_index.intersection(plot_matrix.index)
        ordered.extend(cluster_index(plot_matrix.loc[direction_index]).tolist())

    remaining = plot_matrix.index.difference(pd.Index(ordered), sort=False)
    ordered.extend(cluster_index(plot_matrix.loc[remaining]).tolist())
    return pd.Index(ordered)


def most_common_series_value(series: pd.Series) -> str:
    """Return the most common category value in a display bin."""
    if series.empty:
        return "unknown"
    counts = series.map(clean_value).value_counts()
    if counts.empty:
        return "unknown"
    return str(counts.sort_values(ascending=False).index[0])


def allocate_display_bins(profile_metadata: pd.DataFrame, max_display_rows: int) -> dict[str, int]:
    """Allocate rendered row bins across direction sections."""
    if max_display_rows <= 0 or len(profile_metadata) <= max_display_rows:
        return {
            direction: int((profile_metadata["inc_dec"] == direction).sum())
            for direction in profile_metadata["inc_dec"].dropna().unique()
        }

    directions = [
        direction
        for direction in ["inc", "dec"]
        if int((profile_metadata["inc_dec"] == direction).sum()) > 0
    ]
    directions.extend(
        direction
        for direction in sorted(set(profile_metadata["inc_dec"]) - set(directions))
        if int((profile_metadata["inc_dec"] == direction).sum()) > 0
    )

    remaining_rows = max_display_rows
    allocations: dict[str, int] = {}
    for index, direction in enumerate(directions):
        direction_count = int((profile_metadata["inc_dec"] == direction).sum())
        directions_left = len(directions) - index
        if directions_left == 1:
            allocation = remaining_rows
        else:
            proportional = round(max_display_rows * direction_count / len(profile_metadata))
            allocation = max(1, min(proportional, remaining_rows - directions_left + 1))
        allocation = min(allocation, direction_count)
        allocations[direction] = allocation
        remaining_rows -= allocation

    return allocations


def build_rendered_display_matrix(
    plot_matrix: pd.DataFrame,
    profile_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    max_display_rows: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate adjacent clustered profiles into display bins for a readable PNG."""
    if max_display_rows <= 0 or len(plot_matrix) <= max_display_rows:
        display_matrix = build_display_matrix(plot_matrix, profile_metadata, trait_metadata)
        display_metadata = profile_metadata.copy()
        display_metadata.index = display_matrix.index
        display_metadata["profile_count"] = 1
        display_metadata["study_count_in_bin"] = 1
        return display_matrix, display_metadata

    rows: list[pd.Series] = []
    metadata_rows: list[dict[str, object]] = []
    allocations = allocate_display_bins(profile_metadata, max_display_rows)

    for direction, bin_count in allocations.items():
        direction_index = profile_metadata.index[profile_metadata["inc_dec"] == direction]
        direction_index = direction_index.intersection(plot_matrix.index)
        if len(direction_index) == 0:
            continue

        for bin_number, bin_profile_ids in enumerate(
            np.array_split(direction_index.to_numpy(), bin_count),
            start=1,
        ):
            if len(bin_profile_ids) == 0:
                continue

            bin_index = pd.Index(bin_profile_ids)
            bin_values = plot_matrix.loc[bin_index].mean(axis=0)
            bin_metadata = profile_metadata.loc[bin_index]
            direction_label = DIRECTION_LABELS.get(direction, direction)
            label = (
                f"{direction_label} profile bin {bin_number:03d} | "
                f"n={len(bin_index)} | studies={bin_metadata['study_id'].nunique()}"
            )
            bin_values.name = label
            rows.append(bin_values)
            metadata_rows.append(
                {
                    "display_label": label,
                    "inc_dec": direction,
                    "direction_label": direction_label,
                    "profile_count": len(bin_index),
                    "study_count_in_bin": bin_metadata["study_id"].nunique(),
                    "body_site_group": most_common_series_value(bin_metadata["body_site_group"]),
                    "condition_group": most_common_series_value(bin_metadata["condition_group"]),
                    "first_profile_id": bin_index[0],
                    "last_profile_id": bin_index[-1],
                }
            )

    display_matrix = pd.DataFrame(rows)
    display_matrix.columns = [
        trait_metadata.loc[trait_id, "display_label"] or trait_id
        for trait_id in display_matrix.columns
    ]
    display_metadata = pd.DataFrame.from_records(metadata_rows).set_index("display_label")
    return display_matrix, display_metadata


def build_display_matrix(
    plot_matrix: pd.DataFrame,
    profile_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Apply readable labels for the plotted matrix."""
    display = plot_matrix.copy()
    display.index = [
        " | ".join(
            [
                profile_metadata.loc[profile_id, "study_id"],
                profile_metadata.loc[profile_id, "experiment_id"],
                profile_metadata.loc[profile_id, "direction_label"],
            ]
        )
        for profile_id in display.index
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
    display_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    output_path: Path,
    full_profile_count: int,
    full_study_count: int,
) -> None:
    """Render and save the clustered study profile heatmap."""
    sns.set_theme(style="white")
    display_to_trait_id = dict(zip(display_matrix.columns, trait_metadata.index, strict=True))

    direction_palette = {"inc": "#c0392b", "dec": "#2874a6"}
    body_palette = palette_for_categories(display_metadata["body_site_group"], "tab20")
    condition_palette = palette_for_categories(display_metadata["condition_group"], "Set3")
    trait_palette = palette_for_categories(trait_metadata["display_category"], "Set2")

    row_colors = pd.DataFrame(
        {
            "direction": display_metadata.loc[display_matrix.index, "inc_dec"]
            .map(direction_palette)
            .fillna("#7f8c8d")
            .to_numpy(),
            "body_site": display_metadata.loc[display_matrix.index, "body_site_group"]
            .map(body_palette)
            .to_numpy(),
            "condition": display_metadata.loc[display_matrix.index, "condition_group"]
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
        vmin=-1,
        vmax=1,
        method="average",
        metric="correlation",
        row_cluster=False,
        col_cluster=True,
        row_colors=row_colors,
        col_colors=col_colors,
        linewidths=0.05,
        linecolor="#e5e7eb",
        xticklabels=True,
        yticklabels=len(display_matrix) <= 120,
        figsize=(24, min(40, max(18, len(display_matrix) / 8))),
        cbar_pos=(0.02, 0.80, 0.02, 0.14),
        dendrogram_ratio=(0.12, 0.10),
        colors_ratio=(0.03, 0.03),
        cbar_kws={
            "label": (
                "Column z-score of log1p trait fraction per 1,000 profile counts; "
                "display capped at +/-1"
            )
        },
    )
    cluster.fig.suptitle(
        "Study Trait Profile Clustering\n"
        f"{full_profile_count} BugSigDB study/experiment/direction profiles from "
        f"{full_study_count} studies, aggregated into {len(display_matrix)} "
        "clustered display bins; "
        "increased and decreased sections are separated",
        y=1.08,
        fontsize=16,
        fontweight="bold",
    )
    cluster.ax_heatmap.set_xlabel("Trait")
    cluster.ax_heatmap.set_ylabel("Clustered profile display bin")
    increased_count = int((display_metadata["inc_dec"] == "inc").sum())
    if 0 < increased_count < len(display_matrix):
        cluster.ax_heatmap.axhline(increased_count, color="#111827", linewidth=0.8)
    for label in cluster.ax_heatmap.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(7)
    for label in cluster.ax_heatmap.get_yticklabels():
        label.set_fontsize(7)

    cluster.ax_row_dendrogram.legend(
        handles=legend_handles(direction_palette, display_metadata["inc_dec"]),
        title="Direction",
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
        handles=legend_handles(body_palette, display_metadata["body_site_group"]),
        title="Body site",
        loc="lower left",
        bbox_to_anchor=(0.02, 0.18),
        ncol=1,
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    cluster.fig.legend(
        handles=legend_handles(condition_palette, display_metadata["condition_group"]),
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
    profile_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    rendered_rows: int,
    output_path: Path,
    matrix_output_path: Path,
    profile_metadata_output_path: Path,
    trait_metadata_output_path: Path,
) -> None:
    """Write the plotted matrix and source metadata."""
    matrix_output_path.parent.mkdir(parents=True, exist_ok=True)
    display_matrix.to_csv(matrix_output_path, sep="\t")
    profile_ids = plot_matrix.index
    trait_ids = plot_matrix.columns
    profile_metadata.loc[profile_ids].to_csv(profile_metadata_output_path, sep="\t")
    trait_metadata.loc[trait_ids].to_csv(trait_metadata_output_path, sep="\t")
    print(f"Wrote heatmap: {output_path}")
    print(f"Wrote matrix: {matrix_output_path}")
    print(f"Wrote profile metadata: {profile_metadata_output_path}")
    print(f"Wrote trait metadata: {trait_metadata_output_path}")
    print("Measurement: column z-score of log1p trait fraction per 1,000 profile counts")
    print(f"Profiles analyzed and written: {len(profile_ids)}")
    print(f"Rendered heatmap display rows: {rendered_rows}")
    print(f"Studies represented: {profile_metadata.loc[profile_ids, 'study_id'].nunique()}")
    print(
        "Increased/decreased profiles: "
        f"{int((profile_metadata.loc[profile_ids, 'inc_dec'] == 'inc').sum())}/"
        f"{int((profile_metadata.loc[profile_ids, 'inc_dec'] == 'dec').sum())}"
    )
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
        help="Output selected profile metadata TSV path.",
    )
    parser.add_argument(
        "--trait-metadata-output",
        type=Path,
        default=OUTPUT_DIR / "study_trait_profile_trait_metadata.tsv",
        help="Output selected trait metadata TSV path.",
    )
    parser.add_argument(
        "--max-profiles",
        type=int,
        default=0,
        help="Maximum number of profile rows to plot. Use 0 to include every profile.",
    )
    parser.add_argument(
        "--max-studies",
        type=int,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--max-traits", type=int, default=60)
    parser.add_argument(
        "--display-max-rows",
        type=int,
        default=320,
        help=(
            "Maximum rows to render in the PNG. When fewer than all profiles, adjacent "
            "clustered profiles are averaged into direction-specific display bins."
        ),
    )
    parser.add_argument("--min-traits", type=int, default=0)
    parser.add_argument("--min-total-signal", type=int, default=0)
    parser.add_argument(
        "--min-profile-support",
        "--min-study-support",
        dest="min_profile_support",
        type=int,
        default=6,
    )
    parser.add_argument(
        "--max-profile-support-fraction",
        "--max-study-support-fraction",
        dest="max_profile_support_fraction",
        type=float,
        default=0.95,
    )
    parser.add_argument("--max-body-site-categories", type=int, default=8)
    parser.add_argument("--max-condition-categories", type=int, default=8)
    parser.add_argument(
        "--include-media-traits",
        action="store_true",
        help="Include MediaDive growth-medium traits in the selected traits.",
    )
    args = parser.parse_args()

    profile_matrix, profile_metadata, _trait_columns = load_profile_matrix(
        args.study_trait_profiles_input
    )
    profile_metadata = profile_metadata.join(
        load_study_titles(args.studies_input),
        on="study_id",
        how="left",
    )
    profile_metadata = profile_metadata.join(
        load_raw_study_metadata(args.full_dump_input),
        on="study_id",
        how="left",
    )
    for column in [
        "title",
        "publication_year",
        "condition",
        "body_site",
        "location",
        "host_species",
    ]:
        profile_metadata[column] = profile_metadata[column].map(clean_value)

    max_profiles = args.max_profiles
    if args.max_studies is not None and max_profiles == 0:
        max_profiles = args.max_studies
    selected_profiles = select_profiles(
        profile_matrix,
        profile_metadata,
        max_profiles,
        args.min_traits,
        args.min_total_signal,
    )
    selected_profile_matrix = profile_matrix.loc[selected_profiles]
    profile_metadata = profile_metadata.loc[selected_profiles].copy()
    profile_metadata["body_site_group"] = group_top_categories(
        profile_metadata["body_site"],
        args.max_body_site_categories,
    )
    profile_metadata["condition_group"] = group_top_categories(
        profile_metadata["condition"],
        args.max_condition_categories,
    )

    trait_labels = load_trait_labels(args.trait_columns_input)
    selected_traits, trait_metadata = select_traits(
        selected_profile_matrix,
        trait_labels,
        args.max_traits,
        args.min_profile_support,
        args.max_profile_support_fraction,
        args.include_media_traits,
    )
    plot_matrix = build_plot_matrix(selected_profile_matrix, selected_traits)
    ordered_profiles = order_profiles_by_direction(plot_matrix, profile_metadata)
    plot_matrix = plot_matrix.loc[ordered_profiles]
    profile_metadata = profile_metadata.loc[ordered_profiles]
    display_matrix = build_display_matrix(plot_matrix, profile_metadata, trait_metadata)
    rendered_matrix, rendered_metadata = build_rendered_display_matrix(
        plot_matrix,
        profile_metadata,
        trait_metadata,
        args.display_max_rows,
    )
    save_heatmap(
        rendered_matrix,
        rendered_metadata,
        trait_metadata,
        args.output,
        full_profile_count=len(plot_matrix),
        full_study_count=profile_metadata["study_id"].nunique(),
    )
    write_outputs(
        plot_matrix,
        display_matrix,
        profile_metadata,
        trait_metadata,
        len(rendered_matrix),
        args.output,
        args.matrix_output,
        args.study_metadata_output,
        args.trait_metadata_output,
    )


if __name__ == "__main__":
    main()
