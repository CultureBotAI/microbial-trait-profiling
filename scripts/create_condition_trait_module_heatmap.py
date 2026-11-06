"""Create clustered heatmaps for condition/body-site trait modules."""

from __future__ import annotations

import argparse
import math
from collections import Counter
from pathlib import Path
from textwrap import fill

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

from mtp.composition_denominators import (
    DENOMINATOR_DEFINITION,
    retained_trait_counts,
    trait_fractions,
)
from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "trait-visualization"
PROFILE_METADATA_COLUMNS = ["study_id", "experiment_id", "inc_dec"]
HUMAN_HOST_SPECIES = "Homo sapiens"
GUILD_DEFINITION_TABLE_OUTPUT = OUTPUT_DIR / "condition_trait_module_guild_definition_table.png"
GUILD_DEFINITION_TABLE_SVG_OUTPUT = OUTPUT_DIR / "condition_trait_module_guild_definition_table.svg"

SLATE_50 = "#f8fafc"
SLATE_100 = "#f1f5f9"
SLATE_200 = "#e2e8f0"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"
BLUE_50 = "#eff6ff"
BLUE_600 = "#2563eb"
GREEN_50 = "#ecfdf5"
GREEN_700 = "#047857"
AMBER_50 = "#fffbeb"
AMBER_700 = "#b45309"


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


def wrapped(value: str, width: int) -> str:
    """Wrap text for compact diagram labels."""
    return "\n".join(fill(line, width=width) for line in value.split("\n"))


def most_common_label(values: list[str]) -> str:
    """Return the most common non-empty label with deterministic ties."""
    labels = [clean_value(value) for value in values if clean_value(value) != "unknown"]
    if not labels:
        return "unknown"
    counts = Counter(labels)
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]


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


def load_aggregated_study_profiles(
    profile_path: Path,
    human_experiment_ids: set[str] | None = None,
    direction: str = "all",
    experiment_metadata: pd.DataFrame | None = None,
    group_by: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate experiment profiles to study-level trait profiles.

    By default (`direction="all"`), increased and decreased profiles within each
    study are summed to measure direction-agnostic study trait load. Pass
    `direction="inc"` or `direction="dec"` to restrict aggregation to a single
    direction (#134).

    With `experiment_metadata` (labels indexed by study and experiment), each profile
    keeps its own experiment's labels and profiles are summed per study and `group_by`
    label combination. A study whose experiments differ in condition or body site then
    contributes one row to each of those modules instead of pooling every experiment
    under the study's majority label (#555).
    """
    if direction not in {"all", "inc", "dec"}:
        msg = f"Unsupported direction filter: {direction!r}; expected 'all', 'inc', or 'dec'."
        raise ValueError(msg)
    if experiment_metadata is not None and not group_by:
        # Without grouping fields every study collapses to one row, which hands each
        # experiment its study's majority label again (#555, #646).
        msg = "experiment_metadata needs group_by, the experiment labels to sum profiles by."
        raise ValueError(msg)
    profiles = pd.read_csv(profile_path, sep="\t")
    unfiltered_profile_count = len(profiles)
    unfiltered_study_count = profiles["study_id"].nunique()
    if human_experiment_ids is not None:
        profiles = profiles[profiles["experiment_id"].isin(human_experiment_ids)].copy()
        if profiles.empty:
            msg = "No Homo sapiens profiles were available after applying the human-host filter."
            raise ValueError(msg)
        print(
            "Human-host corpus filter retained "
            f"{len(profiles):,}/{unfiltered_profile_count:,} profiles from "
            f"{profiles['study_id'].nunique():,}/{unfiltered_study_count:,} studies."
        )
    if direction != "all":
        profiles = profiles[profiles["inc_dec"] == direction].copy()
        if profiles.empty:
            msg = f"No profiles remained after filtering for direction={direction!r}."
            raise ValueError(msg)

    trait_columns = [
        column for column in profiles.columns if column not in PROFILE_METADATA_COLUMNS
    ]
    values = retained_trait_counts(
        profiles[trait_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    )
    aggregations: dict[str, tuple[str, object]] = {
        "study_id": ("study_id", "first"),
        "profile_count": ("experiment_id", "size"),
        "experiment_count": ("experiment_id", "nunique"),
    }
    if experiment_metadata is None:
        row_ids = profiles["study_id"]
    else:
        label_columns = experiment_metadata.columns.tolist()
        profiles = profiles.join(experiment_metadata, on=["study_id", "experiment_id"], how="left")
        for column in label_columns:
            profiles[column] = profiles[column].map(clean_value)
            aggregations[column] = (column, most_common_label)
        row_ids = profiles[["study_id", *(group_by or [])]].agg(" | ".join, axis=1)
    # study_id stays a column; an index of the same name breaks reset_index() and
    # join(on="study_id").
    row_ids = row_ids.rename("row_id")
    study_matrix = values.groupby(row_ids, sort=False).sum()

    profile_metadata = profiles.groupby(row_ids, sort=False).agg(**aggregations)
    direction_counts = pd.crosstab(row_ids, profiles["inc_dec"])
    for dir_label in ["inc", "dec"]:
        if dir_label not in direction_counts.columns:
            direction_counts[dir_label] = 0
    profile_metadata = profile_metadata.join(direction_counts[["inc", "dec"]], how="left").fillna(0)
    profile_metadata = profile_metadata.rename(
        columns={"inc": "increased_profile_count", "dec": "decreased_profile_count"}
    )
    profile_metadata["total_trait_signal"] = study_matrix.sum(axis=1)
    profile_metadata["trait_count"] = study_matrix.gt(0).sum(axis=1)
    profile_metadata["denominator_definition"] = DENOMINATOR_DEFINITION
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
    study_matrix = retained_trait_counts(study_matrix)
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
        study_count=("study_id", "nunique"),
        profile_count=("profile_count", "sum"),
        experiment_count=("experiment_count", "sum"),
        condition=("condition", most_common_label),
        body_site=("body_site", most_common_label),
        location=("location", most_common_label),
        host_species=("host_species", most_common_label),
    )
    group_metadata["total_trait_signal"] = group_matrix.sum(axis=1)
    group_metadata["trait_count"] = group_matrix.gt(0).sum(axis=1)
    group_metadata["denominator_definition"] = DENOMINATOR_DEFINITION

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
    if include_media_traits:
        raise ValueError("MediaDive cannot be included in harmonized composition analyses")
    group_matrix = retained_trait_counts(group_matrix)
    composition = trait_fractions(group_matrix)
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
    composition = trait_fractions(group_matrix).loc[:, selected_traits]
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
        "Condition/Disease Trait Modules\n"
        "BugSigDB Condition Groups by Aggregate Study Trait Load (Increased + Decreased Combined)",
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


def summarize_trait_categories(trait_metadata: pd.DataFrame, max_categories: int = 4) -> str:
    """Return a compact category summary for the selected trait set."""
    counts = trait_metadata["display_category"].map(clean_value).value_counts()
    parts = [f"{category}: {count}" for category, count in counts.head(max_categories).items()]
    remaining = len(counts) - len(parts)
    if remaining > 0:
        parts.append(f"+{remaining} more")
    return "; ".join(parts)


def add_table_cell(
    ax: plt.Axes,
    *,
    x: float,
    y: float,
    width: float,
    height: float,
    text: str,
    facecolor: str,
    edgecolor: str,
    color: str = SLATE_700,
    fontsize: float = 7.2,
    fontweight: str = "normal",
    ha: str = "left",
    wrap_width: int = 34,
) -> None:
    rect = Rectangle(
        (x, y),
        width,
        height,
        linewidth=0.65,
        edgecolor=edgecolor,
        facecolor=facecolor,
        zorder=2,
    )
    ax.add_patch(rect)
    text_x = x + width / 2 if ha == "center" else x + 0.75
    ax.text(
        text_x,
        y + height / 2,
        wrapped(text, wrap_width),
        ha=ha,
        va="center",
        fontsize=fontsize,
        fontweight=fontweight,
        color=color,
        linespacing=1.1,
        zorder=3,
    )


def save_guild_definition_table(
    group_metadata: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    output_path: Path,
    svg_output_path: Path | None,
) -> None:
    """Render a diagram table that defines guild logic for condition modules."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = [
        {
            "layer": "Taxon support",
            "definition": "A taxon supports a guild when it carries X OR Y OR Z.",
            "use": "Candidate universe and supporting-taxa evidence.",
            "logic": "OR",
        },
        {
            "layer": "Profile score",
            "definition": "Score each profile as mean z(X, Y, Z) over selected trait columns.",
            "use": "Continuous trait-enrichment signal for profile and direction summaries.",
            "logic": "mean",
        },
        {
            "layer": "Disease module",
            "definition": (
                "Aggregate profiles by condition/body site, then summarize trait-enrichment score."
            ),
            "use": "Condition-level link to a trait-defined guild.",
            "logic": "group mean",
        },
        {
            "layer": "Strict all-trait rule",
            "definition": "Use X AND Y AND Z only when every trait is a required criterion.",
            "use": "Not the default for broad KG-Microbe phenotype guilds.",
            "logic": "avoid AND",
        },
    ]

    fig, ax = plt.subplots(figsize=(15.2, 5.25), dpi=180)
    fig.patch.set_facecolor(SLATE_50)
    ax.set_facecolor(SLATE_50)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 38)
    ax.axis("off")

    ax.text(
        3.5,
        35.2,
        "Condition/disease module guild definition logic",
        fontsize=18,
        fontweight="bold",
        color=SLATE_900,
        ha="left",
        va="top",
    )
    ax.text(
        3.5,
        32.7,
        (
            f"Current context: {len(group_metadata):,} condition/body-site modules x "
            f"{len(trait_metadata):,} selected trait columns."
        ),
        fontsize=8.9,
        color=SLATE_500,
        ha="left",
        va="top",
    )

    badge = FancyBboxPatch(
        (78.5, 31.4),
        17.0,
        2.5,
        boxstyle="round,pad=0.25,rounding_size=1.2",
        linewidth=1.0,
        edgecolor="#a7f3d0",
        facecolor=GREEN_50,
        zorder=3,
    )
    ax.add_patch(badge)
    ax.text(
        87.0,
        32.65,
        "use OR for support",
        ha="center",
        va="center",
        fontsize=7.6,
        fontweight="bold",
        color=GREEN_700,
        zorder=4,
    )

    table_x = 3.5
    table_y = 6.0
    table_w = 93.0
    header_h = 3.0
    row_h = 5.3
    table_h = header_h + row_h * len(rows)
    columns = [
        ("Layer", 17.0, 16),
        ("Definition", 35.0, 44),
        ("Use in diagram", 29.0, 37),
        ("Logic", 12.0, 12),
    ]

    outer = FancyBboxPatch(
        (table_x, table_y),
        table_w,
        table_h,
        boxstyle="round,pad=0.35,rounding_size=0.9",
        linewidth=1.1,
        edgecolor=SLATE_200,
        facecolor="white",
        zorder=1,
    )
    ax.add_patch(outer)

    x = table_x
    for header, width, wrap_width in columns:
        add_table_cell(
            ax,
            x=x,
            y=table_y + row_h * len(rows),
            width=width,
            height=header_h,
            text=header,
            facecolor=BLUE_50,
            edgecolor=SLATE_200,
            color=BLUE_600,
            fontsize=7.4,
            fontweight="bold",
            ha="center" if header == "Logic" else "left",
            wrap_width=wrap_width,
        )
        x += width

    for row_index, row in enumerate(rows):
        y = table_y + row_h * (len(rows) - 1 - row_index)
        facecolor = "white" if row_index % 2 == 0 else SLATE_50
        logic_color = GREEN_700 if row["logic"] == "OR" else AMBER_700
        if row["logic"] == "avoid AND":
            logic_color = SLATE_500
        x = table_x
        for key, width, wrap_width in [
            ("layer", columns[0][1], columns[0][2]),
            ("definition", columns[1][1], columns[1][2]),
            ("use", columns[2][1], columns[2][2]),
            ("logic", columns[3][1], columns[3][2]),
        ]:
            add_table_cell(
                ax,
                x=x,
                y=y,
                width=width,
                height=row_h,
                text=row[key],
                facecolor=facecolor,
                edgecolor=SLATE_200,
                color=logic_color if key == "logic" else SLATE_700,
                fontsize=7.2,
                fontweight="bold" if key in {"layer", "logic"} else "normal",
                ha="center" if key == "logic" else "left",
                wrap_width=wrap_width,
            )
            x += width

    ax.text(
        3.5,
        2.6,
        (
            "Use OR because KG-Microbe annotations are partially overlapping and incomplete. "
            "Requiring every selected trait would measure annotation completeness, "
            "not guild support."
        ),
        fontsize=8.2,
        color=SLATE_500,
        ha="left",
        va="bottom",
    )

    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor=fig.get_facecolor())
    if svg_output_path is not None:
        fig.savefig(svg_output_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


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


def build_parser(
    parser_class: type[argparse.ArgumentParser] = argparse.ArgumentParser,
) -> argparse.ArgumentParser:
    """Build the execution parser shared with arm-plan verification (#973)."""
    parser = parser_class(
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
        "--guild-definition-table-output",
        type=Path,
        default=GUILD_DEFINITION_TABLE_OUTPUT,
        help="Output guild definition table diagram PNG path.",
    )
    parser.add_argument(
        "--guild-definition-table-svg-output",
        type=Path,
        default=GUILD_DEFINITION_TABLE_SVG_OUTPUT,
        help="Output guild definition table diagram SVG path.",
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
        "--direction",
        choices=("all", "inc", "dec"),
        default="all",
        help=(
            "Profile direction filter before study aggregation: "
            "all (direction-agnostic study trait load), inc, or dec (#134)."
        ),
    )
    parser.add_argument(
        "--include-media-traits",
        action="store_true",
        help="Rejected: MediaDive growth-medium traits are excluded before normalization.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.include_media_traits:
        parser.error("MediaDive is excluded before composition normalization (#174).")

    # This sibling-script import belongs to execution. The shared parser must be
    # loadable without a script directory on sys.path or any scientific data access.
    from create_study_trait_profile_heatmap import load_raw_experiment_metadata

    group_by = normalize_group_fields(args.group_by)
    experiment_metadata = load_raw_experiment_metadata(args.full_dump_input)
    # A non-default profile input is held to the 16S and shotgun rule too (#742, #754)
    human_experiment_ids = set(
        experiment_metadata.index.get_level_values("experiment_id")[
            (experiment_metadata["host_species"].map(clean_value) == HUMAN_HOST_SPECIES)
            & experiment_metadata["sequencing_type"].map(clean_value).isin(CORPUS_SEQUENCING_TYPES)
        ]
    )
    study_matrix, profile_metadata = load_aggregated_study_profiles(
        args.study_trait_profiles_input,
        human_experiment_ids,
        direction=args.direction,
        experiment_metadata=experiment_metadata,
        group_by=group_by,
    )
    study_metadata = profile_metadata.join(
        load_study_titles(args.studies_input),
        on="study_id",
        how="left",
    )
    study_metadata["host_species"] = HUMAN_HOST_SPECIES

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
    save_guild_definition_table(
        group_metadata.loc[plot_matrix.index],
        trait_metadata.loc[plot_matrix.columns],
        args.guild_definition_table_output,
        args.guild_definition_table_svg_output,
    )
    print(f"Wrote guild definition table: {args.guild_definition_table_output}")
    print(f"Wrote guild definition table SVG: {args.guild_definition_table_svg_output}")
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
