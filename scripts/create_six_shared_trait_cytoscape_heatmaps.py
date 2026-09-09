"""Create inferential heatmaps from the six-shared-trait Cytoscape inputs."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "trait-visualization" / "trait-profile-investigation" / "six-shared-traits"

GRAPH_MANIFEST_FILENAME = "six_shared_trait_graph_manifest.tsv"
TRAIT_SELECTION_FILENAME = "six_shared_trait_selection.tsv"
MANIFEST_FILENAME = "six_shared_trait_cytoscape_heatmap_manifest.tsv"
PATTERN_TABLE_FILENAME = "six_shared_trait_condition_body_site_trait_pattern_matrix.tsv"
CONTRAST_TABLE_FILENAME = "six_shared_trait_condition_body_site_direction_contrasts.tsv"

MIN_POSITIVE_PROFILES = 3
MIN_STUDIES = 2
P_VALUE_CUTOFF = 0.05
MAX_PATTERN_ROWS = 24
MAX_CONTRAST_ROWS = 24

SLATE_50 = "#f8fafc"
SLATE_100 = "#f1f5f9"
SLATE_200 = "#e2e8f0"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"

LEGACY_CORRELATION_OUTPUTS = (
    "six_shared_trait_condition_body_site_trait_correlations.tsv",
    "six_shared_trait_condition_body_site_trait_vectors.tsv",
    "six_shared_trait_correlation_heatmap_manifest.tsv",
    "six_shared_traits_increased_grouped_condition_body_site_trait_correlation_heatmap.png",
    "six_shared_traits_decreased_grouped_condition_body_site_trait_correlation_heatmap.png",
    (
        "six_shared_traits_increased_individual_conditions_condition_body_site_trait_"
        "correlation_heatmap.png"
    ),
    (
        "six_shared_traits_decreased_individual_conditions_condition_body_site_trait_"
        "correlation_heatmap.png"
    ),
)


@dataclass(frozen=True)
class ContrastMode:
    """One condition-label mode for signed increased-vs-decreased contrasts."""

    key: str
    label: str
    condition_mode: str
    increased_view_key: str
    decreased_view_key: str


CONTRAST_MODES = (
    ContrastMode(
        key="grouped",
        label="Grouped conditions",
        condition_mode="grouped condition labels",
        increased_view_key="increased_grouped",
        decreased_view_key="decreased_grouped",
    ),
    ContrastMode(
        key="individual_conditions",
        label="Individual conditions",
        condition_mode="raw individual condition labels",
        increased_view_key="increased_individual_conditions",
        decreased_view_key="decreased_individual_conditions",
    ),
)


def read_tsv(path: Path) -> pd.DataFrame:
    """Read a TSV with local defaults."""
    return pd.read_csv(path, sep="\t")


def relative_path(path: Path, repo_root: Path = REPO_ROOT) -> str:
    """Return a stable repo-relative path."""
    return path.relative_to(repo_root).as_posix()


def resolve_repo_path(path_text: str, repo_root: Path = REPO_ROOT) -> Path:
    """Resolve a manifest path that may be repo-relative or absolute."""
    path = Path(path_text)
    return path if path.is_absolute() else repo_root / path


def wrap_label(value: object, width: int = 24) -> str:
    """Wrap labels for compact heatmap axes."""
    text = str(value)
    return "\n".join(wrap(text, width=width, break_long_words=False)) or text


def format_p_value(value: float) -> str:
    """Format p-values for cell annotations."""
    if not math.isfinite(value):
        return "NA"
    if value < 0.001:
        return "<0.001"
    return f"{value:.3f}"


def clean_legacy_outputs(output_dir: Path = OUTPUT_DIR) -> None:
    """Remove the superseded correlation heatmap outputs from the prior workflow."""
    for filename in LEGACY_CORRELATION_OUTPUTS:
        path = output_dir / filename
        if path.exists():
            path.unlink()


def read_trait_order(output_dir: Path = OUTPUT_DIR) -> pd.DataFrame:
    """Read the six shared traits in display order."""
    traits = read_tsv(output_dir / TRAIT_SELECTION_FILENAME)
    traits = traits.rename(columns={"trait_name": "trait_label"})
    traits["trait_rank"] = np.arange(len(traits), dtype=int)
    return traits[["trait_id", "trait_label", "slug", "base_color", "trait_rank"]].copy()


def coerce_display_edge_numbers(edges: pd.DataFrame) -> pd.DataFrame:
    """Coerce display-edge numeric columns used by both heatmap types."""
    output = edges.copy()
    numeric_columns = [
        "trait_rank",
        "profile_count",
        "study_count",
        "focus_profile_count",
        "focus_study_count",
        "focus_enriched_profile_count",
        "focus_enriched_study_count",
        "focus_no_enriched_trait_profile_count",
        "comparator_profile_count",
        "comparator_study_count",
        "comparator_enriched_profile_count",
        "background_no_enriched_profile_count",
        "focus_enriched_mean_guild_z",
        "focus_enriched_minus_comparator_guild_z",
        "rank_score",
    ]
    for column in numeric_columns:
        if column in output:
            output[column] = pd.to_numeric(output[column], errors="coerce").fillna(0.0)
    return output


def signal_value(row: pd.Series) -> float:
    """Return the nonnegative signal used by the Cytoscape color scale."""
    contrast = row.get("focus_enriched_minus_comparator_guild_z", math.nan)
    mean = row.get("focus_enriched_mean_guild_z", math.nan)
    value = contrast if pd.notna(contrast) else mean
    if pd.isna(value):
        return 0.0
    return max(0.0, float(value))


def focus_vs_comparator_p_value(row: pd.Series) -> float:
    """Test whether focus positive-z support exceeds comparator positive-z support."""
    focus_positive = int(row["focus_enriched_profile_count"])
    focus_negative = int(row["focus_no_enriched_trait_profile_count"])
    comparator_positive = int(row["comparator_enriched_profile_count"])
    comparator_negative = max(0, int(row["comparator_profile_count"]) - comparator_positive)
    if focus_positive + focus_negative == 0 or comparator_positive + comparator_negative == 0:
        return math.nan
    result = fisher_exact(
        [[focus_positive, focus_negative], [comparator_positive, comparator_negative]],
        alternative="greater",
    )
    return float(result.pvalue)


def signed_direction_p_value(
    increased_positive: int,
    increased_total: int,
    decreased_positive: int,
    decreased_total: int,
) -> float:
    """Test whether positive-z support differs between increased and decreased profiles."""
    if increased_total <= 0 or decreased_total <= 0:
        return math.nan
    result = fisher_exact(
        [
            [increased_positive, max(0, increased_total - increased_positive)],
            [decreased_positive, max(0, decreased_total - decreased_positive)],
        ],
        alternative="two-sided",
    )
    return float(result.pvalue)


def pair_label(condition: object, body_site_group: object) -> str:
    """Return the row label used by condition/body-site heatmaps."""
    return f"{condition} | {body_site_group}"


def prepare_pattern_cells(
    display_edges: pd.DataFrame,
    view: pd.Series,
    trait_order: pd.DataFrame,
) -> pd.DataFrame:
    """Build long-form cells for one condition/body-site x trait pattern heatmap."""
    trait_rank = dict(zip(trait_order["trait_id"], trait_order["trait_rank"], strict=True))
    edges = coerce_display_edge_numbers(display_edges)
    records: list[dict[str, object]] = []
    for row in edges.itertuples(index=False):
        row_series = pd.Series(row._asdict())
        p_value = focus_vs_comparator_p_value(row_series)
        signal = signal_value(row_series)
        strength = math.log1p(float(row.focus_enriched_profile_count)) * signal
        supported = (
            int(row.focus_enriched_profile_count) >= MIN_POSITIVE_PROFILES
            and int(row.focus_enriched_study_count) >= MIN_STUDIES
        )
        significant = supported and math.isfinite(p_value) and p_value <= P_VALUE_CUTOFF
        reason = ""
        if not supported:
            reason = "low positive-z profile or study support"
        elif not significant:
            reason = f"Fisher p>{P_VALUE_CUTOFF:g}"
        records.append(
            {
                "view_key": view.view_key,
                "view_label": view.view_label,
                "focus_direction": view.focus_direction,
                "condition_mode": view.condition_mode,
                "condition": row.condition,
                "body_site_group": row.body_site_group,
                "condition_body_site_pair": pair_label(row.condition, row.body_site_group),
                "trait_id": row.trait_id,
                "trait_label": row.trait_label,
                "trait_rank": trait_rank.get(row.trait_id, 999),
                "focus_positive_profile_count": int(row.focus_enriched_profile_count),
                "focus_positive_study_count": int(row.focus_enriched_study_count),
                "focus_profile_count": int(row.focus_profile_count),
                "comparator_positive_profile_count": int(row.comparator_enriched_profile_count),
                "comparator_profile_count": int(row.comparator_profile_count),
                "signal": signal,
                "strength": strength,
                "fisher_p_value": p_value,
                "displayed_in_heatmap": significant,
                "blank_reason": reason,
            }
        )
    return pd.DataFrame.from_records(records)


def select_pattern_rows(cells: pd.DataFrame, *, max_rows: int = MAX_PATTERN_ROWS) -> list[str]:
    """Select condition/body-site rows with the strongest significant trait recurrence."""
    stats = (
        cells.groupby("condition_body_site_pair", sort=False)
        .agg(
            significant_trait_count=("displayed_in_heatmap", "sum"),
            total_positive_profiles=("focus_positive_profile_count", "sum"),
            max_strength=("strength", "max"),
            total_strength=("strength", "sum"),
        )
        .sort_values(
            [
                "significant_trait_count",
                "total_positive_profiles",
                "total_strength",
                "max_strength",
            ],
            ascending=[False, False, False, False],
        )
    )
    stats = stats[stats["significant_trait_count"] > 0]
    if stats.empty:
        stats = (
            cells.groupby("condition_body_site_pair", sort=False)
            .agg(
                significant_trait_count=("displayed_in_heatmap", "sum"),
                total_positive_profiles=("focus_positive_profile_count", "sum"),
                max_strength=("strength", "max"),
                total_strength=("strength", "sum"),
            )
            .sort_values(
                ["total_positive_profiles", "total_strength", "max_strength"],
                ascending=[False, False, False],
            )
        )
    return stats.head(max_rows).index.astype(str).tolist()


def prepare_contrast_cells(
    increased_edges: pd.DataFrame,
    decreased_edges: pd.DataFrame,
    mode: ContrastMode,
    trait_order: pd.DataFrame,
) -> pd.DataFrame:
    """Build signed increased-vs-decreased cells for one condition-label mode."""
    trait_rank = dict(zip(trait_order["trait_id"], trait_order["trait_rank"], strict=True))
    records: list[dict[str, object]] = []
    for edges in (increased_edges, decreased_edges):
        edges = coerce_display_edge_numbers(edges)
        for row in edges.itertuples(index=False):
            if row.focus_direction == "increased":
                increased_positive = int(row.focus_enriched_profile_count)
                increased_total = int(row.focus_profile_count)
                decreased_positive = int(row.comparator_enriched_profile_count)
                decreased_total = int(row.comparator_profile_count)
            else:
                decreased_positive = int(row.focus_enriched_profile_count)
                decreased_total = int(row.focus_profile_count)
                increased_positive = int(row.comparator_enriched_profile_count)
                increased_total = int(row.comparator_profile_count)
            records.append(
                {
                    "condition": row.condition,
                    "body_site_group": row.body_site_group,
                    "trait_id": row.trait_id,
                    "trait_label": row.trait_label,
                    "trait_rank": trait_rank.get(row.trait_id, 999),
                    "study_count": int(row.study_count),
                    "increased_positive_profile_count": increased_positive,
                    "increased_profile_count": increased_total,
                    "decreased_positive_profile_count": decreased_positive,
                    "decreased_profile_count": decreased_total,
                }
            )
    combined = (
        pd.DataFrame.from_records(records)
        .groupby(["condition", "body_site_group", "trait_id", "trait_label"], as_index=False)
        .agg(
            trait_rank=("trait_rank", "min"),
            study_count=("study_count", "max"),
            increased_positive_profile_count=("increased_positive_profile_count", "max"),
            increased_profile_count=("increased_profile_count", "max"),
            decreased_positive_profile_count=("decreased_positive_profile_count", "max"),
            decreased_profile_count=("decreased_profile_count", "max"),
        )
    )

    output_records: list[dict[str, object]] = []
    for row in combined.itertuples(index=False):
        increased_fraction = (
            row.increased_positive_profile_count / row.increased_profile_count
            if row.increased_profile_count
            else math.nan
        )
        decreased_fraction = (
            row.decreased_positive_profile_count / row.decreased_profile_count
            if row.decreased_profile_count
            else math.nan
        )
        delta = increased_fraction - decreased_fraction
        p_value = signed_direction_p_value(
            int(row.increased_positive_profile_count),
            int(row.increased_profile_count),
            int(row.decreased_positive_profile_count),
            int(row.decreased_profile_count),
        )
        supported = (
            max(
                int(row.increased_positive_profile_count),
                int(row.decreased_positive_profile_count),
            )
            >= MIN_POSITIVE_PROFILES
            and int(row.study_count) >= MIN_STUDIES
        )
        significant = supported and math.isfinite(p_value) and p_value <= P_VALUE_CUTOFF
        reason = ""
        if not supported:
            reason = "low positive-z profile or study support"
        elif not significant:
            reason = f"Fisher p>{P_VALUE_CUTOFF:g}"
        output_records.append(
            {
                "view_key": f"{mode.key}_direction_contrast",
                "view_label": f"{mode.label}: signed increased-vs-decreased contrast",
                "focus_direction": "increased_vs_decreased",
                "condition_mode": mode.condition_mode,
                "condition": row.condition,
                "body_site_group": row.body_site_group,
                "condition_body_site_pair": pair_label(row.condition, row.body_site_group),
                "trait_id": row.trait_id,
                "trait_label": row.trait_label,
                "trait_rank": int(row.trait_rank),
                "study_count": int(row.study_count),
                "increased_positive_profile_count": int(row.increased_positive_profile_count),
                "increased_profile_count": int(row.increased_profile_count),
                "decreased_positive_profile_count": int(row.decreased_positive_profile_count),
                "decreased_profile_count": int(row.decreased_profile_count),
                "increased_positive_fraction": increased_fraction,
                "decreased_positive_fraction": decreased_fraction,
                "increased_minus_decreased_positive_fraction": delta,
                "fisher_p_value": p_value,
                "displayed_in_heatmap": significant,
                "blank_reason": reason,
            }
        )
    return pd.DataFrame.from_records(output_records)


def select_contrast_rows(cells: pd.DataFrame, *, max_rows: int = MAX_CONTRAST_ROWS) -> list[str]:
    """Select condition/body-site rows with the strongest signed direction effects."""
    cells = cells.copy()
    cells["abs_delta"] = cells["increased_minus_decreased_positive_fraction"].abs()
    cells["positive_support"] = (
        cells["increased_positive_profile_count"] + cells["decreased_positive_profile_count"]
    )
    stats = (
        cells.groupby("condition_body_site_pair", sort=False)
        .agg(
            significant_trait_count=("displayed_in_heatmap", "sum"),
            positive_support=("positive_support", "sum"),
            max_abs_delta=("abs_delta", "max"),
            total_abs_delta=("abs_delta", "sum"),
        )
        .sort_values(
            ["significant_trait_count", "positive_support", "total_abs_delta", "max_abs_delta"],
            ascending=[False, False, False, False],
        )
    )
    stats = stats[stats["significant_trait_count"] > 0]
    if stats.empty:
        stats = (
            cells.groupby("condition_body_site_pair", sort=False)
            .agg(
                significant_trait_count=("displayed_in_heatmap", "sum"),
                positive_support=("positive_support", "sum"),
                max_abs_delta=("abs_delta", "max"),
                total_abs_delta=("abs_delta", "sum"),
            )
            .sort_values(
                ["positive_support", "total_abs_delta", "max_abs_delta"],
                ascending=[False, False, False],
            )
        )
    return stats.head(max_rows).index.astype(str).tolist()


def cell_matrix(
    cells: pd.DataFrame,
    row_labels: list[str],
    trait_order: pd.DataFrame,
    *,
    value_column: str,
) -> pd.DataFrame:
    """Pivot significant heatmap cells into a row x trait matrix."""
    trait_ids = trait_order["trait_id"].tolist()
    significant = cells[cells["displayed_in_heatmap"]].copy()
    matrix = significant.pivot_table(
        index="condition_body_site_pair",
        columns="trait_id",
        values=value_column,
        aggfunc="max",
    )
    return matrix.reindex(index=row_labels, columns=trait_ids)


def annotation_matrix(
    cells: pd.DataFrame,
    row_labels: list[str],
    trait_order: pd.DataFrame,
    *,
    mode: str,
) -> pd.DataFrame:
    """Return compact annotations for significant cells."""
    trait_ids = trait_order["trait_id"].tolist()
    annotations = pd.DataFrame("", index=row_labels, columns=trait_ids)
    significant = cells[cells["displayed_in_heatmap"]].copy()
    for row in significant.itertuples(index=False):
        if row.condition_body_site_pair not in annotations.index:
            continue
        if row.trait_id not in annotations.columns:
            continue
        if mode == "pattern":
            annotations.loc[row.condition_body_site_pair, row.trait_id] = (
                f"n={row.focus_positive_profile_count}\np={format_p_value(row.fisher_p_value)}"
            )
        else:
            delta = row.increased_minus_decreased_positive_fraction
            annotations.loc[row.condition_body_site_pair, row.trait_id] = (
                f"{delta:+.2f}\np={format_p_value(row.fisher_p_value)}"
            )
    return annotations


def plot_matrix_heatmap(
    matrix: pd.DataFrame,
    annotations: pd.DataFrame,
    trait_order: pd.DataFrame,
    output_path: Path,
    *,
    title: str,
    subtitle: str,
    colorbar_label: str,
    cmap_name: str,
    vmin: float | None = None,
    vmax: float | None = None,
) -> None:
    """Render one row x six-trait heatmap with blank nonsignificant cells."""
    row_count = max(1, len(matrix.index))
    fig_height = max(6.2, 2.8 + 0.36 * row_count)
    fig, ax = plt.subplots(figsize=(13.8, fig_height))
    fig.subplots_adjust(left=0.34, right=0.90, top=0.82, bottom=0.18)
    fig.patch.set_facecolor("white")
    ax.set_facecolor(SLATE_100)

    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(SLATE_100)
    values = matrix.to_numpy(dtype=float)
    image = ax.imshow(values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)

    ax.set_xticks(
        np.arange(len(trait_order)),
        [wrap_label(label, width=16) for label in trait_order["trait_label"]],
    )
    ax.set_yticks(
        np.arange(len(matrix.index)),
        [wrap_label(label, width=42) for label in matrix.index],
    )
    ax.tick_params(axis="x", labelrotation=35, labelsize=8.8)
    ax.tick_params(axis="y", labelsize=8.6)
    ax.set_xlabel("Six shared traits", color=SLATE_700)
    ax.set_ylabel("Condition | body-site pair", color=SLATE_700)

    ax.set_xticks(np.arange(-0.5, len(trait_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(matrix.index), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.3)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for y_index, row_label in enumerate(matrix.index):
        for x_index, trait_id in enumerate(trait_order["trait_id"]):
            text = annotations.loc[row_label, trait_id]
            if not text:
                continue
            value = matrix.loc[row_label, trait_id]
            text_color = SLATE_900
            if pd.notna(value):
                red, green, blue, _alpha = image.cmap(image.norm(float(value)))
                luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
                text_color = "white" if luminance < 0.56 else SLATE_900
            ax.text(
                x_index,
                y_index,
                text,
                ha="center",
                va="center",
                fontsize=7.2,
                color=text_color,
            )

    fig.text(
        0.04,
        0.965,
        title,
        ha="left",
        va="top",
        fontsize=15.5,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.04,
        0.925,
        "\n".join(wrap(subtitle, width=130)),
        ha="left",
        va="top",
        fontsize=9.2,
        color=SLATE_500,
    )
    cbar = fig.colorbar(image, ax=ax, shrink=0.78, pad=0.02)
    cbar.set_label(colorbar_label, color=SLATE_700)
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_edgecolor(SLATE_200)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_pattern_heatmap(
    cells: pd.DataFrame,
    view: pd.Series,
    trait_order: pd.DataFrame,
    output_dir: Path,
) -> tuple[Path, int, int, int]:
    """Write one approach-2 condition/body-site x trait pattern heatmap."""
    row_labels = select_pattern_rows(cells)
    matrix = cell_matrix(cells, row_labels, trait_order, value_column="strength")
    annotations = annotation_matrix(cells, row_labels, trait_order, mode="pattern")
    output_path = output_dir / (
        f"six_shared_traits_{view.view_key}_condition_body_site_trait_pattern_heatmap.png"
    )
    displayed_cells = int(cells["displayed_in_heatmap"].sum())
    subtitle = (
        "Approach 2: rows are condition/body-site pairs and columns are the six "
        "traits. Cell color is log1p(positive-z focus profiles) multiplied by the "
        "Cytoscape signal. Cells are shown only when Fisher p<=0.05, n>=3, and s>=2."
    )
    plot_matrix_heatmap(
        matrix,
        annotations,
        trait_order,
        output_path,
        title=f"Condition/body-site trait pattern: {view.view_label}",
        subtitle=subtitle,
        colorbar_label="Trait support strength",
        cmap_name="YlGnBu",
    )
    return output_path, len(row_labels), displayed_cells, len(cells)


def write_contrast_heatmap(
    cells: pd.DataFrame,
    mode: ContrastMode,
    trait_order: pd.DataFrame,
    output_dir: Path,
) -> tuple[Path, int, int, int]:
    """Write one approach-4 signed direction contrast heatmap."""
    row_labels = select_contrast_rows(cells)
    matrix = cell_matrix(
        cells,
        row_labels,
        trait_order,
        value_column="increased_minus_decreased_positive_fraction",
    )
    annotations = annotation_matrix(cells, row_labels, trait_order, mode="contrast")
    output_path = (
        output_dir
        / f"six_shared_traits_{mode.key}_condition_body_site_signed_direction_contrast_heatmap.png"
    )
    displayed_cells = int(cells["displayed_in_heatmap"].sum())
    subtitle = (
        "Approach 4: color is the increased-profile positive-z fraction minus the "
        "decreased-profile positive-z fraction. Red favors increased profiles; blue "
        "favors decreased profiles. Cells are shown only when Fisher p<=0.05, "
        "positive-z support is n>=3, and total study support is s>=2."
    )
    plot_matrix_heatmap(
        matrix,
        annotations,
        trait_order,
        output_path,
        title=f"Signed direction contrast: {mode.label}",
        subtitle=subtitle,
        colorbar_label="Increased - decreased positive-z fraction",
        cmap_name="RdBu_r",
        vmin=-1.0,
        vmax=1.0,
    )
    return output_path, len(row_labels), displayed_cells, len(cells)


def append_readme_section(output_dir: Path, manifest: pd.DataFrame) -> None:
    """Append or replace the Cytoscape heatmap README section."""
    readme_path = output_dir / "README.md"
    if not readme_path.exists():
        return
    new_marker = "## Condition-Body-Site Cytoscape Heatmaps"
    old_marker = "## Condition-Body-Site Trait Correlation Heatmaps"
    content = readme_path.read_text(encoding="utf-8")
    split_marker = old_marker if old_marker in content else new_marker
    prefix = content.split(split_marker, maxsplit=1)[0].rstrip()
    lines = [
        "",
        "",
        new_marker,
        "",
        f"- Manifest: `{MANIFEST_FILENAME}`",
        f"- Pattern matrix table: `{PATTERN_TABLE_FILENAME}`",
        f"- Signed contrast table: `{CONTRAST_TABLE_FILENAME}`",
        "",
        "These heatmaps reuse the same six-shared-trait display-edge TSVs that feed "
        "the Cytoscape graph views. They preserve condition/body-site/trait records "
        "instead of correlating already-collapsed condition and body-site nodes.",
        "",
        "- `trait_pattern`: Approach 2. Rows are condition/body-site pairs and "
        "columns are the six traits. Cell color is `log1p(n) * signal`, where `n` "
        "is positive-z focus-profile support and `signal` is the Cytoscape "
        "trait-specific signal, rendered with a yellow-to-blue sequential "
        "palette. Cells require Fisher p<=0.05, n>=3, and s>=2.",
        "- `signed_direction_contrast`: Approach 4. Rows are condition/body-site "
        "pairs and columns are traits. Cell color is the increased-profile "
        "positive-z fraction minus the decreased-profile positive-z fraction. Cells "
        "require Fisher p<=0.05 plus the same support floor.",
        "",
        "| heatmap | type | PNG | displayed cells | tested cells | rows |",
        "| --- | --- | --- | ---: | ---: | ---: |",
    ]
    for row in manifest.itertuples(index=False):
        lines.append(
            f"| {row.view_label} | {row.heatmap_type} | `{row.png_path}` | "
            f"{row.displayed_cell_count} | {row.tested_cell_count} | {row.row_count} |"
        )
    lines.append("")
    readme_path.write_text(prefix + "\n".join(lines), encoding="utf-8")


def write_six_shared_trait_cytoscape_heatmaps(
    output_dir: Path = OUTPUT_DIR,
    repo_root: Path = REPO_ROOT,
) -> pd.DataFrame:
    """Write approach-2 and approach-4 Cytoscape-derived heatmaps."""
    clean_legacy_outputs(output_dir)
    trait_order = read_trait_order(output_dir)
    graph_manifest = read_tsv(output_dir / GRAPH_MANIFEST_FILENAME)
    views_by_key = {row.view_key: row for row in graph_manifest.itertuples(index=False)}

    pattern_tables: list[pd.DataFrame] = []
    contrast_tables: list[pd.DataFrame] = []
    manifest_records: list[dict[str, object]] = []

    display_edges_by_view: dict[str, pd.DataFrame] = {}
    for view in graph_manifest.itertuples(index=False):
        display_edges = read_tsv(resolve_repo_path(view.source_display_edges_path, repo_root))
        display_edges_by_view[view.view_key] = display_edges
        cells = prepare_pattern_cells(display_edges, view, trait_order)
        pattern_tables.append(cells)
        output_path, row_count, displayed_count, tested_count = write_pattern_heatmap(
            cells,
            view,
            trait_order,
            output_dir,
        )
        manifest_records.append(
            {
                "view_key": view.view_key,
                "view_label": view.view_label,
                "heatmap_type": "trait_pattern",
                "focus_direction": view.focus_direction,
                "condition_mode": view.condition_mode,
                "row_count": row_count,
                "displayed_cell_count": displayed_count,
                "tested_cell_count": tested_count,
                "min_positive_profiles": MIN_POSITIVE_PROFILES,
                "min_studies": MIN_STUDIES,
                "p_value_cutoff": P_VALUE_CUTOFF,
                "png_path": relative_path(output_path, repo_root),
                "source_display_edges_path": view.source_display_edges_path,
            }
        )

    for mode in CONTRAST_MODES:
        increased_view = views_by_key[mode.increased_view_key]
        cells = prepare_contrast_cells(
            display_edges_by_view[mode.increased_view_key],
            display_edges_by_view[mode.decreased_view_key],
            mode,
            trait_order,
        )
        contrast_tables.append(cells)
        output_path, row_count, displayed_count, tested_count = write_contrast_heatmap(
            cells,
            mode,
            trait_order,
            output_dir,
        )
        manifest_records.append(
            {
                "view_key": f"{mode.key}_direction_contrast",
                "view_label": f"{mode.label}: signed direction contrast",
                "heatmap_type": "signed_direction_contrast",
                "focus_direction": "increased_vs_decreased",
                "condition_mode": mode.condition_mode,
                "row_count": row_count,
                "displayed_cell_count": displayed_count,
                "tested_cell_count": tested_count,
                "min_positive_profiles": MIN_POSITIVE_PROFILES,
                "min_studies": MIN_STUDIES,
                "p_value_cutoff": P_VALUE_CUTOFF,
                "png_path": relative_path(output_path, repo_root),
                "source_display_edges_path": increased_view.source_display_edges_path,
            }
        )

    pattern_table = pd.concat(pattern_tables, ignore_index=True)
    contrast_table = pd.concat(contrast_tables, ignore_index=True)
    manifest = pd.DataFrame.from_records(manifest_records)

    pattern_table.to_csv(output_dir / PATTERN_TABLE_FILENAME, sep="\t", index=False)
    contrast_table.to_csv(output_dir / CONTRAST_TABLE_FILENAME, sep="\t", index=False)
    manifest.to_csv(output_dir / MANIFEST_FILENAME, sep="\t", index=False)
    append_readme_section(output_dir, manifest)
    return manifest


def main() -> None:
    """Write the Cytoscape-derived heatmap outputs."""
    manifest = write_six_shared_trait_cytoscape_heatmaps()
    print(f"Wrote Cytoscape heatmaps: {len(manifest)}")
    print(f"Wrote manifest: {OUTPUT_DIR / MANIFEST_FILENAME}")


if __name__ == "__main__":
    main()
