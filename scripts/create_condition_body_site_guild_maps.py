"""Create condition-to-body-site guild z-score maps."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from textwrap import shorten, wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path as MplPath

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = VIS_DIR / "trait-profile-investigation"

MATRIX_PATH = VIS_DIR / "study_trait_profile_clustered_heatmap_matrix.tsv"
METADATA_PATH = VIS_DIR / "study_trait_profile_study_metadata.tsv"
TRAIT_METADATA_PATH = VIS_DIR / "study_trait_profile_trait_metadata.tsv"

MIN_FOCUS_PROFILES = 3
MIN_STUDIES = 2
MAX_DISPLAY_EDGES = 18
MAX_DISPLAY_CONDITIONS = 12
MAX_LABEL_WIDTH = 36

SLATE_50 = "#f8fafc"
SLATE_200 = "#e2e8f0"
SLATE_400 = "#94a3b8"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"
EDGE_NEUTRAL = "#475569"


@dataclass(frozen=True)
class MapConfig:
    """Configuration for one condition/body-site guild map."""

    key: str
    title: str
    focus_label: str
    comparator_label: str
    focus_direction: str
    comparator_direction: str
    output_dir: Path
    output_stem: str
    excluded_condition_labels: frozenset[str]

    @property
    def guild_selection_path(self) -> Path:
        return self.output_dir / "guild_trait_selection.tsv"

    @property
    def edge_tsv_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}_edges.tsv"

    @property
    def png_path(self) -> Path:
        return self.output_dir / f"{self.output_stem}.png"


FACULTATIVE_CONFIG = MapConfig(
    key="facultative_aerobe_increased",
    title="Condition-to-Body-Site Map: Facultative-Aerobe Trait Enrichment in Increased Profiles",
    focus_label="increased",
    comparator_label="decreased",
    focus_direction="inc",
    comparator_direction="dec",
    output_dir=INVESTIGATION_DIR / "facultative-aerobe-v2",
    output_stem="condition_body_site_guild_map",
    excluded_condition_labels=frozenset(
        {
            "Body fat distribution",
            "Exposure",
            "Health trait",
            "Microbiome",
            "Population",
            "Physical activity",
        }
    ),
)

ANAEROBE_CONFIG = MapConfig(
    key="anaerobe_decreased",
    title="Condition-to-Body-Site Map: Anaerobe Trait Enrichment in Decreased Profiles",
    focus_label="decreased",
    comparator_label="increased",
    focus_direction="dec",
    comparator_direction="inc",
    output_dir=INVESTIGATION_DIR / "anaerobe-depletion-v2",
    output_stem="condition_body_site_guild_map",
    excluded_condition_labels=frozenset(
        {
            "Body fat distribution",
            "Clinical treatment",
            "Exposure",
            "Health study participation",
            "Health trait",
            "Microbiome",
            "Population",
            "Physical activity",
            "Timepoint",
        }
    ),
)


def clean_text(value: object) -> str:
    """Normalize missing or empty display values."""
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def compact_label(value: object, width: int = MAX_LABEL_WIDTH) -> str:
    """Return a compact label for node text."""
    return shorten(clean_text(value), width=width, placeholder="...")


def wrap_label(value: object, width: int = MAX_LABEL_WIDTH) -> str:
    """Wrap labels without letting them dominate the diagram."""
    label = compact_label(value, width=width)
    return "\n".join(wrap(label, width=max(12, width // 2))) or "unknown"


def is_excluded_condition_label(value: object, excluded: frozenset[str]) -> bool:
    """Return whether a condition label is too broad for display."""
    label = clean_text(value)
    lower_label = label.lower()
    broad_terms = ("exposure", "measurement", "population")
    return (
        label in excluded
        or label == "unknown"
        or any(term in lower_label for term in broad_terms)
    )


def is_display_body_site(value: object) -> bool:
    """Return whether a body-site group should be shown as an anatomical node."""
    label = clean_text(value).lower()
    return label not in {"unknown", "other"}


def read_tsv(path: Path, *, index_col: int | str | None = None) -> pd.DataFrame:
    """Read a TSV file with stable defaults."""
    return pd.read_csv(path, sep="\t", index_col=index_col)


def load_study_heatmap() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load profile-level z-score matrix and metadata with readable trait names."""
    matrix = read_tsv(MATRIX_PATH, index_col=0)
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    metadata = read_tsv(METADATA_PATH, index_col=0)
    trait_metadata = read_tsv(TRAIT_METADATA_PATH, index_col=0)

    matrix.index = metadata.index
    column_labels = {
        row["display_label"]: row["object_name"] for _, row in trait_metadata.iterrows()
    }
    matrix = matrix.rename(columns=column_labels)
    return matrix, metadata


def load_guild_traits(config: MapConfig) -> list[str]:
    """Load the selected traits for one analysis."""
    selection = read_tsv(config.guild_selection_path)
    trait_column = (
        "profile_matrix_object_name"
        if "profile_matrix_object_name" in selection.columns
        else "trait_name"
    )
    traits = [clean_text(value) for value in selection[trait_column]]
    return [trait for trait in traits if trait != "unknown"]


def add_guild_scores(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    config: MapConfig,
) -> pd.DataFrame:
    """Attach row-wise guild z-score means to profile metadata."""
    guild_traits = load_guild_traits(config)
    available_traits = [trait for trait in guild_traits if trait in matrix.columns]
    if not available_traits:
        msg = f"No selected traits for {config.key} were found in the profile matrix."
        raise ValueError(msg)

    scored = metadata.copy()
    scored["guild_z_mean"] = matrix[available_traits].mean(axis=1)
    scored["condition"] = scored["condition"].map(clean_text)
    scored["body_site_group"] = scored["body_site_group"].map(clean_text)
    return scored


def nan_mean(series: pd.Series) -> float:
    """Return a float mean, preserving NaN for empty groups."""
    if series.empty:
        return math.nan
    return float(series.mean())


def summarize_condition_body_site_edges(
    scored: pd.DataFrame,
    config: MapConfig,
) -> pd.DataFrame:
    """Summarize guild z-scores for each condition/body-site pair."""
    condition_is_allowed = ~scored["condition"].map(
        lambda value: is_excluded_condition_label(value, config.excluded_condition_labels)
    )
    source = scored[
        condition_is_allowed & scored["body_site_group"].map(is_display_body_site)
    ].copy()

    records: list[dict[str, object]] = []
    for (condition, body_site_group), group in source.groupby(
        ["condition", "body_site_group"],
        sort=False,
        dropna=False,
    ):
        focus = group[group["inc_dec"] == config.focus_direction]
        comparator = group[group["inc_dec"] == config.comparator_direction]
        if focus.empty:
            continue

        focus_mean = nan_mean(focus["guild_z_mean"])
        comparator_mean = nan_mean(comparator["guild_z_mean"])
        contrast = focus_mean - comparator_mean if np.isfinite(comparator_mean) else math.nan
        records.append(
            {
                "analysis": config.key,
                "condition": clean_text(condition),
                "body_site_group": clean_text(body_site_group),
                "profile_count": len(group),
                "study_count": group["study_id"].nunique(),
                "focus_direction": config.focus_label,
                "comparator_direction": config.comparator_label,
                "focus_profile_count": len(focus),
                "comparator_profile_count": len(comparator),
                "increased_profile_count": int((group["inc_dec"] == "inc").sum()),
                "decreased_profile_count": int((group["inc_dec"] == "dec").sum()),
                "mean_guild_z": nan_mean(group["guild_z_mean"]),
                "increased_mean_guild_z": nan_mean(
                    group.loc[group["inc_dec"] == "inc", "guild_z_mean"]
                ),
                "decreased_mean_guild_z": nan_mean(
                    group.loc[group["inc_dec"] == "dec", "guild_z_mean"]
                ),
                "focus_mean_guild_z": focus_mean,
                "comparator_mean_guild_z": comparator_mean,
                "focus_minus_comparator_guild_z": contrast,
            }
        )

    if not records:
        return pd.DataFrame()

    edges = pd.DataFrame.from_records(records)
    edges["rank_signal"] = edges["focus_minus_comparator_guild_z"].fillna(
        edges["focus_mean_guild_z"]
    )
    edges["rank_score"] = edges["rank_signal"] * np.log1p(edges["focus_profile_count"])
    edges["displayed_in_map"] = False
    return edges.sort_values(
        ["rank_score", "rank_signal", "focus_profile_count", "study_count"],
        ascending=[False, False, False, False],
    )


def select_display_edges(edges: pd.DataFrame) -> pd.DataFrame:
    """Choose a readable set of high-signal supported edges for the PNG."""
    candidates = edges[
        (edges["focus_profile_count"] >= MIN_FOCUS_PROFILES)
        & (edges["study_count"] >= MIN_STUDIES)
        & (edges["focus_mean_guild_z"] > 0)
    ].copy()
    preferred = candidates[candidates["focus_minus_comparator_guild_z"] > 0].copy()
    if len(preferred) >= 6:
        candidates = preferred
    if candidates.empty:
        candidates = edges.head(MAX_DISPLAY_EDGES).copy()

    candidates = candidates.sort_values(
        ["rank_score", "rank_signal", "focus_profile_count", "study_count"],
        ascending=[False, False, False, False],
    )
    per_body_site = candidates.groupby("body_site_group", sort=False).head(2)
    selected_indices: list[int] = []
    selected_conditions: set[str] = set()

    def add_candidate(index: int, row: pd.Series) -> None:
        if index in selected_indices or len(selected_indices) >= MAX_DISPLAY_EDGES:
            return
        condition = str(row["condition"])
        condition_limit_reached = len(selected_conditions) >= MAX_DISPLAY_CONDITIONS
        if condition not in selected_conditions and condition_limit_reached:
            return
        selected_indices.append(index)
        selected_conditions.add(condition)

    for index, row in per_body_site.iterrows():
        add_candidate(index, row)
    for index, row in candidates.iterrows():
        add_candidate(index, row)

    display = candidates.loc[selected_indices].sort_values(
        ["rank_score", "rank_signal", "focus_profile_count", "study_count"],
        ascending=[False, False, False, False],
    )
    return display.copy()


def scaled_line_width(values: pd.Series) -> pd.Series:
    """Scale profile counts to readable line widths."""
    counts = values.astype(float)
    if counts.max() == counts.min():
        return pd.Series(4.0, index=values.index)
    scaled = (np.sqrt(counts) - np.sqrt(counts.min())) / (
        np.sqrt(counts.max()) - np.sqrt(counts.min())
    )
    return 1.2 + scaled * 8.0


def node_positions(labels: list[str]) -> dict[str, float]:
    """Return evenly spaced vertical positions for node labels."""
    if len(labels) == 1:
        return {labels[0]: 0.52}
    positions = np.linspace(0.88, 0.16, len(labels))
    return dict(zip(labels, positions, strict=True))


def weighted_node_stats(edges: pd.DataFrame, node_column: str) -> pd.DataFrame:
    """Summarize displayed node support and focus mean z-score."""
    rows: list[dict[str, object]] = []
    for label, group in edges.groupby(node_column, sort=False):
        weights = group["focus_profile_count"].to_numpy(dtype=float)
        focus_values = group["focus_mean_guild_z"].to_numpy(dtype=float)
        contrast_values = group["focus_minus_comparator_guild_z"].fillna(
            group["focus_mean_guild_z"]
        ).to_numpy(dtype=float)
        rows.append(
            {
                node_column: label,
                "focus_profile_count": int(group["focus_profile_count"].sum()),
                "study_count": int(group["study_count"].max()),
                "focus_mean_guild_z": float(np.average(focus_values, weights=weights)),
                "focus_minus_comparator_guild_z": float(
                    np.average(contrast_values, weights=weights)
                ),
            }
        )
    return pd.DataFrame.from_records(rows).set_index(node_column)


def draw_curve(
    ax: plt.Axes,
    *,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    linewidth: float,
) -> None:
    """Draw one left-to-right Bezier edge."""
    x0, y0 = start
    x1, y1 = end
    control_dx = (x1 - x0) * 0.46
    path = MplPath(
        [
            (x0, y0),
            (x0 + control_dx, y0),
            (x1 - control_dx, y1),
            (x1, y1),
        ],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
    )
    patch = PathPatch(
        path,
        facecolor="none",
        edgecolor=color,
        linewidth=linewidth,
        alpha=0.72,
        capstyle="round",
        joinstyle="round",
    )
    ax.add_patch(patch)


def draw_node(
    ax: plt.Axes,
    *,
    x: float,
    y: float,
    label: str,
    stats_row: pd.Series,
    align: str,
) -> None:
    """Draw a node marker and its text label."""
    ax.add_patch(
        Circle(
            (x, y),
            0.012,
            facecolor="white",
            edgecolor=SLATE_700,
            linewidth=1.0,
            zorder=5,
        )
    )
    text_x = x - 0.025 if align == "right" else x + 0.025
    ax.text(
        text_x,
        y,
        (
            f"{wrap_label(label)}\n"
            f"n={int(stats_row.focus_profile_count)}, z={stats_row.focus_mean_guild_z:+.2f}"
        ),
        ha=align,
        va="center",
        fontsize=8.6,
        color=SLATE_900,
        linespacing=1.08,
    )


def plot_condition_body_site_map(
    display_edges: pd.DataFrame,
    config: MapConfig,
    output_path: Path,
) -> None:
    """Render a condition-to-body-site bipartite guild map."""
    if display_edges.empty:
        msg = f"No display edges available for {config.key}."
        raise ValueError(msg)

    condition_order = (
        display_edges.groupby("condition", sort=False)["rank_score"]
        .max()
        .sort_values(ascending=False)
        .index.tolist()
    )
    body_site_order = (
        display_edges.groupby("body_site_group", sort=False)["focus_profile_count"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    )
    condition_y = node_positions(condition_order)
    body_site_y = node_positions(body_site_order)
    condition_stats = weighted_node_stats(display_edges, "condition")
    body_site_stats = weighted_node_stats(display_edges, "body_site_group")

    height = max(8.0, 1.0 + 0.58 * max(len(condition_order), len(body_site_order)))
    fig, ax = plt.subplots(figsize=(13.5, height))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    values = display_edges["focus_minus_comparator_guild_z"].fillna(
        display_edges["focus_mean_guild_z"]
    )
    vmax = max(0.4, float(values.abs().max()))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    cmap = plt.get_cmap("RdBu_r")
    widths = scaled_line_width(display_edges["focus_profile_count"])

    left_x = 0.34
    right_x = 0.68
    for edge, linewidth in zip(display_edges.itertuples(), widths, strict=True):
        signal = (
            edge.focus_minus_comparator_guild_z
            if np.isfinite(edge.focus_minus_comparator_guild_z)
            else edge.focus_mean_guild_z
        )
        draw_curve(
            ax,
            start=(left_x, condition_y[edge.condition]),
            end=(right_x, body_site_y[edge.body_site_group]),
            color=cmap(norm(signal)),
            linewidth=float(linewidth),
        )

    for condition in condition_order:
        draw_node(
            ax,
            x=left_x,
            y=condition_y[condition],
            label=condition,
            stats_row=condition_stats.loc[condition],
            align="right",
        )
    for body_site in body_site_order:
        draw_node(
            ax,
            x=right_x,
            y=body_site_y[body_site],
            label=body_site,
            stats_row=body_site_stats.loc[body_site],
            align="left",
        )

    fig.text(0.05, 0.965, config.title, ha="left", va="top", fontsize=17, fontweight="bold")
    fig.text(
        0.05,
        0.925,
        (
            "Edges connect BugSigDB conditions to anatomical body-site groups. "
            f"Width = {config.focus_label}-profile support; color = "
            f"{config.focus_label} minus {config.comparator_label} mean trait z-score."
        ),
        ha="left",
        va="top",
        fontsize=10,
        color=SLATE_700,
    )
    ax.text(0.21, 0.93, "Conditions", ha="center", va="center", fontsize=12, color=SLATE_700)
    ax.text(0.79, 0.93, "Body-site groups", ha="center", va="center", fontsize=12, color=SLATE_700)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.030, pad=0.02, shrink=0.55)
    cbar.set_label(
        f"{config.focus_label} - {config.comparator_label} mean trait z-score",
        fontsize=9,
    )
    cbar.ax.tick_params(labelsize=8)

    representative_counts = sorted(
        {
            int(display_edges["focus_profile_count"].min()),
            int(display_edges["focus_profile_count"].median()),
            int(display_edges["focus_profile_count"].max()),
        }
    )
    legend_lines = [
        Line2D(
            [0],
            [0],
            color=EDGE_NEUTRAL,
            linewidth=float(scaled_line_width(pd.Series([count, max(representative_counts)]))[0]),
            alpha=0.72,
            solid_capstyle="round",
            label=f"n={count}",
        )
        for count in representative_counts
    ]
    legend = ax.legend(
        handles=legend_lines,
        title=f"{config.focus_label} profiles",
        loc="lower center",
        bbox_to_anchor=(0.50, -0.015),
        ncol=len(legend_lines),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    ax.add_artist(legend)

    fig.text(
        0.05,
        0.035,
        (
            f"Displayed edges require at least {MIN_FOCUS_PROFILES} {config.focus_label} profiles, "
            f"{MIN_STUDIES} studies, positive focus mean z-score, and are ranked by "
            "contrast-weighted support. Other/unknown body-site groups and broad metadata-like "
            "conditions are excluded."
        ),
        ha="left",
        va="bottom",
        fontsize=8.2,
        color=SLATE_500,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_map_outputs(matrix: pd.DataFrame, metadata: pd.DataFrame, config: MapConfig) -> None:
    """Write the edge table and PNG for one analysis."""
    scored = add_guild_scores(matrix, metadata, config)
    edges = summarize_condition_body_site_edges(scored, config)
    display_edges = select_display_edges(edges)
    display_keys = set(
        zip(display_edges["condition"], display_edges["body_site_group"], strict=True)
    )
    edges["displayed_in_map"] = [
        (condition, body_site) in display_keys
        for condition, body_site in zip(edges["condition"], edges["body_site_group"], strict=True)
    ]
    edges.to_csv(config.edge_tsv_path, sep="\t", index=False)
    plot_condition_body_site_map(display_edges, config, config.png_path)
    print(f"Wrote: {config.edge_tsv_path}")
    print(f"Wrote: {config.png_path}")


def main() -> None:
    matrix, metadata = load_study_heatmap()
    for config in [FACULTATIVE_CONFIG, ANAEROBE_CONFIG]:
        write_map_outputs(matrix, metadata, config)


if __name__ == "__main__":
    main()
