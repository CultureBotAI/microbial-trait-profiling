"""Create an infographic summary of paired trait-profile contrasts."""

from __future__ import annotations

import math
from pathlib import Path
from textwrap import fill, shorten

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle, FancyBboxPatch, PathPatch
from matplotlib.path import Path as MplPath
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = VIS_DIR / "trait-profile-investigation"

PAIRED_MATRIX_PATH = VIS_DIR / "study_trait_profile_paired_enrichment_matrix.tsv"
PAIRED_METADATA_PATH = VIS_DIR / "study_trait_profile_paired_enrichment_metadata.tsv"

OUTPUT_PNG = INVESTIGATION_DIR / "trait_profile_paired_contrast_infographic.png"
TRAIT_STATS_TSV = INVESTIGATION_DIR / "trait_profile_paired_contrast_trait_stats.tsv"
CONDITION_SUMMARY_TSV = INVESTIGATION_DIR / "trait_profile_paired_contrast_condition_summary.tsv"
BODY_SITE_SUMMARY_TSV = INVESTIGATION_DIR / "trait_profile_paired_contrast_body_site_summary.tsv"

TOP_N_TRAITS = 10
TOP_N_CONDITIONS = 10
TOP_N_BODY_SITES = 9
MIN_CONDITION_PAIRS = 10
MIN_BODY_SITE_PAIRS = 3

EXCLUDED_CONDITION_LABELS = {
    "Body fat distribution",
    "Clinical treatment",
    "Exposure",
    "Health trait",
    "Microbiome",
    "Population",
    "Physical activity",
}

SLATE_50 = "#f8fafc"
SLATE_200 = "#e2e8f0"
SLATE_400 = "#94a3b8"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"
INCREASED_BLUE = "#3b82f6"
DECREASED_RED = "#ef4444"


def clean_text(value: object) -> str:
    """Normalize empty metadata and display values."""
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def compact_label(value: object, max_length: int = 42) -> str:
    """Return a compact single-line label for axes."""
    return shorten(clean_text(value).replace("\n", " "), width=max_length, placeholder="...")


def trait_label(value: object, max_length: int = 46) -> str:
    """Return the readable part of a trait display label."""
    label = clean_text(value).replace("\n", " ")
    return shorten(label, width=max_length, placeholder="...")


def is_excluded_condition_label(value: object) -> bool:
    """Return whether a condition label is too broad for condition-panel ranking."""
    label = clean_text(value)
    lower_label = label.lower()
    broad_terms = ("exposure", "measurement", "population")
    return (
        label in EXCLUDED_CONDITION_LABELS
        or label == "unknown"
        or any(term in lower_label for term in broad_terms)
    )


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """Calculate Benjamini-Hochberg adjusted q-values."""
    values = p_values.to_numpy(dtype=float)
    adjusted = np.ones_like(values, dtype=float)
    valid = np.isfinite(values)
    valid_indices = np.flatnonzero(valid)
    if len(valid_indices) == 0:
        return pd.Series(adjusted, index=p_values.index)

    sorted_indices = valid_indices[np.argsort(values[valid])]
    sorted_p = values[sorted_indices]
    ranks = np.arange(1, len(sorted_p) + 1, dtype=float)
    sorted_adjusted = np.minimum.accumulate((sorted_p * len(sorted_p) / ranks)[::-1])[::-1]
    adjusted[sorted_indices] = np.clip(sorted_adjusted, 0.0, 1.0)
    return pd.Series(adjusted, index=p_values.index)


def load_paired_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load paired contrast matrix and metadata, aligned by row order."""
    matrix = pd.read_csv(PAIRED_MATRIX_PATH, sep="\t", index_col=0)
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    metadata = pd.read_csv(PAIRED_METADATA_PATH, sep="\t", index_col=0)
    if len(matrix) != len(metadata):
        msg = (
            f"Paired matrix row count ({len(matrix)}) does not match metadata "
            f"row count ({len(metadata)})."
        )
        raise ValueError(msg)

    # The displayed paired matrix uses readable labels as its index, while the
    # metadata keeps stable pair IDs. The files are written in the same order.
    matrix.index = metadata.index
    return matrix, metadata


def compute_trait_stats(matrix: pd.DataFrame) -> pd.DataFrame:
    """Calculate paired one-sample t-test summaries for each trait column."""
    records: list[dict[str, object]] = []
    for column in matrix.columns:
        values = matrix[column].to_numpy(dtype=float)
        if len(values) < 3 or np.isclose(values.std(ddof=1), 0.0):
            t_stat = float("nan")
            p_value = float("nan")
        else:
            result = stats.ttest_1samp(values, 0.0)
            t_stat = float(result.statistic)
            p_value = float(result.pvalue)

        mean_contrast = float(values.mean())
        records.append(
            {
                "trait": str(column),
                "trait_label": trait_label(column, 80),
                "pair_count": len(values),
                "mean_increased_minus_decreased": mean_contrast,
                "median_increased_minus_decreased": float(np.median(values)),
                "std_increased_minus_decreased": float(values.std(ddof=1)),
                "positive_pair_fraction": float(np.mean(values > 0)),
                "negative_pair_fraction": float(np.mean(values < 0)),
                "nonzero_pair_fraction": float(np.mean(values != 0)),
                "t_statistic": t_stat,
                "p_value": p_value,
                "direction": (
                    "increased-profile higher" if mean_contrast >= 0 else "decreased-profile higher"
                ),
            }
        )

    stats_df = pd.DataFrame.from_records(records)
    stats_df["q_value"] = benjamini_hochberg(stats_df["p_value"].fillna(1.0))
    clipped_p = stats_df["p_value"].clip(lower=np.finfo(float).tiny, upper=1.0)
    stats_df["minus_log10_p"] = -np.log10(clipped_p)
    stats_df["abs_t_statistic"] = stats_df["t_statistic"].abs()
    return stats_df.sort_values(
        ["p_value", "abs_t_statistic", "trait_label"],
        ascending=[True, False, True],
    )


def add_pair_scores(matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Attach per-pair directional trait contrast loads to metadata."""
    scored = metadata.copy()
    scored["increased_profile_higher_load"] = matrix.clip(lower=0).mean(axis=1)
    scored["decreased_profile_higher_load"] = (-matrix.clip(upper=0)).mean(axis=1)
    scored["absolute_contrast_load"] = matrix.abs().mean(axis=1)
    scored["net_increased_minus_decreased"] = matrix.mean(axis=1)
    return scored


def summarize_groups(
    scored: pd.DataFrame,
    group_column: str,
    *,
    min_pairs: int,
) -> pd.DataFrame:
    """Summarize paired profile contrast load by condition or body-site group."""
    source = scored.copy()
    if group_column == "condition":
        source = source[~source[group_column].map(is_excluded_condition_label)]
    else:
        labels = source[group_column].map(clean_text).str.lower()
        source = source[(labels != "unknown") & (labels != "other")]

    summary = (
        source.groupby(group_column, dropna=False)
        .agg(
            pair_count=("study_id", "size"),
            study_count=("study_id", "nunique"),
            increased_profile_higher_load=("increased_profile_higher_load", "mean"),
            decreased_profile_higher_load=("decreased_profile_higher_load", "mean"),
            absolute_contrast_load=("absolute_contrast_load", "mean"),
            net_increased_minus_decreased=("net_increased_minus_decreased", "mean"),
        )
        .reset_index()
    )
    summary = summary[summary["pair_count"] >= min_pairs].copy()
    summary["dominant_direction"] = np.where(
        summary["net_increased_minus_decreased"] >= 0,
        "increased-profile higher",
        "decreased-profile higher",
    )
    return summary.sort_values(
        ["absolute_contrast_load", "pair_count", group_column],
        ascending=[False, False, True],
    )


def anatomical_site_base(value: object) -> tuple[float, float]:
    """Map body-site group labels to approximate body-diagram coordinates."""
    label = clean_text(value).lower()
    if "nasopharynx" in label or "airway" in label:
        return 0.30, 0.86
    if "subgingival" in label or "dental" in label:
        return 0.275, 0.825
    if "saliva" in label:
        return 0.325, 0.825
    if "mouth" in label or "oral" in label:
        return 0.30, 0.825
    if "colon" in label:
        return 0.285, 0.47
    if "caecum" in label or "cecum" in label:
        return 0.33, 0.46
    if "feces" in label or "faeces" in label or "stool" in label:
        return 0.31, 0.37
    if "vagina" in label or "urogenital" in label:
        return 0.31, 0.25
    return 0.16, 0.47


def add_anatomical_outline(ax: plt.Axes) -> None:
    """Draw an original front-facing body-template outline for site callouts."""
    face = "#ffffff"
    edge = SLATE_400

    def add_path(
        start: tuple[float, float],
        segments: list[tuple[tuple[float, float], ...]],
    ) -> None:
        vertices = [start]
        codes = [MplPath.MOVETO]
        for segment in segments:
            if len(segment) == 1:
                vertices.append(segment[0])
                codes.append(MplPath.LINETO)
            else:
                vertices.extend(segment)
                codes.extend([MplPath.CURVE4] * 3)
        vertices.append(start)
        codes.append(MplPath.CLOSEPOLY)
        ax.add_patch(
            PathPatch(
                MplPath(vertices, codes),
                facecolor=face,
                edgecolor=edge,
                lw=1.2,
                joinstyle="round",
                zorder=1,
            )
        )

    def mirror_point(point: tuple[float, float]) -> tuple[float, float]:
        return 0.60 - point[0], point[1]

    def mirror_segments(
        segments: list[tuple[tuple[float, float], ...]],
    ) -> list[tuple[tuple[float, float], ...]]:
        return [tuple(mirror_point(point) for point in segment) for segment in segments]

    ax.add_patch(Circle((0.30, 0.87), 0.045, facecolor=face, edgecolor=edge, lw=1.2, zorder=1))
    torso_segments = [
        ((0.224, 0.800), (0.207, 0.775), (0.196, 0.735)),
        ((0.178, 0.660), (0.178, 0.565), (0.195, 0.480)),
        ((0.205, 0.420), (0.220, 0.355), (0.240, 0.320)),
        ((0.265, 0.295), (0.335, 0.295), (0.360, 0.320)),
        ((0.380, 0.355), (0.395, 0.420), (0.405, 0.480)),
        ((0.422, 0.565), (0.422, 0.660), (0.404, 0.735)),
        ((0.393, 0.775), (0.376, 0.800), (0.355, 0.805)),
        ((0.335, 0.812), (0.265, 0.812), (0.245, 0.805)),
    ]
    left_arm_segments = [
        ((0.175, 0.665), (0.160, 0.585), (0.142, 0.505)),
        ((0.127, 0.445), (0.112, 0.395), (0.105, 0.365)),
        ((0.113, 0.350), (0.135, 0.346), (0.150, 0.358)),
        ((0.162, 0.395), (0.178, 0.455), (0.198, 0.535)),
        ((0.207, 0.600), (0.213, 0.665), (0.224, 0.715)),
    ]
    left_leg_segments = [
        ((0.258, 0.245), (0.250, 0.170), (0.238, 0.090)),
        ((0.233, 0.060), (0.250, 0.045), (0.278, 0.050)),
        ((0.292, 0.100), (0.295, 0.225), (0.295, 0.300)),
    ]
    add_path((0.245, 0.805), torso_segments)
    add_path((0.205, 0.725), left_arm_segments)
    add_path(mirror_point((0.205, 0.725)), mirror_segments(left_arm_segments))
    add_path((0.265, 0.300), left_leg_segments)
    add_path(mirror_point((0.265, 0.300)), mirror_segments(left_leg_segments))
    ax.plot([0.280, 0.280], [0.805, 0.770], color=edge, lw=1.0, zorder=2)
    ax.plot([0.320, 0.320], [0.805, 0.770], color=edge, lw=1.0, zorder=2)


def plot_top_traits(ax: plt.Axes, trait_stats: pd.DataFrame) -> None:
    """Plot panel A: traits with strongest paired contrast evidence."""
    top_traits = trait_stats.head(TOP_N_TRAITS).iloc[::-1].copy()
    colors = [
        INCREASED_BLUE if value >= 0 else DECREASED_RED
        for value in top_traits["mean_increased_minus_decreased"]
    ]
    labels = [trait_label(label, 38) for label in top_traits["trait_label"]]
    bars = ax.barh(labels, top_traits["minus_log10_p"], color=colors, edgecolor="white")
    ax.set_title("A. Top Profile-Differentiating Traits", fontsize=14, fontweight="bold", pad=8)
    ax.set_xlabel("-log10(p) from paired inc-minus-dec profile contrast")
    ax.grid(axis="x", color=SLATE_200, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=9)
    ax.tick_params(axis="x", labelsize=9)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(SLATE_200)
    ax.legend(
        handles=[
            plt.Rectangle((0, 0), 1, 1, color=INCREASED_BLUE, label="Higher in increased"),
            plt.Rectangle((0, 0), 1, 1, color=DECREASED_RED, label="Higher in decreased"),
        ],
        loc="lower right",
        frameon=False,
        fontsize=9,
    )

    x_offset = max(top_traits["minus_log10_p"]) * 0.015
    for bar, row in zip(bars, top_traits.itertuples(index=False), strict=True):
        label = f"mean {row.mean_increased_minus_decreased:+.2f}"
        ax.text(
            bar.get_width() + x_offset,
            bar.get_y() + bar.get_height() / 2,
            label,
            va="center",
            ha="left",
            fontsize=8.2,
            color=SLATE_700,
        )


def plot_condition_summary(ax: plt.Axes, condition_summary: pd.DataFrame) -> None:
    """Plot panel B: condition groups with strongest directional contrast load."""
    plot_data = condition_summary.head(TOP_N_CONDITIONS).iloc[::-1].copy()
    labels = [compact_label(label, 40) for label in plot_data["condition"]]
    decreased = plot_data["decreased_profile_higher_load"]
    increased = plot_data["increased_profile_higher_load"]
    ax.barh(labels, decreased, color=DECREASED_RED, alpha=0.82, label="Decreased-profile higher")
    ax.barh(
        labels,
        increased,
        left=decreased,
        color=INCREASED_BLUE,
        alpha=0.82,
        label="Increased-profile higher",
    )
    ax.set_title(
        "B. Avg Trait-Profile Contrast per Condition", fontsize=14, fontweight="bold", pad=8
    )
    ax.set_xlabel("Mean absolute paired trait-profile contrast")
    ax.grid(axis="x", color=SLATE_200, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=8.8)
    ax.tick_params(axis="x", labelsize=9)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(SLATE_200)
    ax.legend(
        loc="lower right",
        bbox_to_anchor=(1.0, 1.03),
        frameon=False,
        fontsize=9,
        ncol=2,
    )

    max_total = float(plot_data["absolute_contrast_load"].max())
    for y_index, row in enumerate(plot_data.itertuples(index=False)):
        ax.text(
            row.absolute_contrast_load + max_total * 0.015,
            y_index,
            f"{row.absolute_contrast_load:.2f}",
            va="center",
            ha="left",
            fontsize=8.2,
            color=SLATE_700,
            fontweight="bold",
        )


def plot_body_site_summary(ax: plt.Axes, body_site_summary: pd.DataFrame) -> None:
    """Plot panel C: anatomical body-site contrast hotspots."""
    plot_sites = body_site_summary.head(TOP_N_BODY_SITES).copy()
    plot_sites["anchor"] = [anatomical_site_base(site) for site in plot_sites["body_site_group"]]
    plot_sites["anchor_x"] = [point[0] for point in plot_sites["anchor"]]
    plot_sites["anchor_y"] = [point[1] for point in plot_sites["anchor"]]
    plot_sites = plot_sites.sort_values(["anchor_y", "absolute_contrast_load"], ascending=False)

    ax.set_title("C. Anatomical Profile-Contrast Hotspots", fontsize=14, fontweight="bold", pad=8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    add_anatomical_outline(ax)

    max_count = max(float(plot_sites["pair_count"].max()), 1.0)
    max_abs_net = max(float(plot_sites["net_increased_minus_decreased"].abs().max()), 0.05)
    norm = TwoSlopeNorm(vmin=-max_abs_net, vcenter=0, vmax=max_abs_net)
    cmap = plt.get_cmap("RdBu")
    label_ys = np.linspace(0.84, 0.18, len(plot_sites)) if len(plot_sites) > 1 else [0.52]

    for label_y, row in zip(label_ys, plot_sites.itertuples(index=False), strict=False):
        signal = float(row.net_increased_minus_decreased)
        color = cmap(norm(signal))
        size = 130 + 950 * math.sqrt(float(row.pair_count) / max_count)
        ax.scatter(
            [row.anchor_x],
            [row.anchor_y],
            s=size,
            c=[color],
            edgecolors=SLATE_900,
            linewidths=0.8,
            zorder=4,
        )
        ax.plot(
            [row.anchor_x, 0.55, 0.60],
            [row.anchor_y, label_y, label_y],
            color=SLATE_400,
            lw=0.9,
            zorder=2,
        )
        ax.text(
            0.61,
            label_y + 0.018,
            compact_label(row.body_site_group, 28),
            ha="left",
            va="bottom",
            fontsize=9.2,
            fontweight="bold",
            color=SLATE_900,
        )
        ax.text(
            0.61,
            label_y - 0.018,
            (
                f"{int(row.pair_count):,} pairs; load "
                f"{float(row.absolute_contrast_load):.2f}; net {signal:+.2f}"
            ),
            ha="left",
            va="top",
            fontsize=8.0,
            color=SLATE_700,
        )

    ax.text(
        0.02,
        0.04,
        "Bubble size = matched pairs; color = net inc-minus-dec trait contrast.",
        ha="left",
        va="bottom",
        fontsize=8.2,
        color=SLATE_500,
    )


def format_q(value: float) -> str:
    """Format q-values compactly for annotation text."""
    if not np.isfinite(value):
        return "NA"
    if value < 0.001:
        return f"{value:.1e}"
    return f"{value:.3f}"


def plot_story_box(
    ax: plt.Axes,
    trait_stats: pd.DataFrame,
    condition_summary: pd.DataFrame,
    body_site_summary: pd.DataFrame,
    *,
    pair_count: int,
    study_count: int,
    trait_count: int,
) -> None:
    """Plot panel D: interpretation text connecting this to enrichment analysis."""
    ax.set_title("D. How it Fits the Story", fontsize=14, fontweight="bold", pad=8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    box = FancyBboxPatch(
        (0.02, 0.02),
        0.96,
        0.90,
        boxstyle="round,pad=0.025,rounding_size=0.025",
        linewidth=1.0,
        edgecolor=SLATE_400,
        facecolor=SLATE_50,
    )
    ax.add_patch(box)

    top_increased = trait_stats[trait_stats["mean_increased_minus_decreased"] > 0].iloc[0]
    top_decreased = trait_stats[trait_stats["mean_increased_minus_decreased"] < 0].iloc[0]
    top_condition = condition_summary.iloc[0]
    top_body_site = body_site_summary.iloc[0]
    bullets = [
        (
            "Scope",
            (
                f"{trait_count:,} selected profile traits across {pair_count:,} matched "
                f"study/experiment inc-minus-dec pairs from {study_count:,} studies."
            ),
        ),
        (
            "Strongest split",
            (
                f"{trait_label(top_increased.trait_label, 42)} is higher in increased "
                f"profiles (mean {top_increased.mean_increased_minus_decreased:+.2f}, "
                f"q={format_q(float(top_increased.q_value))}); "
                f"{trait_label(top_decreased.trait_label, 42)} is higher in decreased "
                f"profiles (mean {top_decreased.mean_increased_minus_decreased:+.2f}, "
                f"q={format_q(float(top_decreased.q_value))})."
            ),
        ),
        (
            "Where it concentrates",
            (
                f"Highest condition burden: {compact_label(top_condition.condition, 42)}. "
                "Highest explicit body-site load: "
                f"{compact_label(top_body_site.body_site_group, 32)}."
            ),
        ),
        (
            "Relationship to Ilyas",
            (
                "His enrichment infographic ranks common enriched traits across profiles; "
                "this ranks matched directional profile contrasts. Differences are expected "
                "and help separate prevalence from direction-specific structure."
            ),
        ),
        (
            "Caveat",
            (
                "P-values and q-values are descriptive screens over selected profile traits, "
                "not new differential abundance tests."
            ),
        ),
    ]

    ax.text(
        0.06,
        0.86,
        "PROFILE CONTRAST SUMMARY",
        fontsize=11.5,
        fontweight="bold",
        color=SLATE_900,
        va="top",
    )
    y = 0.78
    for heading, body in bullets:
        wrapped = fill(f"{heading}: {body}", width=76)
        ax.text(
            0.06,
            y,
            f"- {wrapped}",
            fontsize=9.3,
            color=SLATE_700,
            va="top",
            linespacing=1.18,
        )
        y -= 0.042 + 0.038 * (wrapped.count("\n") + 1)


def create_infographic(
    trait_stats: pd.DataFrame,
    condition_summary: pd.DataFrame,
    body_site_summary: pd.DataFrame,
    *,
    pair_count: int,
    study_count: int,
    trait_count: int,
) -> None:
    """Render and save the four-panel paired trait-profile contrast infographic."""
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 14,
            "axes.labelsize": 10,
            "figure.facecolor": "white",
        }
    )
    fig = plt.figure(figsize=(18.5, 12.5), dpi=180, facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        left=0.07,
        right=0.97,
        top=0.84,
        bottom=0.08,
        hspace=0.42,
        wspace=0.38,
    )
    ax_traits = fig.add_subplot(grid[0, 0])
    ax_conditions = fig.add_subplot(grid[0, 1])
    ax_body = fig.add_subplot(grid[1, 0])
    ax_story = fig.add_subplot(grid[1, 1])

    fig.suptitle(
        "CultureBotAI: Trait-Profile Contrast Infographic",
        y=0.965,
        fontsize=24,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.5,
        0.925,
        (
            "Analog to enrichment infographic using matched increased-minus-decreased "
            "profile trait z-score contrasts"
        ),
        ha="center",
        va="center",
        fontsize=11.5,
        color=SLATE_500,
    )

    plot_top_traits(ax_traits, trait_stats)
    plot_condition_summary(ax_conditions, condition_summary)
    plot_body_site_summary(ax_body, body_site_summary)
    plot_story_box(
        ax_story,
        trait_stats,
        condition_summary,
        body_site_summary,
        pair_count=pair_count,
        study_count=study_count,
        trait_count=trait_count,
    )
    fig.text(
        0.07,
        0.035,
        (
            "Source: study_trait_profile_paired_enrichment_matrix.tsv and metadata. "
            "Rows are matched BugSigDB study/experiment increased-minus-decreased profile pairs."
        ),
        ha="left",
        va="bottom",
        fontsize=8.5,
        color=SLATE_500,
    )
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PNG, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    matrix, metadata = load_paired_inputs()
    trait_stats = compute_trait_stats(matrix)
    scored = add_pair_scores(matrix, metadata)
    condition_summary = summarize_groups(
        scored,
        "condition",
        min_pairs=MIN_CONDITION_PAIRS,
    )
    body_site_summary = summarize_groups(
        scored,
        "body_site_group",
        min_pairs=MIN_BODY_SITE_PAIRS,
    )

    INVESTIGATION_DIR.mkdir(parents=True, exist_ok=True)
    trait_stats.to_csv(TRAIT_STATS_TSV, sep="\t", index=False)
    condition_summary.to_csv(CONDITION_SUMMARY_TSV, sep="\t", index=False)
    body_site_summary.to_csv(BODY_SITE_SUMMARY_TSV, sep="\t", index=False)

    create_infographic(
        trait_stats,
        condition_summary,
        body_site_summary,
        pair_count=len(matrix),
        study_count=metadata["study_id"].nunique(),
        trait_count=matrix.shape[1],
    )
    print(f"Wrote: {OUTPUT_PNG}")
    print(f"Wrote: {TRAIT_STATS_TSV}")
    print(f"Wrote: {CONDITION_SUMMARY_TSV}")
    print(f"Wrote: {BODY_SITE_SUMMARY_TSV}")


if __name__ == "__main__":
    main()
