"""Create inferential heatmaps from the six-shared-trait Cytoscape inputs."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import fisher_exact

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "trait-visualization" / "trait-profile-investigation" / "six-shared-traits"

GRAPH_MANIFEST_FILENAME = "six_shared_trait_graph_manifest.tsv"
TRAIT_SELECTION_FILENAME = "six_shared_trait_selection.tsv"
MANIFEST_FILENAME = "six_shared_trait_cytoscape_heatmap_manifest.tsv"
PATTERN_TABLE_FILENAME = "six_shared_trait_condition_body_site_trait_pattern_matrix.tsv"
CONTRAST_TABLE_FILENAME = "six_shared_trait_condition_body_site_direction_contrasts.tsv"
SELECTED_CELLS_FILENAME = "six_shared_trait_condition_body_site_bh_selected_cells.tsv"
CONTRAST_PAGES_FILENAME = "six_shared_trait_signed_direction_heatmap_pages.tsv"
GROUPED_CONTRAST_CLUSTER_ROWS_FILENAME = (
    "six_shared_traits_grouped_signed_direction_cluster_rows.tsv"
)
GROUPED_CONTRAST_CLUSTER_SUMMARY_FILENAME = (
    "six_shared_traits_grouped_signed_direction_cluster_summary.tsv"
)
GROUPED_CONTRAST_BRANCH_TESTS_FILENAME = (
    "six_shared_traits_grouped_signed_direction_branch_label_tests.tsv"
)
GROUPED_CONTRAST_PATTERN_TESTS_FILENAME = (
    "six_shared_traits_grouped_signed_direction_branch_pattern_tests.tsv"
)
CONTRAST_ROWS_PER_PAGE = 24

MIN_POSITIVE_PROFILES = 3
MIN_STUDIES = 2
P_VALUE_CUTOFF = 0.01
Q_VALUE_CUTOFF = 0.05
SIGNED_EFFECT_ABS_CUTOFF = 0.01
MIN_SIGNIFICANT_BRANCH_ROWS = 3
MIN_SIGNIFICANT_BRANCH_LABEL_COUNT = 2
MIN_SIGNIFICANT_BRANCH_LABEL_PURITY = 2.0 / 3.0
MIN_PATTERN_DIRECTION_ROWS = 3
MIN_PATTERN_DIRECTION_PURITY = 2.0 / 3.0
MIN_PATTERN_BH_SELECTED_CELLS = 1
MAX_MINUS_LOG10_P = 10.0
MAX_PATTERN_ROWS = 24

SLATE_50 = "#f8fafc"
SLATE_100 = "#f1f5f9"
SLATE_200 = "#e2e8f0"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"
SIGNIFICANCE_BORDER = "#111827"

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


@dataclass(frozen=True)
class ContrastLegendSummary:
    """Counts and display rules shown in a signed-contrast figure legend."""

    condition_label_kind: str
    condition_label_count: int
    anatomy_group_count: int
    row_count: int
    displayed_cell_count: int
    nominal_p_value_cell_count: int
    displayed_bh_selected_cell_count: int
    effect_cutoff_hidden_cell_count: int


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


def write_tsv(table: pd.DataFrame, path: Path) -> None:
    """Write TSVs with explicit missing values instead of trailing empty fields."""
    table.replace("", pd.NA).to_csv(path, sep="\t", index=False, na_rep="NA")


def relative_path(path: Path, repo_root: Path = REPO_ROOT) -> str:
    """Return a stable repo-relative path."""
    return path.relative_to(repo_root).as_posix()


def resolve_repo_path(path_text: str, repo_root: Path = REPO_ROOT) -> Path:
    """Resolve a manifest path that may be repo-relative or absolute."""
    path = Path(path_text)
    return path if path.is_absolute() else repo_root / path


def full_edges_path(display_path: str, repo_root: Path = REPO_ROOT) -> Path:
    """Resolve the complete edge table paired with a graph display table."""
    path = resolve_repo_path(display_path, repo_root)
    suffix = "_display_edges.tsv"
    if not path.name.endswith(suffix):
        raise ValueError(f"Expected a {suffix} source, got {path.name}")
    return path.with_name(path.name.removesuffix(suffix) + "_edges.tsv")


def wrap_label(value: object, width: int = 24) -> str:
    """Wrap labels for compact heatmap axes."""
    text = str(value)
    return "\n".join(wrap(text, width=width, break_long_words=False)) or text


def compact_dendrogram_label(value: object, *, width: int) -> str:
    """Return a compact branch title for placement inside the dendrogram."""
    text = str(value)
    matches = re.findall(
        r"(Increased|Decreased)-associated (.*?) "
        r"\((\d+/\d+) rows; q=([^;]+); (\d+) BH cells?\)",
        text,
    )
    if not matches:
        return wrap_label(text, width=width)
    lines = []
    for direction, trait, row_fraction, q_value, bh_count in matches:
        lines.append(wrap_label(f"{direction}: {trait}", width=width))
        lines.append(f"{row_fraction} | q={q_value} | {bh_count} BH")
    return "\n".join(lines)


def format_p_value(value: float) -> str:
    """Format p-values for cell annotations."""
    if not math.isfinite(value):
        return "NA"
    if value < 0.001:
        return "<0.001"
    if value < P_VALUE_CUTOFF <= round(value, 3):
        return f"<{P_VALUE_CUTOFF:g}"
    return f"{value:.3f}"


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """Calculate Benjamini-Hochberg adjusted q-values."""
    values = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
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


def minus_log10_p_values(p_values: pd.Series) -> pd.Series:
    """Convert finite p-values to positive evidence scores for heatmap coloring."""
    values = pd.to_numeric(p_values, errors="coerce")
    finite = np.isfinite(values)
    output = pd.Series(np.nan, index=p_values.index, dtype=float)
    clipped = values[finite].clip(lower=np.finfo(float).tiny, upper=1.0)
    output.loc[finite] = -np.log10(clipped)
    return output


def finalize_inferential_cells(cells: pd.DataFrame) -> pd.DataFrame:
    """Add p/q-value decisions and heatmap color values to tested cells."""
    output = cells.copy()
    p_values = pd.to_numeric(output["fisher_p_value"], errors="coerce")
    finite_p = np.isfinite(p_values)
    support = output["passes_support"].astype(bool)

    output["fisher_q_value"] = benjamini_hochberg(p_values.fillna(1.0))
    q_values = pd.to_numeric(output["fisher_q_value"], errors="coerce")
    finite_q = np.isfinite(q_values)
    output["minus_log10_fisher_p_value"] = minus_log10_p_values(p_values)
    output["plot_minus_log10_fisher_p_value"] = output["minus_log10_fisher_p_value"].clip(
        upper=MAX_MINUS_LOG10_P
    )
    output["passes_nominal_p_value_cutoff"] = finite_p & (p_values < P_VALUE_CUTOFF)
    output["passes_p_value_cutoff"] = output["passes_nominal_p_value_cutoff"]
    output["nominal_p_value_selected"] = support & output["passes_nominal_p_value_cutoff"]
    output["displayed_in_heatmap"] = output["nominal_p_value_selected"]
    # Compatibility alias for existing consumers. This is a nominal-p flag, not
    # the multiplicity-adjusted selection used for signed-heatmap borders.
    output["significant_in_heatmap"] = output["nominal_p_value_selected"]
    output["selected_after_bh"] = support & finite_q & (q_values <= Q_VALUE_CUTOFF)

    output["blank_reason"] = ""
    output.loc[~support, "blank_reason"] = "low positive-z profile or study support"
    output.loc[support & ~finite_p, "blank_reason"] = "Fisher p=NA"
    output.loc[
        support & finite_p & ~output["passes_nominal_p_value_cutoff"],
        "blank_reason",
    ] = f"Fisher p>={P_VALUE_CUTOFF:g}"

    if "increased_minus_decreased_positive_fraction" in output:
        delta = pd.to_numeric(
            output["increased_minus_decreased_positive_fraction"],
            errors="coerce",
        )
        direction = np.sign(delta.fillna(0.0))
        output["signed_minus_log10_fisher_p_value"] = (
            direction * output["minus_log10_fisher_p_value"]
        )
        output["plot_signed_minus_log10_fisher_p_value"] = (
            direction * output["plot_minus_log10_fisher_p_value"]
        )
        output["passes_absolute_effect_cutoff"] = np.isfinite(delta) & (
            delta.abs() > SIGNED_EFFECT_ABS_CUTOFF
        )
        output["displayed_in_heatmap"] = (
            support & finite_p & output["passes_absolute_effect_cutoff"]
        )
        output.loc[
            support & finite_p & ~output["passes_absolute_effect_cutoff"],
            "blank_reason",
        ] = f"absolute effect <={SIGNED_EFFECT_ABS_CUTOFF:g}"
        output.loc[output["displayed_in_heatmap"], "blank_reason"] = ""

    return output


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


def summarize_contrast_legend(
    cells: pd.DataFrame,
    row_labels: list[str],
    mode: ContrastMode,
) -> ContrastLegendSummary:
    """Summarize the visible condition/anatomy scope and signed-cell rules."""
    selected = cells[cells["condition_body_site_pair"].isin(row_labels)].copy()
    rows = selected.drop_duplicates("condition_body_site_pair")
    displayed = selected["displayed_in_heatmap"].astype(bool)
    bh_selected = selected["selected_after_bh"].astype(bool)
    cutoff_hidden = (
        selected["passes_support"].astype(bool)
        & np.isfinite(pd.to_numeric(selected["fisher_p_value"], errors="coerce"))
        & ~selected["passes_absolute_effect_cutoff"].astype(bool)
    )
    label_kind = "grouped condition labels" if mode.key == "grouped" else "individual labels"
    return ContrastLegendSummary(
        condition_label_kind=label_kind,
        condition_label_count=int(rows["condition"].nunique()),
        anatomy_group_count=int(rows["body_site_group"].nunique()),
        row_count=len(rows),
        displayed_cell_count=int(displayed.sum()),
        nominal_p_value_cell_count=int(selected["nominal_p_value_selected"].sum()),
        displayed_bh_selected_cell_count=int((displayed & bh_selected).sum()),
        effect_cutoff_hidden_cell_count=int(cutoff_hidden.sum()),
    )


def add_contrast_legend(
    fig: Figure,
    summary: ContrastLegendSummary,
    *,
    clustered: bool,
    significant_branch_count: int = 0,
    significant_metadata_label_count: int = 0,
    significant_trait_pattern_count: int = 0,
) -> None:
    """Add a count-aware signed-contrast legend below a heatmap."""
    handles: list[Patch | Line2D] = [
        Patch(
            facecolor="#d6604d",
            edgecolor=SLATE_200,
            label=(
                f"Colored: |effect| > {SIGNED_EFFECT_ABS_CUTOFF:g} "
                f"({summary.displayed_cell_count} cells)"
            ),
        ),
        Patch(
            facecolor="white",
            edgecolor=SIGNIFICANCE_BORDER,
            linewidth=1.8,
            label=(
                f"Dark border: BH-adjusted q <= {Q_VALUE_CUTOFF:g} "
                f"({summary.displayed_bh_selected_cell_count} cells)"
            ),
        ),
        Patch(
            facecolor=SLATE_100,
            edgecolor=SLATE_200,
            label=(
                "Gray: absent, under-supported, or below effect cutoff "
                f"({summary.effect_cutoff_hidden_cell_count} shown-row cells at cutoff)"
            ),
        ),
        Line2D(
            [0],
            [0],
            color=SLATE_500,
            linewidth=1.2,
            linestyle=":",
            label=(
                f"Nominal Fisher p < {P_VALUE_CUTOFF:g}: "
                f"{summary.nominal_p_value_cell_count} cells (audit only)"
            ),
        ),
    ]
    if clustered:
        handles.append(
            Line2D(
                [0],
                [0],
                color="black",
                linewidth=1.2,
                label=(
                    f"Labeled branches: {significant_branch_count} "
                    f"({significant_metadata_label_count} metadata labels; "
                    f"{significant_trait_pattern_count} trait patterns)"
                ),
            )
        )
    legend = fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        ncol=2,
        frameon=True,
        fontsize=7.8,
        title=(
            f"Shown: {summary.condition_label_count} {summary.condition_label_kind} | "
            f"{summary.anatomy_group_count} anatomy groups | "
            f"{summary.row_count} condition x anatomy rows"
        ),
        title_fontsize=8.2,
    )
    legend.get_frame().set_edgecolor(SLATE_200)
    legend.get_frame().set_facecolor("white")


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
                "passes_support": supported,
            }
        )
    return finalize_inferential_cells(pd.DataFrame.from_records(records))


def select_pattern_rows(cells: pd.DataFrame, *, max_rows: int = MAX_PATTERN_ROWS) -> list[str]:
    """Select condition/body-site rows with the strongest significant trait recurrence."""
    cells = cells.copy()
    cells["display_evidence"] = cells["plot_minus_log10_fisher_p_value"].where(
        cells["displayed_in_heatmap"],
        0.0,
    )
    cells["display_positive_profiles"] = cells["focus_positive_profile_count"].where(
        cells["displayed_in_heatmap"],
        0,
    )
    cells["display_strength"] = cells["strength"].where(cells["displayed_in_heatmap"], 0.0)
    stats = (
        cells.groupby("condition_body_site_pair", sort=False)
        .agg(
            significant_trait_count=("displayed_in_heatmap", "sum"),
            bh_selected_trait_count=("selected_after_bh", "sum"),
            total_positive_profiles=("display_positive_profiles", "sum"),
            max_evidence=("display_evidence", "max"),
            total_evidence=("display_evidence", "sum"),
            total_strength=("display_strength", "sum"),
        )
        .sort_values(
            [
                "significant_trait_count",
                "bh_selected_trait_count",
                "total_evidence",
                "total_positive_profiles",
                "total_strength",
                "max_evidence",
            ],
            ascending=[False, False, False, False, False, False],
        )
    )
    stats = stats[stats["significant_trait_count"] > 0]
    if stats.empty:
        stats = (
            cells.groupby("condition_body_site_pair", sort=False)
            .agg(
                significant_trait_count=("displayed_in_heatmap", "sum"),
                total_positive_profiles=("focus_positive_profile_count", "sum"),
                max_evidence=("minus_log10_fisher_p_value", "max"),
                total_evidence=("minus_log10_fisher_p_value", "sum"),
                max_strength=("strength", "max"),
                total_strength=("strength", "sum"),
            )
            .sort_values(
                ["total_evidence", "total_positive_profiles", "total_strength", "max_strength"],
                ascending=[False, False, False, False],
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
                "passes_support": supported,
            }
        )
    return finalize_inferential_cells(pd.DataFrame.from_records(output_records))


def select_contrast_rows(cells: pd.DataFrame) -> list[str]:
    """Group all supported rows by body site, then rank by evidence and coverage."""
    cells = cells.copy()
    cells["abs_delta"] = cells["increased_minus_decreased_positive_fraction"].abs()
    cells["abs_display_evidence"] = (
        cells["plot_signed_minus_log10_fisher_p_value"]
        .abs()
        .where(
            cells["displayed_in_heatmap"],
            0.0,
        )
    )
    cells["positive_support"] = (
        cells["increased_positive_profile_count"] + cells["decreased_positive_profile_count"]
    )
    cells["display_positive_support"] = cells["positive_support"].where(
        cells["displayed_in_heatmap"],
        0,
    )
    cells["display_abs_delta"] = cells["abs_delta"].where(cells["displayed_in_heatmap"], 0.0)
    stats = (
        cells.groupby("condition_body_site_pair", sort=False)
        .agg(
            body_site_group=("body_site_group", "first"),
            displayed_trait_count=("displayed_in_heatmap", "sum"),
            significant_trait_count=("significant_in_heatmap", "sum"),
            bh_selected_trait_count=("selected_after_bh", "sum"),
            positive_support=("display_positive_support", "sum"),
            max_abs_delta=("display_abs_delta", "max"),
            total_abs_delta=("display_abs_delta", "sum"),
            max_evidence=("abs_display_evidence", "max"),
            total_evidence=("abs_display_evidence", "sum"),
        )
        .sort_values(
            [
                "body_site_group",
                "significant_trait_count",
                "displayed_trait_count",
                "bh_selected_trait_count",
                "total_evidence",
                "positive_support",
                "total_abs_delta",
                "max_abs_delta",
            ],
            ascending=[True, False, False, False, False, False, False, False],
        )
    )
    stats = stats[stats["displayed_trait_count"] > 0]
    return stats.index.astype(str).tolist()


def top_values(values: pd.Series, *, max_items: int = 3) -> str:
    """Return a semicolon-delimited descriptive summary for cluster audits."""
    counts = values.astype(str).value_counts()
    if counts.empty:
        return "NA"
    return "; ".join(f"{label} ({count})" for label, count in counts.head(max_items).items())


def spearman_row_linkage(matrix: pd.DataFrame) -> np.ndarray | None:
    """Cluster rows with average linkage on Spearman-correlation distance."""
    if len(matrix.index) < 2:
        return None
    values = matrix.fillna(0.0).astype(float)
    correlations = values.T.corr(method="spearman").fillna(0.0).clip(-1.0, 1.0)
    correlation_values = correlations.to_numpy(dtype=float, copy=True)
    np.fill_diagonal(correlation_values, 1.0)
    distances = np.clip(1.0 - correlation_values, 0.0, 2.0)
    condensed = squareform(distances, checks=False)
    return linkage(condensed, method="average")


def contrast_row_stats(cells: pd.DataFrame) -> pd.DataFrame:
    """Summarize displayed contrast rows for ordering and cluster audits."""
    frame = cells.copy()
    frame["positive_support"] = (
        frame["increased_positive_profile_count"] + frame["decreased_positive_profile_count"]
    ).where(frame["displayed_in_heatmap"], 0)
    frame["abs_delta"] = (
        frame["increased_minus_decreased_positive_fraction"]
        .abs()
        .where(
            frame["displayed_in_heatmap"],
            0.0,
        )
    )
    frame["abs_evidence"] = (
        frame["plot_signed_minus_log10_fisher_p_value"]
        .abs()
        .where(
            frame["displayed_in_heatmap"],
            0.0,
        )
    )
    return frame.groupby("condition_body_site_pair", sort=False).agg(
        condition=("condition", "first"),
        body_site_group=("body_site_group", "first"),
        displayed_trait_count=("displayed_in_heatmap", "sum"),
        significant_trait_count=("significant_in_heatmap", "sum"),
        bh_selected_trait_count=("selected_after_bh", "sum"),
        positive_support=("positive_support", "sum"),
        total_abs_delta=("abs_delta", "sum"),
        max_abs_delta=("abs_delta", "max"),
        total_evidence=("abs_evidence", "sum"),
        max_evidence=("abs_evidence", "max"),
    )


def clustered_contrast_rows(
    matrix: pd.DataFrame,
    cells: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, np.ndarray | None]:
    """Return row order and supported metadata or trait-pattern branch labels."""
    linkage_matrix = spearman_row_linkage(matrix)
    stats = contrast_row_stats(cells)
    if linkage_matrix is None:
        row_order = list(matrix.index)
        branch_tests = pd.DataFrame()
        pattern_tests = pd.DataFrame()
        selected_branches: list[dict[str, object]] = []
    else:
        dendrogram_data = dendrogram(linkage_matrix, labels=list(matrix.index), no_plot=True)
        row_order = [str(label) for label in dendrogram_data["ivl"]]
        row_count = len(matrix.index)
        leaves_by_node: dict[int, frozenset[int]] = {
            index: frozenset({index}) for index in range(row_count)
        }
        for offset, linkage_row in enumerate(linkage_matrix):
            node_id = row_count + offset
            left_id = int(linkage_row[0])
            right_id = int(linkage_row[1])
            leaves_by_node[node_id] = leaves_by_node[left_id] | leaves_by_node[right_id]

        metadata = stats.reindex(matrix.index)
        dimensions = {
            "anatomy": "body_site_group",
            "condition": "condition",
        }
        test_records: list[dict[str, object]] = []
        root_node_id = row_count + len(linkage_matrix) - 1
        for node_id, leaf_indices in leaves_by_node.items():
            branch_size = len(leaf_indices)
            if (
                node_id < row_count
                or node_id == root_node_id
                or branch_size < MIN_SIGNIFICANT_BRANCH_ROWS
            ):
                continue
            branch_positions = sorted(leaf_indices)
            outside_size = row_count - branch_size
            for dimension, column in dimensions.items():
                values = metadata[column].astype(str)
                totals = values.value_counts()
                branch_values = values.iloc[branch_positions]
                for label_value, total_label_count in totals.items():
                    if int(total_label_count) < MIN_SIGNIFICANT_BRANCH_LABEL_COUNT:
                        continue
                    in_branch_count = int(branch_values.eq(label_value).sum())
                    outside_label_count = int(total_label_count) - in_branch_count
                    result = fisher_exact(
                        [
                            [in_branch_count, branch_size - in_branch_count],
                            [outside_label_count, outside_size - outside_label_count],
                        ],
                        alternative="greater",
                    )
                    test_records.append(
                        {
                            "branch_node_id": node_id,
                            "label_dimension": dimension,
                            "label_value": label_value,
                            "branch_row_count": branch_size,
                            "in_branch_label_count": in_branch_count,
                            "outside_branch_label_count": outside_label_count,
                            "total_label_count": int(total_label_count),
                            "branch_label_purity": in_branch_count / branch_size,
                            "label_coverage": in_branch_count / int(total_label_count),
                            "fisher_odds_ratio": float(result.statistic),
                            "fisher_p_value": float(result.pvalue),
                        }
                    )

        branch_tests = pd.DataFrame.from_records(test_records)
        selected_branches = []
        leaf_positions = {
            leaf_index: position for position, leaf_index in enumerate(dendrogram_data["leaves"])
        }
        metadata_candidates: list[dict[str, object]] = []
        if not branch_tests.empty:
            branch_tests["fisher_q_value"] = benjamini_hochberg(branch_tests["fisher_p_value"])
            branch_tests["passes_bh"] = branch_tests["fisher_q_value"] <= Q_VALUE_CUTOFF
            branch_tests["passes_label_count"] = (
                branch_tests["in_branch_label_count"] >= MIN_SIGNIFICANT_BRANCH_LABEL_COUNT
            )
            branch_tests["passes_label_purity"] = (
                branch_tests["branch_label_purity"] >= MIN_SIGNIFICANT_BRANCH_LABEL_PURITY
            )
            branch_tests["passes_display_criteria"] = (
                branch_tests["passes_bh"]
                & branch_tests["passes_label_count"]
                & branch_tests["passes_label_purity"]
            )
            branch_tests["branch_selected_for_display"] = False
            branch_tests["included_in_branch_label"] = False

            qualifying = branch_tests[branch_tests["passes_display_criteria"]].copy()
            qualifying = qualifying.sort_values(
                [
                    "branch_node_id",
                    "label_dimension",
                    "fisher_q_value",
                    "branch_label_purity",
                    "label_coverage",
                    "label_value",
                ],
                ascending=[True, True, True, False, False, True],
            )
            best_by_dimension = qualifying.drop_duplicates(
                ["branch_node_id", "label_dimension"],
                keep="first",
            )
            for node_id, group in best_by_dimension.groupby("branch_node_id", sort=True):
                leaf_indices = leaves_by_node[int(node_id)]
                positions = sorted(leaf_positions[index] for index in leaf_indices)
                label_parts = []
                for row in group.sort_values("label_dimension").itertuples(index=False):
                    dimension_label = "Anatomy" if row.label_dimension == "anatomy" else "Condition"
                    label_parts.append(
                        f"{dimension_label}: {row.label_value} "
                        f"({row.in_branch_label_count}/{row.branch_row_count} rows, "
                        f"q={row.fisher_q_value:.2g})"
                    )
                metadata_candidates.append(
                    {
                        "branch_node_id": int(node_id),
                        "leaf_indices": leaf_indices,
                        "leaf_start": min(positions),
                        "leaf_end": max(positions),
                        "row_count": len(leaf_indices),
                        "cluster_label": "; ".join(label_parts),
                        "label_source": "metadata_enrichment",
                        "minimum_q_value": float(group["fisher_q_value"].min()),
                        "maximum_label_purity": float(group["branch_label_purity"].max()),
                        "metadata_test_indices": group.index.tolist(),
                        "pattern_test_indices": [],
                        "pattern_bh_selected_cell_count": 0,
                    }
                )

            selected_branches.extend(metadata_candidates)

        visible_cells = cells[
            cells["displayed_in_heatmap"] & cells["condition_body_site_pair"].isin(matrix.index)
        ].copy()
        visible_delta = pd.to_numeric(
            visible_cells["increased_minus_decreased_positive_fraction"],
            errors="coerce",
        )
        visible_cells["pattern_direction"] = np.where(
            visible_delta > 0,
            "increased",
            "decreased",
        )
        global_trait_cells = {
            str(trait_id): trait_cells
            for trait_id, trait_cells in visible_cells.groupby("trait_id", sort=True)
        }
        pattern_records: list[dict[str, object]] = []
        for node_id, leaf_indices in leaves_by_node.items():
            branch_size = len(leaf_indices)
            if (
                node_id < row_count
                or node_id == root_node_id
                or branch_size < MIN_SIGNIFICANT_BRANCH_ROWS
            ):
                continue
            branch_labels = {str(matrix.index[index]) for index in leaf_indices}
            branch_cells = visible_cells[
                visible_cells["condition_body_site_pair"].isin(branch_labels)
            ]
            for trait_id, trait_cells in branch_cells.groupby("trait_id", sort=True):
                all_trait_cells = global_trait_cells[str(trait_id)]
                displayed_count = len(trait_cells)
                outside_displayed_count = len(all_trait_cells) - displayed_count
                for direction in ("increased", "decreased"):
                    direction_cells = trait_cells[trait_cells["pattern_direction"].eq(direction)]
                    direction_count = len(direction_cells)
                    global_direction_count = int(
                        all_trait_cells["pattern_direction"].eq(direction).sum()
                    )
                    outside_direction_count = global_direction_count - direction_count
                    result = fisher_exact(
                        [
                            [direction_count, displayed_count - direction_count],
                            [
                                outside_direction_count,
                                outside_displayed_count - outside_direction_count,
                            ],
                        ],
                        alternative="greater",
                    )
                    bh_count = int(direction_cells["selected_after_bh"].sum())
                    nominal_count = int(direction_cells["nominal_p_value_selected"].sum())
                    pattern_records.append(
                        {
                            "branch_node_id": node_id,
                            "trait_id": trait_id,
                            "trait_label": str(trait_cells["trait_label"].iloc[0]),
                            "pattern_direction": direction,
                            "branch_row_count": branch_size,
                            "displayed_trait_row_count": displayed_count,
                            "direction_row_count": direction_count,
                            "direction_purity": direction_count / displayed_count,
                            "outside_displayed_trait_row_count": outside_displayed_count,
                            "outside_direction_row_count": outside_direction_count,
                            "fisher_odds_ratio": float(result.statistic),
                            "fisher_p_value": float(result.pvalue),
                            "bh_selected_direction_cell_count": bh_count,
                            "nominal_direction_cell_count": nominal_count,
                            "mean_signed_effect": float(
                                pd.to_numeric(
                                    direction_cells["increased_minus_decreased_positive_fraction"],
                                    errors="coerce",
                                ).mean()
                            ),
                        }
                    )

        pattern_tests = pd.DataFrame.from_records(pattern_records)
        if not pattern_tests.empty:
            pattern_tests["fisher_q_value"] = benjamini_hochberg(pattern_tests["fisher_p_value"])
            pattern_tests["passes_bh"] = pattern_tests["fisher_q_value"] <= Q_VALUE_CUTOFF
            pattern_tests["passes_direction_rows"] = (
                pattern_tests["direction_row_count"] >= MIN_PATTERN_DIRECTION_ROWS
            )
            pattern_tests["passes_direction_purity"] = (
                pattern_tests["direction_purity"] >= MIN_PATTERN_DIRECTION_PURITY
            )
            pattern_tests["passes_bh_cell_support"] = (
                pattern_tests["bh_selected_direction_cell_count"] >= MIN_PATTERN_BH_SELECTED_CELLS
            )
            pattern_tests["passes_display_criteria"] = (
                pattern_tests["passes_bh"]
                & pattern_tests["passes_direction_rows"]
                & pattern_tests["passes_direction_purity"]
                & pattern_tests["passes_bh_cell_support"]
            )
            pattern_tests["branch_selected_for_display"] = False
            pattern_tests["included_in_branch_label"] = False
            qualifying_patterns = pattern_tests[pattern_tests["passes_display_criteria"]].copy()
            qualifying_patterns["absolute_mean_effect"] = qualifying_patterns[
                "mean_signed_effect"
            ].abs()
            qualifying_patterns = qualifying_patterns.sort_values(
                [
                    "branch_node_id",
                    "fisher_q_value",
                    "bh_selected_direction_cell_count",
                    "direction_purity",
                    "direction_row_count",
                    "absolute_mean_effect",
                    "trait_label",
                ],
                ascending=[True, True, False, False, False, False, True],
            )
            pattern_candidates: list[dict[str, object]] = []
            for node_id, group in qualifying_patterns.groupby("branch_node_id", sort=True):
                leaf_indices = leaves_by_node[int(node_id)]
                positions = sorted(leaf_positions[index] for index in leaf_indices)
                label_parts = []
                for row in group.itertuples(index=False):
                    direction_title = (
                        "Increased-associated"
                        if row.pattern_direction == "increased"
                        else "Decreased-associated"
                    )
                    cell_word = "cell" if row.bh_selected_direction_cell_count == 1 else "cells"
                    label_parts.append(
                        f"{direction_title} {row.trait_label} "
                        f"({row.direction_row_count}/{row.displayed_trait_row_count} rows; "
                        f"q={row.fisher_q_value:.2g}; "
                        f"{row.bh_selected_direction_cell_count} BH {cell_word})"
                    )
                pattern_candidates.append(
                    {
                        "branch_node_id": int(node_id),
                        "leaf_indices": leaf_indices,
                        "leaf_start": min(positions),
                        "leaf_end": max(positions),
                        "row_count": len(leaf_indices),
                        "cluster_label": "Pattern: " + "; ".join(label_parts),
                        "label_source": "bh_supported_trait_pattern",
                        "minimum_q_value": float(group["fisher_q_value"].min()),
                        "maximum_label_purity": float(group["direction_purity"].max()),
                        "metadata_test_indices": [],
                        "pattern_test_indices": group.index.tolist(),
                        "pattern_bh_selected_cell_count": int(
                            group["bh_selected_direction_cell_count"].sum()
                        ),
                    }
                )

            selected_branches.extend(pattern_candidates)

        candidates_by_node: dict[int, dict[str, object]] = {}
        for candidate in selected_branches:
            node_id = int(candidate["branch_node_id"])
            if node_id not in candidates_by_node:
                candidates_by_node[node_id] = candidate.copy()
                continue
            combined = candidates_by_node[node_id]
            combined["cluster_label"] = f"{combined['cluster_label']}; {candidate['cluster_label']}"
            combined["label_source"] = "+".join(
                sorted(
                    {
                        str(combined["label_source"]),
                        str(candidate["label_source"]),
                    }
                )
            )
            combined["minimum_q_value"] = min(
                float(combined["minimum_q_value"]),
                float(candidate["minimum_q_value"]),
            )
            combined["maximum_label_purity"] = max(
                float(combined["maximum_label_purity"]),
                float(candidate["maximum_label_purity"]),
            )
            combined["metadata_test_indices"] = list(combined["metadata_test_indices"]) + list(
                candidate["metadata_test_indices"]
            )
            combined["pattern_test_indices"] = list(combined["pattern_test_indices"]) + list(
                candidate["pattern_test_indices"]
            )
            combined["pattern_bh_selected_cell_count"] = int(
                combined["pattern_bh_selected_cell_count"]
            ) + int(candidate["pattern_bh_selected_cell_count"])

        selected_branches = list(candidates_by_node.values())
        for branch in selected_branches:
            branch_leaves = set(branch["leaf_indices"])
            branch["nesting_depth"] = sum(
                branch_leaves < set(other["leaf_indices"])
                for other in selected_branches
                if other is not branch
            )
        selected_branches.sort(
            key=lambda row: (
                -int(row["row_count"]),
                int(row["leaf_start"]),
                int(row["branch_node_id"]),
            )
        )
        for cluster_id, branch in enumerate(selected_branches, start=1):
            branch["spearman_cluster"] = cluster_id
            if branch["metadata_test_indices"]:
                in_branch = branch_tests["branch_node_id"].eq(branch["branch_node_id"])
                branch_tests.loc[in_branch, "branch_selected_for_display"] = True
                branch_tests.loc[branch["metadata_test_indices"], "included_in_branch_label"] = True
            if branch["pattern_test_indices"]:
                in_branch = pattern_tests["branch_node_id"].eq(branch["branch_node_id"])
                pattern_tests.loc[in_branch, "branch_selected_for_display"] = True
                pattern_tests.loc[branch["pattern_test_indices"], "included_in_branch_label"] = True

    ordered_stats = stats.reindex(row_order)
    row_records = ordered_stats.reset_index().rename(
        columns={"condition_body_site_pair": "row_label"}
    )
    row_records.insert(0, "row_order", np.arange(1, len(row_records) + 1))
    row_records.insert(
        1, "spearman_cluster", pd.Series(pd.NA, index=row_records.index, dtype="Int64")
    )
    row_records.insert(
        2, "branch_node_id", pd.Series(pd.NA, index=row_records.index, dtype="Int64")
    )
    row_records["cluster_label"] = pd.NA
    memberships: dict[str, list[dict[str, object]]] = {
        str(label): [] for label in row_records["row_label"]
    }
    summary_records: list[dict[str, object]] = []
    original_labels = list(matrix.index)
    for branch in selected_branches:
        cluster_id = int(branch["spearman_cluster"])
        labels = [original_labels[index] for index in branch["leaf_indices"]]
        in_branch = row_records["row_label"].isin(labels)
        row_records.loc[in_branch, "spearman_cluster"] = cluster_id
        row_records.loc[in_branch, "branch_node_id"] = int(branch["branch_node_id"])
        row_records.loc[in_branch, "cluster_label"] = branch["cluster_label"]
        for label in labels:
            memberships[str(label)].append(branch)
        group = row_records.loc[in_branch]
        dominant_body_site = group["body_site_group"].astype(str).value_counts().index[0]
        dominant_condition = group["condition"].astype(str).value_counts().index[0]
        metadata_indices = branch["metadata_test_indices"]
        if metadata_indices:
            label_tests = branch_tests.loc[metadata_indices].sort_values("label_dimension")
        else:
            label_tests = pd.DataFrame(columns=["label_dimension"])
        anatomy_test = label_tests[label_tests["label_dimension"].eq("anatomy")]
        condition_test = label_tests[label_tests["label_dimension"].eq("condition")]
        pattern_indices = branch["pattern_test_indices"]
        if pattern_indices:
            selected_patterns = pattern_tests.loc[pattern_indices]
        else:
            selected_patterns = pd.DataFrame()
        if metadata_indices and pattern_indices:
            cluster_label_rule = (
                "all qualifying metadata and signed-trait titles are shown on nested tracks; "
                "metadata uses one-sided Fisher branch enrichment with global BH q<=0.05, "
                ">=3 rows, >=2 matching rows, and purity >=2/3; trait directions use "
                "one-sided Fisher enrichment versus outside displayed rows with global "
                "branch-pattern BH q<=0.05, >=3 matching rows, purity >=2/3, and >=1 "
                "same-direction heatmap-wide BH-selected cell"
            )
        elif metadata_indices:
            cluster_label_rule = (
                "one-sided Fisher branch enrichment; BH q<=0.05; >=3 rows; >=2 matching "
                "rows; label purity >=2/3; all qualifying nested branches shown"
            )
        else:
            cluster_label_rule = (
                "one-sided Fisher enrichment of a signed trait direction versus displayed "
                "rows outside the branch; global branch-pattern BH q<=0.05; >=3 "
                "direction-matching displayed rows; direction purity >=2/3; >=1 "
                "same-direction cell selected at heatmap-wide BH q<=0.05; all qualifying "
                "nested branches shown"
            )
        minimum_q_value = float(branch["minimum_q_value"])
        summary_records.append(
            {
                "spearman_cluster": cluster_id,
                "branch_node_id": int(branch["branch_node_id"]),
                "leaf_start": int(branch["leaf_start"]),
                "leaf_end": int(branch["leaf_end"]),
                "nesting_depth": int(branch["nesting_depth"]),
                "cluster_label": branch["cluster_label"],
                "cluster_label_supported": True,
                "cluster_label_source": branch["label_source"],
                "cluster_label_rule": cluster_label_rule,
                "row_count": len(group),
                "displayed_cell_count": int(group["displayed_trait_count"].sum()),
                "nominal_p_value_cell_count": int(group["significant_trait_count"].sum()),
                "bh_selected_cell_count": int(group["bh_selected_trait_count"].sum()),
                "minimum_branch_label_q_value": (
                    minimum_q_value if math.isfinite(minimum_q_value) else pd.NA
                ),
                "maximum_branch_label_purity": branch["maximum_label_purity"],
                "pattern_bh_selected_cell_count": branch["pattern_bh_selected_cell_count"],
                "metadata_label_count": len(label_tests),
                "trait_pattern_count": len(selected_patterns),
                "pattern_traits": (
                    "; ".join(selected_patterns["trait_label"].astype(str))
                    if not selected_patterns.empty
                    else pd.NA
                ),
                "pattern_directions": (
                    "; ".join(selected_patterns["pattern_direction"].astype(str))
                    if not selected_patterns.empty
                    else pd.NA
                ),
                "anatomy_label": (
                    anatomy_test.iloc[0]["label_value"] if not anatomy_test.empty else pd.NA
                ),
                "anatomy_label_count": (
                    int(anatomy_test.iloc[0]["in_branch_label_count"])
                    if not anatomy_test.empty
                    else pd.NA
                ),
                "anatomy_label_q_value": (
                    anatomy_test.iloc[0]["fisher_q_value"] if not anatomy_test.empty else pd.NA
                ),
                "condition_label": (
                    condition_test.iloc[0]["label_value"] if not condition_test.empty else pd.NA
                ),
                "condition_label_count": (
                    int(condition_test.iloc[0]["in_branch_label_count"])
                    if not condition_test.empty
                    else pd.NA
                ),
                "condition_label_q_value": (
                    condition_test.iloc[0]["fisher_q_value"] if not condition_test.empty else pd.NA
                ),
                "dominant_body_site_group": dominant_body_site,
                "dominant_grouped_condition": dominant_condition,
                "top_body_site_groups": top_values(group["body_site_group"]),
                "top_grouped_conditions": top_values(group["condition"]),
            }
        )
    row_records["supported_branch_count"] = row_records["row_label"].map(
        lambda label: len(memberships[str(label)])
    )
    row_records["supported_branch_node_ids"] = row_records["row_label"].map(
        lambda label: (
            "; ".join(str(branch["branch_node_id"]) for branch in memberships[str(label)]) or pd.NA
        )
    )
    row_records["supported_branch_labels"] = row_records["row_label"].map(
        lambda label: (
            " | ".join(str(branch["cluster_label"]) for branch in memberships[str(label)]) or pd.NA
        )
    )
    summary = pd.DataFrame.from_records(summary_records)
    return row_records, summary, branch_tests, pattern_tests, linkage_matrix


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


def flag_matrix(
    cells: pd.DataFrame,
    row_labels: list[str],
    trait_order: pd.DataFrame,
    *,
    flag_column: str,
) -> pd.DataFrame:
    """Pivot a boolean cell flag into row x trait layout."""
    trait_ids = trait_order["trait_id"].tolist()
    flags = cells[cells[flag_column]].copy()
    matrix = flags.pivot_table(
        index="condition_body_site_pair",
        columns="trait_id",
        values=flag_column,
        aggfunc="max",
    )
    return matrix.reindex(index=row_labels, columns=trait_ids).fillna(False).astype(bool)


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
        p_text = format_p_value(row.fisher_p_value)
        p_label = f"p{p_text}" if p_text.startswith("<") else f"p={p_text}"
        if mode == "pattern":
            annotations.loc[row.condition_body_site_pair, row.trait_id] = (
                f"n={row.focus_positive_profile_count}\n{p_label}"
            )
        else:
            delta = row.increased_minus_decreased_positive_fraction
            annotations.loc[row.condition_body_site_pair, row.trait_id] = f"{delta:+.2f}"
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
    compact: bool = False,
    border_mask: pd.DataFrame | None = None,
    contrast_legend: ContrastLegendSummary | None = None,
) -> None:
    """Render one row x six-trait heatmap with blank nonsignificant cells."""
    row_count = max(1, len(matrix.index))
    fig_height = max(6.2, 2.8 + 0.36 * row_count)
    fig, ax = plt.subplots(figsize=(13.8, fig_height))
    fig.subplots_adjust(left=0.34, right=0.90, top=0.82, bottom=0.18)
    if compact:
        fig.subplots_adjust(top=1 - 1.15 / fig_height, bottom=1.2 / fig_height)
    if contrast_legend is not None:
        fig.subplots_adjust(bottom=max(fig.subplotpars.bottom, 1.9 / fig_height))
    fig.patch.set_facecolor("white")
    ax.set_facecolor(SLATE_100)

    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(SLATE_100)
    values = (
        matrix.to_numpy(dtype=float)
        if len(matrix.index)
        else np.full((1, len(trait_order)), np.nan)
    )
    image = ax.imshow(values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    if matrix.empty:
        ax.text(0.5, 0.5, "No qualifying contrasts", transform=ax.transAxes, ha="center")

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
            if border_mask is not None and bool(border_mask.loc[row_label, trait_id]):
                ax.add_patch(
                    Rectangle(
                        (x_index - 0.5, y_index - 0.5),
                        1,
                        1,
                        fill=False,
                        edgecolor=SIGNIFICANCE_BORDER,
                        linewidth=1.8,
                    )
                )
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
        1 - 0.12 / fig_height if compact else 0.965,
        title,
        ha="left",
        va="top",
        fontsize=15.5,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.04,
        1 - 0.45 / fig_height if compact else 0.925,
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
    if contrast_legend is not None:
        add_contrast_legend(fig, contrast_legend, clustered=False)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_grouped_clustered_contrast_heatmap(
    matrix: pd.DataFrame,
    annotations: pd.DataFrame,
    border_mask: pd.DataFrame,
    trait_order: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    linkage_matrix: np.ndarray | None,
    output_path: Path,
    *,
    legend: ContrastLegendSummary,
    preview: bool = False,
) -> None:
    """Render the grouped signed direction heatmap with a left Spearman dendrogram."""
    row_count = max(1, len(matrix.index))
    significant_branch_count = len(cluster_summary)
    significant_metadata_label_count = int(
        cluster_summary.get("metadata_label_count", pd.Series(dtype=int)).sum()
    )
    significant_trait_pattern_count = int(
        cluster_summary.get("trait_pattern_count", pd.Series(dtype=int)).sum()
    )
    label_track_count = (
        int(cluster_summary["nesting_depth"].max()) + 1 if not cluster_summary.empty else 1
    )
    if preview:
        fig_height = 8.5
        fig_width = 13.8
        left = 0.08
        right = 0.92
        top = 0.86
        bottom = 0.23
        label_width = 0.01
        dendrogram_width = 4.2
        show_row_labels = False
        show_cell_text = False
    else:
        fig_height = max(9.5, 3.2 + 0.34 * row_count)
        fig_width = 22.0
        left = 0.04
        right = 0.965
        top = 1 - 1.55 / fig_height
        bottom = 1.8 / fig_height
        label_width = 3.2
        dendrogram_width = 4.8
        show_row_labels = True
        show_cell_text = True

    fig = plt.figure(figsize=(fig_width, fig_height), facecolor="white")
    grid = fig.add_gridspec(
        1,
        4,
        width_ratios=[dendrogram_width, 6.0, label_width, 0.32],
        left=left,
        right=right,
        top=top,
        bottom=bottom,
        wspace=0.035,
    )
    dendro_ax = fig.add_subplot(grid[0, 0])
    ax = fig.add_subplot(grid[0, 1])
    row_label_spacer_ax = fig.add_subplot(grid[0, 2])
    cbar_ax = fig.add_subplot(grid[0, 3])

    y_extent = row_count * 10
    branch_node_x: dict[int, float] = {}
    branch_node_y: dict[int, float] = {}
    if linkage_matrix is not None:
        dendrogram_data = dendrogram(
            linkage_matrix,
            orientation="left",
            no_labels=True,
            color_threshold=0,
            above_threshold_color=SLATE_500,
            ax=dendro_ax,
        )
        dendro_ax.set_ylim(y_extent, 0)
        branch_node_y.update(
            {
                int(leaf_index): position * 10 + 5
                for position, leaf_index in enumerate(dendrogram_data["leaves"])
            }
        )
        for offset, linkage_row in enumerate(linkage_matrix):
            node_id = row_count + offset
            left_node = int(linkage_row[0])
            right_node = int(linkage_row[1])
            branch_node_x[node_id] = float(linkage_row[2])
            branch_node_y[node_id] = (branch_node_y[left_node] + branch_node_y[right_node]) / 2

    branches_to_plot = (
        cluster_summary.sort_values(["nesting_depth", "leaf_start", "branch_node_id"])
        if not cluster_summary.empty
        else cluster_summary
    )
    for branch in branches_to_plot.itertuples(index=False):
        track = int(branch.nesting_depth)
        node_id = int(branch.branch_node_id)
        node_y = branch_node_y.get(
            node_id,
            (int(branch.leaf_start) + int(branch.leaf_end) + 1) * 5,
        )
        node_x = branch_node_x.get(node_id, 0.0)
        track_spacing = 0.88 / max(1, label_track_count)
        label_x = 0.03 + track * track_spacing
        label_text = compact_dendrogram_label(
            branch.cluster_label,
            width=27 if not preview else 23,
        )
        dendro_ax.annotate(
            label_text,
            xy=(node_x, node_y),
            xycoords="data",
            xytext=(label_x, node_y),
            textcoords=dendro_ax.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=8.0 if not preview else 6.7,
            color="black",
            fontweight="bold",
            bbox={
                "boxstyle": "round,pad=0.14",
                "facecolor": "white",
                "edgecolor": "black",
                "linewidth": 1.2,
                "alpha": 0.96,
            },
            arrowprops={
                "arrowstyle": "-",
                "color": "black",
                "linewidth": 0.9,
                "shrinkA": 2,
                "shrinkB": 1,
            },
            annotation_clip=False,
            clip_on=False,
            zorder=5,
        )
    dendro_ax.axis("off")
    row_label_spacer_ax.patch.set_alpha(0)
    row_label_spacer_ax.set_zorder(-1)
    row_label_spacer_ax.axis("off")

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(SLATE_100)
    values = np.ma.masked_invalid(matrix.to_numpy(dtype=float))
    image = ax.imshow(
        values,
        origin="lower",
        aspect="auto",
        cmap=cmap,
        vmin=-1.0,
        vmax=1.0,
        extent=[-0.5, len(trait_order) - 0.5, 0, y_extent],
    )
    ax.set_ylim(y_extent, 0)
    ax.set_xticks(
        np.arange(len(trait_order)),
        [wrap_label(label, width=16) for label in trait_order["trait_label"]],
    )
    ax.tick_params(axis="x", labelrotation=35, labelsize=8.5 if not preview else 7.2)
    if show_row_labels:
        ax.set_yticks(
            np.arange(len(matrix.index)) * 10 + 5,
            [wrap_label(label, width=48) for label in matrix.index],
        )
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.tick_params(
            axis="y",
            which="major",
            left=False,
            right=False,
            labelleft=False,
            labelright=True,
            labelsize=7.3,
            pad=7,
        )
    else:
        ax.set_yticks([])
    ax.set_xlabel("Six shared traits", color=SLATE_700)
    ax.set_ylabel(
        "" if preview else "Spearman-clustered grouped condition | anatomy rows",
        color=SLATE_700,
        labelpad=18,
    )
    ax.set_xticks(np.arange(-0.5, len(trait_order), 1), minor=True)
    ax.set_yticks(np.arange(0, y_extent + 1, 10), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)
    ax.tick_params(which="minor", bottom=False, left=False, right=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if show_cell_text:
        for y_index, row_label in enumerate(matrix.index):
            y_value = y_index * 10 + 5
            for x_index, trait_id in enumerate(trait_order["trait_id"]):
                if bool(border_mask.loc[row_label, trait_id]):
                    ax.add_patch(
                        Rectangle(
                            (x_index - 0.5, y_index * 10),
                            1,
                            10,
                            fill=False,
                            edgecolor=SIGNIFICANCE_BORDER,
                            linewidth=1.9,
                        )
                    )
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
                    y_value,
                    text,
                    ha="center",
                    va="center",
                    fontsize=6.7,
                    color=text_color,
                )
    else:
        for y_index, row_label in enumerate(matrix.index):
            for x_index, trait_id in enumerate(trait_order["trait_id"]):
                if bool(border_mask.loc[row_label, trait_id]):
                    ax.add_patch(
                        Rectangle(
                            (x_index - 0.5, y_index * 10),
                            1,
                            10,
                            fill=False,
                            edgecolor=SIGNIFICANCE_BORDER,
                            linewidth=0.85,
                        )
                    )

    title = "Grouped signed direction contrast: Spearman-clustered condition/anatomy rows"
    subtitle = (
        "Rows are fixed upstream grouped condition/body-site pairs with at least one displayed "
        f"cell (|effect|>{SIGNED_EFFECT_ABS_CUTOFF:g}). The dendrogram orders rows using average "
        "linkage on Spearman-correlation distance across signed trait deltas; missing cells are "
        "neutral zero for clustering. A boxed title requires either BH-significant metadata "
        "enrichment or BH-significant enrichment of a recurring signed-trait direction versus "
        "the rest of the dendrogram (>=3 rows, >=67% agreement, and >=1 same-direction "
        "BH-selected cell). Titles report their support; unlabeled branches are not interpreted. "
        "Dark cell borders mark BH q<=0.05."
    )
    fig.text(
        0.04,
        0.985 if preview else 1 - 0.12 / fig_height,
        title if not preview else "Grouped signed direction contrast preview",
        ha="left",
        va="top",
        fontsize=15.5 if not preview else 13.0,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.04,
        0.94 if preview else 1 - 0.52 / fig_height,
        "\n".join(wrap(subtitle, width=145 if not preview else 105)),
        ha="left",
        va="top",
        fontsize=8.7 if not preview else 7.6,
        color=SLATE_500,
    )
    cbar = fig.colorbar(image, cax=cbar_ax)
    cbar.set_label("Increased - decreased positive-z fraction", color=SLATE_700)
    cbar.ax.tick_params(labelsize=7.5)
    cbar.outline.set_edgecolor(SLATE_200)
    add_contrast_legend(
        fig,
        legend,
        clustered=True,
        significant_branch_count=significant_branch_count,
        significant_metadata_label_count=significant_metadata_label_count,
        significant_trait_pattern_count=significant_trait_pattern_count,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path, dpi=220 if not preview else 180, bbox_inches="tight", facecolor="white"
    )
    plt.close(fig)


def write_pattern_heatmap(
    cells: pd.DataFrame,
    view: pd.Series,
    trait_order: pd.DataFrame,
    output_dir: Path,
) -> tuple[Path, int, int, int]:
    """Write one approach-2 condition/body-site x trait pattern heatmap."""
    row_labels = select_pattern_rows(cells)
    matrix = cell_matrix(
        cells,
        row_labels,
        trait_order,
        value_column="plot_minus_log10_fisher_p_value",
    )
    annotations = annotation_matrix(cells, row_labels, trait_order, mode="pattern")
    output_path = output_dir / (
        f"six_shared_traits_{view.view_key}_condition_body_site_trait_pattern_heatmap.png"
    )
    displayed_cells = int(cells["displayed_in_heatmap"].sum())
    subtitle = (
        "Approach 2: rows are condition/body-site pairs and columns are the six "
        "traits. Cell color is Fisher evidence, -log10(p), so darker cells have "
        "smaller p-values. Cells are shown for exploratory review only when Fisher "
        f"p<{P_VALUE_CUTOFF:g}, n>=3, and s>=2; BH q-values are written to the TSV."
    )
    plot_matrix_heatmap(
        matrix,
        annotations,
        trait_order,
        output_path,
        title=f"Condition/body-site trait pattern: {view.view_label}",
        subtitle=subtitle,
        colorbar_label=f"Fisher evidence (-log10 p, capped at {MAX_MINUS_LOG10_P:g})",
        cmap_name="YlGnBu",
        vmin=0.0,
        vmax=MAX_MINUS_LOG10_P,
    )
    return output_path, len(row_labels), displayed_cells, len(cells)


def write_contrast_heatmap(
    cells: pd.DataFrame,
    mode: ContrastMode,
    trait_order: pd.DataFrame,
    output_dir: Path,
) -> tuple[Path, int, int, int, dict[str, object]]:
    """Write one approach-4 signed direction contrast heatmap."""
    row_labels = select_contrast_rows(cells)
    matrix = cell_matrix(
        cells,
        row_labels,
        trait_order,
        value_column="increased_minus_decreased_positive_fraction",
    )
    annotations = annotation_matrix(cells, row_labels, trait_order, mode="contrast")
    border_mask = flag_matrix(
        cells,
        row_labels,
        trait_order,
        flag_column="selected_after_bh",
    )
    output_path = (
        output_dir
        / f"six_shared_traits_{mode.key}_condition_body_site_signed_direction_contrast_heatmap.png"
    )
    displayed_cells = int(cells["displayed_in_heatmap"].sum())
    nominal_p_value_cells = int(cells["nominal_p_value_selected"].sum())
    legend = summarize_contrast_legend(cells, row_labels, mode)
    if mode.key == "grouped":
        (
            row_clusters,
            cluster_summary,
            branch_tests,
            pattern_tests,
            row_linkage,
        ) = clustered_contrast_rows(matrix, cells)
        ordered_labels = row_clusters["row_label"].astype(str).tolist()
        matrix = matrix.reindex(ordered_labels)
        annotations = annotations.reindex(ordered_labels)
        border_mask = border_mask.reindex(ordered_labels).fillna(False).astype(bool)

        cluster_by_row = row_clusters.set_index("row_label")["spearman_cluster"].to_dict()
        cluster_label_by_row = row_clusters.set_index("row_label")["cluster_label"].to_dict()
        in_rows = cells["condition_body_site_pair"].map(cluster_by_row)
        cells["spearman_cluster"] = in_rows.astype("Int64")
        cells["spearman_cluster_label"] = cells["condition_body_site_pair"].map(
            cluster_label_by_row
        )
        branch_by_row = row_clusters.set_index("row_label")["branch_node_id"].to_dict()
        cells["spearman_branch_node_id"] = (
            cells["condition_body_site_pair"].map(branch_by_row).astype("Int64")
        )
        order_by_row = row_clusters.set_index("row_label")["row_order"].to_dict()
        cells["spearman_cluster_row_order"] = (
            cells["condition_body_site_pair"].map(order_by_row).astype("Int64")
        )

        cluster_rows_path = output_dir / GROUPED_CONTRAST_CLUSTER_ROWS_FILENAME
        cluster_summary_path = output_dir / GROUPED_CONTRAST_CLUSTER_SUMMARY_FILENAME
        branch_tests_path = output_dir / GROUPED_CONTRAST_BRANCH_TESTS_FILENAME
        pattern_tests_path = output_dir / GROUPED_CONTRAST_PATTERN_TESTS_FILENAME
        write_tsv(row_clusters, cluster_rows_path)
        write_tsv(cluster_summary, cluster_summary_path)
        write_tsv(branch_tests, branch_tests_path)
        write_tsv(pattern_tests, pattern_tests_path)

        preview_path = output_path.with_stem(f"{output_path.stem}_preview")
        plot_grouped_clustered_contrast_heatmap(
            matrix,
            annotations,
            border_mask,
            trait_order,
            cluster_summary,
            row_linkage,
            output_path,
            legend=legend,
        )
        plot_grouped_clustered_contrast_heatmap(
            matrix,
            annotations,
            border_mask,
            trait_order,
            cluster_summary,
            row_linkage,
            preview_path,
            legend=legend,
            preview=True,
        )
        for page_path in output_dir.glob(f"{output_path.stem}_page_*.png"):
            page_path.unlink()
        extra = {
            "display_rule": (
                "supported testable cells with absolute effect above cutoff; BH q cutoff marks "
                "dark borders; nominal p cutoff is audit-only"
            ),
            "page_count": 0,
            "preview_png_path": preview_path,
            "cluster_rows_path": cluster_rows_path,
            "cluster_summary_path": cluster_summary_path,
            "branch_tests_path": branch_tests_path,
            "pattern_tests_path": pattern_tests_path,
            "clustering_method": "row average linkage on Spearman-correlation distance",
            "cluster_label_rule": (
                "metadata titles require one-sided Fisher branch enrichment with global BH "
                "q<=0.05, >=3 branch rows, >=2 matching rows, and purity >=2/3; descriptive "
                "trait-pattern titles require one-sided Fisher enrichment versus outside rows "
                "with global branch-pattern BH q<=0.05, >=3 same-direction displayed rows, "
                "direction purity >=2/3, and >=1 same-direction cell selected at heatmap-wide "
                "BH q<=0.05; every qualifying nested branch and trait pattern displayed"
            ),
            "significant_branch_count": len(cluster_summary),
            "metadata_enrichment_branch_count": int(
                cluster_summary.get("cluster_label_source", pd.Series(dtype=str))
                .astype(str)
                .str.contains("metadata_enrichment", regex=False)
                .sum()
            ),
            "trait_pattern_branch_count": int(
                cluster_summary.get("cluster_label_source", pd.Series(dtype=str))
                .astype(str)
                .str.contains("bh_supported_trait_pattern", regex=False)
                .sum()
            ),
            "metadata_label_count": int(
                cluster_summary.get("metadata_label_count", pd.Series(dtype=int)).sum()
            ),
            "trait_pattern_label_count": int(
                cluster_summary.get("trait_pattern_count", pd.Series(dtype=int)).sum()
            ),
            "legend": legend,
        }
        return output_path, len(row_labels), displayed_cells, len(cells), extra

    subtitle = (
        "Color and numbers show increased minus decreased positive-z fractions. "
        "Red favors increased profiles; blue favors decreased profiles. "
        f"Dark borders mark BH q<={Q_VALUE_CUTOFF:g}; {nominal_p_value_cells} cells meet "
        f"nominal Fisher p<{P_VALUE_CUTOFF:g}. Shown cells require "
        f"|effect|>{SIGNED_EFFECT_ABS_CUTOFF:g}, >=3 positive-z profiles in either direction, "
        "and >=2 total studies. Gray = unavailable, insufficient support, or below cutoff. "
        "Raw p and BH q are in the TSV."
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
        compact=True,
        border_mask=border_mask,
        contrast_legend=legend,
    )
    pages = contrast_row_pages(row_labels)
    for page, labels in enumerate(pages, start=1):
        page_legend = summarize_contrast_legend(cells, labels, mode)
        plot_matrix_heatmap(
            matrix.loc[labels],
            annotations.loc[labels],
            trait_order,
            contrast_page_path(output_path, page),
            title=f"Signed direction contrast: {mode.label} ({page}/{len(pages)})",
            subtitle=subtitle,
            colorbar_label="Increased - decreased positive-z fraction",
            cmap_name="RdBu_r",
            vmin=-1.0,
            vmax=1.0,
            compact=True,
            border_mask=border_mask.loc[labels],
            contrast_legend=page_legend,
        )
    extra = {
        "display_rule": (
            "supported testable cells with absolute effect above cutoff; BH q cutoff marks "
            "dark borders; nominal p cutoff is audit-only"
        ),
        "page_count": len(pages),
        "preview_png_path": "",
        "cluster_rows_path": "",
        "cluster_summary_path": "",
        "branch_tests_path": "",
        "pattern_tests_path": "",
        "clustering_method": "",
        "cluster_label_rule": "not applicable; rows are ranked, not clustered",
        "significant_branch_count": 0,
        "metadata_enrichment_branch_count": 0,
        "trait_pattern_branch_count": 0,
        "metadata_label_count": 0,
        "trait_pattern_label_count": 0,
        "legend": legend,
    }
    return output_path, len(row_labels), displayed_cells, len(cells), extra


def contrast_row_pages(row_labels: list[str]) -> list[list[str]]:
    """Partition the complete row order without dropping or duplicating rows."""
    return [
        row_labels[start : start + CONTRAST_ROWS_PER_PAGE]
        for start in range(0, len(row_labels), CONTRAST_ROWS_PER_PAGE)
    ]


def contrast_page_path(output_path: Path, page: int) -> Path:
    """Return a stable path for a readable page of the full contrast overview."""
    return output_path.with_stem(f"{output_path.stem}_page_{page:02d}")


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
        f"- Grouped signed-direction cluster rows: `{GROUPED_CONTRAST_CLUSTER_ROWS_FILENAME}`",
        "- Grouped signed-direction cluster summary: "
        f"`{GROUPED_CONTRAST_CLUSTER_SUMMARY_FILENAME}`",
        f"- Grouped dendrogram branch-label tests: `{GROUPED_CONTRAST_BRANCH_TESTS_FILENAME}`",
        f"- Grouped dendrogram trait-pattern tests: `{GROUPED_CONTRAST_PATTERN_TESTS_FILENAME}`",
        f"- BH-selected cell table: `{SELECTED_CELLS_FILENAME}`",
        "",
        "Trait-pattern heatmaps use the six-shared-trait display-edge TSVs that feed "
        "the Cytoscape graph views. Signed contrasts use the full increased and decreased "
        "edge TSVs before graph display selection. They preserve condition/body-site/trait records "
        "instead of correlating already-collapsed condition and body-site nodes.",
        "",
        "Signed contrasts display supported, testable effects whose absolute magnitude is "
        f"greater than {SIGNED_EFFECT_ABS_CUTOFF:g}, including nonsignificant results. Gray "
        "cells are unavailable, under-supported, or at or below the effect cutoff. "
        "Trait-pattern heatmaps continue to mask cells that fail the nominal p-value threshold.",
        "",
        "### How Conditions and Cluster Labels Are Chosen",
        "",
        "The signed heatmaps do not choose a top-N disease list. They read the full increased "
        "and decreased Sankey edge tables before graph display selection and include every "
        "condition/body-site row with at least one supported, finite cell above the absolute-"
        "effect cutoff. The grouped view uses the fixed upstream v3 clinical rollup; the "
        "individual view retains the original curated condition label. Both use the upstream "
        "grouped v3 anatomy label. These condition and anatomy names are not inferred from the "
        "dendrogram.",
        "",
        "Grouped conditions are assigned by the ordered, case-insensitive substring rules in "
        "`v3_fully_grouped_condition_display_label`: the first matching rule wins and unmatched "
        "labels become `Other specific clinical labels`. The per-trait "
        "`*_profile_sankey_condition_rollups.tsv` files audit each raw label, displayed label, "
        "rule, profile count, and study count. Every internal dendrogram branch is tested for "
        "one-sided enrichment of each eligible anatomy and grouped-condition label versus all "
        "rows outside that branch. BH adjustment covers the complete tested branch-by-label "
        "family. A branch is named only when q<=0.05, it has at least 3 rows and 2 matching "
        "rows, and the label describes at least two-thirds of the branch. Every qualifying "
        "nested label is placed directly inside the dendrogram and connected to its branch node. "
        "Each visible label reports its matching-row fraction and adjusted q-value; unlabeled "
        "branches are not interpreted.",
        "A branch without significant anatomy/condition enrichment may instead receive a "
        "descriptive trait-pattern title. Each branch/trait/direction is tested for one-sided "
        "enrichment against displayed rows outside the branch, with BH adjustment across the "
        "complete branch-pattern family. A displayed pattern requires q<=0.05, at least three "
        "same-direction displayed cells, at least two-thirds directional agreement, and at least "
        "one same-direction cell already selected at heatmap-wide BH q<=0.05. Pattern titles "
        "report the matching-row fraction, branch-pattern q-value, and BH-cell count. The "
        "pattern-test TSV records every test and selection flag.",
        "",
        "- `trait_pattern`: Approach 2. Rows are condition/body-site pairs and "
        "columns are the six traits. Cell color is Fisher evidence, `-log10(p)`, "
        "rendered with a yellow-to-blue sequential palette so darker cells have "
        f"smaller p-values. Cells require Fisher p<{P_VALUE_CUTOFF:g}, n>=3, "
        "and s>=2 for exploratory display.",
        "- `signed_direction_contrast`: Approach 4. Rows are condition/body-site "
        "pairs and columns are traits. Cell color and numbers are "
        "`increased_fraction - decreased_fraction` on a fixed -1 to +1 scale. Red "
        "favors increased profiles and blue favors decreased profiles. "
        f"A dark border marks BH-adjusted q<={Q_VALUE_CUTOFF:g}; the nominal Fisher "
        f"p<{P_VALUE_CUTOFF:g} flag is reported separately. Unbordered cells are descriptive "
        f"effects. Cells also require absolute effect >{SIGNED_EFFECT_ABS_CUTOFF:g}.",
        "- Condition/body-site rows are chosen from the full increased and decreased edge TSVs, "
        "not from the graph display subset. A row is retained when at least one trait cell has "
        "a finite increased-vs-decreased Fisher test, absolute effect above the cutoff, at "
        "least 3 positive-z profiles in either direction, and at least 2 total studies.",
        "- The grouped signed-direction heatmap is a single clustered diagram. Rows use average "
        "linkage on Spearman-correlation distance across the six signed trait deltas; missing "
        "or unsupported cells are neutral zero for clustering. The left dendrogram orders "
        "condition/anatomy rows. Only branches passing the multiplicity-adjusted enrichment, "
        "minimum-size, and minimum-purity rules receive an anatomy/condition title. Other branches "
        "can receive a descriptive trait-pattern title only under branch-pattern BH significance, "
        "BH-cell support, recurrence, and directional-purity rules. All qualifying nested labels "
        "are overlaid compactly inside the dendrogram and connected to their tested nodes; the "
        "two audit TSVs preserve every rejected and accepted metadata or trait-pattern candidate.",
        "- Individual-condition signed-direction heatmaps still receive readable pages containing "
        f"at most {CONTRAST_ROWS_PER_PAGE} rows each. `{CONTRAST_PAGES_FILENAME}` records only "
        "those individual-condition page paths and counts.",
        f"- BH q-values are computed within each heatmap view. "
        f"`selected_after_bh=True` marks supported cells with q<={Q_VALUE_CUTOFF:g}; visible "
        "selected cells receive the dark border and selected candidates are collected in the "
        "BH-selected cell table. `nominal_p_value_selected` is the explicit nominal flag; "
        "`significant_in_heatmap` remains a compatibility alias for that nominal flag.",
        "- Issue #256 / PR #257 adds pooled, denominator-harmonized, study-weighted fraction "
        "and log-ratio quantities, but no condition-by-anatomy replacement for this heatmap. "
        "Those pooled quantities are not mixed into this profile-level positive-z prevalence "
        "contrast because the estimands and sampling units differ.",
        "",
    ]
    grouped = manifest[manifest["view_key"].eq("grouped_direction_contrast")]
    if not grouped.empty and str(grouped.iloc[0].get("preview_png_path", "")):
        preview_path = Path(str(grouped.iloc[0]["preview_png_path"])).name
        full_path = Path(str(grouped.iloc[0]["png_path"])).name
        lines.extend(
            [
                "### Grouped Signed-Direction Cluster Preview",
                "",
                f"[![Grouped signed-direction cluster preview]({preview_path})]({full_path})",
                "",
            ]
        )
    lines.extend(
        [
            "| heatmap | type | PNG | displayed cells | nominal p<0.01 cells | "
            "BH-selected cells | tested cells | rows | condition labels | anatomy groups | "
            "supported dendrogram branches | metadata labels | trait-pattern labels |",
            "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | "
            "---: | ---: |",
        ]
    )
    for row in manifest.itertuples(index=False):
        nominal_count = (
            "" if pd.isna(row.nominal_p_value_cell_count) else int(row.nominal_p_value_cell_count)
        )
        condition_count = (
            "" if pd.isna(row.condition_label_count) else int(row.condition_label_count)
        )
        anatomy_count = "" if pd.isna(row.anatomy_group_count) else int(row.anatomy_group_count)
        branch_count = (
            "" if pd.isna(row.significant_branch_count) else int(row.significant_branch_count)
        )
        metadata_label_count = (
            "" if pd.isna(row.metadata_label_count) else int(row.metadata_label_count)
        )
        pattern_label_count = (
            "" if pd.isna(row.trait_pattern_label_count) else int(row.trait_pattern_label_count)
        )
        lines.append(
            f"| {row.view_label} | {row.heatmap_type} | `{row.png_path}` | "
            f"{row.displayed_cell_count} | {nominal_count} | "
            f"{row.bh_selected_cell_count} | {row.tested_cell_count} | {row.row_count} | "
            f"{condition_count} | {anatomy_count} | {branch_count} | "
            f"{metadata_label_count} | {pattern_label_count} |"
        )
    lines.append("")
    pages_path = output_dir / CONTRAST_PAGES_FILENAME
    if pages_path.exists():
        lines.extend(["### Signed Direction Pages", ""])
        for page in read_tsv(pages_path).itertuples(index=False):
            lines.append(
                f"- [{page.view_label}, page {page.page}]({Path(page.png_path).name}) "
                f"({page.row_count} rows; {page.condition_label_count} condition labels; "
                f"{page.anatomy_group_count} anatomy groups; "
                f"{page.displayed_cell_count} displayed cells; "
                f"{page.bh_selected_cell_count} BH-bordered cells)"
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
    page_records: list[dict[str, object]] = []

    for view in graph_manifest.itertuples(index=False):
        display_edges = read_tsv(resolve_repo_path(view.source_display_edges_path, repo_root))
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
                "bh_selected_cell_count": int(cells["selected_after_bh"].sum()),
                "tested_cell_count": tested_count,
                "min_positive_profiles": MIN_POSITIVE_PROFILES,
                "min_studies": MIN_STUDIES,
                "p_value_cutoff": P_VALUE_CUTOFF,
                "q_value_cutoff": Q_VALUE_CUTOFF,
                "color_value": "-log10(fisher_p_value)",
                "png_path": relative_path(output_path, repo_root),
                "source_display_edges_path": view.source_display_edges_path,
            }
        )

    for mode in CONTRAST_MODES:
        increased_view = views_by_key[mode.increased_view_key]
        decreased_view = views_by_key[mode.decreased_view_key]
        increased_path = full_edges_path(increased_view.source_display_edges_path, repo_root)
        decreased_path = full_edges_path(decreased_view.source_display_edges_path, repo_root)
        cells = prepare_contrast_cells(
            read_tsv(increased_path),
            read_tsv(decreased_path),
            mode,
            trait_order,
        )
        output_path, row_count, displayed_count, tested_count, extra = write_contrast_heatmap(
            cells,
            mode,
            trait_order,
            output_dir,
        )
        cells["heatmap_page"] = pd.Series(pd.NA, index=cells.index, dtype="Int64")
        if mode.key != "grouped":
            for page, labels in enumerate(
                contrast_row_pages(select_contrast_rows(cells)),
                start=1,
            ):
                in_page = cells["condition_body_site_pair"].isin(labels)
                page_cells = cells.loc[in_page]
                cells.loc[in_page, "heatmap_page"] = page
                page_records.append(
                    {
                        "view_key": f"{mode.key}_direction_contrast",
                        "view_label": mode.label,
                        "page": page,
                        "row_count": len(labels),
                        "condition_label_count": int(page_cells["condition"].nunique()),
                        "anatomy_group_count": int(page_cells["body_site_group"].nunique()),
                        "displayed_cell_count": int(page_cells["displayed_in_heatmap"].sum()),
                        "nominal_p_value_cell_count": int(
                            page_cells["nominal_p_value_selected"].sum()
                        ),
                        "bh_selected_cell_count": int(
                            (
                                page_cells["displayed_in_heatmap"] & page_cells["selected_after_bh"]
                            ).sum()
                        ),
                        "significant_cell_count": int(page_cells["nominal_p_value_selected"].sum()),
                        "absolute_effect_cutoff": SIGNED_EFFECT_ABS_CUTOFF,
                        "png_path": relative_path(contrast_page_path(output_path, page), repo_root),
                    }
                )
        contrast_tables.append(cells)
        preview_path = extra.get("preview_png_path", "")
        cluster_rows_path = extra.get("cluster_rows_path", "")
        cluster_summary_path = extra.get("cluster_summary_path", "")
        branch_tests_path = extra.get("branch_tests_path", "")
        pattern_tests_path = extra.get("pattern_tests_path", "")
        legend = extra["legend"]
        if not isinstance(legend, ContrastLegendSummary):
            raise TypeError("Signed-contrast output is missing its legend summary")
        manifest_records.append(
            {
                "view_key": f"{mode.key}_direction_contrast",
                "view_label": f"{mode.label}: signed direction contrast",
                "heatmap_type": "signed_direction_contrast",
                "focus_direction": "increased_vs_decreased",
                "condition_mode": mode.condition_mode,
                "row_count": row_count,
                "condition_label_count": legend.condition_label_count,
                "anatomy_group_count": legend.anatomy_group_count,
                "displayed_cell_count": displayed_count,
                "nominal_p_value_cell_count": legend.nominal_p_value_cell_count,
                "bh_selected_cell_count": legend.displayed_bh_selected_cell_count,
                "bh_selected_candidate_cell_count": int(cells["selected_after_bh"].sum()),
                "tested_cell_count": tested_count,
                "min_positive_profiles": MIN_POSITIVE_PROFILES,
                "min_studies": MIN_STUDIES,
                "p_value_cutoff": P_VALUE_CUTOFF,
                "q_value_cutoff": Q_VALUE_CUTOFF,
                "absolute_effect_cutoff": SIGNED_EFFECT_ABS_CUTOFF,
                "color_value": "increased_minus_decreased_positive_fraction",
                "display_rule": extra["display_rule"],
                "significant_cell_count": legend.nominal_p_value_cell_count,
                "filled_fraction": displayed_count / (row_count * len(trait_order))
                if row_count
                else 0,
                "page_count": extra["page_count"],
                "png_path": relative_path(output_path, repo_root),
                "preview_png_path": relative_path(preview_path, repo_root)
                if isinstance(preview_path, Path)
                else "",
                "cluster_rows_path": relative_path(cluster_rows_path, repo_root)
                if isinstance(cluster_rows_path, Path)
                else "",
                "cluster_summary_path": relative_path(cluster_summary_path, repo_root)
                if isinstance(cluster_summary_path, Path)
                else "",
                "branch_label_tests_path": relative_path(branch_tests_path, repo_root)
                if isinstance(branch_tests_path, Path)
                else "",
                "branch_pattern_tests_path": relative_path(pattern_tests_path, repo_root)
                if isinstance(pattern_tests_path, Path)
                else "",
                "clustering_method": extra["clustering_method"],
                "cluster_label_rule": extra["cluster_label_rule"],
                "significant_branch_count": extra["significant_branch_count"],
                "metadata_enrichment_branch_count": extra["metadata_enrichment_branch_count"],
                "trait_pattern_branch_count": extra["trait_pattern_branch_count"],
                "metadata_label_count": extra["metadata_label_count"],
                "trait_pattern_label_count": extra["trait_pattern_label_count"],
                "source_display_edges_path": "",
                "source_increased_edges_path": relative_path(increased_path, repo_root),
                "source_decreased_edges_path": relative_path(decreased_path, repo_root),
            }
        )

    pattern_table = pd.concat(pattern_tables, ignore_index=True)
    contrast_table = pd.concat(contrast_tables, ignore_index=True)
    manifest = pd.DataFrame.from_records(manifest_records)
    selected_cells = pd.concat(
        [
            pattern_table.assign(heatmap_type="trait_pattern"),
            contrast_table.assign(heatmap_type="signed_direction_contrast"),
        ],
        ignore_index=True,
        sort=False,
    )
    selected_cells = selected_cells[selected_cells["selected_after_bh"]].sort_values(
        ["view_key", "fisher_q_value", "fisher_p_value", "condition_body_site_pair", "trait_rank"],
        ascending=[True, True, True, True, True],
    )

    write_tsv(pattern_table, output_dir / PATTERN_TABLE_FILENAME)
    write_tsv(contrast_table, output_dir / CONTRAST_TABLE_FILENAME)
    write_tsv(selected_cells, output_dir / SELECTED_CELLS_FILENAME)
    write_tsv(manifest, output_dir / MANIFEST_FILENAME)
    write_tsv(
        pd.DataFrame.from_records(
            page_records,
            columns=[
                "view_key",
                "view_label",
                "page",
                "row_count",
                "condition_label_count",
                "anatomy_group_count",
                "displayed_cell_count",
                "nominal_p_value_cell_count",
                "bh_selected_cell_count",
                "significant_cell_count",
                "absolute_effect_cutoff",
                "png_path",
            ],
        ),
        output_dir / CONTRAST_PAGES_FILENAME,
    )
    append_readme_section(output_dir, manifest)
    return manifest


def main() -> None:
    """Write the Cytoscape-derived heatmap outputs."""
    manifest = write_six_shared_trait_cytoscape_heatmaps()
    print(f"Wrote Cytoscape heatmaps: {len(manifest)}")
    print(f"Wrote manifest: {OUTPUT_DIR / MANIFEST_FILENAME}")


if __name__ == "__main__":
    main()
