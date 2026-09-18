"""Create inferential heatmaps from the six-shared-trait Cytoscape inputs."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
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
CLUSTER_TREE_NODES_FILENAME = "six_shared_trait_signed_direction_cluster_tree_nodes.tsv"
CLUSTER_TREE_LEAVES_FILENAME = "six_shared_trait_signed_direction_cluster_tree_leaves.tsv"
CLUSTER_COMPARISON_FILENAME = "six_shared_trait_signed_direction_cluster_comparison.tsv"
CONTRAST_ROWS_PER_PAGE = 24
EXPORT_FLOAT_SIGNIFICANT_DIGITS = 12

MIN_POSITIVE_PROFILES = 3
MIN_STUDIES = 2
P_VALUE_CUTOFF = 0.01
Q_VALUE_CUTOFF = 0.05
MAX_MINUS_LOG10_P = 10.0
MAX_PATTERN_ROWS = 24
MIN_SPEARMAN_OVERLAP = 3
MIN_DIAGNOSTIC_SHARED_TRAITS = 5
UNDEFINED_SPEARMAN_DISTANCE = 1.0
CLUSTER_LINKAGE_METHOD = "average"
OUTLIER_IQR_MULTIPLIER = 1.5

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

OBSOLETE_CLUSTER_OUTPUTS = (
    "six_shared_trait_grouped_signed_direction_cluster_membership.tsv",
    "six_shared_trait_grouped_signed_direction_cluster_enrichment.tsv",
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
class SpearmanClustering:
    """Auditable row clustering for one signed direction matrix."""

    matrix: pd.DataFrame
    distances: pd.DataFrame
    overlaps: pd.DataFrame
    defined_correlations: pd.DataFrame
    clusterable_matrix: pd.DataFrame
    linkage_matrix: np.ndarray
    clustered_row_labels: list[str]
    insufficient_row_labels: list[str]
    tree_nodes: pd.DataFrame
    tree_leaves: pd.DataFrame
    pair_count: int
    undefined_pair_count: int
    undefined_pair_fraction: float
    distinct_distance_count: int
    distance_tie_excess_count: int
    distance_tie_fraction: float
    linkage_pair_count: int
    linkage_undefined_pair_count: int
    linkage_undefined_pair_fraction: float
    linkage_distinct_distance_count: int
    linkage_distance_tie_excess_count: int
    linkage_distance_tie_fraction: float


@dataclass(frozen=True)
class ContrastHeatmapOutputs:
    """Files, counts, and clustering emitted for one contrast mode."""

    overview_path: Path
    clustered_path: Path
    row_labels: list[str]
    displayed_cell_count: int
    tested_cell_count: int
    clustering: SpearmanClustering


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
    """Write TSVs with explicit missing values and stable float precision."""
    table.replace("", pd.NA).to_csv(
        path,
        sep="\t",
        index=False,
        na_rep="NA",
        float_format=f"%.{EXPORT_FLOAT_SIGNIFICANT_DIGITS}g",
    )


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
    output["passes_p_value_cutoff"] = finite_p & (p_values < P_VALUE_CUTOFF)
    output["displayed_in_heatmap"] = support & output["passes_p_value_cutoff"]
    output["significant_in_heatmap"] = output["displayed_in_heatmap"]
    output["selected_after_bh"] = support & finite_q & (q_values <= Q_VALUE_CUTOFF)

    output["blank_reason"] = ""
    output.loc[~support, "blank_reason"] = "low positive-z profile or study support"
    output.loc[support & ~finite_p, "blank_reason"] = "Fisher p=NA"
    output.loc[
        support & finite_p & ~output["passes_p_value_cutoff"],
        "blank_reason",
    ] = f"Fisher p>={P_VALUE_CUTOFF:g}"

    if "increased_minus_decreased_positive_fraction" in output:
        delta = pd.to_numeric(
            output["increased_minus_decreased_positive_fraction"],
            errors="coerce",
        ).fillna(0.0)
        direction = np.sign(delta)
        output["signed_minus_log10_fisher_p_value"] = (
            direction * output["minus_log10_fisher_p_value"]
        )
        output["plot_signed_minus_log10_fisher_p_value"] = (
            direction * output["plot_minus_log10_fisher_p_value"]
        )
        output["displayed_in_heatmap"] = support & finite_p
        output.loc[output["displayed_in_heatmap"], "blank_reason"] = ""

    return output


def clean_legacy_outputs(output_dir: Path = OUTPUT_DIR) -> None:
    """Remove outputs superseded by the current heatmap workflow."""
    for filename in (*LEGACY_CORRELATION_OUTPUTS, *OBSOLETE_CLUSTER_OUTPUTS):
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


def pairwise_spearman(
    left: np.ndarray,
    right: np.ndarray,
    *,
    min_overlap: int = MIN_SPEARMAN_OVERLAP,
) -> tuple[float, int, str]:
    """Return pairwise-complete Spearman rho and an audit status."""
    finite = np.isfinite(left) & np.isfinite(right)
    overlap = int(finite.sum())
    if overlap < min_overlap:
        return math.nan, overlap, f"fewer than {min_overlap} shared supported traits"

    left_ranks = pd.Series(left[finite]).rank(method="average").to_numpy(dtype=float)
    right_ranks = pd.Series(right[finite]).rank(method="average").to_numpy(dtype=float)
    if np.ptp(left_ranks) == 0.0 or np.ptp(right_ranks) == 0.0:
        return math.nan, overlap, "constant ranks across shared supported traits"

    rho = float(np.corrcoef(left_ranks, right_ranks)[0, 1])
    return float(np.clip(rho, -1.0, 1.0)), overlap, "comparable"


def spearman_distance_matrices(
    matrix: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Calculate pairwise-complete Spearman distance without filling blank cells."""
    labels = matrix.index.astype(str).tolist()
    values = matrix.to_numpy(dtype=float)
    row_count = len(labels)
    distances = np.zeros((row_count, row_count), dtype=float)
    overlaps = np.zeros((row_count, row_count), dtype=int)
    defined = np.zeros((row_count, row_count), dtype=bool)

    for left_index in range(row_count):
        overlaps[left_index, left_index] = int(np.isfinite(values[left_index]).sum())
        defined[left_index, left_index] = True
        for right_index in range(left_index + 1, row_count):
            rho, overlap, status = pairwise_spearman(
                values[left_index],
                values[right_index],
            )
            distance = 1.0 - rho if status == "comparable" else UNDEFINED_SPEARMAN_DISTANCE
            distances[left_index, right_index] = distance
            distances[right_index, left_index] = distance
            overlaps[left_index, right_index] = overlap
            overlaps[right_index, left_index] = overlap
            is_defined = status == "comparable"
            defined[left_index, right_index] = is_defined
            defined[right_index, left_index] = is_defined

    return (
        pd.DataFrame(distances, index=labels, columns=labels),
        pd.DataFrame(overlaps, index=labels, columns=labels),
        pd.DataFrame(defined, index=labels, columns=labels),
    )


def build_spearman_clustering(
    matrix: pd.DataFrame,
    cells: pd.DataFrame,
    *,
    view_key: str,
    view_label: str,
) -> SpearmanClustering:
    """Cluster informative rows after canonicalizing the linkage input order."""
    ranked_order = {
        str(label): order for order, label in enumerate(matrix.index.astype(str), start=1)
    }
    canonical_matrix = matrix.copy()
    canonical_matrix.index = canonical_matrix.index.astype(str)
    canonical_matrix = canonical_matrix.sort_index()
    distances, overlaps, defined = spearman_distance_matrices(canonical_matrix)

    row_count = len(canonical_matrix)
    peer_defined_counts = defined.sum(axis=1).astype(int) - 1
    clusterable_labels = peer_defined_counts[peer_defined_counts > 0].index.astype(str).tolist()
    insufficient_labels = peer_defined_counts[peer_defined_counts == 0].index.astype(str).tolist()
    clusterable_matrix = canonical_matrix.loc[clusterable_labels]
    clusterable_distances = distances.loc[clusterable_labels, clusterable_labels]
    clusterable_defined = defined.loc[clusterable_labels, clusterable_labels]
    clusterable_count = len(clusterable_matrix)

    if clusterable_count > 1:
        condensed = squareform(clusterable_distances.to_numpy(dtype=float), checks=False)
        linkage_matrix = linkage(
            condensed,
            method=CLUSTER_LINKAGE_METHOD,
            optimal_ordering=True,
        )
        leaf_ids = leaves_list(linkage_matrix).astype(int).tolist()
    else:
        linkage_matrix = np.empty((0, 4), dtype=float)
        leaf_ids = list(range(clusterable_count))

    clustered_labels = [str(clusterable_matrix.index[leaf_id]) for leaf_id in leaf_ids]
    final_labels = clustered_labels + insufficient_labels
    linkage_order = {label: order for order, label in enumerate(clustered_labels, start=1)}
    final_order = {label: order for order, label in enumerate(final_labels, start=1)}
    clusterable_leaf_id = {
        str(label): leaf_id for leaf_id, label in enumerate(clusterable_matrix.index.astype(str))
    }

    row_metadata = (
        cells.drop_duplicates("condition_body_site_pair")
        .set_index("condition_body_site_pair")[["condition", "body_site_group"]]
        .reindex(canonical_matrix.index)
    )
    leaf_records: list[dict[str, object]] = []
    for canonical_order, label in enumerate(canonical_matrix.index.astype(str), start=1):
        row_position = canonical_order - 1
        peer_defined = defined.iloc[row_position].to_numpy(dtype=bool).copy()
        peer_defined[row_position] = False
        if peer_defined.any():
            peer_distances = distances.iloc[row_position].to_numpy(dtype=float).copy()
            peer_distances[~peer_defined] = np.inf
            nearest_position = int(np.argmin(peer_distances))
            nearest_label: object = str(canonical_matrix.index[nearest_position])
            nearest_distance: object = float(peer_distances[nearest_position])
        else:
            nearest_label = pd.NA
            nearest_distance = math.nan
        included = label in clusterable_leaf_id
        leaf_records.append(
            {
                "view_key": view_key,
                "view_label": view_label,
                "leaf_id": clusterable_leaf_id.get(label, pd.NA),
                "ranked_heatmap_order": ranked_order[label],
                "canonical_input_order": canonical_order,
                "linkage_leaf_order": linkage_order.get(label, pd.NA),
                "clustered_heatmap_order": final_order[label],
                "condition_body_site_pair": label,
                "condition": row_metadata.loc[label, "condition"],
                "body_site_group": row_metadata.loc[label, "body_site_group"],
                "supported_trait_count": int(canonical_matrix.loc[label].notna().sum()),
                "defined_spearman_neighbor_count": int(peer_defined.sum()),
                "included_in_clustering": included,
                "low_information_for_clustering": not included,
                "clustering_exclusion_reason": ""
                if included
                else "no defined Spearman distance to any other row",
                "nearest_defined_neighbor": nearest_label,
                "nearest_defined_neighbor_distance": nearest_distance,
            }
        )

    node_records = [
        {
            "view_key": view_key,
            "view_label": view_label,
            "node_id": clusterable_count + merge_index,
            "left_child_id": int(merge[0]),
            "right_child_id": int(merge[1]),
            "merge_distance": float(merge[2]),
            "leaf_count": int(merge[3]),
        }
        for merge_index, merge in enumerate(linkage_matrix)
    ]
    tree_nodes = pd.DataFrame.from_records(
        node_records,
        columns=[
            "view_key",
            "view_label",
            "node_id",
            "left_child_id",
            "right_child_id",
            "merge_distance",
            "leaf_count",
        ],
    )
    tree_leaves = (
        pd.DataFrame.from_records(leaf_records)
        .sort_values("clustered_heatmap_order")
        .reset_index(drop=True)
    )

    pair_count = row_count * (row_count - 1) // 2
    if pair_count:
        row_triangle = np.triu_indices(row_count, k=1)
        undefined_pair_count = int((~defined.to_numpy(dtype=bool)[row_triangle]).sum())
        distance_values = distances.to_numpy(dtype=float)[row_triangle]
        distinct_distance_count = int(np.unique(distance_values).size)
    else:
        undefined_pair_count = 0
        distinct_distance_count = 0
    undefined_pair_fraction = undefined_pair_count / pair_count if pair_count else 0.0
    distance_tie_excess_count = max(0, pair_count - distinct_distance_count)
    distance_tie_fraction = distance_tie_excess_count / pair_count if pair_count else 0.0

    linkage_pair_count = clusterable_count * (clusterable_count - 1) // 2
    if linkage_pair_count:
        linkage_triangle = np.triu_indices(clusterable_count, k=1)
        linkage_defined_values = clusterable_defined.to_numpy(dtype=bool)[linkage_triangle]
        linkage_distance_values = clusterable_distances.to_numpy(dtype=float)[linkage_triangle]
        linkage_undefined_pair_count = int((~linkage_defined_values).sum())
        linkage_distinct_distance_count = int(np.unique(linkage_distance_values).size)
    else:
        linkage_undefined_pair_count = 0
        linkage_distinct_distance_count = 0
    linkage_undefined_pair_fraction = (
        linkage_undefined_pair_count / linkage_pair_count if linkage_pair_count else 0.0
    )
    linkage_distance_tie_excess_count = max(
        0,
        linkage_pair_count - linkage_distinct_distance_count,
    )
    linkage_distance_tie_fraction = (
        linkage_distance_tie_excess_count / linkage_pair_count if linkage_pair_count else 0.0
    )

    return SpearmanClustering(
        matrix=canonical_matrix,
        distances=distances,
        overlaps=overlaps,
        defined_correlations=defined,
        clusterable_matrix=clusterable_matrix,
        linkage_matrix=linkage_matrix,
        clustered_row_labels=clustered_labels,
        insufficient_row_labels=insufficient_labels,
        tree_nodes=tree_nodes,
        tree_leaves=tree_leaves,
        pair_count=pair_count,
        undefined_pair_count=undefined_pair_count,
        undefined_pair_fraction=undefined_pair_fraction,
        distinct_distance_count=distinct_distance_count,
        distance_tie_excess_count=distance_tie_excess_count,
        distance_tie_fraction=distance_tie_fraction,
        linkage_pair_count=linkage_pair_count,
        linkage_undefined_pair_count=linkage_undefined_pair_count,
        linkage_undefined_pair_fraction=linkage_undefined_pair_fraction,
        linkage_distinct_distance_count=linkage_distinct_distance_count,
        linkage_distance_tie_excess_count=linkage_distance_tie_excess_count,
        linkage_distance_tie_fraction=linkage_distance_tie_fraction,
    )


def split_joined_labels(value: object) -> list[str]:
    """Split semicolon-delimited rollup labels while preserving literal commas."""
    if pd.isna(value):
        return []
    return [part.strip() for part in str(value).split(";") if part.strip()]


def build_condition_rollup_map(grouped_edges: pd.DataFrame) -> pd.DataFrame:
    """Map raw conditions to the grouped label used at the same body site."""
    records: list[dict[str, str]] = []
    has_raw_conditions = "raw_conditions" in grouped_edges
    for row in grouped_edges.itertuples(index=False):
        raw_value = row.raw_conditions if has_raw_conditions else row.condition
        for raw_condition in split_joined_labels(raw_value):
            records.append(
                {
                    "raw_condition": raw_condition,
                    "body_site_group": str(row.body_site_group),
                    "grouped_condition": str(row.condition),
                    "grouped_condition_body_site_pair": pair_label(
                        row.condition,
                        row.body_site_group,
                    ),
                }
            )
    mapping = pd.DataFrame.from_records(records).drop_duplicates()
    if mapping.empty:
        return mapping
    ambiguity = mapping.groupby(["raw_condition", "body_site_group"])[
        "grouped_condition_body_site_pair"
    ].nunique()
    if (ambiguity > 1).any():
        examples = ", ".join(
            f"{condition} | {body_site}"
            for condition, body_site in ambiguity[ambiguity > 1].index[:5]
        )
        raise ValueError(f"Raw condition/body-site pairs map to multiple grouped rows: {examples}")
    return mapping.drop_duplicates(["raw_condition", "body_site_group"])


def compare_grouped_and_individual_clusters(
    grouped: SpearmanClustering,
    individual: SpearmanClustering,
    rollup_map: pd.DataFrame,
) -> pd.DataFrame:
    """Compare each raw-condition vector with its grouped aggregate and nearest peer."""
    mapping_by_key = {
        (row.raw_condition, row.body_site_group): (
            row.grouped_condition,
            row.grouped_condition_body_site_pair,
        )
        for row in rollup_map.itertuples(index=False)
    }
    grouped_order = dict(
        zip(
            grouped.tree_leaves["condition_body_site_pair"],
            grouped.tree_leaves["clustered_heatmap_order"],
            strict=True,
        )
    )
    individual_order = dict(
        zip(
            individual.tree_leaves["condition_body_site_pair"],
            individual.tree_leaves["clustered_heatmap_order"],
            strict=True,
        )
    )
    individual_metadata = individual.tree_leaves.set_index("condition_body_site_pair")
    records: list[dict[str, object]] = []
    for individual_pair in individual.matrix.index.astype(str):
        condition = str(individual_metadata.loc[individual_pair, "condition"])
        body_site = str(individual_metadata.loc[individual_pair, "body_site_group"])
        grouped_condition, grouped_pair = mapping_by_key.get(
            (condition, body_site),
            (pd.NA, pd.NA),
        )

        individual_position = individual.matrix.index.get_loc(individual_pair)
        peer_defined = (
            individual.defined_correlations.iloc[individual_position].to_numpy(dtype=bool).copy()
        )
        peer_defined[individual_position] = False
        if peer_defined.any():
            peer_distances = (
                individual.distances.iloc[individual_position].to_numpy(dtype=float).copy()
            )
            peer_distances[~peer_defined] = np.inf
            nearest_position = int(np.argmin(peer_distances))
            nearest_pair: object = str(individual.matrix.index[nearest_position])
            nearest_distance: object = float(peer_distances[nearest_position])
        else:
            nearest_pair = pd.NA
            nearest_distance = math.nan

        if pd.isna(grouped_pair):
            rho = math.nan
            shared_traits = 0
            status = "no grouped condition mapping"
        elif grouped_pair not in grouped.matrix.index:
            rho = math.nan
            shared_traits = 0
            status = "grouped row did not meet display support"
        else:
            rho, shared_traits, status = pairwise_spearman(
                individual.matrix.loc[individual_pair].to_numpy(dtype=float),
                grouped.matrix.loc[grouped_pair].to_numpy(dtype=float),
            )
        records.append(
            {
                "individual_condition": condition,
                "body_site_group": body_site,
                "individual_condition_body_site_pair": individual_pair,
                "grouped_condition": grouped_condition,
                "grouped_condition_body_site_pair": grouped_pair,
                "shared_supported_trait_count": shared_traits,
                "individual_vs_grouped_spearman_rho": rho,
                "individual_vs_grouped_spearman_distance": (
                    1.0 - rho if math.isfinite(rho) else math.nan
                ),
                "comparison_status": status,
                "individual_clustered_heatmap_order": individual_order[individual_pair],
                "grouped_clustered_heatmap_order": grouped_order.get(grouped_pair, pd.NA),
                "nearest_individual_condition_body_site_pair": nearest_pair,
                "nearest_individual_spearman_distance": nearest_distance,
            }
        )

    comparison = pd.DataFrame.from_records(records)
    grouped_by_individual_pair = dict(
        zip(
            comparison["individual_condition_body_site_pair"],
            comparison["grouped_condition_body_site_pair"],
            strict=True,
        )
    )
    comparison["nearest_individual_has_same_grouped_condition"] = [
        (
            grouped_by_individual_pair.get(nearest_pair) == grouped_pair
            if pd.notna(nearest_pair) and pd.notna(grouped_pair)
            else pd.NA
        )
        for nearest_pair, grouped_pair in zip(
            comparison["nearest_individual_condition_body_site_pair"],
            comparison["grouped_condition_body_site_pair"],
            strict=True,
        )
    ]
    comparison["grouped_contributing_individual_condition_count"] = (
        comparison.groupby(
            "grouped_condition_body_site_pair",
            dropna=False,
        )["individual_condition_body_site_pair"]
        .transform("nunique")
        .where(comparison["grouped_condition_body_site_pair"].notna(), pd.NA)
        .astype("Int64")
    )
    comparison["pattern_diagnostic_eligible"] = comparison["comparison_status"].eq(
        "comparable"
    ) & comparison["shared_supported_trait_count"].ge(MIN_DIAGNOSTIC_SHARED_TRAITS)
    comparison["_eligible_distance"] = comparison["individual_vs_grouped_spearman_distance"].where(
        comparison["pattern_diagnostic_eligible"]
    )
    comparison["global_difference_rank"] = (
        comparison["_eligible_distance"]
        .rank(
            method="min",
            ascending=False,
        )
        .astype("Int64")
    )
    comparison["within_group_difference_rank"] = (
        comparison.groupby(
            "grouped_condition_body_site_pair",
            dropna=False,
        )["_eligible_distance"]
        .rank(
            method="min",
            ascending=False,
        )
        .astype("Int64")
    )
    comparison["within_group_outlier_threshold"] = comparison.groupby(
        "grouped_condition_body_site_pair",
        dropna=False,
    )[
        "individual_vs_grouped_spearman_distance"
    ].transform(
        lambda values: (
            values.quantile(0.75)
            + OUTLIER_IQR_MULTIPLIER * (values.quantile(0.75) - values.quantile(0.25))
            if values.notna().sum() >= 4
            else math.nan
        )
    )
    comparison["within_group_iqr_outlier"] = (
        comparison["pattern_diagnostic_eligible"]
        & (
            comparison["individual_vs_grouped_spearman_distance"]
            > comparison["within_group_outlier_threshold"]
        )
    ).fillna(False)
    comparison["discordant_or_inverse_pattern"] = (
        comparison["pattern_diagnostic_eligible"]
        & (comparison["individual_vs_grouped_spearman_rho"] <= 0.0)
    ).fillna(False)
    return comparison.drop(columns="_eligible_distance").sort_values(
        ["global_difference_rank", "individual_condition_body_site_pair"],
        na_position="last",
    )


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
    """Pivot a Boolean cell flag into a row x trait matrix."""
    trait_ids = trait_order["trait_id"].tolist()
    flagged = cells[cells[flag_column]].copy()
    matrix = flagged.pivot_table(
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
    significance_borders: pd.DataFrame | None = None,
) -> None:
    """Render one row x six-trait heatmap with blank nonsignificant cells."""
    row_count = max(1, len(matrix.index))
    fig_height = max(6.2, 2.8 + 0.36 * row_count)
    fig, ax = plt.subplots(figsize=(13.8, fig_height))
    fig.subplots_adjust(left=0.34, right=0.90, top=0.82, bottom=0.18)
    if compact:
        fig.subplots_adjust(top=1 - 1.15 / fig_height, bottom=1.2 / fig_height)
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

    if significance_borders is not None:
        for y_index, row_label in enumerate(matrix.index):
            for x_index, trait_id in enumerate(trait_order["trait_id"]):
                if not significance_borders.loc[row_label, trait_id]:
                    continue
                ax.add_patch(
                    Rectangle(
                        (x_index - 0.5, y_index - 0.5),
                        1,
                        1,
                        fill=False,
                        edgecolor=SIGNIFICANCE_BORDER,
                        linewidth=1.8,
                        zorder=4,
                    )
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

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_clustered_contrast_heatmap(
    clustering: SpearmanClustering,
    significance_borders: pd.DataFrame,
    trait_order: pd.DataFrame,
    output_path: Path,
    *,
    title: str,
    subtitle: str,
) -> None:
    """Render a compact signed heatmap with its Spearman row dendrogram."""
    row_count = max(1, len(clustering.matrix))
    clusterable_count = len(clustering.clusterable_matrix)
    insufficient_count = len(clustering.insufficient_row_labels)
    fig_height = min(18.0, max(8.0, 4.0 + 0.13 * row_count))
    fig = plt.figure(figsize=(13.8, fig_height), facecolor="white")
    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.7, 6.0],
        left=0.04,
        right=0.69,
        top=0.84,
        bottom=0.18,
        wspace=0.02,
    )
    dendrogram_axis = fig.add_subplot(grid[0, 0])
    heatmap_axis = fig.add_subplot(grid[0, 1])

    if clusterable_count > 1:
        tree = dendrogram(
            clustering.linkage_matrix,
            orientation="left",
            no_labels=True,
            color_threshold=0,
            above_threshold_color=SLATE_500,
            link_color_func=lambda _node_id: SLATE_500,
            ax=dendrogram_axis,
        )
        row_positions = [int(value) for value in tree["leaves"]]
    else:
        row_positions = list(range(clusterable_count))
    dendrogram_axis.set_ylim(0, 10 * row_count)
    dendrogram_axis.set_xlabel("1 - Spearman rho", fontsize=7.5, color=SLATE_700)
    dendrogram_axis.tick_params(axis="x", labelsize=6.5, colors=SLATE_500)
    dendrogram_axis.tick_params(axis="y", left=False, labelleft=False)
    for spine in dendrogram_axis.spines.values():
        spine.set_visible(False)

    clustered_labels = [str(clustering.clusterable_matrix.index[index]) for index in row_positions]
    ordered_labels = clustered_labels + clustering.insufficient_row_labels
    ordered_matrix = clustering.matrix.loc[ordered_labels]
    ordered_borders = significance_borders.reindex(ordered_matrix.index)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(SLATE_100)
    values = ordered_matrix.to_numpy(dtype=float)
    if ordered_matrix.empty:
        values = np.full((1, len(trait_order)), np.nan)
    image = heatmap_axis.imshow(
        values,
        aspect="auto",
        cmap=cmap,
        vmin=-1.0,
        vmax=1.0,
        origin="lower",
        extent=(-0.5, len(trait_order) - 0.5, 0, 10 * row_count),
    )
    heatmap_axis.set_xticks(
        np.arange(len(trait_order)),
        [wrap_label(label, width=16) for label in trait_order["trait_label"]],
    )
    heatmap_axis.set_yticks(
        5 + 10 * np.arange(len(ordered_matrix)),
        [wrap_label(label, width=42) for label in ordered_matrix.index],
    )
    heatmap_axis.yaxis.tick_right()
    heatmap_axis.yaxis.set_label_position("right")
    heatmap_axis.tick_params(axis="x", labelrotation=35, labelsize=8.0)
    heatmap_axis.tick_params(axis="y", labelsize=4.8, pad=2)
    heatmap_axis.set_xlabel("Six shared traits", color=SLATE_700)
    heatmap_axis.set_ylabel("Condition | body-site pair", color=SLATE_700, labelpad=8)
    heatmap_axis.set_xticks(np.arange(-0.5, len(trait_order), 1), minor=True)
    heatmap_axis.set_yticks(np.arange(0, 10 * row_count + 1, 10), minor=True)
    heatmap_axis.grid(which="minor", color="white", linewidth=0.7)
    heatmap_axis.tick_params(which="minor", bottom=False, left=False)
    for spine in heatmap_axis.spines.values():
        spine.set_visible(False)

    if insufficient_count:
        boundary = 10 * clusterable_count
        dendrogram_axis.axhspan(boundary, 10 * row_count, color=SLATE_100, zorder=-1)
        dendrogram_axis.axhline(boundary, color=SLATE_500, linewidth=0.9)
        dendrogram_axis.text(
            0.5,
            boundary + 5 * insufficient_count,
            "Insufficient overlap\n(not clustered)",
            transform=dendrogram_axis.get_yaxis_transform(),
            ha="center",
            va="center",
            fontsize=7.2,
            color=SLATE_700,
        )
        heatmap_axis.axhline(boundary, color=SLATE_500, linewidth=1.0)

    for y_index, row_label in enumerate(ordered_matrix.index):
        for x_index, trait_id in enumerate(trait_order["trait_id"]):
            if not ordered_borders.loc[row_label, trait_id]:
                continue
            heatmap_axis.add_patch(
                Rectangle(
                    (x_index - 0.5, y_index * 10),
                    1,
                    10,
                    fill=False,
                    edgecolor=SIGNIFICANCE_BORDER,
                    linewidth=1.1,
                    zorder=4,
                )
            )

    fig.text(
        0.04,
        0.975,
        title,
        ha="left",
        va="top",
        fontsize=15.5,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.04,
        0.94,
        "\n".join(wrap(subtitle, width=130)),
        ha="left",
        va="top",
        fontsize=9.0,
        color=SLATE_500,
    )
    colorbar_axis = fig.add_axes((0.30, 0.065, 0.30, 0.014))
    colorbar = fig.colorbar(image, cax=colorbar_axis, orientation="horizontal")
    colorbar.set_label("Increased - decreased positive-z fraction", color=SLATE_700)
    colorbar.ax.tick_params(labelsize=7)
    colorbar.outline.set_edgecolor(SLATE_200)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=240, bbox_inches="tight", facecolor="white")
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
) -> ContrastHeatmapOutputs:
    """Write one approach-4 signed direction contrast heatmap."""
    row_labels = select_contrast_rows(cells)
    matrix = cell_matrix(
        cells,
        row_labels,
        trait_order,
        value_column="increased_minus_decreased_positive_fraction",
    )
    annotations = annotation_matrix(cells, row_labels, trait_order, mode="contrast")
    bh_borders = flag_matrix(
        cells,
        row_labels,
        trait_order,
        flag_column="selected_after_bh",
    )
    output_path = (
        output_dir
        / f"six_shared_traits_{mode.key}_condition_body_site_signed_direction_contrast_heatmap.png"
    )
    clustered_output_path = output_path.with_name(
        output_path.name.replace("_contrast_heatmap.png", "_contrast_clustered_heatmap.png")
    )
    displayed_cells = int(cells["displayed_in_heatmap"].sum())
    raw_p_value_cell_count = int(cells["significant_in_heatmap"].sum())
    bh_selected_cell_count = int(cells["selected_after_bh"].sum())
    subtitle = (
        "Color and numbers show increased minus decreased positive-z fractions. "
        "Red favors increased profiles; blue favors decreased profiles. "
        f"A dark border marks BH q<={Q_VALUE_CUTOFF:g} ({bh_selected_cell_count} cells); "
        f"{raw_p_value_cell_count} cells met uncorrected Fisher p<{P_VALUE_CUTOFF:g}. "
        "All testable cells with >=3 positive-z profiles in either direction and >=2 total "
        "studies are shown. Gray = unavailable or insufficient support, not nonsignificance. "
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
        significance_borders=bh_borders,
    )
    pages = contrast_row_pages(row_labels)
    for page, labels in enumerate(pages, start=1):
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
            significance_borders=bh_borders.loc[labels],
        )
    clustering = build_spearman_clustering(
        matrix,
        cells,
        view_key=f"{mode.key}_direction_contrast",
        view_label=mode.label,
    )
    clustered_subtitle = (
        "Average linkage uses pairwise-complete Spearman distance after ascending label sort "
        "for deterministic tie breaking. "
        f"Before exclusion, {clustering.undefined_pair_count}/{clustering.pair_count} "
        f"({clustering.undefined_pair_fraction:.1%}) row-pair correlations were undefined; "
        f"{len(clustering.insufficient_row_labels)} rows with no defined neighbor are shown "
        "separately and not linked. Among linked rows, "
        f"{clustering.linkage_undefined_pair_count}/{clustering.linkage_pair_count} "
        f"({clustering.linkage_undefined_pair_fraction:.1%}) distances use neutral "
        f"{UNDEFINED_SPEARMAN_DISTANCE:g}; exact-distance tie excess is "
        f"{clustering.linkage_distance_tie_fraction:.1%}. No dendrogram cut is interpreted. "
        f"Gray cells remain missing; dark borders mark BH q<={Q_VALUE_CUTOFF:g} "
        f"({bh_selected_cell_count} cells), versus {raw_p_value_cell_count} cells at "
        f"uncorrected Fisher p<{P_VALUE_CUTOFF:g}."
    )
    plot_clustered_contrast_heatmap(
        clustering,
        bh_borders,
        trait_order,
        clustered_output_path,
        title=f"Clustered signed direction contrast: {mode.label}",
        subtitle=clustered_subtitle,
    )
    return ContrastHeatmapOutputs(
        overview_path=output_path,
        clustered_path=clustered_output_path,
        row_labels=row_labels,
        displayed_cell_count=displayed_cells,
        tested_cell_count=len(cells),
        clustering=clustering,
    )


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
    prefix, prior_section = content.split(split_marker, maxsplit=1)
    next_section_start = prior_section.find("\n## ")
    suffix = ""
    if next_section_start >= 0:
        suffix = prior_section[next_section_start + 1 :].strip()
    lines = [
        new_marker,
        "",
        f"- Manifest: `{MANIFEST_FILENAME}`",
        f"- Pattern matrix table: `{PATTERN_TABLE_FILENAME}`",
        f"- Signed contrast table: `{CONTRAST_TABLE_FILENAME}`",
        f"- BH-selected cell table: `{SELECTED_CELLS_FILENAME}`",
        f"- Cluster-tree nodes: `{CLUSTER_TREE_NODES_FILENAME}`",
        f"- Cluster-tree leaves and row order: `{CLUSTER_TREE_LEAVES_FILENAME}`",
        f"- Individual-versus-grouped comparison: `{CLUSTER_COMPARISON_FILENAME}`",
        "",
        "Trait-pattern heatmaps use the six-shared-trait display-edge TSVs that feed "
        "the Cytoscape graph views. Signed contrasts use the full increased and decreased "
        "edge TSVs before graph display selection. They preserve condition/body-site/trait records "
        "instead of correlating already-collapsed condition and body-site nodes.",
        "",
        "Signed contrasts display all supported, testable effects, including nonsignificant "
        "results. Gray cells are unavailable or under-supported. Trait-pattern heatmaps "
        "continue to mask cells that fail the raw p-value threshold.",
        "",
        "### How Conditions and Blank Cells Are Chosen",
        "",
        "The signed heatmaps do not choose a top-N list of diseases. They read the full "
        "increased and decreased Sankey edge tables before graph display selection, discard "
        "unknown raw condition labels upstream, and include every condition/body-site row "
        "with at least one supported, finite cell. The grouped view applies the v3 clinical "
        "condition rollup rules; the individual view keeps the original curated condition "
        "label. Both views use the same grouped v3 body-site labels.",
        "",
        "The signed contrast table contains candidate edge records, not a complete Cartesian "
        "grid. Its exported `blank_reason` is `low positive-z profile or study support` for "
        "a present candidate that fails display support and is empty (`NA` in the TSV) for a "
        "displayed candidate. A condition/body-site/trait combination for which neither "
        "direction emitted a candidate edge is absent from the table rather than assigned a "
        "blank reason. In the image, both an absent record and a present under-supported "
        "candidate appear gray. A direction can have zero positive-z profiles while the cell "
        "remains colored when total profiles exist and the other direction meets support. "
        "Finite tests that miss either statistical cutoff are still colored and numbered. A "
        "white cell containing `+0.00` is a measured null contrast, not missing data. In "
        "trait-pattern heatmaps only, gray can also mean that the raw one-sided Fisher test "
        f"did not pass p<{P_VALUE_CUTOFF:g}.",
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
        f"A dark border marks BH q<={Q_VALUE_CUTOFF:g}; the uncorrected Fisher "
        f"p<{P_VALUE_CUTOFF:g} count is reported separately below. Unbordered cells are "
        "descriptive effects. Cells require a finite test, "
        "at least 3 positive-z profiles in "
        "either direction, and at least 2 total studies. Every condition/body-site row "
        "with a supported cell is included, without a row cap. Rows are grouped by body "
        "site and ranked by significant-trait count, coverage, and evidence within each site.",
        "- `signed_direction_contrast_clustered`: A separate compact overview clusters rows "
        "with average linkage on pairwise-complete Spearman distance (`1 - rho`). A pairwise "
        f"correlation requires at least {MIN_SPEARMAN_OVERLAP} shared supported traits. Blank "
        "cells are never filled with zero. Rows are sorted by exact condition/body-site pair "
        "label before distance construction so SciPy's tied-distance choices are reproducible. "
        "Rows with no defined Spearman distance to any other row are excluded from linkage and "
        "shown in a separate `Insufficient overlap (not clustered)` block. Remaining undefined "
        f"row-pair distances use neutral distance {UNDEFINED_SPEARMAN_DISTANCE:g}. Figure "
        "subtitles and the manifest report undefined-pair and exact-distance tie fractions. "
        "No dendrogram cut or cluster IDs are interpreted. Trait columns remain in the fixed "
        "six-trait order.",
        f"- Readable pages contain at most {CONTRAST_ROWS_PER_PAGE} rows each. "
        f"`{CONTRAST_PAGES_FILENAME}` records their paths and counts; `heatmap_page` in "
        "the contrast table locates each displayed row. The original PNGs show full overviews.",
        f"- BH q-values are computed within each heatmap view. "
        f"`selected_after_bh=True` marks supported cells with q<={Q_VALUE_CUTOFF:g}; "
        "those cells receive the dark border and are also collected in the BH-selected cell "
        "table. `significant_in_heatmap` retains the uncorrected raw-p flag for audit only.",
        "- The Fisher tables count profiles, not independent studies. A BH-bordered cell is "
        "therefore a profile-level diagnostic and does not establish study-level replication.",
        "",
        "| heatmap | type | full/ranked PNG | clustered PNG | displayed cells | "
        "raw p<0.01 cells | BH-selected cells | tested cells | rows |",
        "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in manifest.itertuples(index=False):
        clustered_path = getattr(row, "clustered_png_path", pd.NA)
        clustered_text = "" if pd.isna(clustered_path) else f"`{clustered_path}`"
        lines.append(
            f"| {row.view_label} | {row.heatmap_type} | `{row.png_path}` | "
            f"{clustered_text} | "
            f"{row.displayed_cell_count} | {int(row.significant_cell_count)} | "
            f"{row.bh_selected_cell_count} | "
            f"{row.tested_cell_count} | {row.row_count} |"
        )
    clustered_rows = manifest[manifest["heatmap_type"] == "signed_direction_contrast"]
    if not clustered_rows.empty:
        lines.extend(
            [
                "",
                "### Clustered Direction Overviews",
                "",
                "The grouped view is previewed below. The 125-row individual view is linked "
                "directly because a page-width thumbnail is illegible. Use the ranked pages "
                "below when numeric annotations need to be read.",
                "",
            ]
        )
        for row in clustered_rows.itertuples(index=False):
            image_name = Path(row.clustered_png_path).name
            lines.extend([f"#### {row.view_label}", ""])
            if row.view_key == "grouped_direction_contrast":
                lines.extend(
                    [
                        f'<a href="{image_name}"><img src="{image_name}" '
                        f'alt="{row.view_label}" width="760"></a>',
                        "",
                    ]
                )
            else:
                lines.extend(
                    [
                        f"[Open the full-resolution clustered heatmap]({image_name}).",
                        "",
                    ]
                )
        lines.extend(
            [
                "### Clustering Diagnostics",
                "",
                "The undefined fraction covers every candidate row pair before low-information "
                "rows are excluded. Tie excess is `(linkage pairs - distinct distances) / "
                "linkage pairs`; it is zero only when every linkage distance is unique.",
                "",
                "| view | rows linked | insufficient-overlap rows | undefined candidate pairs | "
                "candidate tie excess | undefined linkage pairs | linkage tie excess |",
                "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
            ]
        )
        for row in clustered_rows.itertuples(index=False):
            lines.append(
                f"| {row.view_label} | {int(row.clustered_row_count)} | "
                f"{int(row.insufficient_overlap_row_count)} | "
                f"{int(row.undefined_spearman_pair_count)}/{int(row.spearman_pair_count)} "
                f"({row.undefined_spearman_pair_fraction:.1%}) | "
                f"{int(row.distance_tie_excess_count)}/{int(row.spearman_pair_count)} "
                f"({row.distance_tie_fraction:.1%}) | "
                f"{int(row.linkage_undefined_spearman_pair_count)}/"
                f"{int(row.linkage_pair_count)} "
                f"({row.linkage_undefined_spearman_pair_fraction:.1%}) | "
                f"{int(row.linkage_distance_tie_excess_count)}/{int(row.linkage_pair_count)} "
                f"({row.linkage_distance_tie_fraction:.1%}) |"
            )
    lines.append("")
    pages_path = output_dir / CONTRAST_PAGES_FILENAME
    if pages_path.exists():
        lines.extend(["### Signed Direction Pages", ""])
        for page in read_tsv(pages_path).itertuples(index=False):
            lines.append(
                f"- [{page.view_label}, page {page.page}]({Path(page.png_path).name}) "
                f"({page.row_count} rows, {page.displayed_cell_count} filled cells)"
            )
        lines.append("")
    comparison_path = output_dir / CLUSTER_COMPARISON_FILENAME
    if comparison_path.exists():
        comparison = read_tsv(comparison_path)
        comparable = comparison[comparison["pattern_diagnostic_eligible"].astype(bool)].head(10)
        if not comparable.empty:
            lines.extend(
                [
                    "### Largest Individual-versus-Grouped Pattern Differences",
                    "",
                    f"These are the ten largest finite `1 - rho` distances among rows with at "
                    f"least {MIN_DIAGNOSTIC_SHARED_TRAITS} shared supported traits, not a "
                    "formal biological outlier test. Lower-overlap rows remain in the complete "
                    "table but are ineligible for ranks and diagnostic flags. Each individual "
                    "row is compared with a grouped aggregate that includes that row, so the "
                    "comparison is not independent; a one-contributor aggregate is the same "
                    "source row at the broader label. The complete table records contributor "
                    "counts, within-group IQR flags, and each row's nearest Spearman neighbor.",
                    "",
                    "| individual condition | grouped condition | body site | shared traits | "
                    "group contributors | rho | distance | discordant | IQR outlier |",
                    "| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |",
                ]
            )
            for row in comparable.itertuples(index=False):
                lines.append(
                    f"| {row.individual_condition} | {row.grouped_condition} | "
                    f"{row.body_site_group} | {row.shared_supported_trait_count} | "
                    f"{row.grouped_contributing_individual_condition_count} | "
                    f"{row.individual_vs_grouped_spearman_rho:.3f} | "
                    f"{row.individual_vs_grouped_spearman_distance:.3f} | "
                    f"{row.discordant_or_inverse_pattern} | "
                    f"{row.within_group_iqr_outlier} |"
                )
            lines.append("")

    updated = prefix.rstrip() + "\n\n" + "\n".join(lines).rstrip() + "\n"
    if suffix:
        updated += "\n" + suffix + "\n"
    readme_path.write_text(updated, encoding="utf-8")


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
    cluster_tree_nodes: list[pd.DataFrame] = []
    cluster_tree_leaves: list[pd.DataFrame] = []
    clusterings: dict[str, SpearmanClustering] = {}
    contrast_source_edges: dict[str, pd.DataFrame] = {}

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
                "significant_cell_count": int(cells["significant_in_heatmap"].sum()),
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
        increased_edges = read_tsv(increased_path)
        decreased_edges = read_tsv(decreased_path)
        contrast_source_edges[mode.key] = pd.concat(
            [increased_edges, decreased_edges],
            ignore_index=True,
        )
        cells = prepare_contrast_cells(
            increased_edges,
            decreased_edges,
            mode,
            trait_order,
        )
        outputs = write_contrast_heatmap(
            cells,
            mode,
            trait_order,
            output_dir,
        )
        clusterings[mode.key] = outputs.clustering
        cluster_tree_nodes.append(outputs.clustering.tree_nodes)
        cluster_tree_leaves.append(outputs.clustering.tree_leaves)
        clustered_order = dict(
            zip(
                outputs.clustering.tree_leaves["condition_body_site_pair"],
                outputs.clustering.tree_leaves["clustered_heatmap_order"],
                strict=True,
            )
        )
        cells["clustered_heatmap_order"] = (
            cells["condition_body_site_pair"].map(clustered_order).astype("Int64")
        )
        cells["heatmap_page"] = pd.Series(pd.NA, index=cells.index, dtype="Int64")
        for page, labels in enumerate(contrast_row_pages(outputs.row_labels), start=1):
            in_page = cells["condition_body_site_pair"].isin(labels)
            cells.loc[in_page, "heatmap_page"] = page
            page_records.append(
                {
                    "view_key": f"{mode.key}_direction_contrast",
                    "view_label": mode.label,
                    "page": page,
                    "row_count": len(labels),
                    "displayed_cell_count": int(cells.loc[in_page, "displayed_in_heatmap"].sum()),
                    "significant_cell_count": int(
                        cells.loc[in_page, "significant_in_heatmap"].sum()
                    ),
                    "png_path": relative_path(
                        contrast_page_path(outputs.overview_path, page),
                        repo_root,
                    ),
                }
            )
        contrast_tables.append(cells)
        manifest_records.append(
            {
                "view_key": f"{mode.key}_direction_contrast",
                "view_label": f"{mode.label}: signed direction contrast",
                "heatmap_type": "signed_direction_contrast",
                "focus_direction": "increased_vs_decreased",
                "condition_mode": mode.condition_mode,
                "row_count": len(outputs.row_labels),
                "displayed_cell_count": outputs.displayed_cell_count,
                "bh_selected_cell_count": int(cells["selected_after_bh"].sum()),
                "tested_cell_count": outputs.tested_cell_count,
                "min_positive_profiles": MIN_POSITIVE_PROFILES,
                "min_studies": MIN_STUDIES,
                "p_value_cutoff": P_VALUE_CUTOFF,
                "q_value_cutoff": Q_VALUE_CUTOFF,
                "color_value": "increased_minus_decreased_positive_fraction",
                "display_rule": (
                    "all supported testable cells; BH q cutoff marks dark borders; "
                    "raw p cutoff is audit-only"
                ),
                "significant_cell_count": int(cells["significant_in_heatmap"].sum()),
                "filled_fraction": outputs.displayed_cell_count
                / (len(outputs.row_labels) * len(trait_order))
                if outputs.row_labels
                else 0,
                "page_count": len(contrast_row_pages(outputs.row_labels)),
                "png_path": relative_path(outputs.overview_path, repo_root),
                "clustered_png_path": relative_path(outputs.clustered_path, repo_root),
                "cluster_metric": "pairwise-complete Spearman distance (1-rho)",
                "cluster_linkage": CLUSTER_LINKAGE_METHOD,
                "cluster_min_shared_traits": MIN_SPEARMAN_OVERLAP,
                "cluster_input_order": (
                    "condition_body_site_pair ascending before distance construction"
                ),
                "undefined_cluster_distance": UNDEFINED_SPEARMAN_DISTANCE,
                "clustered_row_count": len(outputs.clustering.clustered_row_labels),
                "insufficient_overlap_row_count": len(outputs.clustering.insufficient_row_labels),
                "spearman_pair_count": outputs.clustering.pair_count,
                "undefined_spearman_pair_count": outputs.clustering.undefined_pair_count,
                "undefined_spearman_pair_fraction": (outputs.clustering.undefined_pair_fraction),
                "distinct_spearman_distance_count": (outputs.clustering.distinct_distance_count),
                "distance_tie_excess_count": outputs.clustering.distance_tie_excess_count,
                "distance_tie_fraction": outputs.clustering.distance_tie_fraction,
                "linkage_pair_count": outputs.clustering.linkage_pair_count,
                "linkage_undefined_spearman_pair_count": (
                    outputs.clustering.linkage_undefined_pair_count
                ),
                "linkage_undefined_spearman_pair_fraction": (
                    outputs.clustering.linkage_undefined_pair_fraction
                ),
                "linkage_distinct_distance_count": (
                    outputs.clustering.linkage_distinct_distance_count
                ),
                "linkage_distance_tie_excess_count": (
                    outputs.clustering.linkage_distance_tie_excess_count
                ),
                "linkage_distance_tie_fraction": (outputs.clustering.linkage_distance_tie_fraction),
                "dendrogram_cut_applied": False,
                "cluster_tree_nodes_path": relative_path(
                    output_dir / CLUSTER_TREE_NODES_FILENAME,
                    repo_root,
                ),
                "cluster_tree_leaves_path": relative_path(
                    output_dir / CLUSTER_TREE_LEAVES_FILENAME,
                    repo_root,
                ),
                "cluster_comparison_path": relative_path(
                    output_dir / CLUSTER_COMPARISON_FILENAME,
                    repo_root,
                ),
                "source_display_edges_path": "",
                "source_increased_edges_path": relative_path(increased_path, repo_root),
                "source_decreased_edges_path": relative_path(decreased_path, repo_root),
            }
        )

    pattern_table = pd.concat(pattern_tables, ignore_index=True)
    contrast_table = pd.concat(contrast_tables, ignore_index=True)
    tree_nodes = pd.concat(cluster_tree_nodes, ignore_index=True)
    tree_leaves = pd.concat(cluster_tree_leaves, ignore_index=True)
    rollup_map = build_condition_rollup_map(contrast_source_edges["grouped"])
    cluster_comparison = compare_grouped_and_individual_clusters(
        clusterings["grouped"],
        clusterings["individual_conditions"],
        rollup_map,
    )
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
    write_tsv(tree_nodes, output_dir / CLUSTER_TREE_NODES_FILENAME)
    write_tsv(tree_leaves, output_dir / CLUSTER_TREE_LEAVES_FILENAME)
    write_tsv(cluster_comparison, output_dir / CLUSTER_COMPARISON_FILENAME)
    write_tsv(
        pd.DataFrame.from_records(
            page_records,
            columns=[
                "view_key",
                "view_label",
                "page",
                "row_count",
                "displayed_cell_count",
                "significant_cell_count",
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
