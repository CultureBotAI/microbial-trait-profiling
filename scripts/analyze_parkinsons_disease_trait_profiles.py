from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
TRAIT_VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = TRAIT_VIS_DIR / "trait-profile-investigation"
OUTPUT_DIR = INVESTIGATION_DIR / "parkinsons-disease"

MATRIX_PATH = TRAIT_VIS_DIR / "study_trait_profile_clustered_heatmap_matrix.tsv"
METADATA_PATH = TRAIT_VIS_DIR / "study_trait_profile_study_metadata.tsv"
TRAIT_METADATA_PATH = TRAIT_VIS_DIR / "study_trait_profile_trait_metadata.tsv"

PD_CONDITION = "Parkinson's disease"
BODY_SITE_GROUP = "Feces"

DIRECTION_LABELS = {"inc": "Increased", "dec": "Decreased"}
DIRECTION_COLORS = {"inc": "#c2410c", "dec": "#2563eb"}


@dataclass(frozen=True)
class AnalysisConfig:
    key: str
    label: str
    short_label: str
    score_column: str
    positive_count_column: str
    trait_selection_path: Path
    candidate_taxa_path: Path
    focus_direction: str
    focus_label: str
    comparator_direction: str
    comparator_label: str
    bias_column: str


AEROBE_CONFIG = AnalysisConfig(
    key="aerobe_enrichment_v2",
    label="Aerobe enrichment v2",
    short_label="Aerobe enrichment",
    score_column="aerobe_v2_score",
    positive_count_column="aerobe_v2_positive_trait_count",
    trait_selection_path=INVESTIGATION_DIR / "aerobe-enrichment-v2" / "guild_trait_selection.tsv",
    candidate_taxa_path=INVESTIGATION_DIR
    / "aerobe-enrichment-v2"
    / "candidate_taxa_condition_direction.tsv",
    focus_direction="inc",
    focus_label="increased",
    comparator_direction="dec",
    comparator_label="decreased",
    bias_column="increased_profile_tendency",
)

ANAEROBE_CONFIG = AnalysisConfig(
    key="anaerobe_depletion_v2",
    label="Anaerobe depletion v2",
    short_label="Anaerobe depletion",
    score_column="anaerobe_v2_score",
    positive_count_column="anaerobe_v2_positive_trait_count",
    trait_selection_path=INVESTIGATION_DIR / "anaerobe-depletion-v2" / "guild_trait_selection.tsv",
    candidate_taxa_path=INVESTIGATION_DIR
    / "anaerobe-depletion-v2"
    / "candidate_taxa_condition_direction.tsv",
    focus_direction="dec",
    focus_label="decreased",
    comparator_direction="inc",
    comparator_label="increased",
    bias_column="depletion_bias",
)

ANALYSIS_CONFIGS = (AEROBE_CONFIG, ANAEROBE_CONFIG)


def clean_text(value: object) -> str:
    if pd.isna(value):
        return "unknown"
    text = str(value).strip()
    return text if text else "unknown"


def read_tsv(path: Path, *, index_col: int | str | None = None) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=index_col)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False)


def mean_or_nan(values: pd.Series) -> float:
    if values.empty:
        return math.nan
    return float(values.mean())


def median_or_nan(values: pd.Series) -> float:
    if values.empty:
        return math.nan
    return float(values.median())


def join_unique(values: pd.Series) -> str:
    seen: list[str] = []
    for value in values:
        text = clean_text(value)
        if text != "unknown" and text not in seen:
            seen.append(text)
    return "; ".join(seen)


def fmt_float(value: object, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}f}"


def fmt_p_value(value: object) -> str:
    if pd.isna(value):
        return "NA"
    p_value = float(value)
    if p_value < 0.0001:
        return f"{p_value:.2e}"
    return f"{p_value:.4f}"


def load_study_heatmap() -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix = read_tsv(MATRIX_PATH, index_col=0)
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    metadata = read_tsv(METADATA_PATH, index_col=0)
    trait_metadata = read_tsv(TRAIT_METADATA_PATH, index_col=0)

    if len(matrix) != len(metadata):
        msg = f"Matrix/profile metadata row mismatch: {len(matrix)} vs {len(metadata)}"
        raise ValueError(msg)

    matrix.index = metadata.index
    column_labels = {
        row["display_label"]: row["object_name"] for _, row in trait_metadata.iterrows()
    }
    matrix = matrix.rename(columns=column_labels)
    return matrix, metadata


def load_guild_traits(config: AnalysisConfig) -> list[str]:
    selection = read_tsv(config.trait_selection_path)
    trait_column = (
        "profile_matrix_object_name"
        if "profile_matrix_object_name" in selection.columns
        else "trait_name"
    )
    traits = [clean_text(value) for value in selection[trait_column]]
    return [trait for trait in traits if trait != "unknown"]


def add_trait_scores(matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    scored = metadata.copy()
    scored["profile_id"] = scored.index

    for config in ANALYSIS_CONFIGS:
        guild_traits = load_guild_traits(config)
        available_traits = [trait for trait in guild_traits if trait in matrix.columns]
        if not available_traits:
            msg = f"No selected traits for {config.label} were found in the profile matrix."
            raise ValueError(msg)

        selected = matrix[available_traits]
        scored[config.score_column] = selected.mean(axis=1)
        scored[config.positive_count_column] = selected.gt(0.0).sum(axis=1)
        scored[f"{config.key}_selected_trait_count"] = len(available_traits)

    for column in ("condition", "body_site", "body_site_group", "condition_group"):
        scored[column] = scored[column].map(clean_text)

    scored["direction_label"] = scored["inc_dec"].map(DIRECTION_LABELS).fillna(scored["inc_dec"])
    return scored


def parkinsons_feces_scores(scored: pd.DataFrame) -> pd.DataFrame:
    source = scored[
        (scored["condition"] == PD_CONDITION) & (scored["body_site_group"] == BODY_SITE_GROUP)
    ].copy()
    if source.empty:
        msg = "No Parkinson's disease feces profiles were found."
        raise ValueError(msg)

    columns = [
        "profile_id",
        "study_id",
        "experiment_id",
        "inc_dec",
        "direction_label",
        "condition",
        "body_site",
        "body_site_group",
        "condition_group",
        "publication_year",
        "title",
        "signature_count",
        "taxon_count",
        "aerobe_v2_score",
        "aerobe_v2_positive_trait_count",
        "aerobe_enrichment_v2_selected_trait_count",
        "anaerobe_v2_score",
        "anaerobe_v2_positive_trait_count",
        "anaerobe_depletion_v2_selected_trait_count",
    ]
    return source[columns].sort_values(["study_id", "experiment_id", "inc_dec", "profile_id"])


def summarize_direction_scores(scores: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for direction, direction_label in DIRECTION_LABELS.items():
        group = scores[scores["inc_dec"] == direction]
        records.append(
            {
                "condition": PD_CONDITION,
                "body_site_group": BODY_SITE_GROUP,
                "direction": direction,
                "direction_label": direction_label,
                "profile_count": len(group),
                "study_count": group["study_id"].nunique(),
                "experiment_count": group["experiment_id"].nunique(),
                "aerobe_v2_mean_score": mean_or_nan(group["aerobe_v2_score"]),
                "aerobe_v2_median_score": median_or_nan(group["aerobe_v2_score"]),
                "aerobe_v2_positive_profile_count": int(
                    group["aerobe_v2_positive_trait_count"].gt(0).sum()
                ),
                "anaerobe_v2_mean_score": mean_or_nan(group["anaerobe_v2_score"]),
                "anaerobe_v2_median_score": median_or_nan(group["anaerobe_v2_score"]),
                "anaerobe_v2_positive_profile_count": int(
                    group["anaerobe_v2_positive_trait_count"].gt(0).sum()
                ),
            }
        )
    return pd.DataFrame.from_records(records)


def build_paired_contrasts(scores: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []

    for (study_id, experiment_id), group in scores.groupby(
        ["study_id", "experiment_id"],
        sort=False,
        dropna=False,
    ):
        increased = group[group["inc_dec"] == "inc"]
        decreased = group[group["inc_dec"] == "dec"]
        if increased.empty or decreased.empty:
            continue

        aerobe_inc = mean_or_nan(increased["aerobe_v2_score"])
        aerobe_dec = mean_or_nan(decreased["aerobe_v2_score"])
        anaerobe_inc = mean_or_nan(increased["anaerobe_v2_score"])
        anaerobe_dec = mean_or_nan(decreased["anaerobe_v2_score"])
        aerobe_contrast = aerobe_inc - aerobe_dec
        anaerobe_contrast = anaerobe_dec - anaerobe_inc

        records.append(
            {
                "study_id": study_id,
                "experiment_id": experiment_id,
                "condition": PD_CONDITION,
                "body_site_group": BODY_SITE_GROUP,
                "title": join_unique(group["title"]),
                "publication_year": join_unique(group["publication_year"]),
                "increased_profile_count": len(increased),
                "decreased_profile_count": len(decreased),
                "aerobe_v2_increased_mean": aerobe_inc,
                "aerobe_v2_decreased_mean": aerobe_dec,
                "aerobe_v2_increased_minus_decreased": aerobe_contrast,
                "anaerobe_v2_increased_mean": anaerobe_inc,
                "anaerobe_v2_decreased_mean": anaerobe_dec,
                "anaerobe_v2_decreased_minus_increased": anaerobe_contrast,
                "combined_reciprocal_signal": aerobe_contrast + anaerobe_contrast,
            }
        )

    paired = pd.DataFrame.from_records(records)
    if paired.empty:
        return paired

    return paired.sort_values(
        ["combined_reciprocal_signal", "study_id", "experiment_id"],
        ascending=[False, True, True],
    )


def safe_mannwhitney(focus: pd.Series, comparator: pd.Series) -> float:
    if len(focus) == 0 or len(comparator) == 0:
        return math.nan
    try:
        return float(stats.mannwhitneyu(focus, comparator, alternative="two-sided").pvalue)
    except ValueError:
        return math.nan


def safe_wilcoxon(contrast: pd.Series) -> float:
    values = contrast.dropna()
    if len(values) == 0:
        return math.nan
    try:
        return float(stats.wilcoxon(values, alternative="two-sided").pvalue)
    except ValueError:
        return math.nan


def summarize_paired_tests(scores: pd.DataFrame, paired: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    test_specs = [
        {
            "analysis": AEROBE_CONFIG.label,
            "score_column": "aerobe_v2_score",
            "contrast_column": "aerobe_v2_increased_minus_decreased",
            "focus_direction": "inc",
            "comparator_direction": "dec",
            "positive_effect": "increased profiles have higher aerobe enrichment scores",
        },
        {
            "analysis": ANAEROBE_CONFIG.label,
            "score_column": "anaerobe_v2_score",
            "contrast_column": "anaerobe_v2_decreased_minus_increased",
            "focus_direction": "dec",
            "comparator_direction": "inc",
            "positive_effect": "decreased profiles have higher anaerobe depletion scores",
        },
    ]

    for spec in test_specs:
        focus = scores.loc[scores["inc_dec"] == spec["focus_direction"], spec["score_column"]]
        comparator = scores.loc[
            scores["inc_dec"] == spec["comparator_direction"],
            spec["score_column"],
        ]
        records.append(
            {
                "analysis": spec["analysis"],
                "comparison": "profile-level focus direction vs comparator direction",
                "positive_effect_interpretation": spec["positive_effect"],
                "focus_direction": spec["focus_direction"],
                "comparator_direction": spec["comparator_direction"],
                "focus_count": len(focus),
                "comparator_count": len(comparator),
                "focus_mean": mean_or_nan(focus),
                "comparator_mean": mean_or_nan(comparator),
                "focus_minus_comparator": mean_or_nan(focus) - mean_or_nan(comparator),
                "test": "Mann-Whitney U",
                "p_value": safe_mannwhitney(focus, comparator),
            }
        )

        if not paired.empty:
            contrast = paired[spec["contrast_column"]]
            records.append(
                {
                    "analysis": spec["analysis"],
                    "comparison": "within-experiment paired contrast",
                    "positive_effect_interpretation": spec["positive_effect"],
                    "focus_direction": spec["focus_direction"],
                    "comparator_direction": spec["comparator_direction"],
                    "focus_count": len(contrast),
                    "comparator_count": len(contrast),
                    "focus_mean": mean_or_nan(contrast),
                    "comparator_mean": 0.0,
                    "focus_minus_comparator": mean_or_nan(contrast),
                    "test": "Wilcoxon signed-rank vs zero",
                    "p_value": safe_wilcoxon(contrast),
                }
            )

    return pd.DataFrame.from_records(records)


def load_target_taxa_table(config: AnalysisConfig) -> pd.DataFrame:
    table = read_tsv(config.candidate_taxa_path)
    subset = table[
        (table["condition"] == PD_CONDITION) & (table["body_site_group"] == BODY_SITE_GROUP)
    ].copy()
    if subset.empty:
        return subset

    subset["analysis"] = config.key
    subset["analysis_label"] = config.label
    subset["focus_direction"] = config.focus_label
    subset["comparator_direction"] = config.comparator_label
    subset["focus_profile_count"] = subset[f"{config.focus_label}_profile_count"]
    subset["comparator_profile_count"] = subset[f"{config.comparator_label}_profile_count"]
    subset["focus_study_count"] = subset[f"{config.focus_label}_study_count"]
    subset["comparator_study_count"] = subset[f"{config.comparator_label}_study_count"]
    subset["focus_bias"] = subset[config.bias_column]
    subset["net_increased_minus_decreased"] = (
        subset["increased_profile_count"] - subset["decreased_profile_count"]
    )
    subset["dominant_direction"] = np.select(
        [
            subset["increased_profile_count"] > subset["decreased_profile_count"],
            subset["decreased_profile_count"] > subset["increased_profile_count"],
        ],
        ["increased", "decreased"],
        default="balanced",
    )
    subset["dominant_profile_count"] = subset[
        ["increased_profile_count", "decreased_profile_count"]
    ].max(axis=1)
    subset["dominant_study_count"] = np.where(
        subset["dominant_direction"] == "increased",
        subset["increased_study_count"],
        np.where(
            subset["dominant_direction"] == "decreased",
            subset["decreased_study_count"],
            subset[["increased_study_count", "decreased_study_count"]].max(axis=1),
        ),
    )

    columns = [
        "analysis",
        "analysis_label",
        "target_taxon_name",
        "condition",
        "body_site_group",
        "decreased_profile_count",
        "increased_profile_count",
        "decreased_study_count",
        "increased_study_count",
        "total_profile_count",
        "focus_direction",
        "focus_profile_count",
        "focus_study_count",
        "comparator_direction",
        "comparator_profile_count",
        "comparator_study_count",
        "focus_bias",
        "net_increased_minus_decreased",
        "dominant_direction",
        "dominant_profile_count",
        "dominant_study_count",
    ]
    return subset[columns].sort_values(
        ["analysis", "focus_profile_count", "focus_study_count", "target_taxon_name"],
        ascending=[True, False, False, True],
    )


def build_target_taxa_direction() -> pd.DataFrame:
    taxa_tables = [load_target_taxa_table(config) for config in ANALYSIS_CONFIGS]
    taxa = pd.concat(taxa_tables, ignore_index=True)
    return taxa.sort_values(
        ["analysis", "focus_profile_count", "focus_study_count", "target_taxon_name"],
        ascending=[True, False, False, True],
    )


def summarize_target_taxa(taxa: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for analysis, group in taxa.groupby("analysis_label", sort=False):
        records.append(
            {
                "analysis_label": analysis,
                "taxa_count": len(group),
                "total_profile_count": int(group["total_profile_count"].sum()),
                "focus_profile_count": int(group["focus_profile_count"].sum()),
                "comparator_profile_count": int(group["comparator_profile_count"].sum()),
                "focus_taxa_count": int((group["focus_profile_count"] > 0).sum()),
                "perfect_focus_bias_taxa_count": int((group["focus_bias"] == 1.0).sum()),
                "top_focus_taxon": group.iloc[0]["target_taxon_name"],
                "top_focus_profile_count": int(group.iloc[0]["focus_profile_count"]),
                "top_focus_study_count": int(group.iloc[0]["focus_study_count"]),
            }
        )
    return pd.DataFrame.from_records(records)


def save_figure(fig: plt.Figure, output_path: Path) -> None:
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(output_path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_profile_score_contrast(scores: pd.DataFrame, output_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), constrained_layout=True)
    rng = np.random.default_rng(42)
    specs = [
        ("aerobe_v2_score", "Aerobe enrichment v2 score"),
        ("anaerobe_v2_score", "Anaerobe depletion v2 score"),
    ]

    for ax, (score_column, title) in zip(axes, specs, strict=True):
        data = [
            scores.loc[scores["inc_dec"] == "dec", score_column].to_numpy(),
            scores.loc[scores["inc_dec"] == "inc", score_column].to_numpy(),
        ]
        box = ax.boxplot(
            data,
            positions=[0, 1],
            widths=0.48,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "#111827", "linewidth": 1.4},
        )
        for patch, color in zip(
            box["boxes"],
            [DIRECTION_COLORS["dec"], DIRECTION_COLORS["inc"]],
            strict=True,
        ):
            patch.set_facecolor(color)
            patch.set_alpha(0.18)
            patch.set_edgecolor(color)
            patch.set_linewidth(1.2)

        for position, direction in enumerate(("dec", "inc")):
            values = scores.loc[scores["inc_dec"] == direction, score_column].to_numpy()
            jitter = rng.normal(0, 0.045, size=len(values))
            ax.scatter(
                np.full(len(values), position) + jitter,
                values,
                s=14,
                alpha=0.55,
                color=DIRECTION_COLORS[direction],
                edgecolors="none",
            )
            ax.scatter(
                [position],
                [np.mean(values)],
                marker="D",
                s=44,
                color="#111827",
                zorder=4,
            )

        ax.axhline(0, color="#94a3b8", linewidth=1.0, linestyle="--")
        ax.set_xticks([0, 1], ["Decreased", "Increased"])
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_ylabel("Mean selected-trait z-score")
        ax.grid(axis="y", color="#e5e7eb", linewidth=0.8)
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Parkinson's disease feces profile-score contrast", fontsize=14, fontweight="bold")
    save_figure(fig, output_path)


def plot_paired_contrasts(paired: pd.DataFrame, output_path: Path) -> None:
    if paired.empty:
        return

    fig, ax = plt.subplots(figsize=(9.4, 4.8), constrained_layout=True)
    rng = np.random.default_rng(84)
    specs = [
        (
            "Aerobe: increased - decreased",
            "aerobe_v2_increased_minus_decreased",
            "#c2410c",
        ),
        (
            "Anaerobe: decreased - increased",
            "anaerobe_v2_decreased_minus_increased",
            "#2563eb",
        ),
        ("Combined reciprocal signal", "combined_reciprocal_signal", "#334155"),
    ]

    for y_position, (_label, column, color) in enumerate(specs):
        values = paired[column].dropna().to_numpy()
        jitter = rng.normal(0, 0.035, size=len(values))
        ax.scatter(
            values,
            np.full(len(values), y_position) + jitter,
            s=22,
            alpha=0.58,
            color=color,
            edgecolors="none",
        )
        ax.scatter(
            [np.mean(values)],
            [y_position],
            marker="D",
            s=70,
            color="#111827",
            zorder=4,
        )
        ax.text(
            np.mean(values),
            y_position + 0.2,
            f"mean {np.mean(values):.3f}",
            ha="center",
            va="bottom",
            fontsize=9,
            color="#111827",
        )

    ax.axvline(0, color="#111827", linewidth=1.1)
    ax.set_yticks(range(len(specs)), [spec[0] for spec in specs])
    ax.set_xlabel("Within-experiment score contrast")
    ax.set_title("Parkinson's disease feces paired experiment contrasts", fontweight="bold")
    ax.grid(axis="x", color="#e5e7eb", linewidth=0.8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    save_figure(fig, output_path)


def plot_target_taxa_direction(taxa: pd.DataFrame, output_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.6), constrained_layout=True)
    analysis_order = [AEROBE_CONFIG.label, ANAEROBE_CONFIG.label]

    for ax, analysis_label in zip(axes, analysis_order, strict=True):
        group = taxa[taxa["analysis_label"] == analysis_label].copy()
        group = group.sort_values(
            ["total_profile_count", "focus_profile_count", "target_taxon_name"],
            ascending=[True, True, False],
        )
        y_positions = np.arange(len(group))
        ax.barh(
            y_positions,
            -group["decreased_profile_count"],
            color="#2563eb",
            alpha=0.72,
            label="Decreased",
        )
        ax.barh(
            y_positions,
            group["increased_profile_count"],
            color="#c2410c",
            alpha=0.72,
            label="Increased",
        )
        ax.set_yticks(y_positions, group["target_taxon_name"])
        ax.axvline(0, color="#111827", linewidth=1.1)
        ax.set_title(analysis_label, fontsize=12, fontweight="bold")
        ax.set_xlabel("Profile count")
        ax.grid(axis="x", color="#e5e7eb", linewidth=0.8)
        ax.spines[["top", "right"]].set_visible(False)
        max_count = int(
            max(group["decreased_profile_count"].max(), group["increased_profile_count"].max())
        )
        ax.set_xlim(-(max_count + 2), max_count + 2)

    axes[1].legend(loc="lower right", frameon=False)
    fig.suptitle("Parkinson's disease feces target-taxa directions", fontsize=14, fontweight="bold")
    save_figure(fig, output_path)


def plot_anaerobe_taxa_support(taxa: pd.DataFrame, output_path: Path) -> None:
    group = taxa[taxa["analysis"] == ANAEROBE_CONFIG.key].copy()
    group = group.sort_values(
        ["decreased_profile_count", "decreased_study_count", "target_taxon_name"],
        ascending=[True, True, False],
    )

    fig, ax = plt.subplots(figsize=(8.8, 4.9), constrained_layout=True)
    y_positions = np.arange(len(group))
    ax.barh(
        y_positions,
        group["decreased_profile_count"],
        color="#2563eb",
        alpha=0.74,
        label="Decreased profiles",
    )
    ax.barh(
        y_positions,
        group["increased_profile_count"],
        color="#c2410c",
        alpha=0.74,
        label="Increased profiles",
    )
    ax.scatter(
        group["decreased_study_count"],
        y_positions,
        marker="D",
        s=52,
        color="#111827",
        label="Decreased studies",
        zorder=4,
    )
    ax.set_yticks(y_positions, group["target_taxon_name"])
    ax.set_xlabel("Profile count; black markers show decreased-study count")
    ax.set_title("PD feces anaerobe-depletion target taxa", fontweight="bold")
    ax.grid(axis="x", color="#e5e7eb", linewidth=0.8)
    ax.legend(loc="lower right", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, output_path)


def top_taxa_lines(taxa: pd.DataFrame, analysis_key: str, *, limit: int = 5) -> list[str]:
    group = taxa[taxa["analysis"] == analysis_key].head(limit)
    return [
        (
            f"- {row.target_taxon_name}: {row.decreased_profile_count} decreased / "
            f"{row.increased_profile_count} increased profiles; "
            f"{row.decreased_study_count} decreased / "
            f"{row.increased_study_count} increased studies."
        )
        for row in group.itertuples(index=False)
    ]


def write_readme(
    scores: pd.DataFrame,
    direction_summary: pd.DataFrame,
    tests: pd.DataFrame,
    paired: pd.DataFrame,
    taxa: pd.DataFrame,
    taxa_summary: pd.DataFrame,
) -> None:
    inc_count = int((scores["inc_dec"] == "inc").sum())
    dec_count = int((scores["inc_dec"] == "dec").sum())
    study_count = scores["study_id"].nunique()
    experiment_count = scores["experiment_id"].nunique()
    paired_studies = paired["study_id"].nunique() if not paired.empty else 0

    aerobe_test = tests[tests["analysis"] == AEROBE_CONFIG.label].iloc[0]
    anaerobe_test = tests[tests["analysis"] == ANAEROBE_CONFIG.label].iloc[0]
    aerobe_pair = tests[
        (tests["analysis"] == AEROBE_CONFIG.label)
        & (tests["comparison"] == "within-experiment paired contrast")
    ]
    anaerobe_pair = tests[
        (tests["analysis"] == ANAEROBE_CONFIG.label)
        & (tests["comparison"] == "within-experiment paired contrast")
    ]

    lines = [
        "# Parkinson's Disease Trait-Profile Analysis",
        "",
        "Generated by `scripts/analyze_parkinsons_disease_trait_profiles.py`.",
        "",
        "## Scope",
        "",
        f"- Filter: exact condition `{PD_CONDITION}` and body-site group `{BODY_SITE_GROUP}`.",
        f"- Profiles: {len(scores)} total; {inc_count} increased and {dec_count} decreased.",
        f"- Studies: {study_count}; experiments: {experiment_count}.",
        f"- Paired increased/decreased experiment contrasts: {len(paired)} experiments "
        f"from {paired_studies} studies.",
        "",
        "## Approach 2: PD Feces-Only Trait-Score Contrast",
        "",
        (
            "- Aerobe enrichment profile-level effect "
            f"(increased minus decreased mean): "
            f"{fmt_float(aerobe_test.focus_minus_comparator)}; "
            f"Mann-Whitney p={fmt_p_value(aerobe_test.p_value)}."
        ),
        (
            "- Anaerobe depletion profile-level effect "
            f"(decreased minus increased mean): "
            f"{fmt_float(anaerobe_test.focus_minus_comparator)}; "
            f"Mann-Whitney p={fmt_p_value(anaerobe_test.p_value)}."
        ),
    ]

    if not aerobe_pair.empty and not anaerobe_pair.empty:
        lines.extend(
            [
                (
                    "- Aerobe within-experiment contrast mean "
                    f"(increased minus decreased): "
                    f"{fmt_float(aerobe_pair.iloc[0].focus_minus_comparator)}; "
                    f"Wilcoxon p={fmt_p_value(aerobe_pair.iloc[0].p_value)}."
                ),
                (
                    "- Anaerobe within-experiment contrast mean "
                    f"(decreased minus increased): "
                    f"{fmt_float(anaerobe_pair.iloc[0].focus_minus_comparator)}; "
                    f"Wilcoxon p={fmt_p_value(anaerobe_pair.iloc[0].p_value)}."
                ),
            ]
        )

    lines.extend(
        [
            "",
            "See `parkinsons_feces_profile_trait_score_contrast.png` and "
            "`parkinsons_feces_paired_trait_contrast.png` for the corresponding diagrams.",
            "",
            "## Approach 3: PD Target-Taxa Direction Consistency",
            "",
            "Top aerobe-enrichment target taxa:",
            *top_taxa_lines(taxa, AEROBE_CONFIG.key),
            "",
            "Top anaerobe-depletion target taxa:",
            *top_taxa_lines(taxa, ANAEROBE_CONFIG.key),
            "",
            "See `parkinsons_target_taxa_direction_bar.png` and "
            "`parkinsons_anaerobe_depletion_taxa_support.png` for the corresponding diagrams.",
            "",
            "## Diagram Guide",
            "",
            (
                "- `parkinsons_feces_profile_trait_score_contrast.png`: compares "
                "profile-level v2 scores between decreased and increased PD feces "
                "profiles. Each dot is one BugSigDB study/experiment/direction "
                "profile, boxes show the distribution, black diamonds show means, "
                "and the dashed line marks zero. The aerobe panel is read as "
                "higher scores in increased profiles; the anaerobe panel is read as "
                "higher scores in decreased profiles."
            ),
            (
                "- `parkinsons_feces_paired_trait_contrast.png`: keeps only "
                "study/experiment pairs that have both increased and decreased PD "
                "feces profiles. Dots are within-experiment contrasts, black "
                "diamonds are mean contrasts, and positive values mean the expected "
                "reciprocal pattern: aerobe signal higher in increased profiles and "
                "anaerobe signal higher in decreased profiles."
            ),
            (
                "- `parkinsons_target_taxa_direction_bar.png`: shows PD feces "
                "target taxa from the aerobe-enrichment and anaerobe-depletion "
                "candidate-taxa tables. Orange bars to the right count increased "
                "profiles; blue bars to the left count decreased profiles. Longer "
                "one-sided bars indicate stronger direction consistency."
            ),
            (
                "- `parkinsons_anaerobe_depletion_taxa_support.png`: focuses on "
                "anaerobe-depletion target taxa. Horizontal bars compare decreased "
                "and increased profile support, while black diamonds mark the number "
                "of decreased-supporting studies, separating repeated-study support "
                "from raw profile count."
            ),
            "",
            "## Figure Methods",
            "",
            (
                "- Source data: profile z-score matrix "
                "`study_trait_profile_clustered_heatmap_matrix.tsv`, profile metadata "
                "`study_trait_profile_study_metadata.tsv`, trait metadata "
                "`study_trait_profile_trait_metadata.tsv`, v2 guild trait selections, "
                "and v2 candidate-taxa direction tables."
            ),
            (
                f"- Filter: exact condition `{PD_CONDITION}` and body-site group "
                f"`{BODY_SITE_GROUP}`. Increased and decreased directions are the "
                "BugSigDB profile directions already present in the study-profile "
                "metadata."
            ),
            (
                "- Profile scores: for each profile, the aerobe and anaerobe v2 "
                "scores are the row-wise mean z-score across the selected trait "
                "columns from the corresponding `guild_trait_selection.tsv` file. "
                "Positive-trait counts record how many selected trait z-scores are "
                "greater than zero."
            ),
            (
                "- Paired contrasts: matched study/experiment pairs are retained only "
                "when both increased and decreased PD feces profiles are present. "
                "Aerobe contrast is increased minus decreased; anaerobe contrast is "
                "decreased minus increased; combined reciprocal signal is their sum."
            ),
            (
                "- Statistical summaries: profile-level comparisons use two-sided "
                "Mann-Whitney U tests. Paired contrasts use two-sided Wilcoxon "
                "signed-rank tests against zero. These tests summarize recurrence "
                "patterns and are not causal models."
            ),
            (
                "- Target-taxa diagrams: candidate taxa are read from the aerobe and "
                "anaerobe v2 `candidate_taxa_condition_direction.tsv` files and "
                "filtered to PD feces rows. Counts are profile and unique-study "
                "support for increased versus decreased BugSigDB directions."
            ),
            (
                "- Rendering: diagrams are generated with Matplotlib and saved as both "
                "PNG and SVG. Orange encodes increased-profile support, blue encodes "
                "decreased-profile support, and black diamonds mark means or "
                "decreased-study support depending on the panel."
            ),
            "",
            "The consolidated figure-methods entry is also documented in "
            "`../../figure_methods.md`.",
            "",
            "## Tables",
            "",
            "- `parkinsons_feces_profile_trait_scores.tsv`: profile-level v2 scores.",
            "- `parkinsons_feces_trait_score_summary.tsv`: increased/decreased score summary.",
            "- `parkinsons_feces_paired_trait_contrasts.tsv`: within-experiment contrasts.",
            "- `parkinsons_feces_trait_score_tests.tsv`: profile-level and paired tests.",
            "- `parkinsons_target_taxa_direction.tsv`: PD target taxa by direction.",
            "- `parkinsons_target_taxa_direction_summary.tsv`: target taxa rollup.",
            "",
            "## Notes",
            "",
            (
                "These outputs quantify recurrence and direction consistency in BugSigDB "
                "signatures. They are not causal estimates and should be interpreted as "
                "study-supported signals to prioritize follow-up taxa and traits."
            ),
            "",
            "Direction summary:",
            "",
            direction_summary.to_string(index=False),
            "",
            "Target-taxa summary:",
            "",
            taxa_summary.to_string(index=False),
        ]
    )

    (OUTPUT_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    matrix, metadata = load_study_heatmap()
    scored = add_trait_scores(matrix, metadata)
    scores = parkinsons_feces_scores(scored)
    direction_summary = summarize_direction_scores(scores)
    paired = build_paired_contrasts(scores)
    tests = summarize_paired_tests(scores, paired)
    taxa = build_target_taxa_direction()
    taxa_summary = summarize_target_taxa(taxa)

    write_tsv(scores, OUTPUT_DIR / "parkinsons_feces_profile_trait_scores.tsv")
    write_tsv(direction_summary, OUTPUT_DIR / "parkinsons_feces_trait_score_summary.tsv")
    write_tsv(paired, OUTPUT_DIR / "parkinsons_feces_paired_trait_contrasts.tsv")
    write_tsv(tests, OUTPUT_DIR / "parkinsons_feces_trait_score_tests.tsv")
    write_tsv(taxa, OUTPUT_DIR / "parkinsons_target_taxa_direction.tsv")
    write_tsv(taxa_summary, OUTPUT_DIR / "parkinsons_target_taxa_direction_summary.tsv")

    plot_profile_score_contrast(
        scores,
        OUTPUT_DIR / "parkinsons_feces_profile_trait_score_contrast.png",
    )
    plot_paired_contrasts(
        paired,
        OUTPUT_DIR / "parkinsons_feces_paired_trait_contrast.png",
    )
    plot_target_taxa_direction(
        taxa,
        OUTPUT_DIR / "parkinsons_target_taxa_direction_bar.png",
    )
    plot_anaerobe_taxa_support(
        taxa,
        OUTPUT_DIR / "parkinsons_anaerobe_depletion_taxa_support.png",
    )
    write_readme(scores, direction_summary, tests, paired, taxa, taxa_summary)

    print(f"Wrote Parkinson's disease analysis outputs to {OUTPUT_DIR.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
