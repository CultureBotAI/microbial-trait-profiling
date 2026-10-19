"""Investigate anaerobe-only v2 trait enrichment in decreased profiles."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle, Patch, PathPatch
from matplotlib.path import Path as MplPath

from mtp.body_site_diagram import anatomical_site_base
from mtp.composition_denominators import (
    DENOMINATOR_DEFINITION,
    load_retained_profiles,
    trait_fractions,
)
from mtp.human_host_corpus import (
    excluded_host_note,
    excluded_host_summary,
    join_human_host_profiles,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
OUTPUT_DIR = VIS_DIR / "trait-profile-investigation" / "anaerobe-depletion-v2"
PROFILE_MATRIX_PATH = REPO_ROOT / "data/derived/study_trait_profiles_wide.tsv"

GUILD_TRAITS = {
    "METPO:1000603": "anaerobic",
    "METPO:1000607": "obligately anaerobic",
}
GUILD_TRAIT_SELECTION_JUSTIFICATIONS = {
    "METPO:1000603": (
        "Anaerobic oxygen-tolerance counterpart to the aerobe v2 trait "
        "set; selected for the decreased-profile comparison."
    ),
    "METPO:1000607": (
        "Stricter anaerobe oxygen-tolerance concept; selected to keep the decreased-profile "
        "trait set parallel to oxygen-tolerance rather than fermentation metabolism."
    ),
}
KG_MICROBE_TAXON_UNIVERSE_TRAITS = {
    "METPO:1000603": "anaerobic",
    "METPO:1000607": "obligately anaerobic",
}
KG_MICROBE_TAXON_UNIVERSE_ASSOCIATED_TRAITS: dict[str, str] = {}

TARGET_TAXA = {
    "NCBITaxon:853": "Faecalibacterium prausnitzii",
    "NCBITaxon:39488": "Anaerobutyricum hallii",
    "NCBITaxon:39485": "Lachnospira eligens",
    "NCBITaxon:39491": "Agathobacter rectalis",
    "NCBITaxon:166486": "Roseburia intestinalis",
    "NCBITaxon:40518": "Ruminococcus bromii",
    "NCBITaxon:360807": "Roseburia inulinivorans",
}

EXCLUDED_CONDITION_SUMMARY_LABELS = {
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
# A paired d (and its paired t-test p) from fewer matched profile pairs is reported as NA.
MIN_PAIRED_SMD_PROFILE_PAIRS = 5
INFERRED_GUILD_TOP_TRAITS_PER_GROUP = 12
INFERRED_GUILD_TOP_TAXA_PER_TRAIT = 12
INFERRED_GUILD_TOP_TAXA_PER_GROUP = 12
MIN_KG_MICROBE_UNIVERSE_CORE_TRAITS = 1
GUILD_SPECIFIC_NOTE = (
    "Guild-specific: yes - uses the anaerobe v2 trait-enrichment score "
    "(anaerobic, obligately anaerobic)."
)
NON_GUILD_SPECIFIC_TARGET_TAXA_NOTE = (
    "Guild-specific: no - fixed target-taxon count view; "
    "not recomputed from the trait-enrichment score."
)
TARGET_TAXON_COLORS = {
    "Faecalibacterium prausnitzii": "#2563eb",
    "Lachnospira eligens": "#d97706",
    "Agathobacter rectalis": "#059669",
    "Anaerobutyricum hallii": "#dc2626",
    "Roseburia intestinalis": "#7c3aed",
    "Roseburia inulinivorans": "#0891b2",
    "Ruminococcus bromii": "#64748b",
}
TARGET_TAXON_LEGEND_LABELS = {
    "Faecalibacterium prausnitzii": "F. prausnitzii",
    "Lachnospira eligens": "L. eligens",
    "Agathobacter rectalis": "A. rectalis",
    "Anaerobutyricum hallii": "A. hallii",
    "Roseburia intestinalis": "R. intestinalis",
    "Roseburia inulinivorans": "R. inulinivorans",
    "Ruminococcus bromii": "R. bromii",
}
DIRECTION_LABELS = {
    "inc": "increased profiles",
    "dec": "decreased profiles",
}
DIRECTION_COLORS = {
    "inc": "#dc2626",
    "dec": "#2563eb",
}
PROFILE_COUNT_BIN_DEFINITIONS = (
    (1, 1, "1"),
    (2, 2, "2"),
    (3, 4, "3-4"),
    (5, 9, "5-9"),
    (10, 19, "10-19"),
    (20, 49, "20-49"),
    (50, 99, "50-99"),
    (100, None, "100+"),
)
TRAIT_PROFILE_COUNT_BIN_DEFINITIONS = (
    (1, 99, "1-99"),
    (100, 499, "100-499"),
    (500, 999, "500-999"),
    (1000, 1499, "1,000-1,499"),
    (1500, 1999, "1,500-1,999"),
    (2000, 2499, "2,000-2,499"),
    (2500, 2999, "2,500-2,999"),
    (3000, None, "3,000+"),
)


def add_specificity_note(ax: plt.Axes, note: str) -> None:
    """Add a visible guild-specificity note below the plot title."""
    ax.text(
        0,
        1.02,
        note,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        color="#475569",
    )


def read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def clean_text(value: object) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan" or text.upper() == "NA":
        return "unknown"
    return text


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return math.nan
    return numerator / denominator


def _incomplete_beta_continued_fraction(a: float, b: float, x: float) -> float:
    max_iterations = 200
    epsilon = 3e-14
    tiny = 1e-300

    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d

    for iteration in range(1, max_iterations + 1):
        m2 = 2 * iteration
        numerator = iteration * (b - iteration) * x
        denominator = (qam + m2) * (a + m2)
        d = 1.0 + numerator / denominator * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + numerator / denominator / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c

        numerator = -(a + iteration) * (qab + iteration) * x
        denominator = (a + m2) * (qap + m2)
        d = 1.0 + numerator / denominator * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + numerator / denominator / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < epsilon:
            break

    return h


def regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0

    log_beta_factor = (
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x)
    )
    beta_factor = math.exp(log_beta_factor)
    if x < (a + 1.0) / (a + b + 2.0):
        return beta_factor * _incomplete_beta_continued_fraction(a, b, x) / a
    return 1.0 - beta_factor * _incomplete_beta_continued_fraction(b, a, 1.0 - x) / b


def paired_t_test_p_value(
    mean_difference: float,
    difference_sd: float,
    pair_count: int,
) -> float:
    if pair_count < MIN_PAIRED_SMD_PROFILE_PAIRS or pd.isna(difference_sd) or difference_sd <= 0:
        return math.nan
    t_statistic = mean_difference / (difference_sd / math.sqrt(pair_count))
    degrees_of_freedom = pair_count - 1
    x = degrees_of_freedom / (degrees_of_freedom + t_statistic**2)
    p_value = regularized_incomplete_beta(degrees_of_freedom / 2.0, 0.5, x)
    return max(0.0, min(1.0, p_value))


def paired_guild_enrichment_stats(group: pd.DataFrame) -> dict[str, float | int]:
    """Calculate paired decreased-vs-increased trait enrichment within study/experiment."""
    paired = group.pivot_table(
        index=["study_id", "experiment_id"],
        columns="inc_dec",
        values="guild_z_mean",
        aggfunc="mean",
    )
    if not {"inc", "dec"}.issubset(paired.columns):
        return {
            "paired_profile_pair_count": 0,
            "paired_mean_guild_z_difference": math.nan,
            "paired_guild_z_difference_sd": math.nan,
            "paired_standardized_mean_difference": math.nan,
            "paired_t_test_p_value": math.nan,
        }

    paired = paired.dropna(subset=["inc", "dec"])
    differences = paired["dec"] - paired["inc"]
    pair_count = len(differences)
    mean_difference = differences.mean() if pair_count else math.nan
    difference_sd = differences.std(ddof=1) if pair_count > 1 else math.nan
    smd = math.nan
    p_value = math.nan
    if (
        pair_count >= MIN_PAIRED_SMD_PROFILE_PAIRS
        and not pd.isna(difference_sd)
        and difference_sd > 0
    ):
        smd = mean_difference / difference_sd
        p_value = paired_t_test_p_value(mean_difference, difference_sd, pair_count)

    return {
        "paired_profile_pair_count": pair_count,
        "paired_mean_guild_z_difference": mean_difference,
        "paired_guild_z_difference_sd": difference_sd,
        "paired_standardized_mean_difference": smd,
        "paired_t_test_p_value": p_value,
    }


def profile_id_from_taxa_row(row: pd.Series) -> str:
    direction = "inc" if row["direction"] == "increased" else "dec"
    return f"{row['study_id']}|{row['experiment_id']}|{direction}"


def load_study_heatmap() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    matrix = read_tsv(VIS_DIR / "study_trait_profile_clustered_heatmap_matrix.tsv", index_col=0)
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    metadata = read_tsv(VIS_DIR / "study_trait_profile_study_metadata.tsv", index_col=0)
    trait_metadata = read_tsv(VIS_DIR / "study_trait_profile_trait_metadata.tsv", index_col=0)

    if len(matrix) != len(metadata):
        raise ValueError("Heatmap matrix and metadata have different profile row counts")
    if (
        "denominator_definition" not in metadata.columns
        or not metadata["denominator_definition"].eq(DENOMINATOR_DEFINITION).all()
    ):
        raise ValueError(
            "Heatmap metadata is missing the harmonized non-MediaDive denominator marker"
        )
    matrix_keys = matrix.index.to_series().map(lambda value: str(value).split(" | ", 1)[0])
    matrix_pair_keys = matrix.index.to_series().map(
        lambda value: "/".join(part.strip() for part in str(value).split(" | ")[:2])
    )
    expected_keys = (
        metadata[["study_id", "experiment_id"]].astype(str).agg("/".join, axis=1).to_numpy()
    )
    if not (
        matrix.index.equals(metadata.index)
        or np.array_equal(matrix_keys.to_numpy(), expected_keys)
        or np.array_equal(matrix_pair_keys.to_numpy(), expected_keys)
    ):
        raise ValueError("Heatmap matrix row keys do not match profile metadata")
    matrix.index = metadata.index

    if "trait_id" in trait_metadata.columns:
        trait_metadata = trait_metadata.set_index("trait_id", drop=False)
    column_to_id = {
        row["display_label"]: str(trait_id) for trait_id, row in trait_metadata.iterrows()
    }
    column_labels = {
        row["display_label"]: row["object_name"] for _, row in trait_metadata.iterrows()
    }
    if not set(matrix.columns).issubset(column_to_id):
        raise ValueError("Heatmap matrix has columns missing from trait metadata")
    counts = load_retained_profiles(PROFILE_MATRIX_PATH, metadata)
    trait_ids = [column_to_id[column] for column in matrix.columns]
    expected_comp = trait_fractions(counts).loc[metadata.index, trait_ids]
    expected_scaled = (expected_comp * 1000).apply(lambda column: column.map(math.log1p))
    expected_z = (
        (expected_scaled - expected_scaled.mean(axis=0))
        / expected_scaled.std(axis=0, ddof=0).replace(0, float("nan"))
    ).fillna(0.0)
    expected_z.columns = matrix.columns
    if not np.allclose(matrix.to_numpy(dtype=float), expected_z.to_numpy(dtype=float), atol=1e-9):
        raise ValueError(
            "Stale or incompatible study_trait_profile_clustered_heatmap_matrix.tsv; "
            "regenerate from the non-MediaDive raw profile matrix."
        )
    matrix = matrix.rename(columns=column_labels)
    return matrix, metadata, trait_metadata


def add_guild_scores(matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    available_traits = [trait for trait in GUILD_TRAITS.values() if trait in matrix.columns]
    if not available_traits:
        msg = "No requested anaerobe v2 trait-enrichment traits were found in the matrix."
        raise ValueError(msg)
    scored = metadata.copy()
    scored["guild_z_mean"] = matrix[available_traits].mean(axis=1)
    scored["anaerobe_z_mean"] = scored["guild_z_mean"]
    scored["enriched_guild_trait_count"] = matrix[available_traits].gt(0).sum(axis=1)
    scored["has_enriched_guild_trait"] = scored["enriched_guild_trait_count"] > 0
    for trait_name in available_traits:
        scored[f"z_{trait_name}"] = matrix[trait_name]
    return scored


def summarize_direction(metadata: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for trait_name in GUILD_TRAITS.values():
        column = f"z_{trait_name}"
        inc_mean = metadata.loc[metadata["inc_dec"] == "inc", column].mean()
        dec_mean = metadata.loc[metadata["inc_dec"] == "dec", column].mean()
        rows.append(
            {
                "trait": trait_name,
                "increased_profile_mean_z": inc_mean,
                "decreased_profile_mean_z": dec_mean,
                "decreased_minus_increased_z": dec_mean - inc_mean,
            }
        )
    for score_name in [
        "guild_z_mean",
    ]:
        inc_mean = metadata.loc[metadata["inc_dec"] == "inc", score_name].mean()
        dec_mean = metadata.loc[metadata["inc_dec"] == "dec", score_name].mean()
        rows.append(
            {
                "trait": score_name,
                "increased_profile_mean_z": inc_mean,
                "decreased_profile_mean_z": dec_mean,
                "decreased_minus_increased_z": dec_mean - inc_mean,
            }
        )
    return pd.DataFrame(rows)


def summarize_group(metadata: pd.DataFrame, group_column: str, min_profiles: int) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for group_value, group in metadata.groupby(group_column, dropna=False):
        if len(group) < min_profiles:
            continue
        inc = group[group["inc_dec"] == "inc"]
        dec = group[group["inc_dec"] == "dec"]
        record = {
            group_column: clean_text(group_value),
            "profile_count": len(group),
            "study_count": group["study_id"].nunique(),
            "increased_profile_count": len(inc),
            "decreased_profile_count": len(dec),
            "mean_guild_z": group["guild_z_mean"].mean(),
            "increased_mean_guild_z": inc["guild_z_mean"].mean(),
            "decreased_mean_guild_z": dec["guild_z_mean"].mean(),
            "decreased_minus_increased_guild_z": dec["guild_z_mean"].mean()
            - inc["guild_z_mean"].mean(),
            "decreased_mean_anaerobe_z": dec["anaerobe_z_mean"].mean(),
        }
        if group_column == "condition":
            record.update(paired_guild_enrichment_stats(group))
        records.append(record)
    return pd.DataFrame(records).sort_values(
        ["decreased_mean_guild_z", "profile_count"],
        ascending=[False, False],
    )


def remove_excluded_condition_labels(condition_summary: pd.DataFrame) -> pd.DataFrame:
    """Remove broad study-design labels from disease/injury condition summaries."""
    return condition_summary[
        ~condition_summary["condition"].map(is_excluded_condition_label)
    ].copy()


def is_excluded_condition_label(value: object) -> bool:
    """Return whether a condition label is too broad for disease-context summaries."""
    label = clean_text(value)
    broad_terms = ("exposure", "measurement", "population")
    lower_label = label.lower()
    return label in EXCLUDED_CONDITION_SUMMARY_LABELS or any(
        term in lower_label for term in broad_terms
    )


def select_top_conditions(condition_summary: pd.DataFrame, limit: int) -> pd.DataFrame:
    """Select top condition rows after removing broad excluded labels."""
    condition_summary = remove_excluded_condition_labels(condition_summary)
    top_conditions = (
        condition_summary[
            (condition_summary["decreased_profile_count"] >= 3)
            & (condition_summary["study_count"] >= 2)
            & (condition_summary["decreased_minus_increased_guild_z"] > 0)
            & (condition_summary["decreased_mean_guild_z"] > 0)
        ]
        .sort_values(
            ["decreased_minus_increased_guild_z", "decreased_mean_guild_z"],
            ascending=[False, False],
        )
        .head(limit)
    )
    if top_conditions.empty:
        top_conditions = condition_summary[condition_summary["decreased_profile_count"] >= 3].head(
            limit
        )
    return top_conditions


def paired_smd_label(row: Any) -> str:
    pair_count = int(row.paired_profile_pair_count)
    smd = row.paired_standardized_mean_difference
    p_value = row.paired_t_test_p_value
    p_text = "p=NA" if pd.isna(p_value) else f"p={p_value:.3g}"
    if pd.isna(smd):
        return f"paired d=NA, {p_text}, n={pair_count}"
    return f"paired d={smd:.2f}, {p_text}, n={pair_count}"


def limited_join(values: list[str], limit: int) -> str:
    selected = [value for value in values if clean_text(value) != "unknown"]
    if len(selected) <= limit:
        return "; ".join(selected)
    omitted = len(selected) - limit
    return "; ".join([*selected[:limit], f"+{omitted} more"])


def pluralize(count: int, singular: str, plural: str | None = None) -> str:
    if count == 1:
        return singular
    return plural or f"{singular}s"


def direction_label(inc_dec: object) -> str:
    """Return a display label for compact profile directions."""
    return DIRECTION_LABELS.get(clean_text(inc_dec), clean_text(inc_dec))


def count_zero_enriched_traits_by_direction(profiles: pd.DataFrame) -> dict[str, int]:
    """Count profiles with no positive selected guild-trait z-scores by direction."""
    zero_enriched = profiles["enriched_guild_trait_count"].eq(0)
    counts = zero_enriched.groupby(profiles["inc_dec"]).sum()
    return {direction: int(counts.get(direction, 0)) for direction in DIRECTION_LABELS}


def format_zero_enriched_trait_counts(counts: dict[str, int]) -> str:
    """Return a compact zero-enriched profile-count note."""
    return (
        "Zero enriched selected traits: "
        f"increased {counts.get('inc', 0):,}; decreased {counts.get('dec', 0):,}"
    )


def profile_count_bin_label(profile_count: int) -> str:
    """Return the histogram bin label for an enriched profile count."""
    for lower, upper, label in PROFILE_COUNT_BIN_DEFINITIONS:
        if upper is None and profile_count >= lower:
            return label
        if upper is not None and lower <= profile_count <= upper:
            return label
    return "0"


def trait_profile_count_bin_label(profile_count: int) -> str:
    """Return the histogram bin label for a trait's enriched-profile support."""
    for lower, upper, label in TRAIT_PROFILE_COUNT_BIN_DEFINITIONS:
        if upper is None and profile_count >= lower:
            return label
        if upper is not None and lower <= profile_count <= upper:
            return label
    return "0"


def summarize_condition_enriched_profile_counts(profiles: pd.DataFrame) -> pd.DataFrame:
    """Summarize condition-direction profile counts after zero-enrichment dropouts."""
    selected = profiles.copy()
    selected["condition"] = selected["condition"].map(clean_text)
    selected = selected[~selected["condition"].map(is_excluded_condition_label)].copy()

    rows: list[dict[str, object]] = []
    for (condition, inc_dec), group in selected.groupby(["condition", "inc_dec"], dropna=False):
        enriched_group = group[group["has_enriched_guild_trait"]].copy()
        profile_count = len(group)
        enriched_profile_count = len(enriched_group)
        zero_enriched_profile_count = profile_count - enriched_profile_count
        rows.append(
            {
                "condition": condition,
                "inc_dec": inc_dec,
                "direction_label": direction_label(inc_dec),
                "profile_count": profile_count,
                "enriched_profile_count": enriched_profile_count,
                "zero_enriched_profile_count": zero_enriched_profile_count,
                "enriched_profile_fraction": safe_ratio(
                    float(enriched_profile_count),
                    float(profile_count),
                ),
                "study_count": group["study_id"].nunique(),
                "enriched_study_count": enriched_group["study_id"].nunique(),
                "max_enriched_guild_trait_count": (
                    int(enriched_group["enriched_guild_trait_count"].max())
                    if enriched_profile_count
                    else 0
                ),
            }
        )

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(
        ["enriched_profile_count", "profile_count", "condition", "inc_dec"],
        ascending=[False, False, True, True],
    )


def summarize_body_site_enriched_profile_counts(
    profiles: pd.DataFrame,
    min_profiles: int,
) -> pd.DataFrame:
    """Summarize body-site profile counts after zero-enrichment dropouts."""
    rows: list[dict[str, object]] = []
    for body_site_group, group in profiles.groupby("body_site_group", dropna=False):
        if len(group) < min_profiles:
            continue

        enriched_group = group[group["has_enriched_guild_trait"]].copy()
        inc = group[group["inc_dec"] == "inc"]
        dec = group[group["inc_dec"] == "dec"]
        enriched_inc = enriched_group[enriched_group["inc_dec"] == "inc"]
        enriched_dec = enriched_group[enriched_group["inc_dec"] == "dec"]
        enriched_profile_count = len(enriched_group)
        rows.append(
            {
                "body_site_group": clean_text(body_site_group),
                "profile_count": len(group),
                "enriched_profile_count": enriched_profile_count,
                "zero_enriched_profile_count": len(group) - enriched_profile_count,
                "increased_profile_count": len(inc),
                "decreased_profile_count": len(dec),
                "increased_enriched_profile_count": len(enriched_inc),
                "decreased_enriched_profile_count": len(enriched_dec),
                "increased_zero_enriched_profile_count": len(inc) - len(enriched_inc),
                "decreased_zero_enriched_profile_count": len(dec) - len(enriched_dec),
                "enriched_profile_fraction": safe_ratio(
                    float(enriched_profile_count),
                    float(len(group)),
                ),
                "decreased_enriched_profile_fraction": safe_ratio(
                    float(len(enriched_dec)),
                    float(enriched_profile_count),
                ),
                "mean_enriched_guild_trait_count": enriched_group[
                    "enriched_guild_trait_count"
                ].mean(),
            }
        )

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(
        ["enriched_profile_count", "profile_count", "body_site_group"],
        ascending=[False, False, True],
    )


def summarize_trait_enriched_profile_counts(
    study_matrix: pd.DataFrame,
    profiles: pd.DataFrame,
    study_trait_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize positive trait z-scores within enrichment-filtered profiles."""
    retained_profiles = profiles[profiles["has_enriched_guild_trait"]].copy()
    retained_index = retained_profiles.index.intersection(study_matrix.index)
    retained_profiles = retained_profiles.loc[retained_index]
    if retained_profiles.empty:
        return pd.DataFrame()

    metadata = study_trait_metadata.copy()
    metadata["trait_id"] = metadata["trait_id"].fillna(
        pd.Series(metadata.index, index=metadata.index)
    )
    metadata["trait_name"] = metadata["object_name"].map(clean_text)
    metadata = metadata.drop_duplicates(subset=["trait_name"]).set_index("trait_name", drop=False)

    rows: list[dict[str, object]] = []
    for trait_name in study_matrix.columns:
        trait_values = pd.to_numeric(
            study_matrix.loc[retained_index, trait_name],
            errors="coerce",
        ).fillna(0)
        metadata_row = (
            metadata.loc[trait_name] if trait_name in metadata.index else pd.Series(dtype=object)
        )
        trait_id = str(metadata_row.get("trait_id", ""))
        for inc_dec, group in retained_profiles.groupby("inc_dec", dropna=False):
            direction_values = trait_values.loc[group.index]
            positive_profile_ids = direction_values.index[direction_values > 0]
            positive_group = group.loc[positive_profile_ids]
            trait_enriched_profile_count = len(positive_group)
            retained_profile_count = len(group)
            rows.append(
                {
                    "trait_id": trait_id,
                    "trait_name": trait_name,
                    "trait_category": metadata_row.get("display_category", ""),
                    "is_guild_trait": trait_id in GUILD_TRAITS,
                    "inc_dec": inc_dec,
                    "direction_label": direction_label(inc_dec),
                    "retained_profile_count": retained_profile_count,
                    "trait_enriched_profile_count": trait_enriched_profile_count,
                    "trait_zero_profile_count": (
                        retained_profile_count - trait_enriched_profile_count
                    ),
                    "trait_enriched_profile_fraction": safe_ratio(
                        float(trait_enriched_profile_count),
                        float(retained_profile_count),
                    ),
                    "retained_study_count": group["study_id"].nunique(),
                    "trait_enriched_study_count": positive_group["study_id"].nunique(),
                    "retained_condition_count": group["condition"].map(clean_text).nunique(),
                    "trait_enriched_condition_count": (
                        positive_group["condition"].map(clean_text).nunique()
                    ),
                }
            )

    return pd.DataFrame(rows).sort_values(
        [
            "trait_enriched_profile_count",
            "trait_enriched_study_count",
            "trait_name",
            "inc_dec",
        ],
        ascending=[False, False, True, True],
    )


def build_kg_microbe_taxon_universe() -> pd.DataFrame:
    """Define the KG-Microbe-wide candidate taxon universe from trait criteria."""
    trait_ids = {
        *KG_MICROBE_TAXON_UNIVERSE_TRAITS,
        *KG_MICROBE_TAXON_UNIVERSE_ASSOCIATED_TRAITS,
    }
    matrix = read_tsv(
        REPO_ROOT / "data/derived/taxa_trait_matrix.tsv",
        usecols=lambda column: column in {"taxon_id", "taxon_name", "taxonomic_rank", *trait_ids},
    )
    for trait_id in trait_ids:
        matrix[trait_id] = pd.to_numeric(matrix[trait_id], errors="coerce").fillna(0)

    core_columns = list(KG_MICROBE_TAXON_UNIVERSE_TRAITS)
    associated_columns = list(KG_MICROBE_TAXON_UNIVERSE_ASSOCIATED_TRAITS)
    matrix["kg_microbe_universe_core_trait_count"] = (matrix[core_columns] > 0).sum(axis=1)
    if associated_columns:
        matrix["kg_microbe_universe_associated_trait_count"] = (matrix[associated_columns] > 0).sum(
            axis=1
        )
    else:
        matrix["kg_microbe_universe_associated_trait_count"] = 0
    universe = matrix[
        matrix["kg_microbe_universe_core_trait_count"] >= MIN_KG_MICROBE_UNIVERSE_CORE_TRAITS
    ].copy()
    universe["kg_microbe_universe_core_traits"] = [
        "; ".join(
            trait_name
            for trait_id, trait_name in KG_MICROBE_TAXON_UNIVERSE_TRAITS.items()
            if row[trait_id] > 0
        )
        for _, row in universe.iterrows()
    ]
    if KG_MICROBE_TAXON_UNIVERSE_ASSOCIATED_TRAITS:
        universe["kg_microbe_universe_associated_traits"] = [
            "; ".join(
                trait_name
                for trait_id, trait_name in KG_MICROBE_TAXON_UNIVERSE_ASSOCIATED_TRAITS.items()
                if row[trait_id] > 0
            )
            for _, row in universe.iterrows()
        ]
    else:
        universe["kg_microbe_universe_associated_traits"] = ""
    return universe.sort_values(
        [
            "kg_microbe_universe_core_trait_count",
            "kg_microbe_universe_associated_trait_count",
            "taxon_name",
        ],
        ascending=[False, False, True],
    )


def build_guild_trait_selection(study_trait_metadata: pd.DataFrame) -> pd.DataFrame:
    """Document why each selected guild trait is in the score."""
    trait_ids = set(GUILD_TRAITS)
    taxa_trait_matrix = read_tsv(
        REPO_ROOT / "data/derived/taxa_trait_matrix.tsv",
        usecols=lambda column: column in {"taxon_id", *trait_ids},
    )
    for trait_id in trait_ids:
        taxa_trait_matrix[trait_id] = pd.to_numeric(
            taxa_trait_matrix[trait_id], errors="coerce"
        ).fillna(0)

    metadata = study_trait_metadata.copy()
    metadata["trait_id"] = metadata["trait_id"].fillna(
        pd.Series(metadata.index, index=metadata.index)
    )
    metadata = metadata.drop_duplicates(subset=["trait_id"]).set_index("trait_id", drop=False)

    records: list[dict[str, object]] = []
    for trait_id, trait_name in GUILD_TRAITS.items():
        metadata_row = (
            metadata.loc[trait_id] if trait_id in metadata.index else pd.Series(dtype=object)
        )
        records.append(
            {
                "guild": "anaerobe v2",
                "focus_direction": "decreased",
                "comparator_direction": "increased",
                "candidate_universe_role": "core oxygen-tolerance trait",
                "trait_id": trait_id,
                "trait_name": trait_name,
                "profile_matrix_object_name": metadata_row.get("object_name", trait_name),
                "profile_support": metadata_row.get("profile_support", pd.NA),
                "total_signal": metadata_row.get("total_signal", pd.NA),
                "composition_variance": metadata_row.get("composition_variance", pd.NA),
                "kg_microbe_taxa_with_trait_count": int((taxa_trait_matrix[trait_id] > 0).sum()),
                "selection_justification": GUILD_TRAIT_SELECTION_JUSTIFICATIONS[trait_id],
            }
        )
    return pd.DataFrame(records)


def load_profile_taxa_occurrences(scored_profiles: pd.DataFrame) -> pd.DataFrame:
    """Load all taxon/profile occurrences and join scored profile metadata."""
    taxa = read_tsv(
        REPO_ROOT / "data/derived/bugsigdb_taxa.tsv",
        usecols=[
            "study_id",
            "experiment_id",
            "signature_id",
            "direction",
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
        ],
    )
    taxa["profile_id"] = taxa.apply(profile_id_from_taxa_row, axis=1)
    taxa = taxa.drop_duplicates(subset=["taxon_id", "profile_id"])
    profile_columns = [
        "study_id",
        "experiment_id",
        "inc_dec",
        "direction_label",
        "title",
        "publication_year",
        "condition",
        "body_site",
        "body_site_group",
        "condition_group",
        "location",
        "host_species",
        "guild_z_mean",
        "anaerobe_z_mean",
    ]
    profiles = scored_profiles[profile_columns].copy()
    return taxa.merge(
        profiles,
        left_on="profile_id",
        right_index=True,
        how="inner",
        suffixes=("_taxa", ""),
    )


def trait_metadata_for_inference(
    study_matrix: pd.DataFrame,
    study_trait_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Return trait metadata rows that can be matched across profile and taxa matrices."""
    taxa_trait_path = REPO_ROOT / "data/derived/taxa_trait_matrix.tsv"
    with taxa_trait_path.open(encoding="utf-8") as handle:
        taxa_trait_columns = set(handle.readline().rstrip("\n").split("\t"))

    metadata = study_trait_metadata.copy()
    metadata["trait_id"] = metadata["trait_id"].fillna(
        pd.Series(metadata.index, index=metadata.index)
    )
    metadata = metadata[
        metadata["trait_id"].isin(taxa_trait_columns)
        & metadata["object_name"].isin(study_matrix.columns)
    ].copy()
    metadata = metadata.drop_duplicates(subset=["trait_id"])
    metadata["trait_name"] = metadata["object_name"].map(clean_text)
    metadata = metadata.drop_duplicates(subset=["trait_name"])
    return metadata


def load_taxa_trait_subset(taxon_ids: set[str], trait_ids: set[str]) -> pd.DataFrame:
    """Load taxa-trait matrix columns needed for inferred guild evidence."""
    wanted_columns = {"taxon_id", "taxon_name", "taxonomic_rank", *trait_ids}
    matrix = read_tsv(
        REPO_ROOT / "data/derived/taxa_trait_matrix.tsv",
        usecols=lambda column: column in wanted_columns,
    )
    matrix = matrix[matrix["taxon_id"].isin(taxon_ids)].copy()
    for trait_id in sorted(trait_ids.intersection(matrix.columns)):
        matrix[trait_id] = pd.to_numeric(matrix[trait_id], errors="coerce").fillna(0)
    return matrix


def profile_group_mask(
    profiles: pd.DataFrame,
    group_row: pd.Series,
    group_columns: list[str],
) -> pd.Series:
    mask = pd.Series(True, index=profiles.index)
    for column in group_columns:
        mask &= profiles[column] == group_row[column]
    return mask


def select_enriched_traits_for_group(
    study_matrix: pd.DataFrame,
    trait_metadata: pd.DataFrame,
    focus_profiles: pd.DataFrame,
    comparator_profiles: pd.DataFrame,
) -> pd.DataFrame:
    """Select profile traits enriched in the focus direction for one group."""
    trait_rows = trait_metadata.set_index("trait_name", drop=False)
    trait_names = [trait_name for trait_name in trait_rows.index if trait_name in study_matrix]
    if not trait_names:
        return pd.DataFrame()

    focus_means = study_matrix.loc[focus_profiles.index, trait_names].mean()
    if comparator_profiles.empty:
        comparator_means = pd.Series(math.nan, index=trait_names)
    else:
        comparator_means = study_matrix.loc[comparator_profiles.index, trait_names].mean()

    records: list[dict[str, object]] = []
    for trait_name in trait_names:
        focus_mean = focus_means[trait_name]
        comparator_mean = comparator_means[trait_name]
        delta = focus_mean - comparator_mean
        if pd.isna(focus_mean) or focus_mean <= 0:
            continue
        if not pd.isna(delta) and delta <= 0:
            continue
        metadata_row = trait_rows.loc[trait_name]
        selection_score = focus_mean if pd.isna(delta) else delta * math.sqrt(len(focus_profiles))
        records.append(
            {
                "trait_id": metadata_row["trait_id"],
                "trait_name": trait_name,
                "trait_category": metadata_row["display_category"],
                "focus_mean_z": focus_mean,
                "comparator_mean_z": comparator_mean,
                "focus_minus_comparator_z": delta,
                "selection_score": selection_score,
            }
        )

    if not records:
        return pd.DataFrame()
    return (
        pd.DataFrame(records)
        .sort_values(
            ["focus_minus_comparator_z", "focus_mean_z", "selection_score"],
            ascending=[False, False, False],
        )
        .head(INFERRED_GUILD_TOP_TRAITS_PER_GROUP)
    )


def summarize_inferred_guild_evidence(
    study_matrix: pd.DataFrame,
    profiles: pd.DataFrame,
    study_trait_metadata: pd.DataFrame,
    profile_taxa_occurrences: pd.DataFrame,
    group_rows: pd.DataFrame,
    group_columns: list[str],
    *,
    focus_direction: str,
    taxon_universe_ids: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Infer condition-specific guild evidence from enriched traits and observed taxa."""
    trait_metadata = trait_metadata_for_inference(study_matrix, study_trait_metadata)
    trait_metadata = trait_metadata[trait_metadata["trait_id"].isin(GUILD_TRAITS)].copy()
    if trait_metadata.empty or group_rows.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    selected_profile_ids: set[str] = set()
    group_profile_rows: list[tuple[pd.Series, pd.DataFrame, pd.DataFrame]] = []
    for _, group_row in group_rows.iterrows():
        group_profiles = profiles[profile_group_mask(profiles, group_row, group_columns)].copy()
        focus_profiles = group_profiles[group_profiles["direction_label"] == focus_direction]
        comparator_profiles = group_profiles[group_profiles["direction_label"] != focus_direction]
        if focus_profiles.empty:
            continue
        selected_profile_ids.update(map(str, focus_profiles.index))
        group_profile_rows.append((group_row, focus_profiles, comparator_profiles))

    observed_taxa = profile_taxa_occurrences[
        profile_taxa_occurrences["profile_id"].isin(selected_profile_ids)
    ].copy()
    if taxon_universe_ids is not None:
        observed_taxa = observed_taxa[observed_taxa["taxon_id"].isin(taxon_universe_ids)].copy()
    if observed_taxa.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    taxa_trait_matrix = load_taxa_trait_subset(
        set(observed_taxa["taxon_id"]),
        set(trait_metadata["trait_id"]),
    )
    taxa_trait_by_id = taxa_trait_matrix.set_index("taxon_id", drop=False)

    group_summary_rows: list[dict[str, object]] = []
    trait_summary_rows: list[dict[str, object]] = []
    taxon_summary_rows: list[dict[str, object]] = []

    for group_row, focus_profiles, comparator_profiles in group_profile_rows:
        enriched_traits = select_enriched_traits_for_group(
            study_matrix,
            trait_metadata,
            focus_profiles,
            comparator_profiles,
        )
        if enriched_traits.empty:
            continue

        focus_profile_ids = set(map(str, focus_profiles.index))
        focus_occurrences = observed_taxa[
            observed_taxa["profile_id"].isin(focus_profile_ids)
        ].copy()
        observed_taxon_count = focus_occurrences["taxon_id"].nunique()
        total_focus_taxon_count = profile_taxa_occurrences[
            profile_taxa_occurrences["profile_id"].isin(focus_profile_ids)
        ]["taxon_id"].nunique()
        group_values = {column: group_row[column] for column in group_columns}

        for trait_row in enriched_traits.itertuples(index=False):
            trait_id = str(trait_row.trait_id)
            if trait_id not in taxa_trait_by_id:
                continue
            carrier_ids = set(taxa_trait_by_id.index[pd.to_numeric(taxa_trait_by_id[trait_id]) > 0])
            carrier_occurrences = focus_occurrences[
                focus_occurrences["taxon_id"].isin(carrier_ids)
            ].copy()
            carrier_occurrences = carrier_occurrences.drop_duplicates(
                subset=["taxon_id", "profile_id"]
            )
            if carrier_occurrences.empty:
                continue

            taxa_for_trait = (
                carrier_occurrences.groupby(
                    ["taxon_id", "taxon_name", "taxonomic_rank"],
                    dropna=False,
                )
                .agg(
                    taxon_profile_count=("profile_id", "nunique"),
                    study_count=("study_id", "nunique"),
                )
                .reset_index()
            )
            taxa_for_trait["trait_annotation_count"] = taxa_for_trait["taxon_id"].map(
                taxa_trait_by_id[trait_id]
            )
            taxa_for_trait = taxa_for_trait.sort_values(
                ["taxon_profile_count", "study_count", "taxon_name"],
                ascending=[False, False, True],
            )
            top_taxa = [
                f"{row.taxon_name} ({int(row.taxon_profile_count)})"
                for row in taxa_for_trait.head(INFERRED_GUILD_TOP_TAXA_PER_TRAIT).itertuples(
                    index=False
                )
            ]

            trait_summary_rows.append(
                {
                    **group_values,
                    "focus_direction": focus_direction,
                    "focus_profile_count": len(focus_profiles),
                    "focus_study_count": focus_profiles["study_id"].nunique(),
                    "comparator_profile_count": len(comparator_profiles),
                    "comparator_study_count": comparator_profiles["study_id"].nunique(),
                    "total_observed_taxon_count": total_focus_taxon_count,
                    "trait_id": trait_id,
                    "trait_name": trait_row.trait_name,
                    "trait_category": trait_row.trait_category,
                    "focus_mean_z": trait_row.focus_mean_z,
                    "comparator_mean_z": trait_row.comparator_mean_z,
                    "focus_minus_comparator_z": trait_row.focus_minus_comparator_z,
                    "observed_taxon_count": observed_taxon_count,
                    "carrier_taxon_count": taxa_for_trait["taxon_id"].nunique(),
                    "carrier_taxon_fraction": safe_ratio(
                        float(taxa_for_trait["taxon_id"].nunique()),
                        float(observed_taxon_count),
                    ),
                    "carrier_taxon_profile_count": carrier_occurrences["profile_id"].nunique(),
                    "carrier_study_count": carrier_occurrences["study_id"].nunique(),
                    "top_carrier_taxa": limited_join(top_taxa, INFERRED_GUILD_TOP_TAXA_PER_TRAIT),
                }
            )

            for taxon_row in taxa_for_trait.itertuples(index=False):
                taxon_summary_rows.append(
                    {
                        **group_values,
                        "focus_direction": focus_direction,
                        "trait_id": trait_id,
                        "trait_name": trait_row.trait_name,
                        "trait_category": trait_row.trait_category,
                        "taxon_id": taxon_row.taxon_id,
                        "taxon_name": taxon_row.taxon_name,
                        "taxonomic_rank": taxon_row.taxonomic_rank,
                        "taxon_profile_count": taxon_row.taxon_profile_count,
                        "study_count": taxon_row.study_count,
                        "trait_annotation_count": taxon_row.trait_annotation_count,
                    }
                )

    trait_summary = pd.DataFrame(trait_summary_rows)
    taxon_summary = pd.DataFrame(taxon_summary_rows)
    if trait_summary.empty or taxon_summary.empty:
        return pd.DataFrame(), trait_summary, taxon_summary

    for group_values, group_traits in trait_summary.groupby(group_columns, dropna=False):
        if not isinstance(group_values, tuple):
            group_values = (group_values,)
        group_dict = dict(zip(group_columns, group_values, strict=True))
        selected_taxa = taxon_summary.copy()
        for column, value in group_dict.items():
            selected_taxa = selected_taxa[selected_taxa[column] == value]

        taxa_rank = (
            selected_taxa.groupby(["taxon_id", "taxon_name", "taxonomic_rank"], dropna=False)
            .agg(
                supported_enriched_trait_count=("trait_id", "nunique"),
                max_taxon_profile_count=("taxon_profile_count", "max"),
                max_study_count=("study_count", "max"),
            )
            .reset_index()
            .sort_values(
                ["supported_enriched_trait_count", "max_taxon_profile_count", "max_study_count"],
                ascending=[False, False, False],
            )
        )
        top_trait_text = [
            f"{row.trait_name} ({row.focus_minus_comparator_z:+.2f})"
            for row in group_traits.sort_values(
                ["focus_minus_comparator_z", "focus_mean_z"],
                ascending=[False, False],
            )
            .head(5)
            .itertuples(index=False)
        ]
        top_taxa_text = [
            f"{row.taxon_name} ("
            f"{int(row.supported_enriched_trait_count)} "
            f"{pluralize(int(row.supported_enriched_trait_count), 'trait')}, "
            f"{int(row.max_taxon_profile_count)} "
            f"{pluralize(int(row.max_taxon_profile_count), 'profile')})"
            for row in taxa_rank.head(INFERRED_GUILD_TOP_TAXA_PER_GROUP).itertuples(index=False)
        ]
        first_trait = group_traits.iloc[0]
        group_summary_rows.append(
            {
                **group_dict,
                "focus_direction": focus_direction,
                "focus_profile_count": int(first_trait["focus_profile_count"]),
                "focus_study_count": int(first_trait["focus_study_count"]),
                "comparator_profile_count": int(first_trait["comparator_profile_count"]),
                "comparator_study_count": int(first_trait["comparator_study_count"]),
                "total_observed_taxon_count": int(first_trait["total_observed_taxon_count"]),
                "kg_microbe_universe_taxon_count": int(group_traits["observed_taxon_count"].max()),
                "inferred_enriched_trait_count": group_traits["trait_id"].nunique(),
                "supporting_taxon_count": selected_taxa["taxon_id"].nunique(),
                "top_enriched_traits": limited_join(top_trait_text, 5),
                "top_supporting_taxa": limited_join(
                    top_taxa_text, INFERRED_GUILD_TOP_TAXA_PER_GROUP
                ),
            }
        )

    group_summary = pd.DataFrame(group_summary_rows)
    return group_summary, trait_summary, taxon_summary


def load_taxa_direction_summary() -> pd.DataFrame:
    taxa_metadata = read_tsv(VIS_DIR / "taxa_guild_discovery_taxa_metadata.tsv", index_col=0)
    summary = taxa_metadata.loc[taxa_metadata.index.intersection(TARGET_TAXA.keys())].copy()
    summary["target_taxon_name"] = summary.index.map(TARGET_TAXA)
    summary["decreased_fraction"] = [
        safe_ratio(float(row["decreased"]), float(row["bugsigdb_profile_count"]))
        for _, row in summary.iterrows()
    ]
    summary["depletion_bias"] = [
        safe_ratio(
            float(row["decreased"]) - float(row["increased"]),
            float(row["bugsigdb_profile_count"]),
        )
        for _, row in summary.iterrows()
    ]
    summary["decreased_to_increased_ratio"] = [
        safe_ratio(float(row["decreased"]), float(row["increased"]))
        for _, row in summary.iterrows()
    ]
    return summary[
        [
            "target_taxon_name",
            "taxon_name",
            "taxonomic_rank",
            "bugsigdb_study_count",
            "bugsigdb_profile_count",
            "decreased",
            "increased",
            "dominant_direction",
            "decreased_fraction",
            "depletion_bias",
            "decreased_to_increased_ratio",
            "trait_count",
        ]
    ].sort_values("depletion_bias", ascending=False)


def load_target_taxa_trait_counts() -> pd.DataFrame:
    usecols = ["taxon_id", "taxon_name", "taxonomic_rank", *GUILD_TRAITS.keys()]
    matrix = read_tsv(REPO_ROOT / "data/derived/taxa_trait_matrix.tsv", usecols=usecols)
    matrix = matrix[matrix["taxon_id"].isin(TARGET_TAXA)].copy()
    taxa_metadata = read_tsv(VIS_DIR / "taxa_guild_discovery_taxa_metadata.tsv", index_col=0)
    matrix["target_taxon_name"] = matrix["taxon_id"].map(TARGET_TAXA)
    matrix["heatmap_taxon_name"] = matrix["taxon_id"].map(taxa_metadata["taxon_name"])
    matrix["heatmap_taxonomic_rank"] = matrix["taxon_id"].map(taxa_metadata["taxonomic_rank"])

    for trait_id, trait_name in GUILD_TRAITS.items():
        matrix[f"{trait_name}_count"] = pd.to_numeric(matrix[trait_id], errors="coerce").fillna(0)
        matrix[f"{trait_name}_present"] = matrix[f"{trait_name}_count"] > 0

    present_columns = [f"{trait_name}_present" for trait_name in GUILD_TRAITS.values()]
    matrix["guild_trait_presence_count"] = matrix[present_columns].sum(axis=1)
    output_columns = [
        "taxon_id",
        "target_taxon_name",
        "heatmap_taxonomic_rank",
        "guild_trait_presence_count",
    ]
    for trait_name in GUILD_TRAITS.values():
        output_columns.extend([f"{trait_name}_present", f"{trait_name}_count"])
    output = matrix[output_columns].rename(columns={"heatmap_taxonomic_rank": "taxonomic_rank"})
    return output.sort_values(
        ["guild_trait_presence_count", "target_taxon_name"],
        ascending=[False, True],
    )


def load_target_taxa_occurrences(
    scored_profiles: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return target-taxon occurrences in human-host profiles and those set aside (#587)."""
    taxa = read_tsv(
        REPO_ROOT / "data/derived/bugsigdb_taxa.tsv",
        usecols=[
            "study_id",
            "experiment_id",
            "signature_id",
            "direction",
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
        ],
    )
    taxa = taxa[taxa["taxon_id"].isin(TARGET_TAXA)].copy()
    taxa["profile_id"] = taxa.apply(profile_id_from_taxa_row, axis=1)
    taxa = taxa.drop_duplicates(subset=["taxon_id", "profile_id"])
    taxa["target_taxon_name"] = taxa["taxon_id"].map(TARGET_TAXA)
    experiment_hosts = read_tsv(
        REPO_ROOT / "data/derived/experiments.tsv",
        usecols=["study_id", "experiment_id", "host_species"],
    )
    profile_columns = [
        "study_id",
        "experiment_id",
        "inc_dec",
        "direction_label",
        "title",
        "publication_year",
        "condition",
        "body_site",
        "body_site_group",
        "condition_group",
        "location",
        "host_species",
        "guild_z_mean",
        "anaerobe_z_mean",
    ]
    return join_human_host_profiles(taxa, scored_profiles[profile_columns], experiment_hosts)


def summarize_taxa_conditions(occurrences: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        occurrences.groupby(["target_taxon_name", "condition", "body_site_group", "direction"])[
            ["profile_id", "study_id"]
        ]
        .nunique()
        .rename(columns={"profile_id": "profile_count", "study_id": "study_count"})
        .reset_index()
    )
    if grouped.empty:
        return grouped

    profile_pivot = grouped.pivot_table(
        index=["target_taxon_name", "condition", "body_site_group"],
        columns="direction",
        values="profile_count",
        aggfunc="sum",
        fill_value=0,
    )
    study_pivot = grouped.pivot_table(
        index=["target_taxon_name", "condition", "body_site_group"],
        columns="direction",
        values="study_count",
        aggfunc="sum",
        fill_value=0,
    )
    for direction in ["decreased", "increased"]:
        if direction not in profile_pivot.columns:
            profile_pivot[direction] = 0
        if direction not in study_pivot.columns:
            study_pivot[direction] = 0

    output = profile_pivot.rename(
        columns={
            "decreased": "decreased_profile_count",
            "increased": "increased_profile_count",
        }
    ).join(
        study_pivot.rename(
            columns={
                "decreased": "decreased_study_count",
                "increased": "increased_study_count",
            }
        )
    )
    output = output.reset_index()
    output["total_profile_count"] = (
        output["decreased_profile_count"] + output["increased_profile_count"]
    )
    output["depletion_bias"] = (
        output["decreased_profile_count"] - output["increased_profile_count"]
    ) / output["total_profile_count"].replace(0, pd.NA)
    return output.sort_values(
        ["depletion_bias", "decreased_study_count", "total_profile_count"],
        ascending=[False, False, False],
    )


def count_target_taxa_conditions(
    occurrences: pd.DataFrame,
    direction: str,
) -> pd.Series:
    """Count condition labels for target taxon-profile occurrences."""
    selected = occurrences[occurrences["direction"] == direction].copy()
    selected["condition"] = selected["condition"].map(clean_text)
    selected = selected[~selected["condition"].map(is_excluded_condition_label)].copy()
    return selected["condition"].value_counts()


def save_condition_breakdown_bar_chart(
    counts: pd.Series,
    *,
    output_path: Path,
    title_prefix: str,
    color: str,
    specificity_note: str,
    note: str,
) -> None:
    """Save a horizontal bar chart for condition counts, with `note` as its footer."""
    if counts.empty:
        return

    bar_counts = counts.head(20)
    fig, ax = plt.subplots(figsize=(10, max(5, len(bar_counts) * 0.35)))
    ax.barh(bar_counts.index, bar_counts.to_numpy(), color=color)
    ax.invert_yaxis()
    ax.set_xlabel("Target taxon-profile occurrences")
    ax.set_title(f"{title_prefix}: top conditions", pad=28)
    add_specificity_note(ax, specificity_note)
    for index, value in enumerate(bar_counts.to_numpy()):
        ax.text(value, index, f" {int(value)}", va="center", fontsize=8)
    fig.text(0.01, 0.01, note, fontsize=8, color="#475569")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def dominant_body_site_label(group: pd.DataFrame) -> str:
    """Return a concise dominant body-site label for a condition group."""
    body_site_counts = group["body_site_group"].map(display_body_site).value_counts()
    if body_site_counts.empty:
        return "site unknown"
    dominant_site = str(body_site_counts.index[0])
    dominant_count = int(body_site_counts.iloc[0])
    total_count = int(body_site_counts.sum())
    if dominant_count == total_count:
        return dominant_site
    return f"{dominant_site} {dominant_count}/{total_count}"


def guild_delta_label(condition: str, condition_summary: pd.DataFrame) -> str:
    """Return the condition's trait-enrichment contrast label when available."""
    if condition_summary.empty:
        return ""
    matches = condition_summary[condition_summary["condition"] == condition]
    if matches.empty:
        return ""
    value = matches.iloc[0]["decreased_minus_increased_guild_z"]
    if pd.isna(value):
        return ""
    if abs(value) < 0.005:
        value = 0
    return f"; trait dz={value:.2f}"


def save_condition_target_taxa_stacked_bar_chart(
    occurrences: pd.DataFrame,
    condition_summary: pd.DataFrame,
    *,
    output_path: Path,
    note: str,
    limit: int = 12,
) -> None:
    """Save a stacked bar chart of decreased-profile focal-taxon recurrence by condition.

    `note` follows the occurrence-count caveat in the footer.
    """
    selected = occurrences[occurrences["direction"] == "decreased"].copy()
    selected["condition"] = selected["condition"].map(clean_text)
    selected = selected[~selected["condition"].map(is_excluded_condition_label)].copy()
    if selected.empty:
        return

    counts = selected.pivot_table(
        index="condition",
        columns="target_taxon_name",
        values="profile_id",
        aggfunc="count",
        fill_value=0,
    )
    condition_totals = counts.sum(axis=1).sort_values(ascending=False)
    top_conditions = condition_totals.head(limit).index
    counts = counts.loc[top_conditions]
    taxon_order = counts.sum(axis=0).sort_values(ascending=False)
    taxon_order = [taxon for taxon in taxon_order.index if taxon_order[taxon] > 0]

    selected_top = selected[selected["condition"].isin(top_conditions)]
    study_counts = selected_top.groupby("condition")["study_id"].nunique()
    body_site_labels = {
        condition: dominant_body_site_label(group)
        for condition, group in selected_top.groupby("condition")
    }

    fig_height = max(5.5, len(counts) * 0.45)
    fig, ax = plt.subplots(figsize=(13, fig_height))
    y_positions = list(range(len(counts)))
    left = pd.Series(0, index=counts.index, dtype=float)
    for taxon in taxon_order:
        values = counts[taxon]
        ax.barh(
            y_positions,
            values,
            left=left,
            color=TARGET_TAXON_COLORS.get(taxon, "#475569"),
            label=TARGET_TAXON_LEGEND_LABELS.get(taxon, taxon),
        )
        left = left + values

    totals = counts.sum(axis=1)
    x_max = float(totals.max())
    ax.set_xlim(0, x_max * 1.82)
    x_pad = max(x_max * 0.025, 1)
    for y_pos, condition in enumerate(counts.index):
        total = int(totals.loc[condition])
        study_count = int(study_counts.get(condition, 0))
        body_site = body_site_labels.get(condition, "site unknown")
        label = (
            f"{total} occ.; {study_count} "
            f"{pluralize(study_count, 'study', 'studies')}; {body_site}"
            f"{guild_delta_label(condition, condition_summary)}"
        )
        ax.text(total + x_pad, y_pos, label, va="center", fontsize=8, color="#334155")

    ax.set_yticks(y_positions)
    ax.set_yticklabels(counts.index)
    ax.invert_yaxis()
    ax.set_xlabel("Decreased-profile target-taxon occurrences")
    ax.set_title(
        "Decreased-Profile Anaerobe Target-Taxon Recurrence by Condition",
        pad=28,
    )
    add_specificity_note(ax, NON_GUILD_SPECIFIC_TARGET_TAXA_NOTE)
    ax.legend(
        title="Focal taxon",
        loc="lower right",
        bbox_to_anchor=(0.98, 0.03),
        ncol=1,
        frameon=False,
        fontsize=8,
        title_fontsize=8,
    )
    fig.text(
        0.01,
        0.01,
        "Counts are target-taxon/profile occurrences; a profile can contribute "
        f"once per focal taxon. {note}",
        fontsize=8,
        color="#475569",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def save_trait_enriched_profile_histogram(
    trait_enriched_profile_counts: pd.DataFrame,
    *,
    output_path: Path,
) -> None:
    """Save a histogram of trait counts by enriched profile support."""
    if trait_enriched_profile_counts.empty:
        return

    selected = trait_enriched_profile_counts[
        trait_enriched_profile_counts["trait_enriched_profile_count"] > 0
    ].copy()
    if selected.empty:
        return

    bin_labels = [label for _, _, label in TRAIT_PROFILE_COUNT_BIN_DEFINITIONS]
    selected["profile_count_bin"] = selected["trait_enriched_profile_count"].map(
        lambda count: trait_profile_count_bin_label(int(count))
    )
    histogram = (
        selected.groupby(["profile_count_bin", "inc_dec"], dropna=False)["trait_name"]
        .nunique()
        .unstack(fill_value=0)
        .reindex(bin_labels, fill_value=0)
    )

    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    x_positions = list(range(len(bin_labels)))
    bar_width = 0.38
    offsets = {"inc": -bar_width / 2, "dec": bar_width / 2}

    for direction in ["inc", "dec"]:
        values = (
            histogram[direction].to_numpy()
            if direction in histogram.columns
            else [0] * len(bin_labels)
        )
        ax.bar(
            [x + offsets[direction] for x in x_positions],
            values,
            width=bar_width,
            color=DIRECTION_COLORS[direction],
            label=DIRECTION_LABELS[direction],
        )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(bin_labels, rotation=20, ha="right")
    ax.set_xlabel("Positive-z retained profiles per trait")
    ax.set_ylabel("Number of traits")
    ax.set_title(
        "Trait Enrichment Histogram Within Anaerobe v2 Retained Profiles",
        pad=28,
    )
    add_specificity_note(ax, GUILD_SPECIFIC_NOTE)
    ax.legend(title="BugSigDB profile direction", frameon=False)
    fig.text(
        0.01,
        0.01,
        "Universe: profiles retained by the anaerobe v2 enrichment filter. "
        "Bars count selected heatmap traits by the number of retained profiles "
        "where that trait has a positive z-score.",
        fontsize=8,
        color="#475569",
    )
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def save_condition_enriched_profile_histogram(
    condition_enriched_profile_counts: pd.DataFrame,
    *,
    output_path: Path,
) -> None:
    """Save a histogram of condition counts by enriched profile count."""
    if condition_enriched_profile_counts.empty:
        return

    selected = condition_enriched_profile_counts[
        condition_enriched_profile_counts["enriched_profile_count"] > 0
    ].copy()
    if selected.empty:
        return

    bin_labels = [label for _, _, label in PROFILE_COUNT_BIN_DEFINITIONS]
    selected["profile_count_bin"] = selected["enriched_profile_count"].map(
        lambda count: profile_count_bin_label(int(count))
    )
    histogram = (
        selected.groupby(["profile_count_bin", "inc_dec"], dropna=False)["condition"]
        .nunique()
        .unstack(fill_value=0)
        .reindex(bin_labels, fill_value=0)
    )
    zero_enriched_trait_counts = {
        direction: int(
            condition_enriched_profile_counts.loc[
                condition_enriched_profile_counts["inc_dec"] == direction,
                "zero_enriched_profile_count",
            ].sum()
        )
        for direction in DIRECTION_LABELS
    }

    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    x_positions = list(range(len(bin_labels)))
    bar_width = 0.38
    offsets = {"inc": -bar_width / 2, "dec": bar_width / 2}

    for direction in ["inc", "dec"]:
        values = (
            histogram[direction].to_numpy()
            if direction in histogram.columns
            else [0] * len(bin_labels)
        )
        ax.bar(
            [x + offsets[direction] for x in x_positions],
            values,
            width=bar_width,
            color=DIRECTION_COLORS[direction],
            label=(
                f"{DIRECTION_LABELS[direction]} "
                f"(zero enriched traits: {zero_enriched_trait_counts.get(direction, 0):,})"
            ),
        )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(bin_labels)
    ax.set_xlabel("Enriched profiles per condition")
    ax.set_ylabel("Number of conditions")
    ax.set_title(
        "Condition Histogram After Anaerobe v2 Trait-Enrichment Filtering",
        pad=28,
    )
    add_specificity_note(ax, GUILD_SPECIFIC_NOTE)
    ax.legend(title="BugSigDB profile direction", frameon=False)
    fig.text(
        0.01,
        0.01,
        "A profile is retained when at least one selected anaerobe v2 trait has "
        "a positive z-score; zero-enriched profiles are excluded and counted "
        "in the legend.",
        fontsize=8,
        color="#475569",
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def save_anatomical_enriched_profile_diagram(
    body_site_enriched_profile_summary: pd.DataFrame,
    zero_enriched_trait_counts: dict[str, int],
    *,
    output_path: Path,
    min_enriched_profiles: int = 10,
    max_sites: int = 12,
) -> None:
    """Save an anatomical-site diagram for enrichment-filtered profile counts."""
    if body_site_enriched_profile_summary.empty:
        return

    plot_sites = (
        body_site_enriched_profile_summary[
            ~body_site_enriched_profile_summary["body_site_group"].map(is_other_body_site)
            & (
                body_site_enriched_profile_summary["enriched_profile_count"]
                >= min_enriched_profiles
            )
        ]
        .head(max_sites)
        .copy()
    )
    if plot_sites.empty:
        return

    numeric_columns = [
        "enriched_profile_count",
        "increased_enriched_profile_count",
        "decreased_enriched_profile_count",
    ]
    for column in numeric_columns:
        plot_sites[column] = pd.to_numeric(plot_sites[column], errors="coerce").fillna(0)

    points = [anatomical_site_base(site) for site in plot_sites["body_site_group"]]
    plot_sites["anchor_x"] = [point[0] for point in points]
    plot_sites["anchor_y"] = [point[1] for point in points]
    plot_sites = plot_sites.sort_values(
        ["anchor_y", "enriched_profile_count"],
        ascending=[False, False],
    )

    max_count = max(float(plot_sites["enriched_profile_count"].max()), 1.0)
    fig_height = max(7.2, 1.0 + len(plot_sites) * 0.82)
    fig, ax = plt.subplots(figsize=(11, fig_height))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    add_anatomical_outline(ax)

    label_ys = [0.84]
    if len(plot_sites) > 1:
        step = 0.60 / (len(plot_sites) - 1)
        label_ys = [0.84 - index * step for index in range(len(plot_sites))]

    for label_y, row in zip(label_ys, plot_sites.itertuples(index=False), strict=False):
        enriched_count = int(row.enriched_profile_count)
        increased_count = int(row.increased_enriched_profile_count)
        decreased_count = int(row.decreased_enriched_profile_count)
        size = 95 + 900 * math.sqrt(enriched_count / max_count)
        ax.scatter(
            [row.anchor_x],
            [row.anchor_y],
            s=size,
            c=["#e2e8f0"],
            edgecolors="#111827",
            linewidths=0.8,
            zorder=4,
        )
        ax.plot(
            [row.anchor_x, 0.52, 0.57],
            [row.anchor_y, label_y, label_y],
            color="#6b7280",
            lw=0.9,
            zorder=2,
        )

        bar_width = 0.30 * enriched_count / max_count
        increased_width = bar_width * safe_ratio(float(increased_count), float(enriched_count))
        decreased_width = bar_width * safe_ratio(float(decreased_count), float(enriched_count))
        if increased_width > 0:
            ax.plot(
                [0.58, 0.58 + increased_width],
                [label_y, label_y],
                color=DIRECTION_COLORS["inc"],
                lw=7,
                solid_capstyle="round",
                zorder=3,
            )
        if decreased_width > 0:
            ax.plot(
                [0.58 + increased_width, 0.58 + increased_width + decreased_width],
                [label_y, label_y],
                color=DIRECTION_COLORS["dec"],
                lw=7,
                solid_capstyle="round",
                zorder=3,
            )

        ax.text(
            0.58,
            label_y + 0.018,
            display_body_site(row.body_site_group),
            ha="left",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )
        ax.text(
            0.58,
            label_y - 0.018,
            f"{enriched_count:,} enriched profiles; increased {increased_count:,}, "
            f"decreased {decreased_count:,}",
            ha="left",
            va="top",
            fontsize=8.5,
            color="#374151",
        )

    ax.text(
        0.02,
        0.965,
        "Anaerobe v2 Trait-Enriched Profiles: Anatomical Sites",
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
    )
    ax.text(
        0.02,
        0.925,
        GUILD_SPECIFIC_NOTE,
        ha="left",
        va="top",
        fontsize=8.5,
        color="#475569",
    )
    ax.text(
        0.02,
        0.035,
        "Retained = at least one positive selected-trait z-score.",
        ha="left",
        va="bottom",
        fontsize=8.5,
        color="#374151",
    )
    legend_handles = [
        Patch(
            color=DIRECTION_COLORS["inc"],
            label=(
                f"Increased retained (zero enriched: {zero_enriched_trait_counts.get('inc', 0):,})"
            ),
        ),
        Patch(
            color=DIRECTION_COLORS["dec"],
            label=(
                f"Decreased retained (zero enriched: {zero_enriched_trait_counts.get('dec', 0):,})"
            ),
        ),
    ]
    ax.legend(
        handles=legend_handles,
        title="Profile direction",
        loc="lower right",
        bbox_to_anchor=(0.98, 0.035),
        frameon=True,
        framealpha=1,
        edgecolor="#cbd5e1",
        fontsize=8,
        title_fontsize=8.5,
    )

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def display_body_site(value: object) -> str:
    """Return concise display labels for body-site groups."""
    label = clean_text(value)
    if label.lower() == "unknown":
        return "unknown body site"
    return label


def is_other_body_site(value: object) -> bool:
    return clean_text(value).lower() == "other"


def add_anatomical_outline(ax: plt.Axes) -> None:
    """Draw an original front-facing body-template outline for site callouts."""
    face = "#f8fafc"
    edge = "#94a3b8"

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
                lw=1.4,
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

    ax.add_patch(
        Circle(
            (0.30, 0.87),
            0.045,
            facecolor=face,
            edgecolor=edge,
            lw=1.4,
            zorder=1,
        )
    )
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
    ax.plot([0.280, 0.280], [0.805, 0.770], color=edge, lw=1.2, zorder=2)
    ax.plot([0.320, 0.320], [0.805, 0.770], color=edge, lw=1.2, zorder=2)


def save_anatomical_sites_diagram(
    body_site_summary: pd.DataFrame,
    *,
    output_path: Path,
    title: str,
    signal_column: str,
    signal_label: str,
    zero_enriched_trait_counts: dict[str, int] | None = None,
    max_sites: int = 12,
) -> None:
    """Save an anatomical-site callout diagram for body-site guild summaries."""
    if body_site_summary.empty:
        return

    plot_sites = (
        body_site_summary[~body_site_summary["body_site_group"].map(is_other_body_site)]
        .head(max_sites)
        .copy()
    )
    if plot_sites.empty:
        return
    plot_sites["profile_count"] = pd.to_numeric(
        plot_sites["profile_count"],
        errors="coerce",
    ).fillna(0)
    plot_sites[signal_column] = pd.to_numeric(
        plot_sites[signal_column],
        errors="coerce",
    ).fillna(0)
    points = [anatomical_site_base(site) for site in plot_sites["body_site_group"]]
    plot_sites["anchor_x"] = [point[0] for point in points]
    plot_sites["anchor_y"] = [point[1] for point in points]
    plot_sites = plot_sites.sort_values(
        ["anchor_y", "profile_count"],
        ascending=[False, False],
    )

    max_count = max(float(plot_sites["profile_count"].max()), 1.0)
    max_signal = max(float(plot_sites[signal_column].abs().max()), 0.1)
    norm = TwoSlopeNorm(vmin=-max_signal, vcenter=0, vmax=max_signal)
    cmap = plt.get_cmap("RdBu_r")

    fig_height = max(7.2, 1.0 + len(plot_sites) * 0.82)
    fig, ax = plt.subplots(figsize=(11, fig_height))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    add_anatomical_outline(ax)

    label_ys = [0.84]
    if len(plot_sites) > 1:
        step = 0.70 / (len(plot_sites) - 1)
        label_ys = [0.84 - index * step for index in range(len(plot_sites))]

    for label_y, row in zip(label_ys, plot_sites.itertuples(index=False), strict=False):
        count = int(row.profile_count)
        signal = float(getattr(row, signal_column))
        color = cmap(norm(signal))
        size = 95 + 900 * math.sqrt(count / max_count)
        ax.scatter(
            [row.anchor_x],
            [row.anchor_y],
            s=size,
            c=[color],
            edgecolors="#111827",
            linewidths=0.8,
            zorder=4,
        )
        ax.plot(
            [row.anchor_x, 0.52, 0.57],
            [row.anchor_y, label_y, label_y],
            color="#6b7280",
            lw=0.9,
            zorder=2,
        )
        bar_width = 0.30 * count / max_count
        ax.plot(
            [0.58, 0.58 + bar_width],
            [label_y, label_y],
            color=color,
            lw=7,
            solid_capstyle="round",
            zorder=3,
        )
        ax.text(
            0.58,
            label_y + 0.018,
            display_body_site(row.body_site_group),
            ha="left",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )
        ax.text(
            0.58,
            label_y - 0.018,
            f"{count:,} profiles; {signal_label} {signal:+.2f}",
            ha="left",
            va="top",
            fontsize=8.5,
            color="#374151",
        )

    ax.text(0.02, 0.965, title, ha="left", va="top", fontsize=17, fontweight="bold")
    ax.text(
        0.02,
        0.925,
        GUILD_SPECIFIC_NOTE,
        ha="left",
        va="top",
        fontsize=8.5,
        color="#475569",
    )
    ax.text(
        0.02,
        0.035,
        "Bubble and bar length show profile count; color shows "
        "direction-specific trait enrichment.",
        ha="left",
        va="bottom",
        fontsize=8.5,
        color="#374151",
    )
    if zero_enriched_trait_counts is not None:
        ax.text(
            0.72,
            0.035,
            "Zero enriched selected traits\n"
            f"increased profiles: {zero_enriched_trait_counts.get('inc', 0):,}\n"
            f"decreased profiles: {zero_enriched_trait_counts.get('dec', 0):,}",
            ha="left",
            va="bottom",
            fontsize=8.3,
            color="#374151",
            bbox={
                "boxstyle": "round,pad=0.35",
                "facecolor": "#ffffff",
                "edgecolor": "#cbd5e1",
                "linewidth": 0.8,
            },
        )
    scalar_mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_mappable.set_array([])
    colorbar = fig.colorbar(scalar_mappable, ax=ax, fraction=0.03, pad=0.02)
    colorbar.set_label(signal_label)

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def make_plots(
    taxa_summary: pd.DataFrame,
    condition_summary: pd.DataFrame,
    direction_summary: pd.DataFrame,
    body_site_summary: pd.DataFrame,
    body_site_enriched_profile_summary: pd.DataFrame,
    occurrences: pd.DataFrame,
    excluded_occurrences: pd.DataFrame,
    trait_enriched_profile_counts: pd.DataFrame,
    condition_enriched_profile_counts: pd.DataFrame,
    zero_enriched_trait_counts: dict[str, int],
) -> None:
    plt.rcParams.update({"figure.dpi": 160, "savefig.dpi": 200})

    fig, ax = plt.subplots(figsize=(9, 4.8))
    plot_taxa = taxa_summary.sort_values("depletion_bias")
    ax.barh(plot_taxa["target_taxon_name"], plot_taxa["depletion_bias"], color="#8c2d04")
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("(decreased profiles - increased profiles) / total profiles")
    ax.set_title("Anaerobe Candidate Taxa Are Decreased-Profile Biased", pad=28)
    add_specificity_note(ax, NON_GUILD_SPECIFIC_TARGET_TAXA_NOTE)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "candidate_taxa_depletion_bias.png")
    plt.close(fig)

    top_conditions = select_top_conditions(condition_summary, 12).sort_values(
        "decreased_minus_increased_guild_z"
    )
    fig, ax = plt.subplots(figsize=(11, 5.5))
    condition_values = top_conditions["decreased_minus_increased_guild_z"]
    ax.barh(
        top_conditions["condition"],
        condition_values,
        color="#2563eb",
    )
    if not top_conditions.empty:
        x_max = condition_values.max()
        x_pad = max(x_max * 0.025, 0.02)
        ax.set_xlim(0, x_max * 1.75)
        for y_pos, row in enumerate(top_conditions.itertuples(index=False)):
            ax.text(
                row.decreased_minus_increased_guild_z + x_pad,
                y_pos,
                paired_smd_label(row),
                va="center",
                fontsize=8,
                color="#334155",
            )
    ax.set_xlabel("Anaerobe v2 trait enrichment: decreased mean z - increased mean z")
    ax.set_title(
        "Conditions With Anaerobe v2 Trait Enrichment in Decreased Profiles",
        pad=28,
    )
    add_specificity_note(ax, GUILD_SPECIFIC_NOTE)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "condition_decreased_guild_signal.png")
    plt.close(fig)

    plot_direction = direction_summary[
        direction_summary["trait"].isin(GUILD_TRAITS.values())
    ].sort_values("decreased_minus_increased_z")
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.barh(
        plot_direction["trait"],
        plot_direction["decreased_minus_increased_z"],
        color="#047857",
    )
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("decreased profile mean z - increased profile mean z")
    ax.set_title(
        "Anaerobe v2 Traits Enriched in Decreased Profiles",
        pad=28,
    )
    add_specificity_note(ax, GUILD_SPECIFIC_NOTE)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "guild_trait_direction_shift.png")
    plt.close(fig)

    condition_counts = count_target_taxa_conditions(
        occurrences,
        "decreased",
    )
    excluded_note = excluded_host_note(excluded_occurrences, "decreased")
    save_condition_breakdown_bar_chart(
        condition_counts,
        output_path=OUTPUT_DIR / "target_taxa_condition_breakdown_bar.png",
        title_prefix="Dec-profile anaerobe target taxa",
        color="#2563eb",
        specificity_note=NON_GUILD_SPECIFIC_TARGET_TAXA_NOTE,
        note=excluded_note,
    )
    save_condition_target_taxa_stacked_bar_chart(
        occurrences,
        condition_summary,
        output_path=OUTPUT_DIR / "condition_target_taxa_stacked_bar.png",
        note=excluded_note,
    )
    save_trait_enriched_profile_histogram(
        trait_enriched_profile_counts,
        output_path=OUTPUT_DIR / "trait_enriched_profile_histogram.png",
    )
    save_condition_enriched_profile_histogram(
        condition_enriched_profile_counts,
        output_path=OUTPUT_DIR / "condition_enriched_profile_histogram.png",
    )
    save_anatomical_enriched_profile_diagram(
        body_site_enriched_profile_summary,
        zero_enriched_trait_counts,
        output_path=OUTPUT_DIR / "anatomical_sites_enriched_profile_diagram.png",
    )
    save_anatomical_sites_diagram(
        body_site_summary,
        output_path=OUTPUT_DIR / "anatomical_sites_diagram.png",
        title="Anaerobe v2 Trait Enrichment in Decreased Profiles: Anatomical Sites",
        signal_column="decreased_mean_guild_z",
        signal_label="decreased-profile trait z",
        zero_enriched_trait_counts=zero_enriched_trait_counts,
    )
    save_anatomical_sites_diagram(
        body_site_summary,
        output_path=OUTPUT_DIR / "anatomical_sites_increased_profile_diagram.png",
        title="Anaerobe v2 Trait Enrichment in Increased Profiles: Anatomical Sites",
        signal_column="increased_mean_guild_z",
        signal_label="increased-profile trait z",
        zero_enriched_trait_counts=zero_enriched_trait_counts,
    )


def write_summary_report(
    direction_summary: pd.DataFrame,
    body_site_summary: pd.DataFrame,
    condition_summary: pd.DataFrame,
    taxa_summary: pd.DataFrame,
    condition_taxa_summary: pd.DataFrame,
    trait_enriched_profile_counts: pd.DataFrame,
    condition_enriched_profile_counts: pd.DataFrame,
    body_site_enriched_profile_summary: pd.DataFrame,
    inferred_group_summary: pd.DataFrame,
    kg_microbe_taxon_universe: pd.DataFrame,
    top_depletion_profiles: pd.DataFrame,
    excluded_occurrences: pd.DataFrame,
) -> None:
    direction_rows = direction_summary[
        direction_summary["trait"].isin(GUILD_TRAITS.values())
    ].sort_values("decreased_minus_increased_z", ascending=False)
    top_conditions = select_top_conditions(condition_summary, 8)
    top_taxa = taxa_summary.head(7)
    top_taxon_conditions = condition_taxa_summary[
        (condition_taxa_summary["decreased_profile_count"] >= 3)
        & (condition_taxa_summary["decreased_study_count"] >= 2)
    ].head(10)
    enriched_direction_rows: list[dict[str, object]] = []
    if not condition_enriched_profile_counts.empty:
        for inc_dec, group in condition_enriched_profile_counts.groupby("inc_dec", dropna=False):
            nonzero = group[group["enriched_profile_count"] > 0]
            enriched_direction_rows.append(
                {
                    "direction": direction_label(inc_dec),
                    "condition_count": nonzero["condition"].nunique(),
                    "profile_count": int(group["profile_count"].sum()),
                    "enriched_profile_count": int(group["enriched_profile_count"].sum()),
                    "zero_enriched_profile_count": int(group["zero_enriched_profile_count"].sum()),
                    "enriched_profile_fraction": safe_ratio(
                        float(group["enriched_profile_count"].sum()),
                        float(group["profile_count"].sum()),
                    ),
                }
            )
    top_enriched_trait_rows = pd.DataFrame()
    if not trait_enriched_profile_counts.empty:
        top_enriched_trait_rows = trait_enriched_profile_counts.sort_values(
            [
                "trait_enriched_profile_count",
                "trait_enriched_study_count",
                "trait_name",
            ],
            ascending=[False, False, True],
        ).head(10)
    plotted_body_site_enriched_rows = pd.DataFrame()
    if not body_site_enriched_profile_summary.empty:
        plotted_body_site_enriched_rows = body_site_enriched_profile_summary[
            ~body_site_enriched_profile_summary["body_site_group"].map(is_other_body_site)
            & (body_site_enriched_profile_summary["enriched_profile_count"] >= 10)
        ].head(10)

    lines = [
        "# Anaerobe v2 Depletion Investigation",
        "",
        "This investigation uses the existing study trait profile clustered heatmap matrix, "
        "study metadata, taxa guild heatmap metadata, and derived taxa/trait matrices.",
        "",
        "## Scope",
        "",
        "- Selected traits: "
        + ", ".join(f"{name} ({trait_id})" for trait_id, name in GUILD_TRAITS.items()),
        "- Target taxa: " + ", ".join(TARGET_TAXA.values()),
        "- Formal KG-Microbe taxon universe: "
        f"{len(kg_microbe_taxon_universe):,} taxa with at least one core "
        "oxygen-tolerance annotation among " + ", ".join(KG_MICROBE_TAXON_UNIVERSE_TRAITS.values()),
        "- Trait selection audit: `guild_trait_selection.tsv` records the "
        "selected trait IDs, short rationale, and KG-Microbe taxon counts for "
        "each selected trait.",
        "- Context: all body sites, with body-site and condition summaries to separate "
        "gut, oral, airway, urogenital, and other patterns.",
        "",
        "## Plot Guild Specificity",
        "",
        "| plot | guild-specific? | basis |",
        "| --- | --- | --- |",
        "| `condition_decreased_guild_signal.png` | yes | v2 trait-enrichment score by condition |",
        "| `guild_trait_direction_shift.png` | yes | the two selected v2 traits |",
        "| `trait_enriched_profile_histogram.png` | yes | trait counts inside "
        "profile-level v2 trait-enrichment filtering |",
        "| `condition_enriched_profile_histogram.png` | yes | condition counts after "
        "profile-level v2 trait-enrichment filtering |",
        "| `anatomical_sites_enriched_profile_diagram.png` | yes | anatomical "
        "body-site counts after profile-level v2 trait-enrichment filtering |",
        "| `anatomical_sites_increased_profile_diagram.png` | yes | increased-profile "
        "v2 trait-enrichment score by body-site group |",
        "| `anatomical_sites_diagram.png` | yes | v2 trait-enrichment score by body-site group |",
        "| `candidate_taxa_depletion_bias.png` | no | fixed target-taxon direction counts |",
        "| `condition_target_taxa_stacked_bar.png` | no | fixed target-taxon occurrence "
        "counts stacked by taxon |",
        "| `target_taxa_condition_breakdown_bar.png` | no | fixed target-taxon occurrence counts |",
        "",
        "## Main Signal",
        "",
        "The selected anaerobe-only oxygen-tolerance traits are consistently enriched "
        "in decreased profiles relative to increased profiles. "
        "This is the decreased-profile counterpart to the aerobe v2 "
        "comparison, not a fermentation-defined gut trait set.",
        "",
        "Guardrail: the increased-profile anaerobe score is retained in every "
        "condition and body-site summary. It checks whether the anaerobe signal is "
        "decreased-profile-associated rather than simply high in both directions.",
        "",
        "| trait | decreased - increased z |",
        "| --- | ---: |",
    ]
    for row in direction_rows.itertuples(index=False):
        lines.append(f"| {row.trait} | {row.decreased_minus_increased_z:.3f} |")

    lines.extend(
        [
            "",
            "## Strongest Decreased-Profile-Biased Taxa",
            "",
            "| taxon | decreased / increased profiles | decreased-profile bias |",
            "| --- | ---: | ---: |",
        ]
    )
    for row in top_taxa.itertuples(index=False):
        lines.append(
            f"| {row.target_taxon_name} | {int(row.decreased)} / {int(row.increased)} | "
            f"{row.depletion_bias:.3f} |"
        )

    lines.extend(
        [
            "",
            "## Recurrent Target-Taxon Decreased Profiles By Condition",
            "",
            "| taxon | condition | body site | "
            "decreased / increased profiles | decreased studies |",
            "| --- | --- | --- | ---: | ---: |",
        ]
    )
    for row in top_taxon_conditions.itertuples(index=False):
        lines.append(
            f"| {row.target_taxon_name} | {clean_text(row.condition)} | "
            f"{clean_text(row.body_site_group)} | "
            f"{int(row.decreased_profile_count)} / {int(row.increased_profile_count)} | "
            f"{int(row.decreased_study_count)} |"
        )
    lines.extend(
        [
            "",
            "Target-taxon occurrence views (`target_taxa_condition_breakdown_bar.png`, "
            "`condition_target_taxa_stacked_bar.png`, `candidate_taxa_profile_occurrences.tsv` "
            "and `candidate_taxa_condition_direction.tsv`) count human-host profiles only, "
            "the corpus of the profile metadata. "
            f"{excluded_host_summary(excluded_occurrences)} "
            "They have no profile metadata, so they are not counted under any condition; "
            "`candidate_taxa_excluded_host_occurrences.tsv` lists them.",
        ]
    )

    lines.extend(
        [
            "",
            "## Conditions With Anaerobe v2 Trait Enrichment in Decreased Profiles",
            "",
            "| condition | profiles | studies | decreased - increased z | paired d | "
            "paired p | matched pairs | decreased mean trait z |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in top_conditions.itertuples(index=False):
        diff = row.decreased_minus_increased_guild_z
        diff_text = "" if pd.isna(diff) else f"{diff:.3f}"
        smd = row.paired_standardized_mean_difference
        smd_text = "" if pd.isna(smd) else f"{smd:.3f}"
        p_value = row.paired_t_test_p_value
        p_text = "" if pd.isna(p_value) else f"{p_value:.3g}"
        lines.append(
            f"| {row.condition} | {int(row.profile_count)} | {int(row.study_count)} | "
            f"{diff_text} | {smd_text} | {p_text} | {int(row.paired_profile_pair_count)} | "
            f"{row.decreased_mean_guild_z:.3f} |"
        )

    lines.extend(
        [
            "",
            "Paired d and paired p are blank for conditions with fewer than "
            f"{MIN_PAIRED_SMD_PROFILE_PAIRS} matched decreased/increased profile pairs, or with "
            "zero variance among the paired differences.",
        ]
    )

    if enriched_direction_rows:
        lines.extend(
            [
                "",
                "## Condition Enriched-Profile Histogram",
                "",
                "The enrichment-filtered histogram counts conditions by how many profiles "
                "remain after dropping profiles with zero positive anaerobe v2 selected-trait "
                "z-scores. The dropped profiles are still retained in the main trait-enrichment "
                "summaries and are reported here by profile direction.",
                "",
                "| direction | conditions with retained profiles | retained profiles | "
                "zero-enriched profiles | retained / total profiles |",
                "| --- | ---: | ---: | ---: | ---: |",
            ]
        )
        for row in sorted(enriched_direction_rows, key=lambda item: str(item["direction"])):
            lines.append(
                f"| {row['direction']} | {int(row['condition_count'])} | "
                f"{int(row['enriched_profile_count'])} | "
                f"{int(row['zero_enriched_profile_count'])} | "
                f"{row['enriched_profile_fraction']:.3f} |"
            )

    if not top_enriched_trait_rows.empty:
        lines.extend(
            [
                "",
                "## Trait Enriched-Profile Histogram",
                "",
                "The trait-side companion histogram counts selected heatmap traits by "
                "how many retained enriched profiles have a positive z-score for each "
                "trait. The retained-profile universe is the same one used by the "
                "condition histogram: profiles with at least one positive anaerobe v2 "
                "selected-trait z-score.",
                "",
                "| direction | trait | positive retained profiles | retained profile fraction |",
                "| --- | --- | ---: | ---: |",
            ]
        )
        for row in top_enriched_trait_rows.itertuples(index=False):
            lines.append(
                f"| {direction_label(row.inc_dec)} | {clean_text(row.trait_name)} | "
                f"{int(row.trait_enriched_profile_count)} | "
                f"{row.trait_enriched_profile_fraction:.3f} |"
            )

    if not plotted_body_site_enriched_rows.empty:
        lines.extend(
            [
                "",
                "## Enrichment-Filtered Anatomical Sites",
                "",
                "The anatomical companion keeps body-site groups with at least 10 "
                "profiles containing a positive anaerobe v2 selected-trait z-score. "
                "Counts below are retained profiles after zero-enriched profiles "
                "drop out.",
                "",
                "| body_site_group | retained profiles | increased retained | decreased retained | "
                "zero-enriched profiles |",
                "| --- | ---: | ---: | ---: | ---: |",
            ]
        )
        for row in plotted_body_site_enriched_rows.itertuples(index=False):
            lines.append(
                f"| {row.body_site_group} | {int(row.enriched_profile_count)} | "
                f"{int(row.increased_enriched_profile_count)} | "
                f"{int(row.decreased_enriched_profile_count)} | "
                f"{int(row.zero_enriched_profile_count)} |"
            )

    if not inferred_group_summary.empty:
        lines.extend(
            [
                "",
                "## Data-Inferred Guild Evidence",
                "",
                "For the top condition signals, this table uses all taxa "
                "observed in the decreased profiles for that condition that also "
                "belong to the KG-Microbe-wide anaerobe candidate taxon universe, "
                "then reports the profile v2 traits that are enriched in "
                "those same decreased profiles and the universe taxa that carry "
                "those trait annotations. This is a data-derived guild summary, "
                "not a fixed focal-taxon list.",
                "",
                "| condition | decreased profiles | universe taxa / observed taxa | "
                "enriched traits | supporting taxa |",
                "| --- | ---: | ---: | --- | --- |",
            ]
        )
        for row in inferred_group_summary.head(8).itertuples(index=False):
            lines.append(
                f"| {clean_text(row.condition)} | {int(row.focus_profile_count)} | "
                f"{int(row.kg_microbe_universe_taxon_count)} / "
                f"{int(row.total_observed_taxon_count)} | "
                f"{clean_text(row.top_enriched_traits)} | "
                f"{clean_text(row.top_supporting_taxa)} |"
            )

    lines.extend(
        [
            "",
            "## Body-Site Check",
            "",
            "| body_site_group | profiles | increased mean trait z | decreased mean trait z | "
            "decreased - increased z |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in body_site_summary.head(10).itertuples(index=False):
        diff = row.decreased_minus_increased_guild_z
        diff_text = "" if pd.isna(diff) else f"{diff:.3f}"
        lines.append(
            f"| {row.body_site_group} | {int(row.profile_count)} | "
            f"{row.increased_mean_guild_z:.3f} | {row.decreased_mean_guild_z:.3f} | "
            f"{diff_text} |"
        )

    lines.extend(
        [
            "",
            "## Top Decreased Profiles To Inspect",
            "",
            "| profile | condition | body site | trait z | title |",
            "| --- | --- | --- | ---: | --- |",
        ]
    )
    for row in top_depletion_profiles.head(10).itertuples(index=False):
        title = clean_text(row.title).replace("|", "/")
        if len(title) > 90:
            title = f"{title[:87]}..."
        profile_id = str(row.profile_id).replace("|", "\\|")
        lines.append(
            f"| {profile_id} | {clean_text(row.condition)} | {clean_text(row.body_site)} | "
            f"{row.guild_z_mean:.3f} | {title} |"
        )

    lines.extend(
        [
            "",
            "## Suggested Follow-Up",
            "",
            "1. Treat target-taxon recurrence as an audit view and the selected "
            "oxygen-tolerance traits as the formal trait-enrichment score.",
            "2. Prioritize conditions where decreased profiles have high trait-enrichment scores "
            "and target taxa "
            "are repeatedly reported in decreased profiles across multiple studies.",
            "3. Compare this anaerobe-only v2 signal against the original "
            "five-trait set to separate oxygen-tolerance signal from "
            "metabolic-capability signal.",
            "",
            "## Files",
            "",
            "- `guild_trait_direction_shift.tsv`",
            "- `guild_trait_selection.tsv`",
            "- `body_site_guild_summary.tsv`",
            "- `body_site_enriched_profile_summary.tsv`",
            "- `condition_guild_summary.tsv`",
            "- `trait_enriched_profile_counts.tsv`",
            "- `condition_enriched_profile_counts.tsv`",
            "- `candidate_taxa_direction_summary.tsv`",
            "- `candidate_taxa_trait_presence.tsv`",
            "- `candidate_taxa_condition_direction.tsv`",
            "- `candidate_taxa_profile_occurrences.tsv`",
            "- `candidate_taxa_excluded_host_occurrences.tsv`",
            "- `kg_microbe_candidate_taxon_universe.tsv`",
            "- `condition_inferred_guild_summary.tsv`",
            "- `condition_inferred_guild_trait_evidence.tsv`",
            "- `condition_inferred_guild_taxa_evidence.tsv`",
            "- `condition_body_site_guild_map_edges.tsv`",
            "- `condition_body_site_guild_map_sampled_body_sites.tsv`",
            "- `condition_body_site_guild_map_v2_edges.tsv`",
            "- `condition_body_site_guild_map_v2_sampled_body_sites.tsv`",
            "- `condition_body_site_guild_map_v3_individual_edges.tsv`",
            "- `condition_body_site_guild_map_v3_grouped_edges.tsv`",
            "- `condition_body_site_guild_map_v3_grouped_condition_rollups.tsv`",
            "- `condition_body_site_guild_map_v3_grouped_body_site_rollups.tsv`",
            "- `top_decreased_guild_profiles.tsv`",
            "- Methods are consolidated in `../../figure_methods.md`.",
            "- `candidate_taxa_depletion_bias.png`",
            "- `condition_decreased_guild_signal.png`",
            "- `condition_body_site_guild_map.png`",
            "- `condition_body_site_guild_map_v2.png`",
            "- `condition_body_site_guild_map_v3_individual.png`",
            "- `condition_body_site_guild_map_v3_grouped.png`",
            "- `trait_enriched_profile_histogram.png`",
            "- `condition_enriched_profile_histogram.png`",
            "- `condition_target_taxa_stacked_bar.png`",
            "- `guild_trait_direction_shift.png`",
            "- `anatomical_sites_diagram.png`",
            "- `anatomical_sites_increased_profile_diagram.png`",
            "- `anatomical_sites_enriched_profile_diagram.png`",
            "- `target_taxa_condition_breakdown_bar.png`",
            "",
        ]
    )
    (OUTPUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    study_matrix, study_metadata, study_trait_metadata = load_study_heatmap()
    scored_profiles = add_guild_scores(study_matrix, study_metadata)
    analysis_profiles = scored_profiles.copy()

    direction_summary = summarize_direction(analysis_profiles)
    body_site_summary = summarize_group(analysis_profiles, "body_site_group", min_profiles=10)
    body_site_enriched_profile_summary = summarize_body_site_enriched_profile_counts(
        analysis_profiles,
        min_profiles=10,
    )
    condition_summary = remove_excluded_condition_labels(
        summarize_group(analysis_profiles, "condition", min_profiles=5)
    )
    trait_enriched_profile_counts = summarize_trait_enriched_profile_counts(
        study_matrix,
        analysis_profiles,
        study_trait_metadata,
    )
    condition_enriched_profile_counts = summarize_condition_enriched_profile_counts(
        analysis_profiles
    )
    zero_enriched_trait_counts = count_zero_enriched_traits_by_direction(analysis_profiles)
    taxa_summary = load_taxa_direction_summary()
    trait_presence = load_target_taxa_trait_counts()
    occurrences, excluded_occurrences = load_target_taxa_occurrences(scored_profiles)
    condition_taxa_summary = summarize_taxa_conditions(occurrences)
    profile_taxa_occurrences = load_profile_taxa_occurrences(scored_profiles)
    kg_microbe_taxon_universe = build_kg_microbe_taxon_universe()
    kg_microbe_taxon_universe_ids = set(kg_microbe_taxon_universe["taxon_id"])

    inferred_group_summary, inferred_trait_evidence, inferred_taxa_evidence = (
        summarize_inferred_guild_evidence(
            study_matrix,
            analysis_profiles,
            study_trait_metadata,
            profile_taxa_occurrences,
            select_top_conditions(condition_summary, 8),
            ["condition"],
            focus_direction="decreased",
            taxon_universe_ids=kg_microbe_taxon_universe_ids,
        )
    )

    trait_definitions = study_trait_metadata.loc[
        study_trait_metadata.index.intersection(GUILD_TRAITS.keys())
    ].copy()
    trait_definitions["investigation_trait_name"] = trait_definitions.index.map(GUILD_TRAITS)
    guild_trait_selection = build_guild_trait_selection(study_trait_metadata)

    top_depletion_profiles = analysis_profiles[analysis_profiles["inc_dec"] == "dec"].copy()
    top_depletion_profiles = top_depletion_profiles.sort_values(
        ["guild_z_mean", "total_trait_signal"],
        ascending=[False, False],
    ).head(100)
    top_depletion_profiles = top_depletion_profiles.reset_index(names="profile_id")

    direction_summary.to_csv(OUTPUT_DIR / "guild_trait_direction_shift.tsv", sep="\t", index=False)
    body_site_summary.to_csv(OUTPUT_DIR / "body_site_guild_summary.tsv", sep="\t", index=False)
    body_site_enriched_profile_summary.to_csv(
        OUTPUT_DIR / "body_site_enriched_profile_summary.tsv",
        sep="\t",
        index=False,
    )
    condition_summary.to_csv(OUTPUT_DIR / "condition_guild_summary.tsv", sep="\t", index=False)
    trait_enriched_profile_counts.to_csv(
        OUTPUT_DIR / "trait_enriched_profile_counts.tsv",
        sep="\t",
        index=False,
    )
    condition_enriched_profile_counts.to_csv(
        OUTPUT_DIR / "condition_enriched_profile_counts.tsv",
        sep="\t",
        index=False,
    )
    taxa_summary.to_csv(OUTPUT_DIR / "candidate_taxa_direction_summary.tsv", sep="\t")
    trait_presence.to_csv(OUTPUT_DIR / "candidate_taxa_trait_presence.tsv", sep="\t", index=False)
    condition_taxa_summary.to_csv(
        OUTPUT_DIR / "candidate_taxa_condition_direction.tsv", sep="\t", index=False
    )
    occurrences.sort_values(["target_taxon_name", "profile_id"]).to_csv(
        OUTPUT_DIR / "candidate_taxa_profile_occurrences.tsv",
        sep="\t",
        index=False,
    )
    excluded_occurrences.sort_values(["target_taxon_name", "profile_id"]).to_csv(
        OUTPUT_DIR / "candidate_taxa_excluded_host_occurrences.tsv",
        sep="\t",
        index=False,
    )
    kg_microbe_taxon_universe.to_csv(
        OUTPUT_DIR / "kg_microbe_candidate_taxon_universe.tsv",
        sep="\t",
        index=False,
    )
    inferred_group_summary.to_csv(
        OUTPUT_DIR / "condition_inferred_guild_summary.tsv",
        sep="\t",
        index=False,
    )
    inferred_trait_evidence.to_csv(
        OUTPUT_DIR / "condition_inferred_guild_trait_evidence.tsv",
        sep="\t",
        index=False,
    )
    inferred_taxa_evidence.to_csv(
        OUTPUT_DIR / "condition_inferred_guild_taxa_evidence.tsv",
        sep="\t",
        index=False,
    )
    top_depletion_profiles.to_csv(
        OUTPUT_DIR / "top_decreased_guild_profiles.tsv",
        sep="\t",
        index=False,
    )
    trait_definitions.to_csv(OUTPUT_DIR / "guild_trait_definitions.tsv", sep="\t")
    guild_trait_selection.to_csv(
        OUTPUT_DIR / "guild_trait_selection.tsv",
        sep="\t",
        index=False,
    )

    make_plots(
        taxa_summary,
        condition_summary,
        direction_summary,
        body_site_summary,
        body_site_enriched_profile_summary,
        occurrences,
        excluded_occurrences,
        trait_enriched_profile_counts,
        condition_enriched_profile_counts,
        zero_enriched_trait_counts,
    )
    write_summary_report(
        direction_summary,
        body_site_summary,
        condition_summary,
        taxa_summary,
        condition_taxa_summary,
        trait_enriched_profile_counts,
        condition_enriched_profile_counts,
        body_site_enriched_profile_summary,
        inferred_group_summary,
        kg_microbe_taxon_universe,
        top_depletion_profiles,
        excluded_occurrences,
    )

    print(f"Wrote investigation outputs to: {OUTPUT_DIR}")
    print(f"Study profiles analyzed: {len(scored_profiles)}")
    print(f"Analysis profiles analyzed: {len(analysis_profiles)}")
    print(format_zero_enriched_trait_counts(zero_enriched_trait_counts))
    print(f"Target taxa occurrence rows: {len(occurrences)}")
    print(f"Target taxa occurrences outside the human-host corpus: {len(excluded_occurrences)}")
    print(f"KG-Microbe candidate taxon universe rows: {len(kg_microbe_taxon_universe)}")
    print(f"Inferred guild condition rows: {len(inferred_group_summary)}")
    print(f"Condition/taxon rows: {len(condition_taxa_summary)}")


if __name__ == "__main__":
    main()
