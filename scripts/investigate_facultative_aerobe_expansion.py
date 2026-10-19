"""Investigate aerobe v2 trait enrichment in increased profiles."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path as MplPath

from mtp.body_site_diagram import (
    add_colorbar_beside,
    anatomical_site_base,
    bubble_radius,
    leader_line,
)
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
OUTPUT_DIR = VIS_DIR / "trait-profile-investigation" / "aerobe-enrichment-v2"
PROFILE_MATRIX_PATH = REPO_ROOT / "data/derived/study_trait_profiles_wide.tsv"

GUILD_TRAITS = {
    "METPO:1000608": "facultatively aerobic",
    "METPO:1000602": "aerobic",
    "METPO:1000604": "microaerophilic",
}
GUILD_TRAIT_SELECTION_JUSTIFICATIONS = {
    "METPO:1000608": (
        "Primary facultatively aerobic concept; selected as one core trait for the "
        "increased-profile aerobe v2 oxygen-tolerance trait set."
    ),
    "METPO:1000602": (
        "Oxygen-positive growth trait; selected to capture taxa whose profile signal is "
        "consistent with aerobic oxygen tolerance."
    ),
    "METPO:1000604": (
        "Low-oxygen growth trait; selected to include microaerophilic oxygen tolerance "
        "without broadening the score to unrelated metabolic traits."
    ),
}
KG_MICROBE_TAXON_UNIVERSE_TRAITS = {
    "METPO:1000608": "facultatively aerobic",
    "METPO:1000602": "aerobic",
    "METPO:1000604": "microaerophilic",
}
KG_MICROBE_TAXON_UNIVERSE_ASSOCIATED_TRAITS: dict[str, str] = {}

TARGET_TAXA = {
    "NCBITaxon:33033": "Parvimonas micra",
    "NCBITaxon:573": "Klebsiella pneumoniae",
    "NCBITaxon:1328": "Streptococcus anginosus",
    "NCBITaxon:562": "Escherichia coli",
    "NCBITaxon:851": "Fusobacterium nucleatum",
    "NCBITaxon:208479": "Enterocloster bolteae",
    "NCBITaxon:84112": "Eggerthella lenta",
    "NCBITaxon:28131": "Prevotella intermedia",
    "NCBITaxon:837": "Porphyromonas gingivalis",
    "NCBITaxon:39778": "Veillonella dispar",
}
EXCLUDED_CONDITION_SUMMARY_LABELS = {
    "Body fat distribution",
    "Exposure",
    "Health trait",
    "Microbiome",
    "Population",
    "Physical activity",
    "Time",
    "Timepoint",
}
# A paired d (and its paired t-test p) from fewer matched profile pairs is reported as NA.
MIN_PAIRED_SMD_PROFILE_PAIRS = 5
INFERRED_GUILD_TOP_TRAITS_PER_GROUP = 12
INFERRED_GUILD_TOP_TAXA_PER_TRAIT = 12
INFERRED_GUILD_TOP_TAXA_PER_GROUP = 12
MIN_KG_MICROBE_UNIVERSE_CORE_TRAITS = 1
GUILD_SPECIFIC_NOTE = (
    "Guild-specific: yes - uses the aerobe v2 trait-enrichment score "
    "(facultatively aerobic, aerobic, microaerophilic)."
)
NON_GUILD_SPECIFIC_TARGET_TAXA_NOTE = (
    "Guild-specific: no - fixed target-taxon count view; "
    "not recomputed from the trait-enrichment score."
)
TARGET_TAXON_COLORS = {
    "Escherichia coli": "#dc2626",
    "Fusobacterium nucleatum": "#2563eb",
    "Veillonella dispar": "#0891b2",
    "Klebsiella pneumoniae": "#d97706",
    "Enterocloster bolteae": "#7c3aed",
    "Parvimonas micra": "#059669",
    "Streptococcus anginosus": "#be123c",
    "Eggerthella lenta": "#64748b",
    "Porphyromonas gingivalis": "#9333ea",
    "Prevotella intermedia": "#65a30d",
}
TARGET_TAXON_LEGEND_LABELS = {
    "Escherichia coli": "E. coli",
    "Fusobacterium nucleatum": "F. nucleatum",
    "Veillonella dispar": "V. dispar",
    "Klebsiella pneumoniae": "K. pneumoniae",
    "Enterocloster bolteae": "E. bolteae",
    "Parvimonas micra": "P. micra",
    "Streptococcus anginosus": "S. anginosus",
    "Eggerthella lenta": "E. lenta",
    "Porphyromonas gingivalis": "P. gingivalis",
    "Prevotella intermedia": "P. intermedia",
}


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
    """Calculate paired inc-vs-dec trait enrichment within study/experiment."""
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
    differences = paired["inc"] - paired["dec"]
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
        msg = "No requested aerobe v2 trait-enrichment traits were found in the matrix."
        raise ValueError(msg)
    scored = metadata.copy()
    scored["guild_z_mean"] = matrix[available_traits].mean(axis=1)
    scored["oxygen_tolerance_z_mean"] = scored["guild_z_mean"]
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
                "increased_minus_decreased_z": inc_mean - dec_mean,
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
                "increased_minus_decreased_z": inc_mean - dec_mean,
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
            "increased_minus_decreased_guild_z": inc["guild_z_mean"].mean()
            - dec["guild_z_mean"].mean(),
            "increased_mean_oxygen_tolerance_z": inc["oxygen_tolerance_z_mean"].mean(),
        }
        if group_column == "condition":
            record.update(paired_guild_enrichment_stats(group))
        records.append(record)
    return pd.DataFrame(records).sort_values(
        ["increased_mean_guild_z", "profile_count"],
        ascending=[False, False],
    )


def remove_excluded_condition_labels(condition_summary: pd.DataFrame) -> pd.DataFrame:
    """Remove broad study-design labels from condition summaries."""
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
            (condition_summary["increased_profile_count"] >= 3)
            & (condition_summary["study_count"] >= 2)
            & (condition_summary["increased_minus_decreased_guild_z"] > 0)
            & (condition_summary["increased_mean_guild_z"] > 0)
        ]
        .sort_values(
            ["increased_minus_decreased_guild_z", "increased_mean_guild_z"],
            ascending=[False, False],
        )
        .head(limit)
    )
    if top_conditions.empty:
        top_conditions = condition_summary[condition_summary["increased_profile_count"] >= 3].head(
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
    """Document why each selected trait is in the score."""
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
                "guild": "aerobe v2",
                "focus_direction": "increased",
                "comparator_direction": "decreased",
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
        "oxygen_tolerance_z_mean",
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
                "kg_microbe_universe_taxon_count": group_traits["observed_taxon_count"].max(),
                "inferred_enriched_trait_count": group_traits["trait_id"].nunique(),
                "supporting_taxon_count": selected_taxa["taxon_id"].nunique(),
                "top_enriched_traits": limited_join(top_trait_text, 5),
                "top_supporting_taxa": limited_join(
                    top_taxa_text,
                    INFERRED_GUILD_TOP_TAXA_PER_GROUP,
                ),
            }
        )

    group_summary = pd.DataFrame(group_summary_rows)
    return group_summary, trait_summary, taxon_summary


def load_taxa_direction_summary() -> pd.DataFrame:
    taxa_metadata = read_tsv(VIS_DIR / "taxa_guild_discovery_taxa_metadata.tsv", index_col=0)
    summary = taxa_metadata.loc[taxa_metadata.index.intersection(TARGET_TAXA.keys())].copy()
    summary["target_taxon_name"] = summary.index.map(TARGET_TAXA)
    summary["increased_fraction"] = [
        safe_ratio(float(row["increased"]), float(row["bugsigdb_profile_count"]))
        for _, row in summary.iterrows()
    ]
    summary["increased_profile_tendency"] = [
        safe_ratio(
            float(row["increased"]) - float(row["decreased"]),
            float(row["bugsigdb_profile_count"]),
        )
        for _, row in summary.iterrows()
    ]
    summary["increased_to_decreased_ratio"] = [
        safe_ratio(float(row["increased"]), float(row["decreased"]))
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
            "increased_fraction",
            "increased_profile_tendency",
            "increased_to_decreased_ratio",
            "trait_count",
        ]
    ].sort_values("increased_profile_tendency", ascending=False)


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
        "oxygen_tolerance_z_mean",
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
    output["increased_profile_tendency"] = (
        output["increased_profile_count"] - output["decreased_profile_count"]
    ) / output["total_profile_count"].replace(0, pd.NA)
    return output.sort_values(
        ["increased_profile_tendency", "increased_study_count", "total_profile_count"],
        ascending=[False, False, False],
    )


def count_target_taxa_conditions(occurrences: pd.DataFrame, direction: str) -> pd.Series:
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
    add_specificity_note(ax, NON_GUILD_SPECIFIC_TARGET_TAXA_NOTE)
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
    """Return the condition's guild-score contrast label when available."""
    if condition_summary.empty:
        return ""
    matches = condition_summary[condition_summary["condition"] == condition]
    if matches.empty:
        return ""
    value = matches.iloc[0]["increased_minus_decreased_guild_z"]
    if pd.isna(value):
        return ""
    if abs(value) < 0.005:
        value = 0
    return f"; guild dz={value:.2f}"


def save_condition_target_taxa_stacked_bar_chart(
    occurrences: pd.DataFrame,
    condition_summary: pd.DataFrame,
    *,
    output_path: Path,
    note: str,
    limit: int = 12,
) -> None:
    """Save a stacked bar chart of increased-profile focal-taxon recurrence by condition.

    `note` follows the occurrence-count caveat in the footer.
    """
    selected = occurrences[occurrences["direction"] == "increased"].copy()
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
    ax.set_xlabel("Increased-profile target-taxon occurrences")
    ax.set_title(
        "Increased-Profile Aerobe v2 Target-Taxon Recurrence by Condition",
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
    legend_note: str | None = None,
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

    label_texts = []
    for label_y, row in zip(label_ys, plot_sites.itertuples(index=False), strict=False):
        count = int(row.profile_count)
        signal = float(getattr(row, signal_column))
        color = cmap(norm(signal))
        radius = bubble_radius(count, max_count)
        ax.add_patch(
            Circle(
                (row.anchor_x, row.anchor_y),
                radius,
                facecolor=color,
                edgecolor="#111827",
                linewidth=0.8,
                zorder=4,
            )
        )
        ax.plot(
            *leader_line((row.anchor_x, row.anchor_y), radius, label_y, 0.52, 0.57),
            color="#6b7280",
            lw=0.9,
            zorder=5,
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
        label_texts.append(
            ax.text(
                0.58,
                label_y + 0.018,
                display_body_site(row.body_site_group),
                ha="left",
                va="bottom",
                fontsize=10,
                fontweight="bold",
            )
        )
        label_texts.append(
            ax.text(
                0.58,
                label_y - 0.018,
                f"{count:,} profiles; {signal_label} {signal:+.2f}",
                ha="left",
                va="top",
                fontsize=8.5,
                color="#374151",
            )
        )

    # Title, note and footer sit outside the outline, clear of the head and feet.
    ax.text(0.02, 1.005, title, ha="left", va="bottom", fontsize=17, fontweight="bold")
    ax.text(
        0.02,
        0.995,
        GUILD_SPECIFIC_NOTE,
        ha="left",
        va="top",
        fontsize=8.5,
        color="#475569",
    )
    footer = (
        "Bubble and bar length show profile count; color shows direction-specific trait enrichment."
    )
    if legend_note:
        footer = f"{footer}\n{legend_note}"
    ax.text(
        0.02,
        0.02,
        footer,
        ha="left",
        va="top",
        fontsize=8.5,
        color="#374151",
    )
    scalar_mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_mappable.set_array([])

    fig.tight_layout()
    add_colorbar_beside(fig, ax, scalar_mappable, label_texts, signal_label)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def make_plots(
    taxa_summary: pd.DataFrame,
    condition_summary: pd.DataFrame,
    direction_summary: pd.DataFrame,
    body_site_summary: pd.DataFrame,
    occurrences: pd.DataFrame,
    excluded_occurrences: pd.DataFrame,
) -> None:
    plt.rcParams.update({"figure.dpi": 160, "savefig.dpi": 200})

    fig, ax = plt.subplots(figsize=(9, 5.4))
    plot_taxa = taxa_summary.sort_values("increased_profile_tendency")
    ax.barh(
        plot_taxa["target_taxon_name"],
        plot_taxa["increased_profile_tendency"],
        color="#b91c1c",
    )
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("(increased profiles - decreased profiles) / total profiles")
    ax.set_title("Aerobe v2 Target Taxa: Increased-Profile Tendency", pad=28)
    add_specificity_note(ax, NON_GUILD_SPECIFIC_TARGET_TAXA_NOTE)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "candidate_taxa_increased_profile_tendency.png")
    plt.close(fig)

    top_conditions = select_top_conditions(condition_summary, 12).sort_values(
        "increased_minus_decreased_guild_z"
    )
    fig, ax = plt.subplots(figsize=(11, 5.5))
    condition_values = top_conditions["increased_minus_decreased_guild_z"]
    ax.barh(
        top_conditions["condition"],
        condition_values,
        color="#dc2626",
    )
    if not top_conditions.empty:
        x_max = condition_values.max()
        x_pad = max(x_max * 0.025, 0.02)
        ax.set_xlim(0, x_max * 1.75)
        for y_pos, row in enumerate(top_conditions.itertuples(index=False)):
            ax.text(
                row.increased_minus_decreased_guild_z + x_pad,
                y_pos,
                paired_smd_label(row),
                va="center",
                fontsize=8,
                color="#334155",
            )
    ax.set_xlabel("Aerobe v2 trait enrichment: increased mean z - decreased mean z")
    ax.set_title(
        "Conditions With Aerobe v2 Trait Enrichment in Increased Profiles",
        pad=28,
    )
    add_specificity_note(ax, GUILD_SPECIFIC_NOTE)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "condition_increased_guild_signal.png")
    plt.close(fig)

    plot_direction = direction_summary[
        direction_summary["trait"].isin(GUILD_TRAITS.values())
    ].sort_values("increased_minus_decreased_z")
    fig, ax = plt.subplots(figsize=(8, 4.4))
    ax.barh(
        plot_direction["trait"],
        plot_direction["increased_minus_decreased_z"],
        color="#be123c",
    )
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("increased profile mean z - decreased profile mean z")
    ax.set_title("Aerobe v2 Traits Enriched in Increased Profiles", pad=28)
    add_specificity_note(ax, GUILD_SPECIFIC_NOTE)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "guild_trait_direction_shift.png")
    plt.close(fig)

    condition_counts = count_target_taxa_conditions(occurrences, "increased")
    excluded_note = excluded_host_note(excluded_occurrences, "increased")
    save_condition_breakdown_bar_chart(
        condition_counts,
        output_path=OUTPUT_DIR / "target_taxa_condition_breakdown_bar.png",
        title_prefix="Inc-profile aerobe v2 target taxa",
        color="#dc2626",
        note=excluded_note,
    )
    save_condition_target_taxa_stacked_bar_chart(
        occurrences,
        condition_summary,
        output_path=OUTPUT_DIR / "condition_target_taxa_stacked_bar.png",
        note=excluded_note,
    )
    save_anatomical_sites_diagram(
        body_site_summary,
        output_path=OUTPUT_DIR / "anatomical_sites_diagram.png",
        title="Aerobe v2 Trait Enrichment in Increased Profiles: Anatomical Sites",
        signal_column="increased_mean_guild_z",
        signal_label="increased-profile trait z",
        legend_note=(
            "Unknown body site means no specific BugSigDB body-site label was available; "
            "it is plotted off-body."
        ),
    )


def write_summary_report(
    direction_summary: pd.DataFrame,
    body_site_summary: pd.DataFrame,
    condition_summary: pd.DataFrame,
    taxa_summary: pd.DataFrame,
    condition_taxa_summary: pd.DataFrame,
    inferred_group_summary: pd.DataFrame,
    kg_microbe_taxon_universe: pd.DataFrame,
    top_increased_profiles: pd.DataFrame,
    excluded_occurrences: pd.DataFrame,
) -> None:
    direction_rows = direction_summary[
        direction_summary["trait"].isin(GUILD_TRAITS.values())
    ].sort_values("increased_minus_decreased_z", ascending=False)
    top_conditions = select_top_conditions(condition_summary, 8)
    top_taxa = taxa_summary.head(10)
    top_taxon_conditions = condition_taxa_summary[
        (condition_taxa_summary["increased_profile_count"] >= 3)
        & (condition_taxa_summary["increased_study_count"] >= 2)
    ].head(12)

    lines = [
        "# Aerobe v2 Trait Enrichment in Increased Profiles",
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
        "gut, oral, airway, and other patterns.",
        "",
        "## Plot Guild Specificity",
        "",
        "| plot | guild-specific? | basis |",
        "| --- | --- | --- |",
        "| `condition_increased_guild_signal.png` | yes | v2 trait-enrichment score by condition |",
        "| `guild_trait_direction_shift.png` | yes | the three selected v2 traits |",
        "| `anatomical_sites_diagram.png` | yes | v2 trait-enrichment score by body-site group |",
        "| `candidate_taxa_increased_profile_tendency.png` | no | "
        "fixed target-taxon direction counts |",
        "| `condition_target_taxa_stacked_bar.png` | no | fixed target-taxon occurrence "
        "counts stacked by taxon |",
        "| `target_taxa_condition_breakdown_bar.png` | no | fixed target-taxon occurrence counts |",
        "",
        "## Main Signal",
        "",
        "The selected traits are enriched in increased profiles relative to decreased "
        "profiles, supporting an aerobe v2 increased-profile trait-enrichment signal.",
        "",
        "| trait | increased - decreased z |",
        "| --- | ---: |",
    ]
    for row in direction_rows.itertuples(index=False):
        lines.append(f"| {row.trait} | {row.increased_minus_decreased_z:.3f} |")

    lines.extend(
        [
            "",
            "## Target Taxa Most Often Reported In Increased Profiles",
            "",
            "| taxon | increased / decreased profiles | increased-profile tendency |",
            "| --- | ---: | ---: |",
        ]
    )
    for row in top_taxa.itertuples(index=False):
        lines.append(
            f"| {row.target_taxon_name} | {int(row.increased)} / {int(row.decreased)} | "
            f"{row.increased_profile_tendency:.3f} |"
        )

    lines.extend(
        [
            "",
            "## Recurrent Target-Taxon Increased Profiles By Condition",
            "",
            "| taxon | condition | body site | "
            "increased / decreased profiles | increased studies |",
            "| --- | --- | --- | ---: | ---: |",
        ]
    )
    for row in top_taxon_conditions.itertuples(index=False):
        lines.append(
            f"| {row.target_taxon_name} | {clean_text(row.condition)} | "
            f"{clean_text(row.body_site_group)} | "
            f"{int(row.increased_profile_count)} / {int(row.decreased_profile_count)} | "
            f"{int(row.increased_study_count)} |"
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
            "## Contexts With Aerobe v2 Trait Enrichment in Increased Profiles",
            "",
            "| condition | profiles | studies | increased - decreased z | paired d | "
            "paired p | matched pairs | increased mean trait z |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in top_conditions.itertuples(index=False):
        diff = row.increased_minus_decreased_guild_z
        diff_text = "" if pd.isna(diff) else f"{diff:.3f}"
        smd = row.paired_standardized_mean_difference
        smd_text = "" if pd.isna(smd) else f"{smd:.3f}"
        p_value = row.paired_t_test_p_value
        p_text = "" if pd.isna(p_value) else f"{p_value:.3g}"
        lines.append(
            f"| {row.condition} | {int(row.profile_count)} | {int(row.study_count)} | "
            f"{diff_text} | {smd_text} | {p_text} | {int(row.paired_profile_pair_count)} | "
            f"{row.increased_mean_guild_z:.3f} |"
        )

    lines.extend(
        [
            "",
            "Paired d and paired p are blank for conditions with fewer than "
            f"{MIN_PAIRED_SMD_PROFILE_PAIRS} matched increased/decreased profile pairs, or with "
            "zero variance among the paired differences.",
            "",
            "Broad labels such as `Exposure`, `Population`, `Time`, `Timepoint`, and labels "
            "containing `measurement` come from BugSigDB condition metadata but usually "
            "describe study design, cohort, sampling time, exposure, or measurement context "
            "rather than a disease/condition context. They are excluded from displayed "
            "condition bar charts and condition summaries.",
            "",
            "Oral-site increased-profile trait signals are expected in this analysis because "
            "the oral cavity is oxygen-exposed and several focal taxa are common oral "
            "or airway-associated aerotolerant, facultatively aerobic, or "
            "microaerophilic lineages. These signals should be interpreted as "
            "body-site-specific oral microbial ecology, not as a gut "
            "increased-profile trait-enrichment signal.",
        ]
    )

    if not inferred_group_summary.empty:
        lines.extend(
            [
                "",
                "## Data-Inferred Guild Evidence",
                "",
                "For the top condition signals, this table uses all taxa observed in the "
                "increased profiles for that condition that also belong to the "
                "KG-Microbe-wide candidate taxon universe, then reports the profile "
                "v2 traits that are enriched in those same increased profiles and the "
                "universe taxa that carry those trait annotations. This is a data-derived "
                "guild summary, not a literature-defined aerobe/aerotolerant "
                "guild.",
                "",
                "| condition | increased profiles | universe taxa / observed taxa | "
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
            "| body_site_group | profiles | increased mean trait z | increased - decreased z |",
            "| --- | ---: | ---: | ---: |",
        ]
    )
    for row in body_site_summary.head(10).itertuples(index=False):
        diff = row.increased_minus_decreased_guild_z
        diff_text = "" if pd.isna(diff) else f"{diff:.3f}"
        lines.append(
            f"| {row.body_site_group} | {int(row.profile_count)} | "
            f"{row.increased_mean_guild_z:.3f} | {diff_text} |"
        )

    lines.extend(
        [
            "",
            "## Top Increased Profiles To Inspect",
            "",
            "| profile | condition | body site | trait z | title |",
            "| --- | --- | --- | ---: | --- |",
        ]
    )
    for row in top_increased_profiles.head(10).itertuples(index=False):
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
            "1. Split the signal by body site before interpreting it as an "
            "aerobe v2 increased-profile trait-enrichment signal; "
            "oral and airway profiles can dominate aerobic/microaerophilic traits.",
            "2. In gut-only analyses, test the target taxa against inflammation markers, "
            "oxygen tolerance, virulence genes, and resistome features.",
            "3. Compare increased-profile taxa against the anaerobe v2 trait set to see "
            "whether the same conditions show reciprocal decreased- and "
            "increased-profile patterns.",
            "",
            "## Files",
            "",
            "- `guild_trait_direction_shift.tsv`",
            "- `guild_trait_selection.tsv`",
            "- `body_site_guild_summary.tsv`",
            "- `condition_guild_summary.tsv`",
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
            "- `top_increased_guild_profiles.tsv`",
            "- Methods are consolidated in `../../figure_methods.md`.",
            "- `candidate_taxa_increased_profile_tendency.png`",
            "- `condition_increased_guild_signal.png`",
            "- `condition_body_site_guild_map.png`",
            "- `condition_body_site_guild_map_v2.png`",
            "- `condition_body_site_guild_map_v3_individual.png`",
            "- `condition_body_site_guild_map_v3_grouped.png`",
            "- `condition_target_taxa_stacked_bar.png`",
            "- `guild_trait_direction_shift.png`",
            "- `anatomical_sites_diagram.png`",
            "- `target_taxa_condition_breakdown_bar.png`",
            "",
        ]
    )
    (OUTPUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    study_matrix, study_metadata, study_trait_metadata = load_study_heatmap()
    scored_profiles = add_guild_scores(study_matrix, study_metadata)

    direction_summary = summarize_direction(scored_profiles)
    body_site_summary = summarize_group(scored_profiles, "body_site_group", min_profiles=10)
    condition_summary = summarize_group(scored_profiles, "condition", min_profiles=5)
    condition_summary_output = remove_excluded_condition_labels(condition_summary)
    taxa_summary = load_taxa_direction_summary()
    trait_presence = load_target_taxa_trait_counts()
    occurrences, excluded_occurrences = load_target_taxa_occurrences(scored_profiles)
    condition_taxa_summary = summarize_taxa_conditions(occurrences)
    profile_taxa_occurrences = load_profile_taxa_occurrences(scored_profiles)
    kg_microbe_taxon_universe = build_kg_microbe_taxon_universe()
    kg_microbe_taxon_universe_ids = set(kg_microbe_taxon_universe["taxon_id"])
    guild_trait_selection = build_guild_trait_selection(study_trait_metadata)

    inferred_group_summary, inferred_trait_evidence, inferred_taxa_evidence = (
        summarize_inferred_guild_evidence(
            study_matrix,
            scored_profiles,
            study_trait_metadata,
            profile_taxa_occurrences,
            select_top_conditions(condition_summary, 8),
            ["condition"],
            focus_direction="increased",
            taxon_universe_ids=kg_microbe_taxon_universe_ids,
        )
    )

    trait_definitions = study_trait_metadata.loc[
        study_trait_metadata.index.intersection(GUILD_TRAITS.keys())
    ].copy()
    trait_definitions["investigation_trait_name"] = trait_definitions.index.map(GUILD_TRAITS)

    top_increased_profiles = scored_profiles[scored_profiles["inc_dec"] == "inc"].copy()
    top_increased_profiles = top_increased_profiles.sort_values(
        ["guild_z_mean", "total_trait_signal"],
        ascending=[False, False],
    ).head(100)
    top_increased_profiles = top_increased_profiles.reset_index(names="profile_id")

    direction_summary.to_csv(OUTPUT_DIR / "guild_trait_direction_shift.tsv", sep="\t", index=False)
    body_site_summary.to_csv(OUTPUT_DIR / "body_site_guild_summary.tsv", sep="\t", index=False)
    condition_summary_output.to_csv(
        OUTPUT_DIR / "condition_guild_summary.tsv",
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
    top_increased_profiles.to_csv(
        OUTPUT_DIR / "top_increased_guild_profiles.tsv",
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
        occurrences,
        excluded_occurrences,
    )
    write_summary_report(
        direction_summary,
        body_site_summary,
        condition_summary,
        taxa_summary,
        condition_taxa_summary,
        inferred_group_summary,
        kg_microbe_taxon_universe,
        top_increased_profiles,
        excluded_occurrences,
    )

    print(f"Wrote investigation outputs to: {OUTPUT_DIR}")
    print(f"Study profiles analyzed: {len(scored_profiles)}")
    print(f"Target taxa occurrence rows: {len(occurrences)}")
    print(f"Target taxa occurrences outside the human-host corpus: {len(excluded_occurrences)}")
    print(f"KG-Microbe candidate taxon universe rows: {len(kg_microbe_taxon_universe)}")
    print(f"Inferred guild condition rows: {len(inferred_group_summary)}")
    print(f"Condition/taxon rows: {len(condition_taxa_summary)}")


if __name__ == "__main__":
    main()
