"""Create an FDR-tiered fecal condition-taxon-trait evidence landscape."""

from __future__ import annotations

import math
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from scipy import stats

from mtp.composition_denominators import load_retained_profiles, trait_fractions

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "trait-visualization"
INVESTIGATION_DIR = VIS_DIR / "trait-profile-investigation"
SIX_SHARED_DIR = INVESTIGATION_DIR / "six-shared-traits"
OUTPUT_DIR = INVESTIGATION_DIR / "PD-cellulose"

PROFILE_MATRIX_PATH = REPO_ROOT / "data/derived/study_trait_profiles_wide.tsv"
PROFILE_METADATA_PATH = VIS_DIR / "study_trait_profile_study_metadata.tsv"
TAXA_PATH = REPO_ROOT / "data/derived/bugsigdb_taxa.tsv"
TAXA_TRAIT_MATRIX_PATH = REPO_ROOT / "data/derived/taxa_trait_matrix.tsv"
TRAIT_SUMMARY_PATH = SIX_SHARED_DIR / "six_shared_trait_summary.tsv"

OUTPUT_STEM = "pd_condition_taxon_trait_network"
FECES = "Feces"
HUMAN_HOST = "Homo sapiens"
PD_CONDITION = "Parkinson's disease"
Q_VALUE_CUTOFF = 0.05
EXPLORATORY_Q_CUTOFF = 0.20
MIN_PAIRED_STUDIES = 3
MIN_TAXON_STUDIES = 3
ALLOWED_TAXON_RANKS = {"genus", "species"}
DIRECTION_LABELS = {"inc": "increased", "dec": "decreased"}

INCREASED_COLOR = "#c2410c"
DECREASED_COLOR = "#2563eb"
ANNOTATION_COLOR = "#94a3b8"
PD_COLOR = "#d55e00"
SLATE_50 = "#f8fafc"
SLATE_100 = "#f1f5f9"
SLATE_200 = "#e2e8f0"
SLATE_300 = "#cbd5e1"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"
TRAIT_COLORS = {
    "facultatively aerobic": "#0072b2",
    "cellulose catabolic process": "#d55e00",
    "coccus shaped": "#009e73",
    "sphere shaped": "#cc79a7",
    "arsenate reduction": "#e69f00",
    "moderately halophilic": "#6a3d9a",
}


def read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    """Read a local TSV file."""
    return pd.read_csv(path, sep="\t", **kwargs)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    """Write a stable audit TSV with explicit missing values."""
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.8g")


def clean_text(value: object) -> str:
    """Return a stripped text value for identifiers and labels."""
    if pd.isna(value):
        return ""
    return str(value).strip()


def wrap_label(value: object, width: int = 20) -> str:
    """Wrap a display label without breaking long words."""
    text = clean_text(value)
    return "\n".join(wrap(text, width=width, break_long_words=False)) or "unknown"


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """Return Benjamini-Hochberg adjusted q-values for finite p-values."""
    values = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    adjusted = np.full_like(values, np.nan, dtype=float)
    valid_indices = np.flatnonzero(np.isfinite(values))
    if len(valid_indices) == 0:
        return pd.Series(adjusted, index=p_values.index, dtype=float)
    sorted_indices = valid_indices[np.argsort(values[valid_indices])]
    sorted_p = values[sorted_indices]
    ranks = np.arange(1, len(sorted_p) + 1, dtype=float)
    sorted_adjusted = np.minimum.accumulate((sorted_p * len(sorted_p) / ranks)[::-1])[::-1]
    adjusted[sorted_indices] = np.clip(sorted_adjusted, 0.0, 1.0)
    return pd.Series(adjusted, index=p_values.index, dtype=float)


def exact_mcnemar_p_value(increased_only: int, decreased_only: int) -> float:
    """Test directional asymmetry among discordant paired studies."""
    discordant = increased_only + decreased_only
    if discordant == 0:
        return 1.0
    return float(
        stats.binomtest(
            increased_only,
            n=discordant,
            p=0.5,
            alternative="two-sided",
        ).pvalue
    )


def safe_wilcoxon(values: pd.Series) -> float:
    """Return a two-sided paired Wilcoxon p-value, including all-zero inputs."""
    clean = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if len(clean) == 0 or np.allclose(clean, 0.0):
        return 1.0
    return float(stats.wilcoxon(clean, alternative="two-sided", method="auto").pvalue)


def evidence_tier(q_value: object) -> str:
    """Classify one adjusted p-value without calling exploratory results significant."""
    numeric = pd.to_numeric(pd.Series([q_value]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return "ineligible"
    if float(numeric) <= Q_VALUE_CUTOFF:
        return "FDR-significant"
    if float(numeric) <= EXPLORATORY_Q_CUTOFF:
        return "exploratory"
    return "other tested"


def load_traits() -> pd.DataFrame:
    """Load and validate the six shared traits in their display order."""
    traits = read_tsv(TRAIT_SUMMARY_PATH)[["trait_id", "trait_label"]].copy()
    if traits["trait_id"].duplicated().any() or len(traits) != 6:
        raise ValueError("Expected six unique shared traits")
    traits["trait_id"] = traits["trait_id"].map(clean_text)
    traits["trait_label"] = traits["trait_label"].map(clean_text)
    traits["trait_order"] = np.arange(len(traits))
    return traits


def load_metadata() -> pd.DataFrame:
    """Load profile metadata and retain normalized direction labels."""
    metadata = read_tsv(PROFILE_METADATA_PATH, index_col=0)
    metadata = metadata.copy()
    metadata.index = metadata.index.astype(str)
    metadata["profile_id"] = metadata.index
    text_columns = [
        "study_id",
        "experiment_id",
        "inc_dec",
        "condition",
        "condition_group",
        "body_site_group",
        "host_species",
    ]
    for column in text_columns:
        metadata[column] = metadata[column].map(clean_text)
    metadata = metadata[metadata["inc_dec"].isin(DIRECTION_LABELS)].copy()
    metadata["direction"] = metadata["inc_dec"].map(DIRECTION_LABELS)
    return metadata


def build_condition_eligibility(metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Identify human fecal condition-study pairs containing both directions."""
    feces = metadata[
        metadata["body_site_group"].eq(FECES) & metadata["host_species"].eq(HUMAN_HOST)
    ].copy()
    direction_counts = (
        feces.groupby(["condition", "study_id"], as_index=False)["inc_dec"]
        .nunique()
        .rename(columns={"inc_dec": "direction_count"})
    )
    paired = direction_counts[direction_counts["direction_count"].eq(2)][
        ["condition", "study_id"]
    ].copy()
    eligibility = (
        paired.groupby("condition", as_index=False)
        .agg(paired_study_count=("study_id", "nunique"))
        .sort_values(["paired_study_count", "condition"], ascending=[False, True])
    )
    eligibility["passes_min_paired_studies"] = eligibility["paired_study_count"].ge(
        MIN_PAIRED_STUDIES
    )
    return feces, eligibility


def load_annotated_taxa(traits: pd.DataFrame) -> pd.DataFrame:
    """Load exact genus/species taxa with at least one of the six annotations."""
    trait_ids = traits["trait_id"].tolist()
    matrix = read_tsv(
        TAXA_TRAIT_MATRIX_PATH,
        usecols=["taxon_id", "taxon_name", "taxonomic_rank", *trait_ids],
    )
    for column in ["taxon_id", "taxon_name", "taxonomic_rank"]:
        matrix[column] = matrix[column].map(clean_text)
    matrix = matrix[matrix["taxonomic_rank"].isin(ALLOWED_TAXON_RANKS)].copy()
    annotations = matrix[trait_ids].apply(pd.to_numeric, errors="coerce").fillna(0)
    matrix[trait_ids] = annotations
    return matrix[annotations.gt(0).any(axis=1)].copy()


def build_condition_taxon_tests(
    feces_metadata: pd.DataFrame,
    eligibility: pd.DataFrame,
    annotated_taxa: pd.DataFrame,
) -> pd.DataFrame:
    """Test exact-taxon increased/decreased reporting in paired fecal studies."""
    eligible_conditions = set(
        eligibility.loc[eligibility["passes_min_paired_studies"], "condition"]
    )
    eligible_pairs = (
        feces_metadata[feces_metadata["condition"].isin(eligible_conditions)]
        .groupby(["condition", "study_id"])["inc_dec"]
        .nunique()
        .loc[lambda values: values.eq(2)]
        .reset_index()[["condition", "study_id"]]
    )
    profile_lookup = feces_metadata[
        ["study_id", "experiment_id", "direction", "condition"]
    ].drop_duplicates()
    from mtp.population_tables import is_canonical_derived_path, orient_and_filter_taxa_df

    taxa = read_tsv(
        TAXA_PATH,
        usecols=["study_id", "experiment_id", "direction", "taxon_id"],
    )
    if is_canonical_derived_path(TAXA_PATH, "bugsigdb_taxa.tsv"):
        taxa = orient_and_filter_taxa_df(taxa, filter_to_selection=True)
    for column in ["study_id", "experiment_id", "direction", "taxon_id"]:
        taxa[column] = taxa[column].map(clean_text)
    occurrences = taxa.merge(
        profile_lookup,
        on=["study_id", "experiment_id", "direction"],
        how="inner",
        validate="many_to_many",
    )
    occurrences = occurrences.merge(
        eligible_pairs,
        on=["condition", "study_id"],
        how="inner",
        validate="many_to_one",
    )
    occurrences = occurrences.merge(
        annotated_taxa[["taxon_id", "taxon_name", "taxonomic_rank"]],
        on="taxon_id",
        how="inner",
        validate="many_to_one",
    ).drop_duplicates(["condition", "study_id", "direction", "taxon_id"])

    condition_studies = {
        condition: set(group["study_id"])
        for condition, group in eligible_pairs.groupby("condition", sort=False)
    }
    records: list[dict[str, object]] = []
    grouping = ["condition", "taxon_id", "taxon_name", "taxonomic_rank"]
    for keys, group in occurrences.groupby(grouping, sort=False):
        condition, taxon_id, taxon_name, taxonomic_rank = keys
        increased_studies = set(group.loc[group["direction"].eq("increased"), "study_id"])
        decreased_studies = set(group.loc[group["direction"].eq("decreased"), "study_id"])
        all_studies = condition_studies[str(condition)]
        increased_only = increased_studies - decreased_studies
        decreased_only = decreased_studies - increased_studies
        both = increased_studies & decreased_studies
        reported = increased_studies | decreased_studies
        passes_support = len(reported) >= MIN_TAXON_STUDIES
        p_value = exact_mcnemar_p_value(len(increased_only), len(decreased_only))
        difference = (len(increased_only) - len(decreased_only)) / len(all_studies)
        direction = "increased" if difference > 0 else "decreased" if difference < 0 else "tie"
        records.append(
            {
                "condition": condition,
                "taxon_id": taxon_id,
                "taxon_name": taxon_name,
                "taxonomic_rank": taxonomic_rank,
                "paired_study_count": len(all_studies),
                "reported_study_count": len(reported),
                "increased_reported_studies": len(increased_studies),
                "decreased_reported_studies": len(decreased_studies),
                "increased_only_studies": len(increased_only),
                "decreased_only_studies": len(decreased_only),
                "both_directions_studies": len(both),
                "neither_direction_studies": len(all_studies - reported),
                "discordant_study_count": len(increased_only) + len(decreased_only),
                "directional_difference_fraction": difference,
                "direction": direction,
                "p_value": p_value,
                "passes_min_taxon_studies": passes_support,
                "test": "exact McNemar/binomial test on discordant paired studies",
            }
        )
    tests = pd.DataFrame.from_records(records)
    if tests.empty:
        raise ValueError("No annotated genus/species taxa occur in eligible fecal conditions")
    tests["q_value"] = np.nan
    support = tests["passes_min_taxon_studies"].astype(bool)
    tests.loc[support, "q_value"] = benjamini_hochberg(tests.loc[support, "p_value"])
    tests["selected_after_bh"] = support & tests["q_value"].le(Q_VALUE_CUTOFF)
    return tests.sort_values(
        ["selected_after_bh", "q_value", "condition", "taxonomic_rank", "taxon_name"],
        ascending=[False, True, True, True, True],
        na_position="last",
    ).reset_index(drop=True)


def build_zscore_matrix(metadata: pd.DataFrame, traits: pd.DataFrame) -> pd.DataFrame:
    """Recreate the globally normalized six-trait profile z-scores."""
    trait_ids = traits["trait_id"].tolist()
    values = load_retained_profiles(PROFILE_MATRIX_PATH, metadata)
    composition = trait_fractions(values).loc[:, trait_ids]
    scaled = (composition * 1000).apply(lambda column: column.map(math.log1p))
    zscores = (scaled - scaled.mean(axis=0)) / scaled.std(axis=0, ddof=0).replace(0, float("nan"))
    return zscores.fillna(0.0)


def build_condition_trait_tests(
    metadata: pd.DataFrame,
    zscores: pd.DataFrame,
    traits: pd.DataFrame,
    eligibility: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Test paired study-level increased-minus-decreased trait contrasts."""
    trait_ids = traits["trait_id"].tolist()
    trait_names = dict(zip(traits["trait_id"], traits["trait_label"], strict=True))
    eligible_conditions = set(
        eligibility.loc[eligibility["passes_min_paired_studies"], "condition"]
    )
    profiles = zscores[trait_ids].join(
        metadata[
            [
                "study_id",
                "inc_dec",
                "condition",
                "body_site_group",
                "host_species",
            ]
        ]
    )
    profiles = profiles[
        profiles["condition"].isin(eligible_conditions)
        & profiles["body_site_group"].eq(FECES)
        & profiles["host_species"].eq(HUMAN_HOST)
    ].copy()
    contrast_records: list[dict[str, object]] = []
    for (condition, study_id), group in profiles.groupby(["condition", "study_id"], sort=False):
        direction_means = group.groupby("inc_dec")[trait_ids].mean()
        if not {"inc", "dec"}.issubset(direction_means.index):
            continue
        for trait_id in trait_ids:
            contrast_records.append(
                {
                    "condition": condition,
                    "study_id": study_id,
                    "trait_id": trait_id,
                    "trait": trait_names[trait_id],
                    "contrast": float(
                        direction_means.loc["inc", trait_id] - direction_means.loc["dec", trait_id]
                    ),
                    "contrast_definition": "mean z increased - mean z decreased",
                }
            )
    study_contrasts = pd.DataFrame.from_records(contrast_records)
    records: list[dict[str, object]] = []
    for (condition, trait_id, trait), group in study_contrasts.groupby(
        ["condition", "trait_id", "trait"], sort=False
    ):
        values = pd.to_numeric(group["contrast"], errors="coerce").dropna()
        mean = float(values.mean())
        records.append(
            {
                "condition": condition,
                "trait_id": trait_id,
                "trait": trait,
                "paired_study_count": int(group["study_id"].nunique()),
                "mean_contrast": mean,
                "median_contrast": float(values.median()),
                "positive_study_count": int(values.gt(0).sum()),
                "negative_study_count": int(values.lt(0).sum()),
                "zero_study_count": int(values.eq(0).sum()),
                "direction": "increased" if mean > 0 else "decreased" if mean < 0 else "tie",
                "p_value": safe_wilcoxon(values),
                "test": "two-sided Wilcoxon signed-rank test vs zero",
            }
        )
    tests = pd.DataFrame.from_records(records)
    if tests.empty:
        raise ValueError("No eligible paired fecal condition-trait contrasts")
    tests["passes_min_paired_studies"] = tests["paired_study_count"].ge(MIN_PAIRED_STUDIES)
    tests["q_value"] = np.nan
    support = tests["passes_min_paired_studies"].astype(bool)
    tests.loc[support, "q_value"] = benjamini_hochberg(tests.loc[support, "p_value"])
    tests["selected_after_bh"] = support & tests["q_value"].le(Q_VALUE_CUTOFF)
    tests = tests.sort_values(
        ["selected_after_bh", "q_value", "condition", "trait"],
        ascending=[False, True, True, True],
        na_position="last",
    ).reset_index(drop=True)
    return tests, study_contrasts


NODE_COLUMNS = ["node_id", "label", "layer", "node_type", "taxonomic_rank"]
EDGE_COLUMNS = [
    "source",
    "target",
    "source_layer",
    "target_layer",
    "edge_class",
    "evidence_type",
    "direction",
    "effect",
    "p_value",
    "q_value",
    "study_support",
    "annotation_count",
    "statistically_significant",
    "line_style",
]


def build_network_tables(
    taxon_tests: pd.DataFrame,
    trait_tests: pd.DataFrame,
    annotated_taxa: pd.DataFrame,
    traits: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build nodes and inferential/provenance edges from FDR-selected results."""
    selected_taxa = taxon_tests[taxon_tests["selected_after_bh"]].copy()
    selected_traits = trait_tests[trait_tests["selected_after_bh"]].copy()
    displayed_conditions = sorted(
        set(selected_taxa["condition"]) | set(selected_traits["condition"]),
        key=lambda value: (value != PD_CONDITION, value),
    )
    if not displayed_conditions:
        return (
            pd.DataFrame(columns=NODE_COLUMNS),
            pd.DataFrame(columns=EDGE_COLUMNS),
            selected_taxa,
        )

    trait_names = dict(zip(traits["trait_id"], traits["trait_label"], strict=True))
    trait_ids = traits["trait_id"].tolist()
    selected_taxon_ids = set(selected_taxa["taxon_id"])
    annotation_rows = annotated_taxa[annotated_taxa["taxon_id"].isin(selected_taxon_ids)].copy()

    nodes: list[dict[str, object]] = []
    edges: list[dict[str, object]] = []
    for condition in displayed_conditions:
        nodes.append(
            {
                "node_id": f"condition::{condition}",
                "label": condition,
                "layer": "condition",
                "node_type": "human fecal condition with an FDR-selected association",
                "taxonomic_rank": pd.NA,
            }
        )
    for row in selected_taxa.itertuples(index=False):
        edges.append(
            {
                "source": f"condition::{row.condition}",
                "target": f"taxon::{row.taxon_id}",
                "source_layer": "condition",
                "target_layer": row.taxonomic_rank,
                "edge_class": "FDR-significant condition-taxon direction",
                "evidence_type": row.test,
                "direction": row.direction,
                "effect": row.directional_difference_fraction,
                "p_value": row.p_value,
                "q_value": row.q_value,
                "study_support": row.reported_study_count,
                "annotation_count": pd.NA,
                "statistically_significant": True,
                "line_style": "solid",
            }
        )
    unique_taxa = selected_taxa[["taxon_id", "taxon_name", "taxonomic_rank"]].drop_duplicates(
        "taxon_id"
    )
    for row in unique_taxa.itertuples(index=False):
        nodes.append(
            {
                "node_id": f"taxon::{row.taxon_id}",
                "label": row.taxon_name,
                "layer": row.taxonomic_rank,
                "node_type": "exact BugSigDB taxon with an FDR-selected directional association",
                "taxonomic_rank": row.taxonomic_rank,
            }
        )
    for row in selected_traits.itertuples(index=False):
        edges.append(
            {
                "source": f"condition::{row.condition}",
                "target": f"trait::{row.trait_id}",
                "source_layer": "condition",
                "target_layer": "trait",
                "edge_class": "FDR-significant condition-trait contrast",
                "evidence_type": row.test,
                "direction": row.direction,
                "effect": row.mean_contrast,
                "p_value": row.p_value,
                "q_value": row.q_value,
                "study_support": row.paired_study_count,
                "annotation_count": pd.NA,
                "statistically_significant": True,
                "line_style": "solid",
            }
        )
    annotated_trait_ids: set[str] = set()
    for _, row in annotation_rows.iterrows():
        for trait_id in trait_ids:
            annotation_count = int(pd.to_numeric(row[trait_id], errors="coerce") or 0)
            if annotation_count <= 0:
                continue
            annotated_trait_ids.add(trait_id)
            edges.append(
                {
                    "source": f"taxon::{row['taxon_id']}",
                    "target": f"trait::{trait_id}",
                    "source_layer": row["taxonomic_rank"],
                    "target_layer": "trait",
                    "edge_class": "KG-Microbe annotation provenance",
                    "evidence_type": "KG-Microbe taxon-trait annotation; no inferential test",
                    "direction": pd.NA,
                    "effect": pd.NA,
                    "p_value": pd.NA,
                    "q_value": pd.NA,
                    "study_support": pd.NA,
                    "annotation_count": annotation_count,
                    "statistically_significant": False,
                    "line_style": "dashed",
                }
            )
    displayed_trait_ids = set(selected_traits["trait_id"]) | annotated_trait_ids
    for trait_id in trait_ids:
        if trait_id not in displayed_trait_ids:
            continue
        nodes.append(
            {
                "node_id": f"trait::{trait_id}",
                "label": trait_names[trait_id],
                "layer": "trait",
                "node_type": "six-shared-trait result or KG-Microbe annotation target",
                "taxonomic_rank": pd.NA,
            }
        )
    nodes_frame = pd.DataFrame.from_records(nodes, columns=NODE_COLUMNS).drop_duplicates("node_id")
    edges_frame = pd.DataFrame.from_records(edges, columns=EDGE_COLUMNS)
    return nodes_frame, edges_frame, selected_taxa


def augment_condition_eligibility(
    eligibility: pd.DataFrame,
    taxon_tests: pd.DataFrame,
    trait_tests: pd.DataFrame,
) -> pd.DataFrame:
    """Record why each condition is or is not displayed."""
    taxon_counts = (
        taxon_tests[taxon_tests["selected_after_bh"]].groupby("condition")["taxon_id"].nunique()
    )
    trait_counts = (
        trait_tests[trait_tests["selected_after_bh"]].groupby("condition")["trait_id"].nunique()
    )
    output = eligibility.copy()
    output["fdr_selected_taxon_count"] = output["condition"].map(taxon_counts).fillna(0).astype(int)
    output["fdr_selected_trait_count"] = output["condition"].map(trait_counts).fillna(0).astype(int)
    output["displayed_in_network"] = (
        output[["fdr_selected_taxon_count", "fdr_selected_trait_count"]].sum(axis=1).gt(0)
    )
    output["exclusion_reason"] = pd.NA
    output.loc[~output["passes_min_paired_studies"], "exclusion_reason"] = (
        f"fewer than {MIN_PAIRED_STUDIES} paired human fecal studies"
    )
    no_selected = output["passes_min_paired_studies"] & ~output["displayed_in_network"]
    output.loc[no_selected, "exclusion_reason"] = "no association passes global BH q<=0.05"
    return output.sort_values(
        ["displayed_in_network", "paired_study_count", "condition"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def _even_positions(labels: list[str], top: float = 0.9, bottom: float = 0.08) -> dict[str, float]:
    """Return evenly spaced vertical axes coordinates."""
    if not labels:
        return {}
    if len(labels) == 1:
        return {labels[0]: (top + bottom) / 2}
    values = np.linspace(top, bottom, len(labels))
    return dict(zip(labels, values, strict=True))


def _curved_edge(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str,
    linewidth: float,
    linestyle: str | tuple[int, tuple[int, ...]],
    alpha: float,
    bend: float = 0.0,
) -> None:
    """Draw a cubic Bezier network edge in axes coordinates."""
    x0, y0 = start
    x1, y1 = end
    control_x0 = x0 + (x1 - x0) * 0.42
    control_x1 = x0 + (x1 - x0) * 0.58
    path = MplPath(
        [(x0, y0), (control_x0, y0 + bend), (control_x1, y1 + bend), (x1, y1)],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
    )
    ax.add_patch(
        PathPatch(
            path,
            transform=ax.transAxes,
            facecolor="none",
            edgecolor=color,
            linewidth=linewidth,
            linestyle=linestyle,
            alpha=alpha,
            zorder=1,
        )
    )


def _draw_node(
    ax: plt.Axes,
    x: float,
    y: float,
    label: str,
    *,
    color: str,
    width: float,
    height: float,
    text_color: str = "white",
    fontsize: float = 8.0,
) -> None:
    """Draw one rounded network node."""
    patch = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width,
        height,
        boxstyle="round,pad=0.006,rounding_size=0.008",
        transform=ax.transAxes,
        facecolor=color,
        edgecolor="white",
        linewidth=0.9,
        zorder=3,
    )
    ax.add_patch(patch)
    ax.text(
        x,
        y,
        wrap_label(label, width=18),
        transform=ax.transAxes,
        ha="center",
        va="center",
        color=text_color,
        fontsize=fontsize,
        fontweight="bold",
        zorder=4,
    )


def plot_significant_network(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    taxon_tests: pd.DataFrame,
    trait_tests: pd.DataFrame,
    traits: pd.DataFrame,
    output_path: Path,
) -> None:
    """Render the significant network plus the significant condition-trait matrix."""
    taxon_edges = edges[edges["edge_class"].eq("FDR-significant condition-taxon direction")]
    annotation_edges = edges[edges["edge_class"].eq("KG-Microbe annotation provenance")]
    selected_trait_tests = trait_tests[trait_tests["selected_after_bh"]].copy()
    selected_taxon_count = len(taxon_edges)
    selected_trait_count = len(selected_trait_tests)
    eligible_condition_count = taxon_tests["condition"].nunique()
    eligible_taxon_test_count = int(taxon_tests["passes_min_taxon_studies"].sum())

    condition_ids = taxon_edges["source"].drop_duplicates().tolist()
    condition_metrics = (
        taxon_edges.groupby("source").size().reindex(condition_ids, fill_value=0).to_dict()
    )
    node_index = nodes.set_index("node_id")
    condition_ids = sorted(
        condition_ids,
        key=lambda node_id: (
            node_index.loc[node_id, "label"] != PD_CONDITION,
            -int(condition_metrics.get(node_id, 0)),
            str(node_index.loc[node_id, "label"]),
        ),
    )
    genus_ids = nodes.loc[nodes["layer"].eq("genus"), "node_id"].tolist()
    species_ids = nodes.loc[nodes["layer"].eq("species"), "node_id"].tolist()
    taxon_metrics = taxon_edges.groupby("target").size().to_dict()
    name_lookup = node_index["label"].to_dict()
    genus_ids.sort(key=lambda node_id: (-int(taxon_metrics.get(node_id, 0)), name_lookup[node_id]))
    species_ids.sort(
        key=lambda node_id: (-int(taxon_metrics.get(node_id, 0)), name_lookup[node_id])
    )
    trait_order = traits.set_index("trait_id")["trait_order"].to_dict()
    trait_ids = annotation_edges["target"].drop_duplicates().tolist()
    trait_ids.sort(key=lambda node_id: int(trait_order[node_id.removeprefix("trait::")]))

    matrix_condition_labels = sorted(
        selected_trait_tests["condition"].unique(),
        key=lambda condition: (condition != PD_CONDITION, condition),
    )

    row_pressure = max(
        len(condition_ids),
        len(genus_ids) + len(species_ids),
        len(trait_ids),
        len(matrix_condition_labels),
        6,
    )
    figure_height = max(7.5, min(24.0, 0.62 * row_pressure + 4.0))
    fig = plt.figure(figsize=(18, figure_height), facecolor="white")
    grid = fig.add_gridspec(1, 2, width_ratios=[1.75, 1.0], wspace=0.18)
    ax_network = fig.add_subplot(grid[0, 0])
    ax_matrix = fig.add_subplot(grid[0, 1])
    ax_network.set_axis_off()

    x_positions = {"condition": 0.08, "genus": 0.45, "species": 0.64, "trait": 0.93}
    positions: dict[str, tuple[float, float]] = {}
    for node_id, y in _even_positions(condition_ids).items():
        positions[node_id] = (x_positions["condition"], y)
    for node_id, y in _even_positions(genus_ids, top=0.88, bottom=0.53).items():
        positions[node_id] = (x_positions["genus"], y)
    for node_id, y in _even_positions(species_ids, top=0.45, bottom=0.1).items():
        positions[node_id] = (x_positions["species"], y)
    if genus_ids and not species_ids:
        for node_id, y in _even_positions(genus_ids).items():
            positions[node_id] = (x_positions["genus"], y)
    if species_ids and not genus_ids:
        for node_id, y in _even_positions(species_ids).items():
            positions[node_id] = (x_positions["species"], y)
    for node_id, y in _even_positions(trait_ids).items():
        positions[node_id] = (x_positions["trait"], y)

    for row in taxon_edges.itertuples(index=False):
        q_value = max(float(row.q_value), np.finfo(float).tiny)
        _curved_edge(
            ax_network,
            positions[row.source],
            positions[row.target],
            color=INCREASED_COLOR if row.direction == "increased" else DECREASED_COLOR,
            linewidth=min(6.0, 1.2 + -math.log10(q_value)),
            linestyle="solid",
            alpha=0.72,
        )
    for row in annotation_edges.itertuples(index=False):
        _curved_edge(
            ax_network,
            positions[row.source],
            positions[row.target],
            color=TRAIT_COLORS.get(name_lookup[row.target], ANNOTATION_COLOR),
            linewidth=min(4.0, 0.8 + math.log1p(float(row.annotation_count))),
            linestyle=(0, (3, 3)),
            alpha=0.38,
        )

    for node_id in condition_ids:
        label = name_lookup[node_id]
        _draw_node(
            ax_network,
            *positions[node_id],
            label,
            color=PD_COLOR if label == PD_CONDITION else SLATE_700,
            width=0.2,
            height=0.048,
        )
    for node_id in genus_ids:
        _draw_node(
            ax_network,
            *positions[node_id],
            f"{name_lookup[node_id]} [G]",
            color="#475569",
            width=0.17,
            height=0.043,
            fontsize=7.4,
        )
    for node_id in species_ids:
        _draw_node(
            ax_network,
            *positions[node_id],
            f"{name_lookup[node_id]} [S]",
            color="#64748b",
            width=0.18,
            height=0.043,
            fontsize=7.2,
        )
    for node_id in trait_ids:
        label = name_lookup[node_id]
        _draw_node(
            ax_network,
            *positions[node_id],
            label,
            color=TRAIT_COLORS.get(label, SLATE_500),
            width=0.19,
            height=0.052,
            fontsize=7.5,
        )

    if selected_taxon_count:
        for x, title in [
            (x_positions["condition"], "Condition"),
            (x_positions["genus"], "FDR-selected genus"),
            (
                x_positions["species"],
                "FDR-selected species" if species_ids else "Species: none pass q<=0.05",
            ),
            (x_positions["trait"], "KG-Microbe trait"),
        ]:
            ax_network.text(
                x,
                0.985,
                title,
                transform=ax_network.transAxes,
                ha="center",
                va="top",
                fontsize=10,
                fontweight="bold",
                color=SLATE_900 if "none" not in title else SLATE_500,
            )
        ax_network.legend(
            handles=[
                Line2D(
                    [0],
                    [0],
                    color=INCREASED_COLOR,
                    lw=3,
                    label="solid: taxon favors increased; q<=0.05",
                ),
                Line2D(
                    [0],
                    [0],
                    color=DECREASED_COLOR,
                    lw=3,
                    label="solid: taxon favors decreased; q<=0.05",
                ),
                Line2D(
                    [0],
                    [0],
                    color=ANNOTATION_COLOR,
                    lw=2,
                    linestyle=(0, (3, 3)),
                    label="dashed: KG-Microbe annotation (not a significance test)",
                ),
            ],
            loc="lower left",
            frameon=False,
            fontsize=8.5,
        )
        panel_a_title = "A  Significant condition-taxon directions and trait provenance"
    else:
        stages = [
            (
                0.145,
                f"{eligible_condition_count}\neligible conditions",
                "paired human fecal\nstudies",
                SLATE_700,
            ),
            (
                0.50,
                f"{eligible_taxon_test_count}\neligible tests",
                "genus and species\nkept separate",
                "#475569",
            ),
            (
                0.855,
                "0 FDR-selected\ntaxon results",
                "global BH q<=0.05",
                "#991b1b",
            ),
        ]
        for x, heading, detail, color in stages:
            box = FancyBboxPatch(
                (x - 0.115, 0.48),
                0.23,
                0.19,
                boxstyle="round,pad=0.012,rounding_size=0.015",
                transform=ax_network.transAxes,
                facecolor="white",
                edgecolor=color,
                linewidth=2.2,
            )
            ax_network.add_patch(box)
            ax_network.text(
                x,
                0.59,
                heading,
                transform=ax_network.transAxes,
                ha="center",
                va="center",
                fontsize=12,
                fontweight="bold",
                color=color,
            )
            ax_network.text(
                x,
                0.515,
                detail,
                transform=ax_network.transAxes,
                ha="center",
                va="center",
                fontsize=8.5,
                color=SLATE_700,
            )
        for start, end in [(0.26, 0.385), (0.615, 0.74)]:
            ax_network.annotate(
                "",
                xy=(end, 0.575),
                xytext=(start, 0.575),
                xycoords=ax_network.transAxes,
                arrowprops={"arrowstyle": "-|>", "color": SLATE_300, "lw": 2},
            )
        ax_network.text(
            0.5,
            0.37,
            (
                "No genus or species node is displayed because none passed the "
                "pre-specified FDR threshold.\nThis is a null screen, not evidence "
                "that biological associations are absent."
            ),
            transform=ax_network.transAxes,
            ha="center",
            va="top",
            fontsize=10,
            color=SLATE_700,
            linespacing=1.5,
        )
        ax_network.text(
            0.5,
            0.20,
            "KG-Microbe provenance links are omitted because no taxon was selected.",
            transform=ax_network.transAxes,
            ha="center",
            va="top",
            fontsize=9,
            color=SLATE_500,
        )
        panel_a_title = "A  Condition-taxon screen: no association survives FDR"
    ax_network.set_title(
        panel_a_title,
        loc="left",
        fontsize=12,
        fontweight="bold",
        color=SLATE_900,
        pad=12,
        y=0.70 if not selected_taxon_count else 1.0,
    )

    condition_labels = matrix_condition_labels
    trait_labels = traits["trait_label"].tolist()
    effects = selected_trait_tests.pivot(index="condition", columns="trait", values="mean_contrast")
    q_values = selected_trait_tests.pivot(index="condition", columns="trait", values="q_value")
    counts = selected_trait_tests.pivot(
        index="condition", columns="trait", values="paired_study_count"
    )
    effects = effects.reindex(index=condition_labels, columns=trait_labels)
    q_values = q_values.reindex(index=condition_labels, columns=trait_labels)
    counts = counts.reindex(index=condition_labels, columns=trait_labels)
    finite_effects = effects.to_numpy(dtype=float)
    finite_values = finite_effects[np.isfinite(finite_effects)]
    limit = max(float(np.max(np.abs(finite_values))) if len(finite_values) else 1.0, 0.1)
    color_norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    cmap = plt.get_cmap("RdBu_r")
    for row_index, condition in enumerate(condition_labels):
        for column_index, trait in enumerate(trait_labels):
            value = effects.loc[condition, trait]
            facecolor = "white" if pd.isna(value) else cmap(color_norm(float(value)))
            ax_matrix.add_patch(
                Rectangle(
                    (column_index - 0.5, row_index - 0.5),
                    1,
                    1,
                    facecolor=facecolor,
                    edgecolor=SLATE_200,
                    linewidth=0.8,
                )
            )
            if pd.isna(value):
                continue
            q_value = float(q_values.loc[condition, trait])
            study_count = int(counts.loc[condition, trait])
            text_color = "white" if abs(float(value)) > 0.62 * limit else SLATE_900
            ax_matrix.text(
                column_index,
                row_index,
                f"{float(value):+.2f}\nq={q_value:.3g}; n={study_count}",
                ha="center",
                va="center",
                fontsize=7.0,
                color=text_color,
                fontweight="bold",
            )
    ax_matrix.set_xlim(-0.5, len(trait_labels) - 0.5)
    ax_matrix.set_ylim(len(condition_labels) - 0.5, -0.5)
    ax_matrix.set_xticks(
        range(len(trait_labels)), [wrap_label(label, 14) for label in trait_labels]
    )
    ax_matrix.set_yticks(
        range(len(condition_labels)), [wrap_label(label, 18) for label in condition_labels]
    )
    ax_matrix.set_aspect("equal", adjustable="box")
    ax_matrix.set_anchor("N")
    ax_matrix.tick_params(axis="x", labelrotation=45, labelsize=8)
    ax_matrix.tick_params(axis="y", labelsize=8)
    for spine in ax_matrix.spines.values():
        spine.set_visible(False)
    ax_matrix.set_title(
        "B  Significant paired condition-trait contrasts",
        loc="left",
        fontsize=12,
        fontweight="bold",
        color=SLATE_900,
        pad=12,
    )
    ax_matrix.text(
        1.0,
        1.02,
        "White cells: tested, but did not pass global BH",
        transform=ax_matrix.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color=SLATE_500,
    )
    scalar = plt.cm.ScalarMappable(norm=Normalize(vmin=-limit, vmax=limit), cmap=cmap)
    colorbar = fig.colorbar(
        scalar,
        ax=ax_matrix,
        orientation="horizontal",
        fraction=0.07,
        pad=0.45,
        aspect=28,
    )
    colorbar.set_label("Mean study contrast (increased - decreased z)", fontsize=8)
    colorbar.ax.tick_params(labelsize=7)

    fig.suptitle(
        "Human fecal ecological-reporter network: globally FDR-filtered associations",
        x=0.03,
        y=0.985,
        ha="left",
        fontsize=18,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.03,
        0.925,
        (
            f"All conditions with >={MIN_PAIRED_STUDIES} paired studies and all annotated "
            f"genus/species taxa with >={MIN_TAXON_STUDIES} reporting studies were "
            "tested before selection.\nSolid links/cells pass global "
            f"BH q<={Q_VALUE_CUTOFF:.2f} ({selected_taxon_count} taxon and "
            f"{selected_trait_count} trait results); "
            "dashed links are annotation provenance only."
        ),
        ha="left",
        va="top",
        fontsize=9.2,
        color=SLATE_700,
    )
    fig.subplots_adjust(top=0.87, bottom=0.055, left=0.025, right=0.975)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    svg_path = output_path.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    svg_lines = svg_path.read_text(encoding="utf-8").splitlines()
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_lines) + "\n",
        encoding="utf-8",
    )


def _plot_ready_tests(frame: pd.DataFrame, effect_column: str) -> pd.DataFrame:
    """Return tested rows with numeric effects, q-values, tiers, and plot heights."""
    output = frame.copy()
    output["q_value"] = pd.to_numeric(output["q_value"], errors="coerce")
    output[effect_column] = pd.to_numeric(output[effect_column], errors="coerce")
    output = output.dropna(subset=["q_value", effect_column]).copy()
    output["evidence_tier"] = output["q_value"].map(evidence_tier)
    output["negative_log10_q"] = -np.log10(output["q_value"].clip(lower=np.finfo(float).tiny))
    return output


def _scatter_evidence_tiers(
    ax: plt.Axes,
    frame: pd.DataFrame,
    *,
    effect_column: str,
    support_column: str,
    marker_column: str | None = None,
) -> None:
    """Plot all tests while reserving strong styling for FDR-significant results."""
    marker_groups: list[tuple[str | None, str]] = [(None, "o")]
    if marker_column is not None:
        marker_groups = [("genus", "o"), ("species", "^")]
    tier_order = ["other tested", "exploratory", "FDR-significant"]
    for tier in tier_order:
        for marker_value, marker in marker_groups:
            subset = frame[frame["evidence_tier"].eq(tier)]
            if marker_value is not None and marker_column is not None:
                subset = subset[subset[marker_column].eq(marker_value)]
            if subset.empty:
                continue
            support = pd.to_numeric(subset[support_column], errors="coerce").fillna(0)
            sizes = np.clip(18.0 + 1.6 * support.to_numpy(dtype=float), 20.0, 90.0)
            if tier == "FDR-significant":
                facecolors = np.where(
                    subset[effect_column].to_numpy(dtype=float) >= 0,
                    INCREASED_COLOR,
                    DECREASED_COLOR,
                )
                edgecolors = SLATE_900
                linewidths = 0.8
                alpha = 0.95
                zorder = 4
            elif tier == "exploratory":
                facecolors = "white"
                edgecolors = "#d97706"
                linewidths = 1.6
                alpha = 0.95
                zorder = 3
            else:
                facecolors = SLATE_300
                edgecolors = "white"
                linewidths = 0.35
                alpha = 0.52
                zorder = 2
            ax.scatter(
                subset[effect_column],
                subset["negative_log10_q"],
                s=sizes,
                marker=marker,
                facecolors=facecolors,
                edgecolors=edgecolors,
                linewidths=linewidths,
                alpha=alpha,
                zorder=zorder,
            )


def _format_probability(value: object) -> str:
    """Format a p- or q-value compactly for figure annotations."""
    numeric = float(value)
    return f"{numeric:.3g}"


def _label_priority_results(
    ax: plt.Axes,
    frame: pd.DataFrame,
    *,
    effect_column: str,
    second_label_column: str,
    rank_column: str | None = None,
    limit: int = 8,
) -> None:
    """Label significant and exploratory results without labeling gray background tests."""
    priority = frame[frame["evidence_tier"].isin(["FDR-significant", "exploratory"])].nsmallest(
        limit, "q_value"
    )
    offsets = [(8, 18), (8, -28), (-8, 20), (-8, -30), (10, 34), (10, -42)]
    for index, row in enumerate(priority.itertuples(index=False)):
        label = f"{row.condition}\n{getattr(row, second_label_column)}"
        if rank_column is not None:
            rank = str(getattr(row, rank_column))
            label += f" [{rank[0].upper()}]"
        label += f"\nq={_format_probability(row.q_value)}"
        offset = offsets[index % len(offsets)]
        ax.annotate(
            label,
            xy=(getattr(row, effect_column), row.negative_log10_q),
            xytext=offset,
            textcoords="offset points",
            ha="left" if offset[0] >= 0 else "right",
            va="bottom" if offset[1] >= 0 else "top",
            fontsize=7.5,
            color=SLATE_700,
            arrowprops={"arrowstyle": "-", "color": SLATE_500, "lw": 0.7},
            zorder=5,
        )


def _style_evidence_axis(
    ax: plt.Axes,
    frame: pd.DataFrame,
    *,
    effect_column: str,
    x_label: str,
) -> None:
    """Add shared thresholds and axis styling to an evidence landscape."""
    significant_height = -math.log10(Q_VALUE_CUTOFF)
    exploratory_height = -math.log10(EXPLORATORY_Q_CUTOFF)
    ax.axhline(
        significant_height,
        color=SLATE_900,
        linestyle=(0, (5, 3)),
        linewidth=1.2,
        zorder=1,
    )
    ax.axhline(
        exploratory_height,
        color="#d97706",
        linestyle=(0, (2, 3)),
        linewidth=1.1,
        zorder=1,
    )
    ax.axvline(0.0, color=SLATE_500, linewidth=0.8, zorder=1)
    ax.text(
        0.99,
        significant_height,
        " q=0.05",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=7.5,
        color=SLATE_700,
    )
    ax.text(
        0.99,
        exploratory_height,
        " q=0.20 exploratory",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="#b45309",
    )
    max_effect = max(float(frame[effect_column].abs().max()), 0.1)
    max_height = max(float(frame["negative_log10_q"].max()), significant_height)
    ax.set_xlim(-1.15 * max_effect, 1.15 * max_effect)
    ax.set_ylim(-0.03, max_height + 0.42)
    ax.set_xlabel(x_label, fontsize=9)
    ax.set_ylabel("Statistical evidence (-log10 global q)", fontsize=9)
    ax.grid(axis="y", color=SLATE_200, linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)


def plot_fdr_evidence_landscape(
    taxon_tests: pd.DataFrame,
    trait_tests: pd.DataFrame,
    traits: pd.DataFrame,
    output_path: Path,
) -> None:
    """Show all eligible tests with significant, exploratory, and background tiers."""
    taxon_plot = _plot_ready_tests(
        taxon_tests[taxon_tests["passes_min_taxon_studies"]],
        "directional_difference_fraction",
    )
    trait_plot = _plot_ready_tests(trait_tests, "mean_contrast")
    taxon_counts = taxon_plot["evidence_tier"].value_counts()
    trait_counts = trait_plot["evidence_tier"].value_counts()

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(21, 8.2),
        gridspec_kw={"width_ratios": [1.15, 1.15, 1.05]},
        facecolor="white",
    )
    ax_taxon, ax_trait, ax_pd = axes

    _scatter_evidence_tiers(
        ax_taxon,
        taxon_plot,
        effect_column="directional_difference_fraction",
        support_column="reported_study_count",
        marker_column="taxonomic_rank",
    )
    _style_evidence_axis(
        ax_taxon,
        taxon_plot,
        effect_column="directional_difference_fraction",
        x_label="Directional study fraction (increased - decreased)",
    )
    _label_priority_results(
        ax_taxon,
        taxon_plot,
        effect_column="directional_difference_fraction",
        second_label_column="taxon_name",
        rank_column="taxonomic_rank",
    )
    ax_taxon.set_title(
        f"A  Condition-taxon evidence ({len(taxon_plot)} eligible tests)",
        loc="left",
        fontsize=12,
        fontweight="bold",
        color=SLATE_900,
        pad=12,
    )
    ax_taxon.text(
        0.0,
        -0.17,
        (
            f"FDR-significant: {int(taxon_counts.get('FDR-significant', 0))}  |  "
            f"Exploratory: {int(taxon_counts.get('exploratory', 0))}  |  "
            "Point area: reporting studies"
        ),
        transform=ax_taxon.transAxes,
        fontsize=8,
        color=SLATE_700,
    )

    _scatter_evidence_tiers(
        ax_trait,
        trait_plot,
        effect_column="mean_contrast",
        support_column="paired_study_count",
    )
    _style_evidence_axis(
        ax_trait,
        trait_plot,
        effect_column="mean_contrast",
        x_label="Mean study contrast (increased - decreased z)",
    )
    _label_priority_results(
        ax_trait,
        trait_plot,
        effect_column="mean_contrast",
        second_label_column="trait",
    )
    ax_trait.set_title(
        f"B  Condition-trait evidence ({len(trait_plot)} tests)",
        loc="left",
        fontsize=12,
        fontweight="bold",
        color=SLATE_900,
        pad=12,
    )
    ax_trait.text(
        0.0,
        -0.17,
        (
            f"FDR-significant: {int(trait_counts.get('FDR-significant', 0))}  |  "
            f"Exploratory: {int(trait_counts.get('exploratory', 0))}  |  "
            "Point area: paired studies"
        ),
        transform=ax_trait.transAxes,
        fontsize=8,
        color=SLATE_700,
    )

    trait_order = traits.set_index("trait_label")["trait_order"].to_dict()
    pd_rows = trait_plot[trait_plot["condition"].eq(PD_CONDITION)].copy()
    pd_rows["trait_order"] = pd_rows["trait"].map(trait_order)
    pd_rows = pd_rows.sort_values("trait_order", ascending=False).reset_index(drop=True)
    y_positions = np.arange(len(pd_rows))
    for y_position, row in zip(y_positions, pd_rows.itertuples(index=False), strict=True):
        effect = float(row.mean_contrast)
        tier = str(row.evidence_tier)
        if tier == "FDR-significant":
            facecolor = TRAIT_COLORS.get(row.trait, INCREASED_COLOR)
            edgecolor = SLATE_900
            linewidth = 1.0
        elif tier == "exploratory":
            facecolor = "white"
            edgecolor = "#d97706"
            linewidth = 1.8
        else:
            facecolor = SLATE_300
            edgecolor = "white"
            linewidth = 0.7
        ax_pd.hlines(y_position, 0, effect, color=SLATE_300, linewidth=2.0, zorder=1)
        ax_pd.scatter(
            [effect],
            [y_position],
            s=90,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=3,
        )
    max_pd_effect = max(float(pd_rows["mean_contrast"].abs().max()), 0.1)
    label_x = max_pd_effect * 1.08
    for y_position, row in zip(y_positions, pd_rows.itertuples(index=False), strict=True):
        ax_pd.text(
            label_x,
            y_position,
            (
                f"effect={float(row.mean_contrast):+.2f}; "
                f"p={_format_probability(row.p_value)}; "
                f"q={_format_probability(row.q_value)}; "
                f"n={int(row.paired_study_count)}"
            ),
            ha="left",
            va="center",
            fontsize=7.5,
            color=SLATE_700,
        )
    ax_pd.axvline(0.0, color=SLATE_500, linewidth=0.9)
    ax_pd.set_xlim(-1.25 * max_pd_effect, 2.45 * max_pd_effect)
    ax_pd.set_yticks(y_positions, [wrap_label(value, 18) for value in pd_rows["trait"]])
    ax_pd.set_xlabel("PD mean study contrast", fontsize=9)
    ax_pd.set_title(
        "C  Parkinson's disease: all six traits",
        loc="left",
        fontsize=12,
        fontweight="bold",
        color=SLATE_900,
        pad=12,
    )
    ax_pd.grid(axis="x", color=SLATE_200, linewidth=0.7)
    ax_pd.set_axisbelow(True)
    ax_pd.spines[["top", "right", "left"]].set_visible(False)
    ax_pd.tick_params(axis="y", length=0, labelsize=8)
    ax_pd.tick_params(axis="x", labelsize=8)
    ax_pd.text(
        0.0,
        -0.17,
        "p is unadjusted; q is global across all condition-trait tests.",
        transform=ax_pd.transAxes,
        fontsize=8,
        color=SLATE_700,
    )

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=INCREASED_COLOR,
            markeredgecolor=SLATE_900,
            markersize=7,
            label="FDR-significant (q<=0.05)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="white",
            markeredgecolor="#d97706",
            markeredgewidth=1.6,
            markersize=7,
            label="exploratory (0.05<q<=0.20; not significant)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=SLATE_300,
            markeredgecolor="white",
            markersize=7,
            label="other tested result (q>0.20)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=SLATE_300,
            markeredgecolor=SLATE_500,
            markersize=7,
            label="genus",
        ),
        Line2D(
            [0],
            [0],
            marker="^",
            color="none",
            markerfacecolor=SLATE_300,
            markeredgecolor=SLATE_500,
            markersize=7,
            label="species",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.035),
        ncol=5,
        frameon=False,
        fontsize=8.5,
    )
    fig.suptitle(
        "Human fecal ecological-reporter evidence landscape",
        x=0.035,
        y=0.985,
        ha="left",
        fontsize=19,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.035,
        0.935,
        (
            "All eligible tests are shown. Global BH correction is applied separately "
            "to the condition-taxon and condition-trait families.\nOnly filled points "
            "are FDR-significant; hollow amber points are exploratory and are not "
            "claimed as discoveries."
        ),
        ha="left",
        va="top",
        fontsize=9.5,
        color=SLATE_700,
    )
    fig.subplots_adjust(top=0.82, bottom=0.19, left=0.055, right=0.985, wspace=0.34)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    svg_path = output_path.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    svg_lines = svg_path.read_text(encoding="utf-8").splitlines()
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_lines) + "\n",
        encoding="utf-8",
    )


def write_methods_summary(
    eligibility: pd.DataFrame,
    taxon_tests: pd.DataFrame,
    trait_tests: pd.DataFrame,
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    output_path: Path,
) -> None:
    """Write a compact interpretation and methods companion for the evidence figure."""
    eligible_count = int(eligibility["passes_min_paired_studies"].sum())
    displayed_count = int(eligibility["displayed_in_network"].sum())
    eligible_taxon_test_count = int(taxon_tests["passes_min_taxon_studies"].sum())
    selected_taxa = taxon_tests[taxon_tests["selected_after_bh"]]
    selected_traits = trait_tests[trait_tests["selected_after_bh"]]
    exploratory_taxa = taxon_tests[
        taxon_tests["passes_min_taxon_studies"]
        & taxon_tests["q_value"].gt(Q_VALUE_CUTOFF)
        & taxon_tests["q_value"].le(EXPLORATORY_Q_CUTOFF)
    ]
    exploratory_traits = trait_tests[
        trait_tests["q_value"].gt(Q_VALUE_CUTOFF) & trait_tests["q_value"].le(EXPLORATORY_Q_CUTOFF)
    ]
    genus_count = int(
        selected_taxa.loc[selected_taxa["taxonomic_rank"].eq("genus"), "taxon_id"].nunique()
    )
    species_count = int(
        selected_taxa.loc[selected_taxa["taxonomic_rank"].eq("species"), "taxon_id"].nunique()
    )
    annotation_count = int(edges["edge_class"].eq("KG-Microbe annotation provenance").sum())
    lines = [
        "# FDR-tiered condition-taxon-trait evidence landscape",
        "",
        "## Statistical universe",
        "",
        (
            f"- Human fecal conditions with at least `{MIN_PAIRED_STUDIES}` studies "
            "containing both increased and decreased signatures: "
            f"`{eligible_count}`."
        ),
        (
            "- Exact genus/species taxa required reporting in at least "
            f"`{MIN_TAXON_STUDIES}` paired studies and at least one of the six "
            "KG-Microbe trait annotations."
        ),
        (
            f"- The resulting condition-taxon test family contained "
            f"`{eligible_taxon_test_count}` eligible comparisons."
        ),
        (
            "- Taxon direction was tested with an exact McNemar/binomial test "
            "on discordant paired studies."
        ),
        (
            "- Trait contrasts were tested across paired study-level "
            "increased-minus-decreased z scores with a two-sided Wilcoxon "
            "signed-rank test."
        ),
        (
            "- Benjamini-Hochberg correction was applied globally within the "
            "complete condition-taxon family and within the complete "
            "condition-trait family."
        ),
        f"- FDR-significant results require `q <= {Q_VALUE_CUTOFF:.2f}`.",
        (
            f"- Results with `{Q_VALUE_CUTOFF:.2f} < q <= "
            f"{EXPLORATORY_Q_CUTOFF:.2f}` are shown as exploratory, not significant."
        ),
        "",
        "## Results shown",
        "",
        f"- Eligible conditions represented in the evidence landscape: `{eligible_count}`.",
        (
            f"- FDR-selected condition-taxon results: `{len(selected_taxa)}` "
            f"involving `{genus_count}` genera and `{species_count}` species."
        ),
        f"- FDR-selected condition-trait results: `{len(selected_traits)}`.",
        f"- Exploratory condition-taxon results: `{len(exploratory_taxa)}`.",
        f"- Exploratory condition-trait results: `{len(exploratory_traits)}`.",
        (f"- Conditions in the significant-only node/edge tables: `{displayed_count}`."),
        f"- Dashed KG-Microbe annotation links: `{annotation_count}`.",
        f"- Network nodes/edges: `{len(nodes)}` / `{len(edges)}`.",
        "",
        "## Interpretation boundary",
        "",
        (
            "Filled points are the only results called FDR-significant. Hollow amber "
            "points expose potentially useful follow-up hypotheses but are not "
            "discoveries. Taxon tests measure directional asymmetry of reporting in "
            "curated BugSigDB signatures; they are not abundance-effect meta-analyses "
            "because harmonized study effect sizes and standard errors are unavailable. "
            "KG-Microbe annotations do not establish activity, expression, strain-level "
            "co-occurrence, or causality."
        ),
        "",
        (
            f"[All condition-taxon tests]({OUTPUT_STEM}_condition_taxon_tests.tsv) | "
            f"[All condition-trait tests]({OUTPUT_STEM}_condition_trait_tests.tsv) | "
            f"[Condition eligibility]({OUTPUT_STEM}_condition_eligibility.tsv) | "
            f"[Selected taxa]({OUTPUT_STEM}_taxon_selection.tsv) | "
            f"[Edges]({OUTPUT_STEM}_edges.tsv) | "
            f"[Nodes]({OUTPUT_STEM}_nodes.tsv) | "
            f"[Vector figure]({OUTPUT_STEM}.svg)"
        ),
        "",
    ]
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run the complete FDR-tiered analysis and write audited outputs."""
    traits = load_traits()
    metadata = load_metadata()
    feces_metadata, eligibility = build_condition_eligibility(metadata)
    annotated_taxa = load_annotated_taxa(traits)
    taxon_tests = build_condition_taxon_tests(feces_metadata, eligibility, annotated_taxa)
    zscores = build_zscore_matrix(metadata, traits)
    trait_tests, study_contrasts = build_condition_trait_tests(
        metadata,
        zscores,
        traits,
        eligibility,
    )
    nodes, edges, selected_taxa = build_network_tables(
        taxon_tests,
        trait_tests,
        annotated_taxa,
        traits,
    )
    eligibility = augment_condition_eligibility(eligibility, taxon_tests, trait_tests)
    output_path = OUTPUT_DIR / f"{OUTPUT_STEM}.png"
    plot_fdr_evidence_landscape(taxon_tests, trait_tests, traits, output_path)
    write_tsv(nodes, OUTPUT_DIR / f"{OUTPUT_STEM}_nodes.tsv")
    write_tsv(edges, OUTPUT_DIR / f"{OUTPUT_STEM}_edges.tsv")
    write_tsv(selected_taxa, OUTPUT_DIR / f"{OUTPUT_STEM}_taxon_selection.tsv")
    write_tsv(taxon_tests, OUTPUT_DIR / f"{OUTPUT_STEM}_condition_taxon_tests.tsv")
    write_tsv(trait_tests, OUTPUT_DIR / f"{OUTPUT_STEM}_condition_trait_tests.tsv")
    write_tsv(study_contrasts, OUTPUT_DIR / f"{OUTPUT_STEM}_study_trait_contrasts.tsv")
    write_tsv(eligibility, OUTPUT_DIR / f"{OUTPUT_STEM}_condition_eligibility.tsv")
    write_methods_summary(
        eligibility,
        taxon_tests,
        trait_tests,
        nodes,
        edges,
        OUTPUT_DIR / f"{OUTPUT_STEM}_methods.md",
    )
    print(
        "Wrote FDR-tiered evidence landscape: "
        f"{int(taxon_tests['selected_after_bh'].sum())} taxon results, "
        f"{int(trait_tests['selected_after_bh'].sum())} trait results, "
        f"{int(eligibility['displayed_in_network'].sum())} conditions"
    )


if __name__ == "__main__":
    main()
