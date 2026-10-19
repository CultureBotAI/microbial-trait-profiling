"""Create Parkinson's disease cellulose-trait follow-up plots and audit tables."""

from __future__ import annotations

import math
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
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

PD_CONDITION = "Parkinson's disease"
FECES = "Feces"
CELLULOSE_ID = "METPO:2000103|GO:0030245"
CELLULOSE_LABEL = "cellulose catabolic process"
DIRECTION_LABELS = {"inc": "increased", "dec": "decreased"}
TAXA_DIRECTION_TO_PROFILE = {"increased": "inc", "decreased": "dec"}
# Writes the individual-condition Sankey edge tables this script reads.
SANKEY_RECIPE = "just six-shared-trait-sankeys"

INCREASED_COLOR = "#c2410c"
DECREASED_COLOR = "#2563eb"
PD_COLOR = "#d55e00"
OTHER_COLOR = "#64748b"
CELLULOSE_COLOR = "#d55e00"
SLATE_50 = "#f8fafc"
SLATE_100 = "#f1f5f9"
SLATE_200 = "#e2e8f0"
SLATE_300 = "#cbd5e1"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"
TRAIT_COLORS = {
    "facultatively aerobic": "#0072b2",
    "cellulose catabolic process": CELLULOSE_COLOR,
    "coccus shaped": "#009e73",
    "sphere shaped": "#cc79a7",
    "arsenate reduction": "#e69f00",
    "moderately halophilic": "#6a3d9a",
}


def read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    """Read a local TSV file."""
    return pd.read_csv(path, sep="\t", **kwargs)


def clean_text(value: object) -> str:
    """Return a compact string value suitable for grouping."""
    return str(value or "").strip()


def wrap_label(value: object, width: int = 18) -> str:
    """Wrap labels used inside figures without splitting words."""
    return "\n".join(wrap(clean_text(value), width=width, break_long_words=False)) or "unknown"


def _listing(items: list[str]) -> str:
    """Join names as prose ("a, b and c"), or return "none"."""
    if not items:
        return "none"
    if len(items) == 1:
        return items[0]
    return f"{', '.join(items[:-1])} and {items[-1]}"


def save_figure(fig: plt.Figure, output_path: Path) -> None:
    """Save a figure as PNG and SVG."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(output_path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_tsv(frame: pd.DataFrame, output_path: Path) -> None:
    """Write a TSV file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, sep="\t", index=False)


def load_trait_summary() -> pd.DataFrame:
    """Load the six selected shared traits."""
    summary = read_tsv(SIX_SHARED_DIR / "six_shared_trait_summary.tsv")
    summary["trait_rank"] = np.arange(len(summary))
    return summary


def load_profile_metadata() -> pd.DataFrame:
    """Load profile metadata aligned to the profile matrix."""
    metadata = read_tsv(PROFILE_METADATA_PATH, index_col=0)
    metadata = metadata.copy()
    metadata["profile_id"] = metadata.index
    for column in ["study_id", "experiment_id", "inc_dec", "condition", "body_site_group"]:
        metadata[column] = metadata[column].map(clean_text)
    metadata["body_site"] = metadata.get("body_site", "unknown")
    metadata["body_site"] = metadata["body_site"].map(clean_text)
    metadata["total_trait_signal"] = pd.to_numeric(
        metadata["total_trait_signal"],
        errors="coerce",
    ).fillna(0.0)
    return metadata[metadata["inc_dec"].isin(DIRECTION_LABELS)].copy()


def build_zscore_matrix(metadata: pd.DataFrame, trait_ids: list[str]) -> pd.DataFrame:
    """Recreate the six-shared-traits z-score transform."""
    values = load_retained_profiles(PROFILE_MATRIX_PATH, metadata)
    composition = trait_fractions(values).loc[:, trait_ids]
    scaled = (composition * 1000).apply(lambda column: column.map(math.log1p))
    zscores = (scaled - scaled.mean(axis=0)) / scaled.std(axis=0, ddof=0).replace(0, float("nan"))
    return zscores.fillna(0.0)


def load_individual_display_edges(direction_label: str) -> pd.DataFrame:
    """Load individual-condition Cytoscape source edges for one direction."""
    filename = (
        f"six_shared_traits_{direction_label}_individual_conditions_"
        "profile_sankey_display_edges.tsv"
    )
    path = SIX_SHARED_DIR / filename
    return read_tsv(path)


def load_individual_all_edges(direction_label: str) -> pd.DataFrame:
    """Load all individual-condition source edges for one direction."""
    path = (
        SIX_SHARED_DIR
        / f"six_shared_traits_{direction_label}_individual_conditions_profile_sankey_edges.tsv"
    )
    return read_tsv(path)


def build_pd_display_summary() -> pd.DataFrame:
    """Summarize PD edges shown in the individual-condition Cytoscape diagrams."""
    records: list[dict[str, object]] = []
    for direction_label in ["increased", "decreased"]:
        edges = load_individual_display_edges(direction_label)
        subset = edges[edges["condition"].eq(PD_CONDITION)].copy()
        for row in subset.itertuples(index=False):
            records.append(
                {
                    "direction": direction_label,
                    "trait": row.trait_label,
                    "body_site": row.body_site_group,
                    "focus_enriched_profiles": int(row.focus_enriched_profile_count),
                    "focus_studies": int(row.focus_enriched_study_count),
                    "comparator_enriched_profiles": int(row.comparator_enriched_profile_count),
                    "background_no_enriched_profiles": int(
                        row.background_no_enriched_profile_count
                    ),
                    "focus_enriched_mean_z": float(row.focus_enriched_mean_guild_z),
                    "comparator_mean_z": float(row.comparator_mean_guild_z),
                    "focus_minus_comparator_z": float(row.focus_enriched_minus_comparator_guild_z),
                    "rank_score": float(row.rank_score),
                    "displayed_in_map": bool(row.displayed_in_map),
                }
            )
    return pd.DataFrame.from_records(records)


def build_pd_all_edge_summary() -> pd.DataFrame:
    """Summarize all PD individual-condition edges, including non-displayed small edges."""
    frames: list[pd.DataFrame] = []
    for direction_label in ["increased", "decreased"]:
        edges = load_individual_all_edges(direction_label)
        subset = edges[edges["condition"].eq(PD_CONDITION)].copy()
        subset["direction"] = direction_label
        frames.append(subset)
    all_edges = pd.concat(frames, ignore_index=True)
    summary = (
        all_edges.groupby(["direction", "trait_label"], sort=False)
        .agg(
            focus_enriched_profiles=("focus_enriched_profile_count", "sum"),
            focus_studies=("focus_enriched_study_count", "sum"),
            body_sites=("body_site_group", "nunique"),
            edge_rows=("trait_label", "count"),
            max_focus_minus_comparator_z=("focus_enriched_minus_comparator_guild_z", "max"),
        )
        .reset_index()
        .rename(columns={"trait_label": "trait"})
    )
    return summary.sort_values(
        ["direction", "focus_enriched_profiles"],
        ascending=[True, False],
    )


def recount_pd_feces_cellulose_edges() -> pd.DataFrame:
    """Recount the PD fecal cellulose edges as create_six_shared_trait_sankeys.py counts them.

    The Sankey groups body sites with its own rollups (it counts "Feces,Blood" as feces), so
    its counts can differ from ones taken over body_site_group.
    """
    import create_six_shared_trait_sankeys as sankeys

    metadata = sankeys.load_profile_metadata()
    zscores = sankeys.build_zscore_matrix(
        sankeys.load_raw_trait_counts(sankeys.SHARED_TRAITS, metadata), metadata
    )
    trait = next(trait for trait in sankeys.SHARED_TRAITS if trait.trait_id == CELLULOSE_ID)
    scored = sankeys.score_profiles(metadata, zscores, trait)
    records: list[dict[str, object]] = []
    for inc_dec, direction in DIRECTION_LABELS.items():
        config = sankeys.single_trait_map_config(trait, inc_dec)
        source, _, _ = sankeys.individual_condition_grouped_body_source(scored, config)
        edges = sankeys.summarize_edges_from_source(source, config)
        edge = edges[edges["condition"].eq(PD_CONDITION) & edges["body_site_group"].eq(FECES)]
        records.append(
            {
                "direction": direction,
                "focus_enriched_profiles": int(edge["focus_enriched_profile_count"].sum()),
                "focus_studies": int(edge["focus_enriched_study_count"].sum()),
            }
        )
    return pd.DataFrame.from_records(records)


def check_display_edges_current(pd_display: pd.DataFrame, recounted: pd.DataFrame) -> None:
    """Refuse saved Sankey edge tables that the current profile inputs no longer give (#687)."""
    saved = pd_display[
        pd_display["trait"].eq(CELLULOSE_LABEL) & pd_display["body_site"].eq(FECES)
    ].set_index("direction")
    for row in recounted.itertuples(index=False):
        expected = (row.focus_enriched_profiles, row.focus_studies)
        found = (
            (
                int(saved.at[row.direction, "focus_enriched_profiles"]),
                int(saved.at[row.direction, "focus_studies"]),
            )
            if row.direction in saved.index
            else None
        )
        if found != expected:
            raise ValueError(
                "The six-shared-trait Sankey edge tables are stale: the saved PD fecal cellulose "
                f"{row.direction} edge has (profiles, studies) {found}, the current inputs give "
                f"{expected}. Run `{SANKEY_RECIPE}` first."
            )


def build_pd_trait_direction_summary(
    metadata: pd.DataFrame,
    zscores: pd.DataFrame,
    traits: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize the six traits in PD feces profiles by direction."""
    records: list[dict[str, object]] = []
    pd_feces = metadata["condition"].eq(PD_CONDITION) & metadata["body_site_group"].eq(FECES)
    for row in traits.itertuples(index=False):
        trait_id = str(row.trait_id)
        for direction in ["inc", "dec"]:
            mask = pd_feces & metadata["inc_dec"].eq(direction)
            values = zscores.loc[mask, trait_id]
            records.append(
                {
                    "trait_id": trait_id,
                    "trait": row.trait_label,
                    "direction": DIRECTION_LABELS[direction],
                    "profiles": int(mask.sum()),
                    "studies": int(metadata.loc[mask, "study_id"].nunique()),
                    "positive_z_profiles": int((values > 0).sum()),
                    "positive_z_rate": float((values > 0).mean()),
                    "mean_z": float(values.mean()),
                    "median_z": float(values.median()),
                }
            )

    summary = pd.DataFrame.from_records(records)
    pivot = summary.pivot(index=["trait_id", "trait"], columns="direction")
    flattened = pd.DataFrame(index=pivot.index)
    for metric in ["profiles", "studies", "positive_z_profiles", "positive_z_rate", "mean_z"]:
        for direction in ["increased", "decreased"]:
            flattened[f"{direction}_{metric}"] = pivot[(metric, direction)]
    flattened["increased_minus_decreased_mean_z"] = (
        flattened["increased_mean_z"] - flattened["decreased_mean_z"]
    )
    flattened["increased_minus_decreased_positive_z_rate"] = (
        flattened["increased_positive_z_rate"] - flattened["decreased_positive_z_rate"]
    )
    return flattened.reset_index()


def fisher_positive_test(pd_values: pd.Series, other_values: pd.Series) -> tuple[float, float]:
    """Run Fisher's exact test for positive z-rate enrichment."""
    table = [
        [int((pd_values > 0).sum()), int((pd_values <= 0).sum())],
        [int((other_values > 0).sum()), int((other_values <= 0).sum())],
    ]
    odds_ratio, p_value = stats.fisher_exact(table)
    return float(odds_ratio), float(p_value)


def build_pd_vs_other_cellulose_summary(
    metadata: pd.DataFrame,
    zscores: pd.DataFrame,
) -> pd.DataFrame:
    """Compare PD feces to other feces for cellulose z-scores."""
    records: list[dict[str, object]] = []
    for direction in ["inc", "dec"]:
        feces_mask = metadata["body_site_group"].eq(FECES) & metadata["inc_dec"].eq(direction)
        pd_mask = feces_mask & metadata["condition"].eq(PD_CONDITION)
        other_mask = feces_mask & ~metadata["condition"].eq(PD_CONDITION)
        pd_values = zscores.loc[pd_mask, CELLULOSE_ID]
        other_values = zscores.loc[other_mask, CELLULOSE_ID]
        odds_ratio, fisher_p = fisher_positive_test(pd_values, other_values)
        mann_p = stats.mannwhitneyu(pd_values, other_values, alternative="two-sided").pvalue
        for group_label, mask, values in [
            ("PD feces", pd_mask, pd_values),
            ("Other feces", other_mask, other_values),
        ]:
            records.append(
                {
                    "direction": DIRECTION_LABELS[direction],
                    "group": group_label,
                    "profiles": int(mask.sum()),
                    "studies": int(metadata.loc[mask, "study_id"].nunique()),
                    "positive_z_profiles": int((values > 0).sum()),
                    "positive_z_rate": float((values > 0).mean()),
                    "mean_z": float(values.mean()),
                    "median_z": float(values.median()),
                    "pd_vs_other_positive_z_odds_ratio": odds_ratio,
                    "pd_vs_other_fisher_p": float(fisher_p),
                    "pd_vs_other_mannwhitney_p": float(mann_p),
                }
            )
    return pd.DataFrame.from_records(records)


def fit_cluster_ols(
    frame: pd.DataFrame,
    outcome: str,
    terms: list[str],
    cluster_column: str,
) -> pd.DataFrame:
    """Fit OLS with one-way cluster-robust standard errors."""
    y = frame[outcome].to_numpy(dtype=float)
    x = np.column_stack(
        [np.ones(len(frame)), *[frame[term].to_numpy(dtype=float) for term in terms]]
    )
    term_names = ["intercept", *terms]
    xtx_inv = np.linalg.pinv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    residuals = y - x @ beta
    meat = np.zeros((len(term_names), len(term_names)))
    for _cluster, group in frame.groupby(cluster_column, sort=False):
        indices = group.index.to_numpy()
        xg = x[indices, :]
        ug = residuals[indices][:, None]
        meat += xg.T @ ug @ ug.T @ xg
    groups = frame[cluster_column].nunique()
    observations = len(frame)
    parameters = len(term_names)
    if groups > 1 and observations > parameters:
        finite_sample = (groups / (groups - 1)) * ((observations - 1) / (observations - parameters))
    else:
        finite_sample = 1.0
    cov = finite_sample * xtx_inv @ meat @ xtx_inv
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    t_stat = np.divide(beta, se, out=np.full_like(beta, np.nan), where=se > 0)
    df = max(groups - 1, 1)
    p_values = 2 * stats.t.sf(np.abs(t_stat), df=df)
    ci_low = beta - stats.t.ppf(0.975, df=df) * se
    ci_high = beta + stats.t.ppf(0.975, df=df) * se
    return pd.DataFrame(
        {
            "outcome": outcome,
            "term": term_names,
            "coefficient": beta,
            "std_error_clustered_by_study": se,
            "t_stat": t_stat,
            "df": df,
            "p_value": p_values,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "cluster_count": groups,
            "observation_count": observations,
        }
    )


def build_cellulose_model_tables(
    metadata: pd.DataFrame,
    zscores: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit feces-only cellulose enrichment models and group summaries."""
    frame = metadata[metadata["body_site_group"].eq(FECES)].copy()
    frame = frame.reset_index(drop=True)
    frame["cellulose_z"] = zscores.loc[
        metadata["body_site_group"].eq(FECES), CELLULOSE_ID
    ].to_numpy()
    frame["cellulose_positive_z"] = (frame["cellulose_z"] > 0).astype(float)
    frame["is_pd"] = frame["condition"].eq(PD_CONDITION).astype(float)
    frame["is_decreased"] = frame["inc_dec"].eq("dec").astype(float)
    frame["pd_x_decreased"] = frame["is_pd"] * frame["is_decreased"]

    terms = ["is_pd", "is_decreased", "pd_x_decreased"]
    models = pd.concat(
        [
            fit_cluster_ols(frame, "cellulose_z", terms, "study_id"),
            fit_cluster_ols(frame, "cellulose_positive_z", terms, "study_id"),
        ],
        ignore_index=True,
    )
    models["term_label"] = models["term"].map(
        {
            "intercept": "Other feces, increased baseline",
            "is_pd": "PD effect in increased feces",
            "is_decreased": "Decreased effect in other feces",
            "pd_x_decreased": "Extra PD effect in decreased feces",
        }
    )

    summary = (
        frame.assign(group=np.where(frame["is_pd"].eq(1), "PD feces", "Other feces"))
        .groupby(["inc_dec", "group"], sort=False)
        .agg(
            profiles=("cellulose_z", "size"),
            studies=("study_id", "nunique"),
            positive_z_profiles=("cellulose_positive_z", "sum"),
            positive_z_rate=("cellulose_positive_z", "mean"),
            mean_z=("cellulose_z", "mean"),
            median_z=("cellulose_z", "median"),
        )
        .reset_index()
    )
    summary["direction"] = summary["inc_dec"].map(DIRECTION_LABELS)
    return models, summary


def build_taxon_driver_tables(traits: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Trace observed PD feces taxa that carry the cellulose annotation."""
    trait_ids = traits["trait_id"].astype(str).tolist()
    metadata = load_profile_metadata()
    pd_feces = metadata[
        metadata["condition"].eq(PD_CONDITION)
        & metadata["body_site_group"].eq(FECES)
        & metadata["inc_dec"].isin(["inc", "dec"])
    ].copy()
    pd_feces["taxa_direction"] = pd_feces["inc_dec"].map({"inc": "increased", "dec": "decreased"})

    taxa = read_tsv(TAXA_PATH)
    occurrences = taxa.merge(
        pd_feces[["study_id", "experiment_id", "inc_dec", "taxa_direction", "profile_id"]],
        left_on=["study_id", "experiment_id", "direction"],
        right_on=["study_id", "experiment_id", "taxa_direction"],
        how="inner",
    )
    occurrences = occurrences.drop_duplicates(["profile_id", "inc_dec", "taxon_id"]).copy()

    trait_matrix = read_tsv(
        TAXA_TRAIT_MATRIX_PATH,
        usecols=["taxon_id", "taxon_name", "taxonomic_rank", *trait_ids],
    )
    for trait_id in trait_ids:
        trait_matrix[trait_id] = pd.to_numeric(trait_matrix[trait_id], errors="coerce").fillna(0)
    cellulose_carriers = trait_matrix[trait_matrix[CELLULOSE_ID] > 0].copy()
    occurrence_carriers = occurrences.merge(
        cellulose_carriers,
        on="taxon_id",
        how="inner",
        suffixes=("_reported", "_matrix"),
    )
    occurrence_carriers["taxon_name"] = occurrence_carriers["taxon_name_matrix"].fillna(
        occurrence_carriers["taxon_name_reported"]
    )
    occurrence_carriers["taxonomic_rank"] = occurrence_carriers["taxonomic_rank_matrix"].fillna(
        occurrence_carriers["taxonomic_rank_reported"]
    )

    records: list[dict[str, object]] = []
    for (taxon_id, taxon_name, rank), group in occurrence_carriers.groupby(
        ["taxon_id", "taxon_name", "taxonomic_rank"],
        sort=False,
    ):
        inc = group[group["inc_dec"].eq("inc")]
        dec = group[group["inc_dec"].eq("dec")]
        trait_row = cellulose_carriers[cellulose_carriers["taxon_id"].eq(taxon_id)].iloc[0]
        record = {
            "taxon_id": taxon_id,
            "taxon_name": taxon_name,
            "taxonomic_rank": rank,
            "increased_profile_count": int(inc["profile_id"].nunique()),
            "decreased_profile_count": int(dec["profile_id"].nunique()),
            "increased_study_count": int(inc["study_id"].nunique()),
            "decreased_study_count": int(dec["study_id"].nunique()),
            "total_profile_count": int(group["profile_id"].nunique()),
            "total_study_count": int(group["study_id"].nunique()),
            "net_decreased_minus_increased_profiles": int(
                dec["profile_id"].nunique() - inc["profile_id"].nunique()
            ),
            "cellulose_annotation_count": int(trait_row[CELLULOSE_ID]),
        }
        for trait in traits.itertuples(index=False):
            trait_id = str(trait.trait_id)
            record[f"{trait.slug}_annotation_count"] = int(trait_row[trait_id])
            record[f"{trait.slug}_present"] = bool(trait_row[trait_id] > 0)
        records.append(record)

    drivers = pd.DataFrame.from_records(records)
    drivers = drivers.sort_values(
        ["total_profile_count", "decreased_profile_count", "increased_profile_count"],
        ascending=[False, False, False],
    ).reset_index(drop=True)

    top = drivers.head(18).copy()
    trait_presence_records: list[dict[str, object]] = []
    for row in top.itertuples(index=False):
        for trait in traits.itertuples(index=False):
            slug = str(trait.slug)
            trait_presence_records.append(
                {
                    "taxon_id": row.taxon_id,
                    "taxon_name": row.taxon_name,
                    "taxonomic_rank": row.taxonomic_rank,
                    "trait_id": trait.trait_id,
                    "trait": trait.trait_label,
                    "annotation_count": int(getattr(row, f"{slug}_annotation_count")),
                    "present": bool(getattr(row, f"{slug}_present")),
                    "increased_profile_count": int(row.increased_profile_count),
                    "decreased_profile_count": int(row.decreased_profile_count),
                    "net_decreased_minus_increased_profiles": int(
                        row.net_decreased_minus_increased_profiles
                    ),
                }
            )
    presence = pd.DataFrame.from_records(trait_presence_records)
    return drivers, presence


def plot_direction_specificity(
    pd_display: pd.DataFrame,
    pd_vs_other: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot PD cellulose direction specificity in two panels."""
    cellulose = pd_display[
        pd_display["trait"].eq(CELLULOSE_LABEL) & pd_display["body_site"].eq(FECES)
    ].copy()
    cellulose["direction_order"] = cellulose["direction"].map({"increased": 0, "decreased": 1})
    cellulose = cellulose.sort_values("direction_order")

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2), constrained_layout=True)

    ax = axes[0]
    x = np.arange(len(cellulose))
    width = 0.36
    colors = [INCREASED_COLOR, DECREASED_COLOR]
    bars_focus = ax.bar(
        x - width / 2,
        cellulose["focus_enriched_profiles"],
        width=width,
        color=colors,
        alpha=0.84,
    )
    bars_comparator = ax.bar(
        x + width / 2,
        cellulose["comparator_enriched_profiles"],
        width=width,
        color=SLATE_300,
        alpha=0.94,
    )
    ax.bar_label(bars_focus, padding=3, fontsize=9)
    ax.bar_label(bars_comparator, padding=3, fontsize=9)
    ax.set_xticks(x, ["Increased\nPD feces", "Decreased\nPD feces"])
    ax.set_ylabel("Profiles with cellulose z > 0")
    ax.set_title("PD fecal profiles with cellulose z > 0", fontweight="bold")
    ax.legend(
        handles=[
            Patch(facecolor=INCREASED_COLOR, label="Increased focus edge"),
            Patch(facecolor=DECREASED_COLOR, label="Decreased focus edge"),
            Patch(facecolor=SLATE_300, label="Opposite-direction comparator"),
        ],
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=3,
        fontsize=8.6,
    )
    ax.grid(axis="y", color=SLATE_200, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(
        0,
        max(
            cellulose["focus_enriched_profiles"].max(),
            cellulose["comparator_enriched_profiles"].max(),
        )
        + 18,
    )
    for index, row in enumerate(cellulose.itertuples(index=False)):
        ax.text(
            index,
            max(row.focus_enriched_profiles, row.comparator_enriched_profiles) + 8,
            f"s={row.focus_studies}\nmean z={row.focus_enriched_mean_z:.2f}\n"
            f"delta z={row.focus_minus_comparator_z:.2f}",
            ha="center",
            va="bottom",
            fontsize=8.4,
            color=SLATE_700,
        )

    ax = axes[1]
    rate = pd_vs_other.copy()
    rate["direction_order"] = rate["direction"].map({"increased": 0, "decreased": 1})
    rate["group_order"] = rate["group"].map({"Other feces": 0, "PD feces": 1})
    rate = rate.sort_values(["direction_order", "group_order"])
    direction_positions = {"increased": 0, "decreased": 1}
    group_offsets = {"Other feces": -width / 2, "PD feces": width / 2}
    group_colors = {"Other feces": OTHER_COLOR, "PD feces": PD_COLOR}
    for group_label, group in rate.groupby("group", sort=False):
        xpos = [direction_positions[d] + group_offsets[group_label] for d in group["direction"]]
        bars = ax.bar(
            xpos,
            group["positive_z_rate"] * 100,
            width=width,
            color=group_colors[group_label],
            alpha=0.84,
            label=group_label,
        )
        ax.bar_label(bars, labels=[f"{v:.1f}%" for v in group["positive_z_rate"] * 100], padding=3)
    ax.set_xticks([0, 1], ["Increased", "Decreased"])
    ax.set_ylabel("Profiles with cellulose z > 0 (%)")
    ax.set_title("PD vs other fecal profiles", fontweight="bold")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(axis="y", color=SLATE_200, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    dec_pd = pd_vs_other[
        pd_vs_other["direction"].eq("decreased") & pd_vs_other["group"].eq("PD feces")
    ].iloc[0]
    ax.text(
        0.98,
        0.05,
        f"Decreased PD vs other feces:\nOR={dec_pd.pd_vs_other_positive_z_odds_ratio:.2f}, "
        f"Fisher p={dec_pd.pd_vs_other_fisher_p:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color=SLATE_700,
        bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": SLATE_200},
    )

    fig.suptitle(
        "Cellulose catabolic process annotation in PD fecal profiles, by signature direction",
        fontsize=15,
        fontweight="bold",
    )
    save_figure(fig, output_path)


def plot_pd_trait_lollipop(summary: pd.DataFrame, output_path: Path) -> None:
    """Plot PD feces increased-minus-decreased z contrasts for all six traits."""
    frame = summary.sort_values("increased_minus_decreased_mean_z")
    y = np.arange(len(frame))
    colors = [
        CELLULOSE_COLOR if trait == CELLULOSE_LABEL else TRAIT_COLORS.get(trait, SLATE_500)
        for trait in frame["trait"]
    ]
    fig, ax = plt.subplots(figsize=(9.6, 5.3), constrained_layout=True)
    ax.axvline(0, color=SLATE_900, linewidth=1.0)
    for y_pos, row, color in zip(y, frame.itertuples(index=False), colors, strict=True):
        ax.plot(
            [0, row.increased_minus_decreased_mean_z],
            [y_pos, y_pos],
            color=color,
            linewidth=3.2 if row.trait == CELLULOSE_LABEL else 2.2,
            alpha=0.82,
            solid_capstyle="round",
        )
        ax.scatter(
            row.increased_minus_decreased_mean_z,
            y_pos,
            s=96 if row.trait == CELLULOSE_LABEL else 64,
            color=color,
            edgecolor=SLATE_900 if row.trait == CELLULOSE_LABEL else "white",
            linewidth=1.0,
            zorder=3,
        )
        ax.text(
            row.increased_minus_decreased_mean_z,
            y_pos + 0.22,
            f"{row.increased_minus_decreased_mean_z:+.2f}",
            ha="center",
            va="bottom",
            fontsize=8.4,
            color=SLATE_700,
        )
    ax.set_yticks(y, frame["trait"])
    ax.set_xlabel(
        "PD feces mean z difference: increased - decreased\n"
        "Left of zero = decreased-profile skew; right of zero = increased-profile skew"
    )
    ax.set_title("Six-trait PD feces direction contrast", fontweight="bold")
    ax.grid(axis="x", color=SLATE_200, linewidth=0.8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    save_figure(fig, output_path)


def plot_cytoscape_highlight(pd_display: pd.DataFrame, output_path: Path) -> None:
    """Plot focused Cytoscape-style PD-to-trait-to-feces highlights."""
    frame = pd_display[pd_display["body_site"].eq(FECES)].copy()
    directions = ["increased", "decreased"]
    trait_order = [
        "facultatively aerobic",
        "cellulose catabolic process",
        "coccus shaped",
        "sphere shaped",
        "arsenate reduction",
        "moderately halophilic",
    ]
    min_n = float(frame["focus_enriched_profiles"].min())
    max_n = float(frame["focus_enriched_profiles"].max())

    def line_width(value: float) -> float:
        if max_n <= min_n:
            return 4.0
        return (
            1.6
            + (math.sqrt(value) - math.sqrt(min_n)) / (math.sqrt(max_n) - math.sqrt(min_n)) * 7.4
        )

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.8), constrained_layout=True)
    for ax, direction in zip(axes, directions, strict=True):
        subset = frame[frame["direction"].eq(direction)].set_index("trait")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        condition_pos = (0.12, 0.50)
        body_pos = (0.90, 0.50)
        trait_y = np.linspace(0.82, 0.18, len(trait_order))
        trait_positions = dict(zip(trait_order, trait_y, strict=True))

        for trait in trait_order:
            if trait not in subset.index:
                continue
            row = subset.loc[trait]
            color = TRAIT_COLORS.get(trait, SLATE_500)
            is_cellulose = trait == CELLULOSE_LABEL
            alpha = 0.93 if is_cellulose else 0.22
            width = line_width(float(row.focus_enriched_profiles))
            if is_cellulose:
                width += 2.2
            trait_pos = (0.50, trait_positions[trait])
            bundle_pos = (0.76, 0.50 + (trait_positions[trait] - 0.50) * 0.18)
            ax.plot(
                [condition_pos[0], trait_pos[0]],
                [condition_pos[1], trait_pos[1]],
                color=color,
                linewidth=width,
                alpha=alpha,
                solid_capstyle="round",
                zorder=1,
            )
            ax.plot(
                [trait_pos[0], bundle_pos[0], body_pos[0]],
                [trait_pos[1], bundle_pos[1], body_pos[1]],
                color=color,
                linewidth=width,
                alpha=alpha,
                solid_capstyle="round",
                zorder=1,
            )
            ax.scatter(
                [trait_pos[0]],
                [trait_pos[1]],
                s=760 if is_cellulose else 430,
                color=color,
                edgecolor=SLATE_900 if is_cellulose else "white",
                linewidth=1.2 if is_cellulose else 0.8,
                zorder=3,
            )
            ax.scatter(
                [bundle_pos[0]],
                [bundle_pos[1]],
                s=170 if is_cellulose else 90,
                marker="D",
                color=color,
                edgecolor="none",
                alpha=0.95 if is_cellulose else 0.42,
                zorder=3,
            )
            if is_cellulose:
                ax.text(
                    trait_pos[0],
                    trait_pos[1],
                    "cellulose\ncatabolic\nprocess",
                    ha="center",
                    va="center",
                    fontsize=8.6,
                    color="white",
                    fontweight="bold",
                    linespacing=0.92,
                    zorder=4,
                )
                ax.text(
                    0.50,
                    0.07,
                    f"cellulose edge: n={int(row.focus_enriched_profiles)}, "
                    f"s={int(row.focus_studies)}, delta z={row.focus_minus_comparator_z:.2f}",
                    ha="center",
                    va="bottom",
                    fontsize=9.2,
                    color=SLATE_700,
                )

        ax.scatter(
            [condition_pos[0]],
            [condition_pos[1]],
            s=900,
            color="white",
            edgecolor=SLATE_700,
            linewidth=1.4,
            zorder=4,
        )
        ax.text(
            condition_pos[0],
            condition_pos[1],
            "Parkinson's\ndisease",
            ha="center",
            va="center",
            fontsize=9.6,
            color=SLATE_900,
            fontweight="bold",
            zorder=5,
        )
        ax.scatter(
            [body_pos[0]],
            [body_pos[1]],
            s=900,
            color="white",
            edgecolor=SLATE_700,
            linewidth=1.4,
            zorder=4,
        )
        ax.text(
            body_pos[0],
            body_pos[1],
            "Feces",
            ha="center",
            va="center",
            fontsize=10,
            color=SLATE_900,
            fontweight="bold",
            zorder=5,
        )
        ax.set_title(f"{direction.title()} individual Cytoscape view", fontweight="bold")
        ax.text(
            0.50,
            0.95,
            "Other PD six-trait fecal links are faded for context",
            ha="center",
            va="top",
            fontsize=9,
            color=SLATE_700,
        )

    fig.suptitle("Focused Cytoscape highlight: PD to cellulose to feces", fontweight="bold")
    save_figure(fig, output_path)


def plot_enrichment_test(
    model_table: pd.DataFrame,
    group_summary: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot the cellulose PD-by-direction enrichment test."""
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.4), constrained_layout=True)

    ax = axes[0]
    ordered = group_summary.copy()
    ordered["order"] = ordered["direction"].map({"increased": 0, "decreased": 1}) * 2 + ordered[
        "group"
    ].map({"Other feces": 0, "PD feces": 1})
    ordered = ordered.sort_values("order")
    labels = [f"{row.direction}\n{row.group}".replace(" ", "\n") for row in ordered.itertuples()]
    colors = [
        OTHER_COLOR if row.group == "Other feces" else PD_COLOR for row in ordered.itertuples()
    ]
    bars = ax.bar(
        np.arange(len(ordered)),
        ordered["mean_z"],
        color=colors,
        alpha=0.84,
    )
    ax.axhline(0, color=SLATE_900, linewidth=1.0)
    ax.bar_label(bars, labels=[f"{v:.2f}" for v in ordered["mean_z"]], padding=3, fontsize=8.8)
    ax.set_xticks(np.arange(len(ordered)), labels)
    ax.set_ylabel("Mean cellulose z-score")
    ax.set_title("Feces-only group means", fontweight="bold")
    ax.grid(axis="y", color=SLATE_200, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    terms = ["is_pd", "is_decreased", "pd_x_decreased"]
    plot_frame = model_table[
        model_table["outcome"].eq("cellulose_z") & model_table["term"].isin(terms)
    ].copy()
    plot_frame["order"] = plot_frame["term"].map(
        {"is_pd": 2, "is_decreased": 1, "pd_x_decreased": 0}
    )
    plot_frame = plot_frame.sort_values("order")
    y = np.arange(len(plot_frame))
    colors = [
        CELLULOSE_COLOR if term == "pd_x_decreased" else SLATE_500 for term in plot_frame["term"]
    ]
    ax.axvline(0, color=SLATE_900, linewidth=1.0)
    for y_pos, row, color in zip(y, plot_frame.itertuples(index=False), colors, strict=True):
        ax.plot([row.ci_low, row.ci_high], [y_pos, y_pos], color=color, linewidth=2.2)
        ax.scatter(row.coefficient, y_pos, s=90, color=color, edgecolor="white", zorder=3)
        p_label = "p<0.001" if row.p_value < 0.001 else f"p={row.p_value:.3f}"
        ax.text(
            row.coefficient,
            y_pos + 0.20,
            f"{row.coefficient:+.2f}, {p_label}",
            ha="center",
            va="bottom",
            fontsize=8.4,
            color=SLATE_700,
        )
    # Headroom keeps the top estimate's label below the panel title.
    ax.set_ylim(-0.5, len(plot_frame) - 0.2)
    ax.set_yticks(y, plot_frame["term_label"])
    ax.set_xlabel("OLS coefficient on cellulose z-score")
    ax.set_title("Study-clustered interaction model", fontweight="bold")
    ax.grid(axis="x", color=SLATE_200, linewidth=0.8)
    ax.spines[["top", "right", "left"]].set_visible(False)

    fig.suptitle(
        "Cellulose annotation z-score in fecal profiles by PD status and signature direction",
        fontweight="bold",
    )
    save_figure(fig, output_path)


def plot_taxon_driver_bar(drivers: pd.DataFrame, output_path: Path) -> None:
    """Plot how often PD fecal signatures report each cellulose-annotated taxon, by direction."""
    frame = drivers.head(18).sort_values("total_profile_count", ascending=True)
    labels = [f"{row.taxon_name} ({row.taxonomic_rank})" for row in frame.itertuples()]
    y = np.arange(len(frame))
    fig, ax = plt.subplots(figsize=(10.8, 7.2), constrained_layout=True)
    ax.barh(
        y,
        -frame["decreased_profile_count"],
        color=DECREASED_COLOR,
        alpha=0.78,
        label="Decreased PD feces",
    )
    ax.barh(
        y,
        frame["increased_profile_count"],
        color=INCREASED_COLOR,
        alpha=0.78,
        label="Increased PD feces",
    )
    for y_pos, row in zip(y, frame.itertuples(index=False), strict=True):
        if row.decreased_profile_count:
            ax.text(
                -row.decreased_profile_count - 0.8,
                y_pos,
                str(row.decreased_profile_count),
                ha="right",
                va="center",
                fontsize=8.4,
                color=SLATE_700,
            )
        if row.increased_profile_count:
            ax.text(
                row.increased_profile_count + 0.8,
                y_pos,
                str(row.increased_profile_count),
                ha="left",
                va="center",
                fontsize=8.4,
                color=SLATE_700,
            )
    ax.axvline(0, color=SLATE_900, linewidth=1.0)
    max_count = int(
        max(frame["decreased_profile_count"].max(), frame["increased_profile_count"].max())
    )
    ax.set_xlim(-max_count - 8, max_count + 8)
    ax.set_yticks(y, labels)
    ax.set_xlabel(
        "PD fecal profiles reporting the taxon "
        "(annotated with cellulose catabolic process, itself or via descendants)"
    )
    ax.set_title(
        "PD fecal profiles reporting cellulose-annotated taxa, by signature direction",
        fontweight="bold",
    )
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="x", color=SLATE_200, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, output_path)


def plot_taxon_trait_heatmap(
    drivers: pd.DataFrame,
    presence: pd.DataFrame,
    traits: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot which other shared traits the top cellulose-annotated taxa are annotated with."""
    top_taxa = drivers.head(18).copy()
    top_ids = top_taxa["taxon_id"].tolist()
    trait_labels = traits["trait_label"].tolist()
    matrix = (
        presence[presence["taxon_id"].isin(top_ids)]
        .pivot(index="taxon_id", columns="trait", values="annotation_count")
        .reindex(index=top_ids, columns=trait_labels)
        .fillna(0)
        .astype(int)
    )
    binary = (matrix > 0).astype(int)
    labels = [
        f"{row.taxon_name} ({row.taxonomic_rank})\n"
        f"dec={row.decreased_profile_count}, inc={row.increased_profile_count}"
        for row in top_taxa.itertuples(index=False)
    ]

    fig, ax = plt.subplots(figsize=(11.8, 8.4), constrained_layout=True)
    cmap = ListedColormap([SLATE_100, "#f59e0b"])
    ax.imshow(binary.to_numpy(), aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(trait_labels)), [wrap_label(label, 12) for label in trait_labels])
    ax.set_yticks(np.arange(len(labels)), labels)
    ax.tick_params(axis="x", labelrotation=45)
    for y_pos in range(matrix.shape[0]):
        for x_pos in range(matrix.shape[1]):
            value = int(matrix.iloc[y_pos, x_pos])
            if value > 0:
                ax.text(
                    x_pos,
                    y_pos,
                    str(value),
                    ha="center",
                    va="center",
                    fontsize=7.8,
                    color=SLATE_900,
                    fontweight="bold",
                )
    ax.set_title(
        "Six-trait annotations of the cellulose-annotated taxa most often reported in PD feces",
        fontweight="bold",
    )
    ax.set_xlabel("Six shared traits")
    ax.set_ylabel("Reported taxon label")
    ax.set_xticks(np.arange(-0.5, len(trait_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)
    save_figure(fig, output_path)


def write_interpretations(
    output_path: Path,
    pd_display: pd.DataFrame,
    pd_vs_other: pd.DataFrame,
    trait_summary: pd.DataFrame,
    model_table: pd.DataFrame,
    drivers: pd.DataFrame,
) -> None:
    """Write a Markdown interpretation note for the generated outputs."""
    cell_display = pd_display[
        pd_display["trait"].eq(CELLULOSE_LABEL) & pd_display["body_site"].eq(FECES)
    ].set_index("direction")
    inc_cell = cell_display.loc["increased"]
    dec_cell = cell_display.loc["decreased"]
    dec_vs_other = pd_vs_other[
        pd_vs_other["direction"].eq("decreased") & pd_vs_other["group"].eq("PD feces")
    ].iloc[0]
    dec_other = pd_vs_other[
        pd_vs_other["direction"].eq("decreased") & pd_vs_other["group"].eq("Other feces")
    ].iloc[0]
    inc_vs_other = pd_vs_other[
        pd_vs_other["direction"].eq("increased") & pd_vs_other["group"].eq("PD feces")
    ].iloc[0]
    inc_other = pd_vs_other[
        pd_vs_other["direction"].eq("increased") & pd_vs_other["group"].eq("Other feces")
    ].iloc[0]
    contrasts = trait_summary.set_index("trait")["increased_minus_decreased_mean_z"].sort_values()
    cellulose_contrast = float(contrasts[CELLULOSE_LABEL])
    decreased_skewed = [trait for trait, value in contrasts.items() if value < 0]
    increased_skewed = [trait for trait, value in contrasts[::-1].items() if value > 0]
    interaction = model_table[
        model_table["outcome"].eq("cellulose_z") & model_table["term"].eq("pd_x_decreased")
    ].iloc[0]
    interval = "excludes" if interaction.ci_low > 0 or interaction.ci_high < 0 else "includes"
    top_taxa = drivers.head(8)
    plotted_taxa = drivers.head(18)
    decreased_labels = plotted_taxa.loc[
        plotted_taxa["decreased_profile_count"] > plotted_taxa["increased_profile_count"],
        "taxon_name",
    ].tolist()

    lines = [
        "# Parkinson's Disease and Cellulose Catabolic Process",
        "",
        "## Summary",
        "",
        (
            "This folder follows the cellulose catabolic process trait "
            f"(`{CELLULOSE_ID}`) of the six-shared-traits analysis in Parkinson's disease "
            "(PD) fecal profiles. Each profile is one curated BugSigDB signature. Its trait "
            "values summarize KG-Microbe annotations of the reported taxa, including "
            "annotations rolled up from their descendants; they are not measurements of "
            "cellulose degradation, gene expression or metabolites."
        ),
        "",
        (
            "In PD fecal profiles the mean cellulose annotation z-score is "
            f"{'lower' if cellulose_contrast < 0 else 'higher'} in increased than in "
            f"decreased profiles (increased minus decreased `{cellulose_contrast:.3f}`). "
            "A study-clustered model asks whether the decreased-minus-increased difference "
            "is larger in PD than in other fecal profiles: the extra PD effect in decreased "
            "profiles is "
            f"`{interaction.coefficient:+.3f}` cellulose z (95% CI "
            f"`{interaction.ci_low:+.3f}` to `{interaction.ci_high:+.3f}`, "
            f"`p={interaction.p_value:.3g}`, {int(interaction.cluster_count)} study "
            f"clusters), and the interval {interval} zero. This is one nominal test, not "
            "corrected for the other conditions and traits examined; the corrected screen "
            "of every eligible fecal condition-trait contrast is described in "
            "`pd_condition_taxon_trait_network_methods.md`."
        ),
        "",
        "## Cellulose Edges by Signature Direction",
        "",
        (
            f"In the individual-condition Cytoscape display edges, PD feces has "
            f"`{int(inc_cell.focus_enriched_profiles)}` increased profiles with cellulose "
            f"z > 0 from `{int(inc_cell.focus_studies)}` studies, and "
            f"`{int(dec_cell.focus_enriched_profiles)}` decreased profiles with cellulose "
            f"z > 0 from `{int(dec_cell.focus_studies)}` studies."
        ),
        "",
        (
            "The focus-minus-comparator z contrast is "
            f"`{dec_cell.focus_minus_comparator_z:.3f}` for the decreased edge and "
            f"`{inc_cell.focus_minus_comparator_z:.3f}` for the increased edge."
        ),
        "",
        "See `pd_cellulose_direction_specificity.png` and `pd_cellulose_cytoscape_highlight.png`.",
        "",
        "## PD Six-Trait Context",
        "",
        (
            f"Among the six traits in PD feces, cellulose has an increased-minus-"
            f"decreased mean z difference of `{cellulose_contrast:.3f}`. Negative "
            "values indicate decreased-profile skew."
        ),
        "",
        (
            "In the companion lollipop plot the decreased-skewed traits, most negative "
            f"first, are {_listing(decreased_skewed)}; the increased-skewed traits, most "
            f"positive first, are {_listing(increased_skewed)}."
        ),
        "",
        "See `pd_six_trait_direction_lollipop.png`.",
        "",
        "## PD vs Other Fecal Profiles",
        "",
        (
            "Among decreased fecal profiles, cellulose z > 0 in "
            f"`{dec_vs_other.positive_z_rate:.1%}` of PD profiles and "
            f"`{dec_other.positive_z_rate:.1%}` of other profiles; among increased fecal "
            f"profiles, in `{inc_vs_other.positive_z_rate:.1%}` and "
            f"`{inc_other.positive_z_rate:.1%}`."
        ),
        "",
        (
            f"The decreased-direction PD-vs-other Fisher test gives an odds ratio "
            f"of `{dec_vs_other.pd_vs_other_positive_z_odds_ratio:.2f}` and "
            f"`p={dec_vs_other.pd_vs_other_fisher_p:.3g}`. It treats the profiles of one "
            "study as independent, so the study-clustered model below is the better guide."
        ),
        "",
        (
            f"The feces-only linear model includes PD, decreased direction, and the "
            f"PD-by-decreased interaction, with standard errors clustered by study. The "
            f"cellulose z-score interaction term is `{interaction.coefficient:.3f}` "
            f"(95% CI `{interaction.ci_low:.3f}` to `{interaction.ci_high:.3f}`, "
            f"`p={interaction.p_value:.3g}`). It tests whether the decreased direction has "
            "an extra PD-associated shift in cellulose annotation beyond the main PD and "
            "direction effects."
        ),
        "",
        "See `pd_cellulose_enrichment_test.png` and `pd_cellulose_interaction_model.tsv`.",
        "",
        "## Reported Taxa With the Cellulose Annotation",
        "",
        (
            "Observed PD feces taxa were intersected with the KG-Microbe "
            "taxa-trait matrix, then limited to taxa annotated with cellulose "
            "catabolic process, themselves or through annotated descendants. Counts are "
            "profile-level occurrences after deduplicating each taxon within each profile."
        ),
        "",
        "Cellulose-annotated taxon labels most often reported in PD feces:",
        "",
        "| taxon | rank | decreased profiles | increased profiles |",
        "| --- | --- | ---: | ---: |",
    ]
    for row in top_taxa.itertuples(index=False):
        lines.append(
            f"| {row.taxon_name} | {row.taxonomic_rank} | "
            f"{int(row.decreased_profile_count)} | {int(row.increased_profile_count)} |"
        )
    lines.extend(
        [
            "",
            (
                f"Of the {len(plotted_taxa)} labels plotted, those reported in more decreased "
                f"than increased PD fecal profiles are {_listing(decreased_labels)}. The "
                "labels mix genus, family and higher ranks, so they are BugSigDB reporting "
                "units rather than species-resolved organisms, and an annotation reached "
                "through descendants does not mean every member of the taxon carries it. "
                "The counts describe how often PD signatures report these labels; they do "
                "not show that the taxa produce the cellulose contrast."
            ),
            "",
            (
                "The taxon-trait heatmap shows which of the other five shared traits the "
                "same taxa are annotated with. A higher-rank label can collect different "
                "annotations from different descendants, so co-annotation does not show "
                "that the traits occur together in one organism."
            ),
            "",
            "See `pd_cellulose_taxon_driver_bar.png`, "
            "`pd_cellulose_taxon_six_trait_heatmap.png`, and "
            "`pd_cellulose_taxon_drivers.tsv`.",
            "",
            "## Generated Files",
            "",
            "The file names predate the current figure titles; the titles describe what "
            "each figure shows.",
            "",
            "- `pd_cellulose_direction_specificity.png`: PD fecal profiles with cellulose "
            "z > 0 by signature direction, and PD vs other feces.",
            "- `pd_six_trait_direction_lollipop.png`: all-six-trait PD feces direction contrast.",
            "- `pd_cellulose_cytoscape_highlight.png`: focused Cytoscape-style "
            "PD-cellulose highlight.",
            "- `pd_cellulose_taxon_driver_bar.png`: PD fecal profiles reporting each "
            "cellulose-annotated taxon, by direction.",
            "- `pd_cellulose_taxon_six_trait_heatmap.png`: six-trait annotations of the "
            "most often reported cellulose-annotated taxa.",
            "- `pd_cellulose_enrichment_test.png`: feces-only model and group mean summary.",
            "- `pd_cellulose_display_edge_summary.tsv`: PD display-edge audit data.",
            "- `pd_cellulose_all_edge_summary.tsv`: PD all-edge audit data.",
            "- `pd_six_trait_direction_summary.tsv`: PD feces six-trait z-score summary.",
            "- `pd_cellulose_vs_other_feces_summary.tsv`: PD vs other feces comparison.",
            "- `pd_cellulose_interaction_model.tsv`: study-clustered OLS model table.",
            "- `pd_cellulose_model_group_summary.tsv`: model input group summaries.",
            "- `pd_cellulose_taxon_drivers.tsv`: audit table of the cellulose-annotated "
            "taxa reported in PD feces.",
            "- `pd_cellulose_taxon_trait_presence.tsv`: long-form top-taxon trait matrix.",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Create all PD-cellulose outputs."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    traits = load_trait_summary()
    trait_ids = traits["trait_id"].astype(str).tolist()
    metadata = load_profile_metadata()
    zscores = build_zscore_matrix(metadata, trait_ids)

    pd_display = build_pd_display_summary()
    check_display_edges_current(pd_display, recount_pd_feces_cellulose_edges())
    pd_all_edges = build_pd_all_edge_summary()
    pd_trait_summary = build_pd_trait_direction_summary(metadata, zscores, traits)
    pd_vs_other = build_pd_vs_other_cellulose_summary(metadata, zscores)
    model_table, model_group_summary = build_cellulose_model_tables(metadata, zscores)
    taxon_drivers, taxon_presence = build_taxon_driver_tables(traits)

    write_tsv(pd_display, OUTPUT_DIR / "pd_cellulose_display_edge_summary.tsv")
    write_tsv(pd_all_edges, OUTPUT_DIR / "pd_cellulose_all_edge_summary.tsv")
    write_tsv(pd_trait_summary, OUTPUT_DIR / "pd_six_trait_direction_summary.tsv")
    write_tsv(pd_vs_other, OUTPUT_DIR / "pd_cellulose_vs_other_feces_summary.tsv")
    write_tsv(model_table, OUTPUT_DIR / "pd_cellulose_interaction_model.tsv")
    write_tsv(model_group_summary, OUTPUT_DIR / "pd_cellulose_model_group_summary.tsv")
    write_tsv(taxon_drivers, OUTPUT_DIR / "pd_cellulose_taxon_drivers.tsv")
    write_tsv(taxon_presence, OUTPUT_DIR / "pd_cellulose_taxon_trait_presence.tsv")

    plot_direction_specificity(
        pd_display,
        pd_vs_other,
        OUTPUT_DIR / "pd_cellulose_direction_specificity.png",
    )
    plot_pd_trait_lollipop(
        pd_trait_summary,
        OUTPUT_DIR / "pd_six_trait_direction_lollipop.png",
    )
    plot_cytoscape_highlight(
        pd_display,
        OUTPUT_DIR / "pd_cellulose_cytoscape_highlight.png",
    )
    plot_enrichment_test(
        model_table,
        model_group_summary,
        OUTPUT_DIR / "pd_cellulose_enrichment_test.png",
    )
    plot_taxon_driver_bar(
        taxon_drivers,
        OUTPUT_DIR / "pd_cellulose_taxon_driver_bar.png",
    )
    plot_taxon_trait_heatmap(
        taxon_drivers,
        taxon_presence,
        traits,
        OUTPUT_DIR / "pd_cellulose_taxon_six_trait_heatmap.png",
    )
    write_interpretations(
        OUTPUT_DIR / "interpretations.md",
        pd_display,
        pd_vs_other,
        pd_trait_summary,
        model_table,
        taxon_drivers,
    )

    print(f"Wrote PD-cellulose outputs: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
