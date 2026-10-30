"""Build the PD ecological-reporter tables from the profile inputs and draw the figure.

The condition/trait summary and study contrasts are paired within-study contrasts of the
six-trait profile z-scores that the other PD-cellulose outputs use. The taxon panel crosses
genus and species labels reported in PD fecal signatures with their KG-Microbe annotations
for the six traits. Each profile takes its own experiment's condition and body site.
"""

from __future__ import annotations

import math
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox
from scipy import stats

from create_pd_cellulose_analysis import (
    FECES,
    TAXA_PATH,
    TAXA_TRAIT_MATRIX_PATH,
    build_pd_trait_direction_summary,
    build_zscore_matrix,
    load_profile_metadata,
    load_trait_summary,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "trait-visualization" / "trait-profile-investigation" / "PD-cellulose"
STEM = "pd_ecological_reporter_combined"
SUMMARY_PATH = OUTPUT_DIR / f"{STEM}_condition_trait_summary.tsv"
STUDY_CONTRASTS_PATH = OUTPUT_DIR / f"{STEM}_study_contrasts.tsv"
TAXON_PANEL_PATH = OUTPUT_DIR / f"{STEM}_taxon_panel.tsv"
SELECTION_PATH = OUTPUT_DIR / f"{STEM}_condition_selection.tsv"
OUTPUT_PATH = OUTPUT_DIR / f"{STEM}.png"
FIGURES_NOTE_PATH = OUTPUT_DIR / "pd_ecological_reporter_figures.md"
NUMBERS_BEGIN = "<!-- BEGIN numbers written by `just pd-ecological-reporter-combined` -->"
NUMBERS_END = "<!-- END numbers written by `just pd-ecological-reporter-combined` -->"

PD_CONDITION = "Parkinson's disease"
DIRECTION_LABELS = {"inc": "increased", "dec": "decreased"}
ALLOWED_TAXON_RANKS = {"genus", "species"}
MAX_PANEL_TAXA = 11
CONTRAST_DEFINITION = "mean z increased - mean z decreased"
MIN_PAIRED_STUDIES = 3
SIX_TRAITS = [
    "moderately halophilic",
    "facultatively aerobic",
    "arsenate reduction",
    "coccus shaped",
    "sphere shaped",
    "cellulose catabolic process",
]
INCREASED_SIDE_TRAITS = [
    "moderately halophilic",
    "facultatively aerobic",
    "arsenate reduction",
]
CELLULOSE = "cellulose catabolic process"
THREE_TRAIT_COMPOSITE = "three-trait increased-side composite"
REPORTER_COMPOSITE = "ecological-reporter composite"
PANEL_A_TRAITS = [*SIX_TRAITS, REPORTER_COMPOSITE]
# Panel B label offsets in points, nearest first: the four diagonals, then the four sides.
LABEL_OFFSETS = [
    (dx * radius, dy * radius)
    for radius in (5.0, 10.0, 16.0, 24.0)
    for dx, dy in [(1, 1), (1, -1), (-1, 1), (-1, -1), (1.3, 0), (-1.3, 0), (0, 1.3), (0, -1.3)]
]

PD_COLOR = "#d55e00"
SLATE_100 = "#f1f5f9"
SLATE_200 = "#e2e8f0"
SLATE_300 = "#cbd5e1"
SLATE_500 = "#64748b"
SLATE_700 = "#334155"
SLATE_900 = "#0f172a"


def read_tsv(path: Path) -> pd.DataFrame:
    """Read an audited tab-separated table."""
    return pd.read_csv(path, sep="\t")


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    """Write a stable audit table with explicit missing values."""
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.8g")


def _as_bool(values: pd.Series) -> pd.Series:
    """Normalize boolean-like audit values."""
    return values.astype(str).str.casefold().eq("true")


def _confidence_interval(values: pd.Series) -> tuple[float, float]:
    """Return a descriptive two-sided 95% t interval for study-level contrasts."""
    clean = pd.to_numeric(values, errors="coerce").dropna()
    if len(clean) < 2:
        return math.nan, math.nan
    standard_error = float(clean.std(ddof=1) / math.sqrt(len(clean)))
    critical = float(stats.t.ppf(0.975, len(clean) - 1))
    mean = float(clean.mean())
    return mean - critical * standard_error, mean + critical * standard_error


def build_condition_trait_summaries(
    metadata: pd.DataFrame,
    zscores: pd.DataFrame,
    traits: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build paired within-study direction contrasts and condition-level summaries.

    A study contributes to a fecal condition when it has both increased and decreased
    profiles for it; its contrast is the mean increased z minus the mean decreased z.
    """
    trait_ids = traits["trait_id"].tolist()
    trait_names = dict(zip(traits["trait_id"], traits["trait_label"], strict=True))
    missing_composite = (set(INCREASED_SIDE_TRAITS) | {CELLULOSE}) - set(trait_names.values())
    if missing_composite:
        raise ValueError(f"Missing composite traits: {sorted(missing_composite)}")

    profiles = zscores[trait_ids].join(
        metadata[
            [
                "study_id",
                "experiment_id",
                "inc_dec",
                "condition",
                "condition_group",
                "body_site_group",
            ]
        ]
    )
    profiles = profiles[profiles["body_site_group"].eq(FECES)].copy()
    records: list[dict[str, object]] = []
    for (condition, study_id), group in profiles.groupby(["condition", "study_id"], sort=False):
        direction_means = group.groupby("inc_dec")[trait_ids].mean()
        if not {"inc", "dec"}.issubset(direction_means.index):
            continue
        counts = {
            "increased_profile_count": int(group["inc_dec"].eq("inc").sum()),
            "decreased_profile_count": int(group["inc_dec"].eq("dec").sum()),
        }
        deltas: dict[str, float] = {}
        for trait_id in trait_ids:
            delta = float(
                direction_means.loc["inc", trait_id] - direction_means.loc["dec", trait_id]
            )
            deltas[trait_names[trait_id]] = delta
            records.append(
                {
                    "condition": condition,
                    "condition_group": group["condition_group"].iloc[0],
                    "study_id": study_id,
                    "trait_id": trait_id,
                    "trait": trait_names[trait_id],
                    "contrast": delta,
                    "contrast_definition": CONTRAST_DEFINITION,
                    **counts,
                }
            )
        increased_side = float(np.mean([deltas[name] for name in INCREASED_SIDE_TRAITS]))
        for trait, contrast, definition in [
            (
                THREE_TRAIT_COMPOSITE,
                increased_side,
                "mean increased-minus-decreased z for three increased-side traits",
            ),
            (
                REPORTER_COMPOSITE,
                increased_side - deltas[CELLULOSE],
                "three-trait composite minus cellulose increased-minus-decreased z",
            ),
        ]:
            records.append(
                {
                    "condition": condition,
                    "condition_group": group["condition_group"].iloc[0],
                    "study_id": study_id,
                    "trait_id": "composite",
                    "trait": trait,
                    "contrast": contrast,
                    "contrast_definition": definition,
                    **counts,
                }
            )
    study_contrasts = pd.DataFrame.from_records(records)
    if study_contrasts.empty:
        raise ValueError("No fecal studies contain both increased and decreased profiles")

    summary_records: list[dict[str, object]] = []
    for (condition, trait), group in study_contrasts.groupby(["condition", "trait"]):
        values = group["contrast"]
        ci_low, ci_high = _confidence_interval(values)
        positive_count = int(values.gt(0).sum())
        negative_count = int(values.lt(0).sum())
        summary_records.append(
            {
                "condition": condition,
                "condition_group": group["condition_group"].iloc[0],
                "trait": trait,
                "study_weighted_mean_contrast": float(values.mean()),
                "study_weighted_median_contrast": float(values.median()),
                "ci_low": ci_low,
                "ci_high": ci_high,
                "ci_excludes_zero": bool(ci_low > 0 or ci_high < 0),
                "paired_study_count": int(group["study_id"].nunique()),
                "positive_study_count": positive_count,
                "negative_study_count": negative_count,
                "zero_study_count": int(values.eq(0).sum()),
                "majority_direction_fraction": max(positive_count, negative_count) / len(group),
                "increased_profile_count": int(group["increased_profile_count"].sum()),
                "decreased_profile_count": int(group["decreased_profile_count"].sum()),
                "contrast_definition": group["contrast_definition"].iloc[0],
            }
        )
    return pd.DataFrame.from_records(summary_records), study_contrasts


def build_pd_taxon_carriers(
    metadata: pd.DataFrame,
    zscores: pd.DataFrame,
    traits: pd.DataFrame,
    taxa: pd.DataFrame,
    taxa_trait_matrix: pd.DataFrame,
) -> pd.DataFrame:
    """Count PD fecal profiles and studies reporting each taxon annotated with each trait.

    A taxon carries a trait when the taxa-trait matrix gives it a positive count, which
    includes annotations rolled up from its descendants. Each taxon counts once per
    profile. ``positive_side_in_lollipop`` marks the traits whose PD fecal mean z is higher
    in increased than in decreased profiles.
    """
    trait_ids = traits["trait_id"].tolist()
    direction_summary = build_pd_trait_direction_summary(metadata, zscores, traits)
    positive_ids = set(
        direction_summary.loc[direction_summary["increased_minus_decreased_mean_z"] > 0, "trait_id"]
    )
    profiles = metadata[
        metadata["condition"].eq(PD_CONDITION) & metadata["body_site_group"].eq(FECES)
    ].copy()
    profiles["direction"] = profiles["inc_dec"].map(DIRECTION_LABELS)
    occurrences = taxa.merge(
        profiles[["study_id", "experiment_id", "direction", "profile_id"]],
        on=["study_id", "experiment_id", "direction"],
        how="inner",
        validate="many_to_one",
    ).drop_duplicates(["profile_id", "taxon_id"])
    observed = occurrences.merge(
        taxa_trait_matrix[["taxon_id", "taxon_name", "taxonomic_rank", *trait_ids]],
        on="taxon_id",
        how="left",
        suffixes=("_reported", ""),
        validate="many_to_one",
        indicator=True,
    )
    missing_ids = observed.loc[observed["_merge"].ne("both"), "taxon_id"].unique()
    if len(missing_ids):
        raise ValueError(f"Observed taxa missing from trait matrix: {missing_ids.tolist()}")
    for column in ["taxon_name", "taxonomic_rank"]:
        observed[column] = observed[column].fillna(observed[f"{column}_reported"])
    records: list[dict[str, object]] = []
    for trait in traits.itertuples(index=False):
        carriers = observed[pd.to_numeric(observed[trait.trait_id]) > 0].copy()
        carriers["positive_z"] = carriers["profile_id"].map(zscores[trait.trait_id]) > 0
        for (taxon_id, name, rank), group in carriers.groupby(
            ["taxon_id", "taxon_name", "taxonomic_rank"],
            dropna=False,
        ):
            record: dict[str, object] = {
                "trait_id": trait.trait_id,
                "trait": trait.trait_label,
                "positive_side_in_lollipop": trait.trait_id in positive_ids,
                "taxon_id": taxon_id,
                "taxon_name": name,
                "taxonomic_rank": rank,
                "annotation_count": int(group[trait.trait_id].iloc[0]),
            }
            for direction in ["increased", "decreased"]:
                subset = group[group["direction"].eq(direction)]
                positive = subset[subset["positive_z"]]
                record[f"{direction}_profiles"] = subset["profile_id"].nunique()
                record[f"{direction}_studies"] = subset["study_id"].nunique()
                record[f"{direction}_positive_z_profiles"] = positive["profile_id"].nunique()
                record[f"{direction}_positive_z_studies"] = positive["study_id"].nunique()
            records.append(record)
    return pd.DataFrame(records).sort_values(
        ["trait", "increased_positive_z_profiles", "increased_profiles", "taxon_name"],
        ascending=[True, False, False, True],
    )


def prepare_pd_taxon_panel(
    carriers: pd.DataFrame,
    trait_order: list[str],
    max_taxa: int = MAX_PANEL_TAXA,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select supported genus/species labels and create a complete taxon-trait grid.

    ``total_studies`` adds the increased and decreased study counts, so a study that
    reports a taxon in both directions counts twice; it sizes and ranks the panel.
    """
    selected = carriers[
        carriers["trait"].isin(trait_order) & carriers["taxonomic_rank"].isin(ALLOWED_TAXON_RANKS)
    ].copy()
    if selected.duplicated(["taxon_id", "trait"]).any():
        raise ValueError("Duplicate taxon/trait rows in PD carrier evidence")
    selected["total_studies"] = selected["increased_studies"] + selected["decreased_studies"]
    ranking = (
        selected.groupby(["taxon_id", "taxon_name", "taxonomic_rank"], as_index=False)
        .agg(
            annotated_trait_count=("trait", "nunique"),
            total_studies=("total_studies", "max"),
            increased_studies=("increased_studies", "max"),
            decreased_studies=("decreased_studies", "max"),
        )
        .assign(
            selection_score=lambda frame: (
                frame["total_studies"] * np.sqrt(frame["annotated_trait_count"])
            )
        )
        .sort_values(
            ["selection_score", "total_studies", "taxon_name"],
            ascending=[False, False, True],
        )
        .head(max_taxa)
        .reset_index(drop=True)
    )
    grid = pd.MultiIndex.from_product(
        [ranking["taxon_id"], trait_order], names=["taxon_id", "trait"]
    )
    cells = selected.set_index(["taxon_id", "trait"]).reindex(grid).reset_index()
    lookup = ranking.set_index("taxon_id")
    for column in ["taxon_name", "taxonomic_rank"]:
        cells[column] = cells["taxon_id"].map(lookup[column])
    cells["trait_annotated"] = cells["annotation_count"].notna()
    cells["direction_log2_ratio"] = np.nan
    annotated = cells["trait_annotated"]
    cells.loc[annotated, "direction_log2_ratio"] = np.log2(
        (cells.loc[annotated, "increased_studies"].astype(float) + 0.5)
        / (cells.loc[annotated, "decreased_studies"].astype(float) + 0.5)
    )
    cells["inc_dec_study_counts"] = ""
    cells.loc[annotated, "inc_dec_study_counts"] = (
        cells.loc[annotated, "increased_studies"].astype(int).astype(str)
        + "/"
        + cells.loc[annotated, "decreased_studies"].astype(int).astype(str)
    )
    cells["blank_reason"] = np.where(
        cells["trait_annotated"], "", "no annotation in local taxa-trait matrix"
    )
    return cells, ranking


def qualifying_conditions(summary: pd.DataFrame) -> pd.DataFrame:
    """Select conditions with a nonzero 95% interval for one of the six traits."""
    required = {
        "condition",
        "trait",
        "paired_study_count",
        "ci_excludes_zero",
    }
    missing = required.difference(summary.columns)
    if missing:
        raise ValueError(f"Condition/trait summary is missing columns: {sorted(missing)}")

    six_trait_rows = summary[summary["trait"].isin(SIX_TRAITS)].copy()
    six_trait_rows["paired_study_count"] = pd.to_numeric(
        six_trait_rows["paired_study_count"], errors="coerce"
    )
    six_trait_rows["ci_excludes_zero"] = _as_bool(six_trait_rows["ci_excludes_zero"])
    six_trait_rows = six_trait_rows[six_trait_rows["paired_study_count"].ge(MIN_PAIRED_STUDIES)]

    records: list[dict[str, object]] = []
    for condition, rows in six_trait_rows.groupby("condition", sort=True):
        significant = rows[rows["ci_excludes_zero"]]
        if significant.empty:
            continue
        significant_traits = [trait for trait in SIX_TRAITS if trait in set(significant["trait"])]
        records.append(
            {
                "condition": condition,
                "qualifying_significant_trait_count": len(significant_traits),
                "qualifying_significant_traits": "; ".join(significant_traits),
                "selection_rule": (
                    "at least one of six individual traits has a descriptive 95% "
                    "study-level t interval excluding zero"
                ),
            }
        )
    output = pd.DataFrame.from_records(records)
    if output.empty:
        raise ValueError("No condition has a qualifying six-trait interval")
    output["is_parkinsons_disease"] = output["condition"].eq(PD_CONDITION)
    return output.sort_values(
        ["is_parkinsons_disease", "qualifying_significant_trait_count", "condition"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def _wrap_label(value: object, width: int = 18) -> str:
    """Wrap one plot label without breaking words."""
    return "\n".join(wrap(str(value), width=width, break_long_words=False))


def _summary_matrix(
    summary: pd.DataFrame,
    selected_conditions: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return aligned effect, support, and interval-significance matrices."""
    selected = summary[
        summary["condition"].isin(selected_conditions) & summary["trait"].isin(PANEL_A_TRAITS)
    ].copy()
    selected["study_weighted_mean_contrast"] = pd.to_numeric(
        selected["study_weighted_mean_contrast"], errors="coerce"
    )
    selected["paired_study_count"] = pd.to_numeric(selected["paired_study_count"], errors="coerce")
    selected["ci_excludes_zero"] = _as_bool(selected["ci_excludes_zero"])
    effects = selected.pivot(
        index="condition",
        columns="trait",
        values="study_weighted_mean_contrast",
    ).reindex(index=selected_conditions, columns=PANEL_A_TRAITS)
    support = selected.pivot(
        index="condition",
        columns="trait",
        values="paired_study_count",
    ).reindex(index=selected_conditions, columns=PANEL_A_TRAITS)
    significant = selected.pivot(
        index="condition",
        columns="trait",
        values="ci_excludes_zero",
    ).reindex(index=selected_conditions, columns=PANEL_A_TRAITS)
    return effects, support, significant.fillna(False).astype(bool)


def _plot_condition_heatmap(
    ax: Axes,
    summary: pd.DataFrame,
    selected_conditions: list[str],
) -> None:
    """Plot qualifying condition fingerprints across the six traits and composite."""
    effects, support, significant = _summary_matrix(summary, selected_conditions)
    limit = max(float(np.nanmax(np.abs(effects.to_numpy(dtype=float)))), 0.1)
    image = ax.imshow(
        effects.to_numpy(dtype=float),
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit),
        aspect="auto",
    )
    for y, condition in enumerate(selected_conditions):
        for x, trait in enumerate(PANEL_A_TRAITS):
            effect = effects.loc[condition, trait]
            study_count = support.loc[condition, trait]
            if pd.isna(effect):
                continue
            ax.text(
                x,
                y,
                f"{float(str(effect)):+.2f}\nn={int(float(str(study_count)))}",
                ha="center",
                va="center",
                fontsize=6.5,
                color="white" if abs(float(str(effect))) >= 0.6 * limit else SLATE_900,
            )
            if bool(significant.loc[condition, trait]):
                ax.add_patch(
                    Rectangle(
                        (x - 0.48, y - 0.48),
                        0.96,
                        0.96,
                        fill=False,
                        edgecolor=SLATE_900,
                        linewidth=1.35,
                    )
                )
        if condition == PD_CONDITION:
            ax.add_patch(
                Rectangle(
                    (-0.54, y - 0.53),
                    len(PANEL_A_TRAITS) + 0.08,
                    1.06,
                    fill=False,
                    edgecolor=PD_COLOR,
                    linewidth=2.0,
                    clip_on=False,
                )
            )
    ax.set_xticks(
        np.arange(len(PANEL_A_TRAITS)),
        [_wrap_label(trait, 14) for trait in PANEL_A_TRAITS],
    )
    ax.set_yticks(
        np.arange(len(selected_conditions)),
        [_wrap_label(condition, 24) for condition in selected_conditions],
    )
    ax.tick_params(axis="x", labelsize=7, length=0)
    ax.tick_params(axis="y", labelsize=7.5, length=0)
    ax.set_xticks(np.arange(-0.5, len(PANEL_A_TRAITS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(selected_conditions), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.spines[:].set_visible(False)
    ax.set_title(
        "A  Paired trait contrasts of the qualifying conditions",
        loc="left",
        fontsize=12,
        fontweight="bold",
        color=SLATE_900,
        pad=28,
    )
    ax.text(
        0.0,
        1.02,
        "Cell: mean study contrast and paired-study count; dark border: descriptive 95% t "
        "interval excludes 0",
        transform=ax.transAxes,
        fontsize=7.5,
        color=SLATE_500,
        va="bottom",
    )
    colorbar = ax.figure.colorbar(image, ax=ax, fraction=0.035, pad=0.025)
    colorbar.set_label("Mean increased - decreased z", fontsize=7.5)
    colorbar.ax.tick_params(labelsize=7)


def _specificity_coordinates(
    summary: pd.DataFrame,
    selected_conditions: list[str],
) -> pd.DataFrame:
    """Build the two composite coordinates for qualifying conditions."""
    subset = summary[summary["condition"].isin(selected_conditions)].copy()
    cellulose = subset[subset["trait"].eq(CELLULOSE)][
        ["condition", "study_weighted_mean_contrast", "paired_study_count"]
    ].rename(
        columns={
            "study_weighted_mean_contrast": "cellulose_contrast",
            "paired_study_count": "paired_study_count",
        }
    )
    increased_side = subset[subset["trait"].eq(THREE_TRAIT_COMPOSITE)][
        ["condition", "study_weighted_mean_contrast"]
    ].rename(columns={"study_weighted_mean_contrast": "increased_side_composite"})
    coordinates = cellulose.merge(
        increased_side,
        on="condition",
        how="inner",
        validate="one_to_one",
    )
    for column in ["cellulose_contrast", "paired_study_count", "increased_side_composite"]:
        coordinates[column] = pd.to_numeric(coordinates[column], errors="coerce")
    coordinates["cellulose_decreased_skew"] = -coordinates["cellulose_contrast"]
    return coordinates.dropna().reset_index(drop=True)


def reporter_numbers(summary: pd.DataFrame, selection: pd.DataFrame) -> str:
    """Render the figure note's numbers from the summary and condition-selection tables."""
    selected = selection["condition"].astype(str).tolist()
    coordinates = _specificity_coordinates(summary, selected).set_index("condition")
    upper_right = coordinates[
        coordinates["cellulose_decreased_skew"].gt(0)
        & coordinates["increased_side_composite"].gt(0)
    ]
    rows = summary.set_index(["condition", "trait"])

    def contrast(condition: str, trait: str) -> float:
        return float(rows.at[(condition, trait), "study_weighted_mean_contrast"])

    def interval(trait: str) -> str:
        low = float(rows.at[(PD_CONDITION, trait), "ci_low"])
        high = float(rows.at[(PD_CONDITION, trait), "ci_high"])
        return f"`{low:+.3f}` to `{high:+.3f}`"

    others = [condition for condition in upper_right.index if condition != PD_CONDITION]
    larger_reporter = sum(
        contrast(condition, REPORTER_COMPOSITE) > contrast(PD_CONDITION, REPORTER_COMPOSITE)
        for condition in others
    )
    larger_three_trait = sum(
        contrast(condition, THREE_TRAIT_COMPOSITE) > contrast(PD_CONDITION, THREE_TRAIT_COMPOSITE)
        for condition in others
    )
    studies = int(rows.at[(PD_CONDITION, CELLULOSE), "paired_study_count"])
    pd_in_quadrant = "PD among them" if PD_CONDITION in upper_right.index else "PD not among them"
    return "\n".join(
        [
            f"- Qualifying fecal conditions in panels A and B: `{len(selected)}`.",
            f"- In panel B's upper-right quadrant: `{len(upper_right)}` of the `{len(selected)}`, "
            f"{pd_in_quadrant}. Of the other `{len(others)}`, `{larger_reporter}` have a larger "
            f"ecological-reporter composite than PD and `{larger_three_trait}` a larger "
            "three-trait composite.",
            f"- PD, across `{studies}` paired studies: cellulose contrast "
            f"`{contrast(PD_CONDITION, CELLULOSE):+.3f}` (descriptive 95% t interval "
            f"{interval(CELLULOSE)}); three-trait increased-side composite "
            f"`{contrast(PD_CONDITION, THREE_TRAIT_COMPOSITE):+.3f}` "
            f"({interval(THREE_TRAIT_COMPOSITE)}); ecological-reporter composite "
            f"`{contrast(PD_CONDITION, REPORTER_COMPOSITE):+.3f}` "
            f"({interval(REPORTER_COMPOSITE)}).",
        ]
    )


def write_reporter_numbers(note_path: Path, numbers: str) -> None:
    """Replace the generated numbers block of the figure note."""
    text = note_path.read_text(encoding="utf-8")
    start, end = text.find(NUMBERS_BEGIN), text.find(NUMBERS_END)
    if start < 0 or end < start:
        raise ValueError(f"{note_path.name} has no generated numbers block")
    head = text[: start + len(NUMBERS_BEGIN)]
    note_path.write_text(f"{head}\n{numbers}\n{text[end:]}", encoding="utf-8")


def _label_overlap(box: Bbox, obstacles: list[Bbox], frame: Bbox) -> float:
    """Return the area a label shares with the obstacles plus its area outside the frame."""
    shared = [Bbox.intersection(box, other) for other in obstacles]
    overlap = sum(part.width * part.height for part in shared if part is not None)
    inside = Bbox.intersection(box, frame)
    inside_area = inside.width * inside.height if inside is not None else 0.0
    return overlap + box.width * box.height - inside_area


def _place_point_labels(ax: Axes, coordinates: pd.DataFrame, sizes: np.ndarray) -> None:
    """Label each point at the first offset that overlaps no other label or point.

    PD is placed first, then conditions with more paired studies. When every offset
    collides, the label goes where it overlaps least; far labels get a leader line.
    """
    ax.set_xlim(ax.get_xlim())
    ax.set_ylim(ax.get_ylim())
    points = coordinates[["cellulose_decreased_skew", "increased_side_composite"]].to_numpy(
        dtype=float
    )
    centers = ax.transData.transform(points)
    radii = np.sqrt(sizes) / 2 * ax.figure.dpi / 72
    markers = [
        Bbox.from_extents(x - radius, y - radius, x + radius, y + radius)
        for (x, y), radius in zip(centers, radii, strict=True)
    ]
    frame = ax.get_window_extent()
    conditions = coordinates["condition"].astype(str).tolist()
    order = sorted(
        range(len(conditions)),
        key=lambda index: (
            conditions[index] != PD_CONDITION,
            -float(coordinates["paired_study_count"].iloc[index]),
        ),
    )
    placed: list[Bbox] = []
    for index in order:
        condition = conditions[index]
        style = {
            "textcoords": "offset points",
            "fontsize": 6.8,
            "color": PD_COLOR if condition == PD_CONDITION else SLATE_700,
            "fontweight": "bold" if condition == PD_CONDITION else "normal",
        }
        obstacles = [*placed, *(box for other, box in enumerate(markers) if other != index)]
        candidates = []
        for dx, dy in LABEL_OFFSETS:
            alignment = {
                "ha": "left" if dx > 0 else "right" if dx < 0 else "center",
                "va": "bottom" if dy > 0 else "top" if dy < 0 else "center",
            }
            label = ax.annotate(
                condition, xy=tuple(points[index]), xytext=(dx, dy), **alignment, **style
            )
            box = label.get_window_extent().padded(1.5)
            label.remove()
            candidates.append((_label_overlap(box, obstacles, frame), (dx, dy), alignment, box))
        _overlap, (dx, dy), alignment, box = min(candidates, key=lambda candidate: candidate[0])
        leader = (
            {"arrowstyle": "-", "color": SLATE_300, "linewidth": 0.6, "shrinkA": 0, "shrinkB": 3}
            if math.hypot(dx, dy) > 12
            else None
        )
        ax.annotate(
            condition,
            xy=tuple(points[index]),
            xytext=(dx, dy),
            arrowprops=leader,
            **alignment,
            **style,
        )
        placed.append(box)


def _plot_specificity(
    ax: Axes,
    summary: pd.DataFrame,
    selected_conditions: list[str],
) -> None:
    """Plot the cellulose and three-trait composites of the qualifying conditions."""
    coordinates = _specificity_coordinates(summary, selected_conditions)
    sizes = np.clip(24.0 + 2.2 * coordinates["paired_study_count"], 30.0, 135.0)
    colors = [
        PD_COLOR if condition == PD_CONDITION else "#94a3b8"
        for condition in coordinates["condition"].astype(str)
    ]
    edges = [
        SLATE_900 if condition == PD_CONDITION else "white"
        for condition in coordinates["condition"].astype(str)
    ]
    ax.scatter(
        coordinates["cellulose_decreased_skew"],
        coordinates["increased_side_composite"],
        s=sizes,
        c=colors,
        edgecolors=edges,
        linewidths=1.0,
        alpha=0.95,
        zorder=3,
    )
    ax.axhline(0, color=SLATE_500, linewidth=0.9)
    ax.axvline(0, color=SLATE_500, linewidth=0.9)
    _place_point_labels(ax, coordinates, sizes.to_numpy(dtype=float))
    ax.grid(color=SLATE_200, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.set_xlabel(
        "Cellulose decreased-profile skew\n-(increased - decreased mean z)",
        fontsize=8.5,
    )
    ax.set_ylabel(
        "Three-trait increased-side composite\nmean(increased - decreased z)",
        fontsize=8.5,
    )
    ax.set_title(
        "B  PD among the qualifying fecal conditions",
        loc="left",
        fontsize=12,
        fontweight="bold",
        color=SLATE_900,
        pad=28,
    )
    ax.text(
        0.0,
        1.02,
        (
            f"{len(coordinates)} conditions, each with a descriptive 95% interval excluding 0 "
            "for >=1 of the six traits; point area: paired studies"
        ),
        transform=ax.transAxes,
        fontsize=7.5,
        color=SLATE_500,
        va="bottom",
    )


def _first_nonmissing(values: pd.Series) -> object:
    """Return the first nonempty value from an audit column."""
    for value in values:
        if pd.notna(value) and str(value).strip():
            return value
    return pd.NA


def _plot_taxon_panel(ax: Axes, taxon_panel: pd.DataFrame) -> None:
    """Plot PD genus/species annotation and directional-support provenance."""
    frame = taxon_panel.copy()
    frame["trait_annotated"] = _as_bool(frame["trait_annotated"])
    frame["total_studies"] = pd.to_numeric(frame["total_studies"], errors="coerce")
    frame["direction_log2_ratio"] = pd.to_numeric(frame["direction_log2_ratio"], errors="coerce")
    taxon_ids = frame["taxon_id"].drop_duplicates().tolist()
    metadata = (
        frame.groupby("taxon_id", sort=False)
        .agg(
            taxon_name=("taxon_name", "first"),
            taxonomic_rank=("taxonomic_rank", "first"),
            inc_dec_study_counts=("inc_dec_study_counts", _first_nonmissing),
        )
        .reindex(taxon_ids)
    )
    maximum_ratio = max(float(frame["direction_log2_ratio"].abs().max()), 0.1)
    norm = TwoSlopeNorm(vmin=-maximum_ratio, vcenter=0.0, vmax=maximum_ratio)
    cmap = plt.get_cmap("RdBu_r")
    for y, taxon_id in enumerate(taxon_ids):
        rows = frame[frame["taxon_id"].eq(taxon_id)].set_index("trait")
        for x, trait in enumerate(SIX_TRAITS):
            if trait not in rows.index or not bool(rows.loc[trait, "trait_annotated"]):
                ax.scatter(x, y, marker="x", s=20, color=SLATE_300, linewidth=0.8)
                continue
            ratio = float(str(rows.at[trait, "direction_log2_ratio"]))
            studies = float(str(rows.at[trait, "total_studies"]))
            ax.scatter(
                x,
                y,
                s=np.clip(22.0 + 6.5 * studies, 28.0, 140.0),
                color=cmap(norm(ratio)),
                edgecolor="white",
                linewidth=0.7,
                zorder=3,
            )
    labels = []
    for metadata_row in metadata.itertuples():
        rank = str(metadata_row.taxonomic_rank)[:1].upper()
        support = str(metadata_row.inc_dec_study_counts)
        labels.append(f"{metadata_row.taxon_name} [{rank}]  {support} studies")
    ax.set_xticks(np.arange(len(SIX_TRAITS)), [_wrap_label(trait, 17) for trait in SIX_TRAITS])
    ax.set_yticks(np.arange(len(taxon_ids)), labels)
    ax.invert_yaxis()
    ax.set_xlim(-0.6, len(SIX_TRAITS) - 0.4)
    ax.set_xticks(np.arange(-0.5, len(SIX_TRAITS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(taxon_ids), 1), minor=True)
    ax.grid(which="minor", color=SLATE_200, linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="x", labelsize=7.5, length=0)
    ax.tick_params(axis="y", labelsize=7, length=0)
    ax.spines[:].set_visible(False)
    ax.set_title(
        "C  PD reported taxa: KG-Microbe annotations x BugSigDB direction support",
        loc="left",
        fontsize=12,
        fontweight="bold",
        color=SLATE_900,
        pad=28,
    )
    ax.text(
        0.0,
        1.02,
        "Circle: annotated trait; size: increased + decreased study reports; "
        "x: no local annotation",
        transform=ax.transAxes,
        fontsize=7.5,
        color=SLATE_500,
        va="bottom",
    )
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    colorbar = ax.figure.colorbar(scalar, ax=ax, fraction=0.018, pad=0.018)
    colorbar.set_label(
        "log2((increased studies + 0.5) /\n(decreased studies + 0.5))",
        fontsize=7,
    )
    colorbar.ax.tick_params(labelsize=7)


def plot_combined_figure(
    summary: pd.DataFrame,
    taxon_panel: pd.DataFrame,
    selection: pd.DataFrame,
    output_path: Path,
) -> None:
    """Render the condition-filtered three-panel ecological-reporter figure."""
    selected_conditions = selection["condition"].astype(str).tolist()
    figure_height = max(13.5, 8.0 + 0.36 * len(selected_conditions))
    fig = plt.figure(figsize=(18, figure_height), facecolor="white")
    # Fix the layout first: panel B places its labels against its final position.
    fig.subplots_adjust(top=0.9, bottom=0.055, left=0.135, right=0.965)
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.15, 0.9],
        width_ratios=[1.32, 1.0],
        hspace=0.38,
        wspace=0.27,
    )
    ax_heatmap = fig.add_subplot(grid[0, 0])
    ax_specificity = fig.add_subplot(grid[0, 1])
    ax_taxa = fig.add_subplot(grid[1, :])
    _plot_condition_heatmap(ax_heatmap, summary, selected_conditions)
    _plot_specificity(ax_specificity, summary, selected_conditions)
    _plot_taxon_panel(ax_taxa, taxon_panel)
    fig.suptitle(
        "Parkinson's disease ecological-reporter pattern across fecal conditions, "
        "with taxon provenance",
        x=0.04,
        y=0.985,
        ha="left",
        fontsize=18,
        fontweight="bold",
        color=SLATE_900,
    )
    fig.text(
        0.04,
        0.958,
        (
            "Panels A-B retain only conditions with a descriptive 95% study-level interval "
            "excluding zero for >=1 of the six traits; composites do not qualify a condition. "
            "Contrasts summarize KG-Microbe annotations of taxa in curated BugSigDB "
            "signatures, not measured activity."
        ),
        ha="left",
        va="top",
        fontsize=9,
        color=SLATE_700,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    svg_path = output_path.with_suffix(".svg")
    # A fixed salt gives the SVG's clip-path and image ids the same values on every run.
    with plt.rc_context({"svg.hashsalt": STEM}):
        fig.savefig(
            svg_path,
            bbox_inches="tight",
            facecolor="white",
            metadata={"Date": None},
        )
    plt.close(fig)
    svg_lines = svg_path.read_text(encoding="utf-8").splitlines()
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_lines) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    """Build the audit tables from the profile inputs and draw the combined figure."""
    traits = load_trait_summary()
    trait_ids = traits["trait_id"].tolist()
    missing_traits = set(SIX_TRAITS) - set(traits["trait_label"])
    if missing_traits:
        raise ValueError(f"Figure traits missing from the trait summary: {sorted(missing_traits)}")
    metadata = load_profile_metadata()
    zscores = build_zscore_matrix(metadata, trait_ids)
    from mtp.population_tables import is_canonical_derived_path, orient_and_filter_taxa_df

    summary, study_contrasts = build_condition_trait_summaries(metadata, zscores, traits)
    taxa = read_tsv(TAXA_PATH)
    if is_canonical_derived_path(TAXA_PATH, "bugsigdb_taxa.tsv"):
        taxa = orient_and_filter_taxa_df(taxa, filter_to_selection=True)
    taxa_trait_matrix = pd.read_csv(
        TAXA_TRAIT_MATRIX_PATH,
        sep="\t",
        usecols=["taxon_id", "taxon_name", "taxonomic_rank", *trait_ids],
    )
    carriers = build_pd_taxon_carriers(metadata, zscores, traits, taxa, taxa_trait_matrix)
    taxon_panel, _ranking = prepare_pd_taxon_panel(carriers, SIX_TRAITS)
    write_tsv(summary, SUMMARY_PATH)
    write_tsv(study_contrasts, STUDY_CONTRASTS_PATH)
    write_tsv(taxon_panel, TAXON_PANEL_PATH)
    # Select and draw from the tables as written, so the figure shows what they record.
    summary = read_tsv(SUMMARY_PATH)
    selection = qualifying_conditions(summary)
    write_tsv(selection, SELECTION_PATH)
    plot_combined_figure(summary, read_tsv(TAXON_PANEL_PATH), selection, OUTPUT_PATH)
    write_reporter_numbers(FIGURES_NOTE_PATH, reporter_numbers(summary, read_tsv(SELECTION_PATH)))
    print(
        f"Condition summaries: {len(summary)}; study contrasts: {len(study_contrasts)}; "
        f"taxon-panel cells: {len(taxon_panel)}"
    )
    print(
        f"Wrote {OUTPUT_PATH.name} with {len(selection)} qualifying conditions "
        f"and {len(SIX_TRAITS)} individual traits"
    )


if __name__ == "__main__":
    main()
