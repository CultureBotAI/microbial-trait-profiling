"""Render the audited PD ecological-reporter specificity figure."""

from __future__ import annotations

from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "trait-visualization" / "trait-profile-investigation" / "PD-cellulose"
STEM = "pd_ecological_reporter_combined"
SUMMARY_PATH = OUTPUT_DIR / f"{STEM}_condition_trait_summary.tsv"
TAXON_PANEL_PATH = OUTPUT_DIR / f"{STEM}_taxon_panel.tsv"
SELECTION_PATH = OUTPUT_DIR / f"{STEM}_condition_selection.tsv"
OUTPUT_PATH = OUTPUT_DIR / f"{STEM}.png"

PD_CONDITION = "Parkinson's disease"
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
THREE_TRAIT_COMPOSITE = "three-trait increased-side composite"
REPORTER_COMPOSITE = "ecological-reporter composite"
PANEL_A_TRAITS = [*SIX_TRAITS, REPORTER_COMPOSITE]

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
                color=SLATE_900,
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
        "A  Disease fingerprints for conditions with a six-trait change",
        loc="left",
        fontsize=12,
        fontweight="bold",
        color=SLATE_900,
        pad=28,
    )
    ax.text(
        0.0,
        1.02,
        "Cell: mean study contrast and paired-study count; dark border: 95% interval excludes 0",
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
    cellulose = subset[subset["trait"].eq("cellulose catabolic process")][
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


def _plot_specificity(
    ax: Axes,
    summary: pd.DataFrame,
    selected_conditions: list[str],
) -> None:
    """Plot specificity using only conditions selected by the six-trait rule."""
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
    offsets = [(5, 5), (5, -10), (-5, 6), (-5, -11)]
    ordered = coordinates.sort_values(
        ["condition"],
        key=lambda values: values.ne(PD_CONDITION),
    )
    for index, row in enumerate(ordered.itertuples(index=False)):
        condition = str(row.condition)
        offset = offsets[index % len(offsets)]
        ax.annotate(
            condition,
            xy=(row.cellulose_decreased_skew, row.increased_side_composite),
            xytext=offset,
            textcoords="offset points",
            ha="left" if offset[0] > 0 else "right",
            va="bottom" if offset[1] > 0 else "top",
            fontsize=6.8,
            color=PD_COLOR if condition == PD_CONDITION else SLATE_700,
            fontweight="bold" if condition == PD_CONDITION else "normal",
        )
    ax.axhline(0, color=SLATE_500, linewidth=0.9)
    ax.axvline(0, color=SLATE_500, linewidth=0.9)
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
        "B  Specificity among qualifying fecal conditions",
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
            f"{len(coordinates)} conditions; each has >=1 significant change among the six traits; "
            "point area: paired studies"
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
        "Circle: annotated trait; size: distinct-study support; x: no local annotation",
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
        "Parkinson's disease ecological-reporter hypothesis: specificity and provenance",
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
            "excluding zero for >=1 of the six traits; composites do not qualify a condition."
        ),
        ha="left",
        va="top",
        fontsize=9,
        color=SLATE_700,
    )
    fig.subplots_adjust(top=0.9, bottom=0.055, left=0.135, right=0.965)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    svg_path = output_path.with_suffix(".svg")
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
    """Filter qualifying conditions and redraw the audited combined figure."""
    summary = read_tsv(SUMMARY_PATH)
    taxon_panel = read_tsv(TAXON_PANEL_PATH)
    selection = qualifying_conditions(summary)
    write_tsv(selection, SELECTION_PATH)
    plot_combined_figure(summary, taxon_panel, selection, OUTPUT_PATH)
    print(
        f"Wrote {OUTPUT_PATH.name} with {len(selection)} qualifying conditions "
        f"and {len(SIX_TRAITS)} individual traits"
    )


if __name__ == "__main__":
    main()
