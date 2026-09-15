"""Plot species-level influence contrasts for all six Parkinson's disease traits."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle
from scipy import stats

from create_pd_cellulose_analysis import (
    FECES,
    OUTPUT_DIR,
    PD_CONDITION,
    PROFILE_MATRIX_PATH,
    TAXA_TRAIT_MATRIX_PATH,
    load_profile_metadata,
    read_tsv,
    write_tsv,
)

PROFILE_COLUMNS = ["study_id", "experiment_id", "inc_dec"]
STEM = "pd_species_six_traits_influence_heatmap"
MIN_SUPPORT_STUDIES = 3
TRANSFORM_SCALE = 1000.0
ROW_COLOR_LIMIT = 3.0
INCREASED_COLOR = "#b91c1c"
DECREASED_COLOR = "#1d4ed8"
MISSING_COLOR = "#e5e7eb"


def profile_id(study_id: object, experiment_id: object, inc_dec: object) -> str:
    """Return the profile identifier shared by the matrix and metadata exports."""
    return f"{study_id}|{experiment_id}|{inc_dec}"


def read_profile_counts(path: Path) -> tuple[pd.DataFrame, list[str]]:
    """Read the complete non-MediaDive count matrix with one row per profile."""
    columns = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    trait_columns = [
        column
        for column in columns
        if column not in PROFILE_COLUMNS and "mediadive.medium" not in column
    ]
    profiles = pd.read_csv(path, sep="\t", usecols=[*PROFILE_COLUMNS, *trait_columns])
    ids = [profile_id(*values) for values in profiles[PROFILE_COLUMNS].itertuples(index=False)]
    values = profiles[trait_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    values.index = ids
    return values.groupby(level=0, sort=False).sum(), trait_columns


def read_taxon_vectors(path: Path, trait_columns: list[str]) -> pd.DataFrame:
    """Read reported-taxon vectors aligned to the non-MediaDive profile columns."""
    metadata = ["taxon_id", "taxon_name", "taxonomic_rank"]
    matrix = pd.read_csv(path, sep="\t", usecols=[*metadata, *trait_columns])
    if matrix["taxon_id"].duplicated().any():
        raise ValueError("Duplicate taxon IDs in taxa-trait matrix")
    matrix[trait_columns] = matrix[trait_columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return matrix.set_index("taxon_id")


def fixed_reference(
    profile_counts: pd.DataFrame,
    trait_ids: list[str],
    reference_ids: pd.Index,
) -> tuple[pd.DataFrame, pd.Series, pd.Series, pd.Series]:
    """Transform profile fractions and fit one fixed human-corpus reference."""
    missing = [trait for trait in trait_ids if trait not in profile_counts]
    if missing:
        raise ValueError(f"Missing selected trait columns: {missing}")
    totals = profile_counts.sum(axis=1).astype(float)
    fractions = (
        profile_counts[trait_ids]
        .astype(float)
        .div(totals.where(totals.ne(0), np.nan), axis=0)
        .fillna(0.0)
    )
    transformed = np.log1p(TRANSFORM_SCALE * fractions)
    reference = transformed.loc[reference_ids.intersection(transformed.index)]
    means = reference.mean(axis=0)
    standard_deviations = reference.std(axis=0, ddof=0)
    if reference.empty or standard_deviations.isna().any() or (standard_deviations <= 0).any():
        raise ValueError("The fixed human reference is empty or has a zero-variance selected trait")
    return (transformed - means) / standard_deviations, totals, means, standard_deviations


def complete_pd_pairs(metadata: pd.DataFrame, valid_profiles: pd.Index) -> pd.DataFrame:
    """Return exactly one increased and decreased PD fecal profile per retained experiment."""
    selected = metadata.loc[metadata.index.intersection(valid_profiles)].copy()
    selected = selected[
        selected["condition"].eq(PD_CONDITION)
        & selected["body_site_group"].eq(FECES)
        & selected["inc_dec"].isin(["inc", "dec"])
    ]
    if "profile_id" not in selected:
        selected["profile_id"] = selected.index
    selected = selected.reset_index(drop=True)
    keys = ["study_id", "experiment_id", "inc_dec"]
    if selected.duplicated(keys).any():
        raise ValueError("Duplicate PD study/experiment/direction profiles")
    complete = selected.groupby(["study_id", "experiment_id"])["inc_dec"].agg(set)
    complete_keys = complete[complete.map(lambda directions: directions == {"inc", "dec"})].index
    indexed = selected.set_index(["study_id", "experiment_id"])
    return indexed.loc[indexed.index.isin(complete_keys)].reset_index()


def paired_study_contrasts(scores: pd.Series, pairs: pd.DataFrame) -> pd.Series:
    """Compute increased-minus-decreased experiment contrasts and weight studies equally."""
    frame = pairs[["study_id", "experiment_id", "inc_dec", "profile_id"]].copy()
    frame["score"] = frame["profile_id"].map(scores)
    if frame["score"].isna().any():
        raise ValueError("A retained PD profile is missing its influence score")
    experiments = frame.pivot(
        index=["study_id", "experiment_id"], columns="inc_dec", values="score"
    )
    if not {"inc", "dec"}.issubset(experiments.columns):
        raise ValueError("No complete increased/decreased PD pairs were available")
    contrasts = experiments["inc"] - experiments["dec"]
    return contrasts.groupby(level="study_id").mean()


def mean_interval(values: pd.Series) -> tuple[float, float, float]:
    """Return the mean and two-sided 95% Student t interval."""
    finite = values[np.isfinite(values.to_numpy(float))].astype(float)
    if finite.empty:
        return math.nan, math.nan, math.nan
    mean = float(finite.mean())
    if len(finite) < 2:
        return mean, math.nan, math.nan
    standard_error = float(finite.std(ddof=1) / math.sqrt(len(finite)))
    if standard_error == 0:
        return mean, mean, mean
    radius = float(stats.t.ppf(0.975, len(finite) - 1) * standard_error)
    return mean, mean - radius, mean + radius


def classify_effect(mean: float, low: float, high: float, support_studies: int) -> str:
    """Classify a signed influence contrast without forcing sparse or uncertain cells."""
    if support_studies < MIN_SUPPORT_STUDIES:
        return "low support"
    if not all(math.isfinite(value) for value in [mean, low, high]):
        return "mixed"
    if low > 0:
        return "inc"
    if high < 0:
        return "dec"
    return "mixed"


def taxon_influence(
    profile_counts: pd.DataFrame,
    totals: pd.Series,
    taxon_vector: pd.Series,
    presence_ids: pd.Index,
    trait_ids: list[str],
    reference_standard_deviations: pd.Series,
) -> pd.DataFrame:
    """Recompute selected scores after removing one reported species from each observed profile."""
    present = presence_ids.intersection(profile_counts.index)
    after_counts = profile_counts.loc[present, trait_ids].copy()
    after_counts = after_counts.subtract(taxon_vector[trait_ids], axis=1)
    removed_total = float(taxon_vector[profile_counts.columns].sum())
    after_totals = totals.loc[present] - removed_total
    if (after_counts.to_numpy(float) < 0).any() or (after_totals < 0).any():
        raise ValueError("Taxon removal produced a negative profile count")
    before_fractions = (
        profile_counts.loc[present, trait_ids]
        .astype(float)
        .div(totals.loc[present].where(totals.loc[present].ne(0), np.nan), axis=0)
        .fillna(0.0)
    )
    after_fractions = (
        after_counts.astype(float)
        .div(after_totals.where(after_totals.ne(0), np.nan), axis=0)
        .fillna(0.0)
    )
    before_transformed = np.log1p(TRANSFORM_SCALE * before_fractions)
    after_transformed = np.log1p(TRANSFORM_SCALE * after_fractions)
    influence = pd.DataFrame(0.0, index=profile_counts.index, columns=trait_ids)
    influence.loc[present] = (before_transformed - after_transformed).div(
        reference_standard_deviations[trait_ids], axis=1
    )
    return influence


def build_influence_cells(
    profile_counts: pd.DataFrame,
    metadata: pd.DataFrame,
    taxon_vectors: pd.DataFrame,
    evidence: pd.DataFrame,
    trait_summary: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build six-trait species cells, row annotations, and trait metadata."""
    traits = trait_summary.sort_values("increased_minus_decreased_mean_z", ascending=False).copy()
    trait_ids = traits["trait_id"].astype(str).tolist()
    human_ids = metadata.index[metadata["host_species"].eq("Homo sapiens")]
    baseline_z, totals, _, reference_sd = fixed_reference(profile_counts, trait_ids, human_ids)
    pairs = complete_pd_pairs(metadata, totals[totals > 0].index)
    if pairs.empty:
        raise ValueError("No complete nonzero-total PD fecal pairs were available")

    overall_records = []
    for trait in traits.itertuples(index=False):
        contrasts = paired_study_contrasts(baseline_z[str(trait.trait_id)], pairs)
        mean, low, high = mean_interval(contrasts)
        overall_records.append(
            {
                "trait_id": str(trait.trait_id),
                "trait": trait.trait,
                "overall_mean_z_contrast": mean,
                "overall_ci95_low": low,
                "overall_ci95_high": high,
                "overall_direction": "inc" if mean > 0 else "dec" if mean < 0 else "mixed",
            }
        )
    trait_metadata = pd.DataFrame(overall_records)

    pd_profiles = metadata[
        metadata["condition"].eq(PD_CONDITION)
        & metadata["body_site_group"].eq(FECES)
        & metadata["inc_dec"].isin(["inc", "dec"])
    ].index
    evidence = evidence[
        evidence["trait_id"].isin(trait_ids)
        & evidence["taxonomic_rank"].eq("species")
        & evidence["profile_id"].isin(pd_profiles)
    ].copy()
    occurrence = evidence.drop_duplicates(["taxon_id", "profile_id"])
    species_ids = pd.Index(sorted(occurrence["taxon_id"].unique()))
    missing_taxa = species_ids.difference(taxon_vectors.index)
    if len(missing_taxa):
        raise ValueError(
            f"Observed species missing from taxa-trait matrix: {missing_taxa.tolist()}"
        )

    pair_profile_ids = pd.Index(pairs["profile_id"])
    cells = []
    rows = []
    for taxon_id in species_ids:
        vector = taxon_vectors.loc[taxon_id]
        taxon_occurrence = occurrence[occurrence["taxon_id"].eq(taxon_id)]
        presence_ids = pd.Index(taxon_occurrence["profile_id"].unique())
        paired_presence = taxon_occurrence[taxon_occurrence["profile_id"].isin(pair_profile_ids)]
        support_studies = int(paired_presence["study_id"].nunique())
        influence = taxon_influence(
            profile_counts,
            totals,
            vector,
            presence_ids,
            trait_ids,
            reference_sd,
        )
        effect_magnitudes = []
        for trait in trait_metadata.itertuples(index=False):
            trait_id = str(trait.trait_id)
            annotation_count = int(vector[trait_id])
            contrasts = paired_study_contrasts(influence[trait_id], pairs)
            mean, low, high = mean_interval(contrasts)
            annotated = annotation_count > 0
            status = (
                classify_effect(mean, low, high, support_studies) if annotated else "not annotated"
            )
            if annotated:
                effect_magnitudes.append(abs(mean))
            cells.append(
                {
                    "taxon_id": taxon_id,
                    "taxon_name": vector["taxon_name"],
                    "taxonomic_rank": vector["taxonomic_rank"],
                    "trait_id": trait_id,
                    "trait": trait.trait,
                    "trait_overall_direction": trait.overall_direction,
                    "trait_overall_mean_z_contrast": trait.overall_mean_z_contrast,
                    "annotation_count": annotation_count,
                    "has_trait_annotation": annotated,
                    "paired_profiles": len(pairs),
                    "paired_studies": int(pairs["study_id"].nunique()),
                    "taxon_support_profiles": int(paired_presence["profile_id"].nunique()),
                    "taxon_support_studies": support_studies,
                    "mean_influence_contrast_z": mean if annotated else math.nan,
                    "ci95_low": low if annotated else math.nan,
                    "ci95_high": high if annotated else math.nan,
                    "status": status,
                }
            )

        direction_counts = {
            direction: int(
                taxon_occurrence.loc[
                    taxon_occurrence["direction"].eq(label), "profile_id"
                ].nunique()
            )
            for direction, label in [("inc", "increased"), ("dec", "decreased")]
        }
        inc, dec = direction_counts["inc"], direction_counts["dec"]
        ratio = inc / dec if dec else math.inf
        log_ratio = math.log2(ratio) if inc and dec else math.inf if inc else -math.inf
        rows.append(
            {
                "taxon_id": taxon_id,
                "taxon_name": vector["taxon_name"],
                "taxonomic_rank": vector["taxonomic_rank"],
                "annotated_trait_count": int((vector[trait_ids].astype(float) > 0).sum()),
                "increased_profiles": inc,
                "decreased_profiles": dec,
                "inc_dec_counts": f"{inc}/{dec}",
                "inc_dec_ratio": ratio,
                "log2_inc_dec_ratio": log_ratio,
                "paired_support_studies": support_studies,
                "max_absolute_influence_z": max(effect_magnitudes, default=0.0),
            }
        )

    row_metadata = pd.DataFrame(rows).sort_values(
        [
            "annotated_trait_count",
            "paired_support_studies",
            "max_absolute_influence_z",
            "taxon_name",
        ],
        ascending=[False, False, False, True],
    )
    return pd.DataFrame(cells), row_metadata, trait_metadata


def plot_heatmap(
    cells: pd.DataFrame,
    rows: pd.DataFrame,
    traits: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot trait-specific influence cells with a separate taxon occurrence annotation."""
    row_ids = rows["taxon_id"].tolist()
    trait_ids = traits["trait_id"].tolist()
    values = cells.pivot(index="taxon_id", columns="trait_id", values="mean_influence_contrast_z")
    labels = cells.pivot(index="taxon_id", columns="trait_id", values="status")
    values = values.reindex(index=row_ids, columns=trait_ids)
    labels = labels.reindex(index=row_ids, columns=trait_ids)
    finite = np.abs(values.to_numpy(float)[np.isfinite(values.to_numpy(float))])
    color_limit = max(float(np.quantile(finite, 0.95)) if finite.size else 0.0, 0.02)

    height = max(8.5, len(rows) * 0.32 + 3.4)
    fig = plt.figure(figsize=(15.5, height), facecolor="white")
    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=[6.0, 0.85],
        left=0.29,
        right=0.92,
        top=1 - 1.65 / height,
        bottom=1.15 / height,
        wspace=0.08,
    )
    ax = fig.add_subplot(grid[0, 0])
    row_ax = fig.add_subplot(grid[0, 1], sharey=ax)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(MISSING_COLOR)
    image = ax.imshow(
        values.to_numpy(float), cmap=cmap, vmin=-color_limit, vmax=color_limit, aspect="auto"
    )
    ax.set_xticks(
        range(len(traits)),
        [
            f"{'↑' if row.overall_mean_z_contrast > 0 else '↓'} "
            + "\n".join(wrap(str(row.trait), width=17))
            + f"\n({row.overall_mean_z_contrast:+.2f})"
            for row in traits.itertuples(index=False)
        ],
        fontsize=9,
    )
    ax.set_yticks(range(len(rows)), rows["taxon_name"], fontsize=8.5)
    ax.tick_params(length=0, pad=6)
    ax.set_xticks(np.arange(-0.5, len(traits), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for y, taxon_id in enumerate(row_ids):
        for x, trait_id in enumerate(trait_ids):
            status = labels.loc[taxon_id, trait_id]
            if status == "not annotated":
                continue
            value = values.loc[taxon_id, trait_id]
            if status == "low support":
                ax.add_patch(
                    Rectangle(
                        (x - 0.5, y - 0.5),
                        1,
                        1,
                        facecolor="none",
                        edgecolor="#64748b",
                        hatch="///",
                        linewidth=0,
                    )
                )
            ax.text(
                x,
                y,
                "low\nsupport" if status == "low support" else status,
                ha="center",
                va="center",
                fontsize=7.2,
                color="white" if abs(value) >= 0.58 * color_limit else "#111827",
            )

    row_values = rows["log2_inc_dec_ratio"].clip(-ROW_COLOR_LIMIT, ROW_COLOR_LIMIT).to_numpy(float)
    row_ax.imshow(
        row_values[:, None], cmap=cmap, vmin=-ROW_COLOR_LIMIT, vmax=ROW_COLOR_LIMIT, aspect="auto"
    )
    row_ax.set_xticks([0], ["Taxon occurrence\ninc/dec"], fontsize=9)
    row_ax.tick_params(left=False, labelleft=False, length=0, pad=6)
    row_ax.set_xticks([-0.5, 0.5], minor=True)
    row_ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    row_ax.grid(which="minor", color="white", linewidth=1)
    row_ax.tick_params(which="minor", length=0)
    for y, label in enumerate(rows["inc_dec_counts"]):
        value = row_values[y]
        row_ax.text(
            0,
            y,
            label,
            ha="center",
            va="center",
            fontsize=8,
            color="white" if abs(value) >= 2 else "#111827",
        )
    for spine in row_ax.spines.values():
        spine.set_visible(False)

    bar = fig.colorbar(image, ax=[ax, row_ax], pad=0.02, fraction=0.035, shrink=0.58, extend="both")
    bar.set_label("Species influence on paired PD inc-dec trait contrast (fixed-reference z)")
    fig.text(
        0.03,
        1 - 0.11 / height,
        "PD species contributions across all six traits",
        fontsize=17,
        fontweight="bold",
        va="top",
    )
    subtitle = (
        "Cells show the study-weighted change in each paired increased-decreased trait contrast "
        "attributable to removing one reported species and recomputing its full non-MediaDive "
        f"profile denominator. Inc/dec requires a 95% interval entirely above/below zero and at "
        f"least {MIN_SUPPORT_STUDIES} supporting studies; otherwise cells are mixed or hatched "
        "low support. Gray = no trait annotation at or below the reported species label. Right "
        "column = raw species occurrence "
        "counts in all PD fecal profiles, shown separately from trait-specific influence."
    )
    fig.text(
        0.03,
        1 - 0.48 / height,
        "\n".join(wrap(subtitle, width=140)),
        fontsize=9,
        color="#4b5563",
        va="top",
    )
    fig.legend(
        handles=[
            Patch(facecolor=INCREASED_COLOR, label="positive influence"),
            Patch(facecolor=DECREASED_COLOR, label="negative influence"),
            Patch(facecolor=MISSING_COLOR, label="not annotated"),
            Patch(facecolor="white", edgecolor="#64748b", hatch="///", label="low support"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.60, 1 - 0.94 / height),
        ncol=4,
        frameon=False,
        fontsize=8,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({"svg.hashsalt": STEM}):
        fig.savefig(output_path, dpi=240, bbox_inches="tight", facecolor="white")
        fig.savefig(
            output_path.with_suffix(".svg"),
            bbox_inches="tight",
            facecolor="white",
            metadata={"Date": None},
        )
    plt.close(fig)


def sha256(path: Path) -> str:
    """Return a streamed SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    """Generate the six-trait influence heatmap and its auditable tables."""
    counts, trait_columns = read_profile_counts(PROFILE_MATRIX_PATH)
    metadata = load_profile_metadata()
    taxon_vectors = read_taxon_vectors(TAXA_TRAIT_MATRIX_PATH, trait_columns)
    evidence_path = OUTPUT_DIR / "pd_six_trait_taxon_profile_evidence.tsv"
    evidence = read_tsv(evidence_path)
    trait_summary = read_tsv(OUTPUT_DIR / "pd_six_trait_direction_summary.tsv")
    cells, rows, traits = build_influence_cells(
        counts,
        metadata,
        taxon_vectors,
        evidence,
        trait_summary,
    )
    plot_heatmap(cells, rows, traits, OUTPUT_DIR / f"{STEM}.png")
    write_tsv(cells, OUTPUT_DIR / f"{STEM}_cells.tsv")
    write_tsv(rows, OUTPUT_DIR / f"{STEM}_taxa.tsv")
    write_tsv(traits, OUTPUT_DIR / f"{STEM}_traits.tsv")
    manifest = {
        "method": "fixed-reference leave-one-reported-species-out influence",
        "traits": len(traits),
        "species": len(rows),
        "cells": len(cells),
        "non_mediadive_columns": len(trait_columns),
        "minimum_support_studies": MIN_SUPPORT_STUDIES,
        "inputs": {
            str(PROFILE_MATRIX_PATH.relative_to(PROFILE_MATRIX_PATH.parents[2])): sha256(
                PROFILE_MATRIX_PATH
            ),
            str(TAXA_TRAIT_MATRIX_PATH.relative_to(TAXA_TRAIT_MATRIX_PATH.parents[2])): sha256(
                TAXA_TRAIT_MATRIX_PATH
            ),
            str(evidence_path.relative_to(evidence_path.parents[3])): sha256(evidence_path),
        },
    }
    (OUTPUT_DIR / f"{STEM}_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    note = (
        "# PD Species-by-Six-Trait Influence Heatmap\n\n"
        f"![Heatmap]({STEM}.png)\n\n"
        f"Includes {len(rows)} reported species - the lowest reliable rank in the PD carrier "
        "data - and all six lollipop traits. Each trait cell is now specific to the species-trait "
        "combination: the analysis removes the species' complete non-MediaDive annotation vector "
        "from profiles where it was reported, recomputes the six log-fraction scores using the "
        "unchanged human-corpus reference, forms paired increased-minus-decreased experiment "
        "contrasts, averages experiments within study, and weights studies equally.\n\n"
        f"`inc` and `dec` require at least {MIN_SUPPORT_STUDIES} supporting studies and a 95% "
        "Student t interval entirely above or below zero. `mixed` retains supported estimates "
        "whose interval crosses zero; hatched `low support` cells do not receive a forced "
        "direction. Gray means the trait matrix has no annotation at or below that reported "
        "species label. The separate right-hand annotation retains the raw increased/decreased "
        "profile occurrence counts.\n\n"
        "This is a reported-species influence diagnostic, not a causal attribution or strain-level "
        "measurement. Removing a reported species changes both trait numerators and the nonlinear "
        "full denominator. Parent/child propagation and unresolved source assertions remain the "
        "limitations documented in #190; contrast orientation and trait semantics remain bounded "
        "by #188 and #189.\n\n"
        f"[Cell statistics]({STEM}_cells.tsv) | [Taxon row annotations]({STEM}_taxa.tsv) | "
        f"[Trait summary]({STEM}_traits.tsv) | [Manifest]({STEM}_manifest.json) | "
        f"[Vector figure]({STEM}.svg)\n"
    )
    (OUTPUT_DIR / f"{STEM}.md").write_text(note, encoding="utf-8")
    print(f"Wrote {len(rows)} species x {len(traits)} traits ({len(cells)} cells)")
    print(cells["status"].value_counts().to_string())
    print(OUTPUT_DIR / f"{STEM}.png")


if __name__ == "__main__":
    main()
