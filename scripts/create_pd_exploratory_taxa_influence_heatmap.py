"""Plot an exploratory PD influence heatmap across species, genus, and family taxa."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from create_pd_cellulose_analysis import (
    OUTPUT_DIR,
    PROFILE_MATRIX_PATH,
    TAXA_PATH,
    TAXA_TRAIT_MATRIX_PATH,
    load_profile_metadata,
    read_tsv,
    write_tsv,
)
from create_pd_six_trait_taxon_heatmap import (
    build_influence_cells,
    build_pd_taxa_reference_table,
    literature_reference_table,
    plot_heatmap,
    read_profile_counts,
    read_taxon_vectors,
    sha256_chunks,
    write_reference_interpretation,
)

STEM = "pd_exploratory_taxa_six_traits_influence_heatmap"
REFERENCE_TABLE_STEM = "pd_exploratory_taxa_six_traits_reference"
SELECTED_RANKS = ("species", "genus", "family")
EXPLORATORY_MIN_SUPPORT_STUDIES = 2
MAX_ROWS = 30
RANK_PRIORITY = {"family": 3, "genus": 2, "species": 1}


def rank_priority(rank: object) -> int:
    """Prefer roll-up rows when equally supported."""
    return RANK_PRIORITY.get(str(rank), 0)


def select_exploratory_rows(
    rows: pd.DataFrame,
    reference_table: pd.DataFrame,
    max_rows: int = MAX_ROWS,
) -> pd.DataFrame:
    """Keep literature-seeded taxa first, then top-supported exploratory roll-ups."""
    annotations = reference_table[
        [
            "taxon_id",
            "literature_taxon",
            "literature_direction",
            "literature_match_term",
            "direction_agreement",
        ]
    ].copy()
    annotated = rows.merge(annotations, on="taxon_id", how="left")
    annotated["literature_direction"] = annotated["literature_direction"].fillna("not listed")
    annotated["literature_matched"] = annotated["literature_direction"].ne("not listed")
    annotated["occurrence_profiles"] = annotated["increased_profiles"].astype(int) + annotated[
        "decreased_profiles"
    ].astype(int)
    annotated["rank_priority"] = annotated["taxonomic_rank"].map(rank_priority)
    annotated["selection_reason"] = "top supported exploratory roll-up"
    annotated.loc[annotated["literature_matched"], "selection_reason"] = "PD literature seed"

    sort_columns = [
        "literature_matched",
        "paired_support_studies",
        "occurrence_profiles",
        "rank_priority",
        "annotated_trait_count",
        "max_absolute_influence_z",
        "taxon_name",
    ]
    sorted_rows = annotated.sort_values(
        sort_columns,
        ascending=[False, False, False, False, False, False, True],
    )
    literature_rows = sorted_rows[sorted_rows["literature_matched"]]
    non_literature_rows = sorted_rows[~sorted_rows["literature_matched"]]

    selected = pd.concat(
        [
            literature_rows.head(max_rows),
            non_literature_rows.head(max(0, max_rows - min(len(literature_rows), max_rows))),
        ],
        ignore_index=True,
    ).head(max_rows)
    return selected.drop(columns=["literature_matched", "rank_priority"])


def reference_table_for_rows(
    reference_table: pd.DataFrame, selected_rows: pd.DataFrame
) -> pd.DataFrame:
    """Return the joined reference table in the displayed row order."""
    selected = selected_rows[["taxon_id", "selection_reason"]].copy()
    return selected.merge(reference_table, on="taxon_id", how="left")


def write_exploratory_note(
    rows: pd.DataFrame,
    cells: pd.DataFrame,
    traits: pd.DataFrame,
    reference_table: pd.DataFrame,
    source_rows: int,
    output_path: Path,
) -> None:
    """Write a short Markdown companion for the exploratory heatmap."""
    literature_rows = int(reference_table["literature_direction"].ne("not listed").sum())
    status_counts = cells["status"].value_counts().to_dict()
    note = (
        "# Exploratory PD Taxa-by-Six-Trait Influence Heatmap\n\n"
        f"![Heatmap]({STEM}.png)\n\n"
        "This is the expanded companion to the strict species-only influence heatmap. It includes "
        "reported PD fecal taxa at species, genus, and family rank, uses a lower exploratory "
        f"support threshold of {EXPLORATORY_MIN_SUPPORT_STUDIES} studies, and caps the display at "
        f"{MAX_ROWS} rows so the figure stays readable. Rows are selected by keeping DOI-backed "
        "PD literature matches first, then adding the best-supported non-literature rank "
        "roll-ups.\n\n"
        f"The displayed figure contains {len(rows)} taxa out of {source_rows} eligible "
        f"species/genus/family rows, {literature_rows} of which match the curated PD taxa "
        f"reference. It has {len(traits)} traits and {len(cells)} taxon-trait cells. Cell status "
        f"counts: {status_counts}.\n\n"
        "Use this view to find candidate taxa and rank-level patterns that the conservative "
        "species-only heatmap may miss. Treat it as hypothesis-generating: genus and family rows "
        "can mix multiple species, and the relaxed support threshold is meant to reveal leads, "
        "not final directional claims.\n\n"
        f"[Cell statistics]({STEM}_cells.tsv) | [Taxon row annotations]({STEM}_taxa.tsv) | "
        f"[PD specificity]({STEM}_specificity.tsv) | [Trait summary]({STEM}_traits.tsv) | "
        f"[PD taxa reference]({REFERENCE_TABLE_STEM}_table.tsv) | "
        f"[Reference interpretation]({REFERENCE_TABLE_STEM}_interpretation.md) | "
        f"[Manifest]({STEM}_manifest.json) | [Vector figure]({STEM}.svg)\n"
    )
    output_path.write_text(note, encoding="utf-8")


def main() -> None:
    """Generate the exploratory genus/family/species influence heatmap and tables."""
    counts, trait_columns = read_profile_counts(PROFILE_MATRIX_PATH)
    metadata = load_profile_metadata()
    taxon_vectors = read_taxon_vectors(TAXA_TRAIT_MATRIX_PATH, trait_columns)
    evidence_path = OUTPUT_DIR / "pd_six_trait_taxon_profile_evidence.tsv"
    evidence = read_tsv(evidence_path)
    trait_summary = read_tsv(OUTPUT_DIR / "pd_six_trait_direction_summary.tsv")

    cells, rows, traits, removed_low_support_taxa = build_influence_cells(
        counts,
        metadata,
        taxon_vectors,
        evidence,
        trait_summary,
        taxonomic_ranks=SELECTED_RANKS,
        min_support_studies=EXPLORATORY_MIN_SUPPORT_STUDIES,
        remove_all_low_support_taxa=False,
    )
    source_rows = len(rows)
    literature_reference = literature_reference_table()
    full_reference_table = build_pd_taxa_reference_table(rows, cells, literature_reference)
    selected_rows = select_exploratory_rows(rows, full_reference_table)
    selected_ids = selected_rows["taxon_id"].astype(str).tolist()
    selected_cells = cells[cells["taxon_id"].astype(str).isin(selected_ids)].copy()
    selected_reference_table = reference_table_for_rows(full_reference_table, selected_rows)

    plot_heatmap(
        selected_cells,
        selected_rows,
        traits,
        OUTPUT_DIR / f"{STEM}.png",
        title="PD species + genus + family contributions across all six traits",
        taxon_unit="taxon",
        method_note=(
            "Exploratory rows include reported species, genus, and family taxa; literature-matched "
            "PD taxa are prioritized before top-supported non-literature roll-ups."
        ),
        top_margin_inches=2.85,
        legend_y_inches=2.42,
    )
    write_tsv(selected_cells, OUTPUT_DIR / f"{STEM}_cells.tsv")
    write_tsv(selected_rows, OUTPUT_DIR / f"{STEM}_taxa.tsv")
    specificity_columns = [
        "taxon_id",
        "taxon_name",
        "taxonomic_rank",
        "pd_feces_profiles",
        "other_feces_profiles",
        "pd_other_feces_counts",
        "pd_feces_studies",
        "other_feces_studies",
        "pd_specificity_fraction",
        "pd_vs_other_feces_log2_ratio",
    ]
    write_tsv(selected_rows[specificity_columns], OUTPUT_DIR / f"{STEM}_specificity.tsv")
    write_tsv(traits, OUTPUT_DIR / f"{STEM}_traits.tsv")
    write_tsv(selected_reference_table, OUTPUT_DIR / f"{REFERENCE_TABLE_STEM}_table.tsv")
    write_reference_interpretation(
        selected_reference_table,
        OUTPUT_DIR / f"{REFERENCE_TABLE_STEM}_interpretation.md",
    )

    manifest = {
        "method": "exploratory fixed-reference leave-one-reported-taxon-out influence",
        "taxonomic_ranks": list(SELECTED_RANKS),
        "minimum_support_studies": EXPLORATORY_MIN_SUPPORT_STUDIES,
        "max_displayed_rows": MAX_ROWS,
        "eligible_rows_before_display_cap": source_rows,
        "displayed_rows": len(selected_rows),
        "traits": len(traits),
        "cells": len(selected_cells),
        "filled_cells": int(selected_cells["mean_influence_contrast_z"].notna().sum()),
        "direct_annotation_cells": int(selected_cells["has_trait_annotation"].sum()),
        "denominator_only_cells": int((~selected_cells["has_trait_annotation"]).sum()),
        "rows_removed_all_low_support": removed_low_support_taxa,
        "literature_reference_taxa": len(literature_reference),
        "displayed_rows_with_literature_reference": int(
            selected_reference_table["literature_direction"].ne("not listed").sum()
        ),
        "literature_bug_sig_direction_agreement_counts": {
            str(label): int(count)
            for label, count in selected_reference_table["direction_agreement"]
            .value_counts()
            .items()
        },
        "status_counts": {
            str(label): int(count)
            for label, count in selected_cells["status"].value_counts().items()
        },
        "selection": (
            "DOI-backed PD literature matches first, then top-supported non-literature "
            "species/genus/family roll-ups"
        ),
        "input_sha256_chunks": {
            PROFILE_MATRIX_PATH.relative_to(PROFILE_MATRIX_PATH.parents[2]).as_posix(): (
                sha256_chunks(PROFILE_MATRIX_PATH)
            ),
            TAXA_TRAIT_MATRIX_PATH.relative_to(
                TAXA_TRAIT_MATRIX_PATH.parents[2]
            ).as_posix(): sha256_chunks(TAXA_TRAIT_MATRIX_PATH),
            TAXA_PATH.relative_to(TAXA_PATH.parents[2]).as_posix(): sha256_chunks(TAXA_PATH),
            evidence_path.relative_to(evidence_path.parents[3]).as_posix(): sha256_chunks(
                evidence_path
            ),
        },
        "sha256_chunk_encoding": "concatenate each path's eight chunks in order",
    }
    (OUTPUT_DIR / f"{STEM}_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    write_exploratory_note(
        selected_rows,
        selected_cells,
        traits,
        selected_reference_table,
        source_rows,
        OUTPUT_DIR / f"{STEM}.md",
    )
    print(
        f"Wrote {len(selected_rows)} exploratory taxa x {len(traits)} traits "
        f"from {source_rows} eligible rows"
    )
    print(selected_cells["status"].value_counts().to_string())
    print(OUTPUT_DIR / f"{STEM}.png")


if __name__ == "__main__":
    main()
