"""Trace the PD lollipop traits to observed taxa and profile-level evidence."""

from __future__ import annotations

import pandas as pd

from create_pd_cellulose_analysis import (
    FECES,
    OUTPUT_DIR,
    PD_CONDITION,
    TAXA_PATH,
    TAXA_TRAIT_MATRIX_PATH,
    build_pd_trait_direction_summary,
    build_zscore_matrix,
    load_profile_metadata,
    load_trait_summary,
    read_tsv,
    write_tsv,
)


def main() -> None:
    """Write reproducible carrier lists, shared carriers, and source occurrences."""
    traits = load_trait_summary()
    trait_ids = traits["trait_id"].tolist()
    metadata = load_profile_metadata()
    zscores = build_zscore_matrix(metadata, trait_ids)
    summary = build_pd_trait_direction_summary(metadata, zscores, traits)
    saved = read_tsv(OUTPUT_DIR / "pd_six_trait_direction_summary.tsv")
    pd.testing.assert_frame_equal(
        summary.sort_values("trait_id").reset_index(drop=True),
        saved.sort_values("trait_id").reset_index(drop=True),
        check_dtype=False,
    )
    selected = summary[summary["increased_minus_decreased_mean_z"] > 0]
    selected_ids = selected["trait_id"].tolist()
    if len(selected_ids) != 4:
        raise ValueError("The source lollipop no longer has four positive-side traits")

    profiles = metadata[
        metadata["condition"].eq(PD_CONDITION) & metadata["body_site_group"].eq(FECES)
    ].copy()
    profiles["direction"] = profiles["inc_dec"].map({"inc": "increased", "dec": "decreased"})
    taxa = read_tsv(TAXA_PATH)
    occurrences = taxa.merge(
        profiles[["study_id", "experiment_id", "direction", "profile_id"]],
        on=["study_id", "experiment_id", "direction"],
        how="inner",
        validate="many_to_one",
    ).drop_duplicates(["profile_id", "taxon_id"])
    matrix = read_tsv(
        TAXA_TRAIT_MATRIX_PATH,
        usecols=["taxon_id", "taxon_name", "taxonomic_rank", *trait_ids],
    )
    observed = occurrences.merge(
        matrix,
        on="taxon_id",
        how="left",
        suffixes=("_reported", ""),
        validate="many_to_one",
        indicator=True,
    )
    missing_ids = observed.loc[observed["_merge"].ne("both"), "taxon_id"].unique()
    if len(missing_ids):
        raise ValueError(f"Observed taxa missing from trait matrix: {missing_ids.tolist()}")
    observed["taxon_name"] = observed["taxon_name"].fillna(observed["taxon_name_reported"])
    observed["taxonomic_rank"] = observed["taxonomic_rank"].fillna(
        observed["taxonomic_rank_reported"]
    )
    evidence_frames = []
    records = []
    for trait in traits.itertuples(index=False):
        carriers = observed[pd.to_numeric(observed[trait.trait_id]) > 0].copy()
        carriers["trait_id"] = trait.trait_id
        carriers["trait"] = trait.trait_label
        carriers["annotation_count"] = carriers[trait.trait_id]
        carriers["profile_trait_z"] = carriers["profile_id"].map(zscores[trait.trait_id])
        carriers["profile_trait_positive_z"] = carriers["profile_trait_z"] > 0
        evidence_frames.append(
            carriers[
                [
                    "trait_id",
                    "trait",
                    "taxon_id",
                    "taxon_name",
                    "taxonomic_rank",
                    "study_id",
                    "experiment_id",
                    "profile_id",
                    "direction",
                    "annotation_count",
                    "profile_trait_z",
                    "profile_trait_positive_z",
                ]
            ]
        )
        for (taxon_id, name, rank), group in carriers.groupby(
            ["taxon_id", "taxon_name", "taxonomic_rank"],
            dropna=False,
        ):
            record = {
                "trait_id": trait.trait_id,
                "trait": trait.trait_label,
                "positive_side_in_lollipop": trait.trait_id in selected_ids,
                "taxon_id": taxon_id,
                "taxon_name": name,
                "taxonomic_rank": rank,
                "annotation_count": int(group["annotation_count"].iloc[0]),
            }
            for direction in ["increased", "decreased"]:
                subset = group[group["direction"].eq(direction)]
                enriched = subset[subset["profile_trait_positive_z"]]
                record[f"{direction}_profiles"] = subset["profile_id"].nunique()
                record[f"{direction}_studies"] = subset["study_id"].nunique()
                record[f"{direction}_positive_z_profiles"] = enriched["profile_id"].nunique()
                record[f"{direction}_positive_z_studies"] = enriched["study_id"].nunique()
            records.append(record)
    drivers = pd.DataFrame(records).sort_values(
        ["trait", "increased_positive_z_profiles", "increased_profiles", "taxon_name"],
        ascending=[True, False, False, True],
    )
    evidence = pd.concat(evidence_frames, ignore_index=True)
    if evidence.duplicated(["trait_id", "profile_id", "taxon_id"]).any():
        raise ValueError("Duplicated taxon/profile evidence")

    all_four = observed[observed[selected_ids].gt(0).all(axis=1)]
    overlap_records = []
    for (taxon_id, name, rank), group in all_four.groupby(
        ["taxon_id", "taxon_name", "taxonomic_rank"],
        dropna=False,
    ):
        record = {"taxon_id": taxon_id, "taxon_name": name, "taxonomic_rank": rank}
        for direction in ["increased", "decreased"]:
            subset = group[group["direction"].eq(direction)]
            record[f"{direction}_profiles"] = subset["profile_id"].nunique()
            record[f"{direction}_studies"] = subset["study_id"].nunique()
            simultaneous = zscores.loc[subset["profile_id"], selected_ids].gt(0).all(axis=1)
            record[f"{direction}_all_four_positive_z_profiles"] = int(simultaneous.sum())
        overlap_records.append(record)
    overlap = pd.DataFrame(overlap_records).sort_values("increased_profiles", ascending=False)
    write_tsv(drivers, OUTPUT_DIR / "pd_six_trait_taxon_carriers.tsv")
    write_tsv(evidence, OUTPUT_DIR / "pd_six_trait_taxon_profile_evidence.tsv")
    write_tsv(overlap, OUTPUT_DIR / "pd_four_positive_trait_shared_taxa.tsv")

    lines = [
        "# Taxa Associated With the PD Lollipop Traits",
        "",
        "Source: `pd_six_trait_direction_lollipop.png` and its companion summary. "
        "The recomputed summary matches the saved source table.",
        "",
        f"Scope: {len(profiles[profiles['inc_dec'].eq('inc')])} increased and "
        f"{len(profiles[profiles['inc_dec'].eq('dec')])} decreased PD fecal profiles. "
        "Occurrences are deduplicated by taxon and profile; names and traits come "
        "from the local KG-Microbe taxa-trait matrix.",
        "",
        "The four positive-side traits are moderately halophilic, facultatively aerobic, "
        "arsenate reduction, and coccus shaped. Coccus shaped is nearly neutral (+0.017 mean-z "
        "difference). A positive lollipop value is not evidence of statistical significance "
        "or gene expression. Cellulose breakdown is the strongest negative-side trait.",
        "",
        "These are observed annotation carriers, not estimated causal contributions to "
        "the z-score difference. Profile z>0 describes relative trait enrichment. "
        "The matrix counts trait-bearing taxa at or below each reported taxon "
        "(src/bsdbng/taxa_trait_matrix.py, compute_trait_counts). Different descendant "
        "taxa can supply different annotations, including for species with strain-level "
        "descendants. An all-four matrix intersection therefore does not establish "
        "co-occurrence of four traits in one organism. Counts across taxonomic ranks "
        "must not be added. Raw KG edges are not present locally, so original "
        "assertion and strain provenance were not independently verified.",
        "",
    ]
    for trait in summary.itertuples(index=False):
        rows = drivers[drivers["trait_id"].eq(trait.trait_id)]
        inc = rows[rows["increased_profiles"] > 0]
        lines.extend(
            [
                f"## {trait.trait}",
                "",
                f"Mean-z difference (increased - decreased): "
                f"{trait.increased_minus_decreased_mean_z:+.3f}. "
                f"{len(inc)} annotated taxon labels observed in increased profiles; "
                f"{len(rows)} across both directions.",
                "",
                "Top genus/species labels by occurrence in positive-z increased profiles:",
                "",
                "| Taxon | Rank | Increased profiles | Decreased profiles | "
                "Increased positive-z profiles | Increased positive-z studies |",
                "| --- | --- | ---: | ---: | ---: | ---: |",
            ]
        )
        for row in inc[inc["taxonomic_rank"].isin(["genus", "species"])].head(12).itertuples():
            lines.append(
                f"| {row.taxon_name} | {row.taxonomic_rank} | {row.increased_profiles} | "
                f"{row.decreased_profiles} | {row.increased_positive_z_profiles} | "
                f"{row.increased_positive_z_studies} |"
            )
        lines.append("")
    lines.extend(
        [
            "## Taxa Mapped to All Four Positive-Side Annotations",
            "",
            "Annotation intersection is separate from all four traits having z>0 in a profile.",
            "",
            "| Taxon | Rank | Increased profiles | Decreased profiles | "
            "Increased profiles with all four z>0 |",
            "| --- | --- | ---: | ---: | ---: |",
        ]
    )
    for row in overlap.itertuples():
        lines.append(
            f"| {row.taxon_name} | {row.taxonomic_rank} | {row.increased_profiles} | "
            f"{row.decreased_profiles} | {row.increased_all_four_positive_z_profiles} |"
        )
    lines.extend(
        [
            "",
            "## Complete Evidence",
            "",
            "- [All six trait carrier lists](pd_six_trait_taxon_carriers.tsv)",
            "- [Taxon/profile/study evidence](pd_six_trait_taxon_profile_evidence.tsv)",
            "- [All-four annotation intersection](pd_four_positive_trait_shared_taxa.tsv)",
            "",
        ]
    )
    (OUTPUT_DIR / "pd_trait_taxa_lookup.md").write_text("\n".join(lines), encoding="utf-8")
    print(drivers.groupby("trait").size().to_string())
    print(overlap.to_string(index=False))
    print(f"Wrote {len(evidence)} taxon/profile/trait evidence records to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
