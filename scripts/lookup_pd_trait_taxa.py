"""Trace the six PD lollipop traits to observed taxa and profile-level evidence.

Every trait is traced whatever its sign. Its observed direction, the sign of
increased_minus_decreased_mean_z in pd_six_trait_direction_summary.tsv, is carried as
data in every table (#584).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from create_pd_cellulose_analysis import (
    FECES,
    OUTPUT_DIR,
    PD_CONDITION,
    PROFILE_MATRIX_PATH,
    PROFILE_METADATA_PATH,
    TAXA_PATH,
    TAXA_TRAIT_MATRIX_PATH,
    build_pd_trait_direction_summary,
    build_zscore_matrix,
    load_profile_metadata,
    load_trait_summary,
    read_tsv,
    write_tsv,
)

SUMMARY_NAME = "pd_six_trait_direction_summary.tsv"
CARRIERS_NAME = "pd_six_trait_taxon_carriers.tsv"
EVIDENCE_NAME = "pd_six_trait_taxon_profile_evidence.tsv"
PROVENANCE_NAME = "pd_six_trait_taxon_evidence_provenance.json"
NOTE_NAME = "pd_trait_taxa_lookup.md"
DIRECTIONS = ("increased", "decreased")
DIRECTION_MEANINGS = {
    "increased": "its mean z is higher in increased than in decreased PD fecal profiles",
    "decreased": "its mean z is higher in decreased than in increased PD fecal profiles",
    "neither": "its mean z is the same in both directions",
}
TOP_TAXA_PER_TRAIT = 12
SUMMARY_REBUILD = "`just pd-cellulose-analysis`"


def trait_direction(increased_minus_decreased_mean_z: float) -> str:
    """Name the profile direction in which the trait's mean z is higher."""
    if pd.isna(increased_minus_decreased_mean_z):
        raise ValueError("A lollipop trait has no increased-minus-decreased mean z")
    if increased_minus_decreased_mean_z > 0:
        return "increased"
    if increased_minus_decreased_mean_z < 0:
        return "decreased"
    return "neither"


def trait_directions(summary: pd.DataFrame) -> pd.DataFrame:
    """One row per lollipop trait with its mean-z difference and observed direction."""
    directions = summary[["trait_id", "trait", "increased_minus_decreased_mean_z"]].copy()
    directions["trait_id"] = directions["trait_id"].astype(str)
    directions["trait_direction"] = directions["increased_minus_decreased_mean_z"].map(
        trait_direction
    )
    return directions.sort_values(
        "increased_minus_decreased_mean_z", ascending=False, kind="stable"
    ).reset_index(drop=True)


def pd_feces_mask(metadata: pd.DataFrame) -> pd.Series:
    """PD fecal profiles, the scope of the evidence."""
    return metadata["condition"].eq(PD_CONDITION) & metadata["body_site_group"].eq(FECES)


def pd_feces_profile_ids(metadata: pd.DataFrame) -> list[str]:
    """Sorted ids of the PD fecal profiles that the evidence is built from."""
    return sorted(metadata.loc[pd_feces_mask(metadata), "profile_id"].astype(str))


def sha256_chunks(path: Path) -> list[str]:
    """Return a streamed SHA-256 digest in unambiguous eight-character chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    value = digest.hexdigest()
    return [value[offset : offset + 8] for offset in range(0, len(value), 8)]


def tracked_inputs(output_dir: Path) -> dict[str, Path]:
    """The files the evidence is built from, by name: inputs, profile metadata and lollipop."""
    paths = [
        TAXA_PATH,
        TAXA_TRAIT_MATRIX_PATH,
        PROFILE_MATRIX_PATH,
        PROFILE_METADATA_PATH,
        output_dir / SUMMARY_NAME,
    ]
    return {path.name: path for path in paths}


def input_fingerprints(inputs: dict[str, Path]) -> dict[str, list[str]]:
    """Each file's sha256 in eight-character chunks; a missing file has no chunks."""
    return {name: sha256_chunks(path) if path.exists() else [] for name, path in inputs.items()}


def table_fingerprints(output_dir: Path) -> dict[str, list[str]]:
    """The carrier and evidence tables' hashes, which bind them to their provenance (#703)."""
    return input_fingerprints({name: output_dir / name for name in (CARRIERS_NAME, EVIDENCE_NAME)})


def changed_files(recorded: dict[str, list[str]], current: dict[str, list[str]]) -> list[str]:
    """Names whose hashes differ, counting a file that only one side has."""
    return sorted(
        name for name in set(recorded) | set(current) if recorded.get(name) != current.get(name)
    )


def check_inputs_unchanged(read: dict[str, list[str]], current: dict[str, list[str]]) -> None:
    """Refuse to write evidence from inputs that changed while the lookup read them (#703)."""
    changed = changed_files(read, current)
    if changed:
        msg = (
            f"Inputs changed while `just pd-trait-taxa` read them: {', '.join(changed)}. "
            "Rerun it once they stop changing."
        )
        raise ValueError(msg)


def check_saved_summary(recomputed: pd.DataFrame, saved: pd.DataFrame) -> None:
    """Refuse a saved lollipop summary that the current profiles no longer reproduce (#693)."""
    try:
        pd.testing.assert_frame_equal(
            recomputed.sort_values("trait_id").reset_index(drop=True),
            saved.sort_values("trait_id").reset_index(drop=True),
            check_dtype=False,
        )
    except AssertionError as error:
        detail = str(error).strip().splitlines()[0]
        msg = (
            f"{SUMMARY_NAME} no longer matches the summary recomputed from the current PD fecal "
            f"profiles ({detail}). Rebuild it with {SUMMARY_REBUILD}, then rerun "
            "`just pd-trait-taxa`."
        )
        raise ValueError(msg) from error


def evidence_provenance(
    profile_ids: list[str],
    directions: pd.DataFrame,
    input_hashes: dict[str, list[str]],
    table_hashes: dict[str, list[str]],
) -> dict[str, object]:
    """What the readers check: profile set, lollipop traits, inputs read and tables written."""
    return {
        "evidence": EVIDENCE_NAME,
        "carriers": CARRIERS_NAME,
        "trait_summary": SUMMARY_NAME,
        **profile_set_fingerprint(profile_ids),
        "pd_feces_profile_ids_encoding": "sha256 of the sorted profile_id values, one per line",
        "traits": directions.to_dict(orient="records"),
        "input_sha256_chunks": input_hashes,
        "output_sha256_chunks": table_hashes,
        "sha256_chunk_encoding": "concatenate each file's eight chunks in order",
    }


def profile_set_fingerprint(profile_ids: list[str]) -> dict[str, object]:
    """Count and sha256 of a profile set, so that a later run can tell if it changed."""
    ids = sorted(profile_ids)
    digest = hashlib.sha256("\n".join(ids).encode("utf-8")).hexdigest()
    return {"pd_feces_profiles": len(ids), "pd_feces_profile_ids_sha256": digest}


def observed_side(direction: str) -> str | None:
    """Profile direction whose counts rank carriers of a trait, or None for no difference."""
    return direction if direction in DIRECTIONS else None


def build_carriers_and_evidence(
    observed: pd.DataFrame,
    zscores: pd.DataFrame,
    directions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """List each trait's annotated taxa and their taxon/profile evidence."""
    evidence_frames = []
    records = []
    for trait in directions.itertuples(index=False):
        carriers = observed[pd.to_numeric(observed[trait.trait_id]) > 0].copy()
        carriers["trait_id"] = trait.trait_id
        carriers["trait"] = trait.trait
        carriers["trait_direction"] = trait.trait_direction
        carriers["annotation_count"] = carriers[trait.trait_id]
        carriers["profile_trait_z"] = carriers["profile_id"].map(zscores[trait.trait_id])
        carriers["profile_trait_positive_z"] = carriers["profile_trait_z"] > 0
        evidence_frames.append(
            carriers[
                [
                    "trait_id",
                    "trait",
                    "trait_direction",
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
                "trait": trait.trait,
                "trait_direction": trait.trait_direction,
                "trait_increased_minus_decreased_mean_z": trait.increased_minus_decreased_mean_z,
                "taxon_id": taxon_id,
                "taxon_name": name,
                "taxonomic_rank": rank,
                "annotation_count": int(group["annotation_count"].iloc[0]),
            }
            for direction in DIRECTIONS:
                subset = group[group["direction"].eq(direction)]
                enriched = subset[subset["profile_trait_positive_z"]]
                record[f"{direction}_profiles"] = subset["profile_id"].nunique()
                record[f"{direction}_studies"] = subset["study_id"].nunique()
                record[f"{direction}_positive_z_profiles"] = enriched["profile_id"].nunique()
                record[f"{direction}_positive_z_studies"] = enriched["study_id"].nunique()
            records.append(record)
    carriers_table = sort_carriers(pd.DataFrame(records), directions)
    evidence = pd.concat(evidence_frames, ignore_index=True)
    if evidence.duplicated(["trait_id", "profile_id", "taxon_id"]).any():
        raise ValueError("Duplicated taxon/profile evidence")
    return carriers_table, evidence


def sort_carriers(carriers: pd.DataFrame, directions: pd.DataFrame) -> pd.DataFrame:
    """Order traits as the lollipop does and rank each trait's carriers on its observed side."""
    trait_order = {trait_id: order for order, trait_id in enumerate(directions["trait_id"])}
    ranked = carriers.copy()
    ranked["_trait_order"] = ranked["trait_id"].map(trait_order)
    ranked["_side_positive_z_profiles"] = 0
    ranked["_side_profiles"] = 0
    for direction in DIRECTIONS:
        on_side = ranked["trait_direction"].eq(direction)
        ranked.loc[on_side, "_side_positive_z_profiles"] = ranked.loc[
            on_side, f"{direction}_positive_z_profiles"
        ]
        ranked.loc[on_side, "_side_profiles"] = ranked.loc[on_side, f"{direction}_profiles"]
    neither = ~ranked["trait_direction"].isin(DIRECTIONS)
    ranked.loc[neither, "_side_positive_z_profiles"] = (
        ranked.loc[neither, "increased_positive_z_profiles"]
        + ranked.loc[neither, "decreased_positive_z_profiles"]
    )
    ranked.loc[neither, "_side_profiles"] = (
        ranked.loc[neither, "increased_profiles"] + ranked.loc[neither, "decreased_profiles"]
    )
    ranked = ranked.sort_values(
        ["_trait_order", "_side_positive_z_profiles", "_side_profiles", "taxon_name"],
        ascending=[True, False, False, True],
        kind="stable",
    )
    return ranked.drop(
        columns=["_trait_order", "_side_positive_z_profiles", "_side_profiles"]
    ).reset_index(drop=True)


def direction_phrase(direction: str, difference: float) -> str:
    """Describe a trait's observed direction with its mean-z difference."""
    if direction == "neither":
        return f"no difference in mean z between directions ({difference:+.3f})"
    return f"{direction} ({difference:+.3f})"


def lookup_note(
    directions: pd.DataFrame,
    carriers: pd.DataFrame,
    profiles: pd.DataFrame,
) -> str:
    """Markdown companion listing, per trait, its direction and top carrier labels."""
    phrases = {
        trait.trait_id: direction_phrase(
            str(trait.trait_direction), float(trait.increased_minus_decreased_mean_z)
        )
        for trait in directions.itertuples(index=False)
    }
    trait_list = "; ".join(
        f"{trait.trait} {phrases[trait.trait_id]}" for trait in directions.itertuples(index=False)
    )
    lines = [
        "# Taxa Associated With the PD Lollipop Traits",
        "",
        "Source: `pd_six_trait_direction_lollipop.png` and its companion summary. "
        "The recomputed summary matches the saved source table.",
        "",
        f"Scope: {int(profiles['inc_dec'].eq('inc').sum())} increased and "
        f"{int(profiles['inc_dec'].eq('dec').sum())} decreased PD fecal profiles. "
        "Occurrences are deduplicated by taxon and profile; names and traits come "
        "from the local KG-Microbe taxa-trait matrix.",
        "",
        "All six lollipop traits are traced whatever their sign. Each trait's observed "
        "direction is the profile direction in which its mean z is higher, from the sign of "
        "`increased_minus_decreased_mean_z` in `pd_six_trait_direction_summary.tsv`; the "
        "tables carry it in `trait_direction`. Observed direction (increased minus decreased "
        f"mean z): {trait_list}. A lollipop value is a descriptive difference in mean z, not "
        "evidence of statistical significance or of gene expression.",
        "",
        "These are observed annotation carriers, not estimated causal contributions to "
        "the z-score difference. Profile z>0 describes relative trait enrichment. "
        "The matrix counts trait-bearing taxa at or below each reported taxon "
        "(src/bsdbng/taxa_trait_matrix.py, compute_trait_counts). Different descendant "
        "taxa can supply different annotations, including for species with strain-level "
        "descendants, so one reported label carrying several trait annotations does not "
        "establish that one organism has those traits together. Counts across taxonomic "
        "ranks must not be added. Raw KG edges are not present locally, so original "
        "assertion and strain provenance were not independently verified.",
        "",
    ]
    for trait in directions.itertuples(index=False):
        rows = carriers[carriers["trait_id"].eq(trait.trait_id)]
        side = observed_side(str(trait.trait_direction))
        counts = ", ".join(
            f"{int(rows[f'{direction}_profiles'].gt(0).sum())} in {direction} profiles"
            for direction in DIRECTIONS
        )
        if side is None:
            shown = rows[(rows["increased_profiles"] + rows["decreased_profiles"]) > 0]
            ranking = "Top genus/species labels by occurrence in positive-z profiles:"
        else:
            shown = rows[rows[f"{side}_profiles"] > 0]
            ranking = (
                f"Top genus/species labels by occurrence in positive-z {side} profiles, "
                "the trait's observed direction:"
            )
        lines.extend(
            [
                f"## {trait.trait}",
                "",
                f"Observed direction: {phrases[trait.trait_id]}; "
                f"{DIRECTION_MEANINGS[str(trait.trait_direction)]}. Annotated taxon labels "
                f"observed: {counts}; {len(rows)} across both directions.",
                "",
                ranking,
                "",
                "| Taxon | Rank | Increased profiles | Decreased profiles | "
                "Increased positive-z profiles | Decreased positive-z profiles |",
                "| --- | --- | ---: | ---: | ---: | ---: |",
            ]
        )
        top = shown[shown["taxonomic_rank"].isin(["genus", "species"])].head(TOP_TAXA_PER_TRAIT)
        for row in top.itertuples(index=False):
            lines.append(
                f"| {row.taxon_name} | {row.taxonomic_rank} | {row.increased_profiles} | "
                f"{row.decreased_profiles} | {row.increased_positive_z_profiles} | "
                f"{row.decreased_positive_z_profiles} |"
            )
        lines.append("")
    lines.extend(
        [
            "## Complete Evidence",
            "",
            f"- [All six trait carrier lists]({CARRIERS_NAME})",
            f"- [Taxon/profile/study evidence]({EVIDENCE_NAME})",
            f"- [PD fecal profile set and trait directions used]({PROVENANCE_NAME})",
            "",
        ]
    )
    return "\n".join(lines)


def trace_traits() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read the inputs and trace the lollipop traits.

    Returns the trait directions, the PD fecal profiles, the carriers and the evidence.
    """
    traits = load_trait_summary()
    trait_ids = traits["trait_id"].tolist()
    metadata = load_profile_metadata()
    zscores = build_zscore_matrix(metadata, trait_ids)
    summary = build_pd_trait_direction_summary(metadata, zscores, traits)
    saved = read_tsv(OUTPUT_DIR / SUMMARY_NAME)
    check_saved_summary(summary, saved)
    directions = trait_directions(saved)

    profiles = metadata[pd_feces_mask(metadata)].copy()
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
    carriers, evidence = build_carriers_and_evidence(observed, zscores, directions)
    return directions, profiles, carriers, evidence


def main() -> None:
    """Write carrier lists, taxon/profile evidence, and the provenance the readers check."""
    inputs = tracked_inputs(OUTPUT_DIR)
    # Hashed before reading and checked after, so the recorded hashes are of what was read.
    read_hashes = input_fingerprints(inputs)
    directions, profiles, carriers, evidence = trace_traits()
    check_inputs_unchanged(read_hashes, input_fingerprints(inputs))
    write_tsv(carriers, OUTPUT_DIR / CARRIERS_NAME)
    write_tsv(evidence, OUTPUT_DIR / EVIDENCE_NAME)
    provenance = evidence_provenance(
        pd_feces_profile_ids(profiles), directions, read_hashes, table_fingerprints(OUTPUT_DIR)
    )
    (OUTPUT_DIR / PROVENANCE_NAME).write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    (OUTPUT_DIR / NOTE_NAME).write_text(
        lookup_note(directions, carriers, profiles), encoding="utf-8"
    )
    print(carriers.groupby(["trait", "trait_direction"]).size().to_string())
    print(f"Wrote {len(evidence)} taxon/profile/trait evidence records to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
