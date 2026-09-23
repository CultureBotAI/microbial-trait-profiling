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
    TAXA_DIRECTION_TO_PROFILE,
    TAXA_PATH,
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
UNAVAILABLE_COLOR = "#e5e7eb"
REFERENCE_TABLE_STEM = "pd_species_six_traits_reference"
LITERATURE_REFERENCE_NAME = "pd_taxa_literature_reference"
LITERATURE_REVIEW_DATE = "2026-09-22"
LITERATURE_COMPARISON_COLORS = {
    "agrees": "#dcfce7",
    "differs": "#fee2e2",
    "balanced in BugSigDB": "#fef3c7",
    "literature non-directional": "#ede9fe",
    "not in literature reference": "#f3f4f6",
}
LITERATURE_LABELS = {
    "increased": "lit inc",
    "decreased": "lit dec",
    "variable": "var",
    "associated": "assoc",
    "not listed": "not ref",
}
LITERATURE_COMPARISON_LABELS = {
    "agrees": "agrees",
    "differs": "DIFFERS",
    "balanced in BugSigDB": "balanced",
    "literature non-directional": "non-dir",
    "not in literature reference": "no ref",
}
PD_TAXA_LITERATURE_REFERENCE = [
    {
        "literature_taxon": "Akkermansia",
        "literature_rank": "genus",
        "literature_direction": "increased",
        "literature_evidence": "Often reported as enriched in PD fecal microbiome studies.",
        "doi_citations": "10.1186/s13073-017-0428-y;10.1002/mds.28119;10.1038/s41467-022-34667-x",
        "match_terms": "Akkermansia;Akkermansia muciniphila",
    },
    {
        "literature_taxon": "Bifidobacterium",
        "literature_rank": "genus",
        "literature_direction": "increased",
        "literature_evidence": (
            "Reported as increased at genus/family level, although one cohort showed "
            "attenuation after excluding selected medication users."
        ),
        "doi_citations": "10.1002/mds.26942;10.1016/j.parkreldis.2016.08.019",
        "match_terms": "Bifidobacterium;Bifidobacteriaceae",
    },
    {
        "literature_taxon": "Lactobacillus",
        "literature_rank": "genus",
        "literature_direction": "increased",
        "literature_evidence": "Reported as increased in PD and related motor phenotype studies.",
        "doi_citations": "10.1371/journal.pone.0142164;10.1002/mds.26942",
        "match_terms": "Lactobacillus;Lactobacillaceae",
    },
    {
        "literature_taxon": "Enterobacteriaceae",
        "literature_rank": "family",
        "literature_direction": "increased",
        "literature_evidence": (
            "Reported as enriched in PD and linked with postural instability/gait phenotype."
        ),
        "doi_citations": "10.1002/mds.26069;10.1016/j.parkreldis.2016.08.019",
        "match_terms": (
            "Enterobacteriaceae;Escherichia;Escherichia coli;Shigella;Escherichia-Shigella"
        ),
    },
    {
        "literature_taxon": "Escherichia-Shigella",
        "literature_rank": "genus group",
        "literature_direction": "increased",
        "literature_evidence": "Reported as increased in PD metagenomic profiles.",
        "doi_citations": "10.1007/s11427-016-9001-4",
        "match_terms": "Escherichia-Shigella;Escherichia;Escherichia coli;Shigella",
    },
    {
        "literature_taxon": "Streptococcus",
        "literature_rank": "genus",
        "literature_direction": "increased",
        "literature_evidence": "Often increased, with species-level heterogeneity across cohorts.",
        "doi_citations": "10.1007/s11427-016-9001-4;10.1038/s41467-022-34667-x",
        "match_terms": "Streptococcus",
    },
    {
        "literature_taxon": "Enterococcus",
        "literature_rank": "genus",
        "literature_direction": "increased",
        "literature_evidence": "Reported as increased in PD metagenomic profiles.",
        "doi_citations": "10.1007/s11427-016-9001-4",
        "match_terms": "Enterococcus",
    },
    {
        "literature_taxon": "Proteus",
        "literature_rank": "genus",
        "literature_direction": "increased",
        "literature_evidence": "Reported as increased in PD metagenomic profiles.",
        "doi_citations": "10.1007/s11427-016-9001-4",
        "match_terms": "Proteus",
    },
    {
        "literature_taxon": "Alistipes",
        "literature_rank": "genus/species",
        "literature_direction": "increased",
        "literature_evidence": "Reported as increased in PD, including species-level signals.",
        "doi_citations": "10.1093/brain/awaa201;10.1038/s41467-022-34667-x",
        "match_terms": "Alistipes;Alistipes finegoldii;Alistipes indistinctus",
    },
    {
        "literature_taxon": "Methanobrevibacter",
        "literature_rank": "genus",
        "literature_direction": "increased",
        "literature_evidence": "Reported as increased in PD metagenomic species-level analyses.",
        "doi_citations": "10.1038/s41467-022-34667-x",
        "match_terms": "Methanobrevibacter;Methanobrevibacter smithii",
    },
    {
        "literature_taxon": "Prevotella",
        "literature_rank": "genus/family",
        "literature_direction": "decreased",
        "literature_evidence": (
            "Prevotella/Prevotellaceae depletion is a repeated PD microbiome finding."
        ),
        "doi_citations": "10.1002/mds.26069;10.1186/s13073-017-0428-y",
        "match_terms": "Prevotella;Prevotella copri;Prevotellaceae",
    },
    {
        "literature_taxon": "Faecalibacterium",
        "literature_rank": "genus/species",
        "literature_direction": "decreased",
        "literature_evidence": "Reported as depleted in PD, including F. prausnitzii.",
        "doi_citations": (
            "10.1016/j.parkreldis.2016.08.019;10.1007/s11427-016-9001-4;10.1002/mds.28119"
        ),
        "match_terms": "Faecalibacterium;Faecalibacterium prausnitzii",
    },
    {
        "literature_taxon": "Roseburia",
        "literature_rank": "genus",
        "literature_direction": "decreased",
        "literature_evidence": (
            "Reported as depleted at genus level, with R. faecis also decreased in a "
            "large species-level PD metagenomic study."
        ),
        "doi_citations": ("10.1002/mds.26307;10.1002/mds.28119;10.1038/s41467-022-34667-x"),
        "match_terms": "Roseburia",
    },
    {
        "literature_taxon": "Blautia",
        "literature_rank": "genus",
        "literature_direction": "decreased",
        "literature_evidence": "Reported as depleted in PD gut microbiome cohorts.",
        "doi_citations": "10.1002/mds.26307;10.1007/s11427-016-9001-4",
        "match_terms": "Blautia",
    },
    {
        "literature_taxon": "Ruminococcus",
        "literature_rank": "genus",
        "literature_direction": "decreased",
        "literature_evidence": (
            "Reported as decreased at genus level; species assignments may vary by taxonomy."
        ),
        "doi_citations": "10.1007/s11427-016-9001-4",
        "match_terms": (
            "Ruminococcus;Ruminococcus bromii;Ruminococcus torques;Mediterraneibacter torques"
        ),
    },
    {
        "literature_taxon": "Lachnospiraceae",
        "literature_rank": "family",
        "literature_direction": "decreased",
        "literature_evidence": "Reported as depleted in PD gut microbiome cohorts.",
        "doi_citations": "10.1002/mds.26307;10.1002/mds.26942",
        "match_terms": "Lachnospiraceae",
    },
    {
        "literature_taxon": "Coprococcus",
        "literature_rank": "genus",
        "literature_direction": "decreased",
        "literature_evidence": "Reported as depleted in PD gut microbiome cohorts.",
        "doi_citations": "10.1002/mds.26307",
        "match_terms": "Coprococcus",
    },
    {
        "literature_taxon": "Dorea",
        "literature_rank": "genus",
        "literature_direction": "decreased",
        "literature_evidence": "Reported as depleted in PD gut microbiome cohorts.",
        "doi_citations": "10.1002/mds.26307",
        "match_terms": "Dorea",
    },
    {
        "literature_taxon": "Bacteroides fragilis group",
        "literature_rank": "species group",
        "literature_direction": "decreased",
        "literature_evidence": "Reported as decreased in a PD gut microbiome cohort.",
        "doi_citations": "10.1371/journal.pone.0142164",
        "match_terms": "Bacteroides fragilis;Bacteroides fragilis group",
    },
    {
        "literature_taxon": "Helicobacter pylori",
        "literature_rank": "species",
        "literature_direction": "associated",
        "literature_evidence": (
            "Associated with PD risk, severity, or levodopa response; not expected to behave "
            "like a fecal trait-row species."
        ),
        "doi_citations": "10.1111/hel.12398;10.1016/j.clineuro.2018.09.039",
        "match_terms": "Helicobacter pylori",
    },
]


def format_effect_value(value: float) -> str:
    """Format signed heatmap effects without distracting negative zeros."""
    if not math.isfinite(value):
        return "NA"
    if abs(value) < 0.005:
        return "0.00"
    return f"{value:+.2f}"


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


def taxon_specificity_from_occurrences(
    taxa: pd.DataFrame,
    metadata: pd.DataFrame,
    taxon_ids: pd.Index,
) -> pd.DataFrame:
    """Count each selected reported taxon label in PD feces versus other feces."""
    required_taxa = {"study_id", "experiment_id", "direction", "taxon_id"}
    required_metadata = {
        "profile_id",
        "study_id",
        "experiment_id",
        "inc_dec",
        "condition",
        "body_site_group",
    }
    missing_taxa = sorted(required_taxa.difference(taxa.columns))
    missing_metadata = sorted(required_metadata.difference(metadata.columns))
    if missing_taxa:
        raise ValueError(f"Missing taxa columns for specificity: {missing_taxa}")
    if missing_metadata:
        raise ValueError(f"Missing metadata columns for specificity: {missing_metadata}")

    feces = metadata[
        metadata["body_site_group"].eq(FECES) & metadata["inc_dec"].isin(["inc", "dec"])
    ][["profile_id", "study_id", "experiment_id", "inc_dec", "condition"]].copy()

    reported = taxa[taxa["taxon_id"].isin(taxon_ids)].copy()
    reported["inc_dec"] = reported["direction"].map(TAXA_DIRECTION_TO_PROFILE)
    occurrences = reported.merge(
        feces,
        on=["study_id", "experiment_id", "inc_dec"],
        how="inner",
        suffixes=("_reported", ""),
    )
    occurrences = occurrences.drop_duplicates(["profile_id", "taxon_id"])

    records: list[dict[str, object]] = []
    for taxon_id in taxon_ids:
        group = occurrences[occurrences["taxon_id"].eq(taxon_id)]
        pd_group = group[group["condition"].eq(PD_CONDITION)]
        other_group = group[~group["condition"].eq(PD_CONDITION)]
        pd_profiles = int(pd_group["profile_id"].nunique())
        other_profiles = int(other_group["profile_id"].nunique())
        total_profiles = pd_profiles + other_profiles
        log2_ratio = math.log2((pd_profiles + 0.5) / (other_profiles + 0.5))
        records.append(
            {
                "taxon_id": taxon_id,
                "pd_feces_profiles": pd_profiles,
                "other_feces_profiles": other_profiles,
                "pd_other_feces_counts": f"{pd_profiles}/{other_profiles}",
                "pd_feces_studies": int(pd_group["study_id"].nunique()),
                "other_feces_studies": int(other_group["study_id"].nunique()),
                "pd_specificity_fraction": (
                    pd_profiles / total_profiles if total_profiles else math.nan
                ),
                "pd_vs_other_feces_log2_ratio": log2_ratio,
            }
        )
    return pd.DataFrame.from_records(records)


def load_taxon_specificity(taxon_ids: pd.Index, metadata: pd.DataFrame) -> pd.DataFrame:
    """Load the full taxon occurrence table and build PD-vs-other feces specificity."""
    taxa = read_tsv(
        TAXA_PATH,
        usecols=[
            "study_id",
            "experiment_id",
            "direction",
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
        ],
    )
    return taxon_specificity_from_occurrences(taxa, metadata, taxon_ids)


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


def classify_effect(
    mean: float,
    low: float,
    high: float,
    support_studies: int,
    min_support_studies: int = MIN_SUPPORT_STUDIES,
) -> str:
    """Classify a signed influence contrast without forcing sparse or uncertain cells."""
    if support_studies < min_support_studies:
        return "low support"
    if not all(math.isfinite(value) for value in [mean, low, high]):
        return "mixed"
    if low > 0:
        return "inc"
    if high < 0:
        return "dec"
    return "mixed"


def filter_all_low_support_taxa(
    cells: pd.DataFrame,
    rows: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Drop species whose six trait cells are all low-support estimates."""
    if cells.empty or rows.empty:
        return cells.copy(), rows.copy(), 0
    all_low_support = cells.groupby("taxon_id", sort=False)["status"].apply(
        lambda statuses: statuses.eq("low support").all()
    )
    removed_taxa = set(all_low_support[all_low_support].index.astype(str))
    if not removed_taxa:
        return cells.copy(), rows.copy(), 0
    retained_cells = cells[~cells["taxon_id"].astype(str).isin(removed_taxa)].copy()
    retained_rows = rows[~rows["taxon_id"].astype(str).isin(removed_taxa)].copy()
    return retained_cells, retained_rows, len(removed_taxa)


def literature_reference_table() -> pd.DataFrame:
    """Return the curated DOI-backed PD taxa reference as a stable table."""
    return pd.DataFrame(PD_TAXA_LITERATURE_REFERENCE).copy()


def doi_markdown_links(value: object) -> str:
    """Format a semicolon-delimited DOI field as stable Markdown links."""
    dois = [doi.strip() for doi in str(value).split(";") if doi.strip()]
    if not dois:
        return "none"
    return "; ".join(f"[{doi}](https://doi.org/{doi})" for doi in dois)


def literature_methods_markdown() -> str:
    """Describe the targeted literature review and direction-comparison rules."""
    return "\n".join(
        [
            "## Methods: Literature Review and Direction Comparison",
            "",
            f"The literature direction table is a targeted narrative review updated "
            f"{LITERATURE_REVIEW_DATE}; it is not a systematic review or meta-analysis. "
            "PubMed, PubMed Central, and publisher records were searched using `Parkinson "
            "disease`, `gut microbiome`, each taxon name or synonym, and direction terms "
            "(`increased`, `decreased`, `enriched`, or `depleted`). Existing DOI-backed PD "
            "microbiome references were checked backward for additional direction evidence.",
            "",
            "Eligible evidence comprised human PD-versus-control fecal or intestinal "
            "microbiome studies and cross-cohort human meta-analyses that reported an explicit "
            "taxon direction. Animal-only, intervention-only, non-gut, and severity-only "
            "findings were excluded from increased/decreased calls. Association-only evidence "
            "is retained as `associated` and is not treated as directional.",
            "",
            "Literature direction is assigned at the most specific supported taxonomic level. "
            "An exact species result is preferred; otherwise a genus or family result is used "
            "as a flagged parent-level proxy. `variable` is used when eligible reports conflict "
            "at the matched level. The reference table records the matched term, reference rank, "
            "match type, evidence statement, and DOI list so broader-rank comparisons remain "
            "visible.",
            "",
            "BugSigDB direction is a recurrence summary, not a pooled abundance effect: "
            "`increased` when increased-profile count is greater than decreased-profile count, "
            "`decreased` for the reverse, and `balanced` for a tie. `agrees` and `differs` are "
            "assigned only against directional literature calls; balanced, variable, associated, "
            "and unlisted cases remain separate. These labels do not test statistical agreement "
            "and do not resolve cohort, medication, geography, sequencing, or taxonomic-version "
            "heterogeneity.",
        ]
    )


def normalize_taxon_name(value: object) -> str:
    """Normalize a taxon label for exact/prefix reference matching."""
    normalized = str(value).replace("[", "").replace("]", "").lower()
    return " ".join(normalized.split())


def reference_match_terms(row: pd.Series) -> list[str]:
    """Return all usable literature match terms for one reference row."""
    raw_terms = [row["literature_taxon"]]
    if pd.notna(row.get("match_terms")):
        raw_terms.extend(str(row["match_terms"]).split(";"))
    return [term.strip() for term in raw_terms if term and term.strip()]


def literature_rank_specificity(rank: object) -> int:
    """Rank literature rows so species-level matches beat broad family matches."""
    rank_text = str(rank).lower()
    if "species" in rank_text:
        return 5
    if "genus group" in rank_text:
        return 4
    if "genus" in rank_text:
        return 3
    if "family" in rank_text:
        return 2
    return 1


def match_literature_reference(
    taxon_name: object,
    reference: pd.DataFrame,
) -> dict[str, object]:
    """Find the most specific curated PD literature match for a displayed taxon."""
    name = normalize_taxon_name(taxon_name)
    best: tuple[tuple[int, int, int, int], dict[str, object]] | None = None
    for record in reference.to_dict("records"):
        row = pd.Series(record)
        for term in reference_match_terms(row):
            normalized_term = normalize_taxon_name(term)
            if not normalized_term:
                continue
            exact = name == normalized_term
            prefix = name.startswith(f"{normalized_term} ")
            synonym_contains_taxon = normalized_term.startswith(f"{name} ")
            if not (exact or prefix or synonym_contains_taxon):
                continue
            score = (
                3 if exact else 2 if synonym_contains_taxon else 1,
                len(normalized_term.split()),
                len(normalized_term),
                literature_rank_specificity(record["literature_rank"]),
            )
            match_type = (
                "exact" if exact else "synonym" if synonym_contains_taxon else "genus_or_parent"
            )
            candidate = (
                score,
                {
                    **record,
                    "literature_match_term": term,
                    "literature_match_type": match_type,
                },
            )
            if best is None or candidate[0] > best[0]:
                best = candidate

    if best is None:
        return {
            "literature_taxon": "not listed",
            "literature_rank": "not listed",
            "literature_direction": "not listed",
            "literature_evidence": (
                "No curated PD literature reference matched this displayed species."
            ),
            "doi_citations": "",
            "match_terms": "",
            "literature_match_term": "",
            "literature_match_type": "none",
        }
    return best[1]


def bugsigdb_pd_direction(row: pd.Series) -> str:
    """Classify the displayed taxon's BugSigDB PD occurrence balance."""
    inc = int(row.get("increased_profiles", 0))
    dec = int(row.get("decreased_profiles", 0))
    if inc > dec:
        return "increased"
    if dec > inc:
        return "decreased"
    return "balanced"


def sort_rows_by_bugsigdb_inc_dec(rows: pd.DataFrame) -> pd.DataFrame:
    """Sort from most increased-skewed to most decreased-skewed BugSigDB recurrence."""
    ordered = rows.copy()
    ordered["_inc_minus_dec"] = ordered["increased_profiles"].astype(int) - ordered[
        "decreased_profiles"
    ].astype(int)
    ordered["_total_direction_profiles"] = ordered["increased_profiles"].astype(int) + ordered[
        "decreased_profiles"
    ].astype(int)
    ordered = ordered.sort_values(
        [
            "log2_inc_dec_ratio",
            "_inc_minus_dec",
            "_total_direction_profiles",
            "taxon_name",
        ],
        ascending=[False, False, False, True],
        kind="stable",
    )
    return ordered.drop(columns=["_inc_minus_dec", "_total_direction_profiles"])


def direction_agreement(literature_direction: object, bugsigdb_direction: object) -> str:
    """Compare literature direction to the BugSigDB PD occurrence direction."""
    literature = str(literature_direction)
    bugsigdb = str(bugsigdb_direction)
    if literature == "not listed":
        return "not in literature reference"
    if literature not in {"increased", "decreased"}:
        return "literature non-directional"
    if bugsigdb == "balanced":
        return "balanced in BugSigDB"
    if literature == bugsigdb:
        return "agrees"
    return "differs"


def literature_comparison_label(literature_direction: object, agreement: object) -> str:
    """Format literature direction and BugSigDB agreement for a heatmap row."""
    direction = str(literature_direction) if literature_direction else "not listed"
    comparison = str(agreement) if agreement else "not in literature reference"
    direction_label = LITERATURE_LABELS.get(direction, LITERATURE_LABELS["not listed"])
    comparison_label = LITERATURE_COMPARISON_LABELS.get(comparison, comparison)
    return f"{direction_label}\n{comparison_label}"


def interpret_reference_row(row: pd.Series) -> str:
    """Summarize how to use one literature/BugSigDB/KG-Microbe comparison."""
    agreement = row["direction_agreement"]
    literature = row["literature_direction"]
    bugsigdb = row["bugsigdb_pd_direction"]
    if agreement == "agrees":
        lead = (
            f"BugSigDB PD occurrence agrees with literature ({literature}); treat this as a "
            "higher-priority taxon for trait-level follow-up."
        )
    elif agreement == "differs":
        lead = (
            f"BugSigDB PD occurrence is {bugsigdb}, but the literature reference is "
            f"{literature}; review cohort, taxonomy, and direction labels before interpreting."
        )
    elif agreement == "balanced in BugSigDB":
        lead = (
            f"Literature suggests {literature}, while BugSigDB PD occurrence is balanced; "
            "use as a sensitivity case rather than a directional claim."
        )
    elif agreement == "literature non-directional":
        lead = (
            "The literature match is association-oriented rather than cleanly increased or "
            "decreased; use it as supporting context, not a direction check."
        )
    else:
        lead = (
            "No curated PD literature match was found; treat the row as BugSigDB/KG-Microbe only."
        )

    direct_traits = str(row.get("kg_microbe_direct_six_traits", "none"))
    if direct_traits and direct_traits != "none":
        trait_text = f" KG-Microbe direct six-trait annotations: {direct_traits}."
    else:
        trait_text = " No direct KG-Microbe annotation among the six displayed traits."

    influence = row.get("strongest_influence_z", math.nan)
    if pd.notna(influence) and math.isfinite(float(influence)):
        strongest = (
            f" Strongest current influence: {row['strongest_influence_trait']} "
            f"({float(influence):+.2f} fixed-reference z, {row['strongest_influence_status']})."
        )
    else:
        strongest = ""
    return lead + trait_text + strongest


def build_pd_taxa_reference_table(
    rows: pd.DataFrame,
    cells: pd.DataFrame,
    literature_reference: pd.DataFrame,
) -> pd.DataFrame:
    """Join heatmap taxa to literature, BugSigDB direction, traits, and interpretation."""
    reference_matches = pd.DataFrame(
        [
            match_literature_reference(row.taxon_name, literature_reference)
            for row in rows.itertuples(index=False)
        ]
    )
    table = pd.concat([rows.reset_index(drop=True).copy(), reference_matches], axis=1)
    table["bugsigdb_pd_direction"] = table.apply(bugsigdb_pd_direction, axis=1)
    table["direction_agreement"] = table.apply(
        lambda row: direction_agreement(row["literature_direction"], row["bugsigdb_pd_direction"]),
        axis=1,
    )

    direct_traits = (
        cells[cells["has_trait_annotation"].fillna(False)]
        .groupby("taxon_id")["trait"]
        .apply(lambda traits: "; ".join(sorted({str(trait) for trait in traits})))
        .rename("kg_microbe_direct_six_traits")
    )
    direct_trait_count = (
        cells[cells["has_trait_annotation"].fillna(False)]
        .groupby("taxon_id")["trait_id"]
        .nunique()
        .rename("kg_microbe_direct_trait_count")
    )
    supported_cells = (
        cells[~cells["status"].eq("low support")]
        .groupby("taxon_id")["trait_id"]
        .nunique()
        .rename("supported_influence_cells")
    )
    finite_cells = cells[pd.notna(cells["mean_influence_contrast_z"])].copy()
    finite_cells["absolute_influence_z"] = finite_cells["mean_influence_contrast_z"].abs()
    strongest = finite_cells.loc[
        finite_cells.groupby("taxon_id")["absolute_influence_z"].idxmax(),
        ["taxon_id", "trait", "status", "mean_influence_contrast_z"],
    ].rename(
        columns={
            "trait": "strongest_influence_trait",
            "status": "strongest_influence_status",
            "mean_influence_contrast_z": "strongest_influence_z",
        }
    )

    table = (
        table.merge(direct_traits, on="taxon_id", how="left")
        .merge(direct_trait_count, on="taxon_id", how="left")
        .merge(supported_cells, on="taxon_id", how="left")
        .merge(strongest, on="taxon_id", how="left")
    )
    table["kg_microbe_direct_six_traits"] = table["kg_microbe_direct_six_traits"].fillna("none")
    table["kg_microbe_direct_trait_count"] = (
        table["kg_microbe_direct_trait_count"].fillna(0).astype(int)
    )
    table["supported_influence_cells"] = table["supported_influence_cells"].fillna(0).astype(int)
    table["interpretation"] = table.apply(interpret_reference_row, axis=1)

    columns = [
        "taxon_id",
        "taxon_name",
        "taxonomic_rank",
        "literature_taxon",
        "literature_rank",
        "literature_direction",
        "literature_match_term",
        "literature_match_type",
        "literature_evidence",
        "doi_citations",
        "bugsigdb_pd_direction",
        "increased_profiles",
        "decreased_profiles",
        "inc_dec_counts",
        "direction_agreement",
        "pd_other_feces_counts",
        "pd_feces_profiles",
        "other_feces_profiles",
        "pd_specificity_fraction",
        "pd_vs_other_feces_log2_ratio",
        "kg_microbe_direct_six_traits",
        "kg_microbe_direct_trait_count",
        "annotated_trait_count",
        "supported_influence_cells",
        "strongest_influence_trait",
        "strongest_influence_status",
        "strongest_influence_z",
        "interpretation",
    ]
    return table[columns]


def write_reference_interpretation(reference_table: pd.DataFrame, output_path: Path) -> None:
    """Write a compact Markdown interpretation of the joined PD taxa reference."""
    lines = [
        "# PD Taxa Reference Interpretation",
        "",
        "This table joins the displayed species influence rows to DOI-backed PD taxa literature, "
        "BugSigDB PD occurrence direction, KG-Microbe six-trait annotations, and the strongest "
        "current influence value.",
        "",
        "Rows follow the heatmap order: descending BugSigDB PD increased/decreased log2 "
        "profile-count ratio, with count-based and alphabetical tie-breakers.",
        "",
        "## Agreement Summary",
        "",
        "| Taxon | BugSigDB inc/dec | BugSigDB direction | Literature value | Match scope | "
        "Agreement | Strongest six-trait influence |",
        "| --- | ---: | --- | --- | --- | --- | --- |",
    ]
    for row in reference_table.itertuples(index=False):
        strongest = (
            f"{row.strongest_influence_trait} ({float(row.strongest_influence_z):+.2f}; "
            f"{row.strongest_influence_status})"
            if pd.notna(row.strongest_influence_z)
            else "not available"
        )
        match_scope = f"{row.literature_taxon}; {row.literature_rank}; {row.literature_match_type}"
        lines.append(
            f"| {row.taxon_name} | {row.inc_dec_counts} | {row.bugsigdb_pd_direction} | "
            f"{row.literature_direction} | {match_scope} | {row.direction_agreement} | "
            f"{strongest} |"
        )
    lines.extend(["", literature_methods_markdown(), "", "## Taxon Details", ""])
    for row in reference_table.itertuples(index=False):
        lines.extend(
            [
                f"### {row.taxon_name}",
                f"- Literature match: {row.literature_taxon} ({row.literature_direction}; "
                f"{row.literature_rank}; {row.literature_match_type}); "
                f"source(s): {doi_markdown_links(row.doi_citations)}",
                f"- BugSigDB PD direction: {row.bugsigdb_pd_direction} "
                f"({row.inc_dec_counts} increased/decreased profiles); "
                f"agreement: {row.direction_agreement}",
                f"- KG-Microbe direct six-trait annotations: {row.kg_microbe_direct_six_traits}",
                f"- Interpretation: {row.interpretation}",
                "",
            ]
        )
    output_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


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
    taxonomic_ranks: tuple[str, ...] = ("species",),
    min_support_studies: int = MIN_SUPPORT_STUDIES,
    remove_all_low_support_taxa: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, int]:
    """Build six-trait taxon cells, row annotations, and trait metadata."""
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
        & evidence["taxonomic_rank"].isin(taxonomic_ranks)
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
            status = classify_effect(mean, low, high, support_studies, min_support_studies)
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
                    "influence_basis": (
                        "direct + denominator" if annotated else "denominator only"
                    ),
                    "paired_profiles": len(pairs),
                    "paired_studies": int(pairs["study_id"].nunique()),
                    "taxon_support_profiles": int(paired_presence["profile_id"].nunique()),
                    "taxon_support_studies": support_studies,
                    "mean_influence_contrast_z": mean,
                    "ci95_low": low,
                    "ci95_high": high,
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

    row_metadata = pd.DataFrame(rows)
    cells_table = pd.DataFrame(cells)
    removed_low_support_taxa = 0
    if remove_all_low_support_taxa:
        cells_table, row_metadata, removed_low_support_taxa = filter_all_low_support_taxa(
            cells_table,
            row_metadata,
        )
    row_metadata = sort_rows_by_bugsigdb_inc_dec(row_metadata)
    specificity = load_taxon_specificity(pd.Index(row_metadata["taxon_id"]), metadata)
    row_metadata = row_metadata.merge(specificity, on="taxon_id", how="left")
    return cells_table, row_metadata, trait_metadata, removed_low_support_taxa


def plot_heatmap(
    cells: pd.DataFrame,
    rows: pd.DataFrame,
    traits: pd.DataFrame,
    output_path: Path,
    title: str = "PD species contributions across all six traits",
    taxon_unit: str = "species",
    method_note: str | None = None,
    top_margin_inches: float = 1.65,
    legend_y_inches: float = 1.38,
) -> None:
    """Plot trait influence with BugSigDB recurrence and literature-comparison annotations."""
    row_ids = rows["taxon_id"].tolist()
    trait_ids = traits["trait_id"].tolist()
    values = cells.pivot(index="taxon_id", columns="trait_id", values="mean_influence_contrast_z")
    labels = cells.pivot(index="taxon_id", columns="trait_id", values="status")
    direct_annotations = cells.pivot(
        index="taxon_id", columns="trait_id", values="has_trait_annotation"
    )
    values = values.reindex(index=row_ids, columns=trait_ids)
    labels = labels.reindex(index=row_ids, columns=trait_ids)
    direct_annotations = direct_annotations.reindex(index=row_ids, columns=trait_ids)
    finite = np.abs(values.to_numpy(float)[np.isfinite(values.to_numpy(float))])
    color_limit = max(float(np.quantile(finite, 0.95)) if finite.size else 0.0, 0.02)

    height = max(8.5, len(rows) * 0.32 + 3.4)
    fig = plt.figure(figsize=(15.6, height), facecolor="white")
    grid = fig.add_gridspec(
        1,
        3,
        width_ratios=[6.0, 0.85, 1.10],
        left=0.29,
        right=0.965,
        top=1 - top_margin_inches / height,
        bottom=1.15 / height,
        wspace=0.10,
    )
    ax = fig.add_subplot(grid[0, 0])
    row_ax = fig.add_subplot(grid[0, 1], sharey=ax)
    literature_ax = fig.add_subplot(grid[0, 2], sharey=ax)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(UNAVAILABLE_COLOR)
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
    row_labels = [
        f"{row.taxon_name} ({row.inc_dec_counts})" for row in rows.itertuples(index=False)
    ]
    ax.set_yticks(range(len(rows)), row_labels, fontsize=8.3)
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
            value = values.loc[taxon_id, trait_id]
            if direct_annotations.loc[taxon_id, trait_id]:
                ax.add_patch(
                    Rectangle(
                        (x - 0.49, y - 0.49),
                        0.98,
                        0.98,
                        facecolor="none",
                        edgecolor="#111827",
                        linewidth=0.85,
                    )
                )
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
            cell_label = (
                "low\nsupport" if status == "low support" else format_effect_value(float(value))
            )
            ax.text(
                x,
                y,
                cell_label,
                ha="center",
                va="center",
                fontsize=7.2 if status == "low support" else 7.6,
                color="white" if abs(value) >= 0.58 * color_limit else "#111827",
            )

    row_values = rows["log2_inc_dec_ratio"].clip(-ROW_COLOR_LIMIT, ROW_COLOR_LIMIT).to_numpy(float)
    row_ax.imshow(
        row_values[:, None], cmap=cmap, vmin=-ROW_COLOR_LIMIT, vmax=ROW_COLOR_LIMIT, aspect="auto"
    )
    row_ax.set_xticks([0], ["BugSigDB PD\ninc/dec"], fontsize=9)
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

    literature_directions = rows.get(
        "literature_direction",
        pd.Series(["not listed"] * len(rows), index=rows.index),
    ).fillna("not listed")
    literature_agreements = rows.get(
        "direction_agreement",
        pd.Series(["not in literature reference"] * len(rows), index=rows.index),
    ).fillna("not in literature reference")
    literature_ax.set_xlim(-0.5, 0.5)
    literature_ax.set_ylim(len(rows) - 0.5, -0.5)
    literature_ax.set_xticks([0], ["Literature vs\nBugSigDB"], fontsize=9)
    literature_ax.tick_params(left=False, labelleft=False, length=0, pad=6)
    literature_ax.set_xticks([-0.5, 0.5], minor=True)
    literature_ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    literature_ax.grid(which="minor", color="white", linewidth=1)
    literature_ax.tick_params(which="minor", length=0)
    for y, (direction, agreement) in enumerate(
        zip(literature_directions, literature_agreements, strict=True)
    ):
        direction = str(direction) if direction else "not listed"
        agreement = str(agreement) if agreement else "not in literature reference"
        color = LITERATURE_COMPARISON_COLORS.get(agreement, "#f3f4f6")
        label = literature_comparison_label(direction, agreement)
        literature_ax.add_patch(
            Rectangle(
                (-0.5, y - 0.5),
                1,
                1,
                facecolor=color,
                edgecolor="#b91c1c" if agreement == "differs" else "white",
                linewidth=2.2 if agreement == "differs" else 1,
            )
        )
        literature_ax.text(
            0,
            y,
            label,
            ha="center",
            va="center",
            fontsize=7.3,
            fontweight="bold" if agreement == "differs" else "normal",
            color="#111827",
        )
    for spine in literature_ax.spines.values():
        spine.set_visible(False)

    bar = fig.colorbar(
        image,
        ax=[ax, row_ax, literature_ax],
        pad=0.02,
        fraction=0.032,
        shrink=0.58,
        extend="both",
    )
    bar.set_label(
        f"{taxon_unit.capitalize()} influence on paired PD inc-dec trait contrast "
        "(fixed-reference z)"
    )
    fig.text(
        0.03,
        1 - 0.11 / height,
        title,
        fontsize=17,
        fontweight="bold",
        va="top",
    )
    subtitle = (
        "Cells show the study-weighted change in each paired increased-decreased trait contrast "
        f"attributable to removing one reported {taxon_unit} and recomputing its full "
        "non-MediaDive profile denominator. Taxon row labels include BugSigDB PD "
        "increased/decreased profile "
        "counts in parentheses. Outlined cells include direct removal of that trait; unoutlined "
        f"cells retain the denominator-mediated effect when the {taxon_unit} has no direct "
        "annotation. "
        "Cell text shows the signed fixed-reference z influence value; low-support cells remain "
        f"hatched. Inc/dec calls in the TSV require a 95% interval entirely above/below zero and "
        f"at least {MIN_SUPPORT_STUDIES} supporting studies. The right columns show raw PD "
        "inc/dec occurrence and the matched DOI-backed literature direction together with its "
        "BugSigDB agreement status. Rows marked DIFFERS have opposite recurrence directions."
    )
    if method_note:
        subtitle = f"{subtitle} {method_note}"
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
            Patch(facecolor="white", edgecolor="#111827", label="direct trait annotation"),
            Patch(facecolor="white", edgecolor="#64748b", hatch="///", label="low support"),
            Patch(facecolor=LITERATURE_COMPARISON_COLORS["agrees"], label="direction agrees"),
            Patch(
                facecolor=LITERATURE_COMPARISON_COLORS["differs"],
                edgecolor="#b91c1c",
                linewidth=1.5,
                label="direction differs",
            ),
            Patch(
                facecolor=LITERATURE_COMPARISON_COLORS["balanced in BugSigDB"],
                label="BugSigDB balanced",
            ),
        ],
        loc="upper center",
        bbox_to_anchor=(0.60, 1 - legend_y_inches / height),
        ncol=7,
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


def sha256_chunks(path: Path) -> list[str]:
    """Return a streamed SHA-256 digest in unambiguous eight-character chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    value = digest.hexdigest()
    return [value[offset : offset + 8] for offset in range(0, len(value), 8)]


def main() -> None:
    """Generate the six-trait influence heatmap and its auditable tables."""
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
    )
    literature_reference = literature_reference_table()
    reference_table = build_pd_taxa_reference_table(rows, cells, literature_reference)
    rows = rows.merge(
        reference_table[
            [
                "taxon_id",
                "literature_taxon",
                "literature_direction",
                "literature_match_term",
                "direction_agreement",
            ]
        ],
        on="taxon_id",
        how="left",
    )
    plot_heatmap(cells, rows, traits, OUTPUT_DIR / f"{STEM}.png")
    write_tsv(cells, OUTPUT_DIR / f"{STEM}_cells.tsv")
    write_tsv(rows, OUTPUT_DIR / f"{STEM}_taxa.tsv")
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
    write_tsv(rows[specificity_columns], OUTPUT_DIR / f"{STEM}_specificity.tsv")
    write_tsv(traits, OUTPUT_DIR / f"{STEM}_traits.tsv")
    write_tsv(literature_reference, OUTPUT_DIR / f"{LITERATURE_REFERENCE_NAME}.tsv")
    write_tsv(reference_table, OUTPUT_DIR / f"{REFERENCE_TABLE_STEM}_table.tsv")
    write_reference_interpretation(
        reference_table,
        OUTPUT_DIR / f"{REFERENCE_TABLE_STEM}_interpretation.md",
    )
    manifest = {
        "method": "fixed-reference leave-one-reported-species-out influence",
        "traits": len(traits),
        "species": len(rows),
        "cells": len(cells),
        "filled_cells": int(cells["mean_influence_contrast_z"].notna().sum()),
        "direct_annotation_cells": int(cells["has_trait_annotation"].sum()),
        "denominator_only_cells": int((~cells["has_trait_annotation"]).sum()),
        "species_removed_all_low_support": removed_low_support_taxa,
        "non_mediadive_columns": len(trait_columns),
        "minimum_support_studies": MIN_SUPPORT_STUDIES,
        "row_filter": "exclude species whose six trait cells are all low support",
        "row_order": (
            "descending BugSigDB PD log2(increased/decreased) profile-count ratio; "
            "ties by increased-minus-decreased count, total direction profiles, and taxon name"
        ),
        "literature_review_date": LITERATURE_REVIEW_DATE,
        "literature_review_method": "targeted narrative direction review",
        "displayed_row_annotations": [
            "BugSigDB PD increased/decreased recurrence",
            "literature direction and BugSigDB agreement",
        ],
        "pd_specificity_displayed": False,
        "pd_literature_reference_taxa": len(literature_reference),
        "displayed_species_with_literature_reference": int(
            reference_table["literature_direction"].ne("not listed").sum()
        ),
        "literature_bug_sig_direction_agreement_counts": {
            str(label): int(count)
            for label, count in reference_table["direction_agreement"].value_counts().items()
        },
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
    note = (
        "# PD Species-by-Six-Trait Influence Heatmap\n\n"
        f"![Heatmap]({STEM}.png)\n\n"
        f"Includes {len(rows)} reported species - the lowest reliable rank in the PD carrier "
        "data - and all six lollipop traits. Each trait cell is now specific to the species-trait "
        "combination: the analysis removes the species' complete non-MediaDive annotation vector "
        "from profiles where it was reported, recomputes the six log-fraction scores using the "
        "unchanged human-corpus reference, forms paired increased-minus-decreased experiment "
        "contrasts, averages experiments within study, and weights studies equally. Species "
        "whose six cells are all low support are omitted. Outlined cells include direct removal "
        "of the displayed trait; unoutlined cells show the denominator-mediated effect when the "
        "species has no direct annotation for that trait. The taxon row label also includes the "
        "BugSigDB PD increased/decreased profile count, such as `(0/13)`.\n\n"
        "Rows are sorted from the most increased-skewed to the most decreased-skewed "
        "BugSigDB PD recurrence using the log2 increased/decreased profile-count ratio; "
        "ties use the signed count difference, total support, and taxon name.\n\n"
        f"`inc` and `dec` require at least {MIN_SUPPORT_STUDIES} supporting studies and a 95% "
        "Student t interval entirely above or below zero. The heatmap prints signed numeric "
        "fixed-reference z influence values in supported cells; `mixed` is retained in the TSV "
        "when the interval crosses zero. Hatched `low support` cells do not receive a forced "
        "direction. The first right-hand annotation retains the raw increased/decreased profile "
        "occurrence counts in PD feces. The second annotation shows the DOI-backed literature "
        "direction and explicitly labels whether BugSigDB recurrence agrees, differs, or is "
        "balanced. `DIFFERS` marks *Mediterraneibacter torques* (BugSigDB 3/1 increased versus "
        "literature decreased) and *Escherichia coli* (BugSigDB 1/2 decreased versus literature "
        "increased); *Ruminococcus bromii* is marked balanced. PD-versus-other-feces specificity "
        "is retained in the companion audit TSV but is no longer plotted.\n\n"
        f"The generated reference table joins literature direction, DOI citations, BugSigDB PD "
        "occurrence direction, KG-Microbe six-trait annotations, strongest influence value, and "
        "an interpretation of agreement, disagreement, or balanced evidence.\n\n"
        f"{literature_methods_markdown()}\n\n"
        "This is a reported-species influence diagnostic, not a causal attribution or strain-level "
        "measurement. Removing a reported species changes both trait numerators and the nonlinear "
        "full denominator. Parent/child propagation and unresolved source assertions remain the "
        "limitations documented in #190; contrast orientation and trait semantics remain bounded "
        "by #188 and #189.\n\n"
        f"[Cell statistics]({STEM}_cells.tsv) | [Taxon row annotations]({STEM}_taxa.tsv) | "
        f"[PD specificity audit (not plotted)]({STEM}_specificity.tsv) | "
        f"[Trait summary]({STEM}_traits.tsv) | "
        f"[PD taxa reference]({REFERENCE_TABLE_STEM}_table.tsv) | "
        f"[Reference interpretation]({REFERENCE_TABLE_STEM}_interpretation.md) | "
        f"[Literature DOI list]({LITERATURE_REFERENCE_NAME}.tsv) | "
        f"[Manifest]({STEM}_manifest.json) | "
        f"[Vector figure]({STEM}.svg)\n"
    )
    (OUTPUT_DIR / f"{STEM}.md").write_text(note, encoding="utf-8")
    print(f"Wrote {len(rows)} species x {len(traits)} traits ({len(cells)} cells)")
    print(cells["status"].value_counts().to_string())
    print(OUTPUT_DIR / f"{STEM}.png")


if __name__ == "__main__":
    main()
