#!/usr/bin/env python3
"""Build the selected-figure slide deck (v6) from the figures and tables of a checkout.

Reads the slide figures, the earlier deck manifests (``slide_manifest_v1.tsv``,
``slide_manifest_v2.tsv``) and ``svg_sources/`` under
``data/derived/figure_review_2026-09-21/slides/``, and the rerun tables the slide
text quotes. Writes the HTML deck, the slide manifest, the copied slide assets and
``deck_validation_v6.json`` to the output directory, and updates the figure
provenance YAML in place. Many of these inputs are gitignored local outputs, so a
fresh clone cannot build the deck (#573).

The methods and Tufte plates and the cellulose_decr figures are read from this checkout,
where ``just figure-review-plates`` and ``just cellulose-decr-figures`` draw them; everything
else comes from the input root. Before writing, the build checks every committed file a
slide links or quotes against the revision it cites, and refuses when one differs (#637).

``just figure-review-deck`` runs this script, validates the provenance YAML against its LinkML
schema and prints the PDF with headless Chrome or Chromium.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import re
import shutil
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from statistics import median
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
DECK_DIR = Path("data/derived/figure_review_2026-09-21")
CELLULOSE_DIR = Path("data/analysis/cellulose_decr")
# The folders `just cellulose-decr-figures` draws from the committed cellulose_decr tables.
CELLULOSE_FIGURES = (f"{CELLULOSE_DIR}/figures/", f"{CELLULOSE_DIR}/figures_graphical/")
# Regenerated or untracked figure folders: their slides carry no GitHub figure link.
UNLINKED_FIGURES = ("data/derived/", *CELLULOSE_FIGURES, f"{CELLULOSE_DIR}/figures_tufte/")
REPO_URL = "https://github.com/CultureBotAI/microbial-trait-profiling"
# Where `just figure-review-plates` draws the methods and Tufte plates.
PLATES_DIR = "data/derived/tufte_plates/"
PROVENANCE_NAME = "figure_provenance_and_results.yaml"
# The deck's typed Table 1 numbers (cover note, slide captions, Tufte Plates III and VII)
# were written for this contrast classifier (569 oriented studies). Inputs from another
# version would mix with them, so the build refuses until #579 reads them from the tables.
TYPED_CONTRAST_VERSION = "contrast-orientation-v2"
CONTRAST_PROVENANCE = Path("data/derived/contrast_orientation/contrast_orientation_provenance.json")
COMMIT = "402ef5c25b7105c9aa73398085956dcf1f372bb4"
SHORT_COMMIT = COMMIT[:12]
BRANCH = "main"
GITHUB_BLOB = f"{REPO_URL}/blob/{COMMIT}"
# The cover's figure and code links, at COMMIT
COVER_FIGURE = (
    "trait-visualization/trait-profile-investigation/six-shared-traits/"
    "six_shared_traits_individual_conditions_condition_body_site_signed_direction_contrast_"
    "heatmap_page_01.png"
)
COVER_CODE = "scripts/create_six_shared_trait_cytoscape_heatmaps.py"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# The tables the slide text quotes; build_deck checks the committed ones against COMMIT (#637).
TABLES_READ: set[Path] = set()


def read_tsv(path: Path, quoted: bool = True) -> list[dict[str, str]]:
    """A table the slides read. A missing one stops the build rather than leaving the text
    that quotes it empty or zero (#639). A ``quoted`` table is checked against COMMIT;
    one read only for slides that link it as evidence is checked with those slides."""
    if not path.exists():
        raise FileNotFoundError(f"missing {path}")
    if quoted:
        TABLES_READ.add(path)
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def listing(names: list[str]) -> str:
    """'A', 'A and B' or 'A, B and C'."""
    if len(names) < 2:
        return "".join(names)
    return f"{', '.join(names[:-1])} and {names[-1]}"


def arrow(delta: float) -> str:
    return "↑" if delta > 0 else "↓" if delta < 0 else "="


def fragile_sentence(metrics: dict) -> str:
    """' In the cellulose deep-dive, the HIV result depends on ...' or '' (#590)."""
    fragile = metrics["cellulose_statements"]["fragile"]
    return f" In the cellulose deep-dive, {fragile}." if fragile else ""


def taxon_family_text(metrics: dict) -> str:
    survivors = metrics["pd_taxon_survivors"]
    family = metrics["pd_taxon_family"]
    if survivors == 0:
        return f"no taxon survives the separate {family}-test family"
    return f"{survivors} of the separate {family}-test taxon family survive BH"


def check_contrast_version(repo: Path, allow_mismatch: bool) -> str | None:
    """Refuse contrast-orientation inputs from another classifier than the typed text (#579).

    A mismatch that is allowed returns the reason the deck is a review build (#637).
    """
    path = repo / CONTRAST_PROVENANCE
    version = (
        json.loads(path.read_text(encoding="utf-8")).get("classifier_version")
        if path.exists()
        else None
    )
    if version != TYPED_CONTRAST_VERSION and not allow_mismatch:
        raise SystemExit(
            f"{path} records {version}, but the deck's typed Table 1 numbers were written for "
            f"{TYPED_CONTRAST_VERSION}; a rebuild would mix the two (#579). Pass "
            "--allow-contrast-version-mismatch only for a review build."
        )
    if version != TYPED_CONTRAST_VERSION:
        return f"contrast-orientation inputs from {version}, not {TYPED_CONTRAST_VERSION}"
    return None


def stale_figure_ids(prov_doc: Mapping[str, Any], in_deck: set[str]) -> list[str]:
    """Figure ids that statistical results cite but whose figures are not in the deck.

    Renumbering would otherwise attach those results to whichever figure takes the old
    id's slide number (#638).
    """
    kept = {
        rec["figure_id"]
        for rec in prov_doc.get("figures", [])
        if rec["figure_file_path"] in in_deck
    }
    return sorted(
        {
            fid
            for result in prov_doc.get("statistical_results", [])
            for fid in result.get("figure_ids") or []
            if fid not in kept
        }
    )


def cited_refs() -> list[str]:
    """The revisions whose committed files the slides link or quote (#637): the rerun,
    the commit that drew the cellulose_decr figures and the one that drew the plates."""
    refs = [COMMIT]
    for ref in (cellulose_revision(REPO_ROOT), plate_revision(REPO_ROOT)):
        if ref and ref not in refs:
            refs.append(ref)
    return refs


def load_rerun_metrics(repo: Path) -> dict:
    """Read fresh rerun TSVs to ensure slide titles, takeaways, and cover table match rerun data."""
    metrics: dict = {}

    # 1. Six-shared-traits direction contrasts
    contrasts_path = (
        repo / "trait-visualization/trait-profile-investigation/six-shared-traits/"
        "six_shared_trait_condition_body_site_direction_contrasts.tsv"
    )
    contrasts = read_tsv(contrasts_path)
    metrics["contrasts"] = contrasts

    # Index by (analysis_level, condition_label, body_site_group, trait_label)
    by_key = {}
    for row in contrasts:
        key = (
            row.get("view_key", ""),
            row.get("condition", ""),
            row.get("body_site_group", ""),
            row.get("trait_label", ""),
        )
        by_key[key] = row
    metrics["contrasts_by_key"] = by_key

    # Count individual and grouped FDR survivors
    ind_rows = [
        r for r in contrasts if r.get("view_key") == "individual_conditions_direction_contrast"
    ]
    grp_rows = [r for r in contrasts if r.get("view_key") == "grouped_direction_contrast"]
    ind_sig = [r for r in ind_rows if r.get("selected_after_bh") == "True"]
    grp_sig = [r for r in grp_rows if r.get("selected_after_bh") == "True"]
    metrics["ind_total_tests"] = len(ind_rows)
    metrics["grp_total_tests"] = len(grp_rows)
    metrics["ind_sig_count"] = len(ind_sig)
    metrics["grp_sig_count"] = len(grp_sig)
    metrics["ind_sig_rows"] = sorted(ind_sig, key=lambda r: float(r.get("fisher_q_value", "1")))

    # 2. PD condition-taxon-trait network tests
    pd_net_path = (
        repo / "trait-visualization/trait-profile-investigation/PD-cellulose/"
        "pd_condition_taxon_trait_network_condition_trait_tests.tsv"
    )
    pd_net_rows = read_tsv(pd_net_path)
    metrics["pd_net_total"] = len(pd_net_rows)
    for r in pd_net_rows:
        if (
            r.get("condition") == "Parkinson's disease"
            and r.get("trait") == "moderately halophilic"
        ):
            metrics["pd_halo_wilcoxon"] = r
            break

    for r in pd_net_rows:
        if (
            r.get("condition") == "Parkinson's disease"
            and r.get("trait") == "cellulose catabolic process"
        ):
            metrics["pd_cellulose_wilcoxon"] = r
            break
    if pd_net_rows:
        metrics["pd_net_lead"] = min(pd_net_rows, key=lambda r: float(r["q_value"]))
    pd_dir = repo / "trait-visualization/trait-profile-investigation/PD-cellulose"
    family = [
        r
        for r in read_tsv(pd_dir / "pd_condition_taxon_trait_network_condition_taxon_tests.tsv")
        if r.get("passes_min_taxon_studies") == "True"
    ]
    metrics["pd_taxon_family"] = len(family)
    metrics["pd_taxon_survivors"] = sum(r.get("selected_after_bh") == "True" for r in family)
    drivers = sorted(
        (
            r
            for r in read_tsv(pd_dir / "pd_cellulose_taxon_drivers.tsv")
            if r.get("cellulose_catabolic_process_present") == "True"
            and float(r["net_decreased_minus_increased_profiles"]) > 0
        ),
        key=lambda r: (-float(r["net_decreased_minus_increased_profiles"]), r["taxon_name"]),
    )
    metrics["pd_cellulose_drivers"] = [r["taxon_name"] for r in drivers[:3]]
    metrics["pd_six_trait_summary"] = read_tsv(pd_dir / "pd_six_trait_direction_summary.tsv")
    model_rows = read_tsv(pd_dir / "pd_cellulose_interaction_model.tsv")
    for outcome, key in (
        ("cellulose_z", "pd_interaction"),
        ("cellulose_positive_z", "pd_interaction_positive"),
    ):
        interaction = [
            r for r in model_rows if r["outcome"] == outcome and r["term"] == "pd_x_decreased"
        ]
        if len(interaction) != 1:
            raise ValueError(
                f"pd_cellulose_interaction_model.tsv needs one {outcome} pd_x_decreased row"
            )
        metrics[key] = interaction[0]

    # 3. PD ecological reporter combined
    pd_eco_path = (
        repo / "trait-visualization/trait-profile-investigation/PD-cellulose/"
        "pd_ecological_reporter_combined_condition_trait_summary.tsv"
    )
    pd_eco_rows = read_tsv(pd_eco_path)
    for r in pd_eco_rows:
        if r.get("condition") == "Parkinson's disease":
            if r.get("trait") == "cellulose catabolic process":
                metrics["pd_eco_cellulose"] = r
            elif r.get("trait") == "ecological-reporter composite":
                metrics["pd_eco_composite"] = r

    # 4. Sankey display edges stats
    sankey_stats = {}
    sankey_dir = repo / "trait-visualization/trait-profile-investigation/six-shared-traits"
    for edges_file in sorted(sankey_dir.glob("*_display_edges.tsv")):
        # Only the slides that link a display-edges table as evidence quote it (#637).
        rows = read_tsv(edges_file, quoted=False)
        if not rows:
            continue
        total_support = sum(
            int(float(r.get("focus_enriched_profile_count") or r.get("profile_count") or "0"))
            for r in rows
        )
        feces_support = sum(
            int(float(r.get("focus_enriched_profile_count") or r.get("profile_count") or "0"))
            for r in rows
            if r.get("body_site_group") == "Feces"
        )
        conds = {r.get("condition", "") for r in rows if r.get("condition")}
        sites = {r.get("body_site_group", "") for r in rows if r.get("body_site_group")}
        pct = round(100.0 * feces_support / total_support) if total_support else 0
        sankey_stats[edges_file.name] = {
            "links": len(rows),
            "total_support": total_support,
            "feces_support": feces_support,
            "feces_pct": pct,
            "n_conds": len(conds),
            "n_sites": len(sites),
        }
    metrics["sankey_stats"] = sankey_stats

    # 5. Study-level paired contrasts by stratum, quoted by Tufte Plate II's slide
    paired_dir = repo / "trait-visualization/trait-profile-investigation"
    metrics["paired_contrasts"] = {
        stratum: read_tsv(paired_dir / f"trait_profile_paired_contrast_trait_stats_{stratum}.tsv")
        for stratum in ("pooled", "feces", "feces_excluded", "saliva")
    }
    # 6. Oxygen-switch scores, quoted by Tufte Plate V's slide
    metrics["oxygen_switch"] = read_tsv(
        paired_dir / "taxa-oxygen-switches" / "condition_oxygen_switch_scores.tsv"
    )
    # 7. Tufte Plate I's scales and Plate VI's zero-trait audit, quoted by their slides (#636,
    # #640)
    paired = metrics["paired_contrasts"]
    metrics["pooled_pairs"] = max(int(r["pair_count"]) for r in paired["pooled"])
    metrics["pooled_studies"] = max(int(float(r["study_count"])) for r in paired["pooled"])
    metrics["feces_studies"] = max(int(float(r["study_count"])) for r in paired["feces"])
    metrics["ind_conditions"] = len({r["condition"] for r in ind_rows})
    metrics["ind_sites"] = len({r["body_site_group"] for r in ind_rows})
    metrics["ind_median_studies"] = median(int(r["study_count"]) for r in ind_rows)
    zero_trait = read_tsv(paired_dir / "taxa-coverage" / "top_zero_trait_bugsigdb_taxa.tsv")
    metrics["zero_trait_taxa"] = len(zero_trait)
    metrics["zero_trait_occurrences"] = sum(int(r["signature_occurrences"]) for r in zero_trait)
    # 8. Corpus studies, to say when a slide's prior-result link is one of them (#641)
    metrics["studies"] = {
        r["pmid"]: r for r in read_tsv(repo / "data/derived/studies.tsv") if r.get("pmid")
    }

    return metrics


def p_text(value: float) -> str:
    """A p-value to two significant digits, or more where fewer would round across 0.05."""
    for digits in range(2, 8):
        text = f"{value:#.{digits}g}"
        if (float(text) < 0.05) == (value < 0.05):
            return text
    return repr(value)


def signed(value: float, digits: int = 1) -> str:
    """A signed number that never prints as "-0.0"."""
    text = f"{value:+.{digits}f}"
    return f"{0:.{digits}f}" if float(text) == 0 else text


def plural(count: int, noun: str) -> str:
    return f"{count} {noun}" if count == 1 else f"{count} {noun}s"


def tested_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Traits with a study-level q-value, the family the stratum's BH ran over."""

    def has_q(row: dict[str, str]) -> bool:
        try:
            return not math.isnan(float(row.get("study_q_value") or "nan"))
        except ValueError:
            return False

    return [r for r in rows if has_q(r)]


def study_bh_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Traits passing study-level BH q ≤ 0.05, by mean Δz (increased minus decreased)."""
    passing = [r for r in tested_rows(rows) if float(r["study_q_value"]) <= 0.05]
    return sorted(passing, key=mean_dz)


def mean_dz(row: dict[str, str]) -> float:
    return float(row["study_mean_increased_minus_decreased"])


def trait_names(rows: list[dict[str, str]]) -> str:
    """Trait names as Plate II shows them, e.g. "raffinose (carbon source)" (#607, #617)."""

    def name(row: dict[str, str]) -> str:
        # trait_label ends with the trait's ID, e.g. "xylan catabolic process GO:0045493"
        label = row["trait_label"].rsplit(" ", 1)[0]
        carbon_source = row.get("trait", "").startswith("METPO:2000006|")
        return f"{label} (carbon source)" if carbon_source else label

    return listing([name(r) for r in rows])


def paired_contrast_title(metrics: dict) -> str:
    """Tufte Plate II's slide title, read from the pooled table (#607)."""
    passing = study_bh_rows(metrics["paired_contrasts"]["pooled"])
    verb = "passes" if len(passing) == 1 else "pass"
    return (
        f"Tufte Plate II: {plural(len(passing), 'trait')} {verb} study-level BH pooled, by "
        "direction and body site"
    )


def paired_contrast_text(metrics: dict) -> str:
    """Tufte Plate II's takeaway, read from the stratum tables (#556, #608)."""
    contrasts = metrics["paired_contrasts"]
    passing = study_bh_rows(contrasts["pooled"])
    lower = [r for r in passing if mean_dz(r) < 0]
    higher = [r for r in passing if mean_dz(r) > 0][::-1]
    level = [r for r in passing if mean_dz(r) == 0]
    subject = "the annotation z-score"
    sides = []
    if lower:
        sides.append(
            f"in {len(lower)} {subject} is higher in decreased signatures "
            f"(mean Δz < 0: {trait_names(lower)})"
        )
    if higher:
        sides.append(
            f"in {len(higher)} {'it' if lower else subject} is higher in increased signatures "
            f"(Δz > 0), led by {trait_names(higher[:3])}"
        )
    if level:
        sides.append(f"in {len(level)} the mean Δz is 0")
    tested = len(tested_rows(contrasts["pooled"]))
    feces = len(study_bh_rows(contrasts["feces"]))
    head = (
        f"Of {plural(tested, 'trait')} tested pooled, {len(passing)} "
        f"{'passes' if len(passing) == 1 else 'pass'} study-level BH q ≤ 0.05 ({feces} in Feces)"
    )
    return f"{head}: {'; '.join(sides)}." if sides else f"{head}."


def paired_contrast_caveat(metrics: dict) -> str:
    """Tufte Plate II's caveat: the strata without Feces, and saliva's nominal nitrate sign."""
    contrasts = metrics["paired_contrasts"]

    def studies(stratum: str) -> int:
        return max(int(float(r["study_count"])) for r in contrasts[stratum])

    def nitrate(stratum: str) -> dict[str, str]:
        trait = "METPO:2000103|pathways:nitrate_reduction"
        rows = [r for r in contrasts[stratum] if r["trait"] == trait]
        if len(rows) != 1:
            raise ValueError(f"the {stratum} paired-contrast table needs one {trait} row")
        return rows[0]

    saliva, feces = nitrate("saliva"), nitrate("feces")
    t_saliva, t_feces = float(saliva["study_t_statistic"]), float(feces["study_t_statistic"])
    sign = "the opposite sign to" if (t_saliva > 0) != (t_feces > 0) else "the same sign as"
    if float(saliva["study_q_value"]) <= 0.05:
        support = "at BH q ≤ 0.05"
    elif float(saliva["study_p_value"]) <= 0.05:
        support = "at a nominal p only"
    else:
        support = "without reaching p ≤ 0.05"
    excluded = contrasts["feces_excluded"]
    return (
        f"Without Feces ({studies('feces_excluded')} studies), {len(study_bh_rows(excluded))} "
        f"of {plural(len(tested_rows(excluded)), 'trait')} pass study BH q ≤ 0.05; in "
        f"Saliva ({studies('saliva')} studies), {len(study_bh_rows(contrasts['saliva']))} of "
        f"{len(tested_rows(contrasts['saliva']))}. Saliva's nitrate reduction has {sign} "
        f"Feces' {support} (t = {t_saliva:+.2f}, p = "
        f"{p_text(float(saliva['study_p_value']))}, q = {p_text(float(saliva['study_q_value']))}; "
        f"Feces t = {t_feces:+.2f}, q = {p_text(float(feces['study_q_value']))}), and no test "
        "compares the two sites."
    )


def oxygen_switch_strata(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Plate V's strata, as build_plate5 selects them: the 10 highest-ranked scoring above
    0.20 (rank_score weighs in study support) and the 10 lowest-scoring below -0.15."""

    def score(row: dict[str, str]) -> float:
        return float(row["aerobe_up_anaerobe_down_score"])

    ranked = sorted(rows, key=lambda r: -float(r["rank_score"]))
    up = [r for r in ranked if score(r) > 0.20][:10]
    down = sorted((r for r in rows if score(r) < -0.15), key=score)[:10]
    return up + down


def oxygen_switch_text(metrics: dict) -> str:
    """Tufte Plate V's takeaway: how the strata it shows split by site (#610)."""
    strata = oxygen_switch_strata(metrics["oxygen_switch"])

    def score(row: dict[str, str]) -> float:
        return float(row["aerobe_up_anaerobe_down_score"])

    fecal = [r for r in strata if r["body_site_display"] == "Feces"]
    other = [r for r in strata if r["body_site_display"] != "Feces"]
    exceptions = [r for r in fecal if score(r) < 0] + [r for r in other if score(r) > 0]

    def group(count: int, total: int, kind: str) -> str:
        noun, verb = ("stratum", "scores") if total == 1 else ("strata", "score")
        return f"{count} of {total} {kind} {noun} {verb}"

    counts = []
    if fecal:
        counts.append(
            f"{group(sum(score(r) > 0 for r in fecal), len(fecal), 'fecal')} aerobe-up / "
            "anaerobe-down"
        )
    if other:
        counts.append(
            f"{group(sum(score(r) < 0 for r in other), len(other), 'other')} anaerobe-up / "
            "aerobe-down"
        )
    shown = f"{len(strata)} {'stratum' if len(strata) == 1 else 'strata'}"
    text = f"Of the {shown} shown, {' and '.join(counts)}."
    if exceptions:
        listed = listing(
            [f"{r['condition']} · {r['body_site_display']} {score(r):+.2f}" for r in exceptions]
        )
        text += f" The exceptions: {listed}."
    return text


def plate_one_text(metrics: dict) -> tuple[str, str]:
    """Tufte Plate I's takeaway and caveat, read from the tables the plate draws (#636)."""
    pooled = metrics["paired_contrasts"]["pooled"]
    study_level = sum(
        float(r["study_q_value"]) <= 0.05
        for r in pooled
        if r["trait"] in TABLE1_TRAIT_IDS and r.get("study_q_value")
    )
    every = "All 6" if study_level == 6 else f"{study_level} of 6"
    halophily = metrics["pd_halo_wilcoxon"]
    cells = metrics["ind_sig_count"]
    takeaway = (
        "Traces the 6 Table 1 traits across 6 analytical scales "
        f"({metrics['pooled_pairs']:,} paired experiments → {metrics['pooled_studies']:,} "
        f"studies → {metrics['feces_studies']:,} fecal studies → study sign tests → "
        f"{metrics['ind_total_tests']:,} individual condition\N{MULTIPLICATION SIGN}site profile "
        f"tests → {metrics['pd_net_total']} fecal study Wilcoxon tests)."
    )
    caveat = (
        f"{every} traits {'passes' if study_level == 1 else 'pass'} study-level BH q ≤ 0.05 "
        f"across {metrics['pooled_studies']:,} "
        f"studies; across {metrics['ind_conditions']} conditions \N{MULTIPLICATION SIGN} "
        f"{metrics['ind_sites']} sites (median {metrics['ind_median_studies']:g} "
        f"{'study' if metrics['ind_median_studies'] == 1 else 'studies'} per cell), "
        f"{plural(cells, 'individual profile cell')} {'is' if cells == 1 else 'are'} selected "
        "after BH (q ≤ 0.05 and the study-support rule), and PD fecal halophily reaches "
        "Wilcoxon p = "
        f"{p_text(float(halophily['p_value']))} (BH q = {p_text(float(halophily['q_value']))})."
    )
    return takeaway, caveat


def crc_oxygen_text(rows: list[dict[str, str]]) -> str:
    """'Colorectal cancer splits ... (Feces +0.68 vs Lower GI tract -0.79). ', from the
    oxygen-switch table, or '' when either stratum is missing (#641)."""
    scores = {
        r["body_site_display"]: float(r["aerobe_up_anaerobe_down_score"])
        for r in rows
        if r["condition"] == "Colorectal cancer"
    }
    if "Feces" not in scores or "Lower GI tract" not in scores:
        return ""
    return (
        "Colorectal cancer splits by sampling compartment in O₂ switch score (Feces "
        f"{scores['Feces']:+.2f} vs Lower GI tract {scores['Lower GI tract']:+.2f}). "
    )


def crc_mucosa_shown(metrics: dict) -> str:
    """Whether Plate V shows the CRC Intestinal / GI mucosa stratum (#641)."""
    shown = any(
        r["condition"] == "Colorectal cancer" and r["body_site_display"] == "Intestinal / GI mucosa"
        for r in oxygen_switch_strata(metrics["oxygen_switch"])
    )
    return f"the second is {'' if shown else 'not '}among the strata shown"


def prior_label(url: str, studies: Mapping[str, dict[str, str]]) -> str:
    """A prior-result link's label, saying when the study is in the analysed corpus (#641)."""
    match = re.search(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", url)
    if match is None:
        return url
    pmid = match.group(1)
    study = studies.get(pmid)
    if study is None:
        return f"PubMed {pmid}"
    doi = f" · DOI {study['doi']}" if study.get("doi") else ""
    return (
        f"PubMed {pmid} ({study.get('publication_year', '')}){doi} · a study in the analysed "
        f"corpus (bsdb:{pmid})"
    )


def update_previous_slide_text(row: dict[str, str], metrics: dict) -> tuple[str, str, str]:
    """Update title, takeaway, and caveat for previous slides where rerun metrics changed."""
    fig_path = row["figure_file_path"]
    title = row["title"]
    takeaway = row["takeaway"]
    caveat = row["caveat"]
    by_key = metrics["contrasts_by_key"]

    # Check if this is one of the 18 appendix Sankeys with a display_edges.tsv
    ev_path = row.get("evidence_path", "")
    if ev_path.endswith("_display_edges.tsv"):
        st = metrics["sankey_stats"].get(Path(ev_path).name)
        if st:
            # Update numbers if needed
            if "six_shared_traits_increased_profile_sankey.png" in fig_path:
                title = (
                    f"Feces accounts for {st['feces_pct']}% of displayed increased-profile support"
                )
                takeaway = (
                    f"{st['feces_support']:,} of {st['total_support']:,} displayed positive-z "
                    "memberships "
                    f"connect to Feces ({st['feces_pct']}%); {st['links']} links span "
                    f"{st['n_conds']} condition groups and {st['n_sites']} sites."
                )
            elif "six_shared_traits_decreased_profile_sankey.png" in fig_path:
                title = (
                    f"Feces accounts for {st['feces_pct']}% of displayed decreased-profile support"
                )
                takeaway = (
                    f"{st['feces_support']:,} of {st['total_support']:,} displayed positive-z "
                    "memberships "
                    f"connect to Feces ({st['feces_pct']}%); {st['links']} links span "
                    f"{st['n_conds']} condition groups and {st['n_sites']} sites."
                )
            elif "six_shared_traits_increased_trait_lanes_profile_sankey.png" in fig_path:
                title = f"Six trait lanes retain {st['links']} increased-profile links"
                takeaway = (
                    f"{st['feces_support']:,} of {st['total_support']:,} displayed positive-z "
                    "memberships "
                    f"connect to Feces ({st['feces_pct']}%); {st['links']} links span "
                    f"{st['n_conds']} condition groups and {st['n_sites']} sites."
                )
            elif "six_shared_traits_decreased_trait_lanes_profile_sankey.png" in fig_path:
                title = f"Six trait lanes retain {st['links']} decreased-profile links"
                takeaway = (
                    f"{st['feces_support']:,} of {st['total_support']:,} displayed positive-z "
                    "memberships "
                    f"connect to Feces ({st['feces_pct']}%); {st['links']} links span "
                    f"{st['n_conds']} condition groups and {st['n_sites']} sites."
                )
            elif (
                "six_shared_traits_increased_individual_conditions_profile_sankey.png" in fig_path
                or "six_shared_traits_decreased_individual_conditions_profile_sankey.png"
                in fig_path
            ):
                title = (
                    f"{st['n_conds']} conditions connect to {st['n_sites']} body-site groups "
                    "across six traits"
                )
                takeaway = (
                    f"{st['links']} selected trait links connect {st['n_conds']} individual "
                    "conditions to "
                    f"{st['n_sites']} body-site groups; {st['feces_pct']}% of displayed positive-z "
                    "memberships are fecal."
                )
            else:
                # Single-trait sankeys
                trait_short_map = {
                    "facultatively_aerobic": "facultative-aerobe",
                    "cellulose_catabolic_process": "cellulose",
                    "coccus_shaped": "coccus",
                    "sphere_shaped": "sphere-shaped",
                    "arsenate_reduction": "arsenate-reduction",
                    "moderately_halophilic": "halophily",
                }
                for prefix, t_short in trait_short_map.items():
                    if Path(fig_path).name.startswith(prefix):
                        title = (
                            f"Feces accounts for {st['feces_pct']}% of displayed {t_short} support"
                        )
                        takeaway = (
                            f"{st['feces_support']:,} of {st['total_support']:,} displayed "
                            "positive-z memberships "
                            f"connect to Feces ({st['feces_pct']}%); {st['links']} links span "
                            f"{st['n_conds']} condition groups and {st['n_sites']} sites."
                        )
                        break

    # Specific Main slides with rerun numbers
    if fig_path.endswith(
        "six_shared_traits_individual_conditions_condition_body_site_signed_direction_contrast_heatmap_page_01.png"
    ):
        pd_cell = by_key.get(
            (
                "individual_conditions_direction_contrast",
                "Parkinson's disease",
                "Feces",
                "cellulose catabolic process",
            )
        )
        hiv_cell = by_key.get(
            (
                "individual_conditions_direction_contrast",
                "HIV infection",
                "Feces",
                "cellulose catabolic process",
            )
        )
        crc_cell = by_key.get(
            (
                "individual_conditions_direction_contrast",
                "Colorectal cancer",
                "Feces",
                "cellulose catabolic process",
            )
        )
        cells = [("PD", pd_cell), ("HIV", hiv_cell), ("colorectal cancer", crc_cell)]
        if pd_cell and hiv_cell and crc_cell:
            favored = [
                name
                for name, cell in cells
                if cell["selected_after_bh"] == "True"
                and float(cell["increased_minus_decreased_positive_fraction"]) < 0
            ]
            other_way = [
                name
                for name, cell in cells
                if cell["selected_after_bh"] == "True"
                and float(cell["increased_minus_decreased_positive_fraction"]) >= 0
            ]
            if favored:
                title = (
                    "More cellulose-positive profiles among decreased fecal signatures in "
                    f"{listing(favored)}"
                )
            elif other_way:
                title = (
                    "No cellulose fecal profile contrast passes BH toward decreased signatures; "
                    f"{listing(other_way)} {'passes' if len(other_way) == 1 else 'pass'} the "
                    "other way"
                )
            else:
                title = (
                    "No cellulose fecal profile contrast in PD, HIV or colorectal cancer passes BH"
                )
            takeaway = "Profile-level Fisher BH q for cellulose-positive profiles: " + ", ".join(
                f"{name} {float(cell['fisher_q_value']):.5f} "
                f"({arrow(float(cell['increased_minus_decreased_positive_fraction']))})"
                for name, cell in cells
            )
            caveat = (
                "PD = Parkinson\N{RIGHT SINGLE QUOTATION MARK}s disease. Positive z-scores "
                "describe trait annotations; "
                f"profile-level Fisher tests use a {metrics['ind_total_tests']:,}-row BH family."
                f"{fragile_sentence(metrics)}"
            )

    elif fig_path.endswith(
        "six_shared_traits_grouped_condition_body_site_signed_direction_contrast_heatmap.png"
    ):
        n_grp_sig = metrics["grp_sig_count"]
        n_grp_tot = metrics["grp_total_tests"]
        num_words = {
            1: "One",
            2: "Two",
            3: "Three",
            4: "Four",
            5: "Five",
            6: "Six",
            7: "Seven",
            8: "Eight",
            9: "Nine",
            10: "Ten",
            11: "Eleven",
            12: "Twelve",
            13: "Thirteen",
            14: "Fourteen",
            15: "Fifteen",
        }
        w = num_words.get(n_grp_sig, str(n_grp_sig))
        title = f"{w} grouped condition\N{EN DASH}site trait associations pass FDR"
        takeaway = (
            f"{w} cells pass the {n_grp_tot:,}-row Benjamini\N{EN DASH}Hochberg false "
            "discovery rate (FDR) correction; grouped labels summarize related conditions."
        )

    elif fig_path.endswith("six_shared_traits_combined_shared_high_support_profile_sankey.png"):
        r = by_key.get(
            (
                "grouped_direction_contrast",
                "Neurologic / cognitive disease",
                "Feces",
                "moderately halophilic",
            )
        )
        if r:
            inc_p = round(float(r["increased_positive_fraction"]) * 100)
            dec_p = round(float(r["decreased_positive_fraction"]) * 100)
            qv = float(r["fisher_q_value"])
            takeaway = (
                f"Neurologic/cognitive group: halophily z > 0 in {inc_p}% of increased and "
                f"{dec_p}% of decreased fecal profiles (grouped-screen q = {qv:.4f})."
            )

    elif fig_path.endswith("six_shared_traits_increased_shared_high_support_profile_sankey.png"):
        r = by_key.get(
            ("grouped_direction_contrast", "Infectious disease", "Feces", "coccus shaped")
        )
        if r:
            inc_p = round(float(r["increased_positive_fraction"]) * 100)
            dec_p = round(float(r["decreased_positive_fraction"]) * 100)
            qv = float(r["fisher_q_value"])
            takeaway = (
                f"Infectious-disease group: coccus z > 0 in {inc_p}% of increased and "
                f"{dec_p}% of decreased fecal profiles (grouped-screen q = {qv:.4f})."
            )

    elif fig_path.endswith("six_shared_traits_decreased_shared_high_support_profile_sankey.png"):
        r = by_key.get(
            (
                "grouped_direction_contrast",
                "Neurologic / cognitive disease",
                "Feces",
                "cellulose catabolic process",
            )
        )
        if r:
            inc_p = round(float(r["increased_positive_fraction"]) * 100)
            dec_p = round(float(r["decreased_positive_fraction"]) * 100)
            qv = float(r["fisher_q_value"])
            takeaway = (
                f"Neurologic/cognitive group: cellulose z > 0 in {dec_p}% of decreased and "
                f"{inc_p}% of increased fecal profiles (grouped-screen q = {qv:.5f})."
            )

    elif fig_path.endswith("pd_ecological_reporter_combined.png"):
        c_row = metrics.get("pd_eco_cellulose")
        r_row = metrics.get("pd_eco_composite")
        if c_row and r_row:
            c_val = float(c_row["study_weighted_mean_contrast"])
            r_val = float(r_row["study_weighted_mean_contrast"])
            n_st = c_row["paired_study_count"]
            takeaway = (
                "PD's increased-minus-decreased paired contrast is "
                f"{c_val:+.3f} for cellulose and {r_val:+.3f} for the reporter composite across "
                f"{n_st} paired studies; other conditions share this pattern."
            )

    elif fig_path.endswith("pd_condition_taxon_trait_network.png"):
        w_row = metrics.get("pd_halo_wilcoxon")
        if w_row:
            mean_c = float(w_row["mean_contrast"])
            q_val = float(w_row["q_value"])
            p_val = float(w_row["p_value"])
            n_st = w_row["paired_study_count"]
            if q_val <= 0.05:
                title = "Study-level testing links PD fecal signatures to halophily"
                takeaway = (
                    f"PD moderately halophilic scores rise by {mean_c:.3f} across {n_st} paired "
                    f"studies (BH q = {q_val:.5f}); {taxon_family_text(metrics)}."
                )
            else:
                leads = "leads" if metrics.get("pd_net_lead") is w_row else "in"
                title = (
                    f"PD halophily {leads} the {metrics['pd_net_total']}-test study-level fecal "
                    f"screen (BH q = {q_val:.3g})"
                )
                takeaway = (
                    f"PD moderately halophilic scores change by {mean_c:+.3f} across {n_st} "
                    f"paired studies ({w_row['positive_study_count']} up vs "
                    f"{w_row['negative_study_count']} down, p = {p_val:.2e}, BH q = "
                    f"{q_val:.4f}); {taxon_family_text(metrics)}."
                )

    elif fig_path.endswith("pd_exploratory_taxa_six_traits_influence_heatmap.png"):
        caveat += (
            " The literature reference includes studies that are also in the analysed corpus, "
            "such as Wallen 2020 and 2022."
        )

    elif fig_path == "data/derived/enrichment_infographic.png":
        title = "66,240 ORA enrichment calls describe individual microbial profiles"
        takeaway = (
            "66,240 profile \N{MULTIPLICATION SIGN} trait enrichment calls pass BH q ≤ 0.001 and "
            "fold ≥ 1.5 across "
            "1,480,390 testable pairs (4,670 human profiles \N{MULTIPLICATION SIGN} 317 traits; "
            "35,977 inc vs 30,263 dec)."
        )
        caveat = (
            "Regenerated on current taxon-level hypergeometric ORA model (mean 9.63 hits/profile; "
            "top density: Dental plaque at 224.06 calls/100 samples). Profile ORA is not a paired "
            "condition test."
        )

    elif fig_path == "data/derived/depletion_infographic.png":
        title = "23 ORA depletion calls pass strict fold ≤ 0.5 and BH q ≤ 0.001"
        takeaway = (
            "Under the taxon-level hypergeometric model, only 23 profile \N{MULTIPLICATION SIGN} "
            "trait pairs pass BH q ≤ 0.001 "
            "and fold ≤ 0.5 across 1,480,390 tests (1 increased vs 22 decreased; 13 in Feces, 9 in "
            "Stomach)."
        )
        caveat = (
            "Replaces the legacy assertion-level 100,999-call figure; finite signature sizes limit "
            "one-sided lower-tail power for low-prevalence traits."
        )

    elif fig_path == "data/derived/trait_enrichment_clustered_heatmap.png":
        title = "Trait enrichment varies across 850 profiles and 278 traits"
        takeaway = (
            "Clustered ORA significance matrix (-log10 q, capped at 10) across 850 filtered human "
            "profiles "
            "and 278 traits (min trait prevalence = 14, min profile hits = 12)."
        )
        caveat = (
            "Rows are individual BugSigDB profiles rather than independent cohorts; "
            "the atlas summarizes profile-level ORA structure without testing paired condition "
            "contrasts."
        )

    elif fig_path == "data/derived/Sankey_Diagrams/facultative-aerobe_sankey_inc.png":
        title = "Facultative-aerobe ORA enrichment recurs in increased signatures"
        takeaway = (
            "Flows trace condition → body site → facultatively aerobic ORA enrichment (BH q ≤ "
            "0.001, "
            "fold ≥ 1.5) among increased taxon signatures on the rerun corpus."
        )
        caveat = (
            "Width summarizes enriched profile counts and color summarizes mean -log10(q); "
            "flows are descriptive summaries of profile ORA hits."
        )

    elif fig_path == "data/derived/Sankey_Diagrams/anaerobe_sankey_dec.png":
        title = "Anaerobe ORA enrichment recurs within decreased signatures"
        takeaway = (
            "Flows trace condition → body site → anaerobic/obligately anaerobic ORA enrichment "
            "(BH q ≤ 0.001, fold ≥ 1.5) within decreased taxon signatures."
        )
        caveat = (
            "Enrichment of anaerobe annotations inside a decreased signature means decreased taxa "
            "are "
            "disproportionately anaerobic relative to the background database."
        )

    elif fig_path.endswith("table1_traits_direction_shift.png"):
        title = (
            "Six Table 1 traits show consistent pooled and disease-oriented study direction shifts"
        )
        takeaway = (
            "Across 959 human studies (post-drop_zero_total_pairs), the facultatively aerobic "
            "share is larger in the increased signature in 532/915 non-tied studies (58.1%, BH "
            "q = 5.66e-6) and the cellulose catabolic process share is larger in the decreased "
            "signature in 441/774 (43.0%↑ / 57.0%↓, BH q = 2.35e-4); in 569 disease-vs-healthy "
            "oriented "
            "studies, "
            "the cellulose shift toward decreased signatures strengthens to 283/465 down (60.8%↓, "
            "p = 4.21e-6, q_317 = "
            "2.09e-4)."
        )
        caveat = (
            "Pooled study-level binomial sign tests (BH over the 6 Table 1 traits; q_317 over all "
            "317 traits); the permutation screened_fraction_sign_test = 0.091 reflects "
            "post-selection conditioning on the 6 pre-screened traits."
        )

    return title, takeaway, caveat


# The six Table 1 traits, as the paired-contrast tables name them.
TABLE1_TRAIT_IDS = (
    "METPO:1000608",
    "METPO:2000103|GO:0030245",
    "METPO:1000668",
    "METPO:2000103|pathways:arsenate_reduction",
    "METPO:1000683",
    "METPO:1000623",
)
CONDITION_SHORT = {
    "Parkinson's disease": "PD",
    "HIV infection": "HIV",
    "Colorectal cancer": "CRC",
}
TRAIT_SHORT = {"cellulose catabolic process": "cellulose", "moderately halophilic": "halophily"}


def pd_direction_caveat(model: dict[str, str], positive: dict[str, str]) -> str:
    """The study-clustered test of a PD-specific cellulose direction, from its model rows.

    ``model`` is the cellulose z-score outcome; ``positive`` the share of profiles with a
    positive z-score, which the figure's right panel plots (#590, #608).
    """
    rows = (model, positive)
    shown = [float(r["p_value"]) < 0.05 and float(r["coefficient"]) > 0 for r in rows]
    against = [float(r["p_value"]) < 0.05 and float(r["coefficient"]) < 0 for r in rows]
    if any(against):
        verdict = "at least one measure points the other way, against the figure's title"
    elif all(shown):
        verdict = "a PD-specific direction is shown at p < 0.05"
    elif any(shown):
        verdict = "a PD-specific direction is shown at p < 0.05 for one measure only"
    else:
        verdict = "a PD-specific direction is not shown at p < 0.05"
    return (
        "Exploratory PD diagnostic. Its profile-level Fisher test treats the profiles of one "
        "study as independent. In the study-clustered model "
        "(pd_cellulose_interaction_model.tsv), the extra PD effect in decreased feces is "
        f"{signed(float(model['coefficient']), 2)} cellulose z (95% CI "
        f"{signed(float(model['ci_low']), 2)} to {signed(float(model['ci_high']), 2)}, p = "
        f"{p_text(float(model['p_value']))}), and "
        f"{signed(100 * float(positive['coefficient']))} percentage points (95% CI "
        f"{signed(100 * float(positive['ci_low']))} to {signed(100 * float(positive['ci_high']))}"
        f", p = {p_text(float(positive['p_value']))}) in the share of profiles with cellulose "
        f"z > 0 that the right panel plots, so {verdict}."
    )


def pd_profiles_text(metrics: dict) -> str:
    """'130 increased vs 123 decreased fecal profiles', the lollipop's population (#641)."""
    summary = metrics["pd_six_trait_summary"]

    def counts(column: str) -> str:
        values = sorted({int(r[column]) for r in summary})
        return str(values[0]) if len(values) == 1 else f"{values[0]}\N{EN DASH}{values[-1]}"

    return (
        f"{counts('increased_profiles')} increased vs {counts('decreased_profiles')} decreased "
        "fecal profiles"
    )


def survivor_text(metrics: dict) -> str:
    """The profile-level Fisher FDR survivors with their direction, read from the table."""
    rows = metrics["ind_sig_rows"]
    listed = ", ".join(
        f"{CONDITION_SHORT.get(r['condition'], r['condition'])}"
        f"{'' if r['body_site_group'] == 'Feces' else ' ' + r['body_site_group']} "
        f"{TRAIT_SHORT.get(r['trait_label'], r['trait_label'])} "
        f"{arrow(float(r['increased_minus_decreased_positive_fraction']))}"
        for r in rows
    )
    noun = "survivor" if len(rows) == 1 else "survivors"
    return f"the {len(rows)} profile Fisher FDR {noun} (★ q ≤ 0.05: {listed})"


def get_new_slides_metadata(metrics: dict) -> list[dict[str, str]]:
    """Metadata for the Methods plates, the Tufte plates and the newly generated figures."""
    new_slides = [
        # --- Section 0A: Two High-Level Visual Methods Architecture Plates ---
        {
            "owner": "Marcin / Noel / Ilyas / Jed",
            "title": "Methods Plate I: BugSigDB \N{MULTIPLICATION SIGN} KG-Microbe data "
            "architecture & taxonomic coverage funnel",
            "tag": "Methods architecture · Plate I",
            "figure_file_path": "data/derived/tufte_plates/"
            "methods_plate1_data_architecture_funnel.png",
            "code_path": "scripts/generate_tufte_plates.py",
            "takeaway": "High-level visual overview of the data ingestion funnel: 1,039 BugSigDB "
            "v1.3.1 human studies (2,335 paired experiments, 4,670 profiles, 48,370 "
            "taxon occurrences) joined to 5,829 KG-Microbe organisms "
            "\N{MULTIPLICATION SIGN} 317 normalized traits via NCBI Taxonomy subtree "
            "rollup.",
            "caveat": "Highlights the rank-stratified coverage cliff (Genus/Family 96.5% vs "
            "Species 78.0% occurrences / 54.1% unique species) and defines the OR-rule "
            "oxygen/Gram guilds and 6 Table 1 focal traits.",
            "methods_path": "methods/README.md",
            "evidence_path": "trait-visualization/trait-profile-investigation/taxa-coverage/"
            "taxa_coverage_by_rank.tsv",
            "deck_section": "Main",
            "analysis": "methods_plate_1",
            "selection_reason": "High-level visual Methods diagram replacing dense text-heavy "
            "methodological tables.",
        },
        {
            "owner": "Marcin / Noel / Ilyas / Jed",
            "title": "Methods Plate II: Three-arm statistical architecture & multi-scale "
            "convergence funnel",
            "tag": "Methods architecture · Plate II",
            "figure_file_path": "data/derived/tufte_plates/"
            "methods_plate2_three_arm_architecture.png",
            "code_path": "scripts/generate_tufte_plates.py",
            "takeaway": "Compares the statistical unit, test formulation, and multiplicity control "
            "across Arm 1 (continuous Δz paired contrasts & positive-z Fisher screen), "
            "Arm 2 (taxon-calibrated hypergeometric ORA: 66,240 enrichments vs 23 "
            "depletions), and Arm 3 (study binomial sign tests, McNemar O₂-guild test, "
            "and 15-pair compositional log-ratios across all 959 studies and 569 "
            "disease-vs-healthy oriented studies).",
            "caveat": "Bottom convergence bar shows how 317 KG-Microbe traits narrow to "
            f"{metrics['plate_counts'].get('three_arm_concordant', '?')} three-arm concordant "
            f"study-BH traits ({metrics['plate_counts'].get('three_arm_strict_by', '?')} at "
            "strict BY q ≤ 0.05), explains "
            "why screened_fraction_sign_test = 0.091 is a post-selection artifact, and "
            "highlights 4 BY-significant log-ratios in the 569 oriented "
            "disease-vs-healthy studies.",
            "methods_path": "methods/three_arm_statistics.md",
            "evidence_path": "data/derived/trait_screen/trait_direction_screen.tsv",
            "deck_section": "Main",
            "analysis": "methods_plate_2",
            "selection_reason": "High-level visual Methods architecture diagram explaining how the "
            "three arms complement and validate one another.",
        },
        # --- Section 0B: Seven Simplified High-Density Tufte Analytical Plates ---
        {
            "owner": "Marcin / Noel / Ilyas",
            "title": "Tufte Plate I: Sample partitioning funnel & multi-scale FDR slopegraph",
            "tag": "Tufte synthesis · Plate I",
            "figure_file_path": "data/derived/tufte_plates/"
            "tufte_plate1_scale_attenuation_slopegraph.png",
            "code_path": "scripts/generate_tufte_plates.py",
            "takeaway": plate_one_text(metrics)[0],
            "caveat": plate_one_text(metrics)[1],
            "methods_path": "methods/table1_trait_condition_body_site_maps.md",
            "evidence_path": "trait-visualization/trait-profile-investigation/"
            "trait_profile_paired_contrast_trait_stats_pooled.tsv",
            "deck_section": "Main",
            "analysis": "tufte_plate_1",
            "selection_reason": "Simplified visual Tufte slopegraph resolving why global vs "
            "stratified tests differ in statistical power.",
        },
        {
            "owner": "Marcin / Noel",
            "title": paired_contrast_title(metrics),
            "tag": "Tufte synthesis · Plate II",
            "figure_file_path": "data/derived/tufte_plates/"
            "tufte_plate2_317_trait_guild_dichotomy.png",
            "code_path": "scripts/generate_tufte_plates.py",
            "takeaway": paired_contrast_text(metrics),
            "caveat": paired_contrast_caveat(metrics),
            "methods_path": "trait-visualization/trait-profile-investigation/"
            "paired_contrast_interpretation_notes.md",
            "evidence_path": "trait-visualization/trait-profile-investigation/"
            "trait_profile_paired_contrast_trait_stats_feces.tsv",
            "deck_section": "Main",
            "analysis": "tufte_plate_2",
            "selection_reason": "Tufte small-multiples forest plot expanding beyond Table 1 to all "
            "22 study-significant traits across 4 anatomical strata.",
        },
        {
            "owner": "Noel / Marcin",
            "title": "Tufte Plate III: Compositional closure & reference-explicit paired "
            "log-ratios",
            "tag": "Tufte synthesis · Plate III",
            "figure_file_path": "data/derived/tufte_plates/"
            "tufte_plate3_compositional_log_ratio_matrix.png",
            "code_path": "scripts/generate_tufte_plates.py",
            "takeaway": "Visualizes all 15 pairwise within-study trait log-ratios across all 959 "
            "studies (●) and 569 disease-vs-healthy oriented studies (◆): 3 pairs pass "
            "15-test BY q ≤ 0.05 across all 959 studies (FA/CEL, CEL/ARS, CEL/SPH) and "
            "4 pass in the 569 oriented studies, adding CEL/HAL (91↑:141↓ [39.2%], BY "
            "q = 0.0155★).",
            "caveat": "None of the 10 pairs among the 5 increased-side traits passes BY "
            "q ≤ 0.05 in either cohort. That is a failure to detect a difference among them, "
            "not evidence that they co-move.",
            "methods_path": "methods/table1_paired_log_ratios.md",
            "evidence_path": "data/derived/table1_composition/trait_log_ratio_stats.tsv",
            "deck_section": "Main",
            "analysis": "tufte_plate_3",
            "selection_reason": "Clean 2-panel visual representation of the reference-explicit "
            "compositional log-ratio matrix and forest plot (PR #181) with "
            "569-study oriented sensitivity overlay.",
        },
        {
            "owner": "Noel / Marcin / Ilyas",
            "title": "Tufte Plate IV: Fecal cohort \N{MULTIPLICATION SIGN} trait paired-shift "
            "heatmap & lead disease forest plot",
            "tag": "Tufte synthesis · Plate IV",
            "figure_file_path": "data/derived/tufte_plates/"
            "tufte_plate4_diverging_condition_trait_triangulation.png",
            "code_path": "scripts/generate_tufte_plates.py",
            "takeaway": "Pairs a clean 14-condition \N{MULTIPLICATION SIGN} 7-trait study-weighted "
            "Δz heatmap (left) with a focused forest plot of the 12 lead Fecal Disease "
            "\N{MULTIPLICATION SIGN} Trait contrasts in PD, CRC, HIV, and COVID-19 "
            "(right).",
            "caveat": f"Highlights {survivor_text(metrics)}.{fragile_sentence(metrics)} The "
            "cellulose deep-dive slides test the CRC stool-versus-mucosa difference and the "
            "overlap of decreased taxa across PD, HIV and CRC.",
            "methods_path": "methods/table1_trait_condition_body_site_maps.md",
            "evidence_path": "trait-visualization/trait-profile-investigation/PD-cellulose/"
            "pd_condition_taxon_trait_network_condition_trait_tests.tsv",
            "deck_section": "Main",
            "analysis": "tufte_plate_4",
            "selection_reason": "Simplified visual triangulation across 14 major fecal cohorts.",
        },
        {
            "owner": "Marcin / Noel",
            "title": "Tufte Plate V: Oxygen-switch scores by body site & cross-disease genus "
            "switchers",
            "tag": "Tufte synthesis · Plate V",
            "figure_file_path": "data/derived/tufte_plates/"
            "tufte_plate5_anatomical_oxygen_polarity_and_switchers.png",
            "code_path": "scripts/generate_tufte_plates.py",
            "takeaway": oxygen_switch_text(metrics),
            "caveat": crc_oxygen_text(metrics["oxygen_switch"])
            + "The cellulose deep-dive slides compare CRC stool with non-fecal GI samples, "
            "pooling the colon and intestinal mucosa sites that the oxygen-switch table splits "
            f"into Lower GI tract and Intestinal / GI mucosa ({crc_mucosa_shown(metrics)}); "
            "most of those mucosal contrasts are not patients against controls.",
            "methods_path": "trait-visualization/trait-profile-investigation/taxa-oxygen-switches/"
            "README.md",
            "evidence_path": "trait-visualization/trait-profile-investigation/taxa-oxygen-switches/"
            "condition_oxygen_switch_scores.tsv",
            "deck_section": "Main",
            "analysis": "tufte_plate_5",
            "selection_reason": "Tufte synthesis plate of the lumenal-vs-mucosal/oral oxygen "
            "polarity inversion and genus switchers.",
        },
        {
            "owner": "Marcin / Ilyas / Noel",
            "title": "Tufte Plate VI: Taxonomic resolution cliff & zero-trait root-cause audit",
            "tag": "Tufte synthesis · Plate VI",
            "figure_file_path": "data/derived/tufte_plates/"
            "tufte_plate6_taxonomic_coverage_blindspots.png",
            "code_path": "scripts/generate_tufte_plates.py",
            "takeaway": "Genus (96.5%) and Family (96.5%) occurrences in BugSigDB have high "
            "KG-Microbe coverage, whereas 45.9% of distinct species (1,764/3,846) and "
            "22.0% of species occurrences (4,205/19,113) have zero traits.",
            "caveat": f"Auditing the top {metrics['zero_trait_taxa']} zero-trait taxa "
            f"({metrics['zero_trait_occurrences']:,} occurrences) reveals actionable TaxID "
            "synonyms: Oscillospiraceae (NCBITaxon:541000, 82 occs, 0 traits vs "
            "canonical 216572 with 162 traits), Rothia (508215 vs 32207), and "
            "Escherichia/Shigella slash bins (139 occs).",
            "methods_path": "trait-visualization/trait-profile-investigation/taxa-coverage/"
            "README.md",
            "evidence_path": "trait-visualization/trait-profile-investigation/taxa-coverage/"
            "top_zero_trait_bugsigdb_taxa.tsv",
            "deck_section": "Main",
            "analysis": "tufte_plate_6",
            "selection_reason": "Tufte data-quality plate quantifying the species coverage cliff "
            "and high-leverage TaxID synonym fixes.",
        },
        {
            "owner": "Marcin / Noel / Ilyas",
            "title": "Tufte Plate VII: Three-arm concordance scatter & 19-trait study-FDR forest "
            "plot",
            "tag": "Tufte synthesis · Plate VII",
            "figure_file_path": "data/derived/tufte_plates/"
            "tufte_plate7_three_arm_317_trait_concordance.png",
            "code_path": "scripts/generate_tufte_plates.py",
            "takeaway": "All 6 original Table 1 traits (★) pass all 3 arms with concordant "
            "direction; screening all 317 traits reveals 11 additional "
            "3-arm-concordant traits (17 total at study BH q ≤ 0.05; 8 at strict BY q "
            "≤ 0.05) plus 2 paired-down / ORA-inverted guilds (GC low & Anaerobic).",
            "caveat": "Panel B overlays the 569 disease-vs-healthy oriented subset (○) on the "
            "959-study sign test (●), showing Cellulose catabolism strengthening from "
            "43.0%↑ (57.0%↓, q=2.35e-4) to 39.2%↑ (60.8%↓, p=4.21e-6, q_317=2.09e-4, "
            "BY_317=0.0013).",
            "methods_path": "methods/three_arm_statistics.md",
            "evidence_path": "data/derived/trait_screen/trait_direction_screen.tsv",
            "deck_section": "Main",
            "analysis": "tufte_plate_7",
            "selection_reason": "Visual forest-and-scatter comparison of the 6 Table 1 traits and "
            "all 317 traits across all three analytical arms and oriented "
            "disease-vs-healthy subset.",
        },
        # --- Section A: Unique Publication Deck Figure (PC1 Annotation Density Confounding) ---
        {
            "owner": "Marcin / Ilyas",
            "title": "Annotation density confounds PC1 of raw profile trait scores",
            "tag": "Deck synthesis · PC1 annotation-density diagnostic",
            "figure_file_path": "data/derived/deck_figures/"
            "figure4_annotation_density_pc1_landscape.png",
            "code_path": "scripts/generate_deck_visualizations.py",
            "takeaway": "Demonstrates how per-taxon annotation density in KG-Microbe drives "
            "primary variance across uncalibrated trait profiles.",
            "caveat": "Motivates paired within-study contrasts (Arm 1/3) and background-calibrated "
            "taxon ORA (Arm 2).",
            "methods_path": "methods/README.md",
            "evidence_path": "data/derived/study_trait_profiles_wide.tsv",
            "deck_section": "Main",
            "analysis": "deck_figure_4",
            "selection_reason": "Uniquely visualizes PC1 annotation-density confounding (other "
            "deck_figures pruned as redundant with core plates).",
        },
        # --- Section B: Distinct Noel Figures (Table 1 Feces vs Other & PR148 Gram Sankeys) ---
        {
            "owner": "Noel",
            "title": "Feces vs pooled Other body-site focus for the six Table 1 traits",
            "tag": "Noel · Table 1 anatomical focus",
            "figure_file_path": "data/derived/table1_traits/table1_traits_focus_feces_other.png",
            "code_path": "src/bsdbng/table1_trait_maps.py",
            "takeaway": "Side-by-side comparison of Feces and pooled Other body sites across the "
            "51-condition Table 1 atlas.",
            "caveat": "Solid cell outlines indicate profile-level ORA support; no individual "
            "condition \N{MULTIPLICATION SIGN} site cell passes the 18,152-test sign "
            "screen.",
            "methods_path": "methods/README.md",
            "evidence_path": "data/derived/table1_traits/"
            "table1_traits_condition_body_site_matrix.tsv",
            "deck_section": "Main",
            "analysis": "table1_traits_focus_feces_other",
            "selection_reason": "Retained concise Feces-vs-Other comparison alongside the full "
            "Table 1 condition \N{MULTIPLICATION SIGN} body-site heatmap.",
        },
        {
            "owner": "Noel",
            "title": "Gram-positive trait support across increased microbial signatures",
            "tag": "Noel · Gram Sankey · increased",
            "figure_file_path": "data/derived/noel_analysis_rerun_2026-09-25/pr148/gram_sankey/"
            "gram_positive_inc_sankey.png",
            "code_path": "src/bsdbng/gram_sankey.py",
            "takeaway": "Traces condition → body site → Gram-positive trait support in increased "
            "signatures on the current wide trait matrix.",
            "caveat": "Descriptive flow summary of positive trait shares; does not test paired "
            "within-study sign consistency.",
            "methods_path": "",
            "evidence_path": "data/derived/study_trait_profiles_wide.tsv",
            "deck_section": "Main",
            "analysis": "gram_positive_inc_sankey",
            "selection_reason": "PR #148 Gram-stain Sankey suite regenerated on current wide "
            "matrix.",
        },
        {
            "owner": "Noel",
            "title": "Gram-positive trait support across decreased microbial signatures",
            "tag": "Noel · Gram Sankey · decreased",
            "figure_file_path": "data/derived/noel_analysis_rerun_2026-09-25/pr148/gram_sankey/"
            "gram_positive_dec_sankey.png",
            "code_path": "src/bsdbng/gram_sankey.py",
            "takeaway": "Traces condition → body site → Gram-positive trait support in decreased "
            "signatures on the current wide trait matrix.",
            "caveat": "Descriptive flow summary of positive trait shares; compare alongside "
            "Gram-positive increased flows.",
            "methods_path": "",
            "evidence_path": "data/derived/study_trait_profiles_wide.tsv",
            "deck_section": "Main",
            "analysis": "gram_positive_dec_sankey",
            "selection_reason": "PR #148 Gram-stain Sankey suite regenerated on current wide "
            "matrix.",
        },
        {
            "owner": "Noel",
            "title": "Gram-negative trait support across increased microbial signatures",
            "tag": "Noel · Gram Sankey · increased",
            "figure_file_path": "data/derived/noel_analysis_rerun_2026-09-25/pr148/gram_sankey/"
            "gram_negative_inc_sankey.png",
            "code_path": "src/bsdbng/gram_sankey.py",
            "takeaway": "Traces condition → body site → Gram-negative trait support in increased "
            "signatures on the current wide trait matrix.",
            "caveat": "Descriptive flow summary of positive trait shares; does not isolate "
            "condition-specific FDR effects.",
            "methods_path": "",
            "evidence_path": "data/derived/study_trait_profiles_wide.tsv",
            "deck_section": "Main",
            "analysis": "gram_negative_inc_sankey",
            "selection_reason": "PR #148 Gram-stain Sankey suite regenerated on current wide "
            "matrix.",
        },
        {
            "owner": "Noel",
            "title": "Gram-negative trait support across decreased microbial signatures",
            "tag": "Noel · Gram Sankey · decreased",
            "figure_file_path": "data/derived/noel_analysis_rerun_2026-09-25/pr148/gram_sankey/"
            "gram_negative_dec_sankey.png",
            "code_path": "src/bsdbng/gram_sankey.py",
            "takeaway": "Traces condition → body site → Gram-negative trait support in decreased "
            "signatures on the current wide trait matrix.",
            "caveat": "Descriptive flow summary of positive trait shares; compare alongside "
            "Gram-negative increased flows.",
            "methods_path": "",
            "evidence_path": "data/derived/study_trait_profiles_wide.tsv",
            "deck_section": "Main",
            "analysis": "gram_negative_dec_sankey",
            "selection_reason": "PR #148 Gram-stain Sankey suite regenerated on current wide "
            "matrix.",
        },
        # --- Section C: Distinct Ilyas Figures (Combined Guild Sankey, Six-Trait ORA Sankeys,
        # Cytoscape Networks & Category Distance) ---
        {
            "owner": "Ilyas",
            "title": "Combined facultative-aerobe (inc) and anaerobe (dec) ORA Sankey",
            "tag": "Ilyas · Combined guild ORA Sankey",
            "figure_file_path": "data/derived/Sankey_Diagrams/combined_aerobe_anaerobe_sankey.png",
            "code_path": "scripts/generate_guild_sankey.py",
            "takeaway": "Jointly visualizes facultative-aerobe ORA enrichment in increased "
            "signatures and anaerobe ORA enrichment in decreased signatures.",
            "caveat": "Consolidates the single-guild aerobe-inc and anaerobe-dec ORA Sankeys into "
            "one unified diagram (BH q ≤ 0.001, fold ≥ 1.5).",
            "methods_path": "",
            "evidence_path": "data/derived/trait_enrichment_results.tsv",
            "deck_section": "Main",
            "analysis": "combined_aerobe_anaerobe_sankey",
            "selection_reason": "Consolidates facultative-aerobe (inc) and anaerobe (dec) ORA "
            "Sankeys into one slide.",
        },
        {
            "owner": "Ilyas",
            "title": "Six-trait ORA Sankey across increased microbial signatures",
            "tag": "Ilyas · Six-trait ORA Sankey · increased",
            "figure_file_path": "data/derived/Sankey_Diagrams/six_trait_sankey_increased.png",
            "code_path": "scripts/generate_six_trait_sankey.py",
            "takeaway": "Maps condition → body site → six-shared-trait ORA enrichments within "
            "increased signatures on the rerun corpus.",
            "caveat": "Uses Arm 2 profile ORA calls (BH q ≤ 0.001, fold ≥ 1.5) rather than Arm 1 "
            "positive-z profile shares.",
            "methods_path": "",
            "evidence_path": "data/derived/trait_enrichment_results.tsv",
            "deck_section": "Main",
            "analysis": "six_trait_sankey_increased",
            "selection_reason": "Six-trait ORA Sankey from Ilyas's visualization suite.",
        },
        {
            "owner": "Ilyas",
            "title": "Six-trait ORA Sankey across decreased microbial signatures",
            "tag": "Ilyas · Six-trait ORA Sankey · decreased",
            "figure_file_path": "data/derived/Sankey_Diagrams/six_trait_sankey_decreased.png",
            "code_path": "scripts/generate_six_trait_sankey.py",
            "takeaway": "Maps condition → body site → six-shared-trait ORA enrichments within "
            "decreased signatures on the rerun corpus.",
            "caveat": "Uses Arm 2 profile ORA calls (BH q ≤ 0.001, fold ≥ 1.5) rather than Arm 1 "
            "positive-z profile shares.",
            "methods_path": "",
            "evidence_path": "data/derived/trait_enrichment_results.tsv",
            "deck_section": "Main",
            "analysis": "six_trait_sankey_decreased",
            "selection_reason": "Six-trait ORA Sankey from Ilyas's visualization suite.",
        },
        {
            "owner": "Ilyas",
            "title": "Multipartite ORA network for six shared traits (increased signatures)",
            "tag": "Ilyas · Six-trait ORA network · increased",
            "figure_file_path": "data/derived/Cytoscape_Diagrams/"
            "cytoscape_six_traits_increased.png",
            "code_path": "scripts/generate_cytoscape_diagram.py",
            "takeaway": "Multipartite graph linking conditions, body sites, and the six shared "
            "traits across increased-signature ORA hits.",
            "caveat": "Descriptive network topology of ORA co-enrichment; edge weights reflect "
            "profile hit counts.",
            "methods_path": "",
            "evidence_path": "data/derived/trait_enrichment_results.tsv",
            "deck_section": "Main",
            "analysis": "cytoscape_six_traits_increased",
            "selection_reason": "Cytoscape network from Ilyas's ORA visualization suite.",
        },
        {
            "owner": "Ilyas",
            "title": "Multipartite ORA network for six shared traits (decreased signatures)",
            "tag": "Ilyas · Six-trait ORA network · decreased",
            "figure_file_path": "data/derived/Cytoscape_Diagrams/"
            "cytoscape_six_traits_decreased.png",
            "code_path": "scripts/generate_cytoscape_diagram.py",
            "takeaway": "Multipartite graph linking conditions, body sites, and the six shared "
            "traits across decreased-signature ORA hits.",
            "caveat": "Descriptive network topology of ORA co-enrichment; edge weights reflect "
            "profile hit counts.",
            "methods_path": "",
            "evidence_path": "data/derived/trait_enrichment_results.tsv",
            "deck_section": "Main",
            "analysis": "cytoscape_six_traits_decreased",
            "selection_reason": "Cytoscape network from Ilyas's ORA visualization suite.",
        },
        {
            "owner": "Ilyas",
            "title": "Trait-category directional distance heatmap across ORA calls",
            "tag": "Ilyas · ORA category directional distance",
            "figure_file_path": "data/derived/"
            "trait_enrichment_category_direction_distance_heatmap.png",
            "code_path": "scripts/plot_enrichment_histograms.py",
            "takeaway": "Quantifies directional separation between increased and decreased ORA "
            "enrichment profiles across trait categories.",
            "caveat": "Aggregates profile-level ORA hits across conditions; reflects corpus-wide "
            "directional asymmetry.",
            "methods_path": "",
            "evidence_path": "data/derived/trait_enrichment_results.tsv",
            "deck_section": "Main",
            "analysis": "trait_enrichment_category_direction_distance_heatmap",
            "selection_reason": "ORA category directional distance figure from "
            "plot_enrichment_histograms.py.",
        },
        # --- Section D: Distinct Jed Main Figures (Global Heatmaps, PD Cellulose, Oral Atlas) ---
        {
            "owner": "Jed",
            "title": "Global condition \N{MULTIPLICATION SIGN} trait module clustered heatmap (317 "
            "traits)",
            "tag": "Jed · Global trait-module atlas",
            "figure_file_path": "trait-visualization/condition_trait_module_clustered_heatmap.png",
            "code_path": "scripts/create_trait_Heatmaps.py",
            "takeaway": "Clusters human conditions against co-varying microbial trait modules "
            "across the 4,670-profile wide matrix.",
            "caveat": "Exploratory hierarchical clustering of condition-level trait profiles; "
            "module labels are data-derived.",
            "methods_path": "trait-visualization/trait-profile-investigation/README.md",
            "evidence_path": "trait-visualization/condition_trait_module_membership.tsv",
            "deck_section": "Main",
            "analysis": "condition_trait_module_clustered_heatmap",
            "selection_reason": "Core global heatmap from just trait-heatmaps.",
        },
        {
            "owner": "Jed",
            "title": "Study-level paired trait-profile clustered heatmap",
            "tag": "Jed · Study-level trait atlas",
            "figure_file_path": "trait-visualization/study_trait_profile_clustered_heatmap.png",
            "code_path": "scripts/create_trait_Heatmaps.py",
            "takeaway": "Displays study-level trait profile structure across BugSigDB experiments "
            "joined to KG-Microbe traits.",
            "caveat": "Unadjusted z-score visualization across studies; use paired sign/Wilcoxon "
            "screens for formal inference.",
            "methods_path": "trait-visualization/trait-profile-investigation/README.md",
            "evidence_path": "data/derived/study_trait_profiles_wide.tsv",
            "deck_section": "Main",
            "analysis": "study_trait_profile_clustered_heatmap",
            "selection_reason": "Core study-level heatmap from just trait-heatmaps.",
        },
        {
            "owner": "Jed",
            "title": "Study-level paired trait enrichment contrast heatmap",
            "tag": "Jed · Paired enrichment contrast atlas",
            "figure_file_path": "trait-visualization/"
            "study_trait_profile_paired_enrichment_heatmap.png",
            "code_path": "scripts/create_trait_Heatmaps.py",
            "takeaway": "Summarizes paired increased-minus-decreased trait shifts across studies "
            "and trait modules.",
            "caveat": "Exploratory clustering view; see Table 1 and condition "
            "\N{MULTIPLICATION SIGN} site screens for multiplicity-adjusted tests.",
            "methods_path": "trait-visualization/trait-profile-investigation/README.md",
            "evidence_path": "data/derived/study_trait_profiles_wide.tsv",
            "deck_section": "Main",
            "analysis": "study_trait_profile_paired_enrichment_heatmap",
            "selection_reason": "Core paired-enrichment heatmap from just trait-heatmaps.",
        },
        {
            "owner": "Jed",
            "title": "Data-driven taxa \N{MULTIPLICATION SIGN} guild discovery clustered heatmap",
            "tag": "Jed · Taxa-guild discovery",
            "figure_file_path": "trait-visualization/taxa_guild_discovery_clustered_heatmap.png",
            "code_path": "scripts/create_trait_Heatmaps.py",
            "takeaway": "Clusters recurring BugSigDB taxa by their shared KG-Microbe trait "
            "annotations to define empirical functional guilds.",
            "caveat": "Guild assignments reflect database annotation co-occurrence rather than in "
            "situ transcriptomic activity.",
            "methods_path": "trait-visualization/trait-profile-investigation/README.md",
            "evidence_path": "trait-visualization/taxa_guild_From_data_membership.tsv",
            "deck_section": "Main",
            "analysis": "taxa_guild_discovery_clustered_heatmap",
            "selection_reason": "Core taxa-guild discovery figure from just trait-heatmaps.",
        },
        {
            "owner": "Jed",
            "title": f"PD six-trait direction lollipop: {pd_profiles_text(metrics)}",
            "tag": "Jed · PD six-trait direction summary",
            "figure_file_path": "trait-visualization/trait-profile-investigation/PD-cellulose/"
            "pd_six_trait_direction_lollipop.png",
            "code_path": "scripts/create_pd_cellulose_analysis.py",
            "takeaway": "Mean annotation z-score, increased minus decreased PD fecal profiles: "
            + ", ".join(
                f"{r['trait']} {float(r['increased_minus_decreased_mean_z']):+.2f}"
                for r in metrics["pd_six_trait_summary"]
            )
            + ".",
            "caveat": "Profile-level summary. In the "
            f"{metrics['pd_net_total']}-test study-level Wilcoxon screen, PD halophily has BH "
            f"q = {float(metrics['pd_halo_wilcoxon']['q_value']):.3g} and cellulose q = "
            f"{float(metrics['pd_cellulose_wilcoxon']['q_value']):.3g}.",
            "methods_path": "trait-visualization/trait-profile-investigation/PD-cellulose/"
            "README.md",
            "evidence_path": "trait-visualization/trait-profile-investigation/PD-cellulose/"
            "pd_condition_taxon_trait_network_condition_trait_tests.tsv",
            "deck_section": "Main",
            "analysis": "pd_six_trait_direction_lollipop",
            "selection_reason": "Key PD six-trait summary figure from "
            "create_pd_cellulose_analysis.py.",
        },
        {
            "owner": "Jed",
            "title": "Cellulose-annotated taxa in PD fecal profiles (decreased vs increased)",
            "tag": "Jed · PD cellulose-annotated taxa",
            "figure_file_path": "trait-visualization/trait-profile-investigation/PD-cellulose/"
            "pd_cellulose_taxon_driver_bar.png",
            "code_path": "scripts/create_pd_cellulose_analysis.py",
            "takeaway": "The cellulose-annotated taxa with the most net-decreased PD fecal "
            f"profiles are {listing(metrics['pd_cellulose_drivers'])}.",
            "caveat": "Descriptive contributor breakdown across PD fecal signatures; "
            f"{taxon_family_text(metrics)}. A taxon counts as cellulose-annotated when it or a "
            "taxon below it carries the annotation.",
            "methods_path": "trait-visualization/trait-profile-investigation/PD-cellulose/"
            "README.md",
            "evidence_path": "trait-visualization/trait-profile-investigation/PD-cellulose/"
            "pd_cellulose_taxon_drivers.tsv",
            "deck_section": "Main",
            "analysis": "pd_cellulose_taxon_driver_bar",
            "selection_reason": "Key PD cellulose driver figure from "
            "create_pd_cellulose_analysis.py.",
        },
        {
            "owner": "Jed",
            "title": "PD trait co-carriage network across contributor taxa",
            "tag": "Jed · PD trait co-carriage",
            "figure_file_path": "trait-visualization/trait-profile-investigation/PD-cellulose/"
            "pd_trait_contributors_trait_cocarriage_network.png",
            "code_path": "scripts/create_pd_trait_contributors.py",
            "takeaway": "Quantifies how often the six shared traits are co-annotated on the same "
            "PD contributor taxa in KG-Microbe.",
            "caveat": "Co-carriage explains why facultative aerobiosis, halophily, coccus shape, "
            "and arsenate reduction covary across signatures.",
            "methods_path": "trait-visualization/trait-profile-investigation/PD-cellulose/"
            "README.md",
            "evidence_path": "trait-visualization/trait-profile-investigation/PD-cellulose/"
            "pd_trait_contributors_trait_cocarriage.tsv",
            "deck_section": "Main",
            "analysis": "pd_trait_contributors_trait_cocarriage_network",
            "selection_reason": "Companion co-carriage network to "
            "pd_trait_contributors_species_family_dot_heatmap.png.",
        },
        {
            "owner": "Jed",
            "title": "Oral & dental six-trait signed direction contrast heatmap",
            "tag": "Jed · Oral six-trait atlas",
            "figure_file_path": "trait-visualization/trait-profile-investigation/six-shared-traits/"
            "oral/oral_six_trait_heatmap.png",
            "code_path": "scripts/create_oral_six_trait_heatmap.py",
            "takeaway": "Full overview of signed positive-z direction contrasts across oral and "
            "dental condition \N{MULTIPLICATION SIGN} site strata for the six shared "
            "traits.",
            "caveat": "No oral condition \N{MULTIPLICATION SIGN} site cell passes global FDR "
            "correction; single full figure replaces the paginated oral slices.",
            "methods_path": "trait-visualization/trait-profile-investigation/six-shared-traits/"
            "oral/README.md",
            "evidence_path": "trait-visualization/trait-profile-investigation/six-shared-traits/"
            "oral/oral_six_trait_condition_body_site_contrasts.tsv",
            "deck_section": "Main",
            "analysis": "oral_six_trait_heatmap",
            "selection_reason": "Dedicated oral six-trait heatmap generated by just "
            "oral-six-trait-heatmap.",
        },
    ]
    return new_slides


def get_remaining_appendix_slides(
    already_included: set[str], repo: Path, metrics: dict
) -> list[dict[str, str]]:
    """Collect all remaining distinct regenerated PNGs after pruning redundant/superseded views."""
    # Exclude exact duplicate audit/mirror directories and low-res preview
    excluded_prefixes = (
        "data/derived/noel_analysis_rerun_2026-09-25/pr173/",
        "data/derived/noel_analysis_rerun_2026-09-25/pr148/audit/",
        "data/derived/figure_review_2026-09-21/",
        "methods/figures/",  # identical copies of data/derived/table1_traits/
        PLATES_DIR,  # the plates come from this checkout (#639)
    )
    excluded_exact = {
        "trait-visualization/trait-profile-investigation/six-shared-traits/"
        "six_shared_traits_grouped_condition_body_site_signed_direction_contrast_heatmap_preview.png",
        "trait-visualization/trait-profile-investigation/enrichment-infographic-followup/"
        "my_top3_conditions_relevant_traits_traits_vs_organisms_heatmap.png",
        # Exact SHA-256 byte duplicates
        "data/derived/Sankey_Diagrams/anaerobe_sankey_dec_with_unenriched.png",
        "data/derived/Sankey_Diagrams/facultative-aerobe_sankey_inc_with_unenriched.png",
        # Pruned redundant Main figures (consolidated into composite figures or Tufte plates)
        "data/derived/Sankey_Diagrams/facultative-aerobe_sankey_inc.png",
        "data/derived/Sankey_Diagrams/anaerobe_sankey_dec.png",
        "data/derived/deck_figures/figure1_condition_anatomy_ora_heatmap.png",
        "data/derived/deck_figures/figure2_taxon_calibrated_ora_infographic.png",
        "data/derived/deck_figures/figure3_lead_findings_ora_evidence.png",
        "data/derived/deck_figures/figure5_pd_crc_taxon_attribution_bubble.png",
        "data/derived/table1_traits/table1_traits_condition_heatmap.png",
        "data/derived/table1_traits/table1_traits_focus_top_conditions.png",
        "data/derived/table1_traits/table1_traits_body_site_facets.png",
        "data/derived/table1_traits/table1_traits_first_last_body_site_panels.png",
        "trait-visualization/trait-profile-investigation/taxa-oxygen-switches/"
        "condition_oxygen_switch_scores.png",
        "trait-visualization/trait-profile-investigation/taxa-coverage/"
        "taxa_coverage_histogram_rank_stratified.png",
        "trait-visualization/trait-profile-investigation/six-shared-traits/"
        "individual_conditions_x_body_sites_positive_z_direction_contrasts_page_04.png",
        # Superseded v1/v2 drafts (v3 is retained)
        "trait-visualization/aerobe-enrichment-v2/condition_body_site_guild_map.png",
        "trait-visualization/aerobe-enrichment-v2/condition_body_site_guild_map_v2.png",
        "trait-visualization/all-bugsigdb-condition-body-site-maps/"
        "all_bugsigdb_decreased_profile_condition_body_site_map.png",
        "trait-visualization/all-bugsigdb-condition-body-site-maps/"
        "all_bugsigdb_decreased_profile_condition_body_site_map_v2.png",
        "trait-visualization/all-bugsigdb-condition-body-site-maps/"
        "all_bugsigdb_increased_profile_condition_body_site_map.png",
        "trait-visualization/all-bugsigdb-condition-body-site-maps/"
        "all_bugsigdb_increased_profile_condition_body_site_map_v2.png",
        "trait-visualization/anaerobe-depletion-v2/condition_body_site_guild_map.png",
        "trait-visualization/anaerobe-depletion-v2/condition_body_site_guild_map_v2.png",
        "trait-visualization/trait-profile-investigation/enrichment-infographic-followup/"
        "trait-enrichment-exploration/traits_vs_organisms_heatmap.png",
        # Stretched aspect-ratio duplicates of standard ORA condition heatmaps
        "data/derived/trait_enrichment_x_condition_all_stretched.png",
        "data/derived/trait_enrichment_x_condition_dec_stretched.png",
        "data/derived/trait_enrichment_x_condition_inc_stretched.png",
        # Redundant paginated slices where the full unpaginated heatmap is already in the deck
        "trait-visualization/trait-profile-investigation/six-shared-traits/oral/"
        "oral_six_trait_heatmap_page_01.png",
        "trait-visualization/trait-profile-investigation/six-shared-traits/oral/"
        "oral_six_trait_heatmap_page_02.png",
        "trait-visualization/trait-profile-investigation/six-shared-traits/oral/"
        "oral_six_trait_heatmap_page_03.png",
        "trait-visualization/trait-profile-investigation/six-shared-traits/oral/"
        "oral_six_trait_heatmap_page_04.png",
        "trait-visualization/trait-profile-investigation/six-shared-traits/"
        "six_shared_traits_grouped_condition_body_site_signed_direction_contrast_heatmap_page_01.png",
        "trait-visualization/trait-profile-investigation/six-shared-traits/"
        "six_shared_traits_grouped_condition_body_site_signed_direction_contrast_heatmap_page_02.png",
        "trait-visualization/trait-profile-investigation/six-shared-traits/"
        "six_shared_traits_grouped_condition_body_site_signed_direction_contrast_heatmap_page_03.png",
        "trait-visualization/trait-profile-investigation/six-shared-traits/"
        "six_shared_traits_grouped_condition_body_site_signed_direction_contrast_heatmap_page_04.png",
        "trait-visualization/trait-profile-investigation/six-shared-traits/"
        "six_shared_traits_individual_conditions_condition_body_site_signed_direction_contrast_heatmap_page_02.png",
        "trait-visualization/trait-profile-investigation/six-shared-traits/"
        "six_shared_traits_individual_conditions_condition_body_site_signed_direction_contrast_heatmap_page_03.png",
        "trait-visualization/trait-profile-investigation/six-shared-traits/"
        "six_shared_traits_individual_conditions_condition_body_site_signed_direction_contrast_heatmap_page_05.png",
        "trait-visualization/trait-profile-investigation/six-shared-traits/"
        "six_shared_traits_individual_conditions_condition_body_site_signed_direction_contrast_heatmap_page_06.png",
        # Pure text-table PNGs now represented visually in Methods Plate 1
        "trait-visualization/condition_trait_module_guild_definition_table.png",
        "trait-visualization/anaerobe-depletion-v2/anaerobe_v2_guild_definition_table.png",
        "trait-visualization/trait-profile-investigation/profile-condition-bars/"
        "profile_enrichment_guild_assignment_summary.png",
    }

    all_pngs: list[Path] = []
    for root_sub in ("data/derived", "trait-visualization"):
        root_path = repo / root_sub
        if root_path.exists():
            all_pngs.extend(sorted(root_path.rglob("*.png")))

    appendix_slides: list[dict[str, str]] = []
    for png_path in all_pngs:
        rel = png_path.relative_to(repo).as_posix()
        if rel in already_included or rel in excluded_exact:
            continue
        if any(rel.startswith(p) for p in excluded_prefixes):
            continue
        # The full oral heatmap is a Main slide, so any number of its pages is redundant (#641)
        if re.search(r"/oral_six_trait_heatmap_page_\d+\.png$", rel):
            continue

        # Also skip data/derived/table1_traits/* that were already included via methods/figures/* in
        # v1/v2
        if rel.startswith("data/derived/table1_traits/"):
            equiv = "methods/figures/" + png_path.name
            if equiv in already_included:
                continue

        # Determine owner, code_path, tag, title, takeaway, caveat based on path
        stem = png_path.stem.replace("_", " ").replace("-", " ")
        parent_name = png_path.parent.name

        if rel.startswith("data/derived/Sankey_Diagrams/"):
            owner = "Ilyas"
            code_path = "scripts/generate_guild_sankey.py"
            tag = "Ilyas · Appendix · ORA Sankey"
            title = f"Guild ORA Sankey ({stem})"
            takeaway = (
                "Includes unenriched profile flows alongside ORA-enriched flows "
                f"({png_path.name}) on the rerun corpus."
            )
            caveat = (
                "Descriptive ORA flow summary (BH q ≤ 0.001, fold ≥ 1.5); unenriched bands "
                "show total signature volume."
            )
            evidence_path = "data/derived/trait_enrichment_results.tsv"
        elif rel.startswith("data/derived/"):
            owner = "Ilyas"
            code_path = "scripts/plot_enrichment_histograms.py"
            tag = "Ilyas · Appendix · ORA distribution"
            title = f"ORA enrichment distribution: {stem}"
            takeaway = (
                f"Distribution of taxon-calibrated ORA calls ({png_path.name}) across the "
                "4,670-profile human rerun."
            )
            caveat = (
                "Summarizes marginal ORA hit counts across traits, profiles, or conditions; "
                "not a paired study-level test."
            )
            evidence_path = "data/derived/trait_enrichment_results.tsv"
        elif "/PD-cellulose/" in rel:
            owner = "Jed"
            code_path = "scripts/create_pd_cellulose_analysis.py"
            tag = "Jed · Appendix · PD cellulose investigation"
            title = f"PD cellulose investigation: {stem}"
            takeaway = (
                "Focused Parkinson\N{RIGHT SINGLE QUOTATION MARK}s disease cellulose and "
                f"six-trait diagnostic ({png_path.name}) regenerated on the current wide "
                "matrix."
            )
            caveat = (
                "Exploratory PD diagnostic; evaluate alongside the "
                f"{metrics['pd_net_total']}-test study-level fecal Wilcoxon screen."
            )
            if png_path.name == "pd_cellulose_direction_specificity.png":
                # The figure's own title calls the signal direction-specific (#590).
                caveat = pd_direction_caveat(
                    metrics["pd_interaction"], metrics["pd_interaction_positive"]
                )
            elif png_path.name == "pd_cellulose_taxon_six_trait_heatmap.png":
                caveat = (
                    "Exploratory PD diagnostic. A \N{LEFT DOUBLE QUOTATION MARK}cellulose-"
                    "carrier\N{RIGHT DOUBLE QUOTATION MARK} taxon is one whose KG-Microbe "
                    "annotations, its own or its descendants\N{RIGHT SINGLE QUOTATION MARK}, "
                    "include cellulose catabolic process; cell numbers count the annotated taxa "
                    "at or below each row. It is an annotation, not measured cellulose "
                    "degradation."
                )
            evidence_path = (
                "trait-visualization/trait-profile-investigation/PD-cellulose/"
                "pd_condition_taxon_trait_network_condition_trait_tests.tsv"
            )
        elif "/parkinsons-disease/" in rel:
            owner = "Jed"
            code_path = "scripts/regenerate_trait_profile_investigation.py"
            tag = "Jed · Appendix · Parkinson\N{RIGHT SINGLE QUOTATION MARK}s disease"
            title = f"Parkinson\N{RIGHT SINGLE QUOTATION MARK}s disease trait & taxa view: {stem}"
            takeaway = (
                "Paired fecal trait and target-taxa contrast in "
                f"Parkinson\N{RIGHT SINGLE QUOTATION MARK}s disease ({png_path.name}) from "
                "the rerun corpus."
            )
            caveat = (
                "Descriptive within-condition summary of annotation scores and taxon recurrence."
            )
            evidence_path = "data/derived/study_trait_profiles_wide.tsv"
        elif "/aerobe-enrichment-v2/" in rel:
            owner = "Jed"
            code_path = "scripts/regenerate_trait_profile_investigation.py"
            tag = "Jed · Appendix · Aerobe enrichment v2"
            title = f"Aerobe guild investigation: {stem}"
            takeaway = (
                "Aerobe/facultative-aerobe guild profile and condition-site distribution "
                f"({png_path.name}) on the rerun corpus."
            )
            caveat = (
                "Exploratory guild-level visualization; compare with "
                "Noel\N{RIGHT SINGLE QUOTATION MARK}s study-level McNemar oxygen-guild test."
            )
            evidence_path = "data/derived/study_trait_profiles_wide.tsv"
        elif "/anaerobe-depletion-v2/" in rel:
            owner = "Jed"
            code_path = "scripts/regenerate_trait_profile_investigation.py"
            tag = "Jed · Appendix · Anaerobe depletion v2"
            title = f"Anaerobe guild investigation: {stem}"
            takeaway = (
                "Anaerobe guild decreased-signature distribution and taxon breakdown "
                f"({png_path.name}) on the rerun corpus."
            )
            caveat = (
                "Exploratory guild-level visualization; reflects annotation shares within "
                "decreased signatures."
            )
            evidence_path = "data/derived/study_trait_profiles_wide.tsv"
        elif "/all-bugsigdb-condition-body-site-maps/" in rel:
            owner = "Jed"
            code_path = "scripts/regenerate_trait_profile_investigation.py"
            tag = "Jed · Appendix · All-BugSigDB condition \N{MULTIPLICATION SIGN} site map"
            title = f"Corpus-wide condition \N{MULTIPLICATION SIGN} body-site map: {stem}"
            takeaway = (
                "Corpus-wide condition \N{MULTIPLICATION SIGN} anatomy signature map "
                f"({png_path.name}) across the regenerated BugSigDB human corpus."
            )
            caveat = "Descriptive coverage and guild-signal overview across all human studies."
            evidence_path = "data/derived/study_trait_profiles_wide.tsv"
        elif "/profile-condition-bars/" in rel:
            owner = "Jed"
            code_path = "scripts/create_profile_condition_enrichment_bars.py"
            tag = "Jed · Appendix · Condition profile signal"
            title = f"Condition profile trait signal: {stem}"
            takeaway = (
                "Condition-stratified profile trait signal and coverage diagnostic "
                f"({png_path.name}) on the rerun corpus."
            )
            caveat = (
                "Profile-level summary; coverage differences across conditions can influence "
                "raw signal totals."
            )
            evidence_path = "data/derived/study_trait_profiles_wide.tsv"
        elif "/enrichment-infographic-followup/" in rel:
            owner = "Jed / Ilyas"
            code_path = "scripts/analyze_enrichment_infographic_followup.py"
            tag = "Jed · Appendix · ORA organism-trait followup"
            title = f"Top-enriched traits vs organisms heatmap ({parent_name})"
            takeaway = (
                "Maps top ORA-enriched traits against contributing organisms "
                f"({png_path.name}) from the fresh ORA rerun."
            )
            caveat = "Shows KG-Microbe organism-trait annotations underlying top ORA hits."
            evidence_path = "data/derived/trait_enrichment_results.tsv"
        elif "/six-shared-traits/oral/" in rel:
            owner = "Jed"
            code_path = "scripts/create_oral_six_trait_heatmap.py"
            tag = "Jed · Appendix · Oral six-trait paginated view"
            title = f"Oral six-trait signed direction contrast ({png_path.stem})"
            takeaway = (
                f"Paginated slice ({png_path.name}) of the oral/dental six-shared-trait "
                "signed direction contrast heatmap."
            )
            caveat = (
                "No oral condition \N{MULTIPLICATION SIGN} site cell passes global FDR; shown "
                "for complete inspection of oral strata."
            )
            evidence_path = (
                "trait-visualization/trait-profile-investigation/six-shared-traits/"
                "oral/oral_six_trait_condition_body_site_contrasts.tsv"
            )
        elif "/six-shared-traits/" in rel:
            owner = "Jed"
            code_path = (
                "scripts/create_six_shared_trait_static_graphs.py"
                if "cytoscape_hub_graph" in rel
                else "scripts/create_six_shared_trait_cytoscape_heatmaps.py"
            )
            tag = "Jed · Appendix · Six-shared-traits view"
            title = f"Six-shared-traits view: {stem}"
            takeaway = (
                f"Complete/paginated or hub-network view ({png_path.name}) of the six "
                "shared traits on the rerun corpus."
            )
            caveat = "Companion view to the main-sequence six-shared-traits heatmaps and Sankeys."
            evidence_path = (
                "trait-visualization/trait-profile-investigation/six-shared-traits/"
                "six_shared_trait_condition_body_site_direction_contrasts.tsv"
            )
        elif "/taxa-coverage/" in rel:
            owner = "Jed"
            code_path = "scripts/create_taxa_coverage_histogram.py"
            tag = "Jed · Appendix · Taxa coverage"
            title = f"Taxa trait coverage diagnostic: {stem}"
            takeaway = (
                "Evaluates KG-Microbe trait count distribution and unique vs occurrence "
                f"coverage ({png_path.name})."
            )
            caveat = "Diagnostic of database annotation completeness across BugSigDB taxa."
            evidence_path = (
                "trait-visualization/trait-profile-investigation/taxa-coverage/README.md"
            )
        elif "/bsdb-taxa-only/" in rel:
            owner = "Jed"
            code_path = "scripts/create_trait_Heatmaps.py"
            tag = "Jed · Appendix · Taxa-only clustering"
            title = f"BugSigDB taxa-only clustered heatmap: {stem}"
            takeaway = (
                "Clusters top recurring genera/taxa across BugSigDB profiles "
                f"({png_path.name}) prior to trait projection."
            )
            caveat = "Taxon-occurrence clustering without KG-Microbe trait weighting."
            evidence_path = "data/derived/study_trait_profiles_wide.tsv"
        else:
            owner = "Jed"
            code_path = "scripts/regenerate_trait_profile_investigation.py"
            tag = "Jed · Appendix · Guild assignment"
            title = f"Guild definition & assignment summary: {stem}"
            takeaway = (
                "Reference table/summary of empirical guild definitions and profile "
                f"assignments ({png_path.name})."
            )
            caveat = "Documentation visual for guild membership rules used across Arm 1 analyses."
            evidence_path = "data/derived/study_trait_profiles_wide.tsv"

        appendix_slides.append(
            {
                "owner": owner,
                "title": title,
                "tag": tag,
                "figure_file_path": rel,
                "code_path": code_path,
                "takeaway": takeaway.strip(),
                "caveat": caveat.strip(),
                "methods_path": "",
                "evidence_path": evidence_path,
                "deck_section": "Appendix",
                "analysis": png_path.stem,
                "selection_reason": "Included in v3 Appendix after deduplication of "
                "redundant/superseded views.",
            }
        )

    return appendix_slides


# The shared pipeline steps every regenerated figure depends on.
UPSTREAM_COMMANDS = (
    "just download  # mtp.download.download_exports()",
    "just ingest  # mtp.ingest.ingest()",
    "just study-tsvs  # uv run python -m mtp.study_tsv",
    "just download-kg  # uv run --extra kg python scripts/download_kg_microbe.py",
    "just taxa-trait-matrix  # uv run python scripts/create_taxa_trait_matrix.py",
    "just study-trait-profiles  # uv run python scripts/create_study_trait_profiles.py",
)
# KG-Microbe files `just taxa-trait-matrix` reads (src/mtp/taxa_trait_matrix.py).
KG_INPUTS = ("data/kg-microbe/merged-kg_nodes.tsv", "data/kg-microbe/merged-kg_edges.tsv")
UPSTREAM_STEPS = (
    "step_01_download",
    "step_02_ingest",
    "step_03_study_tsvs",
    "step_04_taxa_trait_matrix",
    "step_05_study_trait_profiles",
)


def current_commands(series: Sequence[str]) -> list[str]:
    """A record's command series with each upstream recipe described as it runs today.

    The source YAML annotated them with a nonexistent `mtp.cli` and left out
    `just download-kg` (#615, #621).
    """
    current = {command.split("  #")[0]: command for command in UPSTREAM_COMMANDS}

    def described(command: str) -> str:
        recipe = command.split("  #")[0]
        if recipe in current:
            return current[recipe]
        # The recipe name is the instruction; a hint naming mtp.cli is wrong.
        return recipe if "mtp.cli" in command else command

    commands = [described(command) for command in series]
    recipes = [command.split("  #")[0] for command in commands]
    if "just taxa-trait-matrix" in recipes and "just download-kg" not in recipes:
        commands.insert(recipes.index("just taxa-trait-matrix"), current["just download-kg"])
    return commands


PLATE_COMMANDS = ("just cellulose-decr-figures", "just figure-review-plates")


def plate_provenance(rec: dict[str, Any], methods_path: str, inputs: Mapping[str, str]) -> None:
    """Point a plate's provenance record at the plate script and the inputs it read.

    The source YAML named data/derived/figure_review_2026-09-21/slides/generate_tufte_plates.py,
    which is not in the repository (#613, #621). Plates IV and V read the captions file that
    `just cellulose-decr-figures` writes, so that recipe runs first. A rebuild reads the
    previous build's record, so the two recipes are appended only once (#639), and the
    inputs are the ones `just figure-review-plates` recorded for this plate (#637).
    """
    rec["generator_script_path"] = "scripts/generate_tufte_plates.py"
    rec["just_recipe"] = "just figure-review-plates"
    # The slide's methods link, which #613 pointed at committed documents (#623)
    rec["methods_doc_path"] = methods_path
    series = [
        c
        for c in rec.get("command_series", [])
        if "generate_tufte_plates.py" not in c and c not in PLATE_COMMANDS
    ]
    rec["command_series"] = [*series, *PLATE_COMMANDS]
    rec["input_data_files"] = list(inputs)
    rec["input_data_sha256s"] = list(inputs.values())


def cellulose_figure_provenance(rec: dict[str, Any], tables: dict[str, str]) -> None:
    """Point a cellulose_decr figure's provenance record at what draws it (#556).

    `just cellulose-decr-figures` draws the figures from the committed tables and records
    the tables it read with their hashes. A record carried over from the source YAML still
    names the generator of AGY's deck, so every field that says how the figure was made is
    replaced, whichever branch built the record.
    """
    rec["generator_script_path"] = "src/mtp/cellulose_decr/figures.py"
    rec["just_recipe"] = "just cellulose-decr-figures"
    rec["command_series"] = [
        *UPSTREAM_COMMANDS,
        "just contrast-orientation",
        "just cellulose-decr",
        "just cellulose-decr-figures",
    ]
    rec["pipeline_step_ids"] = list(UPSTREAM_STEPS)
    rec["upstream_raw_inputs"] = ["data/raw/full_dump.csv", *KG_INPUTS]
    rec["input_data_files"] = [f"{CELLULOSE_DIR}/{name}" for name in tables]
    rec["input_data_sha256s"] = list(tables.values())
    rec["methods_doc_path"] = "methods/cellulose_decr.md"


def load_cellulose_captions(repo: Path) -> dict[str, dict[str, Any]]:
    """Slide text for the cellulose_decr figures, worded from its committed tables (#556).

    ``just cellulose-decr-figures`` writes it beside the figures, so the slides change with
    the tables instead of carrying hand-copied numbers.
    """
    path = repo / CELLULOSE_DIR / "figures" / "figure_captions.json"
    if not path.exists():
        raise FileNotFoundError(f"missing {path}; run `just cellulose-decr-figures` first")
    captions: dict[str, dict[str, Any]] = json.loads(path.read_text(encoding="utf-8"))
    return captions


def recorded_commit(path: Path) -> str | None:
    """A provenance file's code_commit when it is a clean commit, else None."""
    if not path.exists():
        return None
    commit = json.loads(path.read_text(encoding="utf-8")).get("code_commit") or ""
    return commit if re.fullmatch(r"[0-9a-f]{40}", commit) else None


def cellulose_revision(repo: Path) -> str | None:
    """The commit that drew the cellulose_decr figures, or None if it was modified or unknown.

    The cellulose slides link their code, methods and tables at that commit: those files
    are newer than the rerun revision the other slides cite, and a branch name moves (#556).
    """
    return recorded_commit(repo / CELLULOSE_DIR / "figures" / "figures_provenance.json")


def plate_revision(repo: Path) -> str | None:
    """The commit that drew the methods and Tufte plates, where their code is linked (#613).

    scripts/generate_tufte_plates.py is not in the rerun revision the plate slides' tables
    come from, so their CODE link needs the commit `just figure-review-plates` records.
    """
    return recorded_commit(repo / PLATES_DIR / "plates_provenance.json")


def figure_root(repo: Path, rel: str) -> Path:
    """Where a slide file is read: this checkout for the plates and the cellulose_decr files,
    which its own recipes draw, and the input root for everything else (#639)."""
    return REPO_ROOT if rel.startswith((PLATES_DIR, f"{CELLULOSE_DIR}/")) else repo


def git_blob_id(path: Path) -> str:
    """The object id git gives a file's content, to compare it with a committed blob."""
    data = path.read_bytes()
    return hashlib.sha1(b"blob %d\0" % len(data) + data, usedforsecurity=False).hexdigest()


def load_git_trees(folder: Path) -> dict[str, dict[str, str]]:
    """The blob ids of every committed file at each ref, from the `git ls-tree -r -z
    --full-tree <ref>` listings `just figure-review-deck` writes as ``<ref>.tree`` (#637)."""
    trees = {}
    for listing in sorted(folder.glob("*.tree")):
        tree = {}
        for entry in listing.read_text(encoding="utf-8").split("\0"):
            if entry:
                meta, name = entry.split("\t", 1)
                tree[name] = meta.split()[2]
        trees[listing.stem] = tree
    return trees


def revision_problem(
    trees: Mapping[str, Mapping[str, str]], ref: str | None, rel: str, path: Path
) -> str | None:
    """Why ``path`` is not the committed file ``rel`` at ``ref``, or None when it is (#637)."""
    if not ref:
        return "drawn without a clean commit"
    if ref not in trees:
        return f"no listing of {ref[:12]}"
    blob = trees[ref].get(rel)
    if blob is None:
        return f"not committed at {ref[:12]}"
    if not path.exists():
        return "missing"
    if git_blob_id(path) != blob:
        return f"differs from {ref[:12]}"
    return None


# Folders of gitignored outputs, which only their recorded hashes identify.
GITIGNORED = ("data/derived/", *CELLULOSE_FIGURES)


def link_problem(trees: Mapping[str, Mapping[str, str]], ref: str | None, rel: str) -> str | None:
    """Why a link to ``rel`` at ``ref`` would not open, or None when the file is there."""
    if not ref:
        return "no clean revision recorded"
    if ref not in trees:
        return f"no listing of {ref[:12]}"
    return None if rel in trees[ref] else f"not committed at {ref[:12]}"


def committed(rel: str, tree: Mapping[str, str] | None = None) -> bool:
    """Whether a slide file is checked against a revision: it is in the revision's tree, or
    it lies outside the gitignored folders, where a file missing from the tree is a dead
    link rather than an untracked output (#637)."""
    return (tree is not None and rel in tree) or not rel.startswith(GITIGNORED)


def plate_problems(
    repo: Path,
    record: Mapping[str, Any],
    trees: Mapping[str, Mapping[str, str]],
    cellulose: str | None,
    plates_dir: Path = REPO_ROOT / PLATES_DIR,
) -> dict[str, list[str]]:
    """Why each plate does not show the rerun its slide stamps, by plate file (#637).

    A plate must be drawn for COMMIT, from inputs that have not changed since, and its
    committed inputs must be the files at COMMIT, or at the commit that drew the
    cellulose_decr figures for the cellulose_decr tables.
    """
    problems: dict[str, list[str]] = {}
    if record.get("rerun_commit") != COMMIT:
        reason = (
            f"plates drawn for rerun {str(record.get('rerun_commit'))[:12]}, not {SHORT_COMMIT}"
        )
        return {plate: [reason] for plate in record.get("plates", [])}
    outputs = record.get("outputs", {})
    for plate in record.get("plates", []):
        found = []
        inputs = record.get("inputs", {}).get(plate)
        if inputs is None:
            found.append("no inputs recorded")
        for rel, sha in (inputs or {}).items():
            path = figure_root(repo, rel) / rel
            ref = cellulose if rel.startswith(f"{CELLULOSE_DIR}/") else COMMIT
            if not path.exists():
                found.append(f"{rel}: missing")
            elif sha256_file(path) != sha:
                found.append(f"{rel}: changed since the plate was drawn")
            elif committed(rel, trees.get(ref or "")):
                reason = revision_problem(trees, ref, rel, path)
                if reason:
                    found.append(f"{rel}: {reason}")
        # The plate the deck places must be the one drawn, not a later overwrite.
        for drawn in (plate, str(Path(plate).with_suffix(".svg"))):
            path = plates_dir / drawn
            if drawn not in outputs:
                found.append(f"{drawn}: no hash recorded")
            elif not path.exists() or sha256_file(path) != outputs[drawn]:
                found.append(f"{drawn}: not the plate that was drawn")
        if found:
            problems[plate] = found
    return problems


def cellulose_problems(
    record: Mapping[str, Any],
    figures: Sequence[str],
    trees: Mapping[str, Mapping[str, str]],
    revision: str | None,
    root: Path = REPO_ROOT,
) -> list[str]:
    """Why the cellulose_decr figures are not the committed tables drawn at ``revision``:
    a table changed since, differs from the revision, or a figure is not the one drawn
    (#637)."""
    found = []
    for name, sha in record.get("inputs", {}).items():
        rel = f"{CELLULOSE_DIR}/{name}"
        path = root / rel
        reason = (
            "changed since the figures were drawn"
            if sha256_file(path) != sha
            else revision_problem(trees, revision, rel, path)
        )
        if reason:
            found.append(f"{rel}: {reason}")
    outputs = record.get("outputs")
    if not isinstance(outputs, dict):
        return [*found, "the figures' provenance records no output hashes"]
    for rel in figures:
        name = rel.removeprefix(f"{CELLULOSE_DIR}/")
        path = root / rel
        if name not in outputs or not path.exists() or sha256_file(path) != outputs[name]:
            found.append(f"{rel}: not the figure that was drawn")
    return found


def check_slides(
    repo: Path,
    entries: Sequence[dict[str, Any]],
    trees: Mapping[str, Mapping[str, str]],
    tables: Iterable[Path],
) -> tuple[dict[str, list[tuple[str, str]]], list[str]]:
    """Each slide's linked files that are not what the slide claims, by slide figure; and
    the quoted tables that differ from COMMIT (#637).

    The figure and evidence table a slide shows and quotes must be the committed files at
    its revision. Its code and methods documents are only linked, so they must exist there.
    """
    unverified: dict[str, list[tuple[str, str]]] = {}
    for entry in entries:
        fig_rel = entry["figure_file_path"]
        links = slide_links(entry, fig_rel)
        linked = (
            (fig_rel, links["figure_url"] != "NA", links["commit"], True),
            (entry.get("evidence_path", ""), links["evidence_url"] != "NA", links["commit"], True),
            (entry.get("code_path", ""), bool(links["code_url"]), links["code_commit"], False),
            (entry.get("methods_path", ""), bool(links["methods_url"]), links["commit"], False),
        )
        found: list[tuple[str, str]] = []
        for rel, is_linked, ref, shown in linked:
            if rel and is_linked and rel not in {seen for seen, _ in found}:
                reason = (
                    revision_problem(trees, ref or None, rel, figure_root(repo, rel) / rel)
                    if shown
                    else link_problem(trees, ref or None, rel)
                )
                if reason:
                    found.append((rel, reason))
        if found:
            unverified[fig_rel] = found
    quoted = []
    for path in sorted(tables):
        rel = path.relative_to(repo).as_posix() if path.is_relative_to(repo) else str(path)
        if committed(rel, trees.get(COMMIT)):
            reason = revision_problem(trees, COMMIT, rel, path)
            if reason:
                quoted.append(f"{rel}: {reason}")
    return unverified, quoted


def refuse_or_review(problems: Sequence[str], allow_unverified: bool) -> list[str]:
    """Stop the build when files differ from the revisions their slides cite, unless a
    review build is allowed; then return the reason it is one (#637)."""
    if not problems:
        return []
    if not allow_unverified:
        shown = "\n  ".join(problems[:25])
        more = f"\n  ... and {len(problems) - 25} more" if len(problems) > 25 else ""
        raise SystemExit(
            f"{len(problems)} files differ from the revision their slides cite (#637):\n  "
            f"{shown}{more}\nCommit or restore them, redraw the plates and figures, or pass "
            "--allow-unverified for a review build."
        )
    return [f"{len(problems)} files differ from the revision their slides cite"]


def verified_links(
    links: Mapping[str, str], found: Sequence[tuple[str, str]], entry: Mapping[str, Any]
) -> dict[str, str]:
    """A slide's links without the files that differ from their revision (#637)."""
    failed = {rel for rel, _ in found}
    checked = dict(links)
    for key, rel, empty in (
        ("figure_url", entry["figure_file_path"], "NA"),
        ("evidence_url", entry.get("evidence_path", ""), "NA"),
        ("code_url", entry.get("code_path", ""), ""),
        ("methods_url", entry.get("methods_path", ""), ""),
    ):
        if rel in failed:
            checked[key] = empty
    return checked


def check_results_cite_the_deck(
    prov_doc: Mapping[str, Any], in_deck: set[str], source: Path
) -> None:
    """Refuse statistical results that cite a figure outside the deck (#638)."""
    stale = stale_figure_ids(prov_doc, in_deck)
    if stale:
        results = sorted(
            str(result.get("result_id"))
            for result in prov_doc.get("statistical_results", [])
            if set(result.get("figure_ids") or []) & set(stale)
        )
        raise ValueError(
            f"statistical results {results} cite figures that are not in the deck ({stale}). "
            f"Edit them in {source}, which has no generator, or restore the figures."
        )


def refreshed_hashes(
    repo: Path, files: Sequence[str], recorded: Sequence[str]
) -> tuple[list[str], list[str]]:
    """Each input file's current SHA-256, and the files whose recorded hash differs (#637).

    Carried-over records hold "path:sha256" entries or bare hashes in the order of their
    files; both are replaced by bare hashes. A missing file keeps its recorded hash.
    """
    old: dict[str, str] = {}
    for index, value in enumerate(recorded):
        path, sep, sha = value.rpartition(":")
        if sep and re.fullmatch(r"[0-9a-f]{64}", sha):
            old[path] = sha
        elif index < len(files):
            old[files[index]] = value
    current, changed = [], []
    for rel in files:
        path = figure_root(repo, rel) / rel
        sha = sha256_file(path) if path.exists() else old.get(rel, "")
        current.append(sha)
        if (rel in old and old[rel] != sha) or not path.exists():
            changed.append(rel)
    return current, changed


def slide_links(entry: dict[str, Any], fig_rel: str) -> dict[str, str]:
    """A slide's commit and its figure, code, methods and evidence links (#593, #613).

    ``link_ref`` is the revision of the slide's figure and tables: the rerun's, or the commit
    that drew a cellulose_decr figure. ``code_ref``, where set, is its code's. A ref of None
    links nothing.
    """
    ref = entry.get("link_ref", COMMIT)
    code_ref = entry.get("code_ref", ref)
    blob = f"{REPO_URL}/blob/{ref}" if ref else ""
    code_blob = f"{REPO_URL}/blob/{code_ref}" if code_ref else ""
    code_rel = entry.get("code_path", "")
    methods_rel = entry.get("methods_path", "")
    evidence_rel = entry.get("evidence_path", "")
    linked_figure = blob and not fig_rel.startswith(UNLINKED_FIGURES)
    linked_evidence = blob and evidence_rel and not evidence_rel.startswith("data/derived/")
    return {
        "commit": ref or "",
        "code_commit": code_ref or "",
        "figure_url": f"{blob}/{fig_rel}" if linked_figure else "NA",
        "code_url": f"{code_blob}/{code_rel}" if code_blob and code_rel else "",
        "methods_url": f"{blob}/{methods_rel}" if blob and methods_rel else "",
        "evidence_url": f"{blob}/{evidence_rel}" if linked_evidence else "NA",
    }


def revision_stamp(row: dict[str, str]) -> str:
    """A slide's revision: the rerun's, plus the plates' drawing commit, or the commit that
    drew a cellulose_decr figure (#623). A slide whose files differ from it says so (#637)."""
    stamp = base_revision_stamp(row)
    return f"{stamp} · NOT VERIFIED" if row.get("unverified") else stamp


def base_revision_stamp(row: dict[str, str]) -> str:
    if row["figure_file_path"].startswith(PLATES_DIR):
        drawn = row.get("code_commit", "")
        plates = f"plates drawn at {drawn[:12]}" if drawn else "plates drawn without a clean commit"
        return f"Rerun revision {row['commit'][:12]} · {plates}"
    if not row["figure_file_path"].startswith(CELLULOSE_FIGURES):
        return f"Rerun revision {row['commit'][:12]}"
    if row["commit"]:
        return f"Figures drawn at {row['commit'][:12]}"
    return "Figures drawn without a clean commit"


def cellulose_decr_slides(
    captions: dict[str, dict[str, Any]],
    folder: str,
    section: str,
    owner: str,
    revision: str | None,
) -> list[dict[str, Any]]:
    """One slide per cellulose_decr figure in ``folder``, in the order the figures are drawn.

    The presentation figures (figures_graphical/) replace the v5 Tufte plates, which had no
    generator and quoted superseded numbers; the statistical reference figures (figures/) go
    to the Appendix (#556, #570).
    """
    kind = "Figure" if section == "Main" else "Reference figure"
    slides = []
    for figure, caption in captions["figures"].items():
        if not figure.startswith(f"{folder}/"):
            continue
        stem = Path(figure).stem
        number = stem.split("_", 1)[0].removeprefix("fig")
        slides.append(
            {
                "owner": owner,
                "title": caption["title"],
                "tag": f"Cellulose deep-dive · {kind} {number}",
                "figure_file_path": f"{CELLULOSE_DIR}/{figure}",
                "code_path": "src/mtp/cellulose_decr/figures.py",
                "takeaway": caption["takeaway"],
                "caveat": caption["caveat"],
                "methods_path": "methods/cellulose_decr.md",
                "evidence_path": f"{CELLULOSE_DIR}/{caption['evidence']}",
                "deck_section": section,
                "analysis": stem,
                "selection_reason": "cellulose_decr figure drawn by `just "
                "cellulose-decr-figures` from the committed tables (#556).",
                # None when the figures were drawn from modified files: nothing to link.
                "link_ref": revision,
                "source_note": (
                    f"Drawn by `just cellulose-decr-figures` at {revision[:12]} from the "
                    "committed data/analysis/cellulose_decr tables."
                    if revision
                    else "Drawn without a clean recorded commit (modified code or tables, or "
                    "no commit recorded), so no revision is linked."
                ),
            }
        )
    if not slides:
        raise ValueError(f"figure_captions.json lists no figures under {folder}/")
    return slides


def placed_slides(
    items: Sequence[dict[str, Any]], root: Path, included: set[str], what: str
) -> list[dict[str, Any]]:
    """Curated slides in order. A missing figure stops the build, rather than dropping its
    slide and re-attaching its statistical results to another figure (#638)."""
    placed = []
    for item in items:
        fp = item["figure_file_path"]
        if not (root / fp).exists():
            raise FileNotFoundError(f"missing {what}: {root / fp} (#638)")
        if fp not in included:
            item["v1_slide"] = "NA"
            item["v2_slide"] = "NA"
            item["prior_result_url"] = ""
            placed.append(item)
            included.add(fp)
    return placed


def build_deck(
    repo: Path,
    slides_dir: Path,
    provenance_dir: Path,
    trees: Mapping[str, Mapping[str, str]] | None = None,
    allow_unverified: bool = False,
    review: Sequence[str] = (),
) -> None:
    """Build the deck from ``repo``, checking its committed files against ``trees`` (#637).

    ``trees`` holds the blob ids of the revisions the slides cite. A file that differs from
    its revision stops the build unless ``allow_unverified``; the deck is then a review
    build, like one with a ``review`` reason, and the files that differ are not linked.
    """
    if trees is None and not allow_unverified:
        raise SystemExit("no git listings to check the slides against; pass --allow-unverified")
    TABLES_READ.clear()
    metrics = load_rerun_metrics(repo)
    cellulose = load_cellulose_captions(REPO_ROOT)
    metrics["cellulose_statements"] = cellulose["statements"]
    revision = cellulose_revision(REPO_ROOT)
    plate_record = json.loads(
        (REPO_ROOT / PLATES_DIR / "plates_provenance.json").read_text(encoding="utf-8")
    )
    cellulose_record = json.loads(
        (REPO_ROOT / CELLULOSE_DIR / "figures" / "figures_provenance.json").read_text(
            encoding="utf-8"
        )
    )
    # Counts Methods Plate II draws, which its slide quotes (#640)
    metrics["plate_counts"] = plate_record.get("counts", {})

    # The September decks' manifests are inputs; slide_manifest.tsv itself is now this
    # script's output, so it must not stand in for a missing v1 manifest.
    previous = repo / DECK_DIR / "slides"
    v1_manifest_path = previous / "slide_manifest_v1.tsv"
    v2_manifest_path = previous / "slide_manifest_v2.tsv"
    for manifest_path in (v1_manifest_path, v2_manifest_path):
        if not manifest_path.exists():
            raise FileNotFoundError(f"missing earlier deck manifest {manifest_path}")
    v1_rows = read_tsv(v1_manifest_path)
    v2_rows = read_tsv(v2_manifest_path)

    v1_by_fig = {r["figure_file_path"]: r for r in v1_rows}
    v2_by_fig = {r["figure_file_path"]: r for r in v2_rows}

    pruned_prev_main = {
        "data/derived/Sankey_Diagrams/facultative-aerobe_sankey_inc.png",
        "data/derived/Sankey_Diagrams/anaerobe_sankey_dec.png",
        "methods/figures/table1_traits_condition_heatmap.png",
        "methods/figures/table1_traits_focus_top_conditions.png",
        "methods/figures/table1_traits_body_site_facets.png",
        "trait-visualization/trait-profile-investigation/six-shared-traits/"
        "individual_conditions_x_body_sites_positive_z_direction_contrasts_page_04.png",
    }
    demoted_to_appendix_prev_main = {
        "trait-visualization/trait-profile-investigation/parkinsons-disease/"
        "pd_contributor_taxa_fecal.png",
        "trait-visualization/trait-profile-investigation/parkinsons-disease/"
        "pd_leave_one_taxon_out_fecal.png",
        "trait-visualization/trait-profile-investigation/PD-cellulose/"
        "pd_trait_contributors_species_family_heatmap.png",
        "trait-visualization/trait-profile-investigation/PD-cellulose/"
        "pd_trait_contributors_trait_cocarriage.png",
    }

    prev_main_figs: list[str] = []
    demoted_main_figs: list[str] = []
    for r in v2_rows:
        fp = r["figure_file_path"]
        if r.get("deck_section") == "Main":
            if fp in pruned_prev_main:
                continue
            elif fp in demoted_to_appendix_prev_main:
                demoted_main_figs.append(fp)
            else:
                prev_main_figs.append(fp)
    for r in v1_rows:
        fp = r["figure_file_path"]
        if (
            r.get("deck_section") == "Main"
            and fp not in prev_main_figs
            and fp not in demoted_main_figs
        ):
            if fp in pruned_prev_main:
                continue
            elif fp in demoted_to_appendix_prev_main:
                demoted_main_figs.append(fp)
            else:
                prev_main_figs.append(fp)

    prev_appendix_figs: list[str] = []
    for r in v2_rows:
        if r.get("deck_section") == "Appendix":
            prev_appendix_figs.append(r["figure_file_path"])
    for r in v1_rows:
        if r.get("deck_section") == "Appendix" and r["figure_file_path"] not in prev_appendix_figs:
            prev_appendix_figs.append(r["figure_file_path"])

    def make_prev_entry(fig_path: str, section_override: str) -> dict[str, str]:
        base = dict(v2_by_fig.get(fig_path) or v1_by_fig[fig_path])
        v1_num = v1_by_fig[fig_path]["slide"] if fig_path in v1_by_fig else "NA"
        v2_num = v2_by_fig[fig_path]["slide"] if fig_path in v2_by_fig else "NA"
        title, takeaway, caveat = update_previous_slide_text(base, metrics)
        base["title"] = title
        base["takeaway"] = takeaway
        base["caveat"] = caveat
        base["deck_section"] = section_override
        base["v1_slide"] = v1_num
        base["v2_slide"] = v2_num
        if v1_num != "NA" and v2_num != "NA":
            base["selection_reason"] = (
                f"Retained from v1 (slide {v1_num}) and v2 (slide {v2_num}); updated with fresh "
                "rerun data."
            )
        elif v1_num != "NA":
            base["selection_reason"] = (
                f"Restored from v1 (slide {v1_num}, omitted in v2); updated with fresh rerun data."
            )
        else:
            base["selection_reason"] = (
                f"Retained from v2 (slide {v2_num}); updated with fresh rerun data."
            )
        return base

    ordered_entries: list[dict[str, str]] = []
    already_included: set[str] = set()

    new_slides_all = get_new_slides_metadata(metrics)
    tufte_and_methods_slides = [
        item
        for item in new_slides_all
        if item["figure_file_path"].startswith("data/derived/tufte_plates/")
    ]
    other_new_main_slides = [
        item
        for item in new_slides_all
        if not item["figure_file_path"].startswith("data/derived/tufte_plates/")
    ]

    # 1. Opening Methods & Tufte Synthesis Plates (9 slides: 2 Methods + 7 Tufte plates)
    plates_commit = plate_revision(REPO_ROOT)
    for item in tufte_and_methods_slides:
        item["code_ref"] = plates_commit
    ordered_entries += placed_slides(
        tufte_and_methods_slides,
        REPO_ROOT,
        already_included,
        "plate (run `just figure-review-plates`)",
    )

    # 1b. Cellulose deep-dive presentation figures (4 slides from `just cellulose-decr-figures`;
    # a missing figure stops the build below rather than dropping the slide)
    cellulose_main = cellulose_decr_slides(
        cellulose, "figures_graphical", "Main", "Marcin / Noel / Ilyas / Jed", revision
    )
    for item in cellulose_main:
        item["v1_slide"] = "NA"
        item["v2_slide"] = "NA"
        item["prior_result_url"] = ""
        ordered_entries.append(item)
        already_included.add(item["figure_file_path"])

    # 2. Deduplicated Previous Main slides (15 canonical Main slides from v1/v2)
    for fp in prev_main_figs:
        ordered_entries.append(make_prev_entry(fp, "Main"))
        already_included.add(fp)

    # 3. Remaining Curated New Main slides (7 slides across Publication Deck, Noel, Ilyas, Jed)
    ordered_entries += placed_slides(
        other_new_main_slides, repo, already_included, "curated Main figure"
    )

    # 4. Previous Appendix Sankey slides (18 slides) + 4 exploratory single-metric PD variants moved
    # to Appendix
    for fp in prev_appendix_figs + demoted_main_figs:
        ordered_entries.append(make_prev_entry(fp, "Appendix"))
        already_included.add(fp)

    # 5. All remaining distinct non-redundant figures across the rerun (Appendix)
    for item in get_remaining_appendix_slides(already_included, repo, metrics):
        fp = item["figure_file_path"]
        if fp not in already_included and (repo / fp).exists():
            item["v1_slide"] = "NA"
            item["v2_slide"] = "NA"
            item["prior_result_url"] = ""
            ordered_entries.append(item)
            already_included.add(fp)

    # 6. Cellulose deep-dive statistical reference figures (3 slides)
    for item in cellulose_decr_slides(cellulose, "figures", "Appendix", "Marcin", revision):
        item["v1_slide"] = "NA"
        item["v2_slide"] = "NA"
        item["prior_result_url"] = ""
        ordered_entries.append(item)
        already_included.add(item["figure_file_path"])

    # Check every committed file the slides link or quote against the revision it cites,
    # before anything is written (#637).
    # Every statistical result must stay with the figure it describes (#638).
    prov_source = repo / DECK_DIR / "provenance" / PROVENANCE_NAME
    if not prov_source.exists():
        raise FileNotFoundError(f"missing figure provenance {prov_source}")
    with prov_source.open("r", encoding="utf-8") as f:
        prov_doc = yaml.safe_load(f)
    check_results_cite_the_deck(
        prov_doc, {entry["figure_file_path"] for entry in ordered_entries}, prov_source
    )

    checks = trees or {}
    unverified, quoted = check_slides(repo, ordered_entries, checks, TABLES_READ)
    for plate, found in plate_problems(repo, plate_record, checks, revision).items():
        unverified.setdefault(f"{PLATES_DIR}{plate}", []).extend(
            (f"{PLATES_DIR}{plate}", reason) for reason in found
        )
    cellulose_figures_shown = [
        entry["figure_file_path"]
        for entry in ordered_entries
        if entry["figure_file_path"].startswith(CELLULOSE_FIGURES)
    ]
    for reason in cellulose_problems(cellulose_record, cellulose_figures_shown, checks, revision):
        for figure in cellulose_figures_shown:
            unverified.setdefault(figure, []).append(("cellulose_decr figures", reason))
    cover_found = [
        (rel, reason)
        for rel, reason in (
            (COVER_FIGURE, revision_problem(checks, COMMIT, COVER_FIGURE, repo / COVER_FIGURE)),
            (COVER_CODE, link_problem(checks, COMMIT, COVER_CODE)),
        )
        if reason
    ]
    problems = [
        f"{rel}: {reason}" + ("" if rel == fig else f" (slide of {fig})")
        for fig, found in unverified.items()
        for rel, reason in found
    ]
    problems += [f"cover link {rel}: {reason}" for rel, reason in cover_found]
    problems += [f"quoted table {item}" for item in quoted]
    review = list(review)
    if trees is None:
        review.append("nothing was checked against its revision")
    review += refuse_or_review(problems, allow_unverified)

    # Finalize manifest rows, copy fresh PNGs to assets/ and vector SVGs to assets_svg/, and compute
    # SHA-256
    (slides_dir / "assets").mkdir(parents=True, exist_ok=True)
    (slides_dir / "assets_svg").mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, str]] = []
    used_asset_names: set[str] = set()
    figure_states: list[str] = []

    for idx, entry in enumerate(ordered_entries, start=2):
        fig_rel = entry["figure_file_path"]
        src_png = figure_root(repo, fig_rel) / fig_rel
        if not src_png.exists():
            raise FileNotFoundError(f"Missing rerun figure on disk: {src_png}")

        # Resolve true vector .svg for every slide
        cand_svg1 = src_png.with_suffix(".svg")
        cand_svg2 = previous / "svg_sources" / Path(fig_rel).with_suffix(".svg")
        src_svg = cand_svg1 if cand_svg1.exists() else (cand_svg2 if cand_svg2.exists() else None)
        if src_svg is None:
            raise FileNotFoundError(
                f"Missing vector SVG for slide {idx}: checked {cand_svg1} and {cand_svg2}"
            )

        svg_asset_rel = f"assets_svg/v6_slide_{idx:03d}_{src_png.stem}.svg"
        dst_svg = slides_dir / svg_asset_rel
        shutil.copy2(src_svg, dst_svg)

        # Also copy primary file to assets/ (or assets_svg/ if primary is .svg)
        asset_rel = svg_asset_rel
        used_asset_names.add(asset_rel)
        if src_png.suffix.lower() == ".png":
            dst_png = slides_dir / f"assets/v6_slide_{idx:03d}_{src_png.name}"
            dst_png.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src_png, dst_png)

        src_sha = sha256_file(src_png)

        code_rel = entry.get("code_path", "")
        methods_rel = entry.get("methods_path", "")
        evidence_rel = entry.get("evidence_path", "")

        links = slide_links(entry, fig_rel)
        ref = links["commit"]
        fig_url = links["figure_url"]
        code_url = links["code_url"]
        methods_url = links["methods_url"]
        evidence_url = links["evidence_url"]
        found = unverified.get(fig_rel, [])
        checked = verified_links(links, found, entry)
        fig_url = checked["figure_url"]
        evidence_url = checked["evidence_url"]
        code_url = checked["code_url"]
        methods_url = checked["methods_url"]
        figure_state = (
            "verified"
            if fig_url != "NA"
            else "differs"
            if any(rel == fig_rel for rel, _ in found)
            else "untracked"
        )
        figure_states.append(figure_state)
        # The cellulose slides' own note says how they were drawn; the September
        # manifests' notes are kept only as history (#637).
        legacy = entry.get("source_note", "")
        if found:
            listed = "; ".join(f"{rel}: {reason}" for rel, reason in found)
            source_note = (
                f"Not verified against {ref[:12] if ref else 'a revision'}: {listed}. Files "
                "that differ are not linked (review build)."
            )
        elif fig_rel.startswith(CELLULOSE_FIGURES):
            source_note = legacy
        elif figure_state == "verified":
            source_note = f"The figure matches the committed file at {ref[:12]}." + (
                f" Earlier note: {legacy}" if legacy else ""
            )
        else:
            source_note = (
                f"From the rerun at {SHORT_COMMIT}; a gitignored output, so it cannot be "
                "checked against that revision. Its SHA-256 is recorded in the provenance."
                + (f" Earlier note: {legacy}" if legacy else "")
            )

        manifest_rows.append(
            {
                "slide": str(idx),
                "owner": entry["owner"],
                "title": entry["title"],
                "tag": entry["tag"],
                "figure_file_path": fig_rel,
                "code_path": code_rel,
                "commit": ref,
                "code_commit": links["code_commit"],
                "branch": BRANCH,
                "figure_url": fig_url,
                "code_url": code_url,
                "asset": asset_rel,
                "takeaway": entry["takeaway"],
                "caveat": entry["caveat"],
                "methods_path": methods_rel,
                "methods_url": methods_url,
                "source_note": source_note,
                "unverified": "; ".join(f"{rel}: {reason}" for rel, reason in found),
                "figure_sha256": src_sha,
                "deck_section": entry.get("deck_section", "Main"),
                "analysis": entry.get("analysis", ""),
                "evidence_path": evidence_rel,
                "evidence_url": evidence_url,
                "v1_slide": entry.get("v1_slide", "NA"),
                "v2_slide": entry.get("v2_slide", "NA"),
                "prior_result_url": entry.get("prior_result_url", ""),
                "selection_reason": entry.get("selection_reason", ""),
            }
        )

    # Write slide_manifest_v3.tsv and update slide_manifest.tsv
    fieldnames = [
        "slide",
        "owner",
        "title",
        "tag",
        "figure_file_path",
        "code_path",
        "commit",
        "code_commit",
        "branch",
        "figure_url",
        "code_url",
        "asset",
        "takeaway",
        "caveat",
        "methods_path",
        "methods_url",
        "source_note",
        "unverified",
        "figure_sha256",
        "deck_section",
        "analysis",
        "evidence_path",
        "evidence_url",
        "v1_slide",
        "v2_slide",
        "prior_result_url",
        "selection_reason",
    ]
    for manifest_name in ("slide_manifest_v6.tsv", "slide_manifest.tsv"):
        out_tsv = slides_dir / manifest_name
        with out_tsv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(manifest_rows)

    # Build HTML (selected_figures_v6.html and selected_figures.html)
    total_slides = len(manifest_rows) + 1
    n_prev = sum(1 for r in manifest_rows if r["v1_slide"] != "NA" or r["v2_slide"] != "NA")
    n_methods = sum(1 for r in manifest_rows if r["analysis"].startswith("methods_plate"))
    n_tufte = sum(1 for r in manifest_rows if r["analysis"].startswith("tufte_plate"))
    n_new = len(manifest_rows) - n_prev
    n_cellulose = sum(
        1 for r in manifest_rows if r["figure_file_path"].startswith(CELLULOSE_FIGURES)
    )
    n_reference = n_cellulose - len(cellulose_main)
    n_checked = figure_states.count("verified")
    n_differs = figure_states.count("differs")
    n_untracked = figure_states.count("untracked")
    # One line, so the cover header keeps its height; deck_validation_v6.json lists them all.
    review_html = (
        '<div class="tag" style="color:#b8332a;max-width:1520px;white-space:nowrap;'
        'overflow:hidden;text-overflow:ellipsis">Review build, not for distribution: '
        f"{html.escape('; '.join(review))}</div>"
        if review
        else ""
    )
    failed_cover = {rel for rel, _ in cover_found}
    cover_href = "#" if COVER_FIGURE in failed_cover else f"{GITHUB_BLOB}/{COVER_FIGURE}"
    cover_code_href = "#" if COVER_CODE in failed_cover else f"{GITHUB_BLOB}/{COVER_CODE}"

    # Extract cover table rows from rerun contrasts
    by_key = metrics["contrasts_by_key"]
    cover_leads = [
        (
            "Parkinson\N{RIGHT SINGLE QUOTATION MARK}s disease \N{MULTIPLICATION SIGN} Feces",
            "Cellulose catabolism",
            by_key.get(
                (
                    "individual_conditions_direction_contrast",
                    "Parkinson's disease",
                    "Feces",
                    "cellulose catabolic process",
                )
            ),
        ),
        (
            "HIV infection \N{MULTIPLICATION SIGN} Feces",
            "Cellulose catabolism",
            by_key.get(
                (
                    "individual_conditions_direction_contrast",
                    "HIV infection",
                    "Feces",
                    "cellulose catabolic process",
                )
            ),
        ),
        (
            "Colorectal cancer \N{MULTIPLICATION SIGN} Feces",
            "Cellulose catabolism",
            by_key.get(
                (
                    "individual_conditions_direction_contrast",
                    "Colorectal cancer",
                    "Feces",
                    "cellulose catabolic process",
                )
            ),
        ),
        (
            "Parkinson\N{RIGHT SINGLE QUOTATION MARK}s disease \N{MULTIPLICATION SIGN} Feces",
            "Moderately halophilic",
            by_key.get(
                (
                    "individual_conditions_direction_contrast",
                    "Parkinson's disease",
                    "Feces",
                    "moderately halophilic",
                )
            ),
        ),
    ]
    cover_tr_html = []
    for label, trait_disp, r in cover_leads:
        if r:
            delta_pp = float(r["increased_minus_decreased_positive_fraction"]) * 100.0
            sign_str = "\N{MINUS SIGN}" if delta_pp < 0 else "+"
            delta_str = f"{sign_str}{abs(delta_pp):.2f} pp"
            studies_str = r["study_count"]
            q_str = f"{float(r['fisher_q_value']):.5f}"
        else:
            delta_str, studies_str, q_str = "NA", "NA", "NA"
        cover_tr_html.append(
            f"<tr><td>{html.escape(label)}</td><td>{html.escape(trait_disp)}</td>"
            f"<td>{delta_str}</td><td>{studies_str}</td><td>{q_str}</td></tr>"
        )

    pd_w = metrics.get("pd_halo_wilcoxon", {})
    pd_w_q = float(pd_w.get("q_value", 0.05609))
    pd_w_p = float(pd_w.get("p_value", 0.000167))

    statements = cellulose["statements"]
    cellulose_note = statements["overlap"] + (
        f"; {statements['fragile']}" if statements["fragile"] else ""
    )
    html_lines = [
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        "<title>Microbial trait patterns — v6 complete vector deduplicated rerun figure "
        "deck</title><style>",
        "*{box-sizing:border-box}html,body{margin:0;padding:0;background:#edf1f5;color:#142335;font-family:Arial,Helvetica,sans-serif}",
        "@page{size:1600px 1000px;margin:0}",
        ".slide{width:1600px;height:1000px;padding:28px "
        "40px;background:white;display:grid;grid-template-rows:minmax(70px,auto) minmax(0,1fr) "
        "minmax(70px,auto) "
        "86px;gap:10px;break-after:page;page-break-after:always;overflow:hidden}",
        ".slide:last-child{break-after:auto;page-break-after:auto}",
        "header{display:flex;flex-direction:column;justify-content:center;border-bottom:2px solid "
        "#dde5ec}",
        ".tag{font-size:13px;font-weight:700;letter-spacing:1.4px;text-transform:uppercase;color:#17617d;margin-bottom:7px}",
        "h1{font-size:32px;line-height:1.13;margin:0;font-weight:700;color:#12263b}",
        ".visual{min-height:0;min-width:0;display:flex;justify-content:center;align-items:center;overflow:hidden}",
        ".visual "
        "img{display:block;max-height:100%;max-width:100%;width:100%;height:100%;object-fit:contain}",
        ".note{font-size:17px;line-height:1.35;border-left:4px solid #16798e;padding:9px "
        "15px;background:#f2f7fa;display:flex;flex-direction:column;justify-content:center}",
        ".note p{margin:0}.note p+p{margin-top:3px;color:#4c5d6c}",
        'footer{border-top:1px solid #c9d4df;padding-top:9px;font:12px/1.3 "Courier '
        'New",monospace;color:#405466;display:flex;flex-direction:column;justify-content:flex-start;gap:3px}',
        ".path{display:grid;grid-template-columns:63px 1fr;column-gap:8px;align-items:start}",
        ".path b{color:#1e6178;font-size:11px;padding-top:1px}",
        "a{color:inherit;text-decoration:none;overflow-wrap:anywhere;word-break:break-word}",
        ".stamp{margin-top:auto;display:flex;justify-content:space-between;gap:20px;font:12px "
        "Arial,sans-serif;color:#617589;padding-top:4px}",
        ".cover{grid-template-rows:minmax(125px,auto) minmax(0,1fr) minmax(86px,auto) 86px}",
        # One line for the title, so a review build's extra tag line still fits (#637)
        ".cover h1{font-size:40px;line-height:1.08}.cover .tag{font-size:15px}",
        ".cover-body{display:flex;flex-direction:column;justify-content:center}",
        ".deck-subtitle{font-size:22px;line-height:1.45;color:#405568;max-width:1420px;margin:0 0 "
        "22px}",
        "table{border-collapse:collapse;width:100%;font-size:23px;line-height:1.38}",
        "th{font-size:16px;text-transform:uppercase;letter-spacing:.8px;text-align:left;"
        "color:#436477;border-bottom:2px "
        "solid #abc2d0;padding:12px 15px}",
        "td{padding:16px 15px;border-bottom:1px solid "
        "#dce6eb}td:last-child,th:last-child{text-align:right}",
        "tbody "
        "tr:nth-child(odd){background:#f3f7fa}.lead{font-size:19px;line-height:1.45;color:#32536b;margin-top:20px}",
        "@media screen{.slide{margin:20px auto;box-shadow:0 3px 20px #0002}}@media "
        "print{html,body{background:white;-webkit-print-color-adjust:exact;print-color-adjust:exact}.slide{margin:0}}",
        "</style></head><body>",
        f'<section class="slide cover"><header>{review_html}<div class="tag">Version 6 · '
        "Complete Vector Deduplicated Rerun Deck · Marcin / Jed / Noel / Ilyas · 17 October "
        "2026</div>"
        f"<h1>Three-arm microbial trait profiling deck ({len(manifest_rows)} figures across "
        f"{total_slides} slides)</h1></header>"
        '<main class="cover-body"><p class="deck-subtitle">Opens with '
        f"{n_methods} high-level visual Methods "
        f"architectures, {n_tufte} simplified Tufte analytical plates, and "
        f"{len(cellulose_main)} cellulose deep-dive figures, followed by {n_prev} deduplicated "
        f"v1/v2 slides and {n_new - n_methods - n_tufte - n_cellulose} distinct new rerun "
        "figures across all three arms on <code>main</code> "
        f"(<code>{SHORT_COMMIT}</code>, BugSigDB vintage <code>2026-04-24_00:41_UTC</code> "
        "\N{MULTIPLICATION SIGN} KG-Microbe <code>20260120</code>, 4,670 human profiles across "
        f"2,335 paired experiments and 959 studies), and {n_reference} cellulose reference "
        "figures drawn from the committed cellulose_decr tables.</p>"
        "<table><thead><tr><th>Leading condition \N{MULTIPLICATION SIGN} "
        "anatomy</th><th>Trait</th><th>Δ positive-z profiles</th><th>Studies</th><th>BH "
        "q</th></tr></thead>"
        f"<tbody>{''.join(cover_tr_html)}</tbody></table>"
        '<p class="lead">Δ = increased minus decreased proportion of profiles with a positive '
        f"trait z-score ({metrics['ind_total_tests']:,}-row profile-level Fisher family).<br>"
        "Arm 2 taxon-level ORA: 66,240 enrichments (fold ≥ 1.5, BH q ≤ 0.001) and 23 depletions "
        "(fold ≤ 0.5, BH q ≤ 0.001) across 1,480,390 tests. Arm 3 study-level Table 1: all 6 "
        "traits pass BH q &lt; 0.01 across 959 studies (and in 569 disease-vs-healthy oriented "
        "studies, where Cellulose ↓ reaches 60.8%, p = 4.21 \N{MULTIPLICATION SIGN} 10⁻⁶, with 4 "
        "BY-significant log-ratios including CEL/HAL BY q = 0.0155); oxygen-guild McNemar OR = "
        "1.448 (420 vs 290 discordant studies, p = 1.21 \N{MULTIPLICATION SIGN} 10⁻⁶).</p></main>"
        '<aside class="note"><p><strong>Evidence across scales:</strong> In the '
        f"{metrics['pd_net_total']}-test fecal study Wilcoxon screen, PD moderately halophilic "
        f"{'leads' if metrics.get('pd_net_lead') is metrics['pd_halo_wilcoxon'] else 'has'} "
        f"({float(metrics['pd_halo_wilcoxon']['mean_contrast']):+.3f} mean contrast across "
        f"{metrics['pd_halo_wilcoxon']['paired_study_count']} studies, "
        f"{metrics['pd_halo_wilcoxon']['positive_study_count']}↑:"
        f"{metrics['pd_halo_wilcoxon']['negative_study_count']}↓, p = {pd_w_p:.2e}, BH q = "
        f"{pd_w_q:.4f}). In the cellulose deep-dive, "
        f"{html.escape(cellulose_note)}.</p>"
        f"<p>{n_checked} linked figures match the committed files at the revisions their "
        f"slides cite; {n_untracked} are gitignored outputs, identified by SHA-256 in the "
        "provenance"
        + (f"; {n_differs} differ from their revision and are not linked" if n_differs else "")
        + ". Trait annotations describe organism potential, not measured in vivo metabolic "
        "flux.</p></aside>"
        '<footer><div class="path"><b>FIGURE</b><a '
        f'href="{cover_href}">{COVER_FIGURE}</a></div>'
        f'<div class="path"><b>CODE</b><a href="{cover_code_href}">{COVER_CODE}</a></div>'
        '<div class="stamp"><span>CultureBotAI/microbial-trait-profiling · Jed / Noel / Ilyas · '
        f"Rerun revision {SHORT_COMMIT}{' · REVIEW BUILD' if review else ''} · v6 vector"
        "</span><span>01 / "
        f"{total_slides:02d}</span></div></footer></section>",
    ]

    for row in manifest_rows:
        s_num = int(row["slide"])
        owner = row["owner"]
        tag_full = f"{owner} · {row['tag']}" if not row["tag"].startswith(owner) else row["tag"]
        fig_href = row["figure_url"] if row["figure_url"] != "NA" else "#"
        code_href = row["code_url"] if row["code_url"] else "#"

        prior_div = ""
        if row.get("prior_result_url"):
            prior_div = (
                f'<div class="path"><b>PRIOR</b><a href="{html.escape(row["prior_result_url"])}">'
                f"{html.escape(prior_label(row['prior_result_url'], metrics['studies']))}</a></div>"
            )

        lineage_tag = []
        if row["v1_slide"] != "NA":
            lineage_tag.append(f"v1:{row['v1_slide']}")
        if row["v2_slide"] != "NA":
            lineage_tag.append(f"v2:{row['v2_slide']}")
        lineage_str = f" ({', '.join(lineage_tag)})" if lineage_tag else " (new in v3\N{EN DASH}v6)"

        html_lines.append(
            f'<section class="slide"><header><div class="tag">{html.escape(tag_full)}</div>'
            f"<h1>{html.escape(row['title'])}</h1></header>"
            f'<main class="visual"><img src="{html.escape(row["asset"])}" '
            f'alt="{html.escape(row["title"])}"></main>'
            "<aside "
            f'class="note"><p>{html.escape(row["takeaway"])}</p><p>{html.escape(row["caveat"])}</'
            "p></aside>"
            '<footer><div class="path"><b>FIGURE</b><a '
            f'href="{html.escape(fig_href)}">{html.escape(row["figure_file_path"])}</a></div>'
            '<div class="path"><b>CODE</b><a '
            f'href="{html.escape(code_href)}">{html.escape(row["code_path"])}</a></div>'
            f"{prior_div}"
            '<div class="stamp"><span>CultureBotAI/microbial-trait-profiling · '
            f"{html.escape(owner)} · {revision_stamp(row)}"
            f"{' · REVIEW BUILD' if review else ''} · v6 vector"
            f"{lineage_str}</span>"
            f"<span>{s_num:02d} / {total_slides:02d}</span></div></footer></section>"
        )

    html_lines.append("</body></html>")
    html_content = "\n".join(html_lines) + "\n"

    for html_name in ("selected_figures_v6.html", "selected_figures.html"):
        (slides_dir / html_name).write_text(html_content, encoding="utf-8")

    # Synchronize figure_provenance_and_results.yaml with all 153 manifest_rows;
    # `just figure-review-deck` then validates it against its LinkML schema.
    prov_yaml_path = provenance_dir / PROVENANCE_NAME
    cellulose_figures = cellulose_record
    prov_doc["dataset_id"] = "mtp_v6_three_arm_figure_provenance_2026_10_17"
    prov_doc["title"] = (
        "Microbial Trait Profiling v6 Complete Three-Arm Figure Provenance & Standardized "
        "Statistical Results"
    )

    existing_by_fp = {rec["figure_file_path"]: rec for rec in prov_doc.get("figures", [])}
    changed_inputs: list[dict[str, str]] = []
    old_id_to_new_id: dict[str, str] = {}
    new_figures_list = []

    for r in manifest_rows:
        s_num = int(r["slide"])
        new_fid = f"fig_slide_{s_num:03d}"
        fp = r["figure_file_path"]
        if fp in existing_by_fp:
            rec = dict(existing_by_fp[fp])
            old_id_to_new_id[rec["figure_id"]] = new_fid
            rec["figure_id"] = new_fid
            rec["slide_number"] = s_num
            rec["deck_section"] = r["deck_section"]
            rec["owner"] = r["owner"]
            rec["tag"] = r["tag"]
            rec["title"] = r["title"]
            rec["figure_file_path"] = fp
            rec["figure_sha256"] = r["figure_sha256"]
            rec["slide_asset_path"] = f"data/derived/figure_review_2026-09-21/slides/{r['asset']}"
            rec["takeaway"] = r["takeaway"]
            rec["caveat"] = r["caveat"]
        else:
            ev_rel = r.get("evidence_path", "")
            ev_path = figure_root(repo, ev_rel) / ev_rel
            ev_sha = sha256_file(ev_path) if ev_rel and ev_path.exists() else ""
            rec = {
                "figure_id": new_fid,
                "slide_number": s_num,
                "deck_section": r["deck_section"],
                "analytical_arm": "cross_arm_synthesis",
                "owner": r["owner"],
                "tag": r["tag"],
                "title": r["title"],
                "figure_file_path": fp,
                "figure_sha256": r["figure_sha256"],
                "slide_asset_path": f"data/derived/figure_review_2026-09-21/slides/{r['asset']}",
                "generator_script_path": r["code_path"],
                "just_recipe": f"uv run python {r['code_path']}",
                "command_series": [*UPSTREAM_COMMANDS, f"uv run python {r['code_path']}"],
                "pipeline_step_ids": list(UPSTREAM_STEPS),
                "upstream_raw_inputs": ["data/raw/full_dump.csv", *KG_INPUTS],
                "input_data_files": [ev_rel] if ev_rel else [],
                "input_data_sha256s": [ev_sha] if ev_sha else [],
                "methods_doc_path": r.get("methods_path", ""),
                "takeaway": r["takeaway"],
                "caveat": r["caveat"],
                "v1_slide_number": r.get("v1_slide", "NA"),
                "v2_slide_number": r.get("v2_slide", "NA"),
                "statistical_result_ids": [],
            }
        # The source YAML placed the KG files under data/raw, which holds none (#609), and
        # described the upstream recipes with a nonexistent mtp.cli (#621).
        rec["upstream_raw_inputs"] = [
            path.replace("data/raw/merged-kg_", "data/kg-microbe/merged-kg_")
            for path in rec.get("upstream_raw_inputs", [])
        ]
        rec["command_series"] = current_commands(rec.get("command_series", []))
        if fp.startswith(PLATES_DIR):
            plate_provenance(
                rec,
                r.get("methods_path", ""),
                plate_record.get("inputs", {}).get(fp.removeprefix(PLATES_DIR), {}),
            )
        if fp.startswith(CELLULOSE_FIGURES):
            cellulose_figure_provenance(rec, cellulose_figures["inputs"])
        if not fp.startswith((PLATES_DIR, *CELLULOSE_FIGURES)):
            # Today's hashes, in one format, and the inputs that changed since they were
            # recorded (#637)
            rec["input_data_sha256s"], changed = refreshed_hashes(
                repo, rec.get("input_data_files", []), rec.get("input_data_sha256s", [])
            )
            changed_inputs.extend({"figure": fp, "input": rel} for rel in changed)
        new_figures_list.append(rec)

    prov_doc["figures"] = new_figures_list
    for stat_rec in prov_doc.get("statistical_results", []):
        if "figure_ids" in stat_rec and isinstance(stat_rec["figure_ids"], list):
            stat_rec["figure_ids"] = [old_id_to_new_id[fid] for fid in stat_rec["figure_ids"]]

    provenance_dir.mkdir(parents=True, exist_ok=True)
    with prov_yaml_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(prov_doc, f, sort_keys=False, allow_unicode=True, width=120)

    validation = {
        "version": "v6",
        "commit": COMMIT,
        "branch": BRANCH,
        "total_slides_including_cover": total_slides,
        "figure_slides": len(manifest_rows),
        "previous_unique_slides_retained": n_prev,
        "new_or_missed_slides_added": n_new,
        "main_section_slides": sum(1 for r in manifest_rows if r["deck_section"] == "Main"),
        "appendix_section_slides": sum(1 for r in manifest_rows if r["deck_section"] == "Appendix"),
        # What was checked, rather than claims set to True (#637)
        "review_build": bool(review),
        "review_reasons": review,
        "checked_against_revisions": trees is not None,
        "linked_figures_matching_their_revision": n_checked,
        "files_differing_from_their_revision": problems,
        "inputs_changed_since_recorded": changed_inputs,
        "all_slides_have_vector_svg": all(
            r["asset"].startswith("assets_svg/") for r in manifest_rows
        ),
        "provenance_yaml_written": prov_yaml_path.exists(),
        "outputs": {
            "manifest_v6": str(slides_dir / "slide_manifest_v6.tsv"),
            "manifest_default": str(slides_dir / "slide_manifest.tsv"),
            "html_v6": str(slides_dir / "selected_figures_v6.html"),
            "html_default": str(slides_dir / "selected_figures.html"),
            "provenance_yaml": str(prov_yaml_path),
        },
    }
    (slides_dir / "deck_validation_v6.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(validation, indent=2))


def main(argv: Sequence[str] | None = None) -> int:
    # Whole option names only, so the recipe's guards see every option (#617).
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0] if __doc__ else None, allow_abbrev=False
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=REPO_ROOT,
        help="Checkout holding the figures, tables and earlier deck manifests (default: this one).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / DECK_DIR / "slides",
        help="Where the deck, manifest and assets go (default: this checkout's "
        f"{DECK_DIR}/slides, whatever --input-root is).",
    )
    parser.add_argument(
        "--provenance-dir",
        type=Path,
        default=REPO_ROOT / DECK_DIR / "provenance",
        help=f"Where the updated {PROVENANCE_NAME} goes; it is read from the input root "
        f"(default: this checkout's {DECK_DIR}/provenance).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing deck. The deck is a gitignored output with no other copy, so "
        "back up the slides folder first.",
    )
    parser.add_argument(
        "--allow-contrast-version-mismatch",
        action="store_true",
        help=f"Build from contrast-orientation inputs other than {TYPED_CONTRAST_VERSION}, "
        "whose numbers then disagree with the typed text (review builds only; #579).",
    )
    parser.add_argument(
        "--git-trees",
        type=Path,
        help="Folder of `git ls-tree -r -z --full-tree <ref>` listings named <ref>.tree, one per "
        "revision --print-refs names; `just figure-review-deck` writes it (#637).",
    )
    parser.add_argument(
        "--allow-unverified",
        action="store_true",
        help="Build even when a committed file a slide links or quotes differs from the "
        "revision it cites, or without --git-trees; such files are not linked and the deck is "
        "marked as a review build (#637).",
    )
    parser.add_argument(
        "--print-refs",
        action="store_true",
        help="Print the revisions whose files the slides link or quote, one per line, and exit.",
    )
    args = parser.parse_args(argv)
    if args.print_refs:
        print("\n".join(cited_refs()))
        return 0
    repo = args.input_root.resolve()
    mismatch = check_contrast_version(repo, args.allow_contrast_version_mismatch)
    existing = [
        path
        for path in (
            args.output_dir / "selected_figures_v6.html",
            args.provenance_dir / PROVENANCE_NAME,
        )
        if path.exists()
    ]
    if existing and not args.overwrite:
        raise SystemExit(
            f"{existing[0]} exists; back it up and pass --overwrite to replace the deck"
        )
    build_deck(
        repo,
        args.output_dir,
        args.provenance_dir,
        load_git_trees(args.git_trees) if args.git_trees else None,
        args.allow_unverified,
        [mismatch] if mismatch else [],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
