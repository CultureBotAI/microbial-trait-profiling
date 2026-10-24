"""Recount the cellulose_decr numbers that #799 re-derived for the docs (2026-10-23).

Purpose: recompute, from the committed cellulose_decr tables, every number that #799
changed in data/analysis/cellulose_decr/README.md, SIGNIFICANT_RESULTS_SUMMARY.md and
FOUR_HIGHLIGHTED_PAIRS_LITERATURE_SUPPORT.md, format it as the docs do, and check that
the docs quote it (#799, #1004).

Input: --analysis-dir, default data/analysis/cellulose_decr, as regenerated in the final
rerun at code commit 511c891. Expected sha256 prefixes: cellulose_decr_summary.json
75ac4a00, significant_results_summary.tsv 788cbaeb,
compartment_jaccard_and_direction_stats.tsv 4730c140,
compartment_direction_interaction.tsv a9177930, condition_definition_interaction.tsv
ec078ce6, condition_definition_sensitivity.tsv 7f53c297,
sensitivity_definitions_jaccard.tsv 44092c7d, cellulose_decreased_taxa_catalog.tsv
1e326706, profile_and_study_pair_jaccard.tsv 88bd83a4,
four_highlighted_pairs_summary.json 0c776be8.

Output: one tab-separated line per check on stdout (label, quoted text, source, OK or
MISSING, the docs that quote it) and a count. Exits 1 if any check is missing.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections.abc import Sequence
from pathlib import Path

EN_DASH = chr(0x2013)  # the docs write ranges with an en dash
LE = chr(0x2264)  # less-than-or-equal sign


def p2(value: str | float) -> str:
    """Two significant figures in the README style: 0.0086, 3.2e-4, 0.10."""
    x = float(value)
    if x >= 1e-3:
        text = f"{x:.2g}"
        if "e" in text:
            text = f"{x:.4f}"
        if len(text.replace("0.", "").lstrip("0")) < 2:
            text += "0"
        return text
    mantissa, exponent = f"{x:.1e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class Checker:
    def __init__(self, docs: dict[str, str]) -> None:
        self.docs = docs
        self.rows: list[tuple[str, str, str, str, str]] = []

    def check(self, label: str, value: str, names: Sequence[str], source: str) -> None:
        found = [name for name in names if value in self.docs[name]]
        status = "OK" if found else "MISSING"
        self.rows.append((label, value, source, status, ",".join(found)))


def run(analysis: Path) -> int:
    docs = {
        "README": (analysis / "README.md").read_text(encoding="utf-8"),
        "SUMMARY": (analysis / "SIGNIFICANT_RESULTS_SUMMARY.md").read_text(encoding="utf-8"),
        "PAIRS": (analysis / "FOUR_HIGHLIGHTED_PAIRS_LITERATURE_SUPPORT.md").read_text(
            encoding="utf-8"
        ),
    }
    c = Checker(docs)
    rs, r_only = ("README", "SUMMARY"), ("README",)
    summary = json.loads((analysis / "cellulose_decr_summary.json").read_text(encoding="utf-8"))

    comp = {
        (r["contrast_population"], r["compartment"], r["short_condition"]): r
        for r in read_tsv(analysis / "compartment_jaccard_and_direction_stats.tsv")
    }
    populations = (("all_contrasts", "all"), ("disease_vs_healthy_oriented", "disease vs healthy"))
    conditions = (("PD", "PD"), ("HIV", "HIV"), ("CRC", "CRC"), ("Pooled", "All three"))
    for pop, pop_label in populations:
        for cond, cond_label in conditions:
            r = comp[(pop, "Feces", cond)]
            row = (
                f"| {cond_label} | {pop_label} | {r['paired_exp_down_dec_gt_inc']} : "
                f"{r['paired_exp_up_inc_gt_dec']} : {r['paired_exp_tied']} | "
                f"{p2(r['binom_p_two_sided'])} | {r['study_down']} : {r['study_up']} : "
                f"{r['study_tied']} | {p2(r['study_binom_p_two_sided'])} |"
            )
            source = "compartment_jaccard_and_direction_stats.tsv"
            c.check(f"section 4 {cond_label} {pop_label}", row, rs, source)
            wilcoxon = f"| {p2(r['study_wilcoxon_p_two_sided'])} |"
            c.check(f"section 4 {cond_label} {pop_label} Wilcoxon", wilcoxon, rs, source)

    sig = {r["test_id"]: r for r in read_tsv(analysis / "significant_results_summary.tsv")}
    for cond in ("PD", "HIV", "CRC", "Pooled"):
        for pop in ("disease_vs_healthy_oriented", "all_contrasts"):
            test_id = f"feces_{cond}_study_direction_{pop}"
            q = p2(sig[test_id]["bh_q_within_table"])
            c.check(f"q {test_id}", q, rs, "significant_results_summary.tsv: bh_q_within_table")

    nulls = summary["null_models"]
    for key in ("U4", "U5", "U6", "U6_size_matched", "within"):
        mean = nulls[key]["mean_pairwise_jaccard_distance"]
        text = f"{mean['expected']:.3f} ({mean['p2_5']:.3f}{EN_DASH}{mean['p97_5']:.3f})"
        c.check(f"{key} mean d_J", text, r_only, "summary.json: null_models")
        for endpoint in ("mean_pairwise_jaccard_distance", "three_way_core"):
            p = nulls[key][endpoint]["p"]
            text = p2(p) if p > 1e-5 else f"{LE} 1e-5"
            c.check(f"{key} {endpoint} p", text, rs, "summary.json: null_models")
    iterations = nulls["U6_size_matched"]["iterations"]
    c.check("U6 size-matched draws", f"{iterations:,}", rs, "summary.json")
    universes = summary["universes"]
    union = universes["all_feces_decreased_taxa_PD_HIV_CRC_union"]
    c.check("U4 universe", f"on the {union} decreased", r_only, "summary.json: universes")
    units = universes["U6_units_per_condition"]
    text = f"draws {units['PD']}, {units['HIV']} and {units['CRC']} whole"
    c.check("U6 units", text, r_only, "summary.json: universes")
    c.check("U6 studies", f"from all {universes['U6_studies']} feces", r_only, "summary.json")

    total = summary["trait_benchmark_total_traits_tested"]
    biological = summary["trait_benchmark_non_growth_medium_traits_tested"]
    median = summary["benchmark_median_mean_pairwise_jaccard_distance_excluding_growth_media"]
    c.check("benchmark traits", f"{total} matrix traits have", r_only, "summary.json")
    c.check("benchmark non-MediaDive", f"second of the {biological}", r_only, "summary.json")
    c.check("benchmark median", f"d_J {median:.3f}", r_only, "summary.json")
    c.check("summary benchmark", f"of {biological} non-MediaDive", ("SUMMARY",), "summary.json")
    for rank, trait in enumerate(summary["benchmark_top_10_excluding_growth_media"], 1):
        source = "summary.json: benchmark_top_10_excluding_growth_media"
        c.check(f"top 10 #{rank}", f"| {rank} | {trait['trait_label']} |", r_only, source)
        text = (
            f"| {trait['mean_pairwise_jaccard_distance']:.3f} | "
            f"{trait['three_way_intersection']} / {trait['three_way_union']} |"
        )
        c.check(f"top 10 #{rank} d_J", text, r_only, source)

    gi = "Non-fecal GI (mucosa/lumen)"
    inter = {
        (r["contrast_population"], r["short_condition"], r["unit"], r["compartment_b"]): r
        for r in read_tsv(analysis / "compartment_direction_interaction.tsv")
    }
    stratified = "Pooled, stratified by condition"

    def fisher(r: dict[str, str], estimate: str, high: str) -> str:
        return (
            f"| {r['a_down']} : {r['a_up']} | {r['b_down']} : {r['b_up']} | {estimate} "
            f"({float(r['conditional_or_ci95_low']):.3g}{EN_DASH}{high}) | "
            f"{p2(r['fisher_p_two_sided'])} |"
        )

    r = inter[("all_contrasts", stratified, "study", gi)]
    text = fisher(
        r, f"{float(r['conditional_or']):.3g}", f"{float(r['conditional_or_ci95_high']):.3g}"
    )
    c.check(
        "Fisher, all three stratified, studies", text, rs, "compartment_direction_interaction.tsv"
    )
    r = inter[("disease_vs_healthy_oriented", "CRC", "study", gi)]
    text = fisher(
        r, f"{float(r['conditional_or']):.3g}", f"{float(r['conditional_or_ci95_high']):.0f}"
    )
    c.check(
        "Fisher, CRC disease vs healthy, studies", text, rs, "compartment_direction_interaction.tsv"
    )
    r = inter[("disease_vs_healthy_oriented", stratified, "study", gi)]
    text = fisher(
        r, f"{float(r['conditional_or']):.2f}", f"{float(r['conditional_or_ci95_high']):.3g}"
    )
    c.check(
        "Fisher, all three, disease vs healthy",
        text,
        r_only,
        "compartment_direction_interaction.tsv",
    )

    definitions = {
        (
            r["scenario"],
            r["contrast_population"],
            r["short_condition"],
            r["unit"],
            r["compartment_b"],
        ): r
        for r in read_tsv(analysis / "condition_definition_interaction.tsv")
    }
    oriented = "disease_vs_healthy_oriented"
    r = definitions[("paired_tissue_counted", oriented, "CRC", "study", gi)]
    c.check(
        "paired tissue p", p2(r["fisher_p_two_sided"]), rs, "condition_definition_interaction.tsv"
    )
    low, high = float(r["conditional_or_ci95_low"]), float(r["conditional_or_ci95_high"])
    ci = f"({low:.3g}{EN_DASH}{high:.0f}"
    c.check("paired tissue CI", ci, rs, "condition_definition_interaction.tsv")
    r = definitions[("control_tiers_1_2", oriented, "CRC", "study", gi)]
    c.check(
        "Tier 1-2 GI p", p2(r["fisher_p_two_sided"]), rs, "condition_definition_interaction.tsv"
    )
    text = (
        f"{float(r['conditional_or']):.3g} ({float(r['conditional_or_ci95_low']):.3g}{EN_DASH}"
        f"{float(r['conditional_or_ci95_high']):.0f}"
    )
    c.check("Tier 1-2 GI odds ratio", text, ("SUMMARY",), "condition_definition_interaction.tsv")

    sensitivity = {
        (r["scenario"], r["contrast_population"], r["short_condition"], r["unit"]): r
        for r in read_tsv(analysis / "condition_definition_sensitivity.tsv")
    }
    r = sensitivity[("synonyms_merged", oriented, "CRC", "study")]
    text = f"{r['down']} : {r['up']} : {r['tied']} studies (p = {p2(r['binom_p_two_sided'])})"
    c.check(
        "CRC synonyms, disease vs healthy", text, r_only, "condition_definition_sensitivity.tsv"
    )
    pooled = [
        float(v["binom_p_two_sided"])
        for k, v in sensitivity.items()
        if k[2:] == ("Pooled", "study")
    ]
    mantissa, exponent = f"{max(pooled):.2e}".split("e")
    bound = f"p {LE} {mantissa}e{int(exponent)}"
    c.check("pooled sign test, worst scenario", bound, rs, "condition_definition_sensitivity.tsv")
    r = sensitivity[("control_tiers_1_2", oriented, "CRC", "study")]
    untied = int(r["down"]) + int(r["up"])
    text = f"{r['down']} of {untied}\n     (p = {p2(r['binom_p_two_sided'])})"
    c.check("CRC Tier 1-2 studies", text, ("SUMMARY",), "condition_definition_sensitivity.tsv")

    jaccard = {
        r["definition_slug"]: r
        for r in read_tsv(analysis / "sensitivity_definitions_jaccard.tsv")
        if r["comparison"].startswith("Mean")
    }
    r = jaccard["6_baseline_all_feces_dec_taxa"]
    text = (
        f"| {r['n_a']} / {r['n_b']} / {r['n_c']} | {float(r['jaccard_distance']):.3f} | "
        f"{r['intersection']} / {r['union']} | {p2(r['matrix_hypergeom_p'])} |"
    )
    c.check("every decreased feces taxon", text, r_only, "sensitivity_definitions_jaccard.tsv")

    catalog = {
        r["taxon_name"]: r for r in read_tsv(analysis / "cellulose_decreased_taxa_catalog.tsv")
    }
    b = catalog["Bacteroidia"]
    text = (
        f"| Bacteroidia | class | {b['PD_dec_studies']} / {b['PD_inc_studies']} | "
        f"{b['HIV_dec_studies']} / {b['HIV_inc_studies']} | "
        f"{b['CRC_dec_studies']} / {b['CRC_inc_studies']} |"
    )
    c.check("Bacteroidia", text, r_only, "cellulose_decreased_taxa_catalog.tsv")

    pairs = {r["comparison"]: r for r in read_tsv(analysis / "profile_and_study_pair_jaccard.tsv")}
    text = (
        f"| {pairs['PD vs. HIV']['n_studies_a']} / {pairs['PD vs. HIV']['n_studies_b']} / "
        f"{pairs['PD vs. CRC']['n_studies_b']} |"
    )
    c.check(
        "section 8, studies with an annotated taxon",
        text,
        r_only,
        "profile_and_study_pair_jaccard.tsv",
    )
    decreased = [comp[("all_contrasts", "Feces", k)]["n_studies_dec"] for k in ("PD", "HIV", "CRC")]
    text = f"| {' / '.join(decreased)} |"
    c.check("section 8, studies with a decreased signature", text, r_only, "n_studies_dec")

    four = json.loads(
        (analysis / "four_highlighted_pairs_summary.json").read_text(encoding="utf-8")
    )
    for key, label in (
        ("1_PD_Feces_Cellulose_dec", "1. PD, cellulose, decreased"),
        ("4_PD_Feces_ModeratelyHalophilic_inc", "4. PD, moderately halophilic, increased"),
    ):
        text = f"| {label} | {four[key]['n_supporting_studies_in_target_dir']} |"
        c.check(f"four pairs, {label}", text, ("PAIRS",), "four_highlighted_pairs_summary.json")

    for row in c.rows:
        print("\t".join(row))
    missing = [row for row in c.rows if row[3] != "OK"]
    print(f"{len(c.rows)} checks, {len(missing)} missing")
    return 1 if missing else 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=Path("data/analysis/cellulose_decr"))
    args = parser.parse_args()
    sys.exit(run(args.analysis_dir))


if __name__ == "__main__":
    main()
