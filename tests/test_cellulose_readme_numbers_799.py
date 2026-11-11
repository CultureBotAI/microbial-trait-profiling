"""Bind every README result to the committed cellulose tables (#799).

This is a document contract, not a rerun of the statistical analysis. Exact block
identity and complete numeric coverage prevent a value elsewhere in the README
from satisfying the wrong claim. Sources are resolved beside the copied README.
"""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Any

import pytest

DATA = Path(__file__).resolve().parents[1] / "data" / "analysis" / "cellulose_decr"
README = DATA / "README.md"
ORIENTED = "disease_vs_healthy_oriented"
ALL = "all_contrasts"
GI = "Non-fecal GI (mucosa/lumen)"
CONDITIONS = ("PD", "HIV", "CRC")
CELLULOSE = "METPO:2000103|GO:0030245"
SOURCES = (
    "cellulose_decr_summary.json",
    "cellulose_decr_provenance.json",
    "pairwise_jaccard_significance.tsv",
    "cellulose_decreased_taxa_catalog.tsv",
    "cellulose_decreased_taxa_catalog_tier1_2.tsv",
    "compartment_jaccard_and_direction_stats.tsv",
    "significant_results_summary.tsv",
    "condition_definition_sensitivity.tsv",
    "condition_definition_overlap.tsv",
    "compartment_direction_interaction.tsv",
    "condition_definition_interaction.tsv",
    "sensitivity_definitions_jaccard.tsv",
    "trait_convergence_benchmark.tsv",
    "profile_and_study_pair_jaccard.tsv",
    "label_source_corrections.tsv",
)


@dataclass(frozen=True)
class Binding:
    section: str
    text: str


def _rows(name: str, root: Path = DATA) -> list[dict[str, str]]:
    assert name in SOURCES, f"Undeclared quantitative source: {name}"
    with (root / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _one(rows: list[dict[str, str]], **key: str) -> dict[str, str]:
    found = [row for row in rows if all(row[field] == value for field, value in key.items())]
    assert len(found) == 1, f"Source row must be unique: {key}; found {len(found)}"
    return found[0]


def _p(value: str | float) -> str:
    number = float(value)
    if number >= 0.001:
        rendered = f"{number:.2g}"
        if "e" in rendered:
            rendered = f"{number:.4f}"
        if len(rendered.removeprefix("0.").lstrip("0")) < 2:
            rendered += "0" if "." in rendered else ".0"
        return rendered
    mantissa, exponent = f"{number:.1e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def _g(value: str | float) -> str:
    return f"{float(value):.3g}"


def _triple(row: dict[str, Any], names: tuple[str, str, str]) -> str:
    return " : ".join(str(row[name]) for name in names)


def _table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    lines = [headers, tuple("---" for _ in headers), *rows]
    return "\n".join("| " + " | ".join(row) + " |" for row in lines)


def _expected_blocks(root: Path = DATA) -> dict[str, Binding]:
    summary = json.loads((root / "cellulose_decr_summary.json").read_text(encoding="utf-8"))
    provenance = json.loads((root / "cellulose_decr_provenance.json").read_text(encoding="utf-8"))
    blocks: dict[str, Binding] = {}

    def add(key: str, section: str, text: str) -> None:
        assert key not in blocks, f"Duplicate quantitative binding: {key}"
        blocks[key] = Binding(section, text)

    add(
        "producer",
        "Introduction",
        f"The main analysis tables record code commit `{provenance['code_commit'][:7]}` "
        f"and classifier `{provenance['contrast_classifier_version']}`. "
        f"The provenance pins {len(provenance['pinned_inputs'])} input files.",
    )

    history_rows = []
    for condition in CONDITIONS:
        previous = summary["label_source_comparison"]["study_modal"]["feces_direction"][condition]
        current = summary["label_source_comparison"]["experiment"]["feces_direction"][condition]
        names = ("paired_down", "paired_up", "paired_tied")
        history_rows.append((condition, _triple(previous, names), _triple(current, names)))
    add(
        "label_mode_comparison",
        "What changed on 2026-10-17 (#552)",
        "The committed rerun compares study-modal and experiment labels on its own inputs; "
        "these are not preserved results from the earlier dated runs.\n\n"
        + _table(
            ("Condition", "Study-modal down : up : tied", "Experiment down : up : tied"),
            history_rows,
        )
        + "\n\nThe label-source corrections table records "
        f"{len(_rows('label_source_corrections.tsv', root))} "
        "experiments in the target conditions whose labels differ.",
    )

    pairwise = _rows("pairwise_jaccard_significance.tsv", root)
    pairs = [row for row in pairwise if row["endpoint"] == "pairwise_jaccard_distance"]
    core = _one(pairwise, endpoint="three_way_core")
    mean = _one(pairwise, endpoint="mean_pairwise_jaccard_distance")
    overlap_rows = [
        (
            row["comparison"],
            f"{row['n_a']}, {row['n_b']}",
            row["intersection"],
            row["union"],
            f"{float(row['jaccard_similarity']):.3f}",
            f"{float(row['jaccard_distance']):.3f}",
        )
        for row in pairs
    ]
    overlap_rows.extend(
        [
            (
                "Mean of the pairs",
                "",
                "",
                "",
                f"{float(mean['jaccard_similarity']):.3f}",
                f"{float(mean['jaccard_distance']):.3f}",
            ),
            ("All conditions", "", core["intersection"], core["union"], "", ""),
        ]
    )
    add(
        "overlap_table",
        "1. Overlap of the taxon sets",
        _table(("Pair", "Set sizes", "Shared", "Union", "J", "d_J"), overlap_rows),
    )
    catalog = _rows("cellulose_decreased_taxa_catalog.tsv", root)
    core_taxa = [row for row in catalog if row["is_3way_core"] == "True"]
    add(
        "overlap_support",
        "1. Overlap of the taxon sets",
        f"{sum(int(row['n_conditions_feces_dec']) >= 2 for row in catalog)} of the "
        f"{len(catalog)} catalog taxa occur in at least 2 condition sets.",
    )

    nulls = summary["null_models"]
    universes = summary["universes"]
    draws = summary["monte_carlo"]["permutations"]
    cutoff = 1 / draws

    def mc_p(value: float) -> str:
        return f"≤ {cutoff:.0e}".replace("e-0", "e-") if value <= cutoff else _p(value)

    add(
        "null_draws",
        "2. Is the overlap unusual?",
        f"Monte Carlo nulls use {draws:,} draws. The U6 size-window subset retains "
        f"{nulls['U6_size_matched']['iterations']:,} draws under the recorded tolerance: "
        f"`{nulls['U6_size_matched']['tolerance']}`. "
        "Its set sizes need not equal the observed sizes. "
        f"The displayed central intervals use the saved 2.5th and 97.5th percentiles. "
        f"Values at the simulation floor are displayed as {mc_p(cutoff)}, a rounded "
        "display convention rather than the exact finite-simulation minimum.",
    )
    labels = (
        (
            "U1",
            f"U1: uniform, {universes['U1_all_cellulose_taxa_in_matrix']} annotated matrix taxa",
        ),
        (
            "U2",
            f"U2: uniform, {universes['U2_all_cellulose_taxa_in_retained_feces_profiles']} "
            "annotated taxa seen in feces",
        ),
        (
            "U3",
            f"U3: uniform, {universes['U3_all_cellulose_taxa_in_feces_decreased_profiles']} "
            "annotated taxa ever decreased in feces",
        ),
        (
            "U4",
            f"U4: random trait on {universes['all_feces_decreased_taxa_PD_HIV_CRC_union']} "
            "decreased feces taxa",
        ),
        ("U5", "U5: reporting-weighted U3 taxa"),
        ("U6", "U6: whole study signatures"),
        ("U6_size_matched", "U6: retained size-window draws"),
        ("within", "Within-condition reference, any trait"),
    )
    null_rows = []
    for key, label in labels:
        distance = nulls[key]["mean_pairwise_jaccard_distance"]
        ncore = nulls[key]["three_way_core"]
        distance_text = (
            f"{distance['expected']:.3f} ({distance['p2_5']:.3f}\u2013{distance['p97_5']:.3f})"
        )
        if key in ("U1", "U2", "U3"):
            prefix = key.lower()
            core_text = f"{float(core[prefix + '_expected_intersection']):.2f}"
            probability = "exact " + _p(core[prefix + "_hypergeom_p"])
        else:
            core_text = f"{ncore['expected']:.3g} ({ncore['p2_5']:.0f}\u2013{ncore['p97_5']:.0f})"
            probability = mc_p(ncore["p"])
        null_rows.append((label, distance_text, mc_p(distance["p"]), core_text, probability))
    add(
        "null_table",
        "2. Is the overlap unusual?",
        _table(
            (
                "Null",
                "Expected mean d_J (central interval)",
                "p",
                "Expected core (central interval)",
                "p",
            ),
            null_rows,
        ),
    )
    units = universes["U6_units_per_condition"]
    add(
        "null_units",
        "2. Is the overlap unusual?",
        "U6 samples whole signatures: "
        + ", ".join(f"{name} {units[name]}" for name in CONDITIONS)
        + f", from {universes['U6_studies']} feces studies with a decreased signature. "
        + "Its pairwise p-values are "
        + ", ".join(f"{row['comparison']}: {_p(row['u6_signature_perm_p'])}" for row in pairs)
        + ".",
    )
    scenarios = _rows("condition_definition_overlap.tsv", root)
    scenario_rows = []
    for scenario in ("synonyms_merged", "without_site_miscurations"):
        row = _one(scenarios, scenario=scenario, contrast_population=ALL)
        scenario_rows.append(
            (
                row["scenario_label"],
                f"{row['n_PD']} / {row['n_HIV']} / {row['n_CRC']}",
                f"{row['three_way_core']} / {row['three_way_union']}",
                _p(row["u6_mean_jaccard_distance_p"]),
                _p(row["u6_core_p"]),
            )
        )
    add(
        "null_scenarios",
        "2. Is the overlap unusual?",
        _table(
            (
                "All-contrasts condition definition",
                "PD / HIV / CRC",
                "Core / union",
                "U6 mean p",
                "U6 core p",
            ),
            scenario_rows,
        ),
    )

    ordered_core = sorted(
        core_taxa,
        key=lambda row: (
            row["taxonomic_rank"] != "genus",
            -int(row["PD_dec_studies"]),
            row["taxon_name"],
        ),
    )
    catalog_rows = [
        (
            row["taxon_name"],
            row["taxonomic_rank"],
            *(
                f"{row[name + '_dec_studies']} / {row[name + '_inc_studies']}"
                for name in CONDITIONS
            ),
            "yes" if row["is_strictly_net_decreased_all3"] == "True" else "no",
        )
        for row in ordered_core
    ]
    add(
        "core_table",
        "3. Core taxa",
        _table(("Taxon", "Rank", "PD", "HIV", "CRC", "Net decreased"), catalog_rows),
    )
    genera = sum(row["taxonomic_rank"] == "genus" for row in core_taxa)
    direct = sum(row["has_direct_trait_edge"] == "True" for row in core_taxa)
    rank_edges = summary["direct_trait_edges_by_rank"]
    direct_union = [row for row in catalog if row["has_direct_trait_edge"] == "True"]
    add(
        "core_annotations",
        "3. Core taxa",
        f"The core contains {len(core_taxa)} taxa: {genera} genera and {len(core_taxa) - genera} "
        f"family-to-phylum nodes. {direct} core taxa have direct annotation edges. "
        f"The KG has {summary['kg_direct_trait_subjects']} NCBITaxon subjects carrying the edge. "
        f"Among matrix trait-bearing taxa, {rank_edges['genus']['with_direct_edge']} of "
        f"{rank_edges['genus']['trait_bearing_matrix_taxa']} genera and "
        f"{rank_edges['species']['with_direct_edge']} of "
        f"{rank_edges['species']['trait_bearing_matrix_taxa']} species have direct edges. "
        f"The catalog union contains {len(direct_union)} taxa with direct edges, "
        f"{sum(row['taxonomic_rank'] == 'species' for row in direct_union)} of them species.",
    )
    net_study = [row for row in core_taxa if row["is_strictly_net_decreased_all3"] == "True"]
    net_profile = [row for row in core_taxa if row["is_net_decreased_all3_by_profiles"] == "True"]
    lineage_counts = "; ".join(
        f"{name}: {_one(catalog, taxon_name=name)['matrix_trait_bearing_taxa_at_or_below']}"
        for name in ("Lachnospiraceae", "Bacillota")
    )
    add(
        "core_net",
        "3. Core taxa",
        f"Annotated matrix taxa at or below the named nodes are {lineage_counts}. "
        f"Net decreased across all conditions: {len(net_study)} taxa by studies "
        f"({sum(row['taxonomic_rank'] == 'genus' for row in net_study)} genera), "
        f"and {len(net_profile)} by profiles.",
    )
    eubacterium = _one(catalog, taxon_name="Eubacterium")
    without_rectal = _one(scenarios, scenario="without_site_miscurations", contrast_population=ALL)
    add(
        "core_rectal",
        "3. Core taxa",
        f"Eubacterium has {eubacterium['HIV_dec_studies']} HIV decreased-study support. "
        f"Excluding the rectal-mucosa study leaves {without_rectal['three_way_core']} core taxa.",
    )

    comp = _rows("compartment_jaccard_and_direction_stats.tsv", root)
    exp_names = ("paired_exp_down_dec_gt_inc", "paired_exp_up_inc_gt_dec", "paired_exp_tied")
    study_names = ("study_down", "study_up", "study_tied")
    direction_rows = []
    for condition in (*CONDITIONS, "Pooled"):
        for population, label in ((ALL, "all"), (ORIENTED, "disease vs healthy")):
            row = _one(
                comp, short_condition=condition, compartment="Feces", contrast_population=population
            )
            direction_rows.append(
                (
                    condition,
                    label,
                    _triple(row, exp_names),
                    _p(row["binom_p_two_sided"]),
                    _triple(row, study_names),
                    _p(row["study_binom_p_two_sided"]),
                    _p(row["study_wilcoxon_p_two_sided"]),
                )
            )
    add(
        "direction_table",
        "4. Direction within each condition",
        _table(
            (
                "Condition",
                "Contrasts",
                "Experiments down : up : tied",
                "Sign test p",
                "Studies down : up : tied",
                "Sign test p",
                "Study Wilcoxon p",
            ),
            direction_rows,
        ),
    )
    significant = _rows("significant_results_summary.tsv", root)
    q_rows = []
    for condition in (*CONDITIONS, "Pooled"):
        values = [
            _p(
                _one(significant, test_id=f"feces_{condition}_study_direction_{population}")[
                    "bh_q_within_table"
                ]
            )
            for population in (ALL, ORIENTED)
        ]
        q_rows.append((condition, *values))
    add(
        "direction_fdr",
        "4. Direction within each condition",
        f"Benjamini\u2013Hochberg correction covers the {len(significant)} rows in "
        "`significant_results_summary.tsv`; these are the study-direction q-values.\n\n"
        + _table(("Condition", "All contrasts q", "Disease vs healthy q"), q_rows),
    )
    majority_rows = []
    for condition, population in (("PD", ALL), ("PD", ORIENTED), ("CRC", ORIENTED)):
        row = _one(
            comp, short_condition=condition, compartment="Feces", contrast_population=population
        )
        majority_rows.append(
            (
                condition,
                population,
                _triple(row, ("study_majority_down", "study_majority_up", "study_majority_tied")),
                _p(row["study_majority_binom_p_two_sided"]),
            )
        )
    add(
        "direction_majority",
        "4. Direction within each condition",
        _table(("Condition", "Population", "Study majority down : up : tied", "p"), majority_rows),
    )
    filtered_rows = []
    for population in (
        "disease_vs_healthy_strong_evidence",
        "disease_vs_healthy_oriented_excl_msm_mention",
    ):
        row = _one(comp, contrast_population=population, short_condition="HIV", compartment="Feces")
        filtered_rows.append(
            (
                population,
                _triple(row, exp_names),
                _p(row["binom_p_two_sided"]),
                _triple(row, study_names),
                _p(row["study_binom_p_two_sided"]),
            )
        )
    add(
        "direction_additional_filters",
        "4. Direction within each condition",
        _table(("HIV population", "Experiments", "p", "Studies", "p"), filtered_rows),
    )
    definitions = _rows("condition_definition_sensitivity.tsv", root)
    hiv_rows = []
    hiv_scenarios = list(
        dict.fromkeys(
            row["scenario"] for row in definitions if row["contrast_population"] == ORIENTED
        )
    )
    for scenario in hiv_scenarios:
        exp = _one(
            definitions,
            scenario=scenario,
            contrast_population=ORIENTED,
            short_condition="HIV",
            unit="experiment",
        )
        study = _one(
            definitions,
            scenario=scenario,
            contrast_population=ORIENTED,
            short_condition="HIV",
            unit="study",
        )
        hiv_rows.append(
            (
                exp["scenario_label"],
                _triple(exp, ("down", "up", "tied")),
                _p(exp["binom_p_two_sided"]),
                _triple(study, ("down", "up", "tied")),
                _p(study["binom_p_two_sided"]),
            )
        )
    add(
        "direction_hiv",
        "4. Direction within each condition",
        _table(("HIV disease-vs-healthy definition", "Experiments", "p", "Studies", "p"), hiv_rows),
    )
    synonym_rows = []
    for condition, population in (("HIV", ALL), ("CRC", ALL), ("CRC", ORIENTED)):
        exp = _one(
            definitions,
            scenario="synonyms_merged",
            contrast_population=population,
            short_condition=condition,
            unit="experiment",
        )
        study = _one(
            definitions,
            scenario="synonyms_merged",
            contrast_population=population,
            short_condition=condition,
            unit="study",
        )
        synonym_rows.append(
            (
                condition,
                population,
                _triple(exp, ("down", "up", "tied")),
                _p(exp["binom_p_two_sided"]),
                _triple(study, ("down", "up", "tied")),
                _p(study["binom_p_two_sided"]),
            )
        )
    pooled = [
        float(row["binom_p_two_sided"])
        for row in definitions
        if row["short_condition"] == "Pooled" and row["unit"] == "study"
    ]
    add(
        "direction_synonyms",
        "4. Direction within each condition",
        _table(("Synonym labels", "Population", "Experiments", "p", "Studies", "p"), synonym_rows)
        + f"\n\nThe largest pooled study sign-test p across the committed condition-definition "
        f"rows is {_g(max(pooled))}.",
    )

    compartment_rows = []
    compartments = list(
        dict.fromkeys(row["compartment"] for row in comp if row["contrast_population"] == ALL)
    )
    for compartment in compartments:
        studies = [
            _one(comp, contrast_population=ALL, compartment=compartment, short_condition=name)[
                "n_studies_dec"
            ]
            for name in CONDITIONS
        ]
        row = _one(comp, contrast_population=ALL, compartment=compartment, short_condition="Pooled")
        compartment_rows.append(
            (compartment, " / ".join(studies), f"{row[exp_names[0]]} : {row[exp_names[1]]}")
        )
    add(
        "compartment_table",
        "5. Outside feces",
        _table(
            (
                "Compartment",
                "Decreased-signature studies PD / HIV / CRC",
                "Pooled experiments down : up",
            ),
            compartment_rows,
        ),
    )
    outside = [
        row for row in comp if row["contrast_population"] == ALL and row["compartment"] != "Feces"
    ]
    crc_gi = _one(outside, short_condition="CRC", compartment=GI)
    sparse = [row for row in outside if row["short_condition"] != "Pooled" and row != crc_gi]
    max_sparse = max(int(row[exp_names[0]]) + int(row[exp_names[1]]) for row in sparse)
    pooled_counts = [
        int(row[exp_names[0]]) + int(row[exp_names[1]])
        for row in outside
        if row["short_condition"] == "Pooled"
    ]
    add(
        "compartment_power",
        "5. Outside feces",
        f"Outside feces, individual-condition strata other than CRC GI have at most "
        f"{max_sparse} untied experiments. CRC GI has "
        f"{int(crc_gi[exp_names[0]]) + int(crc_gi[exp_names[1]])}, with "
        f"{crc_gi[exp_names[0]]} pointing down. The non-fecal pooled strata have "
        + ", ".join(str(value) for value in pooled_counts)
        + " untied experiments, in the table order.",
    )
    interactions = _rows("compartment_direction_interaction.tsv", root)
    interaction_definitions = _rows("condition_definition_interaction.tsv", root)
    selected_interactions = []
    for population, condition, unit in (
        (ALL, "CRC", "experiment"),
        (ALL, "CRC", "study"),
        (ALL, "Pooled, stratified by condition", "experiment"),
        (ALL, "Pooled, stratified by condition", "study"),
        (ORIENTED, "CRC", "experiment"),
        (ORIENTED, "CRC", "study"),
        (ORIENTED, "Pooled, stratified by condition", "study"),
    ):
        row = _one(
            interactions,
            contrast_population=population,
            short_condition=condition,
            unit=unit,
            compartment_a="Feces",
            compartment_b=GI,
        )
        selected_interactions.append((f"{condition}, {population}", unit, row))
    for scenario in ("paired_tissue_counted", "control_tiers_1_2"):
        row = _one(
            interaction_definitions,
            scenario=scenario,
            contrast_population=ORIENTED,
            short_condition="CRC",
            unit="study",
            compartment_a="Feces",
            compartment_b=GI,
        )
        selected_interactions.append((f"CRC, {scenario}, {ORIENTED}", "study", row))
    interaction_rows = [
        (
            label,
            unit,
            f"{row['a_down']} : {row['a_up']}",
            f"{row['b_down']} : {row['b_up']}",
            f"{_g(row['conditional_or'])} "
            f"({_g(row['conditional_or_ci95_low'])}\u2013{_g(row['conditional_or_ci95_high'])})",
            _p(row["fisher_p_two_sided"]),
        )
        for label, unit, row in selected_interactions
    ]
    add(
        "interaction_table",
        "5. Outside feces",
        _table(
            (
                "Comparison",
                "Unit",
                "Feces",
                "GI mucosa / lumen",
                "Conditional OR (exact interval)",
                "p",
            ),
            interaction_rows,
        ),
    )
    crc = _one(
        interactions,
        contrast_population=ORIENTED,
        short_condition="CRC",
        unit="study",
        compartment_a="Feces",
        compartment_b=GI,
    )
    add(
        "interaction_shared_studies",
        "5. Outside feces",
        f"The CRC disease-vs-healthy comparison includes {crc['studies_in_both_compartments']} "
        "studies sampled in both compartments; the tests nevertheless treat units as independent.",
    )

    sensitivity = [
        row
        for row in _rows("sensitivity_definitions_jaccard.tsv", root)
        if row["comparison"].startswith("Mean")
    ]
    sensitivity_rows = [
        (
            row["definition_label"],
            f"{row['n_a']} / {row['n_b']} / {row['n_c']}",
            f"{float(row['jaccard_distance']):.3f}",
            f"{row['intersection']} / {row['union']}",
            _p(row["matrix_hypergeom_p"]),
        )
        for row in sensitivity
    ]
    add(
        "sensitivity_table",
        "6. Sensitivity to the set definition",
        _table(
            ("Definition", "PD / HIV / CRC", "Mean d_J", "Core / union", "U1 exact core p"),
            sensitivity_rows,
        ),
    )
    benchmark = _rows("trait_convergence_benchmark.tsv", root)
    biological = sorted(
        [row for row in benchmark if row["is_growth_medium"] == "False"],
        key=lambda row: int(row["rank_excluding_growth_media"]),
    )
    cellulose = _one(benchmark, trait_id=CELLULOSE)
    add(
        "benchmark_summary",
        "7. Trait benchmark",
        f"The benchmark contains {len(benchmark)} traits, "
        f"{len(biological)} excluding growth media. "
        f"Cellulose ranks {cellulose['rank_excluding_growth_media']} among those traits and "
        f"{cellulose['rank_by_convergence']} overall. The non-media median mean d_J is "
        f"{median(float(row['mean_pairwise_jaccard_distance']) for row in biological):.3f}. "
        f"The smallest per-condition set size among all benchmarked traits is "
        f"{min(int(row[name + '_dec_taxa']) for row in benchmark for name in CONDITIONS)}.",
    )
    benchmark_rows = [
        (
            row["rank_excluding_growth_media"],
            row["trait_label"],
            "`" + row["trait_id"].replace("|", "\\|") + "`",
            f"{float(row['mean_pairwise_jaccard_distance']):.3f}",
            f"{row['three_way_intersection']} / {row['three_way_union']}",
        )
        for row in biological[:10]
    ]
    add(
        "benchmark_table",
        "7. Trait benchmark",
        _table(("Rank", "Trait", "ID", "Mean d_J", "Core / union"), benchmark_rows),
    )
    pair_studies = _rows("profile_and_study_pair_jaccard.tsv", root)
    pd_hiv = _one(pair_studies, comparison="PD vs. HIV")
    pd_crc = _one(pair_studies, comparison="PD vs. CRC")
    decreased = [
        _one(comp, contrast_population=ALL, compartment="Feces", short_condition=name)[
            "n_studies_dec"
        ]
        for name in CONDITIONS
    ]
    roseburia = _one(catalog, taxon_name="Roseburia")
    add(
        "counting_table",
        "8. Counting rules",
        _table(
            ("Count", "Current pinned run"),
            [
                (
                    "Studies with any decreased feces signature PD / HIV / CRC",
                    " / ".join(decreased),
                ),
                (
                    "Studies reporting an annotated taxon decreased PD / HIV / CRC",
                    f"{pd_hiv['n_studies_a']} / {pd_hiv['n_studies_b']} / {pd_crc['n_studies_b']}",
                ),
                (
                    "Roseburia decreased / increased studies PD; HIV; CRC",
                    "; ".join(
                        f"{roseburia[name + '_dec_studies']} / {roseburia[name + '_inc_studies']}"
                        for name in CONDITIONS
                    ),
                ),
            ],
        ),
    )
    tier_catalog = _rows("cellulose_decreased_taxa_catalog_tier1_2.tsv", root)
    tier_core = [row for row in tier_catalog if row["is_3way_core"] == "True"]
    add(
        "tier_catalog",
        "9. Files",
        "The committed Tier 1\u20132 catalog contains "
        f"{len(tier_catalog)} taxa, {len(tier_core)} in the core, and "
        f"{sum(row['has_direct_trait_edge'] == 'True' for row in tier_core)} "
        "core taxa with a direct edge. "
        "These counts bind directly to the catalog and are cross-checked against its "
        "entries in the committed summary.",
    )
    return blocks


BLOCK = re.compile(r"<!-- numbers: ([a-z_]+) -->\s*(.*?)\s*<!-- /numbers -->", re.DOTALL)
NUMBER_WORD = re.compile(
    r"\b(?:zero|none|one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|"
    r"thirteen|fourteen|fifteen|sixteen|seventeen|eighteen|nineteen|twenty|thirty|forty|"
    r"fifty|sixty|seventy|eighty|ninety|hundred|thousand|million|first|second|third)\b",
    re.IGNORECASE,
)


def _normal(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def _sections(text: str) -> dict[str, str]:
    headings = list(re.finditer(r"^## (.+)$", text, flags=re.MULTILINE))
    names = [match.group(1) for match in headings]
    assert len(names) == len(set(names)), "Duplicate README section"
    sections = {"Introduction": text[: headings[0].start()] if headings else text}
    for index, heading in enumerate(headings):
        end = headings[index + 1].start() if index + 1 < len(headings) else len(text)
        sections[heading.group(1)] = text[heading.end() : end]
    return sections


def _read_blocks(text: str) -> dict[str, tuple[str, str]]:
    found: dict[str, tuple[str, str]] = {}
    for section, content in _sections(text).items():
        for match in BLOCK.finditer(content):
            key, value = match.groups()
            assert key not in found, f"Duplicate quantitative block: {key}"
            found[key] = (section, value)
        leftover = BLOCK.sub("", content)
        assert "<!-- numbers:" not in leftover, f"Malformed quantitative block in {section}"
        assert "<!-- /numbers -->" not in leftover, f"Unmatched quantitative block in {section}"
    return found


def _assert_bound_blocks(text: str, root: Path = DATA) -> None:
    actual = _read_blocks(text)
    expected = _expected_blocks(root)
    assert set(actual) == set(expected), (
        f"Quantitative block inventory differs: missing={set(expected) - set(actual)}, "
        f"unknown={set(actual) - set(expected)}"
    )
    for key, binding in expected.items():
        section, value = actual[key]
        assert section == binding.section, f"Wrong section for {key}: {section}"
        assert _normal(value) == _normal(binding.text), f"Quantitative binding differs: {key}"


def _unbound_numbers(text: str) -> list[str]:
    """Census displayed result tokens, with narrow identifier/structure exemptions."""
    leftovers = BLOCK.sub("", text)
    leftovers = re.sub(r"^## [1-9][0-9]*\. ", "## ", leftovers, flags=re.MULTILINE)
    leftovers = re.sub(r"\b2026-\d{2}-\d{2}\b", "DATE", leftovers)
    leftovers = re.sub(
        r"\b(?:METPO|GO|NCBITaxon|PMID|bsdb):\d+(?:-\d+)?\b", "IDENTIFIER", leftovers
    )
    leftovers = re.sub(r"#\d+\b", "ISSUE", leftovers)
    leftovers = re.sub(r"\bU[1-6](?:\u2013U[1-6])?\b", "NULL", leftovers)
    leftovers = re.sub(r"\bTier [1-3](?:[\u2013-][1-3])?\b", "TIER", leftovers)
    leftovers = re.sub(r"\bclassifier v[3-6]\b", "CLASSIFIER", leftovers)
    leftovers = re.sub(
        r"\bsections? [1-9][0-9]*(?:[\u2013-][1-9][0-9]*)?\b",
        "SECTION",
        leftovers,
        flags=re.IGNORECASE,
    )
    leftovers = re.sub(r"\b(?:three|3)-way\b|\b(?:two|one)-sided\b", "METHOD", leftovers)
    # Paths/commands are recognized only as identifier-like tokens with an extension;
    # numeric prose, link labels, table cells and `17` are never discarded wholesale.
    leftovers = re.sub(
        r"(?<![\w.])(?:[A-Za-z_.-][A-Za-z0-9_./-]*\.(?:py|md|tsv|json|svg|png|pdf|html))\b",
        "FILE",
        leftovers,
    )
    leftovers = re.sub(
        r"https://github\.com/CultureBotAI/microbial-trait-profiling/issues/\d+",
        "ISSUE_URL",
        leftovers,
    )
    leftovers = re.sub(r"(?<!\w)fig[1-9](?:/fig[1-9])?(?!\w)", "FIGURE_ID", leftovers)
    leftovers = re.sub(r"`fig[1-9]_[A-Za-z0-9_]+`", "FIGURE_BASENAME", leftovers)
    found = re.findall(
        r"(?<![A-Za-z])[-+]?\d+(?:[.,]\d+)*(?:e[-+]?\d+)?(?:%|st|nd|rd|th)?", leftovers
    )
    return found + NUMBER_WORD.findall(leftovers)


def test_readme_quantitative_blocks_match_committed_sources() -> None:
    _assert_bound_blocks(README.read_text(encoding="utf-8"))


def test_readme_has_no_unbound_quantitative_claims() -> None:
    text = README.read_text(encoding="utf-8")
    _assert_bound_blocks(text)
    assert not _unbound_numbers(text), (
        f"Unbound quantitative README tokens: {_unbound_numbers(text)}"
    )


@pytest.mark.parametrize(
    "claim",
    [
        "17 taxa",
        "`17` taxa",
        "[17 taxa](source.tsv)",
        "| p | 0.05 |",
        "seven taxa",
        "12%",
        "2.4-fold",
        "z = \u22123.2",
        "p = 1e-5",
    ],
)
def test_numeric_census_refuses_new_unbound_results(claim: str) -> None:
    assert _unbound_numbers(claim), f"Census ignored an unbound result: {claim}"


@pytest.mark.parametrize(
    ("value", "expected"), [(1.0, "1.0"), (0.1, "0.10"), (0.001, "0.0010"), (0.000032, "3.2e-5")]
)
def test_probability_format_preserves_value_at_display_precision(
    value: float, expected: str
) -> None:
    assert _p(value) == expected


def test_duplicate_sections_are_refused() -> None:
    with pytest.raises(AssertionError, match="Duplicate README section"):
        _sections("## Overlap\ntext\n## Overlap\nmore text\n")


def test_every_quantitative_source_is_local_and_declared() -> None:
    assert len(SOURCES) == len(set(SOURCES))
    assert all((DATA / name).is_file() for name in SOURCES)
    assert all(Path(name).name == name for name in SOURCES)


def test_catalog_and_overlap_counts_agree() -> None:
    catalog = _rows("cellulose_decreased_taxa_catalog.tsv")
    core = _one(_rows("pairwise_jaccard_significance.tsv"), endpoint="three_way_core")
    assert len(catalog) == int(core["union"])
    assert sum(row["is_3way_core"] == "True" for row in catalog) == int(core["intersection"])
    ids = Counter(row["taxon_id"] for row in catalog)
    assert all(count == 1 for count in ids.values()), "Duplicate catalog taxon"
    tier = _rows("cellulose_decreased_taxa_catalog_tier1_2.tsv")
    summary = json.loads((DATA / "cellulose_decr_summary.json").read_text(encoding="utf-8"))
    assert (
        dict(Counter(row["core_class"] for row in tier)) == summary["catalog_tier1_2_core_classes"]
    )
    assert (
        sum(
            row["is_3way_core"] == "True" and row["has_direct_trait_edge"] == "True" for row in tier
        )
        == summary["core_taxa_with_direct_trait_edge_tier1_2"]
    )
