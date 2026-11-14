"""Guard SIGNIFICANT_RESULTS_SUMMARY.md pooled sign-test bound and table counts (#833)."""

from __future__ import annotations

import csv
import importlib.util
import re
import shutil
from pathlib import Path
from types import ModuleType

import pytest
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "data" / "analysis" / "cellulose_decr"
RECOUNT_SCRIPT = REPO_ROOT / "scripts" / "adhoc" / "recount_cellulose_readme_numbers_799.py"
LE = "\u2264"


def _load_recount_module(script_path: Path = RECOUNT_SCRIPT) -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "recount_cellulose_readme_numbers_799", script_path
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _format_p2(value: str | float) -> str:
    """Format a p-value to two significant figures in the cellulose_decr doc style."""
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


def _format_bound_2dec_sci(value: float) -> str:
    """Format the worst-case p-value bound as M.MMe-E (e.g. 1.21e-4)."""
    mantissa, exponent = f"{value:.2e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def assert_significant_results_summary_matches_tables(
    analysis_dir: Path = ANALYSIS_DIR,
    recount_script: Path = RECOUNT_SCRIPT,
) -> None:
    """Verify SIGNIFICANT_RESULTS_SUMMARY.md and README.md match the committed TSV tables (#833)."""
    summary_md = (analysis_dir / "SIGNIFICANT_RESULTS_SUMMARY.md").read_text(encoding="utf-8")
    readme_md = (analysis_dir / "README.md").read_text(encoding="utf-8")
    sensitivity_rows = _read_tsv(analysis_dir / "condition_definition_sensitivity.tsv")
    compartment_rows = _read_tsv(analysis_dir / "compartment_jaccard_and_direction_stats.tsv")
    sig_summary_rows = _read_tsv(analysis_dir / "significant_results_summary.tsv")

    pooled_study_rows = [
        row
        for row in sensitivity_rows
        if row["short_condition"] == "Pooled" and row["unit"] == "study"
    ]
    assert len(pooled_study_rows) == 10

    for row in pooled_study_rows:
        down = int(row["down"])
        up = int(row["up"])
        tied = int(row["tied"])
        assert tied >= 0
        assert down / (down + up) == pytest.approx(float(row["share_down_untied"]), rel=1e-5)
        exact_p = float(stats.binomtest(up, down + up, 0.5).pvalue)
        assert float(row["binom_p_two_sided"]) == pytest.approx(exact_p, rel=1e-5)

    by_key = {
        (row["scenario"], row["contrast_population"], row["short_condition"], row["unit"]): row
        for row in sensitivity_rows
    }
    exact_all = by_key[("exact_labels", "all_contrasts", "Pooled", "study")]
    assert (int(exact_all["down"]), int(exact_all["up"]), int(exact_all["tied"])) == (56, 21, 12)
    exact_oriented = by_key[("exact_labels", "disease_vs_healthy_oriented", "Pooled", "study")]
    assert (
        int(exact_oriented["down"]),
        int(exact_oriented["up"]),
        int(exact_oriented["tied"]),
    ) == (49, 14, 12)

    synonyms_all = by_key[("synonyms_merged", "all_contrasts", "Pooled", "study")]
    assert (int(synonyms_all["down"]), int(synonyms_all["up"]), int(synonyms_all["tied"])) == (
        60,
        24,
        12,
    )
    assert stats.binomtest(24, 84, 0.5).pvalue == pytest.approx(
        float(synonyms_all["binom_p_two_sided"]), rel=1e-5
    )

    without_rectal_all = by_key[("without_site_miscurations", "all_contrasts", "Pooled", "study")]
    assert (
        int(without_rectal_all["down"]),
        int(without_rectal_all["up"]),
        int(without_rectal_all["tied"]),
    ) == (55, 21, 12)
    assert stats.binomtest(21, 76, 0.5).pvalue == pytest.approx(
        float(without_rectal_all["binom_p_two_sided"]), rel=1e-5
    )

    worst_all = max(float(row["binom_p_two_sided"]) for row in pooled_study_rows)
    worst_oriented = max(
        float(row["binom_p_two_sided"])
        for row in pooled_study_rows
        if row["contrast_population"] == "disease_vs_healthy_oriented"
    )
    assert worst_all == pytest.approx(float(without_rectal_all["binom_p_two_sided"]), rel=1e-9)
    assert worst_oriented == pytest.approx(
        float(
            by_key[("synonyms_merged", "disease_vs_healthy_oriented", "Pooled", "study")][
                "binom_p_two_sided"
            ]
        ),
        rel=1e-9,
    )
    assert worst_all > 8.8e-5
    assert worst_all > worst_oriented

    expected_bound_str = _format_bound_2dec_sci(worst_all)
    assert expected_bound_str == "1.21e-4"

    summary_bound_match = re.search(
        rf"The pooled sign test stays at\s+p {LE}\s*([0-9.eE+-]+)\s+"
        r"under every label and contrast definition tried",
        summary_md,
    )
    assert summary_bound_match is not None, "Missing pooled sign-test bound sentence in SUMMARY"
    quoted_summary_bound = summary_bound_match.group(1)
    assert quoted_summary_bound == expected_bound_str
    assert float(quoted_summary_bound) >= worst_all

    readme_bound_match = re.search(
        rf"The pooled study-level sign test stays at p {LE}\s*([0-9.eE+-]+)\s+in every scenario\.",
        readme_md,
    )
    assert readme_bound_match is not None, "Missing pooled sign-test bound sentence in README"
    quoted_readme_bound = readme_bound_match.group(1)
    assert quoted_readme_bound == expected_bound_str
    assert float(quoted_readme_bound) >= worst_all

    all_down = int(exact_all["down"])
    all_untied = all_down + int(exact_all["up"])
    all_p_str = _format_p2(exact_all["binom_p_two_sided"])
    oriented_down = int(exact_oriented["down"])
    oriented_untied = oriented_down + int(exact_oriented["up"])
    oriented_p_str = _format_p2(exact_oriented["binom_p_two_sided"])
    expected_pooled_bullet = (
        f"- **Pooled:** {all_down} of {all_untied} untied studies (p = {all_p_str}) and "
        f"{oriented_down} of {oriented_untied} among\n"
        f"     disease-versus-healthy contrasts (p = {oriented_p_str})."
    )
    normalized_summary = " ".join(summary_md.split())
    normalized_bullet = " ".join(expected_pooled_bullet.split())
    assert normalized_bullet in normalized_summary

    sig_by_id = {row["test_id"]: row for row in sig_summary_rows}
    for pop in ("all_contrasts", "disease_vs_healthy_oriented"):
        pop_p_vals = [
            float(row["binom_p_two_sided"])
            for row in pooled_study_rows
            if row["contrast_population"] == pop
        ]
        interp = sig_by_id[f"feces_Pooled_study_direction_{pop}"]["summary_interpretation"]
        range_match = re.search(
            r"Across the condition definitions of condition_definition_sensitivity\.tsv "
            r"this p runs from ([0-9.eE+-]+) to ([0-9.eE+-]+)\.",
            interp,
        )
        assert range_match is not None
        low_str, high_str = range_match.groups()
        assert float(low_str) == pytest.approx(min(pop_p_vals), rel=1e-2)
        assert float(high_str) == pytest.approx(max(pop_p_vals), rel=1e-2)

    comp_by_key = {
        (row["contrast_population"], row["compartment"], row["short_condition"]): row
        for row in compartment_rows
    }
    populations = (("all_contrasts", "all"), ("disease_vs_healthy_oriented", "disease vs healthy"))
    conditions = (("PD", "PD"), ("HIV", "HIV"), ("CRC", "CRC"), ("Pooled", "All three"))
    for pop, pop_label in populations:
        for cond, cond_label in conditions:
            comp_row = comp_by_key[(pop, "Feces", cond)]
            q_val = _format_p2(
                sig_by_id[f"feces_{cond}_study_direction_{pop}"]["bh_q_within_table"]
            )
            expected_table_a_row = (
                f"| {cond_label} | {pop_label} | "
                f"{comp_row['paired_exp_down_dec_gt_inc']} : "
                f"{comp_row['paired_exp_up_inc_gt_dec']} : {comp_row['paired_exp_tied']} | "
                f"{_format_p2(comp_row['binom_p_two_sided'])} | "
                f"{comp_row['study_down']} : {comp_row['study_up']} : {comp_row['study_tied']} | "
                f"{_format_p2(comp_row['study_binom_p_two_sided'])} | {q_val} | "
                f"{_format_p2(comp_row['study_wilcoxon_p_two_sided'])} |"
            )
            assert expected_table_a_row in summary_md

    hiv_table_a2_scenarios = (
        ("exact_labels", "Exact label (primary)"),
        ("synonyms_merged", "+ synonym labels"),
        (
            "named_treated_cases_excluded",
            "Contrasts whose cases are named as treated set aside",
        ),
        ("control_tiers_1_2", "Tier 1\u20132 controls only"),
        ("paired_tissue_counted", "Paired-tissue contrasts also counted"),
        ("comorbidity_counted", "Comorbidity contrasts also counted"),
        (
            "without_site_miscurations",
            "Exact label, without the rectal-mucosa study",
        ),
    )
    for scenario_slug, label in hiv_table_a2_scenarios:
        exp_row = by_key[(scenario_slug, "disease_vs_healthy_oriented", "HIV", "experiment")]
        study_row = by_key[(scenario_slug, "disease_vs_healthy_oriented", "HIV", "study")]
        expected_a2_row = (
            f"| {label} | {exp_row['down']} : {exp_row['up']} : {exp_row['tied']} | "
            f"{_format_p2(exp_row['binom_p_two_sided'])} | "
            f"{study_row['down']} : {study_row['up']} : {study_row['tied']} | "
            f"{_format_p2(study_row['binom_p_two_sided'])} |"
        )
        assert expected_a2_row in summary_md

    recount_mod = _load_recount_module(recount_script)
    assert recount_mod.run(analysis_dir) == 0


def test_pooled_sign_test_bound_matches_worst_scenario_across_all_populations() -> None:
    """Pooled sign-test bound in SIGNIFICANT_RESULTS_SUMMARY.md holds across scenarios (#833)."""
    assert_significant_results_summary_matches_tables(ANALYSIS_DIR, RECOUNT_SCRIPT)


@pytest.mark.parametrize(
    ("doc_name", "target", "replacement"),
    [
        (
            "SIGNIFICANT_RESULTS_SUMMARY.md",
            f"p {LE} 1.21e-4 under every label",
            f"p {LE} 8.8e-5 under every label",
        ),
        (
            "SIGNIFICANT_RESULTS_SUMMARY.md",
            f"p {LE} 1.21e-4 under every label",
            f"p {LE} 1.50e-4 under every label",
        ),
        (
            "SIGNIFICANT_RESULTS_SUMMARY.md",
            "56 of 77 untied studies (p = 8.2e-5)",
            "58 of 78 untied studies (p = 1.9e-5)",
        ),
        (
            "SIGNIFICANT_RESULTS_SUMMARY.md",
            "49 of 63 among\n     disease-versus-healthy contrasts (p = 1.1e-5)",
            "49 of 64 among\n     disease-versus-healthy contrasts (p = 2.4e-5)",
        ),
        (
            "SIGNIFICANT_RESULTS_SUMMARY.md",
            "| All three | all | 107 : 37 : 30 | 4.4e-9 | 56 : 21 : 12 | 8.2e-5 |",
            "| All three | all | 107 : 37 : 30 | 4.4e-9 | 58 : 20 : 11 | 1.9e-5 |",
        ),
        (
            "SIGNIFICANT_RESULTS_SUMMARY.md",
            "| + synonym labels | 13 : 5 : 3 | 0.096 | 10 : 4 : 2 | 0.18 |",
            "| + synonym labels | 13 : 5 : 3 | 0.096 | 10 : 2 : 2 | 0.039 |",
        ),
        (
            "README.md",
            f"stays at p {LE} 1.21e-4 in every scenario.",
            f"stays at p {LE} 8.8e-5 in every scenario.",
        ),
    ],
)
def test_significant_results_summary_validator_rejects_perturbed_bounds_and_counts(
    tmp_path: Path,
    doc_name: str,
    target: str,
    replacement: str,
) -> None:
    """Perturbing the pooled sign-test bound or study counts in the summary docs fails (#833)."""
    scratch_analysis = tmp_path / "cellulose_decr"
    shutil.copytree(ANALYSIS_DIR, scratch_analysis)
    doc_path = scratch_analysis / doc_name
    original = doc_path.read_text(encoding="utf-8")
    assert target in original
    doc_path.write_text(original.replace(target, replacement, 1), encoding="utf-8")

    with pytest.raises(AssertionError):
        assert_significant_results_summary_matches_tables(scratch_analysis, RECOUNT_SCRIPT)
