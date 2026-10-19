"""The cellulose_decr figures and the slide captions they write (#556, #570, #631).

They need matplotlib from the viz extra; conftest.py leaves this module out when it
is missing, and the analysis tests in test_cellulose_decr.py still run."""

from __future__ import annotations

import json
import shutil
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from mtp.cellulose_decr import figures
from mtp.cellulose_decr.analysis import (
    ANALYSIS_DIR,
    ORIENTED,
)

# --------------------------------------------------------------------------- figures


def test_figures_draw_from_the_committed_tables(tmp_path: Path) -> None:
    """Every figure renders from the committed tables, with a slide caption (#556).

    The file names describe the figures rather than claims the analysis rejects (#570),
    and the captions quote the tables.
    """
    for name in figures.FIGURE_INPUTS:
        shutil.copy(ANALYSIS_DIR / name, tmp_path / name)
    commit = "0123456789" * 4
    assert figures.main(["--analysis-dir", str(tmp_path), "--code-commit", commit]) == 0
    provenance = json.loads((tmp_path / "figures" / "figures_provenance.json").read_text())
    assert provenance["script"] == "src/mtp/cellulose_decr/figures.py"
    # Decks link the figures' code and tables at this commit (#556).
    assert provenance["code_commit"] == commit
    assert len(provenance["outputs"]) == 7 * 3 + 2 + 1
    for output in provenance["outputs"]:
        assert (tmp_path / output).stat().st_size > 0
    for claim in ("surprising", "convergence", "guild", "specificity", "podium"):
        assert not [output for output in provenance["outputs"] if claim in output], claim

    captions = json.loads((tmp_path / "figures" / "figure_captions.json").read_text())
    svgs = {output for output in provenance["outputs"] if output.endswith(".svg")}
    assert set(captions["figures"]) == svgs
    for caption in captions["figures"].values():
        assert all(caption[field] for field in ("title", "takeaway", "caveat"))
        assert (ANALYSIS_DIR / caption["evidence"]).exists()

    summary = json.loads((ANALYSIS_DIR / "cellulose_decr_summary.json").read_text())
    observed = summary["observed"]
    overlap = captions["figures"]["figures_graphical/fig1_overlap_and_null_models.svg"]
    assert overlap["title"].startswith(
        f"{observed['three_way_core']} of {observed['three_way_union']} "
    )
    assert overlap["title"] == captions["statements"]["overlap"]
    # PD, CRC and the pooled result keep p < 0.05 under every condition definition; HIV
    # does not (SIGNIFICANT_RESULTS_SUMMARY.md, Table A2).
    direction = captions["figures"]["figures_graphical/fig3_direction_by_body_site.svg"]
    assert direction["title"].endswith(
        "in PD, CRC and all three pooled under every condition definition; in HIV only under some"
    )
    assert "HIV 9 : 0 : 2" in direction["takeaway"]
    assert captions["statements"]["fragile"].startswith("the HIV result depends")
    assert captions["statements"]["crc_compartments"] in direction["caveat"]


# --------------------------------------------------------------------------- captions


@pytest.fixture(scope="module")
def results() -> figures.Results:
    return figures.load(ANALYSIS_DIR)


def _with_scenarios(
    results: figures.Results, edit: Callable[[pd.DataFrame], None]
) -> figures.Results:
    table = results.scenarios.copy()
    edit(table)
    return replace(results, scenarios=table)


def _healthy_studies(table: pd.DataFrame, short: str) -> pd.Series:
    return (
        table["contrast_population"].eq(ORIENTED)
        & table["short_condition"].eq(short)
        & table["unit"].eq("study")
    )


def test_direction_caption_follows_the_direction_not_only_p(results: figures.Results) -> None:
    """A condition significant the other way is named as such, not as a decrease (#589)."""

    def swap_crc(table: pd.DataFrame) -> None:
        rows = _healthy_studies(table, "CRC")
        table.loc[rows, ["down", "up"]] = table.loc[rows, ["up", "down"]].to_numpy()

    title = figures.caption_direction(_with_scenarios(results, swap_crc)).title
    assert "the increased signature does in CRC under every definition" in title
    assert "in PD and all three pooled under every condition definition" in title

    def crc_up_under_one(table: pd.DataFrame) -> None:
        swap_crc(table)
        rows = table.index[_healthy_studies(table, "CRC")]
        table.loc[rows[1:], "binom_p_two_sided"] = 0.5

    edited = _with_scenarios(results, crc_up_under_one)
    assert "the increased signature does in CRC only under some" in (
        figures.caption_direction(edited).title
    )
    assert "CRC: p from " in figures.caption_direction(edited).caveat
    fragile = figures.fragile_statement(edited)
    assert fragile.startswith("the HIV and CRC results depend")
    # Quoted beside decreased-direction claims, it must say which way CRC goes (#608).
    assert fragile.endswith(" toward the increased signature)")


def test_mixed_direction_is_named_not_only_its_p_range(results: figures.Results) -> None:
    """Small p-values in opposite directions are not a robust result (#589)."""

    def flip_one_hiv(table: pd.DataFrame) -> None:
        rows = table.index[_healthy_studies(table, "HIV") & table["binom_p_two_sided"].lt(0.05)]
        down, up = table.loc[rows[0], ["down", "up"]]
        table.loc[rows[0], ["down", "up"]] = [up, down]

    edited = _with_scenarios(results, flip_one_hiv)
    assert figures.direction_verdicts(edited)["HIV"].kind == "mixed"
    caption = figures.caption_direction(edited)
    assert "HIV changes direction with the definition" in caption.title
    assert "HIV: p from " in caption.caveat
    assert " in opposite directions across " in caption.caveat
    assert figures.fragile_statement(edited).endswith(" in opposite directions)")


def test_fragile_pooled_result_starts_its_sentence_and_is_named(
    results: figures.Results,
) -> None:
    def pooled_under_one(table: pd.DataFrame) -> None:
        rows = table.index[_healthy_studies(table, "Pooled")]
        table.loc[rows[1:], "binom_p_two_sided"] = 0.5

    edited = _with_scenarios(results, pooled_under_one)
    assert "All three pooled: p from " in figures.caption_direction(edited).caveat
    statement = figures.fragile_statement(edited)
    assert statement.startswith("the HIV and pooled results depend on the condition definition")
    assert "(study-level p HIV " in statement
    assert "; pooled " in statement


def test_direction_caption_without_significance_or_data(results: figures.Results) -> None:
    def no_significance(table: pd.DataFrame) -> None:
        table["binom_p_two_sided"] = 0.5

    title = figures.caption_direction(_with_scenarios(results, no_significance)).title
    assert title.endswith("in no condition at p < 0.05")
    assert " not " not in title

    def no_hiv(table: pd.DataFrame) -> None:
        table.drop(table.index[_healthy_studies(table, "HIV")], inplace=True)

    title = figures.caption_direction(_with_scenarios(results, no_hiv)).title
    assert title.endswith("; no data for HIV")
    assert figures.fragile_statement(_with_scenarios(results, no_hiv)) == ""


@pytest.mark.parametrize(
    ("mean_p", "core_p", "ending"),
    [
        (0.19, 0.73, "random sets of whole stool signatures overlap about as much"),
        (0.01, 0.02, "they overlap more than random sets of whole stool signatures"),
        (0.99, 0.98, "they overlap less than random sets of whole stool signatures"),
        (0.04, 0.73, "mean-distance p = 0.040 and core p = 0.73"),
    ],
)
def test_overlap_statement_reads_both_u6_endpoints(
    results: figures.Results, mean_p: float, core_p: float, ending: str
) -> None:
    summary = json.loads(json.dumps(results.summary))
    summary["null_models"]["U6"]["mean_pairwise_jaccard_distance"]["p"] = mean_p
    summary["null_models"]["U6"]["three_way_core"]["p"] = core_p
    assert figures.overlap_statement(replace(results, summary=summary)).endswith(ending)


def test_shared_taxa_caption_counts_direct_edges(results: figures.Results) -> None:
    catalog = results.catalog.copy()
    core = catalog.index[catalog["is_3way_core"].astype(bool)]
    catalog.loc[core[0], "has_direct_trait_edge"] = True
    summary = {**results.summary, "core_taxa_with_direct_trait_edge": 1}
    edited = replace(results, catalog=catalog, summary=summary)
    caption = figures.caption_shared_taxa(edited)
    assert caption.title.endswith(", 1 annotated directly")
    assert caption.takeaway.startswith("1 carries its own cellulose annotation edge")
    assert figures.direct_core_text(edited).startswith("1 of the 17 core taxa has its own")


def test_crc_statement_names_the_excluded_experiment_only_when_it_matters(
    results: figures.Results,
) -> None:
    assert "without the tumour-vs-adjacent-tissue experiment" in (
        figures.crc_compartment_statement(results)
    )
    table = results.scenario_interaction.copy()
    primary = figures.interaction_row(results, ORIENTED, "CRC", "study")
    excluded = table["scenario"].eq("crc_within_patient_excluded")
    table.loc[excluded, "fisher_p_two_sided"] = primary["fisher_p_two_sided"]
    statement = figures.crc_compartment_statement(replace(results, scenario_interaction=table))
    assert "without" not in statement
    assert "does not show anatomical specificity" in statement


def test_p_values_keep_two_significant_digits() -> None:
    assert figures.sci(0.0097) == "0.0097"
    assert figures.sci(0.0596) == "0.060"
    # More digits where two would round across 0.05 (#641)
    assert figures.sci(0.0496) == "0.0496"
    assert figures.sci(0.04996) == "0.04996"
    assert figures.sci(0.19) == "0.19"
    assert figures.sci(1.2e-5) == "1.2" + figures.TIMES + "10⁻⁵"


def test_a_missing_benchmark_row_is_named() -> None:
    """A bare StopIteration said nothing about the missing row (#639)."""
    with pytest.raises(ValueError, match="no non-growth-medium row for METPO:x"):
        figures.trait_row([{"trait_id": "METPO:y"}], "METPO:x")


def test_a_failed_figure_run_drops_the_old_captions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """New figures must not sit beside the previous run's captions and commit (#639)."""
    for name in figures.FIGURE_INPUTS:
        shutil.copy(ANALYSIS_DIR / name, tmp_path / name)
    (tmp_path / "figures").mkdir()
    for name in ("figure_captions.json", "figures_provenance.json"):
        (tmp_path / "figures" / name).write_text("{}")

    def fail(_results: figures.Results) -> None:
        raise RuntimeError("figure failed")

    monkeypatch.setattr(figures, "graphical_core", fail)
    with pytest.raises(RuntimeError, match="figure failed"):
        figures.main(["--analysis-dir", str(tmp_path), "--code-commit", "0" * 40])
    assert not (tmp_path / "figures" / "figure_captions.json").exists()
    assert not (tmp_path / "figures" / "figures_provenance.json").exists()
