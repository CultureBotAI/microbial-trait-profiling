"""Pin every current cellulose p-label and the scenario scope (#1869, #1870).

The committed tables cover all seven current figures; the small saved-table
fixture independently exercises the four p-bearing figures. These checks build
the real text artists without writing plots or recalculating analysis results.
"""

from __future__ import annotations

import inspect
import math
import re
from collections.abc import Callable

import matplotlib
import pytest
from matplotlib.figure import Figure
from matplotlib.text import Text

from mtp.cellulose_decr import figures
from mtp.cellulose_decr.analysis import ANALYSIS_DIR, ORIENTED
from test_cellulose_nominal_labels import synthetic_results

DRAWERS = (
    figures.graphical_convergence,
    figures.graphical_core,
    figures.graphical_compartments,
    figures.graphical_benchmark,
    figures.detailed_nulls,
    figures.detailed_compartments,
    figures.detailed_core_and_benchmark,
)
SYNTHETIC_DRAWERS = (
    figures.graphical_convergence,
    figures.graphical_compartments,
    figures.detailed_nulls,
    figures.detailed_compartments,
)
CAPTIONS = (
    figures.caption_overlap,
    figures.caption_shared_taxa,
    figures.caption_direction,
    figures.caption_ranking,
    figures.caption_nulls,
    figures.caption_body_sites,
    figures.caption_shared_by_two,
)
SYNTHETIC_CAPTIONS = (
    figures.caption_overlap,
    figures.caption_direction,
    figures.caption_nulls,
    figures.caption_body_sites,
)
STATEMENTS = (
    figures.overlap_statement,
    figures.crc_compartment_statement,
    figures.fragile_statement,
)
P_MENTION = re.compile(r"(?<![A-Za-z])p(?:-values?)?\s*(?:=|≤|<|>|from\b|\d)")


def _results(source: str) -> figures.Results:
    if source == "committed":
        return figures.load(ANALYSIS_DIR)
    assert source == "synthetic"
    results = synthetic_results()
    # Exercise the actual fragile statement as well as the stable verdicts.
    sensitive = results.scenarios["short_condition"].eq("HIV") & results.scenarios["scenario"].eq(
        "named_treated_cases_excluded"
    )
    results.scenarios.loc[sensitive, "binom_p_two_sided"] = 0.2
    return results


def _assert_nominal(text: str) -> None:
    text = " ".join(text.split())
    unlabelled = [
        text[max(0, match.start() - 45) : match.end() + 30]
        for match in P_MENTION.finditer(text)
        if re.search(r"\bnominal\s+$", text[: match.start()], re.IGNORECASE) is None
    ]
    assert not unlabelled, f"Unlabelled probability display: {unlabelled}"


def test_text_scan_covers_every_current_drawer_caption_and_statement() -> None:
    functions = {
        name
        for name, function in inspect.getmembers(figures, inspect.isfunction)
        if function.__module__ == figures.__name__
    }
    assert {draw.__name__ for draw in DRAWERS} == {
        name for name in functions if name.startswith(("graphical_", "detailed_"))
    }
    assert {caption.__name__ for caption in CAPTIONS} == {
        name for name in functions if name.startswith("caption_")
    }
    assert {statement.__name__ for statement in STATEMENTS} == {
        name for name in functions if name.endswith("_statement")
    }


@pytest.mark.parametrize(
    ("source", "draw"),
    [
        pytest.param(source, draw, id=f"{source}-{draw.__name__}")
        for source, drawers in (("committed", DRAWERS), ("synthetic", SYNTHETIC_DRAWERS))
        for draw in drawers
    ],
)
def test_every_current_figure_probability_is_labelled_nominal(
    source: str, draw: Callable[[figures.Results], Figure]
) -> None:
    with matplotlib.rc_context(figures.STYLE):
        figure = draw(_results(source))
        try:
            for artist in figure.findobj(Text):
                _assert_nominal(artist.get_text())
        finally:
            figure.clear()


@pytest.mark.parametrize(
    ("source", "caption"),
    [
        pytest.param(source, caption, id=f"{source}-{caption.__name__}")
        for source, captions in (("committed", CAPTIONS), ("synthetic", SYNTHETIC_CAPTIONS))
        for caption in captions
    ],
)
def test_every_current_caption_probability_is_labelled_nominal(
    source: str, caption: Callable[[figures.Results], figures.Caption]
) -> None:
    rendered = caption(_results(source))
    for text in (rendered.title, rendered.takeaway, rendered.caveat):
        _assert_nominal(text)


@pytest.mark.parametrize("source", ["committed", "synthetic"])
@pytest.mark.parametrize("statement", STATEMENTS, ids=lambda statement: statement.__name__)
def test_every_current_statement_probability_is_labelled_nominal(
    source: str, statement: Callable[[figures.Results], str]
) -> None:
    text = statement(_results(source))
    assert text
    _assert_nominal(text)


def test_mixed_null_statement_displays_each_probability_as_nominal() -> None:
    results = _results("synthetic")
    null = results.summary["null_models"]["U6"]
    null["mean_pairwise_jaccard_distance"]["p"] = 0.01
    null["three_way_core"]["p"] = 0.4
    text = figures.overlap_statement(results)
    assert len(P_MENTION.findall(text)) == 2
    _assert_nominal(text)


@pytest.mark.parametrize(
    ("site", "expected"),
    [
        ("pooled-experiment", "untied paired experiments (nominal p = 0.040)"),
        ("healthy-study", "disease vs healthy studies 8↓ 2↑, nominal p = 0.040"),
        ("body-site-odds", "2.00 (nominal p = 0.040)"),
    ],
)
def test_previously_unpinned_probability_sites_keep_their_nominal_label(
    site: str, expected: str
) -> None:
    results = _results("synthetic")
    if site == "body-site-odds":
        text = figures.caption_body_sites(results).takeaway
        assert text.count(expected) == 3
    else:
        with matplotlib.rc_context(figures.STYLE):
            figure = figures.graphical_compartments(results)
            try:
                text = " ".join(artist.get_text() for artist in figure.findobj(Text))
            finally:
                figure.clear()
    assert expected in " ".join(text.split())
    _assert_nominal(text)


def test_all_unavailable_scenario_ranges_have_no_invented_numeric_bounds() -> None:
    results = _results("synthetic")
    results.scenarios.loc[:, "binom_p_two_sided"] = math.nan
    verdicts = figures.direction_verdicts(results)
    assert set(verdicts) == {"PD", "HIV", "CRC", "Pooled"}
    for verdict in verdicts.values():
        assert verdict.definitions == 3
        assert math.isnan(verdict.low)
        assert math.isnan(verdict.high)
        assert figures.p_span(verdict, "nominal p from ", " to ") == "nominal p unavailable"


@pytest.mark.parametrize(("source", "endpoint_count"), [("committed", 32), ("synthetic", 16)])
def test_direction_caption_names_the_verdict_endpoints_without_changing_the_count(
    source: str, endpoint_count: int
) -> None:
    results = _results(source)
    scenario_rows = results.scenarios[
        results.scenarios["contrast_population"].eq(ORIENTED)
        & results.scenarios["short_condition"].isin(("PD", "HIV", "CRC", "Pooled"))
        & results.scenarios["unit"].eq("study")
    ]
    # The four primary study probabilities repeat the exact-label scenarios.
    # Three distinct CRC interactions and the counted-tissue sensitivity add four.
    assert len(scenario_rows) + 4 == endpoint_count
    caveat = " ".join(figures.caption_direction(results).caveat.split())
    assert f"Nominal p-values: {endpoint_count} finite endpoints; 0 unavailable." in caveat
    assert "including the scenario endpoints behind the verdicts" in caveat
    assert "not a correction family" in caveat
