"""Current cellulose figure/caption nominal-p display contracts (#1430).

The synthetic tables are deliberately small; repeated p values belong to distinct
endpoints. No committed analysis tables or retired figure routes are consumed.
"""

from __future__ import annotations

import importlib.util
import inspect
import itertools
import json
import math
import re
import runpy
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from importlib.machinery import ModuleSpec
from pathlib import Path
from typing import Any

import matplotlib
import pandas as pd
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.text import Text

from mtp.cellulose_decr import figures
from mtp.cellulose_decr.analysis import (
    ALL_CONTRASTS,
    NON_FECAL_GI,
    ORAL,
    ORIENTED,
    OTHER,
    POOLED_STRATIFIED,
)

NULL_KEYS = ("U1", "U2", "U3", "U4", "U5", "U6", "U6_size_matched")
NULL_COLUMNS = (
    "u1",
    "u2",
    "u3",
    "u4_shared_trait",
    "u5_emp_weighted",
    "u6_signature_perm",
    "u6_size_matched",
)
SCENARIOS = ("exact_labels", "synonyms_merged", "named_treated_cases_excluded")


def synthetic_results() -> figures.Results:
    """Small saved-table-shaped inputs for all four current p-bearing canvases."""
    pairs = ("PD vs. HIV", "PD vs. CRC", "HIV vs. CRC")
    nulls: dict[str, Any] = {}
    for key in NULL_KEYS:
        nulls[key] = {
            "iterations": 99,
            "mean_pairwise_jaccard_distance": {
                "p": 0.04,
                "expected": 0.7,
                "p2_5": 0.4,
                "p97_5": 0.9,
            },
            "three_way_core": {"p": 0.04, "expected": 1.5, "p2_5": 0, "p97_5": 3},
        }
    summary = {
        "bugsigdb_vintage": "synthetic fixture",
        "condition_labels": "exact",
        "core_taxa_with_direct_trait_edge": 1,
        "universes": {"U1_all_cellulose_taxa_in_matrix": 7},
        "observed": {
            "mean_pairwise_jaccard_distance": 2 / 3,
            "three_way_core": 1,
            "three_way_union": 7,
            "two_or_more_way": 4,
            "sizes": dict.fromkeys(figures.SHORT, 4),
            "jaccard_distance": dict.fromkeys(pairs, 2 / 3),
            "intersections": dict.fromkeys(pairs, 2),
            "unions": dict.fromkeys(pairs, 6),
        },
        "four_highlighted_pairs": {
            f"{index}_{short}_Feces_Cellulose_dec": {
                "n_supporting_studies_in_target_dir": 10,
            }
            for index, short in enumerate(figures.SHORT, start=1)
        },
        "null_models": nulls,
    }
    pairwise = []
    for name in (*pairs, "Mean of 3 pairs"):
        row: dict[str, Any] = {
            "comparison": name,
            "endpoint": "mean_pairwise_jaccard_distance" if name.startswith("Mean") else "pair",
            "jaccard_distance": 2 / 3,
        }
        for column in NULL_COLUMNS:
            row[f"{column}_expected_jaccard_distance"] = 0.7
            if column in ("u1", "u2", "u3"):
                row[f"{column}_hypergeom_p"] = 0.04
                row[f"{column}_mc_p"] = 0.04
            else:
                row[f"{column}_p"] = 0.04
        pairwise.append(row)
    memberships = (
        (True, True, True),
        (True, False, True),
        (True, True, False),
        (False, True, True),
        (True, False, False),
        (False, False, True),
        (False, True, False),
    )
    catalog = pd.DataFrame(
        [
            {
                "taxon_name": f"Synthetic taxon {index}",
                "in_PD_feces_dec": pd_member,
                "in_HIV_feces_dec": hiv_member,
                "in_CRC_feces_dec": crc_member,
            }
            for index, (pd_member, hiv_member, crc_member) in enumerate(memberships)
        ]
    )
    compartments = []
    for population in (ALL_CONTRASTS, ORIENTED):
        for compartment in ("Feces", NON_FECAL_GI, ORAL, OTHER):
            for short in (*figures.SHORT, "Pooled"):
                compartments.append(
                    {
                        "contrast_population": population,
                        "compartment": compartment,
                        "short_condition": short,
                        "paired_exp_down_dec_gt_inc": 8,
                        "paired_exp_up_inc_gt_dec": 2,
                        "study_down": 8,
                        "study_up": 2,
                        "study_tied": 0,
                        "paired_experiments": 10,
                        "n_studies_dec": 10,
                        "share_down_untied": 0.8,
                        "binom_p_two_sided": 0.04,
                        "study_binom_p_two_sided": 0.04,
                        "min_attainable_binom_p_two_sided": 2 / 2**10,
                    }
                )
    interaction = []
    for population in (ALL_CONTRASTS, ORIENTED):
        for short in ("CRC", POOLED_STRATIFIED):
            for compartment in (NON_FECAL_GI, ORAL, OTHER):
                for unit in ("experiment", "study"):
                    interaction.append(
                        {
                            "contrast_population": population,
                            "short_condition": short,
                            "compartment_b": compartment,
                            "unit": unit,
                            "fisher_p_two_sided": 0.04,
                            "conditional_or": 2.0,
                            "conditional_or_ci95_low": 0.5,
                            "conditional_or_ci95_high": 8.0,
                            "a_down": 8,
                            "a_up": 2,
                            "b_down": 4,
                            "b_up": 2,
                        }
                    )
    scenarios = pd.DataFrame(
        [
            {
                "scenario": scenario,
                "contrast_population": ORIENTED,
                "short_condition": short,
                "unit": "study",
                "binom_p_two_sided": 0.04,
                "down": 8,
                "up": 2,
                "tied": 0,
            }
            for short in (*figures.SHORT, "Pooled")
            for scenario in SCENARIOS
        ]
    )
    counted = dict(interaction[0])
    counted.update(
        {
            "scenario": "paired_tissue_counted",
            "contrast_population": ORIENTED,
            "short_condition": "CRC",
            "unit": "study",
            "fisher_p_two_sided": 0.12,
        }
    )
    sensitivity = pd.DataFrame(
        [
            {
                "definition_slug": slug,
                "comparison": name,
                "jaccard_distance": 2 / 3,
                "n_a": 4,
                "n_b": 4,
                "intersection": 1 if name.startswith("Mean") else 2,
                "union": 7 if name.startswith("Mean") else 6,
            }
            for slug in figures.SENSITIVITY_LABELS
            for name in (*pairs, "Mean of 3 pairs")
        ]
    )
    return figures.Results(
        summary=summary,
        pairwise=pd.DataFrame(pairwise),
        sensitivity=sensitivity,
        compartments=pd.DataFrame(compartments),
        interaction=pd.DataFrame(interaction),
        catalog=catalog,
        benchmark=pd.DataFrame(),
        scenarios=scenarios,
        scenario_interaction=pd.DataFrame([counted]),
    )


def _figure_text(figure: Figure) -> str:
    return " ".join(artist.get_text() for artist in figure.findobj(Text))


def _caption_text(caption: figures.Caption) -> str:
    return " ".join((caption.title, caption.takeaway, caption.caveat))


def _display(
    consumer: Callable[[figures.Results], Figure | figures.Caption], results: figures.Results
) -> str:
    try:
        rendered = consumer(results)
    except (ValueError, OverflowError, TypeError) as error:
        raise AssertionError(
            f"A saved p-value's availability must not break its display: {error}"
        ) from error
    return _figure_text(rendered) if isinstance(rendered, Figure) else _caption_text(rendered)


def _scope(text: str, finite: int, unavailable: int) -> None:
    compact = " ".join(text.split())
    assert f"Nominal p-values: {finite} finite endpoints; {unavailable} unavailable." in compact
    assert "not a correction family" in compact


NULL_VIEWS = [
    pytest.param(figures.graphical_convergence, 7, id="graphical-7"),
    pytest.param(figures.detailed_nulls, 35, id="detailed-35"),
    pytest.param(figures.caption_overlap, 6, id="overlap-caption-6"),
    pytest.param(figures.caption_nulls, 10, id="null-caption-10"),
]


@pytest.mark.parametrize(("consumer", "endpoint_count"), NULL_VIEWS)
def test_null_views_count_endpoint_identities_not_values_or_repeated_labels(
    consumer: Callable[[figures.Results], Figure | figures.Caption], endpoint_count: int
) -> None:
    """All p are equal: seven nulls are seven tests, and a repeated U6 label is not eight."""
    text = _display(consumer, synthetic_results())
    _scope(text, endpoint_count, 0)
    assert "nominal p" in text.lower()


@pytest.mark.parametrize(
    ("consumer", "finite", "unavailable"),
    [
        pytest.param(figures.graphical_convergence, 5, 2, id="graphical-5-plus-2"),
        pytest.param(figures.detailed_nulls, 32, 3, id="detailed-32-plus-3"),
        pytest.param(figures.caption_overlap, 4, 2, id="overlap-caption-4-plus-2"),
        pytest.param(figures.caption_nulls, 7, 3, id="null-caption-7-plus-3"),
    ],
)
def test_null_views_report_nonfinite_endpoints_without_universal_p_claims(
    consumer: Callable[[figures.Results], Figure | figures.Caption],
    finite: int,
    unavailable: int,
) -> None:
    results = synthetic_results()
    for key, missing in (("U2", math.nan), ("U3", math.inf)):
        results.summary["null_models"][key]["mean_pairwise_jaccard_distance"]["p"] = missing
        mean = results.pairwise["endpoint"].eq("mean_pairwise_jaccard_distance")
        results.pairwise.loc[mean, f"{key.lower()}_mc_p"] = missing
    results.summary["null_models"]["U5"]["three_way_core"]["p"] = math.nan
    text = _display(consumer, results)
    _scope(text, finite, unavailable)
    assert re.search(r"nominal p[^.;\n]{0,24}unavailable", text, re.IGNORECASE)
    assert not re.search(r"\b(?:nan|inf|n/a)\b", text, re.IGNORECASE)
    assert "every p" not in text
    assert "every mean" not in text
    assert "not tested" not in text.lower()


def test_exact_zero_remains_a_finite_nominal_endpoint_in_the_real_null_table() -> None:
    results = synthetic_results()
    selected = results.pairwise["comparison"].eq("PD vs. HIV")
    results.pairwise.loc[selected, "u1_hypergeom_p"] = 0.0
    text = _display(figures.detailed_nulls, results)
    _scope(text, 35, 0)
    assert re.search(r"exact\s+nominal p\s*(?:=\s*)?0(?:\D|$)", text)


def test_monte_carlo_zero_keeps_the_existing_resolution_floor_and_endpoint_count() -> None:
    results = synthetic_results()
    results.summary["null_models"]["U6"]["mean_pairwise_jaccard_distance"]["p"] = 0.0
    text = _display(figures.graphical_convergence, results)
    _scope(text, 7, 0)
    # The same U6 endpoint appears in the subtitle and null panel.
    assert text.count("nominal p ≤ 0.010") == 2


def test_missing_saved_null_value_is_unavailable_without_a_new_test_claim() -> None:
    results = synthetic_results()
    results.summary["null_models"]["U6"]["mean_pairwise_jaccard_distance"]["p"] = None
    text = _display(figures.graphical_convergence, results)
    _scope(text, 6, 1)
    assert text.count("nominal p unavailable") == 2
    assert "not tested" not in text.lower()


COMPARTMENT_VIEWS = [
    pytest.param(figures.graphical_compartments, 33, id="graphical-33"),
    pytest.param(figures.detailed_compartments, 12, id="detailed-12"),
    pytest.param(figures.caption_direction, 16, id="direction-caption-16"),
    pytest.param(figures.caption_body_sites, 3, id="body-sites-caption-3"),
]


@pytest.mark.parametrize(("consumer", "endpoint_count"), COMPARTMENT_VIEWS)
def test_compartment_views_count_only_the_selected_nominal_endpoints(
    consumer: Callable[[figures.Results], Figure | figures.Caption], endpoint_count: int
) -> None:
    results = synthetic_results()
    text = _display(consumer, results)
    _scope(text, endpoint_count, 0)
    assert "nominal p" in text.lower()
    if consumer is figures.graphical_compartments:
        assert "p-values are two-sided sign tests" not in text
        assert "interaction" in text.lower()


def test_graphical_compartment_count_follows_the_actual_sign_test_display_branches() -> None:
    results = synthetic_results()
    all_stool = results.compartments["contrast_population"].eq(ALL_CONTRASTS) & (
        results.compartments["compartment"].eq("Feces")
    )
    pd_row = all_stool & results.compartments["short_condition"].eq("PD")
    results.compartments.loc[pd_row, "min_attainable_binom_p_two_sided"] = 0.25
    hiv_row = all_stool & results.compartments["short_condition"].eq("HIV")
    results.compartments.loc[hiv_row, ["study_down", "study_up"]] = [1, 1]
    # PD contributes neither all-contrast p; HIV contributes its experiment p only.
    # Their oriented study p-values, and every displayed footnote endpoint, remain.
    text = _display(figures.graphical_compartments, results)
    _scope(text, 30, 0)
    assert "too few pairs" in text
    assert "too few studies" in text


@pytest.mark.parametrize(("label_basis", "finite"), [("exact", 16), ("synonyms", 15)])
def test_caption_primary_scenario_alias_uses_the_run_condition_label_basis(
    label_basis: str, finite: int
) -> None:
    results = synthetic_results()
    results.summary["condition_labels"] = label_basis
    # Without the exact CRC scenario, the exact primary adds its own endpoint;
    # the synonym primary still repeats the retained synonyms_merged scenario.
    drop = results.scenarios["short_condition"].eq("CRC") & results.scenarios["scenario"].eq(
        "exact_labels"
    )
    results.scenarios.drop(index=results.scenarios.index[drop], inplace=True)
    _scope(_display(figures.caption_direction, results), finite, 0)


@pytest.mark.parametrize(
    ("consumer", "finite"),
    [
        pytest.param(figures.graphical_compartments, 32, id="graphical"),
        pytest.param(figures.detailed_compartments, 11, id="detailed"),
        pytest.param(figures.caption_direction, 15, id="direction-caption"),
    ],
)
def test_unavailable_interaction_p_does_not_remove_the_finite_odds_ratio(
    consumer: Callable[[figures.Results], Figure | figures.Caption], finite: int
) -> None:
    results = synthetic_results()
    chosen = (
        results.interaction["contrast_population"].eq(ALL_CONTRASTS)
        & results.interaction["short_condition"].eq("CRC")
        & results.interaction["compartment_b"].eq(NON_FECAL_GI)
        & results.interaction["unit"].eq("experiment")
    )
    results.interaction.loc[chosen, "fisher_p_two_sided"] = math.nan
    text = _display(consumer, results)
    _scope(text, finite, 1)
    assert "nominal p unavailable" in text
    assert "not tested" not in text.lower()
    assert results.interaction.loc[chosen, "conditional_or"].item() == 2.0


@pytest.mark.parametrize(
    ("consumer", "finite"),
    [
        pytest.param(figures.graphical_compartments, 33, id="graphical"),
        pytest.param(figures.detailed_compartments, 12, id="detailed"),
    ],
)
def test_finite_interaction_p_counts_even_when_the_odds_ratio_is_not_estimable(
    consumer: Callable[[figures.Results], Figure | figures.Caption], finite: int
) -> None:
    results = synthetic_results()
    chosen = (
        results.interaction["contrast_population"].eq(ALL_CONTRASTS)
        & results.interaction["short_condition"].eq("CRC")
        & results.interaction["compartment_b"].eq(NON_FECAL_GI)
        & results.interaction["unit"].eq("experiment")
    )
    results.interaction.loc[chosen, "conditional_or"] = math.nan
    _scope(_display(consumer, results), finite, 0)


def test_sensitivity_range_retains_an_unavailable_interior_scenario_endpoint() -> None:
    results = synthetic_results()
    hiv = results.scenarios["short_condition"].eq("HIV")
    results.scenarios.loc[hiv, "binom_p_two_sided"] = [0.02, math.nan, 0.04]
    extra = results.scenarios.loc[hiv].iloc[0].to_dict()
    extra.update({"scenario": "extra_synthetic_definition", "binom_p_two_sided": 0.2})
    results.scenarios.loc[len(results.scenarios)] = extra
    primary_hiv = (
        results.compartments["short_condition"].eq("HIV")
        & results.compartments["compartment"].eq("Feces")
        & results.compartments["contrast_population"].eq(ORIENTED)
    )
    results.compartments.loc[primary_hiv, "study_binom_p_two_sided"] = 0.02
    text = _display(figures.caption_direction, results)
    _scope(text, 16, 1)
    assert "nominal p from 0.020 to 0.20" in text
    assert "unavailable" in text
    # The standalone exported statement summarizes all four scenario endpoints,
    # including the missing interior member; it must not count just its extrema.
    assert "include_scope" in inspect.signature(figures.fragile_statement).parameters
    fragile = figures.fragile_statement(results, include_scope=True)
    _scope(fragile, 3, 1)
    assert "nominal p" in fragile


@pytest.mark.parametrize("available", [False, True], ids=["without-viz", "with-viz"])
def test_collection_guard_treats_both_cellulose_figure_suites_together(
    monkeypatch: pytest.MonkeyPatch, available: bool
) -> None:
    find_spec = importlib.util.find_spec

    def simulated_find_spec(name: str, package: str | None = None) -> ModuleSpec | None:
        if name == "matplotlib" and not available:
            return None
        return find_spec(name, package)

    monkeypatch.setattr(importlib.util, "find_spec", simulated_find_spec)
    config = runpy.run_path(str(Path(__file__).with_name("conftest.py")))
    ignored = config["collect_ignore"]
    for name in ("test_cellulose_decr_figures.py", "test_cellulose_nominal_labels.py"):
        assert (name not in ignored) is available


@pytest.mark.parametrize(
    ("consumer", "finite"),
    [
        pytest.param(figures.overlap_statement, 2, id="overlap-two"),
        pytest.param(figures.crc_compartment_statement, 4, id="crc-four"),
    ],
)
def test_standalone_scope_is_separate_from_the_embedded_scientific_sentence(
    consumer: Callable[..., str], finite: int
) -> None:
    assert "include_scope" in inspect.signature(consumer).parameters
    results = synthetic_results()
    sentence = consumer(results)
    assert "Nominal p-values:" not in sentence
    exported = consumer(results, include_scope=True)
    assert exported.startswith(sentence)
    _scope(exported, finite, 0)


def test_standalone_scope_omits_a_conditional_crc_probability_that_is_not_shown() -> None:
    assert "include_scope" in inspect.signature(figures.crc_compartment_statement).parameters
    results = synthetic_results()
    results.scenario_interaction.loc[:, "fisher_p_two_sided"] = 0.04
    exported = figures.crc_compartment_statement(results, include_scope=True)
    _scope(exported, 3, 0)
    assert "counting the tumour-vs-adjacent-tissue contrasts too" not in exported


@pytest.mark.parametrize(
    ("statement_name", "finite"), [("overlap", 2), ("crc_compartments", 4), ("fragile", 3)]
)
def test_main_exports_each_standalone_scope_without_nesting_caption_counts(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, statement_name: str, finite: int
) -> None:
    results = synthetic_results()
    sensitive = results.scenarios["short_condition"].eq("HIV") & results.scenarios["scenario"].eq(
        "named_treated_cases_excluded"
    )
    results.scenarios.loc[sensitive, "binom_p_two_sided"] = 0.2
    for name in figures.FIGURE_INPUTS:
        (tmp_path / name).write_text("synthetic input binding\n", encoding="utf-8")
    monkeypatch.setattr(figures, "load", lambda _: results)
    for name in (
        "graphical_convergence",
        "graphical_core",
        "graphical_compartments",
        "graphical_benchmark",
        "detailed_nulls",
        "detailed_compartments",
        "detailed_core_and_benchmark",
    ):
        monkeypatch.setattr(figures, name, lambda _: Figure(figsize=(1, 1)))
    for name in ("caption_shared_taxa", "caption_ranking", "caption_shared_by_two"):
        monkeypatch.setattr(
            figures,
            name,
            lambda _: figures.Caption(
                "Descriptive", "Synthetic", "No test", figures.FIGURE_INPUTS[0]
            ),
        )

    @contextmanager
    def saved_book(path: Path, **_kwargs: Any) -> Iterator[None]:
        yield None
        path.write_bytes(b"synthetic book save boundary\n")

    def saved_figure(_figure: Figure, directory: Path, stem: str, _book: Any) -> list[Path]:
        directory.mkdir(parents=True, exist_ok=True)
        paths = [directory / f"{stem}{suffix}" for suffix in (".svg", ".pdf", ".png")]
        for path in paths:
            path.write_bytes(b"synthetic figure save boundary\n")
        return paths

    monkeypatch.setattr(figures, "PdfPages", saved_book)
    monkeypatch.setattr(figures, "save", saved_figure)
    assert figures.main(["--analysis-dir", str(tmp_path)]) == 0
    saved = json.loads((tmp_path / "figures" / "figure_captions.json").read_text())
    assert "Nominal p-values:" not in saved["statements"][statement_name]
    _scope(saved["statement_scopes"][statement_name], finite, 0)
    caption_counts = {
        "figures_graphical/fig1_overlap_and_null_models.svg": 6,
        "figures_graphical/fig3_direction_by_body_site.svg": 16,
        "figures/fig1_pairwise_and_3way_jaccard_nulls.svg": 10,
        "figures/fig2_direction_by_body_site_and_definition.svg": 3,
    }
    for name, expected in caption_counts.items():
        caption = saved["figures"][name]
        text = " ".join((caption["title"], caption["takeaway"], caption["caveat"]))
        _scope(text, expected, 0)
        assert text.count("Nominal p-values:") == 1


def test_detailed_null_probability_labels_fit_between_neighboring_cells() -> None:
    """Measure actual text, so restoring the long single-line label fails."""
    with matplotlib.rc_context(figures.STYLE):
        figure = figures.detailed_nulls(synthetic_results())
        canvas: Any = FigureCanvasAgg(figure)
        canvas.draw()
        renderer = canvas.get_renderer()
        table = figure.axes[0]
        labels = [text for text in table.texts if "nominal p" in text.get_text()]
        assert len(labels) == 28
        rows: dict[float, list[Text]] = {}
        for label in labels:
            x, y = label.get_position()
            rows.setdefault(y, []).append(label)
            expected = [
                text
                for text in table.texts
                if text.get_position()[0] == x
                and text.get_text() == "0.700"
                and text.get_position()[1] > y
            ]
            nearest = min(expected, key=lambda text: text.get_position()[1] - y)
            assert label.get_window_extent(renderer).y1 < nearest.get_window_extent(renderer).y0
        assert len(rows) == 4
        for row in rows.values():
            ordered = sorted(row, key=lambda text: text.get_position()[0])
            for left, right in itertools.pairwise(ordered):
                assert left.get_window_extent(renderer).x1 < right.get_window_extent(renderer).x0
