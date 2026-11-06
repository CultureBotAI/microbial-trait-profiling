"""Bind historical overlap/direction text to the saved study-modal source (#1551, #1555).

The fixture is the verbatim cellulose_decr_summary.json from commit 7970cba;
the tests use its label_source_comparison.study_modal object. This preserves
the historical population; current experiment-labelled outputs are not substituted.
These tests construct text artists and SVG strings, never render or save a figure.
"""

from __future__ import annotations

import html
import json
import re
import xml.etree.ElementTree as ET
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pytest
from matplotlib.figure import Figure

from mtp.cellulose_decr import tufte_plates

SOURCE = (
    Path(__file__).resolve().parent
    / "fixtures/cellulose_tufte_historical_sources/study_modal_7970cba.json"
)
PAIR_LABELS = ("PD vs. CRC", "PD vs. HIV", "HIV vs. CRC")
CONDITIONS = ("PD", "HIV", "CRC")


def _plain(value: str) -> str:
    return " ".join(html.unescape(value).split())


def _texts(figure: Figure) -> list[str]:
    return [
        _plain(artist.get_text())
        for artist in [*figure.texts, *(text for axis in figure.axes for text in axis.texts)]
    ]


def _one(values: list[str], prefix: str) -> str:
    found = [value for value in values if value.startswith(prefix)]
    assert len(found) == 1, (prefix, found)
    return found[0]


@pytest.fixture(scope="module")
def reference() -> dict[str, Any]:
    result: dict[str, Any] = json.loads(SOURCE.read_text(encoding="utf-8"))[
        "label_source_comparison"
    ]["study_modal"]
    return result


@pytest.fixture(scope="module")
def plates() -> Iterator[dict[str, Figure]]:
    figures = {
        "I": tufte_plates.build_plate1_convergence_slope_and_overlap(),
        "II": tufte_plates.build_plate2_core_guild_dot_table(),
        "III": tufte_plates.build_plate3_why_only_stool(),
        "IV": tufte_plates.build_plate4_trait_landscape_and_four_pairs(),
    }
    try:
        yield figures
    finally:
        for figure in figures.values():
            plt.close(figure)


@pytest.fixture(scope="module")
def seventh() -> ET.Element:
    return ET.fromstring(tufte_plates.render_fig7_svg())


@pytest.fixture(scope="module")
def eighth() -> ET.Element:
    return ET.fromstring(tufte_plates.render_fig8_svg())


@pytest.fixture(scope="module")
def ninth() -> ET.Element:
    return ET.fromstring(tufte_plates.render_fig9_svg())


def _printed_p(value: str) -> float:
    value = value.translate(str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-"))
    return float(value.replace(" \u00d7 10", "e"))


def _svg_texts(root: ET.Element) -> list[str]:
    return [_plain("".join(node.itertext())) for node in root.findall(".//{*}text")]


def _group(root: ET.Element, identity: str) -> ET.Element:
    node = root.find(f".//{{*}}g[@id='{identity}']")
    assert node is not None, identity
    return node


@pytest.mark.parametrize("pair", PAIR_LABELS)
def test_plate_one_pair_labels_read_historical_intersections_and_unions(
    plates: dict[str, Figure], reference: dict[str, Any], pair: str
) -> None:
    sets = reference["sets"]
    shared, union = sets["intersections"][pair], sets["unions"][pair]
    percent = 100 * (1 - sets["jaccard_distance"][pair])
    texts = _texts(plates["I"])
    assert f"{percent:.1f}% {pair} ({shared}/{union} shared)" in texts
    first, second = pair.split(" vs. ")
    extra = shared - sets["three_way_core"]
    assert (
        f"= {shared} shared of {union} ({percent:.1f}%) (+{extra} in {first} & {second} only)"
    ) in texts


def test_plate_one_partition_is_derived_from_the_historical_sets(
    plates: dict[str, Figure], reference: dict[str, Any]
) -> None:
    sets = reference["sets"]
    core, union = sets["three_way_core"], sets["three_way_union"]
    shared = sets["two_or_more_way"]
    texts = _texts(plates["I"])
    assert texts.count(f"{core} taxa ({100 * core / union:.1f}%)") == 2
    only_two = shared - core
    assert f"{only_two} taxa ({100 * only_two / union:.1f}%)" in texts
    expected = (
        f"Shared by at least 2 diseases: {shared} of {union} taxa ({100 * shared / union:.1f}%)"
    )
    assert _one(texts, "Shared by at least 2 diseases:").startswith(expected)
    assert f"B. How the {union} Lost Fiber-Digesting Taxa Are Shared Across Diseases" in texts
    assert (
        f"1. Partition of all {union} lost fiber-digesting taxa "
        "into three mutually exclusive tiers:"
    ) in texts
    pair_only = [sets["intersections"][pair] - core for pair in PAIR_LABELS]
    assert _one(texts, "Shared by at least 2 diseases:").endswith(
        f"= {core} in all three + {pair_only[0]} in PD & CRC only "
        f"+ {pair_only[1]} in PD & HIV only + {pair_only[2]} in HIV & CRC only"
    )
    mean_percent = 100 * (1 - sets["mean_pairwise_jaccard_distance"])
    mean_labels = [value for value in texts if "3-pair mean" in value]
    assert len(mean_labels) == 1
    assert mean_labels[0].startswith(f"{mean_percent:.1f}% ")
    assert (
        _one(texts, "(p-values nominal)")
        .split("Observed overlap is ")[1]
        .startswith(f"{mean_percent:.1f}%")
    )


def test_plate_one_study_counts_use_the_original_decreased_study_population(
    plates: dict[str, Figure], reference: dict[str, Any]
) -> None:
    counts = [
        reference["feces_direction"][condition]["feces_studies_dec"] for condition in CONDITIONS
    ]
    texts = _texts(plates["I"])
    source = _one(texts, "Source: BugSigDB human stool signatures")
    assert source.startswith(
        f"Source: BugSigDB human stool signatures ({sum(counts)} studies: "
        f"{counts[0]} PD, {counts[1]} HIV, {counts[2]} CRC)"
    )
    subtitle = _plain(plates["I"].texts[2].get_text())
    assert f"across {sum(counts)} stool studies," in subtitle


@pytest.mark.parametrize("pair", [*PAIR_LABELS, "3-Pair Mean"])
def test_plate_two_all_rank_comparator_uses_the_same_historical_sets(
    plates: dict[str, Figure], reference: dict[str, Any], pair: str
) -> None:
    # The genus/species panel has unresolved historical provenance. This check
    # covers only its explicitly all-rank comparator, also shown in Plate I.
    axis = plates["II"].axes[1]
    points = {"PD vs. CRC": 0.58, "PD vs. HIV": 0.48, "HIV vs. CRC": 0.44, "3-Pair Mean": 0.50}
    labels = [text.get_text() for text in axis.texts if text.get_position() == (0.15, points[pair])]
    distance = (
        reference["sets"]["mean_pairwise_jaccard_distance"]
        if pair == "3-Pair Mean"
        else reference["sets"]["jaccard_distance"][pair]
    )
    percent = 100 * (1 - distance)
    assert labels == [f"{percent:.1f}%"]
    assert f"All Ranks ({reference['sets']['three_way_union']} taxa)" in _texts(plates["II"])


@pytest.mark.parametrize("condition", CONDITIONS)
def test_plate_three_feces_direction_counts_percent_and_p_match_the_historical_source(
    plates: dict[str, Figure], reference: dict[str, Any], condition: str
) -> None:
    record = reference["feces_direction"][condition]
    down, up = record["paired_down"], record["paired_up"]
    prefix = f"{100 * down / (down + up):.1f}% ({down} dec vs. {up} inc, nominal p = "
    text = _one(_texts(plates["III"]), prefix)
    displayed = text.removeprefix(prefix).removesuffix(")")
    assert _printed_p(displayed) == float(f"{record['binom_p_two_sided']:.2g}")


def test_plate_three_pooled_counts_and_coverage_follow_the_original_condition_rows(
    plates: dict[str, Figure], reference: dict[str, Any]
) -> None:
    direction = reference["feces_direction"]
    down = sum(direction[condition]["paired_down"] for condition in CONDITIONS)
    up = sum(direction[condition]["paired_up"] for condition in CONDITIONS)
    studies = [direction[condition]["feces_studies_dec"] for condition in CONDITIONS]
    texts = _texts(plates["III"])
    assert _one(texts, f"{100 * down / (down + up):.1f}% (").startswith(
        f"{100 * down / (down + up):.1f}% ({down} dec vs. {up} inc, nominal p = "
    )
    assert f"{100 * down / (down + up):.1f}% of experiments" in _plain(
        plates["III"].texts[2].get_text()
    )
    coverage = _one(texts, "Stool is the primary sample type")
    assert f"({studies[0]} PD, {studies[1]} HIV, {studies[2]} CRC studies)" in coverage
    sets = reference["sets"]
    assert (
        f"{100 * (1 - sets['mean_pairwise_jaccard_distance']):.1f}% shared · "
        f"{sets['three_way_core']} core taxa in all 3 diseases ({sum(studies)} studies)"
    ) in texts


@pytest.mark.parametrize("condition", CONDITIONS)
def test_fig_seven_taxon_and_study_node_counts_read_the_historical_source(
    seventh: ET.Element, reference: dict[str, Any], condition: str
) -> None:
    node = _group(seventh, f"node_{condition.lower()}")
    assert (
        f"{reference['sets']['sizes'][condition]} taxa · "
        f"{reference['feces_direction'][condition]['feces_studies_dec']} studies"
    ) in _svg_texts(node)


@pytest.mark.parametrize("pair", PAIR_LABELS)
def test_fig_seven_edge_counts_and_distances_read_the_historical_sets(
    seventh: ET.Element, reference: dict[str, Any], pair: str
) -> None:
    first, second = pair.split(" vs. ")
    edge = _group(seventh, f"edge_{first.lower()}_{second.lower()}")
    sets = reference["sets"]
    shared, union = sets["intersections"][pair], sets["unions"][pair]
    distance = sets["jaccard_distance"][pair]
    core = sets["three_way_core"]
    texts = _svg_texts(edge)
    assert f"{first} ↔ {second}: {shared} / {union} shared" in texts
    assert (
        f"J = {100 * (1 - distance):.1f}% (d_J = {distance:.3f} · {core} + {shared - core})"
    ) in texts


def test_fig_seven_venn_partition_text_follows_historical_set_arithmetic(
    seventh: ET.Element, reference: dict[str, Any]
) -> None:
    sets = reference["sets"]
    core, union = sets["three_way_core"], sets["three_way_union"]
    shared = sets["two_or_more_way"]
    texts = _svg_texts(seventh)
    assert f"3-Way Core: {core} ({100 * core / union:.1f}%)" in texts
    assert f"2-Way Only: {shared - core} ({100 * (shared - core) / union:.1f}%)" in texts
    for condition in CONDITIONS:
        overlapping = sum(
            value for pair, value in sets["intersections"].items() if condition in pair
        )
        only = sets["sizes"][condition] - overlapping + core
        label = f"{condition} only ({100 * only / union:.1f}%)"
        assert label in texts
        assert texts[texts.index(label) - 1] == str(only)
    for pair in PAIR_LABELS:
        first, second = pair.split(" vs. ")
        extra = sets["intersections"][pair] - core
        label = f"{first} ∩ {second} ({100 * extra / union:.1f}%)"
        assert label in texts
        assert texts[texts.index(label) - 1] == str(extra)
    mean_distance = sets["mean_pairwise_jaccard_distance"]
    assert f"J = {100 * (1 - mean_distance):.1f}% (d_J = {mean_distance:.3f})" in texts
    assert texts[texts.index("3-Way Core") - 1] == str(core)
    assert f"Cellulose Guild (N = {union})" in texts
    assert _one(_svg_texts(_group(seventh, "card_3way")), f"{100 * core / union:.1f}%").startswith(
        f"{100 * core / union:.1f}% ({core} / {union}) "
    )
    assert _one(
        _svg_texts(_group(seventh, "card_2way")), f"{100 * shared / union:.1f}%"
    ).startswith(f"{100 * shared / union:.1f}% ({shared} / {union}) ")
    assert _one(_svg_texts(_group(seventh, "card_2way")), "Only ").startswith(
        f"Only {union - shared} / {union} taxa ({100 * (union - shared) / union:.1f}%) "
        "are unique to a single disease"
    )


@pytest.mark.parametrize("condition", CONDITIONS)
def test_plate_four_cellulose_direction_fragments_read_the_historical_source(
    plates: dict[str, Figure], reference: dict[str, Any], condition: str
) -> None:
    # Effect sizes, q-values and total study counts in these lines have separate
    # unresolved provenance. Only the scoped sign-test fragment is asserted.
    lines = [
        _plain(text.get_text())
        for text in plates["IV"].axes[1].texts
        if " Down vs " in text.get_text()
    ]
    assert len(lines) == len(CONDITIONS)
    line = lines[CONDITIONS.index(condition)]
    fragment = re.search(r"; (\d+) Down vs (\d+) Up, nominal p = ([^)]+)\)", line)
    assert fragment is not None, line
    record = reference["feces_direction"][condition]
    assert int(fragment.group(1)) == record["paired_down"]
    assert int(fragment.group(2)) == record["paired_up"]
    precision = 3 if condition == "HIV" else 2
    assert _printed_p(fragment.group(3)) == float(f"{record['binom_p_two_sided']:.{precision}g}")


@pytest.mark.parametrize("condition", CONDITIONS)
def test_fig_eight_cellulose_direction_and_untied_counts_read_the_historical_source(
    eighth: ET.Element, reference: dict[str, Any], condition: str
) -> None:
    cards = eighth.findall("{*}g")
    assert len(cards) == 4
    texts = _svg_texts(cards[CONDITIONS.index(condition)])
    record = reference["feces_direction"][condition]
    down, up = record["paired_down"], record["paired_up"]
    direction = _one(texts, "• Within-study experiment sign test:")
    prefix = (
        f"• Within-study experiment sign test: {down} Decreased vs. {up} Increased "
        f"({100 * down / (down + up):.1f}% decrease, binomial nominal p = "
    )
    assert direction.startswith(prefix)
    assert _printed_p(direction.removeprefix(prefix).removesuffix(")")) == float(
        f"{record['binom_p_two_sided']:.3g}"
    )
    coverage = _one(texts, "MTP Stool Signal (")
    counts = re.fullmatch(
        r"MTP Stool Signal \((\d+) Decreased Studies / \d+ Total, "
        r"(\d+) Paired Signature Experiments\):",
        coverage,
    )
    assert counts is not None, coverage
    assert int(counts.group(1)) == record["feces_studies_dec"]
    # The historical card prints the untied count here; no claim that this is
    # the total paired population (including tied experiments) is introduced.
    assert int(counts.group(2)) == down + up


def test_fig_nine_disease_totals_and_union_read_the_historical_sets(
    ninth: ET.Element, reference: dict[str, Any]
) -> None:
    texts = _svg_texts(ninth)
    sets = reference["sets"]
    for condition in CONDITIONS:
        assert f"Total |{condition}| = {sets['sizes'][condition]} taxa" in texts
    union = sets["three_way_union"]
    assert f"A. The 7 Mutually Exclusive Regions (Sum = {union} Taxa)" in texts
    assert (
        f"When a Venn diagram partitions the {union} union taxa into non-overlapping subsets:"
    ) in texts


def test_fig_nine_seven_region_counts_follow_historical_set_arithmetic(
    ninth: ET.Element, reference: dict[str, Any]
) -> None:
    texts = _svg_texts(ninth)
    sets = reference["sets"]
    core = sets["three_way_core"]
    regions = {"ALL 3": core}
    for condition in CONDITIONS:
        overlapping = sum(
            value for pair, value in sets["intersections"].items() if condition in pair
        )
        regions[f"{condition} only"] = sets["sizes"][condition] - overlapping + core
    for pair in PAIR_LABELS:
        regions[f"{pair.replace(' vs. ', ' & ')} only"] = sets["intersections"][pair] - core
    for label, count in regions.items():
        assert texts[texts.index(label) - 1] == str(count)
    assert sum(regions.values()) == sets["three_way_union"]


@pytest.mark.parametrize("pair", PAIR_LABELS)
def test_fig_nine_pairwise_counts_and_percentages_read_the_historical_sets(
    ninth: ET.Element, reference: dict[str, Any], pair: str
) -> None:
    cards = ninth.findall("{*}g")
    assert len(cards) == len(PAIR_LABELS)
    texts = _svg_texts(cards[PAIR_LABELS.index(pair)])
    sets = reference["sets"]
    first, second = pair.split(" vs. ")
    shared, union = sets["intersections"][pair], sets["unions"][pair]
    core = sets["three_way_core"]
    heading = texts[0]
    assert heading.endswith(
        f"(|{first} ∩ {second}| = {shared} of {union} union, "
        f"{100 * (1 - sets['jaccard_distance'][pair]):.1f}%)"
    )
    assert f"{core} Core Taxa (in all 3 diseases)" in texts
    increment = _one(texts, "+")
    count = re.match(r"\+\s*(\d+)", increment)
    assert count is not None
    assert int(count.group(1)) == shared - core
    assert f"= {shared} shared (> {core})" in texts


def test_fig_nine_summary_counts_and_percentage_follow_historical_set_arithmetic(
    ninth: ET.Element, reference: dict[str, Any]
) -> None:
    texts = _svg_texts(ninth)
    sets = reference["sets"]
    core, shared = sets["three_way_core"], sets["two_or_more_way"]
    only_two = shared - core
    slices = [sets["intersections"][pair] - core for pair in PAIR_LABELS]
    assert (
        f"Why the 2-way-only slices ({slices[0]}, {slices[1]}, {slices[2]}) "
        f"are smaller than {core}:"
    ) in texts
    assert (
        f"Each 2-way-only slice excludes the {core} core taxa. Total 2-way-only = "
        f"{slices[0]} + {slices[1]} + {slices[2]} = {only_two} taxa."
    ) in texts
    assert (
        f"Summary: Comparing the 3-Way Core ({core}) and the 2-Way-Only Subset ({only_two})"
    ) in texts
    assert (
        f"disease. As a result, {core} of the {shared} shared taxa "
        f"({100 * core / shared:.0f}%) fall into the 3-disease intersection, while {only_two}"
    ) in texts
    assert (
        f"taxa ({slices[0]} + {slices[1]} + {slices[2]}) "
        "appear in exactly two of the three diseases."
    ) in texts
