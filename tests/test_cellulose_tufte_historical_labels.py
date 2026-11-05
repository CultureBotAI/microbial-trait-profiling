"""Pin the visible historical and correction labels without rerunning analyses (#1416)."""

from __future__ import annotations

import html
import re
import xml.etree.ElementTree as ET
from collections.abc import Iterator

import matplotlib.pyplot as plt
import pytest
from matplotlib.figure import Figure

from mtp.cellulose_decr import tufte_plates


def _plain(text: str) -> str:
    return " ".join(html.unescape(re.sub(r"<[^>]+>", "", text)).split())


def _artist_texts(fig: Figure) -> list[str]:
    return [_plain(text.get_text()) for ax in fig.axes for text in ax.texts]


def _svg_texts(svg: str) -> list[str]:
    return [_plain("".join(node.itertext())) for node in ET.fromstring(svg).findall(".//{*}text")]


@pytest.fixture(scope="module")
def label_plates() -> Iterator[dict[str, Figure]]:
    """Construct artists only; output rendering is checked separately."""
    figures: dict[str, Figure] = {}
    builders = {
        "I": tufte_plates.build_plate1_convergence_slope_and_overlap,
        "II": tufte_plates.build_plate2_core_guild_dot_table,
        "III": tufte_plates.build_plate3_why_only_stool,
        "IV": tufte_plates.build_plate4_trait_landscape_and_four_pairs,
        "V": tufte_plates.build_plate5_three_statistical_arms,
        "VI": tufte_plates.build_plate6_combining_results_integrated_evidence,
    }
    try:
        for name, builder in builders.items():
            figures[name] = builder()
        yield figures
    finally:
        for fig in figures.values():
            plt.close(fig)


@pytest.mark.parametrize(
    ("plate", "original", "label"),
    [
        ("I", "PLATE I · CROSS-DISEASE TAXONOMIC CONVERGENCE", "October 16 worked example"),
        (
            "II",
            "PLATE II · THE SHARED FIBER-DEGRADING GUILD AT GENUS & SPECIES RESOLUTION",
            "historical (#1416)",
        ),
        ("III", "PLATE III · ANATOMICAL & ECOLOGICAL SPECIFICITY", "October 16 worked example"),
        (
            "IV",
            "PLATE IV · FUNCTIONAL PATHWAY BENCHMARK & LITERATURE-SUPPORTED SIGNATURES",
            "historical (#1416)",
        ),
        (
            "V",
            "PLATE V · STATISTICAL ARCHITECTURE (METHODS/THREE_ARM_STATISTICS.MD)",
            "historical (#1416)",
        ),
        (
            "VI",
            "PLATE VI · RESULT INTEGRATION & EVIDENCE CONTRACT (METHODS/INTEGRATED_EVIDENCE.MD)",
            "historical (#1416)",
        ),
    ],
)
def test_each_plate_kicker_names_its_own_historical_basis(
    label_plates: dict[str, Figure], plate: str, original: str, label: str
) -> None:
    """An incidental historical/specification string elsewhere cannot satisfy this."""
    assert _plain(label_plates[plate].texts[0].get_text()) == f"{original} · {label}"


@pytest.mark.parametrize(
    "stats",
    [
        "106 dec vs. 39 inc, nominal p = 2.4 \u00d7 10⁻⁸",
        "30 dec vs. 7 inc, nominal p = 1.9 \u00d7 10⁻⁴",
        "19 dec vs. 6 inc, nominal p = 0.015",
        "57 dec vs. 26 inc, nominal p = 8.8 \u00d7 10⁻⁴",
        "13 dec vs. 14 inc, nominal p = 1.00",
        "6 dec vs. 11 inc, nominal p = 0.33",
        "5 dec vs. 6 inc, nominal p = 1.00",
    ],
)
def test_each_plate_three_direction_row_identifies_nominal_p(
    label_plates: dict[str, Figure], stats: str
) -> None:
    matches = [text for text in _artist_texts(label_plates["III"]) if stats in text]
    assert len(matches) == 1
    assert matches[0].endswith(f"({stats})")


def test_plate_four_family_note_is_on_the_four_pair_panel(
    label_plates: dict[str, Figure],
) -> None:
    panel = label_plates["IV"].axes[1]
    subtitle = next(
        text.get_text()
        for text in panel.texts
        if text.get_text().startswith("Difference in share of stool profiles")
    )
    assert _plain(subtitle) == (
        "Difference in share of stool profiles with positive trait annotation "
        "(Increased minus Decreased; percentage points, pp) "
        "BH q: family not recorded (#1525)"
    )


def test_plate_six_qualifies_the_existing_drawn_cell_claim(
    label_plates: dict[str, Figure],
) -> None:
    drawn_status = [
        text for text in _artist_texts(label_plates["VI"]) if text.startswith("drawn (")
    ]
    assert drawn_status == ["drawn (316 cells, 16.9%) [test_status = untested]"]
    claim = next(
        text
        for text in _artist_texts(label_plates["VI"])
        if text.startswith("≥ 3 studies with trait signal; rendered on map")
    )
    assert claim == (
        "≥ 3 studies with trait signal; rendered on map (26 cells have BH q ≤ "
        "0.05 '*', 4 have q ≤ 0.001)"
    )
    note = next(
        text for text in label_plates["VI"].axes[2].texts if text.get_position() == (0.595, 0.42)
    )
    assert note.get_text() == (
        "historical claim (#1416); cells untested; family not recorded (#1526)"
    )


def test_plate_five_preserves_the_original_displayed_family_sizes(
    label_plates: dict[str, Figure],
) -> None:
    text = "\n".join(_artist_texts(label_plates["V"]))
    assert "4,670 profiles \u00d7 317 traits (1.48M tests; #131)" in text
    assert "BH across 1,480,390 tests; 66,240 enrichments" in text
    assert "Six Table 1 Single-Trait Fractions (BH across 6 traits):" in text
    assert "Key Table 1 Paired Log-Ratios (BY across 15 pairs; cancels T):" in text
    assert "BH (6 traits) & BY (15 ratios, 317 traits)" in text


def test_fig_seven_pins_the_kicker_and_unknown_family_note() -> None:
    texts = _svg_texts(tufte_plates.render_fig7_svg())
    assert texts[0] == (
        "PLATE VII · CROSS-DISEASE TAXONOMIC CONVERGENCE & 3-WAY VENN PARTITION "
        "(METPO:2000103 | GO:0030245) · October 16 worked example"
    )
    assert (
        "BH q: family not recorded (#1525); p-values nominal; "
        "benchmark/percentile: historical (#1416)."
    ) in texts
    assert "46 studies · Δ = \u221227.4 pp (q = 0.0077)" in texts


def test_fig_eight_pins_the_title_family_note() -> None:
    texts = _svg_texts(tufte_plates.render_fig8_svg())
    assert texts[0] == (
        "OPPOSITE STOOL TRAIT DIRECTIONS IN HUMAN SIGNATURES (BugSigDB \u00d7 KG-Microbe) "
        "· historical (#1416); BH q: family not recorded (#1525)"
    )


@pytest.mark.parametrize("card_index", range(4))
def test_each_fig_eight_card_identifies_its_fisher_p_as_nominal(card_index: int) -> None:
    root = ET.fromstring(tufte_plates.render_fig8_svg())
    cards = root.findall("{*}g")
    assert len(cards) == 4
    fisher_lines = [
        _plain("".join(node.itertext()))
        for node in cards[card_index].findall("{*}text")
        if "Fisher" in "".join(node.itertext())
    ]
    assert len(fisher_lines) == 1
    assert "Fisher nominal p = " in fisher_lines[0]


def test_fig_nine_has_an_explicit_visible_historical_label() -> None:
    assert "October 16 worked example" in _svg_texts(tufte_plates.render_fig9_svg())


@pytest.fixture(scope="module")
def deck_sections() -> list[str]:
    sections = re.findall(
        r"<section\b[^>]*>(.*?)</section>", tufte_plates.render_slide_deck_html(), re.DOTALL
    )
    assert len(sections) == 16
    return sections


def test_deck_cover_declares_its_historical_basis_in_the_lead(deck_sections: list[str]) -> None:
    lead = re.search(r'<p class="lead">(.*?)</p>', deck_sections[0], re.DOTALL)
    assert lead is not None
    assert "October 16 worked example" in _plain(lead.group(1))
    assert "90-trait benchmark (historical (#1416))" in _plain(deck_sections[0])
    assert "90-trait empirical distributions (historical (#1416))" in _plain(deck_sections[0])


@pytest.mark.parametrize("plate", ["I", "II", "III", "IV", "V", "VI"])
def test_plate_p_values_are_labelled_nominal(label_plates: dict[str, Figure], plate: str) -> None:
    fig = label_plates[plate]
    texts = [text.get_text() for text in fig.texts]
    texts.extend(text.get_text() for ax in fig.axes for text in ax.texts)
    for text in texts:
        for value in re.finditer(r"\bp\s*(?:=|≤|≥|<|>)", text):
            if plate == "I":
                assert text.startswith("(p-values nominal) "), text
            else:
                assert re.search(r"\bnominal\s+$", text[: value.start()]), text


def test_fig_eight_labels_every_nominal_p() -> None:
    for text in _svg_texts(tufte_plates.render_fig8_svg()):
        for value in re.finditer(r"\bp\s*(?:=|≤|≥|<|>)", text):
            assert re.search(r"\bnominal\s+$", text[: value.start()]), text


def test_plate_one_nominal_note_is_above_the_existing_footer_rule(
    label_plates: dict[str, Figure],
) -> None:
    captions = [
        text.get_text()
        for ax in label_plates["I"].axes
        for text in ax.texts
        if "If fiber-digesting bacteria were lost" in text.get_text()
    ]
    assert captions == [
        "(p-values nominal) If fiber-digesting bacteria were lost at the same rate "
        "as other stool bacteria,\n"
        "only ~20% would overlap between diseases. Observed overlap is 50.4%\n"
        "(exact hypergeometric p < 10⁻¹⁰ for every pair; random-guild permutation p = 0.011)."
    ]


@pytest.mark.parametrize("plate", ["I", "II", "III", "IV", "V", "VI"])
def test_new_plate_labels_fit_the_fixed_packet_canvas(
    label_plates: dict[str, Figure], plate: str
) -> None:
    fig = label_plates[plate]
    fig.draw_without_rendering()
    artists = list(fig.texts) + [text for ax in fig.axes for text in ax.texts]
    changed = [
        text
        for text in artists
        if any(
            label in text.get_text()
            for label in (
                "October 16 worked example",
                "historical (#1416)",
                "nominal p",
                "p-values nominal",
                "family not recorded (#1525)",
                "family not recorded (#1526)",
                "test_status = untested",
            )
        )
    ]
    for text in changed:
        bounds = text.get_window_extent()
        assert bounds.x0 >= 0, text.get_text()
        assert bounds.y0 >= 0, text.get_text()
        assert bounds.x1 <= fig.bbox.width, text.get_text()
        assert bounds.y1 <= fig.bbox.height, text.get_text()
        for other in artists:
            if other is text or not other.get_text():
                continue
            # A multiline artist's box includes blank space between lines.
            # Actual word collisions are checked on the rendered PDFs; here
            # test single-line labels, plus the full canvas bounds above.
            if "\n" in text.get_text() or "\n" in other.get_text():
                continue
            assert not bounds.overlaps(other.get_window_extent()), (
                text.get_text(),
                other.get_text(),
            )


@pytest.mark.parametrize("slide_number", range(2, 17))
def test_each_deck_slide_subtitle_declares_its_own_basis(
    deck_sections: list[str], slide_number: int
) -> None:
    subtitle = re.search(
        r'<p class="subtitle"[^>]*>(.*?)</p>', deck_sections[slide_number - 1], re.DOTALL
    )
    assert subtitle is not None
    if 8 <= slide_number <= 14:
        text_basis = (
            "historical (#1416)" if slide_number in {11, 14} else "October 16 worked example"
        )
        label = f"text: {text_basis}; figure drawn from current tables"
    else:
        label = (
            "historical (#1416)"
            if slide_number in {3, 5, 6, 7, 15}
            else "October 16 worked example"
        )
    assert _plain(subtitle.group(1)).endswith(f"· {label}")


@pytest.mark.parametrize("slide_number", range(2, 17))
def test_each_labelled_deck_subtitle_wraps_instead_of_truncating(
    deck_sections: list[str], slide_number: int
) -> None:
    """The deck CSS truncates subtitles with an ellipsis, which would hide a trailing label."""
    tag = re.search(r'<p class="subtitle"[^>]*>', deck_sections[slide_number - 1])
    assert tag is not None
    assert 'style="white-space:normal"' in tag.group(0)


@pytest.mark.parametrize(
    ("slide_number", "image_src"),
    [
        (8, "figures_graphical/fig1_overlap_and_null_models.svg"),
        (9, "figures_graphical/fig2_shared_taxa_and_set_definitions.svg"),
        (10, "figures_graphical/fig3_direction_by_body_site.svg"),
        (11, "figures_graphical/fig4_trait_overlap_ranking.svg"),
        (12, "figures/fig1_pairwise_and_3way_jaccard_nulls.svg"),
        (13, "figures/fig2_direction_by_body_site_and_definition.svg"),
        (14, "figures/fig3_core_taxa_and_trait_benchmark.svg"),
    ],
)
def test_current_table_image_labels_are_scoped_to_the_correct_output_route(
    deck_sections: list[str], slide_number: int, image_src: str
) -> None:
    """Pin the deck's public output routes from figures.main (#1544).

    The table-driven generator and its outputs are documented in
    methods/cellulose_decr.md (the figures section) and the analysis README.
    This checks each labelled slide's link, without executing or duplicating
    the generator's table loading, calculations or rendering.
    """
    section = deck_sections[slide_number - 1]
    assert re.findall(r'<img\s+src="([^"]+)"', section) == [image_src]
    subtitle = re.search(r'<p class="subtitle"[^>]*>(.*?)</p>', section, re.DOTALL)
    assert subtitle is not None
    assert "; figure drawn from current tables" in _plain(subtitle.group(1))


def test_deck_four_pair_footer_identifies_the_unknown_bh_family(deck_sections: list[str]) -> None:
    footer = re.search(r'<footer class="slide-footer">(.*?)</footer>', deck_sections[4], re.DOTALL)
    assert footer is not None
    pair_note = _plain(footer.group(1)).split("Four highlighted pairs:", 1)[1]
    assert "BH q: family not recorded (#1525)" in pair_note


@pytest.mark.parametrize("slide_number", range(1, 17))
def test_deck_p_labels_are_nominal_in_each_slide(
    deck_sections: list[str], slide_number: int
) -> None:
    own_html = re.sub(r"<svg\b.*?</svg>", "", deck_sections[slide_number - 1], flags=re.DOTALL)
    text = _plain(own_html)
    for value in re.finditer(r"\bp\s*(?:=|≤|≥|<|>)", text):
        context = text[max(0, value.start() - 45) : value.end() + 20]
        assert re.search(r"\bnominal\s+$", text[: value.start()]), (
            f"Slide {slide_number}: unlabelled p near {context!r}"
        )


@pytest.mark.parametrize(
    "invented",
    ["48-test family", "184 tested", "10 compartment rows", "three-framework-comparison-v1 run"],
)
def test_no_rejected_family_or_run_name_is_rendered(
    label_plates: dict[str, Figure], deck_sections: list[str], invented: str
) -> None:
    for fig in label_plates.values():
        texts = [_plain(text.get_text()) for text in fig.texts] + _artist_texts(fig)
        assert all(invented not in text for text in texts)
    svgs = (
        tufte_plates.render_fig7_svg(),
        tufte_plates.render_fig8_svg(),
        tufte_plates.render_fig9_svg(),
    )
    for svg in svgs:
        assert all(invented not in text for text in _svg_texts(svg))
    assert all(invented not in _plain(section) for section in deck_sections)
