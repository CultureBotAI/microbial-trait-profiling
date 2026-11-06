"""Keep the retained audit display's historical numbers fixed (#1709).

The fixture's baseline_commit records original review provenance at 92539eb.
Its capture helper executes only the pinned source blob beginning 0672ee7d,
which survives squash merges; current producer code cannot rebaseline it.
It preserves historical numeric tokens, caveats and warnings while excluding
only the exact retirement banner. It does not certify scientific correctness.
Artists are constructed without raster rendering or writing figure products.
"""

from __future__ import annotations

import json
import re
import unicodedata
import xml.etree.ElementTree as ET
from collections.abc import Callable, Iterable
from html.parser import HTMLParser
from pathlib import Path
from typing import Any, cast

import matplotlib.pyplot as plt
import pytest
from matplotlib.figure import Figure
from matplotlib.text import Text

from mtp.cellulose_decr import tufte_plates

BUILDERS = (
    "build_plate1_convergence_slope_and_overlap",
    "build_plate2_core_guild_dot_table",
    "build_plate3_why_only_stool",
    "build_plate4_trait_landscape_and_four_pairs",
    "build_plate5_three_statistical_arms",
    "build_plate6_combining_results_integrated_evidence",
)
SVG_RENDERERS = ("render_fig7_svg", "render_fig8_svg", "render_fig9_svg")
SLIDES = tuple(f"slide-{number:02d}" for number in range(1, 17))
FIXTURE = Path(__file__).parent / "fixtures/cellulose_historical_numeric_snapshot.json"
NUMBER = re.compile(r"[+-]?(?:\d+(?:,\d{3})*(?:\.\d+)?|\.\d+)(?:[eE][+-]?\d+)?")
EXCLUDED_RETIREMENT_BANNER = "RETIRED: historical, unverified numbers (#1572)"
Snapshot = dict[str, dict[str, list[str]]]


@pytest.fixture(autouse=True)
def _local_import() -> None:
    assert Path(tufte_plates.__file__).resolve() == (
        Path(__file__).resolve().parents[1] / "src/mtp/cellulose_decr/tufte_plates.py"
    )


def numeric_tokens(texts: Iterable[str]) -> list[str]:
    """Exclude the exact retirement banner, then preserve ordered normalized numbers."""
    tokens: list[str] = []
    for text in texts:
        historical = text.replace(EXCLUDED_RETIREMENT_BANNER, "")
        plain = unicodedata.normalize("NFKC", historical).replace("\N{MINUS SIGN}", "-")
        tokens.extend(NUMBER.findall(plain))
    return tokens


def figure_numeric_tokens(figure: Figure) -> list[str]:
    """Read figure text, axes text/titles, then labels/ticks on enabled axes."""
    artists: list[Text] = list(figure.texts)
    for axes in figure.axes:
        if not axes.get_visible():
            continue
        artists.extend(axes.texts)
        artists.append(axes.title)
        if axes.axison:
            for axis in (axes.xaxis, axes.yaxis):
                if axis.get_visible():
                    artists.extend((axis.label, axis.get_offset_text()))
                    artists.extend(axis.get_ticklabels(minor=False))
                    artists.extend(axis.get_ticklabels(minor=True))
    return numeric_tokens(artist.get_text() for artist in artists if artist.get_visible())


def svg_numeric_tokens(svg: str) -> list[str]:
    """Read SVG text/tspan content, excluding geometry and style attributes."""
    root = ET.fromstring(svg)
    return numeric_tokens("".join(node.itertext()) for node in root.findall(".//{*}text"))


class SlideTextParser(HTMLParser):
    """Separate all slide text, including inline SVG, from CSS and image routes."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.slides: dict[str, list[str]] = {}
        self.parts: list[str] | None = None
        self.ignored_depth = 0

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        classes = (dict(attrs).get("class") or "").split()
        if tag == "section" and "slide" in classes:
            assert self.parts is None, "Unexpected nested slide"
            self.parts = []
        if self.parts is not None and tag in {"script", "style", "noscript", "title"}:
            self.ignored_depth += 1

    def handle_endtag(self, tag: str) -> None:
        if self.parts is not None and tag in {"script", "style", "noscript", "title"}:
            assert self.ignored_depth > 0
            self.ignored_depth -= 1
        if tag == "section" and self.parts is not None:
            assert self.ignored_depth == 0
            self.slides[f"slide-{len(self.slides) + 1:02d}"] = numeric_tokens(self.parts)
            self.parts = None

    def handle_data(self, data: str) -> None:
        if self.parts is not None and self.ignored_depth == 0:
            self.parts.append(data)


def html_numeric_tokens(html: str) -> dict[str, list[str]]:
    parser = SlideTextParser()
    parser.feed(html)
    parser.close()
    assert parser.parts is None, "Unclosed historical slide"
    return parser.slides


@pytest.fixture(scope="module")
def golden_numeric_tokens() -> Snapshot:
    payload: dict[str, Any] = json.loads(FIXTURE.read_text(encoding="utf-8"))
    assert payload["schema_version"] == 1
    assert payload["baseline_commit"] == "92539eb"
    assert payload["source"] == "src/mtp/cellulose_decr/tufte_plates.py"
    snapshot = payload["numeric_tokens"]
    assert isinstance(snapshot, dict)
    assert set(snapshot) == {"builders", "svgs", "slides"}
    for group, names in (("builders", BUILDERS), ("svgs", SVG_RENDERERS), ("slides", SLIDES)):
        assert set(snapshot[group]) == set(names)
        for name in names:
            values = snapshot[group][name]
            assert isinstance(values, list)
            assert values, f"Empty numeric snapshot: {group}/{name}"
            assert all(isinstance(value, str) and NUMBER.fullmatch(value) for value in values)
    return cast(Snapshot, snapshot)


@pytest.fixture(scope="module")
def historical_html_tokens(golden_numeric_tokens: Snapshot) -> dict[str, list[str]]:
    # Require the reviewed fixture before constructing any display for comparison.
    assert set(golden_numeric_tokens["slides"]) == set(SLIDES)
    return html_numeric_tokens(tufte_plates.render_slide_deck_html())


@pytest.mark.parametrize("builder_name", BUILDERS)
def test_historical_builder_numeric_snapshot(
    golden_numeric_tokens: Snapshot, builder_name: str
) -> None:
    builder: Callable[[], Figure] = getattr(tufte_plates, builder_name)
    figure = builder()
    try:
        assert figure_numeric_tokens(figure) == golden_numeric_tokens["builders"][builder_name]
    finally:
        plt.close(figure)


@pytest.mark.parametrize("renderer_name", SVG_RENDERERS)
def test_historical_svg_numeric_snapshot(
    golden_numeric_tokens: Snapshot, renderer_name: str
) -> None:
    renderer: Callable[[], str] = getattr(tufte_plates, renderer_name)
    assert svg_numeric_tokens(renderer()) == golden_numeric_tokens["svgs"][renderer_name]


@pytest.mark.parametrize("slide", SLIDES)
def test_historical_html_slide_numeric_snapshot(
    golden_numeric_tokens: Snapshot, historical_html_tokens: dict[str, list[str]], slide: str
) -> None:
    assert set(historical_html_tokens) == set(SLIDES)
    assert historical_html_tokens[slide] == golden_numeric_tokens["slides"][slide]


def test_numeric_tokenization_keeps_signs_exponents_and_display_order() -> None:
    assert numeric_tokens(
        ["\N{MINUS SIGN}53.7 pp; 1,872 tests", "2.4 \N{MULTIPLICATION SIGN} 10⁻⁸; +.5; 1e-06"]
    ) == [
        "-53.7",
        "1,872",
        "2.4",
        "10",
        "-8",
        "+.5",
        "1e-06",
    ]


def test_markup_numeric_tokens_exclude_layout_and_keep_nested_visible_text() -> None:
    svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1400 740">'
        '<text x="99">n=17; <tspan>p=0.05</tspan></text></svg>'
    )
    assert svg_numeric_tokens(svg) == ["17", "0.05"]
    html = (
        '<style>.slide {height: 1036px}</style><section class="slide">'
        '<img src="fig9.svg"><p>17 taxa; p=2.4 \N{MULTIPLICATION SIGN} '
        "10<sup>\N{MINUS SIGN}8</sup></p>"
        "<script>const count = 999;</script></section>"
    )
    assert html_numeric_tokens(html) == {"slide-01": ["17", "2.4", "10", "-8"]}


def test_numeric_tokens_exclude_only_the_exact_retirement_banner() -> None:
    banner = "RETIRED: historical, unverified numbers (#1572)"
    assert numeric_tokens(
        [
            f"n=1572; {banner}; 16 slides",
            banner + banner,
            "RETIRED: historical, unverified numbers (#1708)",
            "historical sample count: 1572",
        ]
    ) == ["1572", "16", "1708", "1572"]


def test_figure_numeric_tokens_exclude_banner_but_retain_same_number_in_data() -> None:
    figure = Figure()
    figure.text(0.1, 0.9, "RETIRED: historical, unverified numbers (#1572)")
    figure.text(0.1, 0.5, "n=1572; p=0.05")
    assert figure_numeric_tokens(figure) == ["1572", "0.05"]


def test_markup_numeric_tokens_ignore_single_and_duplicate_exact_banners() -> None:
    banner = "RETIRED: historical, unverified numbers (#1572)"
    svg = (
        '<svg xmlns="http://www.w3.org/2000/svg">'
        f"<text>{banner}</text><text>n=1572; p=0.05</text></svg>"
    )
    assert svg_numeric_tokens(svg) == ["1572", "0.05"]
    html = (
        f'<section class="slide"><div>{banner}</div><p>n=1572</p></section>'
        f'<section class="slide"><div>{banner}</div>{svg}</section>'
    )
    assert html_numeric_tokens(html) == {
        "slide-01": ["1572"],
        "slide-02": ["1572", "0.05"],
    }
