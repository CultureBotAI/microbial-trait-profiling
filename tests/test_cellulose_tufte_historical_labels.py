"""Pin the visible historical and correction labels without rerunning analyses (#1416)."""

from __future__ import annotations

import ast
import html
import re
import xml.etree.ElementTree as ET
from collections.abc import Iterator
from pathlib import Path
from typing import Any

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


def _assigned_expression(name: str, statements: list[ast.stmt]) -> ast.expr | None:
    values: list[ast.expr] = []
    for statement in statements:
        targets = (
            statement.targets
            if isinstance(statement, ast.Assign)
            else [statement.target]
            if isinstance(statement, ast.AnnAssign)
            else []
        )
        if any(isinstance(target, ast.Name) and target.id == name for target in targets):
            assert isinstance(statement, ast.Assign | ast.AnnAssign)
            assert statement.value is not None
            values.append(statement.value)
    assert len(values) <= 1, f"Ambiguous route binding: {name}"
    return values[0] if values else None


def _saved_image_suffixes(tree: ast.Module) -> set[str]:
    saves = [
        node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "save"
    ]
    assert len(saves) == 1
    save = saves[0]
    returns = [node.value for node in save.body if isinstance(node, ast.Return)]
    assert len(returns) == 1
    assert isinstance(returns[0], ast.Name)
    written = returns[0].id
    expression = _assigned_expression(written, save.body)
    assert isinstance(expression, ast.ListComp)
    assert len(expression.generators) == 1
    generator = expression.generators[0]
    assert not generator.ifs
    assert not generator.is_async
    assert isinstance(generator.target, ast.Name)
    iterable: ast.expr | None = generator.iter
    if isinstance(iterable, ast.Name):
        name = iterable.id
        iterable = _assigned_expression(name, save.body)
        if iterable is None:
            iterable = _assigned_expression(name, tree.body)
    assert isinstance(iterable, ast.Tuple | ast.List)
    suffixes = []
    for suffix in iterable.elts:
        assert isinstance(suffix, ast.Constant)
        assert isinstance(suffix.value, str)
        assert suffix.value.startswith(".")
        suffixes.append(suffix.value)
    assert suffixes
    assert len(set(suffixes)) == len(suffixes)
    # Only claim routes for the supported directory / f"{stem}{suffix}" form,
    # and require an actual fig.savefig for every returned path index.
    args = save.args.args
    assert len(args) >= 3
    value = expression.elt
    assert isinstance(value, ast.BinOp)
    assert isinstance(value.op, ast.Div)
    assert isinstance(value.left, ast.Name)
    assert value.left.id == args[1].arg
    assert isinstance(value.right, ast.JoinedStr)
    names = []
    for part in value.right.values:
        assert isinstance(part, ast.FormattedValue)
        assert isinstance(part.value, ast.Name)
        assert part.conversion == -1
        assert part.format_spec is None
        names.append(part.value.id)
    assert names == [args[2].arg, generator.target.id]
    indices = []
    for node in ast.walk(save):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == args[0].arg
            and node.func.attr == "savefig"
        ):
            assert node.args
            assert isinstance(node.args[0], ast.Subscript)
            path = node.args[0]
            assert isinstance(path.value, ast.Name)
            assert path.value.id == written
            assert isinstance(path.slice, ast.Constant)
            assert type(path.slice.value) is int
            indices.append(path.slice.value)
    assert sorted(indices) == list(range(len(suffixes)))
    return set(suffixes)


def _producer_image_routes() -> set[str]:
    """Resolve the table used by main and suffixes actually written by save.

    This bounded AST reader supports a single literal table or suffix sequence
    assigned locally or at module scope, under any name. Dynamic table/suffix
    expressions are unsupported. This checks literal routes, not arbitrary
    control/data flow through save call arguments. No producer loading or
    rendering is run.
    """
    producer = Path(tufte_plates.__file__).with_name("figures.py")
    assert producer.stat().st_size <= 262_144
    source = producer.read_text(encoding="utf-8")
    assert len(source.encode("utf-8")) <= 262_144
    tree = ast.parse(source, filename=str(producer))
    mains = [
        node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "main"
    ]
    assert len(mains) == 1
    main = mains[0]
    loops = [
        node
        for node in main.body
        if isinstance(node, ast.For)
        and isinstance(node.iter, ast.Call)
        and isinstance(node.iter.func, ast.Attribute)
        and node.iter.func.attr == "items"
        and any(
            isinstance(call, ast.Call)
            and isinstance(call.func, ast.Name)
            and call.func.id == "save"
            for call in ast.walk(node)
        )
    ]
    assert len(loops) == 1
    loop = loops[0]
    assert isinstance(loop.iter, ast.Call)
    assert isinstance(loop.iter.func, ast.Attribute)
    assert not loop.iter.args
    assert not loop.iter.keywords
    reference = loop.iter.func.value
    assert isinstance(reference, ast.Name)
    table = _assigned_expression(reference.id, main.body[: main.body.index(loop)])
    if table is None:
        table = _assigned_expression(reference.id, tree.body[: tree.body.index(main)])
    assert isinstance(table, ast.Dict)
    suffixes = _saved_image_suffixes(tree)
    routes: set[str] = set()
    for folder, entries in zip(table.keys, table.values, strict=True):
        assert isinstance(folder, ast.Constant)
        assert isinstance(folder.value, str)
        assert isinstance(entries, ast.List | ast.Tuple)
        for entry in entries.elts:
            assert isinstance(entry, ast.Tuple)
            assert len(entry.elts) == 3
            stem = entry.elts[0]
            assert isinstance(stem, ast.Constant)
            assert isinstance(stem.value, str)
            for suffix in suffixes:
                route = f"{folder.value}/{stem.value}{suffix}"
                assert route not in routes, f"Duplicate producer output: {route}"
                routes.add(route)
    return routes


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
        "BH q: family not recorded (#1525); 4 nominal tests shown"
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
        "BH q: family not recorded (#1525); p-values nominal (14 tests shown); "
        "benchmark/percentile: historical (#1416)."
    ) in texts
    assert "46 studies · Δ = \u221227.4 pp (q = 0.0077)" in texts


def test_fig_eight_pins_the_title_family_note() -> None:
    texts = _svg_texts(tufte_plates.render_fig8_svg())
    assert texts[0] == (
        "OPPOSITE STOOL TRAIT DIRECTIONS IN HUMAN SIGNATURES (BugSigDB \u00d7 KG-Microbe) "
        "· historical (#1416); BH q: family not recorded (#1525); 8 nominal tests (2/card)"
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
        "only ~20% would overlap between diseases. Observed overlap is 50.4% "
        "(4 nominal tests shown:\n"
        "exact hypergeometric p < 10⁻¹⁰ per pair, 3 tests; "
        "random-guild permutation p = 0.011, 1 test)."
    ]


def test_plates_two_and_three_state_nominal_test_counts(
    label_plates: dict[str, Figure],
) -> None:
    """Plates II and III state how many nominal tests their panels show (#1533)."""
    p2_source = _plain(label_plates["II"].texts[3].get_text())
    assert p2_source.endswith("PD\u2013CRC exact nominal p = 5.98 \u00d7 10⁻¹⁰ (1/3 tests).")
    p2_headers = [
        text
        for text in _artist_texts(label_plates["II"])
        if text.startswith("NINE CORE GENERA REPORTED DECREASED IN ALL THREE DISEASES")
    ]
    assert p2_headers == [
        "NINE CORE GENERA REPORTED DECREASED IN ALL THREE DISEASES "
        "(17.1\u00d7 over random expectation, nominal p = 3.5 \u00d7 10⁻¹², 1 core test)"
    ]
    p3_subtitle = next(
        text.get_text()
        for text in label_plates["III"].axes[0].texts
        if text.get_text().startswith("Comparing stool against colon wall biopsies")
    )
    assert _plain(p3_subtitle) == (
        "Comparing stool against colon wall biopsies and saliva relative to a 50% "
        "binomial baseline (7 nominal binomial sign tests)"
    )


@pytest.mark.parametrize("plate", ["I", "II", "III", "IV", "V", "VI"])
def test_new_plate_labels_fit_the_fixed_packet_canvas(
    label_plates: dict[str, Figure], plate: str
) -> None:
    fig = label_plates[plate]
    fig.draw_without_rendering()
    artists = list(fig.texts) + [text for ax in fig.axes for text in ax.texts]
    ticks = [
        tick
        for ax in fig.axes
        for tick in (*ax.get_xticklabels(), *ax.get_yticklabels())
        if tick.get_text()
    ]
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
        for tick in ticks:
            assert not bounds.overlaps(tick.get_window_extent()), (
                text.get_text(),
                tick.get_text(),
            )
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
    """Bind each slide's expected image to figures.main's actual routes (#1556).

    The table-driven generator and its outputs are documented in
    methods/cellulose_decr.md (the figures section) and the analysis README.
    The literals keep the slide-to-figure identity pinned; the independent AST
    read refuses a table-entry rename that would leave those literals stale.
    No producer table loading, calculations or rendering are executed.
    """
    assert image_src in _producer_image_routes()
    section = deck_sections[slide_number - 1]
    assert re.findall(r'<img\s+src="([^"]+)"', section) == [image_src]
    subtitle = re.search(r'<p class="subtitle"[^>]*>(.*?)</p>', section, re.DOTALL)
    assert subtitle is not None
    assert "; figure drawn from current tables" in _plain(subtitle.group(1))


def test_current_table_route_guard_rejects_a_renamed_producer_output(
    deck_sections: list[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    """A current-table label must not survive a producer-only output rename (#1556)."""
    producer = Path(tufte_plates.__file__).with_name("figures.py")
    original = producer.read_text(encoding="utf-8")
    old_stem = "fig1_overlap_and_null_models"
    assert original.count(f'"{old_stem}"') == 1
    changed = original.replace(f'"{old_stem}"', '"renamed_overlap_output"')
    read_text = Path.read_text

    def changed_producer(path: Path, *args: Any, **kwargs: Any) -> str:
        return changed if path == producer else read_text(path, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", changed_producer)
    with pytest.raises(AssertionError):
        test_current_table_image_labels_are_scoped_to_the_correct_output_route(
            deck_sections, 8, f"figures_graphical/{old_stem}.svg"
        )


def test_deck_four_pair_footer_identifies_the_unknown_bh_family(deck_sections: list[str]) -> None:
    footer = re.search(r'<footer class="slide-footer">(.*?)</footer>', deck_sections[4], re.DOTALL)
    assert footer is not None
    pair_note = _plain(footer.group(1)).split("Four highlighted pairs:", 1)[1]
    assert "BH q: family not recorded (#1525)" in pair_note
    assert "(1 nominal test shown)" in pair_note


@pytest.mark.parametrize(
    ("slide_number", "expected_counts"),
    [
        (1, ["(1 nominal test shown)", "(2 nominal tests shown)", "(1 nominal test shown)"]),
        (2, ["(1 nominal test shown)"]),
        (3, ["(2 nominal tests shown)"]),
        (4, ["(3 nominal tests shown)"]),
        (5, ["(1 nominal test shown)"]),
        (8, ["(2 nominal tests shown)"]),
        (9, ["(1 nominal test shown)", "(1 nominal test shown)"]),
        (12, ["(3 nominal tests shown)"]),
        (13, ["(6 nominal tests shown)"]),
    ],
)
def test_deck_cards_and_footers_state_nominal_test_counts(
    deck_sections: list[str], slide_number: int, expected_counts: list[str]
) -> None:
    """Every deck cover card and footer with nominal p states its test count (#1533)."""
    section = deck_sections[slide_number - 1]
    blocks = re.findall(r'<div class="(?:lbl|footer-col)">(.*?)</div>', section, re.DOTALL)
    with_p = [
        _plain(b)
        for b in blocks
        if re.search(r"\bnominal\s+(?:<i>)?p(?:</i>)?\s*(?:=|≤|≥|<|>|&lt;|&gt;)", b)
    ]
    assert len(with_p) == len(expected_counts)
    for block_text, expected in zip(with_p, expected_counts, strict=True):
        assert expected in block_text, (slide_number, block_text)
        assert ") (" not in block_text, (slide_number, block_text)
        if slide_number > 1:
            assert block_text.endswith(f"{expected}."), (slide_number, block_text)


def test_deck_slides_six_and_seven_footers_carry_historical_run_label(
    deck_sections: list[str],
) -> None:
    """Slides 6 and 7 footers quoting 317, 66,240, and 1,912 carry their run label (#1533)."""
    s6_cols = [
        _plain(c)
        for c in re.findall(r'<div class="footer-col">(.*?)</div>', deck_sections[5], re.DOTALL)
    ]
    assert len(s6_cols) == 2
    assert s6_cols[0].startswith("Three statistical arms \u2014 historical (#1416):")
    assert "317 non-MediaDive traits" in s6_cols[0]
    assert "66,240 enrichments" in s6_cols[0]
    assert s6_cols[1].startswith("Table 1 fraction and ratio tests:")
    assert "historical (#1416)" not in s6_cols[1]

    s7_cols = [
        _plain(c)
        for c in re.findall(r'<div class="footer-col">(.*?)</div>', deck_sections[6], re.DOTALL)
    ]
    assert len(s7_cols) == 2
    assert s7_cols[0].startswith("Shared comparison basis \u2014 historical (#1416):")
    assert "317 non-MediaDive trait denominator" in s7_cols[0]
    assert s7_cols[1].startswith("Integrated evidence contract \u2014 historical (#1416):")
    assert "1,912 rows" in s7_cols[1]

    synthetic = (
        '<section class="slide">'
        '<p class="subtitle">Example · historical (#1416)</p>'
        '<footer class="slide-footer">'
        '<div class="footer-col"><b>Unrelated decimal:</b> r = 0.317 across 1,317 pairs.</div>'
        '<div class="footer-col"><b>Screen basis:</b> 317 non-MediaDive traits.</div>'
        "</footer></section>"
    )
    rewritten = tufte_plates._annotate_deck_nominal_counts_and_footer_labels(synthetic)
    syn_cols = [
        _plain(c) for c in re.findall(r'<div class="footer-col">(.*?)</div>', rewritten, re.DOTALL)
    ]
    assert syn_cols == [
        "Unrelated decimal: r = 0.317 across 1,317 pairs.",
        "Screen basis \u2014 historical (#1416): 317 non-MediaDive traits.",
    ]


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


@pytest.mark.parametrize(
    "change",
    [
        "reorder",
        "module-table",
        "rename-local",
        "suffix",
        "save-optional-dpi",
        "module-suffixes",
        "local-suffixes",
    ],
)
def test_route_guard_tracks_outputs_not_table_layout_or_assumed_svg(
    change: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Independent syntax fixture: the real producer may itself already be
    # reordered, renamed, or hoisted without breaking construction of this test.
    tree = ast.parse("""
def save(fig, directory, stem, book):
    paths = [directory / f"{stem}{suffix}" for suffix in (".png", ".svg")]
    fig.savefig(paths[0])
    fig.savefig(paths[1])
    return paths

def main():
    plates = {"plots": [("first", draw, caption), ("second", draw, caption)]}
    for folder, entries in plates.items():
        for stem, draw, caption in entries:
            save(draw(), root / folder, stem, book)
""")
    main = next(
        node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "main"
    )
    save = next(
        node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "save"
    )
    table = main.body[0]
    assert isinstance(table, ast.Assign)
    assert isinstance(table.value, ast.Dict)
    expected_suffixes: tuple[str, ...] = (".png", ".svg")
    if change == "reorder":
        entries = table.value.values[0]
        assert isinstance(entries, ast.List)
        entries.elts.reverse()
    elif change in {"module-table", "rename-local"}:

        class Rename(ast.NodeTransformer):
            def visit_Name(self, node: ast.Name) -> ast.AST:
                if node.id == "plates":
                    node.id = "OUTPUT_ROUTES"
                return node

        Rename().visit(main)
        if change == "module-table":
            main.body.remove(table)
            tree.body.insert(tree.body.index(main), table)
    elif change == "save-optional-dpi":
        save.args.args.append(ast.arg(arg="dpi", annotation=ast.Name(id="int", ctx=ast.Load())))
        save.args.defaults.append(ast.Name(id="PNG_DPI", ctx=ast.Load()))
    elif change in {"module-suffixes", "local-suffixes"}:
        binding = ast.parse('FIGURE_SUFFIXES = (".svg", ".pdf", ".png")').body[0]
        paths_assign = save.body[0]
        assert isinstance(paths_assign, ast.Assign)
        assert isinstance(paths_assign.value, ast.ListComp)
        paths_assign.value.generators[0].iter = ast.Name(id="FIGURE_SUFFIXES", ctx=ast.Load())
        save.body.insert(-1, ast.parse("fig.savefig(paths[2])").body[0])
        if change == "module-suffixes":
            tree.body.insert(0, binding)
        else:
            save.body.insert(0, binding)
        expected_suffixes = (".svg", ".pdf", ".png")
    else:
        suffixes = [
            node
            for node in ast.walk(save)
            if isinstance(node, ast.Constant) and node.value == ".svg"
        ]
        assert len(suffixes) == 1
        suffixes[0].value = ".svgz"
        expected_suffixes = (".png", ".svgz")
    changed = ast.unparse(tree)
    producer = Path(tufte_plates.__file__).with_name("figures.py")
    read_text = Path.read_text

    def changed_producer(path: Path, *args: Any, **kwargs: Any) -> str:
        return changed if path == producer else read_text(path, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", changed_producer)
    assert _producer_image_routes() == {
        f"plots/{stem}{ext}" for stem in ("first", "second") for ext in expected_suffixes
    }


@pytest.mark.parametrize(
    "invalid",
    ["two-args", "unbound-suffixes", "non-sequence-suffixes", "duplicate-suffixes"],
)
def test_route_guard_refuses_underarity_save_or_invalid_suffix_binding(invalid: str) -> None:
    if invalid == "two-args":
        source = """
def save(fig, directory):
    paths = [directory / f"{stem}{suffix}" for suffix in (".png", ".svg")]
    fig.savefig(paths[0])
    fig.savefig(paths[1])
    return paths
"""
    elif invalid == "unbound-suffixes":
        source = """
def save(fig, directory, stem, book):
    paths = [directory / f"{stem}{suffix}" for suffix in UNKNOWN_SUFFIXES]
    fig.savefig(paths[0])
    fig.savefig(paths[1])
    return paths
"""
    elif invalid == "non-sequence-suffixes":
        source = """
FIGURE_SUFFIXES = ".svg"
def save(fig, directory, stem, book):
    paths = [directory / f"{stem}{suffix}" for suffix in FIGURE_SUFFIXES]
    fig.savefig(paths[0])
    return paths
"""
    else:
        source = """
FIGURE_SUFFIXES = (".svg", ".png")
FIGURE_SUFFIXES = (".svg", ".pdf", ".png")
def save(fig, directory, stem, book):
    paths = [directory / f"{stem}{suffix}" for suffix in FIGURE_SUFFIXES]
    fig.savefig(paths[0])
    fig.savefig(paths[1])
    return paths
"""
    tree = ast.parse(source)
    try:
        _saved_image_suffixes(tree)
    except AssertionError:
        return
    except Exception as exc:
        raise AssertionError(
            f"Expected AssertionError for {invalid}, got {type(exc).__name__}"
        ) from exc
    raise AssertionError(f"Route guard accepted invalid save AST: {invalid}")
