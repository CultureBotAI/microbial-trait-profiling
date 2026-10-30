"""The cellulose_decr figures and the slide captions they write (#556, #570, #631).

They need matplotlib from the viz extra; conftest.py leaves this module out when it
is missing, and the analysis tests in test_cellulose_decr.py still run."""

from __future__ import annotations

import hashlib
import importlib
import json
import shutil
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest
from matplotlib.figure import Figure

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
    for output, sha in provenance["outputs"].items():
        # Each output with its hash, so a deck can check it shows these figures (#637)
        assert hashlib.sha256((tmp_path / output).read_bytes()).hexdigest() == sha
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
    assert "HIV 10 : 2 : 2" in direction["takeaway"]
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


def test_crc_statement_names_the_tissue_contrasts_only_when_they_matter(
    results: figures.Results,
) -> None:
    assert "counting the tumour-vs-adjacent-tissue contrasts too" in (
        figures.crc_compartment_statement(results)
    )
    table = results.scenario_interaction.copy()
    primary = figures.interaction_row(results, ORIENTED, "CRC", "study")
    counted = table["scenario"].eq("paired_tissue_counted")
    table.loc[counted, "fisher_p_two_sided"] = primary["fisher_p_two_sided"]
    statement = figures.crc_compartment_statement(replace(results, scenario_interaction=table))
    assert "counting" not in statement
    assert "does not show anatomical specificity" in statement


def test_p_values_keep_two_significant_digits() -> None:
    assert figures.sci(0.0097) == "0.0097"
    assert figures.sci(0.0596) == "0.060"
    # More digits where two would round across 0.05 (#641)
    assert figures.sci(0.0496) == "0.0496"
    assert figures.sci(0.04996) == "0.04996"
    # Not the repr fallback: the loop stops at the first width that stays below 0.05
    assert figures.sci(0.049612) == "0.0496"
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


def test_importing_the_figures_changes_no_matplotlib_setting() -> None:
    """The module used to update rcParams at import, for every later test (#631)."""
    import matplotlib

    before = dict(matplotlib.rcParams)
    importlib.reload(figures)
    assert dict(matplotlib.rcParams) == before
    # The module was imported at collection, so compare with the style too.
    assert matplotlib.rcParams["svg.hashsalt"] != figures.STYLE["svg.hashsalt"]


def test_the_figures_are_drawn_in_their_style(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import matplotlib

    for name in figures.FIGURE_INPUTS:
        shutil.copy(ANALYSIS_DIR / name, tmp_path / name)
    seen = []
    draw = figures.graphical_core

    def recording(results: figures.Results) -> Figure:
        seen.append(matplotlib.rcParams["svg.hashsalt"])
        return draw(results)

    monkeypatch.setattr(figures, "graphical_core", recording)
    figures.main(["--analysis-dir", str(tmp_path), "--code-commit", "0" * 40])
    assert seen == [figures.STYLE["svg.hashsalt"]]


# ----------------------------------------------------------------------- tufte plates & deck (#550)


@pytest.fixture(scope="module")
def generated_tufte_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Generate all Tufte plates, standalone SVGs, and the 16-slide HTML deck once (#550)."""
    from mtp.cellulose_decr import tufte_plates

    out_dir = tmp_path_factory.mktemp("cellulose_decr_tufte")
    assert tufte_plates.main(["--analysis-dir", str(out_dir)]) == 0
    return out_dir


def test_tufte_plates_and_slide_deck_generation(generated_tufte_dir: Path) -> None:
    """The tracked generator writes fig1-fig6 (.svg/.pdf/.png), fig7-fig9 (.svg), PDF, and
    HTML (#550).
    """
    from mtp.cellulose_decr import tufte_plates

    tufte_dir = generated_tufte_dir / "figures_tufte"
    assert tufte_dir.is_dir()
    for stem in tufte_plates.PLATE_STEMS:
        for ext in ("svg", "pdf", "png"):
            path = tufte_dir / f"{stem}.{ext}"
            assert path.is_file(), f"Missing {path.name}"
            assert path.stat().st_size > 0, f"Empty {path.name}"
    for svg_name in tufte_plates.STANDALONE_SVGS:
        path = tufte_dir / svg_name
        assert path.is_file(), f"Missing {svg_name}"
        assert path.stat().st_size > 0, f"Empty {svg_name}"
    combined_pdf = tufte_dir / "cellulose_decr_tufte_figures.pdf"
    assert combined_pdf.is_file()
    assert combined_pdf.stat().st_size > 0

    deck_html = generated_tufte_dir / "cellulose_decr_slide_deck.html"
    assert deck_html.is_file()
    assert deck_html.stat().st_size > 0
    html_text = deck_html.read_text(encoding="utf-8")
    assert html_text.count('<section class="slide') == 16


def test_plate2_species_roles_do_not_overlap_disease_dots(generated_tufte_dir: Path) -> None:
    """Plate II species role annotations must not print over the PD/HIV/CRC dot columns (#550)."""
    from mtp.cellulose_decr import tufte_plates

    svg_text = (
        generated_tufte_dir / "figures_tufte" / "fig2_tufte_core_guild_dot_table.svg"
    ).read_text(encoding="utf-8")
    for redundant in ("(PD &amp; CRC)", "(PD & CRC)", "(CRC)"):
        assert redundant not in svg_text

    fig = tufte_plates.build_plate2_core_guild_dot_table()
    try:
        fig.draw_without_rendering()
        ax_table = fig.axes[0]
        pd_dot_center_px = ax_table.transData.transform((tufte_plates.PLATE2_DOT_X[0], 0.0))[0]
        # Each dot has s=48..52 pt^2 (~4 pt radius = ~5.5 px at 100 dpi); require >= 8 px clearance.
        max_role_right_px = pd_dot_center_px - 8.0
        role_texts = [
            t
            for t in ax_table.texts
            if abs(t.get_position()[0] - tufte_plates.PLATE2_ROLE_X) < 1e-6
            and t.get_text() != "Functional Annotation / Colonic Role"
        ]
        assert len(role_texts) == 14
        for text_obj in role_texts:
            bbox = text_obj.get_window_extent()
            assert bbox.x1 < max_role_right_px, (
                f"Role label {text_obj.get_text()!r} right edge {bbox.x1:.1f}px "
                f"intrudes on PD dot column at {pd_dot_center_px:.1f}px"
            )
    finally:
        import matplotlib.pyplot as plt

        plt.close(fig)


def test_plate3_and_deck_cite_existing_repository_files(generated_tufte_dir: Path) -> None:
    """Plate III and Slide 16 must cite real committed files, never phantom filenames (#550)."""
    import re

    from mtp.paths import REPO_ROOT

    fig3_svg = (generated_tufte_dir / "figures_tufte" / "fig3_tufte_why_only_stool.svg").read_text(
        encoding="utf-8"
    )
    deck_html = (generated_tufte_dir / "cellulose_decr_slide_deck.html").read_text(encoding="utf-8")

    assert "compartment_jaccard_and_direction_stats.tsv" in fig3_svg
    assert "compartment_direction_tests.tsv" not in fig3_svg
    assert "compartment_direction_tests.tsv" not in deck_html
    assert "four_highlighted_pairs_literature_support.tsv" not in deck_html
    assert "cellulose_decr_taxa_presence.tsv" not in deck_html
    assert "cellulose_decreased_taxa_catalog.tsv" in deck_html

    for content in (fig3_svg, deck_html):
        for match in re.findall(
            r"(?:data/analysis/cellulose_decr|methods)/[A-Za-z0-9_./-]+\.(?:tsv|md|json)",
            content,
        ):
            assert (REPO_ROOT / match).is_file(), f"Cited path does not exist: {match}"


def test_plate4_and_deck_cite_nishiwaki_2020_not_2024(generated_tufte_dir: Path) -> None:
    """Plate IV Panel B must cite Nishiwaki 2020 (PMID 32472608), matching Slide 15 (#550)."""
    fig4_svg = (
        generated_tufte_dir / "figures_tufte" / "fig4_tufte_trait_landscape_and_four_pairs.svg"
    ).read_text(encoding="utf-8")
    fig8_svg = (
        generated_tufte_dir / "figures_tufte" / "fig8_literature_support_four_pairs.svg"
    ).read_text(encoding="utf-8")
    deck_html = (generated_tufte_dir / "cellulose_decr_slide_deck.html").read_text(encoding="utf-8")

    assert "Nishiwaki 2020" in fig4_svg
    assert "Nishiwaki 2024" not in fig4_svg
    assert "Nishiwaki et al. (2020)" in fig8_svg
    assert "Nishiwaki 2024" not in fig8_svg
    assert "Nishiwaki et al. (2020)" in deck_html
    assert "Nishiwaki 2024" not in deck_html


def test_slide15_literature_svg_text_fits_inside_boxes(generated_tufte_dir: Path) -> None:
    """Every text line in Slide 15 / fig8 SVG must fit inside its banner or card box (#550)."""
    import re
    import xml.etree.ElementTree as ET

    fig8_path = generated_tufte_dir / "figures_tufte" / "fig8_literature_support_four_pairs.svg"
    fig8_svg = fig8_path.read_text(encoding="utf-8")
    deck_html = (generated_tufte_dir / "cellulose_decr_slide_deck.html").read_text(encoding="utf-8")
    assert fig8_svg.strip() in deck_html

    root = ET.fromstring(fig8_svg)
    ns = {"svg": "http://www.w3.org/2000/svg"}

    def _est_width(text: str, font_size: float, bold: bool = False) -> float:
        factor = 0.56 if bold else 0.52
        return len(text) * font_size * factor

    # Check top-level banner texts (y < 68 spans full 1485px banner; y >= 68 sits in 680px boxes)
    for elem in root.findall("svg:text", ns):
        x = float(elem.attrib["x"])
        y = float(elem.attrib["y"])
        fs = float(elem.attrib.get("font-size", "12"))
        bold = elem.attrib.get("font-weight") == "bold"
        content = "".join(elem.itertext()).strip()
        box_right = 1485.0 if y < 68.0 else 715.0 if x < 750.0 else 1465.0
        est_right = x + _est_width(content, fs, bold)
        assert est_right <= box_right - 6.0, (
            f"Banner text overflows {box_right}: x={x}, est_right={est_right:.1f}, text={content!r}"
        )

    # Check all 4 cards (each card is a <g transform="translate(tx, ty)"> with width=725)
    card_groups = root.findall("svg:g", ns)
    assert len(card_groups) == 4
    for card in card_groups:
        texts = card.findall("svg:text", ns)
        title_elem = texts[0]
        badge_elem = texts[1]
        title_fs = float(title_elem.attrib.get("font-size", "14"))
        badge_fs = float(badge_elem.attrib.get("font-size", "12"))
        title_right = float(title_elem.attrib["x"]) + _est_width(
            "".join(title_elem.itertext()).strip(), title_fs, bold=True
        )
        badge_left = float(badge_elem.attrib["x"]) - _est_width(
            "".join(badge_elem.itertext()).strip(), badge_fs, bold=True
        )
        assert badge_left - title_right >= 15.0, (
            f"Card header title collides with badge: title_right={title_right:.1f}, "
            f"badge_left={badge_left:.1f}"
        )
        for elem in texts[2:]:
            if elem.attrib.get("text-anchor") == "end":
                continue
            x = float(elem.attrib.get("x", "18"))
            fs = float(elem.attrib.get("font-size", "11.5"))
            bold = elem.attrib.get("font-weight") == "bold"
            line = "".join(elem.itertext()).strip()
            est_right = x + _est_width(line, fs, bold)
            assert est_right <= 725.0 - 8.0, (
                f"Card line overflows 725px card: est_right={est_right:.1f}, line={line!r}"
            )
            assert not re.search(r"Nishiwaki\s+2024", line)


def test_tufte_plates_and_deck_are_byte_repeatable(
    generated_tufte_dir: Path, tmp_path: Path
) -> None:
    """Rerunning tufte_plates produces identical bytes across all 23 output files (#1127)."""
    from mtp.cellulose_decr import tufte_plates

    second_dir = tmp_path / "rerun"
    assert tufte_plates.main(["--analysis-dir", str(second_dir), "--code-commit", "deadbeef"]) == 0

    first_files = sorted(
        p.relative_to(generated_tufte_dir) for p in generated_tufte_dir.rglob("*") if p.is_file()
    )
    second_files = sorted(p.relative_to(second_dir) for p in second_dir.rglob("*") if p.is_file())
    assert first_files == second_files
    assert len(first_files) == 6 * 3 + 3 + 1 + 1
    for rel_path in first_files:
        assert (generated_tufte_dir / rel_path).read_bytes() == (
            second_dir / rel_path
        ).read_bytes(), f"Non-repeatable bytes in {rel_path}"


def test_slide_deck_image_paths_match_current_figure_stems(generated_tufte_dir: Path) -> None:
    """Every <img src=...> in the 16-slide deck resolves to a current figure (#1127)."""
    import re

    from mtp.cellulose_decr import tufte_plates

    deck_html = (generated_tufte_dir / "cellulose_decr_slide_deck.html").read_text(encoding="utf-8")
    img_srcs = re.findall(r'<img\s+src="([^"]+)"', deck_html)
    assert len(img_srcs) == 13

    expected_tufte_svgs = {f"figures_tufte/{stem}.svg" for stem in tufte_plates.PLATE_STEMS}
    expected_graphical_svgs = {
        "figures_graphical/fig1_overlap_and_null_models.svg",
        "figures_graphical/fig2_shared_taxa_and_set_definitions.svg",
        "figures_graphical/fig3_direction_by_body_site.svg",
        "figures_graphical/fig4_trait_overlap_ranking.svg",
    }
    expected_statistical_svgs = {
        "figures/fig1_pairwise_and_3way_jaccard_nulls.svg",
        "figures/fig2_direction_by_body_site_and_definition.svg",
        "figures/fig3_core_taxa_and_trait_benchmark.svg",
    }
    assert set(img_srcs) == (
        expected_tufte_svgs | expected_graphical_svgs | expected_statistical_svgs
    )
    for src in expected_tufte_svgs:
        assert (generated_tufte_dir / src).is_file()


def test_tufte_builders_apply_rc_context_without_mutating_global_rcparams() -> None:
    """Standalone build_plate* calls apply _tufte_rcparams() without mutating rcParams (#1127)."""
    import matplotlib
    import matplotlib.pyplot as plt

    from mtp.cellulose_decr import tufte_plates

    before = dict(matplotlib.rcParams)
    importlib.reload(tufte_plates)
    assert dict(matplotlib.rcParams) == before

    fig = tufte_plates.build_plate1_convergence_slope_and_overlap()
    try:
        assert dict(matplotlib.rcParams) == before
        # Text artists constructed inside build_plate1 must carry the Tufte serif stack
        assert fig.texts[0].get_fontfamily() == ["serif"]
    finally:
        plt.close(fig)


def test_tufte_options_are_never_abbreviated() -> None:
    """tufte_plates.parse_args rejects abbreviated flags so Justfile guards hold (#617, #1127)."""
    from mtp.cellulose_decr import tufte_plates

    with pytest.raises(SystemExit) as stopped:
        tufte_plates.parse_args(["--analysis", "/tmp/x"])
    assert stopped.value.code == 2
    with pytest.raises(SystemExit) as stopped:
        tufte_plates.parse_args(["--code", "abc"])
    assert stopped.value.code == 2
