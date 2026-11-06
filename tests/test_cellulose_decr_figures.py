"""The cellulose_decr figures and the slide captions they write (#556, #570, #631).

They need matplotlib from the viz extra; conftest.py leaves this module out when it
is missing, and the analysis tests in test_cellulose_decr.py still run."""

from __future__ import annotations

import hashlib
import importlib
import inspect
import json
import math
import re
import shutil
import xml.etree.ElementTree as ET
from collections.abc import Callable
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Literal, cast

import matplotlib
import pandas as pd
import pytest
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.text import Text
from PIL import Image, ImageChops

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

RETIREMENT_BANNER = "RETIRED: historical, unverified numbers (#1572)"
HISTORICAL_PLATE_STEMS = (
    "fig1_tufte_convergence_slope_and_overlap",
    "fig2_tufte_core_guild_dot_table",
    "fig3_tufte_why_only_stool",
    "fig4_tufte_trait_landscape_and_four_pairs",
    "fig5_tufte_three_statistical_arms",
    "fig6_tufte_combining_results_integrated_evidence",
)
HISTORICAL_SVG_STEMS = (
    *HISTORICAL_PLATE_STEMS,
    "fig7_tufte_convergence_network_and_venn",
    "fig8_literature_support_four_pairs",
    "fig9_set_theory_17_core_vs_pairwise",
)


@pytest.fixture(scope="module", autouse=True)
def retirement_import_origin() -> None:
    from mtp.cellulose_decr import tufte_plates

    assert Path(tufte_plates.__file__).resolve() == (
        Path(__file__).resolve().parents[1] / "src/mtp/cellulose_decr/tufte_plates.py"
    )


@dataclass(frozen=True)
class RetirementSave:
    kind: Literal["standalone", "packet", "packet_figure"]
    name: str
    format: str
    text: str
    visible: bool
    bounds: tuple[float, float, float, float]
    canvas: tuple[float, float, float, float]


@dataclass(frozen=True)
class RetirementExport:
    directory: Path
    saves: tuple[RetirementSave, ...]
    png_without_banner: dict[str, Path]


def _retirement_artist(figure: Figure) -> Text:
    banners = [text for text in figure.texts if text.get_text() == RETIREMENT_BANNER]
    assert len(banners) == 1, "A real save must receive exactly one retirement banner"
    banner = banners[0]
    assert banner.get_text().strip(), "Banner empty before a real save"
    assert banner.get_visible(), "Banner hidden before a real save"
    foreground = to_rgba(banner.get_color(), alpha=banner.get_alpha())
    assert foreground[3] > 0, "Banner has no opaque paint before a real save"
    assert foreground[:3] != figure.get_facecolor()[:3], "Banner paint matches its canvas"
    return banner


def _assert_retirement_save_crop(kwargs: dict[str, Any]) -> None:
    """Pin effective Matplotlib defaults and explicit kwargs before a real export."""
    bbox = kwargs.get("bbox_inches")
    if bbox is None:
        bbox = matplotlib.rcParams["savefig.bbox"]
    if bbox is not None:
        assert isinstance(bbox, str), "Historical export uses an unverified explicit crop"
        assert bbox == "tight", "Historical export must preserve the canvas or tight bounds"
    pad = kwargs.get("pad_inches")
    if pad is None:
        pad = matplotlib.rcParams["savefig.pad_inches"]
    assert isinstance(pad, (int, float)), "Historical export padding must be numeric"
    assert not isinstance(pad, bool)
    assert math.isfinite(pad), "Historical export padding must be finite"
    assert pad >= 0, "Negative export padding can crop away the retirement banner"


def _retirement_save(
    figure: Figure,
    kind: Literal["standalone", "packet", "packet_figure"],
    name: str,
    format_name: str,
) -> RetirementSave:
    banner = _retirement_artist(figure)
    figure.draw_without_rendering()
    bounds = banner.get_window_extent()
    assert bounds.width > 0, f"{name}: empty banner width"
    assert bounds.height > 0, f"{name}: empty banner height"
    assert figure.bbox.contains(bounds.x0, bounds.y0), f"{name}: banner starts outside canvas"
    assert figure.bbox.contains(bounds.x1, bounds.y1), f"{name}: banner ends outside canvas"
    return RetirementSave(
        kind,
        name,
        format_name,
        banner.get_text(),
        banner.get_visible(),
        (bounds.x0, bounds.y0, bounds.x1, bounds.y1),
        (figure.bbox.x0, figure.bbox.y0, figure.bbox.x1, figure.bbox.y1),
    )


@pytest.fixture(scope="module")
def retirement_export(tmp_path_factory: pytest.TempPathFactory) -> RetirementExport:
    """Call the real writer once, checking each save and retaining PNG controls (#1706)."""
    from mtp.cellulose_decr import tufte_plates

    out_dir = tmp_path_factory.mktemp("cellulose_decr_tufte")
    reference_dir = tmp_path_factory.mktemp("cellulose_decr_without_banner")
    saves: list[RetirementSave] = []
    png_without_banner: dict[str, Path] = {}
    # Matplotlib's PDF backend passes PdfPages through Figure.savefig internally;
    # its public filename annotation omits that supported internal destination.
    figure_savefig = cast(Callable[..., None], Figure.savefig)
    packet_savefig = cast(Callable[..., None], PdfPages.savefig)
    pagecount = cast(Callable[[PdfPages], int], PdfPages.get_pagecount)

    def checked_figure_savefig(
        figure: Figure,
        destination: Any,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        _assert_retirement_save_crop(kwargs)
        if isinstance(destination, PdfPages):
            saves.append(_retirement_save(figure, "packet_figure", "packet page", "pdf"))
            figure_savefig(figure, destination, *args, **kwargs)
            return
        assert isinstance(destination, Path), "Unexpected destination from the historical writer"
        format_name = destination.suffix.removeprefix(".")
        saves.append(_retirement_save(figure, "standalone", destination.name, format_name))
        figure_savefig(figure, destination, *args, **kwargs)
        if format_name != "png":
            return

        # Check before changing visibility: never repair an already-hidden banner.
        banner = _retirement_artist(figure)
        original_visibility, original_dpi = banner.get_visible(), figure.dpi
        try:
            figure.set_dpi(float(kwargs["dpi"]))
            figure.draw_without_rendering()
            # Freeze the visible export's tight crop before removing its text.
            # Otherwise tight cropping would move the original kicker into this band.
            crop = figure.get_tightbbox().padded(float(matplotlib.rcParams["savefig.pad_inches"]))
            reference_kwargs = {**kwargs, "bbox_inches": crop, "pad_inches": 0}
            reference = reference_dir / destination.name
            banner.set_visible(False)
            figure_savefig(figure, reference, *args, **reference_kwargs)
            png_without_banner[destination.stem] = reference
        finally:
            banner.set_visible(original_visibility)
            figure.set_dpi(original_dpi)

    def checked_packet_savefig(packet: PdfPages, figure: Figure, **kwargs: Any) -> None:
        _assert_retirement_save_crop(kwargs)
        before = pagecount(packet)
        page = before + 1
        saves.append(_retirement_save(figure, "packet", f"packet page {page}", "pdf"))
        packet_savefig(packet, figure, **kwargs)
        assert pagecount(packet) == page, f"Packet page {page} was not written"

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(Figure, "savefig", checked_figure_savefig)
        patch.setattr(PdfPages, "savefig", checked_packet_savefig)
        assert tufte_plates.main(["--audit-retired", "--analysis-dir", str(out_dir)]) == 0
    return RetirementExport(out_dir, tuple(saves), png_without_banner)


@pytest.fixture(scope="module")
def generated_tufte_dir(retirement_export: RetirementExport) -> Path:
    """Share the one real 23-file export with the existing artifact tests."""
    return retirement_export.directory


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


def test_written_outputs_check_banner_at_every_real_save(
    retirement_export: RetirementExport,
) -> None:
    """All 18 standalone saves and six packet pages reach the actual writers.

    This checks PDF writer inputs and complete saved files. Extracting text from
    PDF bytes is deliberately outside this test's dependency-free scope (#1706).
    """
    standalone = [save for save in retirement_export.saves if save.kind == "standalone"]
    packet = [save for save in retirement_export.saves if save.kind == "packet"]
    packet_figures = [save for save in retirement_export.saves if save.kind == "packet_figure"]
    assert len(standalone) == 18
    assert {save.name for save in standalone} == {
        f"{stem}.{extension}"
        for stem in HISTORICAL_PLATE_STEMS
        for extension in ("svg", "pdf", "png")
    }
    for extension in ("svg", "pdf", "png"):
        assert sum(save.format == extension for save in standalone) == 6
    assert [save.name for save in packet] == [f"packet page {page}" for page in range(1, 7)]
    assert len(packet_figures) == 6
    assert all(save.format == "pdf" for save in (*packet, *packet_figures))
    for save in retirement_export.saves:
        assert save.text == RETIREMENT_BANNER, save.name
        assert save.visible, save.name
        x0, y0, x1, y1 = save.bounds
        left, bottom, right, top = save.canvas
        assert left <= x0 < x1 <= right, save.name
        assert bottom <= y0 < y1 <= top, save.name
    pdf_names = [f"{stem}.pdf" for stem in HISTORICAL_PLATE_STEMS]
    pdf_names.append("cellulose_decr_tufte_figures.pdf")
    for name in pdf_names:
        data = (retirement_export.directory / "figures_tufte" / name).read_bytes()
        assert len(data) > 100, f"{name}: missing PDF content"
        assert data.startswith(b"%PDF-"), f"{name}: missing PDF header"
        assert data.rstrip().endswith(b"%%EOF"), f"{name}: incomplete PDF file"


@pytest.mark.parametrize("stem", HISTORICAL_SVG_STEMS)
def test_written_svg_has_retirement_text(generated_tufte_dir: Path, stem: str) -> None:
    root = ET.fromstring((generated_tufte_dir / "figures_tufte" / f"{stem}.svg").read_text())
    texts = [" ".join("".join(node.itertext()).split()) for node in root.findall(".//{*}text")]
    assert texts.count(RETIREMENT_BANNER) == 1, stem


def _ink_mask(image: Image.Image) -> Image.Image:
    rgb = image.convert("RGB")
    background = Image.new("RGB", rgb.size, "#fbfaf7")
    difference = ImageChops.difference(rgb, background)
    return difference.convert("L").point(lambda value: 255 if value > 32 else 0)


def _first_ink_band_bounds(image: Image.Image) -> tuple[int, int, int, int]:
    """Locate the saved image's first text line, retaining its canvas coordinates."""
    mask = _ink_mask(image)
    rows = mask.getprojection()[1]
    top = next((index for index, occupied in enumerate(rows) if occupied), None)
    assert top is not None, "Saved image has no visible ink"
    bottom = next((index for index in range(top, len(rows)) if not rows[index]), len(rows))
    band = mask.crop((0, top, mask.width, bottom))
    bounds = band.getbbox()
    assert bounds is not None
    return (bounds[0], top + bounds[1], bounds[2], top + bounds[3])


@pytest.mark.parametrize("stem", HISTORICAL_PLATE_STEMS)
def test_written_png_pixels_contain_the_retirement_band(
    retirement_export: RetirementExport,
    stem: str,
) -> None:
    """The saved PNG has ink where the same saved band without the banner is blank."""
    path = retirement_export.directory / "figures_tufte" / f"{stem}.png"
    with (
        Image.open(path) as saved,
        Image.open(retirement_export.png_without_banner[stem]) as control,
    ):
        assert saved.size == control.size, "The no-banner reference must preserve the actual crop"
        bounds = _first_ink_band_bounds(saved)
        visible_band = saved.crop(bounds).convert("RGB")
        removed_band = control.crop(bounds).convert("RGB")
    visible_ink = _ink_mask(visible_band)
    assert sum(visible_ink.histogram()[1:]) > 100, f"{stem}: saved banner has no substantial ink"
    assert _ink_mask(removed_band).getbbox() is None, f"{stem}: control band contains other content"
    assert ImageChops.difference(visible_band, removed_band).getbbox() is not None, stem


def _written_css_properties(html_text: str, selectors: set[str]) -> dict[str, str]:
    """Resolve known flat class rules by importance, class specificity and source order."""
    css = "\n".join(re.findall(r"<style[^>]*>(.*?)</style>", html_text, re.DOTALL))
    properties: dict[str, str] = {}
    priorities: dict[str, tuple[bool, int, int]] = {}
    order = 0
    for selector_list, body in re.findall(r"([^{}]+)\{([^{}]*)\}", css):
        matching = [
            selector.strip()
            for selector in selector_list.split(",")
            if selector.strip() in selectors
        ]
        if not matching:
            continue
        specificity = max(selector.count(".") for selector in matching)
        for declaration in body.split(";"):
            name, separator, value = declaration.partition(":")
            if not separator:
                continue
            name, value = name.strip().lower(), value.strip()
            important = value.lower().endswith("!important")
            if important:
                value = value[: -len("!important")].rstrip()
            order += 1
            priority = (important, specificity, order)
            expanded = {name: value}
            if name == "padding":
                expanded["padding-top"] = value.split()[0]
            if name == "background":
                expanded["background-color"] = value
            for property_name, property_value in expanded.items():
                if property_name not in priorities or priority >= priorities[property_name]:
                    properties[property_name] = property_value
                    priorities[property_name] = priority
    return properties


def _assert_html_banner_style(html_text: str) -> None:
    banner = _written_css_properties(
        html_text,
        {".retirement-banner", ".slide .retirement-banner", ".slide > .retirement-banner"},
    )
    assert not ({"visibility", "display", "opacity"} & banner.keys()), (
        "Retirement banner CSS must not set visibility, display or opacity"
    )
    foreground = to_rgba(banner["color"])
    background = to_rgba(banner["background-color"])
    assert foreground[3] > 0, "HTML retirement text has no opaque paint"
    assert background[3] > 0, "HTML retirement band has no opaque background"
    assert foreground[:3] != background[:3], "HTML retirement paint matches its background"


def test_html_banner_style_cannot_hide_warning() -> None:
    """Check producer CSS directly, without rendering the full 23-file historical export."""
    from mtp.cellulose_decr import tufte_plates

    _assert_html_banner_style(tufte_plates.render_slide_deck_html())


def _css_pixels(value: str) -> float:
    assert re.fullmatch(r"[0-9.]+px", value), f"Expected explicit pixel dimensions, got {value!r}"
    return float(value.removesuffix("px"))


def test_written_html_banners_keep_original_slide_content_clear(generated_tufte_dir: Path) -> None:
    html_text = (generated_tufte_dir / "cellulose_decr_slide_deck.html").read_text()
    _assert_html_banner_style(html_text)
    sections = re.findall(r"<section\b.*?</section>", html_text, re.DOTALL)
    assert len(sections) == 16
    for number, section in enumerate(sections, start=1):
        banners = re.findall(r'<div class="retirement-banner">(.*?)</div>', section)
        references = [
            route
            for route in re.findall(r'<img\b[^>]*\bsrc="([^"]+)"', section)
            if route.startswith("figures_tufte/")
        ]
        expected_references = (
            [f"figures_tufte/{HISTORICAL_PLATE_STEMS[number - 2]}.svg"] if 2 <= number <= 7 else []
        )
        assert references == expected_references, (
            f"Written slide {number}: historical image changed"
        )
        assert banners == ([] if expected_references else [RETIREMENT_BANNER]), number
        embedded_svgs = [
            ET.fromstring((generated_tufte_dir / route).read_text(encoding="utf-8"))
            for route in references
        ]
        embedded_svgs.extend(
            ET.fromstring(raw_svg) for raw_svg in re.findall(r"<svg\b.*?</svg>", section, re.DOTALL)
        )
        embedded_banners = sum(
            " ".join("".join(node.itertext()).split()) == RETIREMENT_BANNER
            for svg in embedded_svgs
            for node in svg.findall(".//{*}text")
        )
        assert len(banners) + embedded_banners == 1, f"Written slide {number}: warning count"
    banner = _written_css_properties(html_text, {".slide > .retirement-banner"})
    slide = _written_css_properties(html_text, {".slide"})
    cover = _written_css_properties(html_text, {".slide", ".cover-slide"})
    page = _written_css_properties(html_text, {"@page"})
    banner_height = _css_pixels(banner["height"])
    assert banner_height == 36
    assert _css_pixels(slide["padding-top"]) >= banner_height + 20
    assert _css_pixels(cover["padding-top"]) >= banner_height + 64
    assert _css_pixels(slide["height"]) >= 1000 + banner_height
    assert _css_pixels(page["size"].split()[1]) >= _css_pixels(slide["height"])


@pytest.mark.parametrize(
    ("slide", "renderer_name"), [(15, "render_fig8_svg"), (16, "render_fig9_svg")]
)
def test_written_inline_svgs_preserve_undecorated_content_and_original_size(
    generated_tufte_dir: Path, slide: int, renderer_name: str
) -> None:
    from mtp.cellulose_decr import tufte_plates

    html = (generated_tufte_dir / "cellulose_decr_slide_deck.html").read_text(encoding="utf-8")
    sections = re.findall(r"<section\b.*?</section>", html, re.DOTALL)
    inline_svgs = re.findall(r"<svg\b.*?</svg>", sections[slide - 1], re.DOTALL)
    assert len(inline_svgs) == 1
    inline = ET.fromstring(inline_svgs[0])
    assert list(map(float, inline.attrib["viewBox"].split())) == [0, 0, 1500, 740]
    assert float(inline.attrib["width"]) == 1500
    assert float(inline.attrib["height"]) == 740
    assert RETIREMENT_BANNER not in "".join(inline.itertext())
    renderer = getattr(tufte_plates, renderer_name)
    undecorated = cast(Callable[[], str], inspect.unwrap(renderer))
    assert inline_svgs[0].strip() == undecorated().strip()


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

    from mtp.cellulose_decr import tufte_plates

    fig8_path = generated_tufte_dir / "figures_tufte" / "fig8_literature_support_four_pairs.svg"
    fig8_svg = fig8_path.read_text(encoding="utf-8")
    deck_html = (generated_tufte_dir / "cellulose_decr_slide_deck.html").read_text(encoding="utf-8")
    undecorated_fig8 = cast(Callable[[], str], inspect.unwrap(tufte_plates.render_fig8_svg))
    assert undecorated_fig8().strip() in deck_html

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
    assert (
        tufte_plates.main(
            ["--audit-retired", "--analysis-dir", str(second_dir), "--code-commit", "deadbeef"]
        )
        == 0
    )

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
        tufte_plates.parse_args(["--audit-retired", "--analysis", "/tmp/x"])
    assert stopped.value.code == 2
    with pytest.raises(SystemExit) as stopped:
        tufte_plates.parse_args(["--audit-retired", "--code", "abc"])
    assert stopped.value.code == 2
    with pytest.raises(SystemExit) as stopped:
        tufte_plates.parse_args(["--audit"])
    assert stopped.value.code == 2
