"""Historical cellulose output is explicit audit-only and visibly retired (#1572)."""

from __future__ import annotations

import inspect
import os
import re
import subprocess
import xml.etree.ElementTree as ET
from collections.abc import Callable
from pathlib import Path
from typing import cast

import matplotlib.pyplot as plt
import pytest
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.text import Text

from mtp.cellulose_decr import tufte_plates

ROOT = Path(__file__).resolve().parents[1]
BANNER = "RETIRED: historical, unverified numbers (#1572)"
BUILDERS = (
    tufte_plates.build_plate1_convergence_slope_and_overlap,
    tufte_plates.build_plate2_core_guild_dot_table,
    tufte_plates.build_plate3_why_only_stool,
    tufte_plates.build_plate4_trait_landscape_and_four_pairs,
    tufte_plates.build_plate5_three_statistical_arms,
    tufte_plates.build_plate6_combining_results_integrated_evidence,
)
PLATE_SLIDE_ROUTES = {
    2: "figures_tufte/fig1_tufte_convergence_slope_and_overlap.svg",
    3: "figures_tufte/fig2_tufte_core_guild_dot_table.svg",
    4: "figures_tufte/fig3_tufte_why_only_stool.svg",
    5: "figures_tufte/fig4_tufte_trait_landscape_and_four_pairs.svg",
    6: "figures_tufte/fig5_tufte_three_statistical_arms.svg",
    7: "figures_tufte/fig6_tufte_combining_results_integrated_evidence.svg",
}


@pytest.fixture(scope="module", autouse=True)
def retirement_import_origin() -> None:
    assert Path(tufte_plates.__file__).resolve() == (
        ROOT / "src/mtp/cellulose_decr/tufte_plates.py"
    )


def test_python_writer_refuses_without_opt_in_before_directory_creation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "absent" / "historical"
    original = Path.mkdir

    def guarded_mkdir(
        path: Path,
        mode: int = 0o777,
        parents: bool = False,
        exist_ok: bool = False,
    ) -> None:
        if path.is_relative_to(output):
            pytest.fail("retired writer attempted output before checking opt-in")
        original(path, mode=mode, parents=parents, exist_ok=exist_ok)

    monkeypatch.setattr(Path, "mkdir", guarded_mkdir)
    with pytest.raises(ValueError, match="--audit-retired"):
        tufte_plates.generate_tufte_plates_and_deck(output)
    assert not output.exists()


def test_cli_refuses_without_opt_in_before_calling_writer(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def forbidden_writer(*_args: object, **_kwargs: object) -> list[Path]:
        pytest.fail("CLI invoked retired writer without explicit opt-in")

    monkeypatch.setattr(tufte_plates, "generate_tufte_plates_and_deck", forbidden_writer)
    with pytest.raises(SystemExit) as stopped:
        tufte_plates.main(["--analysis-dir", str(tmp_path / "absent")])
    assert stopped.value.code == 2
    assert not (tmp_path / "absent").exists()


def test_cli_forwards_explicit_opt_in(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[tuple[Path, int, bool]] = []

    def writer(analysis_dir: Path, *, dpi: int, audit_retired: bool = False) -> list[Path]:
        calls.append((analysis_dir, dpi, audit_retired))
        return []

    monkeypatch.setattr(tufte_plates, "generate_tufte_plates_and_deck", writer)
    assert (
        tufte_plates.main(["--audit-retired", "--analysis-dir", str(tmp_path), "--dpi", "72"]) == 0
    )
    assert calls == [(tmp_path, 72, True)]


def test_help_identifies_retirement_and_pending_replacement(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as stopped:
        tufte_plates.parse_args(["--help"])
    assert stopped.value.code == 0
    text = capsys.readouterr().out
    assert "--audit-retired" in text
    assert "#1572" in text
    assert "#1670" in text
    assert "unverified" in text


@pytest.mark.parametrize("builder", BUILDERS, ids=lambda builder: builder.__name__)
def test_every_historical_canvas_has_an_unclipped_separate_banner(
    builder: Callable[[], Figure],
) -> None:
    figure = builder()
    try:
        banners = [text for text in figure.texts if text.get_text() == BANNER]
        assert len(banners) == 1
        figure.canvas.draw()
        bounds = banners[0].get_window_extent()
        assert figure.bbox.contains(bounds.x0, bounds.y0)
        assert figure.bbox.contains(bounds.x1, bounds.y1)
        for text in figure.findobj(match=Text):
            if text is banners[0] or not text.get_visible() or not text.get_text().strip():
                continue
            assert not bounds.overlaps(text.get_window_extent()), text.get_text()
    finally:
        plt.close(figure)


@pytest.mark.parametrize(
    ("renderer", "width", "height"),
    [
        (tufte_plates.render_fig7_svg, 1400, 770),
        (tufte_plates.render_fig8_svg, 1500, 740),
        (tufte_plates.render_fig9_svg, 1500, 740),
    ],
)
def test_standalone_svg_banner_extends_canvas_above_original_coordinates(
    renderer: Callable[[], str],
    width: int,
    height: int,
) -> None:
    svg = ET.fromstring(renderer())
    banners = [node for node in svg.iter() if node.text == BANNER]
    assert len(banners) == 1
    foreground = to_rgba(banners[0].attrib["fill"])
    backgrounds = [node for node in svg.findall("{*}rect") if float(node.attrib.get("y", "0")) < 0]
    assert len(backgrounds) == 1
    background = to_rgba(backgrounds[0].attrib["fill"])
    assert foreground[3] > 0, "SVG retirement text has no opaque paint"
    assert background[3] > 0, "SVG retirement band has no opaque background"
    assert foreground[:3] != background[:3], "SVG retirement paint matches its background"
    assert float(banners[0].attrib.get("opacity", "1")) > 0
    assert float(banners[0].attrib.get("fill-opacity", "1")) > 0
    x, y, output_width, output_height = map(float, svg.attrib["viewBox"].split())
    assert x == 0
    assert output_width == width
    assert y < 0
    assert y + output_height == height
    assert y < float(banners[0].attrib["y"]) < 0
    assert float(svg.attrib["height"].removesuffix("pt").removesuffix("px")) == output_height


def test_each_of_sixteen_html_slides_has_its_own_retirement_banner() -> None:
    html = tufte_plates.render_slide_deck_html()
    sections = re.findall(r"<section\b.*?</section>", html, re.DOTALL)
    assert len(sections) == 16
    for number, section in enumerate(sections, start=1):
        banners = re.findall(r'<div class="retirement-banner">(.*?)</div>', section)
        references = [
            route
            for route in re.findall(r'<img\b[^>]*\bsrc="([^"]+)"', section)
            if route.startswith("figures_tufte/")
        ]
        expected_references = [PLATE_SLIDE_ROUTES[number]] if number in PLATE_SLIDE_ROUTES else []
        assert references == expected_references, f"Slide {number}: unexpected historical image"
        assert banners == ([] if expected_references else [BANNER]), number
        inline_banners = sum(
            " ".join("".join(node.itertext()).split()) == BANNER
            for raw_svg in re.findall(r"<svg\b.*?</svg>", section, re.DOTALL)
            for node in ET.fromstring(raw_svg).findall(".//{*}text")
        )
        # The six referenced plates have one banner each, independently checked by
        # the canvas tests and by reading their saved SVGs in the export test.
        assert len(banners) + inline_banners + len(expected_references) == 1, number


@pytest.mark.parametrize(
    ("slide", "renderer"), [(15, tufte_plates.render_fig8_svg), (16, tufte_plates.render_fig9_svg)]
)
def test_inline_historical_svgs_keep_original_size_without_a_second_banner(
    slide: int, renderer: Callable[[], str]
) -> None:
    sections = re.findall(
        r"<section\b.*?</section>", tufte_plates.render_slide_deck_html(), re.DOTALL
    )
    inline_svgs = re.findall(r"<svg\b.*?</svg>", sections[slide - 1], re.DOTALL)
    assert len(inline_svgs) == 1
    inline = ET.fromstring(inline_svgs[0])
    assert list(map(float, inline.attrib["viewBox"].split())) == [0, 0, 1500, 740]
    assert float(inline.attrib["width"]) == 1500
    assert float(inline.attrib["height"]) == 740
    assert BANNER not in "".join(inline.itertext())
    undecorated = cast(Callable[[], str], inspect.unwrap(renderer))
    assert inline_svgs[0].strip() == undecorated().strip()


def test_default_figure_recipe_only_invokes_current_renderer(tmp_path: Path) -> None:
    """Run the actual recipe with an inert uv recorder, never an analysis command."""
    executable = tmp_path / "uv"
    log = tmp_path / "calls.txt"
    executable.write_text('#!/bin/sh\nprintf "%s\\n" "$*" >> "$RECORDED_UV_CALLS"\n')
    executable.chmod(0o700)
    result = subprocess.run(
        ["just", "--justfile", str(ROOT / "Justfile"), "cellulose-decr-figures"],
        cwd=ROOT,
        env={
            **os.environ,
            "PATH": f"{tmp_path}{os.pathsep}{os.environ['PATH']}",
            "RECORDED_UV_CALLS": str(log),
        },
        check=False,
        capture_output=True,
        text=True,
        timeout=15,
    )
    assert result.returncode == 0, result.stderr
    calls = log.read_text().splitlines()
    assert len(calls) == 2
    assert calls[0] == "run python -m mtp.cellulose_decr.retired_outputs warn"
    assert "-m mtp.cellulose_decr.figures " in calls[1]


def test_historical_recipe_does_not_supply_audit_permission_implicitly() -> None:
    result = subprocess.run(
        ["just", "--dry-run", "--justfile", str(ROOT / "Justfile"), "cellulose-decr-tufte"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=10,
    )
    command = result.stdout + result.stderr
    assert "-m mtp.cellulose_decr.tufte_plates" in command
    assert "--audit-retired" not in command
