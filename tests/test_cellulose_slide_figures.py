"""Opt-in, lossless slide pages from saved cellulose tables (#1885).

Distinct fixture rows exercise actual artists and pagination, one figure at a
time. The seven report canvases remain a separate, unchanged default export.
"""

from __future__ import annotations

import hashlib
import inspect
import json
from collections.abc import Callable, Iterator
from dataclasses import replace
from pathlib import Path
from typing import Any, cast
from urllib.parse import unquote

import matplotlib
import numpy as np
import pandas as pd
import pytest
from matplotlib.backend_bases import RendererBase
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.text import Text
from numpy.typing import NDArray
from scipy.stats import binomtest

from mtp.cellulose_decr import figures
from mtp.cellulose_decr.analysis import ALL_CONTRASTS, ANALYSIS_DIR
from test_cellulose_nominal_labels import NULL_COLUMNS, NULL_KEYS, synthetic_results


def slide_results() -> figures.Results:
    """Forty-nine region rows cannot fit legibly on one slide."""
    results = synthetic_results()
    catalog = pd.concat([results.catalog for _ in range(7)], ignore_index=True)
    catalog["taxon_name"] = [f"Taxon number {index:02d}" for index in range(len(catalog))]
    members = [f"in_{short}_feces_dec" for short in figures.SHORT]
    catalog["n_conditions_feces_dec"] = catalog[members].sum(axis=1)
    catalog["is_3way_core"] = catalog["n_conditions_feces_dec"].eq(3)
    catalog["has_direct_trait_edge"] = catalog.index % 3 == 0
    catalog["taxonomic_rank"] = "species"
    catalog["is_strictly_net_decreased_all3"] = catalog.index % 2 == 0
    for offset, short in enumerate(figures.SHORT):
        catalog[f"{short}_dec_studies"] = catalog.index + offset + 2
        catalog[f"{short}_inc_studies"] = catalog.index + 2 * offset + 1
    benchmark = pd.DataFrame(
        {
            "trait_id": f"fixture:{index}",
            "trait_label": f"Trait number {index:02d}",
            "object_name": f"Trait number {index:02d}",
            "trait_predicate": "capable of",
            "is_growth_medium": index == 2,
            "mean_pairwise_jaccard_distance": (index + 1) / 25,
        }
        for index in range(18)
    )
    results.summary["trait_id"] = "fixture:0"
    results.summary["observed"]["three_way_core"] = 7
    results.summary["observed"]["three_way_union"] = len(catalog)
    results.pairwise.loc[
        results.pairwise["endpoint"].eq("mean_pairwise_jaccard_distance"), "comparison"
    ] = "Mean across all 3 pairs"
    results.sensitivity.loc[
        results.sensitivity["comparison"].str.startswith("Mean"), "comparison"
    ] = "Mean across all 3 pairs (3-way core)"
    for index, key in enumerate(NULL_KEYS):
        results.summary["null_models"][key]["three_way_core"].update(
            expected=index + 1.2, p2_5=index + 0.2, p97_5=index + 2.2
        )
    for row_index in results.pairwise.index:
        for column_index, column in enumerate(NULL_COLUMNS):
            results.pairwise.loc[row_index, f"{column}_expected_jaccard_distance"] = (
                0.1 + row_index / 10 + column_index / 100
            )
    for index in results.compartments.index:
        results.compartments.loc[index, "paired_exp_down_dec_gt_inc"] = index + 1
        results.compartments.loc[index, "paired_exp_up_inc_gt_dec"] = index + 3
        results.compartments.loc[index, "study_down"] = index + 2
        results.compartments.loc[index, "study_up"] = index + 5
    return replace(results, catalog=catalog, benchmark=benchmark)


def _pages(results: figures.Results) -> Iterator[Any]:
    builder = getattr(figures, "slide_pages", None)
    assert callable(builder), "#1885 requires additive figures.slide_pages(results)"
    pages = builder(results)
    assert inspect.isgenerator(pages), "pages must be lazy so figures do not accumulate"
    return pages


def _texts(figure: Figure) -> dict[str, str]:
    return {
        unquote(str(text.get_gid())): text.get_text()
        for text in figure.findobj(Text)
        if text.get_gid() is not None
    }


def _lines(figure: Figure) -> dict[str, Line2D]:
    return {unquote(str(line.get_gid())): line for line in figure.findobj(Line2D) if line.get_gid()}


def _xdata(line: Line2D) -> NDArray[np.float64]:
    """The tested artists carry numeric coordinates, not categorical plot inputs."""
    return np.asarray(line.get_xdata(), dtype=float)


def test_slide_pages_preserve_null_values_and_every_taxon_region_row() -> None:
    results = slide_results()
    expected: dict[str, str] = {}
    for row in results.pairwise.to_dict("records"):
        for column in NULL_COLUMNS:
            identity = f"null_expected:{row['comparison']}:{column}"
            expected[identity] = f"{row[f'{column}_expected_jaccard_distance']:.3f}"
    names: list[str] = []
    values: dict[str, str] = {}
    core_keys: set[str] = set()
    region_pages = 0
    for page in _pages(results):
        try:
            if page.section == "null_table":
                values.update(
                    {k: v for k, v in _texts(page.figure).items() if k.startswith("null_expected:")}
                )
            elif page.section == "taxon_regions":
                region_pages += 1
                names.extend(
                    v for k, v in _texts(page.figure).items() if k.startswith("taxon_name:")
                )
            elif page.section == "null_core":
                lines = _lines(page.figure)
                for key in NULL_KEYS:
                    row = results.summary["null_models"][key]["three_way_core"]
                    np.testing.assert_allclose(_xdata(lines[f"null_mean:{key}"]), [row["expected"]])
                    np.testing.assert_allclose(
                        _xdata(lines[f"null_interval:{key}"]), [row["p2_5"], row["p97_5"]]
                    )
                    core_keys.add(key)
                np.testing.assert_allclose(_xdata(lines["null_observed"]), [7, 7])
        finally:
            page.figure.clear()
    assert values == expected
    assert sorted(names) == sorted(results.catalog["taxon_name"])
    assert region_pages > 1
    assert core_keys == set(NULL_KEYS)


def test_slide_pages_preserve_share_intervals_sensitivity_and_interactions() -> None:
    results = slide_results()
    expected_shares: dict[str, tuple[float, float, float]] = {}
    for row in results.compartments.to_dict("records"):
        for unit, down_key, up_key in (
            ("experiment", "paired_exp_down_dec_gt_inc", "paired_exp_up_inc_gt_dec"),
            ("study", "study_down", "study_up"),
        ):
            down, up = int(row[down_key]), int(row[up_key])
            interval = binomtest(down, down + up).proportion_ci()
            key = (
                f"{row['contrast_population']}:{row['compartment']}:{row['short_condition']}:{unit}"
            )
            expected_shares[key] = down / (down + up), interval.low, interval.high
    found_shares: dict[str, tuple[float, float, float]] = {}
    sensitivity: dict[str, str] = {}
    interactions: dict[str, tuple[float, float, float]] = {}
    for page in _pages(results):
        try:
            lines = _lines(page.figure)
            for prefix, destination in (("share", found_shares), ("interaction", interactions)):
                for key, line in lines.items():
                    if key.startswith(f"{prefix}_mean:"):
                        identity = key.removeprefix(f"{prefix}_mean:")
                        bounds = _xdata(lines[f"{prefix}_interval:{identity}"])
                        destination[identity] = (
                            float(_xdata(line)[0]),
                            float(bounds[0]),
                            float(bounds[1]),
                        )
            if page.section == "sensitivity":
                sensitivity.update(
                    {
                        k: v
                        for k, v in _texts(page.figure).items()
                        if k.startswith("sensitivity_value:")
                    }
                )
        finally:
            page.figure.clear()
    assert found_shares.keys() == expected_shares.keys()
    for key, value in expected_shares.items():
        np.testing.assert_allclose(found_shares[key], value)
    expected_sensitivity: dict[str, str] = {}
    for row in results.sensitivity.to_dict("records"):
        identity = f"sensitivity_value:{row['definition_slug']}:{row['comparison']}"
        expected_sensitivity[identity] = f"{row['jaccard_distance']:.3f}"
    assert sensitivity == expected_sensitivity
    selected = results.interaction.loc[results.interaction["contrast_population"].eq(ALL_CONTRASTS)]
    assert interactions == {
        f"{row['short_condition']}:{row['compartment_b']}:{row['unit']}": (
            row["conditional_or"],
            row["conditional_or_ci95_low"],
            row["conditional_or_ci95_high"],
        )
        for row in selected.to_dict("records")
    }


def test_slide_pages_preserve_shared_taxon_cells_and_fifteen_benchmark_values() -> None:
    results = slide_results()
    names: list[str] = []
    counts: dict[str, str] = {}
    benchmark: dict[str, float] = {}
    shared_pages = 0
    for page in _pages(results):
        try:
            if page.section == "shared_taxa":
                shared_pages += 1
                text = _texts(page.figure)
                names.extend(v for k, v in text.items() if k.startswith("shared_name:"))
                counts.update({k: v for k, v in text.items() if k.startswith("shared_count:")})
            elif page.section == "trait_benchmark":
                benchmark.update(
                    {
                        k.removeprefix("benchmark_mean:"): float(_xdata(line)[0])
                        for k, line in _lines(page.figure).items()
                        if k.startswith("benchmark_mean:")
                    }
                )
        finally:
            page.figure.clear()
    selected = results.catalog.loc[results.catalog["n_conditions_feces_dec"].ge(2)]
    assert sorted(names) == sorted(selected["taxon_name"])
    assert shared_pages > 1
    expected_counts: dict[str, str] = {}
    for row in selected.to_dict("records"):
        for short in figures.SHORT:
            identity = f"shared_count:{row['taxon_name']}:{short}"
            expected_counts[identity] = (
                f"{row[f'{short}_dec_studies']}/{row[f'{short}_inc_studies']}"
            )
    assert counts == expected_counts
    selected_traits = results.benchmark.loc[~results.benchmark["is_growth_medium"]].head(15)
    assert benchmark == dict(
        zip(
            selected_traits["trait_id"],
            selected_traits["mean_pairwise_jaccard_distance"],
            strict=True,
        )
    )


@pytest.mark.parametrize("source", ["synthetic", "committed"])
def test_every_slide_page_meets_native_slide_size_font_and_canvas_bounds(source: str) -> None:
    results = slide_results() if source == "synthetic" else figures.load(ANALYSIS_DIR)
    sections: set[str] = set()
    stems: set[str] = set()
    sources: set[str] = set()
    with matplotlib.rc_context(figures.STYLE):
        for page in _pages(results):
            try:
                assert page.stem not in stems
                stems.add(page.stem)
                sections.add(page.section)
                sources.add(page.source_figure)
                if page.section == "interactions":
                    assert "Cellulose annotation direction" in page.figure.texts[0].get_text()
                np.testing.assert_allclose(page.figure.get_size_inches(), [13.33, 7.5])
                canvas = FigureCanvasAgg(page.figure)
                draw = cast(Callable[[], None], canvas.draw)
                draw()
                get_renderer = cast(Callable[[], RendererBase], canvas.get_renderer)
                renderer = get_renderer()
                bounds = page.figure.bbox
                # _new creates title, subtitle and footer before section artists.
                # Test the actual wrapped footer, including its complete multiline box.
                footer = page.figure.texts[2]
                footer_box = footer.get_window_extent(renderer)
                for text in page.figure.findobj(Text):
                    if text.get_visible() and text.get_text().strip():
                        assert text.get_fontsize() >= 12, (page.stem, text.get_text())
                        box = text.get_window_extent(renderer)
                        context = (page.stem, text.get_text(), box)
                        assert box.x0 >= bounds.x0 - 1, context
                        assert box.y0 >= bounds.y0 - 1, context
                        assert box.x1 <= bounds.x1 + 1, context
                        assert box.y1 <= bounds.y1 + 1, context
                        if text is not footer:
                            assert not footer_box.overlaps(box), (
                                page.stem,
                                "footer collision",
                                text.get_text(),
                                footer_box,
                                box,
                            )
                for line in page.figure.findobj(Line2D):
                    if line.get_visible() and line.get_linestyle() not in ("None", "none", "", " "):
                        assert line.get_linewidth() >= 0.5
                assert all(
                    getattr(page.caption, field)
                    for field in ("title", "takeaway", "caveat", "evidence")
                )
                vintage = results.summary["bugsigdb_vintage"]
                assert vintage in page.caption.caveat
                assert "Population:" in page.caption.caveat
                assert "PD, Parkinson's disease" in page.caption.caveat
                assert any(
                    f"BugSigDB {vintage}" in text.get_text() for text in page.figure.findobj(Text)
                )
            finally:
                page.figure.clear()
    assert sections == {
        "null_table",
        "null_core",
        "taxon_regions",
        "compartment_shares",
        "sensitivity",
        "interactions",
        "shared_taxa",
        "trait_benchmark",
    }
    assert sources == {
        "fig1_pairwise_and_3way_jaccard_nulls",
        "fig2_direction_by_body_site_and_definition",
        "fig3_core_taxa_and_trait_benchmark",
    }


@pytest.mark.parametrize("enabled", [False, True])
def test_slide_exports_are_opt_in_and_registered_with_captions_and_hashes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, enabled: bool
) -> None:
    builder = getattr(figures, "slide_pages", None)
    page_type = getattr(figures, "SlidePage", None)
    assert callable(builder), "#1885 slide export API is absent"
    assert callable(page_type), "#1885 slide export API is absent"
    caption = figures.Caption(
        "Fixture title", "Fixture takeaway", "Fixture caveat", figures.FIGURE_INPUTS[0]
    )
    seen: list[bool] = []

    def pages(_results: figures.Results) -> Iterator[Any]:
        seen.append(True)
        yield page_type(
            "fixture_page",
            "fixture_source",
            "fixture_section",
            Figure(figsize=(13.33, 7.5)),
            caption,
        )

    def small_figure(_results: figures.Results) -> Figure:
        return Figure(figsize=(1, 1))

    def fixed_caption(_results: figures.Results) -> figures.Caption:
        return caption

    def statement(_results: figures.Results, *, include_scope: bool = False) -> str:
        return "fixture. scope" if include_scope else "fixture"

    def save(figure: Figure, directory: Path, stem: str, book: PdfPages) -> list[Path]:
        directory.mkdir(parents=True, exist_ok=True)
        paths = [directory / f"{stem}{suffix}" for suffix in (".svg", ".pdf", ".png")]
        for path in paths:
            path.write_bytes(path.name.encode())
        save_page = cast(Callable[[Figure], None], book.savefig)
        save_page(figure)
        return paths

    for name, _function in inspect.getmembers(figures, inspect.isfunction):
        if name.startswith(("graphical_", "detailed_")):
            monkeypatch.setattr(figures, name, small_figure)
        elif name.startswith("caption_"):
            monkeypatch.setattr(figures, name, fixed_caption)
        elif name.endswith("_statement"):
            monkeypatch.setattr(figures, name, statement)
    monkeypatch.setattr(figures, "load", lambda _path: slide_results())
    monkeypatch.setattr(figures, "slide_pages", pages)
    monkeypatch.setattr(figures, "save", save)
    for name in figures.FIGURE_INPUTS:
        (tmp_path / name).write_text("fixture\n")
    argv = ["--analysis-dir", str(tmp_path)] + (["--slide-pages"] if enabled else [])
    assert figures.main(argv) == 0
    manifest = json.loads((tmp_path / "figures/figures_provenance.json").read_text())
    captions = json.loads((tmp_path / "figures/figure_captions.json").read_text())
    assert len(manifest["outputs"]) == 24 + (4 if enabled else 0)
    assert seen == ([True] if enabled else [])
    assert (tmp_path / "cellulose_decr_slide_figures.pdf").exists() is enabled
    assert ("figures_slides/fixture_page.svg" in captions["figures"]) is enabled
    if enabled:
        assert (
            captions["figures"]["figures_slides/fixture_page.svg"]["evidence"] == caption.evidence
        )
        assert manifest["slide_pages"] == [
            {
                "stem": "fixture_page",
                "source_figure": "fixture_source",
                "section": "fixture_section",
            }
        ]
        producer_files = manifest["producer_files"]
        assert "src/mtp/cellulose_decr/slide_figures.py" in producer_files
        root = Path(figures.__file__).resolve().parents[3]
        for relative, digest in producer_files.items():
            assert hashlib.sha256((root / relative).read_bytes()).hexdigest() == digest
    for relative, digest in manifest["outputs"].items():
        assert hashlib.sha256((tmp_path / relative).read_bytes()).hexdigest() == digest


@pytest.mark.parametrize(
    ("low", "estimate", "high", "expected_limits"),
    [
        (
            0.0001,
            0.001,
            0.005,
            {"low": (0.01, "<"), "high": (0.01, "<"), "estimate": (0.01, "<")},
        ),
        (
            2000.0,
            4000.0,
            8000.0,
            {"low": (1000.0, ">"), "high": (1000.0, ">"), "estimate": (1000.0, ">")},
        ),
        (0.001, 1.0, 10000.0, {"low": (0.01, "<"), "high": (1000.0, ">")}),
        (0.5, 2.0, 8.0, {}),
        (0.01, 0.01, 1000.0, {}),
        (0.0, 0.0, 0.02, {"low": (0.01, "<"), "estimate": (0.01, "<")}),
        (
            0.5,
            float("inf"),
            float("inf"),
            {"high": (1000.0, ">"), "estimate": (1000.0, ">")},
        ),
    ],
    ids=[
        "wholly-below",
        "wholly-above",
        "spans-both",
        "ordinary",
        "boundaries",
        "zero",
        "infinite",
    ],
)
def test_interaction_slide_marks_axis_overflow_and_preserves_saved_values(
    low: float,
    estimate: float,
    high: float,
    expected_limits: dict[str, tuple[float, str]],
) -> None:
    """CI ends use their axis side; an infinite bound cannot hide the visible CI (S1)."""
    results = slide_results()
    selected = (
        results.interaction.loc[results.interaction["contrast_population"].eq(ALL_CONTRASTS)]
        .head(1)
        .copy()
    )
    selected["conditional_or"] = estimate
    selected["conditional_or_ci95_low"] = low
    selected["conditional_or_ci95_high"] = high
    if estimate == 0:
        selected["a_down"] = 0
    elif estimate == float("inf"):
        selected["a_up"] = 0
    row = selected.iloc[0]
    identity = f"{row['short_condition']}:{row['compartment_b']}:{row['unit']}"
    seen = False
    for page in _pages(replace(results, interaction=selected)):
        try:
            if page.section != "interactions":
                continue
            seen = True
            lines = _lines(page.figure)
            # Both ends clamp independently; wholly off-axis CIs collapse at the
            # appropriate boundary, while a spanning infinite CI has a real segment.
            clipped = [min(max(value, 0.01), 1000.0) for value in (low, high)]
            np.testing.assert_allclose(_xdata(lines[f"interaction_interval:{identity}"]), clipped)
            assert lines[f"interaction_interval:{identity}"].get_visible()
            if 0.01 <= estimate <= 1000.0:
                np.testing.assert_allclose(
                    _xdata(lines[f"interaction_mean:{identity}"]), [estimate]
                )
                assert lines[f"interaction_mean:{identity}"].get_marker() == "o"
            assert f"Saved OR {estimate}; exact 95% interval [{low}, {high}]" in page.caption.caveat
            markers = {
                part: lines[f"interaction_{part}_limit:{identity}"]
                for part in ("low", "high", "estimate")
                if f"interaction_{part}_limit:{identity}" in lines
            }
            assert markers.keys() == expected_limits.keys()
            for part, (boundary, direction) in expected_limits.items():
                np.testing.assert_allclose(_xdata(markers[part]), [boundary])
                assert markers[part].get_marker() == direction
                assert markers[part].get_visible()
                assert not markers[part].get_clip_on()
        finally:
            page.figure.clear()
    assert seen


@pytest.mark.parametrize(
    ("route", "kind"),
    [
        ("figures_slides", "populated_directory"),
        ("figures_slides", "empty_directory"),
        ("figures_slides", "file"),
        ("figures_slides", "dangling_symlink"),
        ("cellulose_decr_slide_figures.pdf", "file"),
        ("cellulose_decr_slide_figures.pdf", "empty_directory"),
        ("cellulose_decr_slide_figures.pdf", "dangling_symlink"),
    ],
)
def test_existing_slide_routes_are_refused_before_any_report_mutation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, route: str, kind: str
) -> None:
    """A shorter rerender cannot strand old paginated pages beside new manifests (S2)."""
    path = tmp_path / route
    if kind in ("populated_directory", "empty_directory"):
        path.mkdir()
        if kind == "populated_directory":
            (path / "old_page_99.svg").write_bytes(b"old slide page")
    elif kind == "dangling_symlink":
        path.symlink_to(tmp_path / "absent-slide-target")
    else:
        path.write_bytes(b"old combined slides")
    original_stat = path.lstat()
    reports = tmp_path / "figures"
    reports.mkdir()
    sentinels = {
        reports / "figure_captions.json": b"old captions",
        reports / "figures_provenance.json": b"old provenance",
        tmp_path / "cellulose_decr_figures.pdf": b"old report PDF",
        tmp_path / "cellulose_decr_graphical_figures.pdf": b"old graphical PDF",
    }
    for file, content in sentinels.items():
        file.write_bytes(content)

    def forbidden_draw(_results: figures.Results) -> Figure:
        raise AssertionError("existing slide route was not refused before report drawing")

    monkeypatch.setattr(figures, "load", lambda _path: slide_results())
    monkeypatch.setattr(figures, "graphical_convergence", forbidden_draw)
    try:
        with pytest.raises(FileExistsError, match=r"[Ss]lide"):
            figures.main(["--analysis-dir", str(tmp_path), "--slide-pages"])
    finally:
        for file, content in sentinels.items():
            assert file.is_file(), f"Prior report was removed: {file.name}"
            assert file.read_bytes() == content
        assert path.lstat().st_ino == original_stat.st_ino
        if kind == "populated_directory":
            assert (path / "old_page_99.svg").read_bytes() == b"old slide page"
        elif kind == "file":
            assert path.read_bytes() == b"old combined slides"
        elif kind == "dangling_symlink":
            assert path.is_symlink()
            assert path.readlink() == tmp_path / "absent-slide-target"
