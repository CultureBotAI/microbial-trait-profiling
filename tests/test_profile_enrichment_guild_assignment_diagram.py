"""Tests for the BugSigDB-to-Sankey workflow stage labels (#748)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import matplotlib.axes
import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"


def _load_workflow_diagram_module() -> ModuleType:
    script_path = SCRIPTS_DIR / "create_profile_enrichment_guild_assignment_diagram.py"
    spec = importlib.util.spec_from_file_location(
        "create_profile_enrichment_guild_assignment_diagram_test",
        script_path,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


workflow_diagram = _load_workflow_diagram_module()


def _sample_counts() -> dict[str, object]:
    return {
        "bugsigdb_version": "v1.3.1",
        "bugsigdb_doi": "10.5281/zenodo.19740674",
        "bugsigdb_snapshot": "2026-04-01",
        "profile_count": 5629,
        "study_count": 1047,
        "experiment_count": 3372,
        "condition_count": 529,
        "body_site_count": 172,
        "body_site_group_count": 9,
        "host_species_count": 1,
        "selected_trait_count": 60,
        "selected_trait_category_count": 4,
        "direction_grid_condition_count": 205,
        "direction_grid_row_count": 820,
        "facultative_universe_count": 1992,
        "facultative_condition_count": 205,
        "anaerobe_condition_count": 205,
        "anaerobe_universe_count": 1722,
        "facultative_target_taxa_count": 10,
        "anaerobe_target_taxa_count": 7,
        "full_dump_row_count": 5400,
        "full_dump_study_count": 1246,
        "full_dump_complete_row_count": 5200,
        "full_dump_homo_sapiens_row_count": 4800,
        "full_dump_homo_sapiens_study_count": 1100,
        "raw_profile_count": 5629,
        "raw_profile_study_count": 1047,
        "raw_profile_experiment_count": 3372,
        "raw_trait_column_count": 2331,
        "paired_profile_count": 2257,
        "paired_study_count": 932,
        "paired_condition_count": 437,
    }


def test_workflow_stage_two_describes_canonical_human_16s_wms_upstream_table() -> None:
    """Box 2 must describe the canonical Homo sapiens 16S/WMS upstream profile matrix (#748)."""
    counts = _sample_counts()
    stages = workflow_diagram.build_workflow_stages(counts)

    assert [stage["number"] for stage in stages] == ["1", "2", "3", "4", "5"]
    stage_one = stages[0]
    stage_two = stages[1]

    assert stage_one["title"] == "BugSigDB data"
    assert "Filters: Homo sapiens host; 16S or WMS sequencing." in stage_one["body"]

    assert stage_two["title"] == "Trait profiles"
    assert stage_two["body"] == (
        "5,629 canonical Homo sapiens 16S or WMS profiles x "
        "2,331 non-empty KG-Microbe trait columns before plotting trait filters."
    )
    assert "before human-host" not in stage_two["body"].lower()
    assert "human-host plotting" not in stage_two["body"].lower()


def test_create_workflow_diagram_renders_stage_two_without_stale_human_host_claim(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Rendered workflow diagram must draw the corrected Stage 2 copy inside its box (#748)."""
    counts = _sample_counts()
    monkeypatch.setattr(workflow_diagram, "get_summary_counts", lambda: counts)
    monkeypatch.setattr(workflow_diagram, "INVESTIGATION_DIR", tmp_path)

    captured_stage_two_Texts: list[tuple[float, float, float, float, str]] = []
    original_add_stage_box = workflow_diagram.add_stage_box

    def _recording_add_stage_box(
        ax: matplotlib.axes.Axes,
        *,
        xy: tuple[float, float],
        width: float,
        height: float,
        number: str,
        title: str,
        body: str,
        facecolor: str = "white",
        edgecolor: str = "#e2e8f0",
        accent: str = "#334155",
        body_width: int = 25,
    ) -> None:
        original_add_stage_box(
            ax,
            xy=xy,
            width=width,
            height=height,
            number=number,
            title=title,
            body=body,
            facecolor=facecolor,
            edgecolor=edgecolor,
            accent=accent,
            body_width=body_width,
        )
        if number == "2":
            body_artist = ax.texts[-1]
            fig = ax.figure
            fig.canvas.draw()
            bbox = ax.transData.inverted().transform_bbox(body_artist.get_window_extent())
            captured_stage_two_Texts.append(
                (bbox.x0, bbox.y0, bbox.x1, bbox.y1, body_artist.get_text())
            )

    monkeypatch.setattr(workflow_diagram, "add_stage_box", _recording_add_stage_box)

    out_png = tmp_path / "summary.png"
    out_svg = tmp_path / "summary.svg"
    workflow_diagram.create_workflow_diagram(output_png=out_png, output_svg=out_svg)

    assert out_png.is_file()
    assert out_png.stat().st_size > 0
    assert out_svg.is_file()
    assert out_svg.stat().st_size > 0

    assert len(captured_stage_two_Texts) == 1
    x0, y0, x1, y1, rendered_text = captured_stage_two_Texts[0]
    assert "5,629 canonical Homo" in rendered_text
    assert "sapiens 16S or WMS" in rendered_text
    assert "2,331 non-empty" in rendered_text
    assert "before plotting trait" in rendered_text
    assert "before human-host" not in rendered_text.lower()
    assert 22.0 <= x0 < x1 <= 22.0 + 17.4
    assert 43.0 <= y0 < y1 <= 43.0 + 12.0

    svg_text = out_svg.read_text(encoding="utf-8")
    assert "before human-host" not in svg_text.lower()
    assert "canonical Homo" in svg_text
