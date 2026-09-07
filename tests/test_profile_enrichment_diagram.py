"""Test suite for the workflow infographic diagram script."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from bsdbng.paths import REPO_ROOT

SCRIPT_PATH = REPO_ROOT / "scripts" / "create_profile_enrichment_guild_assignment_diagram.py"


def test_script_exists() -> None:
    assert SCRIPT_PATH.exists(), f"Expected script at {SCRIPT_PATH}"


def test_generate_workflow_diagram(tmp_path: Path) -> None:
    png_out = tmp_path / "test_workflow.png"
    svg_out = tmp_path / "test_workflow.svg"

    cmd = [
        sys.executable,
        str(SCRIPT_PATH),
        "--output-png",
        str(png_out),
        "--output-svg",
        str(svg_out),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    assert result.returncode == 0
    assert png_out.exists()
    assert png_out.stat().st_size > 0
    assert svg_out.exists()
    assert svg_out.stat().st_size > 0

    svg_content = svg_out.read_text(encoding="utf-8")
    assert "BugSigDB data to trait-enrichment Sankey diagrams" in svg_content
    assert "Trait-Enrichment" in svg_content
    assert "Tripartite Networks" in svg_content
    assert "LinkML Data Architecture" in svg_content
    assert "Taxa-Trait Mapping" in svg_content
    assert "Hypergeometric Pipeline" in svg_content
