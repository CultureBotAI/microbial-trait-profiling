"""Keep follow-up report provenance and variant documentation consistent (#525)."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path("scripts/analyze_enrichment_infographic_followup.py")
REPORTS = Path("trait-visualization/trait-profile-investigation/enrichment-infographic-followup")
NAMES = (
    "README.md",
    "my_top3_README.md",
    "my_top3_conditions_relevant_traits_README.md",
    "my-top3-enriched-traits/README.md",
)


def load_generator(root: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location("followup_readme", root / SCRIPT)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_reports(root: Path, tmp_path: Path) -> None:
    module = load_generator(root)
    empty = pd.DataFrame()
    output = tmp_path / "README.md"
    module.write_report(
        output,
        pd.DataFrame(columns=["condition", "target_condition_rank_score"]),
        empty,
        empty,
        empty,
        empty,
        empty,
        empty,
        empty,
        "test traits",
        "test conditions",
        "test score",
        "",
        {},
        2,
    )
    generated = output.read_text()
    header = generated.split("The target conditions for this run", 1)[0]
    assert "Historical source figure: `data/derived/enrichment_infographic.png`" in header
    assert "CultureBotAI/microbial-trait-profiling" in header
    assert "not a claim that the historical figure is present" in header
    assert "## Report variants" in header
    for name in NAMES:
        assert f"`{name}`" in header
    assert "prefixed artifacts" in header
    assert "separate output directory" in header
    assert "not interchangeable" in header
    for name in NAMES:
        text = (root / REPORTS / name).read_text()
        assert text.startswith(header), name
        assert "\\" not in text, name
        assert "CultureBotAI/bsdbng/blob/" not in text, name
    assert "\\" not in generated
    source = (root / SCRIPT).read_text()
    assert "CultureBotAI/bsdbng/blob/" not in source
    # Native Windows Paths must not leak platform separators into Markdown.
    assert ".relative_to(REPO_ROOT)}`" not in source


def test_report_provenance_and_variants(tmp_path: Path) -> None:
    check_reports(ROOT, tmp_path)
