"""Keep table1-native-results/README.md synchronized with the native BH family wording (#1635)."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import cast

import pytest

from mtp.contribution_inputs import Base
from mtp.native_taxon_influence import native_summary, transformed

REPO_ROOT = Path(__file__).resolve().parents[1]
README_PATH = REPO_ROOT / "methods" / "table1-native-results" / "README.md"

EXPECTED_BH_SENTENCES = (
    "Native p-values test the mean; native BH families contain finite p-values among "
    "selected supported variable traits, not just Table 1. Selected untestable traits "
    "receive q = 1.0 outside the correction vector."
)
EXPECTED_COMPOSITION_PARAGRAPH = (
    "Native effects are experiment-weighted increased-minus-decreased z-score contrasts. "
    "Native p-values test the mean; native BH families contain finite p-values among "
    "selected supported variable traits, not just Table 1. Selected untestable traits "
    "receive q = 1.0 outside the correction vector. "
    "Zero-total profiles remain as zero fractions. "
    "Study-weighted comparisons, sign tests and conditional intervals are separately "
    "labeled in native_composition_comparison.tsv.gz. Refitted scales and fixed baseline "
    "SD comparisons must not be conflated."
)
EXPECTED_POST_1334_COMPOSITION_PARAGRAPH = (
    "Native effects are experiment-weighted increased-minus-decreased z-score contrasts. "
    "Native p-values test the mean; native BH families contain all selected supported "
    "variable traits (family_selected_traits), not just Table 1, padding selected "
    "untestable traits at p = 1.0 (q = 1.0) inside the correction vector. "
    "Zero-total profiles remain as zero fractions. "
    "Study-weighted comparisons, sign tests and conditional intervals are separately "
    "labeled in native_composition_comparison.tsv.gz. Refitted scales and fixed baseline "
    "SD comparisons must not be conflated."
)


def extract_composition_paragraph(markdown_text: str) -> str:
    """Return the first prose paragraph under the ## Composition heading."""
    marker = "## Composition\n\n"
    assert marker in markdown_text, "Missing ## Composition heading in README.md"
    after_heading = markdown_text.split(marker, 1)[1]
    paragraph = after_heading.split("\n\n", 1)[0]
    assert paragraph, "Empty ## Composition paragraph in README.md"
    return paragraph


def generated_published_composition_paragraph() -> str:
    """Render native_summary with a stub and apply the review-bundle filename substitution."""
    stub_base = cast(Base, SimpleNamespace(traits={}))
    raw_summary = native_summary(stub_base, [], [])
    published_summary = raw_summary.replace(
        "native_composition_comparison.tsv",
        "native_composition_comparison.tsv.gz",
    )
    return extract_composition_paragraph(published_summary)


def expected_generator_composition_paragraph() -> str:
    """Determine the required generator paragraph from transformed()'s live BH vector behavior."""
    keys = [
        ("s1", "e1", "inc"),
        ("s1", "e1", "dec"),
        ("s2", "e2", "inc"),
        ("s2", "e2", "dec"),
    ]
    profiles = {
        ("s1", "e1", "inc"): {"T_tested": 3, "T_untestable": 2, "T_filler": 5},
        ("s1", "e1", "dec"): {"T_tested": 1, "T_untestable": 1, "T_filler": 8},
        ("s2", "e2", "inc"): {"T_tested": 4, "T_untestable": 2, "T_filler": 4},
        ("s2", "e2", "dec"): {"T_tested": 1, "T_untestable": 1, "T_filler": 8},
    }
    policy = {
        "minimum_profile_support": 1,
        "maximum_profile_support_fraction": 1.0,
        "transform_multiplier": 100.0,
        "minimum_pairs": 2,
    }
    stub_base = cast(Base, SimpleNamespace(columns=["T_tested", "T_untestable"], traits={}))
    records, _, _ = transformed(stub_base, profiles, keys, policy)
    tested_row = records["T_tested"]
    untestable_row = records["T_untestable"]
    assert tested_row["status"] == "tested"
    assert untestable_row["selected"] is True
    assert untestable_row["status"] == "untestable"
    assert untestable_row["q"] == pytest.approx(1.0)
    tested_p = float(tested_row["p"])
    tested_q = float(tested_row["q"])
    if tested_q == pytest.approx(tested_p):
        return EXPECTED_COMPOSITION_PARAGRAPH
    if tested_q == pytest.approx(min(1.0, 2.0 * tested_p)):
        return EXPECTED_POST_1334_COMPOSITION_PARAGRAPH
    raise AssertionError(
        f"Unexpected BH adjustment ratio q={tested_q} vs p={tested_p} in transformed()"
    )


def check_readme_text(readme_text: str) -> None:
    """Validate the committed README text against the generator and #1635 BH family contract."""
    assert "native BH families contain all supported variable traits" not in readme_text
    assert EXPECTED_BH_SENTENCES in readme_text
    composition_paragraph = extract_composition_paragraph(readme_text)
    assert composition_paragraph == EXPECTED_COMPOSITION_PARAGRAPH
    assert generated_published_composition_paragraph() == expected_generator_composition_paragraph()


def test_table1_native_results_readme_matches_generator_bh_family_wording() -> None:
    readme_text = README_PATH.read_text(encoding="utf-8")
    check_readme_text(readme_text)
