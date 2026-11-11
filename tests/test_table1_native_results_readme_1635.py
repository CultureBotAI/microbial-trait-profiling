"""Keep table1-native-results/README.md synchronized with the native BH family wording (#1635)."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import cast

from mtp.contribution_inputs import Base
from mtp.native_taxon_influence import native_summary

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


def check_readme_text(readme_text: str) -> None:
    """Validate the committed README text against the generator and #1635 BH family contract."""
    assert "native BH families contain all supported variable traits" not in readme_text
    assert EXPECTED_BH_SENTENCES in readme_text
    composition_paragraph = extract_composition_paragraph(readme_text)
    assert composition_paragraph == EXPECTED_COMPOSITION_PARAGRAPH
    assert generated_published_composition_paragraph() in {
        EXPECTED_COMPOSITION_PARAGRAPH,
        EXPECTED_POST_1334_COMPOSITION_PARAGRAPH,
    }


def test_table1_native_results_readme_matches_generator_bh_family_wording() -> None:
    readme_text = README_PATH.read_text(encoding="utf-8")
    check_readme_text(readme_text)
