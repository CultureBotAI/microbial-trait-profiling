"""Tests for project-wide biological interpretation guide and governance docs."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
METHODS_DIR = REPO_ROOT / "methods"
DOCS_DIR = REPO_ROOT / "docs"

BIOLOGICAL_INTERPRETATION_MD = METHODS_DIR / "biological_interpretation.md"
BIOLOGICAL_INTERPRETATION_AUDIT_MD = METHODS_DIR / "biological_interpretation_audit.md"
TABLE1_TRAIT_SEMANTICS_AUDIT_MD = DOCS_DIR / "table1_trait_semantics_audit.md"


def test_biological_interpretation_files_exist() -> None:
    """Verify that all biological interpretation and audit documentation files exist."""
    assert BIOLOGICAL_INTERPRETATION_MD.is_file(), "methods/biological_interpretation.md is missing"
    assert BIOLOGICAL_INTERPRETATION_AUDIT_MD.is_file(), (
        "methods/biological_interpretation_audit.md is missing"
    )
    assert TABLE1_TRAIT_SEMANTICS_AUDIT_MD.is_file(), (
        "docs/table1_trait_semantics_audit.md is missing"
    )


def test_biological_interpretation_structure() -> None:
    """Verify essential sections in methods/biological_interpretation.md."""
    content = BIOLOGICAL_INTERPRETATION_MD.read_text(encoding="utf-8")

    required_sections = [
        "## Evidence chain and measurement units",
        "## Four different meanings of direction",
        "## Wording and the evidence needed to go further",
        "## Minimum result and caption contract",
        "## Worked example: the Rosalind transcript and repository review",
        "## Review responsibilities and open evidence work",
    ]
    for section in required_sections:
        assert section in content, f"Missing section: {section}"

    # Verify resolution of Issue #218 (valid internal reference to
    # three_arm_statistics.md#shared-comparison-basis)
    assert "three_arm_statistics.md#shared-comparison-basis" in content
    assert "#175" in content or "PR #181" in content


def test_biological_interpretation_audit_findings() -> None:
    """Verify that audit findings (C01-C11, P01-P10, R01-R09, G01-G13) are documented."""
    content = BIOLOGICAL_INTERPRETATION_AUDIT_MD.read_text(encoding="utf-8")

    for i in range(1, 12):
        assert f"C{i:02d}" in content, f"Finding C{i:02d} missing from audit"
    for i in range(1, 11):
        assert f"P{i:02d}" in content, f"Finding P{i:02d} missing from audit"
    for i in range(1, 10):
        assert f"R{i:02d}" in content, f"Finding R{i:02d} missing from audit"
    for i in range(1, 14):
        assert f"G{i:02d}" in content, f"Finding G{i:02d} missing from audit"


def test_guild_methods_distinguish_annotation_zero_from_biological_absence() -> None:
    """C01 must describe recorded signal in both row selection and cell states."""
    content = " ".join((METHODS_DIR / "guild_condition_body_site_heatmap.md").read_text().split())
    rows = content.split("## Row selection and clustering", 1)[1].split("##", 1)[0]
    assert "no recorded guild annotation signal" in rows
    assert "no taxon carried" not in rows
    assert "no recorded guild annotation signal under the selected sources and mappings" in content
    assert "does not establish biological absence" in content
    assert "guild is simply absent" not in content


def test_trait_methods_distinguish_recorded_zero_from_biological_absence() -> None:
    """C03 applies to cell states and study ties without changing their policies."""
    content = " ".join(
        (METHODS_DIR / "table1_trait_condition_body_site_maps.md").read_text().split()
    )
    assert "no trait annotation signal was recorded" in content
    assert "zero recorded trait annotations on both sides" in content
    assert "This does not establish biological absence" in content
    assert "trait is absent from both sides" not in content
    assert "no taxon carried the trait" not in content


def test_trait_contrast_correlations_do_not_establish_semantic_equivalence() -> None:
    """C04 must qualify the correlations themselves, not a remote disclaimer."""
    content = " ".join(
        (METHODS_DIR / "table1_trait_condition_body_site_maps.md").read_text().split()
    )
    correlations = content.split("The right panel reports", 1)[1].split("##", 1)[0]
    assert "not proof of semantic equivalence or independent mechanisms" in correlations
    assert "distinct IDs" in correlations


def test_guide_credits_jed_and_links_to_pinned_interpretation_notes() -> None:
    """The historical starting evidence must remain attributable and navigable."""
    content = " ".join(BIOLOGICAL_INTERPRETATION_MD.read_text().split())
    assert "It implements [#179]" in content
    assert "[Jed's paired contrast interpretation notes at `a2b2df98`]" in content
    assert (
        "https://github.com/CultureBotAI/microbial-trait-profiling/blob/"
        "a2b2df98e9beec98b772c9cb5af7aa090d033a98/trait-visualization/"
        "trait-profile-investigation/paired_contrast_interpretation_notes.md"
    ) in content
    assert "../../trait-visualization/" not in content


def test_manuscript_statistical_audit_has_biological_companion() -> None:
    """C10 links the biological limits without implying statistics validate mechanism."""
    content = (REPO_ROOT / "MANUSCRIPT_AUDIT.md").read_text()
    assert "## Biological interpretation companion" in content
    companion = " ".join(content.split("## Biological interpretation companion", 1)[1].split())
    assert "(methods/biological_interpretation.md)" in companion
    assert "(methods/biological_interpretation_audit.md)" in companion
    assert (
        "does not provide patient-level biological measurements or validate a mechanism"
        in companion
    )
    audit = BIOLOGICAL_INTERPRETATION_AUDIT_MD.read_text()
    c10 = next(line for line in audit.splitlines() if line.startswith("| C10 "))
    assert "**Corrected here:**" in c10
    assert "companion" in c10


@pytest.mark.parametrize(
    ("surface", "revision_abbreviation"),
    [
        ("Proposed overlay, PR #173", "bf235552"),
        ("Rosalind archive/review, PR #192", "8268ccce"),
        ("Taxon contribution/influence audit, PR #194", "04aca04c"),
    ],
)
def test_claim_audit_names_the_reviewed_revision(surface: str, revision_abbreviation: str) -> None:
    """Each revision is checked in its own row, not elsewhere in the document."""
    content = BIOLOGICAL_INTERPRETATION_AUDIT_MD.read_text()
    row = next(line for line in content.splitlines() if line.startswith(f"| {surface} |"))
    full_pins = re.findall(r"`([0-9a-f]{40})`", row)
    assert full_pins, f"Missing full Git revision in {surface}"
    assert full_pins[0].startswith(revision_abbreviation)
    if "#173" in surface:
        assert "open" in row


def test_claim_audit_preserves_native_influence_and_later_adoption_limits() -> None:
    """G11's historical panel is not missing its already-delivered native reruns."""
    content = BIOLOGICAL_INTERPRETATION_AUDIT_MD.read_text()
    g11 = next(line for line in content.splitlines() if line.startswith("| G11 "))
    assert "rule-declared 24-scenario panel" in g11
    assert "native composition and fixed/recomputed-background ORA reruns" in g11
    assert "distinct correction families" in g11
    assert "native-owner review" in g11
    adoption = " ".join(content.split("## Adoption record", 1)[1].split())
    assert "pending native-owner review" in adoption
    assert "pending native integrations" not in adoption
    for finding in ("G14", "G15", "G16"):
        assert any(line.startswith(f"| {finding} ") for line in content.splitlines())


def test_table1_trait_semantics_coverage() -> None:
    """Verify that all six Table 1 traits are covered in docs/table1_trait_semantics_audit.md."""
    content = TABLE1_TRAIT_SEMANTICS_AUDIT_MD.read_text(encoding="utf-8")

    expected_traits = [
        "Facultatively aerobic",
        "Cellulose catabolic process",
        "Coccus shaped",
        "Sphere shaped",
        "Arsenate reduction",
        "Moderately halophilic",
    ]
    for trait in expected_traits:
        assert trait in content, f"Trait '{trait}' missing from Table 1 semantics audit"

    expected_trait_ids = [
        "METPO:1000608",
        "GO:0030245",
        "METPO:1000668",
        "METPO:1000683",
        "pathways:arsenate_reduction",
        "METPO:1000623",
    ]
    for tid in expected_trait_ids:
        assert tid in content, f"Trait ID '{tid}' missing from Table 1 semantics audit"


def test_governance_cross_references() -> None:
    """Verify that governance documents reference the biological interpretation standards."""
    agents_content = (REPO_ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "## Biological interpretation and review" in agents_content
    assert "methods/biological_interpretation.md" in agents_content
    assert "docs/table1_trait_semantics_audit.md" in agents_content

    contributing_content = (REPO_ROOT / "CONTRIBUTING.md").read_text(encoding="utf-8")
    assert "methods/biological_interpretation.md" in contributing_content
    assert "docs/table1_trait_semantics_audit.md" in contributing_content

    readme_content = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    assert "methods/biological_interpretation.md" in readme_content

    methods_readme_content = (REPO_ROOT / "methods" / "README.md").read_text(encoding="utf-8")
    assert "biological_interpretation.md" in methods_readme_content
    assert "biological_interpretation_audit.md" in methods_readme_content
    assert "table1_trait_semantics_audit.md" in methods_readme_content

    pr_template_content = (REPO_ROOT / ".github" / "PULL_REQUEST_TEMPLATE.md").read_text(
        encoding="utf-8"
    )
    assert "methods/biological_interpretation.md" in pr_template_content
