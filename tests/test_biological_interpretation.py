"""Tests for project-wide biological interpretation guide and governance docs."""

from __future__ import annotations

from pathlib import Path

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

    section_5 = content.split(
        "## 5. Decisions, Approved Wording, and Unsupported Interpretations", 1
    )[1]
    metpo_608_row = next(
        line for line in section_5.splitlines() if line.startswith("| `METPO:1000608` |")
    )
    assert "pre-#1124" in metpo_608_row
    assert "(capable of growth in both the presence and absence of oxygen)" in metpo_608_row
    assert "facultative (aerobe or anaerobe)" in metpo_608_row
    assert "annotation for growth with or without oxygen" in metpo_608_row
    assert "clean_trait_label" not in content
    delta_z_md = (REPO_ROOT / "methods" / "delta_z_derived.md").read_text(encoding="utf-8")
    assert "facultative (aerobe or anaerobe)" in delta_z_md
    assert "keep the KG node label" in delta_z_md


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
