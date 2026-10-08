from pathlib import Path


def test_biological_interpretation_guide_exists() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    guide_path = repo_root / "methods" / "biological_interpretation.md"
    assert guide_path.exists(), f"Guide missing: {guide_path}"
    content = guide_path.read_text(encoding="utf-8")

    # Core sections and frameworks
    assert "Evidence chain" in content
    assert "Four different meanings of direction" in content
    assert "Signature direction" in content
    assert "Relative annotation contrast" in content
    assert "Enrichment against a background" in content
    assert "Direct biological change" in content
    assert "Minimum result and caption contract" in content
    assert "Defensible wording" in content
    assert "Template for a result" in content
    assert "Rosalind" in content


def test_biological_interpretation_audit_exists() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    audit_path = repo_root / "methods" / "biological_interpretation_audit.md"
    assert audit_path.exists(), f"Audit missing: {audit_path}"
    content = audit_path.read_text(encoding="utf-8")

    # Audit coverage
    for i in range(1, 12):
        assert f"C{i:02d}" in content, f"Missing C{i:02d} in audit"
    for i in range(1, 11):
        assert f"P{i:02d}" in content, f"Missing P{i:02d} in audit"
    for i in range(1, 10):
        assert f"R{i:02d}" in content, f"Missing R{i:02d} in audit"
    for i in range(1, 14):
        assert f"G{i:02d}" in content, f"Missing G{i:02d} in audit"


def test_generator_scripts_no_unsupported_claims() -> None:
    repo_root = Path(__file__).resolve().parent.parent

    enrichment_script = (repo_root / "scripts" / "generate_enrichment_summary.py").read_text(
        encoding="utf-8"
    )
    assert "clinical functional burden" not in enrichment_script
    assert "Hits/100 Samples" not in enrichment_script

    depletion_script = (repo_root / "scripts" / "generate_depletion_summary.py").read_text(
        encoding="utf-8"
    )
    assert "Losses/100 Samples" not in depletion_script
    assert "Losses / 100 Samples" not in depletion_script
    assert "Functional Blind Spots" not in depletion_script


def test_milestone_docs_no_unsupported_claims() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    milestone_doc = (
        repo_root / "docs" / "milestone_taxon_model_before_after_summary.md"
    ).read_text(encoding="utf-8")
    exchangeability_claim = "Restores true statistical exchangeability of biological observations"
    assert exchangeability_claim not in milestone_doc
    assert "genuine biological drivers" not in milestone_doc
    assert "methods/biological_interpretation.md" in milestone_doc


def test_project_documentation_governance() -> None:
    repo_root = Path(__file__).resolve().parent.parent

    readme = (repo_root / "README.md").read_text(encoding="utf-8")
    assert "methods/biological_interpretation.md" in readme
    assert "methods/biological_interpretation_audit.md" in readme

    agents = (repo_root / "AGENTS.md").read_text(encoding="utf-8")
    assert "methods/biological_interpretation.md" in agents

    contributing = (repo_root / "CONTRIBUTING.md").read_text(encoding="utf-8")
    assert "methods/biological_interpretation.md" in contributing

    pr_template = (repo_root / ".github" / "PULL_REQUEST_TEMPLATE.md").read_text(encoding="utf-8")
    assert "methods/biological_interpretation.md" in pr_template
