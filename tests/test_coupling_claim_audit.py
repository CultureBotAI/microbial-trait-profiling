"""Pin source-qualified coupling claims and neutral overlap wording (#1999-#2001)."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METHODS = ROOT / "methods" / "trait_coupling_sensitivity.md"
AUDIT = ROOT / "methods" / "biological_interpretation_audit.md"
REPORT = ROOT / "data" / "derived" / "trait_coupling_report.txt"


def _g16() -> str:
    return next(
        line for line in AUDIT.read_text(encoding="utf-8").splitlines() if line.startswith("| G16 ")
    )


def test_overlap_descriptions_preserve_descriptive_scope_and_caveats() -> None:
    text = METHODS.read_text(encoding="utf-8")
    diagram = text.split("```mermaid", 1)[1].split("```", 1)[0]
    association = text.split("2. **Annotation-count association:**", 1)[1].split("\n---", 1)[0]
    assert "Coccus and sphere are annotated to largely overlapping taxa" in diagram
    assert "Semantic equivalence and the cause of association are not established" in diagram
    assert (
        "Both labels are annotated to largely overlapping taxa in the trait matrix." in association
    )
    assert "do not establish semantic equivalence or the cause of the association" in association


def test_g16_reports_saved_count_profile_values_with_source_and_scope() -> None:
    report = REPORT.read_text(encoding="utf-8")
    g16 = _g16()
    # These are values already printed in the saved report, not recomputed results.
    patterns = (
        r"maximum \|Delta effect\| of ([0-9.]+)",
        r"up to ([0-9.]+) under LCO",
        r"Minimum adjacent row distance:\s+([0-9.]+)",
        r"Delta effect\(coccus\):\s+([+0-9.]+)",
        r"Delta effect\(sphere\):\s+([+0-9.]+)",
        r"maximum \|Delta share_up\| is ([0-9.]+)",
        r"maximum \|Delta mean_fold\| is ([0-9.]+)",
    )
    for pattern in patterns:
        match = re.search(pattern, report)
        assert match is not None, f"Missing saved report value: {pattern}"
        assert match.group(1) in g16, f"G16 omits saved report value: {match.group(1)}"
    assert "../data/derived/trait_coupling_report.txt" in g16
    assert "current saved raw-count report" in g16
    assert "historical binarized-profile" in g16
    assert "not current raw-count results" in g16
    assert "exceed" in g16
    assert "Three LOTO" in g16
    assert "one LCO" in g16


def test_g16_does_not_convert_sensitivity_or_correlation_into_biological_claims() -> None:
    g16 = _g16()
    for unsupported in (
        "due to semantic/ontology redundancy",
        "structurally stable against denominator closure",
        "passing the formal stability decision rule",
    ):
        assert unsupported not in g16
    assert "do not establish semantic equivalence or independent biological mechanisms" in g16
    assert "do not establish unchanged rankings for every leave-out" in g16


def test_methods_state_named_removal_comparisons_without_stability_verdicts() -> None:
    text = METHODS.read_text(encoding="utf-8")
    assert not re.search(r"\b(?:near[- ]twins?|twins?|modestly)\b", text, flags=re.IGNORECASE)
    for framing in (
        "Formal Stability Decision Rule",
        "PASSES CRITERION",
        "passes rule",
        "comparisons pass the historical criterion",
        "the two cannot be reported as separate convergent findings",
    ):
        assert framing not in text
    assert "Named-Removal Reference-Gap Comparison" in text
    assert "historical binarized" in text
    assert "does not establish unchanged rankings for every trait after every exclusion" in text
