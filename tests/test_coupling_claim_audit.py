"""Pin source-qualified coupling claims and neutral overlap wording (#1999-#2001, #2014)."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
METHODS = ROOT / "methods" / "trait_coupling_sensitivity.md"
AUDIT = ROOT / "methods" / "biological_interpretation_audit.md"
REPORT = ROOT / "data" / "derived" / "trait_coupling_report.txt"

HISTORICAL_BINARIZED_BULLETS: tuple[tuple[str, str, tuple[str, ...]], ...] = (
    (
        "**Key LOTO Observations:**",
        (
            "- **Magnitude of perturbation by arm:** "
            "In the historical binarized evaluation, in the composition arm, "
        ),
        (
            r"$|\Delta \text{effect}| = 0.000636$",
            r"$|\Delta \text{share\_up}| = 0.002183$",
            r"$|\Delta \text{mean\_fold}| = 0.024908$",
        ),
    ),
    (
        "**Key LOTO Observations:**",
        (
            "- **Statistical significance (focal-set BH sensitivity check):** "
            "In the historical binarized evaluation, all evaluated traits "
        ),
        (
            r"$q \ll 0.05$ for five traits",
            r"$q \approx 0.15$ for all-profile *facultatively aerobic*",
        ),
    ),
    (
        "### Cross-Arm Comparison",
        (
            "- **Composition Arm:** "
            "In the historical binarized evaluation, maximum absolute shifts are "
        ),
        (
            r"$|\Delta \text{effect}| \le 0.000636$ in LOTO",
            r"$\le 0.001008$ in LCO",
            "0.002040",
        ),
    ),
    (
        "### Cross-Arm Comparison",
        (
            "- **Directional Arm:** "
            "In the historical binarized evaluation, study share up shifts by at most "
        ),
        (
            r"$|\Delta \text{share\_up}| \le 0.002183$ in LOTO",
            r"$\le 0.002736$ in LCO",
            "from $916$ to $917$ under LCO",
        ),
    ),
    (
        "### Cross-Arm Comparison",
        (
            "- **ORA Arm:** "
            "In the historical binarized evaluation, mean fold enrichment shifts by up to "
        ),
        (
            r"$|\Delta \text{mean\_fold}| \le 0.024908$ in LOTO",
            r"$\le 0.022271$ in LCO",
            "from $33$ to $38$ under LCO",
        ),
    ),
)


def _assert_historical_binarized_bullets(text: str) -> None:
    for heading in ("**Key LOTO Observations:**", "### Cross-Arm Comparison"):
        assert text.count(f"\n{heading}\n") == 1, f"Missing or duplicate section heading: {heading}"
        section = text.split(f"\n{heading}\n", 1)[1].split("\n---", 1)[0]
        bullets = [line for line in section.splitlines() if line.startswith("- ")]
        expected_specs = [spec for spec in HISTORICAL_BINARIZED_BULLETS if spec[0] == heading]
        assert len(bullets) == len(expected_specs), f"Unexpected bullet count under {heading}"
        for bullet, (_sec, expected_lead, fragments) in zip(bullets, expected_specs, strict=True):
            assert bullet.startswith(expected_lead), (
                f"Bullet under {heading} must start with {expected_lead!r}: {bullet}"
            )
            for fragment in fragments:
                assert fragment in bullet, f"Bullet under {heading} omits {fragment!r}: {bullet}"


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
    _assert_historical_binarized_bullets(text)


def test_methods_qualify_historical_binarized_loto_and_cross_arm_bullets() -> None:
    """#2014: LOTO observation bullets and Cross-Arm Comparison bullets qualify historical scope."""
    _assert_historical_binarized_bullets(METHODS.read_text(encoding="utf-8"))


@pytest.mark.parametrize(
    ("case_id", "old", "new"),
    [
        pytest.param(
            "loto-magnitude-unqualified",
            (
                "- **Magnitude of perturbation by arm:** "
                "In the historical binarized evaluation, in the composition arm, "
            ),
            "- **Magnitude of perturbation by arm:** In the composition arm, ",
            id="loto-magnitude-unqualified",
        ),
        pytest.param(
            "loto-magnitude-missing-binarized",
            (
                "- **Magnitude of perturbation by arm:** "
                "In the historical binarized evaluation, in the composition arm, "
            ),
            (
                "- **Magnitude of perturbation by arm:** "
                "In the historical evaluation, in the composition arm, "
            ),
            id="loto-magnitude-missing-binarized",
        ),
        pytest.param(
            "loto-significance-unqualified",
            (
                "- **Statistical significance (focal-set BH sensitivity check):** "
                "In the historical binarized evaluation, all evaluated traits "
            ),
            (
                "- **Statistical significance (focal-set BH sensitivity check):** "
                "All evaluated traits "
            ),
            id="loto-significance-unqualified",
        ),
        pytest.param(
            "loto-significance-missing-binarized",
            (
                "- **Statistical significance (focal-set BH sensitivity check):** "
                "In the historical binarized evaluation, all evaluated traits "
            ),
            (
                "- **Statistical significance (focal-set BH sensitivity check):** "
                "In the historical evaluation, all evaluated traits "
            ),
            id="loto-significance-missing-binarized",
        ),
        pytest.param(
            "cross-arm-directional-unqualified",
            (
                "- **Directional Arm:** "
                "In the historical binarized evaluation, study share up shifts by at most "
            ),
            "- **Directional Arm:** Study share up shifts by at most ",
            id="cross-arm-directional-unqualified",
        ),
        pytest.param(
            "cross-arm-directional-missing-binarized",
            (
                "- **Directional Arm:** "
                "In the historical binarized evaluation, study share up shifts by at most "
            ),
            (
                "- **Directional Arm:** "
                "In the historical evaluation, study share up shifts by at most "
            ),
            id="cross-arm-directional-missing-binarized",
        ),
        pytest.param(
            "cross-arm-composition-unqualified",
            (
                "- **Composition Arm:** "
                "In the historical binarized evaluation, maximum absolute shifts are "
            ),
            "- **Composition Arm:** Maximum absolute shifts are ",
            id="cross-arm-composition-unqualified",
        ),
        pytest.param(
            "cross-arm-ora-unqualified",
            (
                "- **ORA Arm:** "
                "In the historical binarized evaluation, mean fold enrichment shifts by up to "
            ),
            "- **ORA Arm:** Mean fold enrichment shifts by up to ",
            id="cross-arm-ora-unqualified",
        ),
        pytest.param(
            "loto-magnitude-missing-stat",
            r"$|\Delta \text{effect}| = 0.000636$",
            r"$|\Delta \text{effect}| = 0.005937$",
            id="loto-magnitude-missing-stat",
        ),
        pytest.param(
            "loto-significance-missing-q",
            r"$q \approx 0.15$ for all-profile *facultatively aerobic*",
            r"$q \approx 1.77 \times 10^{-10}$ for all-profile *facultatively aerobic*",
            id="loto-significance-missing-q",
        ),
        pytest.param(
            "cross-arm-directional-missing-lco-bound",
            r"$\le 0.002736$ in LCO",
            r"$\le 0.006553$ in LCO",
            id="cross-arm-directional-missing-lco-bound",
        ),
        pytest.param(
            "extra-unqualified-loto-bullet",
            "\n---\n\n## 4. Leave-Cluster-Out (LCO) Sensitivity",
            (
                "\n- **Extra Observation:** Unqualified bullet.\n\n"
                "---\n\n## 4. Leave-Cluster-Out (LCO) Sensitivity"
            ),
            id="extra-unqualified-loto-bullet",
        ),
        pytest.param(
            "duplicate-loto-heading",
            "\n**Key LOTO Observations:**\n",
            "\n**Key LOTO Observations:**\n\n**Key LOTO Observations:**\n",
            id="duplicate-loto-heading",
        ),
    ],
)
def test_historical_binarized_bullet_guard_rejects_unqualified_or_incomplete_bullets(
    case_id: str, old: str, new: str
) -> None:
    """#2014: removing historical binarized or pinned stats fails by AssertionError."""
    baseline = METHODS.read_text(encoding="utf-8")
    assert baseline.count(old) == 1, f"Expected unique anchor for {case_id}"
    mutated = baseline.replace(old, new, 1)
    try:
        _assert_historical_binarized_bullets(mutated)
    except AssertionError:
        return
    raise AssertionError(f"Bullet guard accepted invalid methods text for {case_id}")
