"""Pin named-removal comparison wording in scripts/audit_trait_coupling.py (#2015)."""

from __future__ import annotations

import ast
import importlib.util
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_trait_coupling.py"
EXPECTED_BULLET = (
    "- Evaluates the named-removal comparison against the Table 1 minimum adjacent gap."
)


def _load_cli_module(path: Path = SCRIPT) -> ModuleType:
    spec = importlib.util.spec_from_file_location("audit_trait_coupling_2015", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _assert_cli_docstring_and_source(doc: str | None, source: str) -> None:
    assert doc is not None, "scripts/audit_trait_coupling.py must define a module docstring"
    assert EXPECTED_BULLET in doc, (
        f"Missing expected named-removal bullet in docstring: {EXPECTED_BULLET!r}"
    )
    lower_source = source.lower()
    for forbidden in (
        "formal stability decision rule",
        "inter-row spacing",
        "near-twin",
    ):
        assert forbidden not in lower_source, (
            f"Forbidden legacy wording {forbidden!r} found in scripts/audit_trait_coupling.py"
        )


def test_audit_trait_coupling_docstring_uses_named_removal_comparison() -> None:
    """#2015: CLI module docstring states named-removal comparison against minimum adjacent gap."""
    source = SCRIPT.read_text(encoding="utf-8")
    module = _load_cli_module(SCRIPT)
    _assert_cli_docstring_and_source(module.__doc__, source)


@pytest.mark.parametrize(
    ("case_id", "old", "new"),
    [
        pytest.param(
            "legacy-decision-rule-bullet",
            EXPECTED_BULLET,
            "- Evaluates the formal stability decision rule against Table 1 inter-row spacing.",
            id="legacy-decision-rule-bullet",
        ),
        pytest.param(
            "missing-named-removal-bullet",
            f"{EXPECTED_BULLET}\n",
            "",
            id="missing-named-removal-bullet",
        ),
        pytest.param(
            "legacy-inter-row-spacing",
            "the Table 1 minimum adjacent gap",
            "Table 1 inter-row spacing",
            id="legacy-inter-row-spacing",
        ),
        pytest.param(
            "legacy-decision-rule-in-body-comment",
            "def parse_args(",
            "# Evaluates the Formal Stability Decision Rule\ndef parse_args(",
            id="legacy-decision-rule-in-body-comment",
        ),
        pytest.param(
            "legacy-near-twin-in-docstring",
            "Removes coccus + sphere simultaneously",
            "Removes near-twin coccus + sphere simultaneously",
            id="legacy-near-twin-in-docstring",
        ),
        pytest.param(
            "missing-module-docstring",
            '"""Audit trait coupling sensitivity via LOTO, LCO, and multi-arm evaluation.\n\n'
            "Addresses Issue #260:\n"
            "- Leave-One-Trait-Out (LOTO): Computes 6x6 matrix of delta effect and delta q "
            "across Table 1 traits.\n"
            "- Leave-Cluster-Out (LCO): Removes coccus + sphere simultaneously and measures "
            "shifts across 3 arms.\n"
            "- Three-Arm Sensitivity: Compares Paired Composition, Directional Sign-Test,\n"
            "  and ORA Hypergeometric arms.\n"
            "- Evaluates the named-removal comparison against the Table 1 minimum adjacent gap.\n"
            '"""\n',
            "# Module docstring omitted\n",
            id="missing-module-docstring",
        ),
    ],
)
def test_audit_trait_coupling_docstring_guard_rejects_legacy_wording(
    case_id: str, old: str, new: str
) -> None:
    """#2015: guard rejects legacy stability-rule, inter-row-spacing, or near-twin wording."""
    baseline = SCRIPT.read_text(encoding="utf-8")
    assert baseline.count(old) == 1, f"Expected unique anchor for {case_id}"
    mutated = baseline.replace(old, new, 1)
    parsed_doc = ast.get_docstring(ast.parse(mutated), clean=False)
    try:
        _assert_cli_docstring_and_source(parsed_doc, mutated)
    except AssertionError:
        return
    raise AssertionError(f"Docstring guard accepted legacy wording for {case_id}")
