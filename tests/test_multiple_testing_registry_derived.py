"""Require visible same-family BY arithmetic for registry derivations (#1490)."""

from __future__ import annotations

from pathlib import Path
from textwrap import indent
from typing import Any

import pytest

from test_multiple_testing_registry import (
    DEFAULT_LIMITS,
    DiscoveryLimitError,
    Limits,
    SummaryBudget,
    discover,
    example_registry,
    family_procedure_errors,
    registry_errors,
    summarize_source,
)

DERIVATION_ERROR = "derived procedures lack visible arithmetic"
ROOT = Path(__file__).resolve().parents[1]
PREFIX = "m = len(p)\nharmonic_sum = np.sum(1.0 / np.arange(1, m + 1))\nq = benjamini_hochberg(p)\n"


def source(body: str) -> str:
    return (
        "import numpy as np\n"
        "from mtp.enrichment_stats import benjamini_hochberg\n"
        "def screen(p, other):\n" + indent(body, "    ")
    )


def by_registry() -> dict[str, Any]:
    registry = example_registry(("src/x.py", "screen", ["BH", "BY"]))
    registry["families"][0]["anchors"][0]["derived_procedures"] = {
        "BY": "BH multiplied by the same family's harmonic sum and clipped to [0, 1]."
    }
    registry["literal_procedure_counts"][0]["counts"] = {"BH": 1}
    return registry


@pytest.mark.parametrize(
    "body",
    [
        "return benjamini_hochberg(p)\n",
        "return np.clip(benjamini_hochberg(p) * 2, 0, 1)\n",
        "return np.clip(benjamini_hochberg(p) * harmonic(p), 0, 1)\n",
        PREFIX + "return q * harmonic_sum\n",
        PREFIX + "return np.clip(q * harmonic_sum, 0.1, 1)\n",
        PREFIX + "return np.clip(q * harmonic_sum, 0, 2)\n",
        PREFIX + "return np.clip(p * harmonic_sum, 0, 1)\n",
        PREFIX + "q = p\nreturn np.clip(q * harmonic_sum, 0, 1)\n",
        PREFIX + "harmonic_sum = 2\nreturn np.clip(q * harmonic_sum, 0, 1)\n",
        PREFIX + "q += p\nreturn np.clip(q * harmonic_sum, 0, 1)\n",
    ],
    ids=[
        "reason-only",
        "constant-multiplier",
        "opaque-harmonic-helper",
        "no-clip",
        "wrong-lower-clip",
        "wrong-upper-clip",
        "raw-p",
        "overwritten-bh",
        "overwritten-factor",
        "augmented-bh",
    ],
)
def test_derived_by_requires_visible_bh_harmonic_and_clip(body: str) -> None:
    errors = registry_errors({"src/x.py": source(body)}, by_registry())
    assert any(DERIVATION_ERROR in error for error in errors), errors


@pytest.mark.parametrize(
    "harmonic",
    [
        "np.sum(2.0 / np.arange(1, m + 1))",
        "np.sum(1.0 / np.arange(0, m + 1))",
        "np.sum(1.0 / np.arange(1, m))",
        "np.sum(1.0 / np.arange(1, m + 2))",
        "np.sum(1.0 / np.arange(1, m + 1, 2))",
        "sum(2.0 / i for i in range(1, m + 1))",
        "sum(1.0 / i for i in range(0, m + 1))",
        "sum(1.0 / i for i in range(1, m + 1) if i > 1)",
    ],
)
def test_malformed_harmonic_factor_is_not_a_by_derivation(harmonic: str) -> None:
    body = (
        f"m = len(p)\nharmonic_sum = {harmonic}\n"
        "return np.clip(benjamini_hochberg(p) * harmonic_sum, 0, 1)\n"
    )
    errors = registry_errors({"src/x.py": source(body)}, by_registry())
    assert any(DERIVATION_ERROR in error for error in errors), errors


@pytest.mark.parametrize(
    "body",
    [
        "m = len(other)\nharmonic_sum = np.sum(1.0 / np.arange(1, m + 1))\n"
        "return np.clip(benjamini_hochberg(p) * harmonic_sum, 0, 1)\n",
        "m = len(p)\nharmonic_sum = np.sum(1.0 / np.arange(1, m + 1))\n"
        "p = other\nreturn np.clip(benjamini_hochberg(p) * harmonic_sum, 0, 1)\n",
        "m = len(p)\nharmonic_sum = np.sum(1.0 / np.arange(1, m + 1))\n"
        "p_r = np.ones(len(other))\n"
        "return np.clip(benjamini_hochberg(p_r) * harmonic_sum, 0, 1)\n",
        "m = len(p)\nm = len(other)\n"
        "harmonic_sum = np.sum(1.0 / np.arange(1, m + 1))\n"
        "return np.clip(benjamini_hochberg(p) * harmonic_sum, 0, 1)\n",
        "m = len(p)\nm += 1\n"
        "harmonic_sum = np.sum(1.0 / np.arange(1, m + 1))\n"
        "return np.clip(benjamini_hochberg(p) * harmonic_sum, 0, 1)\n",
    ],
    ids=[
        "other-family",
        "rebound-family-vector",
        "different-allocated-length",
        "rebound-scalar-extent",
        "augmented-scalar-extent",
    ],
)
def test_harmonic_extent_must_match_the_bh_vector(body: str) -> None:
    errors = registry_errors({"src/x.py": source(body)}, by_registry())
    assert any(DERIVATION_ERROR in error for error in errors), errors


@pytest.mark.parametrize(
    "tail",
    [
        "alias = q\nalias += p\nreturn np.clip(q * harmonic_sum, 0, 1)\n",
        "alias = q\nq += p\nreturn np.clip(alias * harmonic_sum, 0, 1)\n",
        "alias = q\nalias[0] = 0\nreturn np.clip(q * harmonic_sum, 0, 1)\n",
        "with manager() as q:\n    return np.clip(q * harmonic_sum, 0, 1)\n",
    ],
    ids=["mutated-alias", "mutated-original", "subscript-alias", "with-target"],
)
def test_mutation_or_context_binding_invalidates_bh_fact(tail: str) -> None:
    errors = registry_errors({"src/x.py": source(PREFIX + tail)}, by_registry())
    assert any(DERIVATION_ERROR in error for error in errors), errors


@pytest.mark.parametrize(
    "body",
    [
        PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n",
        PREFIX + "scaled = q * harmonic_sum\nreturn np.clip(scaled, 0.0, 1.0)\n",
        PREFIX + "alias = q\nreturn np.clip(harmonic_sum * alias, 0, 1)\n",
        "m = len(p)\nh = float(np.sum(1.0 / np.arange(1, m + 1)))\n"
        "return np.clip(benjamini_hochberg(p) * h, 0, 1)\n",
        "m = len(p)\nh = sum(1.0 / i for i in range(1, m + 1))\n"
        "return np.clip(benjamini_hochberg(p) * h, 0, 1)\n",
        "m = len(p)\nh = np.sum(1.0 / np.arange(1, m + 1))\n"
        "p_r = np.ones(m, dtype=float)\nq = benjamini_hochberg(p_r)\n"
        "return np.clip(q * h, 0, 1)\n",
        "m = len(p)\nh = np.sum(1.0 / np.arange(1, m + 1))\n"
        "if len(other):\n"
        "    q = benjamini_hochberg(p)\n"
        "    return np.clip(q * h, 0, 1)\n"
        "return p\n",
    ],
    ids=[
        "assigned",
        "separate-clip",
        "alias-and-reversed-product",
        "float-sum",
        "generator-sum",
        "allocated-vector",
        "conditional-derivation",
    ],
)
def test_visible_same_family_by_arithmetic_is_accepted(body: str) -> None:
    assert not registry_errors({"src/x.py": source(body)}, by_registry())


def test_nested_derivation_does_not_credit_the_outer_anchor() -> None:
    body = (
        "def nested(p):\n"
        "    m = len(p)\n"
        "    h = np.sum(1.0 / np.arange(1, m + 1))\n"
        "    return np.clip(benjamini_hochberg(p) * h, 0, 1)\n"
        "return benjamini_hochberg(p)\n"
    )
    registry = by_registry()
    nested_anchor = {
        "path": "src/x.py",
        "function": "screen.nested",
        "procedures": ["BH", "BY"],
        "derived_procedures": {"BY": "Visible nested calculation."},
    }
    registry["families"][0]["anchors"].append(nested_anchor)
    registry["literal_procedure_counts"].append(
        {"path": "src/x.py", "function": "screen.nested", "counts": {"BH": 1}}
    )
    errors = registry_errors({"src/x.py": source(body)}, registry)
    assert any(DERIVATION_ERROR in error and "'screen')" in error for error in errors), errors


def test_derivation_cannot_be_assembled_across_exclusive_branches() -> None:
    body = (
        "m = len(p)\nh = np.sum(1.0 / np.arange(1, m + 1))\n"
        "if len(other):\n"
        "    q = benjamini_hochberg(p)\n"
        "else:\n"
        "    q = p\n"
        "return np.clip(q * h, 0, 1)\n"
    )
    errors = registry_errors({"src/x.py": source(body)}, by_registry())
    assert any(DERIVATION_ERROR in error for error in errors), errors


def test_by_arithmetic_does_not_establish_a_derived_bh_procedure() -> None:
    registry = by_registry()
    registry["families"][0]["anchors"][0]["derived_procedures"] = {"BH": "A reason."}
    registry["literal_procedure_counts"][0]["counts"] = {"BY": 1}
    errors = registry_errors(
        {"src/x.py": "def screen(p):\n    return false_discovery_control(p, method='by')\n"},
        registry,
    )
    assert any(DERIVATION_ERROR in error for error in errors), errors


def test_derived_by_does_not_replace_the_literal_bh_inventory() -> None:
    registry = by_registry()
    registry["literal_procedure_counts"][0]["counts"] = {"BH": 2}
    errors = registry_errors(
        {"src/x.py": source(PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n")}, registry
    )
    assert any("literal procedure count mismatch" in error for error in errors), errors


def test_compact_derivation_fact_is_charged_to_the_summary_budget() -> None:
    valid = source(PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n")
    invalid = valid.replace("harmonic_sum, 0, 1)", "harmonic_sum, 0, 2)")
    valid_budget, invalid_budget = SummaryBudget(DEFAULT_LIMITS), SummaryBudget(DEFAULT_LIMITS)
    summary = summarize_source("src/x.py", valid, valid_budget)
    invalid_summary = summarize_source("src/x.py", invalid, invalid_budget)
    assert summary.derived == ("screen",)
    assert invalid_summary.derived == ()
    assert valid_budget.items == invalid_budget.items + 1
    assert valid_budget.bytes > invalid_budget.bytes
    with pytest.raises(DiscoveryLimitError, match="compact summary item budget"):
        summarize_source(
            "src/x.py", valid, SummaryBudget(Limits(summary_items=invalid_budget.items))
        )


@pytest.mark.parametrize(
    ("path", "function"),
    [
        (
            "src/mtp/composition_calibration/adapters/jed_composition.py",
            "JedCompositionAdapter.run_inference",
        ),
        ("src/mtp/composition_calibration/composition.py", "run_composition_calibration"),
    ],
)
def test_existing_arithmetic_anchors_expose_by_without_importing_modules(
    path: str, function: str
) -> None:
    # Two exact small source files, not repository_sources or a repository scan.
    sites, _ = discover({path: (ROOT / path).read_text(encoding="utf-8")})
    site = sites[(path, function)]
    assert "BY" in site.get("derived_procedures", set())
    registry = example_registry((path, function, ["BH", "BY"]))
    family = registry["families"][0]
    family["anchors"][0]["derived_procedures"] = {"BY": "Visible anchored arithmetic."}
    assert not family_procedure_errors(family, sites)
