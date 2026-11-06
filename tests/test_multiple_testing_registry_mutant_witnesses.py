"""Behavioral witnesses for the 24 independent review mutants (#1578)."""

from __future__ import annotations

import pytest

from test_multiple_testing_registry import registry_errors
from test_multiple_testing_registry_derived import DERIVATION_ERROR, by_registry, source

HEAD = "m = len(p)\nh = np.sum(1.0 / np.arange(1, m + 1))\n"
BH_HEAD = HEAD + "q = benjamini_hochberg(p)\n"
CLIP_Q = "return np.clip(q * h, 0, 1)\n"


@pytest.mark.parametrize(
    "body",
    [
        pytest.param(
            "m = len(p)\nsum = lambda values: 2.0\n"
            "h = sum(1.0 / i for i in range(1, m + 1))\n"
            "return np.clip(benjamini_hochberg(p) * h, 0, 1)\n",
            id="target-ignores-local-rebinding",
        ),
        pytest.param(
            HEAD + "return np.clip(float(benjamini_hochberg(p)) * h, 0, 1)\n",
            id="float-wraps-anything",
        ),
        pytest.param(
            "q = benjamini_hochberg(p)\nm = len(q)\n"
            "h = np.sum(1.0 / np.arange(1, m + 1))\n" + CLIP_Q,
            id="len-of-any-fact",
        ),
        pytest.param(
            HEAD + "p_r = np.ones(m, like=other)\n"
            "return np.clip(benjamini_hochberg(p_r) * h, 0, 1)\n",
            id="ones-any-keywords",
        ),
        pytest.param(
            HEAD + "return np.clip(benjamini_hochberg(p, family_size=2 * m) * h, 0, 1)\n",
            id="bh-call-keywords-allowed",
        ),
        pytest.param(
            HEAD + "return np.clip(benjamini_hochberg(h) * h, 0, 1)\n",
            id="bh-of-any-fact",
        ),
        pytest.param(
            "m = len(p)\nh = np.sum(1.0 / np.arange(1, m + 1, step=2))\n"
            "return np.clip(benjamini_hochberg(p) * h, 0, 1)\n",
            id="arange-keywords-allowed",
        ),
        pytest.param(
            "m = len(p)\nh = np.sum(1.0 / np.arange(1, m - 1))\n"
            "return np.clip(benjamini_hochberg(p) * h, 0, 1)\n",
            id="arange-stop-any-operator",
        ),
        pytest.param(
            "m = len(p)\nh = np.sum(1.0 * np.arange(1, m + 1))\n"
            "return np.clip(benjamini_hochberg(p) * h, 0, 1)\n",
            id="np-harmonic-any-operator",
        ),
        pytest.param(
            "m = len(p)\nh = sum(1.0 * i for i in range(1, m + 1))\n"
            "return np.clip(benjamini_hochberg(p) * h, 0, 1)\n",
            id="generator-any-operator",
        ),
        pytest.param(
            "m = len(p)\nh = sum(1.0 / m for i in range(1, m + 1))\n"
            "return np.clip(benjamini_hochberg(p) * h, 0, 1)\n",
            id="generator-any-divisor",
        ),
        pytest.param(
            HEAD + "return np.clip(benjamini_hochberg(p) / h, 0, 1)\n",
            id="product-any-operator",
        ),
        pytest.param(
            BH_HEAD + "return np.clip(q * h, 0, 1, out=other)\n",
            id="clip-keywords-allowed",
        ),
        pytest.param(
            BH_HEAD + "q, _ = p, None\n" + CLIP_Q,
            id="tuple-targets-not-rebound",
        ),
        pytest.param(
            BH_HEAD + "alias = q\nalias[0] += 1\n" + CLIP_Q,
            id="augassign-subscript-ignored",
        ),
        pytest.param(
            HEAD + "q = p\nif other:\n    q = benjamini_hochberg(p)\n" + CLIP_Q,
            id="if-body-shares-env",
        ),
        pytest.param(
            "m = len(data)\nh = np.sum(1.0 / np.arange(1, m + 1))\n"
            "if other:\n    data = other\n"
            "return np.clip(benjamini_hochberg(data) * h, 0, 1)\n",
            id="merge-drops-one-sided-names",
        ),
        pytest.param(
            HEAD + "q = p\nfor item in other:\n    q = benjamini_hochberg(p)\n" + CLIP_Q,
            id="loop-facts-leak-out",
        ),
        pytest.param(
            BH_HEAD + "for q in other:\n    by = np.clip(q * h, 0, 1)\nreturn by\n",
            id="for-target-not-rebound",
        ),
        pytest.param(
            BH_HEAD + "for item in other:\n    pass\nelse:\n    q = p\n" + CLIP_Q,
            id="loop-else-skipped",
        ),
        pytest.param(
            BH_HEAD + "try:\n    q = p\nexcept ValueError:\n    pass\n" + CLIP_Q,
            id="try-stores-not-invalidated",
        ),
        pytest.param(
            HEAD + "try:\n    q = benjamini_hochberg(p)\n"
            "    by = np.clip(q * h, 0, 1)\nexcept ValueError:\n    pass\nreturn by\n",
            id="try-body-interpreted",
        ),
    ],
)
def test_review_mutant_witnesses_refuse_unproved_derivations(body: str) -> None:
    errors = registry_errors({"src/x.py": source(body)}, by_registry())
    assert any(DERIVATION_ERROR in error for error in errors), errors


def test_with_body_is_interpreted_for_a_valid_bound_derivation() -> None:
    """A real binding in a with block is a witness, even without an as target."""
    body = BH_HEAD + "with manager():\n    by = np.clip(q * h, 0, 1)\nreturn by\n"
    assert registry_errors({"src/x.py": source(body)}, by_registry()) == []


def test_unproved_claim_is_also_absent_from_the_literal_site() -> None:
    """An unsupported claim cannot discharge the independently absent procedure."""
    errors = registry_errors({"src/x.py": source("return benjamini_hochberg(p)\n")}, by_registry())
    assert len(errors) == 2, errors
    assert sum(DERIVATION_ERROR in error for error in errors) == 1, errors
    assert sum("procedures absent from all-literal site" in error for error in errors) == 1, errors
