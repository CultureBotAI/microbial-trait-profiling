"""Control-flow and binding regressions from PR #1570 review (#1575-#1577)."""

from __future__ import annotations

from textwrap import indent

import pytest

from test_multiple_testing_registry import registry_errors
from test_multiple_testing_registry_derived import DERIVATION_ERROR, by_registry, source

HEAD = "m = len(p)\nh = np.sum(1.0 / np.arange(1, m + 1))\n"
BH = "q = benjamini_hochberg(p)\n"
RESULT = "np.clip(q * h, 0, 1)"
DIRECT = "np.clip(benjamini_hochberg(p) * h, 0, 1)"


@pytest.mark.parametrize(
    "body",
    [
        HEAD + f"for chunk in other:\n    by = {DIRECT}\n    p = chunk\nreturn by\n",
        HEAD + f"for chunk in other:\n    by = {DIRECT}\n    h = 2.0\nreturn by\n",
        HEAD + BH + f"for chunk in other:\n    by = {RESULT}\n    q = chunk\nreturn by\n",
        HEAD + f"while other:\n    by = {DIRECT}\n    p = other.pop()\nreturn by\n",
        # The first iteration changes b, the second changes a, and only then does
        # q lose its BH fact. One extra pass is not a fixed-point computation.
        HEAD + BH + "a = q\nb = q\nfor chunk in other:\n"
        f"    by = {RESULT}\n    q = a\n    a = b\n    b = p\nreturn by\n",
    ],
    ids=["for-vector", "for-factor", "for-bh", "while-vector", "multi-pass-fixed-point"],
)
def test_loop_carried_facts_cannot_establish_by(body: str) -> None:
    errors = registry_errors({"src/x.py": source(body)}, by_registry())
    assert any(DERIVATION_ERROR in error for error in errors), errors


@pytest.mark.parametrize(
    "rebinding",
    [
        "(q := p)\n",
        "if (q := p) is None:\n    pass\n",
        "while (q := p) is None:\n    break\n",
        "for item in (q := other):\n    pass\n",
        "from os import sep as q\n",
        "import os as q\n",
        "try:\n    from os import sep as q\nexcept ImportError:\n    pass\n",
        "try:\n    import os.path as q\nexcept ImportError:\n    pass\n",
        "try:\n    raise ValueError()\nexcept ValueError as q:\n    pass\n",
        "match other:\n    case q:\n        pass\n",
        "match other:\n    case [q, *_]:\n        pass\n",
        "match other:\n    case [*q]:\n        pass\n",
        "match other:\n    case {'key': value, **q}:\n        pass\n",
        "del q\n",
        "v = q[:]\nv[:] = p\n",
        "q[:][0] = 5\n",
        "q.flat[:] = p\n",
        "q.flat = p\n",
        "del q[:][0]\n",
    ],
    ids=[
        "walrus-expression",
        "walrus-if-test",
        "walrus-while-test",
        "walrus-for-iterable",
        "from-import",
        "import",
        "from-import-in-try",
        "import-in-try",
        "exception-handler-name",
        "match-as",
        "match-sequence",
        "match-star",
        "match-mapping-rest",
        "del-name",
        "view-write",
        "nested-subscript-write",
        "attribute-subscript-write",
        "attribute-write",
        "nested-subscript-delete",
    ],
)
def test_rebinding_and_writes_invalidate_bh_fact(rebinding: str) -> None:
    errors = registry_errors(
        {"src/x.py": source(HEAD + BH + rebinding + f"return {RESULT}\n")}, by_registry()
    )
    assert any(DERIVATION_ERROR in error for error in errors), errors


@pytest.mark.parametrize(
    "body",
    [
        "m = len(p)\nif (m := len(other)):\n    pass\n"
        "h = np.sum(1.0 / np.arange(1, m + 1))\n" + f"return {DIRECT}\n",
        HEAD + "del p[0]\n" + f"return {DIRECT}\n",
        HEAD + "alias = p\ndel alias[0]\n" + f"return {DIRECT}\n",
    ],
    ids=["walrus-size", "delete-resizes-vector", "delete-resizes-vector-alias"],
)
def test_extent_invalidated_by_rebinding_or_resize(body: str) -> None:
    errors = registry_errors({"src/x.py": source(body)}, by_registry())
    assert any(DERIVATION_ERROR in error for error in errors), errors


@pytest.mark.parametrize(
    ("binding", "body"),
    [
        (
            "sum = lambda values: 2.0",
            "m = len(p)\nh = sum(1.0 / i for i in range(1, m + 1))\n" + f"return {DIRECT}\n",
        ),
        ("len = lambda values: 2", HEAD + f"return {DIRECT}\n"),
        (
            "range = lambda *args: [1, 2]",
            "m = len(p)\nh = sum(1.0 / i for i in range(1, m + 1))\n" + f"return {DIRECT}\n",
        ),
        (
            "float = lambda value: 2.0",
            "m = len(p)\nh = float(np.sum(1.0 / np.arange(1, m + 1)))\n" + f"return {DIRECT}\n",
        ),
        (
            "class sum:\n    pass",
            "m = len(p)\nh = sum(1.0 / i for i in range(1, m + 1))\n" + f"return {DIRECT}\n",
        ),
    ],
    ids=["sum", "len", "range", "float", "class-sum"],
)
def test_module_shadowed_builtins_are_not_trusted(binding: str, body: str) -> None:
    text = source(body).replace("def screen(p, other):", binding + "\ndef screen(p, other):")
    errors = registry_errors({"src/x.py": text}, by_registry())
    assert any(DERIVATION_ERROR in error for error in errors), errors


@pytest.mark.parametrize(
    "body",
    [
        HEAD + BH + f"if False:\n    by = {RESULT}\nreturn q\n",
        HEAD + BH + f"while False:\n    by = {RESULT}\nreturn q\n",
        HEAD + BH + f"if 0:\n    by = {RESULT}\nreturn q\n",
        HEAD + BH + f"while None:\n    by = {RESULT}\nreturn q\n",
        HEAD + BH + f"return q\nby = {RESULT}\n",
        HEAD + BH + f"raise RuntimeError('stop')\nby = {RESULT}\n",
        HEAD + BH + f"{RESULT}\nreturn q\n",
        HEAD + BH + f"for item in other:\n    break\n    by = {RESULT}\nreturn q\n",
        HEAD + BH + f"for item in other:\n    continue\n    by = {RESULT}\nreturn q\n",
        HEAD + BH + f"if True:\n    return q\nelse:\n    by = {RESULT}\nreturn q\n",
        HEAD + BH + f"if other:\n    return q\nelse:\n    raise ValueError()\nby = {RESULT}\n",
        HEAD + BH + f"with manager():\n    return q\nby = {RESULT}\n",
        HEAD + BH + "with outer_manager():\n    with inner_manager():\n        return q\n"
        f"    by = {RESULT}\n",
        HEAD + BH + "if other:\n    with manager():\n        return q\n"
        f"else:\n    return q\nby = {RESULT}\n",
    ],
    ids=[
        "if-false",
        "while-false",
        "if-zero",
        "while-none",
        "after-return",
        "after-raise",
        "discarded-expression",
        "after-break",
        "after-continue",
        "constant-true-else",
        "both-branches-terminate",
        "after-with-return",
        "after-nested-with-return",
        "after-if-with-return",
    ],
)
def test_unreachable_or_discarded_calculation_does_not_establish_by(body: str) -> None:
    errors = registry_errors({"src/x.py": source(body)}, by_registry())
    assert any(DERIVATION_ERROR in error for error in errors), errors


@pytest.mark.parametrize(
    "body",
    [
        HEAD + "for chunk in other:\n"
        "    p_r = np.ones(m, dtype=float)\n    p_r[0] = 0.5\n"
        "    q = benjamini_hochberg(p_r)\n"
        f"    by = {RESULT}\nreturn p\n",
        HEAD + "while other:\n"
        "    p_r = np.ones(m, dtype=float)\n    p_r[0] = 0.5\n"
        "    q = benjamini_hochberg(p_r)\n"
        f"    by = {RESULT}\nreturn p\n",
        HEAD + "with manager():\n" + indent(BH + f"return {RESULT}\n", "    "),
        HEAD + BH + f"if False:\n    q = p\nelse:\n    return {RESULT}\n",
    ],
    ids=["for-rebuilt-vector", "while-rebuilt-vector", "with-body", "reachable-else"],
)
def test_live_supported_derivation_is_accepted(body: str) -> None:
    assert not registry_errors({"src/x.py": source(body)}, by_registry())
