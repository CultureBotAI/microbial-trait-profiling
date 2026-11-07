"""Pin individual derivation charges with independently counted fixtures (#1686)."""

from __future__ import annotations

import ast
import inspect
from dataclasses import replace
from pathlib import Path

import pytest

import test_multiple_testing_registry as guard
from test_multiple_testing_registry_derived import PREFIX, source

ONE_RETURN = "def screen(p):\n    return p\n"
RETURN_CHARGES = {"walk_node": 5, "walk_list": 1, "statement": 1, "expression": 1}
TWO_FACT_BRANCH = "def screen(p, q):\n    if flag:\n        pass\n"
CHARGE_CATEGORIES = {
    "env_copy",
    "walk_node",
    "walk_list",
    "expression",
    "assignment",
    "write_alias_scan",
    "alias_invalidation",
    "env_join",
    "import_names",
    "declared_names",
    "terminates",
    "statement",
    "loop_pass",
}


@pytest.fixture(autouse=True)
def assigned_guard() -> None:
    assert guard.__file__ is not None
    assert (
        Path(guard.__file__).resolve()
        == (Path(__file__).parent / "test_multiple_testing_registry.py").resolve()
    )


def parsed_unit(text: str, limit: int = guard.DEFAULT_LIMITS.derivation_steps) -> guard.Unit:
    unit = guard.Unit(
        "src/counts.py",
        "counts",
        ast.parse(text, filename="src/counts.py"),
        replace(guard.DEFAULT_LIMITS, derivation_steps=limit),
    )
    guard.Visitor(unit).visit(unit.tree)
    return unit


def observed(unit: guard.Unit, counts: dict[str, int]) -> bool:
    parameter = inspect.signature(guard.visible_by_derivation).parameters.get("work_counts")
    assert parameter is not None, "The bounded optional accounting hook must exist"
    assert parameter.kind == inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is None
    result = guard.visible_by_derivation(unit, "screen", work_counts=counts)
    assert counts.keys() <= CHARGE_CATEGORIES
    return result


@pytest.mark.parametrize(
    ("text", "category", "expected"),
    [
        # Both branches copy the two parameter facts: 2 + 2.
        pytest.param(TWO_FACT_BRANCH, "env_copy", 4, id="environment-copy"),
        # Both branches now copy three parameter facts: 3 + 3.
        pytest.param(
            "def screen(p, q, r):\n    if flag:\n        pass\n",
            "env_copy",
            6,
            id="environment-copy-three-facts",
        ),
        # The join reads both two-fact inputs: 2 + 2.
        pytest.param(TWO_FACT_BRANCH, "env_join", 4, id="environment-join"),
        # The join now reads two three-fact inputs: 3 + 3.
        pytest.param(
            "def screen(p, q, r):\n    if flag:\n        pass\n",
            "env_join",
            6,
            id="environment-join-three-facts",
        ),
        # A pass changes no fact, so the first fixed-point pass is stable.
        pytest.param(
            "def screen(p):\n    while flag:\n        pass\n",
            "loop_pass",
            1,
            id="loop-pass",
        ),
        # The first pass changes p; a second pass establishes the fixed point.
        pytest.param(
            "def screen(p):\n    while flag:\n        p = 1\n",
            "loop_pass",
            2,
            id="loop-pass-two-iterations",
        ),
        # Exactly three statements, including the bare return.
        pytest.param(
            "def screen():\n    pass\n    pass\n    return\n",
            "statement",
            3,
            id="statement",
        ),
        # The multiplication and its two operands are three expressions.
        pytest.param(
            "def screen(p, other):\n    return p * other\n",
            "expression",
            3,
            id="expression",
        ),
        # The tuple target and its two names require three assignments.
        pytest.param("def screen(p):\n    a, b = p\n", "assignment", 3, id="assignment"),
        # Two parameters plus alias make the write scan three facts.
        pytest.param(
            "def screen(p, q):\n    alias = p\n    alias[0] = 1\n",
            "write_alias_scan",
            3,
            id="write-alias-scan",
        ),
        # Three parameters plus alias make the write scan four facts.
        pytest.param(
            "def screen(p, q, r):\n    alias = p\n    alias[0] = 1\n",
            "write_alias_scan",
            4,
            id="write-alias-scan-four-facts",
        ),
        # q and alias plus the remembered free origin p are three scan facts.
        pytest.param(
            "def screen(q):\n    alias = p\n    alias[0] = 1\n",
            "write_alias_scan",
            3,
            id="write-alias-scan-free-origin",
        ),
        # The same three facts must be inspected after augmented assignment.
        pytest.param(
            "def screen(p, q):\n    alias = p\n    alias += []\n",
            "alias_invalidation",
            3,
            id="alias-invalidation",
        ),
        # Three parameters and alias require four invalidation inspections.
        pytest.param(
            "def screen(p, q, r):\n    alias = p\n    alias += []\n",
            "alias_invalidation",
            4,
            id="alias-invalidation-four-facts",
        ),
        # The remembered free origin p joins q and alias before this scan.
        pytest.param(
            "def screen(q):\n    alias = p\n    alias += []\n",
            "alias_invalidation",
            3,
            id="alias-invalidation-free-origin",
        ),
        # Both branch statements are inspected to decide whether it terminates.
        pytest.param(
            "def screen():\n    if True:\n        pass\n        pass\n",
            "terminates",
            2,
            id="termination-scan",
        ),
        # Three declared entries cost three even though a is repeated.
        pytest.param(
            "def screen():\n    global a, a, b\n",
            "declared_names",
            3,
            id="declared-names",
        ),
        # The import names are charged independently of binding their names.
        pytest.param(
            "def screen():\n    import a as x, b as y\n",
            "import_names",
            2,
            id="import-names",
        ),
        # Three imported names each incur an import and a binding charge.
        pytest.param(
            "def screen():\n    import a as x, b as y, c as z\n",
            "import_names",
            3,
            id="import-names-three-entries",
        ),
    ],
)
def test_each_overlapping_charge_remains_observable(
    text: str, category: str, expected: int
) -> None:
    counts: dict[str, int] = {}
    assert observed(parsed_unit(text), counts) is False
    assert counts.get(category, 0) == expected
    if category == "import_names":
        # Shared binding accounting includes imports as well as declarations.
        assert counts.get("declared_names", 0) == expected


def test_one_return_has_eight_independently_counted_units() -> None:
    # The module walk visits Module and FunctionDef, plus its one body entry.
    # The return walk visits Return, Name and Load. Then one statement and
    # one expression complete the eight-unit analysis.
    counts: dict[str, int] = {}
    assert observed(parsed_unit(ONE_RETURN), counts) is False
    assert {key: value for key, value in counts.items() if value} == RETURN_CHARGES
    assert sum(counts.values()) == 8


@pytest.mark.parametrize("hooked", [False, True], ids=["unhooked", "hooked"])
@pytest.mark.parametrize(
    "limit", [0, 7, 8], ids=["refuses-at-zero", "refuses-at-seven", "accepts-at-eight"]
)
def test_accounting_hook_preserves_the_exact_budget_boundary(hooked: bool, limit: int) -> None:
    unit = parsed_unit(ONE_RETURN, limit)
    counts: dict[str, int] = {}

    def invoke() -> bool:
        if hooked:
            return observed(unit, counts)
        return guard.visible_by_derivation(unit, "screen")

    if limit < 8:
        attempted = 1 if limit == 0 else 8
        with pytest.raises(
            guard.DiscoveryLimitError,
            match=rf"src/counts[.]py:screen \(attempted {attempted}, limit {limit}\)",
        ):
            invoke()
    else:
        assert invoke() is False
    if hooked:
        # The rejected first node or final expression is still attempted work.
        expected = {"walk_node": 1} if limit == 0 else RETURN_CHARGES
        assert {key: value for key, value in counts.items() if value} == expected
        assert counts.keys() <= CHARGE_CATEGORIES


def test_hook_preserves_a_normal_by_witness_and_discovery_output() -> None:
    text = source(PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n")
    before = guard.discover({"src/counts.py": text})
    assert before[0][("src/counts.py", "screen")]["derived_procedures"] == {"BY"}
    counts: dict[str, int] = {}
    assert observed(parsed_unit(text), counts) is True
    assert guard.visible_by_derivation(parsed_unit(text), "screen") is True
    assert sum(counts.values()) > 8
    assert guard.discover({"src/counts.py": text}) == before


def test_fresh_counters_and_unhooked_calls_retain_no_previous_state() -> None:
    first: dict[str, int] = {}
    assert observed(parsed_unit(ONE_RETURN), first) is False
    snapshot = first.copy()
    second: dict[str, int] = {}
    assert observed(parsed_unit(TWO_FACT_BRANCH), second) is False
    assert second.get("env_copy") == 4
    assert second != first
    assert guard.visible_by_derivation(parsed_unit(TWO_FACT_BRANCH), "screen") is False
    assert first == snapshot
    third: dict[str, int] = {}
    assert observed(parsed_unit(ONE_RETURN), third) is False
    assert {key: value for key, value in third.items() if value} == RETURN_CHARGES
    assert third == first


def test_reused_counter_accumulates_without_sharing_the_per_invocation_limit() -> None:
    counts: dict[str, int] = {}
    assert observed(parsed_unit(ONE_RETURN, limit=8), counts) is False
    try:
        assert observed(parsed_unit(ONE_RETURN, limit=8), counts) is False
    except guard.DiscoveryLimitError as error:
        raise AssertionError("a reused report must not share the budget") from error
    # Each invocation uses eight units; the caller's report accumulates sixteen.
    assert {key: value for key, value in counts.items() if value} == {
        "walk_node": 10,
        "walk_list": 2,
        "statement": 2,
        "expression": 2,
    }
