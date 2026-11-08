"""Compound and external aliases cannot preserve a stale BY extent (#1796, #1798)."""

from __future__ import annotations

import ast
import importlib.util
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
from textwrap import indent
from typing import Literal

import pytest

import test_multiple_testing_registry as guard

Scope = Literal["parameter", "module", "closure"]
SCOPES = ("parameter", "module", "closure")
IMPORTS = (
    "import contextlib\nimport types\nfrom typing import cast\nimport numpy as np\n"
    "from mtp.enrichment_stats import benjamini_hochberg\n"
)
EXTENT = "m = len(p)\nh = np.sum(1.0 / np.arange(1, m + 1))\n"
RESULT = "q = benjamini_hochberg(p)\nreturn np.clip(q * h, 0, 1)\n"
VECTOR = "[0.01, 0.02, 0.03, 0.2]"

# Each pair contains a real write and the same binding/control flow without a
# write. The latter must retain credit: merely encountering an alias is safe.
FORMS = {
    "tuple-unpack": (
        "alias, _unused = p, None\ndel alias[0]\n",
        "alias, _unused = p, None\n",
    ),
    "named-container-unpack": (
        "pair = (p, None)\nalias, _unused = pair\ndel alias[0]\n",
        "pair = (p, None)\nalias, _unused = pair\n",
    ),
    "nested-unpack": (
        "(alias, _inner), _outer = (p, None), None\ndel alias[0]\n",
        "(alias, _inner), _outer = (p, None), None\n",
    ),
    "boolop": ("alias = p or []\ndel alias[0]\n", "alias = p or []\n"),
    "conditional": ("alias = p if m else p\ndel alias[0]\n", "alias = p if m else p\n"),
    "walrus": ("(alias := p)\ndel alias[0]\n", "(alias := p)\n"),
    "branch-join": (
        "if m:\n    alias = p\nelse:\n    alias = []\ndel alias[0]\n",
        "if m:\n    alias = p\nelse:\n    alias = []\n",
    ),
    "branch-no-else": (
        "alias = []\nif m:\n    alias = p\ndel alias[0]\n",
        "alias = []\nif m:\n    alias = p\n",
    ),
    "loop-join": (
        "for _item in (0,):\n    alias = p\ndel alias[0]\n",
        "for _item in (0,):\n    alias = p\n",
    ),
    "for-target": (
        "for alias in (p,):\n    del alias[0]\n",
        "for alias in (p,):\n    pass\n",
    ),
    "later-iteration": (
        "alias = []\nfor item in (0, 1):\n"
        "    if item:\n        del alias[0]\n    else:\n        alias = p\n",
        "alias = []\nfor item in (0, 1):\n"
        "    if item:\n        pass\n    else:\n        alias = p\n",
    ),
    "with-target": (
        "with contextlib.nullcontext(p) as alias:\n    del alias[0]\n",
        "with contextlib.nullcontext(p) as alias:\n    pass\n",
    ),
    "try-alias": (
        "try:\n    alias = p\nexcept Exception:\n    alias = []\ndel alias[0]\n",
        "try:\n    alias = p\nexcept Exception:\n    alias = []\n",
    ),
    "container": ("box = [p]\ndel box[0][0]\n", "box = [p]\n"),
    "container-insert": (
        "box = [None]\nbox[0] = p\ndel box[0][0]\n",
        "box = [None]\nbox[0] = p\n",
    ),
    "namespace": (
        "holder = types.SimpleNamespace(v=p)\ndel holder.v[0]\n",
        "holder = types.SimpleNamespace(v=p)\n",
    ),
    "namespace-insert": (
        "holder = types.SimpleNamespace()\nholder.v = p\ndel holder.v[0]\n",
        "holder = types.SimpleNamespace()\nholder.v = p\n",
    ),
    "cast": ("alias = cast(list, p)\ndel alias[0]\n", "alias = cast(list, p)\n"),
}
SAFE_WRITES = {
    "slice-copy": "alias = p[:]\ndel alias[0]\n",
    "list-copy": "alias = list(p)\ndel alias[0]\n",
    "method-copy": "alias = p.copy()\ndel alias[0]\n",
    "fresh-list": "alias = [0.5, 0.6]\ndel alias[0]\n",
    "local-named-external": "external = [0.5, 0.6]\ndel external[0]\n",
    "fresh-unpacked-list": "alias, _unused = [0.5, 0.6], None\ndel alias[0]\n",
    "sized-array": "arr = np.ones(m)\ndel arr[0]\n",
    "alias-rebound": "alias = p if m else p\nalias = [0.5, 0.6]\ndel alias[0]\n",
    "container-rebound": "box = [p]\nbox = [[0.5, 0.6]]\ndel box[0][0]\n",
    "outer-container-replaced": "box = [p]\nbox[0] = []\n",
    "outer-container-deleted": "box = [p]\ndel box[0]\n",
    "alias-name-deleted": "alias = p if m else p\ndel alias\n",
    "saved-size": "m += 1\n",
}


@pytest.fixture(autouse=True)
def assigned_guard() -> None:
    assert guard.__file__ is not None
    assert Path(guard.__file__).resolve() == (
        Path(__file__).with_name("test_multiple_testing_registry.py").resolve()
    )


def source(scope: Scope, body: str, *, external: bool = False) -> tuple[str, str]:
    if scope == "parameter":
        prefix = f"p = {VECTOR}\nalias = p\n" if external else ""
        return IMPORTS + prefix + "def screen(p):\n" + indent(body, "    "), "screen"
    function = "def screen():\n" + indent(body, "    ")
    prefix = "alias = p\n" if external else ""
    if scope == "module":
        return IMPORTS + f"p = {VECTOR}\n" + prefix + function, "screen"
    return (
        IMPORTS + "def outer(p):\n" + indent(prefix + function + "return screen\n", "    "),
        "outer.screen",
    )


def derived(text: str) -> tuple[str, ...]:
    return guard.summarize_source(
        "src/alias_forms.py", text, guard.SummaryBudget(guard.DEFAULT_LIMITS)
    ).derived


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize("form", FORMS)
def test_untracked_alias_write_invalidates_the_original_extent(scope: Scope, form: str) -> None:
    text, _ = source(scope, EXTENT + FORMS[form][0] + RESULT)
    assert derived(text) == ()


@pytest.mark.parametrize("scope", SCOPES)
def test_loop_aliases_converge_after_arithmetic_facts_stabilize(scope: Scope) -> None:
    # Each iteration advances p by one link. The arithmetic facts become
    # unknown earlier, so their fixed point alone cannot authorize the result.
    chain = (
        "a = []\nb = []\nc = []\nd = []\n"
        "for item in range(4):\n"
        "    d = c if m else c\n"
        "    c = b if m else b\n"
        "    b = a if m else a\n"
        "    a = p if m else p\n"
        "del d[0]\n"
    )
    text, _ = source(scope, EXTENT + chain + RESULT)
    assert derived(text) == ()


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize("form", FORMS)
def test_read_only_compound_alias_preserves_the_witness(scope: Scope, form: str) -> None:
    text, anchor = source(scope, EXTENT + FORMS[form][1] + RESULT)
    assert derived(text) == (anchor,)


@pytest.mark.parametrize("scope", SCOPES)
def test_alias_bound_outside_the_function_invalidates_the_original(scope: Scope) -> None:
    # In the parameter case the caller may pass the module's alias as p. The
    # numerical witness below executes exactly that legal call.
    text, _ = source(scope, EXTENT + "del alias[0]\n" + RESULT, external=True)
    assert derived(text) == ()


@pytest.mark.parametrize("scope", SCOPES)
def test_insertion_updates_an_already_nested_holder_reference(scope: Scope) -> None:
    body = (
        "inner = [None]\nouter = [inner]\nother = inner\nother[0] = p\n"
        + EXTENT
        + "del outer[0][0][0]\n"
        + RESULT
    )
    text, _ = source(scope, body)
    assert derived(text) == ()


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize("nested", [False, True], ids=["direct", "nested"])
def test_holder_identity_survives_a_prior_non_alias_write(scope: Scope, nested: bool) -> None:
    setup = "inner = [None]\nother = inner\n"
    if nested:
        setup += "outer = [inner]\n"
    access = "outer[0][0][0]" if nested else "inner[0][0]"
    body = setup + "inner[0] = None\nother[0] = p\n" + EXTENT + f"del {access}\n" + RESULT
    text, _ = source(scope, body)
    assert derived(text) == ()


@pytest.mark.parametrize("scope", SCOPES)
def test_reading_an_external_alias_length_does_not_prove_it_distinct(scope: Scope) -> None:
    body = EXTENT + "spare = len(alias)\ndel alias[0]\n" + RESULT
    text, _ = source(scope, body, external=True)
    assert derived(text) == ()


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize(
    ("setup", "write"),
    [
        ("box = [None]\nother = box\nother[0] = p\n", "del box[0][0]\n"),
        (
            "holder = types.SimpleNamespace()\nother = holder\nother.v = p\n",
            "del holder.v[0]\n",
        ),
    ],
    ids=["list-holder", "namespace-holder"],
)
def test_insertion_through_a_holder_alias_is_visible_through_the_original(
    scope: Scope, setup: str, write: str
) -> None:
    text, _ = source(scope, setup + EXTENT + write + RESULT)
    assert derived(text) == ()


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize(
    "binding",
    [
        "local_alias = alias\n",
        "local_alias = cast(list, alias)\n",
        "local_alias = alias if alias is not None else alias\n",
    ],
    ids=["name", "cast", "conditional"],
)
def test_local_copy_of_an_external_alias_keeps_its_mutation_risk(
    scope: Scope, binding: str
) -> None:
    body = binding + EXTENT + "del local_alias[0]\n" + RESULT
    text, _ = source(scope, body, external=True)
    assert derived(text) == ()


@pytest.mark.parametrize("scope", SCOPES)
def test_repeated_walrus_binding_preserves_the_alias_mutation_risk(scope: Scope) -> None:
    body = "alias = p\n" + EXTENT + "(alias := p)\ndel alias[0]\n" + RESULT
    text, _ = source(scope, body)
    assert derived(text) == ()


@pytest.mark.parametrize("scope", SCOPES)
def test_repeated_try_binding_preserves_the_alias_mutation_risk(scope: Scope) -> None:
    body = (
        "alias = p\n"
        + EXTENT
        + "try:\n    alias = p\nexcept Exception:\n    pass\ndel alias[0]\n"
        + RESULT
    )
    text, _ = source(scope, body)
    assert derived(text) == ()


def test_caller_alias_does_not_depend_on_the_parameter_name() -> None:
    body = (EXTENT + "del alias[0]\n" + RESULT).replace("len(p)", "len(values)")
    body = body.replace("benjamini_hochberg(p)", "benjamini_hochberg(values)")
    text = IMPORTS + f"p = {VECTOR}\nalias = p\ndef screen(values):\n" + indent(body, "    ")
    # Calling screen(alias) is legal; its parameter need not share the name of
    # either external binding. The write still changes the BH input's extent.
    assert derived(text) == ()


@pytest.mark.parametrize("scope", SCOPES)
def test_fresh_local_alias_shadows_an_external_binding_without_mutating_it(scope: Scope) -> None:
    body = EXTENT + "alias = [0.5, 0.6]\ndel alias[0]\n" + RESULT
    text, anchor = source(scope, body, external=True)
    assert derived(text) == (anchor,)


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize("control", SAFE_WRITES)
def test_copy_allocation_and_rebinding_writes_preserve_the_witness(
    scope: Scope, control: str
) -> None:
    text, anchor = source(scope, EXTENT + SAFE_WRITES[control] + RESULT)
    assert derived(text) == (anchor,)


@pytest.mark.parametrize(
    ("binding", "write"),
    [
        ("alias, _unused = p, None\n", "del alias[0]\n"),
        ("box = [p]\n", "del box[0][0]\n"),
    ],
    ids=["unpacked-alias", "container"],
)
def test_old_compound_alias_cannot_invalidate_a_rebound_parameter(binding: str, write: str) -> None:
    body = binding + "p = other\n" + EXTENT + write + RESULT
    text = IMPORTS + "def screen(p, other):\n" + indent(body, "    ")
    assert derived(text) == ("screen",)


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize("form", ["tuple-unpack", "container", "external"])
def test_real_series_alias_witness_is_numerically_wrong_and_not_credited(
    tmp_path: Path, scope: Scope, form: str
) -> None:
    from mtp.cellulose_decr import analysis

    mutation = "del alias[0]\n" if form == "external" else FORMS[form][0]
    text, _ = source(scope, EXTENT + mutation + RESULT, external=form == "external")
    text = text.replace(
        "from mtp.enrichment_stats import benjamini_hochberg",
        "import pandas as pd\nfrom mtp.cellulose_decr.analysis import benjamini_hochberg",
    ).replace(f"p = {VECTOR}", f"p = pd.Series({VECTOR})")
    fixture_path = tmp_path / "series_alias_forms.py"
    fixture_path.write_text(text, encoding="utf-8")
    spec = importlib.util.spec_from_file_location("series_alias_forms", fixture_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.benjamini_hochberg is analysis.benjamini_hochberg
    assert Path(analysis.__file__).resolve() == (
        Path(__file__).resolve().parents[1] / "src/mtp/cellulose_decr/analysis.py"
    )
    if scope == "module" or (scope == "parameter" and form == "external"):
        vector = module.p
    else:
        vector = module.pd.Series([0.01, 0.02, 0.03, 0.2])
    try:
        if scope == "module":
            result = module.screen()
        elif scope == "parameter":
            result = module.screen(vector)
        else:
            result = module.outer(vector)()
    except (TypeError, ValueError, IndexError, KeyError) as error:
        raise AssertionError("The real Series witness must execute to numeric values") from error
    assert list(vector.index) == [1, 2, 3]
    assert list(vector) == [0.02, 0.03, 0.2]
    assert list(result.index) == [1, 2, 3]

    # Independent BY on the remaining sorted family, with rational arithmetic.
    remaining = (Fraction(2, 100), Fraction(3, 100), Fraction(1, 5))
    count = len(remaining)
    harmonic = sum((Fraction(1, rank) for rank in range(1, count + 1)), Fraction(0))
    expected = [
        float(
            min(
                Fraction(1),
                *(count * harmonic * remaining[j] / (j + 1) for j in range(i, count)),
            )
        )
        for i in range(count)
    ]
    assert expected == pytest.approx([33 / 400, 33 / 400, 11 / 30], rel=1e-12, abs=1e-14)
    assert list(result) == pytest.approx([3 / 32, 3 / 32, 5 / 12], rel=1e-12, abs=1e-14)
    assert list(result) != pytest.approx(expected, rel=1e-12, abs=1e-14)
    assert fixture_path.read_text(encoding="utf-8") == text
    assert derived(text) == ()


def work_unit(parameters: list[str], value: str, limit: int) -> guard.Unit:
    text = f"def screen({', '.join(parameters)}):\n    holder = {value}\n"
    tree = ast.parse(text, filename="src/alias_work.py")
    # Both fixtures have fewer than 500 AST nodes, with no loops or joins.
    assert sum(1 for _ in ast.walk(tree)) < 500
    unit = guard.Unit(
        "src/alias_work.py",
        "alias_work",
        tree,
        replace(guard.DEFAULT_LIMITS, derivation_steps=limit),
    )
    guard.Visitor(unit).visit(unit.tree)
    return unit


@pytest.mark.parametrize("count", [2, 40], ids=["small-control", "wide-union"])
def test_union_cardinality_is_bounded_without_control_flow(count: int) -> None:
    parameters = [f"p{index}" for index in range(count)]
    unit = work_unit(parameters, "[" + ", ".join(parameters) + "]", 1000)
    if count == 40:
        # Cumulative unions inspect 1+...+40=820 origins, with at least two
        # identifier characters per origin, before final binding costs.
        with pytest.raises(
            guard.DiscoveryLimitError,
            match=r"^derivation work budget exceeded: src/alias_work\.py:screen .*limit 1000\)",
        ):
            guard.visible_by_derivation(unit, "screen")
    else:
        counts: dict[str, int] = {}
        try:
            result = guard.visible_by_derivation(unit, "screen", work_counts=counts)
        except guard.DiscoveryLimitError as error:
            raise AssertionError("The two-origin control must fit the fixed allowance") from error
        assert result is False
        # Two-character p0/p1: each read costs dispatch 1, named lookup 1,
        # singleton lookup/construction 1+6, merge visits 2, cache entry 7.
        # The first accumulator merge visits two inputs; the second adds a
        # lookup, three three-unit scans and an entry. Containment costs
        # lookup 1, six-unit scan, six-unit construction and five-unit entry.
        # The arithmetic/list dispatches each cost 1; holder:2:4 costs 15.
        assert counts["expression"] == (
            1 + 1 + 2 * (1 + 1 + 1 + 6 + 2 + 7) + 2 + (2 + 1 + 3 + 3 + 3 + 4) + (1 + 6 + 6 + 5) + 15
        )


@pytest.mark.parametrize("depth", [2, 40], ids=["small-control", "deep-selector"])
def test_nested_selector_paths_are_bounded_without_control_flow(depth: int) -> None:
    unit = work_unit(["p"], "[" * depth + "p" + "]" * depth, 3000)
    counts: dict[str, int] = {}
    if depth == 40:
        # The two traversals of paths of lengths 0,...,39 already inspect
        # 2*780 selector positions; each '[]' selector also has two bytes.
        with pytest.raises(
            guard.DiscoveryLimitError,
            match=r"^derivation work budget exceeded: src/alias_work\.py:screen .*limit 3000\)",
        ):
            guard.visible_by_derivation(unit, "screen", work_counts=counts)
    else:
        try:
            result = guard.visible_by_derivation(unit, "screen", work_counts=counts)
        except guard.DiscoveryLimitError as error:
            raise AssertionError("The depth-two control must fit the fixed allowance") from error
        assert result is False
        # For [[p]]: arithmetic check 1, list dispatches 2, p read 17
        # (dispatch 1 + named lookup 1 + singleton 6 + merge 2 + cache 7),
        # inner/outer merge visits 2 each, containment 11/14 (lookup, scan,
        # path construction and cache entry), holder:2:4 singleton 15.
        # This pins path work independently of the generous small allowance.
        assert counts["expression"] == 1 + 2 + 17 + 2 + 11 + 2 + 14 + 15
