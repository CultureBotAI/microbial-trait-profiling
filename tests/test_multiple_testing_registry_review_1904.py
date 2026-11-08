"""Targeted alias and cache regressions from PR #1904 review (#1911-#1915)."""

from __future__ import annotations

from pathlib import Path
from textwrap import indent
from typing import Literal

import pytest

import test_multiple_testing_registry as guard
import test_multiple_testing_registry_alias_forms as aliases

SCOPES: tuple[aliases.Scope, ...] = ("parameter", "module")
FreeScope = Literal["module", "closure"]
SLOTS = {
    "list-slot": ("box = [p]\n", "box[0]"),
    "dict-entry": ('holder = {"v": p}\n', 'holder["v"]'),
    "namespace-attribute": ("holder = types.SimpleNamespace(v=p)\n", "holder.v"),
}
LOOP_JOINS = {
    "while-not-alias": "alias = []\nwhile not alias:\n    alias = p\n",
    "two-pass-for": "alias = []\nfor _pass in (0, 1):\n    alias = p\n",
}


@pytest.fixture(autouse=True)
def assigned_guard() -> None:
    assert guard.__file__ is not None
    assert Path(guard.__file__).resolve() == (
        Path(__file__).with_name("test_multiple_testing_registry.py").resolve()
    )
    assert aliases.__file__ is not None
    assert Path(aliases.__file__).resolve() == (
        Path(__file__).with_name("test_multiple_testing_registry_alias_forms.py").resolve()
    )


def parameter_write_source(scope: FreeScope, operation: str) -> tuple[str, str]:
    """The caller may pass the free BH vector as the mutation parameter (#1911)."""
    function = "def screen(out):\n" + indent(aliases.EXTENT + operation + aliases.RESULT, "    ")
    if scope == "module":
        return (
            aliases.IMPORTS + f"p = {aliases.VECTOR}\n" + function + "screen(p)\n",
            "screen",
        )
    return (
        aliases.IMPORTS + "def outer(p):\n" + indent(function + "return screen(p)\n", "    "),
        "outer.screen",
    )


@pytest.mark.parametrize("scope", ["module", "closure"])
def test_parameter_write_invalidates_observed_free_vector_extent(scope: FreeScope) -> None:
    """A caller-bound parameter can shorten the free vector after its harmonic factor."""
    text, _ = parameter_write_source(scope, "del out[0]\n")
    assert aliases.derived(text) == ()


@pytest.mark.parametrize("scope", ["module", "closure"])
def test_parameter_read_keeps_observed_free_vector_extent(scope: FreeScope) -> None:
    text, anchor = parameter_write_source(scope, "first = out[0]\n")
    assert aliases.derived(text) == (anchor,)


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize("slot", SLOTS)
def test_augmented_slot_write_invalidates_the_referent_extent(
    scope: aliases.Scope, slot: str
) -> None:
    """For list vectors, += mutates the object before storing it back (#1912)."""
    setup, target = SLOTS[slot]
    text, _ = aliases.source(
        scope, aliases.EXTENT + setup + f"{target} += [0.5]\n" + aliases.RESULT
    )
    assert aliases.derived(text) == ()


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize("slot", SLOTS)
def test_plain_slot_replacement_keeps_the_original_vector_extent(
    scope: aliases.Scope, slot: str
) -> None:
    """Replacing the holder's slot is not an in-place write to its old vector."""
    setup, target = SLOTS[slot]
    text, anchor = aliases.source(
        scope, aliases.EXTENT + setup + f"{target} = [0.5]\n" + aliases.RESULT
    )
    assert aliases.derived(text) == (anchor,)


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize("join", LOOP_JOINS)
def test_repeated_loop_join_retains_the_possible_alias_until_its_write(
    scope: aliases.Scope, join: str
) -> None:
    """A merge-cache hit must retain the union, not return its left input (#1913)."""
    text, _ = aliases.source(
        scope, aliases.EXTENT + LOOP_JOINS[join] + "del alias[0]\n" + aliases.RESULT
    )
    assert aliases.derived(text) == ()


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize("join", LOOP_JOINS)
def test_repeated_loop_join_without_a_write_keeps_the_extent(
    scope: aliases.Scope, join: str
) -> None:
    text, anchor = aliases.source(scope, aliases.EXTENT + LOOP_JOINS[join] + aliases.RESULT)
    assert aliases.derived(text) == (anchor,)


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize(
    "copy_expression", ["copy.copy(p)", "sorted(p)"], ids=["copy-copy", "sorted"]
)
def test_independent_copy_delete_does_not_invalidate_original_extent(
    scope: aliases.Scope, copy_expression: str
) -> None:
    """Both standard-library copies preserve the original vector's length (#1914)."""
    text, anchor = aliases.source(
        scope,
        aliases.EXTENT + f"alias = {copy_expression}\ndel alias[0]\n" + aliases.RESULT,
    )
    assert aliases.derived("import copy\n" + text) == (anchor,)


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize("attribute", ["v", "w"], ids=["vector-attribute", "other-attribute"])
def test_namespace_mutation_respects_the_selected_attribute(
    scope: aliases.Scope, attribute: str
) -> None:
    """Selecting w must not inherit v's possible reference (#1915)."""
    text, anchor = aliases.source(
        scope,
        aliases.EXTENT
        + "holder = types.SimpleNamespace(v=p, w=[0.5, 0.6])\n"
        + f"del holder.{attribute}[0]\n"
        + aliases.RESULT,
    )
    assert aliases.derived(text) == (() if attribute == "v" else (anchor,))


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize("write", ["", "allocated[0] = 0.25\n"], ids=["untouched", "index-write"])
def test_zeros_scalar_shape_preserves_the_matching_harmonic_extent(
    scope: aliases.Scope, write: str
) -> None:
    result = aliases.RESULT.replace("benjamini_hochberg(p)", "benjamini_hochberg(allocated)")
    text, anchor = aliases.source(
        scope, aliases.EXTENT + "allocated = np.zeros(m, dtype=float)\n" + write + result
    )
    assert aliases.derived(text) == (anchor,)


@pytest.mark.parametrize("allocator", ["ones", "zeros"])
@pytest.mark.parametrize("case", ["other-extent", "deleted-element"])
def test_fixed_allocation_requires_matching_and_unmodified_extent(
    allocator: str, case: str
) -> None:
    shape = "len(other)" if case == "other-extent" else "m"
    write = "del allocated[0]\n" if case == "deleted-element" else ""
    result = aliases.RESULT.replace("benjamini_hochberg(p)", "benjamini_hochberg(allocated)")
    body = aliases.EXTENT + f"allocated = np.{allocator}({shape}, dtype=float)\n" + write + result
    text = aliases.IMPORTS + "def screen(p, other):\n" + indent(body, "    ")
    # AST-only: deletion must invalidate the allocated extent even though an
    # actual ndarray rejects deletion. A different family length also cannot match.
    assert aliases.derived(text) == ()


@pytest.mark.parametrize("scope", SCOPES)
@pytest.mark.parametrize("shape", ["p", "(reps, m)"], ids=["shape-vector", "matrix-shape"])
def test_zero_allocation_does_not_retain_shape_or_dtype_parameter_references(
    scope: aliases.Scope, shape: str
) -> None:
    # [1, 1] is a valid NumPy shape and also lies in the nominal-p range. The
    # separate dtype argument denotes a dtype, never the BH vector itself.
    prefix = "p = [1, 1]\n" if scope == "module" else ""
    parameters = "reps, dtype" if scope == "module" else "p, reps, dtype"
    body = (
        aliases.EXTENT
        + f"allocated = np.zeros({shape}, dtype=dtype)\n"
        + "allocated[0] = 0.25\n"
        + aliases.RESULT
    )
    text = aliases.IMPORTS + prefix + f"def screen({parameters}):\n" + indent(body, "    ")
    assert aliases.derived(text) == ("screen",)


@pytest.mark.parametrize("allocator", ["ones", "zeros"])
@pytest.mark.parametrize(
    "arguments", ["m, dtype=float, like=p", "p, float"], ids=["like-dispatch", "two-positional"]
)
def test_allocation_outside_the_narrow_signature_remains_an_opaque_alias(
    allocator: str, arguments: str
) -> None:
    # like= can dispatch through p.__array_function__ and return p. The
    # optimization deliberately recognizes exactly one shape argument and dtype
    # keywords only; two positional arguments stay outside that narrow rule.
    setup = aliases.EXTENT + f"allocated = np.{allocator}({arguments})\n"
    text, anchor = aliases.source("parameter", setup + aliases.RESULT)
    assert aliases.derived(text) == (anchor,)
    changed, _ = aliases.source("parameter", setup + "del allocated[0]\n" + aliases.RESULT)
    assert aliases.derived(changed) == ()


@pytest.mark.parametrize("allocator", ["ones", "zeros"])
@pytest.mark.parametrize("shadow", ["callable-name", "namespace-parameter", "namespace-local"])
def test_shadowed_allocation_call_cannot_claim_fresh_reference_credit(
    allocator: str, shadow: str
) -> None:
    # Keep the harmonic factor and clipping bound to real NumPy. Otherwise
    # shadowing np could remove the BY witness before the allocation is examined.
    imports = aliases.IMPORTS + "import numpy as real_np\n"
    parameters = "p"
    if shadow == "callable-name":
        imports += f"from numpy import {allocator}\n"
        setup = f"{allocator} = lambda values: values\n"
        call = f"{allocator}(p)"
    else:
        parameters = "p, np" if shadow == "namespace-parameter" else "p, holder"
        setup = "" if shadow == "namespace-parameter" else "np = holder\n"
        call = f"np.{allocator}(p)"
    prefix = setup + aliases.EXTENT.replace("np.", "real_np.") + f"allocated = {call}\n"
    result = aliases.RESULT.replace("np.", "real_np.")
    for write, expected in (("", ("screen",)), ("del allocated[0]\n", ())):
        text = imports + f"def screen({parameters}):\n" + indent(prefix + write + result, "    ")
        assert aliases.derived(text) == expected


@pytest.mark.parametrize("allocator", ["ones", "zeros"])
@pytest.mark.parametrize(
    ("dtype_setup", "dtype"),
    [("", "object"), ("dtype = object\n", "dtype"), ("float = object\n", "float")],
)
@pytest.mark.parametrize("through_alias", [False, True])
def test_object_array_keeps_references_inserted_after_allocation(
    allocator: str, dtype_setup: str, dtype: str, through_alias: bool
) -> None:
    setup = (
        aliases.EXTENT
        + dtype_setup
        + f"box = np.{allocator}(m, dtype={dtype})\n"
        + ("holder = box\nholder[0] = p\n" if through_alias else "box[0] = p\n")
    )
    control, anchor = aliases.source("parameter", setup + aliases.RESULT)
    changed, _ = aliases.source("parameter", setup + "del box[0][0]\n" + aliases.RESULT)
    assert aliases.derived(control) == (anchor,)
    assert aliases.derived(changed) == ()
