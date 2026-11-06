"""Pin conservative declared-scope reads in the BY witness guard (#1609).

The resolver does not implement declaration-aware lexical lookup. An explicit
``global`` or ``nonlocal`` therefore invalidates the declared name, including
callables; it is not evidence that the runtime object is the imported primitive.
Undeclared free vectors retain the existing bounded name/extent model. These
regressions pin behavior already present on main, rather than introduce a fix.
"""

from __future__ import annotations

from textwrap import indent
from typing import Literal

import pytest

from test_multiple_testing_registry import DEFAULT_LIMITS, SummaryBudget, summarize_source

Scope = Literal["global", "nonlocal"]
IMPORTS = "import numpy as np\nfrom mtp.enrichment_stats import benjamini_hochberg\n"
EXTENT = "m = len(p)\nh = np.sum(1.0 / np.arange(1, m + 1))\n"
RESULT = "q = benjamini_hochberg(p)\nreturn np.clip(q * h, 0, 1)\n"


def vector_source(scope: Scope, body: str) -> tuple[str, str]:
    """Use valid module or enclosing-function bindings without importing them."""
    function = "def screen(other):\n" + indent(body, "    ")
    if scope == "global":
        return IMPORTS + "p = [0.1, 0.2]\n" + function, "screen"
    return (
        IMPORTS + "def outer(p):\n" + indent(function + "return screen\n", "    "),
        "outer.screen",
    )


def derived(source: str) -> tuple[str, ...]:
    return summarize_source("src/x.py", source, SummaryBudget(DEFAULT_LIMITS)).derived


@pytest.mark.parametrize("scope", ["global", "nonlocal"])
def test_undeclared_free_vector_reads_are_accepted(scope: Scope) -> None:
    source, anchor = vector_source(scope, EXTENT + RESULT)
    assert derived(source) == (anchor,)


@pytest.mark.parametrize("scope", ["global", "nonlocal"])
def test_explicit_scope_declarations_invalidate_vector_reads(scope: Scope) -> None:
    source, _ = vector_source(scope, f"{scope} p\n" + EXTENT + RESULT)
    assert derived(source) == ()


@pytest.mark.parametrize("scope", ["global", "nonlocal"])
def test_declaration_in_dead_branch_still_invalidates_scope_read(scope: Scope) -> None:
    # Python declarations affect the whole scope even when their branch is dead.
    source, _ = vector_source(scope, f"if False:\n    {scope} p\n" + EXTENT + RESULT)
    assert derived(source) == ()


@pytest.mark.parametrize("scope", ["global", "nonlocal"])
def test_declared_vector_rebound_before_whole_calculation_is_supported(scope: Scope) -> None:
    # Invalidation is a fact reset, not a permanent ban on the enclosing anchor.
    source, anchor = vector_source(scope, f"{scope} p\np = other\n" + EXTENT + RESULT)
    assert derived(source) == (anchor,)


@pytest.mark.parametrize("callable_name", ["np", "benjamini_hochberg"])
def test_nonlocal_callable_parameter_is_not_trusted_as_import(callable_name: str) -> None:
    # Unit.resolve sees the module import, but runtime lookup reaches this
    # arbitrary outer parameter. Removing declaration invalidation is unsafe.
    function = (
        "def screen(p):\n"
        + indent(f"nonlocal {callable_name}\n" + EXTENT + RESULT, "    ")
        + "return screen\n"
    )
    source = IMPORTS + f"def outer({callable_name}):\n" + indent(function, "    ")
    assert derived(source) == ()


@pytest.mark.parametrize("scope", ["global", "nonlocal"])
def test_free_vector_resize_invalidates_earlier_extent(scope: Scope) -> None:
    source, _ = vector_source(scope, EXTENT + "del p[0]\n" + RESULT)
    assert derived(source) == ()
