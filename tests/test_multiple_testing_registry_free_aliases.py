"""Free-vector alias writes must not establish a stale BY witness (#1671)."""

from __future__ import annotations

import importlib.util
from fractions import Fraction
from pathlib import Path
from textwrap import indent
from typing import Literal

import pytest

import test_multiple_testing_registry as guard

Scope = Literal["module", "closure"]
IMPORTS = "import numpy as np\nfrom mtp.enrichment_stats import benjamini_hochberg\n"
EXTENT = "m = len(p)\nh = np.sum(1.0 / np.arange(1, m + 1))\n"
RESULT = "q = benjamini_hochberg(p)\nreturn np.clip(q * h, 0, 1)\n"
MUTATIONS = (
    "alias = p\ndel alias[0]\n",
    "alias = p\nalias += [0.5, 0.6]\n",
    "first = p\nsecond = first\ndel second[0]\n",
    "alias = p\nalias[:] = alias[:2]\n",
)


@pytest.fixture(autouse=True)
def assigned_guard() -> None:
    assert guard.__file__ is not None
    assert (
        Path(guard.__file__).resolve()
        == Path(__file__).with_name("test_multiple_testing_registry.py").resolve()
    )


def free_source(scope: Scope, body: str) -> tuple[str, str]:
    function = "def screen():\n" + indent(body, "    ")
    if scope == "module":
        return IMPORTS + "p = [0.01, 0.02, 0.03, 0.2]\n" + function, "screen"
    return (
        IMPORTS + "def outer(p):\n" + indent(function + "return screen\n", "    "),
        "outer.screen",
    )


def derived(source: str) -> tuple[str, ...]:
    return guard.summarize_source(
        "src/alias_witness.py", source, guard.SummaryBudget(guard.DEFAULT_LIMITS)
    ).derived


@pytest.mark.parametrize("scope", ["module", "closure"])
@pytest.mark.parametrize(
    "mutation", MUTATIONS, ids=["delete", "augmented", "chain-delete", "slice-write"]
)
def test_free_vector_alias_mutation_cannot_reuse_its_previous_extent(
    scope: Scope, mutation: str
) -> None:
    text, _ = free_source(scope, EXTENT + mutation + RESULT)
    assert derived(text) == ()


@pytest.mark.parametrize("scope", ["module", "closure"])
@pytest.mark.parametrize("use_alias", [False, True], ids=["original", "alias-chain"])
def test_untouched_free_vector_aliases_keep_the_same_family_witness(
    scope: Scope, use_alias: bool
) -> None:
    result = (
        RESULT.replace("benjamini_hochberg(p)", "benjamini_hochberg(second)")
        if use_alias
        else RESULT
    )
    text, anchor = free_source(scope, EXTENT + "first = p\nsecond = first\n" + result)
    assert derived(text) == (anchor,)


@pytest.mark.parametrize("scope", ["module", "closure"])
@pytest.mark.parametrize("saved_callable", ["saved_np = np\n", "saved_bh = benjamini_hochberg\n"])
def test_saving_a_callable_alias_does_not_rebind_the_original_import(
    scope: Scope, saved_callable: str
) -> None:
    text, anchor = free_source(scope, saved_callable + EXTENT + RESULT)
    assert derived(text) == (anchor,)


@pytest.mark.parametrize("scope", ["module", "closure"])
@pytest.mark.parametrize("read_location", ["branch", "loop"])
def test_free_vector_read_does_not_become_a_mutation_at_a_control_flow_join(
    scope: Scope, read_location: str
) -> None:
    if read_location == "branch":
        body = "if p is not None:\n    observed = len(p)\n" + EXTENT + RESULT
    else:
        body = "for item in (0,):\n" + indent(EXTENT + RESULT, "    ")
    text, anchor = free_source(scope, body)
    assert derived(text) == (anchor,)


@pytest.mark.parametrize("scope", ["module", "closure"])
@pytest.mark.parametrize("alias_action", ["delete-name", "rebind-before-write"])
def test_changing_an_alias_name_without_mutating_the_free_vector_is_safe(
    scope: Scope, alias_action: str
) -> None:
    action = (
        "del alias\n" if alias_action == "delete-name" else "alias = [0.5, 0.6]\ndel alias[0]\n"
    )
    text, anchor = free_source(scope, EXTENT + "alias = p\n" + action + RESULT)
    assert derived(text) == (anchor,)


def test_mutating_an_old_parameter_alias_does_not_invalidate_a_rebound_vector() -> None:
    body = "old = p\np = other\n" + EXTENT + "del old[0]\n" + RESULT
    text = IMPORTS + "def screen(p, other):\n" + indent(body, "    ")
    assert derived(text) == ("screen",)


@pytest.mark.parametrize("scope", ["module", "closure"])
@pytest.mark.parametrize("container", ["branch", "loop"])
def test_free_alias_mutation_inside_control_flow_invalidates_the_previous_extent(
    scope: Scope, container: str
) -> None:
    header = "if p is not None:\n" if container == "branch" else "for item in (0,):\n"
    body = EXTENT + header + indent(MUTATIONS[0], "    ") + RESULT
    text, _ = free_source(scope, body)
    assert derived(text) == ()


@pytest.mark.parametrize(
    "mutation", MUTATIONS, ids=["delete", "augmented", "chain-delete", "slice-write"]
)
def test_parameter_alias_mutation_remains_rejected(mutation: str) -> None:
    text = IMPORTS + "def screen(p):\n" + indent(EXTENT + mutation + RESULT, "    ")
    assert derived(text) == ()


def augmented_flow_source(
    scope: Literal["parameter", "module", "closure"], container: str, binding: str
) -> tuple[str, str]:
    if container == "try":
        block = "try:\n    alias += [0.5]\nexcept Exception:\n    pass\n"
    else:
        block = "match 0:\n    case _:\n        alias += [0.5]\n"
    body = EXTENT + binding + block + RESULT
    if scope == "parameter":
        return IMPORTS + "def screen(p):\n" + indent(body, "    "), "screen"
    return free_source(scope, body)


@pytest.mark.parametrize("scope", ["parameter", "module", "closure"])
@pytest.mark.parametrize("container", ["try", "match"])
def test_alias_augmented_assignment_in_unsupported_flow_invalidates_original(
    scope: Literal["parameter", "module", "closure"], container: str
) -> None:
    text, _ = augmented_flow_source(scope, container, "alias = p\n")
    assert derived(text) == ()


@pytest.mark.parametrize("scope", ["parameter", "module", "closure"])
@pytest.mark.parametrize("container", ["try", "match"])
def test_unsupported_flow_mutating_an_unrelated_allocation_keeps_the_witness(
    scope: Literal["parameter", "module", "closure"], container: str
) -> None:
    # The same augmented write changes a fresh list, not the vector used for BH.
    text, anchor = augmented_flow_source(scope, container, "alias = [0.25]\n")
    assert derived(text) == (anchor,)


@pytest.mark.parametrize("scope", ["module", "closure"])
@pytest.mark.parametrize(
    "mutation",
    ["arr = np.ones(m)\ndel arr[0]\n", "m += 1\n"],
    ids=["separate-allocation", "saved-size"],
)
def test_mutating_a_nonvector_fact_does_not_rebind_the_free_vector(
    scope: Scope, mutation: str
) -> None:
    # Both facts carry p's extent key, but neither write changes p. The saved
    # harmonic sum remains valid even when its earlier size variable is rebound.
    text, anchor = free_source(scope, EXTENT + mutation + RESULT)
    assert derived(text) == (anchor,)


@pytest.mark.parametrize("scope", ["module", "closure"])
def test_direct_free_vector_deletion_remains_rejected(scope: Scope) -> None:
    text, _ = free_source(scope, EXTENT + "del p[0]\n" + RESULT)
    assert derived(text) == ()


@pytest.mark.parametrize("scope", ["module", "closure"])
def test_explicit_scope_declaration_remains_conservative(scope: Scope) -> None:
    declaration = "global" if scope == "module" else "nonlocal"
    text, _ = free_source(scope, f"{declaration} p\n" + EXTENT + MUTATIONS[0] + RESULT)
    assert derived(text) == ()


def test_free_series_alias_deletion_runs_to_the_wrong_by_and_is_not_credited(
    tmp_path: Path,
) -> None:
    """The repository's Series-returning BH makes this a numerical witness, not a TypeError."""
    from mtp.cellulose_decr import analysis

    text = (
        "import numpy as np\nimport pandas as pd\n"
        "from mtp.cellulose_decr.analysis import benjamini_hochberg\n"
        "p = pd.Series([0.01, 0.02, 0.03, 0.2])\n"
        "def screen():\n" + indent(EXTENT + MUTATIONS[0] + RESULT, "    ")
    )
    fixture_path = tmp_path / "free_series_alias_witness.py"
    fixture_path.write_text(text, encoding="utf-8")
    spec = importlib.util.spec_from_file_location("free_series_alias_witness", fixture_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.benjamini_hochberg is analysis.benjamini_hochberg
    assert Path(analysis.__file__).resolve() == (
        Path(__file__).resolve().parents[1] / "src/mtp/cellulose_decr/analysis.py"
    )
    try:
        result = module.screen()
    except (TypeError, ValueError, IndexError, KeyError) as error:
        raise AssertionError(
            "The real Series alias witness must run to a numeric result"
        ) from error
    assert list(module.p.index) == [1, 2, 3]
    assert list(module.p) == [0.02, 0.03, 0.2]
    assert list(result.index) == [1, 2, 3]

    # Independent BY definition on the three remaining, already sorted p-values:
    # min over ranks j >= i of n * H_n * p_j / j, capped at one.
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
    # Guard exactly the bytes which were executed, not a simplified stand-in.
    assert fixture_path.read_text(encoding="utf-8") == text
    assert derived(text) == ()
