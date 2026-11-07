"""Bound derivation analysis and preserve candidate headroom (#1581, #1683, #1700)."""

import re
from dataclasses import replace
from pathlib import Path
from textwrap import indent

import pytest
import yaml

from test_multiple_testing_registry import (
    DEFAULT_LIMITS,
    REGISTRY,
    ROOT,
    DiscoveryLimitError,
    discover,
    repository_sources,
)
from test_multiple_testing_registry_derived import PREFIX, source

DERIVATION_HEADROOM_DIVISOR = 16


def nested_derivation(depth: int) -> str:
    body = PREFIX + "result = np.clip(q * harmonic_sum, 0, 1)\n"
    for index in range(depth):
        body = f"for item_{index} in other:\n" + indent(body, "    ")
    return source(body + "return result\n")


def branching_derivation(bindings: int, branches: int) -> str:
    body = "".join(f"value_{index} = p\n" for index in range(bindings))
    body += "if other:\n    pass\n" * branches
    return source(body + PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n")


def test_default_budget_refuses_small_exponential_loop_tree() -> None:
    """This source is tiny; an AST-size limit alone cannot bound re-analysis."""
    text = nested_derivation(12)
    assert len(text) < 3000
    with pytest.raises(
        DiscoveryLimitError, match=r"^derivation work budget exceeded:.*src/x.py.*screen"
    ):
        discover({"src/x.py": text})


@pytest.mark.parametrize("text", [nested_derivation(2), branching_derivation(3, 2)])
def test_ordinary_derivations_still_provide_a_witness(text: str) -> None:
    discovered, _ = discover({"src/x.py": text})
    assert discovered[("src/x.py", "screen")]["derived_procedures"] == {"BY"}


def check_repository_derivation_headroom(root: Path) -> None:
    """Reserve function and file headroom, including unclaimed candidates."""
    registry = yaml.safe_load((root / REGISTRY.relative_to(ROOT)).read_text(encoding="utf-8"))
    anchors: set[tuple[str, str]] = {
        (anchor["path"], anchor["function"])
        for family in registry["families"]
        for anchor in family["anchors"]
        if "BY" in anchor.get("derived_procedures", {})
    }
    assert anchors, "The registry must supply real derived-BY headroom cases"
    limits = replace(
        DEFAULT_LIMITS,
        derivation_steps=DEFAULT_LIMITS.derivation_steps // DERIVATION_HEADROOM_DIVISOR,
        file_derivation_steps=DEFAULT_LIMITS.file_derivation_steps // DERIVATION_HEADROOM_DIVISOR,
    )
    discovered, _ = discover(repository_sources(root), limits=limits)
    # Claims provide a non-vacuity check, never the source-selection rule.
    for anchor in sorted(anchors):
        name = f"{anchor[0]}:{anchor[1]}"
        assert anchor in discovered, f"Missing derived-BY headroom anchor: {name}"
        assert "BY" in discovered[anchor].get("derived_procedures", set()), (
            f"Uncredited derived-BY headroom anchor: {name}"
        )


def test_all_repository_candidates_use_at_most_sixteenth_of_default_budget() -> None:
    try:
        check_repository_derivation_headroom(ROOT)
    except DiscoveryLimitError as error:
        raise AssertionError(
            "Every repository candidate must retain the declared headroom"
        ) from error


def test_custom_work_budget_refuses_before_a_witness_is_credited() -> None:
    limits = replace(DEFAULT_LIMITS, derivation_steps=1)
    with pytest.raises(DiscoveryLimitError, match=r"^derivation work budget exceeded:") as stopped:
        discover(
            {"src/x.py": source(PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n")}, limits=limits
        )
    message = str(stopped.value)
    assert "src/x.py:screen" in message
    # One module visit plus its three body entries is charged before traversal.
    # Report that attempted charge, not merely the limit or limit + 1.
    assert re.search(r"\battempted 4\b", message)
    assert re.search(r"\blimit 1\b", message)
    assert "Simplify" in message
    assert "Limits.derivation_steps" in message


def test_environment_work_counts_toward_budget() -> None:
    """Copies and joins cost more than visiting the short sequence of if nodes."""
    limits = replace(DEFAULT_LIMITS, derivation_steps=5_000)
    with pytest.raises(DiscoveryLimitError, match=r"^derivation work budget exceeded:"):
        discover({"src/x.py": branching_derivation(70, 30)}, limits=limits)


def test_non_candidate_does_not_use_derivation_budget() -> None:
    limits = replace(DEFAULT_LIMITS, derivation_steps=0)
    discover({"src/x.py": "def ordinary(value):\n    return value\n"}, limits=limits)


def test_module_ast_walk_counts_toward_derivation_budget() -> None:
    limits = replace(DEFAULT_LIMITS, derivation_steps=500)
    text = "0\n" * 300 + source(PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n")
    with pytest.raises(DiscoveryLimitError, match=r"^derivation work budget exceeded:"):
        discover({"src/x.py": text}, limits=limits)


def test_declared_names_are_charged_even_when_they_share_one_fact() -> None:
    """Global names are strings, not child AST nodes; duplicates keep env tiny."""
    limits = replace(DEFAULT_LIMITS, derivation_steps=500)
    names = ", ".join(["sentinel"] * 1000)
    text = source(f"global {names}\n" + PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n")
    with pytest.raises(DiscoveryLimitError, match=r"^derivation work budget exceeded:"):
        discover({"src/x.py": text}, limits=limits)


def test_module_walk_charges_non_ast_list_entries() -> None:
    """ast.iter_child_nodes visits string-list entries without yielding nodes."""
    limits = replace(DEFAULT_LIMITS, derivation_steps=500)
    names = ", ".join(["sentinel"] * 1000)
    text = f"global {names}\n" + source(PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n")
    with pytest.raises(DiscoveryLimitError, match=r"^derivation work budget exceeded:"):
        discover({"src/x.py": text}, limits=limits)


def test_witness_does_not_excuse_later_over_budget_work() -> None:
    limits = replace(DEFAULT_LIMITS, derivation_steps=500)
    body = PREFIX + "result = np.clip(q * harmonic_sum, 0, 1)\n"
    body += "global " + ", ".join(["sentinel"] * 1000) + "\nreturn result\n"
    with pytest.raises(DiscoveryLimitError, match=r"^derivation work budget exceeded:"):
        discover({"src/x.py": source(body)}, limits=limits)
