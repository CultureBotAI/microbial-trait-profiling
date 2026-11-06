"""Bound repeated derivation analysis, including small deeply nested ASTs (#1581)."""

from dataclasses import replace
from textwrap import indent

import pytest

from test_multiple_testing_registry import DEFAULT_LIMITS, DiscoveryLimitError, discover
from test_multiple_testing_registry_derived import PREFIX, source


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
    with pytest.raises(DiscoveryLimitError, match=r"derivation work budget.*src/x.py.*screen"):
        discover({"src/x.py": text})


@pytest.mark.parametrize("text", [nested_derivation(2), branching_derivation(3, 2)])
def test_ordinary_derivations_still_provide_a_witness(text: str) -> None:
    discovered, _ = discover({"src/x.py": text})
    assert discovered[("src/x.py", "screen")]["derived_procedures"] == {"BY"}


def test_custom_work_budget_refuses_before_a_witness_is_credited() -> None:
    limits = replace(DEFAULT_LIMITS, derivation_steps=1)
    with pytest.raises(DiscoveryLimitError, match="derivation work budget"):
        discover(
            {"src/x.py": source(PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n")}, limits=limits
        )


def test_environment_work_counts_toward_budget() -> None:
    """Copies and joins cost more than visiting the short sequence of if nodes."""
    limits = replace(DEFAULT_LIMITS, derivation_steps=5_000)
    with pytest.raises(DiscoveryLimitError, match="derivation work budget"):
        discover({"src/x.py": branching_derivation(70, 30)}, limits=limits)


def test_non_candidate_does_not_use_derivation_budget() -> None:
    limits = replace(DEFAULT_LIMITS, derivation_steps=0)
    discover({"src/x.py": "def ordinary(value):\n    return value\n"}, limits=limits)


def test_module_ast_walk_counts_toward_derivation_budget() -> None:
    limits = replace(DEFAULT_LIMITS, derivation_steps=500)
    text = "0\n" * 300 + source(PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n")
    with pytest.raises(DiscoveryLimitError, match="derivation work budget"):
        discover({"src/x.py": text}, limits=limits)


def test_declared_names_are_charged_even_when_they_share_one_fact() -> None:
    """Global names are strings, not child AST nodes; duplicates keep env tiny."""
    limits = replace(DEFAULT_LIMITS, derivation_steps=500)
    names = ", ".join(["sentinel"] * 1000)
    text = source(f"global {names}\n" + PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n")
    with pytest.raises(DiscoveryLimitError, match="derivation work budget"):
        discover({"src/x.py": text}, limits=limits)


def test_module_walk_charges_non_ast_list_entries() -> None:
    """ast.iter_child_nodes visits string-list entries without yielding nodes."""
    limits = replace(DEFAULT_LIMITS, derivation_steps=500)
    names = ", ".join(["sentinel"] * 1000)
    text = f"global {names}\n" + source(PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n")
    with pytest.raises(DiscoveryLimitError, match="derivation work budget"):
        discover({"src/x.py": text}, limits=limits)


def test_witness_does_not_excuse_later_over_budget_work() -> None:
    limits = replace(DEFAULT_LIMITS, derivation_steps=500)
    body = PREFIX + "result = np.clip(q * harmonic_sum, 0, 1)\n"
    body += "global " + ", ".join(["sentinel"] * 1000) + "\nreturn result\n"
    with pytest.raises(DiscoveryLimitError, match="derivation work budget"):
        discover({"src/x.py": source(body)}, limits=limits)
