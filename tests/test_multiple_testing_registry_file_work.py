"""Bound the combined derivation work of candidates in one file (#1685)."""

from __future__ import annotations

import re
from dataclasses import fields, replace
from pathlib import Path
from typing import Any

import pytest

import test_multiple_testing_registry as guard
from test_multiple_testing_registry_derived import PREFIX, source

TWO_CANDIDATES = (
    "import numpy as np\n"
    "from mtp.enrichment_stats import benjamini_hochberg\n"
    "def first(p):\n"
    "    return np.clip(benjamini_hochberg(p), 0, 1)\n"
    "def second(p):\n"
    "    return np.clip(benjamini_hochberg(p), 0, 1)\n"
)


@pytest.fixture(autouse=True)
def assigned_guard() -> None:
    assert guard.__file__ is not None
    assert (
        Path(guard.__file__).resolve()
        == (Path(__file__).parent / "test_multiple_testing_registry.py").resolve()
    )


def configured(file_steps: int, function_steps: int = 1_000_000) -> guard.Limits:
    assert "file_derivation_steps" in {item.name for item in fields(guard.Limits)}, (
        "A separate per-file derivation allowance must exist"
    )
    return replace(
        guard.DEFAULT_LIMITS,
        file_derivation_steps=file_steps,
        derivation_steps=function_steps,
    )


def accepted(
    sources: dict[str, str], limits: guard.Limits
) -> dict[tuple[str, str], dict[str, Any]]:
    try:
        sites, _ = guard.discover(sources, limits=limits)
    except guard.DiscoveryLimitError as error:
        raise AssertionError(
            "The independently bounded positive control must be accepted"
        ) from error
    return sites


def assert_file_refusal(message: str, attempted: int, limit: int) -> None:
    assert message.startswith("per-file derivation work budget exceeded: ")
    assert "src/two.py:second" in message
    assert re.search(rf"\battempted {attempted}\b", message)
    assert re.search(rf"\blimit {limit}\b", message)
    assert "Limits.file_derivation_steps" in message


def test_default_file_allowance_is_bounded_without_lowering_the_function_allowance() -> None:
    assert "file_derivation_steps" in {item.name for item in fields(guard.Limits)}
    assert guard.DEFAULT_LIMITS.file_derivation_steps == 1_000_000
    assert guard.DEFAULT_LIMITS.derivation_steps == 1_000_000


def test_two_candidates_fit_exactly_sixty_eight_units_in_one_file() -> None:
    # Per candidate: module walk 13; return subtree walk 17; statement 1;
    # clip, BH call and p expression 3. Two candidates therefore cost 68.
    sites = accepted({"src/two.py": TWO_CANDIDATES}, configured(68, 34))
    assert set(sites) == {("src/two.py", "first"), ("src/two.py", "second")}
    assert all(site["procedures"] == {"BH"} for site in sites.values())
    assert all(not site.get("derived_procedures") for site in sites.values())


@pytest.mark.parametrize(("limit", "attempted"), [(34, 35), (67, 68)])
def test_second_candidate_uses_the_first_candidates_remaining_file_allowance(
    limit: int, attempted: int
) -> None:
    with pytest.raises(guard.DiscoveryLimitError) as stopped:
        guard.discover({"src/two.py": TWO_CANDIDATES}, limits=configured(limit, 34))
    assert_file_refusal(str(stopped.value), attempted, limit)


@pytest.mark.parametrize("file_limit", [33, 68])
def test_existing_per_function_refusal_keeps_its_precedence_and_diagnostic(file_limit: int) -> None:
    with pytest.raises(guard.DiscoveryLimitError) as stopped:
        guard.discover({"src/two.py": TWO_CANDIDATES}, limits=configured(file_limit, 33))
    message = str(stopped.value)
    assert message.startswith("derivation work budget exceeded: src/two.py:first ")
    assert "attempted 34, limit 33" in message
    assert "Limits.derivation_steps" in message
    assert "file_derivation_steps" not in message


def test_zero_file_allowance_refuses_the_first_attempted_candidate_node() -> None:
    with pytest.raises(guard.DiscoveryLimitError) as stopped:
        guard.discover({"src/two.py": TWO_CANDIDATES}, limits=configured(0))
    message = str(stopped.value)
    assert message.startswith("per-file derivation work budget exceeded: ")
    assert "src/two.py:first" in message
    assert "attempted 1, limit 0" in message


def test_non_candidates_do_not_consume_the_file_derivation_allowance() -> None:
    text = "def ordinary(value):\n    return value\n"
    assert accepted({"src/ordinary.py": text}, configured(0, 0)) == {}


def test_separate_files_each_receive_their_own_allowance_within_one_discovery() -> None:
    sites = accepted(
        {"src/one.py": TWO_CANDIDATES, "scripts/two.py": TWO_CANDIDATES}, configured(68, 34)
    )
    assert set(sites) == {
        ("src/one.py", "first"),
        ("src/one.py", "second"),
        ("scripts/two.py", "first"),
        ("scripts/two.py", "second"),
    }


def test_later_candidate_refusal_cannot_return_a_partial_compact_summary() -> None:
    budget = guard.SummaryBudget(configured(67, 34))
    with pytest.raises(guard.DiscoveryLimitError) as stopped:
        guard.summarize_source("src/two.py", TWO_CANDIDATES, budget)
    assert_file_refusal(str(stopped.value), 68, 67)


def test_normal_by_witness_remains_accepted_with_both_default_limits() -> None:
    text = source(PREFIX + "return np.clip(q * harmonic_sum, 0, 1)\n")
    sites = accepted({"src/normal.py": text}, configured(1_000_000))
    assert sites[("src/normal.py", "screen")]["derived_procedures"] == {"BY"}
