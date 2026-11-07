"""Headroom includes unclaimed candidates and reserves more than a quarter (#1700)."""

from __future__ import annotations

import re
from dataclasses import replace
from pathlib import Path

import pytest
import yaml

import test_multiple_testing_registry as guard
import test_multiple_testing_registry_work as work


@pytest.fixture(autouse=True)
def assigned_guard_and_headroom_check() -> None:
    for module in (guard, work):
        assert module.__file__ is not None
        assert (
            Path(module.__file__).resolve()
            == (Path(__file__).parent / (module.__name__ + ".py")).resolve()
        )


def repository_fixture(
    root: Path, location: str, *, depth: int, claimed: bool
) -> tuple[dict[str, str], str]:
    """The registry is valid for its declared sources; extra is independently optional."""
    extra_path = f"{location}/extra.py"
    sources = {"src/x.py": work.nested_derivation(0), extra_path: work.nested_derivation(depth)}
    for path, text in sources.items():
        target = root / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")
    # A source omitted from both Git and the registry must still be analysed.
    (root / ".gitignore").write_text(extra_path + "\n", encoding="utf-8")
    anchors = [("src/x.py", "screen", ["BH", "BY"])]
    if claimed:
        anchors.append((extra_path, "screen", ["BH", "BY"]))
    registry = guard.example_registry(*anchors)
    for anchor in registry["families"][0]["anchors"]:
        anchor["derived_procedures"] = {
            "BY": "BH multiplied by the same family's harmonic sum and clipped to [0, 1]."
        }
    for counts in registry["literal_procedure_counts"]:
        counts["counts"] = {"BH": 1}
    registry_path = root / "methods/multiple_testing_families.yaml"
    registry_path.parent.mkdir()
    registry_path.write_text(yaml.safe_dump(registry), encoding="utf-8")
    assert dict(guard.repository_sources(root)) == sources
    return sources, extra_path


@pytest.mark.parametrize("location", ["src", "scripts"])
@pytest.mark.parametrize("claimed", [False, True], ids=["unclaimed", "claimed"])
def test_headroom_refuses_a_candidate_that_fits_the_full_and_quarter_budgets(
    tmp_path: Path, location: str, claimed: bool
) -> None:
    # This small depth-eight candidate fits the ordinary limit and its quarter,
    # but not the required sixteenth. No expected work count comes from the checker.
    sources, extra_path = repository_fixture(tmp_path, location, depth=8, claimed=claimed)
    try:
        discovered, _ = guard.discover(sources, limits=guard.DEFAULT_LIMITS)
    except guard.DiscoveryLimitError as error:
        raise AssertionError("The control candidate must fit the ordinary budget") from error
    assert discovered[(extra_path, "screen")]["derived_procedures"] == {"BY"}
    quarter_limits = replace(
        guard.DEFAULT_LIMITS, derivation_steps=guard.DEFAULT_LIMITS.derivation_steps // 4
    )
    try:
        quarter_discovered, _ = guard.discover(sources, limits=quarter_limits)
    except guard.DiscoveryLimitError as error:
        raise AssertionError("The control candidate must also fit the quarter budget") from error
    assert quarter_discovered[(extra_path, "screen")]["derived_procedures"] == {"BY"}
    with pytest.raises(
        guard.DiscoveryLimitError,
        match=rf"derivation work budget exceeded: {re.escape(extra_path)}:screen ",
    ):
        work.check_repository_derivation_headroom(tmp_path)


@pytest.mark.parametrize("location", ["src", "scripts"])
def test_cheap_unclaimed_candidates_remain_eligible_for_headroom_checking(
    tmp_path: Path, location: str
) -> None:
    sources, extra_path = repository_fixture(tmp_path, location, depth=0, claimed=False)
    discovered, _ = guard.discover(sources, limits=guard.DEFAULT_LIMITS)
    assert discovered[(extra_path, "screen")]["derived_procedures"] == {"BY"}
    try:
        work.check_repository_derivation_headroom(tmp_path)
    except guard.DiscoveryLimitError as error:
        raise AssertionError(
            "An unclaimed candidate with ample headroom must be accepted"
        ) from error
