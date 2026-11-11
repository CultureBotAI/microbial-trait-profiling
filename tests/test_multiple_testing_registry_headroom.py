"""Headroom includes unclaimed candidates and reserves more than a quarter (#1700)."""

from __future__ import annotations

import re
from dataclasses import replace
from pathlib import Path
from typing import Any

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
    root: Path, location: str, *, depth: int, claimed: bool, padding: int = 0
) -> tuple[dict[str, str], str]:
    """The registry is valid for its declared sources; extra is independently optional."""
    extra_path = f"{location}/extra.py"
    sources = {"src/x.py": work.nested_derivation(0), extra_path: work.nested_derivation(depth)}
    if padding:
        header = "def screen(p, other):\n"
        assert sources[extra_path].count(header) == 1
        sources[extra_path] = sources[extra_path].replace(
            header, header + "    pass\n" * padding, 1
        )
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
        match=rf"^derivation work budget exceeded: {re.escape(extra_path)}:screen ",
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


def test_declared_headroom_divisor_is_exactly_sixteen() -> None:
    assert getattr(work, "DERIVATION_HEADROOM_DIVISOR", None) == 16


def test_headroom_discovery_uses_exactly_one_sixteenth(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources, _ = repository_fixture(tmp_path, "scripts", depth=0, claimed=False)
    calls: list[tuple[dict[str, str], guard.Limits]] = []
    original = guard.discover

    def record_discovery(
        inputs: guard.Sources, *, limits: guard.Limits = guard.DEFAULT_LIMITS
    ) -> tuple[dict[tuple[str, str], dict[str, Any]], set[tuple[str, str]]]:
        materialized = dict(inputs)
        calls.append((materialized, limits))
        return original(materialized, limits=limits)

    monkeypatch.setattr(work, "discover", record_discovery)
    work.check_repository_derivation_headroom(tmp_path)
    assert calls == [
        (
            sources,
            replace(
                guard.DEFAULT_LIMITS,
                derivation_steps=guard.DEFAULT_LIMITS.derivation_steps // 16,
                file_derivation_steps=guard.DEFAULT_LIMITS.file_derivation_steps // 16,
            ),
        )
    ]


def test_padded_candidate_fits_fifteenth_but_exceeds_headroom(tmp_path: Path) -> None:
    # This narrower witness distinguishes a sixteenth from even a fifteenth;
    # the original depth-eight fixture cannot distinguish an eighth.
    # A pass costs one statement and one AST visit. The 32,000 passes
    # therefore cost 64,000 units independently of alias-model bookkeeping:
    # above 62,500 (a sixteenth), below 66,666 (a fifteenth), with room for BY.
    sources, extra_path = repository_fixture(
        tmp_path, "scripts", depth=0, claimed=False, padding=32_000
    )
    fifteenth = replace(
        guard.DEFAULT_LIMITS, derivation_steps=guard.DEFAULT_LIMITS.derivation_steps // 15
    )
    try:
        discovered, _ = guard.discover(sources, limits=fifteenth)
    except guard.DiscoveryLimitError as error:
        raise AssertionError("The boundary control must fit the fifteenth budget") from error
    assert discovered[(extra_path, "screen")]["derived_procedures"] == {"BY"}
    with pytest.raises(
        guard.DiscoveryLimitError,
        match=rf"^derivation work budget exceeded: {re.escape(extra_path)}:screen ",
    ):
        work.check_repository_derivation_headroom(tmp_path)


def multiple_candidate_repository_fixture(root: Path, copies: int) -> tuple[dict[str, str], str]:
    """Put independently named candidates in one ignored, unclaimed source file."""
    # Each candidate contributes 24,000 known pass units plus a short BY
    # witness. Two fit 62,500; three exceed it while remaining below 125,000.
    sources, extra_path = repository_fixture(
        root, "scripts", depth=0, claimed=False, padding=12_000
    )
    header, separator, function = sources[extra_path].partition("def screen")
    assert separator == "def screen"
    sources[extra_path] = header + "".join(
        f"def screen_{index}" + function for index in range(copies)
    )
    (root / extra_path).write_text(sources[extra_path], encoding="utf-8")
    assert dict(guard.repository_sources(root)) == sources
    return sources, extra_path


def test_two_candidates_in_one_file_retain_both_sixteenth_allowances(tmp_path: Path) -> None:
    sources, extra_path = multiple_candidate_repository_fixture(tmp_path, 2)
    limits = replace(
        guard.DEFAULT_LIMITS,
        derivation_steps=guard.DEFAULT_LIMITS.derivation_steps // 16,
        file_derivation_steps=guard.DEFAULT_LIMITS.file_derivation_steps // 16,
    )
    try:
        discovered, _ = guard.discover(sources, limits=limits)
    except guard.DiscoveryLimitError as error:
        raise AssertionError("Two candidates must fit both sixteenth allowances") from error
    assert {anchor: site["derived_procedures"] for anchor, site in discovered.items()} == {
        ("src/x.py", "screen"): {"BY"},
        (extra_path, "screen_0"): {"BY"},
        (extra_path, "screen_1"): {"BY"},
    }
    try:
        work.check_repository_derivation_headroom(tmp_path)
    except guard.DiscoveryLimitError as error:
        raise AssertionError("The headroom helper must admit the two-candidate file") from error


def test_file_headroom_refuses_three_individually_small_candidates(tmp_path: Path) -> None:
    # Each candidate fits the function sixteenth; together they fit a file eighth
    # but exceed a file sixteenth. No expected work count is read from the checker.
    sources, extra_path = multiple_candidate_repository_fixture(tmp_path, 3)
    limits = replace(
        guard.DEFAULT_LIMITS,
        derivation_steps=guard.DEFAULT_LIMITS.derivation_steps // 16,
        file_derivation_steps=guard.DEFAULT_LIMITS.file_derivation_steps // 8,
    )
    try:
        discovered, _ = guard.discover(sources, limits=limits)
    except guard.DiscoveryLimitError as error:
        raise AssertionError(
            "Three candidates must fit the function sixteenth and file eighth"
        ) from error
    assert {anchor: site["derived_procedures"] for anchor, site in discovered.items()} == {
        ("src/x.py", "screen"): {"BY"},
        (extra_path, "screen_0"): {"BY"},
        (extra_path, "screen_1"): {"BY"},
        (extra_path, "screen_2"): {"BY"},
    }
    with pytest.raises(
        guard.DiscoveryLimitError,
        match=rf"^per-file derivation work budget exceeded: {re.escape(extra_path)}:screen_2 ",
    ):
        work.check_repository_derivation_headroom(tmp_path)


def test_missing_claimed_anchor_is_reported_by_an_assertion(tmp_path: Path) -> None:
    repository_fixture(tmp_path, "scripts", depth=0, claimed=False)
    (tmp_path / "src/x.py").write_text("def unrelated(value):\n    return value\n")
    discovered, _ = guard.discover(guard.repository_sources(tmp_path))
    assert ("src/x.py", "screen") not in discovered
    try:
        with pytest.raises(
            AssertionError, match=r"Missing derived-BY headroom anchor: src/x\.py:screen"
        ):
            work.check_repository_derivation_headroom(tmp_path)
    except KeyError as error:
        raise AssertionError(
            "A missing claimed anchor must produce a named assertion, not KeyError"
        ) from error


def test_present_but_uncredited_claim_names_its_anchor(tmp_path: Path) -> None:
    repository_fixture(tmp_path, "scripts", depth=0, claimed=False)
    (tmp_path / "src/x.py").write_text(
        "from mtp.enrichment_stats import benjamini_hochberg\n"
        "def screen(p):\n    return benjamini_hochberg(p)\n"
    )
    discovered, _ = guard.discover(guard.repository_sources(tmp_path))
    assert ("src/x.py", "screen") in discovered
    assert "BY" not in discovered[("src/x.py", "screen")].get("derived_procedures", set())
    with pytest.raises(
        AssertionError, match=r"Uncredited derived-BY headroom anchor: src/x\.py:screen"
    ):
        work.check_repository_derivation_headroom(tmp_path)


def test_headroom_requires_at_least_one_declared_derived_anchor(tmp_path: Path) -> None:
    repository_fixture(tmp_path, "scripts", depth=0, claimed=False)
    registry_path = tmp_path / guard.REGISTRY.relative_to(guard.ROOT)
    registry = yaml.safe_load(registry_path.read_text(encoding="utf-8"))
    for family in registry["families"]:
        for anchor in family["anchors"]:
            anchor.pop("derived_procedures", None)
    registry_path.write_text(yaml.safe_dump(registry), encoding="utf-8")
    with pytest.raises(
        AssertionError, match="The registry must supply real derived-BY headroom cases"
    ):
        work.check_repository_derivation_headroom(tmp_path)
