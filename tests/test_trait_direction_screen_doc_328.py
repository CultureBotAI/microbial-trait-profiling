"""Keep #328's historical matrix observation distinct from execution evidence."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

METHODS = Path(__file__).resolve().parents[1] / "methods"


def _section(text: str, title: str) -> str:
    match = re.search(rf"^## {re.escape(title)}\n(.*?)(?=^## |\Z)", text, re.MULTILINE | re.DOTALL)
    assert match is not None, f"Missing methods section: {title}"
    return match.group(1)


def _normalized(text: str) -> str:
    return " ".join(text.replace("**", "").replace("`", "").split())


def _candidate_family() -> str:
    text = (METHODS / "trait_direction_screen.md").read_text(encoding="utf-8")
    return _normalized(_section(text, "The candidate family"))


def _assert_historical_matrix_counts_name_their_population(doc_text: str) -> None:
    text = _normalized(doc_text)
    family = _normalized(_section(doc_text, "The candidate family"))
    assert "historical v1.3.1" in family
    assert "historical profile matrix" in family
    assert "5,829 input profiles" in family
    assert re.search(r"2,290 (?:excluded MediaDive columns|MediaDive columns excluded)", family)
    assert re.search(r"317 (?:retained traits|candidates)", family)
    assert re.search(r"4,670 paired(?: human)? profiles", family)
    assert "2,335 experiments" in family
    assert "959 studies" in family
    assert "2,333" not in family
    assert "2,333" not in text


def test_historical_matrix_counts_name_their_population() -> None:
    doc_text = (METHODS / "trait_direction_screen.md").read_text(encoding="utf-8")
    _assert_historical_matrix_counts_name_their_population(doc_text)


@pytest.mark.parametrize(
    ("anchor", "replacement"),
    [
        pytest.param(
            "## Why the family changed",
            "## Why the family changed\n\nHistorical typo: 2,333 experiments.",
            id="before-candidate-family",
        ),
        pytest.param(
            "## Outputs",
            "## Outputs\n\nHistorical typo: **2,333** experiments.",
            id="after-candidate-family-bold",
        ),
    ],
)
def test_historical_matrix_counts_rejects_2333_outside_candidate_family(
    anchor: str, replacement: str
) -> None:
    """#2196: `"2,333"` must fail anywhere in `trait_direction_screen.md`, not just one section."""
    doc_text = (METHODS / "trait_direction_screen.md").read_text(encoding="utf-8")
    assert doc_text.count(anchor) == 1
    mutated = doc_text.replace(anchor, replacement, 1)
    with pytest.raises(AssertionError):
        _assert_historical_matrix_counts_name_their_population(mutated)


def test_historical_matrix_identity_has_a_pinned_source() -> None:
    orientation = (METHODS / "contrast_orientation.md").read_text(encoding="utf-8")
    historical = _section(orientation, "Changes from `contrast-orientation-v3` to `-v4` (#569)")
    profile_digests = re.findall(r"profile matrix SHA-256\s+`([0-9a-f]{64})`", historical)
    assert len(profile_digests) == 1, "The source must identify one complete profile digest"
    assert profile_digests[0].startswith("c07d589f")
    family = _candidate_family()
    assert profile_digests[0] in family, "Candidate counts must name the same full matrix digest"
    revisions = re.findall(
        r"https://github\.com/CultureBotAI/microbial-trait-profiling/blob/"
        r"([0-9a-f]{40})/methods/table1_log_ratio_validation\.md(?:#[^\s)]*)?",
        family,
    )
    assert len(revisions) == 1, "Use one commit-pinned historical comparison record"
    assert revisions[0].startswith("de1a54b9")


def test_matrix_observation_does_not_reconstruct_the_screen_execution() -> None:
    family = _candidate_family().lower()
    assert re.search(
        r"original directional screen (?:still )?lacks a complete execution manifest", family
    )
    assert "matching these counts does not reconstruct it" in family


def test_historical_counts_are_not_current_gzip_or_future_constants() -> None:
    family = _candidate_family().lower()
    assert "not the currently tracked gzip" in family
    assert "or a future run" in family
    assert "three_arm_statistics.md#vintage" in family


def test_historical_screen_results_and_later_phenotype_context_are_preserved() -> None:
    text = _normalized((METHODS / "trait_direction_screen.md").read_text(encoding="utf-8"))
    assert "historical off screen recorded by the protocol amendment" in text
    assert "m = 317 and 18 traits met BY q <= 0.05" in text
    assert "m = 112 by Li & Ji on the study-contrast basis, 56 traits" in text
    family = _candidate_family()
    assert re.search(r"phenotype.{0,50}\b344\b", family)
    assert "5,629-profile table" in family
    assert "../docs/reviews/2026-10-25-rebuild-of-F-record.md" in family
