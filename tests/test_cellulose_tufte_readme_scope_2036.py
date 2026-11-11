"""Verify README scope and unpinned-token disposition (#1997, #2036, #2037, #2038)."""

from __future__ import annotations

from pathlib import Path

import pytest

README_PATH = (
    Path(__file__).resolve().parent
    / "fixtures"
    / "cellulose_tufte_historical_sources"
    / "README.md"
)

SCOPE_HEADING = "## Scope of the source pins (#1596)"
KEPT_SCOPE_OPENING = (
    "The table below lists what `test_cellulose_tufte_statistics_sources.py` binds to\n"
    "identified historical sources; the other test files are summarised after it."
)
KEPT_SCOPE_LIMITATION = "They do not pin every numeric token or every repeated occurrence,"
REPEATED_SCOPE_SENTENCE = (
    "These tests bind selected displayed values and labels to identified historical sources."
)

OVERLAP_SUMMARY = (
    "`test_cellulose_tufte_overlap_sources.py` pins selected set sizes, direction labels, "
    "percentages, distances and Plate III's nominal p-values for Plates I to IV and "
    "figures 7 to 9; it is summarised here, not tabulated."
)
STALE_OVERLAP_SUMMARY = (
    "`test_cellulose_tufte_overlap_sources.py` pins selected set sizes and direction labels"
)

MUTATION_SUMMARY = (
    "`test_cellulose_tufte_source_pin_mutations.py` checks that every tabulated pin, "
    "and one overlap pin, fails on targeted mutations."
)
STALE_MUTATION_SUMMARY = "checks that the pins above fail on targeted mutations"

STALE_DEFERRED_TAIL = "Broader source pins remain a separate task."
UNPINNED_TOKENS: tuple[str, ...] = (
    "10 PD\u2229CRC",
    "= 40 taxa",
    "46 PD",
    "27 shared (57.4%)",
    "17 &lt; 27 (PD\u2013CRC)",
    "17 of 32",
    "45%\u201357%",
    "50.4% shared",
    "n = 2,335 pairs",
    "Both FA (58.2% up",
)
UNPINNED_TOKEN_GROUPS: tuple[str, ...] = (
    "'10 PD\u2229CRC', '= 40 taxa', '46 PD', and '27 shared (57.4%)' in fig7",
    "'17 &lt; 27 (PD\u2013CRC)' in fig9",
    "'17 of 32' in Plate I",
    "'45%\u201357%' in Plate III",
    "'50.4% shared' in Plate IV",
    "'n = 2,335 pairs' in Plate V",
    "'Both FA (58.2% up' in Plate VI",
)
UNPINNED_DISPOSITION = (
    "won't pin (#1997), because the plates were retired from default generation by #1693 (#1572)."
)


def _normalize_whitespace(text: str) -> str:
    return " ".join(text.split())


def check_historical_sources_readme_scope(text: str) -> None:
    """Assert #1997, #2036, #2037, and #2038 invariants on README content."""
    assert SCOPE_HEADING in text
    scope_section = text.split(SCOPE_HEADING, 1)[1]
    normalized = _normalize_whitespace(scope_section)

    # #2038: drop repeated sentence while keeping the opening and limitation sentence.
    assert KEPT_SCOPE_OPENING in scope_section
    assert KEPT_SCOPE_LIMITATION in scope_section
    assert REPEATED_SCOPE_SENTENCE not in normalized
    assert scope_section.count("identified historical sources") == 1

    # #2036: accurate summary of test_cellulose_tufte_overlap_sources.py.
    assert OVERLAP_SUMMARY in scope_section
    assert OVERLAP_SUMMARY in normalized
    assert STALE_OVERLAP_SUMMARY not in normalized

    # #2037: accurate summary of test_cellulose_tufte_source_pin_mutations.py.
    assert MUTATION_SUMMARY in scope_section
    assert MUTATION_SUMMARY in normalized
    assert STALE_MUTATION_SUMMARY not in normalized

    # #1997: explicit list of 10 sourced-but-unpinned tokens and won't-pin disposition.
    assert STALE_DEFERRED_TAIL not in normalized
    assert len(UNPINNED_TOKENS) == 10
    for token in UNPINNED_TOKENS:
        assert token in scope_section, f"Missing sourced-but-unpinned token: {token}"
    for group in UNPINNED_TOKEN_GROUPS:
        assert group in normalized, f"Missing figure/plate token group: {group}"
    assert UNPINNED_DISPOSITION in scope_section
    assert UNPINNED_DISPOSITION in normalized


def test_historical_sources_readme_scope_and_unpinned_disposition() -> None:
    text = README_PATH.read_text(encoding="utf-8")
    check_historical_sources_readme_scope(text)


@pytest.mark.parametrize(
    ("target", "replacement"),
    [
        (
            KEPT_SCOPE_OPENING,
            f"{KEPT_SCOPE_OPENING} {REPEATED_SCOPE_SENTENCE}",
        ),
        (
            OVERLAP_SUMMARY,
            "`test_cellulose_tufte_overlap_sources.py` pins selected set sizes and direction "
            "labels for Plates I to IV and figures 7 to 9; it is summarised here, not tabulated.",
        ),
        (
            "percentages, distances and Plate III's nominal p-values",
            "percentages and distances",
        ),
        (
            MUTATION_SUMMARY,
            "`test_cellulose_tufte_source_pin_mutations.py` checks that the pins above fail on "
            "targeted mutations.",
        ),
        (
            UNPINNED_DISPOSITION,
            STALE_DEFERRED_TAIL,
        ),
        *((token, "REDACTED_TOKEN") for token in UNPINNED_TOKENS),
    ],
)
def test_check_historical_sources_readme_scope_rejects_regressions(
    target: str, replacement: str
) -> None:
    text = README_PATH.read_text(encoding="utf-8")
    assert target in text
    mutated = text.replace(target, replacement, 1)
    with pytest.raises(AssertionError):
        check_historical_sources_readme_scope(mutated)
