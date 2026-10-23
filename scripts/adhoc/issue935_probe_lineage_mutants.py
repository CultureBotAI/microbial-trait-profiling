"""Probe lineage-review regressions against plausible implementation mutations.

Inputs: one mutation name and the repository's synthetic unit tests.
Outputs: pytest diagnostics and a nonzero exit if the regression survives.
Issue #935 / PR #891; 2026-10-23. Source bytes are restored even on failure.
Run each mutation in a fresh locked uv process; no production inputs are read.
"""

from __future__ import annotations

import argparse
from importlib.util import cache_from_source
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "src/mtp/influence_extension.py"
TEST = "tests/test_influence_extension.py"
MUTATIONS = {
    "ascending_rank": (
        'key=lambda row: (-row["count"], row["taxon"]),',
        'key=lambda row: (row["count"], row["taxon"]),',
        "test_ranked_panel_has_top_three_supported_lineages_with_real_tie_at_cutoff",
    ),
    "top_one": (
        'top = protocol["top_per_trait_direction_rank_scope"]',
        "top = 1",
        "test_ranked_panel_has_top_three_supported_lineages_with_real_tie_at_cutoff",
    ),
    "no_support_threshold": (
        'row["status"] == "resolved" and row["studies"] >= minimum',
        'row["status"] == "resolved"',
        "test_ranked_panel_has_top_three_supported_lineages_with_real_tie_at_cutoff",
    ),
    "summary_minimum": (
        '"minimum_effect": min(effects) if effects else "",',
        '"minimum_effect": max(effects) if effects else "",',
        "test_loso_summary_pins_actual_extrema_baseline_and_reversal_count",
    ),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mutation", choices=sorted(MUTATIONS))
    mutation = parser.parse_args().mutation
    before, after, test = MUTATIONS[mutation]
    original = SOURCE.read_bytes()
    text = original.decode("utf-8")
    if text.count(before) != 1:
        raise ValueError(f"Expected exactly one mutation target: {mutation}")
    try:
        SOURCE.write_text(text.replace(before, after), encoding="utf-8")
        # Same-length edits within one timestamp tick must not reuse cached code.
        Path(cache_from_source(str(SOURCE))).unlink(missing_ok=True)
        code = pytest.main(["-q", str(ROOT / TEST) + "::" + test])
    finally:
        SOURCE.write_bytes(original)
        Path(cache_from_source(str(SOURCE))).unlink(missing_ok=True)
    if SOURCE.read_bytes() != original:
        raise AssertionError("Mutation probe did not restore source bytes")
    if code != pytest.ExitCode.TESTS_FAILED:
        raise AssertionError(f"Mutation {mutation} was not rejected by its regression: {code}")
    print(f"Rejected {mutation}; exact source bytes restored.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
