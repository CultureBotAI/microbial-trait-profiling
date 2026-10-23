"""Probe lineage-review regressions against plausible implementation mutations.

Inputs: one mutation name and the repository's synthetic unit tests.
Outputs: pytest diagnostics and a nonzero exit if the regression survives.
Issues #935/#960 / PR #891; 2026-10-23. Mutates an isolated module copy only.
Run each mutation in a fresh locked uv process; no production inputs are read.
"""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

import mtp
from mtp import influence_extension

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "src/mtp/influence_extension.py"
TEST = "tests/test_influence_extension.py"
JOINT_TEST = "tests/test_influence_extension_joint.py"
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
    "joint_unchanged": (
        "remaining, ambiguous = joint_removal(\n"
        '                base, graph, ledger, json.loads(scenario["lineages"])\n'
        "            )",
        "remaining, ambiguous = base.profiles, {}",
        JOINT_TEST + "::test_joint_effects_match_exact_hand_calculation",
    ),
    "joint_as_single": (
        "remaining, ambiguous = joint_removal(\n"
        '                base, graph, ledger, json.loads(scenario["lineages"])\n'
        "            )",
        'remaining, _, ambiguous = removal(base, graph, ledger, scenario["lineage"])',
        JOINT_TEST + "::test_joint_effects_match_exact_hand_calculation",
    ),
    "effects_q_columns": (
        'row[f"{population}_q_bh_six"] = "" if bh is None else bh\n'
        '                row[f"{population}_q_by_six"] = "" if by is None else by',
        'row[f"{population}_q_bh_six"] = "" if by is None else by\n'
        '                row[f"{population}_q_by_six"] = "" if bh is None else bh',
        JOINT_TEST + "::test_lineage_bh_by_columns_match_independent_family_oracle",
    ),
    "summary_effects_q_columns": (
        "f\"{row['before_q_bh_six']} / {row['after_q_bh_six']}\",\n"
        "            f\"{row['before_q_by_six']} / {row['after_q_by_six']}\",",
        "f\"{row['before_q_by_six']} / {row['after_q_by_six']}\",\n"
        "            f\"{row['before_q_bh_six']} / {row['after_q_bh_six']}\",",
        JOINT_TEST + "::test_summary_bh_by_cells_are_in_exact_columns",
    ),
    "summary_loso_q_columns": (
        "f\"{row['minimum_q_bh_six']} to {row['maximum_q_bh_six']} | \"\n"
        "            f\"{row['minimum_q_by_six']} to {row['maximum_q_by_six']} |\\n\"",
        "f\"{row['minimum_q_by_six']} to {row['maximum_q_by_six']} | \"\n"
        "            f\"{row['minimum_q_bh_six']} to {row['maximum_q_bh_six']} |\\n\"",
        JOINT_TEST + "::test_summary_bh_by_cells_are_in_exact_columns",
    ),
    "no_protocol_preflight": (
        'raise ValueError("Lineage-extension protocol hash differs from the declared rule")',
        "pass",
        "tests/test_influence_extension_core.py::"
        "test_protocol_rule_tampering_is_rejected_before_output_creation",
    ),
}


def reject_mutation(mutation: str, before: str, after: str, test_path: str) -> None:
    """Exercise one synthetic regression with copied source and the real test loader."""
    if Path(influence_extension.__file__).resolve() != SOURCE.resolve():
        raise AssertionError("Imported extension does not belong to this checkout")
    original = SOURCE.read_bytes()
    text = original.decode("utf-8")
    if text.count(before) != 1:
        raise ValueError(f"Expected exactly one mutation target: {mutation}")
    scratch = ROOT / ".tmp"
    scratch.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="lineage-mutation-", dir=scratch) as directory:
        copied = Path(directory) / "mtp"
        shutil.copytree(SOURCE.parent, copied, ignore=shutil.ignore_patterns("__pycache__"))
        mutated = copied / SOURCE.name
        mutated.write_text(text.replace(before, after), encoding="utf-8")
        spec = importlib.util.spec_from_file_location("mtp.influence_extension", mutated)
        if spec is None or spec.loader is None:
            raise AssertionError("Cannot load the isolated mutation")
        module = importlib.util.module_from_spec(spec)
        try:
            sys.modules["mtp.influence_extension"] = module
            mtp.influence_extension = module
            spec.loader.exec_module(module)
            code = pytest.main(["-q", "-p", "no:cacheprovider", str(ROOT / test_path)])
        finally:
            sys.modules["mtp.influence_extension"] = influence_extension
            mtp.influence_extension = influence_extension
    if SOURCE.read_bytes() != original:
        raise AssertionError("Tracked source changed during isolated mutation probe")
    if code != pytest.ExitCode.TESTS_FAILED:
        raise AssertionError(f"Mutation {mutation} was not rejected by its regression: {code}")
    print(f"Rejected {mutation}; tracked source bytes were never edited.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mutation", choices=sorted(MUTATIONS))
    mutation = parser.parse_args().mutation
    before, after, test = MUTATIONS[mutation]
    reject_mutation(mutation, before, after, test if "::" in test else TEST + "::" + test)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
