"""Capture synthetic calibration parity for #784, 2026-10-22.

Input: the fixed synthetic paired profiles below and seed 1729. Output: JSON of
the replicate ledger, criteria and guild outcomes. The baseline was captured on
main 1ce8277 with the historical helper's --legacy-interface option, before
removing the ignored matrix argument. Review #900 removes that obsolete capture
path; this helper now exercises only the supported interface. No production data.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from mtp.calibration import procedures_for, run_calibration
from mtp.trait_guilds import TABLE1_TRAITS, Profile

SEED = 1729


def fixture_profiles() -> list[Profile]:
    """Include repeated experiments and variable, nontrivial oxygen-guild signs."""
    profiles = []
    for study in range(12):
        for experiment in range(1 + study % 2):
            for direction, offset in (("inc", 3), ("dec", 1)):
                counts = {
                    trait: (study * 3 + index * 5 + offset) % 17
                    for index, trait in enumerate(TABLE1_TRAITS)
                }
                counts["aerobe"] = (study + offset * (study % 3 + 1)) % 13 + 1
                counts["anaerobe"] = (study * 2 - offset * (study % 4 + 1)) % 11 + 1
                profiles.append(
                    Profile(
                        study_id=f"s{study}",
                        experiment_id=f"s{study}-{experiment}",
                        inc_dec=direction,
                        condition=f"condition{study % 3}",
                        body_site=f"site{study % 10}",
                        guild_counts=counts,
                        total=sum(counts.values()) + 400,
                    )
                )
    return profiles


def capture() -> dict[str, Any]:
    """Capture scientific parity using the supported matrix-free interface."""
    profiles = fixture_profiles()
    procedures = procedures_for()
    result = run_calibration(profiles, procedures, 4, 2, SEED)
    return {
        "seed": SEED,
        "studies": 12,
        "null_replicates": 4,
        "power_replicates": 2,
        "ledger": result.ledger,
        "criteria": result.criteria,
        "guild_outcomes": [asdict(outcome) for outcome in procedures["guild_mcnemar"](profiles)],
        "completed": result.completed,
    }


def main() -> None:
    """Save a comparison snapshot; fixture size is intentionally below calibration budgets."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            capture(),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
