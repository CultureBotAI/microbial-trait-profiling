"""Capture synthetic calibration parity for #784, 2026-10-22.

Input: the fixed synthetic paired profiles below and seed 1729. Output: JSON of
the replicate ledger, criteria and guild outcomes. The baseline was captured on
main 1ce8277 with --legacy-interface, before removing the ignored matrix argument.
This helper runs no production-data analysis and is retained with its fixture.
"""

from __future__ import annotations

import argparse
import inspect
import json
from collections.abc import Callable
from dataclasses import asdict
from pathlib import Path
from typing import Any

from mtp.calibration import Procedure, procedures_for, run_calibration
from mtp.guild_maps import build_matrix
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


def capture(*, legacy_interface: bool = False) -> dict[str, Any]:
    """The legacy path captures the old ignored argument; normal use has no matrix."""
    profiles = fixture_profiles()
    if legacy_interface:
        # Historical capture only: this eight-site matrix influenced no result.
        if "matrix" not in inspect.signature(procedures_for).parameters:
            raise ValueError("Legacy capture requires the pre-#784 calibration implementation")
        legacy_factory: Callable[..., dict[str, Procedure]] = procedures_for
        procedures = legacy_factory(
            build_matrix(profiles, max_conditions=60, max_body_sites=8, min_studies=3)
        )
    else:
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
    parser.add_argument("--legacy-interface", action="store_true")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            capture(legacy_interface=args.legacy_interface),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
