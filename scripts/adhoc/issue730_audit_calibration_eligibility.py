"""Audit calibration/production guild eligibility for #730 (2026-10-19).

Inputs: explicit count-profile and experiment TSVs. Output: JSON with observed
zero-total exclusions, raw and production guild tests, and calibration parity.
Historical source files are read only. Retain this script with its results.

The unused eight-body-site display matrix was removed for #784 on 2026-10-22.
The original audit's McNemar and rank correlations never read that matrix;
their scientific outputs are unchanged (producer source hashes change).
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from mtp.calibration import guild_mcnemar_with
from mtp.guild_maps import compute_correlation
from mtp.paths import file_sha256
from mtp.trait_guilds import OXYGEN_GUILDS, TABLE1_TRAITS, drop_zero_total_pairs, load_profiles


def audit(profiles_path: Path, experiments_path: Path) -> dict[str, object]:
    profiles, stats = load_profiles(
        profiles_path,
        experiments_path,
        {**TABLE1_TRAITS, **OXYGEN_GUILDS},
        reject_missing_counts=True,
    )
    exclusion = drop_zero_total_pairs(profiles)
    raw = compute_correlation(profiles)
    production = compute_correlation(exclusion.profiles)
    calibrated = guild_mcnemar_with()(profiles)[0]
    return {
        "issue": 730,
        "inputs": {
            str(path.resolve()): file_sha256(path) for path in (profiles_path, experiments_path)
        },
        "source_sha256": {
            str(path): file_sha256(path)
            for path in (
                Path(__file__),
                Path(__file__).parents[2] / "src/mtp/calibration.py",
                Path(__file__).parents[2] / "src/mtp/guild_maps.py",
            )
        },
        "corpus": asdict(stats),
        "eligibility": {
            "zero_total_experiments": exclusion.zero_total_experiments,
            "studies_before": exclusion.studies_before,
            "studies_after": exclusion.studies_after,
            "studies_lost": exclusion.studies_lost,
            "profiles_after": len(exclusion.profiles),
            "null_invariance": "Direction swaps preserve the unordered pair of retained totals, "
            "so these same zero-total experiment pairs occur in every complete-null replicate.",
        },
        "unfiltered_guild": {"p": raw.mcnemar_p, "discordant": list(raw.discordant)},
        "production_guild": {
            "p": production.mcnemar_p,
            "discordant": list(production.discordant),
        },
        "calibration_guild": asdict(calibrated),
        "calibration_matches_production": (
            calibrated.p_value == production.mcnemar_p
            and calibrated.scored == sum(production.discordant)
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--experiments", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for path in (args.profiles, args.experiments, Path(__file__)):
        if args.output.resolve() == path.resolve() or (
            args.output.exists() and path.samefile(args.output)
        ):
            parser.error("Audit output must not overwrite an input or the audit script")
    result = audit(args.profiles, args.experiments)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
