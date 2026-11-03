"""Count the calibration injection's complement columns, for the #1345 review fixes.

Purpose: confirm how many retained trait columns no calibration procedure reads,
and that the injection headroom (``profile.total - sum(profile.guild_counts)``)
equals exactly the counts in those columns once the Table 1 and oxygen-guild term
sets are disjoint (#1237). Before #1237, METPO:1000608 was both a Table 1 trait and
an aerobe-guild member, so its counts were subtracted twice.

The tracked groups are the ones ``mtp.calibration.main`` loads:
``{**TABLE1_TRAITS, **OXYGEN_GUILDS}``.

Inputs: --profiles PATH, a profile matrix (the committed
        data/derived/study_trait_profiles_wide.tsv.gz).
Outputs: a JSON summary on stdout: group and distinct-term counts, pairwise
         overlaps, retained, tracked and complement column counts, and how many
         paired human profiles have a headroom different from their complement sum.
Issue: #1345 (review of PR #1338 for #1237). Date: 2026-11-02.
"""

from __future__ import annotations

import argparse
import csv
import json
from itertools import combinations
from pathlib import Path

from mtp.count_input import open_input_text
from mtp.trait_guilds import (
    OXYGEN_GUILDS,
    TABLE1_TRAITS,
    load_profiles,
    resolve_trait_columns,
    retained_trait_columns,
)


def summarize(profiles_path: Path, experiments_path: Path) -> dict[str, object]:
    groups = {**TABLE1_TRAITS, **OXYGEN_GUILDS}
    overlaps = {
        f"{left}&{right}": sorted(set(left_terms) & set(right_terms))
        for (left, left_terms), (right, right_terms) in combinations(groups.items(), 2)
        if set(left_terms) & set(right_terms)
    }
    distinct_terms = {term for terms in groups.values() for term in terms}
    with open_input_text(profiles_path) as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    _, retained, _ = retained_trait_columns(header)
    resolved = resolve_trait_columns(profiles_path, groups)
    tracked = {column for pairs in resolved.values() for _, column in pairs}
    complement = [column for column in retained if column not in tracked]
    profiles, _ = load_profiles(profiles_path, experiments_path, groups, reject_missing_counts=True)
    positions = {column: index for index, column in enumerate(header)}
    complement_sums: dict[tuple[str, str, str], int] = {}
    with open_input_text(profiles_path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            key = (row[0], row[1], row[2])
            complement_sums[key] = sum(int(row[positions[column]]) for column in complement)
    mismatched = [
        p
        for p in profiles
        if p.total - sum(p.guild_counts.values())
        != complement_sums[(p.study_id, p.experiment_id, p.inc_dec)]
    ]
    return {
        "issue": 1345,
        "groups": len(groups),
        "distinct_terms": len(distinct_terms),
        "pairwise_overlaps": overlaps,
        "retained_columns": len(retained),
        "tracked_columns": len(tracked),
        "complement_columns": len(complement),
        "paired_human_profiles": len(profiles),
        "profiles_headroom_differs_from_complement": len(mismatched),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument(
        "--experiments",
        type=Path,
        required=True,
        help="an experiments table covering the matrix, e.g. the #1237 canary's placeholder",
    )
    args = parser.parse_args()
    print(json.dumps(summarize(args.profiles, args.experiments), indent=2))


if __name__ == "__main__":
    main()
