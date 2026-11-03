"""Measure the oxygen guilds on the committed profile matrix, for the #1237 canary.

Purpose: compare the aerobe and anaerobe guild totals and the study-level guild
McNemar before and after METPO:1000608 left the aerobe guild (#1237). Run
``measure`` once with the base commit's ``src`` first on PYTHONPATH and once on
the branch, then pass both output directories to ``compare``.

The production experiments table and BugSigDB dump are not committed, so the
condition and body site of each experiment are unknown here. ``measure`` writes a
placeholder experiments table giving every experiment of the matrix one condition,
one body site and a human host (the committed matrix holds human profiles only,
#656). The real loader then applies its pair filter, the shared zero-total
exclusion and the guild sums, and the guild maps' own McNemar runs on the result.
Cell statuses need the real metadata and are not measured.

Inputs:
  measure --profiles PATH --output-dir DIR
      PATH is the committed data/derived/study_trait_profiles_wide.tsv.gz;
      DIR is outside the repository.
  compare --base DIR --branch DIR --output PATH
Outputs:
  DIR/placeholder_experiments.tsv, DIR/study_guild_contrasts.tsv and
  DIR/guild_summary.json from ``measure``; the JSON at PATH from ``compare``.
Issue: #1237. Date: 2026-11-02.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import mtp
from mtp.count_input import open_input_text
from mtp.guild_maps import SIGN_TOLERANCE, compute_correlation, study_direction_means
from mtp.paths import file_sha256
from mtp.trait_guilds import (
    OXYGEN_GUILDS,
    drop_zero_total_pairs,
    load_profiles,
    resolve_trait_columns,
)

PLACEHOLDER = "unassigned (#1237 canary)"


def write_placeholder_experiments(profiles: Path, output: Path) -> int:
    """One human row per (study, experiment) of the matrix, condition unknown."""
    keys: set[tuple[str, str]] = set()
    with open_input_text(profiles) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            keys.add((row["study_id"], row["experiment_id"]))
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["study_id", "experiment_id", "condition", "body_site", "host_species"])
        for study_id, experiment_id in sorted(keys):
            writer.writerow([study_id, experiment_id, PLACEHOLDER, PLACEHOLDER, "Homo sapiens"])
    return len(keys)


def measure(profiles_path: Path, output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    experiments_path = output_dir / "placeholder_experiments.tsv"
    experiments = write_placeholder_experiments(profiles_path, experiments_path)
    profiles, stats = load_profiles(
        profiles_path, experiments_path, OXYGEN_GUILDS, reject_missing_counts=True
    )
    exclusion = drop_zero_total_pairs(profiles)
    kept = exclusion.profiles
    correlation = compute_correlation(kept)
    means = study_direction_means(kept)
    with (output_dir / "study_guild_contrasts.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["study_id", "aerobe_contrast", "anaerobe_contrast"])
        for study, per_key in sorted(means.items()):
            writer.writerow(
                [
                    study,
                    repr(per_key["aerobe_inc"] - per_key["aerobe_dec"]),
                    repr(per_key["anaerobe_inc"] - per_key["anaerobe_dec"]),
                ]
            )
    summary: dict[str, object] = {
        "issue": 1237,
        "mtp_module": str(Path(mtp.__file__).resolve().parent),
        "oxygen_guilds": {guild: list(terms) for guild, terms in OXYGEN_GUILDS.items()},
        "resolved_columns": {
            guild: [column for _, column in pairs]
            for guild, pairs in resolve_trait_columns(profiles_path, OXYGEN_GUILDS).items()
        },
        "profiles_sha256": file_sha256(profiles_path),
        "placeholder_experiments": experiments,
        "kept_profiles_before_zero_total": stats.kept_profiles,
        "zero_total_experiments": exclusion.zero_total_experiments,
        "studies": exclusion.studies_after,
        "profiles": len(kept),
        "guild_count_mass": {
            guild: sum(p.guild_counts[guild] for p in kept) for guild in OXYGEN_GUILDS
        },
        "profiles_with_guild_signal": {
            guild: sum(1 for p in kept if p.guild_counts[guild] > 0) for guild in OXYGEN_GUILDS
        },
        "total_count_mass": sum(p.total for p in kept),
        "sign_table": correlation.table,
        "mcnemar_discordant": list(correlation.discordant),
        "mcnemar_p": correlation.mcnemar_p,
        "odds_ratio": correlation.odds_ratio,
        "mcnemar_status": correlation.status,
        "studies_scored": correlation.studies,
        "studies_without_oxygen_annotation": correlation.excluded_no_annotation,
        "spearman": {key: list(value) for key, value in correlation.spearman.items()},
    }
    (output_dir / "guild_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def _contrasts(path: Path) -> dict[str, tuple[float, float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            row["study_id"]: (float(row["aerobe_contrast"]), float(row["anaerobe_contrast"]))
            for row in csv.DictReader(handle, delimiter="\t")
        }


def _sign(value: float) -> int:
    """The guild maps' own study sign, tolerance included, so the counts reconcile."""
    if abs(value) < SIGN_TOLERANCE:
        return 0
    return 1 if value > 0 else -1


def compare(base_dir: Path, branch_dir: Path, output: Path) -> dict[str, object]:
    base = json.loads((base_dir / "guild_summary.json").read_text(encoding="utf-8"))
    branch = json.loads((branch_dir / "guild_summary.json").read_text(encoding="utf-8"))
    base_contrasts = _contrasts(base_dir / "study_guild_contrasts.tsv")
    branch_contrasts = _contrasts(branch_dir / "study_guild_contrasts.tsv")
    if set(base_contrasts) != set(branch_contrasts):
        raise ValueError("The two runs scored different studies")
    flips: dict[str, int] = {}
    for study, (aerobe, _) in base_contrasts.items():
        key = f"{_sign(aerobe):+d}->{_sign(branch_contrasts[study][0]):+d}"
        flips[key] = flips.get(key, 0) + 1
    anaerobe_changed = sum(
        1
        for study, (_, anaerobe) in base_contrasts.items()
        if _sign(anaerobe) != _sign(branch_contrasts[study][1])
    )
    base_mass = base["guild_count_mass"]["aerobe"]
    branch_mass = branch["guild_count_mass"]["aerobe"]
    result: dict[str, object] = {
        "issue": 1237,
        "same_matrix": base["profiles_sha256"] == branch["profiles_sha256"],
        "same_population": all(
            base[key] == branch[key] for key in ("studies", "profiles", "total_count_mass")
        ),
        "aerobe_columns": {
            "base": base["resolved_columns"]["aerobe"],
            "branch": branch["resolved_columns"]["aerobe"],
        },
        "aerobe_count_mass": {
            "base": base_mass,
            "branch": branch_mass,
            "relative_change": (branch_mass - base_mass) / base_mass,
        },
        "anaerobe_count_mass": {
            "base": base["guild_count_mass"]["anaerobe"],
            "branch": branch["guild_count_mass"]["anaerobe"],
        },
        "profiles_with_aerobe_signal": {
            "base": base["profiles_with_guild_signal"]["aerobe"],
            "branch": branch["profiles_with_guild_signal"]["aerobe"],
        },
        "aerobe_study_sign_transitions": dict(sorted(flips.items())),
        "anaerobe_study_signs_changed": anaerobe_changed,
        "sign_table": {"base": base["sign_table"], "branch": branch["sign_table"]},
        "mcnemar": {
            run: {
                "discordant": summary["mcnemar_discordant"],
                "p": summary["mcnemar_p"],
                "odds_ratio": summary["odds_ratio"],
                "status": summary["mcnemar_status"],
                "studies_scored": summary["studies_scored"],
            }
            for run, summary in (("base", base), ("branch", branch))
        },
        "contrast_spearman_aerobe_vs_anaerobe": {
            run: summary["spearman"].get("aerobe_contrast~anaerobe_contrast")
            for run, summary in (("base", base), ("branch", branch))
        },
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("measure")
    run.add_argument("--profiles", type=Path, required=True)
    run.add_argument("--output-dir", type=Path, required=True)
    both = sub.add_parser("compare")
    both.add_argument("--base", type=Path, required=True)
    both.add_argument("--branch", type=Path, required=True)
    both.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "measure":
        summary = measure(args.profiles, args.output_dir)
        print(json.dumps({key: summary[key] for key in ("mtp_module", "oxygen_guilds")}))
    else:
        print(json.dumps(compare(args.base, args.branch, args.output), indent=2))


if __name__ == "__main__":
    main()
