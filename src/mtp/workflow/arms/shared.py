"""Explicit cohort selection for legacy consumers of human count profiles (#658)."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

from mtp.paths import file_sha256
from mtp.workflow.arms.types import ArmContext, Command

HUMAN_HOST = "Homo sapiens"


def human_count_path(context: ArmContext) -> Path:
    """The execution-local human count table, retaining all feature columns."""
    return context.output / "human_counts" / "study_trait_profiles_wide.tsv"


def human_count_command(context: ArmContext) -> Command:
    """Record host selection as an executed, hashed prerequisite to legacy analyses."""
    profiles = context.products["count_profiles"]
    experiments = context.products["experiments"]
    output = human_count_path(context)
    summary = output.parent / "summary.json"
    return Command(
        "human_counts",
        (
            sys.executable,
            "-m",
            "mtp.workflow.arms.shared",
            "--profiles",
            str(profiles),
            "--experiments",
            str(experiments),
            "--output",
            str(output),
            "--summary",
            str(summary),
        ),
        (profiles, experiments),
        (output, summary),
    )


def select_human_count_profiles(
    profiles: Path, experiments: Path, output: Path, summary: Path
) -> None:
    """Keep every human profile verbatim in value and order; account for every exclusion."""
    hosts: dict[tuple[str, str], str] = {}
    with experiments.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"study_id", "experiment_id", "host_species"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Experiment metadata lacks study, experiment, or host fields")
        for metadata in reader:
            key = metadata["study_id"], metadata["experiment_id"]
            if key in hosts or not all(key):
                raise ValueError(f"Missing or duplicate experiment identity: {key}")
            if metadata["host_species"] is None:
                raise ValueError(f"Missing host metadata: {key}")
            hosts[key] = metadata["host_species"]

    retained = 0
    excluded: Counter[str] = Counter()
    with profiles.open(encoding="utf-8", newline="") as handle:
        profile_reader = csv.reader(handle, delimiter="\t")
        header = next(profile_reader, [])
        if not {"study_id", "experiment_id", "inc_dec"} <= set(header):
            raise ValueError("Count profiles lack study, experiment, or direction fields")
        study_index, experiment_index = header.index("study_id"), header.index("experiment_id")
        for row in profile_reader:
            if len(row) != len(header):
                raise ValueError("Count profile row does not match its header")
            key = row[study_index], row[experiment_index]
            if key not in hosts:
                raise ValueError(f"Count profile has no experiment metadata: {key}")
            if hosts[key] == HUMAN_HOST:
                retained += 1
            else:
                excluded[hosts[key] or "<unspecified>"] += 1
    if not retained:
        raise ValueError("Count profiles contain no Homo sapiens rows")
    destinations = {output.resolve(), summary.resolve()}
    if len(destinations) != 2 or destinations.intersection(
        {profiles.resolve(), experiments.resolve()}
    ):
        raise ValueError("Human count outputs must be distinct from each other and their sources")
    if any(path.exists() or path.is_symlink() for path in (output, summary)):
        raise FileExistsError("Human count outputs already exist")
    output.parent.mkdir(parents=True, exist_ok=True)
    with (
        profiles.open(encoding="utf-8", newline="") as source,
        output.open("w", encoding="utf-8", newline="") as handle,
    ):
        reader_rows = csv.reader(source, delimiter="\t")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(next(reader_rows))
        for row in reader_rows:
            if hosts[(row[study_index], row[experiment_index])] == HUMAN_HOST:
                writer.writerow(row)
    summary.parent.mkdir(parents=True, exist_ok=True)
    summary.write_text(
        json.dumps(
            {
                "host_species": HUMAN_HOST,
                "selection": "exact experiments.host_species equality; no other row/feature filter",
                "profiles": {"path": str(profiles), "sha256": file_sha256(profiles)},
                "experiments": {"path": str(experiments), "sha256": file_sha256(experiments)},
                "output": {"path": str(output), "sha256": file_sha256(output)},
                "input_rows": retained + excluded.total(),
                "retained_rows": retained,
                "excluded_rows": excluded.total(),
                "excluded_rows_by_host": dict(sorted(excluded.items())),
                "feature_and_identifier_columns": len(header),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def main(argv: list[str] | None = None) -> None:
    """Materialize a human cohort from explicit verified-core file paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--experiments", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args(argv)
    select_human_count_profiles(args.profiles, args.experiments, args.output, args.summary)


if __name__ == "__main__":
    main()
