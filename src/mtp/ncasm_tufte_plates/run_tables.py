"""Read the run tables an NCASM plate is drawn from, with their hashes (#2264, #776).

A plate that reads tables takes a run root: the folder of one analysis run holding
``comparison/`` and ``arms/``, for example ``.../F-511c891/runs/tier12-off``. Every
loader here checks that the tables describe one population and one test family, and
returns the files it read with their SHA-256, which the build writes into the manifest.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from mtp.ncasm_tufte_plates.data_status import InputFile

__all__ = [
    "ARM1_POOLED_STATS",
    "COMPARISON_RUN",
    "DIRECTION_SCREEN",
    "RunIdentity",
    "RunTableError",
    "TraitCloud",
    "TraitPoint",
    "load_trait_cloud",
]

COMPARISON_RUN = Path("comparison/comparison_run.json")
ARM1_POOLED_STATS = Path(
    "comparison/arm1_pooled_stats/trait_profile_paired_contrast_trait_stats_pooled.tsv"
)
DIRECTION_SCREEN = Path("arms/noel/trait_direction_screen/trait_direction_screen.tsv")
SIGNIFICANCE_LEVEL = 0.05


class RunTableError(ValueError):
    """A run table is missing, or does not describe what the plate assumes."""


@dataclass(frozen=True)
class RunIdentity:
    """Which run the tables belong to, from ``comparison/comparison_run.json``."""

    comparison_id: str
    contrast_population: str
    direction_basis: str
    core_run_id: str
    git_commit: str
    git_dirty: bool


@dataclass(frozen=True)
class TraitPoint:
    """One trait tested by both arms."""

    trait: str
    name: str
    mean_shift: float
    studies_up: int
    studies_down: int
    share_up: float
    by_q: float
    in_table1: bool

    @property
    def percent_up(self) -> float:
        return 100.0 * self.share_up

    @property
    def passes_by(self) -> bool:
        return self.by_q <= SIGNIFICANCE_LEVEL


@dataclass(frozen=True)
class TraitCloud:
    """Per-trait pooled shift (Arm 1) against the study sign screen (Arm 3)."""

    identity: RunIdentity
    points: tuple[TraitPoint, ...]
    comparisons: int
    arm1_studies: int
    screen_traits: int
    screen_studies_max: int
    screen_by_hits: tuple[str, ...]
    correlation: float
    inputs: tuple[InputFile, ...]

    @property
    def table1(self) -> tuple[TraitPoint, ...]:
        return tuple(point for point in self.points if point.in_table1)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _input(run_root: Path, relative: Path) -> InputFile:
    path = run_root / relative
    if not path.is_file():
        raise RunTableError(f"{relative} is missing under the run root {run_root}")
    return {"path": relative.as_posix(), "sha256": _sha256(path)}


def _rows(run_root: Path, relative: Path, required: tuple[str, ...]) -> list[dict[str, str]]:
    with (run_root / relative).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [column for column in required if column not in (reader.fieldnames or [])]
        if missing:
            raise RunTableError(f"{relative} lacks the columns {', '.join(missing)}")
        rows = list(reader)
    if not rows:
        raise RunTableError(f"{relative} has no rows")
    return rows


def _one(rows: list[dict[str, str]], column: str, relative: Path) -> str:
    values = {row[column] for row in rows}
    if len(values) != 1:
        raise RunTableError(f"{relative}: {column} is not constant ({sorted(values)[:4]})")
    return values.pop()


def _number(row: dict[str, str], column: str, relative: Path) -> float:
    try:
        value = float(row[column])
    except ValueError as error:
        raise RunTableError(f"{relative}: {column} of {row['trait']} is not a number") from error
    if not math.isfinite(value):
        raise RunTableError(f"{relative}: {column} of {row['trait']} is not finite")
    return value


def _count(row: dict[str, str], column: str, relative: Path) -> int:
    value = _number(row, column, relative)
    if value != int(value) or value < 0:
        raise RunTableError(f"{relative}: {column} of {row['trait']} is not a count")
    return int(value)


def _display_name(object_name: str) -> str:
    name = object_name.replace("_", " ").strip()
    return name[:1].upper() + name[1:]


def _identity(run_root: Path) -> RunIdentity:
    record = json.loads((run_root / COMPARISON_RUN).read_text(encoding="utf-8"))
    if record.get("status") != "completed":
        raise RunTableError(f"{COMPARISON_RUN}: the comparison is not completed")
    try:
        return RunIdentity(
            comparison_id=str(record["id"]),
            contrast_population=str(record["contrast_population"]),
            direction_basis=str(record["direction_basis"]),
            core_run_id=str(record["core_run_id"]),
            git_commit=str(record["git_commit"]),
            git_dirty=bool(record["git_dirty"]),
        )
    except KeyError as error:
        raise RunTableError(f"{COMPARISON_RUN} lacks {error}") from error


def load_trait_cloud(run_root: Path) -> TraitCloud:
    """The traits tested by both the pooled comparison shift and the study sign screen.

    Refuses tables that are stratified, untested or inconsistent, and a pooled-shift
    trait that the screen does not list, rather than drawing a partial cloud."""
    inputs = tuple(
        _input(run_root, relative)
        for relative in (COMPARISON_RUN, ARM1_POOLED_STATS, DIRECTION_SCREEN)
    )
    identity = _identity(run_root)

    arm1 = _rows(
        run_root,
        ARM1_POOLED_STATS,
        (
            "trait",
            "stratum_kind",
            "pair_count",
            "study_count",
            "pair_test_status",
            "mean_increased_minus_decreased",
        ),
    )
    if _one(arm1, "stratum_kind", ARM1_POOLED_STATS) != "pooled":
        raise RunTableError(f"{ARM1_POOLED_STATS} is not the pooled stratum")
    if _one(arm1, "pair_test_status", ARM1_POOLED_STATS) != "tested":
        raise RunTableError(f"{ARM1_POOLED_STATS} holds untested traits")
    comparisons = int(_one(arm1, "pair_count", ARM1_POOLED_STATS))
    arm1_studies = int(_one(arm1, "study_count", ARM1_POOLED_STATS))

    screen = _rows(
        run_root,
        DIRECTION_SCREEN,
        (
            "trait",
            "object_name",
            "in_table1",
            "studies_up",
            "studies_down",
            "studies_scored",
            "share_up",
            "by_q_nominal",
            "test_status",
            "family_size_nominal",
        ),
    )
    if _one(screen, "test_status", DIRECTION_SCREEN) != "tested":
        raise RunTableError(f"{DIRECTION_SCREEN} holds untested traits")
    if int(_one(screen, "family_size_nominal", DIRECTION_SCREEN)) != len(screen):
        raise RunTableError(f"{DIRECTION_SCREEN}: the family size is not the number of rows")
    by_trait = {row["trait"]: row for row in screen}
    if len(by_trait) != len(screen) or len({row["trait"] for row in arm1}) != len(arm1):
        raise RunTableError("a trait is listed twice in a run table")
    absent = sorted(row["trait"] for row in arm1 if row["trait"] not in by_trait)
    if absent:
        raise RunTableError(f"{DIRECTION_SCREEN} lacks pooled-shift traits: {absent[:4]}")

    points = []
    for row in arm1:
        vote = by_trait[row["trait"]]
        flag = vote["in_table1"].lower()
        if flag not in {"true", "false"}:
            raise RunTableError(f"{DIRECTION_SCREEN}: in_table1 of {row['trait']} is {flag!r}")
        points.append(
            TraitPoint(
                trait=row["trait"],
                name=_display_name(vote["object_name"]),
                mean_shift=_number(row, "mean_increased_minus_decreased", ARM1_POOLED_STATS),
                studies_up=_count(vote, "studies_up", DIRECTION_SCREEN),
                studies_down=_count(vote, "studies_down", DIRECTION_SCREEN),
                share_up=_number(vote, "share_up", DIRECTION_SCREEN),
                by_q=_number(vote, "by_q_nominal", DIRECTION_SCREEN),
                in_table1=flag == "true",
            )
        )
    by_hits = tuple(
        _display_name(row["object_name"])
        for row in screen
        if _number(row, "by_q_nominal", DIRECTION_SCREEN) <= SIGNIFICANCE_LEVEL
    )
    correlation = float(
        np.corrcoef([point.mean_shift for point in points], [point.share_up for point in points])[
            0, 1
        ]
    )
    return TraitCloud(
        identity=identity,
        points=tuple(points),
        comparisons=comparisons,
        arm1_studies=arm1_studies,
        screen_traits=len(screen),
        screen_studies_max=max(_count(row, "studies_scored", DIRECTION_SCREEN) for row in screen),
        screen_by_hits=by_hits,
        correlation=correlation,
        inputs=inputs,
    )
