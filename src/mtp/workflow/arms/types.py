"""Typed command plans shared by the three owner reruns (#658)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ArmContext:
    """Resolved core inputs and explicit settings for one fresh rerun."""

    repository: Path
    manifest: Path
    output: Path
    inputs: Path
    products: dict[str, Path]
    code_commit: str
    replicates: int | None
    seed: int | None
    figures: bool
    contrast_population: str = "all_paired"
    direction_basis: str = "group_1_vs_group_0"


@dataclass(frozen=True)
class Command:
    """An argument-vector invocation with required inputs and output contracts."""

    id: str
    argv: tuple[str, ...]
    inputs: tuple[Path, ...]
    outputs: tuple[Path, ...]
