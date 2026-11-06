"""Explicit population selection for supplied figure tables (#1344).

No filename selects a population and no path-keyed cache is used. Callers keep
existing behavior by not invoking these helpers. A verified core already carries
its population and direction basis; its constructor preserves those rows without
reclassifying or relabelling them. These helpers do not verify a core themselves.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mtp import contrast_orientation as orientation
from mtp.human_host_corpus import (
    CORPUS_SEQUENCING_TYPES,
    HUMAN_HOST_SPECIES,
    experiment_in_profile_population,
)
from mtp.trait_guilds import Profile
from mtp.workflow.population import (
    ALL_PAIRED,
    NATIVE_DIRECTION_BASIS,
    ORIENTED_DIRECTION_BASIS,
    POPULATIONS,
)

ExperimentKey = tuple[str, str]
PROFILE_COLUMNS = ("study_id", "experiment_id", "inc_dec")
PAIR_COLUMNS = ("study_id", "experiment_id")
DIRECTIONS = {"inc": "inc", "dec": "dec", "increased": "inc", "decreased": "dec"}
FLIP = {"inc": "dec", "dec": "inc", "increased": "decreased", "decreased": "increased"}


def _population_basis(population: str, basis: str, *, verified_core: bool = False) -> str:
    if population not in POPULATIONS:
        raise ValueError(f"Unknown contrast population: {population!r}")
    if basis not in {NATIVE_DIRECTION_BASIS, ORIENTED_DIRECTION_BASIS}:
        raise ValueError(f"Unknown input direction basis: {basis!r}")
    output = NATIVE_DIRECTION_BASIS if population == ALL_PAIRED else ORIENTED_DIRECTION_BASIS
    if (population == ALL_PAIRED or verified_core) and basis != output:
        raise ValueError(
            f"Direction basis {basis!r} is incompatible with population {population!r}"
        )
    return output


@dataclass(frozen=True)
class PopulationSelection:
    """Declared membership and the direction changes still required for these inputs."""

    contrast_population: str
    input_direction_basis: str
    output_direction_basis: str
    selected_keys: frozenset[ExperimentKey]
    keys_to_flip: frozenset[ExperimentKey]

    def __post_init__(self) -> None:
        expected = _population_basis(self.contrast_population, self.input_direction_basis)
        if self.output_direction_basis != expected:
            raise ValueError("Selection output direction basis differs from its population")
        if not self.keys_to_flip <= self.selected_keys:
            raise ValueError("Direction changes contain an unselected experiment")
        if self.keys_to_flip and (
            self.contrast_population == ALL_PAIRED
            or self.input_direction_basis == ORIENTED_DIRECTION_BASIS
        ):
            raise ValueError("Identity or disease-oriented inputs cannot be flipped again")


def _columns(frame: pd.DataFrame, required: tuple[str, ...]) -> None:
    if not frame.columns.is_unique or not set(required) <= set(frame.columns):
        raise ValueError(f"Missing or duplicate identity/direction columns: {required}")


def _identity(value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("Missing experiment identity")
    return value


def _pairs(frame: pd.DataFrame) -> list[ExperimentKey]:
    _columns(frame, PAIR_COLUMNS)
    return [
        (_identity(study), _identity(experiment))
        for study, experiment in frame[list(PAIR_COLUMNS)].itertuples(index=False, name=None)
    ]


def _profile_pairs(frame: pd.DataFrame) -> list[ExperimentKey]:
    _columns(frame, PROFILE_COLUMNS)
    pairs = _pairs(frame)
    directions = frame["inc_dec"].tolist()
    if any(not isinstance(value, str) or value not in {"inc", "dec"} for value in directions):
        raise ValueError("Unknown profile direction; expected inc or dec")
    if len(set(zip(pairs, directions, strict=True))) != len(pairs):
        raise ValueError("Duplicate profile identity")
    return pairs


def selection_from_verified_core(
    profiles: pd.DataFrame, *, contrast_population: str, direction_basis: str
) -> PopulationSelection:
    """Preserve an already verified core's declared population and stored directions.

    Verification and binding of the supplied table to its core remain the caller's
    responsibility. No raw export or classification is consulted here.
    """
    output = _population_basis(contrast_population, direction_basis, verified_core=True)
    return PopulationSelection(
        contrast_population,
        direction_basis,
        output,
        frozenset(_profile_pairs(profiles)),
        frozenset(),
    )


def _metadata(path: Path) -> dict[ExperimentKey, tuple[object, object]]:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    _columns(frame, (*PAIR_COLUMNS, "host_species", "sequencing_type"))
    pairs = _pairs(frame)
    if len(set(pairs)) != len(pairs):
        raise ValueError("Duplicate experiment metadata identity")
    return {
        key: (host, method)
        for key, (host, method) in zip(
            pairs,
            frame[["host_species", "sequencing_type"]].itertuples(index=False, name=None),
            strict=True,
        )
    }


def _selected(found: orientation.Classification, population: str) -> bool:
    if found.contrast_class != orientation.DISEASE_VS_HEALTHY:
        return False
    if (
        population != "disease_vs_healthy_oriented"
        and found.contrast_type != orientation.CASE_VS_CONTROL_SUBJECT
    ):
        return False
    if population == "disease_vs_healthy_tier1_2_subject_controls" and found.control_tier not in {
        orientation.TIER_1_HEALTHY_NORMAL,
        orientation.TIER_2_CONDITION_NEGATIVE_SCREENED,
    }:
        return False
    if found.orientation not in {orientation.NATIVE, orientation.REVERSED}:
        raise ValueError("Selected experiment classification has an invalid orientation")
    return True


def get_population_selection(
    profiles_path: Path,
    experiments_path: Path,
    *,
    contrast_population: str,
    input_direction_basis: str,
    raw_dump_path: Path | None = None,
) -> PopulationSelection:
    """Select an explicitly requested population from standalone table identities.

    Counts are neither read nor altered; their scientific validation remains with
    the caller. ``all_paired`` retains workflow identity semantics, including
    unpaired and non-human rows, and does not read experiment/raw metadata.
    """
    output = _population_basis(contrast_population, input_direction_basis)
    profiles = pd.read_csv(
        profiles_path, sep="\t", usecols=list(PROFILE_COLUMNS), dtype=str, keep_default_na=False
    )
    pairs = _profile_pairs(profiles)
    if contrast_population == ALL_PAIRED:
        return PopulationSelection(
            contrast_population, input_direction_basis, output, frozenset(pairs), frozenset()
        )
    if raw_dump_path is None:
        raise ValueError("Explicit subset selection requires the raw export path")
    metadata = _metadata(experiments_path)
    classes = {
        key: orientation.classify(record)
        for key, record in orientation.read_experiment_records(raw_dump_path).items()
    }
    directions: dict[ExperimentKey, set[str]] = {}
    for key, direction in zip(pairs, profiles["inc_dec"], strict=True):
        if key not in metadata:
            raise ValueError(f"Missing experiment metadata for profile: {key}")
        if key not in classes:
            raise ValueError(f"Missing raw classification for profile: {key}")
        directions.setdefault(key, set()).add(direction)
    selected = frozenset(
        key
        for key, values in directions.items()
        if values == {"inc", "dec"}
        and experiment_in_profile_population(
            *metadata[key],
            host_species=HUMAN_HOST_SPECIES,
            sequencing_types=CORPUS_SEQUENCING_TYPES,
        )
        and _selected(classes[key], contrast_population)
    )
    if not selected:
        raise ValueError(f"Contrast population {contrast_population!r} selects no profiles")
    flips: frozenset[ExperimentKey] = frozenset()
    if input_direction_basis == NATIVE_DIRECTION_BASIS:
        probes = [
            Profile(study, experiment, "inc", "", "", {}, 0)
            for study, experiment in sorted(selected)
        ]
        oriented = orientation.oriented_profiles(probes, classes)
        if len(oriented) != len(probes):
            raise ValueError("Classifier dropped an explicitly selected experiment")
        flips = frozenset(
            (profile.study_id, profile.experiment_id)
            for profile in oriented
            if profile.inc_dec == "dec"
        )
    return PopulationSelection(contrast_population, input_direction_basis, output, selected, flips)


def _apply(
    frame: pd.DataFrame,
    pairs: list[ExperimentKey],
    columns: tuple[str, ...],
    selection: PopulationSelection,
) -> pd.DataFrame:
    if selection.contrast_population == ALL_PAIRED:
        return frame.copy()
    keep = [key in selection.selected_keys for key in pairs]
    result = frame.loc[keep].copy()
    flip = [
        key in selection.keys_to_flip for key, retained in zip(pairs, keep, strict=True) if retained
    ]
    if any(flip):
        for column in columns:
            result.loc[flip, column] = result.loc[flip, column].map(FLIP)
    return result


def select_profiles_pandas(frame: pd.DataFrame, selection: PopulationSelection) -> pd.DataFrame:
    """Select rows and, when required, relabel inc/dec without changing count values."""
    return _apply(frame, _profile_pairs(frame), ("inc_dec",), selection)


def select_taxa_pandas(frame: pd.DataFrame, selection: PopulationSelection) -> pd.DataFrame:
    """Preserve occurrence multiplicity/order and consistently relabel direction columns."""
    pairs = _pairs(frame)
    columns = tuple(name for name in ("direction", "inc_dec") if name in frame)
    if not columns:
        raise ValueError("Taxon table lacks a direction column")
    normalized = []
    for column in columns:
        values = frame[column].tolist()
        if any(not isinstance(value, str) or value not in DIRECTIONS for value in values):
            raise ValueError("Unknown taxon direction")
        normalized.append([DIRECTIONS[value] for value in values])
    if len(normalized) == 2 and normalized[0] != normalized[1]:
        raise ValueError("Taxon direction columns disagree")
    return _apply(frame, pairs, columns, selection)
