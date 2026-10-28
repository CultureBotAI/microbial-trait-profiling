"""Validate and select matrix-coverage inputs without silent loss (#1117)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

EXPERIMENT_KEYS = ["study_id", "experiment_id"]
PROFILE_KEYS = [*EXPERIMENT_KEYS, "inc_dec"]
MATRIX_METADATA = {"taxon_id", "taxon_name", "taxonomic_rank"}


def identifiers(frame: pd.DataFrame, keys: list[str], description: str, *, unique: bool) -> None:
    if not frame.columns.is_unique or not set(keys).issubset(frame.columns):
        raise ValueError(f"Missing or duplicate {description} identity columns")
    for key in keys:
        if any(
            not isinstance(value, str)
            or not value
            or value.strip() != value
            or any(character in value for character in "\t\r\n")
            for value in frame[key]
        ):
            raise ValueError(f"Invalid {description} identifier: {key}")
    if unique and frame.duplicated(keys).any():
        raise ValueError(f"Duplicate {description} identity")


def normalize_direction(value: object) -> str:
    if isinstance(value, str):
        raw = value.strip().lower()
        if raw in {"inc", "increased", "up"}:
            return "inc"
        if raw in {"dec", "decreased", "down"}:
            return "dec"
    raise ValueError(f"Unknown profile direction: {value!r}")


def validated_matrix(matrix: pd.DataFrame, columns: pd.DataFrame) -> pd.DataFrame:
    identifiers(matrix, ["taxon_id"], "matrix taxon", unique=True)
    identifiers(columns, ["column_id"], "feature", unique=True)
    if not set(columns.column_id).issubset(matrix.columns):
        raise ValueError("Feature metadata column identity is absent from the matrix")
    result = matrix.copy()
    for name in matrix.columns:
        if name in MATRIX_METADATA:
            continue
        if any(isinstance(value, (bool, np.bool_)) for value in matrix[name]):
            raise ValueError("Matrix cells must contain finite nonnegative numeric values")
        try:
            values = pd.to_numeric(matrix[name], errors="raise")
            numeric = values.to_numpy(dtype=float)
        except (ValueError, TypeError) as error:
            raise ValueError("Matrix contains non-numeric values") from error
        if not np.isfinite(numeric).all() or (numeric < 0).any():
            raise ValueError(
                "Matrix contains negative or non-numeric values (including nonfinite values)"
            )
        result[name] = values
    return result


def _keys(frame: pd.DataFrame, keys: list[str]) -> pd.MultiIndex:
    return pd.MultiIndex.from_frame(frame[keys])


def _covered(frame: pd.DataFrame, metadata: pd.DataFrame, description: str) -> None:
    if not _keys(frame, EXPERIMENT_KEYS).isin(_keys(metadata, EXPERIMENT_KEYS)).all():
        raise ValueError(f"{description} contains an orphan experiment without metadata")


@dataclass(frozen=True)
class PopulationInputs:
    active_keys: pd.DataFrame
    active_experiments: pd.DataFrame
    active_taxa: pd.DataFrame
    source_experiments: pd.DataFrame
    source_taxa: pd.DataFrame
    input_occurrence_rows: int
    excluded_occurrence_rows: int
    excluded_profile_keys: int


def population_inputs(
    tier12_experiments: pd.DataFrame,
    tier12_taxa: pd.DataFrame,
    tier12_profile_keys: pd.DataFrame,
    source_experiments: pd.DataFrame,
    source_taxa: pd.DataFrame,
) -> PopulationInputs:
    """Reject malformed joins; report valid occurrences outside the active keys."""
    for frame in (tier12_experiments, source_experiments):
        identifiers(frame, EXPERIMENT_KEYS, "experiment metadata", unique=True)
        if not {"body_site", "condition", "host_species"}.issubset(frame.columns):
            raise ValueError("Experiment metadata lacks body site, condition or host")
    prepared = []
    for frame in (tier12_taxa, source_taxa):
        identifiers(frame, [*EXPERIMENT_KEYS, "taxon_id"], "reported taxon", unique=False)
        if "direction" not in frame:
            raise ValueError("Reported taxa lack a direction")
        frame = frame.copy()
        frame["inc_dec"] = frame.direction.map(normalize_direction)
        prepared.append(frame)
    tier_taxa, full_taxa = prepared
    # Valid non-active taxa may only have metadata in the full source export.
    _covered(tier_taxa, source_experiments, "Tier 1+2 reported taxa")
    _covered(full_taxa, source_experiments, "Source reported taxa")
    identifiers(tier12_profile_keys, EXPERIMENT_KEYS, "profile", unique=False)
    if "inc_dec" not in tier12_profile_keys:
        raise ValueError("Profile keys lack a direction")
    active = tier12_profile_keys[PROFILE_KEYS].copy()
    active["inc_dec"] = active.inc_dec.map(normalize_direction)
    identifiers(active, PROFILE_KEYS, "profile", unique=True)
    _covered(active, tier12_experiments, "Active profile")
    _covered(active, source_experiments, "Active profile source")
    selected = _keys(tier_taxa, PROFILE_KEYS).isin(_keys(active, PROFILE_KEYS))
    excluded = tier_taxa.loc[~selected]
    selected_experiments = _keys(tier12_experiments, EXPERIMENT_KEYS).isin(
        _keys(active, EXPERIMENT_KEYS)
    )
    return PopulationInputs(
        active,
        tier12_experiments.loc[selected_experiments].copy(),
        tier_taxa.loc[selected].copy(),
        source_experiments.copy(),
        full_taxa,
        len(tier_taxa),
        len(excluded),
        len(excluded[PROFILE_KEYS].drop_duplicates()),
    )


def select_experiments(frame: pd.DataFrame, experiments: pd.DataFrame) -> pd.DataFrame:
    """Select by the complete experiment identity, never by experiment ID alone."""
    return frame.loc[_keys(frame, EXPERIMENT_KEYS).isin(_keys(experiments, EXPERIMENT_KEYS))].copy()
