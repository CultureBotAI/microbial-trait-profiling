"""Shared study-level direction flip null generator and replicate RNG (#263, #182, #343).

Provides the shared exchangeability null where study-level direction labels are
randomly flipped (with 50% probability per study), preserving all intra-study
correlations, taxon co-annotations, signature sizes, and marginal distributions.

Re-exported by ``bsdbng.calibration`` so that downstream directional calibration
modules (such as ``screen_calibration.py`` in #319) and ORA calibration can import
from one uniform interface without module-versus-package collisions.
"""

from __future__ import annotations

import hashlib
from collections.abc import Hashable, Sequence
from dataclasses import replace
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from bsdbng.trait_guilds import Profile

SEED = 182
SHARED_FLIP_REPLICATES = 100


def _stream_key(name: str) -> int:
    """A stable integer stream id for a named scenario, independent of its position.

    ``hash`` is salted per interpreter, and a scenario's index would move whenever
    the tuple is reordered; a digest of the name survives both.
    """
    return int.from_bytes(hashlib.sha256(name.encode("utf-8")).digest()[:8], "big")


def replicate_rng(seed: int, stream: str, replicate: int) -> np.random.Generator:
    """The generator for one replicate of one named stream.

    ``numpy.random.default_rng`` takes the whole triple as seed entropy, so every
    (seed, stream, replicate) is independently reproducible and another arm can
    regenerate any single replicate without replaying the ones before it.
    """
    return np.random.default_rng([seed, _stream_key(stream), replicate])


def unit_keys(profiles: Sequence[Profile], unit: str = "study") -> list[Hashable]:
    """The resampling units of the corpus, in the order the flip draws them."""
    if unit not in ("study", "experiment"):
        raise ValueError(f"unknown flip unit {unit!r}")
    keys: set[Hashable] = (
        {p.study_id for p in profiles}
        if unit == "study"
        else {(p.study_id, p.experiment_id) for p in profiles}
    )
    return sorted(keys, key=repr)


def flipped_units(
    profiles: Sequence[Profile], rng: np.random.Generator, *, unit: str = "study"
) -> set[Hashable]:
    """Which units this replicate exchanges: one uniform draw per unit, in key order."""
    return {key for key in unit_keys(profiles, unit) if rng.random() < 0.5}


def flip_directions(
    profiles: Sequence[Profile], rng: np.random.Generator, *, unit: str = "study"
) -> list[Profile]:
    """Exchange the two direction labels within a random half of the units.

    ``unit`` is ``study`` or ``experiment``. Every profile of a flipped unit is
    relabelled, so the corpus stays paired and each unit keeps its own counts.
    """
    flipped = flipped_units(profiles, rng, unit=unit)
    swap = {"inc": "dec", "dec": "inc"}
    result: list[Profile] = []
    for profile in profiles:
        key = profile.study_id if unit == "study" else (profile.study_id, profile.experiment_id)
        result.append(
            replace(profile, inc_dec=swap[profile.inc_dec]) if key in flipped else profile
        )
    return result


def _unit_columns(key: Hashable, unit: str) -> tuple[str, str | None]:
    """Split a resampling unit into atomic study and experiment columns."""
    if unit == "study":
        return str(key), None
    if isinstance(key, tuple) and len(key) == 2:
        return str(key[0]), str(key[1])
    raise ValueError(f"unit {unit!r} does not describe {key!r}")


def flip_assignments(
    profiles: Sequence[Profile],
    *,
    unit: str = "study",
    replicates: int = SHARED_FLIP_REPLICATES,
    seed: int = SEED,
    stream: str = "complete_null_study_flip",
) -> list[dict[str, object]]:
    """The per-replicate flip set, as a table another arm can consume (#263).

    ``stream`` is the scenario name whose replicates these are, so replicate ``i``
    here is replicate ``i`` of that scenario: an arm that reruns its own pipeline on
    these assignments is measured on the identical null, not merely a similar one.
    """
    keys = unit_keys(profiles, unit)
    rows: list[dict[str, object]] = []
    for index in range(replicates):
        rng = replicate_rng(seed, stream, index)
        flipped = flipped_units(profiles, rng, unit=unit)
        for key in keys:
            study, experiment = _unit_columns(key, unit)
            rows.append(
                {
                    "seed": seed,
                    "stream": stream,
                    "unit": unit,
                    "replicate": index,
                    "study_id": study,
                    "experiment_id": experiment,
                    "flipped": int(key in flipped),
                }
            )
    return rows


def flip_profile_dataframe(
    profiles_df: pd.DataFrame,
    rng: np.random.Generator,
    *,
    unit: str = "study",
) -> pd.DataFrame:
    """Exchange direction labels within a random half of units in a DataFrame (#263).

    Supports both 'inc'/'dec' and 'increased'/'decreased' encodings.
    """
    if unit not in ("study", "experiment"):
        raise ValueError(f"unknown flip unit {unit!r}")
    if profiles_df.empty:
        return profiles_df.copy()

    df = profiles_df.copy()
    if unit == "study":
        units = sorted(df["study_id"].unique(), key=repr)
        flipped_keys = {u for u in units if rng.random() < 0.5}
        flip_mask = df["study_id"].isin(flipped_keys)
    else:
        row_keys = list(zip(df["study_id"], df["experiment_id"], strict=True))
        unique_exps = sorted(set(row_keys), key=repr)
        flipped_keys = {u for u in unique_exps if rng.random() < 0.5}
        flip_mask = pd.Series([k in flipped_keys for k in row_keys], index=df.index)

    swap_map = {
        "inc": "dec",
        "dec": "inc",
        "increased": "decreased",
        "decreased": "increased",
    }
    df.loc[flip_mask, "inc_dec"] = df.loc[flip_mask, "inc_dec"].map(
        lambda val: swap_map.get(str(val), str(val))
    )
    return df
