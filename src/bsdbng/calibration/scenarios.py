"""Simulation scenario generators for statistical calibration.

Signature sizes (hypergeometric draw counts) are read from the profile table's
``taxa_count`` column, which is written by ``bsdbng.study_trait_profiles``. A
``signature_size`` column is accepted as an alias. There is deliberately no
fallback value: a missing or malformed draw count raises ``ValueError``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PROFILE_METADATA_COLUMNS: frozenset[str] = frozenset(
    {
        "study_id",
        "experiment_id",
        "inc_dec",
        "condition",
        "body_site",
        "signature_size",
        "taxa_count",
    }
)
MATRIX_METADATA_COLUMNS: frozenset[str] = frozenset({"taxon_id", "taxon_name", "taxonomic_rank"})
SIGNATURE_SIZE_COLUMNS: tuple[str, ...] = ("taxa_count", "signature_size")
PROFILE_ID_COLUMNS: tuple[str, ...] = ("study_id", "experiment_id", "inc_dec")


def shared_trait_columns(profiles_df: pd.DataFrame, taxa_matrix_df: pd.DataFrame) -> list[str]:
    """Return trait columns present in both the profile table and the taxa matrix."""
    profile_cols = set(profiles_df.columns)
    return [
        col
        for col in taxa_matrix_df.columns
        if col not in MATRIX_METADATA_COLUMNS
        and col not in PROFILE_METADATA_COLUMNS
        and col in profile_cols
    ]


def signature_size_column(profiles_df: pd.DataFrame) -> str:
    """Return the name of the draw-count column, or raise if none is present."""
    for col in SIGNATURE_SIZE_COLUMNS:
        if col in profiles_df.columns:
            return col
    raise ValueError(
        "Profiles table has no signature size column; expected one of "
        f"{list(SIGNATURE_SIZE_COLUMNS)}. Regenerate profiles with `just study-trait-profiles`."
    )


def resolve_signature_sizes(profiles_df: pd.DataFrame) -> np.ndarray:
    """Return integer signature sizes per profile, validating every value."""
    col = signature_size_column(profiles_df)
    numeric = pd.to_numeric(profiles_df[col], errors="coerce")
    if numeric.isna().any():
        bad = int(numeric.isna().sum())
        raise ValueError(f"Column '{col}' has {bad} missing or non-numeric values")
    values = numeric.to_numpy(dtype=float)
    if not np.all(np.equal(np.mod(values, 1.0), 0.0)):
        raise ValueError(f"Column '{col}' contains non-integer values")
    sizes: np.ndarray = values.astype(np.int64)
    if np.any(sizes < 0):
        raise ValueError(f"Column '{col}' contains negative values")
    return sizes


def profile_identifiers(profiles_df: pd.DataFrame) -> np.ndarray:
    """Build ``study:experiment:direction`` identifiers for each profile row."""
    missing = [col for col in PROFILE_ID_COLUMNS if col not in profiles_df.columns]
    if missing:
        raise ValueError(f"Profiles table is missing identifier columns: {missing}")
    parts = [profiles_df[col].astype(str) for col in PROFILE_ID_COLUMNS]
    ids = parts[0].str.cat(parts[1:], sep=":")
    id_array: np.ndarray = ids.to_numpy(dtype=object)
    return id_array


def trait_population_counts(taxa_matrix_df: pd.DataFrame, trait_columns: list[str]) -> np.ndarray:
    """Count taxa carrying each trait (K) in the background universe."""
    numeric = taxa_matrix_df[trait_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    counts: np.ndarray = (numeric.to_numpy() > 0).sum(axis=0).astype(np.int64)
    return counts


def generate_complete_null_ora_profiles(
    empirical_profiles_df: pd.DataFrame,
    taxa_trait_matrix_df: pd.DataFrame,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Generate complete null profiles under the ORA null hypothesis (H0_ORA).

    Each profile's trait counts are drawn independently from
    Hypergeometric(M, K_trait, n_profile), where M is the number of taxa in the
    matrix, K_trait the number carrying the trait, and n_profile the empirical
    signature size.
    """
    trait_columns = shared_trait_columns(empirical_profiles_df, taxa_trait_matrix_df)
    if not trait_columns:
        raise ValueError("No trait columns shared between profiles and taxa matrix")
    k_pop_arr = trait_population_counts(taxa_trait_matrix_df, trait_columns)
    n_universe = len(taxa_trait_matrix_df)
    signature_sizes = resolve_signature_sizes(empirical_profiles_df)
    if signature_sizes.size and int(signature_sizes.max()) > n_universe:
        raise ValueError(
            f"Signature size {int(signature_sizes.max())} exceeds taxa universe size {n_universe}"
        )

    n_profiles = len(empirical_profiles_df)
    n_traits = len(trait_columns)
    simulated_k = np.zeros((n_profiles, n_traits), dtype=np.int64)
    ngood = k_pop_arr[np.newaxis, :]
    nbad = (n_universe - k_pop_arr)[np.newaxis, :]
    for n_draw in np.unique(signature_sizes):
        if n_draw == 0:
            continue
        p_mask = signature_sizes == n_draw
        n_sub = int(np.sum(p_mask))
        simulated_k[p_mask, :] = rng.hypergeometric(
            ngood=ngood, nbad=nbad, nsample=int(n_draw), size=(n_sub, n_traits)
        )

    size_col = signature_size_column(empirical_profiles_df)
    meta_cols = [
        col
        for col in empirical_profiles_df.columns
        if col in PROFILE_METADATA_COLUMNS and col != size_col
    ]
    sim_dict: dict[str, np.ndarray] = {
        col: empirical_profiles_df[col].to_numpy() for col in meta_cols
    }
    sim_dict[size_col] = signature_sizes
    for t_idx, col in enumerate(trait_columns):
        sim_dict[col] = simulated_k[:, t_idx]

    return pd.DataFrame(sim_dict)


def generate_mixed_signal_ora_profiles(
    empirical_profiles_df: pd.DataFrame,
    taxa_trait_matrix_df: pd.DataFrame,
    spike_traits: list[str],
    signal_strength: float = 0.5,
    rng: np.random.Generator | None = None,
) -> tuple[pd.DataFrame, set[tuple[str, str]]]:
    """Generate mixed signal profiles with known ground-truth spikes.

    Starts from a complete null draw, then adds ``round(n * signal_strength)``
    (at least 1) trait-bearing taxa to every testable profile for each spike
    trait, capped at the profile's signature size and the trait's population
    count. Only (profile, trait) pairs whose count actually increased are
    recorded as true signals.
    """
    if not 0.0 < signal_strength <= 1.0:
        raise ValueError(f"signal_strength must be in (0, 1], got {signal_strength}")
    if rng is None:
        rng = np.random.default_rng()

    null_df = generate_complete_null_ora_profiles(empirical_profiles_df, taxa_trait_matrix_df, rng)
    trait_columns = shared_trait_columns(null_df, taxa_trait_matrix_df)
    unknown = [t for t in spike_traits if t not in trait_columns]
    if unknown:
        raise ValueError(f"Spike traits not found in shared trait columns: {unknown}")

    profile_ids = profile_identifiers(null_df)
    signature_sizes = resolve_signature_sizes(null_df)
    k_pop = dict(
        zip(
            trait_columns,
            trait_population_counts(taxa_trait_matrix_df, trait_columns).tolist(),
            strict=True,
        )
    )

    true_signal_pairs: set[tuple[str, str]] = set()
    for trait in spike_traits:
        total_trait_taxa = int(k_pop[trait])
        if total_trait_taxa == 0:
            continue
        spike_counts = np.maximum(1, np.round(signature_sizes * signal_strength).astype(np.int64))
        current_k = null_df[trait].to_numpy(dtype=np.int64)
        cap = np.minimum(signature_sizes, total_trait_taxa)
        new_k = np.minimum(cap, current_k + spike_counts)
        new_k = np.where(signature_sizes > 0, new_k, current_k)
        null_df[trait] = new_k
        for p_id in profile_ids[new_k > current_k]:
            true_signal_pairs.add((str(p_id), trait))

    return null_df, true_signal_pairs
