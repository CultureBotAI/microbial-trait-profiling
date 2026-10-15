"""Simulation scenarios and data generation for statistical calibration (#348, #351, #343).

Generates synthetic null and mixed-signal profiles under:
1. Uniform hypergeometric urn draws without replacement (H0_ORA benchmark check)
2. Shared study-level direction flip permutation null (#263 / #343)
3. Mixed signal profiles with known ground-truth spikes for empirical FDR assessment

All matrix inputs are strictly binarized to presence (> 0) so that observed trait counts
k never exceed draw counts N or population counts K (resolving #348).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


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


PROFILE_METADATA_COLUMNS: list[str] = [
    "study_id",
    "experiment_id",
    "condition",
    "body_site",
    "antibiotics_exclusion",
    "host_species",
    "sequencing_type",
    "sample_type",
    "inc_dec",
    "taxa_count",
    "n_taxa",
    "signature_size",
]

TAXA_METADATA_COLUMNS: list[str] = [
    "taxon_name",
    "taxonomic_rank",
    "ncbi_taxonomy_id",
    "rank",
    "name",
]

PROFILE_ID_COLUMNS: list[str] = ["study_id", "experiment_id", "inc_dec"]


def signature_size_column(profiles_df: pd.DataFrame) -> str:
    """Find the column specifying signature size N per profile."""
    for col in ("taxa_count", "n_taxa", "signature_size"):
        if col in profiles_df.columns:
            return col
    raise ValueError(
        "Profiles table missing signature size column ('taxa_count', 'n_taxa', or 'signature_size')"
    )


def shared_trait_columns(profiles_df: pd.DataFrame, taxa_matrix_df: pd.DataFrame) -> list[str]:
    """Find trait columns present in both profiles and taxa matrix."""
    profile_traits = {
        col for col in profiles_df.columns if col not in set(PROFILE_METADATA_COLUMNS)
    }
    taxa_traits = {col for col in taxa_matrix_df.columns if col not in set(TAXA_METADATA_COLUMNS)}
    shared = sorted(profile_traits & taxa_traits)
    return shared


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
    ids = parts[0]
    for part in parts[1:]:
        ids = ids.str.cat(part, sep=":")
    id_array: np.ndarray = ids.to_numpy(dtype=object)
    return id_array


def _parse_strict_count(val: Any) -> int:
    """Parse count strictly matching production rules."""
    if val is None or val == "" or (isinstance(val, float) and np.isnan(val)):
        return 0
    if isinstance(val, int | np.integer):
        if val < 0:
            raise ValueError(f"Negative count: {val}")
        return int(val)
    val_str = str(val).strip()
    if not val_str or val_str == "0":
        return 0
    try:
        val_float = float(val_str)
    except ValueError as err:
        raise ValueError(f"Invalid count value: {val!r}") from err
    if not val_float.is_integer():
        raise ValueError(f"Count must be an integer: {val!r}")
    if val_float < 0:
        raise ValueError(f"Negative count: {val!r}")
    return int(val_float)


def trait_population_counts(taxa_matrix_df: pd.DataFrame, trait_columns: list[str]) -> np.ndarray:
    """Count taxa carrying each trait (K) in the background universe."""
    counts = np.zeros(len(trait_columns), dtype=np.int64)
    for col_idx, col in enumerate(trait_columns):
        series = taxa_matrix_df[col]
        if pd.api.types.is_integer_dtype(series):
            arr = series.to_numpy()
            if np.any(arr < 0):
                raise ValueError(f"Negative count in column {col}")
            counts[col_idx] = int(np.sum(arr > 0))
        else:
            col_k = 0
            for val in series:
                c = _parse_strict_count(val)
                if c > 0:
                    col_k += 1
            counts[col_idx] = col_k
    return counts


def generate_complete_null_ora_profiles(
    empirical_profiles_df: pd.DataFrame,
    taxa_trait_matrix_df: pd.DataFrame,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Generate synthetic complete null profiles under the ORA null hypothesis (H0_ORA).

    Each profile of signature size N is generated by sampling N unique taxa
    without replacement from the M taxa universe (the rows of taxa_trait_matrix_df).
    Taxa matrix cells are strictly binarized to presence (> 0) so k <= N always (#348).
    Summing the whole taxon vectors across the sampled taxa preserves all
    inter-trait correlations and co-annotation bundles by construction.
    """
    trait_columns = shared_trait_columns(empirical_profiles_df, taxa_trait_matrix_df)
    if not trait_columns:
        raise ValueError("No trait columns shared between profiles and taxa matrix")
    n_universe = len(taxa_trait_matrix_df)
    signature_sizes = resolve_signature_sizes(empirical_profiles_df)
    if signature_sizes.size and int(signature_sizes.max()) > n_universe:
        raise ValueError(
            f"Signature size {int(signature_sizes.max())} exceeds taxa universe size {n_universe}"
        )

    n_profiles = len(empirical_profiles_df)
    n_traits = len(trait_columns)
    simulated_k = np.zeros((n_profiles, n_traits), dtype=np.int64)

    # Convert taxa trait matrix to numeric binary presence matrix (M x T) (#348)
    numeric_matrix = (
        taxa_trait_matrix_df[trait_columns]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
        .to_numpy()
    )
    taxa_traits = (numeric_matrix > 0).astype(np.int64)

    for i, n_draw in enumerate(signature_sizes):
        if n_draw <= 0:
            continue
        sampled_taxa = rng.choice(n_universe, size=int(n_draw), replace=False)
        simulated_k[i, :] = taxa_traits[sampled_taxa].sum(axis=0)

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


def generate_study_flip_null_ora_profiles(
    empirical_profiles_df: pd.DataFrame,
    rng: np.random.Generator,
    *,
    unit: str = "study",
) -> pd.DataFrame:
    """Generate study-level direction flip null profiles (#263, #343).

    Note (#388): While this exchangeability null is valid for cross-profile screen/direction
    statistics (such as PR #319), swapping inc/dec labels leaves individual profile
    taxa compositions unchanged and therefore cannot alter within-profile hypergeometric
    ORA p-values. It is retained here for shared study-flip generation but dropped from
    the active within-profile ORA calibration scenarios in the harness.
    """
    return flip_profile_dataframe(empirical_profiles_df, rng, unit=unit)


def generate_mixed_signal_ora_profiles(
    empirical_profiles_df: pd.DataFrame,
    taxa_trait_matrix_df: pd.DataFrame,
    spike_traits: list[str],
    signal_strength: float = 0.5,
    rng: np.random.Generator | None = None,
    *,
    test_type: str = "enrichment",
    spike_fraction: float = 0.5,
) -> tuple[pd.DataFrame, set[tuple[str, str]]]:
    """Generate mixed signal profiles with known ground-truth spikes (#348, #351, #390).

    Preserves whole taxon vectors by sampling N unique taxa per profile:
    - Leaves unspiked profiles (1 - spike_fraction) as pure uniform null draws (#390).
    - For spiked profiles, draws signal-carrying taxa (enriched or depleted) without replacement.
    - Treats near-twin/co-annotated traits with identical presence vectors in the universe
      as ground-truth true discoveries alongside the spiked trait (#390).
    - Supports both enrichment and depletion test scopes (#351, #390).
    """
    if not 0.0 < signal_strength <= 1.0:
        raise ValueError(f"signal_strength must be in (0, 1], got {signal_strength}")
    if not 0.0 < spike_fraction <= 1.0:
        raise ValueError(f"spike_fraction must be in (0, 1], got {spike_fraction}")
    if rng is None:
        rng = np.random.default_rng()

    trait_columns = shared_trait_columns(empirical_profiles_df, taxa_trait_matrix_df)
    unknown = [t for t in spike_traits if t not in trait_columns]
    if unknown:
        raise ValueError(f"Spike traits not found in shared trait columns: {unknown}")

    n_universe = len(taxa_trait_matrix_df)
    signature_sizes = resolve_signature_sizes(empirical_profiles_df)
    if signature_sizes.size and int(signature_sizes.max()) > n_universe:
        raise ValueError(
            f"Signature size {int(signature_sizes.max())} exceeds taxa universe size {n_universe}"
        )

    profile_ids = profile_identifiers(empirical_profiles_df)
    n_profiles = len(empirical_profiles_df)
    n_traits = len(trait_columns)
    simulated_k = np.zeros((n_profiles, n_traits), dtype=np.int64)

    # Convert taxa trait matrix to numeric binary presence matrix (M x T) (#348)
    numeric_matrix = (
        taxa_trait_matrix_df[trait_columns]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
        .to_numpy()
    )
    taxa_traits = (numeric_matrix > 0).astype(np.int64)

    trait_to_idx = {col: idx for idx, col in enumerate(trait_columns)}
    true_signal_pairs: set[tuple[str, str]] = set()

    # Issue #390: Map each spike trait to all near-twin traits with identical presence vectors
    twin_traits_map: dict[str, list[str]] = {}
    for trait in spike_traits:
        t_idx = trait_to_idx[trait]
        target_col = taxa_traits[:, t_idx]
        twins = [
            other_trait
            for other_trait, other_idx in trait_to_idx.items()
            if np.array_equal(target_col, taxa_traits[:, other_idx])
        ]
        twin_traits_map[trait] = twins

    spike_info = {}
    for trait in spike_traits:
        t_idx = trait_to_idx[trait]
        pos_mask = taxa_traits[:, t_idx] > 0
        pos_indices = np.flatnonzero(pos_mask)
        neg_indices = np.flatnonzero(~pos_mask)
        spike_info[trait] = (t_idx, pos_indices, neg_indices)

    num_spikes = len(spike_traits)
    # Determine which profiles are spiked versus unspiked null (#390)
    # E.g. with spike_fraction=0.5, even-indexed profiles are spiked, odd are pure null.
    stride = max(1, round(1.0 / spike_fraction))

    for i, n_draw in enumerate(signature_sizes):
        if n_draw <= 0:
            continue
        p_id = str(profile_ids[i])
        is_spiked_profile = (i % stride) == 0

        if not is_spiked_profile:
            # Pure null profile: uniform draw from universe (#390)
            sampled_taxa = rng.choice(n_universe, size=int(n_draw), replace=False)
            simulated_k[i, :] = taxa_traits[sampled_taxa].sum(axis=0)
            continue

        spike_trait = spike_traits[(i // stride) % num_spikes]
        t_idx, pos_indices, neg_indices = spike_info[spike_trait]
        k_total = len(pos_indices)

        if k_total == 0:
            sampled_taxa = rng.choice(n_universe, size=int(n_draw), replace=False)
            simulated_k[i, :] = taxa_traits[sampled_taxa].sum(axis=0)
            continue

        if test_type == "depletion":
            # Depletion spike (#351, #390): draw taxa exclusively from negative indices
            # so k_observed = 0, creating a significant depletion signal
            if len(neg_indices) >= int(n_draw):
                sampled_taxa = rng.choice(neg_indices, size=int(n_draw), replace=False)
            else:
                n_neg = len(neg_indices)
                n_pos = int(n_draw) - n_neg
                pos_taxa = rng.choice(pos_indices, size=n_pos, replace=False)
                sampled_taxa = np.concatenate([pos_taxa, neg_indices])

            simulated_k[i, :] = taxa_traits[sampled_taxa].sum(axis=0)
            for twin in twin_traits_map[spike_trait]:
                true_signal_pairs.add((p_id, twin))
        else:
            # Enrichment spike (#348, #390): oversample positive taxa
            desired_pos = max(1, round(n_draw * signal_strength))
            min_pos = max(0, int(n_draw - len(neg_indices)))
            max_pos = min(int(n_draw), k_total)
            n_pos = min(max_pos, max(min_pos, desired_pos))
            n_neg = int(n_draw) - n_pos

            pos_taxa = (
                rng.choice(pos_indices, size=n_pos, replace=False)
                if n_pos > 0
                else np.array([], dtype=int)
            )
            neg_taxa = (
                rng.choice(neg_indices, size=n_neg, replace=False)
                if n_neg > 0
                else np.array([], dtype=int)
            )
            sampled_taxa = np.concatenate([pos_taxa, neg_taxa])
            simulated_k[i, :] = taxa_traits[sampled_taxa].sum(axis=0)

            if n_pos > 0:
                for twin in twin_traits_map[spike_trait]:
                    true_signal_pairs.add((p_id, twin))

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

    return pd.DataFrame(sim_dict), true_signal_pairs
