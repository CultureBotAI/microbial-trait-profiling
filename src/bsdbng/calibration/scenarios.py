"""Simulation scenario generators for statistical calibration."""

import numpy as np
import pandas as pd


def generate_complete_null_ora_profiles(
    empirical_profiles_df: pd.DataFrame,
    taxa_trait_matrix_df: pd.DataFrame,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Generate complete null profiles under the ORA null hypothesis (H0_ORA)."""
    all_taxa = taxa_trait_matrix_df.index.to_numpy()
    trait_columns = taxa_trait_matrix_df.columns.tolist()

    meta_cols = [
        col
        for col in ["study_id", "experiment_id", "inc_dec", "condition", "body_site"]
        if col in empirical_profiles_df.columns
    ]

    simulated_rows = []

    for _, row in empirical_profiles_df.iterrows():
        trait_counts = row[trait_columns].to_numpy(dtype=int)
        max_c = (
            int(np.max(trait_counts))
            if len(trait_counts) > 0 and np.max(trait_counts) > 0
            else 10
        )
        profile_n = int(row.get("signature_size", max_c))
        if profile_n <= 0:
            profile_n = 10

        sampled_taxa_idx = rng.choice(
            len(all_taxa), size=min(profile_n, len(all_taxa)), replace=False
        )
        sampled_matrix = taxa_trait_matrix_df.iloc[sampled_taxa_idx]

        k_counts = (sampled_matrix > 0).sum(axis=0).to_dict()

        row_dict = {col: row[col] for col in meta_cols}
        row_dict["signature_size"] = profile_n
        row_dict.update(k_counts)
        simulated_rows.append(row_dict)

    return pd.DataFrame(simulated_rows)


def generate_mixed_signal_ora_profiles(
    empirical_profiles_df: pd.DataFrame,
    taxa_trait_matrix_df: pd.DataFrame,
    spike_traits: list[str],
    signal_strength: float = 0.5,
    rng: np.random.Generator | None = None,
) -> tuple[pd.DataFrame, set[tuple[str, str]]]:
    """Generate mixed signal profiles with known ground-truth spikes."""
    if rng is None:
        rng = np.random.default_rng()

    null_df = generate_complete_null_ora_profiles(empirical_profiles_df, taxa_trait_matrix_df, rng)
    trait_columns = taxa_trait_matrix_df.columns.tolist()
    true_signal_pairs = set()

    for idx, row in null_df.iterrows():
        study = row.get("study_id", "")
        exp = row.get("experiment_id", "")
        direction = row.get("inc_dec", "")
        profile_id = f"{study}:{exp}:{direction}"
        n = int(row.get("signature_size", 10))

        for trait in spike_traits:
            if trait not in trait_columns:
                continue

            trait_taxa_mask = taxa_trait_matrix_df[trait] > 0
            trait_taxa_indices = np.where(trait_taxa_mask)[0]

            if len(trait_taxa_indices) == 0:
                continue

            spike_count = max(1, round(n * signal_strength))
            chosen_trait_taxa = rng.choice(
                trait_taxa_indices,
                size=min(spike_count, len(trait_taxa_indices)),
                replace=False,
            )

            current_k = row[trait]
            new_k = min(n, current_k + len(chosen_trait_taxa))
            null_df.at[idx, trait] = new_k
            true_signal_pairs.add((profile_id, trait))

    return null_df, true_signal_pairs
