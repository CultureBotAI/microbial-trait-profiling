"""Audit clustering integration without changing scientific inputs or exports.

Purpose: verify inference preservation and diagnostic correctness for #198, #220,
and #221. Inputs: baseline contrast TSV and regenerated six-trait/oral exports.
Output: JSON audit summary on stdout; nonzero exit on any failed invariant.
Date: 2026-10-15. Run with the repository's audit-six-trait-diagnostics Justfile recipe.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_DIR = (
    REPO_ROOT / "trait-visualization" / "trait-profile-investigation" / "six-shared-traits"
)
CONTRAST_FILENAME = "six_shared_trait_condition_body_site_direction_contrasts.tsv"
KEYS = ["view_key", "condition_body_site_pair", "trait_id"]
IDENTITY_COLUMNS = [
    "view_label",
    "focus_direction",
    "condition_mode",
    "condition",
    "body_site_group",
    "trait_label",
    "blank_reason",
]
COUNT_COLUMNS = [
    "trait_rank",
    "study_count",
    "increased_positive_profile_count",
    "increased_profile_count",
    "decreased_positive_profile_count",
    "decreased_profile_count",
]
FLOAT_COLUMNS = [
    "increased_positive_fraction",
    "decreased_positive_fraction",
    "increased_minus_decreased_positive_fraction",
    "fisher_p_value",
    "fisher_q_value",
    "minus_log10_fisher_p_value",
    "plot_minus_log10_fisher_p_value",
    "signed_minus_log10_fisher_p_value",
    "plot_signed_minus_log10_fisher_p_value",
]
FLAG_COLUMNS = [
    "passes_support",
    "passes_nominal_p_value_cutoff",
    "passes_p_value_cutoff",
    "nominal_p_value_selected",
    "displayed_in_heatmap",
    "significant_in_heatmap",
    "selected_after_bh",
    "passes_absolute_effect_cutoff",
]
NUMERIC_RTOL = 1e-11  # Covers rounding to twelve significant digits, including tiny p-values.


def require(condition: bool, message: str) -> None:
    """Fail an audit invariant even when Python assertions are disabled."""
    if not condition:
        raise ValueError(message)


def read_table(path: Path) -> pd.DataFrame:
    """Read an artifact, rejecting unresolved merge markers."""
    text = path.read_text(encoding="utf-8")
    require(
        not any(line.startswith(("<<<<<<< ", ">>>>>>> ")) for line in text.splitlines()),
        f"Unresolved merge markers: {path}",
    )
    return pd.read_csv(path, sep="\t")


def keyed(table: pd.DataFrame, keys: list[str], label: str) -> pd.DataFrame:
    """Require complete unique keys before aligning tables."""
    require(not table[keys].isna().any().any(), f"Missing keys in {label}")
    require(not table.duplicated(keys).any(), f"Duplicate keys in {label}")
    return table.set_index(keys).sort_index()


def flags(values: pd.Series) -> pd.Series:
    """Parse Boolean exports without treating a literal False as truthy."""
    strings = values.astype(str)
    require(strings.isin(["True", "False"]).all(), f"Invalid Boolean field {values.name}")
    return strings.eq("True")


def same_numeric(left: pd.Series, right: pd.Series, *, exact: bool = False) -> bool:
    """Compare numbers, preserving missing values and significant-digit precision."""
    a = pd.to_numeric(left, errors="raise").to_numpy(dtype=float, na_value=np.nan)
    b = pd.to_numeric(right, errors="raise").to_numpy(dtype=float, na_value=np.nan)
    return bool(np.allclose(a, b, rtol=0.0 if exact else NUMERIC_RTOL, atol=0.0, equal_nan=True))


def verify_inference(expected: pd.DataFrame, actual: pd.DataFrame, label: str) -> None:
    """Check complete cell identity, inference, support, and selection preservation."""
    left = keyed(expected, KEYS, f"{label} expected")
    right = keyed(actual, KEYS, f"{label} actual")
    require(left.index.equals(right.index), f"Changed or missing cell keys: {label}")
    for column in IDENTITY_COLUMNS:
        equal = left[column].fillna("").astype(str).eq(right[column].fillna("").astype(str))
        require(equal.all(), f"Changed {column}: {label}")
    for column in COUNT_COLUMNS + FLOAT_COLUMNS:
        require(
            same_numeric(left[column], right[column], exact=column in COUNT_COLUMNS),
            f"Changed numeric field {column}: {label}",
        )
    for column in FLAG_COLUMNS:
        require(flags(left[column]).equals(flags(right[column])), f"Changed {column}: {label}")


def display_matrices(contrasts: pd.DataFrame, traits: list[str]) -> dict[str, pd.DataFrame]:
    """Build supported displayed vectors directly from the audited contrast table."""
    matrices = {}
    for view, cells in contrasts.groupby("view_key"):
        selected = cells.loc[flags(cells["displayed_in_heatmap"])].copy()
        # Reconstruct the calculation precision from integer counts. Ranking the
        # twelve-digit exported effects can merge near-equal floating-point ranks.
        selected["audit_effect"] = (
            selected["increased_positive_profile_count"] / selected["increased_profile_count"]
            - selected["decreased_positive_profile_count"] / selected["decreased_profile_count"]
        )
        matrices[str(view)] = selected.pivot(
            index="condition_body_site_pair",
            columns="trait_id",
            values="audit_effect",
        ).reindex(columns=traits)
    return matrices


def correlation(left: np.ndarray, right: np.ndarray) -> tuple[float, int]:
    """Independently recover pairwise-complete Spearman rho and shared support."""
    finite = np.isfinite(left) & np.isfinite(right)
    shared = int(finite.sum())
    if shared < 3:
        return np.nan, shared
    a = pd.Series(left[finite]).rank().to_numpy()
    b = pd.Series(right[finite]).rank().to_numpy()
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan, shared
    return float(np.clip(np.corrcoef(a, b)[0, 1], -1.0, 1.0)), shared


def source_contributors(
    output_dir: Path,
) -> tuple[dict[str, set[str]], dict[tuple[str, str], str]]:
    """Count distinct raw condition labels across full grouped source edge tables."""
    manifest = read_table(output_dir / "six_shared_trait_graph_manifest.tsv")
    contributors: dict[str, set[str]] = {}
    mappings: dict[tuple[str, str], str] = {}
    for row in manifest.loc[
        manifest["view_key"].isin(["increased_grouped", "decreased_grouped"])
    ].itertuples(index=False):
        display_path = Path(row.source_display_edges_path)
        display_path = display_path if display_path.is_absolute() else REPO_ROOT / display_path
        suffix = "_display_edges.tsv"
        require(display_path.name.endswith(suffix), f"Unexpected source path {display_path}")
        source_path = display_path.with_name(display_path.name.removesuffix(suffix) + "_edges.tsv")
        edges = read_table(source_path)
        for edge in edges.itertuples(index=False):
            pair = f"{edge.condition} | {edge.body_site_group}"
            raw = edge.raw_conditions
            if pd.isna(raw):
                continue
            for condition in str(raw).split(";"):
                condition = condition.strip()
                if not condition:
                    continue
                key = (condition, str(edge.body_site_group))
                require(key not in mappings or mappings[key] == pair, f"Ambiguous rollup: {key}")
                mappings[key] = pair
                contributors.setdefault(pair, set()).add(condition)
    require(bool(contributors), "No grouped source contributors found")
    return contributors, mappings


def audit_comparisons(output_dir: Path, matrices: dict[str, pd.DataFrame]) -> dict[str, int]:
    """Recompute overlap eligibility, diagnostic families, ranks, and contributor counts."""
    table = read_table(output_dir / "six_shared_trait_signed_direction_cluster_comparison.tsv")
    table = keyed(table, ["individual_condition_body_site_pair"], "comparison")
    individual = matrices["individual_conditions_direction_contrast"]
    grouped = matrices["grouped_direction_contrast"]
    require(set(table.index) == set(individual.index), "Comparison rows differ from displayed rows")
    contributors, mappings = source_contributors(output_dir)
    displayed_counts = table.groupby("grouped_condition_body_site_pair").size()
    rhos = []
    for pair, row in table.iterrows():
        group_pair = row["grouped_condition_body_site_pair"]
        source_pair = mappings.get((row["individual_condition"], row["body_site_group"]))
        require(
            (pd.isna(group_pair) and source_pair is None) or group_pair == source_pair,
            f"Incorrect grouped condition mapping: {pair}",
        )
        if group_pair in grouped.index:
            rho, shared = correlation(
                individual.loc[pair].to_numpy(dtype=float),
                grouped.loc[group_pair].to_numpy(dtype=float),
            )
        else:
            rho, shared = np.nan, 0
        require(shared == row["shared_supported_trait_count"], f"Incorrect shared traits: {pair}")
        rhos.append(rho)
        if pd.notna(group_pair):
            require(
                row["grouped_contributing_individual_condition_count"]
                == len(contributors.get(group_pair, set())),
                f"Incorrect source contributors: {pair}",
            )
            require(
                row["grouped_displayed_individual_condition_count"] == displayed_counts[group_pair],
                f"Incorrect displayed contributors: {pair}",
            )
    rho = pd.Series(rhos, index=table.index)
    require(
        same_numeric(rho, table["individual_vs_grouped_spearman_rho"]), "Incorrect comparison rho"
    )
    distances = 1.0 - rho
    require(
        same_numeric(distances, table["individual_vs_grouped_spearman_distance"]),
        "Incorrect comparison distances",
    )
    eligible = rho.notna() & table["shared_supported_trait_count"].ge(5)
    require(flags(table["pattern_diagnostic_eligible"]).equals(eligible), "Incorrect eligibility")
    eligible_distances = distances.where(eligible)
    groups = table["grouped_condition_body_site_pair"]
    require(
        same_numeric(
            eligible_distances.rank(method="min", ascending=False),
            table["global_difference_rank"],
            exact=True,
        ),
        "Incorrect global ranks",
    )
    require(
        same_numeric(
            eligible_distances.groupby(groups).rank(method="min", ascending=False),
            table["within_group_difference_rank"],
            exact=True,
        ),
        "Incorrect group ranks",
    )
    thresholds = pd.Series(np.nan, index=table.index)
    for _group, members in table.groupby("grouped_condition_body_site_pair"):
        values = eligible_distances.loc[members.index].dropna().to_numpy()
        if len(values) >= 4:
            lower, upper = np.quantile(values, [0.25, 0.75])
            thresholds.loc[members.index] = upper + 1.5 * (upper - lower)
    require(
        np.allclose(
            thresholds,
            table["within_group_outlier_threshold"],
            rtol=NUMERIC_RTOL,
            atol=1e-11,
            equal_nan=True,
        ),
        "Incorrect eligible-only IQR threshold",
    )
    outliers = eligible & distances.gt(thresholds)
    discordant = eligible & rho.le(0.0)
    require(flags(table["within_group_iqr_outlier"]).equals(outliers), "Incorrect IQR flags")
    require(
        flags(table["discordant_or_inverse_pattern"]).equals(discordant), "Incorrect inverse flags"
    )
    return {
        "rows": len(table),
        "eligible_rows": int(eligible.sum()),
        "iqr_eligible_rows": int((eligible & thresholds.notna()).sum()),
        "iqr_outliers": int(outliers.sum()),
        "discordant_rows": int(discordant.sum()),
    }


def audit_leaves(output_dir: Path, matrices: dict[str, pd.DataFrame]) -> dict[str, dict[str, int]]:
    """Independently confirm that every linked row has at least one defined neighbor."""
    leaves = read_table(output_dir / "six_shared_trait_signed_direction_cluster_tree_leaves.tsv")
    keyed(leaves, ["view_key", "condition_body_site_pair"], "tree leaves")
    summary = {}
    require(set(leaves["view_key"]) == set(matrices), "Tree leaf views differ from contrast views")
    for view, matrix in matrices.items():
        rows = leaves.loc[leaves["view_key"].eq(view)].set_index("condition_body_site_pair")
        require(set(rows.index) == set(matrix.index), f"Tree leaves omit displayed rows: {view}")
        values = matrix.to_numpy(dtype=float)
        for position, label in enumerate(matrix.index):
            has_neighbor = any(
                np.isfinite(correlation(values[position], other)[0])
                for index, other in enumerate(values)
                if index != position
            )
            included = bool(flags(rows.loc[[label], "included_in_clustering"]).iloc[0])
            require(included == has_neighbor, f"Incorrect linkage inclusion: {view}/{label}")
            require(
                (rows.loc[label, "defined_spearman_neighbor_count"] > 0) == has_neighbor,
                f"Incorrect defined-neighbor support: {view}/{label}",
            )
        included_count = int(flags(rows["included_in_clustering"]).sum())
        summary[view] = {
            "rows": len(rows),
            "linked_rows": included_count,
            "insufficient_overlap_rows": len(rows) - included_count,
        }
    return summary


def audit_oral(output_dir: Path, contrasts: pd.DataFrame) -> dict[str, object]:
    """Validate source hashes, oral selection, and retained full-family inference."""
    oral_dir = output_dir / "oral"
    provenance = json.loads((oral_dir / "oral_six_trait_provenance.json").read_text())
    for name_key, hash_key in [
        ("source_table", "source_table_sha256"),
        ("source_trait_table", "source_trait_table_sha256"),
    ]:
        actual_hash = hashlib.sha256((output_dir / provenance[name_key]).read_bytes()).hexdigest()
        require(actual_hash == provenance[hash_key], f"Stale oral source hash: {name_key}")
    individual = contrasts.loc[contrasts["view_key"].eq("individual_conditions_direction_contrast")]
    oral_site = individual["body_site_group"].isin(provenance["oral_body_site_groups"])
    disease_labels = set(provenance["oral_disease_labels_casefolded"])
    oral_disease = individual["condition"].map(
        lambda value: any(part.strip().casefold() in disease_labels for part in value.split(","))
    )
    expected = individual.loc[oral_site | oral_disease]
    cells = read_table(oral_dir / "oral_six_trait_cells.tsv")
    verify_inference(expected, cells, "oral subset retains source inference")
    rows = keyed(
        read_table(oral_dir / "oral_six_trait_rows.tsv"),
        ["condition_body_site_pair"],
        "oral row summary",
    )
    require(set(rows.index) == set(cells["condition_body_site_pair"]), "Oral row coverage changed")
    displayed = flags(cells["displayed_in_heatmap"])
    bh = flags(cells["selected_after_bh"])
    counts = {
        "source_individual_bh_family_size": len(individual),
        "selected_rows": len(rows),
        "main_figure_rows": int(flags(rows["in_main_figure"]).sum()),
        "selected_cells": len(cells),
        "displayed_cells": int(displayed.sum()),
        "bh_significant_displayed_cells": int((displayed & bh).sum()),
        "bh_selected_cell_count": int(bh.sum()),
        "nominal_p_value_cell_count": int(flags(cells["nominal_p_value_selected"]).sum()),
        "significant_cell_count": int(flags(cells["significant_in_heatmap"]).sum()),
    }
    for key, value in counts.items():
        require(provenance[key] == value, f"Incorrect oral provenance count: {key}")
    for filename in provenance["figures"]:
        require((oral_dir / filename).is_file(), f"Missing oral figure: {filename}")
    return {
        **counts,
        "source_table_sha256": provenance["source_table_sha256"],
        "source_trait_table_sha256": provenance["source_trait_table_sha256"],
    }


def main() -> None:
    """Print an auditable summary or a structured failure and nonzero exit status."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-contrasts", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    try:
        baseline = read_table(args.baseline_contrasts)
        contrasts = read_table(args.output_dir / CONTRAST_FILENAME)
        verify_inference(baseline, contrasts, "baseline contrast preservation")
        traits = read_table(args.output_dir / "six_shared_trait_selection.tsv")["trait_id"].tolist()
        require(len(traits) == len(set(traits)) == 6, "Expected six distinct selected traits")
        matrices = display_matrices(contrasts, traits)
        summary = {
            "status": "passed",
            "baseline_contrasts": str(args.baseline_contrasts),
            "baseline_sha256": hashlib.sha256(args.baseline_contrasts.read_bytes()).hexdigest(),
            "numeric_relative_tolerance": NUMERIC_RTOL,
            "preserved_columns": IDENTITY_COLUMNS + COUNT_COLUMNS + FLOAT_COLUMNS + FLAG_COLUMNS,
            "contrast_cells": len(contrasts),
            "views": {
                str(view): {
                    "cells": len(cells),
                    "displayed_cells": int(flags(cells["displayed_in_heatmap"]).sum()),
                    "bh_selected_cells": int(flags(cells["selected_after_bh"]).sum()),
                }
                for view, cells in contrasts.groupby("view_key")
            },
            "comparisons": audit_comparisons(args.output_dir, matrices),
            "clustering": audit_leaves(args.output_dir, matrices),
            "oral": audit_oral(args.output_dir, contrasts),
        }
    except (ValueError, KeyError, OSError) as error:
        print(json.dumps({"status": "failed", "error": str(error)}, indent=2))
        raise SystemExit(1) from error
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
