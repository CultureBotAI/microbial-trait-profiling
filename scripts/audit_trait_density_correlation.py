"""Audit trait correlation, density adjustment, and near-twin collinearity (#264, #352).

Addresses Issue #264 and Issue #352:
- Computes baseline pairwise correlation distribution and PC1 variance across non-MediaDive traits.
- Evaluates PC1 regression adjustment and profile composition normalization.
- Tests collapse of facultatively aerobic x moderately halophilic vs. survival of convergent six.
- Tests conditional contrast and VIF for coccus shaped vs. sphere shaped.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mtp.density_adjustment import (
    DEFAULT_COLUMNS_PATH,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_PROFILES_PATH,
    run_density_audit,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Audit trait correlation, annotation density adjustment, and conditional near-twins."
        )
    )
    parser.add_argument(
        "--profiles",
        type=Path,
        default=DEFAULT_PROFILES_PATH,
        help=f"Path to study_trait_profiles_wide.tsv (default: {DEFAULT_PROFILES_PATH})",
    )
    parser.add_argument(
        "--columns",
        type=Path,
        default=DEFAULT_COLUMNS_PATH,
        help=f"Path to taxa_trait_matrix_columns.tsv (default: {DEFAULT_COLUMNS_PATH})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory to write summary TSVs and reports (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--code-commit",
        type=str,
        default="",
        help="Code commit hash for provenance tracking",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        results = run_density_audit(
            profiles_path=args.profiles,
            columns_path=args.columns,
            output_dir=args.output_dir,
            code_commit=args.code_commit,
        )
    except Exception as exc:
        sys.stderr.write(f"Error during density adjustment audit: {exc}\n")
        return 1

    print("\n" + "=" * 80)
    print("DENSITY ADJUSTMENT & REDUNDANCY AUDIT COMPLETED")
    print("=" * 80)
    print(f"Summary metrics written to:      {results['summary_tsv']}")
    print(f"Table 1 adjusted effects to:    {results['table1_tsv']}")
    print(f"Comprehensive report to:        {results['report_txt']}")
    print("\nKey Findings:")
    raw = results["raw_pearson"]
    adj = results["adj_pc1_pearson"]
    twin = results["twin_summary"]
    pair_attenuation = adj.test_pair_r - raw.test_pair_r
    print(
        f"  - Baseline PC1 explains {raw.pc1_variance_pct:.2f}% of variance "
        f"across {raw.n_traits} traits."
    )
    print(
        f"  - Post-adjustment pairwise correlation centered at median {adj.median_pairwise_r:+.4f} "
        f"(was {raw.median_pairwise_r:+.4f})."
    )
    print(
        f"  - facultatively aerobic x moderately halophilic r: {raw.test_pair_r:+.4f} "
        f"-> {adj.test_pair_r:+.4f} (attenuated by {pair_attenuation:+.4f})."
    )
    print(
        f"  - Near-twins (coccus x sphere) collinearity: r = {twin.profile_correlation:+.4f}, "
        f"VIF = {twin.vif_profile:.2f}."
    )
    print("=" * 80 + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
