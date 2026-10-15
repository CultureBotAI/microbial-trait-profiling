"""Audit trait coupling sensitivity via LOTO, LCO, and multi-arm evaluation.

Addresses Issue #260:
- Leave-One-Trait-Out (LOTO): Computes 6x6 matrix of delta effect and delta q across Table 1 traits.
- Leave-Cluster-Out (LCO): Removes coccus + sphere simultaneously and measures shifts across 3 arms.
- Three-Arm Sensitivity: Compares Paired Composition, Directional Sign-Test,
  and ORA Hypergeometric arms.
- Evaluates the formal stability decision rule against Table 1 inter-row spacing.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mtp.coupling import DEFAULT_OUTPUT_DIR, run_coupling_audit
from mtp.density_adjustment import DEFAULT_COLUMNS_PATH, DEFAULT_PROFILES_PATH


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit trait coupling sensitivity via LOTO, LCO, and multi-arm evaluation."
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
        help=f"Directory to write output TSVs and reports (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--code-commit",
        type=str,
        default="",
        help="Code commit hash for provenance tracking",
    )
    return parser.parse_args(argv)


def main() -> int:
    args = parse_args()
    try:
        results = run_coupling_audit(
            profiles_path=args.profiles,
            columns_path=args.columns,
            output_dir=args.output_dir,
            code_commit=args.code_commit,
        )
    except Exception as exc:
        sys.stderr.write(f"Error during trait coupling audit: {exc}\n")
        return 1

    rule = results["rule_result"]
    rule_rev = results["rule_result_rev"]

    print("\n" + "=" * 80)
    print("TRAIT COUPLING SENSITIVITY AUDIT COMPLETED")
    print("=" * 80)
    print(f"LOTO matrix written to:       {results['loto_tsv']}")
    print(f"LCO summary written to:      {results['lco_tsv']}")
    print(f"Comprehensive report to:     {results['report_txt']}")
    print("\nDecision Rule Summary:")
    print(f"  - Minimum adjacent Table 1 distance: {rule.min_adjacent_distance:.6f}")
    print(
        f"  - Dropping sphere -> Coccus effect change: {rule.delta_effect:+.6f} "
        f"(|Delta| = {rule.abs_delta_effect:.6f})"
    )
    print(
        f"  - Dropping coccus -> Sphere effect change: {rule_rev.delta_effect:+.6f} "
        f"(|Delta| = {rule_rev.abs_delta_effect:.6f})"
    )
    print(f"  - Criterion Passed: {rule.passes_stability_criterion}")
    print(f"  - Verdict: {rule.verdict}")
    print("=" * 80 + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
