"""Run and verify the standard Ilyas, Noel and Jed analysis suites (#658)."""

from __future__ import annotations

import argparse
from pathlib import Path

from mtp.workflow.arms import ilyas, jed, noel
from mtp.workflow.arms.runtime import execute_arm, verify_arm
from mtp.workflow.comparison import execute_comparison_run, verify_comparison_run


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("ilyas", "noel", "jed"):
        arm = commands.add_parser(
            name, help=f"Execute {name}'s standard suite from a verified core run"
        )
        arm.add_argument("core_run", type=Path, help="Completed core run.json from just core-run")
        arm.add_argument(
            "output", type=Path, help="Fresh destination; existing directories are rejected"
        )
        arm.add_argument(
            "--replicates",
            type=int,
            help="Override stochastic-stage replicate counts; otherwise preserve stage defaults",
        )
        arm.add_argument(
            "--seed",
            type=int,
            help="Override stochastic-stage seeds; otherwise preserve stage defaults",
        )
        arm.add_argument(
            "--no-figures",
            action="store_true",
            help="Run analysis and calibration without optional figure rendering",
        )
    verify = commands.add_parser(
        "verify", help="Verify a completed arm execution and its core bindings"
    )
    verify.add_argument("record", type=Path)
    compare = commands.add_parser(
        "compare",
        help="Execute cross-arm comparison consuming verified core, Noel, Jed, and Ilyas manifests",
    )
    compare.add_argument("core_run", type=Path, help="Verified core run.json")
    compare.add_argument("noel_run", type=Path, help="Verified Noel arm_run.json")
    compare.add_argument("jed_run", type=Path, help="Verified Jed arm_run.json")
    compare.add_argument("ilyas_run", type=Path, help="Verified Ilyas arm_run.json")
    compare.add_argument("output", type=Path, help="Fresh destination for comparison_run.json")
    compare.add_argument(
        "--arm1-pooled-stats",
        type=Path,
        help="Optional precomputed trait_profile_paired_contrast_trait_stats_pooled.tsv",
    )
    verify_cmp = commands.add_parser(
        "verify-comparison",
        help="Verify a completed ComparisonRun manifest and its complete cross-arm outputs",
    )
    verify_cmp.add_argument("record", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "verify":
        record = verify_arm(args.record)
        print(f"Verified {record.id}: executed stages, core binding and output bytes")
        return 0
    if args.command == "compare":
        cmp_record = execute_comparison_run(
            args.core_run,
            args.noel_run,
            args.jed_run,
            args.ilyas_run,
            args.output,
            arm1_pooled_stats=args.arm1_pooled_stats,
        )
        print(cmp_record.manifest_path)
        return 0
    if args.command == "verify-comparison":
        cmp_record = verify_comparison_run(args.record)
        print(f"Verified {cmp_record.id}: cross-arm bindings, analysis bases and complete outputs")
        return 0
    if args.replicates is not None and args.replicates < 1:
        parser.error("--replicates must be positive")
    planners = {
        "ilyas": ilyas.build_commands,
        "noel": noel.build_commands,
        "jed": jed.build_commands,
    }
    record = execute_arm(
        args.command,
        args.core_run,
        args.output,
        planners[args.command],
        replicates=args.replicates,
        seed=args.seed,
        figures=not args.no_figures,
    )
    print(record.manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
