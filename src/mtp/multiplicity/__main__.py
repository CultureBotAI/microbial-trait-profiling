"""CLI for additive verified post-arm multiplicity sidecars."""

import argparse
from pathlib import Path

from mtp.multiplicity.stage import run_post_arm, verify_post_arm


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run")
    run.add_argument("--comparison-run", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True)
    verify = commands.add_parser("verify")
    verify.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    record = (
        run_post_arm(args.comparison_run, args.output)
        if args.command == "run"
        else verify_post_arm(args.run)
    )
    print(f"Verified {len(record.families)} additive families for {record.core_run_id}")


if __name__ == "__main__":
    main()
