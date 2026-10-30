"""CLI for additive verified post-arm multiplicity sidecars."""

import argparse
from pathlib import Path

from mtp.multiplicity.stage import run_post_arm, verify_post_arm
from mtp.workflow.output_roots import add_temporary_output_option, require_durable_output


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run")
    run.add_argument("--comparison-run", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True)
    add_temporary_output_option(run)
    verify = commands.add_parser("verify")
    verify.add_argument("--run", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "run":
        require_durable_output(parser, args.output, allowed=args.allow_temporary_output)
        record = run_post_arm(args.comparison_run, args.output)
    else:
        record = verify_post_arm(args.run)
    print(f"Verified {len(record.families)} additive families for {record.core_run_id}")


if __name__ == "__main__":
    main()
