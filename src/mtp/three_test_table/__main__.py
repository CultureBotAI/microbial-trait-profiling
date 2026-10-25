"""Command-line interface for the declared per-trait sidecar (#1044)."""

import argparse
from pathlib import Path

from mtp.three_test_table.stage import run_table, verify_table


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "verify"))
    parser.add_argument("--comparison-run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scratch-dir", type=Path, required=True)
    parser.add_argument("--code-commit")
    parser.add_argument("--multiplicity-run", type=Path)
    parser.add_argument("--label-support", type=Path)
    args = parser.parse_args(argv)
    if args.action == "run":
        record = run_table(
            args.comparison_run,
            args.output,
            scratch=args.scratch_dir,
            code_commit=args.code_commit,
            multiplicity=args.multiplicity_run,
            label_support=args.label_support,
        )
    else:
        if args.code_commit is not None:
            parser.error("verify reads the producer commit from the recorded manifest")
        record = verify_table(
            args.output,
            args.comparison_run,
            scratch=args.scratch_dir,
            multiplicity=args.multiplicity_run,
            label_support=args.label_support,
        )
    print(f"{record['run_id']}: {record['family_size']} traits; {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
