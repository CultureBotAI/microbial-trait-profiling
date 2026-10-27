"""Compare pinned off and phenotype evidence on a common population and producer."""

from __future__ import annotations

import argparse
from pathlib import Path

from mtp.rollup_comparison.stage import run_comparison


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for scope in ("off", "phenotype"):
        parser.add_argument(f"--{scope}-three-test-run", type=Path, required=True)
        parser.add_argument(f"--{scope}-sha256", required=True)
        parser.add_argument(f"--{scope}-comparison-run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--code-commit")
    args = parser.parse_args(argv)
    record = run_comparison(
        args.off_three_test_run,
        args.off_comparison_run,
        args.phenotype_three_test_run,
        args.phenotype_comparison_run,
        args.output,
        off_sha256=args.off_sha256,
        phenotype_sha256=args.phenotype_sha256,
        code_commit=args.code_commit,
    )
    print(f"{record['run_id']}: {record['counts']['shared_traits']} shared traits; {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
