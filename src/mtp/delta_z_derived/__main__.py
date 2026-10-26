"""Derive labelled summaries from pinned, completed Delta-z tables without recomputing them."""

from __future__ import annotations

import argparse
from pathlib import Path

from mtp.delta_z_derived.stage import run_derived


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--three-test-run", type=Path, required=True)
    parser.add_argument("--three-test-sha256", required=True)
    parser.add_argument("--comparison-run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--trait-list", type=Path)
    selection.add_argument("--summary-only", action="store_true")
    parser.add_argument("--top-diseases", type=int, default=7)
    parser.add_argument("--code-commit")
    args = parser.parse_args(argv)
    if args.top_diseases < 1:
        parser.error("--top-diseases must be positive")
    record = run_derived(
        args.three_test_run,
        args.comparison_run,
        args.output,
        expected_sha256=args.three_test_sha256,
        trait_list=args.trait_list,
        summary_only=args.summary_only,
        code_commit=args.code_commit,
        top_diseases=args.top_diseases,
    )
    print(f"{record['run_id']}: {record['family_size']} traits; {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
