"""Render publication figures from pinned, completed sidecars without replaying analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

from mtp.publication_figures.stage import FIGURES, build_figures
from mtp.workflow.output_roots import add_temporary_output_option


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--three-test-run", type=Path, required=True)
    parser.add_argument("--three-test-sha256", required=True)
    parser.add_argument("--comparison-run", type=Path, required=True)
    parser.add_argument("--derived-run", type=Path, required=True)
    parser.add_argument("--derived-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure", choices=FIGURES, action="append", dest="figures")
    parser.add_argument(
        "--ranking", choices=("comparison_mean", "study_mean"), default="comparison_mean"
    )
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--code-commit")
    add_temporary_output_option(parser)
    args = parser.parse_args(argv)
    record = build_figures(
        args.three_test_run,
        args.comparison_run,
        args.derived_run,
        args.output,
        three_test_sha256=args.three_test_sha256,
        derived_sha256=args.derived_sha256,
        figures=tuple(args.figures) if args.figures else FIGURES,
        ranking=args.ranking,
        top_n=args.top_n,
        dpi=args.dpi,
        code_commit=args.code_commit,
        allow_temporary_output=args.allow_temporary_output,
    )
    print(f"{len(record['figures'])} figures published: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
