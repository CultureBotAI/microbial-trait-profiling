#!/usr/bin/env python3
"""Validate a taxon or profile count matrix against the shared count contract."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from bsdbng.count_inputs import (
    COMPLETE_COUNT_MODE,
    COUNT_INPUT_MODES,
    profile_count_reader,
    read_profile_coverage,
    taxa_trait_count_reader,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("profile", "taxa"))
    parser.add_argument("input", type=Path)
    parser.add_argument(
        "--count-input-mode",
        choices=COUNT_INPUT_MODES,
        default=COMPLETE_COUNT_MODE,
    )
    parser.add_argument(
        "--coverage-input",
        type=Path,
        help="Coverage sidecar to reconcile with a profile matrix.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    reader = (
        profile_count_reader(args.input, mode=args.count_input_mode)
        if args.kind == "profile"
        else taxa_trait_count_reader(args.input, mode=args.count_input_mode)
    )
    keys = {row.key for row in reader}
    if args.coverage_input is not None:
        if args.kind != "profile":
            raise ValueError("--coverage-input is valid only for a profile matrix")
        coverage = read_profile_coverage(args.coverage_input)
        missing = sorted(keys.difference(coverage))
        extra = sorted(set(coverage).difference(keys))
        if missing or extra:
            raise ValueError(
                f"Profile/coverage key mismatch: missing={missing[:5]}, extra={extra[:5]}"
            )
        print(f"Coverage rows: {len(coverage)}")
        print(f"Incomplete coverage rows: {sum(not row.is_complete for row in coverage.values())}")
    print(f"Validated {args.kind} count input: {args.input}")
    print(f"Rows: {reader.summary.rows}")
    print(f"Count columns: {reader.summary.count_columns}")
    print(f"Count input mode: {reader.summary.mode}")
    print(f"Blank counts interpreted as zero: {reader.summary.blank_counts_interpreted_as_zero}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
