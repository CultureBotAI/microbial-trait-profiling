"""Exercise recorded direction arguments with fixture-only files (#843, 2026-10-23).

Inputs: an explicit fixture profile and direction basis. Output: an explicit copy
used by execution/verification tests; no scientific analysis is performed.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def build_parser(
    parser_class: type[argparse.ArgumentParser] = argparse.ArgumentParser,
) -> argparse.ArgumentParser:
    """Supply the same fixture CLI contract to the runtime's targeted tests."""
    parser = parser_class(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--direction-basis",
        choices=("group_1_vs_group_0", "disease_vs_comparison"),
        required=True,
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.output.write_bytes(args.input.read_bytes())


if __name__ == "__main__":
    main()
