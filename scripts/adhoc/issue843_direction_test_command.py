"""Exercise recorded direction arguments with fixture-only files (#843, 2026-10-23).

Inputs: an explicit fixture profile and direction basis. Output: an explicit copy
used by execution/verification tests; no scientific analysis is performed.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--direction-basis",
        choices=("group_1_vs_group_0", "disease_vs_comparison"),
        required=True,
    )
    args = parser.parse_args()
    args.output.write_bytes(args.input.read_bytes())


if __name__ == "__main__":
    main()
