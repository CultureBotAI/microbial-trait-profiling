"""Prove #839 rendered geometry rejects undersized layout constants (2026-10-23).

Inputs: this checkout's guild renderer and regression test, synthetic profiles only.
Outputs: baseline/mutant diagnostic PNGs and JSON verdicts in --output-dir.
Constants are changed in memory; source files and scientific outputs are untouched.
"""

from __future__ import annotations

import argparse
import json
import runpy
from collections.abc import Callable
from pathlib import Path
from typing import cast

import pytest

from mtp import guild_maps

Check = Callable[[Path, pytest.MonkeyPatch, str | None, bool], None]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    test_path = Path(__file__).resolve().parents[2] / "tests" / "test_trait_guilds.py"
    outcomes: list[dict[str, object]] = []
    for constant in (None, "TITLE_BLOCK_HEIGHT", "SUBTITLE_LINE_HEIGHT"):
        name = constant or "baseline"
        destination = args.output_dir / name
        destination.mkdir(exist_ok=True)
        with pytest.MonkeyPatch.context() as patch:
            if constant:
                patch.setattr(guild_maps, constant, 1)
            # Load after mutation: an old formula-repeating assertion would import
            # the same shrunken value, and would still pass.
            namespace = runpy.run_path(str(test_path))
            check = cast(
                Check, namespace["test_the_plotting_area_starts_below_the_last_subtitle_line"]
            )
            try:
                check(destination, patch, "Population note.<br>A second population line.", True)
            except AssertionError as error:
                outcomes.append({"case": name, "rejected": True, "error": str(error)})
            else:
                outcomes.append({"case": name, "rejected": False})
    args.output_dir.joinpath("verdicts.json").write_text(
        json.dumps(outcomes, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(outcomes, indent=2))
    return 0 if [row["rejected"] for row in outcomes] == [False, True, True] else 1


if __name__ == "__main__":
    raise SystemExit(main())
