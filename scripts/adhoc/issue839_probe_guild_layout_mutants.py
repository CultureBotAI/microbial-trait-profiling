"""Probe #839/#929 layout regressions against real rendered tests (2026-10-23).

Inputs: this checkout's guild renderer and regression test, synthetic profiles only.
Outputs: baseline/mutant diagnostic PNGs and JSON verdicts in --output-dir.
Constants are changed in memory; source files and scientific outputs are untouched.
"""

from __future__ import annotations

import argparse
import io
import json
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest

from mtp import guild_maps


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tests = Path(__file__).resolve().parents[2] / "tests"
    outcomes: list[dict[str, object]] = []
    for constant in (
        None,
        "TITLE_BLOCK_HEIGHT",
        "SUBTITLE_LINE_HEIGHT",
        "SITE_HEADER_HEIGHT",
        "entry_counting",
    ):
        name = constant or "baseline"
        destination = args.output_dir / name
        destination.mkdir(exist_ok=True)
        with pytest.MonkeyPatch.context() as patch:
            if constant == "entry_counting":
                patch.setattr(
                    guild_maps,
                    "_title_block_height",
                    lambda lines: (
                        guild_maps.TITLE_BLOCK_HEIGHT + guild_maps.SUBTITLE_LINE_HEIGHT * len(lines)
                    ),
                )
            elif constant:
                patch.setattr(guild_maps, constant, getattr(guild_maps, constant) // 2)
            # Import each test after mutation, just as a fresh checkout run would.
            # Retained imported constants would spuriously reject a formula-only test.
            for module in ("test_trait_guilds", "test_table1_trait_maps"):
                sys.modules.pop(module, None)
            captured = io.StringIO()
            with redirect_stdout(captured), redirect_stderr(captured):
                code = pytest.main(
                    [
                        "-q",
                        "--basetemp",
                        str(destination / "fixtures"),
                        str(tests / "test_trait_guilds.py")
                        + "::test_the_plotting_area_starts_below_the_last_subtitle_line",
                        str(tests / "test_table1_trait_maps.py")
                        + "::test_the_site_names_sit_below_the_last_subtitle_line",
                    ]
                )
            destination.joinpath("pytest.log").write_text(captured.getvalue(), encoding="utf-8")
            if code not in (pytest.ExitCode.OK, pytest.ExitCode.TESTS_FAILED):
                raise RuntimeError(f"{name}: pytest infrastructure failed with exit {code}")
            outcomes.append({"case": name, "rejected": code == pytest.ExitCode.TESTS_FAILED})
    args.output_dir.joinpath("verdicts.json").write_text(
        json.dumps(outcomes, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(outcomes, indent=2))
    return 0 if [row["rejected"] for row in outcomes] == [False, True, True, True, True] else 1


if __name__ == "__main__":
    raise SystemExit(main())
