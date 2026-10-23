"""Reproduce missing direction-option verification (#843/#944, 2026-10-23).

Inputs: an existing synthetic core fixture and the reviewed runtime checkout.
Outputs: a fixture-only arm record, tampered-plan evidence and parser exit status.
No production inputs or scientific analysis are used.
"""

from __future__ import annotations

import argparse
import io
import json
import runpy
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from mtp.paths import REPO_ROOT
from mtp.workflow.arms.runtime import execute_arm, verify_arm
from mtp.workflow.arms.types import ArmContext, Command
from mtp.workflow.runtime import content_hash
from mtp.workflow.verification import verify_core_run


def parser_exit(argv: list[str]) -> tuple[int, str]:
    """Exercise an entrypoint's required-argument rejection in this process."""
    previous_argv = sys.argv
    stderr, stdout = io.StringIO(), io.StringIO()
    try:
        sys.argv = argv
        with redirect_stderr(stderr), redirect_stdout(stdout):
            try:
                runpy.run_path(argv[0], run_name="__main__")
            except SystemExit as error:
                code = error.code if isinstance(error.code, int) else int(error.code is not None)
            else:
                code = 0
    finally:
        sys.argv = previous_argv
    return code, stderr.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-fixture", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--known-heatmap",
        action="store_true",
        help="Probe the registered heatmap entrypoint in the recorded plan.",
    )
    args = parser.parse_args()
    core = verify_core_run(args.core_fixture)
    if not any(p.name == "source_mode" and p.value == "fixture" for p in core.parameters or []):
        raise ValueError("This reproduction accepts only an explicitly synthetic fixture")

    def plan(context: ArmContext) -> list[Command]:
        source = context.products["count_profiles"]
        destination = context.output / "observed.tsv"
        return [
            Command(
                "observe_direction",
                (
                    sys.executable,
                    str(REPO_ROOT / "scripts/adhoc/issue843_direction_test_command.py"),
                    "--input",
                    str(source),
                    "--output",
                    str(destination),
                    "--direction-basis",
                    context.direction_basis,
                ),
                (source,),
                (destination,),
            )
        ]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    record = execute_arm("jed", args.core_fixture, args.output_dir / "arm", plan, figures=False)
    path = Path(record.manifest_path)
    changed = json.loads(path.read_text(encoding="utf-8"))
    changed["plan"][0]["argv"] = changed["plan"][0]["argv"][:-2]
    fixture_argv = list(changed["plan"][0]["argv"])
    if args.known_heatmap:
        # Verification probe only: the synthetic copy operation remains the sole
        # executed command. The actual heatmap parser is checked separately below.
        changed["plan"][0]["argv"][1] = str(REPO_ROOT / "scripts/create_taxa_guild_heatmap.py")
    changed["plan_sha256"] = content_hash(changed["plan"])
    path.write_text(json.dumps(changed), encoding="utf-8")
    rejection: str | None
    try:
        verified = verify_arm(path)
    except ValueError as error:
        accepted, rejection = False, str(error)
    else:
        accepted, rejection = verified.id == record.id, None
    execution_exit, execution_error = parser_exit(fixture_argv[1:])
    actual_exit, actual_error = parser_exit(
        [str(REPO_ROOT / "scripts/create_taxa_guild_heatmap.py")]
    )
    evidence = {
        "tampered_plan_verified": accepted,
        "registered_heatmap_probe": args.known_heatmap,
        "verification_rejection": rejection,
        "fixture_entrypoint_exit": execution_exit,
        "fixture_entrypoint_stderr": execution_error,
        "actual_heatmap_parser_exit": actual_exit,
        "actual_heatmap_parser_requires_direction": "--direction-basis" in actual_error,
    }
    (args.output_dir / "evidence.json").write_text(
        json.dumps(evidence, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(evidence, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
