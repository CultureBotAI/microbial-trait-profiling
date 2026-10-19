"""Public rerun interfaces must retain explicit paths and options (#658)."""

import subprocess
from pathlib import Path

import pytest

from mtp.paths import REPO_ROOT
from mtp.workflow.arms.__main__ import build_parser


@pytest.mark.parametrize("arm", ["ilyas", "noel", "jed"])
def test_owner_cli_requires_core_manifest_and_destination(arm: str) -> None:
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([arm])
    args = parser.parse_args(
        [
            arm,
            "/input with spaces/run.json",
            "/fresh output",
            "--no-figures",
            "--replicates",
            "2",
            "--seed",
            "0",
        ]
    )
    assert args.core_run == Path("/input with spaces/run.json")
    assert args.output == Path("/fresh output")
    assert args.no_figures
    assert args.replicates == 2
    assert args.seed == 0


@pytest.mark.parametrize("recipe", ["rerun-ilyas", "rerun-noel", "rerun-jed"])
def test_just_owner_recipes_have_no_legacy_rebuild_dependencies(recipe: str) -> None:
    completed = subprocess.run(
        ["just", "--dry-run", recipe, "/input/run.json", "/new/output", "--no-figures"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    commands = completed.stdout + completed.stderr
    assert commands.count("uv run") == 1
    assert "mtp.workflow.arms" in commands
    assert '"$@"' in commands


@pytest.mark.parametrize(
    "script",
    [
        "rerun_enrichment_ora_arm_ilyas.sh",
        "rerun_directional_consistency_arm_noel.sh",
        "rerun_composition_arm_jed.sh",
    ],
)
def test_shell_entrypoints_fail_closed_without_explicit_inputs(script: str) -> None:
    completed = subprocess.run(
        ["bash", str(REPO_ROOT / "scripts" / script)], capture_output=True, text=True, check=False
    )
    assert completed.returncode == 2
    assert "CORE_RUN_JSON FRESH_OUTPUT_DIR" in completed.stderr
