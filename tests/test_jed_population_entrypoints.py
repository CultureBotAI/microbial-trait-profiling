"""Run real CLI parsers and the CI-compatible Tier 1+2 recipe (#1344, port #1219)."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from mtp.workflow.__main__ import build_parser as core_parser
from mtp.workflow.arms.__main__ import build_parser as arms_parser

REPOSITORY = Path(__file__).resolve().parents[1]
POPULATION = "disease_vs_healthy_tier1_2_subject_controls"


@pytest.mark.parametrize("arm", ["jed", "noel", "ilyas"])
def test_explicit_figures_retains_the_default_and_refuses_conflicting_switches(arm: str) -> None:
    parser = arms_parser()
    argv = [arm, "core/run.json", "fresh-output"]
    assert parser.parse_args(argv).no_figures is False
    assert parser.parse_args([*argv, "--figures"]).no_figures is False
    assert parser.parse_args([*argv, "--no-figures"]).no_figures is True
    with pytest.raises(SystemExit):
        parser.parse_args([*argv, "--figures", "--no-figures"])


@pytest.mark.parametrize("extra", [[], ["off", "fixture"]])
def test_tier12_recipe_executes_with_explicit_source_paths(
    tmp_path: Path, extra: list[str]
) -> None:
    """CI runs this with just 1.39.0; the shim captures argv without a scientific run."""
    just = shutil.which("just")
    assert just is not None, "just is a required development/CI tool"
    shim = tmp_path / "bin" / "uv"
    shim.parent.mkdir()
    shim.write_text(
        f"#!{sys.executable}\nimport json, sys\nprint(json.dumps(sys.argv[1:]))\n",
        encoding="utf-8",
    )
    shim.chmod(0o755)
    inputs = [
        str(tmp_path / name) for name in ("fresh core", "raw dump.csv", "edges.tsv", "nodes.tsv")
    ]
    result = subprocess.run(
        [just, "--justfile", str(REPOSITORY / "Justfile"), "rerun-core-548", *inputs, *extra],
        env={**os.environ, "PATH": f"{shim.parent}{os.pathsep}{os.environ['PATH']}"},
        capture_output=True,
        text=True,
        check=True,
    )
    captured = json.loads(result.stdout)
    assert captured[:4] == ["run", "python", "-m", "mtp.workflow"]
    parsed = core_parser().parse_args(captured[4:])
    assert parsed.command == "prepare-core"
    assert [
        str(parsed.output_root),
        str(parsed.raw),
        str(parsed.edges),
        str(parsed.nodes),
    ] == inputs
    assert parsed.contrast_population == POPULATION
    assert parsed.strain_rollup == ("off" if extra else "phenotype")
    assert parsed.source_mode == ("fixture" if extra else "pinned")
    assert not any(Path(path).exists() for path in inputs)
