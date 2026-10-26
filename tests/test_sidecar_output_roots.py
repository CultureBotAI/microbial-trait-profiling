"""Run-bound sidecars refuse a temporary output root, as the reruns do (#1087).

Sidecars with their own CLI check their output root themselves. The core-* and trait-list
density recipes run modules that the arms also run as subprocesses, so their recipes check
the root through scripts/guard_output_root.sh and the modules stay unguarded (#1087).
"""

from __future__ import annotations

import importlib
import os
import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from unittest.mock import patch

import pytest

import mtp.workflow.output_roots as output_roots
from mtp import influence_extension, label_support_audit, presence_depth_sensitivity
from mtp.cellulose_decr import selection_scan, selection_scan_scopes, study_split_run
from mtp.paths import REPO_ROOT
from mtp.workflow.arms import ilyas, jed, noel
from mtp.workflow.arms.types import ArmContext, Command

with patch.object(sys, "path", [str(REPO_ROOT / "scripts"), *sys.path]):
    sankeys = importlib.import_module("create_six_shared_trait_sankeys")

DURABLE = Path("/srv/mtp-1087-never-created/output")
REFUSAL = "is under the temporary directory"
OPT_IN = "--allow-temporary-output"


class Reached(Exception):
    """Raised by the first step after the guard, so nothing runs or is written."""


@dataclass(frozen=True)
class Sidecar:
    recipe: str
    module: object
    main: Callable[[list[str]], object]
    first_step: str
    argv: Callable[[Path], list[str]]


SIDECARS = [
    Sidecar(
        "cellulose-decr-selection-scan",
        selection_scan,
        selection_scan.main,
        "open_tier1_2",
        lambda output: ["--core-run", "/in/run.json", "--output-dir", str(output)],
    ),
    Sidecar(
        "cellulose-decr-selection-scan-scopes",
        selection_scan_scopes,
        selection_scan_scopes.main,
        "load_protocol",
        lambda output: ["--scan", "/in/a", "--scan", "/in/b", "--output-dir", str(output)],
    ),
    Sidecar(
        "table1-influence-extension",
        influence_extension,
        influence_extension.main,
        "run_extension",
        lambda output: ["--core-run", "/in/run.json", "--output-dir", str(output)],
    ),
    Sidecar(
        "audit-label-support",
        label_support_audit,
        label_support_audit.main,
        "resolve_kg_file",
        lambda output: [
            *("--taxa", "/in/taxa.tsv", "--matrix", "/in/matrix.tsv"),
            *("--profiles", "/in/profiles.tsv", "--strain-rollup", "off"),
            *("--output-dir", str(output)),
        ],
    ),
    Sidecar(
        "selection-study-splits",
        study_split_run,
        study_split_run.main,
        "execute",
        lambda output: ["--comparison-run", "/in/comparison_run.json", "--output-dir", str(output)],
    ),
    Sidecar(
        "presence-depth-sensitivity",
        presence_depth_sensitivity,
        presence_depth_sensitivity.main,
        "run_sensitivity",
        lambda output: ["--core-run", "/in/run.json", "--output-dir", str(output)],
    ),
    Sidecar(
        "trait-list-direction-contrasts",
        sankeys,
        sankeys.cli,
        "load_trait_list",
        lambda output: ["--trait-list", "/in/traits.txt", "--output-dir", str(output)],
    ),
]

# Recipes whose module the arms also run; the recipe, not the module, checks the root
GUARDED_RECIPES = {
    "core-screen": ["/in/run.json"],
    "core-composition": ["/in/run.json"],
    "core-ora": ["/in/run.json"],
    "core-strata": ["/in/run.json"],
    "density-adjustment-trait-list": ["/in/traits.txt", "/in/profiles.tsv", "/in/columns.tsv"],
}


@pytest.fixture
def temporary(tmp_path: Path) -> Path:
    """An output root in pytest's scratch, which lies in the system temporary directory."""
    assert output_roots.temporary_root(tmp_path) is not None
    return tmp_path / "output"


@pytest.fixture(params=SIDECARS, ids=[sidecar.recipe for sidecar in SIDECARS])
def sidecar(request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch) -> Sidecar:
    selected: Sidecar = request.param

    def reached(*_: object, **__: object) -> object:
        raise Reached

    monkeypatch.setattr(selected.module, selected.first_step, reached)
    return selected


def test_a_sidecar_refuses_a_temporary_output_root_before_any_work(
    sidecar: Sidecar, temporary: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as refused:
        sidecar.main(sidecar.argv(temporary))
    assert refused.value.code == 2
    error = capsys.readouterr().err
    assert REFUSAL in error
    assert OPT_IN in error
    assert not temporary.exists()


def test_the_opt_in_admits_a_sidecar_temporary_root(sidecar: Sidecar, temporary: Path) -> None:
    with pytest.raises(Reached):
        sidecar.main([*sidecar.argv(temporary), OPT_IN])


def test_a_durable_sidecar_root_needs_no_opt_in(sidecar: Sidecar) -> None:
    with pytest.raises(Reached):
        sidecar.main(sidecar.argv(DURABLE))


NO_WRITE_MODES = [
    Sidecar(
        "selection-scan-plan",
        selection_scan,
        selection_scan.main,
        "open_tier1_2",
        lambda output: ["--core-run", "/in/run.json", "--output-dir", str(output), "--plan"],
    ),
    Sidecar(
        "study-splits-verify",
        study_split_run,
        study_split_run.main,
        "execute",
        lambda output: ["--comparison-run", "/in/c.json", "--output-dir", str(output), "--verify"],
    ),
]


@pytest.mark.parametrize("mode", NO_WRITE_MODES, ids=[mode.recipe for mode in NO_WRITE_MODES])
def test_modes_that_write_nothing_are_not_refused(
    mode: Sidecar, temporary: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def reached(*_: object, **__: object) -> object:
        raise Reached

    monkeypatch.setattr(mode.module, mode.first_step, reached)
    with pytest.raises(Reached):
        mode.main(mode.argv(temporary))


def test_the_trait_list_opt_in_needs_a_trait_list(capsys: pytest.CaptureFixture[str]) -> None:
    with patch.object(sankeys, "main") as six_trait_main, pytest.raises(SystemExit):
        sankeys.cli([OPT_IN])
    six_trait_main.assert_not_called()
    assert "apply only with --trait-list" in capsys.readouterr().err


def test_the_check_command_refuses_and_admits(
    temporary: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as refused:
        output_roots.main([str(temporary)])
    assert refused.value.code == 2
    assert REFUSAL in capsys.readouterr().err
    output_roots.main([str(temporary), OPT_IN])
    output_roots.main([str(DURABLE)])
    assert not temporary.exists()


def _guard(*argv: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["bash", str(REPO_ROOT / "scripts" / "guard_output_root.sh"), *argv],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def test_the_recipe_guard_consumes_the_opt_in_and_runs_the_rest(temporary: Path) -> None:
    refused = _guard(str(temporary), "printf", "%s\\n", "a", "b")
    assert refused.returncode == 2
    assert REFUSAL in refused.stderr
    assert refused.stdout == ""

    admitted = _guard(str(temporary), "printf", "%s\\n", "a", OPT_IN, "b")
    assert admitted.returncode == 0, admitted.stderr
    assert admitted.stdout == "a\nb\n"

    durable = _guard(str(DURABLE), "printf", "%s\\n", "a")
    assert (durable.returncode, durable.stdout) == (0, "a\n")


def _recipe(recipe: str, output: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["just", recipe, *GUARDED_RECIPES[recipe], str(output), *extra],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.mark.parametrize("recipe", GUARDED_RECIPES)
def test_an_arm_module_recipe_refuses_a_temporary_root_before_its_module_runs(
    recipe: str, temporary: Path
) -> None:
    completed = _recipe(recipe, temporary)
    assert completed.returncode != 0
    assert REFUSAL in completed.stderr
    assert "Traceback" not in completed.stderr
    assert not temporary.exists()


@pytest.mark.parametrize("recipe", GUARDED_RECIPES)
def test_an_arm_module_recipe_passes_on_everything_but_the_opt_in(
    recipe: str, temporary: Path
) -> None:
    """With the opt-in the module runs, and fails only on the absent inputs."""
    extra = ("upper", OPT_IN) if recipe == "core-ora" else (OPT_IN,)
    completed = _recipe(recipe, temporary, *extra)
    assert completed.returncode != 0
    assert REFUSAL not in completed.stderr
    assert "unrecognized arguments" not in completed.stderr
    assert "/in/" in completed.stderr


def test_core_ora_defaults_to_the_upper_tail(tmp_path: Path) -> None:
    """Without a tail, `just core-ora` runs enrichment, not depletion (#1111).

    A stub `uv` first on PATH prints each command the recipe would run, one argument per line.
    """
    stub = tmp_path / "bin" / "uv"
    stub.parent.mkdir()
    stub.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "$@"\n', encoding="utf-8")
    stub.chmod(0o755)
    environment = {**os.environ, "PATH": f"{stub.parent}{os.pathsep}{os.environ['PATH']}"}

    def tail(*extra: str) -> str:
        completed = subprocess.run(
            ["just", "core-ora", "/in/run.json", str(DURABLE), *extra],
            cwd=REPO_ROOT,
            env=environment,
            capture_output=True,
            text=True,
            check=True,
        )
        arguments = completed.stdout.splitlines()
        assert "scripts/trait_enrichment.py" in arguments
        return arguments[arguments.index("--tail") + 1]

    assert tail() == "upper"
    assert tail("lower") == "lower"


class _Products(dict[str, Path]):
    def __missing__(self, key: str) -> Path:
        return Path("/core") / f"{key}.tsv"


@pytest.mark.parametrize("figures", [True, False])
@pytest.mark.parametrize(
    "build",
    [jed.build_commands, noel.build_commands, ilyas.build_commands],
    ids=["jed", "noel", "ilyas"],
)
def test_no_arm_plan_hands_the_opt_in_to_a_child_module(
    build: Callable[[ArmContext], list[Command]], figures: bool
) -> None:
    """The arm checks its own directory (#1049); its recorded child argv stay as they were."""
    context = ArmContext(
        repository=REPO_ROOT,
        manifest=Path("/core/run.json"),
        output=Path("/arm"),
        inputs=Path("/arm/inputs"),
        products=_Products(),
        code_commit="revision",
        replicates=None,
        seed=None,
        figures=figures,
    )
    argv = [token for command in build(context) for token in command.argv]
    assert argv
    assert OPT_IN not in argv
