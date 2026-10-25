"""Every rerun and archive command refuses a temporary output root (#1049).

Each command's executor is replaced by a recorder, so nothing runs and nothing is
written; the output paths below are never created.
"""

from __future__ import annotations

import tempfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import pytest

import mtp.multiplicity.__main__ as multiplicity_cli
import mtp.workflow.__main__ as workflow_cli
import mtp.workflow.arms.__main__ as arms_cli

TEMPORARY = Path(tempfile.gettempdir()) / "mtp-1049-never-created" / "output"
DURABLE = Path("/srv/mtp-1049-never-created/output")
RECORD = SimpleNamespace(manifest_path="run.json")
INPUTS = [f"/in/{name}.json" for name in ("core", "noel", "jed", "ilyas")]


@dataclass(frozen=True)
class Entrypoint:
    recipe: str
    module: object
    main: Callable[[list[str]], object]
    executor: str
    argv: Callable[[Path], list[str]]
    result: object


def _arm(arm: str) -> Entrypoint:
    return Entrypoint(
        f"rerun-{arm}",
        arms_cli,
        arms_cli.main,
        "execute_arm",
        lambda output: [arm, "/in/run.json", str(output), "--no-figures"],
        RECORD,
    )


ENTRYPOINTS = [
    Entrypoint(
        "rerun-core",
        workflow_cli,
        workflow_cli.main,
        "execute_workflow",
        lambda output: ["prepare-core", str(output), "/in/raw.csv", "/in/e.tsv", "/in/n.tsv"],
        RECORD,
    ),
    Entrypoint(
        "core-run",
        workflow_cli,
        workflow_cli.main,
        "execute_workflow",
        lambda output: ["run", "--output-root", str(output)],
        RECORD,
    ),
    Entrypoint(
        "evidence-archive",
        workflow_cli,
        workflow_cli.main,
        "create_archive",
        lambda output: ["archive", "--output", str(output), "--consumer", "/in/record.json"],
        "archive.json",
    ),
    *(_arm(arm) for arm in ("ilyas", "noel", "jed")),
    Entrypoint(
        "rerun-compare",
        arms_cli,
        arms_cli.main,
        "execute_comparison_run",
        lambda output: ["compare", *INPUTS, str(output)],
        RECORD,
    ),
    Entrypoint(
        "rerun-multiplicity",
        multiplicity_cli,
        multiplicity_cli.main,
        "run_post_arm",
        lambda output: ["run", "--comparison-run", INPUTS[0], "--output", str(output)],
        SimpleNamespace(families=[], core_run_id="core"),
    ),
]


@pytest.fixture(params=ENTRYPOINTS, ids=[entrypoint.recipe for entrypoint in ENTRYPOINTS])
def entrypoint(request: pytest.FixtureRequest) -> Entrypoint:
    selected: Entrypoint = request.param
    return selected


@pytest.fixture
def executed(entrypoint: Entrypoint, monkeypatch: pytest.MonkeyPatch) -> list[tuple[object, ...]]:
    """The positional arguments of every call that reached the command's executor."""
    calls: list[tuple[object, ...]] = []

    def record(*args: object, **_: object) -> object:
        calls.append(args)
        return entrypoint.result

    monkeypatch.setattr(entrypoint.module, entrypoint.executor, record)
    return calls


def test_a_temporary_output_root_is_refused_before_any_work(
    entrypoint: Entrypoint,
    executed: list[tuple[object, ...]],
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as refused:
        entrypoint.main(entrypoint.argv(TEMPORARY))
    assert refused.value.code == 2
    error = capsys.readouterr().err
    assert "is under the temporary directory" in error
    assert "--allow-temporary-output" in error
    assert executed == []
    assert not TEMPORARY.exists()


def test_the_opt_in_admits_a_temporary_output_root(
    entrypoint: Entrypoint, executed: list[tuple[object, ...]]
) -> None:
    entrypoint.main([*entrypoint.argv(TEMPORARY), "--allow-temporary-output"])
    assert [TEMPORARY in call for call in executed] == [True]


def test_a_durable_output_root_needs_no_opt_in(
    entrypoint: Entrypoint, executed: list[tuple[object, ...]]
) -> None:
    entrypoint.main(entrypoint.argv(DURABLE))
    assert [DURABLE in call for call in executed] == [True]
