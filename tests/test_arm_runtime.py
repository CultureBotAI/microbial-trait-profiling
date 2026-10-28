"""Real subprocess and corruption tests for LinkML arm execution records (#658)."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from dataclasses import dataclass, replace
from pathlib import Path

import pytest

from mtp.datamodel.workflow_pydantic import ArmRun
from mtp.paths import REPO_ROOT
from mtp.workflow import runtime as core_runtime
from mtp.workflow.arms import runtime
from mtp.workflow.arms.runtime import ArmExecutionError, Planner, execute_arm, verify_arm
from mtp.workflow.arms.types import ArmContext, Command
from mtp.workflow.operations import OPERATIONS, validate_core_product

HELPER = Path("scripts/adhoc/issue658_arm_test_command.py")


@dataclass(frozen=True)
class Harness:
    repository: Path
    manifest: Path
    output: Path

    def execute(self, planner: Planner) -> ArmRun:
        return execute_arm(
            "jed", self.manifest, self.output, planner, figures=False, repository=self.repository
        )


@pytest.fixture
def harness(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Harness:
    repository = tmp_path / "source"
    (repository / "src").mkdir(parents=True)
    (repository / "src" / "version.py").write_text("VERSION = 1\n", encoding="utf-8")
    # The core uses the real imported operations. Retain their matrix/policy source
    # bytes under the expected inventory paths instead of inventing a producer (#1129).
    for name in ("taxa_trait_matrix.py", "trait_key_policy.py"):
        relative = Path("src/mtp") / name
        (repository / relative).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO_ROOT / relative, repository / relative)
    for name in ("Justfile", "pyproject.toml", "uv.lock"):
        (repository / name).write_text("# test fixture\n", encoding="utf-8")
    (repository / HELPER).parent.mkdir(parents=True)
    shutil.copyfile(REPO_ROOT / HELPER, repository / HELPER)

    def git_state(root: Path) -> tuple[str, bool]:
        assert root == repository
        return "fixture-revision", False

    monkeypatch.setattr(runtime, "repository_state", git_state)
    monkeypatch.setattr(core_runtime, "repository_state", git_state)
    fixture = REPO_ROOT / "tests" / "fixtures" / "workflow"
    core = core_runtime.execute_workflow(
        REPO_ROOT / "workflows" / "bugsigdb_core.yaml",
        {
            "raw_dump": fixture / "full_dump.csv",
            "kg_edges": fixture / "merged-kg_edges.tsv",
            "kg_nodes": fixture / "merged-kg_nodes.tsv",
        },
        tmp_path / "core",
        OPERATIONS,
        parameter_overrides={"source_mode": "fixture"},
        validator=validate_core_product,
        repository_root=repository,
    )
    return Harness(repository, Path(core.manifest_path), tmp_path / "arm outputs")


def _command(
    context: ArmContext,
    *,
    name: str = "copy",
    source: Path | None = None,
    output: Path | None = None,
    mode: str = "copy",
    extra: tuple[str, ...] = (),
) -> Command:
    source = source or context.products["count_profiles"]
    output = output or context.output / "result.txt"
    return Command(
        name,
        (
            sys.executable,
            str(context.repository / HELPER),
            "--input",
            str(source),
            "--output",
            str(output),
            "--mode",
            mode,
            *extra,
        ),
        (source,),
        (output,),
    )


def _copy_plan(context: ArmContext) -> list[Command]:
    return [_command(context)]


def test_successful_real_execution_is_directly_linkml_validated(harness: Harness) -> None:
    record = harness.execute(_copy_plan)
    manifest = Path(record.manifest_path)
    assert verify_arm(manifest) == record
    validated = subprocess.run(
        [
            sys.executable,
            "-m",
            "linkml.validator.cli",
            "--schema",
            str(REPO_ROOT / "schema" / "workflow.yaml"),
            "--target-class",
            "ArmRun",
            str(manifest),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert validated.returncode == 0, validated.stdout + validated.stderr
    assert record.steps
    assert record.steps[0].exit_code == 0
    assert record.steps[0].log
    assert Path(record.steps[0].log.path).is_file()
    assert record.steps[0].inputs
    assert (harness.output / "result.txt").read_bytes() == Path(
        record.steps[0].inputs[0].path
    ).read_bytes()
    assert json.loads(record.runtime_json)["process_settings"]["MPLBACKEND"] == "Agg"


@pytest.mark.parametrize("mode", ["fail", "missing", "mutate_input", "symlink", "hardlink"])
def test_failed_commands_never_verify_and_keep_diagnostics(harness: Harness, mode: str) -> None:
    def plan(context: ArmContext) -> list[Command]:
        return [
            _command(context, source=context.inputs / "study_trait_profiles_wide.tsv", mode=mode)
        ]

    with pytest.raises(ArmExecutionError) as failure:
        harness.execute(plan)
    record = failure.value.record
    assert record.status == "failed"
    assert record.steps
    assert record.steps[0].status == "failed"
    saved = ArmRun.model_validate_json(Path(record.manifest_path).read_text(encoding="utf-8"))
    assert saved.error
    assert saved.ended_at
    with pytest.raises(ValueError, match="not a completed"):
        verify_arm(Path(record.manifest_path))


def test_code_mutation_fails_and_keeps_the_original_code_inventory(harness: Harness) -> None:
    def plan(context: ArmContext) -> list[Command]:
        return [
            _command(
                context,
                extra=("--mutate-file", str(context.repository / "src" / "version.py")),
            )
        ]

    with pytest.raises(ArmExecutionError, match="Code or configuration changed"):
        harness.execute(plan)


def test_collection_input_additions_fail_when_not_declared_outputs(harness: Harness) -> None:
    def plan(context: ArmContext) -> list[Command]:
        return [_command(context, source=context.inputs, mode="add_input")]

    with pytest.raises(ArmExecutionError, match="Recorded bytes changed"):
        harness.execute(plan)


def _shared_directory_plan(context: ArmContext) -> list[Command]:
    directory = context.output / "shared"
    first = _command(context, output=directory / "first.txt")
    second = _command(context, name="extend", source=directory, output=directory / "second.txt")
    return [
        replace(first, outputs=(directory, *first.outputs)),
        replace(second, outputs=(directory, *second.outputs)),
    ]


def test_shared_output_directory_extensions_preserve_stage_snapshots(harness: Harness) -> None:
    record = harness.execute(_shared_directory_plan)
    assert verify_arm(Path(record.manifest_path)) == record
    assert record.steps
    assert record.steps[0].outputs
    assert record.steps[1].outputs
    assert len(record.steps[0].outputs[0].inventory) == 1
    assert len(record.steps[1].outputs[0].inventory) == 2
    (harness.output / "shared" / "unexpected.txt").write_text("extra", encoding="utf-8")
    with pytest.raises(ValueError, match="output inventory changed"):
        verify_arm(Path(record.manifest_path))


def test_rewriting_an_earlier_stage_output_fails_the_attempt(harness: Harness) -> None:
    def plan(context: ArmContext) -> list[Command]:
        first = _command(context, name="first", output=context.output / "first.txt")
        second = _command(
            context,
            name="second",
            output=context.output / "second.txt",
            extra=("--mutate-file", str(context.output / "first.txt")),
        )
        return [first, second]

    with pytest.raises(ArmExecutionError, match="Recorded bytes changed"):
        harness.execute(plan)


def test_argument_vector_preserves_shell_metacharacters_as_literal_data(harness: Harness) -> None:
    marker = harness.output.parent / "shell-was-executed"
    literal = f"$(touch {marker}); `touch {marker}` && 'quoted'"

    def plan(context: ArmContext) -> list[Command]:
        return [_command(context, extra=("--payload", literal))]

    record = harness.execute(plan)
    assert (harness.output / "result.txt").read_text(encoding="utf-8").endswith(literal)
    assert not marker.exists()
    verify_arm(Path(record.manifest_path))


@pytest.mark.parametrize("tampering", ["argv", "missing_step", "output_path", "log_path"])
def test_manifest_plan_and_execution_tampering_are_rejected(
    harness: Harness, tampering: str
) -> None:
    record = harness.execute(_copy_plan)
    path = Path(record.manifest_path)
    changed = json.loads(path.read_text(encoding="utf-8"))
    if tampering == "argv":
        changed["plan"][0]["argv"].append("--fictional-option")
    elif tampering == "missing_step":
        changed["steps"] = []
    elif tampering == "output_path":
        changed["steps"][0]["outputs"][0]["path"] = str(harness.output / "other.txt")
    else:
        changed["steps"][0]["log"]["path"] = str(harness.output / "result.txt")
    path.write_text(json.dumps(changed), encoding="utf-8")
    with pytest.raises(
        ValueError, match=r"plan checksum|ordered plan|evidence disagrees|log location"
    ):
        verify_arm(path)


@pytest.mark.parametrize("collision", ["source", "record", "inputs", "logs", "symlink"])
def test_invalid_plan_outputs_fail_before_starting_commands(
    harness: Harness, collision: str
) -> None:
    targets = {
        "source": harness.repository / "src" / "version.py",
        "record": harness.output / "arm_run.json",
        "inputs": harness.output / "inputs" / "full_dump.csv",
        "logs": harness.output / "logs" / "copy.log",
        "symlink": harness.output.parent / "alias" / "result.txt",
    }
    if collision == "symlink":
        targets[collision].parent.symlink_to(harness.repository / "src", target_is_directory=True)

    def plan(context: ArmContext) -> list[Command]:
        return [_command(context, output=targets[collision])]

    with pytest.raises(ValueError, match=r"outside dedicated|execution record|reserved|aliases"):
        harness.execute(plan)
    assert not harness.output.exists()


def test_unbound_and_future_inputs_are_rejected(harness: Harness) -> None:
    def plan(context: ArmContext) -> list[Command]:
        return [_command(context, source=context.output / "future-result.txt")]

    with pytest.raises(ValueError, match="Input is not core"):
        harness.execute(plan)
    assert not harness.output.exists()
