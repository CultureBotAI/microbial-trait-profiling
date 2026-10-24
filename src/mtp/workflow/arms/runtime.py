"""Execute declared owner commands and persist checked LinkML records (#658)."""

from __future__ import annotations

import argparse
import asyncio
import importlib.metadata
import json
import os
import platform
import re
import runpy
import shutil
import sys
import uuid
from collections.abc import Callable, Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import NoReturn

from mtp import taxon_contributions
from mtp.datamodel.workflow_pydantic import (
    ArmCommandDefinition,
    ArmCommandExecution,
    ArmRun,
    ContrastPopulation,
    DirectionBasis,
    ExecutionStatus,
    FileInventoryEntry,
    ParameterValue,
    RecordedPath,
    StorageKind,
)
from mtp.paths import REPO_ROOT
from mtp.workflow.arms.types import ArmContext, Command
from mtp.workflow.consumers import (
    CoreBundle,
    ProfileKind,
    open_core_bundle,
    preflight_consumer_outputs,
)
from mtp.workflow.operations import run_strain_rollup
from mtp.workflow.runtime import (
    artifact_path,
    content_hash,
    inventory,
    repository_state,
    sha256_file,
)
from mtp.workflow.verification import verify_core_run

Planner = Callable[[ArmContext], list[Command]]
# Steps that rebuild matrix counts from the KG, and so must apply the core's strain rollup
# scope, with the argument parser each module runs (#788, #881)
SCOPED_MODULES: dict[str, Callable[[type[argparse.ArgumentParser]], argparse.ArgumentParser]] = {
    "mtp.taxon_contributions": taxon_contributions.argument_parser
}
# These entrypoints require an explicit basis rather than choosing a default.
DIRECTION_REQUIRED_SCRIPTS = frozenset(
    {
        "scripts/create_taxa_guild_heatmap.py",
        "scripts/create_study_trait_profile_heatmap.py",
    }
)
INPUT_NAMES = {
    "raw_dump": "full_dump.csv",
    "kg_edges": "merged-kg_edges.tsv",
    "kg_nodes": "merged-kg_nodes.tsv",
    "source_identity": "source_identity.json",
    "study_table": "studies.tsv",
    "experiments": "experiments.tsv",
    "experiment_ontology_terms": "experiment_ontology_terms.tsv",
    "signatures": "signatures.tsv",
    "reported_taxa": "bugsigdb_taxa.tsv",
    "taxon_profiles_wide": "taxon_profiles_wide.tsv",
    "taxon_profiles_long": "taxon_profiles_long.tsv",
    "taxa_matrix": "taxa_trait_matrix.tsv",
    "feature_metadata": "taxa_trait_matrix_columns.tsv",
    "feature_catalog": "feature_catalog.tsv",
    "count_profiles": "study_trait_profiles_wide.tsv",
    "count_coverage": "count_coverage.tsv",
    "ora_profiles": "study_trait_profiles_ora_wide.tsv",
    "ora_coverage": "ora_coverage.tsv",
}


class ArmExecutionError(ValueError):
    """A failed attempt whose diagnostic manifest remains available."""

    def __init__(self, message: str, record: ArmRun) -> None:
        super().__init__(f"{message}; execution record: {record.manifest_path}")
        self.record = record


def _snapshot(path: Path) -> RecordedPath:
    if not path.is_absolute() or path.resolve() != path:
        raise ValueError(f"Snapshot path must be absolute and free of symlink aliases: {path}")
    kind = StorageKind.collection if path.is_dir() else StorageKind.file
    digest, size, members = inventory(path, kind)
    return RecordedPath(
        path=str(path.absolute()),
        storage_kind=kind,
        sha256=digest,
        size_bytes=size,
        inventory=members,
    )


def _check_snapshot(record: RecordedPath, *, exact: bool = False) -> None:
    path = Path(record.path)
    if not path.is_absolute() or path.resolve() != path:
        raise ValueError(f"Recorded path is not canonical: {path}")
    if record.storage_kind == "file" or exact:
        if _snapshot(path) != record:
            raise ValueError(f"Recorded bytes changed: {path}")
        return
    if path.is_symlink() or not path.is_dir():
        raise ValueError(f"Recorded collection unavailable: {path}")
    if content_hash([item.model_dump(mode="json") for item in record.inventory]) != record.sha256:
        raise ValueError(f"Invalid collection inventory digest: {path}")
    if sum(item.size_bytes for item in record.inventory) != record.size_bytes:
        raise ValueError(f"Invalid collection inventory size: {path}")
    seen: set[str] = set()
    for member in record.inventory:
        relative = Path(member.path)
        if relative.is_absolute() or ".." in relative.parts or member.path in seen:
            raise ValueError(f"Invalid collection inventory member: {member.path}")
        seen.add(member.path)
        target = path / relative
        if target.resolve() != target.absolute() or sha256_file(target) != member.sha256:
            raise ValueError(f"Recorded collection member changed: {target}")
        if target.stat().st_size != member.size_bytes:
            raise ValueError(f"Recorded collection member size changed: {target}")


def _code_inventory(repository: Path) -> list[FileInventoryEntry]:
    paths = {repository / name for name in ("Justfile", "pyproject.toml", "uv.lock")}
    for directory in ("src", "scripts", "methods", "schema", "workflows"):
        paths.update(
            path
            for path in (repository / directory).rglob("*")
            if path.is_file()
            and "__pycache__" not in path.parts
            and path.suffix in {".py", ".sh", ".json", ".yaml", ".yml", ".tsv", ".md"}
        )
    return [
        FileInventoryEntry(
            path=path.relative_to(repository).as_posix(),
            sha256=sha256_file(path),
            size_bytes=path.stat().st_size,
        )
        for path in sorted(paths)
    ]


def _save(record: ArmRun) -> None:
    checked = ArmRun.model_validate(record.model_dump(mode="json"))
    path = Path(record.manifest_path)
    if path.parent.resolve() != path.parent:
        raise ValueError(f"Execution record directory changed: {path.parent}")
    temporary = path.with_name(f".arm-run-{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("x", encoding="utf-8") as stream:
            stream.write(checked.model_dump_json(indent=2) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _tree_inventory(root: Path) -> list[FileInventoryEntry]:
    entries = []
    for path in sorted(root.rglob("*")):
        if path == root / "arm_run.json":
            continue
        if path.is_symlink() or not (path.is_dir() or path.is_file()):
            raise ValueError(f"Non-regular rerun output: {path}")
        if path.is_file():
            entries.append(
                FileInventoryEntry(
                    path=path.relative_to(root).as_posix(),
                    sha256=sha256_file(path),
                    size_bytes=path.stat().st_size,
                )
            )
    return entries


def _definition(command: Command) -> ArmCommandDefinition:
    return ArmCommandDefinition(
        id=command.id,
        argv=list(command.argv),
        inputs=[str(path.absolute()) for path in command.inputs],
        outputs=[str(path.absolute()) for path in command.outputs],
    )


class _StepArguments(argparse.ArgumentParser):
    """Parse a planned step's arguments as its module does, raising instead of exiting."""

    def error(self, message: str) -> NoReturn:
        raise ValueError(message)

    def exit(self, status: int = 0, message: str | None = None) -> NoReturn:
        raise ValueError(message or f"argument parsing stopped with status {status}")


class _MissingDirectionBasis(ValueError):
    """A required figure basis is absent; historical replay needs its recorded code."""


def _step_strain_rollup(command: Command, *, absent_scope_is_off: bool) -> str:
    """The strain rollup scope a scoped step runs under, read by its module's own parser.

    The parser's rules decide: the equals form and an unambiguous abbreviation count, the
    last of repeated options wins, and a malformed option is refused. A plan must state the
    scope; only a plan recorded before #788 states none, and it ran under off (#881).
    """
    parser = SCOPED_MODULES[command.argv[2]](_StepArguments)
    parser.set_defaults(strain_rollup=None)
    try:
        scope = parser.parse_args(command.argv[3:]).strain_rollup
    except ValueError as error:
        raise ValueError(
            f"{command.id} has invalid {command.argv[2]} arguments: {error}"
        ) from error
    if scope is not None:
        return str(scope)
    if not absent_scope_is_off:
        raise ValueError(
            f"{command.id} does not state --strain-rollup; a new plan must name its core's scope"
        )
    return "off"


def _step_direction_basis(command: Command, repository: Path) -> str | None:
    """Use registered figures' complete execution parsers, including rejection rules."""
    script = Path(command.argv[1])
    relative = (
        script.relative_to(repository).as_posix() if script.is_relative_to(repository) else ""
    )
    registered = next(
        (name for name in DIRECTION_REQUIRED_SCRIPTS if name.casefold() == relative.casefold()),
        None,
    )
    if registered is not None and relative != registered:
        raise ValueError(f"{command.id} must use canonical script path {registered!r}")
    if registered is not None:
        # Load only this verifier's canonical registered implementation, never code
        # selected by an arbitrary recorded repository path. main() uses this builder.
        namespace = runpy.run_path(str(REPO_ROOT / registered))
        parser = namespace["build_parser"](_StepArguments)
        if not isinstance(parser, argparse.ArgumentParser):
            raise TypeError(f"{registered} did not supply an argument parser")
    else:
        # Unrelated entrypoints retain their own options (e.g. --direction inc).
        parser = _StepArguments(add_help=False, allow_abbrev=False)
        parser.add_argument("--direction-basis", choices=[basis.value for basis in DirectionBasis])
    arguments = command.argv[3:] if command.argv[1] == "-m" else command.argv[2:]
    try:
        parsed = (
            parser.parse_args(arguments) if registered else parser.parse_known_args(arguments)[0]
        )
    except ValueError as error:
        if str(error).startswith(
            "the following arguments are required:"
        ) and "--direction-basis" in str(error):
            raise _MissingDirectionBasis(
                f"{command.id} has invalid direction-basis arguments: {error}"
            ) from error
        raise ValueError(f"{command.id} has invalid direction-basis arguments: {error}") from error
    return str(parsed.direction_basis) if parsed.direction_basis is not None else None


def _validate_plan(
    commands: Sequence[Command], context: ArmContext, *, absent_scope_is_off: bool = False
) -> None:
    """Check a plan before execution or on verification.

    ``absent_scope_is_off`` is only for verifying a plan recorded on a core that recorded no
    strain rollup scope, as every core before #788 did (#881).
    """
    seen: set[str] = set()
    available = {context.manifest, *context.products.values()}
    available.update(context.inputs / name for name in INPUT_NAMES.values())
    available.add(context.inputs)
    previous_outputs: set[Path] = set()
    for command in commands:
        _definition(command)
        if (
            command.id in seen
            or len(command.argv) < 2
            or not Path(command.argv[0]).is_absolute()
            or any("\0" in argument for argument in command.argv)
        ):
            raise ValueError(f"Invalid or duplicate Python command: {command.id}")
        if command.argv[1] == "-m":
            if len(command.argv) < 3 or not re.fullmatch(
                r"mtp(?:\.[A-Za-z_][A-Za-z_0-9]*)+", command.argv[2]
            ):
                raise ValueError(f"Command module is outside the declared package: {command.id}")
        else:
            script = Path(command.argv[1])
            if (
                script.suffix != ".py"
                or not script.is_relative_to(context.repository / "scripts")
                or script.resolve() != script
            ):
                raise ValueError(f"Command must execute a tracked Python entrypoint: {command.id}")
        if command.argv[1] == "-m" and command.argv[2] in SCOPED_MODULES:
            scope = _step_strain_rollup(command, absent_scope_is_off=absent_scope_is_off)
            if scope != context.strain_rollup:
                raise ValueError(
                    f"{command.id} uses strain rollup scope {scope!r}, but its core "
                    f"recorded {context.strain_rollup!r}"
                )
        direction_basis = _step_direction_basis(command, context.repository)
        if direction_basis is not None and direction_basis != context.direction_basis:
            raise ValueError(
                f"{command.id} uses direction basis {direction_basis!r}, but its core "
                f"recorded {context.direction_basis!r}"
            )
        seen.add(command.id)
        if not command.outputs:
            raise ValueError(f"No declared outputs: {command.id}")
        for path in (*command.inputs, *command.outputs):
            if not path.is_absolute() or path.resolve() != path:
                raise ValueError(f"Command path must be absolute without aliases: {path}")
        for source in command.inputs:
            static_source = any(
                source.is_relative_to(context.repository / directory)
                for directory in ("src", "scripts", "methods", "schema", "workflows")
            )
            prior_source = any(source.is_relative_to(path) for path in previous_outputs)
            if source not in available and not static_source and not prior_source:
                raise ValueError(f"Input is not core, source, or an earlier output: {source}")
        for output in command.outputs:
            resolved = output.resolve()
            if not resolved.is_relative_to(context.output) or resolved == context.output:
                raise ValueError(f"Output outside dedicated rerun directory: {output}")
            if resolved.is_relative_to(context.inputs) or resolved.is_relative_to(
                context.output / "logs"
            ):
                raise ValueError(f"Output overlaps reserved input/log directory: {output}")
            if resolved in {context.output / "arm_run.json", context.output / "arm_run.tmp"}:
                raise ValueError(f"Output overlaps execution record: {output}")
        previous_outputs.update(command.outputs)
    if not commands:
        raise ValueError("Empty arm command plan")


def _plan_hash(plan: Sequence[ArmCommandDefinition]) -> str:
    return content_hash([item.model_dump(mode="json") for item in plan])


def _can_extend(path: Path, outputs: Sequence[Path]) -> bool:
    """An input collection may gain files only under explicitly declared outputs."""
    return any(output.is_relative_to(path) for output in outputs)


def _check_input_after_command(record: RecordedPath, outputs: Sequence[Path]) -> None:
    path = Path(record.path)
    if record.storage_kind != "collection" or not _can_extend(path, outputs):
        _check_snapshot(record, exact=True)
        return
    _check_snapshot(record)
    current = _snapshot(path)
    old_members = {member.path for member in record.inventory}
    for member in current.inventory:
        target = path / member.path
        if member.path not in old_members and not any(
            target.is_relative_to(output) for output in outputs
        ):
            raise ValueError(f"Input collection gained an undeclared member: {target}")


async def _run_process(command: Command, repository: Path, log: Path) -> int:
    environment = dict(os.environ)
    environment.update({"PYTHONPATH": str(repository / "src"), "MPLBACKEND": "Agg"})
    with log.open("wb") as stream:
        process = await asyncio.create_subprocess_exec(
            *command.argv,
            cwd=repository,
            env=environment,
            stdout=stream,
            stderr=asyncio.subprocess.STDOUT,
        )
        try:
            return await process.wait()
        finally:
            if process.returncode is None:
                process.terminate()
                await process.wait()


def _prepare_inputs(bundle: CoreBundle, directory: Path) -> list[RecordedPath]:
    directory.mkdir()
    prepared = []
    for product, filename in INPUT_NAMES.items():
        source = bundle.path(product)
        destination = directory / filename
        shutil.copyfile(source, destination)
        observed = _snapshot(destination)
        artifact = next(item for item in bundle.run.artifacts or [] if item.product_id == product)
        if observed.sha256 != artifact.sha256 or observed.size_bytes != artifact.size_bytes:
            raise ValueError(f"Core product changed while preparing inputs: {product}")
        prepared.append(observed)
    return prepared


def execute_arm(
    arm: str,
    manifest: Path,
    output: Path,
    planner: Planner,
    *,
    replicates: int | None = None,
    seed: int | None = None,
    figures: bool = True,
    repository: Path = REPO_ROOT,
) -> ArmRun:
    """Verify sources, execute the complete selected plan, and retain failed attempts."""
    if output.is_symlink():
        raise ValueError(f"Output must not be a symlink: {output}")
    repository, manifest, output = repository.resolve(), manifest.resolve(), output.resolve()
    if arm not in {"ilyas", "noel", "jed"} or (replicates is not None and replicates < 1):
        raise ValueError("Invalid arm or replicate count")
    if output.exists():
        raise ValueError(f"Choose a fresh output directory; already exists: {output}")
    for directory in (
        "src",
        "scripts",
        "methods",
        "schema",
        "workflows",
        "data/raw",
        "data/studies",
        "data/kg-microbe",
    ):
        if output.is_relative_to(repository / directory):
            raise ValueError(f"Rerun output overlaps source directory: {output}")
    profile: ProfileKind = "ora_profiles" if arm == "ilyas" else "count_profiles"
    bundle = open_core_bundle(manifest, profile)
    preflight_consumer_outputs(bundle, output / "arm_run.json", output_directory=output)
    commit, dirty = repository_state(repository)
    context = ArmContext(
        repository=repository,
        manifest=manifest,
        output=output,
        inputs=output / "inputs",
        products={
            item.product_id: bundle.path(item.product_id) for item in bundle.run.artifacts or []
        },
        code_commit=commit + ("+dirty" if dirty else ""),
        replicates=replicates,
        seed=seed,
        figures=figures,
        contrast_population=ContrastPopulation(bundle.run.contrast_population).value,
        direction_basis=DirectionBasis(bundle.run.direction_basis).value,
        strain_rollup=run_strain_rollup(bundle.run),
    )
    commands = planner(context)
    _validate_plan(commands, context)
    if any(command.argv[0] != sys.executable for command in commands):
        raise ValueError("Commands must use this Python executable")
    plan = [_definition(command) for command in commands]
    record = ArmRun(
        id=f"arm-{uuid.uuid4().hex}",
        arm=arm,
        scope="standard_analysis_calibration_and_optional_figures",
        status=ExecutionStatus.running,
        started_at=datetime.now(UTC),
        manifest_path=str(output / "arm_run.json"),
        core_manifest_path=str(manifest),
        core_manifest_sha256=sha256_file(manifest),
        core_run_id=bundle.run.run_id,
        profile_product=profile,
        contrast_population=bundle.run.contrast_population,
        direction_basis=bundle.run.direction_basis,
        git_commit=commit,
        git_dirty=dirty,
        repository_path=str(repository),
        code_inventory=_code_inventory(repository),
        runtime_json=json.dumps(
            {
                "python": sys.version,
                "executable": sys.executable,
                "platform": platform.platform(),
                "packages": sorted(
                    (str(item.metadata["Name"]), item.version)
                    for item in importlib.metadata.distributions()
                ),
                "process_settings": {"PYTHONPATH": str(repository / "src"), "MPLBACKEND": "Agg"},
            },
            sort_keys=True,
        ),
        parameters=[
            ParameterValue(name=name, value=str(value), source=source)
            for name, value, source in (
                (
                    "replicates",
                    replicates if replicates is not None else "stage_defaults",
                    "override" if replicates is not None else "default",
                ),
                (
                    "seed",
                    seed if seed is not None else "stage_defaults",
                    "override" if seed is not None else "default",
                ),
                ("figures", figures, "override" if not figures else "default"),
            )
        ],
        plan=plan,
        plan_sha256=_plan_hash(plan),
        steps=[],
        prepared_inputs=[],
        outputs=[],
    )
    output.mkdir(parents=True, exist_ok=False)
    (output / "logs").mkdir()
    _save(record)
    try:
        record.prepared_inputs = _prepare_inputs(bundle, context.inputs)
        _save(record)
        for command in commands:
            step = ArmCommandExecution(
                id=command.id,
                status=ExecutionStatus.running,
                started_at=datetime.now(UTC),
                inputs=[],
                outputs=[],
            )
            record.steps = [*(record.steps or []), step]
            log = output / "logs" / f"{command.id}.log"
            _save(record)
            try:
                if _code_inventory(repository) != record.code_inventory:
                    raise ValueError("Code or configuration changed during rerun")
                step.inputs = [_snapshot(path) for path in dict.fromkeys(command.inputs)]
                for path in command.outputs:
                    path.parent.mkdir(parents=True, exist_ok=True)
                print(f"[{arm}] {command.id}", flush=True)
                step.exit_code = asyncio.run(_run_process(command, repository, log))
                step.log = _snapshot(log)
                if step.exit_code != 0:
                    raise ValueError(f"Command {command.id} exited {step.exit_code}; see {log}")
                preflight_consumer_outputs(bundle, output / "arm_run.json", output_directory=output)
                for item in step.inputs:
                    _check_input_after_command(item, command.outputs)
                step.outputs = [_snapshot(path) for path in dict.fromkeys(command.outputs)]
                if _code_inventory(repository) != record.code_inventory:
                    raise ValueError("Code or configuration changed during rerun")
                step.status = ExecutionStatus.completed
            except BaseException as error:
                step.status, step.error = ExecutionStatus.failed, str(error)
                if log.is_file():
                    step.log = _snapshot(log)
                raise
            finally:
                step.ended_at = datetime.now(UTC)
                _save(record)
        if sha256_file(manifest) != record.core_manifest_sha256:
            raise ValueError("Upstream core manifest changed during rerun")
        verify_core_run(manifest)
        for item in record.prepared_inputs:
            _check_snapshot(item, exact=True)
        for completed_step in record.steps or []:
            for item in completed_step.outputs or []:
                _check_snapshot(item)
        if _code_inventory(repository) != record.code_inventory:
            raise ValueError("Code or configuration changed during rerun")
        record.outputs = _tree_inventory(output)
        record.status = ExecutionStatus.completed
    except BaseException as error:
        record.status, record.error = ExecutionStatus.failed, str(error)
        raise ArmExecutionError(str(error), record) from error
    finally:
        record.ended_at = datetime.now(UTC)
        _save(record)
    return record


def verify_arm(path: Path) -> ArmRun:
    """Reconcile completed stage evidence, immutable core inputs and final output bytes."""
    path = path.resolve()
    record = ArmRun.model_validate_json(path.read_text(encoding="utf-8"))
    if record.status != ExecutionStatus.completed or not record.ended_at or record.error:
        raise ValueError("Arm rerun is not a completed successful execution")
    if Path(record.manifest_path) != path:
        raise ValueError("Arm manifest path does not match recorded location")
    manifest = Path(record.core_manifest_path)
    if sha256_file(manifest) != record.core_manifest_sha256:
        raise ValueError("Core manifest checksum changed")
    core = verify_core_run(manifest)
    if core.run_id != record.core_run_id:
        raise ValueError("Core execution identity changed")
    if (
        core.contrast_population != record.contrast_population
        or core.direction_basis != record.direction_basis
    ):
        raise ValueError("Arm contrast population or direction basis differs from its core")
    expected_profile = "ora_profiles" if record.arm == "ilyas" else "count_profiles"
    if record.profile_product != expected_profile:
        raise ValueError("Arm profile unit is incompatible")
    if _plan_hash(record.plan) != record.plan_sha256:
        raise ValueError("Arm command plan checksum changed")
    products = {
        item.product_id: artifact_path(core, item.product_id) for item in core.artifacts or []
    }
    context = ArmContext(
        repository=Path(record.repository_path),
        manifest=manifest,
        output=path.parent,
        inputs=path.parent / "inputs",
        products=products,
        code_commit=record.git_commit,
        replicates=None,
        seed=None,
        figures=True,
        contrast_population=ContrastPopulation(record.contrast_population).value,
        direction_basis=DirectionBasis(record.direction_basis).value,
        strain_rollup=run_strain_rollup(core),
    )
    commands = [
        Command(
            item.id,
            tuple(item.argv),
            tuple(map(Path, item.inputs or [])),
            tuple(map(Path, item.outputs)),
        )
        for item in record.plan
    ]
    # Only a core from before #788 recorded no scope, and only its plans may omit it (#881)
    # Known heatmap entrypoints always require an explicit direction basis (#843/#944).
    # Pre-#949 plans use their recorded implementation for historical verification;
    # the current verifier does not infer missing direction arguments (#955).
    declared = {item.name for item in core.parameters or []}
    try:
        _validate_plan(commands, context, absent_scope_is_off="strain_rollup" not in declared)
    except _MissingDirectionBasis as error:
        raise ValueError(
            f"{error}. For a historical record, verify with its recorded commit "
            f"{record.git_commit} and locked environment; retain the original checkout "
            "and input/output paths. See docs/workflow.md for historical arm and "
            "comparison verification. The current verifier does not infer a missing basis."
        ) from error
    runtime = json.loads(record.runtime_json)
    if not isinstance(runtime, dict) or any(
        command.argv[0] != runtime.get("executable") for command in commands
    ):
        raise ValueError("Command executable disagrees with recorded runtime")
    if [step.id for step in record.steps or []] != [command.id for command in record.plan]:
        raise ValueError("Arm execution does not cover its complete ordered plan")
    if len({command.id for command in record.plan}) != len(record.plan) or not record.plan:
        raise ValueError("Arm plan has missing or duplicate commands")
    for definition, step in zip(record.plan, record.steps or [], strict=True):
        if (
            step.status != ExecutionStatus.completed
            or step.exit_code != 0
            or not step.ended_at
            or step.error
            or step.log is None
        ):
            raise ValueError(f"Incomplete or failed arm command: {step.id}")
        if [item.path for item in step.inputs or []] != list(
            dict.fromkeys(definition.inputs or [])
        ):
            raise ValueError(f"Command input evidence disagrees: {step.id}")
        if [item.path for item in step.outputs or []] != list(dict.fromkeys(definition.outputs)):
            raise ValueError(f"Command output evidence disagrees: {step.id}")
        expected_log = path.parent / "logs" / f"{step.id}.log"
        if Path(step.log.path) != expected_log:
            raise ValueError(f"Command log location disagrees: {step.id}")
        for item in [*(step.inputs or []), *(step.outputs or []), step.log]:
            _check_snapshot(item)
    expected_names = {str(path.parent / "inputs" / name) for name in INPUT_NAMES.values()}
    if {item.path for item in record.prepared_inputs or []} != expected_names:
        raise ValueError("Prepared core inputs do not match the declared product set")
    artifacts = {item.product_id: item for item in core.artifacts or []}
    prepared = {item.path: item for item in record.prepared_inputs or []}
    for product, filename in INPUT_NAMES.items():
        item = prepared[str(path.parent / "inputs" / filename)]
        if (
            item.sha256 != artifacts[product].sha256
            or item.size_bytes != artifacts[product].size_bytes
        ):
            raise ValueError(f"Prepared input differs from upstream product: {product}")
        _check_snapshot(item, exact=True)
    if _tree_inventory(path.parent) != record.outputs:
        raise ValueError("Arm output inventory changed")
    return record
