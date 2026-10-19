"""Execute registered Python stages with verified, immutable bundles (#629).

This is an explicit-I/O adapter, not a filesystem sandbox or shell scheduler.
Operations must obtain scientific inputs and settings from their arguments.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import inspect
import json
import os
import platform
import sys
import uuid
from collections.abc import Callable, Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path

from dulwich.errors import NotGitRepository
from dulwich.porcelain import status
from dulwich.repo import Repo

from mtp.workflow.definition import read_definition, safe_relative_path, unique, validate_definition
from mtp.workflow.models import (
    ArtifactInstance,
    ExecutionStatus,
    FileInventoryEntry,
    ParameterValue,
    ProductDefinition,
    StepExecution,
    ValidationResult,
    ValidationStatus,
    WorkflowDefinition,
    WorkflowRun,
)

Operation = Callable[[dict[str, Path], dict[str, Path], dict[str, object]], None]
Validator = Callable[[ProductDefinition, Path, Mapping[str, Path]], Sequence[ValidationResult]]


class WorkflowExecutionError(ValueError):
    """Execution failure with the persisted diagnostic run attached."""

    def __init__(self, message: str, run: WorkflowRun) -> None:
        super().__init__(message)
        self.run = run


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def content_hash(value: object) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    """Hash a regular file, checking that it did not change while being read."""
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Expected regular file, without symlinks: {path}")
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    if (before.st_ino, before.st_size, before.st_mtime_ns) != (
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise ValueError(f"File changed while hashing: {path}")
    return digest.hexdigest()


def inventory(path: Path, storage_kind: str) -> tuple[str, int, list[FileInventoryEntry]]:
    """Inventory all collection files, including hidden members; reject symlinks."""
    if path.is_symlink():
        raise ValueError(f"Symlink artifact is not permitted: {path}")
    if storage_kind == "file":
        digest = sha256_file(path)
        member = FileInventoryEntry(path=path.name, sha256=digest, size_bytes=path.stat().st_size)
        return digest, member.size_bytes, [member]
    if storage_kind != "collection" or not path.is_dir():
        raise ValueError(f"Expected collection directory: {path}")
    members: list[FileInventoryEntry] = []
    for item in sorted(path.rglob("*")):
        if item.is_symlink() or not (item.is_file() or item.is_dir()):
            raise ValueError(f"Non-regular collection member: {item}")
        if item.is_file():
            members.append(
                FileInventoryEntry(
                    path=item.relative_to(path).as_posix(),
                    sha256=sha256_file(item),
                    size_bytes=item.stat().st_size,
                )
            )
    if not members:
        raise ValueError(f"Empty collection: {path}")
    digest = content_hash([member.model_dump(mode="json") for member in members])
    return digest, sum(member.size_bytes for member in members), members


def _atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    with temporary.open("x", encoding="utf-8") as stream:
        stream.write(canonical_json(value) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _save_run(run: WorkflowRun) -> None:
    validated = WorkflowRun.model_validate(run.model_dump(mode="json"))
    _atomic_json(Path(run.manifest_path), validated.model_dump(mode="json"))


def _resolve_parameters(
    definition: WorkflowDefinition, overrides: Mapping[str, object]
) -> list[ParameterValue]:
    definitions = {item.name: item for item in definition.parameters or []}
    cohort_parameters = {item.cohort_parameter for item in definition.products}
    if set(overrides) - set(definitions):
        raise ValueError(
            f"Undeclared parameter overrides: {sorted(set(overrides) - set(definitions))}"
        )
    values: list[ParameterValue] = []
    for name, item in sorted(definitions.items()):
        value: object = item.default
        source = "default"
        if item.environment and item.environment in os.environ:
            value, source = os.environ[item.environment], f"environment:{item.environment}"
        if name in overrides:
            value, source = overrides[name], "override"
        if value is None:
            if item.required:
                raise ValueError(f"Required parameter missing: {name}")
            continue
        if not isinstance(value, (str, int, float, bool)):
            raise ValueError(f"Parameter {name} must be a scalar")
        resolved = str(value).lower() if isinstance(value, bool) else str(value)
        if name in cohort_parameters:
            resolved = resolved.strip()
            if not resolved:
                raise ValueError(f"Cohort parameter {name} cannot be empty")
        values.append(ParameterValue(name=name, value=resolved, source=source))
    return values


def _effective_product(
    product: ProductDefinition, parameters: Mapping[str, str]
) -> ProductDefinition:
    value = product.model_dump(mode="json")
    if product.cohort_parameter:
        if product.cohort_parameter not in parameters:
            raise ValueError(f"Missing cohort parameter for {product.id}")
        value["cohort"] = parameters[product.cohort_parameter]
    return ProductDefinition.model_validate(value)


def _code_identity(
    root: Path, registry: Mapping[str, Operation]
) -> tuple[str, list[FileInventoryEntry]]:
    paths = set((root / "src").rglob("*.py")) | set((root / "schema").rglob("*.yaml"))
    for operation in registry.values():
        filename = inspect.getsourcefile(operation)
        if not filename:
            raise ValueError("Operations must have inspectable Python source files")
        paths.add(Path(filename).resolve())
    entries: list[FileInventoryEntry] = []
    for path in sorted(paths):
        path = path.resolve()
        name = path.relative_to(root).as_posix() if path.is_relative_to(root) else str(path)
        entries.append(
            FileInventoryEntry(path=name, sha256=sha256_file(path), size_bytes=path.stat().st_size)
        )
    if not entries:
        raise ValueError("No code files available for fingerprint")
    return content_hash([item.model_dump(mode="json") for item in entries]), entries


def repository_state(root: Path) -> tuple[str, bool]:
    """Read commit and dirty state, including linked worktrees, without spawning Git."""
    try:
        with Repo(str(root)) as repository:
            commit = repository.head().decode("ascii")
            changes = status(
                repository, ignored=False, untracked_files="normal", optional_locks=False
            )
            dirty = bool(changes.unstaged or changes.untracked or any(changes.staged.values()))
    except (NotGitRepository, KeyError, OSError, UnicodeError) as error:
        raise ValueError(f"Cannot read repository provenance from {root}: {error}") from error
    return commit, dirty


def _fingerprint(run: WorkflowRun) -> str:
    external = {product.id for product in run.definition.products if product.external}
    runtime = json.loads(run.runtime_json)
    runtime.pop("executable", None)
    return content_hash(
        {
            "definition": run.definition_sha256,
            "spec": run.spec_sha256,
            "code": run.code_sha256,
            "lock": run.lock_sha256,
            "runtime": runtime,
            "parameters": sorted((item.name, item.value) for item in run.parameters or []),
            "inputs": sorted(
                (item.product_id, item.sha256)
                for item in run.artifacts or []
                if item.product_id in external
            ),
        }
    )


def _make_artifact(
    run: WorkflowRun,
    product: ProductDefinition,
    path: Path,
    producer: str | None,
    validator: Validator | None,
    context: Mapping[str, Path],
) -> ArtifactInstance:
    digest, size, members = inventory(path, str(product.storage_kind))
    validations = [
        ValidationResult(
            validator="artifact_inventory_v1",
            status=ValidationStatus.passed,
            message="All regular files inventoried and hashed; no symlinks or special files.",
        )
    ]
    if validator:
        extra = list(validator(product, path, context))
        if not extra:
            raise ValueError(f"Validator returned no results for {product.id}")
        validations.extend(ValidationResult.model_validate(item) for item in extra)
    elif product.kind != "structural":
        raise ValueError(f"Scientific validator required for product {product.id} ({product.kind})")
    if any(result.status != "passed" for result in validations):
        raise ValueError(f"Product validation failed for {product.id}: {validations}")
    if inventory(path, str(product.storage_kind))[0] != digest:
        raise ValueError(f"Validator changed artifact {product.id}")
    return ArtifactInstance(
        id=f"{run.run_id}:artifact:{product.id}",
        product_id=product.id,
        run_id=run.run_id,
        producer_step=producer,
        path=str(path),
        sha256=digest,
        size_bytes=size,
        inventory=members,
        validations=validations,
        unit=product.unit,
        cohort=product.cohort,
        feature_policy=product.feature_policy,
        coverage_policy=product.coverage_policy,
    )


def artifact_path(run: WorkflowRun, product_id: str) -> Path:
    """Resolve exactly one recorded product binding, never guess legacy paths."""
    matches = [item for item in run.artifacts or [] if item.product_id == product_id]
    if len(matches) != 1:
        raise ValueError(f"Expected one artifact for {product_id}; found {len(matches)}")
    path = Path(matches[0].path)
    return path if path.is_absolute() else Path(run.manifest_path).parent / path


def _verify_artifact(run: WorkflowRun, artifact: ArtifactInstance) -> None:
    products = {item.id: item for item in run.definition.products}
    product = products[artifact.product_id]
    path = artifact_path(run, product.id)
    if (
        product.external
        and artifact.requested_path
        and Path(artifact.requested_path).resolve(strict=True) != path
    ):
        raise ValueError(f"External binding target changed: {product.id}")
    if not product.external:
        run_dir = Path(run.manifest_path).parent
        if not path.resolve().is_relative_to(run_dir.resolve()):
            raise ValueError(f"Output escapes run directory: {product.id}")
        if any(parent.is_symlink() for parent in path.parents if parent.is_relative_to(run_dir)):
            raise ValueError(f"Symlink output directory: {product.id}")
    actual = inventory(path, str(product.storage_kind))
    expected = (
        artifact.sha256,
        artifact.size_bytes,
        [item.model_dump(mode="json") for item in artifact.inventory],
    )
    observed = (actual[0], actual[1], [item.model_dump(mode="json") for item in actual[2]])
    if observed != expected:
        raise ValueError(f"Artifact content or inventory changed: {product.id}")


def execute_workflow(
    spec_path: str | Path,
    input_paths: Mapping[str, Path],
    output_root: str | Path,
    operation_registry: Mapping[str, Operation],
    *,
    parameter_overrides: Mapping[str, object] | None = None,
    validator: Validator | None = None,
    repository_root: Path | None = None,
) -> WorkflowRun:
    """Run explicit-I/O Python operations; failures persist diagnostics and raise."""
    spec_path = Path(spec_path).resolve()
    definition = read_definition(spec_path)
    order = validate_definition(definition)
    products = {item.id: item for item in definition.products}
    steps = {item.id: item for item in definition.steps}
    required = {item.id for item in products.values() if item.external}
    if set(input_paths) != required:
        raise ValueError(f"External bindings must be exactly {sorted(required)}")
    operations = {item.operation for item in definition.steps}
    if operations - set(operation_registry):
        raise ValueError(f"Unregistered operations: {sorted(operations - set(operation_registry))}")
    registry = {name: operation_registry[name] for name in sorted(operations)}
    parameters = _resolve_parameters(definition, parameter_overrides or {})
    resolved = {item.name: item.value for item in parameters}
    root = (repository_root or Path(__file__).resolve().parents[3]).resolve()
    code_hash, code_entries = _code_identity(root, registry)
    lock_hash = sha256_file(root / "uv.lock")
    git_commit, git_dirty = repository_state(root)
    run_id = f"run-{uuid.uuid4().hex}"
    run_dir = Path(output_root).resolve() / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    (run_dir / ".staging").mkdir()
    (run_dir / "steps").mkdir()
    run = WorkflowRun(
        run_id=run_id,
        definition=definition,
        status=ExecutionStatus.running,
        started_at=datetime.now(UTC),
        manifest_path=str(run_dir / "run.json"),
        spec_path=str(spec_path),
        spec_sha256=sha256_file(spec_path),
        definition_sha256=content_hash(definition.model_dump(mode="json")),
        execution_fingerprint="0" * 64,
        code_sha256=code_hash,
        code_inventory=code_entries,
        lock_sha256=lock_hash,
        git_commit=git_commit,
        git_dirty=git_dirty,
        runtime_json=canonical_json(
            {
                "python": sys.version,
                "executable": sys.executable,
                "platform": platform.platform(),
                "packages": sorted(
                    (str(item.metadata["Name"]), item.version)
                    for item in importlib.metadata.distributions()
                ),
                "operations": {
                    name: f"{value.__module__}.{value.__qualname__}"
                    for name, value in registry.items()
                },
            }
        ),
        parameters=parameters,
        steps=[],
        artifacts=[],
        outputs=[],
    )
    _save_run(run)
    context = {key: Path(value).resolve() for key, value in input_paths.items()}
    try:
        for product_id in sorted(required):
            artifact = _make_artifact(
                run,
                _effective_product(products[product_id], resolved),
                context[product_id],
                None,
                validator,
                context,
            )
            artifact.requested_path = str(Path(input_paths[product_id]).absolute())
            run.artifacts = [*(run.artifacts or []), artifact]
        run.execution_fingerprint = _fingerprint(run)
        _save_run(run)
        for step_id in order:
            step = steps[step_id]
            step_run = StepExecution(
                id=f"{run_id}:step:{step_id}",
                step_id=step_id,
                status=ExecutionStatus.running,
                started_at=datetime.now(UTC),
                input_artifacts=[
                    item.id for item in run.artifacts or [] if item.product_id in step.inputs
                ],
                output_artifacts=[],
                parameters=[item for item in parameters if item.name in (step.parameters or [])],
            )
            run.steps = [*(run.steps or []), step_run]
            _save_run(run)
            staging = run_dir / ".staging" / step_id
            staging.mkdir()
            inputs = {key: context[key] for key in step.inputs}
            outputs = {key: staging / products[key].filename for key in step.outputs}
            for key, path in outputs.items():
                path.parent.mkdir(parents=True, exist_ok=True)
                if products[key].storage_kind == "collection":
                    path.mkdir()
            for artifact in run.artifacts or []:
                _verify_artifact(run, artifact)
            registry[step.operation](
                inputs, outputs, {item.name: item.value for item in step_run.parameters or []}
            )
            for artifact in run.artifacts or []:
                _verify_artifact(run, artifact)
            current_context = {**inputs, **outputs}
            produced = [
                _make_artifact(
                    run,
                    _effective_product(products[key], resolved),
                    outputs[key],
                    step_run.id,
                    validator,
                    current_context,
                )
                for key in step.outputs
            ]
            declared_members = {
                path.relative_to(staging).as_posix()
                for key, path in outputs.items()
                if products[key].storage_kind == "file"
            }
            for key, path in outputs.items():
                if products[key].storage_kind == "collection":
                    declared_members.update(
                        item.relative_to(staging).as_posix()
                        for item in path.rglob("*")
                        if item.is_file()
                    )
            actual_members = {
                path.relative_to(staging).as_posix()
                for path in staging.rglob("*")
                if path.is_file()
            }
            if declared_members != actual_members:
                raise ValueError(f"Undeclared output files from {step_id}")
            if (
                sha256_file(spec_path) != run.spec_sha256
                or _code_identity(root, registry)[0] != code_hash
            ):
                raise ValueError("Specification or code changed during execution")
            if sha256_file(root / "uv.lock") != lock_hash:
                raise ValueError("Lockfile changed during execution")
            completed = step_run.model_copy(deep=True)
            completed.output_artifacts = [item.id for item in produced]
            completed.status = ExecutionStatus.completed
            completed.finished_at = datetime.now(UTC)
            _atomic_json(staging / "_step.json", completed.model_dump(mode="json"))
            destination = run_dir / "steps" / step_id
            staging.rename(destination)
            step_run.output_artifacts = completed.output_artifacts
            step_run.status = completed.status
            step_run.finished_at = completed.finished_at
            for artifact in produced:
                artifact.path = (
                    Path("steps") / step_id / products[artifact.product_id].filename
                ).as_posix()
                context[artifact.product_id] = run_dir / artifact.path
            run.artifacts = [*(run.artifacts or []), *produced]
            run.outputs = [*(run.outputs or []), *(item.id for item in produced)]
            _save_run(run)
        run.status = ExecutionStatus.completed
        run.finished_at = datetime.now(UTC)
        _validate_run_record(run, Path(run.manifest_path), verify_inputs=True)
        _save_run(run)
        return run
    except BaseException as error:
        run.status = ExecutionStatus.failed
        run.finished_at = datetime.now(UTC)
        run.error = f"{type(error).__name__}: {error}"
        for step_run in run.steps or []:
            if step_run.status == "running":
                step_run.status = ExecutionStatus.failed
                step_run.finished_at = run.finished_at
                step_run.error = run.error
        _save_run(run)
        if isinstance(error, (KeyboardInterrupt, SystemExit)):
            raise
        raise WorkflowExecutionError(run.error, run) from error


def validate_run(run_path: str | Path, verify_inputs: bool = True) -> WorkflowRun:
    """Require a complete, internally consistent run and verify artifact bytes.

    The embedded definition and recorded execution fingerprint are authoritative;
    current checkout code/spec files need not remain available to read an old run.
    This validates evidence consistency, not cryptographic authorship.
    """
    manifest = Path(run_path).resolve()
    run = WorkflowRun.model_validate_json(manifest.read_text(encoding="utf-8"))
    run.manifest_path = str(manifest)
    return _validate_run_record(run, manifest, verify_inputs)


def _validate_run_record(run: WorkflowRun, manifest: Path, verify_inputs: bool) -> WorkflowRun:
    order = validate_definition(run.definition)
    if run.status != "completed" or run.finished_at is None:
        raise ValueError(f"Run is not completed: {run.status}")
    if run.definition_sha256 != content_hash(run.definition.model_dump(mode="json")):
        raise ValueError("Embedded workflow definition checksum mismatch")
    if run.code_sha256 != content_hash(
        [item.model_dump(mode="json") for item in run.code_inventory]
    ):
        raise ValueError("Code inventory checksum mismatch")
    unique([item.path for item in run.code_inventory], "code inventory paths")
    unique([item.name for item in run.parameters or []], "resolved parameter names")
    parameters = {item.name: item.value for item in run.parameters or []}
    parameter_defs = {item.name: item for item in run.definition.parameters or []}
    if set(parameters) - set(parameter_defs):
        raise ValueError("Undeclared resolved parameter")
    if any(item.required and item.name not in parameters for item in parameter_defs.values()):
        raise ValueError("Missing required resolved parameter")
    unique([item.id for item in run.artifacts or []], "artifact identifiers")
    unique([item.product_id for item in run.artifacts or []], "artifact product bindings")
    unique([item.id for item in run.steps or []], "execution identifiers")
    unique([item.step_id for item in run.steps or []], "executed step bindings")
    products = {item.id: _effective_product(item, parameters) for item in run.definition.products}
    artifacts = {item.id: item for item in run.artifacts or []}
    by_product = {item.product_id: item for item in artifacts.values()}
    executions = {item.step_id: item for item in run.steps or []}
    if set(by_product) != set(products) or set(executions) != set(order):
        raise ValueError("Run lacks required product or step bindings")
    produced_ids: list[str] = []
    for definition in run.definition.steps:
        step = executions[definition.id]
        if step.status != "completed" or step.finished_at is None:
            raise ValueError(f"Uncompleted producer: {step.step_id}")
        if step.id != f"{run.run_id}:step:{step.step_id}":
            raise ValueError("Step execution identity does not match run")
        expected_inputs = {by_product[key].id for key in definition.inputs}
        expected_outputs = {by_product[key].id for key in definition.outputs}
        unique(step.input_artifacts or [], "step input bindings")
        unique(step.output_artifacts or [], "step output bindings")
        if set(step.input_artifacts or []) != expected_inputs:
            raise ValueError(f"Incorrect input lineage for {step.step_id}")
        if set(step.output_artifacts or []) != expected_outputs:
            raise ValueError(f"Incorrect output lineage for {step.step_id}")
        expected_parameters = [
            item.model_dump(mode="json")
            for item in run.parameters or []
            if item.name in (definition.parameters or [])
        ]
        if [item.model_dump(mode="json") for item in step.parameters or []] != expected_parameters:
            raise ValueError(f"Incorrect resolved parameters for {step.step_id}")
        publication = manifest.parent / "steps" / definition.id / "_step.json"
        published = StepExecution.model_validate_json(publication.read_text(encoding="utf-8"))
        if published != step:
            raise ValueError(f"Step publication record mismatch: {step.step_id}")
        produced_ids.extend(step.output_artifacts or [])
        for key in definition.outputs:
            artifact = by_product[key]
            expected_path = (Path("steps") / definition.id / products[key].filename).as_posix()
            if artifact.producer_step != step.id or artifact.path != expected_path:
                raise ValueError(f"Incorrect producer or output location for {key}")
    unique(run.outputs or [], "run output bindings")
    if set(run.outputs or []) != set(produced_ids):
        raise ValueError("Run output bindings do not match producing steps")
    for artifact in artifacts.values():
        product = products[artifact.product_id]
        if artifact.run_id != run.run_id or artifact.id != f"{run.run_id}:artifact:{product.id}":
            raise ValueError(f"Artifact identity does not match run: {product.id}")
        if product.external:
            if artifact.producer_step or not Path(artifact.path).is_absolute():
                raise ValueError(f"External input has invalid lineage: {product.id}")
        else:
            safe_relative_path(artifact.path)
        for field in ("unit", "cohort", "feature_policy", "coverage_policy"):
            if getattr(artifact, field) != getattr(product, field):
                raise ValueError(f"Artifact {field} differs from effective contract: {product.id}")
        if not artifact.validations or any(
            item.status != "passed" for item in artifact.validations
        ):
            raise ValueError(f"Artifact has failed or missing validation: {product.id}")
        if product.kind != "structural" and not any(
            item.validator != "artifact_inventory_v1" for item in artifact.validations
        ):
            raise ValueError(f"Artifact lacks scientific validation: {product.id}")
        unique([item.path for item in artifact.inventory], "artifact inventory paths")
        if verify_inputs or not product.external:
            _verify_artifact(run, artifact)
    if _fingerprint(run) != run.execution_fingerprint:
        raise ValueError("Execution fingerprint mismatch")
    return run
