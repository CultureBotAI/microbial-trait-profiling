"""New arm plans preserve core trait keys without relabelling historical plans (#1359)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path
from typing import NoReturn

import pytest

from mtp.datamodel.workflow_pydantic import (
    ArmRun,
    ExecutionStatus,
    FileInventoryEntry,
    ParameterValue,
    WorkflowRun,
)
from mtp.paths import REPO_ROOT
from mtp.trait_key_policy import TraitKeyPolicy
from mtp.workflow.arms import runtime
from mtp.workflow.arms.noel import build_commands
from mtp.workflow.arms.types import ArmContext, Command
from mtp.workflow.consumers import CoreBundle
from mtp.workflow.runtime import sha256_file
from mtp.workflow.trait_policy import HISTORICAL_MATRIX_POLICIES, MATRIX_SOURCE

RAW: TraitKeyPolicy = "raw_object_v1"
UNION: TraitKeyPolicy = "facultative_union_v1"
HISTORICAL_CORE = REPO_ROOT / "tests/fixtures/trait_key_policy/historical_raw_core/core/run.json"


def _core(policy: TraitKeyPolicy) -> WorkflowRun:
    """Exercise the real audited-producer resolver, not an invented column heuristic."""
    run = WorkflowRun.model_validate_json(HISTORICAL_CORE.read_bytes())
    identity = next(key for key, value in HISTORICAL_MATRIX_POLICIES.items() if value == policy)
    run.code_inventory = [
        FileInventoryEntry(path=MATRIX_SOURCE, sha256=identity[0], size_bytes=identity[1])
    ]
    return run


def _context(tmp_path: Path, policy: TraitKeyPolicy | None) -> ArmContext:
    products = {
        product: tmp_path / "core" / filename for product, filename in runtime.INPUT_NAMES.items()
    }
    return ArmContext(
        repository=REPO_ROOT,
        manifest=tmp_path / "core/run.json",
        output=tmp_path / "arm",
        inputs=tmp_path / "arm/inputs",
        products=products,
        code_commit="fixture",
        replicates=None,
        seed=None,
        figures=False,
        trait_key_policy=policy,
    )


def _ledger(context: ArmContext) -> Command:
    return next(
        command for command in build_commands(context) if command.id == "taxon_contributions"
    )


def _without_policy(command: Command) -> Command:
    args = list(command.argv)
    index = args.index("--trait-key-policy")
    del args[index : index + 2]
    return replace(command, argv=tuple(args))


def _pairs(argv: Sequence[str]) -> tuple[tuple[str, str], ...]:
    return tuple((argv[i], argv[i + 1]) for i in range(len(argv) - 1))


@pytest.mark.parametrize("policy", [RAW, UNION])
def test_noel_ledger_names_the_resolved_core_policy(tmp_path: Path, policy: TraitKeyPolicy) -> None:
    context = _context(tmp_path, policy)
    ledger = _ledger(context)
    assert ("--trait-key-policy", policy) in _pairs(ledger.argv)
    runtime._validate_plan([ledger], context)


@pytest.mark.parametrize(
    "options",
    [
        ("--trait-key-policy", UNION),
        (f"--trait-key-policy={UNION}",),
        ("--trait-key-p", UNION),
        ("--trait-key-policy", RAW, "--trait-key-policy", UNION),
        ("--trait-key-policy", "unknown"),
        ("--trait-key-policy",),
        (),
    ],
)
def test_new_plan_checks_effective_policy_with_the_execution_parser(
    tmp_path: Path, options: tuple[str, ...]
) -> None:
    context = _context(tmp_path, RAW)
    ledger = _without_policy(_ledger(context))
    ledger = replace(ledger, argv=(*ledger.argv, *options))
    with pytest.raises(ValueError, match=r"trait.key.policy|trait-key-policy"):
        runtime._validate_plan([ledger], context)


@pytest.mark.parametrize(
    "options",
    [
        (f"--trait-key-policy={RAW}",),
        ("--trait-key-policy", UNION, "--trait-key-policy", RAW),
    ],
)
def test_plan_accepts_the_policy_the_execution_parser_will_use(
    tmp_path: Path, options: tuple[str, ...]
) -> None:
    context = _context(tmp_path, RAW)
    ledger = _without_policy(_ledger(context))
    runtime._validate_plan([replace(ledger, argv=(*ledger.argv, *options))], context)


def test_historical_plan_absence_remains_unrecorded(tmp_path: Path) -> None:
    context = _context(tmp_path, None)
    ledger = _ledger(context)
    assert "--trait-key-policy" not in ledger.argv
    assert context.trait_key_policy is None
    runtime._validate_plan([ledger], context)


def test_policy_argument_requires_a_matching_arm_declaration(tmp_path: Path) -> None:
    context = _context(tmp_path, None)
    ledger = _ledger(context)
    ledger = replace(ledger, argv=(*ledger.argv, "--trait-key-policy", RAW))
    with pytest.raises(ValueError, match="trait-key-policy"):
        runtime._validate_plan([ledger], context)


class _RecordBoundaryReached(AssertionError):
    def __init__(self) -> None:
        super().__init__("verification guard bypassed before evidence checks")


@pytest.mark.parametrize("policy", [RAW, UNION])
def test_execution_resolves_and_records_the_verified_core_policy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, policy: TraitKeyPolicy
) -> None:
    # Isolate expensive core verification and subprocess execution. Policy resolution,
    # plan construction/validation, and creation of the execution record are real.
    core = _core(policy)
    bundle = CoreBundle(HISTORICAL_CORE, core, "count_profiles", "all")
    monkeypatch.setattr(runtime, "open_core_bundle", lambda *_args: bundle)
    monkeypatch.setattr(runtime, "preflight_consumer_outputs", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(runtime, "repository_state", lambda _root: ("fixture", False))
    monkeypatch.setattr(runtime, "_code_inventory", lambda _root: [])
    contexts: list[ArmContext] = []
    records: list[ArmRun] = []

    def plan(context: ArmContext) -> list[Command]:
        contexts.append(context)
        return [_ledger(context)]

    def stop_at_record(record: ArmRun) -> NoReturn:
        records.append(record)
        raise _RecordBoundaryReached

    monkeypatch.setattr(runtime, "_save", stop_at_record)
    with pytest.raises(_RecordBoundaryReached):
        runtime.execute_arm(
            "noel", HISTORICAL_CORE, tmp_path / "output", plan, repository=REPO_ROOT
        )
    assert contexts[0].trait_key_policy == policy
    assert [
        (item.value, item.source)
        for item in records[0].parameters
        if item.name == "trait_key_policy"
    ] == [(policy, "verified_core")]


def _record(core: WorkflowRun, tmp_path: Path, values: Sequence[str]) -> ArmRun:
    return ArmRun(
        id="fixture-arm",
        arm="noel",
        scope="standard_analysis_calibration_and_optional_figures",
        status=ExecutionStatus.completed,
        started_at=datetime.now(UTC),
        ended_at=datetime.now(UTC),
        manifest_path=str(tmp_path / "arm_run.json"),
        core_manifest_path=core.manifest_path,
        core_manifest_sha256=sha256_file(HISTORICAL_CORE),
        core_run_id=core.run_id,
        profile_product="count_profiles",
        contrast_population=core.contrast_population,
        direction_basis=core.direction_basis,
        git_commit="fixture",
        git_dirty=False,
        repository_path=str(REPO_ROOT),
        code_inventory=[],
        runtime_json="{}",
        parameters=[
            ParameterValue(name="trait_key_policy", value=value, source="verified_core")
            for value in values
        ],
        plan=[],
        plan_sha256=runtime._plan_hash([]),
    )


def _stop_before_plan_evidence(monkeypatch: pytest.MonkeyPatch, contexts: list[ArmContext]) -> None:
    def stop(_commands: Sequence[Command], context: ArmContext, **_kwargs: object) -> NoReturn:
        contexts.append(context)
        raise _RecordBoundaryReached

    monkeypatch.setattr(runtime, "_validate_plan", stop)


@pytest.mark.parametrize("policy", [RAW, UNION])
@pytest.mark.parametrize("declared", [False, True])
def test_verifier_only_names_policy_when_the_arm_record_declares_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, policy: TraitKeyPolicy, declared: bool
) -> None:
    core = _core(policy)
    record = _record(core, tmp_path, [policy] if declared else [])
    contexts: list[ArmContext] = []
    _stop_before_plan_evidence(monkeypatch, contexts)
    original = record.model_dump(mode="json")
    with pytest.raises(_RecordBoundaryReached):
        runtime._verify_arm_record(Path(record.manifest_path), record, core)
    assert contexts[0].trait_key_policy == (policy if declared else None)
    assert record.model_dump(mode="json") == original
    assert core.trait_key_policy is None


@pytest.mark.parametrize("values", [[UNION], [RAW, RAW], [RAW, UNION], ["unknown"]])
def test_verifier_refuses_contradictory_or_nonunique_recorded_policy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, values: list[str]
) -> None:
    core = _core(RAW)
    record = _record(core, tmp_path, values)
    contexts: list[ArmContext] = []
    _stop_before_plan_evidence(monkeypatch, contexts)
    with pytest.raises(ValueError, match=r"trait.key.policy|trait-key-policy"):
        runtime._verify_arm_record(Path(record.manifest_path), record, core)
    assert not contexts


@pytest.mark.parametrize("source", ["default", "override"])
def test_recorded_arm_policy_must_be_derived_from_the_verified_core(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, source: str
) -> None:
    core = _core(RAW)
    record = _record(core, tmp_path, [RAW])
    record.parameters[0].source = source
    contexts: list[ArmContext] = []
    _stop_before_plan_evidence(monkeypatch, contexts)
    with pytest.raises(ValueError, match="verified_core"):
        runtime._verify_arm_record(Path(record.manifest_path), record, core)
    assert not contexts


def test_raw_core_ledger_command_carries_raw_object_v1(tmp_path: Path) -> None:
    from mtp.workflow.trait_policy import resolve_trait_policy

    resolved_policy = resolve_trait_policy(_core(RAW)).policy
    context = _context(tmp_path, resolved_policy)
    ledger = _ledger(context)
    assert resolved_policy == RAW
    assert ("--trait-key-policy", RAW) in _pairs(ledger.argv)


def test_step_policy_parser_returns_none_when_flag_is_omitted(tmp_path: Path) -> None:
    ledger = _without_policy(_ledger(_context(tmp_path, RAW)))
    assert runtime._step_trait_key_policy(ledger) is None
