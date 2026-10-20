"""Adversarial execution and lineage checks for issue #629."""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
import yaml

from mtp.workflow import runtime
from mtp.workflow.definition import read_definition
from mtp.workflow.models import ProductDefinition, ValidationResult, ValidationStatus, WorkflowRun
from mtp.workflow.runtime import (
    Operation,
    WorkflowExecutionError,
    artifact_path,
    execute_workflow,
    validate_run,
)


def _copy(inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]) -> None:
    outputs["result"].write_text(
        inputs["raw"].read_text() + str(parameters.get("suffix", "")), encoding="utf-8"
    )


def _collection(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    _copy(inputs, {"result": outputs["result"] / "part.txt"}, parameters)
    (outputs["result"] / ".hidden").write_text("included", encoding="utf-8")


def _fail(inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]) -> None:
    _copy(inputs, outputs, parameters)
    raise RuntimeError("deliberate failure")


def _mutate_input(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    _copy(inputs, outputs, parameters)
    inputs["raw"].write_text("changed", encoding="utf-8")


def _extra_output(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    _copy(inputs, outputs, parameters)
    outputs["result"].with_name("undeclared.txt").write_text("extra", encoding="utf-8")


def _symlink_output(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    del parameters
    outputs["result"].symlink_to(inputs["raw"])


def _second(
    inputs: dict[str, Path], outputs: dict[str, Path], parameters: dict[str, object]
) -> None:
    _copy(inputs, {"result": outputs["second"]}, parameters)


@dataclass
class Harness:
    source: Path
    spec_path: Path
    raw: Path
    output_root: Path
    spec: dict[str, Any]

    def save(self) -> None:
        self.spec_path.write_text(yaml.safe_dump(self.spec), encoding="utf-8")

    def execute(
        self, operation: Operation = _copy, overrides: dict[str, object] | None = None
    ) -> WorkflowRun:
        self.save()
        return execute_workflow(
            self.spec_path,
            {"raw": self.raw},
            self.output_root,
            {"copy": operation},
            parameter_overrides=overrides,
            repository_root=self.source,
        )


@pytest.fixture
def harness(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Harness:
    source = tmp_path / "source"
    (source / "src").mkdir(parents=True)
    (source / "src" / "version.py").write_text("VERSION = 1\n", encoding="utf-8")
    (source / "uv.lock").write_text("version = 1\n", encoding="utf-8")
    raw = tmp_path / "input.txt"
    raw.write_text("alpha", encoding="utf-8")

    def fake_git(root: Path) -> tuple[str, bool]:
        assert root == source
        return "fixture-commit", False

    monkeypatch.setattr(runtime, "repository_state", fake_git)
    spec: dict[str, Any] = {
        "id": "fixture",
        "version": "1",
        "products": [
            {
                "id": "raw",
                "kind": "structural",
                "storage_kind": "file",
                "format": "txt",
                "filename": "input.txt",
                "external": True,
            },
            {
                "id": "result",
                "kind": "structural",
                "storage_kind": "file",
                "format": "txt",
                "filename": "result.txt",
                "unit": "annotation",
                "cohort": "all",
                "cohort_parameter": "host",
                "feature_policy": "all_features",
            },
        ],
        "steps": [
            {
                "id": "copy",
                "operation": "copy",
                "inputs": ["raw"],
                "outputs": ["result"],
                "parameters": ["suffix", "host"],
            }
        ],
        "parameters": [
            {"name": "suffix", "default": "!", "environment": "MTP_TEST_SUFFIX"},
            {"name": "host", "default": "all", "required": True},
        ],
    }
    monkeypatch.delenv("MTP_TEST_SUFFIX", raising=False)
    result = Harness(source, tmp_path / "workflow.yaml", raw, tmp_path / "runs", spec)
    result.save()
    return result


def _rewrite(run: WorkflowRun, value: dict[str, Any]) -> None:
    Path(run.manifest_path).write_text(json.dumps(value), encoding="utf-8")


def test_attempt_identity_distinct_from_deterministic_fingerprint(harness: Harness) -> None:
    first, second = harness.execute(), harness.execute()
    assert first.run_id != second.run_id
    assert first.execution_fingerprint == second.execution_fingerprint
    assert first.definition_sha256 == second.definition_sha256
    assert artifact_path(first, "result").read_text() == "alpha!"
    assert validate_run(first.manifest_path) == first
    assert not list((Path(first.manifest_path).parent / ".staging").iterdir())


def test_parameter_environment_and_override_change_effective_cohort(
    harness: Harness, monkeypatch: pytest.MonkeyPatch
) -> None:
    baseline = harness.execute()
    monkeypatch.setenv("MTP_TEST_SUFFIX", "environment")
    environment = harness.execute()
    override = harness.execute(overrides={"suffix": "override", "host": "Homo sapiens"})
    assert artifact_path(environment, "result").read_text() == "alphaenvironment"
    assert artifact_path(override, "result").read_text() == "alphaoverride"
    assert len({item.execution_fingerprint for item in (baseline, environment, override)}) == 3
    assert (
        next(item for item in override.artifacts or [] if item.product_id == "result").cohort
        == "Homo sapiens"
    )
    assert any(
        item.source == "environment:MTP_TEST_SUFFIX" for item in environment.parameters or []
    )


@pytest.mark.parametrize("target", ["input", "code", "lock", "definition"])
def test_fingerprint_changes_for_actual_semantic_inputs(harness: Harness, target: str) -> None:
    first = harness.execute()
    if target == "input":
        harness.raw.write_text("beta", encoding="utf-8")
    elif target == "code":
        (harness.source / "src" / "version.py").write_text("VERSION = 2\n", encoding="utf-8")
    elif target == "lock":
        (harness.source / "uv.lock").write_text("version = 2\n", encoding="utf-8")
    else:
        harness.spec["products"][1]["feature_policy"] = "retained_traits"
    assert harness.execute().execution_fingerprint != first.execution_fingerprint


@pytest.mark.parametrize("operation", [_fail, _mutate_input, _extra_output, _symlink_output])
def test_failed_steps_preserve_diagnostics_and_never_publish(
    harness: Harness, operation: Operation
) -> None:
    with pytest.raises(WorkflowExecutionError) as raised:
        harness.execute(operation)
    run = raised.value.run
    assert run.status == "failed"
    assert Path(run.manifest_path).exists()
    assert not list((Path(run.manifest_path).parent / "steps").iterdir())
    assert not run.outputs
    assert (run.steps or [])[0].status == "failed"
    with pytest.raises(ValueError, match="not completed"):
        validate_run(run.manifest_path)


def test_interruption_before_atomic_publication_is_not_consumable(
    harness: Harness, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = Path.rename

    def fail_rename(path: Path, target: str | Path) -> Path:
        if Path(target).parent.name == "steps":
            raise OSError("publication interrupted")
        return original(path, target)

    monkeypatch.setattr(Path, "rename", fail_rename)
    with pytest.raises(WorkflowExecutionError, match="publication interrupted") as raised:
        harness.execute()
    run = raised.value.run
    assert run.status == "failed"
    assert not run.outputs
    assert not list((Path(run.manifest_path).parent / "steps").iterdir())
    with pytest.raises(ValueError, match="not completed"):
        validate_run(run.manifest_path)


@pytest.mark.parametrize("change", ["add", "remove", "rename"])
def test_collection_inventory_detects_member_changes(harness: Harness, change: str) -> None:
    harness.spec["products"][1]["storage_kind"] = "collection"
    run = harness.execute(_collection)
    path = artifact_path(run, "result")
    artifact = next(item for item in run.artifacts or [] if item.product_id == "result")
    assert {item.path for item in artifact.inventory} == {"part.txt", ".hidden"}
    if change == "add":
        (path / "extra.txt").write_text("extra", encoding="utf-8")
    elif change == "remove":
        (path / ".hidden").unlink()
    else:
        (path / "part.txt").rename(path / "other.txt")
    with pytest.raises(ValueError, match="changed"):
        validate_run(run.manifest_path)


def test_external_symlink_resolves_actual_bytes_and_retargeting_is_detected(
    harness: Harness,
) -> None:
    original = harness.raw
    binding = original.with_name("binding.txt")
    binding.symlink_to(original)
    harness.raw = binding
    run = harness.execute()
    source = next(item for item in run.artifacts or [] if item.product_id == "raw")
    assert source.path == str(original.resolve())
    assert source.requested_path == str(binding)
    other = original.with_name("other.txt")
    other.write_bytes(original.read_bytes())
    binding.unlink()
    binding.symlink_to(other)
    with pytest.raises(ValueError, match="binding target changed"):
        validate_run(run.manifest_path)


def test_changed_external_bytes_rejected_unless_verification_explicitly_disabled(
    harness: Harness,
) -> None:
    run = harness.execute()
    harness.raw.write_text("altered", encoding="utf-8")
    with pytest.raises(ValueError, match="changed"):
        validate_run(run.manifest_path)
    assert validate_run(run.manifest_path, verify_inputs=False).status == "completed"


@pytest.mark.parametrize(
    "mutation", ["cohort", "unit", "producer", "definition", "missing", "lineage", "fingerprint"]
)
def test_structurally_valid_false_provenance_is_rejected(harness: Harness, mutation: str) -> None:
    run = harness.execute()
    value = run.model_dump(mode="json")
    output = next(item for item in value["artifacts"] if item["product_id"] == "result")
    if mutation in {"cohort", "unit"}:
        output[mutation] = "incorrect"
    elif mutation == "producer":
        output["producer_step"] = "another-run:copy"
    elif mutation == "definition":
        value["definition"]["products"][1]["kind"] = "different_contract"
    elif mutation == "missing":
        value["artifacts"] = value["artifacts"][:1]
    elif mutation == "lineage":
        value["steps"][0]["input_artifacts"] = []
    else:
        value["execution_fingerprint"] = "0" * 64
    _rewrite(run, value)
    with pytest.raises(ValueError, match=r"differs|Incorrect|checksum|lacks|fingerprint"):
        validate_run(run.manifest_path)


def test_step_publication_and_output_bytes_are_verified(harness: Harness) -> None:
    run = harness.execute()
    step_path = Path(run.manifest_path).parent / "steps" / "copy" / "_step.json"
    step = json.loads(step_path.read_text())
    step["parameters"] = []
    step_path.write_text(json.dumps(step), encoding="utf-8")
    with pytest.raises(ValueError, match="publication record mismatch"):
        validate_run(run.manifest_path)
    step_path.write_text(json.dumps((run.steps or [])[0].model_dump(mode="json")), encoding="utf-8")
    artifact_path(run, "result").write_text("tampered", encoding="utf-8")
    with pytest.raises(ValueError, match="changed"):
        validate_run(run.manifest_path)


@pytest.mark.parametrize(
    "mutation", ["duplicate", "unknown", "cycle", "overlap", "unsafe", "producer", "parameter"]
)
def test_graph_contract_rejects_ambiguous_or_invalid_definitions(
    harness: Harness, mutation: str
) -> None:
    if mutation == "duplicate":
        harness.spec["products"].append(harness.spec["products"][0].copy())
    elif mutation == "unknown":
        harness.spec["steps"][0]["inputs"] = ["unknown"]
    elif mutation == "cycle":
        harness.spec["steps"][0]["depends_on"] = ["copy"]
    elif mutation == "overlap":
        product = {**harness.spec["products"][1], "id": "another", "filename": "result.txt/nested"}
        harness.spec["products"].append(product)
        harness.spec["steps"][0]["outputs"].append("another")
    elif mutation == "unsafe":
        harness.spec["products"][1]["filename"] = "../escape"
    elif mutation == "producer":
        harness.spec["steps"].append({**harness.spec["steps"][0], "id": "other"})
    else:
        harness.spec["steps"][0]["parameters"] = ["unknown"]
    harness.save()
    with pytest.raises(ValueError, match=r"Duplicate|Unknown|cycle|Overlapping|Unsafe"):
        read_definition(harness.spec_path)


def test_scientific_validator_is_required_and_failed_results_prevent_publication(
    harness: Harness,
) -> None:
    harness.spec["products"][1]["kind"] = "scientific"
    with pytest.raises(WorkflowExecutionError, match="Scientific validator required"):
        harness.execute()

    def reject(
        product: ProductDefinition, path: Path, context: Mapping[str, Path]
    ) -> list[ValidationResult]:
        assert path.exists()
        assert context["raw"].exists()
        return [
            ValidationResult(
                validator="fixture_science",
                status=(
                    ValidationStatus.failed if product.id == "result" else ValidationStatus.passed
                ),
                message="fixture",
            )
        ]

    with pytest.raises(WorkflowExecutionError, match="Product validation failed") as raised:
        execute_workflow(
            harness.spec_path,
            {"raw": harness.raw},
            harness.output_root,
            {"copy": _copy},
            validator=reject,
            repository_root=harness.source,
        )
    assert not raised.value.run.outputs


def test_undeclared_binding_and_override_rejected_before_execution(harness: Harness) -> None:
    with pytest.raises(ValueError, match="Undeclared parameter"):
        harness.execute(overrides={"unknown": "value"})
    with pytest.raises(ValueError, match="External bindings"):
        execute_workflow(harness.spec_path, {}, harness.output_root, {"copy": _copy})


def test_effective_cohort_whitespace_is_normalized_before_execution(harness: Harness) -> None:
    run = harness.execute(overrides={"host": " all "})
    assert next(item for item in run.parameters or [] if item.name == "host").value == "all"
    output = next(item for item in run.artifacts or [] if item.product_id == "result")
    assert output.cohort == "all"
    assert validate_run(run.manifest_path).status == "completed"


@pytest.mark.parametrize("filename", ["_step.json", "_step.json/child.tsv"])
def test_publication_record_path_is_reserved(harness: Harness, filename: str) -> None:
    harness.spec["products"][1]["filename"] = filename
    harness.save()
    with pytest.raises(ValueError, match="reserved"):
        read_definition(harness.spec_path)


def test_unused_parameter_cannot_claim_to_change_the_executed_analysis(harness: Harness) -> None:
    harness.spec["parameters"].append({"name": "annotation_rollup", "default": "direct_only"})
    harness.save()
    with pytest.raises(ValueError, match="Unused declared parameters"):
        read_definition(harness.spec_path)


def test_output_validator_receives_only_declared_inputs_and_sibling_outputs(
    harness: Harness,
) -> None:
    harness.spec["products"].append(
        {**harness.spec["products"][1], "id": "second", "filename": "second.txt"}
    )
    harness.spec["steps"].append(
        {
            "id": "second",
            "operation": "second",
            "inputs": ["raw"],
            "outputs": ["second"],
            "depends_on": ["copy"],
            "parameters": ["suffix", "host"],
        }
    )
    harness.save()
    seen: dict[str, set[str]] = {}

    def check(
        product: ProductDefinition, path: Path, context: Mapping[str, Path]
    ) -> list[ValidationResult]:
        assert path.exists()
        seen[product.id] = set(context)
        return [
            ValidationResult(
                validator="context_contract",
                status=ValidationStatus.passed,
                message="Explicit declared validation dependencies.",
            )
        ]

    run = execute_workflow(
        harness.spec_path,
        {"raw": harness.raw},
        harness.output_root,
        {"copy": _copy, "second": _second},
        validator=check,
        repository_root=harness.source,
    )
    assert run.status == "completed"
    assert seen == {"raw": {"raw"}, "result": {"raw", "result"}, "second": {"raw", "second"}}
