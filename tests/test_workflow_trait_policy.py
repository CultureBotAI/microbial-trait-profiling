"""Bind historical and explicit trait-key semantics to core records (#1129)."""

from __future__ import annotations

import json
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest

from mtp.datamodel.workflow_pydantic import ArtifactInstance, WorkflowRun
from mtp.workflow import verification
from mtp.workflow.runtime import content_hash
from mtp.workflow.trait_policy import resolve_trait_policy

ROOT = Path(__file__).resolve().parents[1]
HISTORICAL_FIXTURE = Path("tests/fixtures/trait_key_policy/historical_implementations.json")
MATRIX_SOURCE = "src/mtp/taxa_trait_matrix.py"
RAW = "raw_object_v1"
UNION = "facultative_union_v1"


def _historical_inventory(commit: str) -> dict[str, Any]:
    data = json.loads((ROOT / HISTORICAL_FIXTURE).read_text(encoding="utf-8"))
    assert data["format"] == "historical-trait-policy-inventories-v1"
    matches = [record for record in data["records"] if record["source_reference"] == commit]
    assert len(matches) == 1
    return dict(matches[0]["code_inventory"])


def test_historical_inventory_fixture_does_not_require_git_history(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = (ROOT / HISTORICAL_FIXTURE).read_bytes()
    target = tmp_path / HISTORICAL_FIXTURE
    target.parent.mkdir(parents=True)
    target.write_bytes(payload)
    records = json.loads(payload)["records"]
    monkeypatch.setitem(globals(), "ROOT", tmp_path)
    for record in records:
        assert _historical_inventory(record["source_reference"]) == record["code_inventory"]


def _record(*, commit: str = "511c891", policy: str | None = None) -> dict[str, Any]:
    """A schema-valid unit record; full execution validation has separate integration tests."""
    parameters = [{"name": "strain_rollup", "value": "off", "source": "override"}]
    definitions: list[dict[str, Any]] = [{"name": "strain_rollup", "default": "phenotype"}]
    bound = ["strain_rollup"]
    if policy is not None:
        parameters.append({"name": "trait_key_policy", "value": policy, "source": "default"})
        definitions.append({"name": "trait_key_policy", "default": policy, "required": True})
        bound.append("trait_key_policy")
    products = [
        {
            "id": identity,
            "kind": identity,
            "storage_kind": "file",
            "format": "tsv",
            "filename": f"{identity}.tsv",
            "external": identity == "kg_edges",
        }
        for identity in ("kg_edges", "taxa_matrix", "feature_metadata", "feature_catalog")
    ]
    steps = [
        {
            "id": "matrix",
            "operation": "taxa_matrix",
            "inputs": ["kg_edges"],
            "outputs": ["taxa_matrix", "feature_metadata"],
            "parameters": list(bound),
        },
        {
            "id": "catalog",
            "operation": "feature_catalog",
            "inputs": ["taxa_matrix", "feature_metadata"],
            "outputs": ["feature_catalog"],
            "parameters": list(bound),
        },
    ]
    definition = {
        "id": "policy-test",
        "version": "1",
        "products": products,
        "steps": steps,
        "parameters": definitions,
    }
    inventory = [_historical_inventory(commit)]
    if policy is not None:
        # Explicit-policy producers need not appear in the legacy compatibility registry.
        inventory = [{"path": MATRIX_SOURCE, "sha256": "c" * 64, "size_bytes": 30000}]
    record: dict[str, Any] = {
        "run_id": "run-policy-test",
        "definition": definition,
        "status": "completed",
        "started_at": "2026-10-27T00:00:00Z",
        "finished_at": "2026-10-27T00:00:01Z",
        "manifest_path": "/durable/policy-test/run.json",
        "spec_path": "/source/workflows/bugsigdb_core.yaml",
        "spec_sha256": "a" * 64,
        "definition_sha256": content_hash(definition),
        "execution_fingerprint": "b" * 64,
        "code_sha256": content_hash(inventory),
        "code_inventory": inventory,
        "lock_sha256": "d" * 64,
        "git_commit": commit,
        "git_dirty": False,
        "runtime_json": "{}",
        "parameters": parameters,
        "steps": [
            {
                "id": f"run-policy-test:step:{step['id']}",
                "step_id": step["id"],
                "status": "completed",
                "started_at": "2026-10-27T00:00:00Z",
                "finished_at": "2026-10-27T00:00:01Z",
                "parameters": deepcopy(parameters),
                "input_artifacts": [],
                "output_artifacts": [],
            }
            for step in steps
        ],
        "artifacts": [],
        "outputs": [],
    }
    if policy is not None:
        record["trait_key_policy"] = policy
    return record


@pytest.mark.parametrize(("commit", "policy"), [("511c891", RAW), ("0e49e87", UNION)])
def test_historical_policy_uses_exact_recorded_source_and_preserves_record(
    commit: str, policy: str
) -> None:
    run = WorkflowRun.model_validate(_record(commit=commit))
    before = run.model_dump(mode="json")
    resolved = resolve_trait_policy(run)
    assert resolved.policy == policy
    assert resolved.include_policy_metadata is False
    assert run.model_dump(mode="json") == before


@pytest.mark.parametrize("policy", [RAW, UNION])
def test_explicit_policy_selects_extended_metadata(policy: str) -> None:
    run = WorkflowRun.model_validate(_record(policy=policy))
    resolved = resolve_trait_policy(run)
    assert resolved.policy == policy
    assert resolved.include_policy_metadata is True


@pytest.mark.parametrize("commit", ["511c891", "0e49e87"])
@pytest.mark.parametrize("policy", [RAW, UNION])
def test_explicit_policy_cannot_relabel_an_implicit_policy_producer(
    commit: str, policy: str
) -> None:
    record = _record(policy=policy)
    record["code_inventory"] = [_historical_inventory(commit)]
    record["code_sha256"] = content_hash(record["code_inventory"])
    with pytest.raises(ValueError, match=r"(?i)policy|implementation"):
        resolve_trait_policy(WorkflowRun.model_validate(record))


@pytest.mark.parametrize("damage", ["unknown_hash", "wrong_size", "missing", "duplicate", "path"])
def test_historical_unknown_or_ambiguous_inventory_fails_closed(damage: str) -> None:
    record = _record()
    entry = record["code_inventory"][0]
    if damage == "unknown_hash":
        entry["sha256"] = "f" * 64
    elif damage == "wrong_size":
        entry["size_bytes"] += 1
    elif damage == "missing":
        record["code_inventory"] = []
    elif damage == "duplicate":
        record["code_inventory"].append(deepcopy(entry))
    else:
        entry["path"] = "src/bsdbng/taxa_trait_matrix.py"
    record["code_sha256"] = content_hash(record["code_inventory"])
    with pytest.raises(ValueError, match=r"(?i)policy|inventory|implementation"):
        resolve_trait_policy(WorkflowRun.model_validate(record))


@pytest.mark.parametrize(
    "damage",
    [
        "field_absent",
        "parameter_absent",
        "declaration_absent",
        "parameter_duplicate",
        "parameter_disagrees",
        "matrix_binding_absent",
        "catalog_binding_absent",
        "execution_parameter_absent",
        "execution_parameter_disagrees",
    ],
)
def test_partial_or_inconsistent_explicit_policy_rejected(damage: str) -> None:
    record = _record(policy=UNION)
    if damage == "field_absent":
        del record["trait_key_policy"]
    elif damage == "parameter_absent":
        record["parameters"] = record["parameters"][:1]
    elif damage == "declaration_absent":
        record["definition"]["parameters"] = record["definition"]["parameters"][:1]
    elif damage == "parameter_duplicate":
        record["parameters"].append(deepcopy(record["parameters"][-1]))
    elif damage == "parameter_disagrees":
        record["parameters"][-1]["value"] = RAW
    elif damage in {"matrix_binding_absent", "catalog_binding_absent"}:
        index = 0 if damage == "matrix_binding_absent" else 1
        record["definition"]["steps"][index]["parameters"] = ["strain_rollup"]
    elif damage == "execution_parameter_absent":
        record["steps"][0]["parameters"] = record["steps"][0]["parameters"][:1]
    else:
        record["steps"][0]["parameters"][-1]["value"] = RAW
    with pytest.raises(ValueError, match=r"(?i)policy|parameter|binding"):
        resolve_trait_policy(WorkflowRun.model_validate(record))


@pytest.mark.parametrize("location", ["definition", "execution"])
def test_stray_policy_binding_cannot_be_treated_as_historical(location: str) -> None:
    record = _record()
    if location == "definition":
        record["definition"]["steps"][0]["parameters"].append("trait_key_policy")
    else:
        record["steps"][0]["parameters"].append(
            {"name": "trait_key_policy", "value": UNION, "source": "default"}
        )
    with pytest.raises(ValueError, match=r"(?i)policy|parameter|binding"):
        resolve_trait_policy(WorkflowRun.model_validate(record))


def test_historical_resolution_is_not_commit_date_or_current_checkout_dependent() -> None:
    record = _record()
    record["git_commit"] = "0e49e87"
    record["git_dirty"] = True
    resolved = resolve_trait_policy(WorkflowRun.model_validate(record))
    assert resolved.policy == RAW
    assert resolved.include_policy_metadata is False


def test_verifier_passes_historical_policy_without_skipping_product_checks(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    record = _record()
    run = WorkflowRun.model_validate(record)
    calls: list[tuple[str, Mapping[str, Any]]] = []

    def validated(path: Path, **kwargs: Any) -> WorkflowRun:
        assert path == tmp_path / "run.json"
        calls.append(("validate_run", kwargs))
        return run

    def product_path(source: WorkflowRun, product_id: str, resolver: Any = None) -> Path:
        assert source is run
        assert resolver is None
        return tmp_path / f"{product_id}.tsv"

    def validate_product(product: Any, path: Path, context: Any, **kwargs: Any) -> list[Any]:
        assert path == context[product.id] == tmp_path / f"{product.id}.tsv"
        calls.append((product.id, kwargs))
        return []

    # Structural/inventory validation is deliberately a seam here; integration tests
    # exercise real files. All declared scientific products must still be checked.
    run.artifacts = []
    for product in run.definition.products:
        run.artifacts.append(
            ArtifactInstance(
                id=f"run-policy-test:artifact:{product.id}",
                product_id=product.id,
                run_id=run.run_id,
                path=str(tmp_path / f"{product.id}.tsv"),
                sha256="a" * 64,
                size_bytes=1,
                inventory=[],
                validations=[],
            )
        )
    monkeypatch.setattr(verification, "validate_run", validated)
    monkeypatch.setattr(verification, "artifact_path", product_path)
    monkeypatch.setattr(verification, "validate_core_product", validate_product)
    assert verification.verify_core_run(tmp_path / "run.json") is run
    assert [name for name, _ in calls] == [
        "validate_run",
        *[product.id for product in run.definition.products],
    ]
    for _, kwargs in calls[1:]:
        assert kwargs["trait_key_policy"] == RAW
        assert kwargs["include_policy_metadata"] is False
