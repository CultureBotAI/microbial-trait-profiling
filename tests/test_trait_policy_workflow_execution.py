"""Execution-bound trait policies and non-bypass historical compatibility (#1129)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import artifact_path, execute_workflow, validate_run
from mtp.workflow.verification import verify_core_run

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "workflow"


@pytest.mark.parametrize("spec_name", ["bugsigdb_core.yaml", "bugsigdb_core_subset.yaml"])
@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_execution_binds_policy_to_run_steps_metadata_and_complete_verification(
    tmp_path: Path, spec_name: str, policy: str
) -> None:
    run = execute_workflow(
        ROOT / "workflows" / spec_name,
        {
            "raw_dump": FIXTURE / "full_dump.csv",
            "kg_edges": FIXTURE / "merged-kg_edges.tsv",
            "kg_nodes": FIXTURE / "merged-kg_nodes.tsv",
        },
        tmp_path / "runs",
        OPERATIONS,
        parameter_overrides={"source_mode": "fixture", "trait_key_policy": policy},
        validator=validate_core_product,
        repository_root=ROOT,
    )
    assert run.trait_key_policy == policy
    assert {item.name: item.value for item in run.parameters or []}["trait_key_policy"] == policy
    for step in run.steps or []:
        if step.step_id in {"taxa_matrix", "feature_catalog"}:
            assert {item.name: item.value for item in step.parameters or []}[
                "trait_key_policy"
            ] == policy
    expected = ["METPO:1000605", "METPO:1000608", "METPO:1000612"]
    for product in ("feature_metadata", "feature_catalog"):
        with artifact_path(run, product).open(encoding="utf-8", newline="") as handle:
            records = list(csv.DictReader(handle, delimiter="\t"))
        focal = next(row for row in records if row["column_id"] == "METPO:1000608")
        assert focal["trait_key_policy"] == policy
        assert json.loads(focal["trait_member_ids_json"]) == (
            expected if policy == "facultative_union_v1" else ["METPO:1000608"]
        )
    verified = verify_core_run(Path(run.manifest_path))
    assert verified == run
    original = Path(run.manifest_path).read_bytes()
    raw = json.loads(original)
    raw["trait_key_policy"] = (
        "raw_object_v1" if policy == "facultative_union_v1" else "facultative_union_v1"
    )
    Path(run.manifest_path).write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match=r"(?i)policy"):
        validate_run(Path(run.manifest_path))
    Path(run.manifest_path).write_bytes(original)
    assert verify_core_run(Path(run.manifest_path)) == run
