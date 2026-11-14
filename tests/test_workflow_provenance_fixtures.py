"""Routing fixtures retain the production model and real trait-policy resolver (#1359)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from mtp.datamodel.workflow_pydantic import WorkflowRun
from mtp.workflow.runtime import validate_run
from mtp.workflow.trait_policy import resolve_trait_policy
from workflow_provenance_fixtures import HISTORICAL_CORE, routing_core


def test_routing_fixture_uses_a_real_historical_raw_policy_record(tmp_path: Path) -> None:
    product = tmp_path / "count.tsv"
    product.write_text("synthetic routing bytes\n", encoding="utf-8")
    core = routing_core(tmp_path / "core.json", {"count_profiles": product})
    assert isinstance(core, WorkflowRun)
    assert core.trait_key_policy is None
    resolved = resolve_trait_policy(core)
    assert resolved.policy == "raw_object_v1"
    assert resolved.include_policy_metadata is False
    assert core.artifacts is not None
    assert core.artifacts[0].path == str(product.resolve())


def test_typed_core_load_refuses_missing_code_inventory(tmp_path: Path) -> None:
    payload = json.loads(HISTORICAL_CORE.read_bytes())
    del payload["code_inventory"]
    manifest = tmp_path / "missing-inventory.json"
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    # validate_run parses the required WorkflowRun fields before any scientific
    # verification or policy resolver can see this malformed declaration.
    with pytest.raises(ValidationError, match="code_inventory"):
        validate_run(manifest)
