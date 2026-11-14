"""Record a supplied-root decision without adding a machine locator (#1820)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from mtp.datamodel.workflow_pydantic import ArmRun, ComparisonRun
from mtp.paths import REPO_ROOT
from mtp.workflow.arms.runtime import execute_arm, verify_arm
from mtp.workflow.comparison import execute_comparison_run, verify_comparison_run
from test_arm_runtime import Harness, _copy_plan
from test_arm_runtime import harness as harness
from test_comparison_run import stage as stage

FIELD = "provenance_root_supplied"


@pytest.mark.parametrize("model", [ArmRun, ComparisonRun])
def test_schema_and_model_preserve_unknown_for_historical_records(
    model: type[ArmRun] | type[ComparisonRun],
) -> None:
    schema = yaml.safe_load((REPO_ROOT / "schema/workflow.yaml").read_text(encoding="utf-8"))
    attributes = schema["classes"][model.__name__]["attributes"]
    assert FIELD in attributes, "The portable execution declaration must be authoritative LinkML"
    assert attributes[FIELD]["range"] == "boolean"
    assert not attributes[FIELD].get("required", False)
    assert "ifabsent" not in attributes[FIELD], "Historical omission is unknown, not false"
    assert FIELD in model.model_fields
    assert model.model_fields[FIELD].default is None


def _assert_declaration(record: ArmRun | ComparisonRun, expected: bool) -> None:
    path = Path(record.manifest_path)
    raw = json.loads(path.read_bytes())
    assert FIELD in raw, "New execution records must declare the actual option decision"
    assert raw[FIELD] is expected
    assert record.model_dump(mode="json")[FIELD] is expected
    # The new provenance-root declaration is only a boolean, never a root locator.
    assert {key: value for key, value in raw.items() if key.startswith("provenance_root")} == {
        FIELD: expected
    }
    historical = dict(raw)
    del historical[FIELD]
    restored = type(record).model_validate(historical)
    assert restored.model_dump(mode="json")[FIELD] is None
    assert restored.model_dump(mode="json", exclude={FIELD}) == record.model_dump(
        mode="json", exclude={FIELD}
    )


@pytest.mark.parametrize(
    ("arm", "supplied"), [("noel", False), ("noel", True), ("jed", False), ("ilyas", False)]
)
def test_real_arm_execution_records_supplied_or_omitted_root(
    harness: Harness, arm: str, supplied: bool
) -> None:
    root = harness.output.parent if supplied else None
    record = execute_arm(
        arm,
        harness.manifest,
        harness.output,
        _copy_plan,
        repository=harness.repository,
        figures=False,
        provenance_root=root,
    )
    assert verify_arm(Path(record.manifest_path)) == record
    _assert_declaration(record, supplied)


def test_real_comparison_omission_and_historical_read(stage: dict[str, Path]) -> None:
    record = verify_comparison_run(stage["comparison"])
    _assert_declaration(record, False)


def test_real_comparison_records_explicit_root_even_with_record_relative_basis(
    stage: dict[str, Path],
) -> None:
    root = stage["comparison"].parent
    record = execute_comparison_run(
        stage["core"],
        stage["noel"],
        stage["jed"],
        stage["ilyas"],
        root / "comparison-explicit-root",
        arm1_pooled_stats=(
            stage["comparison"]
            / "arm1_pooled_stats/trait_profile_paired_contrast_trait_stats_pooled.tsv"
        ),
        provenance_root=root,
    )
    assert verify_comparison_run(Path(record.manifest_path)) == record
    _assert_declaration(record, True)
    provenance = json.loads(
        (
            Path(record.manifest_path).parent
            / "integrated_evidence/integrated_evidence_provenance.json"
        ).read_bytes()
    )
    assert provenance["path_basis"]["kind"] == "root-relative-v1"
