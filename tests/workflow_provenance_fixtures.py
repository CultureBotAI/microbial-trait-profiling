"""Typed historical core records for bounded provenance-routing tests (#1359).

These tests replace the expensive complete-core verification boundary. Preserve
the committed historical producer declarations so policy resolution itself stays
real; substituting tiny routing artifacts does not make a scientifically verified
core bundle. Callers retain their mutable product-to-path mapping separately.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from mtp.datamodel.workflow_pydantic import ArtifactInstance, FileInventoryEntry, WorkflowRun
from mtp.paths import file_sha256
from mtp.workflow.trait_policy import resolve_trait_policy

HISTORICAL_CORE = (
    Path(__file__).parent / "fixtures/trait_key_policy/historical_raw_core/core/run.json"
)


def routing_core(
    manifest: Path, products: Mapping[str, Path], *, run_id: str = "synthetic-core"
) -> WorkflowRun:
    """Retain an actual historical declaration model with local routing artifacts."""
    run = WorkflowRun.model_validate_json(HISTORICAL_CORE.read_bytes())
    templates = {artifact.product_id: artifact for artifact in run.artifacts or []}
    artifacts: list[ArtifactInstance] = []
    for name, path in products.items():
        template = templates[name]
        sha256 = file_sha256(path)
        size_bytes = path.stat().st_size
        values = template.model_dump(mode="python")
        values.update(
            id=f"{run_id}:artifact:{name}",
            product_id=name,
            run_id=run_id,
            path=str(path.resolve()),
            requested_path=None,
            sha256=sha256,
            size_bytes=size_bytes,
            inventory=[FileInventoryEntry(path=path.name, sha256=sha256, size_bytes=size_bytes)],
            validations=[],
        )
        artifacts.append(ArtifactInstance.model_validate(values))
    values = run.model_dump(mode="python")
    values.update(run_id=run_id, manifest_path=str(manifest.resolve()), artifacts=artifacts)
    checked = WorkflowRun.model_validate(values)
    assert checked.trait_key_policy is None
    assert resolve_trait_policy(checked).policy == "raw_object_v1"
    return checked
