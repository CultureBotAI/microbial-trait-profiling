"""Revalidate persisted scientific products with their effective execution policies."""

from pathlib import Path

from mtp.datamodel.workflow_pydantic import WorkflowRun
from mtp.workflow.operations import validate_core_product
from mtp.workflow.runtime import artifact_path, validate_run


def verify_core_run(manifest: Path) -> WorkflowRun:
    """Check actual bytes and then scientific relationships across the product bundle."""
    run = validate_run(manifest)
    artifacts = {artifact.product_id: artifact for artifact in run.artifacts or []}
    paths = {key: artifact_path(run, key) for key in artifacts}
    for product in run.definition.products:
        effective = product.model_copy(update={"cohort": artifacts[product.id].cohort})
        validate_core_product(effective, paths[product.id], paths)
    return run
