"""Revalidate persisted scientific products with their effective execution policies."""

from pathlib import Path

from mtp.datamodel.workflow_pydantic import WorkflowRun
from mtp.workflow.locations import LocationResolver
from mtp.workflow.operations import validate_core_product
from mtp.workflow.runtime import artifact_path, validate_run
from mtp.workflow.trait_policy import resolve_trait_policy


def verify_core_run(manifest: Path, *, resolver: LocationResolver | None = None) -> WorkflowRun:
    """Check actual bytes and then scientific relationships across the product bundle."""
    run = validate_run(manifest, resolver=resolver)
    trait_policy = resolve_trait_policy(run)
    artifacts = {artifact.product_id: artifact for artifact in run.artifacts or []}
    paths = {key: artifact_path(run, key, resolver) for key in artifacts}
    recorded_sources = (
        {
            product.id: Path(artifacts[product.id].path)
            for product in run.definition.products
            if product.external
        }
        if resolver is not None
        else None
    )
    parameters = {item.name: item.value for item in run.parameters or []}
    for product in run.definition.products:
        effective = product.model_copy(update={"cohort": artifacts[product.id].cohort})
        validate_core_product(
            effective,
            paths[product.id],
            paths,
            recorded_sources=recorded_sources,
            parameters=parameters,
            trait_key_policy=trait_policy.policy,
            include_policy_metadata=trait_policy.include_policy_metadata,
        )
    return run
