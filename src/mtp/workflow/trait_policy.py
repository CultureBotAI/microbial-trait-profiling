"""Resolve producer-bound trait semantics without executing historical code (#1129).

Old runs have no explicit policy. Only the exact audited matrix implementation
bytes identify their implicit rule; a commit label or observed output vocabulary
cannot establish it. The caller retains full run/inventory and scientific
verification. Resolution never mutates historical manifests or shared state.
"""

from __future__ import annotations

from dataclasses import dataclass

from mtp.datamodel.workflow_pydantic import FileInventoryEntry, WorkflowRun
from mtp.trait_key_policy import TraitKeyPolicy, parse_trait_key_policy

POLICY_PARAMETER = "trait_key_policy"
MATRIX_SOURCE = "src/mtp/taxa_trait_matrix.py"

# Audited Git blobs: F 511c891 (raw IDs), then #1124 0e49e87 (union IDs).
# Match both the complete file digest and size; unfamiliar historical producers
# require an explicit compatibility review, never an inferred default.
HISTORICAL_MATRIX_POLICIES: dict[tuple[str, int], TraitKeyPolicy] = {
    ("e7208fbc26c505cb5b7ad3e35b515b9ec1bd4976307d06046e27788f2ac8402f", 26364): "raw_object_v1",
    ("db02030a8ce6c33a7d57151517962722470f7f0973d9ea49db30ca58c96fc874", 26851): (
        "facultative_union_v1"
    ),
}


@dataclass(frozen=True)
class ResolvedTraitPolicy:
    """Independent scientific key policy and producer's metadata serialization format."""

    policy: TraitKeyPolicy
    include_policy_metadata: bool


def _matrix_implementation(run: WorkflowRun) -> FileInventoryEntry:
    paths = [entry.path for entry in run.code_inventory]
    if len(paths) != len(set(paths)):
        raise ValueError("Duplicate code inventory paths cannot identify a trait policy")
    entries = [entry for entry in run.code_inventory if entry.path == MATRIX_SOURCE]
    if len(entries) != 1:
        raise ValueError("Trait policy requires one recorded matrix implementation inventory")
    return entries[0]


def _validate_producer_bindings(run: WorkflowRun, policy: TraitKeyPolicy) -> None:
    products = {
        product.id
        for product in run.definition.products
        if product.kind in {"taxa_matrix", "feature_catalog"}
    }
    producers = [step for step in run.definition.steps if products.intersection(step.outputs)]
    if not producers:
        raise ValueError("Explicit trait policy has no matrix/catalog producer binding")
    for producer in producers:
        if (producer.parameters or []).count(POLICY_PARAMETER) != 1:
            raise ValueError(f"Producer {producer.id} must bind exactly one trait policy")
        executions = [step for step in run.steps or [] if step.step_id == producer.id]
        if len(executions) != 1:
            raise ValueError(f"Trait policy requires one producer execution: {producer.id}")
        values = [
            item.value for item in executions[0].parameters or [] if item.name == POLICY_PARAMETER
        ]
        if values != [policy]:
            raise ValueError(f"Producer {producer.id} trait policy differs from the core policy")


def resolve_trait_policy(run: WorkflowRun) -> ResolvedTraitPolicy:
    """Resolve the checked core's rule, preserving its original recorded fields.

    Full execution/inventory verification remains the caller's responsibility.
    An explicit policy requires matching manifest, parameter and producer-step
    declarations. Otherwise only audited historical implementation bytes select
    an implicit policy and the historical metadata format.
    """
    implementation = _matrix_implementation(run)
    identity = (implementation.sha256, implementation.size_bytes)
    recorded = [item.value for item in run.parameters or [] if item.name == POLICY_PARAMETER]
    declarations = [
        item for item in run.definition.parameters or [] if item.name == POLICY_PARAMETER
    ]
    step_declaration = any(
        POLICY_PARAMETER in (step.parameters or []) for step in run.definition.steps
    )
    execution_declaration = any(
        item.name == POLICY_PARAMETER for step in run.steps or [] for item in step.parameters or []
    )
    if run.trait_key_policy is None:
        if recorded or declarations or step_declaration or execution_declaration:
            raise ValueError("Partial explicit trait policy cannot be treated as historical")
        if identity not in HISTORICAL_MATRIX_POLICIES:
            raise ValueError("Unknown historical matrix implementation for trait policy")
        return ResolvedTraitPolicy(HISTORICAL_MATRIX_POLICIES[identity], False)

    policy = parse_trait_key_policy(run.trait_key_policy)
    if recorded != [policy] or len(declarations) != 1:
        raise ValueError("Core trait policy differs from its unique declared resolved parameter")
    if identity in HISTORICAL_MATRIX_POLICIES:
        raise ValueError(
            "Historical implicit-policy implementation cannot claim an explicit policy"
        )
    _validate_producer_bindings(run, policy)
    return ResolvedTraitPolicy(policy, True)
