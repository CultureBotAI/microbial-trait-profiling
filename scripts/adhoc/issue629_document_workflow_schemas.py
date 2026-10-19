"""Document LinkML contracts for issue #629 (2026-10-19).

Inputs/outputs: schema/workflow.yaml and schema/core_products.yaml, updated in
place with field titles and descriptions. Existing explicit prose is preserved.
This saved helper records the mechanical documentation expansion, not analysis.
"""

from __future__ import annotations

import re
from pathlib import Path

import yaml

DESCRIPTIONS = {
    "id": "Identifier unique within its containing document.",
    "storage_kind": "Whether the artifact is a regular file or an inventoried directory.",
    "external": "Whether this product is supplied explicitly rather than produced by a step.",
    "name": "Declared parameter identifier.",
    "required": "Whether execution requires an effective parameter value.",
    "value": "Effective resolved parameter value passed to the operation.",
    "source": "Origin of the effective value: default, declared environment variable, or override.",
    "operation": "Key resolving to an explicitly registered Python callable.",
    "inputs": "Product references whose verified artifact paths are passed to the operation.",
    "outputs": "Declared outputs; references are resolved and checked against their producers.",
    "depends_on": "Additional prerequisite step references beyond input-producing dependencies.",
    "parameters": "Declared or effective parameters scoped to this definition or execution.",
    "version": "Version of the declared workflow contract.",
    "products": "Scientific product contracts available to this workflow.",
    "steps": "Ordered or dependency-resolved workflow stages or their attempted executions.",
    "sha256": "Lowercase SHA-256 of file bytes or the canonical collection inventory.",
    "size_bytes": "Total number of regular-file bytes represented by this record.",
    "validator": "Versioned validation procedure identifier.",
    "status": "Lifecycle or validation outcome; only completed verified runs are consumable.",
    "message": "Human-readable validation result.",
    "product_id": "Reference to the embedded product contract.",
    "run_id": "Unique execution-attempt identifier, distinct from the deterministic fingerprint.",
    "producer_step": "Reference to the execution that produced this artifact; absent for inputs.",
    "inventory": "Complete sorted list of regular-file members and their byte hashes.",
    "validations": "Recorded inventory and scientific product validation outcomes.",
    "step_id": "Reference to the declared step contract.",
    "started_at": "UTC timestamp when this execution attempt began.",
    "finished_at": "UTC timestamp when this execution attempt finished or failed.",
    "input_artifacts": "Concrete artifact references actually bound as operation inputs.",
    "output_artifacts": "Concrete artifact references published by this completed operation.",
    "definition": "Embedded validated workflow definition used by this attempt.",
    "manifest_path": "Manifest location; relocation does not change execution identity.",
    "spec_path": "Source location of the workflow specification at execution time.",
    "spec_sha256": "SHA-256 of the exact source specification bytes.",
    "definition_sha256": "SHA-256 of the canonical embedded validated workflow definition.",
    "execution_fingerprint": "Hash of contracts, inputs, code, environment, and resolved values.",
    "code_sha256": "SHA-256 of the complete recorded code inventory.",
    "code_inventory": "Source and schema files, including untracked implementation code.",
    "lock_sha256": "SHA-256 of the dependency lockfile used for this attempt.",
    "git_commit": "Repository commit at execution start; code hashes also capture local changes.",
    "git_dirty": "Whether tracked or untracked repository changes existed at execution start.",
    "runtime_json": "Canonical JSON of Python, package versions, platform, and operations.",
    "artifacts": "Concrete external input and published output artifact instances.",
    "run_manifest_path": "Explicit path to the verified producing run manifest.",
    "run_manifest_sha256": "SHA-256 of the producing run manifest consumed by the analysis.",
    "consumer": "Identifier of the downstream analysis or presentation operation.",
    "product_ids": "Product contracts explicitly consumed by this analysis.",
    "artifact_ids": "Concrete verified artifacts consumed by this analysis.",
    "created_at": "UTC timestamp when this consumer record was created.",
    "source_path": "Location of the historical source whose claims are being preserved.",
    "source_sha256": "SHA-256 of the historical source bytes at import time.",
    "imported_at": "UTC timestamp of legacy import, not the unknown historical execution time.",
    "verification_status": "Unverified historical execution; byte integrity alone is insufficient.",
    "claims_json": "Original historical claims preserved as JSON without inventing executed steps.",
    "features": "Predicate-object feature definitions for the matrix vocabulary.",
    "coverage": "Per-profile accounting of reported, matched, and annotated taxa.",
    "studies": "Canonical normalized BugSigDB study records.",
    "assertions": "Canonical KG-Microbe trait assertion records.",
    "taxon_trait_rows": "Canonical taxon-to-feature matrix rows.",
    "column_id": "Stable predicate-object column identifier matching the matrix header.",
    "column_label": "Human-readable matrix feature label.",
    "predicate_id": "Relation identifier in the KG-Microbe assertion.",
    "predicate_name": "Human-readable relation name.",
    "object_id": "Trait or culture-medium object identifier.",
    "object_name": "Human-readable trait or medium name.",
    "object_category": "KG category of the assertion object.",
    "feature_class": "Whether the feature is a trait annotation or a culture-medium relation.",
    "annotation_evidence": "Direct versus descendant evidence remains unresolved by this matrix.",
    "study_id": "Study component of the unique profile key.",
    "experiment_id": "Experiment component of the unique profile key.",
    "inc_dec": "Reported increased or decreased direction component of the profile key.",
    "reported_distinct_taxa": "Distinct taxon identifiers reported in the source profile.",
    "lookup_matched_taxa": "Reported taxon identifiers present in the explicit taxon-trait lookup.",
    "lookup_unmatched_taxa": "Reported taxon identifiers absent from the taxon-trait lookup.",
    "matched_taxa_with_annotations": "Matched taxa with at least one retained feature annotation.",
    "matched_taxa_without_annotations": "Matched taxa with no retained feature annotations.",
    "coverage_status": "Explicit coverage category determined by the reconciled counts.",
    "trait_filter_id": "SHA-256 identifying the exact ordered retained feature vocabulary.",
    "lineage_resolution": "Unresolved lineage; presence does not establish direct annotation.",
}


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    for name in ("workflow", "core_products"):
        path = root / "schema" / f"{name}.yaml"
        schema = yaml.safe_load(path.read_text(encoding="utf-8"))
        descriptions = dict(DESCRIPTIONS)
        for value in schema.get("classes", {}).values():
            for field, definition in value.get("attributes", {}).items():
                if "description" in definition:
                    descriptions.setdefault(field, definition["description"])
        for section in ("classes", "enums"):
            for key, value in schema.get(section, {}).items():
                title = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", key)
                value.setdefault("title", title)
                value.setdefault("description", f"Permitted {title.lower()} values.")
                for field, definition in value.get("attributes", {}).items():
                    definition.setdefault("title", field.replace("_", " ").capitalize())
                    if "description" not in definition:
                        definition["description"] = descriptions[field]
        path.write_text(yaml.safe_dump(schema, sort_keys=False, width=100), encoding="utf-8")


if __name__ == "__main__":
    main()
