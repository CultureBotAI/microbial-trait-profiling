"""Build the tracked subset graph while preserving the full-source ORA background.

Input: workflows/bugsigdb_core.yaml. Output: workflows/bugsigdb_core_subset.yaml.
Serves #658 / #660; 2026-10-19. No analysis data are read or modified.
"""

from copy import deepcopy
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
EXPORTS = (
    "study_table",
    "experiments",
    "experiment_ontology_terms",
    "signatures",
    "reported_taxa",
    "taxon_profiles_wide",
    "taxon_profiles_long",
)


class IndentedSafeDumper(yaml.SafeDumper):
    """Keep generated sequence indentation consistent with repository YAML."""

    def increase_indent(self, flow: bool = False, indentless: bool = False) -> None:
        super().increase_indent(flow=flow, indentless=indentless if flow else False)


def main() -> None:
    definition = yaml.safe_load((ROOT / "workflows/bugsigdb_core.yaml").read_text())
    definition["id"] = "bugsigdb-core-subset"
    definition["parameters"].append({"name": "contrast_population", "default": "all_paired"})
    definition["description"] = (
        "Verified contrast-population selection and disease orientation before canonical "
        "exports and profiles. Original ingested studies and exports remain preserved. "
        "The trait matrix, vocabulary and ORA background use the full source population."
    )
    products = {product["id"]: product for product in definition["products"]}
    source_studies = deepcopy(products["studies"])
    source_studies.update(id="source_studies", filename="source/studies")
    definition["products"].append(source_studies)
    products["studies"].update(
        kind="contrast_study_collection",
        cohort="all_paired",
        cohort_parameter="contrast_population",
        coverage_policy="contrast_population_v1",
        description=(
            "Exact selected projection of source_studies. Group metadata remains original; "
            "signature directions follow the recorded WorkflowRun direction_basis."
        ),
    )
    for name in EXPORTS:
        original = deepcopy(products[name])
        original.update(id=f"source_{name}", filename=f"source/{original['filename']}")
        definition["products"].append(original)
    definition["products"].append(
        {
            "id": "contrast_selection",
            "kind": "contrast_selection",
            "storage_kind": "file",
            "format": "json",
            "filename": "provenance/contrast_selection.json",
            "unit": "experiment_selection_decision",
            "cohort": "all_paired",
            "cohort_parameter": "contrast_population",
            "coverage_policy": "contrast_population_v1",
        }
    )
    steps = {step["id"]: step for step in definition["steps"]}
    steps["ingest"].update(operation="ingest_source", outputs=["source_studies", "ingest_log"])
    steps["taxa_matrix"].update(
        operation="source_taxa_matrix", inputs=["source_reported_taxa", "kg_edges", "kg_nodes"]
    )
    definition["steps"][2:2] = [
        {
            "id": "source_export_tables",
            "operation": "export_source_tables",
            "inputs": ["source_studies"],
            "outputs": [f"source_{name}" for name in EXPORTS],
        },
        {
            "id": "select_contrasts",
            "operation": "select_contrasts",
            "inputs": ["source_studies", "raw_dump"],
            "outputs": ["studies", "contrast_selection"],
            "parameters": ["contrast_population"],
        },
    ]
    destination = ROOT / "workflows/bugsigdb_core_subset.yaml"
    destination.write_text(
        yaml.dump(definition, Dumper=IndentedSafeDumper, sort_keys=False, width=100)
    )


if __name__ == "__main__":
    main()
