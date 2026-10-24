"""Bind actual #828-produced AnalysisBasis records to executed core fixtures (#853/#869)."""

from pathlib import Path
from typing import Literal

import pytest

from mtp.datamodel.workflow_pydantic import WorkflowRun
from mtp.multiplicity.sources import bind_analysis_basis
from mtp.workflow.comparison import core_analysis_bases
from mtp.workflow.consumers import open_core_bundle
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import execute_workflow
from test_workflow_strain_rollup import FIXTURE, ROOT, SPEC, STRAIN_EDGES, STRAIN_NODES


@pytest.fixture(scope="module")
def actual_basis_runs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, WorkflowRun]:
    root = tmp_path_factory.mktemp("multiplicity-basis")
    edges, nodes = root / "edges.tsv", root / "nodes.tsv"
    edges.write_text(
        (FIXTURE / "merged-kg_edges.tsv").read_text()
        + "".join("\t".join([*row, "fixture:strain"]) + "\n" for row in STRAIN_EDGES)
    )
    nodes.write_text(
        (FIXTURE / "merged-kg_nodes.tsv").read_text()
        + "".join("\t".join(row) + "\n" for row in STRAIN_NODES)
    )
    return {
        scope: execute_workflow(
            SPEC,
            {"raw_dump": FIXTURE / "full_dump.csv", "kg_edges": edges, "kg_nodes": nodes},
            root / scope,
            OPERATIONS,
            parameter_overrides={"source_mode": "fixture", "strain_rollup": scope},
            validator=validate_core_product,
            repository_root=ROOT,
        )
        for scope in ("off", "phenotype", "full")
    }


@pytest.mark.parametrize("scope", ["off", "phenotype"])
def test_actual_producer_bases_bind_to_verified_core(
    actual_basis_runs: dict[str, WorkflowRun], scope: str
) -> None:
    run = actual_basis_runs[scope]
    bases = {basis.profile_unit: basis for basis in core_analysis_bases(run)}
    for product in ("count_profiles", "ora_profiles"):
        bundle = open_core_bundle(Path(run.manifest_path), product)
        basis, features = bind_analysis_basis(bundle, bases[product])
        assert basis == bases[product]
        assert basis.retained_feature_count == len(features)
        assert basis.rollup_policy == f"ncbi_taxonomy_closure_rollup_{scope}"
        assert all("mediadive.medium" not in feature for feature in features)
    off = core_analysis_bases(actual_basis_runs["off"])[0]
    phenotype = core_analysis_bases(actual_basis_runs["phenotype"])[0]
    assert phenotype.retained_feature_count > off.retained_feature_count
    assert phenotype.feature_filter_id != off.feature_filter_id


def test_actual_full_scope_is_explicitly_unsupported(
    actual_basis_runs: dict[str, WorkflowRun],
) -> None:
    run = actual_basis_runs["full"]
    basis = core_analysis_bases(run)[0]
    bundle = open_core_bundle(Path(run.manifest_path))
    with pytest.raises(ValueError, match="Unsupported strain rollup scope"):
        bind_analysis_basis(bundle, basis)


@pytest.mark.parametrize("product", ["count_profiles", "ora_profiles"])
@pytest.mark.parametrize(
    "field",
    ["aggregation_policy", "denominator_policy", "background_inventory_id", "weighting_policy"],
)
def test_actual_basis_rejects_changed_producer_semantics(
    actual_basis_runs: dict[str, WorkflowRun],
    product: Literal["count_profiles", "ora_profiles"],
    field: str,
) -> None:
    run = actual_basis_runs["phenotype"]
    basis = next(basis for basis in core_analysis_bases(run) if basis.profile_unit == product)
    altered = basis.model_copy(update={field: "different-scientific-basis"})
    bundle = open_core_bundle(Path(run.manifest_path), product)
    with pytest.raises(ValueError, match="producer contract"):
        bind_analysis_basis(bundle, altered)
