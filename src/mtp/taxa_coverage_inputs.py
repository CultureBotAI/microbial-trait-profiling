"""Bind the coverage figure's inputs and population to one verified core (#1344)."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from mtp.datamodel.workflow_pydantic import ContrastPopulation, DirectionBasis
from mtp.taxa_trait_matrix import StrainRollup
from mtp.workflow.consumers import CoreBundle, open_core_bundle
from mtp.workflow.operations import run_strain_rollup

CORE_PRODUCTS = {
    "taxa": "reported_taxa",
    "matrix": "taxa_matrix",
    "profiles": "count_profiles",
    "experiments": "experiments",
    "raw_dump": "raw_dump",
    "taxonomy_edges": "kg_edges",
    "kg_nodes": "kg_nodes",
    "source_identity": "source_identity",
}


@dataclass(frozen=True)
class CoverageInputs:
    """Authoritative paths and direction state; no raw contrast classification."""

    bundle: CoreBundle
    paths: dict[str, Path]
    strain_rollup: StrainRollup
    contrast_population: str
    direction_basis: str


def bind_coverage_inputs(
    manifest: Path,
    *,
    supplied: Mapping[str, Path | None],
    strain_rollup: str | None,
    contrast_population: str | None,
    direction_basis: str | None,
) -> CoverageInputs:
    """Verify once and reject conflicting declarations before a figure writes output."""
    bundle = open_core_bundle(manifest)
    paths = {name: bundle.path(product) for name, product in CORE_PRODUCTS.items()}
    for name, provided in supplied.items():
        if name not in paths:
            raise ValueError(f"Unknown core coverage input: {name}")
        if provided is not None and provided.resolve() != paths[name].resolve():
            raise ValueError(f"--core-run conflicts with explicit --{name.replace('_', '-')}")
    scope = run_strain_rollup(bundle.run)
    population = ContrastPopulation(bundle.run.contrast_population).value
    direction = DirectionBasis(bundle.run.direction_basis).value
    for name, declaration, recorded in (
        ("strain-rollup", strain_rollup, scope),
        ("contrast-population", contrast_population, population),
        ("direction-basis", direction_basis, direction),
    ):
        if declaration is not None and declaration != recorded:
            raise ValueError(f"--core-run conflicts with explicit --{name}: {declaration!r}")
    return CoverageInputs(bundle, paths, scope, population, direction)
