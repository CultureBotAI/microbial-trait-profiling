"""Tests verifying mtp package imports and transitional bsdbng compatibility shims (#423)."""

from __future__ import annotations

import importlib
import pkgutil

import pytest

import bsdbng
import mtp


def _get_submodule_names(package_name: str) -> list[str]:
    """Discover all submodule names within a package."""
    pkg = importlib.import_module(package_name)
    submodules: list[str] = []
    if not hasattr(pkg, "__path__"):
        return submodules
    for info in pkgutil.walk_packages(pkg.__path__, prefix=f"{package_name}."):
        # Exclude linkml-generated pydantic modules from direct attribute checks
        if "pydantic" in info.name:
            continue
        submodules.append(info.name)
    return sorted(submodules)


def test_mtp_version_matches_bsdbng() -> None:
    """The version exposed by bsdbng shim must match mtp."""
    assert mtp.__version__ == bsdbng.__version__
    assert mtp.__version__ == "0.1.0"


def test_mtp_modules_can_be_imported_directly() -> None:
    """Every submodule under mtp must be importable without error."""
    mtp_submodules = _get_submodule_names("mtp")
    assert len(mtp_submodules) >= 15
    for mod_name in mtp_submodules:
        mod = importlib.import_module(mod_name)
        assert mod is not None


@pytest.mark.parametrize(
    "submodule",
    [
        "calibration",
        "contrast_orientation",
        "contribution_inputs",
        "composition_calibration",
        "composition_calibration.adapters",
        "composition_calibration.adapters.jed_composition",
        "composition_calibration.composition",
        "composition_calibration.protocol",
        "composition_calibration.scenarios",
        "corpus_identity",
        "datamodel",
        "deck_visualizations",
        "density_adjustment",
        "download",
        "enrichment_results",
        "enrichment_stats",
        "guild_maps",
        "ingest",
        "integrated_evidence",
        "log_ratio_sensitivity",
        "native_influence_publication",
        "native_taxon_influence",
        "ontology_audit",
        "ora_calibration",
        "ora_calibration.adapters",
        "ora_calibration.adapters.ilyas_ora",
        "ora_calibration.harness",
        "ora_calibration.metrics",
        "ora_calibration.protocol",
        "ora_calibration.scenarios",
        "paths",
        "stats",
        "stratification",
        "study_trait_profiles",
        "study_tsv",
        "table1_trait_maps",
        "taxa_trait_matrix",
        "taxon_contributions",
        "taxon_influence",
        "trait_composition",
        "trait_guilds",
        "trait_screen",
        "trait_terms",
        "validate",
        "validation_evidence",
    ],
)
def test_bsdbng_shim_mirrors_mtp_exports(submodule: str) -> None:
    """The bsdbng compatibility shim must re-export all symbols defined in mtp."""
    mtp_mod = importlib.import_module(f"mtp.{submodule}")
    bsdbng_mod = importlib.import_module(f"bsdbng.{submodule}")

    if hasattr(bsdbng_mod, "__all__"):
        for name in bsdbng_mod.__all__:
            assert hasattr(mtp_mod, name), f"mtp.{submodule} missing export {name}"
            assert hasattr(bsdbng_mod, name), f"bsdbng.{submodule} missing export {name}"
            val_mtp = getattr(mtp_mod, name)
            val_bsdbng = getattr(bsdbng_mod, name)
            assert val_mtp is val_bsdbng, f"Export {name} identity mismatch in {submodule}"


def test_datamodel_classes_identical() -> None:
    """Pydantic datamodel classes from bsdbng.datamodel and mtp.datamodel must match."""
    import bsdbng.datamodel as b_dm
    import mtp.datamodel as m_dm

    for model_name in [
        "Study",
        "Signature",
        "Experiment",
        "Taxon",
        "TaxonTraitMatrixRow",
        "DirectionEnum",
    ]:
        b_cls = getattr(b_dm, model_name)
        m_cls = getattr(m_dm, model_name)
        assert b_cls is m_cls, f"Datamodel class {model_name} should be identical"
