"""Tests verifying mtp package imports and transitional bsdbng compatibility shims (#423)."""

from __future__ import annotations

import ast
import importlib
import pkgutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest

import bsdbng
import mtp
from mtp.paths import REPO_ROOT


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
        "coupling",
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
    _assert_mirror_exports(mtp_mod, bsdbng_mod, submodule)
    # import_module resolves sys.modules replacement, so the imported object can
    # be mtp itself. Check the source's explicit export declarations as well.
    location = REPO_ROOT / "src" / "bsdbng" / submodule.replace(".", "/")
    source_path = location.with_suffix(".py")
    if not source_path.is_file():
        source_path = location / "__init__.py"
    _assert_declared_shim_imports(source_path.read_text(encoding="utf-8"), submodule)


def _assert_mirror_exports(mtp_mod: ModuleType, bsdbng_mod: ModuleType, submodule: str) -> None:
    """Canonical exports define the contract, even when the shim omits __all__."""
    exports = getattr(mtp_mod, "__all__", None)
    if exports is None:
        exports = [name for name in dir(mtp_mod) if not name.startswith("_")]
    for name in exports:
        assert hasattr(mtp_mod, name), f"mtp.{submodule} missing export {name}"
        assert hasattr(bsdbng_mod, name), f"bsdbng.{submodule} missing export {name}"
        assert getattr(mtp_mod, name) is getattr(bsdbng_mod, name), (
            f"Export {name} identity mismatch in {submodule}"
        )


def _assert_declared_shim_imports(source: str, submodule: str) -> None:
    """An alias or __getattr__ fallback must not hide a broken explicit import list."""
    tree = ast.parse(source)
    declarations = [
        node.value
        for node in tree.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "__all__" for target in node.targets)
    ]
    # Pure forwarding shims have no explicit export list. Their actual public
    # surface is still checked above against the canonical module.
    if not declarations:
        return
    assert len(declarations) == 1, f"bsdbng.{submodule} repeats its export declaration"
    exported = ast.literal_eval(declarations[0])
    assert isinstance(exported, (list, tuple))
    assert all(isinstance(name, str) for name in exported)
    assert len(set(exported)) == len(exported), f"Duplicate shim exports in {submodule}"
    imported = {
        alias.asname or alias.name
        for node in tree.body
        if isinstance(node, ast.ImportFrom)
        and node.module is not None
        and node.module.startswith("mtp")
        for alias in node.names
    }
    missing = sorted(set(exported) - imported)
    assert not missing, f"bsdbng.{submodule} declares exports without imports: {missing}"


def test_mirror_rejects_missing_export_without_all() -> None:
    """Mutation of a real alias's public surface must not evade the mirror check (#812)."""
    target = importlib.import_module("mtp.guild_maps")
    broken = ModuleType("bsdbng.guild_maps")
    broken.__dict__.update(vars(target))
    broken.__dict__.pop("__all__", None)
    del broken.__dict__["CELL_NO_SIGNAL"]
    with pytest.raises(AssertionError, match="CELL_NO_SIGNAL"):
        _assert_mirror_exports(target, broken, "guild_maps")


def test_mirror_uses_canonical_all_instead_of_a_truncated_shim_all() -> None:
    target = ModuleType("mtp.fixture")
    target.__dict__.update(__all__=["first", "second"], first=object(), second=object())
    broken = ModuleType("bsdbng.fixture")
    broken.__dict__.update(__all__=["first"], first=target.first)
    with pytest.raises(AssertionError, match="second"):
        _assert_mirror_exports(target, broken, "fixture")


def test_mirror_rejects_wrong_export_identity() -> None:
    target = ModuleType("mtp.fixture")
    target.__dict__["constant"] = object()
    broken = ModuleType("bsdbng.fixture")
    broken.__dict__["constant"] = object()
    with pytest.raises(AssertionError, match="identity mismatch"):
        _assert_mirror_exports(target, broken, "fixture")


def test_removing_a_real_shim_import_is_detected_despite_module_alias() -> None:
    path = REPO_ROOT / "src/bsdbng/guild_maps.py"
    source = path.read_text(encoding="utf-8")
    _assert_declared_shim_imports(source, "guild_maps")
    mutation = source.replace("    CELL_NO_SIGNAL,\n", "", 1)
    assert mutation != source
    with pytest.raises(AssertionError, match="CELL_NO_SIGNAL"):
        _assert_declared_shim_imports(mutation, "guild_maps")


def _check_isolated_ora_imports(first_namespace: str) -> None:
    """Exercise the package search path before another test can prepopulate aliases."""
    namespaces = [first_namespace, "mtp" if first_namespace == "bsdbng" else "bsdbng"]
    for submodule in ["adapters.ilyas_ora", "harness", "metrics", "protocol", "scenarios"]:
        modules = [
            importlib.import_module(f"{namespace}.ora_calibration.{submodule}")
            for namespace in namespaces
        ]
        assert modules[0] is modules[1], f"Duplicate implementation of {submodule}"
        assert modules[0].__name__ == f"mtp.ora_calibration.{submodule}"


@pytest.mark.parametrize("first_namespace", ["bsdbng", "mtp"])
def test_ora_children_keep_canonical_identity_in_both_import_orders(first_namespace: str) -> None:
    result = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), first_namespace],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


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


if __name__ == "__main__":
    _check_isolated_ora_imports(sys.argv[1])
