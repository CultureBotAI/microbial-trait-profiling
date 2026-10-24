"""Tests verifying mtp package imports and transitional bsdbng compatibility shims (#423)."""

from __future__ import annotations

import ast
import importlib
import pkgutil
import runpy
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest

import bsdbng
import mtp
from mtp.paths import REPO_ROOT

_EXPORT_AUDIT = runpy.run_path(str(REPO_ROOT / "scripts/adhoc/issue937_audit_exports.py"))
cross_module_uses = _EXPORT_AUDIT["cross_module_uses"]
export_omissions = _EXPORT_AUDIT["export_omissions"]
public_definitions = _EXPORT_AUDIT["public_definitions"]
source_inventory = _EXPORT_AUDIT["source_inventory"]
declared_exports = _EXPORT_AUDIT["declared_exports"]


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


def _shim_sources(root: Path = REPO_ROOT / "src/bsdbng") -> dict[str, Path]:
    """Discover source shims, including packages; generated models belong to #886."""
    paths = {}
    for path in sorted(root.rglob("*.py")):
        if path.name.endswith("_pydantic.py"):
            continue
        relative = path.relative_to(root).with_suffix("")
        parts = relative.parts[:-1] if relative.name == "__init__" else relative.parts
        paths[".".join(parts)] = path
    return paths


@pytest.mark.parametrize("submodule", _shim_sources())
def test_bsdbng_shim_mirrors_mtp_exports(submodule: str) -> None:
    """The bsdbng compatibility shim must re-export all symbols defined in mtp."""
    suffix = f".{submodule}" if submodule else ""
    mtp_mod = importlib.import_module(f"mtp{suffix}")
    bsdbng_mod = importlib.import_module(f"bsdbng{suffix}")
    _assert_mirror_exports(mtp_mod, bsdbng_mod, submodule)
    # import_module resolves sys.modules replacement, so the imported object can
    # be mtp itself. Check the source's explicit export declarations as well.
    source_path = _shim_sources()[submodule]
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
    if mtp_mod is not bsdbng_mod:
        _assert_star_import_surface(mtp_mod, bsdbng_mod)


def _canonical_declarations(module: ModuleType) -> set[str]:
    """Explicit __all__, or public definitions/assignments, excluding incidental imports.

    Runtime identity checks still cover all public runtime names. This separate
    source contract checks handwritten shim lists against what the implementation
    actually defines, so sys.modules aliases cannot make it vacuous.
    """
    explicit = getattr(module, "__all__", None)
    if explicit is not None:
        return set(explicit)
    assert module.__file__ is not None
    tree = ast.parse(Path(module.__file__).read_text(encoding="utf-8"))
    return set(public_definitions(tree))


def _assert_declared_shim_imports(source: str, submodule: str) -> None:
    """An alias or __getattr__ fallback must not hide a broken explicit import list."""
    tree = ast.parse(source)
    exported = declared_exports(tree, f"bsdbng.{submodule}")
    # Pure forwarding shims have no explicit export list. Their actual public
    # surface is still checked above against the canonical module.
    if exported is None:
        return
    suffix = f".{submodule}" if submodule else ""
    canonical = importlib.import_module(f"mtp{suffix}")
    absent = sorted(_canonical_declarations(canonical) - set(exported))
    assert not absent, f"bsdbng.{submodule} omits canonical declarations: {absent}"
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


def test_shared_ledger_argument_parser_is_an_explicit_public_export() -> None:
    """The arm verifier consumes this public parser; both export lists must include it."""
    target = importlib.import_module("mtp.taxon_contributions")
    assert "argument_parser" in target.__all__
    source = _shim_sources()["taxon_contributions"].read_text(encoding="utf-8")
    _assert_declared_shim_imports(source, "taxon_contributions")


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


def test_removing_both_shim_lists_cannot_hide_canonical_drift() -> None:
    source = _shim_sources()["guild_maps"].read_text(encoding="utf-8")
    mutation = source.replace("    CELL_NO_SIGNAL,\n", "").replace('    "CELL_NO_SIGNAL",\n', "")
    assert mutation != source
    with pytest.raises(AssertionError, match=r"omits canonical.*CELL_NO_SIGNAL"):
        _assert_declared_shim_imports(mutation, "guild_maps")


def _assert_star_import_surface(canonical: ModuleType, legacy: ModuleType) -> None:
    """Import-star uses __all__, otherwise public module-dictionary names, never dir()."""

    def names(module: ModuleType) -> set[str]:
        explicit = getattr(module, "__all__", None)
        return (
            set(explicit)
            if explicit is not None
            else {n for n in vars(module) if not n.startswith("_")}
        )

    expected = names(canonical)
    actual = names(legacy)
    assert actual == expected, f"Star-import export drift: {expected ^ actual}"
    for name in expected:
        assert getattr(legacy, name) is getattr(canonical, name), (
            f"Star-import identity drift: {name}"
        )


def test_ora_star_import_surface_includes_power_exports() -> None:
    canonical = importlib.import_module("mtp.ora_calibration")
    legacy = importlib.import_module("bsdbng.ora_calibration")
    _assert_star_import_surface(canonical, legacy)
    mutation = ModuleType("bsdbng.ora_calibration")
    mutation.__dict__.update(vars(legacy))
    mutation.__dict__["__all__"] = [
        name for name in legacy.__all__ if name not in {"PowerMetric", "compute_empirical_power"}
    ]
    with pytest.raises(AssertionError, match=r"Star-import export drift.*PowerMetric"):
        _assert_star_import_surface(canonical, mutation)
    source = _shim_sources()["ora_calibration"].read_text(encoding="utf-8")
    for name in ("PowerMetric", "compute_empirical_power"):
        source = source.replace(f"    {name},\n", "").replace(f'    "{name}",\n', "")
    with pytest.raises(AssertionError, match=r"omits canonical.*PowerMetric"):
        _assert_declared_shim_imports(source, "ora_calibration")


def test_annotated_shim_exports_cannot_hide_missing_import_or_export() -> None:
    submodule = "ora_calibration.adapters"
    source = _shim_sources()[submodule].read_text(encoding="utf-8")
    annotated = source.replace("__all__ = [", "__all__: list[str] = [")
    assert annotated != source
    _assert_declared_shim_imports(annotated, submodule)
    for mutation in (
        annotated.replace('    "IlyasORAAdapter",\n', ""),
        annotated.replace("    IlyasORAAdapter,\n", ""),
    ):
        with pytest.raises(AssertionError, match="IlyasORAAdapter"):
            _assert_declared_shim_imports(mutation, submodule)


@pytest.mark.parametrize(
    ("declaration", "message"),
    [
        ("__all__: list[str]", "uninitialized"),
        ("__all__: list[str] = build_exports()", "nonliteral"),
        ("__all__: list[str] = 'IlyasORAAdapter'", "list or tuple"),
        ("__all__: list[str] = [1]", "must be names"),
        ("__all__: list[str] = ['x', 'x']", "Duplicate"),
        ("__all__: list[str] = []\n__all__ = []", "repeats"),
    ],
)
def test_annotated_shim_declaration_rejects_malformed_values(
    declaration: str, message: str
) -> None:
    with pytest.raises(AssertionError, match=message):
        _assert_declared_shim_imports(declaration, "ora_calibration.adapters")


def test_non_alias_star_omission_is_detected_despite_getattr(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Mutate every discovered non-alias shim, and demonstrate actual import-star loss."""
    importer = tmp_path / "star_import.py"
    exercised = set()
    for submodule in _shim_sources():
        suffix = f".{submodule}" if submodule else ""
        canonical = importlib.import_module(f"mtp{suffix}")
        legacy = importlib.import_module(f"bsdbng{suffix}")
        if canonical is legacy:
            continue
        importer.write_text(f"from bsdbng{suffix} import *\n", encoding="utf-8")
        imported = runpy.run_path(str(importer))
        exports = getattr(canonical, "__all__", None)
        if exports is None:
            exports = [name for name in vars(canonical) if not name.startswith("_")]
        assert exports, f"Non-alias shim {submodule} needs a nonempty control"
        name = sorted(exports)[0]
        assert imported[name] is getattr(canonical, name)
        with monkeypatch.context() as patch:
            patch.setattr(legacy, "__all__", [n for n in exports if n != name], raising=False)
            assert getattr(legacy, name) is getattr(canonical, name)
            assert name not in runpy.run_path(str(importer))
            with pytest.raises(AssertionError, match=f"Star-import export drift.*{name}"):
                _assert_mirror_exports(canonical, legacy, submodule)
        exercised.add(submodule)
    assert {"ora_calibration", "ora_calibration.adapters"} <= exercised


def test_star_surface_without_all_uses_module_dictionary() -> None:
    canonical = ModuleType("mtp.fixture")
    canonical.__dict__.update(shared=object(), _private=object())
    canonical.__dict__["__dir__"] = lambda: ["shared", "invented"]
    legacy = ModuleType("bsdbng.fixture")
    legacy.__dict__["shared"] = canonical.shared
    _assert_star_import_surface(canonical, legacy)
    legacy.__dict__["accidental"] = object()
    with pytest.raises(AssertionError, match=r"Star-import export drift.*accidental"):
        _assert_star_import_surface(canonical, legacy)
    del legacy.__dict__["accidental"]
    del legacy.__dict__["shared"]
    legacy.__dict__["__getattr__"] = lambda name: getattr(canonical, name)
    assert legacy.shared is canonical.shared
    with pytest.raises(AssertionError, match=r"Star-import export drift.*shared"):
        _assert_star_import_surface(canonical, legacy)


def test_canonical_export_lists_cover_cross_module_public_use() -> None:
    sources, packages, _ = source_inventory(REPO_ROOT)
    assert not export_omissions(sources, packages)


def test_jed_constant_is_public_even_if_both_export_lists_omit_it() -> None:
    canonical = importlib.import_module("mtp.composition_calibration.adapters.jed_composition")
    assert "PROFILE_NONTRAIT_COLUMNS" in canonical.__all__
    sources, packages, _ = source_inventory(REPO_ROOT)
    module = canonical.__name__
    sources[module] = sources[module].replace('    "PROFILE_NONTRAIT_COLUMNS",\n', "")
    assert export_omissions(sources, packages)[module] == ["PROFILE_NONTRAIT_COLUMNS"]
    shim = _shim_sources()["composition_calibration.adapters.jed_composition"]
    source = shim.read_text(encoding="utf-8")
    _assert_declared_shim_imports(source, "composition_calibration.adapters.jed_composition")
    source = source.replace('    "PROFILE_NONTRAIT_COLUMNS",\n', "").replace(
        "    PROFILE_NONTRAIT_COLUMNS,\n", ""
    )
    with pytest.raises(AssertionError, match="PROFILE_NONTRAIT_COLUMNS"):
        _assert_declared_shim_imports(source, "composition_calibration.adapters.jed_composition")


@pytest.mark.parametrize(
    "consumer",
    [
        "from mtp.fixture import shared as renamed",
        "from .fixture import shared",
        "import mtp.fixture as alias\nvalue = alias.shared",
        "import mtp.fixture\nvalue = mtp.fixture.shared",
        "from mtp import fixture as alias\nvalue = alias.shared",
    ],
)
def test_cross_module_audit_recognizes_import_forms(consumer: str) -> None:
    sources = {"mtp.fixture": "shared = 1\n__all__ = []", "mtp.consumer": consumer}
    assert export_omissions(sources, set()) == {"mtp.fixture": ["shared"]}
    assert cross_module_uses(sources, set())[("mtp.fixture", "shared")] == {"mtp.consumer"}
    sources["mtp.fixture"] = "shared = 1\n__all__: list[str] = ['shared']"
    assert export_omissions(sources, set()) == {}


def test_cross_module_audit_distinguishes_canonical_reexports_from_dependencies() -> None:
    sources = {
        "mtp.fixture": (
            "import json\nfrom typing import Any\nfrom .origin import shared as public\n"
            "_private = 1\n__all__ = []"
        ),
        "mtp.origin": "shared = 1",
        "mtp.consumer": "from mtp.fixture import json, Any, _private, public",
    }
    assert export_omissions(sources, set()) == {"mtp.fixture": ["public"]}
    sources["mtp.fixture"] += "\n__all__ = ['public']"
    with pytest.raises(AssertionError, match="repeats"):
        export_omissions(sources, set())


def test_canonical_module_reexports_are_audited_without_third_party_imports() -> None:
    sources = {
        "mtp.origin": "shared = 1",
        "mtp.facade": "import mtp.origin as public\nimport json as incidental\n__all__ = []",
        "mtp.consumer": "from mtp.facade import public, incidental",
    }
    assert export_omissions(sources, set()) == {"mtp.facade": ["public"]}
    sources["mtp.facade"] = sources["mtp.facade"].replace("__all__ = []", "__all__ = ['public']")
    assert export_omissions(sources, set()) == {}


def test_cross_module_aliases_are_scoped_to_their_functions() -> None:
    sources = {
        "mtp.first": "first = 1\n__all__ = []",
        "mtp.second": "second = 2\n__all__ = []",
        "mtp.consumer": (
            "def one():\n    import mtp.first as module\n    return module.first\n"
            "def two():\n    import mtp.second as module\n    return module.second\n"
        ),
    }
    assert export_omissions(sources, set()) == {
        "mtp.first": ["first"],
        "mtp.second": ["second"],
    }


def test_public_definitions_include_only_bound_assignment_names() -> None:
    tree = ast.parse(
        "scalar = 1\nfirst, [second, third] = values\nannotated: str = ''\n"
        "cache[key] = 1\nmodule.attribute = 1\n_private = 1\n"
    )
    assert public_definitions(tree) == {"scalar", "first", "second", "third", "annotated"}


def test_cross_module_audit_handles_package_relative_import_and_parameter_shadow() -> None:
    sources = {
        "mtp.package": "from .origin import public as renamed",
        "mtp.package.origin": "public = 1\n__all__ = []",
        "mtp.consumer": (
            "import mtp.package.origin as module\n"
            "def unrelated(module):\n    return module.incidental\n"
        ),
    }
    uses = cross_module_uses(sources, {"mtp.package"})
    assert uses == {("mtp.package.origin", "public"): {"mtp.package"}}
    assert export_omissions(sources, {"mtp.package"}) == {"mtp.package.origin": ["public"]}


def test_class_local_alias_does_not_replace_method_global_binding() -> None:
    sources = {
        "mtp.first": "first = 1\n__all__ = []",
        "mtp.second": "second = 2\n__all__ = []",
        "mtp.consumer": (
            "import mtp.first as module\n"
            "class Example:\n    import mtp.second as module\n    value = module.second\n"
            "    def method(self):\n        return module.first\n"
        ),
    }
    assert export_omissions(sources, set()) == {
        "mtp.first": ["first"],
        "mtp.second": ["second"],
    }


@pytest.mark.parametrize(
    "class_body",
    [
        "    callback = lambda value=module.second: module.first\n",
        "    values = [module.first for item in module.second]\n",
        "    class Nested(module.second):\n        value = module.first\n",
    ],
)
def test_class_expression_bodies_use_non_class_lexical_parent(class_body: str) -> None:
    sources = {
        "mtp.first": "first = 1\n__all__ = []",
        "mtp.second": "second = 2\n__all__ = []",
        "mtp.consumer": (
            "import mtp.first as module\n"
            "class Example:\n    import mtp.second as module\n" + class_body
        ),
    }
    assert cross_module_uses(sources, set()) == {
        ("mtp.first", "first"): {"mtp.consumer"},
        ("mtp.second", "second"): {"mtp.consumer"},
    }
    assert export_omissions(sources, set()) == {
        "mtp.first": ["first"],
        "mtp.second": ["second"],
    }


@pytest.mark.parametrize(
    ("class_body", "uses_class_member"),
    [
        ("    def f(self, x=lambda: module.first):\n        pass\n", False),
        ("    def f(self, x=(module.second, lambda: module.first)):\n        pass\n", True),
        ("    class Nested(factory(module.second, lambda: module.first)):\n        pass\n", True),
        ("    callback = lambda x=(module.second, lambda: module.first): x\n", True),
        ("    values = [item for item in (module.second, lambda: module.first)]\n", True),
    ],
)
def test_class_defining_expressions_preserve_nested_body_context(
    class_body: str, uses_class_member: bool
) -> None:
    sources = {
        "mtp.first": "first = 1\n__all__ = []",
        "mtp.second": "second = 2\n__all__ = []",
        "mtp.consumer": (
            "import mtp.first as module\n"
            "class Example:\n    import mtp.second as module\n" + class_body
        ),
    }
    expected = {("mtp.first", "first"): {"mtp.consumer"}}
    if uses_class_member:
        expected[("mtp.second", "second")] = {"mtp.consumer"}
    assert cross_module_uses(sources, set()) == expected


@pytest.mark.parametrize(
    "consumer",
    [
        "import mtp.fixture as module\ndef f(x=module.public):\n    pass",
        "import mtp.fixture as module\n@module.public\ndef f():\n    pass",
        "import mtp.fixture as module\nclass C(module.public):\n    pass",
        "import mtp.fixture as module\nvalue = lambda: module.public",
        "import mtp.fixture as module\nvalue = [module.public for item in values]",
        "import mtp.fixture as module\nvalue = [item for module in module.public]",
        "import mtp.fixture as module\ndef f(x: module.public):\n    pass",
        "import mtp.fixture as module\ndef f() -> module.public:\n    pass",
        "import mtp.fixture as module\ndef f(x=lambda module: module.public):\n"
        "    return module.public",
    ],
)
def test_qualified_uses_in_deferred_or_nested_expressions(consumer: str) -> None:
    sources = {"mtp.fixture": "public = 1\n__all__ = []", "mtp.consumer": consumer}
    assert export_omissions(sources, set()) == {"mtp.fixture": ["public"]}


@pytest.mark.parametrize(
    "expression",
    ["lambda module: module.public", "[module.public for module in values]"],
)
def test_nested_expression_bindings_do_not_inherit_shadowed_module_alias(expression: str) -> None:
    sources = {
        "mtp.fixture": "public = 1\n__all__ = []",
        "mtp.consumer": f"import mtp.fixture as module\nvalue = {expression}",
    }
    assert export_omissions(sources, set()) == {}


def test_shim_discovery_includes_new_modules_without_a_handwritten_list(tmp_path: Path) -> None:
    assert "screen_calibration" in _shim_sources()
    (tmp_path / "new_package").mkdir()
    for filename in ("__init__.py", "new_package/__init__.py", "new_package/new_module.py"):
        (tmp_path / filename).write_text('"""A newly added shim."""\n')
    assert set(_shim_sources(tmp_path)) == {"", "new_package", "new_package.new_module"}


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
