"""Fresh-process import identity contracts for the legacy datamodel package (#886)."""

from __future__ import annotations

import importlib
import json
import os
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest

MODELS = {
    "mtp_pydantic": "Study",
    "core_products_pydantic": "CountProfileRow",
    "workflow_pydantic": "ProductDefinition",
    "multiple_testing_pydantic": "HypothesisFamily",
}
ALIASES = {**{name: name for name in MODELS}, "bsdbng_pydantic": "mtp_pydantic"}
FIRST_IMPORTS = {
    "legacy-child-mtp": "bsdbng.datamodel.mtp_pydantic",
    "legacy-child-bsdbng": "bsdbng.datamodel.bsdbng_pydantic",
    "canonical-child": "mtp.datamodel.mtp_pydantic",
    "legacy-package": "bsdbng.datamodel",
    "canonical-package": "mtp.datamodel",
}


def _available(name: str) -> ModuleType:
    try:
        return importlib.import_module(name)
    except ImportError as error:
        raise AssertionError(
            f"Datamodel import unavailable: {name} ({type(error).__name__}: {error})"
        ) from error


def _probe(order: str) -> None:
    source = Path(__file__).resolve().parents[1] / "src"
    first = _available(FIRST_IMPORTS[order])
    canonical = _available("mtp.datamodel")
    modules = {name: _available(f"mtp.datamodel.{name}") for name in MODELS}
    classes = {name: getattr(module, MODELS[name]) for name, module in modules.items()}
    study = canonical.Study
    legacy = _available("bsdbng.datamodel")
    for alias, name in ALIASES.items():
        module = _available(f"bsdbng.datamodel.{alias}")
        assert module is modules[name], f"Duplicate legacy datamodel module: {alias}"
        assert sys.modules[f"bsdbng.datamodel.{alias}"] is modules[name]
        assert vars(legacy).get(alias) is modules[name], f"Missing legacy parent attribute: {alias}"
        assert getattr(module, MODELS[name]) is classes[name], f"Duplicate model class: {alias}"
        for child, saved in modules.items():
            assert sys.modules[f"mtp.datamodel.{child}"] is saved, (
                f"Canonical module rebound: {child}"
            )
            assert vars(canonical).get(child) is saved, (
                f"Canonical parent attribute rebound: {child}"
            )
            assert getattr(saved, MODELS[child]) is classes[child], (
                f"Canonical class rebound: {child}"
            )
        assert canonical.Study is study, "Canonical package Study rebound"
    assert legacy is not canonical, "Legacy datamodel package must retain its own search path"
    assert sys.modules[FIRST_IMPORTS[order]] is first, "First imported object was rebound"
    assert vars(_available("bsdbng")).get("datamodel") is legacy
    assert vars(_available("mtp")).get("datamodel") is canonical
    # Exercise actual from/dotted syntax as well as importlib's sys.modules lookup.
    import mtp.datamodel.mtp_pydantic as canonical_child
    from bsdbng.datamodel import Study as legacy_study
    from bsdbng.datamodel import bsdbng_pydantic as historical_child
    from bsdbng.datamodel import mtp_pydantic as modern_child

    assert legacy_study is study is modules["mtp_pydantic"].Study
    assert historical_child is modern_child is canonical_child is modules["mtp_pydantic"]
    origins = {}
    for name, module in modules.items():
        qualified = f"mtp.datamodel.{name}"
        assert module.__name__ == qualified, f"Noncanonical module name: {name}"
        assert classes[name].__module__ == qualified, f"Noncanonical class module: {name}"
        assert module.__file__ is not None
        assert module.__spec__ is not None
        assert module.__spec__.origin is not None
        expected = source / "mtp" / "datamodel" / f"{name}.py"
        assert Path(module.__file__).resolve() == expected, f"Wrong module source: {name}"
        assert Path(module.__spec__.origin).resolve() == expected, f"Wrong spec source: {name}"
        origins[name] = {
            "file": module.__file__,
            "spec_origin": module.__spec__.origin,
            "class_module": classes[name].__module__,
        }
    print(
        json.dumps({"order": order, "origins": origins, "aliases": list(ALIASES)}, sort_keys=True)
    )


@pytest.mark.parametrize("order", FIRST_IMPORTS)
def test_legacy_datamodel_identity_in_fresh_import_orders(order: str) -> None:
    source = Path(__file__).resolve().parents[1] / "src"
    result = subprocess.run(
        [sys.executable, "-B", str(Path(__file__).resolve()), order],
        env={**os.environ, "PYTHONPATH": str(source), "PYTHONDONTWRITEBYTECODE": "1"},
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert result.returncode == 0, f"Fresh-process datamodel assertion failed:\n{result.stderr}"
    receipt = json.loads(result.stdout)
    assert receipt["order"] == order
    assert set(receipt["origins"]) == set(MODELS)
    assert receipt["aliases"] == list(ALIASES)


if __name__ == "__main__":
    _probe(sys.argv[1])
