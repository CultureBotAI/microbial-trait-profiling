"""The #961 mutation probe may write only its exclusive, bounded source copy."""

from __future__ import annotations

import hashlib
import importlib.util
from collections.abc import Callable
from pathlib import Path
from types import ModuleType
from typing import Any, cast

import pytest
from git import Repo

from mtp.paths import REPO_ROOT

SOURCE = "src/mtp/coassertion_audit.py"
FILES = {
    "src/mtp/__init__.py": b"\n",
    SOURCE: b"TARGET = True\n",
    "src/mtp/taxa_trait_matrix.py": b"PARENT = True\n",
    "schema/workflow.yaml": b"id: https://example.org/fixture\n",
    "tests/conftest.py": b"\n",
    "tests/test_coassertion_audit.py": b"def test_fixture():\n    assert True\n",
    "methods/coassertion_audit.md": b"Fixture protocol.\n",
    "methods/calibration_protocol.json": b"{}\n",
    "pyproject.toml": b"[tool.pytest.ini_options]\naddopts = '-q'\n",
    "uv.lock": b"version = 1\n",
}


@pytest.fixture
def verifier() -> ModuleType:
    path = REPO_ROOT / "scripts/adhoc/issue948_coassertion_mutations.py"
    spec = importlib.util.spec_from_file_location("coassertion_isolation_under_test", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def api(verifier: ModuleType, name: str) -> Callable[..., Any]:
    function = getattr(verifier, name, None)
    assert callable(function), f"The mutation runner must expose its {name} boundary"
    return cast(Callable[..., Any], function)


@pytest.fixture
def repository(tmp_path: Path) -> Path:
    root = tmp_path / "canonical"
    root.mkdir()
    repo = Repo.init(root)
    for relative, raw in FILES.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    unrelated = root / "tests/test_unrelated.py"
    unrelated.write_text("UNRELATED = True\n", encoding="utf-8")
    repo.index.add([*FILES, "tests/test_unrelated.py"])
    (root / "src/untracked.py").write_text("UNTRACKED = True\n", encoding="utf-8")
    return root


def test_isolation_copies_only_bound_tracked_inputs_without_shared_inodes(
    verifier: ModuleType, repository: Path, tmp_path: Path
) -> None:
    output = tmp_path / "evidence"
    isolated = api(verifier, "isolate_copy")(repository, output)
    assert isolated == output / "repository"
    for relative, raw in FILES.items():
        assert (isolated / relative).read_bytes() == raw
        assert (repository / relative).read_bytes() == raw
        assert (isolated / relative).stat().st_ino != (repository / relative).stat().st_ino
    assert not (isolated / "src/untracked.py").exists()
    assert not (isolated / "tests/test_unrelated.py").exists()
    assert not (isolated / ".git").exists()
    assert (output / "ISOLATION.json").stat().st_size > 0


def test_isolation_refuses_reusing_an_existing_output(
    verifier: ModuleType, repository: Path, tmp_path: Path
) -> None:
    output = tmp_path / "evidence"
    output.mkdir()
    marker = output / "preserved.txt"
    marker.write_bytes(b"preserve me")
    with pytest.raises((ValueError, FileExistsError), match=r"exist|exclusive|reuse"):
        api(verifier, "isolate_copy")(repository, output)
    assert marker.read_bytes() == b"preserve me"


def test_isolation_refuses_an_output_within_the_canonical_repository(
    verifier: ModuleType, repository: Path
) -> None:
    output = repository / "mutation-output"
    with pytest.raises(ValueError, match=r"canonical|repository|outside"):
        api(verifier, "isolate_copy")(repository, output)
    assert not output.exists()
    assert (repository / SOURCE).read_bytes() == FILES[SOURCE]


def test_isolation_refuses_a_selected_source_symlink(
    verifier: ModuleType, repository: Path, tmp_path: Path
) -> None:
    outside = tmp_path / "outside.py"
    outside.write_bytes(FILES[SOURCE])
    (repository / SOURCE).unlink()
    (repository / SOURCE).symlink_to(outside)
    with pytest.raises(ValueError, match=r"symlink|regular|escape"):
        api(verifier, "isolate_copy")(repository, tmp_path / "evidence")
    assert outside.read_bytes() == FILES[SOURCE]


def test_mutation_changes_only_the_isolated_source_and_records_both_hashes(
    verifier: ModuleType, repository: Path, tmp_path: Path
) -> None:
    isolated = api(verifier, "isolate_copy")(repository, tmp_path / "evidence")
    result = api(verifier, "apply_mutation")(isolated, SOURCE, "TARGET = True", "TARGET = False", 1)
    expected = b"TARGET = False\n"
    assert (repository / SOURCE).read_bytes() == FILES[SOURCE]
    assert (isolated / SOURCE).read_bytes() == expected
    assert result["original_source_sha256"] == hashlib.sha256(FILES[SOURCE]).hexdigest()
    assert result["mutated_source_sha256"] == hashlib.sha256(expected).hexdigest()


@pytest.mark.parametrize("expected", [0, 2])
def test_stale_mutation_anchor_never_changes_either_copy(
    verifier: ModuleType, repository: Path, tmp_path: Path, expected: int
) -> None:
    isolated = api(verifier, "isolate_copy")(repository, tmp_path / "evidence")
    with pytest.raises(ValueError, match=r"anchor|target|count"):
        api(verifier, "apply_mutation")(
            isolated, SOURCE, "TARGET = True", "TARGET = False", expected
        )
    assert (repository / SOURCE).read_bytes() == FILES[SOURCE]
    assert (isolated / SOURCE).read_bytes() == FILES[SOURCE]


@pytest.mark.parametrize("relative", ["../outside.py", "/outside.py", "src/../../outside.py"])
def test_mutation_refuses_paths_outside_the_bound_source_inventory(
    verifier: ModuleType, repository: Path, tmp_path: Path, relative: str
) -> None:
    isolated = api(verifier, "isolate_copy")(repository, tmp_path / "evidence")
    with pytest.raises(ValueError, match=r"path|bound|inventory|relative"):
        api(verifier, "apply_mutation")(isolated, relative, "True", "False", 1)
    assert (repository / SOURCE).read_bytes() == FILES[SOURCE]
    assert (isolated / SOURCE).read_bytes() == FILES[SOURCE]


def test_mutation_refuses_the_original_repository_even_with_a_matching_anchor(
    verifier: ModuleType, repository: Path
) -> None:
    with pytest.raises(ValueError, match=r"isolat|receipt|scratch"):
        api(verifier, "apply_mutation")(repository, SOURCE, "TARGET = True", "TARGET = False", 1)
    assert (repository / SOURCE).read_bytes() == FILES[SOURCE]
