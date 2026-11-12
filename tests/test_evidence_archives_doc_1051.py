"""Pin Linux temporary-directory documentation and verifier isolation (Refs #1051)."""

import ast
import importlib.util
import inspect
import re
import sys
import tempfile
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_linux_temp_environment_documentation() -> None:
    text = (ROOT / "docs/evidence-archives.md").read_text()
    paragraph = " ".join(text.split("On Linux,", 1)[1].split("\n\n", 1)[0].split())
    assert paragraph.startswith("`/tmp` is refused only when it is the process temporary directory")
    sentence = paragraph.split("If ", 1)[1]
    assert "`$TMPDIR`, `$TEMP` or `$TMP` points somewhere else" in sentence
    assert "a Linux `/tmp` output root is **not** refused, so do not use one." in sentence
    tree = ast.parse(inspect.getsource(tempfile))
    functions = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "_candidate_tempdir_list"
    ]
    assert len(functions) == 1
    candidates = [
        node.value
        for node in ast.walk(functions[0])
        if isinstance(node, ast.Constant) and node.value in ("TMPDIR", "TEMP", "TMP")
    ]
    assert candidates == ["TMPDIR", "TEMP", "TMP"]
    assert re.findall(r"`\$(\w+)`", sentence) == candidates


@pytest.fixture
def verifier() -> Iterator[ModuleType]:
    spec = importlib.util.spec_from_file_location(
        "verifier_1051_test", ROOT / "scripts/adhoc/verify_1051_mutants_codex2_20261111.py"
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
        yield module
    finally:
        sys.modules.pop(spec.name, None)


@pytest.mark.parametrize("case", ["root", "nested", "inbound", "outbound"])
def test_scratch_containment_before_creation(
    case: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, verifier: ModuleType
) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    external = tmp_path / "external"
    external.mkdir()
    inbound = external / "alias"
    inbound.symlink_to(repo, target_is_directory=True)
    (repo / "alias").symlink_to(external, target_is_directory=True)
    paths = {
        "root": repo,
        "nested": repo / "new",
        "inbound": inbound / "new",
        "outbound": repo / "alias" / "new",
    }
    monkeypatch.setattr(verifier, "ROOT", repo)

    def forbidden(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("containment must precede mkdir and mkdtemp")

    monkeypatch.setattr(Path, "mkdir", forbidden)
    monkeypatch.setattr(verifier.tempfile, "mkdtemp", forbidden)
    with pytest.raises(ValueError, match="outside the repository root"):
        verifier.main(["--scratch-dir", str(paths[case]), "--canary"])
    assert not (repo / "new").exists()
    assert not (external / "new").exists()


@pytest.mark.parametrize("mode", ["internal", "setup"])
def test_verifier_rejects_invalid_pytest_session(
    mode: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, verifier: ModuleType
) -> None:
    (tmp_path / "tests").mkdir()
    relative = "tests/probe_1051.py"
    body = (
        "import pytest\n@pytest.fixture(autouse=True)\ndef broken():\n"
        "    raise ValueError('AssertionError in setup')\n"
        "def test_probe():\n    pass\n"
        if mode == "setup"
        else "def test_probe():\n    assert False, 'call assertion'\n"
    )
    (tmp_path / relative).write_text(body)
    if mode == "internal":
        (tmp_path / "conftest.py").write_text(
            "def pytest_runtest_logfinish(nodeid, location):\n"
            "    raise RuntimeError('internal error after call')\n"
        )
    monkeypatch.setattr(verifier, "TEST", relative)
    with monkeypatch.context() as context:
        context.chdir(tmp_path)
        context.setattr(sys, "path", list(sys.path))
        context.setenv("MPLBACKEND", "Agg")
        context.setenv("PYTHONPATH", str(tmp_path / "src"))
        verifier.worker(str(tmp_path), "test_probe")
    assert int((tmp_path / "exitcode.txt").read_text()) == (3 if mode == "internal" else 1)
    with pytest.raises(RuntimeError):
        verifier.validate_result(tmp_path, "test_probe", mutant=True)
    sys.modules.pop("probe_1051", None)
