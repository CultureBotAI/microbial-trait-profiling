"""Unavailable historical probe cases still publish honest summaries (#1303)."""

import importlib.util
import json
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest
from git.exc import GitCommandError

ROOT = Path(__file__).resolve().parents[1]
HELPER = Path("scripts/adhoc/issue1295_1300_probe_loader_mutants.py")


def load_probe(root: Path = ROOT) -> ModuleType:
    spec = importlib.util.spec_from_file_location("probe1303", root / HELPER)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("unavailable", [("head",), ("base",), ("head", "base"), ()])
def test_summary_preserves_available_cases(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, unavailable: tuple[str, ...]
) -> None:
    probe = load_probe()
    repo = tmp_path / "repo"
    target = repo / probe.TARGET
    target.parent.mkdir(parents=True)
    target.write_text("source\n")
    monkeypatch.setattr(probe, "ROOT", repo)
    monkeypatch.setattr(
        probe, "REVISIONS", {"head": ("head", ("target",)), "base": ("base", ("target",))}
    )
    monkeypatch.setattr(probe, "CASES", ("fixed", "head", "base"))
    calls = []

    def revision(rev: str) -> str:
        if rev in unavailable:
            raise probe.UnavailableRevision(rev)
        return "historical\n"

    def execute(case_dir: Path) -> tuple[int, SimpleNamespace]:
        calls.append(case_dir.name)
        outcome = "passed" if case_dir.name == "fixed" else "failed"
        return int(outcome == "failed"), SimpleNamespace(results={"target": outcome}, reasons={})

    monkeypatch.setattr(probe, "revision_source", revision)
    monkeypatch.setattr(probe, "execute_tests", execute)
    output = tmp_path / "output"
    code = probe.main(["--output-dir", str(output)])
    assert code == int(bool(unavailable))
    summary = json.loads((output / "summary.json").read_text())
    assert [row["case"] for row in summary] == ["fixed", "head", "base"]
    assert calls == [name for name in ("fixed", "head", "base") if name not in unavailable]
    for row in summary:
        missing = row["case"] in unavailable
        assert row["status"] == ("unavailable" if missing else "verified")
        assert row["verified"] is not missing
        assert json.loads((output / row["case"] / "evidence.json").read_text()) == row
        if missing:
            assert row["revision"] == row["case"]
            assert row["pytest_exit"] is None
            assert row["selected"] == 0
    assert target.read_text() == "source\n"


@pytest.mark.parametrize("mode", ["missing", "git_other", "permission"])
def test_revision_errors_are_narrow(monkeypatch: pytest.MonkeyPatch, mode: str) -> None:
    probe = load_probe()

    class Repository:
        @property
        def git(self) -> "Repository":
            return self

        def __enter__(self) -> "Repository":
            return self

        def __exit__(self, *args: object) -> None:
            pass

        def rev_parse(self, *_args: str) -> None:
            if mode == "missing":
                raise GitCommandError("git rev-parse", 1)
            if mode == "git_other":
                raise GitCommandError("git rev-parse", 128)
            raise PermissionError("cannot read object database")

    monkeypatch.setattr(probe, "Repo", lambda _root: Repository())
    error = None
    try:
        probe.revision_source("revision")
    except Exception as exc:
        error = exc
    expected = {
        "missing": probe.UnavailableRevision,
        "git_other": GitCommandError,
        "permission": PermissionError,
    }[mode]
    assert isinstance(error, expected)


def test_real_missing_revision() -> None:
    probe = load_probe()
    with pytest.raises(probe.UnavailableRevision):
        probe.revision_source("0" * 40)
