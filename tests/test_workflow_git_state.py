"""Compare issue-629 repository provenance with Git's observed repository state."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

from mtp.workflow import runtime


def git(root: Path, *arguments: str) -> str:
    """Use Git as the independent reference implementation in temporary fixtures."""
    environment = {
        **os.environ,
        "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull,
    }
    completed = subprocess.run(
        [
            "git",
            "-c",
            f"core.hooksPath={os.devnull}",
            "-c",
            "commit.gpgsign=false",
            "-C",
            str(root),
            *arguments,
        ],
        capture_output=True,
        text=True,
        check=True,
        env=environment,
    )
    return completed.stdout.strip()


@pytest.fixture
def repository(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", os.devnull)
    monkeypatch.setenv("GIT_CONFIG_SYSTEM", os.devnull)
    root = tmp_path / "repository"
    root.mkdir()
    git(root, "init", "--initial-branch=main")
    git(root, "config", "user.name", "Workflow Fixture")
    git(root, "config", "user.email", "workflow-fixture@example.invalid")
    (root / "tracked.txt").write_text("original\n", encoding="utf-8")
    (root / ".gitignore").write_text("ignored/\n", encoding="utf-8")
    git(root, "add", "tracked.txt", ".gitignore")
    git(root, "commit", "-m", "fixture")
    return root


@pytest.mark.parametrize(
    "state", ["clean", "staged", "unstaged", "untracked", "ignored", "deleted"]
)
def test_repository_state_matches_git(repository: Path, state: str) -> None:
    if state in {"staged", "unstaged"}:
        (repository / "tracked.txt").write_text("modified\n", encoding="utf-8")
        if state == "staged":
            git(repository, "add", "tracked.txt")
    elif state == "untracked":
        (repository / "new.txt").write_text("untracked\n", encoding="utf-8")
    elif state == "ignored":
        ignored = repository / "ignored"
        ignored.mkdir()
        (ignored / "data.txt").write_text("ignored output\n", encoding="utf-8")
    elif state == "deleted":
        (repository / "tracked.txt").unlink()
    expected = (
        git(repository, "rev-parse", "HEAD"),
        bool(git(repository, "status", "--porcelain", "--untracked-files=normal")),
    )
    assert runtime.repository_state(repository) == expected
    assert expected[1] == (state not in {"clean", "ignored"})


def test_repository_state_supports_detached_head(repository: Path) -> None:
    git(repository, "checkout", "--detach", "HEAD")
    assert runtime.repository_state(repository) == (git(repository, "rev-parse", "HEAD"), False)


def test_repository_provenance_does_not_refresh_or_rewrite_the_index(repository: Path) -> None:
    index = repository / ".git" / "index"
    before = index.read_bytes()
    tracked = repository / "tracked.txt"
    original = tracked.stat()
    os.utime(tracked, ns=(original.st_atime_ns, original.st_mtime_ns + 1_000_000_000))
    assert runtime.repository_state(repository)[1] is False
    assert index.read_bytes() == before


def test_linked_worktree_uses_its_own_head_index_and_working_files(repository: Path) -> None:
    linked = repository.parent / "linked"
    git(repository, "worktree", "add", "--detach", str(linked), "HEAD")
    assert (linked / ".git").is_file()
    (repository / "tracked.txt").write_text("main worktree modified\n", encoding="utf-8")
    assert runtime.repository_state(repository)[1]
    assert runtime.repository_state(linked) == (git(linked, "rev-parse", "HEAD"), False)
    (linked / "tracked.txt").write_text("linked worktree modified\n", encoding="utf-8")
    git(linked, "add", "tracked.txt")
    assert runtime.repository_state(linked) == (git(linked, "rev-parse", "HEAD"), True)


def test_missing_repository_or_unborn_head_cannot_invent_provenance(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="repository provenance"):
        runtime.repository_state(tmp_path)
    git(tmp_path, "init", "--initial-branch=main")
    with pytest.raises(ValueError, match="repository provenance"):
        runtime.repository_state(tmp_path)
