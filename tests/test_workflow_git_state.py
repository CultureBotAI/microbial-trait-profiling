"""Compare issue-635 repository provenance with Git's observed repository state."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import pytest
from dulwich.porcelain import status
from dulwich.repo import Repo

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


@pytest.mark.parametrize(
    "case", ["tracked", "changed_target", "real_directory", "untracked", "untracked_empty_target"]
)
def test_directory_symlink_is_classified_by_its_own_index_entry(
    repository: Path, case: str
) -> None:
    target = repository / "target"
    target.mkdir()
    (target / "content.txt").write_text("tracked target\n", encoding="utf-8")
    other = repository / "other-target"
    other.mkdir()
    (other / "other.txt").write_text("another tracked target\n", encoding="utf-8")
    git(repository, "add", "target/content.txt", "other-target/other.txt")
    git(repository, "commit", "-m", "directory target fixture")
    link_target = "target"
    if case == "untracked_empty_target":
        (repository / "empty-target").mkdir()
        link_target = "empty-target"
    (repository / "directory-link").symlink_to(link_target, target_is_directory=True)
    if not case.startswith("untracked"):
        git(repository, "add", "directory-link")
        git(repository, "commit", "-m", "tracked directory symlink fixture")
    if case == "changed_target":
        (repository / "directory-link").unlink()
        (repository / "directory-link").symlink_to("other-target", target_is_directory=True)
    elif case == "real_directory":
        (repository / "directory-link").unlink()
        (repository / "directory-link").mkdir()
        (repository / "directory-link" / "new.txt").write_text("new file\n", encoding="utf-8")
    expected_dirty = bool(git(repository, "status", "--porcelain", "--untracked-files=normal"))
    assert expected_dirty == (case != "tracked")
    observed_by_mode: dict[str, bool] = {}
    with Repo(str(repository)) as repo:
        for mode in ("normal", "all"):
            changes = status(repo, ignored=False, untracked_files=mode, optional_locks=False)
            observed_by_mode[mode] = bool(
                changes.unstaged or changes.untracked or any(changes.staged.values())
            )
    print(json.dumps({"case": case, "git_dirty": expected_dirty, "dulwich": observed_by_mode}))
    assert runtime.repository_state(repository)[1] == expected_dirty


def test_repository_provenance_does_not_refresh_or_rewrite_the_index(repository: Path) -> None:
    git_dir = repository / ".git"
    before = {
        path.relative_to(git_dir): (path.read_bytes(), path.stat().st_mtime_ns)
        for path in git_dir.rglob("*")
        if path.is_file()
    }
    tracked = repository / "tracked.txt"
    original = tracked.stat()
    os.utime(tracked, ns=(original.st_atime_ns, original.st_mtime_ns + 1_000_000_000))
    assert runtime.repository_state(repository)[1] is False
    after = {
        path.relative_to(git_dir): (path.read_bytes(), path.stat().st_mtime_ns)
        for path in git_dir.rglob("*")
        if path.is_file()
    }
    assert after == before


def test_repository_provenance_works_while_another_process_holds_the_index_lock(
    repository: Path,
) -> None:
    lock = repository / ".git" / "index.lock"
    lock.write_bytes(b"another process owns this lock\n")
    before = lock.read_bytes()
    assert runtime.repository_state(repository) == (git(repository, "rev-parse", "HEAD"), False)
    assert lock.read_bytes() == before


def test_repository_provenance_does_not_search_parent_directories(repository: Path) -> None:
    nested = repository / "nested"
    nested.mkdir()
    with pytest.raises(ValueError, match="repository provenance"):
        runtime.repository_state(nested)


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
