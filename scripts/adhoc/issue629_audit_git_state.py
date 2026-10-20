"""Audit the Git metadata regression in issues #629/#635 (2026-10-19).

Inputs: a repository/worktree and optional workflow run JSON.
Output: JSON containing the exact staged, unstaged, and untracked paths observed
by Dulwich, native Git state, selected nonsecret configuration, and run identity.
This diagnostic deliberately does not modify source, schema, lock, or Git index.
Dulwich is a declared QA dependency; install it with the standard `just sync`.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path

from dulwich.porcelain import status
from dulwich.repo import Repo

from mtp.workflow.runtime import repository_state


def text(value: bytes | str) -> str:
    return value.decode("utf-8", "surrogateescape") if isinstance(value, bytes) else value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--run", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    with Repo(str(root)) as repository:
        changes = status(repository, ignored=False, untracked_files="normal", optional_locks=False)
        configuration: dict[str, str | None] = {}
        for key in (b"excludesFile", b"ignoreCase", b"precomposeUnicode", b"fileMode"):
            try:
                configuration[text(key)] = text(repository.get_config_stack().get(b"core", key))
            except KeyError:
                configuration[text(key)] = None
        report: dict[str, object] = {
            "root": str(root),
            "dulwich_version": importlib.metadata.version("dulwich"),
            "commit": repository.head().decode("ascii"),
            "staged": {
                text(key): [text(path) for path in paths] for key, paths in changes.staged.items()
            },
            "unstaged": [text(path) for path in changes.unstaged],
            "untracked": [text(path) for path in changes.untracked],
            "dirty": bool(changes.unstaged or changes.untracked or any(changes.staged.values())),
            "core_config": configuration,
        }
        all_changes = status(repository, ignored=False, untracked_files="all", optional_locks=False)
        report["dulwich_all"] = {
            "dirty": bool(
                all_changes.unstaged or all_changes.untracked or any(all_changes.staged.values())
            ),
            "untracked": [text(path) for path in all_changes.untracked],
        }
    native_commit, native_dirty = repository_state(root)
    report["native_git"] = {
        "gitpython_version": importlib.metadata.version("GitPython"),
        "commit": native_commit,
        "dirty": native_dirty,
    }
    if args.run:
        run = json.loads(args.run.read_text(encoding="utf-8"))
        report["run"] = {
            key: run[key]
            for key in (
                "run_id",
                "git_commit",
                "git_dirty",
                "started_at",
                "status",
                "code_sha256",
            )
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(args.output)
    print(json.dumps({key: report[key] for key in ("commit", "dirty", "core_config", "staged")}))
    print(json.dumps({"unstaged": report["unstaged"], "untracked": report["untracked"]}))
    print(json.dumps({"native_git": report["native_git"], "dulwich_all": report["dulwich_all"]}))


if __name__ == "__main__":
    main()
