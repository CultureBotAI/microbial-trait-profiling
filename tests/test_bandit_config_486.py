"""Regression tests for plain bandit configuration without the no-op [toml] extra (#486)."""

from __future__ import annotations

import re
import tomllib
from pathlib import Path
from typing import Any

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
BANDIT_REPO_URL = "https://github.com/PyCQA/bandit"


def _iter_pyproject_specs(pyproject_data: dict[str, Any]) -> list[str]:
    groups = [
        pyproject_data.get("project", {}).get("dependencies", []),
        *pyproject_data.get("project", {}).get("optional-dependencies", {}).values(),
        *pyproject_data.get("dependency-groups", {}).values(),
    ]
    return [spec for group in groups for spec in group if isinstance(spec, str)]


def validate_bandit_configuration(
    *,
    pyproject_path: Path = ROOT / "pyproject.toml",
    pre_commit_path: Path = ROOT / ".pre-commit-config.yaml",
    uv_lock_path: Path = ROOT / "uv.lock",
) -> None:
    """Verify bandit is configured without the Python <3.11 [toml] extra across repo configs."""
    pyproject_data = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))
    requires_python = str(pyproject_data.get("project", {}).get("requires-python", ""))
    assert requires_python == ">=3.11", (
        f"Expected project.requires-python == '>=3.11', got {requires_python!r}"
    )

    all_specs = _iter_pyproject_specs(pyproject_data)
    bandit_specs = [
        spec for spec in all_specs if re.match(r"^\s*bandit(?:\b|\[|[><=!~])", spec, re.IGNORECASE)
    ]
    assert len(bandit_specs) == 1, (
        f"Expected exactly one bandit dependency specifier in pyproject.toml, got {bandit_specs!r}"
    )
    bandit_spec = bandit_specs[0].strip()
    assert "[" not in bandit_spec, (
        f"pyproject.toml must not specify extras on bandit (#486): {bandit_spec!r}"
    )
    spec_match = re.fullmatch(r"bandit(>=|==)([0-9]+\.[0-9]+\.[0-9]+)", bandit_spec)
    assert spec_match is not None, (
        f"Expected plain bandit version specifier in pyproject.toml, got {bandit_spec!r}"
    )
    pyproject_version = spec_match.group(2)

    qa_specs = [
        spec
        for spec in pyproject_data.get("dependency-groups", {}).get("qa", [])
        if isinstance(spec, str)
    ]
    assert bandit_spec in qa_specs, f"Expected {bandit_spec!r} in dependency-groups.qa"
    assert "bandit" in pyproject_data.get("tool", {}), "Missing [tool.bandit] in pyproject.toml"

    lock_data = tomllib.loads(uv_lock_path.read_text(encoding="utf-8"))
    locked_bandit = [
        pkg
        for pkg in lock_data.get("package", [])
        if isinstance(pkg, dict) and pkg.get("name") == "bandit" and "version" in pkg
    ]
    assert len(locked_bandit) == 1, f"Expected one locked bandit package, got {locked_bandit!r}"
    locked_version = str(locked_bandit[0]["version"])
    assert locked_version == pyproject_version, (
        f"pyproject.toml bandit floor {pyproject_version!r} does not match "
        f"uv.lock version {locked_version!r}"
    )

    project_name = str(pyproject_data.get("project", {}).get("name", ""))
    mtp_pkgs = [
        pkg
        for pkg in lock_data.get("package", [])
        if isinstance(pkg, dict) and pkg.get("name") == project_name
    ]
    assert len(mtp_pkgs) == 1, f"Missing {project_name!r} root package in uv.lock"
    qa_requires = mtp_pkgs[0].get("metadata", {}).get("requires-dev", {}).get("qa", [])
    bandit_requires = [
        entry for entry in qa_requires if isinstance(entry, dict) and entry.get("name") == "bandit"
    ]
    assert len(bandit_requires) == 1, (
        f"Expected one bandit entry in uv.lock requires-dev.qa, got {bandit_requires!r}"
    )
    assert "extras" not in bandit_requires[0], (
        f"uv.lock requires-dev.qa must not retain bandit extras (#486): {bandit_requires[0]!r}"
    )

    pre_commit_data = yaml.safe_load(pre_commit_path.read_text(encoding="utf-8"))
    bandit_repos = [
        repo
        for repo in pre_commit_data.get("repos", [])
        if isinstance(repo, dict)
        and str(repo.get("repo", "")).rstrip("/").removesuffix(".git") == BANDIT_REPO_URL
    ]
    assert len(bandit_repos) == 1, (
        f"Expected one {BANDIT_REPO_URL} repo in .pre-commit-config.yaml, got {bandit_repos!r}"
    )
    bandit_repo = bandit_repos[0]
    rev = str(bandit_repo.get("rev", "")).removeprefix("v")
    assert rev == locked_version, (
        f"Bandit pre-commit rev {rev!r} does not match uv.lock version {locked_version!r}"
    )

    hooks = [
        hook
        for hook in bandit_repo.get("hooks", [])
        if isinstance(hook, dict) and hook.get("id") == "bandit"
    ]
    assert len(hooks) == 1, f"Expected one bandit hook entry, got {hooks!r}"
    bandit_hook = hooks[0]
    assert bandit_hook.get("args") == ["-c", "pyproject.toml"], (
        f"Bandit hook must pass ['-c', 'pyproject.toml'], got {bandit_hook.get('args')!r}"
    )
    additional_deps = bandit_hook.get("additional_dependencies", [])
    assert not additional_deps, (
        "Bandit pre-commit hook must not declare redundant or [toml] "
        f"additional_dependencies (#486): {additional_deps!r}"
    )


def test_repo_bandit_config_drops_noop_toml_extra() -> None:
    """Repository pyproject.toml, uv.lock, and .pre-commit-config.yaml must not use bandit[toml]."""
    validate_bandit_configuration()


def test_validate_bandit_configuration_rejects_regressions(tmp_path: Path) -> None:
    """Validator must fail closed on [toml] extras or missing -c pyproject.toml args (#486)."""
    pyproject_path = tmp_path / "pyproject.toml"
    pre_commit_path = tmp_path / ".pre-commit-config.yaml"
    uv_lock_path = tmp_path / "uv.lock"

    valid_pyproject = (
        '[project]\nname = "microbial-trait-profiling"\nrequires-python = ">=3.11"\n'
        '[dependency-groups]\nqa = ["bandit>=1.9.4"]\n'
        '[tool.bandit]\nexclude_dirs = ["tests"]\n'
    )
    valid_pre_commit = (
        "repos:\n"
        "  - repo: https://github.com/PyCQA/bandit\n"
        "    rev: 1.9.4\n"
        "    hooks:\n"
        "      - id: bandit\n"
        "        args: [-c, pyproject.toml]\n"
    )
    valid_uv_lock = (
        '[[package]]\nname = "bandit"\nversion = "1.9.4"\n\n'
        '[[package]]\nname = "microbial-trait-profiling"\nversion = "0.2.0"\n'
        "[package.metadata.requires-dev]\n"
        'qa = [{ name = "bandit", specifier = ">=1.9.4" }]\n'
    )

    pyproject_path.write_text(valid_pyproject, encoding="utf-8")
    pre_commit_path.write_text(valid_pre_commit, encoding="utf-8")
    uv_lock_path.write_text(valid_uv_lock, encoding="utf-8")
    validate_bandit_configuration(
        pyproject_path=pyproject_path,
        pre_commit_path=pre_commit_path,
        uv_lock_path=uv_lock_path,
    )

    pyproject_path.write_text(
        valid_pyproject.replace("bandit>=1.9.4", "bandit[toml]>=1.9.4"),
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="must not specify extras on bandit"):
        validate_bandit_configuration(
            pyproject_path=pyproject_path,
            pre_commit_path=pre_commit_path,
            uv_lock_path=uv_lock_path,
        )
    pyproject_path.write_text(valid_pyproject, encoding="utf-8")

    uv_lock_path.write_text(
        valid_uv_lock.replace(
            '{ name = "bandit", specifier = ">=1.9.4" }',
            '{ name = "bandit", extras = ["toml"], specifier = ">=1.9.4" }',
        ),
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="must not retain bandit extras"):
        validate_bandit_configuration(
            pyproject_path=pyproject_path,
            pre_commit_path=pre_commit_path,
            uv_lock_path=uv_lock_path,
        )
    uv_lock_path.write_text(valid_uv_lock, encoding="utf-8")

    pre_commit_path.write_text(
        valid_pre_commit + '        additional_dependencies: ["bandit[toml]"]\n',
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="additional_dependencies"):
        validate_bandit_configuration(
            pyproject_path=pyproject_path,
            pre_commit_path=pre_commit_path,
            uv_lock_path=uv_lock_path,
        )

    pre_commit_path.write_text(
        valid_pre_commit.replace("        args: [-c, pyproject.toml]\n", ""),
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="pyproject\\.toml"):
        validate_bandit_configuration(
            pyproject_path=pyproject_path,
            pre_commit_path=pre_commit_path,
            uv_lock_path=uv_lock_path,
        )
