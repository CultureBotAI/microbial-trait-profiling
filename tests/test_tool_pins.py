"""Tools pinned in pyproject.toml must match their pre-commit hook revs (#410, #420, #445)."""

import re
import tomllib
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]

# pyproject dependency name -> pre-commit repo URL whose rev must equal its pin
PINNED_HOOKS = {
    "uv": "https://github.com/astral-sh/uv-pre-commit",
    "detect-secrets": "https://github.com/Yelp/detect-secrets",  # pragma: allowlist secret
}

# Tools with a >= floor in pyproject: the hook rev must equal the version uv.lock resolves,
# so `uv run <tool>` and the pre-commit hook in CI run the same release.
LOCKED_HOOKS = {
    "ruff": "https://github.com/astral-sh/ruff-pre-commit",
    "bandit": "https://github.com/PyCQA/bandit",
}


def _normalize(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def _pyproject_pins() -> dict[str, set[str]]:
    data = tomllib.loads((ROOT / "pyproject.toml").read_text())
    groups = [
        data.get("project", {}).get("dependencies", []),
        *data.get("project", {}).get("optional-dependencies", {}).values(),
        *data.get("dependency-groups", {}).values(),
    ]
    pins: dict[str, set[str]] = {}
    for group in groups:
        for spec in group:
            if not isinstance(spec, str):  # PEP 735 {include-group = ...}
                continue
            match = re.match(r"\s*([A-Za-z0-9_.-]+)\s*(?:\[[^\]]*\])?\s*==\s*([^\s;,]+)", spec)
            if match:
                pins.setdefault(_normalize(match.group(1)), set()).add(match.group(2))
    return pins


def _hook_revs() -> dict[str, str]:
    config = yaml.safe_load((ROOT / ".pre-commit-config.yaml").read_text())
    return {repo["repo"]: repo["rev"] for repo in config["repos"] if "rev" in repo}


@pytest.mark.parametrize("name", sorted(PINNED_HOOKS))
def test_pyproject_pin_matches_hook_rev(name: str) -> None:
    pins = _pyproject_pins().get(name)
    assert pins, f"{name} is no longer pinned with == in pyproject.toml"
    rev = _hook_revs()[PINNED_HOOKS[name]]
    assert pins == {rev.removeprefix("v")}, (
        f"{name}: pyproject.toml pins {sorted(pins)} but .pre-commit-config.yaml uses rev {rev}"
    )


def _locked_versions() -> dict[str, str]:
    lock = tomllib.loads((ROOT / "uv.lock").read_text())
    return {_normalize(pkg["name"]): pkg["version"] for pkg in lock["package"] if "version" in pkg}


@pytest.mark.parametrize("name", sorted(LOCKED_HOOKS))
def test_hook_rev_matches_locked_version(name: str) -> None:
    locked = _locked_versions().get(name)
    assert locked, f"{name} is not in uv.lock"
    rev = _hook_revs()[LOCKED_HOOKS[name]]
    assert rev.removeprefix("v") == locked, (
        f"{name}: uv.lock resolves {locked} but .pre-commit-config.yaml uses rev {rev}"
    )
