"""Tools pinned in pyproject.toml must match their pre-commit hook revs (#410, #420)."""

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
