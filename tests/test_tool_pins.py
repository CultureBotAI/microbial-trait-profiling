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


def _pyproject_pins() -> dict[str, str]:
    data = tomllib.loads((ROOT / "pyproject.toml").read_text())
    pins: dict[str, str] = {}
    for group in data.get("dependency-groups", {}).values():
        for spec in group:
            match = re.fullmatch(r"([A-Za-z0-9_.-]+)(?:\[[^\]]*\])?==([^\s;]+)", spec)
            if match:
                pins[match.group(1)] = match.group(2)
    return pins


def _hook_revs() -> dict[str, str]:
    config = yaml.safe_load((ROOT / ".pre-commit-config.yaml").read_text())
    return {repo["repo"]: repo["rev"] for repo in config["repos"] if "rev" in repo}


@pytest.mark.parametrize("name", sorted(PINNED_HOOKS))
def test_pyproject_pin_matches_hook_rev(name: str) -> None:
    pin = _pyproject_pins().get(name)
    assert pin is not None, f"{name} is no longer pinned with == in pyproject.toml"
    rev = _hook_revs()[PINNED_HOOKS[name]]
    assert rev.removeprefix("v") == pin, (
        f"{name}: pyproject.toml pins {pin} but .pre-commit-config.yaml uses rev {rev}"
    )
