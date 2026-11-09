"""Tools pinned in pyproject.toml must match their pre-commit hook revs (#410, #420, #445, #449)."""

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


def _normalize_repo_url(url: str) -> str:
    normalized = url.strip().casefold().rstrip("/")
    normalized = normalized.removesuffix(".git").rstrip("/")
    return normalized


def _pyproject_pins(pyproject_path: Path = ROOT / "pyproject.toml") -> dict[str, set[str]]:
    data = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))
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
            pkg_match = re.match(
                r"\s*([A-Za-z0-9_.-]+)\s*(?:\[[^\]]*\])?\s*(.*)\Z",
                spec,
                re.DOTALL,
            )
            if not pkg_match:
                continue
            name = _normalize(pkg_match.group(1))
            rest = pkg_match.group(2).strip()
            version_part, sep, _ = rest.partition(";")
            version_part = version_part.strip()
            has_marker = bool(sep)
            exact_match = (
                None
                if has_marker or version_part.startswith("@")
                else re.fullmatch(
                    r"==(?!=)\s*([^\s;,)]+)|\(\s*==(?!=)\s*([^\s;,)]+)\s*\)",
                    version_part,
                )
            )
            has_equals = not version_part.startswith("@") and bool(
                re.search(r"(?<![=!<>~])==(?!=)", version_part)
            )
            if exact_match:
                pin_version = exact_match.group(1) or exact_match.group(2)
                pins.setdefault(name, set()).add(pin_version)
            elif name in PINNED_HOOKS and (has_marker or has_equals):
                raise AssertionError(
                    f"{name}: marker-gated or compound pin in pyproject.toml "
                    f"is not unconditional: {spec!r}"
                )
    return pins


def _hook_revs(config_path: Path = ROOT / ".pre-commit-config.yaml") -> dict[str, list[str]]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    revs: dict[str, list[str]] = {}
    for repo in config["repos"]:
        if "rev" in repo:
            key = _normalize_repo_url(str(repo["repo"]))
            revs.setdefault(key, []).append(str(repo["rev"]))
    return revs


def _single_hook_rev(repo_url: str, config_path: Path = ROOT / ".pre-commit-config.yaml") -> str:
    revs = _hook_revs(config_path)
    key = _normalize_repo_url(repo_url)
    assert key in revs, f"{repo_url} is missing from .pre-commit-config.yaml"
    repo_revs = revs[key]
    assert len(repo_revs) == 1, (
        f"{repo_url} appears {len(repo_revs)} times in .pre-commit-config.yaml: {repo_revs}"
    )
    return repo_revs[0]


def _check_pyproject_pin_matches_hook_rev(
    name: str,
    *,
    pyproject_path: Path = ROOT / "pyproject.toml",
    config_path: Path = ROOT / ".pre-commit-config.yaml",
) -> None:
    pins = _pyproject_pins(pyproject_path).get(name)
    assert pins, f"{name} is no longer pinned with == in pyproject.toml"
    rev = _single_hook_rev(PINNED_HOOKS[name], config_path)
    assert pins == {rev.removeprefix("v")}, (
        f"{name}: pyproject.toml pins {sorted(pins)} but .pre-commit-config.yaml uses rev {rev}"
    )


@pytest.mark.parametrize("name", sorted(PINNED_HOOKS))
def test_pyproject_pin_matches_hook_rev(name: str) -> None:
    _check_pyproject_pin_matches_hook_rev(name)


def _locked_versions(lock_path: Path = ROOT / "uv.lock") -> dict[str, set[str]]:
    lock = tomllib.loads(lock_path.read_text(encoding="utf-8"))
    versions: dict[str, set[str]] = {}
    for pkg in lock["package"]:
        if "version" in pkg:
            versions.setdefault(_normalize(pkg["name"]), set()).add(pkg["version"])
    return versions


def _check_hook_rev_matches_locked_version(
    name: str,
    *,
    lock_path: Path = ROOT / "uv.lock",
    config_path: Path = ROOT / ".pre-commit-config.yaml",
) -> None:
    locked = _locked_versions(lock_path).get(name)
    assert locked, f"{name} is not in uv.lock"
    assert len(locked) == 1, f"{name} resolves to several versions in uv.lock: {sorted(locked)}"
    rev = _single_hook_rev(LOCKED_HOOKS[name], config_path)
    assert {rev.removeprefix("v")} == locked, (
        f"{name}: uv.lock resolves {sorted(locked)} but .pre-commit-config.yaml uses rev {rev}; "
        "on a Dependabot bump, set the hook rev to the new version on the same PR"
    )


@pytest.mark.parametrize("name", sorted(LOCKED_HOOKS))
def test_hook_rev_matches_locked_version(name: str) -> None:
    _check_hook_rev_matches_locked_version(name)


def test_duplicate_hook_repo_fails_closed(tmp_path: Path) -> None:
    """Listing a hook repo twice in .pre-commit-config.yaml must not pass via last-wins (#449)."""
    pyproject_path = tmp_path / "pyproject.toml"
    pyproject_path.write_text(
        '[dependency-groups]\nqa = ["uv==0.13.4"]\n',
        encoding="utf-8",
    )
    lock_path = tmp_path / "uv.lock"
    lock_path.write_text(
        '[[package]]\nname = "ruff"\nversion = "0.17.0"\n',
        encoding="utf-8",
    )
    stale_then_matching = tmp_path / "stale_then_matching.yaml"
    stale_then_matching.write_text(
        "repos:\n"
        "  - repo: https://github.com/astral-sh/uv-pre-commit\n"
        "    rev: 0.12.22\n"
        "  - repo: https://github.com/astral-sh/uv-pre-commit\n"
        "    rev: 0.13.4\n",
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="appears 2 times"):
        _check_pyproject_pin_matches_hook_rev(
            "uv",
            pyproject_path=pyproject_path,
            config_path=stale_then_matching,
        )

    identical_duplicate = tmp_path / "identical_duplicate.yaml"
    identical_duplicate.write_text(
        "repos:\n"
        "  - repo: https://github.com/astral-sh/ruff-pre-commit\n"
        "    rev: v0.17.0\n"
        "  - repo: https://github.com/astral-sh/ruff-pre-commit\n"
        "    rev: v0.17.0\n",
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="appears 2 times"):
        _check_hook_rev_matches_locked_version(
            "ruff",
            lock_path=lock_path,
            config_path=identical_duplicate,
        )

    for alias_url in (
        "https://github.com/astral-sh/uv-pre-commit.git",
        "https://github.com/Astral-SH/uv-pre-commit",
        "https://github.com/astral-sh/uv-pre-commit/",
    ):
        alias_duplicate = tmp_path / "alias_duplicate.yaml"
        alias_duplicate.write_text(
            f"repos:\n"
            f"  - repo: {alias_url}\n"
            f"    rev: 0.12.22\n"
            f"  - repo: https://github.com/astral-sh/uv-pre-commit\n"
            f"    rev: 0.13.4\n",
            encoding="utf-8",
        )
        with pytest.raises(AssertionError, match="appears 2 times"):
            _check_pyproject_pin_matches_hook_rev(
                "uv",
                pyproject_path=pyproject_path,
                config_path=alias_duplicate,
            )

    dot_git_only = tmp_path / "dot_git_only.yaml"
    dot_git_only.write_text(
        "repos:\n  - repo: https://github.com/astral-sh/uv-pre-commit.git\n    rev: 0.13.4\n",
        encoding="utf-8",
    )
    _check_pyproject_pin_matches_hook_rev(
        "uv",
        pyproject_path=pyproject_path,
        config_path=dot_git_only,
    )


@pytest.mark.parametrize(
    "bad_spec",
    [
        "uv==0.13.4; python_version<'3.0'",
        "uv==0.13.4; python_version>='3.11'",
        "uv==0.13.4,!=0.13.4",
        "uv>=0.13.0,==0.13.4",
    ],
)
def test_marker_gated_or_compound_pyproject_pin_fails_closed(
    tmp_path: Path,
    bad_spec: str,
) -> None:
    """Marker-gated or compound == specs must not satisfy unconditional tool pins (#449)."""
    pyproject_path = tmp_path / "pyproject.toml"
    pyproject_path.write_text(
        f'[dependency-groups]\nqa = ["{bad_spec}"]\n',
        encoding="utf-8",
    )
    config_path = tmp_path / ".pre-commit-config.yaml"
    config_path.write_text(
        "repos:\n  - repo: https://github.com/astral-sh/uv-pre-commit\n    rev: 0.13.4\n",
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="marker-gated or compound pin"):
        _check_pyproject_pin_matches_hook_rev(
            "uv",
            pyproject_path=pyproject_path,
            config_path=config_path,
        )

    side_by_side = tmp_path / "side_by_side_pyproject.toml"
    side_by_side.write_text(
        '[dependency-groups]\ndev = ["uv>=0.12.0; python_version<\'3.0\'"]\nqa = ["uv==0.13.4"]\n',
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="marker-gated or compound pin"):
        _check_pyproject_pin_matches_hook_rev(
            "uv",
            pyproject_path=side_by_side,
            config_path=config_path,
        )


def test_unrelated_markers_urls_and_parenthesised_pin_pass(tmp_path: Path) -> None:
    """Unrelated marker specs, URL fragments with ==, and PEP 508 (==...) pins must pass (#2093)."""
    pyproject_path = tmp_path / "pyproject.toml"
    pyproject_path.write_text(
        "[project]\n"
        "dependencies = [\n"
        "  \"pywin32>=306; sys_platform == 'win32'\",\n"
        "  \"colorama>=0.4; platform_system == 'Windows'\",\n"
        "  \"tomli==2.0.1; python_version < '3.11'\",\n"
        '  "pkg @ https://example.org/archive.tar.gz#sha256==abc",\n'
        '  "requests>=2,==2.32.0",\n'
        "]\n"
        "[dependency-groups]\n"
        'qa = ["uv (==0.13.4)"]\n',
        encoding="utf-8",
    )
    config_path = tmp_path / ".pre-commit-config.yaml"
    config_path.write_text(
        "repos:\n  - repo: https://github.com/astral-sh/uv-pre-commit\n    rev: 0.13.4\n",
        encoding="utf-8",
    )
    assert _pyproject_pins(pyproject_path) == {"uv": {"0.13.4"}}
    _check_pyproject_pin_matches_hook_rev(
        "uv",
        pyproject_path=pyproject_path,
        config_path=config_path,
    )
