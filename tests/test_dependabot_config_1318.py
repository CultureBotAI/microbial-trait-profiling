"""Regression tests for Dependabot linkml + linkml-runtime group and pydantic floor (#1318)."""

from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
MIN_PYDANTIC_VERSION: tuple[int, int, int] = (2, 13, 5)
EXPECTED_LINKML_GROUP_PATTERNS: list[str] = ["linkml", "linkml-runtime"]


def _parse_version_tuple(version_str: str) -> tuple[int, int, int]:
    parts = version_str.strip().split(".")
    assert len(parts) == 3, f"Expected three-part numeric version string, got {version_str!r}"
    assert all(part.isdigit() for part in parts), (
        f"Expected numeric version components, got {version_str!r}"
    )
    return (int(parts[0]), int(parts[1]), int(parts[2]))


def validate_dependabot_linkml_group(
    *,
    dependabot_path: Path = ROOT / ".github" / "dependabot.yml",
    uv_lock_path: Path = ROOT / "uv.lock",
) -> None:
    """Verify .github/dependabot.yml groups linkml + linkml-runtime and locks pydantic >= 2.13.5."""
    data = yaml.safe_load(dependabot_path.read_text(encoding="utf-8"))
    assert isinstance(data, dict), "Expected .github/dependabot.yml top-level mapping"
    assert data.get("version") == 2, (
        f"Expected Dependabot version == 2, got {data.get('version')!r}"
    )

    updates = data.get("updates")
    assert isinstance(updates, list), "Expected updates list in .github/dependabot.yml"
    assert bool(updates), "Expected non-empty updates list in .github/dependabot.yml"

    uv_entries = [
        entry
        for entry in updates
        if isinstance(entry, dict)
        and entry.get("package-ecosystem") == "uv"
        and entry.get("directory") == "/"
    ]
    assert len(uv_entries) == 1, (
        f"Expected exactly one uv ecosystem entry for '/', got {uv_entries!r}"
    )
    uv_entry = uv_entries[0]

    for other_entry in updates:
        if isinstance(other_entry, dict) and other_entry is not uv_entry:
            other_groups = other_entry.get("groups", {})
            has_other_linkml = isinstance(other_groups, dict) and "linkml" in other_groups
            assert not has_other_linkml, (
                "groups.linkml must be configured under package-ecosystem 'uv', "
                "not another ecosystem"
            )

    groups = uv_entry.get("groups")
    assert isinstance(groups, dict), (
        "Expected groups.linkml under the 'uv' package-ecosystem entry "
        "in .github/dependabot.yml (#1318)"
    )
    assert "linkml" in groups, (
        "Expected groups.linkml under the 'uv' package-ecosystem entry "
        "in .github/dependabot.yml (#1318)"
    )
    linkml_group = groups["linkml"]
    assert isinstance(linkml_group, dict), (
        f"Expected groups.linkml to be a mapping, got {linkml_group!r}"
    )
    patterns = linkml_group.get("patterns")
    assert isinstance(patterns, list), (
        f"Expected groups.linkml.patterns to be a list of strings, got {patterns!r}"
    )
    assert all(isinstance(item, str) for item in patterns), (
        f"Expected groups.linkml.patterns to be a list of strings, got {patterns!r}"
    )
    assert patterns == EXPECTED_LINKML_GROUP_PATTERNS, (
        f"Expected groups.linkml.patterns == {EXPECTED_LINKML_GROUP_PATTERNS!r}, got {patterns!r}"
    )

    ignore_entries = uv_entry.get("ignore", [])
    assert isinstance(ignore_entries, list), (
        f"Expected uv ignore entry to be a list, got {ignore_entries!r}"
    )
    ignored_names = {
        str(item.get("dependency-name", "")) for item in ignore_entries if isinstance(item, dict)
    }
    for required_pkg in EXPECTED_LINKML_GROUP_PATTERNS:
        assert required_pkg not in ignored_names, (
            f"{required_pkg!r} must not appear in the uv ignore list (#1318): {ignored_names!r}"
        )

    lock_data: dict[str, Any] = tomllib.loads(uv_lock_path.read_text(encoding="utf-8"))
    packages = [pkg for pkg in lock_data.get("package", []) if isinstance(pkg, dict)]
    by_name = {
        str(pkg["name"]): str(pkg["version"])
        for pkg in packages
        if "name" in pkg and "version" in pkg
    }
    for required_pkg in ("linkml", "linkml-runtime", "pydantic"):
        assert required_pkg in by_name, f"Missing {required_pkg!r} in uv.lock"

    assert by_name["linkml"] == by_name["linkml-runtime"], (
        f"Expected locked linkml ({by_name['linkml']!r}) and linkml-runtime "
        f"({by_name['linkml-runtime']!r}) versions to match"
    )
    pydantic_version = _parse_version_tuple(by_name["pydantic"])
    assert pydantic_version >= MIN_PYDANTIC_VERSION, (
        f"Expected locked pydantic >= {MIN_PYDANTIC_VERSION} for linkml compatibility (#1318), "
        f"got {by_name['pydantic']!r}"
    )


def test_repo_dependabot_groups_linkml_and_linkml_runtime() -> None:
    """Repository .github/dependabot.yml must group linkml and linkml-runtime under uv (#1318)."""
    validate_dependabot_linkml_group()


def test_validate_dependabot_linkml_group_rejects_regressions(tmp_path: Path) -> None:
    """Validator must fail closed on missing groups, missing patterns, ignore, or low pydantic."""
    dependabot_path = tmp_path / "dependabot.yml"
    uv_lock_path = tmp_path / "uv.lock"

    valid_dependabot = (
        "version: 2\n"
        "updates:\n"
        '  - package-ecosystem: "github-actions"\n'
        '    directory: "/"\n'
        "    schedule:\n"
        '      interval: "weekly"\n'
        '  - package-ecosystem: "uv"\n'
        '    directory: "/"\n'
        "    schedule:\n"
        '      interval: "weekly"\n'
        "    groups:\n"
        "      linkml:\n"
        "        patterns:\n"
        '          - "linkml"\n'
        '          - "linkml-runtime"\n'
        "    ignore:\n"
        '      - dependency-name: "uv"\n'
    )
    valid_uv_lock = (
        '[[package]]\nname = "linkml"\nversion = "1.12.1"\n\n'
        '[[package]]\nname = "linkml-runtime"\nversion = "1.12.1"\n\n'
        '[[package]]\nname = "pydantic"\nversion = "2.14.2"\n'
    )

    dependabot_path.write_text(valid_dependabot, encoding="utf-8")
    uv_lock_path.write_text(valid_uv_lock, encoding="utf-8")
    validate_dependabot_linkml_group(
        dependabot_path=dependabot_path,
        uv_lock_path=uv_lock_path,
    )

    dependabot_path.write_text(
        valid_dependabot.replace('          - "linkml-runtime"\n', ""),
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="groups\\.linkml\\.patterns"):
        validate_dependabot_linkml_group(
            dependabot_path=dependabot_path,
            uv_lock_path=uv_lock_path,
        )
    dependabot_path.write_text(valid_dependabot, encoding="utf-8")

    dependabot_path.write_text(
        valid_dependabot + '      - dependency-name: "linkml"\n',
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="must not appear in the uv ignore list"):
        validate_dependabot_linkml_group(
            dependabot_path=dependabot_path,
            uv_lock_path=uv_lock_path,
        )
    dependabot_path.write_text(valid_dependabot, encoding="utf-8")

    uv_lock_path.write_text(
        valid_uv_lock.replace('version = "2.14.2"', 'version = "2.13.0"'),
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="Expected locked pydantic >="):
        validate_dependabot_linkml_group(
            dependabot_path=dependabot_path,
            uv_lock_path=uv_lock_path,
        )
