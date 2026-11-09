"""Tests that README.md documents release creation, assets, and verification (#59)."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
README_MD = REPO_ROOT / "README.md"
RELEASE_WORKFLOW_YAML = REPO_ROOT / ".github" / "workflows" / "release.yaml"

FORBIDDEN_RELEASE_PHRASES: tuple[str, ...] = (
    "turbomam/bsdbng",
    "without explicit human approval",
    "without user approval",
    "without asking the user",
    "no approval required",
    "no human approval",
    "skip `just check`",
    "skip just check",
)


def _extract_h2_section(markdown: str, heading: str) -> str:
    """Return the markdown text under `## <heading>` up to the next `## ` heading."""
    pattern = re.compile(
        rf"^## {re.escape(heading)}\s*$\n(.*?)(?=^## |\Z)",
        re.MULTILINE | re.DOTALL,
    )
    match = pattern.search(markdown)
    if match is None:
        raise AssertionError(f"Missing '## {heading}' section in README.md")
    return match.group(1)


def _extract_h3_section(section_markdown: str, subheading: str) -> str:
    """Return the markdown text under `### <subheading>` up to the next `### ` heading."""
    pattern = re.compile(
        rf"^### {re.escape(subheading)}\s*$\n(.*?)(?=^### |\Z)",
        re.MULTILINE | re.DOTALL,
    )
    match = pattern.search(section_markdown)
    if match is None:
        raise AssertionError(f"Missing '### {subheading}' subsection under '## Releases'")
    return match.group(1)


def _load_release_workflow(workflow_path: Path) -> dict[Any, Any]:
    loaded = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    assert isinstance(loaded, dict), "release.yaml must parse to a mapping"
    return loaded


def _workflow_tag_triggers(workflow: dict[Any, Any]) -> list[str]:
    on_block = workflow.get("on") or workflow.get(True)
    assert isinstance(on_block, dict), "release.yaml must define an 'on' trigger mapping"
    push_block = on_block.get("push")
    assert isinstance(push_block, dict), "release.yaml must define 'on.push'"
    tags = push_block.get("tags")
    assert isinstance(tags, list), "release.yaml must define 'on.push.tags' as a list"
    assert tags, "release.yaml must define non-empty 'on.push.tags'"
    return [str(tag) for tag in tags]


def _workflow_expected_asset_templates(workflow: dict[Any, Any]) -> set[str]:
    jobs = workflow.get("jobs")
    assert isinstance(jobs, dict), "release.yaml must define 'jobs'"
    build_job = jobs.get("build-and-release")
    assert isinstance(build_job, dict), "release.yaml must define job 'build-and-release'"
    steps = build_job.get("steps")
    assert isinstance(steps, list), "build-and-release must define 'steps'"

    files_block: str | None = None
    for step in steps:
        if not isinstance(step, dict):
            continue
        uses = str(step.get("uses", ""))
        if uses.startswith("softprops/action-gh-release"):
            with_block = step.get("with")
            assert isinstance(with_block, dict), "action-gh-release step must have 'with'"
            raw_files = with_block.get("files")
            assert isinstance(raw_files, str), "action-gh-release 'with.files' must be a string"
            assert raw_files.strip(), "action-gh-release 'with.files' must be non-empty"
            files_block = raw_files
            break

    assert files_block is not None, "release.yaml must include a softprops/action-gh-release step"
    globs = [line.strip() for line in files_block.splitlines() if line.strip()]
    assert globs, "action-gh-release 'with.files' must list at least one asset pattern"

    templates: set[str] = set()
    for glob_pattern in globs:
        assert "*" in glob_pattern, f"Expected wildcard asset pattern, got {glob_pattern!r}"
        templates.add(glob_pattern.replace("*", "<tag>"))
    return templates


def _workflow_just_commands(workflow: dict[Any, Any]) -> list[str]:
    jobs = workflow.get("jobs", {})
    steps = jobs.get("build-and-release", {}).get("steps", [])
    commands: list[str] = []
    for step in steps:
        if not isinstance(step, dict):
            continue
        run_script = step.get("run")
        if not isinstance(run_script, str):
            continue
        for match in re.finditer(r"\bjust\s+([a-z0-9-]+)", run_script):
            cmd = f"just {match.group(1)}"
            if cmd not in commands:
                commands.append(cmd)
    assert commands, "Expected at least one 'just <recipe>' command in release.yaml"
    return commands


def validate_release_documentation(readme_text: str, workflow: dict[Any, Any]) -> None:
    """Validate that README.md's `## Releases` section matches `.github/workflows/release.yaml`."""
    releases_section = _extract_h2_section(readme_text, "Releases")
    normalized_section = " ".join(releases_section.lower().split())

    for forbidden in FORBIDDEN_RELEASE_PHRASES:
        assert forbidden.lower() not in normalized_section, (
            f"Forbidden or stale phrase {forbidden!r} found in '## Releases' section"
        )

    # Human approval boundary in the intro block (before first '### ' subheading)
    intro_block = releases_section.split("\n### ", 1)[0]
    normalized_intro = " ".join(intro_block.lower().split())
    assert (
        "explicit human approval is required before tagging or publishing a release"
        in normalized_intro
    ), (
        "'## Releases' intro must require explicit human approval before "
        "tagging or publishing a release"
    )
    assert "never tag or publish a release unless a human has approved it" in normalized_intro, (
        "'## Releases' intro must forbid tagging or publishing a release "
        "unless a human has approved it"
    )
    assert "always ask the user before downloading large release assets" in normalized_intro, (
        "'## Releases' intro must require asking the user before downloading large release assets"
    )

    # Trigger tags from release.yaml
    trigger_tags = _workflow_tag_triggers(workflow)
    for trigger_tag in trigger_tags:
        assert f"`{trigger_tag}`" in releases_section, (
            f"'## Releases' must document workflow tag trigger `{trigger_tag}`"
        )
    assert ".github/workflows/release.yaml" in releases_section, (
        "'## Releases' must reference `.github/workflows/release.yaml`"
    )

    # Web UI subsection
    web_ui_section = _extract_h3_section(releases_section, "Creating a release (GitHub Web UI)")
    assert "https://github.com/CultureBotAI/microbial-trait-profiling/releases" in web_ui_section, (
        "Web UI release instructions must link to CultureBotAI/microbial-trait-profiling/releases"
    )
    for required_ui_token in (
        "Draft a new release",
        "Create a new release",
        "Choose a tag",
        "vX.Y.Z",
        "Create new tag on publish",
        "`main`",
        "Generate release notes",
        "Publish release",
    ):
        assert required_ui_token in web_ui_section, (
            f"Web UI release instructions must include {required_ui_token!r}"
        )

    # CLI subsection
    cli_section = _extract_h3_section(releases_section, "Creating a release (CLI)")
    for required_cli_token in (
        "git checkout main && git pull",
        "just check",
        "git tag vX.Y.Z",
        "git push origin vX.Y.Z",
    ):
        assert required_cli_token in cli_section, (
            f"CLI release instructions must include {required_cli_token!r}"
        )

    # Pipeline & attached assets subsection
    assets_section = _extract_h3_section(releases_section, "Release pipeline and attached assets")
    expected_assets = _workflow_expected_asset_templates(workflow)
    documented_assets = set(re.findall(r"`(bsdbng-[a-z0-9-]+-<tag>\.[a-z0-9.]+)`", assets_section))
    assert documented_assets == expected_assets, (
        "Documented release assets in '### Release pipeline and attached assets' "
        f"({sorted(documented_assets)}) do not match .github/workflows/release.yaml "
        f"({sorted(expected_assets)})"
    )

    for just_cmd in _workflow_just_commands(workflow):
        assert f"`{just_cmd}`" in assets_section, (
            "'### Release pipeline and attached assets' must document "
            f"`{just_cmd}` from release.yaml"
        )

    # Post-release verification subsection
    verification_section = _extract_h3_section(releases_section, "Post-release verification")
    for required_verification_token in (
        "gh release list --repo CultureBotAI/microbial-trait-profiling --limit 1",
        "**Ask the user before downloading**",
        "gh release download <tag> --repo CultureBotAI/microbial-trait-profiling",
        "--pattern 'bsdbng-studies-*.tar.gz'",
        "--pattern 'bsdbng-stats-*.txt'",
        "--pattern 'bsdbng-ingest-log-*.json'",
        "--pattern 'bsdbng-raw-data-*.tar.gz'",
        "--pattern 'bsdbng-pipeline-log-*.log'",
        "taxon_<ID>",
        "https://bugsigdb.org/Main_Page",
    ):
        assert required_verification_token in verification_section, (
            f"Post-release verification section must include {required_verification_token!r}"
        )


def test_readme_release_section_matches_release_workflow() -> None:
    """Pin README.md's `## Releases` section against `.github/workflows/release.yaml` (#59)."""
    readme_text = README_MD.read_text(encoding="utf-8")
    workflow = _load_release_workflow(RELEASE_WORKFLOW_YAML)
    validate_release_documentation(readme_text, workflow)


def test_validate_release_documentation_rejects_missing_assets_and_negations() -> None:
    """Verify that drift in assets, missing steps, or negated approval fails validation."""
    readme_text = README_MD.read_text(encoding="utf-8")
    workflow = _load_release_workflow(RELEASE_WORKFLOW_YAML)

    # Dropping bsdbng-stats-<tag>.txt from the assets subsection must fail
    missing_stats = readme_text.replace("`bsdbng-stats-<tag>.txt`", "`bsdbng-other-<tag>.txt`")
    with pytest.raises(AssertionError, match="do not match"):
        validate_release_documentation(missing_stats, workflow)

    # Replacing the canonical repo with legacy turbomam/bsdbng must fail
    legacy_repo = readme_text.replace(
        "https://github.com/CultureBotAI/microbial-trait-profiling/releases",
        "https://github.com/turbomam/bsdbng/releases",
    )
    with pytest.raises(AssertionError, match="turbomam/bsdbng"):
        validate_release_documentation(legacy_repo, workflow)

    # Negating the human approval rule must fail
    negated_approval = readme_text.replace(
        "never tag or publish a release unless a human has\napproved it",
        "you may tag a release without user approval",
    )
    with pytest.raises(AssertionError, match="without user approval"):
        validate_release_documentation(negated_approval, workflow)
