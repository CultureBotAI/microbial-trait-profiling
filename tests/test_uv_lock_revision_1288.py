"""Regression guard for Issue #1288: `uv.lock` revision and `upload-time` metadata."""

from __future__ import annotations

import importlib.util
import sys
import tomllib
from pathlib import Path
from typing import Any, Protocol, cast

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


class _VerifierModule(Protocol):
    def validate_scratch_dir(self, scratch_dir: Path, repo_root: Path) -> Path: ...

    def evaluate_mutant_outcome(
        self,
        exit_code: int,
        inspector: Any,
        junit_xml: Path,
        output: str,
    ) -> tuple[bool, str]: ...

    def _run_targeted_pytest(
        self,
        scratch_repo: Path,
        nodeid: str,
        junit_xml: Path,
        *,
        repo_root: Path | None = None,
    ) -> tuple[bool, str]: ...


def _load_verifier_module() -> _VerifierModule:
    script_path = REPO_ROOT / "scripts/adhoc/verify_1288_mutants_agy_20261111.py"
    spec = importlib.util.spec_from_file_location("verify_1288_mutants_agy_20261111", script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return cast(_VerifierModule, module)


def _has_nonempty_upload_time(artifact: object) -> bool:
    if not isinstance(artifact, dict):
        return False
    upload_time = artifact.get("upload-time")
    return isinstance(upload_time, str) and bool(upload_time.strip())


def validate_uv_lock_format(lock_path: Path) -> dict[str, int]:
    """Parse `uv.lock` and assert `revision >= 5` plus `upload-time` on all registry artifacts."""
    data = tomllib.loads(lock_path.read_text(encoding="utf-8"))
    version = data.get("version")
    if isinstance(version, bool) or not isinstance(version, int) or version != 1:
        raise AssertionError(f"Expected uv.lock version == 1, got {version!r}")

    revision = data.get("revision")
    if isinstance(revision, bool) or not isinstance(revision, int) or revision < 5:
        raise AssertionError(
            f"Expected uv.lock revision >= 5 (uv >= 0.13.4 lock format), got {revision!r}"
        )

    packages = data.get("package")
    if not isinstance(packages, list) or not packages:
        raise AssertionError("Expected non-empty [[package]] list in uv.lock")

    registry_packages = 0
    checked_artifacts = 0
    for pkg in packages:
        if not isinstance(pkg, dict):
            raise AssertionError(f"Expected package table dict, got {type(pkg).__name__}")
        pkg_name = str(pkg.get("name") or "<unnamed>")
        pkg_version = str(pkg.get("version") or "<unversioned>")
        source = pkg.get("source")
        if not isinstance(source, dict) or "registry" not in source:
            continue

        registry_packages += 1
        sdist = pkg.get("sdist")
        wheels = pkg.get("wheels")
        pkg_artifacts = 0

        if sdist is not None:
            if not _has_nonempty_upload_time(sdist):
                raise AssertionError(
                    f"Registry package {pkg_name}=={pkg_version} sdist is missing a "
                    "non-empty 'upload-time'"
                )
            pkg_artifacts += 1
            checked_artifacts += 1

        if wheels is not None:
            if not isinstance(wheels, list) or not wheels:
                raise AssertionError(
                    f"Registry package {pkg_name}=={pkg_version} has invalid wheels={wheels!r}"
                )
            for idx, wheel in enumerate(wheels):
                if not _has_nonempty_upload_time(wheel):
                    wheel_id = (
                        wheel.get("url", f"wheels[{idx}]")
                        if isinstance(wheel, dict)
                        else f"wheels[{idx}]"
                    )
                    raise AssertionError(
                        f"Registry package {pkg_name}=={pkg_version} wheel ({wheel_id}) is "
                        "missing a non-empty 'upload-time'"
                    )
                pkg_artifacts += 1
                checked_artifacts += 1

        if pkg_artifacts == 0:
            raise AssertionError(
                f"Registry package {pkg_name}=={pkg_version} has neither sdist nor wheels"
            )

    if checked_artifacts == 0:
        raise AssertionError("Expected at least one registry distribution artifact in uv.lock")

    return {
        "version": version,
        "revision": revision,
        "registry_packages": registry_packages,
        "checked_artifacts": checked_artifacts,
        "missing_upload_time": 0,
    }


def test_committed_uv_lock_has_revision_5_and_upload_times() -> None:
    """Committed uv.lock keeps revision >= 5 and upload-time on every registry artifact."""
    summary = validate_uv_lock_format(REPO_ROOT / "uv.lock")
    assert summary["version"] == 1
    assert summary["revision"] == 5
    assert summary["registry_packages"] >= 200
    # A dependency update may legitimately add or remove distribution artifacts.
    # Derive the expected count from the lock, while validating every upload-time.
    lock_data = tomllib.loads((REPO_ROOT / "uv.lock").read_text(encoding="utf-8"))
    expected_artifacts = sum(
        int("sdist" in package) + len(package.get("wheels", []))
        for package in lock_data["package"]
        if "registry" in package.get("source", {})
    )
    assert summary["checked_artifacts"] == expected_artifacts
    assert summary["missing_upload_time"] == 0


def test_validate_uv_lock_format_rejects_regressed_or_missing_revision(tmp_path: Path) -> None:
    """revision = 1, revision = 4, or omitted revision fail with AssertionError."""
    valid_body = (
        "[[package]]\n"
        'name = "alabaster"\n'
        'version = "1.0.0"\n'
        'source = { registry = "https://pypi.org/simple" }\n'
        'sdist = { url = "https://example.org/a.tar.gz", upload-time = "2024-07-26T18:15:03Z" }\n'
    )
    for idx, header in enumerate(
        (
            "version = 1\nrevision = 1\n",
            "version = 1\nrevision = 4\n",
            "version = 1\n",
        )
    ):
        lock_file = tmp_path / f"uv_rev_{idx}.lock"
        lock_file.write_text(header + valid_body, encoding="utf-8")
        caught: AssertionError | None = None
        try:
            validate_uv_lock_format(lock_file)
        except AssertionError as exc:
            caught = exc
        assert caught is not None
        assert "revision >= 5" in str(caught)


def test_validate_uv_lock_format_rejects_missing_or_blank_sdist_upload_time(
    tmp_path: Path,
) -> None:
    """A registry package sdist without non-empty upload-time fails with AssertionError."""
    for idx, sdist_line in enumerate(
        (
            'sdist = { url = "https://example.org/alabaster-1.0.0.tar.gz", size = 24210 }\n',
            'sdist = { url = "https://example.org/alabaster-1.0.0.tar.gz", upload-time = "  " }\n',
        )
    ):
        lock_file = tmp_path / f"uv_sdist_{idx}.lock"
        lock_file.write_text(
            "version = 1\n"
            "revision = 5\n"
            "[[package]]\n"
            'name = "alabaster"\n'
            'version = "1.0.0"\n'
            'source = { registry = "https://pypi.org/simple" }\n' + sdist_line,
            encoding="utf-8",
        )
        caught: AssertionError | None = None
        try:
            validate_uv_lock_format(lock_file)
        except AssertionError as exc:
            caught = exc
        assert caught is not None
        assert "alabaster==1.0.0" in str(caught)
        assert "sdist" in str(caught)
        assert "upload-time" in str(caught)


def test_validate_uv_lock_format_rejects_missing_or_blank_wheel_upload_time(
    tmp_path: Path,
) -> None:
    """A registry package wheel without non-empty upload-time fails with AssertionError."""
    for idx, wheel_line in enumerate(
        (
            'wheels = [{ url = "https://example.org/alabaster.whl", size = 13929 }]\n',
            'wheels = [{ url = "https://example.org/alabaster.whl", upload-time = "" }]\n',
        )
    ):
        lock_file = tmp_path / f"uv_wheel_{idx}.lock"
        lock_file.write_text(
            "version = 1\n"
            "revision = 5\n"
            "[[package]]\n"
            'name = "alabaster"\n'
            'version = "1.0.0"\n'
            'source = { registry = "https://pypi.org/simple" }\n' + wheel_line,
            encoding="utf-8",
        )
        caught: AssertionError | None = None
        try:
            validate_uv_lock_format(lock_file)
        except AssertionError as exc:
            caught = exc
        assert caught is not None
        assert "alabaster==1.0.0" in str(caught)
        assert "wheel" in str(caught)
        assert "upload-time" in str(caught)


def test_validate_uv_lock_format_allows_virtual_package_and_rejects_zero_registry_artifacts(
    tmp_path: Path,
) -> None:
    """Non-registry packages are allowed alongside registry ones, but 0 registry artifacts fail."""
    mixed_lock = tmp_path / "uv_mixed.lock"
    mixed_lock.write_text(
        "version = 1\n"
        "revision = 5\n"
        "[[package]]\n"
        'name = "mtp"\n'
        'version = "0.1.0"\n'
        'source = { editable = "." }\n'
        "[[package]]\n"
        'name = "alabaster"\n'
        'version = "1.0.0"\n'
        'source = { registry = "https://pypi.org/simple" }\n'
        'sdist = { url = "https://example.org/a.tar.gz", upload-time = "2024-07-26T18:15:03Z" }\n',
        encoding="utf-8",
    )
    summary = validate_uv_lock_format(mixed_lock)
    assert summary["registry_packages"] == 1
    assert summary["checked_artifacts"] == 1

    virtual_only_lock = tmp_path / "uv_virtual_only.lock"
    virtual_only_lock.write_text(
        "version = 1\n"
        "revision = 5\n"
        "[[package]]\n"
        'name = "mtp"\n'
        'version = "0.1.0"\n'
        'source = { virtual = "." }\n',
        encoding="utf-8",
    )
    caught_zero: AssertionError | None = None
    try:
        validate_uv_lock_format(virtual_only_lock)
    except AssertionError as exc:
        caught_zero = exc
    assert caught_zero is not None
    assert "at least one registry distribution artifact" in str(caught_zero)


def test_verifier_rejects_scratch_dir_inside_repo_before_mkdir(tmp_path: Path) -> None:
    """validate_scratch_dir rejects repo root, child path, and inbound symlink before mkdir."""
    verifier = _load_verifier_module()
    fake_repo = tmp_path / "fake_repo"
    fake_repo.mkdir()

    inside_child = fake_repo / "nested_scratch_that_must_not_be_created"
    caught_child: ValueError | None = None
    try:
        verifier.validate_scratch_dir(inside_child, fake_repo)
    except ValueError as exc:
        caught_child = exc
    assert caught_child is not None
    assert not inside_child.exists()

    inbound_symlink = tmp_path / "symlink_into_repo"
    inbound_symlink.symlink_to(fake_repo / "inside_target")
    caught_symlink: ValueError | None = None
    try:
        verifier.validate_scratch_dir(inbound_symlink, fake_repo)
    except ValueError as exc:
        caught_symlink = exc
    assert caught_symlink is not None
    assert not (fake_repo / "inside_target").exists()


def test_verifier_rejects_setup_error_and_internal_exit_code(tmp_path: Path) -> None:
    """Verifier rejects both setup-phase ValueError('AssertionError') and exit_code == 3."""
    verifier = _load_verifier_module()
    scratch_repo = tmp_path / "repo"
    tests_dir = scratch_repo / "tests"
    tests_dir.mkdir(parents=True)
    test_file = tests_dir / "test_setup_boom.py"
    test_file.write_text(
        "import pytest\n\n"
        "@pytest.fixture\n"
        "def bad_fixture() -> None:\n"
        '    raise ValueError("AssertionError in setup")\n\n'
        "def test_boom(bad_fixture: None) -> None:\n"
        "    assert True\n",
        encoding="utf-8",
    )
    junit_xml = tmp_path / "junit_setup_boom.xml"
    killed, detail = verifier._run_targeted_pytest(
        scratch_repo, "tests/test_setup_boom.py::test_boom", junit_xml
    )
    assert killed is False
    assert "setup" in detail
    assert "ValueError" in detail

    fake_junit = tmp_path / "junit_fake_call_fail.xml"
    fake_junit.write_text(
        '<testsuite tests="1"><testcase classname="t" name="n">'
        '<failure message="AssertionError">AssertionError</failure>'
        "</testcase></testsuite>",
        encoding="utf-8",
    )

    class _FakeReport:
        failed = True

    class _FakeInspector:
        def __init__(self) -> None:
            self.call_reports = [_FakeReport()]
            self.non_call_failures: list[tuple[str, str]] = []
            self.call_exc_types = ["AssertionError"]

    killed_internal, detail_internal = verifier.evaluate_mutant_outcome(
        int(pytest.ExitCode.INTERNAL_ERROR),
        _FakeInspector(),
        fake_junit,
        "INTERNALERROR>",
    )
    assert killed_internal is False
    assert "exit_code=3" in detail_internal
