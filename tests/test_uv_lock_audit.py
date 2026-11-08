"""Keep dependency drift distinct from lock serialization metadata (#1286)."""

from __future__ import annotations

import hashlib
import json
import runpy
from pathlib import Path
from typing import Any

import pytest

from mtp.local_path_patterns import compile_absolute_local_paths
from mtp.provenance_paths import assert_portable_text, resolve_reference, verify_reference

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts/adhoc/issue1286_compare_uv_locks.py"
GUARD_PROBE_SCRIPT = REPO_ROOT / "scripts/adhoc/issue928_guard_isolation.py"
HELPER = runpy.run_path(str(SCRIPT))
GUARD_PROBE = runpy.run_path(str(GUARD_PROBE_SCRIPT))
ABSOLUTE_LOCAL_PATH = compile_absolute_local_paths(flags=0, classified=False)
LOCK = """version = 1
revision = 1
requires-python = ">=3.11"
resolution-markers = ["python_version < '3.12'", "python_version >= '3.12'"]
[[package]]
name = "example"
version = "1.0"
source = { registry = "https://example.org/simple" }
dependencies = [{name = "other", marker = "python_version < '3.12'"}]
sdist = {url = "https://example.org/source.tar.gz", hash = "sha256:original", size = 12}
wheels = [{url = "https://example.org/package.whl", hash = "sha256:wheel", size = 9}]
[[package]]
name = "example"
version = "2.0"
source = { registry = "https://example.org/simple" }
"""


def compare(tmp_path: Path, before: str, after: str) -> dict[str, Any]:
    paths = [tmp_path / "before.lock", tmp_path / "after.lock"]
    for path, content in zip(paths, (before, after), strict=True):
        path.write_text(content)
    result: dict[str, Any] = HELPER["compare_locks"](*paths)
    assert paths[0].read_text() == before
    assert paths[1].read_text() == after
    for role, path in zip(("before", "after"), paths, strict=True):
        raw = path.read_bytes()
        assert result[role]["path"] == path.name
        assert not Path(result[role]["path"]).is_absolute()
        assert result[role]["sha256"] == hashlib.sha256(raw).hexdigest()
        assert result[role]["size_bytes"] == len(raw)
    assert result.get("path_basis") == {"kind": "root-relative-v1", "root": "."}
    assert_portable_text(json.dumps(result))
    return result


def test_only_root_revision_and_artifact_upload_times_are_metadata(tmp_path: Path) -> None:
    after = LOCK.replace("revision = 1", "revision = 5")
    after = after.replace("size = 12}", 'size = 12, upload-time = "2026-01-01T00:00:00Z"}')
    after = after.replace("size = 9}", 'size = 9, upload-time = "2026-01-02T00:00:00Z"}')
    result = compare(tmp_path, LOCK, after)
    assert result["substantive_equal"] is True
    assert result["normalized_diff"] == ""
    assert result["package_version_changes"] == {}
    assert result["format_metadata"]["after"] == [
        {"path": ["revision"], "value": 5},
        {"path": ["package", 0, "sdist", "upload-time"], "value": "2026-01-01T00:00:00Z"},
        {"path": ["package", 0, "wheels", 0, "upload-time"], "value": "2026-01-02T00:00:00Z"},
    ]


@pytest.mark.parametrize(
    ("old", "new"),
    [
        ("example.org/simple", "different.example/simple"),
        ("sha256:original", "sha256:changed"),
        ("size = 12", "size = 13"),
        ("python_version < '3.12'", "python_version < '3.13'"),
        ("version = 1\n", "version = 2\n"),
        ('requires-python = ">=3.11"', 'requires-python = ">=3.12"'),
    ],
)
def test_sources_hashes_sizes_markers_and_other_fields_remain_substantive(
    tmp_path: Path, old: str, new: str
) -> None:
    result = compare(tmp_path, LOCK, LOCK.replace(old, new))
    assert result["substantive_equal"] is False
    assert result["normalized_diff"].startswith("--- before.normalized.json\n")


def test_metadata_named_fields_elsewhere_are_not_discarded(tmp_path: Path) -> None:
    before = LOCK.replace("source = { registry =", 'source = { revision = "first", registry =')
    after = before.replace('revision = "first"', 'revision = "second"')
    result = compare(tmp_path, before, after)
    assert result["substantive_equal"] is False
    before = LOCK + '\n[custom]\nupload-time = "first"\n'
    result = compare(
        tmp_path, before, before.replace('upload-time = "first"', 'upload-time = "next"')
    )
    assert result["substantive_equal"] is False


def test_array_and_package_order_are_not_normalized_away(tmp_path: Path) -> None:
    markers = "[\"python_version < '3.12'\", \"python_version >= '3.12'\"]"
    reversed_markers = "[\"python_version >= '3.12'\", \"python_version < '3.12'\"]"
    assert (
        compare(tmp_path, LOCK, LOCK.replace(markers, reversed_markers))["substantive_equal"]
        is False
    )
    header, first, second = LOCK.split("[[package]]")
    result = compare(tmp_path, LOCK, header + "[[package]]" + second + "[[package]]" + first)
    assert result["substantive_equal"] is False


def test_multiple_versions_are_reported_without_collapsing_a_package(tmp_path: Path) -> None:
    result = compare(tmp_path, LOCK, LOCK.replace('version = "2.0"', 'version = "2.1"'))
    assert result["package_version_changes"] == {
        "example": {"before": ["1.0", "2.0"], "after": ["1.0", "2.1"]}
    }


def test_cli_writes_fresh_report_and_never_overwrites(tmp_path: Path) -> None:
    before, after, output = (
        tmp_path / name for name in ("before.lock", "after.lock", "audit.json")
    )
    before.write_text(LOCK)
    after.write_text(LOCK)
    argv = ["--before", str(before), "--after", str(after), "--output", str(output)]
    assert HELPER["main"](argv) == 0
    saved = output.read_bytes()
    assert json.loads(saved)["substantive_equal"] is True
    with pytest.raises(FileExistsError):
        HELPER["main"](argv)
    assert output.read_bytes() == saved
    alias = tmp_path / "dangling.json"
    alias.symlink_to(tmp_path / "missing.json")
    with pytest.raises(FileExistsError):
        HELPER["main"]([*argv[:-1], str(alias)])
    assert not (tmp_path / "missing.json").exists()


def test_compare_locks_uses_portable_provenance_and_refuses_account_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#1290: lock comparator emits portable paths and refuses account-bearing paths."""
    home = tmp_path / "Users" / "residentacct"
    workspace = home / "workspace"
    workspace.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", lambda: home)
    monkeypatch.setenv("USER", "residentacct")
    monkeypatch.setenv("LOGNAME", "residentacct")

    before = workspace / "before.lock"
    after = workspace / "sub" / "after.lock"
    after.parent.mkdir()
    before.write_text(LOCK, encoding="utf-8")
    after.write_text(LOCK, encoding="utf-8")
    output = workspace / "reports" / "audit.json"
    output.parent.mkdir()

    assert (
        HELPER["main"](["--before", str(before), "--after", str(after), "--output", str(output)])
        == 0
    )
    raw_report = output.read_text(encoding="utf-8")
    assert ABSOLUTE_LOCAL_PATH.search(raw_report) is None
    assert "residentacct" not in raw_report.casefold()
    report = json.loads(raw_report)
    assert report.get("path_basis") == {"kind": "root-relative-v1", "root": ".."}
    assert report["before"]["path"] == "before.lock"
    assert report["after"]["path"] == "sub/after.lock"
    assert_portable_text(raw_report)
    assert (
        verify_reference(output, report["before"], path_basis=report["path_basis"])
        == before.resolve()
    )
    assert (
        verify_reference(output, report["after"], path_basis=report["path_basis"])
        == after.resolve()
    )

    account_input = workspace / "residentacct-locks" / "before.lock"
    account_input.parent.mkdir()
    account_input.write_text(LOCK, encoding="utf-8")
    refused_output = workspace / "reports" / "refused.json"
    with pytest.raises(ValueError, match="Provenance contains a local account name"):
        HELPER["main"](
            [
                "--before",
                str(account_input),
                "--after",
                str(after),
                "--output",
                str(refused_output),
            ]
        )
    assert not refused_output.exists()

    with pytest.raises(ValueError, match=r"Provenance root must be below a home or temporary root"):
        HELPER["main"](
            [
                "--before",
                str(before),
                "--after",
                str(after),
                "--output",
                str(refused_output),
                "--provenance-root",
                str(home),
            ]
        )
    assert not refused_output.exists()


def test_guard_isolation_probe_emits_portable_paths_and_refuses_account_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#1027: site-stats guard probe emits portable paths and refuses account-bearing paths."""
    home = tmp_path / "Users" / "probeacct"
    workspace = home / "workspace"
    workspace.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", lambda: home)
    monkeypatch.setenv("USER", "probeacct")
    monkeypatch.setenv("LOGNAME", "probeacct")

    account_output = workspace / "probeacct-evidence"
    with pytest.raises(ValueError, match="Provenance contains a local account name"):
        GUARD_PROBE["main"](["--output-dir", str(account_output)])
    assert not account_output.exists()

    with pytest.raises(ValueError, match=r"Provenance root must be below a home or temporary root"):
        GUARD_PROBE["main"](
            ["--output-dir", str(home / "evidence"), "--provenance-root", str(home)]
        )
    assert not (home / "evidence").exists()

    probe_output = workspace / "guard-evidence"
    assert GUARD_PROBE["main"](["--output-dir", str(probe_output)]) == 0
    isolation_path = probe_output / "isolation.json"
    raw_isolation = isolation_path.read_text(encoding="utf-8")
    record = json.loads(raw_isolation)
    assert record["verified"] is True
    assert record.get("path_basis") == {"kind": "supplied-root-v1"}
    assert record["canonical_origin"] == "scripts/create_trait_profile_infographic.py"
    assert record["isolated_origin"] == "infographic_without_guard.py"
    assert record["checkout_outputs"] == "trait-visualization/trait-profile-investigation"
    assert not Path(record["canonical_origin"]).is_absolute()
    assert not Path(record["isolated_origin"]).is_absolute()
    assert ABSOLUTE_LOCAL_PATH.search(raw_isolation) is None
    assert "probeacct" not in raw_isolation.casefold()
    assert str(REPO_ROOT) not in raw_isolation
    assert str(probe_output) not in raw_isolation
    assert_portable_text(raw_isolation)
    assert set(record["source_before"]) == {
        "scripts/create_trait_profile_infographic.py",
        "tests/test_trait_profile_infographic.py",
        "scripts/adhoc/issue928_guard_isolation.py",
    }
    assert record["source_before"] == record["source_after"]
    assert (
        resolve_reference(
            isolation_path,
            record["canonical_origin"],
            path_basis=record["path_basis"],
            provenance_root=REPO_ROOT,
        )
        == (REPO_ROOT / "scripts/create_trait_profile_infographic.py").resolve()
    )
