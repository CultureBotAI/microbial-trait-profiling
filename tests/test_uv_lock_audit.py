"""Keep dependency drift distinct from lock serialization metadata (#1286)."""

from __future__ import annotations

import hashlib
import json
import runpy
from pathlib import Path
from typing import Any

import pytest

SCRIPT = Path(__file__).parents[1] / "scripts/adhoc/issue1286_compare_uv_locks.py"
HELPER = runpy.run_path(str(SCRIPT))
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
        assert result[role]["sha256"] == hashlib.sha256(raw).hexdigest()
        assert result[role]["size_bytes"] == len(raw)
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
