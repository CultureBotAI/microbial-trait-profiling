"""Pin the capture source to the surviving Git blob, independent of PR refs (#1717)."""

from __future__ import annotations

import importlib.util
import io
import json
import re
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest
from git import Repo
from git.exc import BadName, InvalidGitRepositoryError, NoSuchPathError

ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "scripts/adhoc/capture_cellulose_numeric_snapshot_1709.py"
REVIEWED_BLOB = "0672ee7d"
DOCUMENTATION = ROOT / "methods/cellulose_decr.md"


def documented_blob() -> str:
    match = re.search(
        r"snapshot source blob is\s+`([0-9a-f]{40})`", DOCUMENTATION.read_text(encoding="utf-8")
    )
    assert match is not None
    return match.group(1)


@pytest.fixture(scope="module")
def capture_helper() -> ModuleType:
    spec = importlib.util.spec_from_file_location("cellulose_snapshot_capture", HELPER)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.__file__ is not None
    assert Path(module.__file__).resolve() == HELPER.resolve()
    return module


class TrackedStream(io.BytesIO):
    def __init__(self, content: bytes) -> None:
        super().__init__(content)
        self.read_sizes: list[int | None] = []

    def read(self, size: int | None = -1) -> bytes:
        self.read_sizes.append(size)
        return super().read(size)


class BlobRepository:
    def __init__(self, content: bytes, object_type: str, object_id: str) -> None:
        self.stream = TrackedStream(content)
        self.object_type = object_type
        self.object_id = object_id
        self.revisions: list[str] = []
        self.closed = False

    def rev_parse(self, revision: str) -> SimpleNamespace:
        self.revisions.append(revision)
        return SimpleNamespace(
            type=self.object_type, hexsha=self.object_id, data_stream=self.stream
        )

    def close(self) -> None:
        self.closed = True


def install_repository(
    capture_helper: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    root: Path,
    content: bytes,
    *,
    object_type: str = "blob",
    object_id: str | None = None,
) -> BlobRepository:
    repository = BlobRepository(content, object_type, object_id or documented_blob())

    def open_repository(path: Path) -> BlobRepository:
        assert path == root
        return repository

    monkeypatch.setattr(capture_helper, "ROOT", root)
    monkeypatch.setattr(capture_helper, "Repo", open_repository)
    return repository


def test_reads_the_full_source_blob_without_resolving_the_squashed_commit(
    capture_helper: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    assert capture_helper.REVIEWED_SOURCE_BLOB == "0672ee7d"
    source = b"HISTORICAL_VALUE = 1572\n"
    repository = install_repository(capture_helper, monkeypatch, tmp_path, source)
    assert capture_helper.read_reviewed_source() == source
    assert repository.revisions == ["0672ee7d"]
    assert repository.stream.read_sizes == [capture_helper.MAX_SOURCE_BYTES + 1]
    assert repository.closed


@pytest.mark.parametrize(
    ("object_type", "object_id"), [("tree", REVIEWED_BLOB), ("blob", "0" * 40)]
)
def test_refuses_any_object_other_than_the_pinned_source_blob(
    capture_helper: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    object_type: str,
    object_id: str,
) -> None:
    repository = install_repository(
        capture_helper,
        monkeypatch,
        tmp_path,
        b"VALUE = 1\n",
        object_type=object_type,
        object_id=object_id,
    )
    with pytest.raises(ValueError, match="pinned source blob"):
        capture_helper.read_reviewed_source()
    assert repository.stream.read_sizes == []
    assert repository.closed


def test_refuses_oversize_source_and_closes_the_repository(
    capture_helper: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    repository = install_repository(
        capture_helper, monkeypatch, tmp_path, b" " * (capture_helper.MAX_SOURCE_BYTES + 1)
    )
    with pytest.raises(ValueError, match="bounded size"):
        capture_helper.read_reviewed_source()
    assert repository.stream.read_sizes == [capture_helper.MAX_SOURCE_BYTES + 1]
    assert repository.closed


def test_reviewed_module_imports_only_the_pinned_payload_and_removes_its_copy(
    capture_helper: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    def pinned_source() -> bytes:
        return b"SOURCE_IDENTITY = 'reviewed-blob'\n"

    monkeypatch.setattr(capture_helper, "read_reviewed_source", pinned_source)
    with capture_helper.reviewed_module() as baseline:
        assert baseline.SOURCE_IDENTITY == "reviewed-blob"
        assert baseline.__file__ is not None
        temporary_source = Path(baseline.__file__)
        assert temporary_source.read_bytes() == pinned_source()
        assert temporary_source != ROOT / capture_helper.SOURCE
        assert sys.modules[baseline.__name__] is baseline
    assert not temporary_source.exists()
    assert baseline.__name__ not in sys.modules


@pytest.mark.parametrize("existing_kind", ["file", "dangling-symlink"])
def test_existing_output_is_refused_before_any_capture(
    capture_helper: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    existing_kind: str,
) -> None:
    output = tmp_path / "snapshot.json"
    if existing_kind == "file":
        output.write_bytes(b"keep this historical snapshot\n")
    else:
        output.symlink_to(tmp_path / "missing-target.json")

    def forbidden_capture() -> ModuleType:
        raise AssertionError("An existing output must be refused before loading producers")

    monkeypatch.setattr(sys, "argv", [str(HELPER), "--output", str(output)])
    monkeypatch.setattr(capture_helper, "load_extractors", forbidden_capture)
    with pytest.raises(ValueError, match=r"already exists|Aliased"):
        capture_helper.main()
    if existing_kind == "file":
        assert output.read_bytes() == b"keep this historical snapshot\n"
    else:
        assert output.is_symlink()
        assert not output.exists()


def test_main_captures_the_reviewed_module_not_the_current_producer(
    capture_helper: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    baseline = ModuleType("reviewed_baseline")
    extractors = ModuleType("numeric_extractors")
    current = ModuleType("current_producer_must_not_be_captured")
    monkeypatch.setattr(extractors, "tufte_plates", current, raising=False)
    snapshot = {"builders": {"plate": ["1572"]}, "svgs": {}, "slides": {}}

    @contextmanager
    def reviewed_module() -> Iterator[ModuleType]:
        yield baseline

    def capture(module: ModuleType, helpers: ModuleType) -> dict[str, dict[str, list[str]]]:
        assert module is baseline
        assert helpers is extractors
        return snapshot

    output = tmp_path / "fresh-snapshot.json"
    monkeypatch.setattr(sys, "argv", [str(HELPER), "--output", str(output)])
    monkeypatch.setattr(capture_helper, "load_extractors", lambda: extractors)
    monkeypatch.setattr(capture_helper, "reviewed_module", reviewed_module)
    monkeypatch.setattr(capture_helper, "capture", capture)
    assert capture_helper.main() == 0
    payload = json.loads(output.read_bytes())
    assert payload["baseline_commit"] == "92539eb"
    assert payload["numeric_tokens"] == snapshot


def test_real_git_blob_matches_documented_full_identity(capture_helper: ModuleType) -> None:
    """Missing/ambiguous object names in a full clone are failures, never skips (#1730)."""
    try:
        repository = Repo(ROOT)
    except (InvalidGitRepositoryError, NoSuchPathError):
        pytest.skip("Historical Git-object check requires a Git checkout")
    try:
        if repository.git.rev_parse("--is-shallow-repository").strip() == "true":
            pytest.skip("Historical Git-object check requires full history")
        assert capture_helper.REVIEWED_SOURCE_BLOB == "0672ee7d"
        blob = repository.rev_parse("0672ee7d")
        assert blob.type == "blob"
        assert blob.hexsha.startswith("0672ee7d")
        assert blob.hexsha == documented_blob()
        assert f"`{blob.hexsha}`" in DOCUMENTATION.read_text(encoding="utf-8")
        assert capture_helper.read_reviewed_source() == blob.data_stream.read(
            capture_helper.MAX_SOURCE_BYTES + 1
        )
    finally:
        repository.close()


def test_real_git_check_propagates_bad_name_instead_of_skipping(
    capture_helper: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    def missing_blob(revision: str) -> None:
        assert revision == "0672ee7d"
        raise BadName(revision)

    repository = SimpleNamespace(
        git=SimpleNamespace(rev_parse=lambda _flag: "false"),
        rev_parse=missing_blob,
        close=lambda: None,
    )
    monkeypatch.setattr(sys.modules[__name__], "Repo", lambda _root: repository)
    # This synthetic full clone must fail if any BadName handler tries to skip.
    monkeypatch.setattr(pytest, "skip", pytest.fail)
    with pytest.raises(BadName):
        test_real_git_blob_matches_documented_full_identity(capture_helper)
