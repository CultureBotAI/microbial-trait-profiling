"""Evidence survives removal of original paths, without weakening verification (#726)."""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import pytest

from mtp.paths import REPO_ROOT, file_sha256
from mtp.workflow.__main__ import main
from mtp.workflow.archive import ArchiveResolver, create_archive, verify_archive
from mtp.workflow.consumers import open_core_bundle, write_consumer_record
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import execute_workflow, inventory


@pytest.fixture
def evidence(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    original = tmp_path / "original"
    inputs = original / "inputs"
    shutil.copytree(REPO_ROOT / "tests/fixtures/workflow", inputs)
    alias = original / "raw-alias.csv"
    alias.symlink_to(inputs / "full_dump.csv")
    run = execute_workflow(
        REPO_ROOT / "workflows/bugsigdb_core.yaml",
        {
            "raw_dump": alias,
            "kg_edges": inputs / "merged-kg_edges.tsv",
            "kg_nodes": inputs / "merged-kg_nodes.tsv",
        },
        original / "runs",
        OPERATIONS,
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    core = Path(run.manifest_path)
    result = original / "result.tsv"
    result.write_text("statistic\tvalue\nexample\t1\n")
    consumer = write_consumer_record(
        open_core_bundle(core),
        original / "consumer.json",
        consumer="archive-test",
        outputs=[result],
    )
    extra = original / ".supporting-evidence.txt"
    extra.write_text("Supplemental evidence is retained, not promoted to verified execution.\n")
    return original, core, consumer, extra


def archive(evidence: tuple[Path, Path, Path, Path], destination: Path) -> Path:
    _, core, consumer, extra = evidence
    return create_archive(
        destination, core_manifests=[core], consumer_records=[consumer], retained_paths=[extra]
    )


def test_relocated_archive_verifies_without_originals(
    evidence: tuple[Path, Path, Path, Path], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    original, core, consumer, extra = evidence
    original_bytes = {path: path.read_bytes() for path in (core, consumer, extra)}
    manifest = archive(evidence, tmp_path / "archive")
    before = json.loads(manifest.read_text())
    moved = tmp_path / "moved"
    manifest.parent.rename(moved)
    shutil.rmtree(original)
    retained = moved / "archive.json"
    record = verify_archive(retained)
    assert record.id == before["id"]
    resolver = ArchiveResolver(retained)
    for path, expected in original_bytes.items():
        assert resolver.resolve(path).read_bytes() == expected
    assert main(["verify-archive", str(retained)]) == 0
    assert "core contracts and consumer chains" in capsys.readouterr().out
    assert main(["archive-resolve", str(retained), str(extra)]) == 0
    assert str(resolver.resolve(extra)) in capsys.readouterr().out


@pytest.mark.parametrize("damage", ["corrupt", "missing", "extra", "symlink", "hardlink"])
def test_retained_members_are_not_replaced_by_live_inputs(
    evidence: tuple[Path, Path, Path, Path], tmp_path: Path, damage: str
) -> None:
    manifest = archive(evidence, tmp_path / "archive")
    resolver = ArchiveResolver(manifest)
    extra = resolver.resolve(evidence[3])
    if damage == "corrupt":
        extra.write_text("changed\n")
    elif damage == "missing":
        extra.unlink()
    elif damage == "extra":
        (manifest.parent / "data/undeclared.txt").write_text("extra\n")
    elif damage == "symlink":
        extra.unlink()
        extra.symlink_to(evidence[3])
    else:
        extra.unlink()
        os.link(evidence[3], extra)
    with pytest.raises(
        (ValueError, FileNotFoundError),
        match=r"changed|No such file|undeclared files|escapes archive|Symlink|Hard-linked",
    ):
        verify_archive(manifest)


@pytest.mark.parametrize(
    ("damage", "message"),
    [
        ("traversal", "Unsafe relative"),
        ("absolute", "Unsafe relative"),
        ("duplicate", "Duplicate archive source mapping"),
        ("conflict", "Conflicting archive source mappings"),
        ("unretained", "Path is not retained"),
    ],
)
def test_archive_mapping_rejects_ambiguous_or_unretained_locations(
    evidence: tuple[Path, Path, Path, Path], tmp_path: Path, damage: str, message: str
) -> None:
    manifest = archive(evidence, tmp_path / "archive")
    data = json.loads(manifest.read_text())
    first = data["locations"][0]
    if damage == "traversal":
        first["archive_path"] = "data/../escape"
    elif damage == "absolute":
        first["archive_path"] = "/absolute/path"
    elif damage == "duplicate":
        data["locations"].append(first.copy())
    elif damage == "conflict":
        second = data["locations"][1]
        second["original_paths"] = [first["original_paths"][0] + "/conflicting-child"]
    else:
        data["consumer_records"].append(str(tmp_path / "unretained.json"))
    manifest.write_text(json.dumps(data))
    with pytest.raises(ValueError, match=message):
        verify_archive(manifest)


def test_cli_infers_core_from_consumer_and_preserves_record_bytes(
    evidence: tuple[Path, Path, Path, Path], tmp_path: Path
) -> None:
    _, core, consumer, _ = evidence
    checksum = file_sha256(consumer)
    destination = tmp_path / "archive"
    # tmp_path is temporary storage, which the archive command refuses without the opt-in
    argv = ["archive", "--output", str(destination), "--consumer", str(consumer)]
    assert main([*argv, "--allow-temporary-output"]) == 0
    record = verify_archive(destination / "archive.json")
    assert record.core_manifests == [str(core)]
    assert file_sha256(consumer) == checksum
    with pytest.raises(FileExistsError, match="fresh directory"):
        archive(evidence, destination)


def test_archive_cannot_write_inside_retained_source(
    evidence: tuple[Path, Path, Path, Path],
) -> None:
    original, core, consumer, _ = evidence
    with pytest.raises(ValueError, match="overlaps"):
        create_archive(
            original / "archive",
            core_manifests=[core],
            consumer_records=[consumer],
            retained_paths=[original],
        )
    assert not (original / "archive").exists()


def test_declared_source_alias_is_retained_as_a_mapping(
    evidence: tuple[Path, Path, Path, Path], tmp_path: Path
) -> None:
    original, _, _, _ = evidence
    manifest = archive(evidence, tmp_path / "archive")
    resolver = ArchiveResolver(manifest)
    assert resolver.resolve(original / "raw-alias.csv") == resolver.resolve(
        original / "inputs/full_dump.csv"
    )


def test_retained_paths_do_not_create_an_execution_claim(tmp_path: Path) -> None:
    extra = tmp_path / "legacy.json"
    extra.write_text('{"status": "passed"}')
    with pytest.raises(ValueError, match="requires a core run"):
        create_archive(
            tmp_path / "archive", core_manifests=[], consumer_records=[], retained_paths=[extra]
        )


def test_source_identity_keeps_recorded_paths_even_after_consistent_rehash(
    evidence: tuple[Path, Path, Path, Path], tmp_path: Path
) -> None:
    manifest = archive(evidence, tmp_path / "archive")
    resolver = ArchiveResolver(manifest)
    core_path = resolver.resolve(evidence[1])
    run = json.loads(core_path.read_text())
    artifact = next(item for item in run["artifacts"] if item["product_id"] == "source_identity")
    identity_path = core_path.parent / artifact["path"]
    identity = json.loads(identity_path.read_text())
    raw = identity["consumed_files"]["raw_dump"]
    raw["path"] = str(resolver.resolve(Path(raw["path"])))
    identity_path.write_text(json.dumps(identity))
    digest, size, members = inventory(identity_path, "file")
    artifact.update(
        sha256=digest, size_bytes=size, inventory=[item.model_dump(mode="json") for item in members]
    )
    core_path.write_text(json.dumps(run))
    archived = json.loads(manifest.read_text())
    for location in archived["locations"]:
        path = manifest.parent / location["archive_path"]
        digest, size, members = inventory(path, location["storage_kind"])
        location.update(
            sha256=digest,
            size_bytes=size,
            inventory=[item.model_dump(mode="json") for item in members],
        )
    manifest.write_text(json.dumps(archived))
    with pytest.raises(ValueError, match="Source identity does not match consumed file: raw_dump"):
        verify_archive(manifest)
