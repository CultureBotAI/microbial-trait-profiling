"""Retain and independently verify immutable execution evidence (#726)."""

from __future__ import annotations

import shutil
import uuid
from datetime import UTC, datetime
from pathlib import Path

from mtp.datamodel.workflow_pydantic import (
    ArchivedLocation,
    ConsumerRecord,
    EvidenceArchive,
    StorageKind,
)
from mtp.workflow.consumers import verify_consumer_record
from mtp.workflow.definition import safe_relative_path, unique
from mtp.workflow.runtime import inventory, validate_run
from mtp.workflow.verification import verify_core_run


def _absolute(path: Path) -> Path:
    if not path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Expected normalized absolute recorded path: {path}")
    return path


def _relative(value: str) -> Path:
    safe_relative_path(value)
    return Path(value)


class ArchiveResolver:
    """Resolve exclusively through a validated archive; never consult original inputs."""

    def __init__(self, manifest: Path) -> None:
        if manifest.is_symlink():
            raise ValueError("Archive manifest cannot be a symlink")
        self.root = manifest.resolve(strict=True).parent
        self.record = EvidenceArchive.model_validate_json(manifest.read_text(encoding="utf-8"))
        if not self.record.core_manifests or not self.record.locations:
            raise ValueError("Archive requires at least one core and retained location")
        self.mappings: dict[Path, tuple[Path, bool]] = {}
        self.files: set[Path] = set()
        for entry in self.record.locations:
            relative = _relative(entry.archive_path)
            if relative.parts[0] != "data":
                raise ValueError("Retained archive locations must be inside data/")
            target = self.root / relative
            if not entry.original_paths:
                raise ValueError("Archived location lacks original paths")
            for original in entry.original_paths:
                source = _absolute(Path(original))
                if source in self.mappings:
                    raise ValueError(f"Duplicate archive source mapping: {source}")
                self.mappings[source] = (target, entry.storage_kind == "collection")
            if entry.storage_kind == "file":
                self.files.add(target)
            else:
                self.files.update(target / _relative(item.path) for item in entry.inventory)
            self._within_archive(target)
            actual = inventory(target, str(entry.storage_kind))
            if actual != (entry.sha256, entry.size_bytes, entry.inventory):
                raise ValueError(f"Archived content or inventory changed: {entry.archive_path}")
        for source, (target, _) in self.mappings.items():
            for parent in source.parents:
                if parent in self.mappings:
                    parent_target, collection = self.mappings[parent]
                    if not collection or target != parent_target / source.relative_to(parent):
                        raise ValueError(f"Conflicting archive source mappings: {source}")
        observed = set()
        for member in (self.root / "data").rglob("*"):
            self._within_archive(member)
            if member.is_file():
                if member.stat().st_nlink != 1:
                    raise ValueError(f"Hard-linked retained file is not independent: {member}")
                observed.add(member)
            elif not member.is_dir():
                raise ValueError(f"Non-regular archive member: {member}")
        if observed != self.files:
            raise ValueError("Archive contains missing or undeclared files")
        for label, values in (
            ("core manifests", self.record.core_manifests),
            ("consumer records", self.record.consumer_records or []),
            ("retained paths", self.record.retained_paths or []),
        ):
            unique(values, label)
            for value in values:
                self.resolve(_absolute(Path(value)))

    def _within_archive(self, path: Path) -> None:
        if not path.resolve(strict=True).is_relative_to(self.root / "data"):
            raise ValueError(f"Retained location escapes archive data directory: {path}")
        if any(
            item.is_symlink() for item in (path, *path.parents) if item.is_relative_to(self.root)
        ):
            raise ValueError(f"Symlink inside evidence archive: {path}")

    def resolve(self, path: Path) -> Path:
        """Map a historical absolute path or check an already-resolved archive member."""
        path = _absolute(path)
        if path.is_relative_to(self.root / "data"):
            target = path
        else:
            target = None
            for original in (path, *path.parents):
                if original in self.mappings:
                    retained, collection = self.mappings[original]
                    if original != path and not collection:
                        raise ValueError(f"Cannot resolve a child of archived file: {path}")
                    target = retained / path.relative_to(original)
                    break
        if target is None:
            raise ValueError(f"Path is not retained in this archive: {path}")
        self._within_archive(target)
        if target.is_file() and target not in self.files:
            raise ValueError(f"Uninventoried archive member: {target}")
        return target


def verify_archive(manifest: Path) -> EvidenceArchive:
    """Revalidate recorded core science and consumer chains using retained bytes only."""
    resolver = ArchiveResolver(manifest)
    for core in resolver.record.core_manifests:
        verify_core_run(Path(core), resolver=resolver.resolve)
    for consumer in resolver.record.consumer_records or []:
        record = verify_consumer_record(Path(consumer), resolver=resolver.resolve)
        if record.run_manifest_path not in resolver.record.core_manifests:
            raise ValueError("Archive consumer references an undeclared core manifest")
    return resolver.record


def create_archive(
    output: Path,
    *,
    core_manifests: list[Path],
    consumer_records: list[Path],
    retained_paths: list[Path],
) -> Path:
    """Copy declared dependencies into a fresh, verified archive without editing originals."""
    output = output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Archive output must be a fresh directory: {output}")
    sources: dict[Path, Path] = {}

    def retain(path: Path) -> Path:
        original = _absolute(path.absolute())
        source = original.resolve(strict=True)
        if not (source.is_file() or source.is_dir()):
            raise ValueError(f"Cannot retain non-regular input: {source}")
        sources[original] = source
        sources[source] = source
        return original

    cores = {retain(path) for path in core_manifests}
    consumers = {retain(path) for path in consumer_records}
    extras = {retain(path) for path in retained_paths}
    for path in consumers:
        record = ConsumerRecord.model_validate_json(path.read_text(encoding="utf-8"))
        cores.add(retain(Path(record.run_manifest_path)))
        for member in record.outputs:
            retain(Path(member.path))
    if not cores:
        raise ValueError("An evidence archive requires a core run or a core-bound consumer")
    for manifest in cores:
        run = validate_run(manifest)
        retain(manifest.parent)
        for product in run.definition.products:
            if product.external:
                artifact = next(
                    item for item in run.artifacts or [] if item.product_id == product.id
                )
                retain(Path(artifact.path))
                if artifact.requested_path:
                    retain(Path(artifact.requested_path))
    roots: list[Path] = []
    for source in sorted(set(sources.values()), key=lambda path: (len(path.parts), str(path))):
        if not any(source.is_relative_to(root) for root in roots):
            roots.append(source)
    destination = output.resolve()
    for source in roots:
        if destination.is_relative_to(source) or source.is_relative_to(destination):
            raise ValueError("Archive destination overlaps retained source paths")
    output.mkdir(parents=True)
    staging = output / ".staging"
    staging.mkdir()
    try:
        locations = []
        destinations = {}
        for index, source in enumerate(roots):
            kind = StorageKind.collection if source.is_dir() else StorageKind.file
            before = inventory(source, kind)
            relative = Path("data") / f"{index:04d}" / source.name
            target = staging / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            if kind == "collection":
                shutil.copytree(source, target)
            else:
                shutil.copyfile(source, target)
            if inventory(target, kind) != before or inventory(source, kind) != before:
                raise ValueError(f"Source changed during archival copy: {source}")
            destinations[source] = relative
            locations.append(
                ArchivedLocation(
                    original_paths=[str(source)],
                    archive_path=relative.as_posix(),
                    storage_kind=kind,
                    sha256=before[0],
                    size_bytes=before[1],
                    inventory=before[2],
                )
            )
        for original, source in sorted(sources.items()):
            if original == source:
                continue
            root = next(root for root in roots if source.is_relative_to(root))
            relative = destinations[root] / source.relative_to(root)
            matched = next(
                (entry for entry in locations if entry.archive_path == relative.as_posix()), None
            )
            if matched is not None:
                matched.original_paths.append(str(original))
            else:
                kind = StorageKind.collection if source.is_dir() else StorageKind.file
                digest, size, members = inventory(staging / relative, kind)
                locations.append(
                    ArchivedLocation(
                        original_paths=[str(original)],
                        archive_path=relative.as_posix(),
                        storage_kind=kind,
                        sha256=digest,
                        size_bytes=size,
                        inventory=members,
                    )
                )
        archive_record = EvidenceArchive(
            id=f"archive:{uuid.uuid4()}",
            archive_version="1",
            created_at=datetime.now(UTC),
            core_manifests=sorted(map(str, cores)),
            consumer_records=sorted(map(str, consumers)),
            retained_paths=sorted(map(str, extras)),
            locations=locations,
        )
        manifest = staging / "archive.json"
        manifest.write_text(archive_record.model_dump_json(indent=2) + "\n", encoding="utf-8")
        verify_archive(manifest)
        (staging / "data").rename(output / "data")
        manifest.rename(output / "archive.json")
        staging.rmdir()
        return output / "archive.json"
    except BaseException:
        shutil.rmtree(output)
        raise
