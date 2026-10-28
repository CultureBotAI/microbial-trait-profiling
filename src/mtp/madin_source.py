"""The pinned Madin et al. condensed trait table read by the label-support audit (#1048).

The audit resolves ``infores:madin_etal`` KG annotations to the datasets Madin et al.
compiled. That needs ``condensed_traits_NCBI.csv``: one row per source record, with its
NCBI ``tax_id`` and ``data_source``. ``condensed_species_NCBI.csv`` merges records per
species and drops both columns, so it cannot resolve sub-sources.

The table is fetched by COMMIT SHA, as the BugSigDB dump is (#151): a branch or tag URL
could serve different bytes later. The repository's only release, v1.0.0, and the
``jmadin`` fork both carry an earlier table, so neither can stand in for this pin.
``methods/label_support_audit.md`` records how the pin was confirmed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import stat
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from uuid import uuid4

import httpx

from mtp.paths import REPO_ROOT, file_sha256
from mtp.workflow.output_roots import require_durable_output, temporary_root

PROVENANCE_FILE = "download_provenance.json"


@dataclass(frozen=True)
class SourcePin:
    """A file fixed to one commit of a public GitHub repository, with its size and hash."""

    repository: str
    commit: str
    commit_date: str
    path: str
    sha256: str
    size_bytes: int

    @property
    def filename(self) -> str:
        return Path(self.path).name

    @property
    def url(self) -> str:
        return f"https://raw.githubusercontent.com/{self.repository}/{self.commit}/{self.path}"

    def descriptor(self) -> dict[str, object]:
        return {**asdict(self), "url": self.url}


# Upstream master has pointed at this commit since 2020-10-16.
MADIN_PIN = SourcePin(
    repository="bacteria-archaea-traits/bacteria-archaea-traits",
    commit="396dcd4d98ca20b613ba280f5514611651b6fb91",
    commit_date="2020-10-16",
    path="output/condensed_traits_NCBI.csv",
    sha256="6100fc11679cd9facae65cad37a2969978bc64c3e57984d896f8147f78db6440",
    size_bytes=46362149,
)


class MadinPinError(ValueError):
    """Raised when a Madin table's bytes differ from the pinned release."""


def _check_size(label: str, size: int) -> None:
    pin = MADIN_PIN
    where = f"{pin.path} at {pin.repository}@{pin.commit[:7]}"
    if size != pin.size_bytes:
        msg = f"{label} is {size:,} bytes, expected {pin.size_bytes:,} for {where}"
        raise MadinPinError(msg)


def _check(label: str, size: int, digest: str) -> None:
    _check_size(label, size)
    pin = MADIN_PIN
    where = f"{pin.path} at {pin.repository}@{pin.commit[:7]}"
    if digest != pin.sha256:
        msg = f"{label} has sha256 {digest}, expected {pin.sha256} for {where}"
        raise MadinPinError(msg)


def verify_madin_table(path: Path) -> dict[str, object]:
    """Return the pin's descriptor if *path* has its size and sha256; refuse it otherwise."""
    information = path.stat()
    if not stat.S_ISREG(information.st_mode):
        raise MadinPinError(f"{path} is not a regular file (a directory cannot be a Madin table)")
    _check_size(path.name, information.st_size)
    _check(path.name, information.st_size, file_sha256(path))
    return MADIN_PIN.descriptor()


def verify_madin_record(label: str, size_bytes: int, sha256: str) -> dict[str, object]:
    """Return the pin's descriptor if a recorded size and sha256 are the pin's; refuse them."""
    _check(label, size_bytes, sha256)
    return MADIN_PIN.descriptor()


def _external_durable_path(path: Path) -> Path:
    """Apply main's temporary-root rule and refuse this working tree before any write."""
    expanded = path.expanduser()
    lexical = Path(os.path.abspath(expanded))
    resolved = expanded.resolve()
    roots = {REPO_ROOT.absolute(), REPO_ROOT.resolve()}
    for spelling in (lexical, resolved):
        if temporary := temporary_root(spelling):
            raise MadinPinError(
                f"{path} is under temporary storage {temporary}; use a durable data directory"
            )
        if any(spelling.is_relative_to(root) for root in roots):
            raise MadinPinError(f"{path} is inside the repository working tree")
    return resolved


def _read_provenance(path: Path) -> dict[str, dict[str, object]]:
    """Refuse malformed or contradictory records before network access or placement."""
    if not path.exists():
        return {}
    if not path.is_file():
        raise MadinPinError(f"Madin provenance record is not a file: {path}")
    try:
        decoded = json.loads(path.read_bytes())
        # JSON's non-finite extensions cannot be preserved in a standards-compliant record.
        json.dumps(decoded, allow_nan=False)
    except (ValueError, UnicodeError) as error:
        raise MadinPinError(f"Invalid Madin provenance JSON: {path}") from error
    if not isinstance(decoded, dict) or any(
        not isinstance(entry, dict) for entry in decoded.values()
    ):
        raise MadinPinError("Madin provenance must be an object of record objects")
    records: dict[str, dict[str, object]] = decoded
    selected = records.get(MADIN_PIN.filename, {})
    if any(
        key in selected and selected[key] != value for key, value in MADIN_PIN.descriptor().items()
    ):
        raise MadinPinError("Existing Madin provenance record contradicts the source pin")
    return records


def _write_provenance(path: Path, records: dict[str, dict[str, object]], *, reused: bool) -> None:
    """Replace a complete record atomically, retaining known download time on reuse."""
    now = datetime.now(tz=UTC).isoformat()
    entry = {
        **records.get(MADIN_PIN.filename, {}),
        **MADIN_PIN.descriptor(),
        "verified_at": now,
        "retrieval": "verified_reuse" if reused else "downloaded",
        "integrity_verified": True,
    }
    if not reused:
        entry["downloaded_at"] = now
    updated = {**records, MADIN_PIN.filename: entry}
    raw = json.dumps(updated, indent=2, allow_nan=False) + "\n"
    partial = path.with_name(f".{path.name}.{uuid4().hex}.part")
    owned = False
    try:
        with partial.open("x", encoding="utf-8") as handle:
            owned = True
            handle.write(raw)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(partial, path)
    finally:
        if owned:
            partial.unlink(missing_ok=True)


def download_madin_table(dest_dir: Path) -> Path:
    """Fetch the pinned table into *dest_dir*, anonymously; never keep bytes that fail the pin.

    An existing copy is checked and reused, never overwritten. Each successful reuse
    records verification without inventing a download time. The request carries no
    credentials or identifying headers. Metadata publication is atomic; if it fails
    after a new verified table is placed, that table may remain for verified reuse.
    """
    pin = MADIN_PIN
    dest_dir = _external_durable_path(dest_dir)
    target = dest_dir / pin.filename
    _external_durable_path(target)
    provenance_path = dest_dir / PROVENANCE_FILE
    provenance = _read_provenance(provenance_path)
    if target.exists() or target.is_symlink():
        try:
            verify_madin_table(target)
        except MadinPinError as error:
            msg = f"{error}; remove it to fetch the pinned file again"
            raise MadinPinError(msg) from error
        _write_provenance(provenance_path, provenance, reused=True)
        return target
    dest_dir.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(f".{target.name}.{uuid4().hex}.part")
    digest = hashlib.sha256()
    size = 0
    owned = False
    try:
        with partial.open("xb") as handle:
            owned = True
            with (
                httpx.Client(follow_redirects=True, timeout=120) as client,
                client.stream("GET", pin.url) as response,
            ):
                response.raise_for_status()
                for block in response.iter_bytes():
                    next_size = size + len(block)
                    if next_size > pin.size_bytes:
                        _check_size(pin.url, next_size)
                    handle.write(block)
                    digest.update(block)
                    size = next_size
        _check(pin.url, size, digest.hexdigest())
        try:
            os.link(partial, target)
        except FileExistsError as error:
            raise MadinPinError(f"Madin destination appeared during download: {target}") from error
    finally:
        if owned:
            partial.unlink(missing_ok=True)
    _write_provenance(provenance_path, provenance, reused=False)
    return target


def main(argv: list[str] | None = None) -> int:
    """Download the pinned Madin et al. condensed trait table and check its sha256."""
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument(
        "--dest-dir",
        type=Path,
        required=True,
        help="a durable data directory outside the repository, not a temporary directory",
    )
    args = parser.parse_args(argv)
    destination = args.dest_dir.expanduser()
    require_durable_output(parser, destination, allowed=False)
    try:
        path = download_madin_table(destination)
    except MadinPinError as error:
        raise SystemExit(str(error)) from error
    pin = MADIN_PIN
    print(f"{path} matches {pin.path} at {pin.repository}@{pin.commit[:7]} (sha256 {pin.sha256})")
    return 0


__all__ = [
    "MADIN_PIN",
    "PROVENANCE_FILE",
    "MadinPinError",
    "SourcePin",
    "download_madin_table",
    "main",
    "verify_madin_record",
    "verify_madin_table",
]


if __name__ == "__main__":
    raise SystemExit(main())
