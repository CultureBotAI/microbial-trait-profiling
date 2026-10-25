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
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

import httpx

from mtp.paths import file_sha256

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


def _check(label: str, size: int, digest: str) -> None:
    pin = MADIN_PIN
    where = f"{pin.path} at {pin.repository}@{pin.commit[:7]}"
    if size != pin.size_bytes:
        msg = f"{label} is {size:,} bytes, expected {pin.size_bytes:,} for {where}"
        raise MadinPinError(msg)
    if digest != pin.sha256:
        msg = f"{label} has sha256 {digest}, expected {pin.sha256} for {where}"
        raise MadinPinError(msg)


def verify_madin_table(path: Path) -> dict[str, object]:
    """Return the pin's descriptor if *path* has its size and sha256; refuse it otherwise."""
    _check(path.name, path.stat().st_size, file_sha256(path))
    return MADIN_PIN.descriptor()


def verify_madin_record(label: str, size_bytes: int, sha256: str) -> dict[str, object]:
    """Return the pin's descriptor if a recorded size and sha256 are the pin's; refuse them."""
    _check(label, size_bytes, sha256)
    return MADIN_PIN.descriptor()


def download_madin_table(dest_dir: Path) -> Path:
    """Fetch the pinned table into *dest_dir*, anonymously; never keep bytes that fail the pin.

    An existing copy is checked and reused, never overwritten. The request carries no
    credentials or identifying headers.
    """
    pin = MADIN_PIN
    dest_dir.mkdir(parents=True, exist_ok=True)
    target = dest_dir / pin.filename
    if target.exists():
        try:
            verify_madin_table(target)
        except MadinPinError as error:
            msg = f"{error}; remove it to fetch the pinned file again"
            raise MadinPinError(msg) from error
        return target
    partial = target.with_name(f"{target.name}.part")
    digest = hashlib.sha256()
    size = 0
    try:
        with (
            httpx.Client(follow_redirects=True, timeout=120) as client,
            client.stream("GET", pin.url) as response,
        ):
            response.raise_for_status()
            with partial.open("wb") as handle:
                for block in response.iter_bytes(1 << 20):
                    handle.write(block)
                    digest.update(block)
                    size += len(block)
        _check(pin.url, size, digest.hexdigest())
        partial.replace(target)
    finally:
        partial.unlink(missing_ok=True)
    provenance_path = dest_dir / PROVENANCE_FILE
    provenance = (
        json.loads(provenance_path.read_text(encoding="utf-8")) if provenance_path.exists() else {}
    )
    provenance[pin.filename] = {
        **pin.descriptor(),
        "downloaded_at": datetime.now(tz=UTC).isoformat(),
        "integrity_verified": True,
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
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
    try:
        path = download_madin_table(args.dest_dir)
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
