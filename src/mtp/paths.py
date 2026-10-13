"""Project path helpers."""

from __future__ import annotations

import hashlib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_DATA_DIR = REPO_ROOT / "data" / "raw"
STUDY_DATA_DIR = REPO_ROOT / "data" / "studies"
DERIVED_DATA_DIR = REPO_ROOT / "data" / "derived"
KG_MICROBE_DIR = REPO_ROOT / "data" / "kg-microbe"
SCHEMA_PATH = REPO_ROOT / "schema" / "mtp.yaml"

__all__ = [
    "DERIVED_DATA_DIR",
    "KG_MICROBE_DIR",
    "RAW_DATA_DIR",
    "REPO_ROOT",
    "SCHEMA_PATH",
    "STUDY_DATA_DIR",
    "file_sha256",
]


def file_sha256(path: Path | str) -> str:
    """The hex SHA-256 of a file, for provenance blocks that identify inputs."""
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()
