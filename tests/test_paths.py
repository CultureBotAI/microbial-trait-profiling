"""Tests for mtp.paths and canonical path/file helpers."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

import bsdbng.paths
import mtp.paths
from mtp.paths import SCHEMA_PATH, file_sha256


def test_file_sha256_empty_file(tmp_path: Path) -> None:
    """Empty file returns known empty-file SHA-256 digest."""
    empty_file = tmp_path / "empty.txt"
    empty_file.write_bytes(b"")
    expected = hashlib.sha256(b"").hexdigest()
    assert expected == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    assert file_sha256(empty_file) == expected
    # Also verify str argument works
    assert file_sha256(str(empty_file)) == expected


def test_file_sha256_known_content(tmp_path: Path) -> None:
    """Content with known bytes returns correct SHA-256 digest."""
    sample_file = tmp_path / "sample.txt"
    content = b"hello world\n"
    sample_file.write_bytes(content)
    expected = hashlib.sha256(content).hexdigest()
    assert file_sha256(sample_file) == expected
    assert file_sha256(str(sample_file)) == expected


def test_file_sha256_schema_path() -> None:
    """Computes correct SHA-256 for an actual existing repository file."""
    assert SCHEMA_PATH.is_file()
    expected = hashlib.sha256(SCHEMA_PATH.read_bytes()).hexdigest()
    assert file_sha256(SCHEMA_PATH) == expected


def test_file_sha256_missing_file_raises(tmp_path: Path) -> None:
    """Missing file raises FileNotFoundError."""
    missing = tmp_path / "does_not_exist.txt"
    with pytest.raises(FileNotFoundError):
        file_sha256(missing)


def test_file_sha256_shim_parity() -> None:
    """Shim bsdbng.paths.file_sha256 is identical object to mtp.paths.file_sha256."""
    assert bsdbng.paths.file_sha256 is mtp.paths.file_sha256

