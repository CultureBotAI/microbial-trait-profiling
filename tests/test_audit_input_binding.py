"""Parsed audit text must be bound to the exact streamed input bytes (#1357)."""

from __future__ import annotations

import gzip
import hashlib
import io
from contextlib import AbstractContextManager
from pathlib import Path
from typing import IO, Any

import pytest


def _verified(
    path: Path, expected: bytes, *, gzip_compressed: bool = False, size_delta: int = 0
) -> AbstractContextManager[IO[str]]:
    # Local import keeps the test-first absence a test failure, not a collection failure.
    from mtp.audit_input_binding import verified_input_text

    return verified_input_text(
        path,
        hashlib.sha256(expected).hexdigest(),
        len(expected) + size_delta,
        gzip_compressed=gzip_compressed,
    )


class _BoundedBytes(io.BytesIO):
    """Refuse unbounded raw reads, even when the consumer asks for all text."""

    def __init__(self, value: bytes) -> None:
        super().__init__(value)
        self.requests: list[int] = []
        self.bytes_read = 0

    def read(self, size: int | None = -1) -> bytes:
        assert size is not None, "unbounded raw read"
        assert 0 <= size <= 65536, "unbounded raw read"
        self.requests.append(size)
        result = super().read(size)
        self.bytes_read += len(result)
        return result

    def readinto(self, buffer: Any) -> int:
        view = memoryview(buffer)
        assert len(view) <= 65536, "oversized raw readinto"
        self.requests.append(len(view))
        count = super().readinto(view)
        self.bytes_read += count
        return count


def _install_raw(monkeypatch: pytest.MonkeyPatch, path: Path, raw: _BoundedBytes) -> list[str]:
    calls: list[str] = []

    def open_once(candidate: Path, mode: str = "r", *_args: object, **_kwargs: object) -> IO[bytes]:
        assert candidate == path, f"unexpected file open: {candidate}"
        assert mode == "rb"
        calls.append(mode)
        assert len(calls) == 1, "the nominal input was reopened"
        return raw

    # All fixture writes precede this patch; the bound reader may open only this input.
    monkeypatch.setattr(Path, "open", open_once)
    return calls


@pytest.mark.parametrize("compressed", [False, True])
def test_bound_reader_hashes_raw_bytes_and_decodes_text(tmp_path: Path, compressed: bool) -> None:
    plain = "name\tcount\r\nÉmile\t2\r\n".encode()
    raw = gzip.compress(plain, mtime=0) if compressed else plain
    path = tmp_path / "input.bin"
    path.write_bytes(raw)
    with _verified(path, raw, gzip_compressed=compressed) as text:
        assert text.read() == plain.decode()
    assert text.closed


def test_same_size_valid_text_flip_is_rejected_on_exit(tmp_path: Path) -> None:
    expected = b"name\tcount\na\t2\n"
    changed = expected.replace(b"2", b"9")
    assert len(changed) == len(expected)
    path = tmp_path / "input.tsv"
    path.write_bytes(changed)
    parsed: list[str] = []
    with (
        pytest.raises(ValueError, match="Audit input bytes differ"),
        _verified(path, expected) as text,
    ):
        parsed.extend(text.readlines())
    assert parsed[-1] == "a\t9\n"


@pytest.mark.parametrize("size_delta", [-1, 1])
def test_correct_digest_with_wrong_raw_size_is_rejected(tmp_path: Path, size_delta: int) -> None:
    raw = b"name\tcount\na\t2\n"
    path = tmp_path / "input.tsv"
    path.write_bytes(raw)
    with (
        pytest.raises(ValueError, match="Audit input bytes differ"),
        _verified(path, raw, size_delta=size_delta) as text,
    ):
        assert text.readline() == "name\tcount\n"


def test_gzip_header_flip_is_detected_even_when_decoded_text_matches(tmp_path: Path) -> None:
    raw = gzip.compress(b"name\tcount\na\t2\n", mtime=0)
    changed = bytearray(raw)
    changed[4] ^= 1  # MTIME differs; both members decode to the same valid text.
    assert gzip.decompress(changed) == gzip.decompress(raw)
    path = tmp_path / "input.gz"
    path.write_bytes(changed)
    with (
        pytest.raises(ValueError, match="Audit input bytes differ"),
        _verified(path, raw, gzip_compressed=True) as text,
    ):
        assert text.read() == "name\tcount\na\t2\n"


@pytest.mark.parametrize("compressed", [False, True])
def test_early_header_return_drains_all_bytes_in_bounded_reads(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, compressed: bool
) -> None:
    plain = b"header\n" + b"body row\n" * 300000
    encoded = gzip.compress(plain, mtime=0) if compressed else plain
    path = tmp_path / "nominal-input"
    raw = _BoundedBytes(encoded)
    calls = _install_raw(monkeypatch, path, raw)

    def header_only() -> str:
        with _verified(path, encoded, gzip_compressed=compressed) as text:
            return text.readline()

    assert header_only() == "header\n"
    assert raw.bytes_read == len(encoded)
    assert raw.requests
    assert max(raw.requests) <= 65536
    assert len(calls) == 1
    assert raw.closed


def test_early_header_return_refuses_a_changed_unread_tail(tmp_path: Path) -> None:
    expected = b"header\n" + b"body row\n" * 30000
    path = tmp_path / "input.tsv"
    path.write_bytes(expected[:-2] + b"X\n")
    with (
        pytest.raises(ValueError, match="Audit input bytes differ"),
        _verified(path, expected) as text,
    ):
        assert text.readline() == "header\n"


def test_early_header_return_still_decodes_the_unread_text(tmp_path: Path) -> None:
    raw = b"header\n" + b"body row\n" * 30000 + b"\xff\n"
    path = tmp_path / "input.tsv"
    path.write_bytes(raw)
    with pytest.raises(UnicodeDecodeError), _verified(path, raw) as text:
        assert text.readline() == "header\n"


def test_consumer_full_read_remains_bounded_and_uses_one_raw_open(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    expected = b"header\n" + b"0123456789\n" * 150000
    raw = _BoundedBytes(expected)
    path = tmp_path / "nominal-input"
    calls = _install_raw(monkeypatch, path, raw)
    with _verified(path, expected) as text:
        assert text.read() == expected.decode()
    assert len(calls) == 1
    assert raw.bytes_read == len(expected)
    assert len(raw.requests) > 1


def test_parser_exception_is_not_replaced_by_hash_failure_or_drain(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    raw = _BoundedBytes(b"different bytes and an invalid tail \xff")
    path = tmp_path / "nominal-input"
    _install_raw(monkeypatch, path, raw)
    error = RuntimeError("original parser failure")
    with pytest.raises(RuntimeError) as caught, _verified(path, b"expected bytes"):
        raise error
    assert caught.value is error
    assert raw.bytes_read == 0
    assert raw.closed
