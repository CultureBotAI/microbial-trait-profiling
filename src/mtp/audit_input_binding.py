"""Bind streamed audit parsing to the raw bytes named by its input receipt.

The digest covers the same raw stream that supplies the parser, including gzip
headers and trailers when compressed. A successful context exit consumes and
checks any unread text. A parser exception skips this final verification so the
original error is preserved. No input snapshot or temporary copy is created.
"""

from __future__ import annotations

import gzip
import hashlib
import io
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Any

_BLOCK_SIZE = 64 * 1024


class _HashingReader(io.RawIOBase):
    """Hash bounded raw reads before their bytes reach a decoding buffer."""

    def __init__(self, source: IO[bytes]) -> None:
        super().__init__()
        self._source = source
        self.digest = hashlib.sha256()
        self.size = 0

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: Any) -> int:
        view = memoryview(buffer).cast("B")
        chunk = self._source.read(min(len(view), _BLOCK_SIZE))
        length = len(chunk)
        view[:length] = chunk
        self.digest.update(chunk)
        self.size += length
        return length


@contextmanager
def verified_input_text(
    path: Path,
    expected_sha256: str,
    expected_size: int,
    *,
    gzip_compressed: bool = False,
) -> Iterator[IO[str]]:
    """Yield UTF-8 text, then bind successful parsing to its expected raw bytes.

    The context owns the returned stream; readers receiving it must not close or
    rewind it. Hashing and the final drain use bounded reads even when a caller
    returns after inspecting only a header. Compression is explicit rather than
    inferred from the nominal path. The path is used only to open the input and
    identify a verification error; this helper produces no provenance record.
    """
    with path.open("rb") as raw:
        hashed = _HashingReader(raw)
        with io.BufferedReader(hashed, buffer_size=_BLOCK_SIZE) as buffered:
            decoded = gzip.GzipFile(fileobj=buffered, mode="rb") if gzip_compressed else buffered
            with io.TextIOWrapper(decoded, encoding="utf-8", newline="") as text:
                yield text
                # This block is reached only after normal caller completion;
                # exceptions thrown through yield retain their original cause.
                while text.read(_BLOCK_SIZE):
                    pass
                if hashed.size != expected_size or hashed.digest.hexdigest() != expected_sha256:
                    raise ValueError(
                        f"Audit input bytes differ for {path}: "
                        f"expected {expected_size} bytes / {expected_sha256}, "
                        f"read {hashed.size} bytes / {hashed.digest.hexdigest()}"
                    )
