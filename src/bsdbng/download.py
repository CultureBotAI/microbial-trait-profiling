"""Compatibility shim forwarding bsdbng.download to mtp.download."""

import sys
from typing import Any

from mtp.download import (
    BUGSIGDB_EXPORTS_DOI,
    BUGSIGDB_EXPORTS_TAG_SHA,
    BUGSIGDB_EXPORTS_VERSION,
    DEFAULT_SOURCE,
    DEVEL_FULL_DUMP_URL,
    EXPORT_BASE,
    EXPORT_FILES,
    FULL_DUMP_NAME,
    FULL_DUMP_URL,
    INITIAL_BACKOFF_SECONDS,
    MAX_RETRIES,
    PINNED_FULL_DUMP_BYTES,
    PINNED_FULL_DUMP_MD5,
    PINNED_FULL_DUMP_URL,
    PROVENANCE_FILE,
    SNAPSHOT_RE,
    SOURCE_ENV_VAR,
    SOURCE_LIVE,
    SOURCE_PINNED,
    IntegrityError,
    SourceResolutionError,
    _fetch_with_retry,
    describe_source,
    download_exports,
    parse_snapshot,
    release_descriptor,
    resolve_source,
    verify_pinned_dump,
)

__all__ = [
    "BUGSIGDB_EXPORTS_DOI",
    "BUGSIGDB_EXPORTS_TAG_SHA",
    "BUGSIGDB_EXPORTS_VERSION",
    "DEFAULT_SOURCE",
    "DEVEL_FULL_DUMP_URL",
    "EXPORT_BASE",
    "EXPORT_FILES",
    "FULL_DUMP_NAME",
    "FULL_DUMP_URL",
    "INITIAL_BACKOFF_SECONDS",
    "MAX_RETRIES",
    "PINNED_FULL_DUMP_BYTES",
    "PINNED_FULL_DUMP_MD5",
    "PINNED_FULL_DUMP_URL",
    "PROVENANCE_FILE",
    "SNAPSHOT_RE",
    "SOURCE_ENV_VAR",
    "SOURCE_LIVE",
    "SOURCE_PINNED",
    "IntegrityError",
    "SourceResolutionError",
    "_fetch_with_retry",
    "describe_source",
    "download_exports",
    "parse_snapshot",
    "release_descriptor",
    "resolve_source",
    "verify_pinned_dump",
]

import mtp.download as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.download", run_name="__main__")
