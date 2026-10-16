"""Compatibility shim forwarding bsdbng.corpus_identity to mtp.corpus_identity."""

import sys
from typing import Any

from mtp.corpus_identity import (
    DEFAULT_KG_PROVENANCE_PATH,
    UNKNOWN_RELEASE,
    CorpusIdentity,
    KgInput,
    MatrixIdentity,
    corpus_identity,
    measure_matrix,
    read_kg_inputs,
    summary_lines,
    verify_provenance_source_hashes,
    write_identity,
)

__all__ = [
    "DEFAULT_KG_PROVENANCE_PATH",
    "UNKNOWN_RELEASE",
    "CorpusIdentity",
    "KgInput",
    "MatrixIdentity",
    "corpus_identity",
    "measure_matrix",
    "read_kg_inputs",
    "summary_lines",
    "verify_provenance_source_hashes",
    "write_identity",
]

import mtp.corpus_identity as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)


if __name__ == "__main__":
    import runpy

    runpy.run_module("mtp.corpus_identity", run_name="__main__")
