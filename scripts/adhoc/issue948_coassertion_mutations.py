"""Check the four #948 audit mutants against synthetic tests (2026-10-23).

Inputs: current coassertion audit source and fixture tests. Outputs: pytest logs
and JSON mutation outcomes. Each invocation mutates one source expression and
restores exact bytes in a finally block; run only in the issue worktree.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest

from mtp.paths import REPO_ROOT

MUTATIONS = {
    "off_strain_evidence": (
        "else admitted.get(subject, set())",
        "else graph.strain_traits.get(subject, set())",
        1,
    ),
    "off_strain_attribution": (
        "if admitted.get(strain, set()).intersection",
        "if graph.strain_traits.get(strain, set()).intersection",
        2,
    ),
    "off_strain_inventory": (
        "for strain, keys in admitted.items()",
        "for strain, keys in graph.strain_traits.items()",
        1,
    ),
    "same_strain_union": (
        "(direct_left and direct_right) or same_strains",
        "(direct_left and direct_right)",
        1,
    ),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mutation", choices=sorted(MUTATIONS), required=True)
    args = parser.parse_args()
    if "mtp.coassertion_audit" in sys.modules:
        raise ValueError("Run each mutation in a fresh interpreter before importing the audit")
    source = REPO_ROOT / "src/mtp/coassertion_audit.py"
    original = source.read_bytes()
    text = original.decode("utf-8")
    old, new, expected = MUTATIONS[args.mutation]
    if text.count(old) != expected:
        raise ValueError("Mutation target changed; inspect the current implementation")
    output = io.StringIO()
    test = "test_scope_counts_preserve_sources_and_do_not_propagate_to_ancestors"
    try:
        source.write_text(text.replace(old, new), encoding="utf-8")
        with redirect_stdout(output), redirect_stderr(output):
            code = pytest.main(["-q", f"{REPO_ROOT}/tests/test_coassertion_audit.py::{test}"])
    finally:
        source.write_bytes(original)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / f"{args.mutation}.log").write_text(output.getvalue(), encoding="utf-8")
    record = {
        "mutation": args.mutation,
        "pytest_exit": int(code),
        "killed": code == pytest.ExitCode.TESTS_FAILED,
        "original_source_sha256": hashlib.sha256(original).hexdigest(),
        "restored_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    }
    (args.output_dir / f"{args.mutation}.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(record, indent=2))
    return 0 if record["killed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
