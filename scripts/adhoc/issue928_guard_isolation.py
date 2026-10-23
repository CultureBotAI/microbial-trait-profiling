"""Verify #928 guard-test isolation (2026-10-23).

Inputs: infographic producer, its synthetic guard tests, and existing checkout
outputs. Outputs: pytest log and before/after hash inventory. Temporarily removes
the argument guard, runs pytest directly, and restores exact producer bytes.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import runpy
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest

from mtp.paths import REPO_ROOT, file_sha256


def inventory(root: Path) -> dict[str, str]:
    """Include ignored and untracked files as well as tracked checkout outputs."""
    return {
        path.relative_to(root).as_posix(): file_sha256(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    producer = REPO_ROOT / "scripts/create_trait_profile_infographic.py"
    original = producer.read_bytes()
    namespace = runpy.run_path(str(producer))
    checkout_outputs = namespace["INVESTIGATION_DIR"]
    if not isinstance(checkout_outputs, Path):
        raise ValueError("The producer must declare a Path output directory")
    before = inventory(checkout_outputs)
    if not before:
        raise ValueError("This isolation audit requires existing checkout outputs")
    guard = (
        "        if not args.stratified_stats_only:\n"
        '            parser.error("--min-stratified-pairs requires --stratified-stats-only")\n'
    )
    text = original.decode("utf-8")
    if text.count(guard) != 1:
        raise ValueError("Guard target changed; inspect the producer before mutation")
    captured = io.StringIO()
    try:
        producer.write_bytes(text.replace(guard, "").encode("utf-8"))
        with redirect_stdout(captured), redirect_stderr(captured):
            code = pytest.main(
                [
                    "-q",
                    f"{REPO_ROOT}/tests/test_trait_profile_infographic.py::"
                    "test_min_stratified_pairs_is_rejected_outside_site_stats",
                ]
            )
    finally:
        producer.write_bytes(original)
    after = inventory(checkout_outputs)
    record = {
        "pytest_exit": int(code),
        "mutation_detected": code == pytest.ExitCode.TESTS_FAILED,
        "checkout_outputs": checkout_outputs.relative_to(REPO_ROOT).as_posix(),
        "before": before,
        "after": after,
        "checkout_unchanged": before == after,
        "producer_original_sha256": hashlib.sha256(original).hexdigest(),
        "producer_restored_sha256": file_sha256(producer),
    }
    (args.output_dir / "pytest.log").write_text(captured.getvalue(), encoding="utf-8")
    (args.output_dir / "isolation.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    print(f"mutation_detected={record['mutation_detected']}; checkout_unchanged={before == after}")
    return 0 if record["mutation_detected"] and record["checkout_unchanged"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
