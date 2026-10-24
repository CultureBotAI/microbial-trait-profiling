"""Reproduce #933 sensitivity mutants with synthetic tests (2026-10-23).

Inputs: this checkout's two presence/depth modules and regression tests.
Outputs: per-case pytest logs and a JSON result/source-restoration inventory.
Run only in an isolated, otherwise idle worktree: source is temporarily mutated
and restored in finally blocks. No production inputs or outputs are used.
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from pathlib import Path

import pytest

import mtp.presence_depth_families as families
import mtp.presence_depth_sensitivity as sensitivity
from mtp.paths import REPO_ROOT, file_sha256

TEST_FILE = "tests/test_presence_depth_sensitivity.py"


@dataclass(frozen=True)
class Mutation:
    """One reviewed regression and the assertions expected to reject it."""

    name: str
    module: str
    original: str
    replacement: str
    tests: tuple[str, ...]


MUTATIONS = (
    Mutation(
        "overwrite_primary_q",
        "presence_depth_families",
        "q_by_six_trait=adjusted.q_by,\n                    q_bh_six_trait=adjusted.q_bh,",
        "q_by=adjusted.q_by,\n                    q_bh=adjusted.q_bh,",
        ("test_complete_family_keeps_untestable_slots_in_authoritative_correction",),
    ),
    Mutation(
        "swap_focal_slug",
        "presence_depth_sensitivity",
        "focal_trait=slug,",
        "focal_trait=next(iter(TABLE1_TRAITS)),",
        ("test_focal_rows_keep_primary_corrections_and_exact_trait_identity",),
    ),
    Mutation(
        "revert_rollup_vocabulary",
        "presence_depth_sensitivity",
        'rollup_policy=f"ncbi_taxonomy_closure_rollup_{rollup}",',
        'rollup_policy=f"ncbi_taxonomy_rollup_{rollup}",',
        ("test_verified_core_run_records_actual_family_and_population[off-False]",),
    ),
    Mutation(
        "allow_extra_ora_population",
        "presence_depth_sensitivity",
        "if set(ora_rows) != expected:",
        "if not expected.issubset(ora_rows):",
        ("test_binary_ora_population_requires_exact_equality[extra]",),
    ),
    Mutation(
        "remove_sign_tally_guard",
        "presence_depth_sensitivity",
        'raise ValueError("Native and reconstructed study sign tallies disagree")',
        "pass",
        ("test_native_sign_tally_mismatch_is_rejected",),
    ),
)


def reload_sources() -> None:
    """Reload byte-exact source and test imports, avoiding cached pre-mutation functions."""
    for module in (families, sensitivity):
        cached = getattr(module, "__cached__", None)
        if cached:
            Path(cached).unlink(missing_ok=True)
        importlib.invalidate_caches()
        importlib.reload(module)
    test_path = (REPO_ROOT / TEST_FILE).resolve()
    for name, module in sys.modules.copy().items():
        filename = getattr(module, "__file__", None)
        if filename and Path(filename).resolve() == test_path:
            sys.modules.pop(name)


def run_tests(output: Path, name: str, tests: tuple[str, ...]) -> int:
    """Run only the specified synthetic regressions and retain their assertion output."""
    reload_sources()
    with (
        (output / f"{name}.log").open("w", encoding="utf-8") as log,
        redirect_stdout(log),
        redirect_stderr(log),
    ):
        return int(
            pytest.main(
                [
                    "-q",
                    "--basetemp",
                    str(output / name),
                    *(f"{TEST_FILE}::{test}" for test in tests),
                ]
            )
        )


def main() -> int:
    """Require a passing control, rejected mutants, and exact source restoration."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    paths = {
        name: REPO_ROOT / "src/mtp" / f"{name}.py"
        for name in ("presence_depth_families", "presence_depth_sensitivity")
    }
    originals = {name: path.read_bytes() for name, path in paths.items()}
    before = {name: file_sha256(path) for name, path in paths.items()}
    results = []
    try:
        controls = tuple(dict.fromkeys(test for item in MUTATIONS for test in item.tests))
        results.append({"name": "baseline", "exit_code": run_tests(output, "baseline", controls)})
        for mutation in MUTATIONS:
            path = paths[mutation.module]
            source = originals[mutation.module].decode("utf-8")
            if mutation.original not in source:
                raise ValueError(f"Mutation target absent: {mutation.name}")
            try:
                path.write_text(
                    source.replace(mutation.original, mutation.replacement), encoding="utf-8"
                )
                code = run_tests(output, mutation.name, mutation.tests)
                results.append({"name": mutation.name, "exit_code": code})
            finally:
                path.write_bytes(originals[mutation.module])
    finally:
        for name, path in paths.items():
            if path.read_bytes() != originals[name]:
                path.write_bytes(originals[name])
        reload_sources()
        evidence = {
            "results": results,
            "source_sha256_before": before,
            "source_sha256_after": {name: file_sha256(path) for name, path in paths.items()},
        }
        (output / "evidence.json").write_text(
            json.dumps(evidence, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
    print(json.dumps(evidence, indent=2))
    return int(
        results[0]["exit_code"] != 0
        or len(results) != len(MUTATIONS) + 1
        or any(row["exit_code"] != 1 for row in results[1:])
        or evidence["source_sha256_before"] != evidence["source_sha256_after"]
    )


if __name__ == "__main__":
    raise SystemExit(main())
