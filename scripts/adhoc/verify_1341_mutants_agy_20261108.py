#!/usr/bin/env python3
"""Scratch-copy mutation verifier for #1341 and #463 contributor doc contracts.

Purpose:
    Verify that single-site revert mutations in ``AGENTS.md``,
    ``.github/copilot-instructions.md``, ``CONTRIBUTING.md``, and ``README.md``
    are killed by ``AssertionError`` in
    ``tests/test_sync_agent_docs.py::test_contributor_docs_agree_on_package_paths_and_check_before_commit``
    without modifying any tracked repository file.

Inputs:
    Tracked ``pyproject.toml``, ``AGENTS.md``, ``.github/copilot-instructions.md``,
    ``CONTRIBUTING.md``, ``README.md``, ``scripts/sync_agent_docs.py``, and
    ``tests/test_sync_agent_docs.py``.

Outputs:
    Per-case ``PASS``/``FAIL`` lines and exit code ``0`` when the unmutated
    control passes and all 8 mutants fail by assertion.

Issues: #1341, #463
Date: 2026-11-08
"""

from __future__ import annotations

import hashlib
import io
import shutil
import sys
import tempfile
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
TRACKED_FILES: tuple[str, ...] = (
    "pyproject.toml",
    "AGENTS.md",
    ".github/copilot-instructions.md",
    "CONTRIBUTING.md",
    "README.md",
    "scripts/sync_agent_docs.py",
    "tests/test_sync_agent_docs.py",
)


@dataclass(frozen=True)
class DocMutant:
    """Single-site revert mutation on a contributor documentation file."""

    name: str
    rel_path: str
    target: str
    replacement: str


MUTANTS: tuple[DocMutant, ...] = (
    DocMutant(
        name="M1_agents_reverts_to_before_requesting_review_1341",
        rel_path="AGENTS.md",
        target="3. **Run `just check`** — all checks must pass before committing",
        replacement="4. **Run `just check`** — all checks must pass before requesting review",
    ),
    DocMutant(
        name="M2_copilot_reverts_to_before_requesting_review_1341",
        rel_path=".github/copilot-instructions.md",
        target="3. **Run `just check`** — all checks must pass before committing",
        replacement="4. **Run `just check`** — all checks must pass before requesting review",
    ),
    DocMutant(
        name="M3_contributing_reverts_to_src_bsdbng_1341",
        rel_path="CONTRIBUTING.md",
        target="- Source code goes in `src/mtp/` (with compatibility shims in `src/bsdbng/`)",
        replacement="- Source code goes in `src/bsdbng/`",
    ),
    DocMutant(
        name="M4_contributing_drops_before_committing_sentence_1341",
        rel_path="CONTRIBUTING.md",
        target="**All checks must pass before committing.**",
        replacement="**All checks must pass.**",
    ),
    DocMutant(
        name="M5_readme_reverts_to_before_opening_a_pr_1341",
        rel_path="README.md",
        target="4. Run `just check` — **all checks must pass before committing**",
        replacement="4. Run `just check` — **all checks must pass before opening a PR**",
    ),
    DocMutant(
        name="M6_readme_reverts_to_package_remains_bsdbng_463",
        rel_path="README.md",
        target="The Python package is `mtp` (`src/mtp/`, `schema/mtp.yaml`).",
        replacement="The Python package remains `bsdbng`.",
    ),
    DocMutant(
        name="M7_readme_reverts_to_unpinned_1900_plus_studies_463",
        rel_path="README.md",
        target=(
            "just ingest      # parse into per-study YAML for the pinned BugSigDB release (~22s)"
        ),
        replacement="just ingest      # parse into per-study YAML (~22s, 1900+ studies)",
    ),
    DocMutant(
        name="M8_readme_mongo_load_reverts_to_1245_463",
        rel_path="README.md",
        target=(
            "just mongo-load         # load study YAMLs into MongoDB "
            "(~45s; for the pinned BugSigDB release)"
        ),
        replacement=(
            "just mongo-load         # load study YAMLs into MongoDB "
            "(~45s; 1,245 studies in pinned v1.3.1)"
        ),
    ),
)


def _tracked_hashes() -> dict[str, str]:
    return {
        rel: hashlib.sha256((REPO_ROOT / rel).read_bytes()).hexdigest() for rel in TRACKED_FILES
    }


def _copy_tree(scratch_root: Path) -> None:
    for rel in TRACKED_FILES:
        src = REPO_ROOT / rel
        dst = scratch_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _run_scratch_pytest(scratch_root: Path) -> tuple[int, str]:
    test_file = scratch_root / "tests" / "test_sync_agent_docs.py"
    buf = io.StringIO()
    sys.modules.pop("test_sync_agent_docs", None)
    with redirect_stdout(buf), redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    str(test_file),
                    "-k",
                    "test_contributor_docs_agree_on_package_paths_and_check_before_commit",
                    "--tb=short",
                    "-q",
                ]
            )
        )
    return rc, buf.getvalue()


def main() -> int:
    before = _tracked_hashes()
    failures = 0

    with tempfile.TemporaryDirectory(prefix="mtp-1341-control-") as tmp:
        scratch = Path(tmp)
        _copy_tree(scratch)
        rc, out = _run_scratch_pytest(scratch)
        if rc != 0:
            print(f"FAIL [CONTROL]: expected rc=0, got rc={rc}\n{out}")
            return 1
        print("PASS [CONTROL]: rc=0, PYTEST_ORIGIN_OK")

    for mutant in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mtp-1341-{mutant.name}-") as tmp:
            scratch = Path(tmp)
            _copy_tree(scratch)
            target_file = scratch / mutant.rel_path
            original_text = target_file.read_text(encoding="utf-8")
            if mutant.target not in original_text:
                print(f"FAIL [{mutant.name}]: target not found in {mutant.rel_path}")
                failures += 1
                continue
            target_file.write_text(
                original_text.replace(mutant.target, mutant.replacement, 1),
                encoding="utf-8",
            )
            rc, out = _run_scratch_pytest(scratch)
            if rc == 1 and "AssertionError" in out:
                print(f"PASS [{mutant.name}]: killed by AssertionError (rc=1, PYTEST_ORIGIN_OK)")
            else:
                print(f"FAIL [{mutant.name}]: expected rc=1 + AssertionError, got rc={rc}\n{out}")
                failures += 1

    after = _tracked_hashes()
    if before != after:
        print("FAIL [INTEGRITY]: tracked repository files changed during mutation check")
        return 1

    if failures:
        print(f"FAILED: {failures}/{len(MUTANTS)} mutants survived")
        return 1
    print(f"ALL {len(MUTANTS)}/{len(MUTANTS)} MUTANTS KILLED BY ASSERTION")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
