#!/usr/bin/env python3
"""Scratch-isolated mutation verifier for issue #2086.

Purpose:
    Verify that tests/test_taxon_influence_date_scope_2086.py kills mutations across:
    - methods/table1_taxon_influence_results.md (reverting to timeless present tense
      'remain unreviewed', dropping 'historical `190-v2`', dropping the status
      checkpoint link, dropping the non-recertification / biological-validation
      caveat, or altering #188 / #189 closure dates)
    - methods/table1_taxon_influence.md (reverting to 'remain subject to #188;
      feature semantics remain subject to #189', dropping 'historical `190-v2`',
      dropping the status checkpoint link, dropping the non-recertification /
      biological-validation caveat, or altering #188 / #189 closure dates)
    - methods/biological_interpretation_audit.md (altering the status checkpoint
      heading, status, or closure dates)
    - tests/test_taxon_influence_date_scope_2086.py (weakening the negative guard)
    without mutating the live worktree in place.

Inputs:
    - methods/table1_taxon_influence_results.md
    - methods/table1_taxon_influence.md
    - methods/biological_interpretation_audit.md
    - tests/test_taxon_influence_date_scope_2086.py

Outputs:
    - <scratch-dir>/summary.json
    - <scratch-dir>/verdicts.json

Issue: #2086
Date: 2026-11-10
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
import shutil
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from types import ModuleType

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_REL = Path("methods/table1_taxon_influence_results.md")
METHODS_REL = Path("methods/table1_taxon_influence.md")
AUDIT_REL = Path("methods/biological_interpretation_audit.md")
TEST_REL = Path("tests/test_taxon_influence_date_scope_2086.py")


@dataclass(frozen=True)
class MutationSpec:
    """Single scratch-isolated mutation specification for #2086."""

    name: str
    rel_path: Path
    needle: str
    replacement: str
    test_node: str
    expected_detail_substring: str


@dataclass(frozen=True)
class MutationVerdict:
    """Recorded outcome for a single scratch-isolated mutant execution."""

    name: str
    target_file: str
    test_node: str
    killed: bool
    exit_code: int
    original_sha256: str
    mutated_sha256: str
    expected_detail_substring: str
    matched_detail: bool
    detail: str


MUTATIONS: tuple[MutationSpec, ...] = (
    MutationSpec(
        name="results_revert_remain_unreviewed",
        rel_path=RESULTS_REL,
        needle="were unreviewed at execution time",
        replacement="remain unreviewed",
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring=(
            "Timeless present-tense phrase 'remain unreviewed' must not appear"
        ),
    ),
    MutationSpec(
        name="results_drop_historical_190_v2",
        rel_path=RESULTS_REL,
        needle="In the historical `190-v2` run, contrast orientation (#188)",
        replacement="In this run, contrast orientation (#188)",
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring=(
            "Expected date-scoping to the historical 190-v2 run in "
            "table1_taxon_influence_results.md"
        ),
    ),
    MutationSpec(
        name="results_drop_checkpoint_link",
        rel_path=RESULTS_REL,
        needle=(
            "[`methods/biological_interpretation_audit.md#status-checkpoint-2026-10-30`]"
            "(biological_interpretation_audit.md#status-checkpoint-2026-10-30)"
        ),
        replacement="methods/biological_interpretation_audit.md",
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring=(
            "Expected relative link to "
            "biological_interpretation_audit.md#status-checkpoint-2026-10-30"
        ),
    ),
    MutationSpec(
        name="results_drop_non_recertification",
        rel_path=RESULTS_REL,
        needle=(
            "does not retroactively recertify the historical `190-v2` outputs or "
            "establish biological validation"
        ),
        replacement="completes the historical tracker record",
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring=(
            "Expected non-recertification statement for historical 190-v2 outputs"
        ),
    ),
    MutationSpec(
        name="results_wrong_188_date",
        rel_path=RESULTS_REL,
        needle="2026-10-11",
        replacement="2026-10-12",
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring=(
            "Expected #188 closure date (2026-10-11 UTC) in table1_taxon_influence_results.md"
        ),
    ),
    MutationSpec(
        name="results_wrong_189_date",
        rel_path=RESULTS_REL,
        needle="2026-10-08",
        replacement="2026-10-09",
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring=(
            "Expected #189 closure date (2026-10-08 UTC) in table1_taxon_influence_results.md"
        ),
    ),
    MutationSpec(
        name="methods_revert_remain_subject",
        rel_path=METHODS_REL,
        needle="were scoped to #188 and feature\nsemantics to #189",
        replacement="remain subject to #188; feature semantics remain subject to #189",
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring=(
            "Timeless present-tense phrase 'remain subject to #188' must not appear"
        ),
    ),
    MutationSpec(
        name="methods_drop_historical_190_v2",
        rel_path=METHODS_REL,
        needle="In the historical\n`190-v2` run, contrast eligibility",
        replacement="In the current\nrun, contrast eligibility",
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring=(
            "Expected date-scoping to the historical 190-v2 run in table1_taxon_influence.md"
        ),
    ),
    MutationSpec(
        name="methods_drop_checkpoint_link",
        rel_path=METHODS_REL,
        needle=(
            "[`methods/biological_interpretation_audit.md#status-checkpoint-2026-10-30`]"
            "(biological_interpretation_audit.md#status-checkpoint-2026-10-30)"
        ),
        replacement="methods/biological_interpretation_audit.md",
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring=(
            "Expected relative link to "
            "biological_interpretation_audit.md#status-checkpoint-2026-10-30"
        ),
    ),
    MutationSpec(
        name="methods_drop_non_recertification",
        rel_path=METHODS_REL,
        needle=(
            "does not retroactively recertify the historical `190-v2`\n"
            "outputs or establish biological validation"
        ),
        replacement="completes the tracker milestone for these prerequisites",
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring=(
            "Expected non-recertification statement for historical 190-v2 outputs"
        ),
    ),
    MutationSpec(
        name="methods_wrong_188_date",
        rel_path=METHODS_REL,
        needle="2026-10-11",
        replacement="2026-10-12",
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring=(
            "Expected #188 closure date (2026-10-11 UTC) in table1_taxon_influence.md"
        ),
    ),
    MutationSpec(
        name="methods_wrong_189_date",
        rel_path=METHODS_REL,
        needle="2026-10-08",
        replacement="2026-10-09",
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring=(
            "Expected #189 closure date (2026-10-08 UTC) in table1_taxon_influence.md"
        ),
    ),
    MutationSpec(
        name="audit_drop_checkpoint_heading",
        rel_path=AUDIT_REL,
        needle="## Status checkpoint (2026-10-30)",
        replacement="## Historical status table",
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring=(
            "Missing status checkpoint section in biological_interpretation_audit.md"
        ),
    ),
    MutationSpec(
        name="audit_mark_188_open",
        rel_path=AUDIT_REL,
        needle=(
            "| [#188](https://github.com/CultureBotAI/microbial-trait-profiling/issues/188)"
            ", contrast orientation | Closed | 2026-10-11 |"
        ),
        replacement=(
            "| [#188](https://github.com/CultureBotAI/microbial-trait-profiling/issues/188)"
            ", contrast orientation | Open | 2026-10-11 |"
        ),
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring="Expected #188 to be Closed in status checkpoint",
    ),
    MutationSpec(
        name="audit_wrong_189_date",
        rel_path=AUDIT_REL,
        needle=(
            "| [#189](https://github.com/CultureBotAI/microbial-trait-profiling/issues/189)"
            ", trait semantics | Closed | 2026-10-08 |"
        ),
        replacement=(
            "| [#189](https://github.com/CultureBotAI/microbial-trait-profiling/issues/189)"
            ", trait semantics | Closed | 2026-10-10 |"
        ),
        test_node="test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
        expected_detail_substring="Expected #189 closure date 2026-10-08 in status checkpoint",
    ),
    MutationSpec(
        name="test_weaken_negative_guard",
        rel_path=TEST_REL,
        needle='        for forbidden in (\n            "remain unreviewed",',
        replacement="        for forbidden in (",
        test_node=(
            "test_taxon_influence_date_scope_guard_rejects_unscoped_or_inconsistent_text"
            "[results-timeless-remain-unreviewed]"
        ),
        expected_detail_substring="DID NOT RAISE",
    ),
)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _copy_baseline_tree(scratch_case_dir: Path) -> None:
    if scratch_case_dir.exists():
        shutil.rmtree(scratch_case_dir)
    for rel in (RESULTS_REL, METHODS_REL, AUDIT_REL, TEST_REL):
        src = REPO_ROOT / rel
        dst = scratch_case_dir / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _load_module_from_path(module_name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module spec from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _run_scratch_test_node(scratch_case_dir: Path, test_node: str) -> tuple[int, str]:
    test_file = scratch_case_dir / TEST_REL
    sys.modules.pop("test_taxon_influence_date_scope_2086", None)
    importlib.invalidate_caches()
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    old_stdout, old_stderr = sys.stdout, sys.stderr
    try:
        sys.stdout = stdout_buf
        sys.stderr = stderr_buf
        exit_code = int(
            pytest.main(
                [
                    "-q",
                    "--tb=short",
                    "--import-mode=importlib",
                    "-o",
                    f"cache_dir={scratch_case_dir / '.pytest_cache'}",
                    f"{test_file}::{test_node}",
                ]
            )
        )
    finally:
        sys.stdout = old_stdout
        sys.stderr = old_stderr
        sys.modules.pop("test_taxon_influence_date_scope_2086", None)
    combined = (stdout_buf.getvalue() + "\n" + stderr_buf.getvalue()).strip()
    if exit_code == 0:
        return 0, combined or "passed"

    try:
        mod = _load_module_from_path(f"_scratch_2086_{scratch_case_dir.name}", test_file)
        if "[" not in test_node:
            fn = getattr(mod, test_node)
            fn()
    except BaseException as exc:
        exc_detail = f"{type(exc).__name__}: {exc}"
        if exc_detail not in combined:
            combined = f"{exc_detail}\n{combined}".strip()
    return exit_code, combined


def _evaluate_mutation(spec: MutationSpec, scratch_root: Path) -> MutationVerdict:
    case_dir = scratch_root / spec.name
    _copy_baseline_tree(case_dir)

    target_path = case_dir / spec.rel_path
    original_bytes = target_path.read_bytes()
    original_text = original_bytes.decode("utf-8")
    count = original_text.count(spec.needle)
    if count != 1:
        raise RuntimeError(
            f"Expected needle for {spec.name} to appear once in {spec.rel_path}, found {count}"
        )

    mutated_text = original_text.replace(spec.needle, spec.replacement, 1)
    mutated_bytes = mutated_text.encode("utf-8")
    target_path.write_bytes(mutated_bytes)

    exit_code, detail = _run_scratch_test_node(case_dir, spec.test_node)
    matched_detail = spec.expected_detail_substring in detail
    killed = exit_code == 1 and matched_detail
    return MutationVerdict(
        name=spec.name,
        target_file=str(spec.rel_path),
        test_node=spec.test_node,
        killed=killed,
        exit_code=exit_code,
        original_sha256=_sha256_bytes(original_bytes),
        mutated_sha256=_sha256_bytes(mutated_bytes),
        expected_detail_substring=spec.expected_detail_substring,
        matched_detail=matched_detail,
        detail=detail,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scratch-dir",
        "--output-dir",
        dest="scratch_dir",
        type=Path,
        required=True,
        help="Scratch directory outside the git worktree for mutant copies and JSON artifacts.",
    )
    parser.add_argument(
        "--canary",
        action="store_true",
        help="Run only the first mutation specification as a single-unit canary.",
    )
    parser.add_argument(
        "--mutation",
        type=str,
        default=None,
        help="Run only the named mutation specification.",
    )
    args = parser.parse_args(argv)

    scratch_root: Path = args.scratch_dir.resolve()
    scratch_root.mkdir(parents=True, exist_ok=True)

    baseline_dir = scratch_root / "_baseline"
    _copy_baseline_tree(baseline_dir)
    baseline_code, baseline_detail = _run_scratch_test_node(
        baseline_dir,
        "test_taxon_influence_docs_date_scope_188_and_189_to_historical_190_v2_run",
    )
    if baseline_code != 0:
        raise RuntimeError(f"Baseline scratch test failed before mutation run: {baseline_detail}")

    selected: list[MutationSpec]
    if args.mutation is not None:
        selected = [m for m in MUTATIONS if m.name == args.mutation]
        if not selected:
            raise ValueError(f"Unknown mutation name: {args.mutation}")
    elif args.canary:
        selected = [MUTATIONS[0]]
    else:
        selected = list(MUTATIONS)

    verdicts: list[MutationVerdict] = []
    for spec in selected:
        verdict = _evaluate_mutation(spec, scratch_root)
        verdicts.append(verdict)

    killed_count = sum(1 for v in verdicts if v.killed)
    survived = [v.name for v in verdicts if not v.killed]

    verdicts_payload = {
        "issue": 2086,
        "mode": "canary" if args.canary else ("single" if args.mutation else "full"),
        "scratch_dir": str(scratch_root),
        "verdicts": [asdict(v) for v in verdicts],
    }
    summary_payload = {
        "issue": 2086,
        "mode": verdicts_payload["mode"],
        "scratch_dir": str(scratch_root),
        "baseline_passed": True,
        "baseline_detail": baseline_detail,
        "total": len(verdicts),
        "killed": killed_count,
        "survived": len(survived),
        "survived_names": survived,
        "all_killed": len(survived) == 0,
        "verdict_summaries": [
            {
                "name": v.name,
                "target_file": v.target_file,
                "test_node": v.test_node,
                "killed": v.killed,
                "matched_detail": v.matched_detail,
                "expected_detail_substring": v.expected_detail_substring,
                "detail_preview": v.detail.splitlines()[0] if v.detail else "",
            }
            for v in verdicts
        ],
    }

    verdicts_path = scratch_root / "verdicts.json"
    summary_path = scratch_root / "summary.json"
    verdicts_path.write_text(json.dumps(verdicts_payload, indent=2) + "\n", encoding="utf-8")
    summary_path.write_text(json.dumps(summary_payload, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(summary_payload, indent=2))
    return 0 if not survived else 1


if __name__ == "__main__":
    raise SystemExit(main())
