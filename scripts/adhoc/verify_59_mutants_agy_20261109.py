"""Isolated scratch-copy mutation verifier for Issue #59 (README.md release docs).

Purpose:
  Verify that tests/test_release_documentation.py kills 12 distinct mutants of
  README.md and .github/workflows/release.yaml in an isolated scratch directory,
  and verify that the test binds to the scratch directory rather than the live
  checkout.
Inputs:
  - README.md
  - .github/workflows/release.yaml
  - tests/test_release_documentation.py
Outputs:
  - Stdout summary table of control + 12 mutant runs; exits 0 iff all 12 mutants
    fail with AssertionError.
Issue: #59
Date: 2026-11-09
"""

from __future__ import annotations

import contextlib
import io
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
REL_FILES: tuple[str, ...] = (
    "README.md",
    ".github/workflows/release.yaml",
    "tests/test_release_documentation.py",
)

MUTANTS: tuple[tuple[str, str, str, str], ...] = (
    (
        "M01_drop_releases_h2_heading",
        "README.md",
        "## Releases\n",
        "## Release notes archive\n",
    ),
    (
        "M02_replace_canonical_repo_with_legacy_turbomam",
        "README.md",
        "https://github.com/CultureBotAI/microbial-trait-profiling/releases",
        "https://github.com/turbomam/bsdbng/releases",
    ),
    (
        "M03_negate_human_approval_boundary",
        "README.md",
        "Explicit human approval is required before tagging\nor publishing a release",
        "Releases may be tagged without explicit human approval",
    ),
    (
        "M04_drop_ask_before_downloading_rule",
        "README.md",
        "and always ask the user before\ndownloading large release assets.\n",
        "and download release assets immediately.\n",
    ),
    (
        "M05_omit_stats_asset_from_readme_table",
        "README.md",
        "| `bsdbng-stats-<tag>.txt` | Summary counts and statistics output from `just stats` |\n",
        "",
    ),
    (
        "M06_omit_pipeline_log_asset_from_readme_table",
        "README.md",
        "| `bsdbng-pipeline-log-<tag>.log` | Pipeline step timing log (`logs/pipeline.log`) |\n",
        "",
    ),
    (
        "M07_drop_web_ui_create_new_tag_on_publish_step",
        "README.md",
        'and click **"Create new tag on publish"**',
        "and save draft",
    ),
    (
        "M08_drop_cli_just_check_step",
        "README.md",
        "just check          # verify all checks pass before tagging\n",
        "",
    ),
    (
        "M09_drop_cli_git_push_tag_step",
        "README.md",
        "git push origin vX.Y.Z\n",
        "git push origin main\n",
    ),
    (
        "M10_drop_workflow_validate_studies_reference",
        "README.md",
        "`just ingest`, `just validate-studies`, and `just stats`",
        "`just ingest`, and `just stats`",
    ),
    (
        "M11_drop_post_release_placeholder_taxon_check",
        "README.md",
        "rather than `taxon_<ID>` placeholders.",
        "in the output directory.",
    ),
    (
        "M12_workflow_adds_undocumented_sixth_asset",
        ".github/workflows/release.yaml",
        "            bsdbng-stats-*.txt\n",
        "            bsdbng-stats-*.txt\n            bsdbng-extra-*.tar.gz\n",
    ),
)


def _populate_scratch(scratch_root: Path) -> None:
    for rel in REL_FILES:
        src = REPO_ROOT / rel
        dst = scratch_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _run_pytest_in_scratch(scratch_root: Path) -> tuple[int, str]:
    for mod_name in list(sys.modules):
        if mod_name == "test_release_documentation" or mod_name.startswith(
            "tests.test_release_documentation"
        ):
            sys.modules.pop(mod_name, None)

    buf = io.StringIO()
    test_file = scratch_root / "tests" / "test_release_documentation.py"
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    "-q",
                    "--tb=native",
                    "-o",
                    "addopts=",
                    str(test_file),
                ]
            )
        )
    return rc, buf.getvalue()


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="verify_59_mutants_") as tmp_dir:
        scratch_root = Path(tmp_dir)
        _populate_scratch(scratch_root)

        # 1. Verify unmodified control passes
        control_rc, control_out = _run_pytest_in_scratch(scratch_root)
        if control_rc != 0:
            print("CONTROL FAILED UNEXPECTEDLY:\n" + control_out)
            return 1

        # 2. Verify scratch origin binding: append a canary test asserting REPO_ROOT == scratch_root
        test_path = scratch_root / "tests" / "test_release_documentation.py"
        original_test_text = test_path.read_text(encoding="utf-8")
        test_path.write_text(
            original_test_text
            + "\n\ndef test_scratch_origin_binding() -> None:\n"
            + f"    assert REPO_ROOT == Path({str(scratch_root)!r}).resolve()\n",
            encoding="utf-8",
        )
        origin_rc, origin_out = _run_pytest_in_scratch(scratch_root)
        test_path.write_text(original_test_text, encoding="utf-8")
        if origin_rc != 0:
            print("SCRATCH ORIGIN BINDING FAILED:\n" + origin_out)
            return 1

        print(f"Control + scratch-origin binding passed in {scratch_root}")

        # 3. Run all 12 mutants
        killed = 0
        for name, rel_path, target, replacement in MUTANTS:
            _populate_scratch(scratch_root)
            file_path = scratch_root / rel_path
            original = file_path.read_text(encoding="utf-8")
            if target not in original:
                raise RuntimeError(f"{name}: target {target!r} not found in {rel_path}")
            mutated = original.replace(target, replacement, 1)
            if mutated == original:
                raise RuntimeError(f"{name}: replacement produced no change in {rel_path}")
            file_path.write_text(mutated, encoding="utf-8")

            rc, out = _run_pytest_in_scratch(scratch_root)
            has_assertion_error = "AssertionError" in out
            has_bad_error = any(
                err in out for err in ("SyntaxError", "ImportError", "ModuleNotFoundError")
            )
            if rc == 1 and has_assertion_error and not has_bad_error:
                killed += 1
                print(f"  [KILLED] {name} (rc={rc}, AssertionError verified)")
            else:
                print(f"  [SURVIVED/INVALID] {name} (rc={rc}):\n{out}")
                return 1

        print(f"Summary: {killed}/{len(MUTANTS)} mutants killed by AssertionError.")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
