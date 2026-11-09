"""Isolated scratch-copy mutant verifier for issue #568 (calibration protocol & ORA urn docs).

Purpose:
  Verify in an isolated temporary copy that the control test suite passes and that
  mutations to methods/calibration_protocol.json, src/mtp/calibration.py,
  methods/calibration.md, methods/three_arm_statistics.md, and
  methods/table1_native_influence.md are each killed by pytest assertions
  (other_fails == 0).
Inputs:
  Worktree files in the parent repository of this script.
Outputs:
  Stdout summary of control pass and 12/12 mutant kills; non-zero exit on any survivor
  or non-assertion failure.
Issue: #568
Date: 2026-11-08
"""

from __future__ import annotations

import contextlib
import io
import os
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

WORKTREE = Path(__file__).resolve().parents[2]

FILES_TO_COPY = (
    "pyproject.toml",
    "methods/calibration_protocol.json",
    "methods/calibration.md",
    "methods/three_arm_statistics.md",
    "methods/table1_native_influence.md",
    "tests/test_calibration.py",
)

DIRS_TO_COPY = (
    "src/mtp",
    "src/bsdbng",
)

MUTANTS: tuple[tuple[str, str, str, str], ...] = (
    (
        "M01_protocol_revert_declared_by_568",
        "methods/calibration_protocol.json",
        (
            "and its scope and not_covered declarations were amended in #568 to record "
            "the three ORA background/urn definitions and the separate arm-calibration "
            "harnesses (#550, #551). "
        ),
        "",
    ),
    (
        "M02_protocol_revert_scope_to_old_131_owner_wording",
        "methods/calibration_protocol.json",
        (
            "procedures owned by this repository's directional-consistency framework "
            "(src/mtp/calibration.py). As amended in #568, the enrichment arm's "
            "production hypergeometric procedure draws distinct reported taxa against "
            "the all-host taxa_trait_matrix.tsv background (M = 5,829 BugSigDB-reported "
            "taxa, including 645 non-human-only and 862 non-paired-human taxa; "
            "src/mtp/three_test_table/ora.py, calibrated separately in "
            "src/mtp/ora_calibration/ under #551), whereas the historical #131 "
            "assertion-count urn (106,993,585 assertions across 4,517 testable profiles) "
            "is retained only in the Table 1 native influence ORA reruns "
            "(src/mtp/native_taxon_influence.py, methods/table1_native_influence.md) and "
            "the harmonized baseline in methods/three_arm_statistics.md defines a "
            "proposed paired-human assertion-weighted background; the composition arm's "
            "native transformed procedure is calibrated separately in "
            "src/mtp/composition_calibration/ (#550)."
        ),
        (
            "procedures owned by this repository's directional-consistency framework. "
            "The composition arm's and the enrichment arm's native procedures are not "
            "evaluated here; their adapters belong to their owners (#182 assigns them), "
            "and the assertion-level urn of #131 is a separate model question."
        ),
    ),
    (
        "M03_protocol_corrupt_production_taxon_count",
        "methods/calibration_protocol.json",
        (
            "M = 5,829 BugSigDB-reported taxa, including 645 non-human-only and 862 "
            "non-paired-human taxa"
        ),
        (
            "M = 4,967 BugSigDB-reported taxa, including 645 non-human-only and 862 "
            "non-paired-human taxa"
        ),
    ),
    (
        "M04_protocol_revert_not_covered_0",
        "methods/calibration_protocol.json",
        (
            "the enrichment arm's hypergeometric procedures — production distinct-taxon "
            "urn over taxa_trait_matrix.tsv (M = 5,829 all-host taxa, including 645 "
            "non-human-only and 862 non-paired-human taxa; calibrated separately in "
            "src/mtp/ora_calibration/, #551), historical #131 assertion-count urn "
            "retained in the Table 1 native influence reruns (106,993,585 assertions "
            "across 4,517 testable profiles; methods/table1_native_influence.md), and "
            "proposed harmonized-baseline paired-human assertion-weighted background "
            "(methods/three_arm_statistics.md) — and the composition arm's native "
            "transformed procedure (calibrated separately in "
            "src/mtp/composition_calibration/, #550)"
        ),
        (
            "the enrichment arm's hypergeometric urn (#131) and the composition arm's "
            "native transform (owners' adapters)"
        ),
    ),
    (
        "M05_protocol_corrupt_amendment_date_or_reason",
        "methods/calibration_protocol.json",
        '"date": "2026-11-08"',
        '"date": "2026-10-19"',
    ),
    (
        "M06_calibration_py_revert_docstring",
        "src/mtp/calibration.py",
        "As recorded in the v3.1 (#568) amendment to methods/calibration_protocol.json,",
        ("The composition and enrichment arms' native procedures need adapters from their owners,"),
    ),
    (
        "M07_calibration_py_revert_not_calibrated_here_tuple",
        "src/mtp/calibration.py",
        (
            '    "the enrichment arm\x27s hypergeometric procedures (production M = 5,829 "\n'
            '    "all-host distinct-taxon urn in src/mtp/ora_calibration/ #551; historical "\n'
            '    "#131 assertion-count urn in Table 1 native influence reruns; proposed "\n'
            '    "harmonized-baseline paired-human assertion-weighted background #568)",'
        ),
        '    "the enrichment arm\x27s hypergeometric urn (#131)",',
    ),
    (
        "M08_calibration_py_status_hardcodes_old_not_calibrated_here",
        "src/mtp/calibration.py",
        '"not_calibrated_here": list(NOT_CALIBRATED_HERE),',
        '"not_calibrated_here": ["the enrichment arm\x27s hypergeometric urn (#131)"],',
    ),
    (
        "M09_calibration_md_revert_ora_urn_bullet",
        "methods/calibration.md",
        (
            "The enrichment arm's hypergeometric procedures across all three documented ORA\n"
            "  background definitions (`methods/three_arm_statistics.md`, #568):"
        ),
        "The enrichment arm's hypergeometric procedure; its adapter belongs to its owner:",
    ),
    (
        "M10_calibration_md_drop_v31_amendment_note",
        "methods/calibration.md",
        (
            "[`calibration_protocol.json`](calibration_protocol.json) records this scope\n"
            "  amendment explicitly (`v3.1`, #568) rather than silently rewriting its earlier\n"
            "  #131 assertion-level and #182 owner-adapter declaration."
        ),
        "[`calibration_protocol.json`](calibration_protocol.json) records the criteria.",
    ),
    (
        "M11_three_arm_statistics_drop_calibration_protocol_crossref",
        "methods/three_arm_statistics.md",
        (
            " [`calibration_protocol.json`](calibration_protocol.json) (`v3.1`, #568) and "
            "[`calibration.md`](calibration.md) explicitly record how these three ORA "
            "background definitions relate to the directional calibration scope and the "
            "separate ORA calibration harness (`src/mtp/ora_calibration/`, #551)."
        ),
        "",
    ),
    (
        "M12_table1_native_influence_drop_calibration_protocol_crossref",
        "methods/table1_native_influence.md",
        (
            " [`calibration_protocol.json`](calibration_protocol.json) (`v3.1`, #568) and "
            "[`calibration.md`](calibration.md) record this distinction between the "
            "historical #131 assertion-count urn retained here, the production `M = 5,829` "
            "all-host distinct-taxon urn calibrated in `src/mtp/ora_calibration/` (#551), "
            "and the proposed harmonized-baseline paired-human assertion-weighted background."
        ),
        "",
    ),
    (
        "M13_bsdbng_calibration_shim_omit_not_calibrated_here",
        "src/bsdbng/calibration.py",
        '    "NOT_CALIBRATED_HERE",\n',
        "",
    ),
)


def _run_pytest(root: Path) -> tuple[int, str]:
    old_cwd = Path.cwd()
    old_path = list(sys.path)
    os.environ["GIT_CONFIG_GLOBAL"] = "/dev/null"
    os.environ["GIT_CONFIG_NOSYSTEM"] = "1"
    try:
        os.chdir(root)
        sys.path.insert(0, str(root / "src"))
        for mod_name in [
            key
            for key in list(sys.modules)
            if key in {"mtp", "bsdbng", "tests.test_calibration"}
            or key.startswith(("mtp.", "bsdbng.", "tests."))
        ]:
            sys.modules.pop(mod_name, None)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            rc = int(
                pytest.main(
                    [
                        "-s",
                        "--tb=short",
                        "--rootdir",
                        str(root),
                        str(root / "tests" / "test_calibration.py"),
                        "-k",
                        (
                            "test_calibration_protocol_and_docs_declare_ora_urn_scope_amendment or "
                            "test_calibration_summary_marks_p_value_columns_as_uncorrected_diagnostics"
                        ),
                    ]
                )
            )
        return rc, buf.getvalue()
    finally:
        sys.path[:] = old_path
        os.chdir(old_cwd)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="verify_568_mutants_") as tmp_dir:
        root = Path(tmp_dir)
        for rel_dir in DIRS_TO_COPY:
            shutil.copytree(WORKTREE / rel_dir, root / rel_dir)
        for rel in FILES_TO_COPY:
            src = WORKTREE / rel
            dst = root / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        conftest = (
            "from pathlib import Path\n"
            "import mtp.calibration as _m\n\n"
            "def pytest_sessionstart(session):\n"
            "    actual = Path(_m.__file__).resolve()\n"
            "    expected = (\n"
            '        Path(__file__).resolve().parents[1] / "src" / "mtp" / "calibration.py"\n'
            "    ).resolve()\n"
            "    if actual != expected:\n"
            '        raise RuntimeError(f"Import origin mismatch: {actual} != {expected}")\n'
            '    print(f"PYTEST_ORIGIN_OK:{actual}")\n'
        )
        (root / "tests" / "conftest.py").write_text(conftest, encoding="utf-8")

        rc, control_out = _run_pytest(root)
        if rc != 0 or "PYTEST_ORIGIN_OK:" not in control_out:
            sys.stderr.write(f"CONTROL FAILED (rc={rc}):\n{control_out}\n")
            return 1
        print("control: PASS (rc=0, origin verified inside pytest)")

        killed = 0
        other_fails = 0
        survivors: list[str] = []

        for name, rel_path, target, replacement in MUTANTS:
            file_path = root / rel_path
            original = file_path.read_text(encoding="utf-8")
            count = original.count(target)
            if count != 1:
                sys.stderr.write(f"{name}: expected target once in {rel_path}, found {count}\n")
                return 2
            try:
                file_path.write_text(original.replace(target, replacement, 1), encoding="utf-8")
                m_rc, output = _run_pytest(root)
                if m_rc == 0:
                    survivors.append(name)
                    print(f"  {name}: SURVIVED")
                elif "PYTEST_ORIGIN_OK:" in output and (
                    "AssertionError" in output or " assert " in output
                ):
                    killed += 1
                    print(f"  {name}: KILLED (AssertionError)")
                else:
                    other_fails += 1
                    print(f"  {name}: OTHER_FAIL (rc={m_rc})\n{output}")
            finally:
                file_path.write_text(original, encoding="utf-8")

        total = len(MUTANTS)
        print(
            f"summary: killed={killed}/{total}, survivors={len(survivors)}, "
            f"other_fails={other_fails}"
        )
        if killed != total or survivors or other_fails:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
