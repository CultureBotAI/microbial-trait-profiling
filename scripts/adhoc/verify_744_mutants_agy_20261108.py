"""Isolated scratch-copy mutant verifier for issue #744 (validation feasibility docs).

Purpose:
  Verify in an isolated temporary copy that the control test suite passes and that
  mutations to methods/validation_feasibility.md and src/mtp/validation_evidence.py
  docstrings (reversions and semantic reversals between the 1,080 all-modality human
  export universe and the 959 paired 16S/WMS discovery universe) are each killed by
  pytest assertions (other_fails == 0).
Inputs:
  Worktree files in the parent repository of this script.
Outputs:
  Stdout summary of control pass and 12/12 mutant kills; non-zero exit on any survivor
  or non-assertion failure.
Issue: #744
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
    "methods/validation_feasibility.md",
    "tests/test_validation_evidence.py",
)

DIRS_TO_COPY = (
    "src/mtp",
    "src/bsdbng",
)

MUTANTS: tuple[tuple[str, str, str, str], ...] = (
    (
        "M01_feasibility_md_drop_study_universe_paragraph",
        "methods/validation_feasibility.md",
        "The inventory reports two nested human study universes side by side (#744):",
        "The inventory reports a single human study universe:",
    ),
    (
        "M02_feasibility_md_swap_1080_and_959_denominators",
        "methods/validation_feasibility.md",
        "all 1,080 human-host studies in `full_dump.csv` (**Human studies** /",
        "all 959 human-host studies in `full_dump.csv` (**Human studies** /",
    ),
    (
        "M03_feasibility_md_drop_non_16s_wms_modalities",
        "methods/validation_feasibility.md",
        "export — including ITS, PCR, 18S or untyped rows — as well as unpaired",
        "export — restricted to 16S and WMS rows — as well as paired",
    ),
    (
        "M04_feasibility_md_drop_denominator_decision_note",
        "methods/validation_feasibility.md",
        (
            "Keeping the whole-export human denominator alongside the\n"
            "discovery-corpus denominator is intentional so validation-feasibility\n"
            "pointers in non-16S/WMS or unpaired human studies remain visible while\n"
            "discovery-corpus counts stay separate; changing the underlying denominator\n"
            "would be a separate scientific decision."
        ),
        "Both columns share the same denominator.",
    ),
    (
        "M05_feasibility_md_drop_marker_split_bucket_paragraph",
        "methods/validation_feasibility.md",
        (
            "As in the study-level summary above, the **Human** column\n"
            '(`bucket["human"]` in `validation_evidence_provenance.json`) counts every'
        ),
        "The **Human** column counts only paired 16S/WMS experiments in",
    ),
    (
        "M06_feasibility_md_corrupt_its_experiment_example",
        "methods/validation_feasibility.md",
        (
            "Non-16S/WMS or\n"
            "unpaired human experiments (such as the ITS experiment `bsdb:684/1` in\n"
            "`marker_named_groups_not_measured`) therefore increment **Human** without\n"
            "being eligible for **In discovery corpus**."
        ),
        "All human experiments in the export are eligible for **In discovery corpus**.",
    ),
    (
        "M07_validation_evidence_py_revert_module_docstring",
        "src/mtp/validation_evidence.py",
        (
            "Both study-level and experiment-level summaries report two nested human\n"
            "evidence universes (#744): all `Homo sapiens` rows in `full_dump.csv`\n"
            "(1,080 human studies across all sequencing modalities — including ITS,\n"
            "PCR, 18S, or untyped rows — and pairing states) alongside the\n"
            "discovery-eligible paired human 16S/WMS subset (959 studies in\n"
            "`study_trait_profiles_wide.tsv`, #742)."
        ),
        "",
    ),
    (
        "M08_validation_evidence_py_revert_evidence_rows_docstring",
        "src/mtp/validation_evidence.py",
        (
            '    """One wide row per human study in `full_dump.csv` (all 1,080 `Homo sapiens`\n'
            "    export studies across all sequencing modalities, including ITS, PCR, 18S or\n"
            "    untyped and unpaired studies), with `in_discovery_corpus` flagging\n"
            "    membership in the 959-study paired human 16S/WMS discovery corpus (#742,\n"
            "    #744).\n"
            '    """'
        ),
        '    """The per-study evidence-availability matrix for human studies."""',
    ),
    (
        "M09_validation_evidence_py_reverse_evidence_rows_universe_meaning",
        "src/mtp/validation_evidence.py",
        "all 1,080 `Homo sapiens`\n    export studies across all sequencing modalities",
        "only the 959 `Homo sapiens`\n    paired 16S/WMS discovery studies",
    ),
    (
        "M10_validation_evidence_py_revert_summary_rows_docstring",
        "src/mtp/validation_evidence.py",
        (
            '    """One count row per evidence type, reporting both `human_studies` /\n'
            "    `human_studies_mentioning` over all `Homo sapiens` export studies (1,080\n"
            "    across all modalities and pairing states) and `discovery_studies` /\n"
            "    `discovery_studies_mentioning` over the paired human 16S/WMS discovery\n"
            "    corpus (959 studies, #742, #744).\n"
            '    """'
        ),
        '    """Studies mentioning each evidence type, overall and within the discovery corpus."""',
    ),
    (
        "M11_validation_evidence_py_revert_marker_split_candidates_docstring",
        "src/mtp/validation_evidence.py",
        (
            " `in_discovery_corpus`\n"
            "    flags paired human 16S/WMS experiments in `study_trait_profiles_wide.tsv`\n"
            "    (#742), whereas whole-export human rows also include non-16S/WMS (e.g. ITS,\n"
            "    PCR, 18S, untyped) and unpaired human experiments (#744)."
        ),
        "",
    ),
    (
        "M12_validation_evidence_py_revert_classification_counts_docstring",
        "src/mtp/validation_evidence.py",
        (
            '    """Per-classification counts for `validation_evidence_provenance.json`,\n'
            "    separating all export experiments (`experiments`), all `Homo sapiens`\n"
            "    export experiments across every sequencing modality and pairing state\n"
            "    (`human`), and paired human 16S/WMS discovery-corpus experiments\n"
            "    (`in_discovery_corpus`, #742, #744).\n"
            '    """'
        ),
        (
            '    """Candidate counts per classification, overall and for human and discovery '
            'subsets."""'
        ),
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
            if key in {"mtp", "bsdbng", "tests.test_validation_evidence"}
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
                        str(root / "tests" / "test_validation_evidence.py"),
                        "-k",
                        (
                            "test_validation_feasibility_docs_and_docstrings_"
                            "distinguish_human_and_discovery_universes"
                        ),
                    ]
                )
            )
        return rc, buf.getvalue()
    finally:
        sys.path[:] = old_path
        os.chdir(old_cwd)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="verify_744_mutants_") as tmp_dir:
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
            "import mtp.validation_evidence as _m\n\n"
            "def pytest_sessionstart(session):\n"
            "    actual = Path(_m.__file__).resolve()\n"
            "    expected = (\n"
            "        Path(__file__).resolve().parents[1]\n"
            '        / "src"\n'
            '        / "mtp"\n'
            '        / "validation_evidence.py"\n'
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
