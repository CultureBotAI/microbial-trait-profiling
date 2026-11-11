#!/usr/bin/env python3
"""Scratch-copy mutation verifier for #1274 and #1022 taxa-coverage follow-ups.

Purpose:
  Verify that tests/test_taxa_coverage_taxonomy.py and tests/test_taxa_coverage_histogram.py
  kill source-level mutations of scripts/create_taxa_coverage_histogram.py (#1274, #1022)
  via AssertionError (rc=1) when executed against an isolated scratch copy with verified
  import origin inside pytest.

Inputs:
  - scripts/create_taxa_coverage_histogram.py
  - tests/test_taxa_coverage_taxonomy.py
  - tests/test_taxa_coverage_histogram.py

Outputs:
  - stdout summary of CONTROL (rc=0) and 6 killed mutants (rc=1, AssertionError,
    PYTEST_ORIGIN_OK).

Issue: #1274, #1022
Date: 2026-11-08
"""

from __future__ import annotations

import contextlib
import io
import multiprocessing
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
TARGET_REL = Path("scripts/create_taxa_coverage_histogram.py")
TESTS_REL = (
    Path("tests/test_taxa_coverage_taxonomy.py"),
    Path("tests/test_taxa_coverage_histogram.py"),
)

ASSERTION_MARKERS = (
    "AssertionError",
    "Failed: DID NOT RAISE",
    "Regex pattern did not match",
)
ERROR_MARKERS = (
    "ImportError",
    "ModuleNotFoundError",
    "SyntaxError",
    "IndentationError",
    "NameError",
    "TypeError",
    "AttributeError",
    "KeyError",
    "INTERNALERROR",
    "ERROR collecting",
)


@dataclass(frozen=True)
class Mutant:
    """One source-level mutation in scripts/create_taxa_coverage_histogram.py."""

    name: str
    old: str
    new: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M1_silently_skip_numeric_subject_pointing_at_rank_root_1022",
        old=(
            '                if subject.startswith("NCBITaxon:"):\n'
            "                    raise ValueError("
            "f\"Invalid taxonomy identifier in edge {row['id']}: {subject}\")\n"
        ),
        new="",
    ),
    Mutant(
        name="M2_silently_skip_malformed_ncbi_identifier_1022",
        old=(
            "            for identifier in (subject, parent):\n"
            '                if identifier.startswith("NCBITaxon:") and not re.fullmatch(\n'
            '                    r"NCBITaxon:[1-9][0-9]*", identifier\n'
            "                ):\n"
            "                    raise ValueError(\n"
            "                        f\"Invalid taxonomy identifier in edge {row['id']}: "
            '{identifier}"\n'
            "                    )\n"
        ),
        new=(
            "            if any(\n"
            '                identifier.startswith("NCBITaxon:")\n'
            '                and not re.fullmatch(r"NCBITaxon:[1-9][0-9]*", identifier)\n'
            "                for identifier in (subject, parent)\n"
            "            ):\n"
            "                continue\n"
        ),
    ),
    Mutant(
        name="M3_accept_rank_class_only_file_1022",
        old="        if not accepted:\n",
        new="        if not accepted and not rank_class_records:\n",
    ),
    Mutant(
        name="M4_omit_excluded_records_note_1274",
        old=(
            "            excluded_records_note=(\n"
            '                "excluded_records counts malformed KG edge records dropped by '
            'the strict "\n'
            '                "edge parser (always 0); taxonomy_rank_class_edge_records_skipped '
            'counts "\n'
            '                "well-formed ontology rank-class subclass_of records skipped '
            'during "\n'
            '                "organism-ancestry graph projection (#1008, #1274)."\n'
            "            ),\n"
        ),
        new='            excluded_records_note="",\n',
    ),
    Mutant(
        name="M5_revert_readme_rank_class_wording_1274",
        old=(
            '            "Ontology rank-class edges are skipped during organism-ancestry '
            'projection",\n'
            '            "(counted in parser_audit as taxonomy_rank_class_edge_records_skipped, '
            'while",\n'
            '            "excluded_records: 0 confirms the strict KG edge parser dropped no '
            'record).",\n'
        ),
        new=(
            '            "Ontology rank-class edges are excluded from organism ancestry and '
            'counted in",\n'
            '            "the provenance audit as taxonomy_rank_class_edge_records_skipped.",\n'
        ),
    ),
    Mutant(
        name="M6_misdiagnose_candidatus_lineage_1274",
        old='            "Candidatus Lineage",\n',
        new='            "Rare / Novel Isolate (Lacks KG Phenotype)",\n',
    ),
)


def _setup_scratch(tmp_dir: Path) -> None:
    """Copy minimal source, scripts, fixtures, and test trees into the scratch directory."""
    for rel_dir in ("src/mtp", "src/bsdbng", "tests/fixtures/coverage"):
        src_dir = REPO_ROOT / rel_dir
        if src_dir.exists():
            shutil.copytree(src_dir, tmp_dir / rel_dir)
    for rel_file in (
        TARGET_REL,
        Path("scripts/generate_tufte_plates.py"),
        *TESTS_REL,
        Path("pyproject.toml"),
        Path("uv.lock"),
        Path(".git"),
    ):
        dst = tmp_dir / rel_file
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO_ROOT / rel_file, dst)
    conftest = tmp_dir / "tests/conftest.py"
    expected_root = tmp_dir.resolve()
    conftest.write_text(
        "from pathlib import Path\n"
        "import mtp.paths as _paths\n"
        "from test_taxa_coverage_histogram import coverage as _cov\n\n"
        "def pytest_runtest_setup(item):\n"
        "    actual = _paths.REPO_ROOT.resolve()\n"
        f"    expected = Path({str(expected_root)!r}).resolve()\n"
        "    if actual != expected:\n"
        "        raise RuntimeError(f'wrong origin: {actual} != {expected}')\n"
        "    cov_path = Path(_cov.__file__).resolve()\n"
        "    if cov_path != expected / 'scripts/create_taxa_coverage_histogram.py':\n"
        "        raise RuntimeError(f'wrong coverage origin: {cov_path}')\n"
        "    print('PYTEST_ORIGIN_OK', flush=True)\n",
        encoding="utf-8",
    )


def _worker(tmp_dir_str: str, out_file_str: str) -> None:
    """Run pytest in a clean child process rooted at tmp_dir."""
    tmp_dir = Path(tmp_dir_str)
    os.chdir(tmp_dir)
    sys.path.insert(0, str(tmp_dir / "tests"))
    sys.path.insert(0, str(tmp_dir / "src"))
    for name in [
        k
        for k in list(sys.modules)
        if k in {"mtp", "test_taxa_coverage_histogram", "taxa_coverage_histogram"}
        or k.startswith(("mtp.", "bsdbng"))
    ]:
        sys.modules.pop(name, None)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    "-q",
                    "-s",
                    "--tb=short",
                    "--rootdir",
                    str(tmp_dir),
                    str(tmp_dir / "tests/test_taxa_coverage_taxonomy.py"),
                    f"{tmp_dir / 'tests/test_taxa_coverage_histogram.py'}"
                    "::test_ncbi_bacterial_zero_trait_categories_and_parser_audit_wording",
                ]
            )
        )
    Path(out_file_str).write_text(f"RC={rc}\n{buf.getvalue()}", encoding="utf-8")


def _run_case(tmp_dir: Path) -> tuple[int, str]:
    out_file = tmp_dir / "pytest_out.txt"
    ctx = multiprocessing.get_context("spawn")
    proc = ctx.Process(target=_worker, args=(str(tmp_dir), str(out_file)))
    proc.start()
    proc.join()
    if proc.exitcode != 0 or not out_file.exists():
        raise RuntimeError(f"worker process failed with exitcode={proc.exitcode}")
    raw = out_file.read_text(encoding="utf-8")
    first, _, rest = raw.partition("\n")
    return int(first.removeprefix("RC=")), rest


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="mtp_1274_control_") as tmp:
        tmp_dir = Path(tmp)
        _setup_scratch(tmp_dir)
        rc, out = _run_case(tmp_dir)
        if rc != 0 or "PYTEST_ORIGIN_OK" not in out:
            print(f"FAIL [CONTROL]: rc={rc}\n{out}")
            return 1
        print("CONTROL: rc=0 PYTEST_ORIGIN_OK")

    for mutant in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mtp_1274_{mutant.name}_") as tmp:
            tmp_dir = Path(tmp)
            _setup_scratch(tmp_dir)
            target_path = tmp_dir / TARGET_REL
            source = target_path.read_text(encoding="utf-8")
            count = source.count(mutant.old)
            if count != 1:
                print(f"FAIL [{mutant.name}]: expected 1 occurrence, found {count}")
                return 1
            target_path.write_text(source.replace(mutant.old, mutant.new, 1), encoding="utf-8")
            rc, out = _run_case(tmp_dir)
            if rc != 1:
                print(f"FAIL [{mutant.name}]: expected rc=1, got rc={rc}\n{out}")
                return 1
            if "PYTEST_ORIGIN_OK" not in out:
                print(f"FAIL [{mutant.name}]: missing PYTEST_ORIGIN_OK\n{out}")
                return 1
            if any(err in out for err in ERROR_MARKERS):
                print(f"FAIL [{mutant.name}]: non-assertion error marker in output\n{out}")
                return 1
            if not any(marker in out for marker in ASSERTION_MARKERS):
                print(f"FAIL [{mutant.name}]: no assertion marker in output\n{out}")
                return 1
            print(f"KILLED {mutant.name}: rc=1 AssertionError PYTEST_ORIGIN_OK")

    print(f"ALL {len(MUTANTS)} #1274/#1022 MUTANTS KILLED CLEANLY.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
