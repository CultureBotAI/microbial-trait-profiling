"""Isolated scratch-copy mutant verifier for #1334/#2304/#2305 (2026-11-11).

Purpose:
    Verify in an isolated temporary copy (never mutating the live worktree) that
    numerical, publication and versioned-protocol regression tests kill mutants
    via AssertionError. Each follow-up mutant uses only the new regression tests.
Inputs:
    Worktree src/mtp package, the three selected test files, and the two selected
    methods/protocol files. No repository or environment is copied.
Outputs:
    Per-mutant status table printed to stdout; non-zero exit if control fails,
    any mutant survives or fails without AssertionError, or any worktree file changes.
Issue:
    #1334, #2304, #2305
Date:
    2026-11-11
"""

from __future__ import annotations

import hashlib
import importlib
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import NamedTuple

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.dont_write_bytecode = True
TARGET_FILES = (
    Path("src/mtp/native_taxon_influence.py"),
    Path("src/mtp/native_influence_publication.py"),
    Path("src/mtp/influence_assessment.py"),
    Path("tests/test_native_taxon_influence.py"),
    Path("tests/test_native_influence_publication.py"),
    Path("tests/test_influence_assessment.py"),
    Path("methods/multiple_testing.md"),
    Path("methods/table1_native_influence_protocol.json"),
)

TEST_NAMES = (
    "test_constant_paired_effect_is_not_a_significant_native_result",
    "test_native_exports_end_to_end_with_an_empty_stratum",
    "test_unselected_but_testable_traits_are_excluded_from_native_bh_family",
    "test_selected_untestable_traits_are_padded_at_one_in_native_bh_family",
    "test_multiple_testing_methods_doc_covers_native_influence_and_pd_network_bh_families",
)


class _TestPlan(NamedTuple):
    name: str
    test_file: Path
    test_names: tuple[str, ...]
    module_name: str
    expected_passes: int


ORIGINAL = _TestPlan(
    "original",
    Path("tests/test_native_taxon_influence.py"),
    TEST_NAMES,
    "mtp.native_taxon_influence",
    len(TEST_NAMES),
)
PUBLICATION = _TestPlan(
    "publication",
    Path("tests/test_native_influence_publication.py"),
    (
        "test_handoff_distinguishes_correction_family_from_tested_count",
        "test_handoff_refuses_undeclared_correction_family",
    ),
    "mtp.native_influence_publication",
    8,
)
PROTOCOL = _TestPlan(
    "protocol",
    Path("tests/test_native_taxon_influence.py"),
    (
        "test_native_v3_manifest_declares_the_executed_padded_family",
        "test_native_run_refuses_historical_v2_before_writing_outputs",
        "test_native_run_refuses_v3_family_contract_contradictions",
    ),
    "mtp.native_taxon_influence",
    4,
)
READER = _TestPlan(
    "reader",
    Path("tests/test_influence_assessment.py"),
    ("test_native_protocol_versions_preserve_the_recorded_historical_or_current_meaning",),
    "mtp.influence_assessment",
    2,
)

MUTANTS: tuple[tuple[str, Path, str, str], ...] = (
    (
        "M01_bh_over_tested_only",
        Path("src/mtp/native_taxon_influence.py"),
        '    selected_rows = [r for r in records.values() if r["selected"]]',
        '    selected_rows = [r for r in records.values() if r["p"] != ""]',
    ),
    (
        "M02_bh_over_all_records",
        Path("src/mtp/native_taxon_influence.py"),
        '    selected_rows = [r for r in records.values() if r["selected"]]',
        "    selected_rows = list(records.values())",
    ),
    (
        "M03_pad_untestable_with_zero",
        Path("src/mtp/native_taxon_influence.py"),
        'adjust([r["p"] if r["p"] != "" else 1.0 for r in selected_rows])',
        'adjust([r["p"] if r["p"] != "" else 0.0 for r in selected_rows])',
    ),
    (
        "M04_pad_untestable_with_half",
        Path("src/mtp/native_taxon_influence.py"),
        'adjust([r["p"] if r["p"] != "" else 1.0 for r in selected_rows])',
        'adjust([r["p"] if r["p"] != "" else 0.5 for r in selected_rows])',
    ),
    (
        "M05_family_tests_uses_selected_len",
        Path("src/mtp/native_taxon_influence.py"),
        '        row["family_tests"] = len(tested)',
        '        row["family_tests"] = len(selected_rows)',
    ),
    (
        "M06_family_selected_uses_tested_len",
        Path("src/mtp/native_taxon_influence.py"),
        '        row["family_selected_traits"] = len(selected_rows)',
        '        row["family_selected_traits"] = len(tested)',
    ),
    (
        "M07_summary_says_outside_vector",
        Path("src/mtp/native_taxon_influence.py"),
        '"untestable traits at p = 1.0 (q = 1.0) inside the correction vector. "',
        '"untestable traits receive q = 1.0 outside the correction vector. "',
    ),
    (
        "M08_methods_doc_omits_pd_network_rule",
        Path("methods/multiple_testing.md"),
        "passes_min_taxon_studies",
        "omitted_taxon_support_rule",
    ),
    (
        "M09_methods_doc_omits_shared_sign_test",
        Path("methods/multiple_testing.md"),
        "correct_shared_comparisons",
        "omitted_shared_comparisons",
    ),
)

FOLLOWUP_MUTANTS: tuple[tuple[str, Path, str, str, _TestPlan], ...] = (
    (
        "M26_known_v1_protocol_refused",
        Path("src/mtp/native_influence_publication.py"),
        'elif version == "190-native-v1":',
        'elif version == "removed-native-v1":',
        PUBLICATION,
    ),
    (
        "M10_current_family_uses_tested_count",
        Path("src/mtp/native_influence_publication.py"),
        '        family_field = "native_before_family_selected_traits"',
        '        family_field = "native_before_family_tests"',
        PUBLICATION,
    ),
    (
        "M11_historical_family_uses_selected_count",
        Path("src/mtp/native_influence_publication.py"),
        '"native_family_size": arsenate_native[family_field],',
        '"native_family_size": arsenate_native.get("native_before_family_selected_traits", '
        "arsenate_native[family_field]),",
        PUBLICATION,
    ),
    (
        "M12_tested_count_uses_correction_family",
        Path("src/mtp/native_influence_publication.py"),
        '"native_tested_count": arsenate_native["native_before_family_tests"],',
        '"native_tested_count": arsenate_native[family_field],',
        PUBLICATION,
    ),
    (
        "M13_missing_family_silently_falls_back",
        Path("src/mtp/native_influence_publication.py"),
        "        if not record.get(family_field):\n"
        '            raise ValueError(f"Missing correction family field for {version}: '
        '{family_field}")',
        "        if not record.get(family_field):\n"
        '            record[family_field] = record["native_before_family_tests"]',
        PUBLICATION,
    ),
    (
        "M14_unknown_protocol_silently_falls_back",
        Path("src/mtp/native_influence_publication.py"),
        '        raise ValueError(f"Unknown native protocol for correction family: {version!r}")',
        '        family_field = "native_before_family_tests"\n'
        '        family_basis = "historical_finite_tests"',
        PUBLICATION,
    ),
    (
        "M15_current_family_mislabels_its_basis",
        Path("src/mtp/native_influence_publication.py"),
        '        family_basis = "all_selected_traits_padded_at_one"',
        '        family_basis = "historical_finite_tests"',
        PUBLICATION,
    ),
    (
        "M16_protocol_version_reverted",
        Path("methods/table1_native_influence_protocol.json"),
        '"version": "190-native-v3"',
        '"version": "190-native-v2"',
        PROTOCOL,
    ),
    (
        "M17_protocol_denominator_reverted",
        Path("methods/table1_native_influence_protocol.json"),
        '"family_size_field": "family_selected_traits"',
        '"family_size_field": "family_tests"',
        PROTOCOL,
    ),
    (
        "M18_protocol_padding_changed",
        Path("methods/table1_native_influence_protocol.json"),
        '"untestable_correction_p": 1.0',
        '"untestable_correction_p": 0.0',
        PROTOCOL,
    ),
    (
        "M19_protocol_prose_reverted",
        Path("methods/table1_native_influence_protocol.json"),
        "BH over all selected traits per scenario/scope, including untestable selected traits "
        "padded at p=1 inside the correction vector; preserve their reported missing p, "
        "native q=1 and explicit status",
        "BH over finite p-values of all selected traits per scenario/scope",
        PROTOCOL,
    ),
    (
        "M20_historical_protocol_executes_current_code",
        Path("src/mtp/native_taxon_influence.py"),
        "def _validate_protocol(protocol: Record) -> None:\n",
        "def _validate_protocol(protocol: Record) -> None:\n    return\n",
        PROTOCOL,
    ),
    (
        "M21_contradictory_family_field_accepted",
        Path("src/mtp/native_taxon_influence.py"),
        'composition.get("family_size_field") != "family_selected_traits"',
        "False",
        PROTOCOL,
    ),
    (
        "M22_contradictory_padding_accepted",
        Path("src/mtp/native_taxon_influence.py"),
        'composition.get("untestable_correction_p") != 1.0',
        "False",
        PROTOCOL,
    ),
    (
        "M23_reader_rejects_current_protocol",
        Path("src/mtp/influence_assessment.py"),
        'if protocol.get("version") not in {NATIVE_VERSION, "190-native-v3"}:',
        'if protocol.get("version") != NATIVE_VERSION:',
        READER,
    ),
    (
        "M24_reader_rejects_historical_protocol",
        Path("src/mtp/influence_assessment.py"),
        'if protocol.get("version") not in {NATIVE_VERSION, "190-native-v3"}:',
        'if protocol.get("version") != "190-native-v3":',
        READER,
    ),
    (
        "M25_reader_rewrites_historical_protocol",
        Path("src/mtp/influence_assessment.py"),
        'detail += _detail("protocol", "native_influence", protocol)',
        'detail += _detail("protocol", "native_influence", '
        '{**protocol, "version": "190-native-v3"})',
        READER,
    ),
)


class _OutcomeCapturePlugin:
    def __init__(
        self, expected_test_file: Path, expected_module_file: Path, module_name: str
    ) -> None:
        self.expected_test_file = expected_test_file.resolve()
        self.expected_module_file = expected_module_file.resolve()
        self.module_name = module_name
        self.passed = 0
        self.assertion_failures = 0
        self.other_failures: list[str] = []
        self.origin_checked = False

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        for item in items:
            item_path = Path(str(item.fspath)).resolve()
            if item_path != self.expected_test_file:
                raise RuntimeError(f"Unexpected test origin: {item_path}")
        mod = sys.modules.get(self.module_name)
        if mod is None or mod.__file__ is None:
            raise RuntimeError(f"{self.module_name} not loaded during collection")
        if Path(mod.__file__).resolve() != self.expected_module_file:
            raise RuntimeError(f"Unexpected module origin: {mod.__file__}")
        self.origin_checked = True

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call" and report.passed:
            self.passed += 1
        elif report.failed:
            longrepr = getattr(report, "longrepr", None)
            crash = getattr(longrepr, "reprcrash", None)
            msg = getattr(crash, "message", str(longrepr))
            if report.when == "call" and (
                msg.startswith(("AssertionError", "Failed: DID NOT RAISE", "assert "))
                or "AssertionError" in msg
            ):
                self.assertion_failures += 1
            else:
                self.other_failures.append(f"{report.nodeid}:{report.when}:{msg}")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _snapshot_hashes(root: Path) -> dict[Path, str]:
    return {rel: _sha256(root / rel) for rel in TARGET_FILES}


def _populate_scratch(scratch: Path) -> None:
    shutil.copytree(
        REPO_ROOT / "src" / "mtp",
        scratch / "src" / "mtp",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    for relative in TARGET_FILES:
        if relative.parts[0] != "src":
            target = scratch / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(REPO_ROOT / relative, target)


def _run_in_scratch(scratch: Path, plan: _TestPlan) -> tuple[int, _OutcomeCapturePlugin]:
    os.environ["GIT_CONFIG_GLOBAL"] = "/dev/null"
    os.environ["GIT_CONFIG_NOSYSTEM"] = "1"
    scratch_src = str((scratch / "src").resolve())
    repo_src = str((REPO_ROOT / "src").resolve())
    sys.path[:] = [p for p in sys.path if Path(p).resolve() != Path(repo_src)]
    if scratch_src not in sys.path:
        sys.path.insert(0, scratch_src)

    test_modules = {path.stem for path in TARGET_FILES if path.parts[0] == "tests"}
    for mod_name in list(sys.modules):
        if (
            mod_name == "mtp"
            or mod_name.startswith("mtp.")
            or mod_name in test_modules
            or any(mod_name.endswith(f".{name}") for name in test_modules)
        ):
            sys.modules.pop(mod_name, None)
    importlib.invalidate_caches()

    test_path = (scratch / plan.test_file).resolve()
    mod_path = (scratch / "src" / Path(*plan.module_name.split("."))).with_suffix(".py")
    plugin = _OutcomeCapturePlugin(test_path, mod_path, plan.module_name)
    rc = int(
        pytest.main(
            [
                "-q",
                "-p",
                "no:cacheprovider",
                "--basetemp",
                str(scratch / "pytest-temp"),
                "-o",
                "addopts=--import-mode=importlib",
                "-k",
                " or ".join(plan.test_names),
                str(test_path),
            ],
            plugins=[plugin],
        )
    )
    return rc, plugin


def main() -> int:
    before = _snapshot_hashes(REPO_ROOT)
    with tempfile.TemporaryDirectory(prefix="verify_1334_mutants_") as tmp_dir:
        scratch = Path(tmp_dir)
        _populate_scratch(scratch)
        for plan in (ORIGINAL, PUBLICATION, PROTOCOL, READER):
            rc, control = _run_in_scratch(scratch, plan)
            if (
                rc != 0
                or not control.origin_checked
                or control.passed != plan.expected_passes
                or control.assertion_failures != 0
                or control.other_failures
            ):
                print(f"CONTROL {plan.name} FAILED:", rc, control.passed, control.other_failures)
                return 1
            print(f"CONTROL {plan.name}: rc=0 PYTEST_ORIGIN_OK ({control.passed} passed)")

        mutants: list[tuple[str, Path, str, str, _TestPlan]] = [
            (*mutant, ORIGINAL) for mutant in MUTANTS
        ]
        mutants.extend(FOLLOWUP_MUTANTS)
        killed = 0
        for name, rel_path, target, replacement, plan in mutants:
            file_path = scratch / rel_path
            original = file_path.read_text(encoding="utf-8")
            if original.count(target) != 1:
                print(f"{name}: target count {original.count(target)} != 1 in {rel_path}")
                return 1
            file_path.write_text(original.replace(target, replacement, 1), encoding="utf-8")
            try:
                m_rc, outcome = _run_in_scratch(scratch, plan)
            finally:
                file_path.write_text(original, encoding="utf-8")
            if (
                m_rc != 0
                and outcome.origin_checked
                and outcome.assertion_failures >= 1
                and not outcome.other_failures
            ):
                killed += 1
                print(
                    f"{name}: KILLED (rc={m_rc}, assertion_failures={outcome.assertion_failures}, "
                    f"test_plan={plan.name})"
                )
            else:
                print(
                    f"{name}: SURVIVED/OTHER (rc={m_rc}, "
                    f"assertions={outcome.assertion_failures}, other={outcome.other_failures})"
                )
                return 1

    after = _snapshot_hashes(REPO_ROOT)
    if before != after:
        print("WORKTREE HASH MISMATCH:", before, after)
        return 1
    print(f"SUMMARY: {killed}/{len(mutants)} mutants killed by AssertionError; worktree untouched.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
