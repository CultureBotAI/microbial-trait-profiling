#!/usr/bin/env python3
"""Scratch-copy mutation verifier for #1265 and #1322 three_arm_statistics.md contracts.

Purpose:
    Verify that single-site revert mutations in ``methods/three_arm_statistics.md``
    are killed by ``AssertionError`` in
    ``tests/test_biological_interpretation.py::test_three_arm_statistics_pins_common_quantity_modules_and_historical_screen_provenance``
    without modifying any tracked repository file.

Inputs:
    Tracked ``methods/three_arm_statistics.md`` and
    ``tests/test_biological_interpretation.py``.

Outputs:
    Per-case ``PASS``/``FAIL`` lines and exit code ``0`` when the unmutated
    control passes and all 8 mutants fail by assertion.

Issues: #1265, #1322
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
    "methods/three_arm_statistics.md",
    "tests/test_biological_interpretation.py",
)


@dataclass(frozen=True)
class DocMutant:
    """Single-site revert mutation on ``methods/three_arm_statistics.md``."""

    name: str
    rel_path: str
    target: str
    replacement: str


MUTANTS: tuple[DocMutant, ...] = (
    DocMutant(
        name="M1_intro_reverts_to_to_implement_in_178_1265",
        rel_path="methods/three_arm_statistics.md",
        target=(
            "This section specifies the harmonized comparison implemented in\n"
            "[`src/mtp/integrated_evidence.py`](../src/mtp/integrated_evidence.py)"
        ),
        replacement=(
            "This section specifies the harmonized comparison to implement in\n"
            "[#178](https://github.com/CultureBotAI/microbial-trait-profiling/issues/178)"
        ),
    ),
    DocMutant(
        name="M2_fraction_denominator_drops_8b280a3_pin_1322",
        rel_path="methods/three_arm_statistics.md",
        target="344 under `phenotype` in the historical scope audit at main 8b280a3.",
        replacement="344 under `phenotype` in the scope audit.",
    ),
    DocMutant(
        name="M3_denominators_bullet_drops_historical_8b280a3_pin_1322",
        rel_path="methods/three_arm_statistics.md",
        target=(
            "- **Denominators.** In that historical 8b280a3 audit, every profile's fraction sums "
            "344 retained columns."
        ),
        replacement="- **Denominators.** Every profile's fraction sums 344 retained columns.",
    ),
    DocMutant(
        name="M4_sensitivity_reverts_to_tracks_planned_183_1265",
        rel_path="methods/three_arm_statistics.md",
        target=(
            "[`src/mtp/log_ratio_sensitivity.py`](../src/mtp/log_ratio_sensitivity.py) "
            "([#183](https://github.com/CultureBotAI/microbial-trait-profiling/issues/183); "
            "[`table1_log_ratio_sensitivity.md`](table1_log_ratio_sensitivity.md)) "
            "implements this\nsensitivity analysis."
        ),
        replacement=(
            "[#183](https://github.com/CultureBotAI/microbial-trait-profiling/issues/183) tracks "
            "this planned\nsensitivity analysis."
        ),
    ),
    DocMutant(
        name="M5_typical_n_drops_historical_screen_and_shift_table_pins_1322",
        rel_path="methods/three_arm_statistics.md",
        target=(
            "across the 317 `off`-scope candidates in the historical v1.3.1 screen "
            "([`trait_direction_screen.md`](trait_direction_screen.md); [Vintage](#vintage)), "
            "606 to 898 for the six of Table 1 in the historical `off`-scope shift table "
            "([`table1_trait_condition_body_site_maps.md`]"
            "(table1_trait_condition_body_site_maps.md), ties excluded)"
        ),
        replacement=(
            "across the 317 `off`-scope candidates, 606 to 898 for the six of Table 1 "
            "(ties excluded)"
        ),
    ),
    DocMutant(
        name="M6_main_limitations_drops_historical_screen_pin_1322",
        rel_path="methods/three_arm_statistics.md",
        target=(
            "across the 317 `off`-scope candidates in that historical screen "
            "([`trait_direction_screen.md`](trait_direction_screen.md); [Vintage](#vintage))"
        ),
        replacement="across the 317 `off`-scope candidates",
    ),
    DocMutant(
        name="M7_harmonization_reverts_to_until_work_is_complete_1265",
        rel_path="methods/three_arm_statistics.md",
        target=(
            "[`src/mtp/composition_denominators.py`](../src/mtp/composition_denominators.py) "
            "([#174](https://github.com/CultureBotAI/microbial-trait-profiling/issues/174)), "
            "`drop_zero_total_pairs` in [`src/mtp/trait_guilds.py`](../src/mtp/trait_guilds.py) "
            "(#228), [`src/mtp/log_ratio_sensitivity.py`](../src/mtp/log_ratio_sensitivity.py) "
            "([#183](https://github.com/CultureBotAI/microbial-trait-profiling/issues/183)), and "
            "[`src/mtp/integrated_evidence.py`](../src/mtp/integrated_evidence.py) "
            "([#178](https://github.com/CultureBotAI/microbial-trait-profiling/issues/178)) "
            "implement the shared non-MediaDive denominator and common paired quantities on "
            "`main`; the historical results below retain their original run-specific\n"
            "denominator definitions from the snapshots in which they were produced."
        ),
        replacement=(
            "Issue\n"
            "[#174](https://github.com/CultureBotAI/microbial-trait-profiling/issues/174) tracks "
            "harmonization and\nthe numerical rerun; the historical results below must retain "
            "their original\ndenominator definitions until that work is complete."
        ),
    ),
    DocMutant(
        name="M8_noel_integrated_evidence_drops_module_pin_1265",
        rel_path="methods/three_arm_statistics.md",
        target=(
            "Noel owns the integrated evidence table in\n"
            "[#178](https://github.com/CultureBotAI/microbial-trait-profiling/issues/178), "
            "implemented in [`src/mtp/integrated_evidence.py`](../src/mtp/integrated_evidence.py) "
            "([`integrated_evidence.md`](integrated_evidence.md))."
        ),
        replacement=(
            "Noel owns the integrated evidence table in\n"
            "[#178](https://github.com/CultureBotAI/microbial-trait-profiling/issues/178)."
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
    test_file = scratch_root / "tests" / "test_biological_interpretation.py"
    buf = io.StringIO()
    sys.modules.pop("test_biological_interpretation", None)
    with redirect_stdout(buf), redirect_stderr(buf):
        rc = int(
            pytest.main(
                [
                    str(test_file),
                    "-k",
                    "test_three_arm_statistics_pins_common_quantity_modules_and_historical_screen_provenance",
                    "--tb=long",
                    "-q",
                ]
            )
        )
    return rc, buf.getvalue()


def main() -> int:
    before = _tracked_hashes()
    failures = 0

    with tempfile.TemporaryDirectory(prefix="mtp-1265-control-") as tmp:
        scratch = Path(tmp)
        _copy_tree(scratch)
        rc, out = _run_scratch_pytest(scratch)
        if rc != 0:
            print(f"FAIL [CONTROL]: expected rc=0, got rc={rc}\n{out}")
            return 1
        print("PASS [CONTROL]: rc=0, PYTEST_ORIGIN_OK")

    for mutant in MUTANTS:
        with tempfile.TemporaryDirectory(prefix=f"mtp-1265-{mutant.name}-") as tmp:
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
