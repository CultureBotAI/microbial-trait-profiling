"""Check the two stool-selection regression mutants (#1142, 2026-10-28).

Inputs: the synthetic exact-stool integration test and current production code.
Outputs: per-worker pytest logs, source hashes and JSON outcomes under --output.
Mutations replace Python callables only in isolated spawned processes; production
files remain unchanged. --expect-survivors records the pre-fix coverage gap.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import multiprocessing
import sys
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
TEST = (
    "tests/test_six_shared_trait_list_contrasts.py::"
    "test_stool_group_selection_is_additive_and_uses_exact_site_after_corpus_zscores"
)
MUTANTS = ("baseline", "old_grouped_stool", "profile_count_ranking")
SOURCES = (
    "src/mtp/stool_group_contrasts.py",
    "scripts/create_six_shared_trait_sankeys.py",
    "scripts/create_six_shared_trait_cytoscape_heatmaps.py",
    "tests/test_six_shared_trait_list_contrasts.py",
    "scripts/adhoc/issue1142_stool_selection_mutants.py",
)


class Reports:
    def __init__(self) -> None:
        self.counts: Counter[str] = Counter()
        self.failures: list[str] = []

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        if report.when == "call":
            self.counts[report.outcome] += 1
        if report.failed:
            self.failures.append(report.nodeid)


def _hashes() -> dict[str, str]:
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in SOURCES}


def worker(mutant: str, output: Path) -> int:
    before = _hashes()
    with pytest.MonkeyPatch.context() as patch:
        patch.syspath_prepend(str(ROOT / "scripts"))
        heatmaps = importlib.import_module("create_six_shared_trait_cytoscape_heatmaps")
        prepare = heatmaps.prepare_contrast_cells
        selector = heatmaps.select_stool_group_contrasts
        grouped_calls = 0
        corpus_cells: pd.DataFrame | None = None

        def old_population(*args: Any, **kwargs: Any) -> pd.DataFrame:
            nonlocal grouped_calls, corpus_cells
            cells: pd.DataFrame = prepare(*args, **kwargs)
            if args[2].key == "grouped":
                grouped_calls += 1
                if grouped_calls % 2:
                    corpus_cells = cells.copy()
                else:
                    if corpus_cells is None:
                        raise AssertionError("Mutant needs the preceding full-corpus cells")
                    return corpus_cells.loc[corpus_cells["body_site_group"].eq("Feces")].copy()
            return cells

        def rank_profiles(
            cells: pd.DataFrame, *, group_count: int = 6
        ) -> tuple[pd.DataFrame, pd.DataFrame]:
            all_selected, ledger = selector(cells, group_count=max(1, len(cells)))
            profiles = (
                cells.groupby("condition", sort=False)[
                    ["increased_profile_count", "decreased_profile_count"]
                ]
                .first()
                .sum(axis=1)
            )
            ledger["ranking_profiles"] = ledger["condition"].map(profiles)
            ledger = ledger.sort_values(
                ["eligible", "ranking_profiles", "condition"], ascending=[False, False, True]
            ).reset_index(drop=True)
            ledger["rank"] = pd.array(
                [
                    index + 1 if eligible else None
                    for index, eligible in enumerate(ledger["eligible"])
                ],
                dtype="Int64",
            )
            ledger["chosen"] = ledger["rank"].le(group_count).fillna(False).astype(bool)
            chosen = ledger.loc[ledger["chosen"], "condition"]
            selected = all_selected.loc[all_selected["condition"].isin(chosen)].copy()
            return selected, ledger.drop(columns="ranking_profiles")

        if mutant == "old_grouped_stool":
            patch.setattr(heatmaps, "prepare_contrast_cells", old_population)
        elif mutant == "profile_count_ranking":
            patch.setattr(heatmaps, "select_stool_group_contrasts", rank_profiles)
        reports = Reports()
        code = int(pytest.main(["-q", "--tb=short", str(ROOT / TEST)], plugins=[reports]))
    if before != _hashes():
        raise AssertionError(
            "Production, test or runner source changed during the mutation campaign"
        )
    record = {
        "mutant": mutant,
        "pytest_exit": code,
        "counts": dict(reports.counts),
        "failures": reports.failures,
        "source_sha256": before,
    }
    (output / "result.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return code


def worker_entrypoint(mutant: str, output: Path) -> None:
    with (
        (output / "pytest.log").open("w", encoding="utf-8") as log,
        redirect_stdout(log),
        redirect_stderr(log),
    ):
        raise SystemExit(worker(mutant, output))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expect-survivors", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    context = multiprocessing.get_context("spawn")
    records = []
    for mutant in MUTANTS:
        output = args.output / mutant
        output.mkdir()
        process = context.Process(target=worker_entrypoint, args=(mutant, output))
        process.start()
        process.join()
        record_path = output / "result.json"
        record = json.loads(record_path.read_text()) if record_path.exists() else {}
        expected = 0 if mutant == "baseline" or args.expect_survivors else 1
        valid = process.exitcode == expected and (expected == 0 or bool(record.get("failures")))
        records.append({**record, "expected_outcome": valid, "worker_exit": process.exitcode})
        print(f"{mutant}: exit={process.exitcode}; {record.get('counts', {})}", flush=True)
    (args.output / "summary.json").write_text(
        json.dumps(records, indent=2) + "\n", encoding="utf-8"
    )
    return int(not all(record["expected_outcome"] for record in records))


if __name__ == "__main__":
    sys.exit(main())
