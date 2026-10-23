"""Check fitted-label mutants without modifying canonical files (#970; 2026-10-23).

Inputs: this checkout's coverage tests, optionally read from a historical Git ref.
Outputs: pytest diagnostics and a JSON verdict under --output-dir. Only an isolated
test copy and in-memory plot wrapper are changed; no production data are used.
"""

from __future__ import annotations

import argparse
import io
import json
from collections.abc import Callable
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest
from git import Repo
from matplotlib.axes import Axes
from matplotlib.text import Text

from mtp.paths import REPO_ROOT, file_sha256

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "scripts/create_taxa_coverage_histogram.py"
TEST = ROOT / "tests/test_taxa_coverage_histogram.py"
Labels = list[tuple[Text, float, float, str, str]]


class InsideMutation:
    def __init__(self, mutation: str) -> None:
        self.mutation = mutation

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        """Alter the loaded figure after fitting; never rewrite the producer."""
        modules = {item.module for item in items if isinstance(item, pytest.Function)}
        for module in modules:
            coverage = module.coverage
            if Path(coverage.__file__).resolve() != SOURCE.resolve():
                raise AssertionError("Coverage producer is not from this checkout")
            original = coverage.fit_coverage_segment_labels
            mutation = self.mutation

            def altered(
                axis: Axes,
                labels: Labels,
                bar_height: float,
                original: Callable[[Axes, Labels, float], None] = original,
                mutation: str = mutation,
            ) -> None:
                original(axis, labels, bar_height)
                for text, _left, _right, _compact, description in labels:
                    if text.axes is axis and description.startswith("0 traits: 7.8%"):
                        if mutation == "drop_label":
                            text.remove()
                        elif mutation == "drop_count":
                            text.set_text(text.get_text().split("\n")[0])

            coverage.fit_coverage_segment_labels = altered


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mutation", choices=["baseline", "drop_label", "drop_count"])
    parser.add_argument("--tests-from-ref")
    parser.add_argument("--expect-survive", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if REPO_ROOT.resolve() != ROOT:
        raise AssertionError("Imported mtp belongs to a different checkout")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    before = {str(path.relative_to(ROOT)): file_sha256(path) for path in (SOURCE, TEST)}
    if args.tests_from_ref:
        with Repo(ROOT, search_parent_directories=False, expand_vars=False) as repository:
            blob = (
                repository.commit(args.tests_from_ref).tree
                / "tests/test_taxa_coverage_histogram.py"
            )
            text = blob.data_stream.read().decode("utf-8")
    else:
        text = TEST.read_text()
    assignment = (
        'SCRIPT = Path(__file__).resolve().parents[1] / "scripts/create_taxa_coverage_histogram.py"'
    )
    if text.count(assignment) != 1:
        raise AssertionError("Expected the test module's single producer-path assignment")
    copied = output / "test_coverage_copy.py"
    copied.write_text(text.replace(assignment, f"SCRIPT = Path({str(SOURCE)!r})"))
    captured = io.StringIO()
    with redirect_stdout(captured), redirect_stderr(captured):
        code = int(
            pytest.main(
                ["-q", "-p", "no:cacheprovider", str(copied), "-k", "narrow_segments"],
                plugins=[InsideMutation(args.mutation)],
            )
        )
    log = captured.getvalue()
    (output / "pytest.log").write_text(log)
    after = {str(path.relative_to(ROOT)): file_sha256(path) for path in (SOURCE, TEST)}
    expected_code = 0 if args.mutation == "baseline" or args.expect_survive else 1
    marker = "AssertionError: Missing fitted segment label:"
    verified = code == expected_code and before == after and (expected_code == 0 or marker in log)
    record = {
        "mutation": args.mutation,
        "tests_from_ref": args.tests_from_ref,
        "exit_code": code,
        "expected_code": expected_code,
        "expected_marker": marker,
        "canonical_before": before,
        "canonical_after": after,
        "verified": verified,
    }
    (output / "evidence.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))
    return int(not verified)


if __name__ == "__main__":
    raise SystemExit(main())
