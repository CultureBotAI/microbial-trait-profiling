"""Real-producer headroom and hash-seed-stable charges (#1909, #1910).

The loop variants reproduce the reviewer's AST-only edits; no composition
analysis is imported or executed. As a script, this file measures the same
sources in a fresh interpreter for the cross-process regression test.
"""

from __future__ import annotations

import ast
import hashlib
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from textwrap import indent
from typing import Any

import pytest

import test_multiple_testing_registry as guard
import test_multiple_testing_registry_work as headroom

ROOT = Path(__file__).resolve().parents[1]
RELATIVE_SOURCE = "src/mtp/composition_calibration/composition.py"
FUNCTION = "run_composition_calibration"
VARIANTS = ("unedited", "three_line_loop", "four_line_loop")
PREFIX_CASES = ("nested_prefixes", "loop_nested_prefixes")


def _prefix_sources() -> dict[str, str]:
    """Unequal-depth aliases force ordered intermediate cache construction."""
    header = (
        "import types\nimport numpy as np\n"
        "from mtp.enrichment_stats import benjamini_hochberg\n"
        "def run_composition_calibration(p, other):\n"
        "    holder = types.SimpleNamespace(a=p, b=[p], c=[[p]])\n"
    )
    tail = (
        "    safe = [0.01, 0.02, 0.03]\n"
        "    m = len(safe)\n"
        "    h = np.sum(1.0 / np.arange(1, m + 1))\n"
        "    q = benjamini_hochberg(safe)\n"
        "    return np.clip(q * h, 0, 1)\n"
    )
    return {
        "nested_prefixes": header + "    p[0] = other\n" + tail,
        "loop_nested_prefixes": header + "    for _ in (0, 1):\n        p[0] = other\n" + tail,
    }


def _assert_assigned_sources() -> None:
    for module in (guard, headroom):
        assert module.__file__ is not None
        assert Path(module.__file__).resolve() == ROOT / "tests" / f"{module.__name__}.py"
    assert guard.ROOT.resolve() == ROOT


@pytest.fixture(autouse=True)
def assigned_guard() -> None:
    _assert_assigned_sources()


def _variants(source: str) -> dict[str, str]:
    tree = ast.parse(source, filename=RELATIVE_SOURCE)
    function = next(
        node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == FUNCTION
    )
    first = function.body[0]
    if (
        isinstance(first, ast.Expr)
        and isinstance(first.value, ast.Constant)
        and isinstance(first.value.value, str)
    ):
        first = function.body[1]
    loop = next(node for node in function.body if isinstance(node, ast.For))
    assert isinstance(loop.target, ast.Name)
    assert loop.target.id == "thresh"
    assert isinstance(loop.iter, ast.Name)
    assert loop.iter.id == "THRESHOLDS"
    segment = ast.get_source_segment(source, loop)
    assert segment is not None
    loop_lines = segment.splitlines()
    assert len(loop_lines) == 3, "Keep the reviewer's three-line loop witness explicit"
    loop_lines = [loop_lines[0]] + [line[loop.col_offset :] for line in loop_lines[1:]]
    blocks = {
        "three_line_loop": "\n".join(loop_lines) + "\n",
        "four_line_loop": (
            "for _key in study_keys:\n"
            "    _row = study_contrasts.loc[_key]\n"
            "    if _row is None:\n"
            "        break\n"
        ),
    }
    lines = source.splitlines(keepends=True)
    insertion = first.lineno - 1
    return {
        "unedited": source,
        **{
            name: (
                "".join(lines[:insertion])
                + indent(block, " " * first.col_offset)
                + "".join(lines[insertion:])
            )
            for name, block in blocks.items()
        },
    }


@dataclass(frozen=True)
class Measurement:
    credited: bool
    categories: dict[str, int]

    @property
    def total(self) -> int:
        return sum(self.categories.values())


def _measure(text: str) -> Measurement:
    module = RELATIVE_SOURCE.removeprefix("src/").removesuffix(".py").replace("/", ".")
    unit = guard.Unit(
        RELATIVE_SOURCE,
        module,
        ast.parse(text, filename=RELATIVE_SOURCE),
        guard.DEFAULT_LIMITS,
    )
    guard.Visitor(unit).visit(unit.tree)
    categories: dict[str, int] = {}
    credited = guard.visible_by_derivation(unit, FUNCTION, work_counts=categories)
    assert categories, "Headroom cannot be established from absent accounting"
    assert all(value >= 0 for value in categories.values())
    assert sum(categories.values()) > 0
    return Measurement(credited, categories)


def test_main_headroom_allowances_remain_unchanged() -> None:
    assert headroom.DERIVATION_HEADROOM_DIVISOR == 16
    assert guard.DEFAULT_LIMITS.derivation_steps // headroom.DERIVATION_HEADROOM_DIVISOR == 62_500
    assert (
        guard.DEFAULT_LIMITS.file_derivation_steps // headroom.DERIVATION_HEADROOM_DIVISOR == 62_500
    )


@pytest.mark.parametrize("variant", VARIANTS)
def test_real_composition_and_ordinary_loop_edits_keep_by_within_the_headroom_allowance(
    variant: str,
) -> None:
    source = (ROOT / RELATIVE_SOURCE).read_text(encoding="utf-8")
    measured = _measure(_variants(source)[variant])
    assert measured.credited, f"Real producer lost its BY witness: {variant}"
    # Exercise ordinary growth against the owner's existing allowance. A new,
    # stricter ceiling would recreate the routine-edit failure being fixed.
    limit = guard.DEFAULT_LIMITS.derivation_steps // headroom.DERIVATION_HEADROOM_DIVISOR
    assert measured.total <= limit, (variant, measured.total, limit, measured.categories)


def _measurement_report() -> dict[str, Any]:
    _assert_assigned_sources()
    assert guard.__file__ is not None
    guard_path = Path(guard.__file__).resolve()
    source_bytes = (ROOT / RELATIVE_SOURCE).read_bytes()
    measurements = {}
    cases = _variants(source_bytes.decode("utf-8")) | _prefix_sources()
    for name, source in cases.items():
        measured = _measure(source)
        measurements[name] = {
            "credited": measured.credited,
            "categories": measured.categories,
            "total": measured.total,
        }
    return {
        "seed": os.environ.get("PYTHONHASHSEED"),
        "guard_path": str(guard_path),
        "guard_sha256": hashlib.sha256(guard_path.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "measurements": measurements,
    }


def test_real_producer_category_charges_match_in_two_fresh_hash_seed_processes() -> None:
    assert guard.__file__ is not None
    guard_path = Path(guard.__file__).resolve()
    expected_guard_sha = hashlib.sha256(guard_path.read_bytes()).hexdigest()
    expected_source_sha = hashlib.sha256((ROOT / RELATIVE_SOURCE).read_bytes()).hexdigest()
    reports = []
    for seed in ("0", "1"):
        environment = os.environ.copy()
        environment["PYTHONHASHSEED"] = seed
        environment["PYTHONDONTWRITEBYTECODE"] = "1"
        environment["PYTHONPATH"] = os.pathsep.join(str(ROOT / part) for part in ("src", "tests"))
        completed = subprocess.run(
            [sys.executable, "-B", str(Path(__file__).resolve()), "--measure"],
            cwd=ROOT,
            env=environment,
            check=True,
            shell=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
        report = json.loads(completed.stdout)
        assert report["seed"] == seed
        assert report["guard_path"] == str(guard_path)
        assert report["guard_sha256"] == expected_guard_sha
        assert report["source_sha256"] == expected_source_sha
        assert set(report["measurements"]) == set(VARIANTS) | set(PREFIX_CASES)
        for measured in report["measurements"].values():
            assert measured["credited"] is True
            assert measured["categories"]
            assert measured["total"] == sum(measured["categories"].values())
        reports.append(report)
    for variant in (*VARIANTS, *PREFIX_CASES):
        first = reports[0]["measurements"][variant]
        second = reports[1]["measurements"][variant]
        assert first["categories"] == second["categories"], (variant, first, second)
        assert first["total"] == second["total"]


if __name__ == "__main__":
    if sys.argv[1:] != ["--measure"]:
        raise SystemExit("Use --measure to collect the hash-seed regression measurements")
    print(json.dumps(_measurement_report(), sort_keys=True))
