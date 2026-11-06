"""Pinned complete-bundle byte parity, independent of native verification (#1169)."""

from __future__ import annotations

import argparse
import hashlib
import json
import runpy
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

SCRIPT = Path(__file__).parents[1] / "scripts/adhoc/issue1169_compare_derived_bytes.py"
HELPER = runpy.run_path(str(SCRIPT))
TABLES = (
    "composition_pair_delta_z.tsv",
    "composition_study_delta_z.tsv",
    "trait_summary.tsv",
    "quadrant_pairs.tsv",
    "quadrant_summary.tsv",
    "stool_disease_shifts.tsv",
    "stool_study_sensitivity.tsv",
    "site_disease_shifts.tsv",
)
COMMIT = "c" * 40


def identity(raw: bytes) -> dict[str, str | int]:
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


@dataclass
class Pair:
    baseline: Path
    current: Path
    records: dict[str, dict[str, Any]]
    pins: dict[str, str]

    def seal(self, role: str) -> None:
        path = getattr(self, role)
        raw = json.dumps(self.records[role], sort_keys=True).encode()
        path.write_bytes(raw)
        self.pins[role] = hashlib.sha256(raw).hexdigest()

    def compare(self, **kwargs: Any) -> dict[str, Any]:
        options = {
            "baseline_sha256": self.pins["baseline"],
            "current_sha256": self.pins["current"],
            "expected_current_commit": COMMIT,
            **kwargs,
        }
        result: dict[str, Any] = HELPER["compare_bundles"](self.baseline, self.current, **options)
        return result


@pytest.fixture
def pair(tmp_path: Path) -> Pair:
    parent_inputs = {"experiments.tsv@sha256:" + "e" * 64: {"sha256": "e" * 64, "size_bytes": 20}}
    record: dict[str, Any] = {
        "version": "delta-z-derived-v1",
        "summary_only": False,
        "run_id": "derived-baseline",
        "recorded_at": "2026-10-28T00:00:00+00:00",
        "code_commit": "b" * 40,
        "code_dirty": False,
        "parent_run_id": "parent",
        "parent_manifest_sha256": "a" * 64,
        "parent_producer": {"code_commit": "d" * 40, "code_dirty": False},
        "core_run_id": "core",
        "comparison_run_id": "comparison",
        "source_run_ids": {arm: arm for arm in ("core", "noel", "jed", "ilyas")},
        "source_commits": dict.fromkeys(("core", "noel", "jed", "ilyas"), "e" * 40),
        "strain_rollup": "phenotype",
        "scope_role": "primary",
        "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
        "direction_basis": "disease_vs_comparison",
        "analysis_bases": {"noel": "count"},
        "family_size": 1,
        "retained_traits": ["METPO:1"],
        "counts": {"traits": 1, "comparisons": 1, "studies": 1},
        "trait_list": {
            "trait_list_file": "headline.txt",
            "trait_list_sha256": "f" * 64,
            "trait_count": 1,
            "trait_ids": ["METPO:1"],
        },
        "policies": {"row3_top_diseases": 7, "row2_zero_policy": "exact_zero"},
        "parent_input_inventory": parent_inputs,
        "inputs": parent_inputs,
    }
    paths = {role: tmp_path / role / "delta_z_derived_run.json" for role in ("baseline", "current")}
    result = Pair(paths["baseline"], paths["current"], {}, {})
    for role, path in paths.items():
        path.parent.mkdir()
        row = deepcopy(record)
        if role == "current":
            row.update(code_commit=COMMIT, run_id="derived-current", recorded_at="later")
        row["outputs"] = {}
        for name in TABLES:
            raw = f"row\tvalue\n{name}\t1\n".encode()
            (path.parent / name).write_bytes(raw)
            row["outputs"][name] = identity(raw)
        result.records[role] = row
        result.seal(role)
    return result


def test_complete_byte_identical_bundles_pass_and_keep_provenance(pair: Pair) -> None:
    result = pair.compare()
    assert result["status"] == "PASS"
    assert set(result["tables"]) == set(TABLES)
    assert result["baseline"]["record"] == pair.records["baseline"]
    assert result["current"]["record"] == pair.records["current"]
    assert result["native_verification_performed"] is False
    for role in ("baseline", "current"):
        assert result[role]["sha256"] == pair.pins[role]


@pytest.mark.parametrize("scope", ["off", "phenotype"])
def test_both_f_scopes_are_supported(pair: Pair, scope: str) -> None:
    for role in pair.records:
        pair.records[role]["strain_rollup"] = scope
        pair.records[role]["scope_role"] = "sensitivity" if scope == "off" else "primary"
        pair.seal(role)
    assert pair.compare()["status"] == "PASS"


@pytest.mark.parametrize("role", ["baseline", "current"])
def test_independent_manifest_pin_is_required(pair: Pair, role: str) -> None:
    with pytest.raises(ValueError, match="manifest pin"):
        pair.compare(**{f"{role}_sha256": "0" * 64})


@pytest.mark.parametrize("mutation", ["changed_bound_byte", "unbound_byte", "missing", "extra"])
def test_table_identity_inventory_and_parity_fail_closed(pair: Pair, mutation: str) -> None:
    name = TABLES[-1]
    path = pair.current.parent / name
    if mutation in {"changed_bound_byte", "unbound_byte"}:
        raw = path.read_bytes().replace(b"\t1\n", b"\t2\n")
        path.write_bytes(raw)
        if mutation == "changed_bound_byte":
            pair.records["current"]["outputs"][name] = identity(raw)
    elif mutation == "missing":
        path.unlink()
        del pair.records["current"]["outputs"][name]
    else:
        path = pair.current.parent / "unexpected.tsv"
        path.write_bytes(b"extra\n")
        pair.records["current"]["outputs"][path.name] = identity(path.read_bytes())
    pair.seal("current")
    with pytest.raises(ValueError, match=r"bytes differ|inventory|identity"):
        pair.compare()


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("parent_manifest_sha256", "9" * 64),
        ("parent_run_id", "other-parent"),
        ("core_run_id", "other-core"),
        ("source_commits", {"core": "f" * 40}),
        ("strain_rollup", "off"),
        ("direction_basis", "as_reported"),
        ("trait_list", {"trait_ids": ["METPO:999"]}),
        ("policies", {"row3_top_diseases": 6, "row2_zero_policy": "exact_zero"}),
        ("counts", {"traits": 99}),
        ("inputs", {}),
    ],
)
def test_changed_contract_refused(pair: Pair, field: str, value: Any) -> None:
    pair.records["current"][field] = value
    pair.seal("current")
    with pytest.raises(ValueError, match=field):
        pair.compare()


def test_missing_contract_field_in_both_is_not_equality(pair: Pair) -> None:
    for role in pair.records:
        del pair.records[role]["policies"]
        pair.seal(role)
    with pytest.raises(ValueError, match="policies"):
        pair.compare()


@pytest.mark.parametrize("mutation", ["commit", "dirty", "summary"])
def test_current_producer_and_full_bundle_required(pair: Pair, mutation: str) -> None:
    key, value = {
        "commit": ("code_commit", "d" * 40),
        "dirty": ("code_dirty", True),
        "summary": ("summary_only", True),
    }[mutation]
    pair.records["current"][key] = value
    pair.seal("current")
    with pytest.raises(ValueError, match=r"producer|summary|complete"):
        pair.compare()


def test_manifest_parsing_uses_hashed_first_buffer(
    pair: Pair, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = Path.read_bytes
    first = pair.current.read_bytes()
    changed = json.dumps({**pair.records["current"], "code_commit": "f" * 40}).encode()
    reads = 0

    def read(path: Path) -> bytes:
        nonlocal reads
        if path == pair.current:
            reads += 1
            return changed if reads == 2 else first
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", read)
    record, _tables = HELPER["load_bundle"](pair.current, pair.pins["current"], {})
    assert record["code_commit"] == COMMIT


def test_cli_emits_json_stdout_without_creating_report_files(
    pair: Pair, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    args = argparse.Namespace(
        baseline_manifest=pair.baseline,
        current_manifest=pair.current,
        baseline_sha256=pair.pins["baseline"],
        current_sha256=pair.pins["current"],
        expected_current_commit=COMMIT,
    )
    monkeypatch.setattr(argparse.ArgumentParser, "parse_args", lambda _self: args)
    before = set(pair.current.parent.parent.rglob("*"))
    assert HELPER["main"]() == 0
    assert json.loads(capsys.readouterr().out)["status"] == "PASS"
    assert set(pair.current.parent.parent.rglob("*")) == before
