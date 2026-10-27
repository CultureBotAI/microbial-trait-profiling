"""Pinned-parent signed-rank sensitivity audit regressions (#1073, 2026-10-26)."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from mtp.three_test_table.wilcoxon_audit import AuditResult, audit_parent

LEGACY = "study-mean-delta-z-wilcoxon-v1"
OUTPUTS = (
    "composition_pair_delta_z.tsv",
    "composition_study_delta_z.tsv",
    "per_trait_evidence.tsv",
)


def _tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    path.write_text(stream.getvalue(), encoding="utf-8")


def _rows(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(path.read_text(encoding="utf-8")), delimiter="\t"))


@dataclass
class Parent:
    root: Path
    record: dict[str, Any]
    values: dict[str, list[float]]

    @property
    def manifest(self) -> Path:
        return self.root / "three_test_run.json"

    def seal(self) -> str:
        self.record["outputs"] = {
            name: {
                "sha256": hashlib.sha256(raw := (self.root / name).read_bytes()).hexdigest(),
                "size_bytes": len(raw),
            }
            for name in OUTPUTS
        }
        self.manifest.write_text(json.dumps(self.record, indent=2) + "\n", encoding="utf-8")
        return hashlib.sha256(self.manifest.read_bytes()).hexdigest()

    def audit(self) -> AuditResult:
        return audit_parent(self.manifest, self.seal(), scope="off")


def _make_parent(tmp_path: Path, values: dict[str, list[float]] | None = None) -> Parent:
    root = tmp_path / "pinned parent"
    root.mkdir()
    values = values or {
        "METPO:positive": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        "METPO:tied": [1.0, 1.0, -1.0, 2.0, -2.0, 0.0],
        "METPO:tiny": [1e-13, -1e-13, 0.0, 0.0, 0.0, 0.0],
        "METPO:zero": [0.0] * 6,
        "METPO:adjacent": [0.3, 0.30000000000000004, -0.1, -0.2, 0.0, 0.0],
    }
    studies = [
        {"study_id": f"s{i}", **{trait: vector[i] for trait, vector in values.items()}}
        for i in range(len(next(iter(values.values()))))
    ]
    pairs = [
        {
            "study_id": row["study_id"],
            "experiment_id": f"e{i}",
            **{trait: row[trait] for trait in values},
        }
        for i, row in enumerate(studies)
    ]
    pairs.append({**pairs[0], "experiment_id": "e-extra"})
    records: list[dict[str, Any]] = []
    probabilities: list[float] = []
    for trait, vector in values.items():
        nonzero = int(np.count_nonzero(vector))
        statistic, probability = "", ""
        if nonzero:
            result = stats.wilcoxon(vector, method="approx", zero_method="wilcox", correction=False)
            statistic, probability = str(float(result.statistic)), str(float(result.pvalue))
        probabilities.append(float(probability) if nonzero else 1.0)
        records.append(
            {
                "trait": trait,
                "trait_label": trait.removeprefix("METPO:"),
                "composition_contract": "full-catalog-study-delta-z-v1",
                "composition_family_size": len(values),
                "composition_pairs": len(pairs),
                "composition_studies": len(studies),
                "composition_test_unit": "study",
                "composition_effect_units": "standardized_log1p_1000_retained_fraction_difference",
                "composition_zero_total_policy": "retained_as_zero_fraction",
                "study_wilcoxon_contract": LEGACY,
                "study_wilcoxon_status": "tested" if nonzero else "untestable",
                "study_wilcoxon_reason": "" if nonzero else "all_tied",
                "study_wilcoxon_nonzero_studies": nonzero,
                "study_wilcoxon_zero_studies": len(studies) - nonzero,
                "study_wilcoxon_statistic": statistic,
                "study_wilcoxon_p": probability,
                "strain_rollup": "off",
                "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
                "direction_basis": "disease_vs_comparison",
                "core_run_id": "core-fixture",
            }
        )
    for correction in ("bh", "by"):
        for row, q in zip(
            records, stats.false_discovery_control(probabilities, method=correction), strict=True
        ):
            row[f"study_wilcoxon_q_{correction}"] = (
                float(q) if row["study_wilcoxon_status"] == "tested" else ""
            )
    _tsv(root / OUTPUTS[0], pairs)
    _tsv(root / OUTPUTS[1], studies)
    _tsv(root / OUTPUTS[2], records)
    return Parent(
        root,
        {
            "version": "three-test-per-trait-v1",
            "run_id": "three-test-fixture",
            "recorded_at": "2026-10-26T00:00:00+00:00",
            "code_commit": "a" * 40,
            "code_dirty": False,
            "code_inventory": {"declared-only.py": "b" * 64},
            "environment": {
                "python": "3.11.5",
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scipy": "1.17.1",
            },
            "strain_rollup": "off",
            "scope_role": "sensitivity",
            "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
            "direction_basis": "disease_vs_comparison",
            "core_run_id": "core-fixture",
            "comparison_run_id": "comparison-fixture",
            "comparison_manifest_sha256": "c" * 64,
            "source_run_ids": {
                "noel": "noel-fixture",
                "jed": "jed-fixture",
                "ilyas": "ilyas-fixture",
            },
            "source_commits": {"core": "a" * 40, "comparison": "a" * 40},
            "retained_traits": list(values),
            "family_size": len(values),
            "test_families": ["secondary-full-retained-study-wilcoxon"],
        },
        values,
    )


@pytest.fixture
def parent(tmp_path: Path) -> Parent:
    return _make_parent(tmp_path)


def test_reproduces_legacy_then_corrects_entire_family_without_source_writes(
    parent: Parent,
) -> None:
    digest = parent.seal()
    before = {path: path.read_bytes() for path in parent.root.iterdir()}
    result = audit_parent(parent.manifest, digest, scope="off")
    table = result.table.set_index("trait")
    assert list(table.index) == list(parent.values)
    assert table.loc["METPO:positive", "new_p"] == pytest.approx(0.03125)
    assert table.loc["METPO:tiny", "archived_status"] == "tested"
    assert table.loc["METPO:tiny", "new_status"] == "untestable"
    assert table.loc["METPO:tiny", "new_nonzero_studies"] == 0
    assert pd.isna(table.loc["METPO:tiny", "new_p"])
    assert pd.isna(table.loc["METPO:zero", "new_q_by"])
    new_p = [
        float(table.loc[trait, "new_p"]) if table.loc[trait, "new_status"] == "tested" else 1.0
        for trait in parent.values
    ]
    for correction in ("by", "bh"):
        expected = stats.false_discovery_control(new_p, method=correction)
        assert table.loc["METPO:positive", f"new_q_{correction}"] == expected[0]
    for metric in ("statistic", "p", "q_by", "q_bh"):
        for trait in ("METPO:positive", "METPO:tied", "METPO:tiny", "METPO:adjacent"):
            assert table.loc[trait, f"reproduced_{metric}"] == pytest.approx(
                table.loc[trait, f"archived_{metric}"], rel=1e-12, abs=1e-15
            )
    assert result.summary["family_size"] == 5
    assert result.summary["n_studies"] == 6
    assert result.summary["n_pairs"] == 7
    assert result.summary["status_changes"] == 1
    assert result.summary["changes"]["p"]["comparable_tested"] == 3
    assert not table.loc["METPO:tiny", "p_changed_exact"]
    assert pd.isna(table.loc["METPO:tiny", "p_absolute_delta"])
    assert not table.loc["METPO:zero", "p_changed_exact"]
    assert (
        result.summary["verification_scope"] == "pinned_parent_outputs_and_legacy_reproduction_only"
    )
    assert result.input_hashes == {
        path.resolve(): hashlib.sha256(raw).hexdigest() for path, raw in before.items()
    }
    assert {path: path.read_bytes() for path in parent.root.iterdir()} == before


def test_round_trip_preserves_adjacent_values_and_absolute_rank_distinctions(
    parent: Parent,
) -> None:
    result = parent.audit().table.set_index("trait")
    assert result.loc["METPO:adjacent", "reproduced_statistic"] == 3.0
    assert result.loc["METPO:adjacent", "new_method"] == "exact"
    expected = stats.wilcoxon([0.3, 0.30000000000000004, -0.1, -0.2], method="exact")
    assert result.loc["METPO:adjacent", "new_p"] == expected.pvalue


def test_delta_counts_compare_new_to_reproduced_not_archived_rounding(tmp_path: Path) -> None:
    parent = _make_parent(
        tmp_path,
        {
            "METPO:positive": [float(i) for i in range(1, 52)],
            "METPO:zero": [0.0] * 51,
        },
    )
    rows = _rows(parent.root / OUTPUTS[2])
    rows[0]["study_wilcoxon_p"] = repr(math.nextafter(float(rows[0]["study_wilcoxon_p"]), 1.0))
    _tsv(parent.root / OUTPUTS[2], rows)
    result = parent.audit()
    table = result.table.set_index("trait")
    assert table.loc["METPO:positive", "archived_p"] != table.loc["METPO:positive", "reproduced_p"]
    assert table.loc["METPO:positive", "new_p"] == table.loc["METPO:positive", "reproduced_p"]
    assert not table.loc["METPO:positive", "p_changed_exact"]
    for metric in ("p", "q_by", "q_bh"):
        assert result.summary["changes"][metric]["exact"] == 0
        assert result.summary["changes"][metric]["material"] == 0
        assert result.summary["changes"][metric]["comparable_tested"] == 1


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("version", "unknown"),
        ("strain_rollup", "phenotype"),
        ("scope_role", "primary"),
        ("family_size", 4),
        ("family_size", True),
        ("retained_traits", ["duplicate", "duplicate"]),
    ],
)
def test_refuses_malformed_declared_parent_contract(parent: Parent, field: str, value: Any) -> None:
    parent.record[field] = value
    with pytest.raises(ValueError, match=r"Parent.*(?:invalid|catalog)"):
        parent.audit()


def test_refuses_wrong_pin_before_any_output_read(parent: Parent) -> None:
    parent.seal()
    for name in OUTPUTS:
        (parent.root / name).unlink()
    with pytest.raises(ValueError, match=r"pin|SHA256"):
        audit_parent(parent.manifest, "f" * 64, scope="off")


def test_refuses_duplicate_manifest_keys(parent: Parent) -> None:
    parent.seal()
    raw = parent.manifest.read_text(encoding="utf-8").replace(
        '"family_size": 5', '"family_size": 4, "family_size": 5'
    )
    parent.manifest.write_text(raw, encoding="utf-8")
    digest = hashlib.sha256(parent.manifest.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="Duplicate manifest key"):
        audit_parent(parent.manifest, digest, scope="off")


@pytest.mark.parametrize("field", ["study_wilcoxon_contract", "study_wilcoxon_p"])
def test_refuses_missing_legacy_contract_or_inference_columns(parent: Parent, field: str) -> None:
    rows = _rows(parent.root / OUTPUTS[2])
    for row in rows:
        del row[field]
    _tsv(parent.root / OUTPUTS[2], rows)
    with pytest.raises(ValueError, match=r"contract differs|Missing legacy"):
        parent.audit()


def test_parent_manifest_alias_uses_actual_bundle_and_records_both_identities(
    parent: Parent,
) -> None:
    digest = parent.seal()
    alias = parent.root.parent / "parent-manifest-alias.json"
    alias.symlink_to(parent.manifest)
    result = audit_parent(alias, digest, scope="off")
    assert result.input_hashes[alias] == digest
    assert result.input_hashes[parent.manifest.resolve()] == digest
    assert len(result.input_hashes) == 5


@pytest.mark.parametrize("damage", ["hash", "size", "missing", "extra", "symlink"])
def test_refuses_output_inventory_damage(parent: Parent, damage: str) -> None:
    digest = parent.seal()
    path = parent.root / OUTPUTS[1]
    if damage == "hash":
        path.write_bytes(path.read_bytes().replace(b"1.0", b"9.0", 1))
    elif damage == "size":
        parent.record["outputs"][OUTPUTS[1]]["size_bytes"] += 1
    elif damage == "missing":
        path.unlink()
    elif damage == "extra":
        (parent.root / "unexpected.tsv").write_text("x", encoding="utf-8")
    else:
        target = parent.root.parent / "aliased-study.tsv"
        path.rename(target)
        path.symlink_to(target)
    if damage == "size":
        parent.manifest.write_text(json.dumps(parent.record), encoding="utf-8")
        digest = hashlib.sha256(parent.manifest.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match=r"Parent output"):
        audit_parent(parent.manifest, digest, scope="off")


@pytest.mark.parametrize(
    "damage",
    [
        "study_duplicate",
        "study_blank",
        "study_missing",
        "study_nonfinite",
        "study_reordered",
        "pair_duplicate",
        "pair_orphan",
        "trait_duplicate",
        "trait_reordered",
    ],
)
def test_refuses_ledger_identity_and_numerical_damage(parent: Parent, damage: str) -> None:
    name = (
        OUTPUTS[0]
        if damage.startswith("pair")
        else OUTPUTS[2]
        if damage.startswith("trait")
        else OUTPUTS[1]
    )
    rows = _rows(parent.root / name)
    if damage.endswith("duplicate"):
        rows.append(rows[0].copy())
    elif damage == "study_blank":
        rows[0]["study_id"] = " "
    elif damage == "study_missing":
        rows.pop()
    elif damage == "study_nonfinite":
        rows[0]["METPO:positive"] = "nan"
    elif damage == "study_reordered":
        rows = [{key: row[key] for key in reversed(row)} for row in rows]
    elif damage == "pair_orphan":
        rows[0]["study_id"] = "unknown-study"
    else:
        rows.reverse()
    _tsv(parent.root / name, rows)
    with pytest.raises(ValueError, match=r"Parent (?:TSV|ledger|numerical)"):
        parent.audit()


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("study_wilcoxon_contract", "study-mean-delta-z-wilcoxon-v2"),
        ("study_wilcoxon_status", "untestable"),
        ("study_wilcoxon_reason", "all_tied"),
        ("study_wilcoxon_nonzero_studies", "5"),
        ("study_wilcoxon_zero_studies", "1"),
        ("study_wilcoxon_statistic", "1.0"),
        ("study_wilcoxon_p", "0.4"),
        ("study_wilcoxon_q_by", "0.4"),
        ("study_wilcoxon_q_bh", "0.4"),
        ("composition_studies", "5"),
        ("composition_pairs", "6"),
        ("study_wilcoxon_p", "nan"),
        ("study_wilcoxon_p", ""),
        ("strain_rollup", "phenotype"),
        ("core_run_id", "other-core"),
    ],
)
def test_refuses_every_legacy_inference_or_support_mismatch(
    parent: Parent, field: str, value: str
) -> None:
    rows = _rows(parent.root / OUTPUTS[2])
    rows[0][field] = value
    _tsv(parent.root / OUTPUTS[2], rows)
    with pytest.raises(ValueError, match=r"Parent|Legacy"):
        parent.audit()


def test_refuses_inferential_values_on_untestable_archived_row(parent: Parent) -> None:
    rows = _rows(parent.root / OUTPUTS[2])
    rows[3]["study_wilcoxon_p"] = "1.0"
    _tsv(parent.root / OUTPUTS[2], rows)
    with pytest.raises(ValueError, match="Untestable legacy inference"):
        parent.audit()
