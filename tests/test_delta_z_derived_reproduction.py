"""Inspect completed F outputs against declared references; no Delta-z replay (#1052).

Inputs: off summary, phenotype derived bundle, historical numerical JSON, producer,
core experiments table and an explicit historical JSON SHA-256 pin.
Output: pytest assertions/log only. Date: 2026-10-25.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Any

import pytest
from git import Blob, Repo

from delta_z_reference_helpers import F_COMMIT
from mtp.paths import file_sha256

pytestmark = pytest.mark.integration
Inputs = tuple[Path, Path, Path, str, Path, str]
HEADLINE_TRAIT_LIST = "methods/trait_lists/headline_ten_phenotype.txt"


def table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.fixture(scope="module")
def inputs() -> Inputs:
    names = (
        "MTP_DZ_DERIVED_OFF",
        "MTP_DZ_DERIVED_PHENOTYPE",
        "MTP_DZ_HISTORICAL_STATS",
        "MTP_DZ_DERIVED_COMMIT",
        "MTP_DZ_CORE_EXPERIMENTS",
        "MTP_DZ_HISTORICAL_SHA256",
    )
    if not any(name in os.environ for name in names):
        pytest.skip("Requires completed F derived bundles and historical numerical references")
    missing = [name for name in names if not os.environ.get(name, "").strip()]
    if missing:
        raise ValueError(f"Missing reproduction inputs: {missing}")
    return (
        Path(os.environ[names[0]]),
        Path(os.environ[names[1]]),
        Path(os.environ[names[2]]),
        os.environ[names[3]],
        Path(os.environ[names[4]]),
        os.environ[names[5]],
    )


def record(root: Path, scope: str, producer: str) -> dict[str, Any]:
    result: dict[str, Any] = json.loads((root / "delta_z_derived_run.json").read_text())
    assert result["code_commit"] == producer
    assert result["code_dirty"] is False
    assert result["strain_rollup"] == scope
    assert result["contrast_population"] == "disease_vs_healthy_tier1_2_subject_controls"
    assert result["direction_basis"] == "disease_vs_comparison"
    assert set(result["source_commits"].values()) == {F_COMMIT}
    assert {path.name for path in root.iterdir()} == {
        "delta_z_derived_run.json",
        *result["outputs"],
    }
    for name, identity in result["outputs"].items():
        assert Path(name).name == name
        assert not (root / name).is_symlink()
        assert file_sha256(root / name) == identity["sha256"]
        assert (root / name).stat().st_size == identity["size_bytes"]
    return result


def test_active_study_references(inputs: Inputs) -> None:
    off, phenotype, _, producer, _, _ = inputs
    for root, scope, traits, active in (
        # Historical F sources and persisted ledgers give 291; the brief's 289 was stale (#1101).
        (off, "off", 317, 291),
        (phenotype, "phenotype", 344, 316),
    ):
        result = record(root, scope, producer)
        assert result["counts"] == {
            "traits": traits,
            "comparisons": 935,
            "studies": 548,
            "traits_with_at_least_50_active_studies": active,
        }
        rows = table(root / "trait_summary.tsv")
        assert len(rows) == traits
        assert len({row["trait"] for row in rows}) == traits
        flags = [row["sdz_active_at_least_50"] for row in rows]
        assert set(flags).issubset({"True", "False"})
        assert (
            sum(flag == "True" for flag in flags)
            == result["counts"]["traits_with_at_least_50_active_studies"]
        )


def test_oxygen_quadrant_references(inputs: Inputs) -> None:
    _, phenotype, _, _, _, _ = inputs
    rows = {row["population"]: row for row in table(phenotype / "quadrant_summary.tsv")}
    assert set(rows) == {"all", "stool"}
    for population, count, q4, q2 in (("all", 935, 328, 253), ("stool", 595, 256, 140)):
        row = rows[population]
        assert row["fa_trait"] == "METPO:1000608"
        assert row["oa_trait"] == "METPO:1000607"
        assert int(row["n_comparisons"]) == count
        assert int(row["fa_up_oa_down"]) == q4
        assert int(row["fa_down_oa_up"]) == q2
        assert (
            sum(
                int(row[name])
                for name in (
                    "fa_up_oa_down",
                    "fa_down_oa_up",
                    "both_up",
                    "both_down",
                    "axis",
                )
            )
            == count
        )
    # Published reference precisions, declared before qualification.
    assert round(float(rows["all"]["sign_test_p"]), 4) == 0.0021
    assert round(float(rows["all"]["pearson_r"]), 2) == -0.45
    assert format(float(rows["stool"]["sign_test_p"]), ".2g") == "5.9e-09"


def committed_headline_traits(result: dict[str, Any], producer: str) -> list[str]:
    """Require the exact committed list and its full provenance, independently of outputs."""
    with Repo(Path(__file__).resolve().parents[1]) as repository:
        blob = repository.commit(producer).tree / HEADLINE_TRAIT_LIST
        assert isinstance(blob, Blob)
        raw = blob.data_stream.read()
    declared = [
        line
        for source in raw.decode("utf-8").splitlines()
        if (line := source.strip()) and not line.startswith("#")
    ]
    assert len(declared) == len(set(declared))
    assert len(declared) in {9, 10}
    assert result["trait_list"] == {
        "trait_list_file": Path(HEADLINE_TRAIT_LIST).name,
        "trait_list_sha256": hashlib.sha256(raw).hexdigest(),
        "trait_count": len(declared),
        "trait_ids": declared,
    }
    return declared


def pinned_historical_reference(path: Path, expected_sha256: str) -> dict[str, Any]:
    """Hash exactly the historical bytes being parsed, using the caller's explicit pin."""
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == expected_sha256
    reference: dict[str, Any] = json.loads(raw)
    assert isinstance(reference, dict)
    return reference


def core_pair_metadata(
    path: Path, result: dict[str, Any]
) -> dict[tuple[str, str], tuple[str, str]]:
    """Independently parse the exact source metadata bytes named by both inventories."""
    bindings = [
        (identity, value)
        for identity, value in result["inputs"].items()
        if identity.partition("@sha256:")[0] == "experiments.tsv"
    ]
    assert len(bindings) == 1
    identity, expected = bindings[0]
    assert identity == f"experiments.tsv@sha256:{expected['sha256']}"
    assert result["parent_input_inventory"][identity] == expected
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == expected["sha256"]
    assert len(raw) == expected["size_bytes"]
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8"), newline=""), delimiter="\t")
    required = ["study_id", "experiment_id", "body_site", "condition"]
    assert reader.fieldnames is not None
    assert len(reader.fieldnames) == len(set(reader.fieldnames))
    assert set(required).issubset(reader.fieldnames)
    metadata = {}
    for row in reader:
        assert all(isinstance(row[name], str) for name in required)
        key = row["study_id"], row["experiment_id"]
        assert all(value.strip() for value in key)
        assert key not in metadata
        metadata[key] = row["body_site"], row["condition"]
    return metadata


def test_seventy_stool_cells_against_historical_numbers(inputs: Inputs) -> None:
    _, phenotype, historical, producer, _, historical_sha256 = inputs
    result = record(phenotype, "phenotype", producer)
    declared = committed_headline_traits(result, producer)
    historical_rows = pinned_historical_reference(historical, historical_sha256)["s14_tier12"]
    # Historical display aliases only (#1102); preserve exact core condition identities.
    aliases = {
        "Irritable bowel": "Irritable bowel syndrome",
        "Polycystic ovary": "Polycystic ovary syndrome",
    }
    reference = {aliases.get(row["name"], row["name"]): row for row in historical_rows}
    assert len(reference) == len(historical_rows)
    rows = table(phenotype / "stool_disease_shifts.tsv")
    assert len(rows) == 7 * len(declared)
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row["condition"]].append(row)
    assert len(groups) == 7
    assert set(groups) == set(reference), "Historical labels differ from exact core conditions"
    assert sorted((int(group[0]["n_comparisons"]) for group in groups.values()), reverse=True) == [
        70,
        45,
        31,
        31,
        21,
        18,
        17,
    ]
    for condition, group in groups.items():
        expected = reference[condition]
        assert [row["trait"] for row in group] == declared
        for row in group:
            original = expected["traits"][row["trait"]]
            assert row["body_site"] == "Feces"
            assert int(row["n_comparisons"]) == expected["n_comp"]
            assert int(row["n_studies"]) == expected["n_stud"]
            assert row["p_adjustment"] == "none_nominal"
            assert float(row["mean_delta_z"]) == pytest.approx(
                original["mean"],
                rel=1e-12,
                abs=1e-14,
            )
            assert float(row["wilcoxon_p"]) == pytest.approx(
                original["p_wilcoxon"],
                rel=1e-12,
                abs=1e-14,
            )


def test_site_condition_cells_against_persisted_pairs(inputs: Inputs) -> None:
    _, phenotype, _, producer, experiments, _ = inputs
    result = record(phenotype, "phenotype", producer)
    pair_rows = table(phenotype / "composition_pair_delta_z.tsv")
    pairs = {(row["study_id"], row["experiment_id"]): row for row in pair_rows}
    assert len(pair_rows) == len(pairs) == 935
    traits = result["retained_traits"]
    assert len(traits) == len(set(traits)) == result["family_size"]
    assert all(set(row) == {"study_id", "experiment_id", *traits} for row in pair_rows)
    source_metadata = core_pair_metadata(experiments, result)
    assert set(pairs).issubset(source_metadata)
    reported_rows = table(phenotype / "quadrant_pairs.tsv")
    reported = {(row["study_id"], row["experiment_id"]): row for row in reported_rows}
    assert len(reported_rows) == len(reported) == len(pairs)
    assert set(reported) == set(pairs)
    for key, row in reported.items():
        assert (row["body_site"], row["condition"]) == source_metadata[key]
    groups: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    for key in pairs:
        groups[source_metadata[key]].append(key)
    cells = table(phenotype / "site_disease_shifts.tsv")
    assert len(cells) == len(groups) * result["family_size"]
    identities = {(row["body_site"], row["condition"], row["trait"]) for row in cells}
    assert len(identities) == len(cells)
    assert identities == {
        (site, condition, trait) for site, condition in groups for trait in traits
    }
    for row in cells:
        keys = groups[(row["body_site"], row["condition"])]
        assert int(row["n_comparisons"]) == len(keys)
        assert int(row["n_studies"]) == len({key[0] for key in keys})
        # Independent scalar aggregation of the persisted Delta-z ledger only.
        mean = math.fsum(float(pairs[key][row["trait"]]) for key in keys) / len(keys)
        assert float(row["mean_delta_z"]) == pytest.approx(mean, rel=1e-12, abs=1e-14)
