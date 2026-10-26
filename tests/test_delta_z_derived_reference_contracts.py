"""Independent qualification-input regressions (#1096, 2026-10-25).

Synthetic completed bundles exercise the actual F qualification assertions. No
upstream pipeline, real-data analysis or Delta-z transformation is executed.
"""

from __future__ import annotations

import csv
import hashlib
import json
import tempfile
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest
from git import Repo

import test_delta_z_derived_reproduction as qualification
from mtp.paths import file_sha256

REPOSITORY = Path(__file__).resolve().parents[1]
TRAIT_LIST = "methods/trait_lists/headline_ten_phenotype.txt"
EXTRA_TRAIT = "METPO:2000103|pathways:trithionate_oxidation"
F_COMMIT = "511c891229b659e207cfe0f4281a97c3fc039a80"
HISTORICAL_DISPLAY_LABELS = {
    "Irritable bowel syndrome": "Irritable bowel",
    "Polycystic ovary syndrome": "Polycystic ovary",
}


@pytest.fixture
def tmp_path() -> Iterator[Path]:
    root = REPOSITORY / ".tmp" / "delta-z-reference-tests"
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="case-", dir=root) as directory:
        yield Path(directory)


def write_table(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


@dataclass
class QualificationFixture:
    root: Path
    historical: Path
    experiments: Path
    producer: str
    historical_sha256: str
    record: dict[str, Any]

    def save(self) -> None:
        self.record["outputs"] = {
            path.name: {"sha256": file_sha256(path), "size_bytes": path.stat().st_size}
            for path in sorted(self.root.glob("*.tsv"))
        }
        (self.root / "delta_z_derived_run.json").write_text(
            json.dumps(self.record), encoding="utf-8"
        )

    def inputs(self) -> qualification.Inputs:
        return (
            self.root,
            self.root,
            self.historical,
            self.producer,
            self.experiments,
            self.historical_sha256,
        )


def base_fixture(tmp_path: Path) -> QualificationFixture:
    root = tmp_path / "phenotype"
    root.mkdir()
    historical = tmp_path / "historical.json"
    historical.write_text("{}\n", encoding="utf-8")
    experiments = tmp_path / "core" / "experiments.tsv"
    experiments.parent.mkdir()
    with Repo(REPOSITORY) as repository:
        producer = repository.head.commit.hexsha
        blob = repository.commit(producer).tree / TRAIT_LIST
        raw_list = blob.data_stream.read()
    traits = [
        line.strip()
        for line in raw_list.decode().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    record: dict[str, Any] = {
        "code_commit": producer,
        "code_dirty": False,
        "strain_rollup": "phenotype",
        "contrast_population": "disease_vs_healthy_tier1_2_subject_controls",
        "direction_basis": "disease_vs_comparison",
        "source_commits": dict.fromkeys(["comparison", "core", "noel", "jed", "ilyas"], F_COMMIT),
        "trait_list": {
            "trait_list_file": Path(TRAIT_LIST).name,
            "trait_list_sha256": hashlib.sha256(raw_list).hexdigest(),
            "trait_count": len(traits),
            "trait_ids": traits,
        },
        "family_size": len(traits),
        "retained_traits": traits,
        "inputs": {},
        "parent_input_inventory": {},
    }
    return QualificationFixture(
        root, historical, experiments, producer, file_sha256(historical), record
    )


def stool_fixture(tmp_path: Path, change: str = "") -> QualificationFixture:
    fixture = base_fixture(tmp_path)
    declared = fixture.record["trait_list"]
    original_traits = list(declared["trait_ids"])
    if change == "extra_trithionate":
        declared["trait_ids"] = [*original_traits[:-1], EXTRA_TRAIT]
    elif change == "reordered":
        declared["trait_ids"] = original_traits[::-1]
    elif change == "filename":
        declared["trait_list_file"] = "unrelated.txt"
    elif change == "hash":
        declared["trait_list_sha256"] = "f" * 64
    elif change == "count":
        declared["trait_count"] = 11
    reference = []
    rows = []
    for number, count in enumerate([70, 45, 31, 31, 21, 18, 17]):
        condition = f"Condition {number}"
        values = {
            trait: {"mean": number + (position + 1) / 10, "p_wilcoxon": (position + 1) / 100}
            for position, trait in enumerate([*original_traits, EXTRA_TRAIT])
        }
        reference.append(
            {"name": condition, "n_comp": count, "n_stud": count // 2, "traits": values}
        )
        for trait in declared["trait_ids"]:
            rows.append(
                {
                    "condition": condition,
                    "trait": trait,
                    "body_site": "Feces",
                    "n_comparisons": count,
                    "n_studies": count // 2,
                    "p_adjustment": "none_nominal",
                    "mean_delta_z": values[trait]["mean"],
                    "wilcoxon_p": values[trait]["p_wilcoxon"],
                }
            )
    fixture.historical.write_text(json.dumps({"s14_tier12": reference}), encoding="utf-8")
    fixture.historical_sha256 = file_sha256(fixture.historical)
    write_table(fixture.root / "stool_disease_shifts.tsv", rows)
    fixture.save()
    return fixture


def named_stool_fixture(tmp_path: Path, *, display_labels: bool) -> QualificationFixture:
    """Use exact consumer condition names and optional historical aliases (#1102)."""
    fixture = stool_fixture(tmp_path)
    names = dict(zip(["Condition 0", "Condition 1"], HISTORICAL_DISPLAY_LABELS, strict=True))
    rows = qualification.table(fixture.root / "stool_disease_shifts.tsv")
    for row in rows:
        row["condition"] = names.get(row["condition"], row["condition"])
    write_table(fixture.root / "stool_disease_shifts.tsv", rows)
    reference = json.loads(fixture.historical.read_text(encoding="utf-8"))
    for original in reference["s14_tier12"]:
        canonical = names.get(original["name"], original["name"])
        original["name"] = (
            HISTORICAL_DISPLAY_LABELS.get(canonical, canonical) if display_labels else canonical
        )
    save_reference(fixture, reference)
    return fixture


def save_reference(fixture: QualificationFixture, reference: dict[str, Any]) -> None:
    """Pin each synthetic reference after deliberate fixture construction, before validation."""
    fixture.historical.write_text(json.dumps(reference), encoding="utf-8")
    fixture.historical_sha256 = file_sha256(fixture.historical)
    fixture.save()


def site_fixture(tmp_path: Path, change: str = "") -> QualificationFixture:
    fixture = base_fixture(tmp_path)
    traits = ["TEST:a", "TEST:b"]
    fixture.record["retained_traits"] = traits
    fixture.record["family_size"] = len(traits)
    pairs: list[dict[str, Any]] = []
    core: list[dict[str, Any]] = []
    consumer: list[dict[str, Any]] = []
    for number in range(935):
        identity = {"study_id": f"s{number // 3}", "experiment_id": f"e{number}"}
        pairs.append({**identity, "TEST:a": float(number % 7 - 3), "TEST:b": float(number % 5 - 2)})
        original = {
            **identity,
            "body_site": "Feces" if number % 2 else "Saliva",
            "body_site_category": "Gut" if number % 2 else "Oral",
            "condition": "Condition A" if number % 3 else "Condition B",
        }
        core.append(original)
        reported = {key: value for key, value in original.items() if key != "body_site_category"}
        if change == "category_field":
            reported["body_site"] = original["body_site_category"]
        elif change == "site_relabel":
            reported["body_site"] = "Stool" if original["body_site"] == "Feces" else "Saliva"
        elif change == "condition_relabel":
            reported["condition"] = "Renamed " + original["condition"]
        consumer.append(reported)
    write_table(fixture.experiments, core)
    input_identity = {
        "sha256": file_sha256(fixture.experiments),
        "size_bytes": fixture.experiments.stat().st_size,
    }
    input_name = f"experiments.tsv@sha256:{input_identity['sha256']}"
    fixture.record["inputs"][input_name] = input_identity
    fixture.record["parent_input_inventory"][input_name] = input_identity
    groups: dict[tuple[str, str], list[int]] = defaultdict(list)
    for number, row in enumerate(consumer):
        groups[(row["body_site"], row["condition"])].append(number)
    cells = []
    for (site, condition), members in groups.items():
        for trait in traits:
            cells.append(
                {
                    "body_site": site,
                    "condition": condition,
                    "trait": trait,
                    "n_comparisons": len(members),
                    "n_studies": len({pairs[index]["study_id"] for index in members}),
                    "mean_delta_z": sum(pairs[index][trait] for index in members) / len(members),
                }
            )
    write_table(fixture.root / "composition_pair_delta_z.tsv", pairs)
    write_table(fixture.root / "quadrant_pairs.tsv", consumer)
    write_table(fixture.root / "site_disease_shifts.tsv", cells)
    fixture.save()
    return fixture


def test_reference_controls_accept_correct_committed_list_and_source_metadata(
    tmp_path: Path,
) -> None:
    stool_root, site_root = tmp_path / "stool", tmp_path / "site"
    stool_root.mkdir()
    site_root.mkdir()
    stool = stool_fixture(stool_root)
    site = site_fixture(site_root)
    qualification.test_seventy_stool_cells_against_historical_numbers(stool.inputs())
    qualification.test_site_condition_cells_against_persisted_pairs(site.inputs())


@pytest.mark.parametrize("change", ["extra_trithionate", "reordered", "filename", "hash", "count"])
def test_qualification_rejects_wrong_trait_list_even_with_seventy_matching_reference_cells(
    tmp_path: Path, change: str
) -> None:
    fixture = stool_fixture(tmp_path, change)
    with pytest.raises(AssertionError):
        qualification.test_seventy_stool_cells_against_historical_numbers(fixture.inputs())


def test_qualification_rejects_changed_historical_reference_bytes(tmp_path: Path) -> None:
    fixture = stool_fixture(tmp_path)
    reference = json.loads(fixture.historical.read_text(encoding="utf-8"))
    reference["changed_after_pin"] = True
    fixture.historical.write_text(json.dumps(reference), encoding="utf-8")
    with pytest.raises(AssertionError):
        qualification.test_seventy_stool_cells_against_historical_numbers(fixture.inputs())


def test_historical_condition_display_labels_match_exact_consumer_names(tmp_path: Path) -> None:
    fixture = named_stool_fixture(tmp_path, display_labels=True)
    qualification.test_seventy_stool_cells_against_historical_numbers(fixture.inputs())


@pytest.mark.parametrize("canonical", HISTORICAL_DISPLAY_LABELS)
def test_historical_condition_alias_collision_is_rejected(tmp_path: Path, canonical: str) -> None:
    fixture = named_stool_fixture(tmp_path, display_labels=False)
    reference = json.loads(fixture.historical.read_text(encoding="utf-8"))
    original = next(row for row in reference["s14_tier12"] if row["name"] == canonical)
    reference["s14_tier12"].append({**original, "name": HISTORICAL_DISPLAY_LABELS[canonical]})
    save_reference(fixture, reference)
    with pytest.raises(AssertionError):
        qualification.test_seventy_stool_cells_against_historical_numbers(fixture.inputs())


@pytest.mark.parametrize("replace_existing", [False, True])
def test_unknown_historical_condition_alias_is_rejected(
    tmp_path: Path, replace_existing: bool
) -> None:
    fixture = named_stool_fixture(tmp_path, display_labels=False)
    reference = json.loads(fixture.historical.read_text(encoding="utf-8"))
    unknown = {**reference["s14_tier12"][0], "name": "Irritable bowel syndrome (IBS)"}
    if replace_existing:
        reference["s14_tier12"][0] = unknown
    else:
        reference["s14_tier12"].append(unknown)
    save_reference(fixture, reference)
    with pytest.raises(AssertionError):
        qualification.test_seventy_stool_cells_against_historical_numbers(fixture.inputs())


def test_historical_condition_aliases_do_not_normalize_consumer_output(tmp_path: Path) -> None:
    fixture = named_stool_fixture(tmp_path, display_labels=True)
    rows = qualification.table(fixture.root / "stool_disease_shifts.tsv")
    for row in rows:
        row["condition"] = HISTORICAL_DISPLAY_LABELS.get(row["condition"], row["condition"])
    write_table(fixture.root / "stool_disease_shifts.tsv", rows)
    fixture.save()
    with pytest.raises(AssertionError):
        qualification.test_seventy_stool_cells_against_historical_numbers(fixture.inputs())


@pytest.mark.parametrize("change", ["category_field", "site_relabel", "condition_relabel"])
def test_qualification_rejects_consistent_wrong_consumer_metadata(
    tmp_path: Path, change: str
) -> None:
    fixture = site_fixture(tmp_path, change)
    with pytest.raises(AssertionError):
        qualification.test_site_condition_cells_against_persisted_pairs(fixture.inputs())


def test_qualification_rejects_changed_core_experiment_bytes(tmp_path: Path) -> None:
    fixture = site_fixture(tmp_path)
    fixture.experiments.write_text(
        fixture.experiments.read_text(encoding="utf-8") + "\n", encoding="utf-8"
    )
    with pytest.raises(AssertionError):
        qualification.test_site_condition_cells_against_persisted_pairs(fixture.inputs())
