"""Selection, direction and source-projection regression tests for #658."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import pytest
import yaml

from mtp import contrast_orientation as orientation
from mtp.ingest import ingest
from mtp.paths import file_sha256
from mtp.workflow.models import ContrastPopulation
from mtp.workflow.population import POPULATIONS, select_contrasts, validate_population

ROOT = Path(__file__).resolve().parents[1]
HEADER = (
    "BSDB ID",
    "Study",
    "Title",
    "Year",
    "Study design",
    "PMID",
    "DOI",
    "Condition",
    "EFO ID",
    "Body site",
    "UBERON ID",
    "Host species",
    "Sequencing type",
    "Group 0 name",
    "Group 1 name",
    "Group 1 definition",
    "Group 0 sample size",
    "Group 1 sample size",
    "Description",
    "Abundance in Group 1",
    "NCBI Taxonomy IDs",
    "MetaPhlAn taxon names",
)


@dataclass
class Corpus:
    raw: Path
    source: Path
    output: Path
    audit: Path

    def select(self, population: str) -> dict[str, Any]:
        select_contrasts(
            {"source_studies": self.source, "raw_dump": self.raw},
            {"studies": self.output, "contrast_selection": self.audit},
            {"contrast_population": population},
        )
        value = json.loads(self.audit.read_text(encoding="utf-8"))
        assert isinstance(value, dict)
        return value

    def validate(self, population: str) -> dict[str, object]:
        return validate_population(self.output, self.audit, self.source, self.raw, population)


def _write_raw(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=HEADER)
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture
def corpus(tmp_path: Path) -> Corpus:
    examples = (
        ("1", "1", "Asthma", "Healthy controls", "Asthma patients", "Homo sapiens", True),
        ("1", "2", "Asthma", "Asthma patients", "Healthy controls", "Homo sapiens", True),
        ("1", "3", "HIV infection", "HIV negative", "HIV positive", "Homo sapiens", True),
        (
            "1",
            "4",
            "Colorectal cancer",
            "Tumor tissue",
            "Adjacent non-cancerous tissue",
            "Homo sapiens",
            True,
        ),
        (
            "1",
            "5",
            "Constipation",
            "Maintenance Hemodialysis without Constipation (MHDNCG)",
            "Maintenance Hemodialysis with Constipation (MHDCG)",
            "Homo sapiens",
            True,
        ),
        (
            "1",
            "6",
            "HIV infection",
            "Controls",
            "HIV Patients treated with Maraviroc",
            "Homo sapiens",
            True,
        ),
        ("2", "1", "Asthma", "Healthy controls", "Asthma patients", "Mus musculus", True),
        ("3", "1", "Asthma", "Healthy controls", "Asthma patients", "Homo sapiens", False),
        ("4", "1", "Asthma", "Placebo group", "Drug treatment group", "Homo sapiens", True),
        ("5", "1", "Asthma", "Healthy controls", "Asthma patients", "Homo sapiens", True),
        ("6", "1", "Asthma", "Healthy controls", "Asthma patients", "Homo sapiens", False),
    )
    # Every other experiment is 16S, which the profiles keep (#742)
    fungal = {("5", "1"): "ITS / ITS2", ("6", "1"): "ITS / ITS2"}
    rows = []
    for study, experiment, condition, group0, group1, host, paired in examples:
        directions = ["increased", "decreased"] if paired else ["increased"]
        if (study, experiment) == ("1", "2"):
            directions.append("increased")
        for index, direction in enumerate(directions, 1):
            row = dict.fromkeys(HEADER, "")
            row.update(
                {
                    "BSDB ID": f"bsdb:{study}/{experiment}/{index}",
                    "Study": f"Study {study}",
                    "Title": "Population fixture",
                    "Year": "2026",
                    "Study design": "case-control",
                    "Condition": condition,
                    "Body site": "stool",
                    "Host species": host,
                    "Sequencing type": fungal.get((study, experiment), "16S"),
                    "Group 0 name": group0,
                    "Group 1 name": group1,
                    "Group 1 definition": group1,
                    "Group 0 sample size": "1" if (study, experiment) == ("1", "1") else "12",
                    "Group 1 sample size": "15",
                    "Description": "Reported taxon occurrences",
                    "Abundance in Group 1": direction,
                    "NCBI Taxonomy IDs": "100;100",
                    "MetaPhlAn taxon names": "g__Examplegenus,g__Examplegenus",
                }
            )
            rows.append(row)
    result = Corpus(
        tmp_path / "raw.csv",
        tmp_path / "source",
        tmp_path / "selected",
        tmp_path / "selection.json",
    )
    _write_raw(result.raw, rows)
    ingest(
        study_dir=result.source,
        derived_dir=tmp_path / "ingest",
        full_dump_path=result.raw,
    )
    return result


def _read(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert isinstance(value, dict)
    return value


def _save(path: Path, value: object) -> None:
    path.write_text(yaml.safe_dump(value, sort_keys=False), encoding="utf-8")


def _selected_ids(audit: dict[str, Any]) -> set[str]:
    return {row["experiment_id"] for row in audit["selected_experiments"]}


def test_population_names_match_linkml_enum() -> None:
    assert set(POPULATIONS) == {population.value for population in ContrastPopulation}


def test_default_is_byte_identity_including_nonhuman_and_incomplete(corpus: Corpus) -> None:
    before = {path.name: path.read_bytes() for path in corpus.source.iterdir()}
    audit = corpus.select("all_paired")
    assert {path.name: path.read_bytes() for path in corpus.output.iterdir()} == before
    assert audit["source_experiments"] == audit["selected_experiments"]
    assert audit["reversed_experiments"] == []
    assert audit["output_direction_basis"] == "group_1_vs_group_0"
    assert audit["classifier"]["applied"] is False
    assert audit["counts"]["excluded_experiments"] == 0
    assert corpus.validate("all_paired")["population"] == "all_paired"


def test_identity_does_not_require_classifier_columns(tmp_path: Path) -> None:
    raw = ROOT / "tests/fixtures/workflow/full_dump.csv"
    source = tmp_path / "source"
    ingest(study_dir=source, derived_dir=tmp_path / "ingest", full_dump_path=raw)
    fixture = Corpus(raw, source, tmp_path / "selected", tmp_path / "selection.json")
    fixture.select("all_paired")
    fixture.validate("all_paired")


@pytest.mark.parametrize(
    ("population", "experiments", "reversed_ids"),
    [
        ("disease_vs_healthy_oriented", {1, 2, 3, 4, 5, 6}, {2, 4}),
        ("case_vs_control_subject", {1, 2, 3, 6}, {2}),
        ("disease_vs_healthy_tier1_2_subject_controls", {1, 2, 3}, {2}),
    ],
)
def test_exact_nested_populations_and_single_reversal(
    corpus: Corpus, population: str, experiments: set[int], reversed_ids: set[int]
) -> None:
    source_bytes = {path.name: path.read_bytes() for path in corpus.source.iterdir()}
    audit = corpus.select(population)
    assert _selected_ids(audit) == {f"bsdb:1-{number}" for number in experiments}
    assert {row["experiment_id"] for row in audit["reversed_experiments"]} == {
        f"bsdb:1-{number}" for number in reversed_ids
    }
    assert audit["input_direction_basis"] == "group_1_vs_group_0"
    assert audit["output_direction_basis"] == "disease_vs_comparison"
    assert audit["policy"]["trait_total_filter"] == "none"
    assert list(corpus.output.iterdir()) == [corpus.output / "bsdb_1.yaml"]
    original = _read(corpus.source / "bsdb_1.yaml")
    selected = _read(corpus.output / "bsdb_1.yaml")
    assert {key: value for key, value in selected.items() if key != "experiments"} == {
        key: value for key, value in original.items() if key != "experiments"
    }
    source_experiments = {row["id"]: row for row in original["experiments"]}
    for experiment in selected["experiments"]:
        before = source_experiments[experiment["id"]]
        assert {key: value for key, value in experiment.items() if key != "signatures"} == {
            key: value for key, value in before.items() if key != "signatures"
        }
        assert len(experiment["signatures"]) == len(before["signatures"])
        for old, new in zip(before["signatures"], experiment["signatures"], strict=True):
            assert {key: value for key, value in old.items() if key != "direction"} == {
                key: value for key, value in new.items() if key != "direction"
            }
            reversed_ = int(experiment["id"].rsplit("-", 1)[1]) in reversed_ids
            assert (old["direction"] != new["direction"]) == reversed_
            assert len(new["taxa"]) == 2
    assert source_bytes == {path.name: path.read_bytes() for path in corpus.source.iterdir()}
    assert audit["counts"]["source"]["experiments"] == 11
    assert audit["counts"]["selected"]["experiments"] == len(experiments)
    assert corpus.validate(population)["population"] == population


def test_tier2_weak_evidence_is_retained_and_treated_flag_is_not_a_separate_filter(
    corpus: Corpus,
) -> None:
    audit = corpus.select("case_vs_control_subject")
    decisions = {row["experiment_id"]: row for row in audit["experiments"]}
    assert "weak_comparison_evidence" in decisions["bsdb:1-3"]["flags"]
    assert decisions["bsdb:1-3"]["selected"]
    assert "treated_disease_group" in decisions["bsdb:1-6"]["flags"]
    assert decisions["bsdb:1-6"]["selected"]
    assert decisions["bsdb:2-1"]["exclusion_reason"] == "non_human"
    assert decisions["bsdb:3-1"]["exclusion_reason"] == "unpaired"
    assert decisions["bsdb:4-1"]["exclusion_reason"] == "contrast_class"
    # A human paired asthma contrast like bsdb:2-1, but sequenced by ITS (#742)
    assert decisions["bsdb:5-1"]["exclusion_reason"] == "sequencing_type"
    # An ITS experiment excluded for another reason keeps that reason (#761)
    assert decisions["bsdb:6-1"]["exclusion_reason"] == "unpaired"
    assert audit["schema_version"] == "contrast-selection-v2"
    assert decisions["bsdb:5-1"]["sequencing_type"] == "ITS / ITS2"
    assert audit["policy"]["sequencing_types"] == ["16S", "WMS"]


@pytest.mark.parametrize(
    "mutation",
    [
        "group_metadata",
        "taxon_loss",
        "direction",
        "signature_loss",
        "experiment_loss",
        "extra_file",
    ],
)
def test_projection_rejects_corruption_even_if_output_hash_is_rewritten(
    corpus: Corpus, mutation: str
) -> None:
    audit = corpus.select("disease_vs_healthy_oriented")
    path = corpus.output / "bsdb_1.yaml"
    study = _read(path)
    experiment = study["experiments"][0]
    if mutation == "group_metadata":
        experiment["group_0_name"] = "Changed group"
    elif mutation == "taxon_loss":
        experiment["signatures"][0]["taxa"].pop()
    elif mutation == "direction":
        experiment["signatures"][0]["direction"] = "decreased"
    elif mutation == "signature_loss":
        experiment["signatures"].pop()
    elif mutation == "experiment_loss":
        study["experiments"].pop()
    else:
        (corpus.output / ".hidden").write_text("unrecorded", encoding="utf-8")
    _save(path, study)
    audit["selected_files"][0]["sha256"] = file_sha256(path)
    corpus.audit.write_text(json.dumps(audit), encoding="utf-8")
    with pytest.raises(ValueError, match=r"projection differs|Unexpected member"):
        corpus.validate("disease_vs_healthy_oriented")


@pytest.mark.parametrize(
    ("field", "value"),
    [("group_0_sample_size", True), ("group_1_sample_size", 15.0)],
)
def test_projection_rejects_equal_value_with_wrong_yaml_type(
    corpus: Corpus, field: str, value: bool | float
) -> None:
    audit = corpus.select("disease_vs_healthy_oriented")
    path = corpus.output / "bsdb_1.yaml"
    study = _read(path)
    experiment = study["experiments"][0]
    assert type(experiment[field]) is int
    assert experiment[field] == value
    experiment[field] = value
    assert study == _read(path)
    _save(path, study)
    audit["selected_files"][0]["sha256"] = file_sha256(path)
    corpus.audit.write_text(json.dumps(audit), encoding="utf-8")
    with pytest.raises(ValueError, match="selected Study projection differs"):
        corpus.validate("disease_vs_healthy_oriented")


@pytest.mark.parametrize(
    "field",
    [
        "population",
        "schema_version",
        "raw_dump_sha256",
        "output_direction_basis",
        "classifier",
        "selected_experiments",
        "reversed_experiments",
        "source_experiments",
        "policy",
        "counts",
        "experiments",
    ],
)
def test_audit_is_recomputed_not_trusted(corpus: Corpus, field: str) -> None:
    audit = corpus.select("disease_vs_healthy_oriented")
    audit[field] = "tampered"
    corpus.audit.write_text(json.dumps(audit), encoding="utf-8")
    with pytest.raises(ValueError, match="audit differs"):
        corpus.validate("disease_vs_healthy_oriented")


@pytest.mark.parametrize(
    ("section", "field", "value"),
    [
        ("classifier", "applied", 1),
        ("policy", "strong_evidence_only", 0),
        ("selected_counts", "experiments", 6.0),
        ("experiment", "selected", 1),
    ],
)
def test_audit_rejects_equal_value_with_wrong_json_type(
    corpus: Corpus, section: str, field: str, value: int | float
) -> None:
    audit = corpus.select("disease_vs_healthy_oriented")
    if section == "selected_counts":
        target = audit["counts"]["selected"]
    elif section == "experiment":
        target = audit["experiments"][0]
    else:
        target = audit[section]
    assert target[field] == value
    assert type(target[field]) is not type(value)
    target[field] = value
    corpus.audit.write_text(json.dumps(audit), encoding="utf-8")
    with pytest.raises(ValueError, match="audit differs"):
        corpus.validate("disease_vs_healthy_oriented")


def test_wrong_expected_population_fails(corpus: Corpus) -> None:
    corpus.select("disease_vs_healthy_oriented")
    with pytest.raises(ValueError, match="projection differs"):
        corpus.validate("case_vs_control_subject")


def test_missing_raw_experiment_is_not_silently_excluded(corpus: Corpus) -> None:
    with corpus.raw.open(encoding="utf-8", newline="") as stream:
        rows = [row for row in csv.DictReader(stream) if not row["BSDB ID"].startswith("bsdb:2/")]
    _write_raw(corpus.raw, rows)
    with pytest.raises(ValueError, match="no raw export record"):
        corpus.select("disease_vs_healthy_oriented")
    assert not corpus.output.exists()
    assert not corpus.audit.exists()


def test_invalid_selected_orientation_fails(
    corpus: Corpus, monkeypatch: pytest.MonkeyPatch
) -> None:
    classify = orientation.classify

    def invalid(record: orientation.ExperimentRecord) -> orientation.Classification:
        result = classify(record)
        if record.experiment_id == "bsdb:1-1":
            return replace(result, orientation=orientation.NOT_APPLICABLE)
        return result

    monkeypatch.setattr(orientation, "classify", invalid)
    with pytest.raises(ValueError, match="invalid orientation"):
        corpus.select("disease_vs_healthy_oriented")
    assert not corpus.audit.exists()


@pytest.mark.parametrize("population", ["unknown", "strong_evidence_only"])
def test_undeclared_populations_fail(corpus: Corpus, population: str) -> None:
    with pytest.raises(ValueError, match="Unknown contrast population"):
        corpus.select(population)
    assert not corpus.output.exists()


def test_source_and_existing_output_are_protected(corpus: Corpus) -> None:
    original = {path.name: path.read_bytes() for path in corpus.source.iterdir()}
    corpus.output = corpus.source
    with pytest.raises(ValueError, match="empty regular directory"):
        corpus.select("all_paired")
    assert original == {path.name: path.read_bytes() for path in corpus.source.iterdir()}
    corpus.output = corpus.source.parent / "selected"
    corpus.select("all_paired")
    with pytest.raises(ValueError, match="new file"):
        corpus.select("all_paired")


def test_empty_population_fails_before_output(corpus: Corpus) -> None:
    for path in corpus.source.iterdir():
        value = _read(path)
        for experiment in value["experiments"]:
            experiment["host_species"] = "Mus musculus"
        _save(path, value)
    with pytest.raises(ValueError, match="selects no experiments"):
        corpus.select("case_vs_control_subject")
    assert not corpus.output.exists()


def test_population_names_are_explicit() -> None:
    assert POPULATIONS == (
        "all_paired",
        "disease_vs_healthy_oriented",
        "case_vs_control_subject",
        "disease_vs_healthy_tier1_2_subject_controls",
    )
