"""Inventory of measured evidence for independent validation (#191)."""

from __future__ import annotations

from pathlib import Path

import pytest

from bsdbng.validation_evidence import (
    EVIDENCE_TYPES,
    ExternalCohort,
    StudyRecord,
    evidence_rows,
    external_overlap,
    find_mentions,
    inferred_function_traits,
    marker_stratified_experiments,
    parse_supplement,
    read_studies,
    summary_rows,
    title_agreement,
)

HEADER = (
    "BSDB ID,Study,Study design,PMID,DOI,URL,Authors list,Title,Journal,Year,Keywords,Experiment,"
    "Location of subjects,Host species,Body site,UBERON ID,Condition,EFO ID,Group 0 name,"
    "Group 1 name,Group 1 definition,Group 0 sample size,Group 1 sample size,Antibiotics exclusion,"
    "Sequencing type,Matched on,Confounders controlled for,Description"
)


def _dump(path: Path, rows: list[dict[str, str]]) -> Path:
    names = HEADER.split(",")
    lines = ["# BugSigDB test", HEADER]
    for row in rows:
        lines.append(",".join('"' + row.get(name, "NA").replace('"', "") + '"' for name in names))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _row(**changes: str) -> dict[str, str]:
    base = {
        "BSDB ID": "bsdb:1/1/1",
        "PMID": "1",
        "DOI": "NA",
        "Authors list": "Li X, Wang Y",
        "Title": "Gut microbiome in disease",
        "Year": "2020",
        "Keywords": "NA",
        "Location of subjects": "China",
        "Host species": "Homo sapiens",
        "Body site": "Feces",
        "Condition": "Colorectal cancer",
        "Group 0 name": "Healthy controls",
        "Group 1 name": "Patients",
        "Group 1 definition": "NA",
        "Sequencing type": "16S",
        "Antibiotics exclusion": "NA",
        "Matched on": "NA",
        "Confounders controlled for": "NA",
        "Description": "NA",
    }
    base.update(changes)
    return base


def test_studies_pool_text_across_rows_and_drop_na(tmp_path: Path) -> None:
    path = _dump(
        tmp_path / "dump.csv",
        [
            _row(Description="faecal calprotectin measured"),
            _row(**{"BSDB ID": "bsdb:1/2/1", "Description": "NA", "Sequencing type": "WMS"}),
        ],
    )
    studies = read_studies(path)
    assert list(studies) == ["bsdb:1"]
    record = studies["bsdb:1"]
    assert record.texts["Description"] == ["faecal calprotectin measured"]
    assert record.sequencing_types == "16S;WMS"
    assert record.antibiotics_exclusion_field == ""


def test_mentions_report_the_field_and_the_word() -> None:
    record = StudyRecord(
        study_id="bsdb:1",
        pmid="1",
        doi="",
        title="",
        authors="",
        year="2020",
        location="China",
        host_species="Homo sapiens",
        sequencing_types="16S",
        antibiotics_exclusion_field="",
        texts={"Description": ["Faecal calprotectin and IL-6 were measured"], "Title": ["A study"]},
    )
    mention = find_mentions(record, EVIDENCE_TYPES["inflammation_marker"])
    assert mention.status == "mentioned_in_curated_metadata"
    assert mention.fields == "Description"
    assert mention.terms == "calprotectin;il-6"
    assert find_mentions(record, EVIDENCE_TYPES["transit"]).status == (
        "not_mentioned_in_curated_metadata"
    )


def test_rows_cover_human_studies_and_flag_the_corpus_and_modality(tmp_path: Path) -> None:
    path = _dump(
        tmp_path / "dump.csv",
        [
            _row(Description="short-chain fatty acids quantified", **{"Sequencing type": "WMS"}),
            _row(**{"BSDB ID": "bsdb:2/1/1", "Host species": "Mus musculus"}),
        ],
    )
    rows = evidence_rows(read_studies(path), {"bsdb:1"})
    assert [row["study_id"] for row in rows] == ["bsdb:1"]
    assert rows[0]["in_discovery_corpus"] is True
    assert rows[0]["metagenomic_reads_possible"] is True
    assert rows[0]["scfa_metabolite_status"] == "mentioned_in_curated_metadata"
    summary = {row["evidence"]: row for row in summary_rows(rows)}
    assert summary["scfa_metabolite"]["discovery_studies_mentioning"] == 1
    assert summary["inflammation_marker"]["discovery_studies_mentioning"] == 0


def test_marker_conditions_are_measurements_not_drugs(tmp_path: Path) -> None:
    path = _dump(
        tmp_path / "dump.csv",
        [
            _row(**{"BSDB ID": "bsdb:1/1/1", "Condition": "Calprotectin measurement"}),
            _row(**{"BSDB ID": "bsdb:2/1/1", "Condition": "Medroxyprogesterone acetate"}),
            _row(**{"BSDB ID": "bsdb:3/1/1", "Condition": "Clostridium butyricum"}),
        ],
    )
    found = marker_stratified_experiments(path, {("bsdb:1", "bsdb:1-1")})
    assert [row["study_id"] for row in found] == ["bsdb:1"]
    assert found[0]["in_discovery_corpus"] is True
    assert "not independent" in str(found[0]["independence"])


def test_function_traits_report_what_the_annotation_asserts(tmp_path: Path) -> None:
    profiles = tmp_path / "profiles.tsv"
    profiles.write_text(
        "study_id\texperiment_id\tinc_dec\tMETPO:2000006|CHEBI:17272\tCHEBI:16865"
        "\tCHEBI:18401\tMETPO:2000006|mediadive.medium:1\n",
        encoding="utf-8",
    )
    nodes = tmp_path / "nodes.tsv"
    nodes.write_text(
        "id\tcategory\tname\n"
        "METPO:2000006\tx\tuses as carbon source\n"
        "CHEBI:17272\tx\tpropionate\n"
        "CHEBI:16865\tx\tgamma-aminobutyric acid\n"
        "CHEBI:18401\tx\tphenylacetate\n",
        encoding="utf-8",
    )
    rows = inferred_function_traits(profiles, nodes)
    assert [row["term"] for row in rows] == ["CHEBI:17272"]
    assert rows[0]["qualifier_label"] == "uses as carbon source"


@pytest.mark.parametrize(
    ("left", "right", "matches"),
    [
        (
            "Ametagenome-wide association study of gut microbiota in type 2 diabetes",
            "A metagenome-wide association study of gut microbiota in type 2 diabetes",
            True,
        ),
        (
            "Alterations of the gut microbiome in hypertension",
            "Alterations of the human gut microbiome in liver cirrhosis",
            False,
        ),
    ],
)
def test_title_agreement_matches_the_same_paper_only(left: str, right: str, matches: bool) -> None:
    assert (title_agreement(left, right) >= 0.9) is matches


TITLE = "Ametagenome-wide association study of gut microbiota in type 2 diabetes"


def _cohort(study_id: str = "QinJ_2012", reference_title: str = TITLE) -> ExternalCohort:
    return ExternalCohort(
        study_id=study_id,
        reference_title=reference_title,
        year="2012",
        disease_name="type 2 diabetes",
        abbreviation="T2D",
        accession="PRJNA422434",
        region="Guangdong",
        antibiotic_usage="N",
        probiotics_usage="no data",
        medication_usage="no data",
        projects="QinJ_2012.T2D",
        samples=369,
        disease_samples=186,
        control_samples=183,
        role="Discover",
        samples_with_age=369,
        samples_with_bmi=1,
    )


def test_overlap_is_decided_by_title_and_states_when_it_cannot_be(tmp_path: Path) -> None:
    path = _dump(
        tmp_path / "dump.csv",
        [
            _row(
                **{
                    "BSDB ID": "bsdb:9/1/1",
                    "Title": (
                        "A metagenome-wide association study of gut microbiota in type 2 diabetes"
                    ),
                }
            )
        ],
    )
    studies = read_studies(path)
    matched, untitled = external_overlap(
        [_cohort(), _cohort(study_id="LiS_2021a", reference_title="")], studies, {"bsdb:9"}
    )
    assert matched["overlap_status"] == "overlaps_discovery_corpus"
    assert matched["matched_bugsigdb_studies"] == "bsdb:9"
    assert matched["match_basis"] == "title"
    assert matched["measured_inflammation_or_scfa_available"] is False
    assert untitled["overlap_status"] == "not_assessable_without_title"
    outside = external_overlap([_cohort()], studies, set())[0]
    assert outside["overlap_status"] == "curated_in_bugsigdb_outside_discovery_corpus"


def test_supplement_parsing_joins_wrapped_titles_and_stops_at_the_next_table() -> None:
    text = (
        "Study ID      Reference\n"
        "QinJ_2012     Ametagenome-wide association study of gut microbiota\n"
        "              in type 2 diabetes\n"
        "Comparison ID               Observed number of species\n"
        "QinJ_2012.T2D    0.123    0.456\n"
    )
    parsed = parse_supplement(text)
    assert parsed["QinJ_2012"]["reference_title"].endswith("in type 2 diabetes")
    assert "0.123" not in parsed["QinJ_2012"]["reference_title"]
