from pathlib import Path

import pytest
import yaml

from bsdbng.ingest import (
    _clean,
    _clean_float,
    _clean_int,
    _clean_list,
    assert_required_exports,
    ingest,
)


def test_assert_required_exports_reports_missing_files(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match=r"experiments\.csv"):
        assert_required_exports(tmp_path)


# --- clean helper tests ---


def test_clean_strips_and_nones() -> None:
    assert _clean("  hello  ") == "hello"
    assert _clean("NA") is None
    assert _clean("") is None
    assert _clean("  ") is None


def test_clean_int() -> None:
    assert _clean_int("42") == 42
    assert _clean_int("NA") is None
    assert _clean_int("") is None
    assert _clean_int("3.5") is None


def test_clean_float() -> None:
    assert _clean_float("0.05") == 0.05
    assert _clean_float("NA") is None
    assert _clean_float("") is None
    assert _clean_float("not a number") is None


def test_clean_list() -> None:
    assert _clean_list("UBERON:0001988,UBERON:0002114") == ["UBERON:0001988", "UBERON:0002114"]
    assert _clean_list("EFO:0001073") == ["EFO:0001073"]
    assert _clean_list("NA") is None
    assert _clean_list("") is None


# --- ingest tests ---

FULL_DUMP_HEADER = (
    "# BugSigDB 2026-04-16, License: CC-BY-4.0\n"
    "BSDB ID,Study,Study design,PMID,DOI,URL,Authors list,Title,Journal,Year,"
    "Keywords,Experiment,Location of subjects,Host species,Body site,UBERON ID,"
    "Condition,EFO ID,Group 0 name,Group 1 name,Group 1 definition,"
    "Group 0 sample size,Group 1 sample size,Antibiotics exclusion,"
    "Sequencing type,16S variable region,Sequencing platform,"
    "Data transformation,Statistical test,Significance threshold,"
    "MHT correction,LDA Score above,Matched on,Confounders controlled for,"
    "Pielou,Shannon,Chao1,Simpson,Inverse Simpson,Richness,"
    "Signature page name,Source,Curated date,Curator,Revision editor,Description,"
    "Abundance in Group 1,NCBI Taxonomy IDs,MetaPhlAn taxon names,State,Reviewer\n"
)

# Row 1: two taxa (semicolon-separated), each with lineage (pipe-separated)
FULL_DUMP_ROW_1 = (
    "bsdb:99/1/1,Study 99,case-control,12345678,10.1234/test,,John Doe,"
    "Test Study Title,Test Journal,2024,keyword1,"
    "Experiment 1,USA,Homo sapiens,stool,UBERON:0001988,"
    "obesity,EFO:0001073,lean,obese,obese subjects,10,10,yes,"
    "16S,V4,Illumina,relative abundances,Mann-Whitney,0.05,"
    "BH,NA,age,BMI,,,,,,NA,"
    "Signature 1,Table S1,2024-01-15,TestCurator,TestEditor,Obesity-related taxa,"
    "increased,1783272|91061|1350;1783272|1239|1578|1485,"
    '"k__Bacillati|c__Bacilli|g__Enterococcus,k__Bacillati|p__Bacillota|c__Clostridia|g__Clostridium",'
    "Complete,TestReviewer\n"
)

# Row 2: single taxon, no lineage
FULL_DUMP_ROW_2 = (
    "bsdb:99/1/2,Study 99,case-control,12345678,10.1234/test,,John Doe,"
    "Test Study Title,Test Journal,2024,keyword1,"
    "Experiment 1,USA,Homo sapiens,stool,UBERON:0001988,"
    "obesity,EFO:0001073,lean,obese,obese subjects,10,10,yes,"
    "16S,V4,Illumina,relative abundances,Mann-Whitney,0.05,"
    "BH,NA,age,BMI,,,,,,NA,"
    "Signature 2,Fig 3,2024-01-15,TestCurator,TestEditor,Decreased in obesity,"
    "decreased,5678,"
    "s__Bacteroides fragilis,"
    "Complete,TestReviewer\n"
)


def _write_full_dump(raw_dir: Path) -> None:
    (raw_dir / "full_dump.csv").write_text(FULL_DUMP_HEADER + FULL_DUMP_ROW_1 + FULL_DUMP_ROW_2)
    for name in ("studies.csv", "experiments.csv", "signatures.csv"):
        (raw_dir / name).touch()


def test_ingest_produces_yaml(tmp_path: Path) -> None:
    raw_dir = tmp_path / "raw"
    study_dir = tmp_path / "studies"
    raw_dir.mkdir()
    _write_full_dump(raw_dir)

    written = ingest(raw_dir, study_dir, tmp_path / "derived")
    assert len(written) == 1

    study = yaml.safe_load(written[0].read_text())

    # Study-level fields
    assert study["id"] == "bsdb:99"
    assert study["pmid"] == 12345678
    assert study["publication_year"] == 2024
    assert study["doi"] == "10.1234/test"
    assert study["study_design"] == "case-control"
    assert study["authors_list"] == "John Doe"
    assert study["journal"] == "Test Journal"
    assert study["keywords"] == "keyword1"
    assert len(study["experiments"]) == 1

    # Experiment-level fields
    exp = study["experiments"][0]
    assert exp["id"] == "bsdb:99-1"
    assert exp["host_species"] == "Homo sapiens"
    assert exp["body_site"] == "stool"
    assert exp["body_site_ontology_id"] == ["UBERON:0001988"]
    assert exp["condition"] == "obesity"
    assert exp["condition_ontology_id"] == ["EFO:0001073"]
    assert exp["location_of_subjects"] == "USA"
    assert exp["group_0_sample_size"] == 10
    assert exp["group_1_sample_size"] == 10
    assert exp["group_1_definition"] == "obese subjects"
    assert exp["antibiotics_exclusion"] == "yes"
    assert exp["sequencing_type"] == "16S"
    assert exp["variable_region_16s"] == "V4"
    assert exp["sequencing_platform"] == "Illumina"
    assert exp["significance_threshold"] == 0.05
    assert exp["statistical_test"] == "Mann-Whitney"
    assert exp["matched_on"] == "age"
    assert exp["confounders_controlled_for"] == "BMI"
    assert exp["data_transformation"] == "relative abundances"
    assert exp["mht_correction"] == "BH"
    assert exp["lda_score_above"] is None
    assert exp["pielou"] is None
    assert exp["shannon"] is None
    assert exp["richness"] is None
    assert len(exp["signatures"]) == 2

    # Signature-level fields
    sig_increased = next(s for s in exp["signatures"] if s["direction"] == "increased")
    assert sig_increased["signature_source"] == "Table S1"
    assert sig_increased["signature_description"] == "Obesity-related taxa"
    assert len(sig_increased["taxa"]) == 2
    assert sig_increased["taxa"][0]["id"] == "NCBITaxon:1350"
    assert sig_increased["taxa"][0]["taxon_name"] == "Enterococcus"
    assert sig_increased["taxa"][0]["taxonomic_rank"] == "genus"
    assert sig_increased["taxa"][1]["id"] == "NCBITaxon:1485"
    assert sig_increased["taxa"][1]["taxon_name"] == "Clostridium"

    sig_decreased = next(s for s in exp["signatures"] if s["direction"] == "decreased")
    assert sig_decreased["signature_source"] == "Fig 3"
    assert sig_decreased["signature_description"] == "Decreased in obesity"
    assert len(sig_decreased["taxa"]) == 1
    assert sig_decreased["taxa"][0]["id"] == "NCBITaxon:5678"
    assert sig_decreased["taxa"][0]["taxon_name"] == "Bacteroides fragilis"
    assert sig_decreased["taxa"][0]["taxonomic_rank"] == "species"


def test_ingest_no_placeholder_taxon_names(tmp_path: Path) -> None:
    raw_dir = tmp_path / "raw"
    study_dir = tmp_path / "studies"
    raw_dir.mkdir()
    _write_full_dump(raw_dir)

    written = ingest(raw_dir, study_dir, tmp_path / "derived")
    assert written, "ingest() should write at least one YAML study file"

    for path in written:
        study = yaml.safe_load(path.read_text(encoding="utf-8"))
        for exp in study.get("experiments", []):
            for sig in exp.get("signatures", []):
                for taxon in sig.get("taxa", []):
                    name = taxon.get("taxon_name", "")
                    assert not name.startswith("taxon_"), (
                        f"{path.name}: placeholder taxon name {name!r} — "
                        f"name delimiter may be wrong"
                    )


def test_ingest_yaml_round_trips_through_pydantic(tmp_path: Path) -> None:
    raw_dir = tmp_path / "raw"
    study_dir = tmp_path / "studies"
    raw_dir.mkdir()
    _write_full_dump(raw_dir)

    written = ingest(raw_dir, study_dir, tmp_path / "derived")
    study_data = yaml.safe_load(written[0].read_text())

    from bsdbng.datamodel import Study

    study = Study.model_validate(study_data)
    assert study.id == "bsdb:99"
    assert study.pmid == 12345678
