import json
from pathlib import Path

import pytest
import yaml

from bsdbng.ingest import (
    _MAX_INFO_EXAMPLES,
    _clean,
    _clean_float,
    _clean_int,
    _clean_list,
    _summarize_ingest_log,
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


def test_ingest_log_summarizes_repetitive_info_entries_by_default(tmp_path: Path) -> None:
    """Default ingest_log.json summarizes repeated info entries and keeps skips (#65)."""
    raw_dir = tmp_path / "raw"
    study_dir = tmp_path / "studies"
    derived_dir = tmp_path / "derived"
    raw_dir.mkdir()

    # Build a row with 7 three-level taxa and 2 four-level taxa, plus a skipped signature row
    three_level_ids = ";".join(f"1783272|91061|{1350 + idx}" for idx in range(7))
    three_level_names = ",".join(f"k__Bacillati|c__Bacilli|g__Genus{idx}" for idx in range(7))
    four_level_ids = ";".join(f"1783272|1239|1578|{1485 + idx}" for idx in range(2))
    four_level_names = ",".join(
        f"k__Bacillati|p__Bacillota|c__Clostridia|g__Clostridium{idx}" for idx in range(2)
    )
    row_many_taxa = (
        "bsdb:99/1/1,Study 99,case-control,12345678,https://doi.org/10.1234/test,,John Doe,"
        "Test Study Title,Test Journal,2024,keyword1,"
        "Experiment 1,USA,Homo sapiens,stool,UBERON:0001988,"
        "obesity,EFO:0001073,lean,obese,obese subjects,10,10,yes,"
        "16S,V4,Illumina,relative abundances,Mann-Whitney,0.05,"
        "BH,NA,age,BMI,,,,,,NA,"
        "Signature 1,Table S1,2024-01-15,TestCurator,TestEditor,Obesity-related taxa,"
        f"increased,{three_level_ids};{four_level_ids},"
        f'"{three_level_names},{four_level_names}",'
        "Complete,TestReviewer\n"
    )
    row_skipped_sig = (
        "bsdb:99/1/3,Study 99,case-control,12345678,https://doi.org/10.1234/test,,John Doe,"
        "Test Study Title,Test Journal,2024,keyword1,"
        "Experiment 1,USA,Homo sapiens,stool,UBERON:0001988,"
        "obesity,EFO:0001073,lean,obese,obese subjects,10,10,yes,"
        "16S,V4,Illumina,relative abundances,Mann-Whitney,0.05,"
        "BH,NA,age,BMI,,,,,,NA,"
        "Signature 3,Table S2,2024-01-15,TestCurator,TestEditor,Invalid direction,"
        "unchanged,5678,s__Bacteroides fragilis,Complete,TestReviewer\n"
    )
    (raw_dir / "full_dump.csv").write_text(
        FULL_DUMP_HEADER + row_many_taxa + FULL_DUMP_ROW_2 + row_skipped_sig,
        encoding="utf-8",
    )

    ingest(raw_dir, study_dir, derived_dir)
    events = json.loads((derived_dir / "ingest_log.json").read_text(encoding="utf-8"))

    assert isinstance(events, list)
    for event in events:
        assert isinstance(event, dict)
        assert all(isinstance(k, str) and isinstance(v, str) for k, v in event.items())

    skips = [e for e in events if e.get("level") == "skip"]
    assert skips == [
        {
            "level": "skip",
            "entity": "signature",
            "id": "bsdb:99-1-3",
            "reason": "unparsable direction: 'unchanged'",
        }
    ]

    info_events = [e for e in events if e.get("level") == "info"]
    assert len(info_events) == 3
    infos = {e["reason"]: e for e in info_events}
    assert set(infos) == {
        "lineage with 3 levels, used leaf",
        "lineage with 4 levels, used leaf",
        "DOI normalized from URL",
    }

    three_level_info = infos["lineage with 3 levels, used leaf"]
    assert three_level_info.get("entity") == "taxon"
    assert three_level_info.get("count") == "7"
    examples_3 = three_level_info.get("examples", "").split(", ")
    assert len(examples_3) == _MAX_INFO_EXAMPLES
    assert examples_3[0] == "bsdb:99-1-1:NCBITaxon:1350"
    assert examples_3[-1] == f"bsdb:99-1-1:NCBITaxon:{1350 + _MAX_INFO_EXAMPLES - 1}"

    four_level_info = infos["lineage with 4 levels, used leaf"]
    assert four_level_info.get("entity") == "taxon"
    assert four_level_info.get("count") == "2"
    assert four_level_info.get("examples") == (
        "bsdb:99-1-1:NCBITaxon:1485, bsdb:99-1-1:NCBITaxon:1486"
    )

    doi_info = infos["DOI normalized from URL"]
    assert doi_info.get("entity") == "study"
    assert doi_info.get("count") == "1"
    assert doi_info.get("examples") == ("bsdb:99: 'https://doi.org/10.1234/test' → '10.1234/test'")


def test_ingest_verbose_log_preserves_individual_info_entries(tmp_path: Path) -> None:
    """Passing verbose_log=True writes uncollapsed per-taxon info entries (#65)."""
    raw_dir = tmp_path / "raw"
    study_dir = tmp_path / "studies"
    derived_dir = tmp_path / "derived"
    raw_dir.mkdir()
    _write_full_dump(raw_dir)

    ingest(raw_dir, study_dir, derived_dir, verbose_log=True)
    events = json.loads((derived_dir / "ingest_log.json").read_text(encoding="utf-8"))
    assert events == [
        {
            "level": "info",
            "entity": "taxon",
            "id": "NCBITaxon:1350",
            "signature": "bsdb:99-1-1",
            "reason": "lineage with 3 levels, used leaf",
        },
        {
            "level": "info",
            "entity": "taxon",
            "id": "NCBITaxon:1485",
            "signature": "bsdb:99-1-1",
            "reason": "lineage with 4 levels, used leaf",
        },
    ]


def test_summarize_ingest_log_validates_max_examples_and_deduplicates_examples() -> None:
    """_summarize_ingest_log rejects invalid max_examples and deduplicates repeated example IDs."""
    with pytest.raises(ValueError, match="max_examples must be at least 1"):
        _summarize_ingest_log([], max_examples=0)

    raw_events = [
        {
            "level": "info",
            "entity": "taxon",
            "id": "NCBITaxon:100",
            "signature": "bsdb:1-1-1",
            "reason": "lineage with 2 levels, used leaf",
        },
        {
            "level": "info",
            "entity": "taxon",
            "id": "NCBITaxon:100",
            "signature": "bsdb:1-1-1",
            "reason": "lineage with 2 levels, used leaf",
        },
        {
            "level": "info",
            "entity": "taxon",
            "id": "NCBITaxon:101",
            "signature": "bsdb:1-1-2",
            "reason": "lineage with 2 levels, used leaf",
        },
    ]
    summarized = _summarize_ingest_log(raw_events, max_examples=2)
    assert summarized == [
        {
            "level": "info",
            "entity": "taxon",
            "reason": "lineage with 2 levels, used leaf",
            "count": "3",
            "examples": "bsdb:1-1-1:NCBITaxon:100, bsdb:1-1-2:NCBITaxon:101",
        }
    ]


@pytest.mark.parametrize(
    ("reason", "values"),
    [
        ("non-numeric publication year dropped", ["'unknown'", "'2020?'", "'in press'"]),
        ("non-numeric PMID dropped", ["'pmid:123'", "'pending'", "'12x'"]),
        (
            "DOI normalized from URL",
            [
                "'https://doi.org/10.1234/one' → '10.1234/one'",
                "'http://doi.org/10.1234/two' → '10.1234/two'",
                "'https://doi.org/10.1234/three' → '10.1234/three'",
            ],
        ),
        ("unrecognized DOI format dropped", ["'pending'", "'doi:123'", "'not, a DOI'"]),
        ("non-HTTP URL dropped", ["'ftp://one'", "'www.two.org'", "'see: supplement'"]),
    ],
)
def test_summarize_study_info_collapses_distinct_raw_values(reason: str, values: list[str]) -> None:
    """Repeated study normalizations retain each ID and value under one reason (#2204)."""
    events = [
        {"level": "info", "entity": "study", "id": f"bsdb:{i}", "reason": f"{reason}: {value}"}
        for i, value in enumerate(values, start=1)
    ]
    original = [dict(event) for event in events]
    assert _summarize_ingest_log(events) == [
        {
            "level": "info",
            "entity": "study",
            "reason": reason,
            "count": "3",
            "examples": ", ".join(f"bsdb:{i}: {value}" for i, value in enumerate(values, start=1)),
        }
    ]
    assert events == original


def test_summarize_study_info_caps_examples_without_losing_distinct_values() -> None:
    """Count every event; cap unique ID/value pairs, including repeats within a study."""
    pairs = [(1, "first"), (1, "first"), (1, "second"), *[(i, f"value{i}") for i in range(2, 8)]]
    events = [
        {
            "level": "info",
            "entity": "study",
            "id": f"bsdb:{study}",
            "reason": f"non-numeric PMID dropped: {value!r}",
        }
        for study, value in pairs
    ]
    summaries = _summarize_ingest_log(events, max_examples=2)
    assert summaries == [
        {
            "level": "info",
            "entity": "study",
            "reason": "non-numeric PMID dropped",
            "count": "9",
            "examples": "bsdb:1: 'first', bsdb:1: 'second'",
        }
    ]
    assert _summarize_ingest_log(events)[0]["examples"] == (
        "bsdb:1: 'first', bsdb:1: 'second', bsdb:2: 'value2', bsdb:3: 'value3', bsdb:4: 'value4'"
    )


def test_summarize_study_info_preserves_other_entities_reasons_and_levels() -> None:
    """The five study-info templates do not generalize to other entries or mutate input."""
    reason = "non-HTTP URL dropped: 'ftp://one'"
    unchanged = [
        {"level": level, "entity": "study", "id": "bsdb:1", "reason": reason}
        for level in ("warning", "error", "skip")
    ]
    other_info = [
        {"level": "info", "entity": "taxon", "id": "NCBITaxon:1", "reason": reason},
        {"level": "info", "entity": "study", "id": "bsdb:2", "reason": "custom reason: one"},
        {"level": "info", "entity": "study", "id": "bsdb:3", "reason": "custom reason: two"},
        {
            "level": "info",
            "entity": "study",
            "id": "bsdb:4",
            "reason": "non-HTTP URL dropped elsewhere: 'ftp://two'",
        },
    ]
    raw = [other_info[0], *unchanged, *other_info[1:]]
    original = [dict(event) for event in raw]
    assert _summarize_ingest_log(raw) == [
        *unchanged,
        *[
            {
                "level": "info",
                "entity": event["entity"],
                "reason": event["reason"],
                "count": "1",
                "examples": event["id"],
            }
            for event in other_info
        ],
    ]
    assert raw == original


@pytest.mark.parametrize("verbose_log", [False, True], ids=["summary", "verbose"])
def test_ingest_multiple_study_doi_normalizations(tmp_path: Path, verbose_log: bool) -> None:
    """Real ingestion collapses distinct DOI events only in the nonverbose log (#2204)."""
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    rows = [
        FULL_DUMP_ROW_1.replace("bsdb:99/", f"bsdb:{i}/")
        .replace("Study 99", f"Study {i}")
        .replace("10.1234/test", f"https://doi.org/10.1234/study{i}")
        for i in range(1, 4)
    ]
    (raw_dir / "full_dump.csv").write_text(FULL_DUMP_HEADER + "".join(rows), encoding="utf-8")
    derived = tmp_path / "derived"
    written = ingest(raw_dir, tmp_path / "studies", derived, verbose_log=verbose_log)
    assert len(written) == 3
    assert [yaml.safe_load(path.read_text())["doi"] for path in written] == [
        "10.1234/study1",
        "10.1234/study2",
        "10.1234/study3",
    ]
    events = json.loads((derived / "ingest_log.json").read_text(encoding="utf-8"))
    study_info = [e for e in events if e["level"] == "info" and e["entity"] == "study"]
    details = [f"'https://doi.org/10.1234/study{i}' → '10.1234/study{i}'" for i in range(1, 4)]
    if verbose_log:
        assert study_info == [
            {
                "level": "info",
                "entity": "study",
                "id": f"bsdb:{i}",
                "reason": f"DOI normalized from URL: {value}",
            }
            for i, value in enumerate(details, start=1)
        ]
    else:
        assert study_info == [
            {
                "level": "info",
                "entity": "study",
                "reason": "DOI normalized from URL",
                "count": "3",
                "examples": ", ".join(
                    f"bsdb:{i}: {value}" for i, value in enumerate(details, start=1)
                ),
            }
        ]


def test_summarize_study_info_keeps_mixed_groups_separate() -> None:
    """Interleaved templates, bare reasons, taxa and warnings retain independent groups."""
    warning = {"level": "warning", "entity": "study", "id": "bsdb:1", "reason": "keep warning"}
    events = [
        {
            "level": "info",
            "entity": "study",
            "id": "bsdb:1",
            "reason": "non-numeric PMID dropped: 'unknown'",
        },
        {
            "level": "info",
            "entity": "taxon",
            "id": "NCBITaxon:1",
            "reason": "lineage with 2 levels, used leaf",
        },
        {
            "level": "info",
            "entity": "study",
            "id": "bsdb:1",
            "reason": "non-HTTP URL dropped: 'ftp://one'",
        },
        warning,
        {
            "level": "info",
            "entity": "study",
            "id": "bsdb:2",
            "reason": "non-numeric PMID dropped: 'unknown'",
        },
        {
            "level": "info",
            "entity": "study",
            "id": "bsdb:3",
            "reason": "non-numeric PMID dropped: 'pending'",
        },
        {
            "level": "info",
            "entity": "study",
            "id": "bsdb:2",
            "reason": "non-HTTP URL dropped: 'ftp://two'",
        },
        {"level": "info", "entity": "study", "id": "bsdb:4", "reason": "non-HTTP URL dropped"},
    ]
    assert _summarize_ingest_log(events) == [
        warning,
        {
            "level": "info",
            "entity": "study",
            "reason": "non-numeric PMID dropped",
            "count": "3",
            "examples": "bsdb:1: 'unknown', bsdb:2: 'unknown', bsdb:3: 'pending'",
        },
        {
            "level": "info",
            "entity": "taxon",
            "reason": "lineage with 2 levels, used leaf",
            "count": "1",
            "examples": "NCBITaxon:1",
        },
        {
            "level": "info",
            "entity": "study",
            "reason": "non-HTTP URL dropped",
            "count": "2",
            "examples": "bsdb:1: 'ftp://one', bsdb:2: 'ftp://two'",
        },
        {
            "level": "info",
            "entity": "study",
            "reason": "non-HTTP URL dropped",
            "count": "1",
            "examples": "bsdb:4",
        },
    ]
