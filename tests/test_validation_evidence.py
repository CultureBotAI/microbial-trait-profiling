"""Inventory of measured evidence for independent validation (#191)."""

from __future__ import annotations

from pathlib import Path

import pytest

from bsdbng.validation_evidence import (
    EVIDENCE_TYPES,
    INTERVENTION_CONTRAST,
    MARKER_NAMED_NOT_MEASURED,
    MARKER_STRATIFIED,
    ExternalCohort,
    StudyRecord,
    classify_marker_split,
    evidence_rows,
    external_overlap,
    find_mentions,
    inferred_function_traits,
    marker_split_candidates,
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


def _source(bsdb_id: str, experiment: str) -> dict[str, str]:
    """A dump row carrying one real experiment's condition and group fields."""
    condition, group_0, group_1, definition = SOURCE_ROWS[experiment]
    return {
        "BSDB ID": bsdb_id,
        "Condition": condition,
        "Group 0 name": group_0,
        "Group 1 name": group_1,
        "Group 1 definition": definition,
    }


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


# Condition, Group 0 name, Group 1 name and Group 1 definition exactly as the pinned
# BugSigDBExports v1.3.1 full_dump.csv carries them, for the experiments named in #292.
# Verified against data/raw/full_dump.csv (BugSigDB 2026-04-24_00:41_UTC) - the typo
# "fecal protectin", the doubled spaces and the misspelt "granulocute" are the source's.
SOURCE_ROWS: dict[str, tuple[str, str, str, str]] = {
    # Missed: the IBD calprotectin threshold split #191 is about. The condition names
    # no marker at all; the group names carry it.
    "bsdb:35888992-6": (
        "Increased inflammatory response",
        "Inflammatory bowel disease patients with fecal calprotectin levels < 1000 mg/kg",
        "Inflammatory bowel disease patients with fecal calprotectin levels > 1000 mg/kg",
        "Inflammatory bowel disease patients with higher fecal protectin levels "
        "(31 with CD and 12 with UC had above 1000 mg/kg).",
    ),
    # Missed: a presence/absence cut on systemic inflammation, in an RCT whose other
    # experiments are trial arms.
    "bsdb:38890302-6": (
        "Inflammation",
        "No inflammation",
        "Inflammation",
        "This group includes children grouped according to  systemic inflammation",
    ),
    # Missed: defined by measured cytokine levels, so the marker is in the definition.
    "bsdb:33220709-1": (
        "Vaginal inflammation",
        "Low Genital Inflammation",
        "High Genital Inflammation",
        "Women with high levels of genital inflammation, defined by elevated "
        "inflammatory cytokine levels in their cervicovaginal secretions.",
    ),
    # Missed: a non-interleukin cytokine row of bsdb:39103909.
    "bsdb:39103909-4": (
        "Monocyte chemotactic protein 1 measurement",
        "Decreased MCP-1 (monocyte chemoattractant protein 1)",
        "Increased MCP-1 (monocyte chemoattractant protein 1)",
        "Increased levels of monocyte chemoattractant protein 1 (MCP-1)",
    ),
    # Missed: and another, a colony-stimulating factor rather than an interleukin.
    "bsdb:39103909-13": (
        "Granulocyte colony-stimulating factor measurement",
        "Decreased GCSF (granulocute colony-stimulating factor)",
        "Increased GCSF (granulocute colony-stimulating factor)",
        "Increased levels of granulocute colony-stimulating factor (GCSF)",
    ),
    # Kept, and must stay kept: the condition names a transplant, but the cut is on
    # measured fecal butyrate, so the intervention test must not read the condition.
    "bsdb:31065565-2": (
        "Response to allogeneic hematopoietic stem cell transplant",
        "low fecal butyrate levels",
        "high fecal butyrate levels",
        "patients undergoing T cell repleted allogeneic hematopoietic stem cell transplantation",
    ),
    # Wrongly admitted: oxygen therapy is an intervention.
    "bsdb:879-4": (
        "Oxygen",
        "No oxygen therapy",
        "Oxygen therapy",
        "COVID-19 patients who had oxygen therapy.",
    ),
    # Wrongly admitted: baseline against a placebo arm, not a measured split.
    "bsdb:39748236-1": (
        "Oxidative Stress",
        "Baseline",
        "Placebo",
        "Healthy male soccer players taking placebo pills blindly for the first period of 12weeks.",
    ),
    # Wrongly admitted: baseline against a probiotic arm.
    "bsdb:39748236-2": (
        "Oxidative Stress",
        "Baseline",
        "Probiotic",
        "Healthy male soccer players taking probiotic pills blindly for the second "
        "period of 12 weeks.",
    ),
    # Must not be admitted: the condition says "measurement", the groups are disease
    # status, so the condition label alone cannot decide.
    "bsdb:37180439-9": (
        "C-reactive protein measurement",
        "Healthy Controls",
        "Patients with Metabolic dysfunction-associated liver disease",
        "Patients diagnosed with Non-alcoholic fatty liver disease (now known as "
        'Metabolic dysfunction-associated fatty liver disease "MAFLD")',
    ),
    # Must not be admitted: the condition names a cytokine, the groups are diabetes
    # status, and the definition's "high levels of blood sugar" is not the marker.
    "bsdb:36687728-2": (
        "Cytokine",
        "Non-Type 2 Diabetes (Prediabetes)",
        "Type 2 Diabetes",
        "Participants with chronic metabolic condition characterized by high levels "
        "of blood sugar (glucose), insulin resistance, and relative insulin "
        "deficiency, after adjusting for age (used as a covariate).",
    ),
    # Must not be admitted: prose describing a disease, not a cut on inflammation.
    "bsdb:39787116-2": (
        "Chronic obstructive pulmonary disease",
        "Healthy Control Group (HC)",
        "AECOPD Group (Acute Exacerbation of Chronic Obstructive Pulmonary Disease)",
        "Patients diagnosed with Acute Exacerbation of Chronic Obstructive Pulmonary "
        "Disease (AECOPD). It is a sudden worsening of COPD symptoms, including "
        "increased breathlessness, cough, mucus production, and inflammation.",
    ),
    # Must not be admitted: a steroid ester, not a measured short-chain fatty acid.
    "bsdb:30866773-2": (
        "Medroxyprogesterone acetate",
        "black naive vaginal microbiome (baseline visit )",
        "white naiive vaginal microbiome (baseline)",
        "white participants at baseline",
    ),
}


@pytest.mark.parametrize(
    ("experiment", "classification", "evidence_type"),
    [
        ("bsdb:35888992-6", MARKER_STRATIFIED, "inflammation_marker"),
        ("bsdb:38890302-6", MARKER_STRATIFIED, "inflammation_marker"),
        ("bsdb:33220709-1", MARKER_STRATIFIED, "inflammation_marker"),
        ("bsdb:39103909-4", MARKER_STRATIFIED, "inflammation_marker"),
        ("bsdb:39103909-13", MARKER_STRATIFIED, "inflammation_marker"),
        ("bsdb:31065565-2", MARKER_STRATIFIED, "scfa_metabolite"),
        ("bsdb:879-4", INTERVENTION_CONTRAST, "nitrate_oxygen_redox"),
        ("bsdb:39748236-1", INTERVENTION_CONTRAST, "nitrate_oxygen_redox"),
        ("bsdb:39748236-2", INTERVENTION_CONTRAST, "nitrate_oxygen_redox"),
        ("bsdb:37180439-9", MARKER_NAMED_NOT_MEASURED, "inflammation_marker"),
        ("bsdb:36687728-2", MARKER_NAMED_NOT_MEASURED, "inflammation_marker"),
        ("bsdb:39787116-2", MARKER_NAMED_NOT_MEASURED, "inflammation_marker"),
    ],
)
def test_marker_split_is_decided_by_the_group_fields_not_the_condition(
    experiment: str, classification: str, evidence_type: str
) -> None:
    split = classify_marker_split(*SOURCE_ROWS[experiment])
    assert split is not None, experiment
    assert split.classification == classification, experiment
    assert split.evidence_type == evidence_type, experiment
    assert (split.exclusion_reason == "") is (classification == MARKER_STRATIFIED)


def test_calprotectin_split_reports_the_group_field_that_decided_it() -> None:
    split = classify_marker_split(*SOURCE_ROWS["bsdb:35888992-6"])
    assert split is not None
    assert split.split_evidence_field == "Group 0 name"
    assert split.split_evidence_term == "levels"
    assert "calprotectin" in split.terms


def test_an_absent_marker_in_one_group_is_a_measured_split() -> None:
    split = classify_marker_split(*SOURCE_ROWS["bsdb:38890302-6"])
    assert split is not None
    assert split.split_evidence_term == "marker named as absent"


def test_intervention_reason_names_the_term_that_disqualified_the_contrast() -> None:
    therapy = classify_marker_split(*SOURCE_ROWS["bsdb:879-4"])
    probiotic = classify_marker_split(*SOURCE_ROWS["bsdb:39748236-2"])
    assert therapy is not None
    assert probiotic is not None
    assert "therapy" in therapy.exclusion_reason
    assert "probiotic" in probiotic.exclusion_reason


def test_a_condition_naming_no_marker_at_all_is_not_a_candidate() -> None:
    assert classify_marker_split("Colorectal cancer", "Healthy controls", "Patients", "") is None
    # A probiotic species carries an SCFA name only inside its epithet.
    assert classify_marker_split("Clostridium butyricum", "control", "treated", "") is None


def test_marker_candidates_keep_the_excluded_rows_with_their_reason(tmp_path: Path) -> None:
    path = _dump(
        tmp_path / "dump.csv",
        [
            _row(**_source("bsdb:1/1/1", "bsdb:35888992-6")),
            _row(**_source("bsdb:2/1/1", "bsdb:879-4")),
            _row(**{"BSDB ID": "bsdb:3/1/1", "Condition": "Colorectal cancer"}),
        ],
    )
    candidates = marker_split_candidates(path, {("bsdb:1", "bsdb:1-1")})
    found = {str(row["experiment_id"]): row for row in candidates}
    assert sorted(found) == ["bsdb:1-1", "bsdb:2-1"]
    assert found["bsdb:1-1"]["classification"] == MARKER_STRATIFIED
    assert found["bsdb:1-1"]["in_discovery_corpus"] is True
    assert "not independent" in str(found["bsdb:1-1"]["independence"])
    assert found["bsdb:2-1"]["classification"] == INTERVENTION_CONTRAST
    assert found["bsdb:2-1"]["in_discovery_corpus"] is False


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
