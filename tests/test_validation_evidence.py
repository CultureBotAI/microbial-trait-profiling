"""Inventory of measured evidence for independent validation (#191)."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import fields
from pathlib import Path

import pytest

from bsdbng import validation_evidence
from bsdbng.validation_evidence import (
    EVIDENCE_EXCLUSIONS,
    EVIDENCE_TYPES,
    INFERRED_FUNCTION_FIELDS,
    INTERVENTION_CONTRAST,
    MARKER_NAMED_NOT_MEASURED,
    MARKER_STRATIFIED,
    NODES_MISSING,
    NODES_VERIFIED,
    SUN2024_SUPPLEMENT_PDF_SHA256,
    ExternalCohort,
    Mention,
    StudyRecord,
    _write,
    classify_marker_split,
    evidence_rows,
    external_overlap,
    find_mentions,
    inferred_function_traits,
    main,
    marker_split_candidates,
    nodes_status,
    parse_supplement,
    read_studies,
    summarize_sun2024,
    summary_rows,
    title_agreement,
    write_external_snapshot,
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


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


# Curated text exactly as the pinned dump carries it, for the #292 SCFA examples.
SCFA_SOURCE_TEXT: dict[str, str] = {
    # A steroid ester, the condition of all four experiments of bsdb:30866773.
    "bsdb:30866773": "Medroxyprogesterone acetate",
    # A corticosteroid ester, in the title of bsdb:27543211.
    "bsdb:27543211": (
        "Cutaneous microbiome effects of fluticasone propionate cream and adjunctive "
        "bleach baths in childhood atopic dermatitis"
    ),
    # An enzyme name, in the figure description of bsdb:24296500. The source spells it
    # "CoAtransferase", with no hyphen.
    "bsdb:24296500": (
        "Asterisks indicate a significant difference (q < 0.05) compared to baseline. "
        "BCoAT: Butyryl-CoA:acetate CoAtransferase, NND: New Nordic Diet, ADD: "
        "Average Danish"
    ),
    # A measured metabolite, for contrast: bsdb:31065565's group name.
    "bsdb:31065565": "low fecal butyrate levels",
    # Untargeted metabolomics, which is its own evidence type.
    "bsdb:36133313": "gut microbiome, microbiome, mycobiome, plasma metabolites, metabolomics",
}


def _text_record(text: str) -> StudyRecord:
    return StudyRecord(
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
        texts={"Description": [text]},
    )


@pytest.mark.parametrize(
    ("study", "scfa", "excluded", "gene", "metabolomics"),
    [
        ("bsdb:30866773", False, "acetate", False, False),
        ("bsdb:27543211", False, "propionate", False, False),
        ("bsdb:24296500", False, "acetate", True, False),
        ("bsdb:31065565", True, "", False, False),
        ("bsdb:36133313", False, "", False, True),
    ],
)
def test_scfa_categories_separate_drugs_genes_and_metabolomics(
    study: str, scfa: bool, excluded: str, gene: bool, metabolomics: bool
) -> None:
    record = _text_record(SCFA_SOURCE_TEXT[study])

    def mentioned(name: str) -> Mention:
        return find_mentions(record, EVIDENCE_TYPES[name], EVIDENCE_EXCLUSIONS.get(name))

    metabolite = mentioned("scfa_metabolite")
    assert (metabolite.status == "mentioned_in_curated_metadata") is scfa, study
    assert metabolite.excluded_terms == excluded, study
    assert (mentioned("scfa_gene_or_enzyme").status == "mentioned_in_curated_metadata") is gene, (
        study
    )
    assert (
        mentioned("untargeted_metabolomics").status == "mentioned_in_curated_metadata"
    ) is metabolomics, study


def test_the_enzyme_name_is_reported_as_an_enzyme_not_as_acetate() -> None:
    record = _text_record(SCFA_SOURCE_TEXT["bsdb:24296500"])
    enzyme = find_mentions(record, EVIDENCE_TYPES["scfa_gene_or_enzyme"])
    assert "bcoat" in enzyme.terms
    assert "butyryl-coa:acetate" in enzyme.terms
    # The second half of the same source name has to be recognised on its own, so the
    # enzyme rule does not rest on the butyryl-CoA prefix. The source spells it with no
    # hyphen before "transferase", which the published enzyme name carries.
    bare = _text_record("acetate CoAtransferase")
    assert find_mentions(bare, EVIDENCE_TYPES["scfa_gene_or_enzyme"]).terms == (
        "acetate coatransferase"
    )
    assert (
        find_mentions(
            bare, EVIDENCE_TYPES["scfa_metabolite"], EVIDENCE_EXCLUSIONS["scfa_metabolite"]
        ).excluded_terms
        == "acetate"
    )


def test_an_ester_does_not_make_a_marker_split_candidate() -> None:
    # All four experiments of bsdb:30866773 are conditioned on the ester.
    assert (
        classify_marker_split(
            "Medroxyprogesterone acetate",
            "black naive vaginal microbiome (baseline visit )",
            "white naiive vaginal microbiome (baseline)",
            "white participants at baseline",
        )
        is None
    )


def test_a_study_with_no_searchable_text_is_not_reported_as_not_mentioned() -> None:
    empty = StudyRecord(
        study_id="bsdb:1",
        pmid="1",
        doi="",
        title="",
        authors="",
        year="",
        location="",
        host_species="Homo sapiens",
        sequencing_types="",
        antibiotics_exclusion_field="",
        texts={},
    )
    assert find_mentions(empty, EVIDENCE_TYPES["scfa_metabolite"]).status == (
        "no_curated_text_to_search"
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


# The two profile columns of data/derived/study_trait_profiles_wide.tsv that the
# "consumption, not production" reading rests on, with their KG-Microbe labels.
PROFILE_HEADER = (
    "study_id\texperiment_id\tinc_dec\tMETPO:2000006|CHEBI:17272\tCHEBI:16865"
    "\tCHEBI:18401\tMETPO:2000006|mediadive.medium:1\n"
)
NODE_TABLE = (
    "id\tcategory\tname\n"
    "METPO:2000006\tx\tuses as carbon source\n"
    "CHEBI:17272\tx\tpropionate\n"
    "CHEBI:16865\tx\tgamma-aminobutyric acid\n"
    "CHEBI:18401\tx\tphenylacetate\n"
)


def _kg(tmp_path: Path, *, md5: str | None = None, sidecar: bool = True) -> Path:
    """A node table beside a provenance sidecar, as `just download-kg` leaves them."""
    directory = tmp_path / "kg-microbe"
    directory.mkdir(parents=True, exist_ok=True)
    nodes = directory / "merged-kg_nodes.tsv"
    nodes.write_text(NODE_TABLE, encoding="utf-8")
    if sidecar:
        digest = md5 or hashlib.md5(nodes.read_bytes(), usedforsecurity=False).hexdigest()
        (directory / "provenance.json").write_text(
            json.dumps({"files": {"merged-kg_nodes.tsv": {"md5": digest}}}),
            encoding="utf-8",
        )
    return nodes


def test_function_traits_report_what_the_annotation_asserts(tmp_path: Path) -> None:
    profiles = tmp_path / "profiles.tsv"
    profiles.write_text(PROFILE_HEADER, encoding="utf-8")
    nodes = _kg(tmp_path)
    assert nodes_status(nodes) == NODES_VERIFIED
    rows = inferred_function_traits(profiles, nodes)
    assert [row["term"] for row in rows] == ["CHEBI:17272"]
    assert rows[0]["qualifier_label"] == "uses as carbon source"
    assert list(rows[0]) == list(INFERRED_FUNCTION_FIELDS)


def test_nodes_status_separates_missing_mismatched_and_unverifiable(tmp_path: Path) -> None:
    assert nodes_status(tmp_path / "merged-kg_nodes.tsv") == NODES_MISSING
    mismatched = nodes_status(_kg(tmp_path, md5="0" * 32))
    assert mismatched.startswith("present but MD5 ")
    assert "0" * 32 in mismatched
    unverifiable = nodes_status(_kg(tmp_path / "bare", sidecar=False))
    assert unverifiable == "present but unverifiable: no provenance.json beside it"
    no_entry = tmp_path / "empty"
    no_entry.mkdir()
    (no_entry / "provenance.json").write_text('{"files": {}}', encoding="utf-8")
    (no_entry / "merged-kg_nodes.tsv").write_text(NODE_TABLE, encoding="utf-8")
    assert nodes_status(no_entry / "merged-kg_nodes.tsv") == (
        "present but unverifiable: provenance.json records no md5 for it"
    )


@pytest.mark.parametrize("status", [NODES_MISSING, "present but MD5 abc does not match"])
def test_unverified_nodes_refuse_to_produce_labels(tmp_path: Path, status: str) -> None:
    profiles = tmp_path / "profiles.tsv"
    profiles.write_text(PROFILE_HEADER, encoding="utf-8")
    with pytest.raises(ValueError, match="cannot be interpreted"):
        inferred_function_traits(profiles, _kg(tmp_path), status)


def test_an_empty_table_still_gets_its_header(tmp_path: Path) -> None:
    written = _write(tmp_path / "inferred_function_traits.tsv", [], INFERRED_FUNCTION_FIELDS)
    assert written.read_text(encoding="utf-8") == "\t".join(INFERRED_FUNCTION_FIELDS) + "\n"


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


# The pinned sample.group's header, and two rows built from its real column values.
SAMPLE_GROUP = (
    "Sample\tGroup\tGroup_1\tGroup_2\tLevelB\tLevelA\tgroup\tProject\tProject_1"
    "\tProject_2\tLevelA2\tLevelA3\tLevelA4\tLevelA5\tAge\tGender\tBMI\tLocation_1"
    "\tLocation\tsequence_type\tDNA_extract_type\n"
    "SRR8845307\tDisease\tDisease\tDisease\tBL\tmetabolic disease\tlow\tWangQ_2021.BL"
    "\tWangQ_2021\tDiscover\tcardiometabolic disease\tcardiometabolic disease\tBL"
    "\tcardiometabolic disease\t55.3\tFemale\t15.81685744\tShenzhen"
    "\tShenzhen (Guangdong)\tBGISEQ-500\tQinJ_2012\n"
    "SRR8845308\tControl\tControl\tControl\tBL\tmetabolic disease\tlow\tWangQ_2021.BL"
    "\tWangQ_2021\tDiscover\tcardiometabolic disease\tcardiometabolic disease\tBL"
    "\tcardiometabolic disease\tNA\tMale\tNA\tShenzhen"
    "\tShenzhen (Guangdong)\tBGISEQ-500\tQinJ_2012\n"
)


def test_summarize_joins_the_sample_counts_to_the_supplement_description(
    tmp_path: Path,
) -> None:
    sample_group = tmp_path / "sample.group"
    sample_group.write_text(SAMPLE_GROUP, encoding="utf-8")
    supplement = tmp_path / "supplement.txt"
    supplement.write_text(
        "Study ID      Reference\nWangQ_2021    Gut microbiota and metabolic disease\n",
        encoding="utf-8",
    )
    (cohort,) = summarize_sun2024(sample_group, supplement)
    assert cohort.study_id == "WangQ_2021"
    assert cohort.year == "2021"
    assert cohort.samples == 2
    assert cohort.disease_samples == 1
    assert cohort.control_samples == 1
    assert cohort.role == "Discover"
    assert cohort.projects == "WangQ_2021.BL"
    # Age and BMI are counted per sample, and "NA" is absent, not a value.
    assert cohort.samples_with_age == 1
    assert cohort.samples_with_bmi == 1
    assert cohort.reference_title == "Gut microbiota and metabolic disease"


def test_the_snapshot_refuses_an_unpinned_sample_table(tmp_path: Path) -> None:
    sample_group = tmp_path / "sample.group"
    sample_group.write_text(SAMPLE_GROUP, encoding="utf-8")
    supplement = tmp_path / "supplement.txt"
    supplement.write_text("Study ID      Reference\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not the pinned table"):
        write_external_snapshot(sample_group, supplement, tmp_path / "out.tsv")


def test_the_snapshot_refuses_a_supplement_text_it_cannot_identify(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sample_group = tmp_path / "sample.group"
    sample_group.write_text(SAMPLE_GROUP, encoding="utf-8")
    monkeypatch.setattr(
        validation_evidence,
        "SUN2024_SAMPLE_GROUP_SHA256",
        hashlib.sha256(sample_group.read_bytes()).hexdigest(),
    )
    supplement = tmp_path / "supplement.txt"
    supplement.write_text("Study ID      Reference\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not the pinned supplement text") as raised:
        write_external_snapshot(sample_group, supplement, tmp_path / "out.tsv")
    # The message has to say how to produce the right bytes, not just that these are wrong.
    assert "pdftotext -layout" in str(raised.value)
    assert SUN2024_SUPPLEMENT_PDF_SHA256 in str(raised.value)


def test_no_title_match_is_reported_as_such_and_not_as_independence(tmp_path: Path) -> None:
    path = _dump(tmp_path / "dump.csv", [_row(Title="Something else entirely")])
    (row,) = external_overlap([_cohort()], read_studies(path), {"bsdb:1"})
    assert row["overlap_status"] == "no_title_match_found"
    assert row["matched_bugsigdb_studies"] == ""
    assert row["match_basis"] == ""
    assert row["best_title_agreement"] == 0.0


def test_a_newline_in_curated_text_does_not_split_the_row(tmp_path: Path) -> None:
    # bsdb:27362264's Antibiotics exclusion field, exactly as the pinned dump carries it.
    path = tmp_path / "dump.csv"
    path.write_text(
        "# BugSigDB test\n"
        + HEADER
        + "\n"
        + ",".join(
            '"' + _row(**{"Antibiotics exclusion": "X"}).get(name, "NA") + '"'
            for name in HEADER.split(",")
        ).replace('"X"', '"2 weeks prior to delivery\n2 weeks after delivery"')
        + "\n",
        encoding="utf-8",
    )
    studies = read_studies(path)
    assert "\n" in studies["bsdb:1"].antibiotics_exclusion_field
    written = _write(tmp_path / "out.tsv", evidence_rows(studies, {"bsdb:1"}))
    body = written.read_text(encoding="utf-8")
    assert body.count("\n") == 2, "one header line and one row line"
    assert "2 weeks prior to delivery 2 weeks after delivery" in body


def test_the_searched_fields_are_reported_beside_the_absence_of_a_mention(
    tmp_path: Path,
) -> None:
    path = _dump(tmp_path / "dump.csv", [_row()])
    (row,) = evidence_rows(read_studies(path), set())
    assert row["transit_status"] == "not_mentioned_in_curated_metadata"
    # The Group 1 definition and Keywords of this row are NA, so they were not searched.
    assert row["searched_fields"] == "Condition;Group 0 name;Group 1 name;Title"
    assert "access" not in row


def test_expression_and_flux_are_inventoried(tmp_path: Path) -> None:
    path = _dump(
        tmp_path / "dump.csv",
        # bsdb:32143621's and bsdb:33220709's keywords, as the dump carries them.
        [
            _row(Keywords="metatranscriptomics"),
            _row(**{"BSDB ID": "bsdb:2/1/1", "Keywords": "metaproteomics"}),
            _row(**{"BSDB ID": "bsdb:3/1/1"}),
        ],
    )
    rows = evidence_rows(read_studies(path), {"bsdb:1", "bsdb:2", "bsdb:3"})
    statuses = {str(row["study_id"]): row["expression_or_flux_status"] for row in rows}
    assert statuses["bsdb:1"] == "mentioned_in_curated_metadata"
    assert statuses["bsdb:2"] == "mentioned_in_curated_metadata"
    assert statuses["bsdb:3"] == "not_mentioned_in_curated_metadata"
    summary = {row["evidence"]: row for row in summary_rows(rows)}
    assert summary["expression_or_flux"]["discovery_studies_mentioning"] == 2


TABLE1_COLUMNS = (
    "METPO:1000608\tGO:0030245\tMETPO:1000668\tMETPO:1000683"
    "\tpathways:arsenate_reduction\tMETPO:1000623\tMETPO:2000006|CHEBI:17272"
)


def _run_inputs(tmp_path: Path) -> dict[str, Path]:
    """The four inputs main() reads, each a minimal copy of the real thing."""
    profile = tmp_path / "study_trait_profiles_wide.tsv"
    profile.write_text(
        "study_id\texperiment_id\tinc_dec\t" + TABLE1_COLUMNS + "\n"
        "bsdb:1\tbsdb:1-1\tinc\t1\t0\t0\t0\t0\t0\t2\n"
        "bsdb:1\tbsdb:1-1\tdec\t0\t1\t0\t0\t0\t0\t1\n",
        encoding="utf-8",
    )
    experiments = tmp_path / "experiments.tsv"
    experiments.write_text(
        "study_id\texperiment_id\texperiment_name\tcondition\tbody_site\thost_species"
        "\tlocation_of_subjects\tgroup_0_name\tgroup_1_name\tsignature_count\ttaxon_count\n"
        "bsdb:1\tbsdb:1-1\tExperiment 1\tIncreased inflammatory response\tFeces"
        "\tHomo sapiens\tSpain\tlow calprotectin\thigh calprotectin\t2\t3\n",
        encoding="utf-8",
    )
    dump = _dump(tmp_path / "dump.csv", [_row(**_source("bsdb:1/1/1", "bsdb:35888992-6"))])
    cohorts = tmp_path / "cohorts.tsv"
    cohorts.write_text(
        "\t".join(f.name for f in fields(ExternalCohort))
        + "\n"
        + "\t".join(str(getattr(_cohort(), f.name)) for f in fields(ExternalCohort))
        + "\n",
        encoding="utf-8",
    )
    return {
        "--profile-input": profile,
        "--experiments-input": experiments,
        "--full-dump-input": dump,
        "--external-cohorts": cohorts,
    }


def _argv(tmp_path: Path, *extra: str) -> list[str]:
    out = tmp_path / "out"
    args = [part for name, path in _run_inputs(tmp_path).items() for part in (name, str(path))]
    return [*args, "--output-dir", str(out), *extra]


def test_main_writes_every_table_and_a_provenance_record(tmp_path: Path) -> None:
    nodes = _kg(tmp_path)
    assert main([*_argv(tmp_path), "--nodes-input", str(nodes)]) == 0
    out = tmp_path / "out"
    written = sorted(path.name for path in out.iterdir())
    assert written == [
        "evidence_availability_by_study.tsv",
        "evidence_availability_summary.tsv",
        "external_cohort_overlap.tsv",
        "inferred_function_traits.tsv",
        "marker_split_candidates.tsv",
        "validation_evidence_provenance.json",
    ]
    provenance = json.loads(
        (out / "validation_evidence_provenance.json").read_text(encoding="utf-8")
    )
    assert provenance["nodes_status"] == NODES_VERIFIED
    assert provenance["inferred_function_traits_interpretable"] is True
    assert provenance["marker_split_candidates"] == {
        MARKER_STRATIFIED: {"experiments": 1, "human": 1, "in_discovery_corpus": 1}
    }
    assert provenance["external_supplement_extraction"].startswith("pdftotext -layout")
    # The calprotectin experiment is the one measured split in these inputs.
    (candidate,) = _read_rows(out / "marker_split_candidates.tsv")
    assert candidate["experiment_id"] == "bsdb:1-1"
    assert candidate["classification"] == MARKER_STRATIFIED
    assert "no per-study access value" in " ".join(provenance["notes"])


def test_main_warns_and_still_writes_a_header_without_the_node_table(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    missing = tmp_path / "nowhere" / "merged-kg_nodes.tsv"
    assert main([*_argv(tmp_path), "--nodes-input", str(missing)]) == 0
    captured = capsys.readouterr()
    assert "is missing" in captured.err
    assert "just download-kg" in captured.err
    table = (tmp_path / "out" / "inferred_function_traits.tsv").read_text(encoding="utf-8")
    assert table == "\t".join(INFERRED_FUNCTION_FIELDS) + "\n"
    provenance = json.loads(
        (tmp_path / "out" / "validation_evidence_provenance.json").read_text(encoding="utf-8")
    )
    assert provenance["nodes_status"] == NODES_MISSING
    assert provenance["inferred_function_traits_interpretable"] is False


def test_main_fails_on_a_missing_node_table_when_asked_to(tmp_path: Path) -> None:
    missing = tmp_path / "nowhere" / "merged-kg_nodes.tsv"
    with pytest.raises(SystemExit) as raised:
        main([*_argv(tmp_path), "--nodes-input", str(missing), "--require-nodes"])
    assert raised.value.code == 2


def test_report_renders_the_methods_tables_from_the_artifacts(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    nodes = _kg(tmp_path)
    assert main([*_argv(tmp_path, "--report"), "--nodes-input", str(nodes)]) == 0
    report = capsys.readouterr().out
    assert "| Evidence | Discovery studies mentioning |" in report
    assert f"| `{MARKER_STRATIFIED}` | 1 | 1 | 1 |" in report
    assert "| Inflammation markers | 1 | 1 |" in report
    # A pipe inside a profile column has to be escaped or it would split the cell.
    assert r"`METPO:2000006\|CHEBI:17272`" in report
    assert "Node table: verified_against_provenance." in report


def test_supplement_parsing_reads_the_disease_and_accession_sections() -> None:
    # Lines 46-53 and the accession table of the pinned supplement text, verbatim.
    text = (
        " Study ID      Disease name                             Disease type"
        "              Abbr.         Number of     # case   # control   Total"
        "   DNA extraction method\n"
        "                                                                      "
        "                          comparisons\n"
        " ZhuJ_2018     breast cancer                            cancer"
        "                    BC            1             62       71          133"
        "     QIAamp DNA Stool Mini Kit\n"
        " WangQ_2021    bone mass loss                           cardiometabolic disease"
        "   BL            1             231      119         350     QinJ_2012\n"
        "Study ID      Sequencing platform     Antibiotic usage   Probiotics"
        "   Medication   Region       Accession     Journal                    Year\n"
        "ZhuJ_2018     Illumina HiSeq X Ten    N                  no data"
        "      no data      Guangxi      PRJNA453965   Microbiome                 2018\n"
        "WangQ_2021    BGISEQ-500              N                  N"
        "            no data      Guangdong    PRJNA530339   GigaByte                   2021\n"
    )
    parsed = parse_supplement(text)
    assert parsed["ZhuJ_2018"]["disease_name"] == "breast cancer"
    assert parsed["ZhuJ_2018"]["abbreviation"] == "BC"
    assert parsed["ZhuJ_2018"]["accession"] == "PRJNA453965"
    assert parsed["ZhuJ_2018"]["region"] == "Guangxi"
    assert parsed["ZhuJ_2018"]["antibiotic_usage"] == "N"
    assert parsed["ZhuJ_2018"]["probiotics_usage"] == "no data"
    assert parsed["ZhuJ_2018"]["medication_usage"] == "no data"
    assert parsed["WangQ_2021"]["disease_name"] == "bone mass loss"
    assert parsed["WangQ_2021"]["accession"] == "PRJNA530339"
    assert parsed["WangQ_2021"]["probiotics_usage"] == "N"


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
