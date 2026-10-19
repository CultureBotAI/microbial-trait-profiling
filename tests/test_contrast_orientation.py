"""Contrast classes and orientation for BugSigDB experiments (#188)."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
from dataclasses import replace
from pathlib import Path

import pytest

from bsdbng.contrast_orientation import (
    ARM_GROUP,
    CASE_VS_CONTROL_SUBJECT,
    CLASSIFIER_VERSION,
    COMORBIDITY_WITHIN_DISEASE,
    DISEASE_VS_DISEASE,
    DISEASE_VS_HEALTHY,
    INTERVENTION_EVIDENCE,
    LONGITUDINAL,
    MARKER,
    NATIVE,
    NOT_APPLICABLE,
    NOT_REVIEWED,
    OTHER,
    PAIRED_TISSUE_AUTOLOGOUS,
    RESOLVED_DISEASE,
    REVERSED,
    REVIEW_FIELDS,
    TIER_1_HEALTHY_NORMAL,
    TIER_2_CONDITION_NEGATIVE_SCREENED,
    TIER_3_HOSPITAL_COMORBID_OR_TISSUE,
    TREATMENT,
    UNRESOLVED,
    Classification,
    ExperimentRecord,
    class_counts,
    classify,
    comparison_rows,
    control_tier_counts,
    main,
    oriented_profiles,
    read_experiment_records,
    read_review,
    review_sample,
    write_manifest,
)
from bsdbng.trait_guilds import TABLE1_TRAITS, Profile

BASE = ExperimentRecord(
    study_id="bsdb:1",
    experiment_id="bsdb:1-1",
    pmid="1",
    doi="",
    study_design="case-control",
    condition="Colorectal cancer",
    efo_id="EFO:1",
    body_site="Feces",
    uberon_id="UBERON:1",
    host_species="Homo sapiens",
    group_0_name="Healthy controls",
    group_1_name="CRC patients",
    group_1_definition="patients with colorectal cancer",
    group_0_sample_size="10",
    group_1_sample_size="10",
    description="",
)


def _classify(**changes: str) -> tuple[str, str]:
    found = classify(replace(BASE, **changes))
    return found.contrast_class, found.orientation


@pytest.mark.parametrize(
    ("changes", "expected"),
    [
        ({}, (DISEASE_VS_HEALTHY, NATIVE)),
        # A reversal needs a plain comparison group and a disease group opposite it.
        (
            {
                "group_0_name": "Women with polycystic ovary syndrome",
                "group_1_name": "Healthy control",
                "condition": "Polycystic ovary syndrome",
            },
            (DISEASE_VS_HEALTHY, REVERSED),
        ),
        # "normal pH" is a level, not a healthy comparison, so nothing is reversed.
        (
            {"group_0_name": "abnormal pH", "group_1_name": "normal pH"},
            (UNRESOLVED, NOT_APPLICABLE),
        ),
        (
            {
                "group_0_name": "Adults",
                "group_1_name": "normal weight children",
                "condition": "Obesity",
            },
            (UNRESOLVED, NOT_APPLICABLE),
        ),
        (
            {
                "group_0_name": "Mild COVID-19 patients",
                "group_1_name": "Severe COVID-19 patients",
                "condition": "COVID-19",
            },
            (DISEASE_VS_DISEASE, NOT_APPLICABLE),
        ),
        ({"group_0_name": "Baseline", "group_1_name": "Follow up"}, (LONGITUDINAL, NOT_APPLICABLE)),
        (
            {"group_0_name": "pre-FMT", "group_1_name": "post-FMT (2 months)"},
            (LONGITUDINAL, NOT_APPLICABLE),
        ),
        (
            {"group_0_name": "Recovered patients", "group_1_name": "Healthy controls"},
            (LONGITUDINAL, NOT_APPLICABLE),
        ),
        ({"group_0_name": "Placebo", "group_1_name": "Probiotic"}, (TREATMENT, NOT_APPLICABLE)),
        ({"group_0_name": "non-smokers", "group_1_name": "smokers"}, (TREATMENT, NOT_APPLICABLE)),
        (
            {
                "group_0_name": "Caucasians",
                "group_1_name": "African Americans",
                "condition": "Ethnic group",
            },
            (OTHER, NOT_APPLICABLE),
        ),
        (
            {
                "group_0_name": "Low level of LDL",
                "group_1_name": "High level of LDL",
                "condition": "LDL",
            },
            (OTHER, NOT_APPLICABLE),
        ),
        ({"group_0_name": "", "group_1_name": "CRC patients"}, (UNRESOLVED, NOT_APPLICABLE)),
    ],
)
def test_contrast_classes(changes: dict[str, str], expected: tuple[str, str]) -> None:
    assert _classify(**changes) == expected


def test_a_disease_whose_name_starts_with_non_is_not_its_own_absence() -> None:
    """The disease "Non-alcoholic fatty liver disease" is not alcohol exposure or its absence."""
    assert _classify(
        group_0_name="healthy controls",
        group_1_name="non-alcoholic fatty liver disease patients",
        condition="Non-alcoholic fatty liver disease",
    ) == (DISEASE_VS_HEALTHY, NATIVE)


def test_absence_of_the_condition_is_a_comparison_but_only_a_weak_one() -> None:
    found = classify(
        replace(
            BASE, group_0_name="Non-IBS", group_1_name="IBS-D", condition="Irritable bowel syndrome"
        )
    )
    assert (found.contrast_class, found.orientation) == (DISEASE_VS_HEALTHY, NATIVE)
    assert "weak_comparison_evidence" in found.flags


def test_test_markers_compare_only_the_marker_that_differs() -> None:
    assert _classify(
        group_0_name="HPV-/HIV-",
        group_1_name="HPV+/HIV-",
        condition="Human papilloma virus infection",
    ) == (DISEASE_VS_HEALTHY, NATIVE)


def test_a_hyphen_used_as_punctuation_is_not_a_negative_test() -> None:
    """ "NAFLD- [fatty liver]" against "NASH" is two liver diseases, not a sign.

    The result is ``unresolved`` rather than ``disease_vs_disease`` because "NASH" is
    not lexically the condition: Group 0 is not read as a comparison, which is what this
    guards, and Group 1 supports no reading at all.
    """
    assert _classify(
        group_0_name="NAFLD- [fatty liver]",
        group_1_name="NASH",
        condition="Non-alcoholic fatty liver disease",
    ) == (UNRESOLVED, NOT_APPLICABLE)


def test_a_marker_sign_is_absence_of_the_condition_only_when_it_is_the_condition() -> None:
    """Her2 status within breast cancer is not breast cancer against health."""
    assert _classify(
        group_0_name="tumor HP-",
        group_1_name="tumor HP+",
        condition="Gastric cancer",
    ) == (DISEASE_VS_DISEASE, NOT_APPLICABLE)


def test_a_two_by_two_marker_design_has_no_unaffected_group() -> None:
    """Both markers flip, so neither group is the other's comparison."""
    assert _classify(
        group_0_name="HIV+ MASLD-",
        group_1_name="HIV- MASLD+",
        condition="Non-alcoholic fatty liver disease",
    ) == (UNRESOLVED, NOT_APPLICABLE)


#: Every contrast that review #291 found wrongly admitted to, or wrongly kept out of,
#: ``disease_vs_healthy``, with the curated strings copied verbatim from the pinned
#: v1.3.1 export. The identifier is PMID and experiment number, as the review cites them.
ELIGIBILITY_FIXTURES: tuple[tuple[str, dict[str, str], tuple[str, str]], ...] = (
    # A bare "negative" whose subject is a marker, an antibody or a drug, not the
    # condition. Everyone in these contrasts has the disease.
    (
        "32468338-1 Her2-/+ within breast carcinoma",
        {
            "condition": "HER2 Positive Breast Carcinoma",
            "group_0_name": "Her2 negative tumor",
            "group_1_name": "Her2 positive tumor",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    (
        "31764438-1 ASCA-/+ within Crohn's disease",
        {
            "condition": "Crohn's disease",
            "group_0_name": "crohn disease ASCA negative patients",
            "group_1_name": "crohn disease ASCA positive patients",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    (
        "29596446-1 a drug whose misspelling defeats the term list",
        {
            "condition": "Type II diabetes mellitus",
            "group_0_name": "Metaformin negative among Opoid negative",
            "group_1_name": "Metaformin positive",
            "description": (
                "the interactive influence of metformin and opioids on gut microbiota in "
                "the subgroup with diabetes."
            ),
        },
        (TREATMENT, NOT_APPLICABLE),
    ),
    (
        "36494405-3 one dementia negated inside another",
        {
            "condition": "Lewy body dementia",
            "group_0_name": "Parkinson\u2019s disease dementia negative (PDD-)",
            "group_1_name": "Dementia with Lewy Bodies (DLB) patients",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    # "non-/no <modifier> <condition>": the negation qualifies the condition instead of
    # denying it, so both groups still have it.
    (
        "29518419-3 non-exacerbated asthma",
        {
            "condition": "Asthma",
            "group_0_name": "non-exacerbated asthma",
            "group_1_name": "exacerbated asthma",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    (
        "34609073-2 non-cachectic cancer",
        {
            "condition": "Cancer",
            "group_0_name": "Non-Cachectic Cancer",
            "group_1_name": "Cachetic Cancer",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    (
        "38006331-2 the negation qualifies the tumour, not the patient",
        {
            "condition": "Intestinal obstruction",
            "group_0_name": "Non-obstruction tumor group",
            "group_1_name": "Obstruction tumor group",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    (
        "38596264-1 no liver metastasis within metastatic colorectal cancer",
        {
            "condition": "Metastatic colorectal cancer",
            "group_0_name": "No liver metastasis (NLM)",
            "group_1_name": "Liver metastasis (LM)",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    (
        "38589368-4 non-inflammatory depression",
        {
            "condition": "Depressive disorder",
            "group_0_name": "Non-inflammatory depression",
            "group_1_name": "Inflammatory depression",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    (
        "40450263-5 non-severe preeclampsia",
        {
            "condition": "Preeclampsia",
            "group_0_name": "Mild/Non-severe Preeclampsia (P)",
            "group_1_name": "Severe Preeclampsia (SP)",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    (
        "26901400-3 non-endometrioma cysts are still cysts",
        {
            "condition": "Endometriosis",
            "group_0_name": "Patients with non-endometrioma cysts",
            "group_1_name": "Patients with ovarian endometrioma cysts",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    (
        "32777190-1 pregnancy is a state, not a disease",
        {
            "condition": "Pregnancy",
            "group_0_name": "non pregnant women without gingivitis",
            "group_1_name": "pregnant women without gingivitis",
        },
        (OTHER, NOT_APPLICABLE),
    ),
    (
        "34631595-1 non-pregnant against clinical pregnancy",
        {
            "condition": "Pregnancy",
            "group_0_name": "Non-pregnant (FN)",
            "group_1_name": "Clinical pregnancy (FP)",
        },
        (OTHER, NOT_APPLICABLE),
    ),
    # "Control" naming an untreated arm rather than a healthy comparison.
    (
        "33505487-1 control against an experimental arm",
        {
            "condition": "Chronic constipation",
            "group_0_name": "Control group",
            "group_1_name": "Experimental group",
            "group_1_definition": (
                "Experimental group of patients before and after AMT (Acupoint Massage Therapy)"
            ),
        },
        (TREATMENT, NOT_APPLICABLE),
    ),
    (
        "36760344-1 control against a meditation arm",
        {
            "condition": "Exercise",
            "group_0_name": "Control group",
            "group_1_name": "Meditation group",
            "group_1_definition": (
                "Tibetan Buddhist monks that performed meditation practices of Samatha and "
                "Vipassana for at least 2 hours a day for 3 to 30 years."
            ),
        },
        (TREATMENT, NOT_APPLICABLE),
    ),
    (
        "31481657-1 controls against a marijuana usage arm",
        {
            "condition": "Cannabis use",
            "group_0_name": "controls (lateral border of tongue)",
            "group_1_name": "marijuana usage group",
        },
        (TREATMENT, NOT_APPLICABLE),
    ),
    (
        "34264976-4 control against a trial arm",
        {
            "condition": "Metabolic process",
            "group_0_name": "Control group",
            "group_1_name": "Trial group; 35d (Intestinal)",
        },
        (TREATMENT, NOT_APPLICABLE),
    ),
    # A healthy-sounding label that names a subgroup of the disease.
    (
        "28457228-2 an IBS subgroup that resembles healthy controls",
        {
            "condition": "Irritable bowel syndrome",
            "group_0_name": "HC-like IBS(Healthy Control like Irritable Bowel Syndrome)",
            "group_1_name": "IBS1 (Irritable Bowel Syndrome 1)",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    (
        "40110195-1 one value of the condition denied, the other kept",
        {
            "condition": "Gastric adenocarcinoma,Epstein-Barr virus infection",
            "group_0_name": "Epstein-Barr virus non-infected (EBV- GC)",
            "group_1_name": "Epstein-Barr virus infected (EBV+ GC)",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    # Age bands read as time points.
    (
        "27474122-7 age bands are not time points",
        {
            "condition": "Age",
            "group_0_name": "2-12 years",
            "group_1_name": "18-62 years",
        },
        (OTHER, NOT_APPLICABLE),
    ),
    (
        "25484919-1 age bands within one control group",
        {
            "condition": "Chlamydia trachomatis",
            "group_0_name": "Controls aged < = 0 years",
            "group_1_name": "Controls aged > 10 years",
        },
        (UNRESOLVED, NOT_APPLICABLE),
    ),
    (
        "29808701-4 a multi-valued condition with only one value denied",
        {
            "condition": "Coinfection,Age",
            "group_0_name": "young HIV-infected group without coinfection",
            "group_1_name": "aging HIV-infected women with coinfections",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    # Contrasts the previous rules excluded too conservatively.
    (
        "39355268-4 a normal group is a comparison group",
        {
            "condition": "Gastric cancer",
            "group_0_name": "Normal group",
            "group_1_name": "Gastric Cancer",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    (
        "31109865-1 noncancer, spelled without a separator and in other words",
        {
            "condition": "Squamous cell carcinoma",
            "group_0_name": "noncancer individuals",
            "group_1_name": "Oral Squamous Cell Cancer patients",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    (
        "24451087-1 a drug negated in the disease group, not the condition",
        {
            "condition": "HIV infection",
            "group_0_name": "Healthy controls",
            "group_1_name": "HIV infected cART negative",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    (
        "34100340-4 the sign comes before its subject",
        {
            "condition": "COVID-19",
            "group_0_name": "Negative for COVID-19 viral RNA in feces",
            "group_1_name": "Positive for COVID-19 viral RNA in feces",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    (
        "30782659-1 negative for a lesion is not negative for the condition",
        {
            "condition": "Human papilloma virus infection",
            "group_0_name": "negative for intraepithelial lesion or malignancy",
            "group_1_name": "high-grade squamous intraepithelial lesions",
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    (
        "25764541-3 a negator chained onto a disease named non-something",
        {
            "condition": "Non-alcoholic fatty liver disease",
            "group_0_name": "obese without non-alcoholic fatty liver disease",
            "group_1_name": "non-alcoholic fatty liver disease patients",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    # Genuine comparisons that the eligibility rules must not take away.
    (
        "30576008-1 a condition named only inside the negation",
        {
            "condition": "Parkinson's disease",
            "group_0_name": "healthy controls without any form of Parkinson's Disease",
            "group_1_name": "patients with idiopathic parkinson's disease",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    (
        "29486796-14 an uninfected group whose name explains the zero count",
        {
            "condition": "Helminthiasis",
            "group_0_name": "Uninfected (Zero STH infection)",
            "group_1_name": "Heavy Infected (Ascaris)",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    (
        "25058515-1 a hyphenated uninfected group",
        {
            "condition": "HIV infection",
            "group_0_name": "HIV-uninfected men",
            "group_1_name": "HIV-infected men",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    (
        "31808612-1 comparison words do not state a condition called Health trait",
        {
            "condition": "Health trait",
            "group_0_name": "healthy controls",
            "group_1_name": "suboptimal health status",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    (
        "39818212-1 a study label in the group name is not a diagnosis",
        {
            "condition": "Chlamydia trachomatis",
            "group_0_name": "Chlamydia trichomonas Controls (t0)",
            "group_1_name": "Chlamydia trichomonas Cases (t0)",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    (
        "38594374-1 a subtype name is not a severity modifier",
        {
            "condition": "Chronic constipation",
            "group_0_name": "Children without Functional constipation",
            "group_1_name": "Children with Functional constipation",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    (
        "37119735-1 prevalent says which case, not how bad",
        {
            "condition": "Atrial fibrillation",
            "group_0_name": "Individuals without Prevalent Atrial Fibrillation (AF)",
            "group_1_name": "Individuals with Prevalent Atrial Fibrillation (AF)",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    (
        "40258140-3 experimental inside a disease name is not an arm",
        {
            "condition": "Experimental autoimmune encephalomyelitis",
            "group_0_name": "Healthy",
            "group_1_name": "EAE (Experimental Autoimmune Encephalomyelitis)",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    (
        "27562258-4 usage at the end of a name is not a usage arm",
        {
            "condition": "Eczema",
            "group_0_name": "Normal control",
            "group_1_name": "AD cases with no reported cream usage",
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
    (
        "33531650-1 a drug in the cohort's history is not the contrast",
        {
            "condition": "Type II diabetes mellitus",
            "group_0_name": "Healthy controls",
            "group_1_name": (
                "Individuals having Type 2 Diabetes Mellitus without Diabetic Retinopathy"
            ),
            "group_1_definition": (
                "Individuals having Type 2 Diabetes Mellitus without Diabetic Retinopathy "
                "(T2D group); the T2DM cohort had history of taking anti-diabetic "
                "medication such as metformin"
            ),
        },
        (UNRESOLVED, NOT_APPLICABLE),
    ),
    (
        "31611994-1 cancer therapy in a description is not a therapy contrast",
        {
            "condition": "Colorectal cancer",
            "group_0_name": "Healthy Controls",
            "group_1_name": "CRC",
            "description": (
                "Patients were sampled before chemotherapy or radiotherapy was started."
            ),
        },
        (DISEASE_VS_HEALTHY, NATIVE),
    ),
)


@pytest.mark.parametrize(
    ("changes", "expected"),
    [pytest.param(changes, expected, id=name) for name, changes, expected in ELIGIBILITY_FIXTURES],
)
def test_disease_versus_healthy_eligibility(
    changes: dict[str, str], expected: tuple[str, str]
) -> None:
    assert _classify(**changes) == expected


def test_cancer_therapy_is_not_read_as_an_intervention_in_free_text() -> None:
    """A case-control description routinely dates its samples against cancer therapy.

    Scanning the Group 1 definition and description for these words moved eight genuine
    healthy-control contrasts out of the population, so they are excluded by name. The
    approximate-match gate in :func:`_stated_intervention` also blocks them, and this
    pins the term list itself so neither guard is the only one.
    """
    for word in ("chemotherapy", "radiotherapy", "immunotherapy", "vaccinated"):
        assert not INTERVENTION_EVIDENCE.search(f"sampled before {word} was started")
    for word in ("placebo", "metformin", "meditation", "acupoint", "massage"):
        assert INTERVENTION_EVIDENCE.search(f"the group received {word} daily")


def test_an_arm_label_needs_the_word_group_to_count_as_an_arm() -> None:
    """Otherwise a disease called "Experimental Autoimmune Encephalomyelitis" is an arm."""
    assert ARM_GROUP.search("Trial group; 35d (Intestinal)")
    assert ARM_GROUP.search("marijuana usage group")
    assert not ARM_GROUP.search("EAE (Experimental Autoimmune Encephalomyelitis)")
    assert not ARM_GROUP.search("Low Trial Making Test (TMT) B Scores")
    assert not ARM_GROUP.search("AD cases with no reported cream usage")


def test_a_shared_timepoint_is_not_a_time_contrast() -> None:
    assert _classify(
        group_0_name="Stable PD patients at baseline",
        group_1_name="Progressed PD patients at baseline",
        condition="Parkinson's disease",
    ) == (DISEASE_VS_DISEASE, NOT_APPLICABLE)


def test_negated_treatment_names_untreated_disease() -> None:
    assert _classify(
        group_0_name="HIV negative patients",
        group_1_name="HIV-positive newly diagnosed without antiretroviral treatment",
        condition="HIV infection",
    ) == (DISEASE_VS_HEALTHY, NATIVE)


def test_missing_definitions_are_flagged_not_filled() -> None:
    found = classify(replace(BASE, group_1_definition="NA"))
    assert "no_group_1_definition" in found.flags
    assert "no_group_0_definition_in_export" in found.flags


def _profile(experiment: str, direction: str) -> Profile:
    return Profile("bsdb:1", experiment, direction, "C", "Feces", {"t": 1}, 2)


def test_reversed_experiments_are_relabelled_and_others_dropped() -> None:
    profiles = [
        _profile(e, d) for e in ("bsdb:1-1", "bsdb:1-2", "bsdb:1-3") for d in ("inc", "dec")
    ]
    classes = {
        ("bsdb:1", "bsdb:1-1"): Classification(DISEASE_VS_HEALTHY, NATIVE, "", ""),
        ("bsdb:1", "bsdb:1-2"): Classification(DISEASE_VS_HEALTHY, REVERSED, "", ""),
        ("bsdb:1", "bsdb:1-3"): Classification(TREATMENT, NOT_APPLICABLE, "", ""),
    }
    kept = oriented_profiles(profiles, classes)
    assert {(p.experiment_id, p.inc_dec) for p in kept} == {
        ("bsdb:1-1", "inc"),
        ("bsdb:1-1", "dec"),
        ("bsdb:1-2", "dec"),
        ("bsdb:1-2", "inc"),
    }
    assert [p.inc_dec for p in kept if p.experiment_id == "bsdb:1-2"] == ["dec", "inc"]


HEADER = (
    "BSDB ID,Study design,PMID,DOI,Condition,EFO ID,Body site,UBERON ID,Host species,"
    "Group 0 name,Group 1 name,Group 1 definition,Group 0 sample size,Group 1 sample size,"
    "Description"
)


def test_export_rows_collapse_to_experiments_and_conflicts_raise(tmp_path: Path) -> None:
    good = tmp_path / "dump.csv"
    good.write_text(
        "# BugSigDB test\n" + HEADER + "\n"
        "bsdb:7/2/1,case-control,1,NA,Asthma,E,Feces,U,Homo sapiens,Controls,Asthma,d,1,1,x\n"
        "bsdb:7/2/2,case-control,1,NA,Asthma,E,Feces,U,Homo sapiens,Controls,Asthma,d,1,1,x\n",
        encoding="utf-8",
    )
    records = read_experiment_records(good)
    assert list(records) == [("bsdb:7", "bsdb:7-2")]
    assert records[("bsdb:7", "bsdb:7-2")].doi == ""
    bad = tmp_path / "bad.csv"
    bad.write_text(
        HEADER + "\n"
        "bsdb:7/2/1,case-control,1,NA,Asthma,E,Feces,U,Homo sapiens,Controls,Asthma,d,1,1,x\n"
        "bsdb:7/2/2,case-control,1,NA,Asthma,E,Feces,U,Homo sapiens,Patients,Asthma,d,1,1,x\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="conflicting"):
        read_experiment_records(bad)


def test_a_curated_line_break_does_not_split_a_manifest_row(tmp_path: Path) -> None:
    """187 of 4,457 curated descriptions hold a line break.

    The csv module quotes them, so the file stays valid CSV, but every row after the
    first such description lands on its own physical line and grep, cut and head all
    read the wrong fields.
    """
    dump = tmp_path / "dump.csv"
    dump.write_text(
        HEADER + "\n"
        "bsdb:7/2/1,case-control,1,NA,Asthma,E,Feces,U,Homo sapiens,Controls,Asthma,"
        '"a definition\nover two lines",1,1,"a description\nover two lines"\n',
        encoding="utf-8",
    )
    records = read_experiment_records(dump)
    record = records[("bsdb:7", "bsdb:7-2")]
    assert record.description == "a description over two lines"
    assert record.group_1_definition == "a definition over two lines"
    manifest = write_manifest(
        tmp_path / "manifest.tsv", records, {k: classify(v) for k, v in records.items()}, [], None
    )
    assert len(manifest.read_text(encoding="utf-8").splitlines()) == 2


def test_the_review_sample_is_reproducible_from_its_seed() -> None:
    records = {
        ("bsdb:1", f"bsdb:1-{n}"): replace(BASE, study_id="bsdb:1", experiment_id=f"bsdb:1-{n}")
        for n in range(1, 21)
    }
    classes = {key: classify(record) for key, record in records.items()}
    corpus = list(records)
    first = review_sample(records, classes, corpus, per_class=5, seed=188)
    assert first == review_sample(records, classes, corpus, per_class=5, seed=188)
    assert first != review_sample(records, classes, corpus, per_class=5, seed=2026)
    # The corpus order must not matter, or the sample depends on how the caller built it.
    assert first == review_sample(records, classes, list(reversed(corpus)), per_class=5, seed=188)
    assert [row["experiment_id"] for row in first] == sorted(row["experiment_id"] for row in first)
    assert len(first) == 5
    assert {row["sampled_class"] for row in first} == {DISEASE_VS_HEALTHY}
    # A class with fewer than per_class experiments contributes all of them.
    assert len(review_sample(records, classes, corpus, per_class=50, seed=188)) == len(records)
    # Each experiment's place depends on its own identifier alone, so a larger sample
    # extends a smaller one instead of replacing it. A per-draw shuffle would not nest.
    nested = [
        {
            row["experiment_id"]
            for row in review_sample(records, classes, corpus, per_class=size, seed=188)
        }
        for size in range(1, 8)
    ]
    for smaller, larger in itertools.pairwise(nested):
        assert smaller < larger


def test_the_committed_review_covers_the_sample_it_claims_and_is_joined_back() -> None:
    """The reported agreement rate is only auditable if these cases are committed."""
    review = Path(__file__).resolve().parents[1] / "methods" / "contrast_orientation_review.tsv"
    rows = list(csv.DictReader(review.open(encoding="utf-8", newline=""), delimiter="\t"))
    assert len(rows) == 60
    assert {row["verdict"] for row in rows} == {"agree", "disagree"}
    assert all(row["reviewer"] for row in rows)
    assert all(row["note"] for row in rows if row["verdict"] == "disagree")
    statuses = read_review(review)
    assert len(statuses) == len(rows)
    assert set(statuses.values()) == {"reviewed_agreed", "reviewed_disagreed"}


def test_an_unreviewed_row_in_the_review_table_is_an_error(tmp_path: Path) -> None:
    """A blank verdict would otherwise be counted as neither agreement nor disagreement."""
    review = tmp_path / "review.tsv"
    review.write_text(
        "\t".join(REVIEW_FIELDS) + "\n" + "\t".join(["bsdb:1", "bsdb:1-1"] + [""] * 11) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="must be reviewed"):
        read_review(review)


def test_a_verdict_reaches_the_manifest_and_an_unsampled_row_stays_unreviewed(
    tmp_path: Path,
) -> None:
    assert read_review(tmp_path / "missing.tsv") == {}
    records = {
        ("bsdb:1", "bsdb:1-1"): replace(BASE, experiment_id="bsdb:1-1"),
        ("bsdb:1", "bsdb:1-2"): replace(BASE, experiment_id="bsdb:1-2"),
    }
    classes = {key: classify(record) for key, record in records.items()}
    manifest = write_manifest(
        tmp_path / "manifest.tsv",
        records,
        classes,
        [],
        None,
        {("bsdb:1", "bsdb:1-1"): "reviewed_disagreed"},
    )
    rows = list(csv.DictReader(manifest.open(encoding="utf-8", newline=""), delimiter="\t"))
    assert [(row["experiment_id"], row["validation_status"]) for row in rows] == [
        ("bsdb:1-1", "reviewed_disagreed"),
        ("bsdb:1-2", NOT_REVIEWED),
    ]


def test_both_groups_reading_as_comparisons_is_unresolved_not_a_guess() -> None:
    found = classify(
        replace(
            BASE,
            group_0_name="Healthy controls",
            group_1_name="Non-IBS volunteers",
            condition="Irritable bowel syndrome",
        )
    )
    assert (found.contrast_class, found.orientation) == (UNRESOLVED, NOT_APPLICABLE)
    assert "both groups read as comparisons" in found.evidence


def test_a_possible_reversal_without_disease_evidence_stays_unresolved() -> None:
    """Group 1 looks like the comparison, but Group 0 shows no disease to reverse toward."""
    found = classify(
        replace(
            BASE,
            group_0_name="Group A",
            group_1_name="Healthy controls",
            condition="Colorectal cancer",
        )
    )
    assert (found.contrast_class, found.orientation) == (UNRESOLVED, NOT_APPLICABLE)
    assert "possible_reversal_unverified" in found.flags.split(",")


def test_a_randomized_trial_without_a_comparison_group_name_is_treatment() -> None:
    found = classify(
        replace(
            BASE,
            study_design="randomized controlled trial",
            group_0_name="Arm A",
            group_1_name="Arm B",
            condition="Colorectal cancer",
        )
    )
    assert (found.contrast_class, found.orientation) == (TREATMENT, NOT_APPLICABLE)
    assert "design_rct" in found.flags.split(",")
    assert "randomized controlled trial" in found.evidence


def test_a_randomized_trial_with_a_healthy_comparison_is_still_disease_versus_healthy() -> None:
    found = classify(replace(BASE, study_design="randomized controlled trial"))
    assert (found.contrast_class, found.orientation) == (DISEASE_VS_HEALTHY, NATIVE)
    assert "design_rct" in found.flags.split(",")


def test_baseline_beside_a_healthy_comparison_only_dates_the_disease_sample() -> None:
    found = classify(
        replace(BASE, group_0_name="Healthy controls", group_1_name="CRC patients at baseline")
    )
    assert (found.contrast_class, found.orientation) == (DISEASE_VS_HEALTHY, NATIVE)
    assert "disease_sample_dated_at_baseline" in found.flags.split(",")


def test_a_variable_condition_beside_a_comparison_group_is_other() -> None:
    found = classify(
        replace(
            BASE,
            group_0_name="Healthy controls",
            group_1_name="Urban residents",
            condition="Population",
        )
    )
    assert (found.contrast_class, found.orientation) == (OTHER, NOT_APPLICABLE)
    assert "the condition is a variable" in found.evidence


def test_one_disease_value_makes_a_multi_valued_condition_more_than_a_variable() -> None:
    """ "Coinfection,Age" is not an age comparison, because coinfection is a disease."""
    assert _classify(
        group_0_name="Younger subjects",
        group_1_name="Older subjects",
        condition="Coinfection,Age",
    ) != (OTHER, NOT_APPLICABLE)
    assert _classify(
        group_0_name="Younger subjects", group_1_name="Older subjects", condition="Age,Sex"
    ) == (OTHER, NOT_APPLICABLE)


def test_a_withdrawn_comparison_reading_is_flagged_not_silently_dropped() -> None:
    found = classify(
        replace(
            BASE,
            group_0_name="HC-like IBS(Healthy Control like Irritable Bowel Syndrome)",
            group_1_name="IBS1 (Irritable Bowel Syndrome 1)",
            condition="Irritable bowel syndrome",
        )
    )
    assert found.contrast_class == DISEASE_VS_DISEASE
    assert "comparison_reading_withdrawn" in found.flags.split(",")


def test_an_age_band_is_flagged_so_the_reading_can_be_checked() -> None:
    found = classify(
        replace(BASE, group_0_name="2-12 years", group_1_name="18-62 years", condition="Age")
    )
    assert found.contrast_class == OTHER
    assert "age_band_not_timepoint" in found.flags.split(",")


def test_the_strong_evidence_population_drops_weakly_evidenced_contrasts() -> None:
    profiles = [_profile(e, d) for e in ("bsdb:1-1", "bsdb:1-2") for d in ("inc", "dec")]
    classes = {
        ("bsdb:1", "bsdb:1-1"): Classification(DISEASE_VS_HEALTHY, NATIVE, "", "a_flag"),
        ("bsdb:1", "bsdb:1-2"): Classification(
            DISEASE_VS_HEALTHY, NATIVE, "", "a_flag,weak_comparison_evidence"
        ),
    }
    assert {p.experiment_id for p in oriented_profiles(profiles, classes)} == {
        "bsdb:1-1",
        "bsdb:1-2",
    }
    assert {
        p.experiment_id for p in oriented_profiles(profiles, classes, strong_evidence_only=True)
    } == {"bsdb:1-1"}


def test_class_counts_count_experiments_and_distinct_studies_in_class_order() -> None:
    classes = {
        ("bsdb:1", "bsdb:1-1"): Classification(DISEASE_VS_HEALTHY, NATIVE, "", ""),
        ("bsdb:1", "bsdb:1-2"): Classification(DISEASE_VS_HEALTHY, NATIVE, "", ""),
        ("bsdb:2", "bsdb:2-1"): Classification(DISEASE_VS_HEALTHY, REVERSED, "", ""),
        ("bsdb:3", "bsdb:3-1"): Classification(TREATMENT, NOT_APPLICABLE, "", ""),
        ("bsdb:4", "bsdb:4-1"): Classification(UNRESOLVED, NOT_APPLICABLE, "", ""),
    }
    # bsdb:4-1 is classified but outside the analysed corpus, so it is not counted.
    corpus = [key for key in classes if key[0] != "bsdb:4"]
    assert class_counts(classes, corpus) == [
        {
            "contrast_class": DISEASE_VS_HEALTHY,
            "orientation": NATIVE,
            "experiments": 2,
            "studies": 1,
        },
        {
            "contrast_class": DISEASE_VS_HEALTHY,
            "orientation": REVERSED,
            "experiments": 1,
            "studies": 1,
        },
        {
            "contrast_class": TREATMENT,
            "orientation": NOT_APPLICABLE,
            "experiments": 1,
            "studies": 1,
        },
    ]


def test_comparison_rows_report_both_analyses_for_each_population() -> None:
    # Two experiments in one study, so a study count cannot pass as an experiment count.
    profiles = [
        Profile("bsdb:1", experiment, d, "C", "Feces", dict.fromkeys(TABLE1_TRAITS, 3), 10)
        for experiment in ("bsdb:1-1", "bsdb:1-2")
        for d in ("inc", "dec")
    ]
    rows = comparison_rows([("one", profiles), ("two", profiles)])
    assert {row["population"] for row in rows} == {"one", "two"}
    assert {row["comparison"] for row in rows} == {"fraction", "log_ratio"}
    for population in ("one", "two"):
        mine = [row for row in rows if row["population"] == population]
        assert sum(row["comparison"] == "fraction" for row in mine) == len(TABLE1_TRAITS)
        assert {row["population_studies"] for row in mine} == {1}
        assert {row["population_experiments"] for row in mine} == {2}
    assert {row["family"] for row in rows} == {
        "Benjamini-Hochberg over 6 fractions",
        "Benjamini-Yekutieli over 15 ratios",
    }


def _write_corpus(tmp_path: Path) -> tuple[Path, Path, Path]:
    """A two-experiment corpus: one disease-versus-healthy, one treatment contrast."""
    dump = tmp_path / "full_dump.csv"
    dump.write_text(
        "# BugSigDB 2026-04-24_00:41_UTC, License: Creative Commons Attribution 4.0\n"
        + HEADER
        + "\n"
        "bsdb:1/1/1,case-control,1,NA,Colorectal cancer,E,Feces,U,Homo sapiens,"
        "Healthy controls,CRC patients,d,10,10,x\n"
        "bsdb:2/1/1,case-control,2,NA,Colorectal cancer,E,Feces,U,Homo sapiens,"
        "Placebo,Probiotic,d,10,10,x\n",
        encoding="utf-8",
    )
    terms = [term for ids in TABLE1_TRAITS.values() for term in ids]
    header = ["study_id", "experiment_id", "inc_dec", *terms]
    rows = [
        "\t".join([study, experiment, direction, *[str(count)] * len(terms)])
        for study, experiment, increased in (("bsdb:1", "bsdb:1-1", 4), ("bsdb:2", "bsdb:2-1", 5))
        for direction, count in (("inc", increased), ("dec", increased + 1))
    ]
    profiles = tmp_path / "profiles.tsv"
    profiles.write_text("\n".join(["\t".join(header), *rows]) + "\n", encoding="utf-8")
    experiments = tmp_path / "experiments.tsv"
    experiments.write_text(
        "\n".join(
            [
                "\t".join(["study_id", "experiment_id", "condition", "body_site", "host_species"]),
                "bsdb:1\tbsdb:1-1\tColorectal cancer\tFeces\tHomo sapiens",
                "bsdb:2\tbsdb:2-1\tColorectal cancer\tFeces\tHomo sapiens",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return dump, profiles, experiments


def test_main_writes_every_artefact_and_accounts_for_every_experiment(tmp_path: Path) -> None:
    dump, profiles, experiments = _write_corpus(tmp_path)
    output = tmp_path / "out"
    assert (
        main(
            [
                "--full-dump-input",
                str(dump),
                "--profile-input",
                str(profiles),
                "--experiments-input",
                str(experiments),
                "--output-dir",
                str(output),
                "--code-commit",
                "abc123",
            ]
        )
        == 0
    )
    manifest = list(
        csv.DictReader((output / "contrast_manifest.tsv").open(encoding="utf-8"), delimiter="\t")
    )
    assert [row["experiment_id"] for row in manifest] == ["bsdb:1-1", "bsdb:2-1"]
    assert [row["contrast_class"] for row in manifest] == [DISEASE_VS_HEALTHY, TREATMENT]
    assert [row["contrast_type"] for row in manifest] == [CASE_VS_CONTROL_SUBJECT, TREATMENT]
    assert [row["control_tier"] for row in manifest] == [TIER_1_HEALTHY_NORMAL, NOT_APPLICABLE]
    assert {row["in_analysed_corpus"] for row in manifest} == {"true"}
    # Spelled out, so that bumping the version is a deliberate edit to this test too:
    # the manifest and provenance are what downstream work keys off.
    assert CLASSIFIER_VERSION == "contrast-orientation-v4"
    assert {row["classifier_version"] for row in manifest} == {CLASSIFIER_VERSION}
    assert {row["bugsigdb_vintage"] for row in manifest} == {"2026-04-24_00:41_UTC"}

    counts = list(
        csv.DictReader(
            (output / "contrast_class_counts.tsv").open(encoding="utf-8"), delimiter="\t"
        )
    )
    # Every analysed experiment appears in exactly one class, and none is lost.
    assert sum(int(row["experiments"]) for row in counts) == len(manifest)

    provenance = json.loads(
        (output / "contrast_orientation_provenance.json").read_text(encoding="utf-8")
    )
    assert provenance["code_commit"] == "abc123"
    assert provenance["classifier_version"] == CLASSIFIER_VERSION
    assert provenance["analysed_experiments"] == len(manifest)
    assert sum(row["experiments"] for row in provenance["class_counts"]) == len(manifest)
    assert provenance["control_tier_counts"] == [
        {
            "contrast_type": CASE_VS_CONTROL_SUBJECT,
            "control_tier": TIER_1_HEALTHY_NORMAL,
            "experiments": 1,
            "studies": 1,
        }
    ]
    assert provenance["full_dump_sha256"] == hashlib.sha256(dump.read_bytes()).hexdigest()

    populations = {
        row["population"]
        for row in csv.DictReader(
            (output / "disease_control_comparison.tsv").open(encoding="utf-8"), delimiter="\t"
        )
    }
    assert populations == {
        "historical_human_paired",
        "disease_vs_healthy_oriented",
        "disease_vs_healthy_strong_evidence",
        "disease_vs_healthy_tier1_healthy_normal",
        "disease_vs_healthy_tier1_2_subject_controls",
    }


def test_main_refuses_a_corpus_the_export_does_not_describe(tmp_path: Path) -> None:
    dump, profiles, experiments = _write_corpus(tmp_path)
    dump.write_text(
        "\n".join(
            line
            for line in dump.read_text(encoding="utf-8").splitlines()
            if not line.startswith("bsdb:2/")
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="no export record"):
        main(
            [
                "--full-dump-input",
                str(dump),
                "--profile-input",
                str(profiles),
                "--experiments-input",
                str(experiments),
                "--output-dir",
                str(tmp_path / "out"),
            ]
        )


def test_v3_contrast_type_and_control_tier_stratification() -> None:
    # Tier 1: explicit healthy/normal subject controls
    tier1 = classify(
        replace(
            BASE,
            condition="Preterm birth",
            group_0_name="Full-term infants",
            group_1_name="Preterm infants",
        )
    )
    assert tier1.contrast_class == DISEASE_VS_HEALTHY
    assert tier1.orientation == NATIVE
    assert tier1.contrast_type == CASE_VS_CONTROL_SUBJECT
    assert tier1.control_tier == TIER_1_HEALTHY_NORMAL

    # Tier 2: condition-negative / screened asymptomatic controls (including self-scoped absence)
    tier2 = classify(
        replace(
            BASE,
            condition="glucose metabolism disease",
            group_0_name="non-diabetic",
            group_1_name="diabetic",
        )
    )
    assert tier2.contrast_class == DISEASE_VS_HEALTHY
    assert tier2.orientation == NATIVE
    assert tier2.contrast_type == CASE_VS_CONTROL_SUBJECT
    assert tier2.control_tier == TIER_2_CONDITION_NEGATIVE_SCREENED

    # Tier 3: paired autologous tissue comparison
    tissue = classify(
        replace(
            BASE,
            condition="Colorectal cancer",
            group_0_name="Tumor tissue",
            group_1_name="Adjacent non-cancerous tissue",
        )
    )
    assert tissue.contrast_class == DISEASE_VS_HEALTHY
    assert tissue.orientation == REVERSED
    assert tissue.contrast_type == PAIRED_TISSUE_AUTOLOGOUS
    assert tissue.control_tier == TIER_3_HOSPITAL_COMORBID_OR_TISSUE
    assert "paired_tissue_comparison" in tissue.flags.split(",")

    # Tier 3: comorbidity within shared background disease
    comorbid = classify(
        replace(
            BASE,
            condition="Constipation",
            group_0_name="Maintenance Hemodialysis without Constipation (MHDNCG)",
            group_1_name="Maintenance Hemodialysis with Constipation (MHDCG)",
        )
    )
    assert comorbid.contrast_class == DISEASE_VS_HEALTHY
    assert comorbid.orientation == NATIVE
    assert comorbid.contrast_type == COMORBIDITY_WITHIN_DISEASE
    assert comorbid.control_tier == TIER_3_HOSPITAL_COMORBID_OR_TISSUE
    assert "comorbidity_control" in comorbid.flags.split(",")

    # Tier 3: treated disease group vs healthy control
    treated_conf = classify(
        replace(
            BASE,
            condition="HIV infection",
            group_0_name="Controls",
            group_1_name="HIV Patients treated with Maraviroc",
        )
    )
    assert treated_conf.contrast_class == DISEASE_VS_HEALTHY
    assert treated_conf.orientation == NATIVE
    assert treated_conf.contrast_type == CASE_VS_CONTROL_SUBJECT
    assert treated_conf.control_tier == TIER_3_HOSPITAL_COMORBID_OR_TISSUE
    assert "treated_disease_group" in treated_conf.flags.split(",")

    # Resolved disease subtype of longitudinal_timepoint
    resolved = classify(
        replace(
            BASE,
            condition="COVID-19",
            group_0_name="Healthy controls",
            group_1_name="Recovered patients",
        )
    )
    assert resolved.contrast_class == LONGITUDINAL
    assert resolved.contrast_type == RESOLVED_DISEASE
    assert resolved.control_tier == NOT_APPLICABLE


def test_v3_review_sample_disagreement_fixes() -> None:
    # Shared treatment term in both groups cancels out when one group is a control
    shared_tx = classify(
        replace(
            BASE,
            condition="Oral squamous cell carcinoma",
            group_0_name="Nonsmoking, nondrinking, and dentate healthy controls",
            group_1_name="Nonsmoking, nondrinking, and dentate pre-treatment OSCC patients",
        )
    )
    assert shared_tx.contrast_class == DISEASE_VS_HEALTHY
    assert shared_tx.orientation == NATIVE
    assert "shared_treatment_term" in shared_tx.flags.split(",")

    # Pre-DM is a clinical stage, not a longitudinal timepoint
    pre_dm = classify(
        replace(
            BASE,
            condition="Diabetes mellitus",
            group_0_name="Prediabetes (Pre-DM)",
            group_1_name="Type 2 diabetes mellitus (T2DM)",
        )
    )
    assert pre_dm.contrast_class == DISEASE_VS_DISEASE

    # Pre-treatment disease baseline vs healthy control (reversed)
    pre_tx_rev = classify(
        replace(
            BASE,
            condition="Colorectal cancer",
            group_0_name="Pre-treatment colorectal cancer patients",
            group_1_name="Healthy controls",
        )
    )
    assert pre_tx_rev.contrast_class == DISEASE_VS_HEALTHY
    assert pre_tx_rev.orientation == REVERSED
    assert "disease_sample_dated_at_baseline" in pre_tx_rev.flags.split(",")

    # Shared "controls" label where only one group has specific condition-negative evidence
    obese_ctrl = classify(
        replace(
            BASE,
            condition="Obesity",
            group_0_name="normal-weight controls",
            group_1_name="obese controls",
        )
    )
    assert obese_ctrl.contrast_class == DISEASE_VS_HEALTHY
    assert obese_ctrl.orientation == NATIVE

    # Both groups are healthy controls differing by demographic cohort or body site
    cohort_diff = classify(
        replace(
            BASE,
            condition="Obesity",
            group_0_name="Lean Ghanaians(GL)",
            group_1_name="Lean Americans (UL)",
        )
    )
    assert cohort_diff.contrast_class == OTHER

    # Calendar years on both groups read as longitudinal
    cal_years = classify(
        replace(
            BASE,
            condition="Helminthiasis",
            group_0_name="2008 - Uninfected",
            group_1_name="2010 - Uninfected",
        )
    )
    assert cal_years.contrast_class == LONGITUDINAL


def test_oriented_profiles_filters_by_control_tiers() -> None:
    p1 = Profile("bsdb:1", "bsdb:1-1", "inc", "Feces", "CRC", {"a": 5}, 10)
    p2 = Profile("bsdb:2", "bsdb:2-1", "inc", "Feces", "CRC", {"a": 5}, 10)
    p3 = Profile("bsdb:3", "bsdb:3-1", "inc", "Feces", "CRC", {"a": 5}, 10)
    classes = {
        ("bsdb:1", "bsdb:1-1"): Classification(
            DISEASE_VS_HEALTHY,
            NATIVE,
            "strong",
            "",
            CASE_VS_CONTROL_SUBJECT,
            TIER_1_HEALTHY_NORMAL,
        ),
        ("bsdb:2", "bsdb:2-1"): Classification(
            DISEASE_VS_HEALTHY,
            NATIVE,
            "weak",
            "weak_comparison_evidence",
            CASE_VS_CONTROL_SUBJECT,
            TIER_2_CONDITION_NEGATIVE_SCREENED,
        ),
        ("bsdb:3", "bsdb:3-1"): Classification(
            DISEASE_VS_HEALTHY,
            REVERSED,
            "tissue",
            "paired_tissue_comparison",
            PAIRED_TISSUE_AUTOLOGOUS,
            TIER_3_HOSPITAL_COMORBID_OR_TISSUE,
        ),
    }
    tier1_only = oriented_profiles([p1, p2, p3], classes, control_tiers={TIER_1_HEALTHY_NORMAL})
    assert [p.experiment_id for p in tier1_only] == ["bsdb:1-1"]

    tier1_2 = oriented_profiles(
        [p1, p2, p3],
        classes,
        control_tiers={TIER_1_HEALTHY_NORMAL, TIER_2_CONDITION_NEGATIVE_SCREENED},
    )
    assert [p.experiment_id for p in tier1_2] == ["bsdb:1-1", "bsdb:2-1"]
    assert control_tier_counts(classes, list(classes)) == [
        {
            "contrast_type": CASE_VS_CONTROL_SUBJECT,
            "control_tier": TIER_1_HEALTHY_NORMAL,
            "experiments": 1,
            "studies": 1,
        },
        {
            "contrast_type": CASE_VS_CONTROL_SUBJECT,
            "control_tier": TIER_2_CONDITION_NEGATIVE_SCREENED,
            "experiments": 1,
            "studies": 1,
        },
        {
            "contrast_type": PAIRED_TISSUE_AUTOLOGOUS,
            "control_tier": TIER_3_HOSPITAL_COMORBID_OR_TISSUE,
            "experiments": 1,
            "studies": 1,
        },
    ]


def _stratified(**changes: str) -> tuple[str, str, str, str]:
    found = classify(replace(BASE, **changes))
    return found.contrast_class, found.orientation, found.contrast_type, found.control_tier


#: HIV-negative comparison groups as the export spells them, and people with HIV on
#: antiretroviral therapy compared against them (#569), with the curated strings copied
#: verbatim from the pinned v1.3.1 export. The identifier is PMID and experiment number.
HIV_STATUS_FIXTURES: tuple[tuple[str, dict[str, str], tuple[str, str, str, str]], ...] = (
    # Spellings of an HIV-negative comparison group.
    (
        "33102244-1 an en dash is a minus sign, before a word as at the end",
        {
            "condition": "HIV infection",
            "group_0_name": "HIV\u2013 individuals",
            "group_1_name": "HIV+ individuals",
        },
        (
            DISEASE_VS_HEALTHY,
            NATIVE,
            CASE_VS_CONTROL_SUBJECT,
            TIER_2_CONDITION_NEGATIVE_SCREENED,
        ),
    ),
    (
        "30782659-3 a hyphen at the end of the name",
        {"condition": "HIV infection", "group_0_name": "HIV-", "group_1_name": "HIV+"},
        (
            DISEASE_VS_HEALTHY,
            NATIVE,
            CASE_VS_CONTROL_SUBJECT,
            TIER_2_CONDITION_NEGATIVE_SCREENED,
        ),
    ),
    (
        "31078141-1 a negative test for HIV-1",
        {
            "condition": "HIV infection",
            "group_0_name": "HIV-1-negative individuals",
            "group_1_name": "Recently HIV-1 infected individuals",
        },
        (
            DISEASE_VS_HEALTHY,
            NATIVE,
            CASE_VS_CONTROL_SUBJECT,
            TIER_2_CONDITION_NEGATIVE_SCREENED,
        ),
    ),
    (
        "39854172-1 people without HIV",
        {
            "condition": "Human immunodeficiency virus",
            "group_0_name": "People Without HIV",
            "group_1_name": "People With HIV",
        },
        (
            DISEASE_VS_HEALTHY,
            NATIVE,
            CASE_VS_CONTROL_SUBJECT,
            TIER_2_CONDITION_NEGATIVE_SCREENED,
        ),
    ),
    (
        "37721699-2 the condition's acronym spelled out",
        {
            "condition": "HIV infection",
            "group_0_name": (
                "Pregnant women \u2265 18 years of age without Human immunodeficiency virus "
                "(HIV) at Postpartum (PPt) visit"
            ),
            "group_1_name": (
                "Pregnant women \u2265 18 years of age with Human immunodeficiency virus "
                "(HIV) at Postpartum (PPt) visit."
            ),
        },
        (
            DISEASE_VS_HEALTHY,
            NATIVE,
            CASE_VS_CONTROL_SUBJECT,
            TIER_2_CONDITION_NEGATIVE_SCREENED,
        ),
    ),
    (
        "24399150-1 uninfected controls",
        {
            "condition": "HIV infection",
            "group_0_name": "uninfected controls",
            "group_1_name": "HIV-infected subjects",
        },
        (DISEASE_VS_HEALTHY, NATIVE, CASE_VS_CONTROL_SUBJECT, TIER_1_HEALTHY_NORMAL),
    ),
    # People with HIV on antiretroviral therapy against an HIV-negative comparison group:
    # disease against health with treated cases, not a treatment contrast.
    (
        "29411528-5 HAART-treated cases against healthy controls",
        {
            "condition": "HIV infection",
            "group_0_name": "healthy control",
            "group_1_name": "HAART-treated HIV+ patients",
        },
        (
            DISEASE_VS_HEALTHY,
            NATIVE,
            CASE_VS_CONTROL_SUBJECT,
            TIER_3_HOSPITAL_COMORBID_OR_TISSUE,
        ),
    ),
    (
        "38383483-3 people living with HIV on ART against people living without it",
        {
            "condition": "HIV infection",
            "group_0_name": "Controls/ People living without HIV (PLWOH)",
            "group_1_name": "People living with HIV on Antiretroviral Therapy (PLWH on ART)",
        },
        (
            DISEASE_VS_HEALTHY,
            NATIVE,
            CASE_VS_CONTROL_SUBJECT,
            TIER_3_HOSPITAL_COMORBID_OR_TISSUE,
        ),
    ),
    (
        "39880958-1 on ART against a bare HIV-",
        {
            "condition": "HIV infection",
            "group_0_name": "HIV-",
            "group_1_name": "People Living With HIV on Antiretroviral therapy (ART+ PLWH)",
        },
        (
            DISEASE_VS_HEALTHY,
            NATIVE,
            CASE_VS_CONTROL_SUBJECT,
            TIER_3_HOSPITAL_COMORBID_OR_TISSUE,
        ),
    ),
    (
        "28118207-1 ART-treated cases against HIV-uninfected men",
        {
            "condition": "HIV infection",
            "group_0_name": "HIV-uninfected MSM",
            "group_1_name": "ART-treated HIV infected MSM",
        },
        (
            DISEASE_VS_HEALTHY,
            NATIVE,
            CASE_VS_CONTROL_SUBJECT,
            TIER_3_HOSPITAL_COMORBID_OR_TISSUE,
        ),
    ),
    # Everyone has HIV, or nobody does: no HIV-negative comparison group.
    (
        "27477587-2 treated against untreated people with HIV",
        {
            "condition": "HIV infection",
            "group_0_name": "HAART-naive HIV-1 patients",
            "group_1_name": "HAART-treated HIV-1 patients",
        },
        (TREATMENT, NOT_APPLICABLE, TREATMENT, NOT_APPLICABLE),
    ),
    (
        "34336719-5 the same people before and on ART",
        {
            "condition": "Human immunodeficiency virus",
            "group_0_name": "Acute HIV-infected participants at baseline (A0)",
            "group_1_name": "Acute HIV-infected participants on Antiretroviral Therapy (A12)",
        },
        (TREATMENT, NOT_APPLICABLE, TREATMENT, NOT_APPLICABLE),
    ),
    (
        "36073815-1 exposure in utero between uninfected infants",
        {
            "condition": "HIV mother to child transmission",
            "group_0_name": "HIV-unexposed infants (HUU) at 6 weeks",
            "group_1_name": "HIV-exposed uninfected infants (HEU) at 6 weeks",
        },
        (TREATMENT, NOT_APPLICABLE, TREATMENT, NOT_APPLICABLE),
    ),
)


@pytest.mark.parametrize(
    ("changes", "expected"),
    [pytest.param(changes, expected, id=name) for name, changes, expected in HIV_STATUS_FIXTURES],
)
def test_hiv_status_and_treated_cases(
    changes: dict[str, str], expected: tuple[str, str, str, str]
) -> None:
    assert _stratified(**changes) == expected


def test_treated_cases_are_flagged_rather_than_made_a_treatment_contrast() -> None:
    """The therapy is recorded on the experiment, and a weak comparison group is enough."""
    for control in ("healthy control", "HIV-"):
        found = classify(
            replace(
                BASE,
                condition="HIV infection",
                group_0_name=control,
                group_1_name="People Living With HIV on Antiretroviral therapy (ART+ PLWH)",
            )
        )
        assert found.contrast_class == DISEASE_VS_HEALTHY
        assert "treated_disease_group" in found.flags.split(",")


def test_a_treated_group_named_by_time_keeps_the_treatment_reading() -> None:
    """29808701-1: "for more than 1 year" reads as a time point, not as a duration.

    Clearing the treatment reading would leave only that time reading, so the experiment
    stays ``treatment_exposure`` instead of turning into a ``longitudinal_timepoint`` it is
    not. It compares HIV-negative and HIV-positive women, and is a known remaining miss.
    """
    assert _classify(
        condition="HIV infection",
        group_0_name="Age-matched HIV-uninfected women",
        group_1_name="HIV-positive women who have been taking ART for more than 1 year",
    ) == (TREATMENT, NOT_APPLICABLE)


def test_a_minus_sign_or_en_dash_before_a_word_is_a_sign_but_a_hyphen_is_not() -> None:
    """A hyphen before a word separates a label as often as it negates a test."""
    assert MARKER.findall("HIV\u2013 individuals") == [("HIV", "\u2013")]
    assert MARKER.findall("HIV\u2212 women") == [("HIV", "\u2212")]
    assert MARKER.findall("HIV\u2013") == [("HIV", "\u2013")]
    assert MARKER.findall("HIV-") == [("HIV", "-")]
    assert not MARKER.findall("Grade C (previously aggressive) Periodontitis- Shallow pocket")
    assert not MARKER.findall("BV- Positive (BV: Bacterial Vaginosis) at Visit 1")
    # An en dash inside a word or a range, or a spaced dash, is not a sign either.
    assert not MARKER.findall("Pre\u2013Collagen-Induced Arthritis (Pre-CIA) mice")
    assert not MARKER.findall("Group 1 (18\u201335 years old)")
    assert not MARKER.findall("HIV - negative")


def test_a_spaced_sign_is_read_without_its_space() -> None:
    """The sign's group holds the sign alone, so "HPV +" reads as positive (#596)."""
    assert MARKER.findall("HIV +") == [("HIV", "+")]
    assert MARKER.findall("HPV +") == [("HPV", "+")]
    assert MARKER.findall("HIV \u2013") == [("HIV", "\u2013")]
    assert MARKER.findall("HIV -, HPV +") == [("HIV", "-"), ("HPV", "+")]
    assert not MARKER.findall("HIV \u2013 individuals")
    assert _classify(
        condition="Human papilloma virus infection",
        group_0_name="healthy control",
        group_1_name="HPV +",
    ) == (DISEASE_VS_HEALTHY, NATIVE)


def test_hiv_status_is_a_comparison_only_where_hiv_is_the_condition() -> None:
    """28328799-1 compares HIV status within a periodontitis study; PLWOH is not a disease."""
    assert _classify(condition="Periodontitis", group_0_name="HIV\u2013", group_1_name="HIV+") == (
        UNRESOLVED,
        NOT_APPLICABLE,
    )
    assert _stratified(condition="HIV infection", group_0_name="PLWOH", group_1_name="PLWH") == (
        DISEASE_VS_HEALTHY,
        NATIVE,
        CASE_VS_CONTROL_SUBJECT,
        TIER_2_CONDITION_NEGATIVE_SCREENED,
    )
    assert _classify(condition="Periodontitis", group_0_name="PLWOH", group_1_name="PLWH") == (
        UNRESOLVED,
        NOT_APPLICABLE,
    )


def test_only_the_condition_s_own_acronym_spelled_out_denies_it() -> None:
    with_hiv = "Pregnant women with Human immunodeficiency virus (HIV)"
    assert _classify(
        condition="HIV infection",
        group_0_name="Pregnant women without Human immunodeficiency virus (HIV)",
        group_1_name=with_hiv,
    ) == (DISEASE_VS_HEALTHY, NATIVE)
    assert _classify(
        condition="HIV infection",
        group_0_name="Pregnant women without Hepatitis B virus (HBV)",
        group_1_name=with_hiv,
    ) == (DISEASE_VS_DISEASE, NOT_APPLICABLE)


def test_tumour_against_adjacent_tissue_can_be_kept_out_of_subject_level_controls() -> None:
    """25699023-2 compares tissue within the same patients, not patients with controls.

    It stays ``disease_vs_healthy``, so a population built from the class alone still holds
    it; its contrast type and control tier are what keep it out of a patients-versus-controls
    population (#569).
    """
    tissue = classify(
        replace(
            BASE,
            condition="Colorectal cancer",
            group_0_name="adjacent non-cancerous tissues",
            group_1_name="cancerous tissues",
        )
    )
    assert (tissue.contrast_class, tissue.orientation) == (DISEASE_VS_HEALTHY, NATIVE)
    assert tissue.contrast_type == PAIRED_TISSUE_AUTOLOGOUS
    assert tissue.control_tier == TIER_3_HOSPITAL_COMORBID_OR_TISSUE
    assert "paired_tissue_comparison" in tissue.flags.split(",")

    profiles = [_profile(e, d) for e in ("bsdb:1-1", "bsdb:1-2") for d in ("inc", "dec")]
    classes = {("bsdb:1", "bsdb:1-1"): tissue, ("bsdb:1", "bsdb:1-2"): classify(BASE)}
    assert {p.experiment_id for p in oriented_profiles(profiles, classes)} == {
        "bsdb:1-1",
        "bsdb:1-2",
    }
    subject_level = oriented_profiles(
        profiles,
        classes,
        control_tiers={TIER_1_HEALTHY_NORMAL, TIER_2_CONDITION_NEGATIVE_SCREENED},
    )
    assert {p.experiment_id for p in subject_level} == {"bsdb:1-2"}
