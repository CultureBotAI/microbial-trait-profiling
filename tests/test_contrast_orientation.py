"""Contrast classes and orientation for BugSigDB experiments (#188)."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import re
from collections import Counter
from collections.abc import Callable
from dataclasses import asdict, fields, make_dataclass, replace
from functools import partial
from pathlib import Path
from typing import Any, cast

import pytest

from bsdbng.contrast_orientation import (
    ARM_GROUP,
    CASE_VS_CONTROL_SUBJECT,
    CLASSES,
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
from mtp.contrast_orientation import _acronym_spelled_out
from mtp.table1_trait_maps import compute_direction_shifts
from mtp.trait_composition import compute_log_ratio_shifts

#: The record that the regression tables and most tests below change a few fields of. The
#: classifier-version guard (#793) classifies the regression tables built on it, so an edit
#: here changes the guard's cases as an edit to those tables would.
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
            {"group_0_name": "COVID-19 patients", "group_1_name": "Recovered COVID-19 patients"},
            (LONGITUDINAL, NOT_APPLICABLE),
        ),
        (
            {"group_0_name": "Recovered patients", "group_1_name": "Healthy controls"},
            (OTHER, NOT_APPLICABLE),
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
#: The classifier-version guard (#793) classifies these cases, so they stay as they are:
#: a new case goes in a test of its own.
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
    assert Counter(row["sampled_class"] for row in rows) == dict.fromkeys(CLASSES, 10)
    assert {row["verdict"] for row in rows} == {"agree"}
    assert all(row["reviewer"] for row in rows)
    statuses = read_review(review)
    assert len(statuses) == len(rows)
    assert set(statuses.values()) == {"reviewed_agreed"}


def test_an_unreviewed_row_in_the_review_table_is_an_error(tmp_path: Path) -> None:
    """A blank verdict or invalid label/consistency in the review table raises ValueError."""
    review = tmp_path / "review.tsv"
    review.write_text(
        "\t".join(REVIEW_FIELDS) + "\n" + "\t".join(["bsdb:1", "bsdb:1-1"] + [""] * 11) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="must be reviewed"):
        read_review(review)

    # Invalid sampled_class raises
    review.write_text(
        "\t".join(REVIEW_FIELDS)
        + "\n"
        + "\t".join(
            [
                "bsdb:1",
                "bsdb:1-1",
                "1",
                "CRC",
                "HC",
                "CRC",
                "invalid_class",
                NATIVE,
                "rev",
                "agree",
                "",
                DISEASE_VS_HEALTHY,
                NATIVE,
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="invalid sampled labels"):
        read_review(review)

    # Disagree without a note raises
    review.write_text(
        "\t".join(REVIEW_FIELDS)
        + "\n"
        + "\t".join(
            [
                "bsdb:1",
                "bsdb:1-1",
                "1",
                "CRC",
                "HC",
                "CRC",
                UNRESOLVED,
                NOT_APPLICABLE,
                "rev",
                "disagree",
                "",
                DISEASE_VS_HEALTHY,
                NATIVE,
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="without a note"):
        read_review(review)

    # Valid disagree row with note passes when classifier disagrees, and raises when it matches
    review.write_text(
        "\t".join(REVIEW_FIELDS)
        + "\n"
        + "\t".join(
            [
                "bsdb:1",
                "bsdb:1-1",
                "1",
                "CRC",
                "HC",
                "CRC",
                UNRESOLVED,
                NOT_APPLICABLE,
                "rev",
                "disagree",
                "reason",
                DISEASE_VS_HEALTHY,
                NATIVE,
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    assert read_review(
        review,
        classes={("bsdb:1", "bsdb:1-1"): Classification(UNRESOLVED, NOT_APPLICABLE, "", "")},
    ) == {("bsdb:1", "bsdb:1-1"): "reviewed_disagreed"}
    with pytest.raises(ValueError, match="classifier matches reviewed labels"):
        read_review(
            review,
            classes={
                ("bsdb:1", "bsdb:1-1"): Classification(DISEASE_VS_HEALTHY, NATIVE, "", ""),
            },
        )


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


@pytest.mark.parametrize(("fraction_count", "ratio_count"), [(3, 4), (4, 3)])
@pytest.mark.parametrize("all_tied", [False, True])
def test_comparison_rows_use_executed_families_including_untestable_members(
    monkeypatch: pytest.MonkeyPatch,
    fraction_count: int,
    ratio_count: int,
    all_tied: bool,
) -> None:
    """Labels follow each real kernel's complete family, including internal p=1s."""
    traits = tuple(TABLE1_TRAITS)
    fractions = traits[:fraction_count]
    ratios = traits[:ratio_count]
    profiles = [
        Profile(
            f"bsdb:{study}",
            f"bsdb:{study}-1",
            direction,
            "C",
            "Feces",
            {
                trait: (
                    0
                    if index == 0
                    else 4
                    if all_tied or (index == 1) == (direction == "inc")
                    else 1
                )
                for index, trait in enumerate(traits)
            },
            10,
        )
        for study in range(1, 5)
        for direction in ("inc", "dec")
    ]
    expected_fractions = compute_direction_shifts(profiles, trait_order=fractions)
    expected_ratios = compute_log_ratio_shifts(profiles, traits=ratios)
    monkeypatch.setattr(
        "mtp.contrast_orientation.compute_direction_shifts",
        partial(compute_direction_shifts, trait_order=fractions),
    )
    monkeypatch.setattr(
        "mtp.contrast_orientation.TABLE1_TRAITS",
        ratios,
    )

    rows = comparison_rows([("synthetic", profiles)])
    fraction_rows = [row for row in rows if row["comparison"] == "fraction"]
    ratio_rows = [row for row in rows if row["comparison"] == "log_ratio"]
    pair_count = len(list(itertools.combinations(ratios, 2)))
    assert len(fraction_rows) == fraction_count
    assert len(ratio_rows) == pair_count
    assert {row["family"] for row in fraction_rows} == {
        f"Benjamini-Hochberg over {fraction_count} fractions"
    }
    assert {row["family"] for row in ratio_rows} == {
        f"Benjamini-Yekutieli over {pair_count} ratios"
    }
    for actual, expected in zip(fraction_rows, expected_fractions.shifts, strict=True):
        assert (actual["numerator"], actual["p_value"], actual["q_value"]) == (
            expected.trait,
            expected.p_value,
            expected.q_value,
        )
    for actual, ratio in zip(ratio_rows, expected_ratios.summaries, strict=True):
        assert (
            actual["numerator"],
            actual["reference"],
            actual["p_value"],
            actual["q_value"],
        ) == (ratio.numerator, ratio.reference, ratio.sign_p, ratio.by_q)
    for family in (fraction_rows, ratio_rows):
        finite_count = sum(row["p_value"] is not None for row in family)
        assert finite_count < len(family)
        assert (finite_count == 0) is all_tied


def test_comparison_rows_emit_untestable_rows_for_an_empty_population_without_changing_others(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """An empty population emits untestable rows across both families without aborting (#1838)."""
    assert comparison_rows.__doc__ is not None
    normalized_doc = " ".join(comparison_rows.__doc__.split())
    assert "An empty population emits untestable rows across both families." in normalized_doc
    assert "no paired profiles" not in comparison_rows.__doc__.lower()
    assert (
        "(``ratio_traits``, defaulting to ``TABLE1_TRAITS``; every named trait must be a "
        "known trait in ``TABLE1_TRAITS``)."
    ) in normalized_doc
    traits = tuple(TABLE1_TRAITS)
    profiles = [
        Profile(
            f"bsdb:{study}",
            f"bsdb:{study}-1",
            direction,
            "C",
            "Feces",
            {
                trait: (4 if (index % 2 == 0) == (direction == "inc") else 1)
                for index, trait in enumerate(traits)
            },
            10,
        )
        for study in range(1, 5)
        for direction in ("inc", "dec")
    ]
    baseline = comparison_rows([("populated", profiles)])
    try:
        combined = comparison_rows(
            [("empty_before", []), ("populated", profiles), ("empty_after", [])]
        )
    except ValueError as error:
        pytest.fail(f"an empty population must not abort the table: {error}")
    assert [row for row in combined if row["population"] == "populated"] == baseline

    expected_pairs = list(itertools.combinations(traits, 2))
    for empty_name in ("empty_before", "empty_after"):
        empty_rows = [row for row in combined if row["population"] == empty_name]
        fraction_rows = [row for row in empty_rows if row["comparison"] == "fraction"]
        ratio_rows = [row for row in empty_rows if row["comparison"] == "log_ratio"]
        assert len(fraction_rows) == len(traits)
        assert len(ratio_rows) == len(expected_pairs)
        assert [row["numerator"] for row in fraction_rows] == list(traits)
        assert {row["reference"] for row in fraction_rows} == {""}
        assert [(row["numerator"], row["reference"]) for row in ratio_rows] == expected_pairs
        assert {row["family"] for row in fraction_rows} == {
            f"Benjamini-Hochberg over {len(traits)} fractions"
        }
        assert {row["family"] for row in ratio_rows} == {
            f"Benjamini-Yekutieli over {len(expected_pairs)} ratios"
        }
        for row in empty_rows:
            assert (
                row["population_studies"],
                row["population_experiments"],
                row["studies_scored"],
                row["studies_up"],
                row["studies_down"],
                row["studies_tied"],
                row["share_up"],
                row["ci_low"],
                row["ci_high"],
                row["mean_effect"],
                row["p_value"],
                row["q_value"],
            ) == (0, 0, 0, 0, 0, 0, None, None, None, None, None, None)

    # Non-default trait list: the empty population's labels and pairs match populated rows (#1888).
    custom_fractions = traits[:3]
    custom_ratios = traits[:4]
    expected_custom_pairs = list(itertools.combinations(custom_ratios, 2))
    with monkeypatch.context() as patch:
        patch.setattr(
            "mtp.contrast_orientation.compute_direction_shifts",
            partial(compute_direction_shifts, trait_order=custom_fractions),
        )
        patch.setattr("mtp.contrast_orientation.TABLE1_TRAITS", custom_ratios)
        try:
            custom_rows = comparison_rows([("custom_populated", profiles), ("custom_empty", [])])
        except ValueError as error:
            pytest.fail(f"an empty population must not abort the table: {error}")
    pop_frac = [
        row
        for row in custom_rows
        if row["population"] == "custom_populated" and row["comparison"] == "fraction"
    ]
    pop_rat = [
        row
        for row in custom_rows
        if row["population"] == "custom_populated" and row["comparison"] == "log_ratio"
    ]
    custom_frac = [
        row
        for row in custom_rows
        if row["population"] == "custom_empty" and row["comparison"] == "fraction"
    ]
    custom_rat = [
        row
        for row in custom_rows
        if row["population"] == "custom_empty" and row["comparison"] == "log_ratio"
    ]
    assert len(custom_frac) == len(pop_frac) == 3
    assert len(custom_rat) == len(pop_rat) == 6
    assert (
        {row["family"] for row in custom_frac}
        == {row["family"] for row in pop_frac}
        == {"Benjamini-Hochberg over 3 fractions"}
    )
    assert (
        {row["family"] for row in custom_rat}
        == {row["family"] for row in pop_rat}
        == {"Benjamini-Yekutieli over 6 ratios"}
    )
    assert (
        [(row["numerator"], row["reference"]) for row in custom_rat]
        == [(row["numerator"], row["reference"]) for row in pop_rat]
        == expected_custom_pairs
    )

    try:
        explicit_rows = comparison_rows(
            [("explicit_populated", profiles), ("explicit_empty", [])],
            ratio_traits=custom_ratios,
        )
    except ValueError as error:
        pytest.fail(f"an empty population must not abort the table: {error}")
    explicit_pop_rat = [
        row
        for row in explicit_rows
        if row["population"] == "explicit_populated" and row["comparison"] == "log_ratio"
    ]
    explicit_emp_rat = [
        row
        for row in explicit_rows
        if row["population"] == "explicit_empty" and row["comparison"] == "log_ratio"
    ]
    assert (
        {row["family"] for row in explicit_emp_rat}
        == {row["family"] for row in explicit_pop_rat}
        == {"Benjamini-Yekutieli over 6 ratios"}
    )
    assert (
        [(row["numerator"], row["reference"]) for row in explicit_emp_rat]
        == [(row["numerator"], row["reference"]) for row in explicit_pop_rat]
        == expected_custom_pairs
    )

    # Invalid ratio trait list still raises on an empty population; unpaired profiles still raise.
    for bad_traits in (("only_one",), (traits[0], traits[0])):
        with monkeypatch.context() as patch:
            patch.setattr("mtp.contrast_orientation.TABLE1_TRAITS", bad_traits)
            with pytest.raises(ValueError, match="Provide at least two distinct traits"):
                comparison_rows([("invalid_traits", [])])
        with pytest.raises(ValueError, match="Provide at least two distinct traits"):
            comparison_rows([("invalid_traits", [])], ratio_traits=bad_traits)
    for unknown_traits in (
        ("no_such_a", "no_such_b"),
        (traits[0], "no_such_b"),
        ("no_such_a", traits[0]),
    ):
        for target_populations in (
            (),
            (("empty", ()),),
            (("populated", profiles),),
            (("populated", profiles), ("empty", ())),
        ):
            with pytest.raises(
                ValueError,
                match=r"Unknown ratio trait no_such_[ab]: ratio_traits must be Table 1 traits",
            ):
                comparison_rows(target_populations, ratio_traits=unknown_traits)
    extra_profiles = [
        replace(profile, guild_counts={**profile.guild_counts, "extra_profile_column": 3})
        for profile in profiles
    ]
    for extra_ratio_traits in (
        (traits[0], "extra_profile_column"),
        ("extra_profile_column", traits[0]),
    ):
        for target_populations in (
            (),
            (("empty", ()),),
            (("populated", extra_profiles),),
            (("populated", extra_profiles), ("empty", ())),
        ):
            with pytest.raises(
                ValueError,
                match=(
                    r"Unknown ratio trait extra_profile_column: "
                    r"ratio_traits must be Table 1 traits"
                ),
            ):
                comparison_rows(target_populations, ratio_traits=extra_ratio_traits)
    with pytest.raises(ValueError, match="1 increased vs 0 decreased profiles"):
        comparison_rows([("unpaired_study", [profiles[0]])])
    with pytest.raises(ValueError, match="Unpaired experiment"):
        comparison_rows(
            [
                (
                    "unpaired_experiment",
                    [profiles[0], replace(profiles[1], experiment_id="bsdb:1-2")],
                )
            ]
        )

    # CLI writes disease_control_comparison.tsv even when Tier 1 and strong subsets are empty.
    dump, profile_tsv, experiments_tsv = _write_corpus(tmp_path)
    dump.write_text(
        dump.read_text(encoding="utf-8").replace("Healthy controls", "Non-colorectal cancer"),
        encoding="utf-8",
    )
    output = tmp_path / "out_tier2_only"
    assert (
        main(
            [
                "--full-dump-input",
                str(dump),
                "--profile-input",
                str(profile_tsv),
                "--experiments-input",
                str(experiments_tsv),
                "--output-dir",
                str(output),
                "--code-commit",
                "abc123",
            ]
        )
        == 0
    )
    written = list(
        csv.DictReader(
            (output / "disease_control_comparison.tsv").open(encoding="utf-8"), delimiter="\t"
        )
    )
    for empty_tier in (
        "disease_vs_healthy_strong_evidence",
        "disease_vs_healthy_tier1_healthy_normal",
    ):
        tier_rows = [row for row in written if row["population"] == empty_tier]
        assert len(tier_rows) == len(traits) + len(expected_pairs)
        assert {row["population_studies"] for row in tier_rows} == {"0"}
        assert {row["population_experiments"] for row in tier_rows} == {"0"}
        assert {row["studies_scored"] for row in tier_rows} == {"0"}
        assert {row["p_value"] for row in tier_rows} == {""}
        assert {row["q_value"] for row in tier_rows} == {""}
    populated_tier = [
        row for row in written if row["population"] == "disease_vs_healthy_tier1_2_subject_controls"
    ]
    assert {row["population_studies"] for row in populated_tier} == {"1"}
    assert {row["population_experiments"] for row in populated_tier} == {"1"}


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
    assert CLASSIFIER_VERSION == "contrast-orientation-v6"
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

    # Tier 3: treated disease group vs healthy control (native and reversed)
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

    treated_conf_rev = classify(
        replace(
            BASE,
            condition="HIV infection",
            group_0_name="HIV Patients treated with Maraviroc",
            group_1_name="Controls",
        )
    )
    assert treated_conf_rev.contrast_class == DISEASE_VS_HEALTHY
    assert treated_conf_rev.orientation == REVERSED
    assert treated_conf_rev.contrast_type == CASE_VS_CONTROL_SUBJECT
    assert treated_conf_rev.control_tier == TIER_3_HOSPITAL_COMORBID_OR_TISSUE
    assert "treated_disease_group" in treated_conf_rev.flags.split(",")

    # Tier 3: hospital or symptomatic non-healthy control
    hosp_ctrl = classify(
        replace(
            BASE,
            condition="Colorectal cancer",
            group_0_name="Hospital controls",
            group_1_name="CRC patients",
        )
    )
    assert hosp_ctrl.contrast_class == DISEASE_VS_HEALTHY
    assert hosp_ctrl.orientation == NATIVE
    assert hosp_ctrl.contrast_type == CASE_VS_CONTROL_SUBJECT
    assert hosp_ctrl.control_tier == TIER_3_HOSPITAL_COMORBID_OR_TISSUE
    assert "hospital_or_symptomatic_control" in hosp_ctrl.flags.split(",")

    # Cross-sectional recovered vs healthy control routes to OTHER with RESOLVED_DISEASE
    resolved_cs = classify(
        replace(
            BASE,
            condition="COVID-19",
            group_0_name="Healthy controls",
            group_1_name="Recovered patients",
        )
    )
    assert resolved_cs.contrast_class == OTHER
    assert resolved_cs.contrast_type == RESOLVED_DISEASE
    assert resolved_cs.control_tier == NOT_APPLICABLE

    # Longitudinal active vs recovered routes to LONGITUDINAL with RESOLVED_DISEASE
    resolved_long = classify(
        replace(
            BASE,
            condition="COVID-19",
            group_0_name="SARS-CoV-2 infected patients",
            group_1_name="SARS-CoV-2 recovered individuals",
        )
    )
    assert resolved_long.contrast_class == LONGITUDINAL
    assert resolved_long.contrast_type == RESOLVED_DISEASE
    assert resolved_long.control_tier == NOT_APPLICABLE


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
#: The classifier-version guard (#793) classifies these cases, so they stay as they are:
#: a new case goes in a test of its own.
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


def test_known_miss_a_treated_group_named_by_time_keeps_the_treatment_reading() -> None:
    """29808701-1: "for more than 1 year" reads as a duration band, not a longitudinal visit (#581).

    With the duration qualifier masked out of :func:`_timepoints`, the treated-disease-group
    rule classifies HIV-uninfected women against HIV-positive women on ART (>1 year) as
    ``disease_vs_healthy`` Tier 3 (``treated_disease_group`` and ``age_band_not_timepoint``),
    resolving the limitation pinned in #581 / #628.
    """
    found = classify(
        replace(
            BASE,
            condition="HIV infection",
            group_0_name="Age-matched HIV-uninfected women",
            group_1_name="HIV-positive women who have been taking ART for more than 1 year",
        )
    )
    assert (found.contrast_class, found.orientation) == (DISEASE_VS_HEALTHY, NATIVE)
    assert found.contrast_type == CASE_VS_CONTROL_SUBJECT
    assert found.control_tier == TIER_3_HOSPITAL_COMORBID_OR_TISSUE
    assert "treated_disease_group" in found.flags.split(",")
    assert "age_band_not_timepoint" in found.flags.split(",")


def test_spaced_pre_art_and_parenthesized_art_naive_read_as_untreated_baseline_cases() -> None:
    """Untreated baseline cases versus controls, and two untreated Fiebig stages.

    33324725-1 and 38310301-1 compare baseline cases with controls; 33324725-4
    compares two HIV-positive groups before treatment (#581, #806).
    """
    pre_art = classify(
        replace(
            BASE,
            condition="HIV infection",
            group_0_name="Human immunodeficiency virus negative (HIV-)Uninfected Controls",
            group_1_name=(
                "Human immunodeficiency virus positive(HIV+) pre - antiretroviral therapy(ART)"
            ),
        )
    )
    assert (pre_art.contrast_class, pre_art.orientation) == (DISEASE_VS_HEALTHY, NATIVE)
    assert pre_art.contrast_type == CASE_VS_CONTROL_SUBJECT
    assert pre_art.control_tier == TIER_1_HEALTHY_NORMAL
    assert "negated_treatment_term" in pre_art.flags.split(",")
    assert "disease_sample_dated_at_baseline" in pre_art.flags.split(",")

    # 33324725-4: verbatim export labels compare Fiebig stages, not frailty.
    pre_art_stages = classify(
        replace(
            BASE,
            condition="HIV infection",
            group_0_name=(
                "Human immunodeficiency virus positive(HIV+) pre- antiretroviral therapy (ART) "
                "in Fiebig stage 1-2(F1-2)"
            ),
            group_1_name=(
                "Human immunodeficiency virus positive (HIV+) pre- antiretroviral therapy(ART) "
                "in Fiebig stage 3(F3)"
            ),
        )
    )
    assert (pre_art_stages.contrast_class, pre_art_stages.orientation) == (
        DISEASE_VS_DISEASE,
        NOT_APPLICABLE,
    )
    assert pre_art_stages.contrast_type == DISEASE_VS_DISEASE
    assert pre_art_stages.control_tier == NOT_APPLICABLE
    assert "negated_treatment_term" in pre_art_stages.flags.split(",")
    assert "treated_disease_group" not in pre_art_stages.flags.split(",")

    art_naive = classify(
        replace(
            BASE,
            condition="HIV infection",
            group_0_name="Healthy controls (HC)",
            group_1_name="Antiretroviral therapy (ART) Na\u00efve at baseline",
        )
    )
    assert (art_naive.contrast_class, art_naive.orientation) == (DISEASE_VS_HEALTHY, NATIVE)
    assert art_naive.contrast_type == CASE_VS_CONTROL_SUBJECT
    assert art_naive.control_tier == TIER_1_HEALTHY_NORMAL
    assert "negated_treatment_term" in art_naive.flags.split(",")
    assert "disease_sample_dated_at_baseline" in art_naive.flags.split(",")


def test_treatment_and_multivalued_conditions_do_not_let_group_names_override_scope() -> None:
    """#626: group-name comparisons under TREATMENT_CONDITION or multi-valued conditions
    stay excluded so FMT recipient mouse trials (38477534-4) and comorbidity-matched cohorts
    (41011452-3) are not misclassified as primary disease-vs-healthy contrasts.
    """
    # 25057045-1: curated under "Response to antiviral drug" (ART in group_1_definition)
    assert _classify(
        condition="Response to antiviral drug",
        group_0_name="HIV-uninfected controls",
        group_1_name="Chronic HIV-infected cases",
    ) == (TREATMENT, NOT_APPLICABLE)

    # 37660163-13: curated under "Response to antiviral drug" (HCV/HIV coinfection)
    assert _classify(
        condition="Response to antiviral drug",
        group_0_name="Healthy controls",
        group_1_name="patients with HCV/HIV coinfection at baseline",
    ) == (TREATMENT, NOT_APPLICABLE)

    # 38477534-4: FMT recipient trial with "Healthy controls" vs "Chronic Constipation"
    assert _classify(
        condition="Response to transplant",
        group_0_name="Healthy controls - Fecal Microbiota Transplantation (HC-FMT)",
        group_1_name="Chronic Constipation - Fecal Microbiota Transplantation (CC-FMT)",
    ) == (TREATMENT, NOT_APPLICABLE)

    # 41011452-3: both groups have oral candidiasis; only general comparison tokens override
    assert _classify(
        condition="HIV infection,Oral candidiasis",
        group_0_name="HIV-uninfected participants with oral candidiasis",
        group_1_name="HIV-infected participants with oral candidiasis",
    ) == (DISEASE_VS_DISEASE, NOT_APPLICABLE)


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


@pytest.mark.parametrize(
    ("group_0", "group_1"),
    [
        ("Healthy controls", "HIV patients on HAART"),
        ("HIV-negative controls", "HIV patients on cART"),
        ("HIV-negative controls", "People with HIV taking ART"),
        ("HIV-negative controls", "People with HIV receiving ART"),
        ("HIV-negative controls", "People with HIV under ART"),
        ("PLWOH", "PLWH on ART"),
        ("HIV-negative controls", "PWH on ART"),
        ("Healthy controls", "Human immunodeficiency virus positive on ART"),
    ],
)
def test_every_antiretroviral_cue_and_name_for_hiv_reaches_the_treated_rule(
    group_0: str, group_1: str
) -> None:
    """ART spellings that TREATMENT_TERMS cannot read, and PLWH, PWH or the spelled-out
    virus as the name of the disease group, all reach the treated-cases rule (#627, #628)."""
    found = classify(
        replace(BASE, condition="HIV infection", group_0_name=group_0, group_1_name=group_1)
    )
    assert (found.contrast_class, found.orientation) == (DISEASE_VS_HEALTHY, NATIVE)
    assert found.control_tier == TIER_3_HOSPITAL_COMORBID_OR_TISSUE
    assert "treated_disease_group" in found.flags.split(",")


@pytest.mark.parametrize(
    ("condition", "group_0", "group_1"),
    [
        # Group 1 does not name the condition
        ("Obesity", "Lean controls", "Patients treated with metformin"),
        # Group 1 is a comparison group itself
        ("Obesity", "Lean controls", "Non-obese controls treated with metformin"),
        # Group 0 carries a treatment word too (bsdb:32733788-1)
        (
            "Response to immunochemotherapy",
            "non-responders (patients treated with conventional dose)",
            "responders (patients treated with metronomic dose of capecitabine)",
        ),
    ],
)
def test_the_treated_rule_needs_each_of_its_guards(
    condition: str, group_0: str, group_1: str
) -> None:
    found = classify(replace(BASE, condition=condition, group_0_name=group_0, group_1_name=group_1))
    assert found.contrast_class == TREATMENT
    assert "treated_disease_group" not in found.flags.split(",")


@pytest.mark.parametrize(
    ("condition", "group_0", "group_1"),
    [
        # The "II" of "Type II" is a number, not an acronym of "insulin injection"
        (
            "Type II diabetes mellitus",
            "T2D without insulin injection",
            "T2D with insulin injection",
        ),
        # PLWOH is a comparison only where HIV is the condition, not any "hi" initialism
        ("Helminth infection", "PLWOH", "PLWH"),
        ("HOMA-IR", "PLWOH", "PLWH"),
    ],
)
def test_acronym_and_hiv_readings_stay_with_their_conditions(
    condition: str, group_0: str, group_1: str
) -> None:
    """Both were unresolved in v3 and must stay so (#627)."""
    assert _classify(condition=condition, group_0_name=group_0, group_1_name=group_1) == (
        UNRESOLVED,
        NOT_APPLICABLE,
    )


def test_acronyms_are_the_capital_words_of_the_condition() -> None:
    """Pins case, length and the rest of the phrase (#628)."""
    # Lower-case words such as the "to" of "Response to …" are never acronyms
    assert _acronym_spelled_out(["toxic", "outcome"], "Response to antiviral drug") == [
        "toxic",
        "outcome",
    ]
    # The rest of the phrase follows the acronym
    assert _acronym_spelled_out(
        ["human", "immunodeficiency", "virus", "positive"], "HIV infection"
    ) == ["hiv", "positive"]
    # Two-letter acronyms count
    assert _acronym_spelled_out(["crohn's", "disease", "free"], "CD flare") == ["cd", "free"]
    # Roman numerals do not
    assert _acronym_spelled_out(["insulin", "injection"], "Type II diabetes mellitus") == [
        "insulin",
        "injection",
    ]


def test_reversed_antiretroviral_cases_require_a_strong_control_group() -> None:
    """Reversed ON_ANTIRETROVIRAL contrasts need strong_1 to reverse (#574, #630)."""
    for g0 in ("People living with HIV on ART", "People living with HIV on HAART"):
        rev_strong = classify(
            replace(
                BASE,
                condition="HIV infection",
                group_0_name=g0,
                group_1_name="Healthy controls",
            )
        )
        assert (rev_strong.contrast_class, rev_strong.orientation) == (
            DISEASE_VS_HEALTHY,
            REVERSED,
        )
        assert rev_strong.contrast_type == CASE_VS_CONTROL_SUBJECT
        assert rev_strong.control_tier == TIER_3_HOSPITAL_COMORBID_OR_TISSUE
        assert "treated_disease_group" in rev_strong.flags.split(",")

    rev_weak_art = classify(
        replace(
            BASE,
            condition="HIV infection",
            group_0_name="People living with HIV on ART",
            group_1_name="HIV-",
        )
    )
    assert (rev_weak_art.contrast_class, rev_weak_art.orientation) == (TREATMENT, NOT_APPLICABLE)
    assert "comparison_group_present" in rev_weak_art.flags.split(",")
    assert "treated_disease_group" not in rev_weak_art.flags.split(",")

    rev_weak_haart = classify(
        replace(
            BASE,
            condition="HIV infection",
            group_0_name="People living with HIV on HAART",
            group_1_name="HIV-",
        )
    )
    assert (rev_weak_haart.contrast_class, rev_weak_haart.orientation) == (
        UNRESOLVED,
        NOT_APPLICABLE,
    )
    assert "possible_reversal_unverified" in rev_weak_haart.flags.split(",")
    assert "treated_disease_group" not in rev_weak_haart.flags.split(",")


def test_precursor_lesions_in_controls_are_tier3_subject_controls_not_comorbidity() -> None:
    """31171880-1..4 (Yachida et al. 2019): CRC controls with polyps are Tier 3 subject
    controls (hospital_or_symptomatic_control), not comorbidity_within_disease (#672).
    Explicitly negated adenoma controls remain Tier 1 or Tier 2.
    """
    yachida = classify(
        replace(
            BASE,
            condition="Colorectal cancer",
            group_0_name="healthy control (normal and a few polyps)",
            group_1_name="stage I/II (SI/II)",
        )
    )
    assert (yachida.contrast_class, yachida.orientation) == (DISEASE_VS_HEALTHY, NATIVE)
    assert yachida.contrast_type == CASE_VS_CONTROL_SUBJECT
    assert yachida.control_tier == TIER_3_HOSPITAL_COMORBID_OR_TISSUE
    assert "hospital_or_symptomatic_control" in yachida.flags.split(",")
    assert "comorbidity_control" not in yachida.flags.split(",")

    # 22622349-1: negated adenoma with "controls" stays Tier 1
    adenoma_free = classify(
        replace(
            BASE,
            condition="Adenoma",
            group_0_name="Adenoma-free (controls)",
            group_1_name="Adenoma cases",
        )
    )
    assert (adenoma_free.contrast_class, adenoma_free.orientation) == (
        DISEASE_VS_HEALTHY,
        NATIVE,
    )
    assert adenoma_free.contrast_type == CASE_VS_CONTROL_SUBJECT
    assert adenoma_free.control_tier == TIER_1_HEALTHY_NORMAL
    assert "hospital_or_symptomatic_control" not in adenoma_free.flags.split(",")

    # 23335968-1: negated adenoma without strong token stays Tier 2
    non_adenoma = classify(
        replace(
            BASE,
            condition="Colorectal adenoma",
            group_0_name="Non-adenoma",
            group_1_name="Colorectal adenoma",
        )
    )
    assert (non_adenoma.contrast_class, non_adenoma.orientation) == (
        DISEASE_VS_HEALTHY,
        NATIVE,
    )
    assert non_adenoma.contrast_type == CASE_VS_CONTROL_SUBJECT
    assert non_adenoma.control_tier == TIER_2_CONDITION_NEGATIVE_SCREENED
    assert "hospital_or_symptomatic_control" not in non_adenoma.flags.split(",")


#: The four experiments that contrast-orientation-v6 reclassifies (#739, #793), with the
#: condition and group names that its version history in methods/contrast_orientation.md
#: lists, and the v6 reading it gives them. v5 read all four as ``treatment_exposure``.
#: The classifier-version guard (#793) classifies these cases, so they stay as they are:
#: a new case goes in a test of its own.
V6_RECLASSIFIED_FIXTURES: tuple[tuple[str, dict[str, str], tuple[str, str, str, str]], ...] = (
    (
        "29808701-1 therapy for more than 1 year is a duration, not a visit",
        {
            "condition": "HIV infection",
            "group_0_name": "Age-matched HIV-uninfected women",
            "group_1_name": "HIV-positive women who have been taking ART for more than 1 year",
        },
        (
            DISEASE_VS_HEALTHY,
            NATIVE,
            CASE_VS_CONTROL_SUBJECT,
            TIER_3_HOSPITAL_COMORBID_OR_TISSUE,
        ),
    ),
    (
        "33324725-1 a spaced pre - ART baseline against uninfected controls",
        {
            "condition": "HIV infection",
            "group_0_name": "Human immunodeficiency virus negative (HIV-)Uninfected Controls",
            "group_1_name": (
                "Human immunodeficiency virus positive(HIV+) pre - antiretroviral therapy(ART)"
            ),
        },
        (DISEASE_VS_HEALTHY, NATIVE, CASE_VS_CONTROL_SUBJECT, TIER_1_HEALTHY_NORMAL),
    ),
    (
        "33324725-4 two Fiebig stages, both before ART",
        {
            "condition": "HIV infection",
            "group_0_name": (
                "Human immunodeficiency virus positive(HIV+) pre- antiretroviral therapy (ART) "
                "in Fiebig stage 1-2(F1-2)"
            ),
            "group_1_name": (
                "Human immunodeficiency virus positive (HIV+) pre- antiretroviral therapy(ART) "
                "in Fiebig stage 3(F3)"
            ),
        },
        (DISEASE_VS_DISEASE, NOT_APPLICABLE, DISEASE_VS_DISEASE, NOT_APPLICABLE),
    ),
    (
        "38310301-1 ART-naive cases at baseline against healthy controls",
        {
            "condition": "HIV infection",
            "group_0_name": "Healthy controls (HC)",
            "group_1_name": "Antiretroviral therapy (ART) Na\u00efve at baseline",
        },
        (DISEASE_VS_HEALTHY, NATIVE, CASE_VS_CONTROL_SUBJECT, TIER_1_HEALTHY_NORMAL),
    ),
)


@pytest.mark.parametrize(
    ("changes", "expected"),
    [
        pytest.param(changes, expected, id=name)
        for name, changes, expected in V6_RECLASSIFIED_FIXTURES
    ],
)
def test_the_experiments_v6_reclassifies_read_as_its_version_history_says(
    changes: dict[str, str], expected: tuple[str, str, str, str]
) -> None:
    assert _stratified(**changes) == expected


#: The classifier-version guard's pins (#793): for each CLASSIFIER_VERSION, SHA-256 digests
#: of the guard's cases and of the classifier's decisions on them. A pinned version's digests
#: never change. A rule change that moves any case needs a new CLASSIFIER_VERSION, an entry in
#: the version history of methods/contrast_orientation.md, and new digests under the new
#: version. The frozen per-case reference below is validated against those original pins.
CLASSIFIER_GUARD_DIGESTS = (
    Path(__file__).parent / "fixtures" / "contrast_orientation_guard_digests.json"
)


# Frozen from the existing committed guard cases only (#803); no new export rows.
_GUARD_REFERENCE_VERSION = "contrast-orientation-v6"
_GUARD_REFERENCE_COMMON: dict[str, str] = {
    "body_site": "",
    "condition": "HIV infection",
    "description": "",
    "doi": "",
    "efo_id": "",
    "experiment_id": "bsdb:1-1",
    "group_0_name": "Healthy Controls",
    "group_0_sample_size": "",
    "group_1_definition": "",
    "group_1_name": "18-62 years",
    "group_1_sample_size": "",
    "host_species": "",
    "pmid": "1",
    "study_design": "",
    "study_id": "bsdb:1",
    "uberon_id": "",
}

_GUARD_REFERENCE_DECISIONS: tuple[tuple[str, str, str, str, str], ...] = (
    (
        "disease_vs_disease",
        "not_applicable",
        "disease_vs_disease",
        "not_applicable",
        "age_band_not_timepoint,negated_treatment_term,no_group_0_definition_in_export,shared_timepoint",
    ),
    (
        "disease_vs_disease",
        "not_applicable",
        "disease_vs_disease",
        "not_applicable",
        "age_band_not_timepoint,no_group_0_definition_in_export,no_group_1_definition",
    ),
    (
        "disease_vs_disease",
        "not_applicable",
        "disease_vs_disease",
        "not_applicable",
        "comparison_reading_withdrawn,no_group_0_definition_in_export",
    ),
    (
        "disease_vs_disease",
        "not_applicable",
        "disease_vs_disease",
        "not_applicable",
        "no_group_0_definition_in_export",
    ),
    (
        "disease_vs_disease",
        "not_applicable",
        "disease_vs_disease",
        "not_applicable",
        "no_group_0_definition_in_export,no_group_1_definition",
    ),
    (
        "disease_vs_healthy",
        "native",
        "case_vs_control_subject",
        "tier_1_healthy_normal",
        "disease_sample_dated_at_baseline,negated_treatment_term,no_group_0_definition_in_export",
    ),
    (
        "disease_vs_healthy",
        "native",
        "case_vs_control_subject",
        "tier_1_healthy_normal",
        (
            "disease_sample_dated_at_baseline,negated_treatment_term,"
            "no_group_0_definition_in_export,no_group_1_definition,shared_treatment_term"
        ),
    ),
    (
        "disease_vs_healthy",
        "native",
        "case_vs_control_subject",
        "tier_1_healthy_normal",
        "no_group_0_definition_in_export",
    ),
    (
        "disease_vs_healthy",
        "native",
        "case_vs_control_subject",
        "tier_1_healthy_normal",
        "no_group_0_definition_in_export,no_group_1_definition",
    ),
    (
        "disease_vs_healthy",
        "native",
        "case_vs_control_subject",
        "tier_1_healthy_normal",
        "no_group_0_definition_in_export,no_group_1_definition,shared_timepoint",
    ),
    (
        "disease_vs_healthy",
        "native",
        "case_vs_control_subject",
        "tier_1_healthy_normal",
        "no_group_0_definition_in_export,shared_timepoint",
    ),
    (
        "disease_vs_healthy",
        "native",
        "case_vs_control_subject",
        "tier_2_condition_negative_screened",
        "no_group_0_definition_in_export,no_group_1_definition,weak_comparison_evidence",
    ),
    (
        "disease_vs_healthy",
        "native",
        "case_vs_control_subject",
        "tier_2_condition_negative_screened",
        "no_group_0_definition_in_export,shared_timepoint,weak_comparison_evidence",
    ),
    (
        "disease_vs_healthy",
        "native",
        "case_vs_control_subject",
        "tier_2_condition_negative_screened",
        "no_group_0_definition_in_export,weak_comparison_evidence",
    ),
    (
        "disease_vs_healthy",
        "native",
        "case_vs_control_subject",
        "tier_3_hospital_comorbid_or_tissue",
        "age_band_not_timepoint,no_group_0_definition_in_export,treated_disease_group",
    ),
    (
        "disease_vs_healthy",
        "native",
        "case_vs_control_subject",
        "tier_3_hospital_comorbid_or_tissue",
        "no_group_0_definition_in_export,no_group_1_definition,treated_disease_group",
    ),
    (
        "disease_vs_healthy",
        "native",
        "case_vs_control_subject",
        "tier_3_hospital_comorbid_or_tissue",
        "no_group_0_definition_in_export,treated_disease_group",
    ),
    (
        "disease_vs_healthy",
        "native",
        "case_vs_control_subject",
        "tier_3_hospital_comorbid_or_tissue",
        "no_group_0_definition_in_export,treated_disease_group,weak_comparison_evidence",
    ),
    (
        "disease_vs_healthy",
        "native",
        "comorbidity_within_disease",
        "tier_3_hospital_comorbid_or_tissue",
        (
            "comorbidity_control,no_group_0_definition_in_export,no_group_1_definition,"
            "weak_comparison_evidence"
        ),
    ),
    (
        "disease_vs_healthy",
        "native",
        "comorbidity_within_disease",
        "tier_3_hospital_comorbid_or_tissue",
        "comorbidity_control,no_group_0_definition_in_export,weak_comparison_evidence",
    ),
    (
        "disease_vs_healthy",
        "native",
        "paired_tissue_autologous",
        "tier_3_hospital_comorbid_or_tissue",
        "no_group_0_definition_in_export,no_group_1_definition,paired_tissue_comparison",
    ),
    (
        "longitudinal_timepoint",
        "not_applicable",
        "longitudinal_timepoint",
        "not_applicable",
        "no_group_0_definition_in_export,no_group_1_definition",
    ),
    (
        "other",
        "not_applicable",
        "other",
        "not_applicable",
        "age_band_not_timepoint,no_group_0_definition_in_export",
    ),
    (
        "other",
        "not_applicable",
        "other",
        "not_applicable",
        "age_band_not_timepoint,no_group_0_definition_in_export,no_group_1_definition",
    ),
    (
        "other",
        "not_applicable",
        "other",
        "not_applicable",
        "comparison_reading_withdrawn,no_group_0_definition_in_export,no_group_1_definition",
    ),
    ("other", "not_applicable", "other", "not_applicable", "no_group_0_definition_in_export"),
    (
        "other",
        "not_applicable",
        "other",
        "not_applicable",
        "no_group_0_definition_in_export,no_group_1_definition",
    ),
    (
        "other",
        "not_applicable",
        "resolved_disease",
        "not_applicable",
        "no_group_0_definition_in_export,no_group_1_definition",
    ),
    (
        "treatment_exposure",
        "not_applicable",
        "treatment_exposure",
        "not_applicable",
        (
            "comparison_group_present,negated_treatment_term,"
            "no_group_0_definition_in_export,shared_timepoint"
        ),
    ),
    (
        "treatment_exposure",
        "not_applicable",
        "treatment_exposure",
        "not_applicable",
        "comparison_group_present,no_group_0_definition_in_export",
    ),
    (
        "treatment_exposure",
        "not_applicable",
        "treatment_exposure",
        "not_applicable",
        "comparison_group_present,no_group_0_definition_in_export,no_group_1_definition",
    ),
    (
        "treatment_exposure",
        "not_applicable",
        "treatment_exposure",
        "not_applicable",
        "negated_treatment_term,no_group_0_definition_in_export",
    ),
    (
        "treatment_exposure",
        "not_applicable",
        "treatment_exposure",
        "not_applicable",
        "no_group_0_definition_in_export",
    ),
    (
        "treatment_exposure",
        "not_applicable",
        "treatment_exposure",
        "not_applicable",
        "no_group_0_definition_in_export,no_group_1_definition",
    ),
    (
        "unresolved",
        "not_applicable",
        "not_applicable",
        "not_applicable",
        "age_band_not_timepoint,no_group_0_definition_in_export",
    ),
    (
        "unresolved",
        "not_applicable",
        "not_applicable",
        "not_applicable",
        "no_group_0_definition_in_export",
    ),
    (
        "unresolved",
        "not_applicable",
        "not_applicable",
        "not_applicable",
        "no_group_0_definition_in_export,no_group_1_definition",
    ),
)

_GUARD_REFERENCE_CASES: tuple[tuple[dict[str, str], int], ...] = (
    (
        {
            "condition": "Chronic kidney disease",
            "experiment_id": "bsdb:30300382-1",
            "group_0_name": "healthy controls",
            "group_1_name": "Chronic kidney disease",
            "pmid": "30300382",
            "study_id": "bsdb:30300382",
        },
        8,
    ),
    (
        {
            "condition": "Arthritis",
            "experiment_id": "bsdb:31842923-1",
            "group_0_name": "healthy controls",
            "group_1_name": "Juvenile Idiopathic Arthritis",
            "pmid": "31842923",
            "study_id": "bsdb:31842923",
        },
        8,
    ),
    (
        {
            "condition": "Colorectal cancer",
            "experiment_id": "bsdb:32083132-1",
            "group_0_name": "Healthy Controls (HC)",
            "group_1_name": "Colorectal Cancer (CC)",
            "pmid": "32083132",
            "study_id": "bsdb:32083132",
        },
        8,
    ),
    (
        {
            "condition": "Human papilloma virus infection",
            "experiment_id": "bsdb:32842982-5",
            "group_0_name": "Healthy control",
            "group_1_name": "HPV+",
            "pmid": "32842982",
            "study_id": "bsdb:32842982",
        },
        8,
    ),
    (
        {
            "condition": "Irritable bowel syndrome",
            "experiment_id": "bsdb:33663411-1",
            "group_0_name": "Healthy Control",
            "group_1_name": "IBS-D",
            "pmid": "33663411",
            "study_id": "bsdb:33663411",
        },
        8,
    ),
    (
        {
            "condition": "Hypoglycemia",
            "experiment_id": "bsdb:38481313-1",
            "group_0_name": "Controls",
            "group_1_name": "Diabetic hypoglycemia",
            "pmid": "38481313",
            "study_id": "bsdb:38481313",
        },
        8,
    ),
    (
        {
            "condition": "Head and neck squamous cell carcinoma",
            "experiment_id": "bsdb:39325441-2",
            "group_1_name": "Head and neck squamous cell cancer (HNSCC) participants",
            "pmid": "39325441",
            "study_id": "bsdb:39325441",
        },
        8,
    ),
    (
        {
            "condition": "Parkinson's disease",
            "experiment_id": "bsdb:39443634-5",
            "group_1_name": "Parkinson\u2019s disease patients (PD patients)",
            "pmid": "39443634",
            "study_id": "bsdb:39443634",
        },
        8,
    ),
    (
        {
            "condition": "Constipation",
            "experiment_id": "bsdb:39588509-4",
            "group_0_name": "Maintenance Hemodialysis without Constipation (MHDNCG)",
            "group_1_name": "Maintenance Hemodialysis with Constipation (MHDCG)",
            "pmid": "39588509",
            "study_id": "bsdb:39588509",
        },
        18,
    ),
    (
        {
            "condition": "Sleep Disorder",
            "experiment_id": "bsdb:40247863-1",
            "group_0_name": "Healthy control group",
            "group_1_name": "Catathrenia group",
            "pmid": "40247863",
            "study_id": "bsdb:40247863",
        },
        8,
    ),
    (
        {
            "condition": "Acute myeloid leukemia",
            "experiment_id": "bsdb:28245856-2",
            "group_0_name": "no infection",
            "group_1_name": "infection",
            "pmid": "28245856",
            "study_id": "bsdb:28245856",
        },
        4,
    ),
    (
        {
            "condition": "Chlamydophila infectious disease,Human papilloma virus infection",
            "experiment_id": "bsdb:29313867-2",
            "group_0_name": "HPV positive women",
            "group_1_name": "patient with chlamydia infection",
            "pmid": "29313867",
            "study_id": "bsdb:29313867",
        },
        4,
    ),
    (
        {
            "condition": "Gastric cancer",
            "experiment_id": "bsdb:30584008-4",
            "group_0_name": "early stage normal",
            "group_1_name": "early stage tumor",
            "pmid": "30584008",
            "study_id": "bsdb:30584008",
        },
        20,
    ),
    (
        {
            "condition": "Oral mucositis",
            "experiment_id": "bsdb:32537095-1",
            "group_0_name": "No/mild oral mucositis (Muc0-1)",
            "group_1_name": "Moderate oral mucositis (Muc2)",
            "pmid": "32537095",
            "study_id": "bsdb:32537095",
        },
        4,
    ),
    (
        {
            "condition": "Schizophrenia",
            "experiment_id": "bsdb:33067643-2",
            "group_0_name": "Schizophrenia patients without violence (NV.SCV)",
            "group_1_name": "Schizophrenia patients with violence (V.SCZ)",
            "pmid": "33067643",
            "study_id": "bsdb:33067643",
        },
        4,
    ),
    (
        {
            "condition": "Psoriasis",
            "experiment_id": "bsdb:33685393-2",
            "group_0_name": "Pustular Psoriasis (PP)",
            "group_1_name": "Psoriasis Vulgaris (PV)",
            "pmid": "33685393",
            "study_id": "bsdb:33685393",
        },
        4,
    ),
    (
        {
            "condition": "Body mass index,Physical activity",
            "experiment_id": "bsdb:36381318-5",
            "group_0_name": "Normal weight/active patients",
            "group_1_name": "Overweight/Obese/active patients",
            "pmid": "36381318",
            "study_id": "bsdb:36381318",
        },
        24,
    ),
    (
        {
            "condition": "Colorectal cancer",
            "experiment_id": "bsdb:38306898-3",
            "group_0_name": "Small-sized colorectal cancer tumors in aoCRC",
            "group_1_name": "Large-sized colorectal cancer tumors in aoCRC",
            "pmid": "38306898",
            "study_id": "bsdb:38306898",
        },
        4,
    ),
    (
        {
            "condition": "Metabolic toxicity",
            "experiment_id": "bsdb:39417117-26",
            "group_0_name": "Non-diabetic",
            "group_1_name": "Diabetic",
            "pmid": "39417117",
            "study_id": "bsdb:39417117",
        },
        11,
    ),
    (
        {
            "condition": "Spontaneous preterm birth",
            "experiment_id": "bsdb:39569373-1",
            "group_0_name": "Term birth-Trimester one (TB-T1)",
            "group_1_name": "Preterm birth-Trimester one (PTB-T1)",
            "pmid": "39569373",
            "study_id": "bsdb:39569373",
        },
        9,
    ),
    (
        {
            "condition": "Nutraceutical",
            "experiment_id": "bsdb:26238090-2",
            "group_0_name": "HGT Group - Healthy controls",
            "group_1_name": "PGT group - colorectal cancer patients perioperative probiotics",
            "pmid": "26238090",
            "study_id": "bsdb:26238090",
        },
        30,
    ),
    (
        {
            "condition": "Oral squamous cell carcinoma",
            "experiment_id": "bsdb:27259999-3",
            "group_0_name": "surgery",
            "group_1_name": "chemotherapy-radiation and surgery",
            "pmid": "27259999",
            "study_id": "bsdb:27259999",
        },
        33,
    ),
    (
        {
            "condition": "Smoking behavior",
            "experiment_id": "bsdb:28844190-1",
            "group_0_name": "non-smokers",
            "group_1_name": "smokers",
            "pmid": "28844190",
            "study_id": "bsdb:28844190",
        },
        33,
    ),
    (
        {
            "condition": "Cesarean section",
            "experiment_id": "bsdb:29538354-4",
            "group_0_name": "Furry pets at home (NO)",
            "group_1_name": "Furry pets at home (YES)",
            "pmid": "29538354",
            "study_id": "bsdb:29538354",
        },
        33,
    ),
    (
        {
            "condition": "Parkinson's disease",
            "experiment_id": "bsdb:31815177-11",
            "group_1_name": "Parkinson's disease Patients treated with entacapone",
            "pmid": "31815177",
            "study_id": "bsdb:31815177",
        },
        15,
    ),
    (
        {
            "condition": "Contraception",
            "experiment_id": "bsdb:33149114-1",
            "group_0_name": "NET-En (norethisterone enanthate ) arm",
            "group_1_name": "COC (combined oral contraceptives) arm",
            "pmid": "33149114",
            "study_id": "bsdb:33149114",
        },
        36,
    ),
    (
        {
            "condition": "Oral squamous cell carcinoma",
            "experiment_id": "bsdb:37542310-2",
            "group_0_name": "Nonsmoking, nondrinking, and dentate healthy controls",
            "group_1_name": "Nonsmoking, nondrinking, and dentate pre-treatment OSCC patients",
            "pmid": "37542310",
            "study_id": "bsdb:37542310",
        },
        6,
    ),
    (
        {
            "condition": "Viral disease",
            "experiment_id": "bsdb:40032727-6",
            "group_0_name": "Anti-NoV GII.4 Sydney 2012 pAb + isotype control nor the specific "
            "pAb + IgG isotype control pAb (NoV+N+C)",
            "group_1_name": "Total bacteria isolated before treatments (T)",
            "pmid": "40032727",
            "study_id": "bsdb:40032727",
        },
        30,
    ),
    (
        {
            "condition": "Diet measurement",
            "experiment_id": "bsdb:40102379-7",
            "group_0_name": "Mother + Infant Control",
            "group_1_name": "Mother + Infant Intervention",
            "pmid": "40102379",
            "study_id": "bsdb:40102379",
        },
        30,
    ),
    (
        {
            "condition": "High altitude adaptation",
            "experiment_id": "bsdb:40257273-5",
            "group_0_name": "G-3(Group 3 prolonged altitude exposure at 4500 m.)",
            "group_1_name": "G-2 (Group 2 acute altitude exposure at 4500 m)",
            "pmid": "40257273",
            "study_id": "bsdb:40257273",
        },
        33,
    ),
    (
        {
            "condition": "Diabetes mellitus",
            "experiment_id": "bsdb:24013136-6",
            "group_0_name": "Prediabetes (Pre-DM)",
            "group_1_name": "Type 2 diabetes mellitus (T2DM)",
            "pmid": "24013136",
            "study_id": "bsdb:24013136",
        },
        4,
    ),
    (
        {
            "condition": "Premature birth",
            "experiment_id": "bsdb:25444008-2",
            "group_0_name": "10-day old Full-term infants(born after uncomplicated pregnancy)",
            "group_1_name": "10-day old Pre-term",
            "pmid": "25444008",
            "study_id": "bsdb:25444008",
        },
        8,
    ),
    (
        {
            "condition": "Chronic constipation",
            "experiment_id": "bsdb:28870194-10",
            "group_0_name": "Placebo moderate subgroup V1 (Visit 1)",
            "group_1_name": "Placebo moderate subgroup V3 (Visit 3)",
            "pmid": "28870194",
            "study_id": "bsdb:28870194",
        },
        21,
    ),
    (
        {
            "condition": "Chronic constipation",
            "experiment_id": "bsdb:28870194-9",
            "group_0_name": "EpiCor Severe subgroup V1(Visit 1)",
            "group_1_name": "EpiCor Severe subgroup V3 (Visit 3)",
            "pmid": "28870194",
            "study_id": "bsdb:28870194",
        },
        21,
    ),
    (
        {
            "condition": "COVID-19",
            "experiment_id": "bsdb:33431578-2",
            "group_0_name": "Healthy controls",
            "group_1_name": "Recovered COVID-19 cases",
            "pmid": "33431578",
            "study_id": "bsdb:33431578",
        },
        27,
    ),
    (
        {
            "condition": "Oral squamous cell carcinoma",
            "experiment_id": "bsdb:37542310-3",
            "group_0_name": "pre-treatment OSCC patients",
            "group_1_name": "post-treatment OSCC patients",
            "pmid": "37542310",
            "study_id": "bsdb:37542310",
        },
        21,
    ),
    (
        {
            "condition": "Gastric cancer",
            "experiment_id": "bsdb:37740058-1",
            "group_0_name": "Gastric cancer patients (Pre-gastrectomy)",
            "group_1_name": "Gastric cancer patients (Post-gastrectomy)",
            "pmid": "37740058",
            "study_id": "bsdb:37740058",
        },
        21,
    ),
    (
        {
            "condition": "Response to diet",
            "experiment_id": "bsdb:39454570-2",
            "group_0_name": "Non-Healthy Low Carbohydrate Diet at Baseline (Non-HLCD)",
            "group_1_name": "Non-Healthy Low Carbohydrate Diet at Endpoint (Non-HLCD)",
            "pmid": "39454570",
            "study_id": "bsdb:39454570",
        },
        21,
    ),
    (
        {
            "condition": "Respiratory Syncytial Virus Infection",
            "experiment_id": "bsdb:39642873-2",
            "group_1_name": "RSV (Respiratory Syncytial Virus) Convalescence",
            "pmid": "39642873",
            "study_id": "bsdb:39642873",
        },
        27,
    ),
    (
        {
            "condition": "Response to allogeneic hematopoietic stem cell transplant",
            "experiment_id": "bsdb:41077329-6",
            "group_0_name": "Baseline samples",
            "group_1_name": "Day 60 post-HSCT samples",
            "pmid": "41077329",
            "study_id": "bsdb:41077329",
        },
        21,
    ),
    (
        {
            "condition": "Ethnic group,Place of residence measurement",
            "experiment_id": "bsdb:20679230-3",
            "group_0_name": "European Children in Italy",
            "group_1_name": "Rural African Children in Burkina Faso",
            "pmid": "20679230",
            "study_id": "bsdb:20679230",
        },
        26,
    ),
    (
        {
            "condition": "Autoantibody measurement",
            "experiment_id": "bsdb:26068542-9",
            "group_0_name": "One autoantibody",
            "group_1_name": "Two to four autoantibodies",
            "pmid": "26068542",
            "study_id": "bsdb:26068542",
        },
        26,
    ),
    (
        {
            "condition": "Age",
            "experiment_id": "bsdb:27474122-7",
            "group_0_name": "2-12 years",
            "pmid": "27474122",
            "study_id": "bsdb:27474122",
        },
        23,
    ),
    (
        {
            "condition": "Stroke",
            "experiment_id": "bsdb:30778376-13",
            "group_0_name": "Low Risk (LR)",
            "group_1_name": "High Risk (HR)",
            "pmid": "30778376",
            "study_id": "bsdb:30778376",
        },
        26,
    ),
    (
        {
            "condition": "Waist circumference",
            "experiment_id": "bsdb:31339887-18",
            "group_0_name": "Low waist circumference (male and female)",
            "group_1_name": "High waist circumference (male and female)",
            "pmid": "31339887",
            "study_id": "bsdb:31339887",
        },
        26,
    ),
    (
        {
            "condition": "Hypertension",
            "experiment_id": "bsdb:33294219-1",
            "group_0_name": "Han Hypertensive at Low Altitude (LHH)",
            "group_1_name": "Tibetan Hypertensive at High Altitude (HTH)",
            "pmid": "33294219",
            "study_id": "bsdb:33294219",
        },
        26,
    ),
    (
        {
            "condition": "Population",
            "experiment_id": "bsdb:37596536-18",
            "group_0_name": "Himba",
            "group_1_name": "!Xun",
            "pmid": "37596536",
            "study_id": "bsdb:37596536",
        },
        26,
    ),
    (
        {
            "condition": "Population",
            "experiment_id": "bsdb:37596536-3",
            "group_0_name": "Kuvale",
            "group_1_name": "Tjimba",
            "pmid": "37596536",
            "study_id": "bsdb:37596536",
        },
        26,
    ),
    (
        {
            "condition": "Environmental exposure measurement",
            "experiment_id": "bsdb:40595471-2",
            "group_0_name": "Western Human Microbiome Project (HMP)",
            "group_1_name": "Adult Yanomami (16\u201350 yrs)",
            "pmid": "40595471",
            "study_id": "bsdb:40595471",
        },
        23,
    ),
    (
        {
            "condition": "Grey matter volume measurement",
            "experiment_id": "bsdb:41056270-6",
            "group_0_name": "Lower grey matter volume (GMV)",
            "group_1_name": "Higher grey matter volume (GMV)",
            "pmid": "41056270",
            "study_id": "bsdb:41056270",
        },
        26,
    ),
    (
        {
            "condition": "Age at assessment",
            "experiment_id": "bsdb:29196767-1",
            "group_0_name": "Young",
            "group_1_name": "Old",
            "pmid": "29196767",
            "study_id": "bsdb:29196767",
        },
        26,
    ),
    (
        {
            "condition": "Helminthiasis",
            "experiment_id": "bsdb:29486796-10",
            "group_0_name": "2008 - Uninfected",
            "group_1_name": "2010 - Uninfected",
            "pmid": "29486796",
            "study_id": "bsdb:29486796",
        },
        21,
    ),
    (
        {
            "condition": "Parkinson's disease",
            "experiment_id": "bsdb:29776865-2",
            "group_0_name": "PD duration less than 5 years",
            "group_1_name": "PD duration more than 5 years",
            "pmid": "29776865",
            "study_id": "bsdb:29776865",
        },
        1,
    ),
    (
        {
            "condition": "Obesity",
            "experiment_id": "bsdb:30459320-5",
            "group_0_name": "Lean Ghanaians(GL)",
            "group_1_name": "Lean Americans (UL)",
            "pmid": "30459320",
            "study_id": "bsdb:30459320",
        },
        26,
    ),
    (
        {
            "condition": "Air pollution",
            "experiment_id": "bsdb:31126003-3",
            "group_0_name": "moderately poluuted region (B)",
            "group_1_name": "heavily polluted region (C)",
            "pmid": "31126003",
            "study_id": "bsdb:31126003",
        },
        33,
    ),
    (
        {
            "condition": "Asthma",
            "experiment_id": "bsdb:31201890-3",
            "group_0_name": "Samples that were RV-A negative-",
            "group_1_name": "Samples with rhinovirus-A (RV-A)",
            "pmid": "31201890",
            "study_id": "bsdb:31201890",
        },
        36,
    ),
    (
        {
            "condition": "Small intestine enteropathy",
            "experiment_id": "bsdb:31674052-3",
            "group_0_name": "IBS.sibo-",
            "group_1_name": "IBS.sibo+",
            "pmid": "31674052",
            "study_id": "bsdb:31674052",
        },
        36,
    ),
    (
        {
            "condition": "Parkinson's disease",
            "experiment_id": "bsdb:33185912-3",
            "group_0_name": "Levodopa (LD)",
            "group_1_name": "Levodopa-carbidopa intestinal gel (LCIG)",
            "pmid": "33185912",
            "study_id": "bsdb:33185912",
        },
        33,
    ),
    (
        {
            "condition": "Non-alcoholic fatty liver disease",
            "experiment_id": "bsdb:38162648-3",
            "group_0_name": "HIV+ MASLD-",
            "group_1_name": "HIV- MASLD+",
            "pmid": "38162648",
            "study_id": "bsdb:38162648",
        },
        36,
    ),
    (
        {
            "condition": "Phenylketonuria",
            "experiment_id": "bsdb:38745393-10",
            "group_0_name": "Control Adults Saliva (Ctrl-Adults-S) patients",
            "group_1_name": "Control Adults Subgingival plaque (Ctrl-Adults-SP) patients",
            "pmid": "38745393",
            "study_id": "bsdb:38745393",
        },
        26,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "HER2 Positive Breast Carcinoma",
            "efo_id": "EFO:1",
            "group_0_name": "Her2 negative tumor",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Her2 positive tumor",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        3,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Crohn's disease",
            "efo_id": "EFO:1",
            "group_0_name": "crohn disease ASCA negative patients",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "crohn disease ASCA positive patients",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        3,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Type II diabetes mellitus",
            "description": "the interactive influence of metformin and opioids on gut "
            "microbiota in the subgroup with diabetes.",
            "efo_id": "EFO:1",
            "group_0_name": "Metaformin negative among Opoid negative",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Metaformin positive",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        32,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Lewy body dementia",
            "efo_id": "EFO:1",
            "group_0_name": "Parkinson\u2019s disease dementia negative (PDD-)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Dementia with Lewy Bodies (DLB) patients",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        3,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Asthma",
            "efo_id": "EFO:1",
            "group_0_name": "non-exacerbated asthma",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "exacerbated asthma",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        3,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Cancer",
            "efo_id": "EFO:1",
            "group_0_name": "Non-Cachectic Cancer",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Cachetic Cancer",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        3,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Intestinal obstruction",
            "efo_id": "EFO:1",
            "group_0_name": "Non-obstruction tumor group",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Obstruction tumor group",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        3,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Metastatic colorectal cancer",
            "efo_id": "EFO:1",
            "group_0_name": "No liver metastasis (NLM)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Liver metastasis (LM)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        3,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Depressive disorder",
            "efo_id": "EFO:1",
            "group_0_name": "Non-inflammatory depression",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Inflammatory depression",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        3,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Preeclampsia",
            "efo_id": "EFO:1",
            "group_0_name": "Mild/Non-severe Preeclampsia (P)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Severe Preeclampsia (SP)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        3,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Endometriosis",
            "efo_id": "EFO:1",
            "group_0_name": "Patients with non-endometrioma cysts",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Patients with ovarian endometrioma cysts",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        3,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Pregnancy",
            "efo_id": "EFO:1",
            "group_0_name": "non pregnant women without gingivitis",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "pregnant women without gingivitis",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        25,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Pregnancy",
            "efo_id": "EFO:1",
            "group_0_name": "Non-pregnant (FN)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Clinical pregnancy (FP)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        25,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Chronic constipation",
            "efo_id": "EFO:1",
            "group_0_name": "Control group",
            "group_0_sample_size": "10",
            "group_1_definition": "Experimental group of patients before and after AMT "
            "(Acupoint Massage Therapy)",
            "group_1_name": "Experimental group",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        29,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Exercise",
            "efo_id": "EFO:1",
            "group_0_name": "Control group",
            "group_0_sample_size": "10",
            "group_1_definition": "Tibetan Buddhist monks that performed meditation practices "
            "of Samatha and Vipassana for at least 2 hours a day for 3 to "
            "30 years.",
            "group_1_name": "Meditation group",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        29,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Cannabis use",
            "efo_id": "EFO:1",
            "group_0_name": "controls (lateral border of tongue)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "marijuana usage group",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        29,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Metabolic process",
            "efo_id": "EFO:1",
            "group_0_name": "Control group",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Trial group; 35d (Intestinal)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        29,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Irritable bowel syndrome",
            "efo_id": "EFO:1",
            "group_0_name": "HC-like IBS(Healthy Control like Irritable Bowel Syndrome)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "IBS1 (Irritable Bowel Syndrome 1)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        2,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Gastric adenocarcinoma,Epstein-Barr virus infection",
            "efo_id": "EFO:1",
            "group_0_name": "Epstein-Barr virus non-infected (EBV- GC)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Epstein-Barr virus infected (EBV+ GC)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        2,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Age",
            "efo_id": "EFO:1",
            "group_0_name": "2-12 years",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        22,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Chlamydia trachomatis",
            "efo_id": "EFO:1",
            "group_0_name": "Controls aged < = 0 years",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Controls aged > 10 years",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        34,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Coinfection,Age",
            "efo_id": "EFO:1",
            "group_0_name": "young HIV-infected group without coinfection",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "aging HIV-infected women with coinfections",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        2,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Gastric cancer",
            "efo_id": "EFO:1",
            "group_0_name": "Normal group",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Gastric Cancer",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        7,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Squamous cell carcinoma",
            "efo_id": "EFO:1",
            "group_0_name": "noncancer individuals",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Oral Squamous Cell Cancer patients",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        13,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "Healthy controls",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "HIV infected cART negative",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        7,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "COVID-19",
            "efo_id": "EFO:1",
            "group_0_name": "Negative for COVID-19 viral RNA in feces",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Positive for COVID-19 viral RNA in feces",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        13,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Human papilloma virus infection",
            "efo_id": "EFO:1",
            "group_0_name": "negative for intraepithelial lesion or malignancy",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "high-grade squamous intraepithelial lesions",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        3,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Non-alcoholic fatty liver disease",
            "efo_id": "EFO:1",
            "group_0_name": "obese without non-alcoholic fatty liver disease",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "non-alcoholic fatty liver disease patients",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        19,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Parkinson's disease",
            "efo_id": "EFO:1",
            "group_0_name": "healthy controls without any form of Parkinson's Disease",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "patients with idiopathic parkinson's disease",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        7,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Helminthiasis",
            "efo_id": "EFO:1",
            "group_0_name": "Uninfected (Zero STH infection)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Heavy Infected (Ascaris)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        7,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "HIV-uninfected men",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "HIV-infected men",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        7,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Health trait",
            "efo_id": "EFO:1",
            "group_0_name": "healthy controls",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "suboptimal health status",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        7,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Chlamydia trachomatis",
            "efo_id": "EFO:1",
            "group_0_name": "Chlamydia trichomonas Controls (t0)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Chlamydia trichomonas Cases (t0)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        10,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Chronic constipation",
            "efo_id": "EFO:1",
            "group_0_name": "Children without Functional constipation",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Children with Functional constipation",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        13,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Atrial fibrillation",
            "efo_id": "EFO:1",
            "group_0_name": "Individuals without Prevalent Atrial Fibrillation (AF)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Individuals with Prevalent Atrial Fibrillation (AF)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        13,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Experimental autoimmune encephalomyelitis",
            "efo_id": "EFO:1",
            "group_0_name": "Healthy",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "EAE (Experimental Autoimmune Encephalomyelitis)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        7,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Eczema",
            "efo_id": "EFO:1",
            "group_0_name": "Normal control",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "AD cases with no reported cream usage",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        7,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Type II diabetes mellitus",
            "efo_id": "EFO:1",
            "group_0_name": "Healthy controls",
            "group_0_sample_size": "10",
            "group_1_definition": "Individuals having Type 2 Diabetes Mellitus without Diabetic "
            "Retinopathy (T2D group); the T2DM cohort had history of "
            "taking anti-diabetic medication such as metformin",
            "group_1_name": "Individuals having Type 2 Diabetes Mellitus without Diabetic "
            "Retinopathy",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        35,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Colorectal cancer",
            "description": "Patients were sampled before chemotherapy or radiotherapy was started.",
            "efo_id": "EFO:1",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "CRC",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        7,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "HIV\u2013 individuals",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "HIV+ individuals",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        13,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "HIV-",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "HIV+",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        13,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "HIV-1-negative individuals",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Recently HIV-1 infected individuals",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        13,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Human immunodeficiency virus",
            "efo_id": "EFO:1",
            "group_0_name": "People Without HIV",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "People With HIV",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        13,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "Pregnant women ≥ 18 years of age without Human immunodeficiency "
            "virus (HIV) at Postpartum (PPt) visit",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Pregnant women ≥ 18 years of age with Human immunodeficiency virus "
            "(HIV) at Postpartum (PPt) visit.",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        12,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "uninfected controls",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "HIV-infected subjects",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        7,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "healthy control",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "HAART-treated HIV+ patients",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        16,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "Controls/ People living without HIV (PLWOH)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "People living with HIV on Antiretroviral Therapy (PLWH on ART)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        16,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "HIV-",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "People Living With HIV on Antiretroviral therapy (ART+ PLWH)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        17,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "HIV-uninfected MSM",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "ART-treated HIV infected MSM",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        16,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "HAART-naive HIV-1 patients",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "HAART-treated HIV-1 patients",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        31,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "Human immunodeficiency virus",
            "efo_id": "EFO:1",
            "group_0_name": "Acute HIV-infected participants at baseline (A0)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Acute HIV-infected participants on Antiretroviral Therapy (A12)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        32,
    ),
    (
        {
            "body_site": "Feces",
            "condition": "HIV mother to child transmission",
            "efo_id": "EFO:1",
            "group_0_name": "HIV-unexposed infants (HUU) at 6 weeks",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "HIV-exposed uninfected infants (HEU) at 6 weeks",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        28,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "Age-matched HIV-uninfected women",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "HIV-positive women who have been taking ART for more than 1 year",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        14,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "Human immunodeficiency virus negative (HIV-)Uninfected Controls",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Human immunodeficiency virus positive(HIV+) pre - antiretroviral "
            "therapy(ART)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        5,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "Human immunodeficiency virus positive(HIV+) pre- antiretroviral "
            "therapy (ART) in Fiebig stage 1-2(F1-2)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Human immunodeficiency virus positive (HIV+) pre- antiretroviral "
            "therapy(ART) in Fiebig stage 3(F3)",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        0,
    ),
    (
        {
            "body_site": "Feces",
            "efo_id": "EFO:1",
            "group_0_name": "Healthy controls (HC)",
            "group_0_sample_size": "10",
            "group_1_definition": "patients with colorectal cancer",
            "group_1_name": "Antiretroviral therapy (ART) Naïve at baseline",
            "group_1_sample_size": "10",
            "host_species": "Homo sapiens",
            "study_design": "case-control",
            "uberon_id": "UBERON:1",
        },
        5,
    ),
)

_GUARD_REFERENCE_SNAPSHOT: list[dict[str, Any]] = [
    {
        "record": {**_GUARD_REFERENCE_COMMON, **changes},
        "decision": list(_GUARD_REFERENCE_DECISIONS[index]),
    }
    for changes, index in _GUARD_REFERENCE_CASES
]


def _guard_cases() -> list[ExperimentRecord]:
    """The fixed cases the classifier-version guard classifies (#793).

    They are the repository's committed tables of experiment strings: the 60 experiments
    of the committed hand review, which records their condition and group names but not
    their design, Group 1 definition or description, and the #291, #569 and v6 regression
    tables above, each case built on BASE. They stay fixed so that a change in the
    decisions on them can only come from the classifier.
    """
    review = Path(__file__).resolve().parents[1] / "methods" / "contrast_orientation_review.tsv"
    with review.open(encoding="utf-8", newline="") as handle:
        reviewed = [
            ExperimentRecord(
                study_id=row["study_id"],
                experiment_id=row["experiment_id"],
                pmid=row["pmid"],
                doi="",
                study_design="",
                condition=row["condition"],
                efo_id="",
                body_site="",
                uberon_id="",
                host_species="",
                group_0_name=row["group_0_name"],
                group_1_name=row["group_1_name"],
                group_1_definition="",
                group_0_sample_size="",
                group_1_sample_size="",
                description="",
            )
            for row in csv.DictReader(handle, delimiter="\t")
        ]
    regression = [
        replace(BASE, **changes)
        for _, changes, _ in (
            *ELIGIBILITY_FIXTURES,
            *HIV_STATUS_FIXTURES,
            *V6_RECLASSIFIED_FIXTURES,
        )
    ]
    return reviewed + regression


def _guard_digests(cases: list[ExperimentRecord]) -> dict[str, str]:
    """SHA-256 of the cases, and of the classifier's decision on each.

    A decision is what a consumer can select on: the class, orientation, contrast type,
    control tier and flags. The evidence string explains a decision and is left out, so
    rewording it needs no new version.
    """
    inputs, decisions = [], []
    for record in cases:
        case = json.dumps(asdict(record), ensure_ascii=False, sort_keys=True)
        found = classify(record)
        flags = ",".join(sorted(flag for flag in found.flags.split(",") if flag))
        inputs.append(case)
        decisions.append(
            "\t".join(
                [
                    case,
                    found.contrast_class,
                    found.orientation,
                    found.contrast_type,
                    found.control_tier,
                    flags,
                ]
            )
        )
    return {
        "inputs": hashlib.sha256("\n".join(sorted(inputs)).encode()).hexdigest(),
        "decisions": hashlib.sha256("\n".join(sorted(decisions)).encode()).hexdigest(),
    }


_GuardKey = tuple[str, str, str, str, str]
_GuardDecision = tuple[str, str, str, str, str]
_GUARD_KEY_FIELDS = ("study_id", "experiment_id", "condition", "group_0_name", "group_1_name")
_GUARD_DECISION_FIELDS = ("contrast_class", "orientation", "contrast_type", "control_tier", "flags")


def _guard_reference_snapshot() -> list[dict[str, Any]]:
    """Return fresh containers while keeping the original record field set (#803)."""
    return [
        {"record": dict(row["record"]), "decision": list(row["decision"])}
        for row in _GUARD_REFERENCE_SNAPSHOT
    ]


def _guard_case_key(record: dict[str, Any]) -> _GuardKey:
    return cast(_GuardKey, tuple(record[field] for field in _GUARD_KEY_FIELDS))


def _guard_decision(record: ExperimentRecord) -> _GuardDecision:
    found = classify(record)
    flags = ",".join(sorted(flag for flag in found.flags.split(",") if flag))
    return (found.contrast_class, found.orientation, found.contrast_type, found.control_tier, flags)


def _guard_decision_rows(
    cases: list[ExperimentRecord],
) -> dict[_GuardKey, tuple[_GuardDecision, ...]]:
    """Group only decisions by the five-field key, retaining duplicate-key rows."""
    rows: dict[_GuardKey, list[_GuardDecision]] = {}
    for record in cases:
        rows.setdefault(_guard_case_key(asdict(record)), []).append(_guard_decision(record))
    return {key: tuple(sorted(values)) for key, values in rows.items()}


def _guard_reference_digests(reference: list[dict[str, Any]]) -> dict[str, str]:
    """Validate the frozen raw records and decisions without calling classify."""
    inputs, decisions = [], []
    for row in reference:
        case = json.dumps(row["record"], ensure_ascii=False, sort_keys=True)
        inputs.append(case)
        decisions.append("\t".join([case, *row["decision"]]))
    return {
        "inputs": hashlib.sha256("\n".join(sorted(inputs)).encode()).hexdigest(),
        "decisions": hashlib.sha256("\n".join(sorted(decisions)).encode()).hexdigest(),
    }


def _guard_case_label(key: _GuardKey) -> str:
    return (
        f"{key[0]} / {key[1]} (condition={key[2]!r}; "
        f"group_0_name={key[3]!r}; group_1_name={key[4]!r})"
    )


def _guard_input_differences(
    cases: list[ExperimentRecord], reference: list[dict[str, Any]]
) -> list[str]:
    """Describe changed raw input rows without treating multiplicity as a rule."""
    previous: dict[_GuardKey, Counter[str]] = {}
    current: dict[_GuardKey, Counter[str]] = {}
    for target, records in (
        (previous, [row["record"] for row in reference]),
        (current, [asdict(record) for record in cases]),
    ):
        for record in records:
            key = _guard_case_key(record)
            raw = json.dumps(record, ensure_ascii=False, sort_keys=True)
            target.setdefault(key, Counter())[raw] += 1
    differences = []
    for key in sorted(previous.keys() | current.keys()):
        old, new = previous.get(key, Counter()), current.get(key, Counter())
        removed, added = list((old - new).elements()), list((new - old).elements())
        paired = min(len(removed), len(added))
        if paired:
            changed_fields: set[str] = set()
            for old_raw, new_raw in zip(
                sorted(removed)[:paired], sorted(added)[:paired], strict=True
            ):
                before, after = json.loads(old_raw), json.loads(new_raw)
                changed_fields.update(
                    field
                    for field in before.keys() | after.keys()
                    if field not in before or field not in after or before[field] != after[field]
                )
            differences.append(
                f"edited/replaced {paired} variant(s) of {_guard_case_label(key)}; "
                f"differing fields: {', '.join(sorted(changed_fields))}"
            )
        if len(removed) > paired:
            differences.append(f"removed {len(removed) - paired} row(s): {_guard_case_label(key)}")
        if len(added) > paired:
            differences.append(f"added {len(added) - paired} row(s): {_guard_case_label(key)}")
    return differences


def _guard_changed_reference_decisions(reference: list[dict[str, Any]]) -> list[str]:
    """Classify every original record, even when evolving fixtures remove it."""
    key_counts = Counter(_guard_case_key(row["record"]) for row in reference)
    seen: Counter[_GuardKey] = Counter()
    record_fields = {field.name for field in fields(ExperimentRecord)}
    changed = []
    for row in reference:
        raw = row["record"]
        key = _guard_case_key(raw)
        seen[key] += 1
        # Preserve raw snapshot bytes above; reconstruct only for calling the
        # current API, whose unused/defaulted record fields may have evolved.
        record = ExperimentRecord(
            **{name: value for name, value in raw.items() if name in record_fields}
        )
        observed, expected = _guard_decision(record), tuple(row["decision"])
        if observed != expected:
            variant = f" variant {seen[key]}" if key_counts[key] > 1 else ""
            fields_changed = ", ".join(
                f"{name}: {old!r} -> {new!r}"
                for name, old, new in zip(_GUARD_DECISION_FIELDS, expected, observed, strict=True)
                if old != new
            )
            changed.append(f"{_guard_case_label(key)}{variant}: {fields_changed}")
    return changed


def _assert_classifier_guard(
    cases: list[ExperimentRecord],
    reference: list[dict[str, Any]],
    pinned: dict[str, str],
    *,
    reference_version: str | None = None,
) -> None:
    if reference_version is None:
        reference_version = _GUARD_REFERENCE_VERSION
    assert reference_version == CLASSIFIER_VERSION, (
        f"reference version {reference_version!r} differs from {CLASSIFIER_VERSION!r}; "
        "capture a reviewed reference for the approved classifier version."
    )
    assert _guard_reference_digests(reference) == pinned, (
        "the frozen reference does not match the pinned digests; restore the independently "
        "verified snapshot before diagnosing the current classifier."
    )
    found = _guard_digests(cases)
    input_changes = _guard_input_differences(cases, reference)
    decision_changes = _guard_changed_reference_decisions(reference)
    messages = []
    if input_changes:
        messages.append("fixture-input drift:\n" + "\n".join(input_changes))
    if decision_changes:
        messages.append(
            "classification decision changes on frozen guard cases:\n"
            + "\n".join(decision_changes)
            + "\nThis requires a CLASSIFIER_VERSION bump, a methods version-history entry, "
            "and new pins/reference under that version; never replace the existing pins."
        )
    elif input_changes:
        frozen_rows: dict[_GuardKey, set[_GuardDecision]] = {}
        for row in reference:
            frozen_rows.setdefault(_guard_case_key(row["record"]), set()).add(
                cast(_GuardDecision, tuple(row["decision"]))
            )
        current_rows = _guard_decision_rows(cases)
        shared = frozen_rows.keys() & current_rows.keys()
        changed_shared = [
            key for key in sorted(shared) if frozen_rows[key] != set(current_rows[key])
        ]
        messages.append(
            "no classifier decisions changed on the frozen guard cases. "
            f"Shared keys: {len(shared)}; "
            f"shared-key decision-value changes: {len(changed_shared)}. "
            "Added/removed keys have no shared counterpart; row multiplicity is fixture data. "
            "Restore the frozen fixture inputs and put new/edited cases in a separate test; "
            "this is not a CLASSIFIER_VERSION bump."
        )
        if changed_shared:
            messages.append(
                "Input-driven decision values differ at: "
                + "; ".join(_guard_case_label(key) for key in changed_shared)
            )
    assert not messages, "\n".join(messages)
    assert found == pinned, (
        "legacy guard digests changed without a diagnosed input or decision change"
    )


def test_the_classifier_s_decisions_change_only_with_its_version() -> None:
    """A classifier that decides differently is a new classifier version (#793).

    #739 changed how duration bands, a spaced "pre - <therapy>" and "(ART) Na\u00efve" are
    read, and moved four experiments out of ``treatment_exposure``, while the classifier
    still reported ``contrast-orientation-v5``. Run records, manifests and provenance then
    gave one version for two classifications.
    """
    cases = _guard_cases()
    found = _guard_digests(cases)
    pins = json.loads(CLASSIFIER_GUARD_DIGESTS.read_text(encoding="utf-8"))
    pinned = pins.get(CLASSIFIER_VERSION)
    assert pinned is not None, (
        f"no guard digests are pinned for {CLASSIFIER_VERSION}. Once the version history in "
        f"methods/contrast_orientation.md describes the change, pin {found} under "
        f"{CLASSIFIER_VERSION!r} in {CLASSIFIER_GUARD_DIGESTS.name} and capture its reviewed "
        "per-case reference. Existing version pins must remain unchanged."
    )
    _assert_classifier_guard(cases, _guard_reference_snapshot(), pinned)


def test_the_classifier_version_has_an_entry_in_the_version_history() -> None:
    """Each version says in the methods note what it changed and which experiments moved."""
    methods = Path(__file__).resolve().parents[1] / "methods" / "contrast_orientation.md"
    version = CLASSIFIER_VERSION.removeprefix("contrast-orientation-")
    heading = re.compile(
        rf"^## Changes from `contrast-orientation-v\d+` to `-{re.escape(version)}`",
        re.MULTILINE,
    )
    assert heading.search(methods.read_text(encoding="utf-8")), (
        f"methods/contrast_orientation.md has no '## Changes from `contrast-orientation-vN` "
        f"to `-{version}`' section for {CLASSIFIER_VERSION}"
    )


_UNIT_GUARD_DECISION = (
    "disease_vs_healthy",
    "native",
    "case_vs_control_subject",
    "tier_1_healthy_normal",
    "no_group_0_definition_in_export",
)


def _guard_api(name: str) -> Callable[..., Any]:
    function = globals().get(name)
    assert callable(function), f"The classifier guard must provide its {name} boundary"
    return cast(Callable[..., Any], function)


def _unit_guard_answer() -> Classification:
    return Classification(
        contrast_class=_UNIT_GUARD_DECISION[0],
        orientation=_UNIT_GUARD_DECISION[1],
        evidence="synthetic diagnostic fixture",
        flags=_UNIT_GUARD_DECISION[4],
        contrast_type=_UNIT_GUARD_DECISION[2],
        control_tier=_UNIT_GUARD_DECISION[3],
    )


def _unit_guard_data(
    monkeypatch: pytest.MonkeyPatch, cases: list[ExperimentRecord] | None = None
) -> tuple[list[ExperimentRecord], list[dict[str, Any]], dict[str, str]]:
    records = cases if cases is not None else [BASE]
    monkeypatch.setitem(globals(), "classify", lambda _record: _unit_guard_answer())
    reference = [
        {"record": asdict(record), "decision": list(_UNIT_GUARD_DECISION)} for record in records
    ]
    return records, reference, _guard_digests(records)


def _guard_failure_message(
    cases: list[ExperimentRecord],
    reference: list[dict[str, Any]],
    pinned: dict[str, str],
    *,
    reference_version: str | None = None,
) -> str:
    """Require refusal before checking its diagnosis, including in mutant proofs."""
    try:
        _guard_api("_assert_classifier_guard")(
            cases, reference, pinned, reference_version=reference_version
        )
    except AssertionError as error:
        return str(error)
    raise AssertionError("The classifier guard did not refuse the changed case")


@pytest.mark.parametrize("edit", ["metadata", "added", "removed", "key", "default-field"])
def test_guard_fixture_drift_does_not_request_a_classifier_version_bump(
    monkeypatch: pytest.MonkeyPatch, edit: str
) -> None:
    records, reference, pinned = _unit_guard_data(monkeypatch)
    if edit == "metadata":
        changed = [replace(records[0], pmid="corrected-pmid")]
        witness = "pmid"
    elif edit == "added":
        changed = [*records, replace(records[0], experiment_id="new-case")]
        witness = "added"
    elif edit == "removed":
        changed = []
        witness = "removed"
    elif edit == "key":
        changed = [replace(records[0], group_0_name="Corrected healthy group")]
        witness = "Corrected healthy group"
    else:
        extended = make_dataclass(
            "ExtendedGuardRecord",
            [("unused_extension", str, "")],
            bases=(ExperimentRecord,),
            frozen=True,
        )
        changed = [extended(**asdict(records[0]))]
        monkeypatch.setitem(globals(), "ExperimentRecord", extended)
        witness = "unused_extension"
    message = _guard_failure_message(changed, reference, pinned)
    assert "fixture-input drift" in message
    assert witness in message
    assert "no classifier decisions changed" in message
    assert "separate test" in message
    assert "not a CLASSIFIER_VERSION bump" in message


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("contrast_class", "other"),
        ("orientation", "reversed"),
        ("contrast_type", "disease_vs_disease"),
        ("control_tier", "tier_2_condition_negative_screened"),
        ("flags", "new_guard_flag"),
    ],
)
def test_guard_reports_each_changed_classification_field(
    monkeypatch: pytest.MonkeyPatch, field: str, value: str
) -> None:
    records, reference, pinned = _unit_guard_data(monkeypatch)
    changed_answer = replace(_unit_guard_answer(), **{field: value})
    monkeypatch.setitem(globals(), "classify", lambda _record: changed_answer)
    message = _guard_failure_message(records, reference, pinned)
    assert "classification decision changes" in message
    assert records[0].experiment_id in message
    assert field in message
    assert "requires a CLASSIFIER_VERSION bump" in message


def test_guard_keeps_decision_changes_visible_during_fixture_drift(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records, reference, pinned = _unit_guard_data(monkeypatch)
    changed_answer = replace(_unit_guard_answer(), orientation="reversed")
    monkeypatch.setitem(globals(), "classify", lambda _record: changed_answer)
    message = _guard_failure_message([replace(records[0], pmid="edited")], reference, pinned)
    assert "fixture-input drift" in message
    assert "classification decision changes" in message
    assert "orientation" in message
    assert "requires a CLASSIFIER_VERSION bump" in message


def test_guard_keeps_changed_rules_visible_after_all_fixture_rows_are_removed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records, reference, pinned = _unit_guard_data(monkeypatch)
    changed_answer = replace(_unit_guard_answer(), orientation="reversed")
    monkeypatch.setitem(globals(), "classify", lambda _record: changed_answer)
    message = _guard_failure_message([], reference, pinned)
    assert "fixture-input drift" in message
    assert "removed" in message
    assert "classification decision changes on frozen guard cases" in message
    assert records[0].experiment_id in message
    assert "orientation: 'native' -> 'reversed'" in message
    assert "requires a CLASSIFIER_VERSION bump" in message
    assert "no classifier decisions changed" not in message


def test_guard_distinguishes_input_driven_decisions_from_changed_rules(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records, reference, pinned = _unit_guard_data(monkeypatch)
    first = _unit_guard_answer()
    reversed_answer = replace(first, orientation="reversed")
    monkeypatch.setitem(
        globals(),
        "classify",
        lambda record: reversed_answer if record.description == "edited" else first,
    )
    message = _guard_failure_message([replace(records[0], description="edited")], reference, pinned)
    assert "fixture-input drift" in message
    assert "Input-driven decision values differ" in message
    assert records[0].experiment_id in message
    assert "no classifier decisions changed on the frozen guard cases" in message
    assert "not a CLASSIFIER_VERSION bump" in message


def test_guard_decision_rows_preserve_duplicate_keys_and_ignore_unused_fields(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records, _reference, _pinned = _unit_guard_data(monkeypatch)
    answer = replace(_unit_guard_answer(), flags="beta,,alpha")
    monkeypatch.setitem(globals(), "classify", lambda _record: answer)
    rows = _guard_api("_guard_decision_rows")([records[0], replace(records[0], pmid="other")])
    key = (
        records[0].study_id,
        records[0].experiment_id,
        records[0].condition,
        records[0].group_0_name,
        records[0].group_1_name,
    )
    expected = (*_UNIT_GUARD_DECISION[:4], "alpha,beta")
    assert rows == {key: (expected, expected)}


def test_guard_reordering_and_evidence_wording_do_not_change_the_pins(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records, reference, pinned = _unit_guard_data(
        monkeypatch, [BASE, replace(BASE, experiment_id="second-case")]
    )
    answer = replace(_unit_guard_answer(), evidence="changed explanation only")
    monkeypatch.setitem(globals(), "classify", lambda _record: answer)
    _guard_api("_assert_classifier_guard")(list(reversed(records)), reference, pinned)


def test_guard_removing_an_identical_duplicate_is_only_fixture_drift(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records, reference, pinned = _unit_guard_data(monkeypatch, [BASE, BASE])
    message = _guard_failure_message(records[:1], reference, pinned)
    assert "removed" in message
    assert "no classifier decisions changed" in message
    assert "not a CLASSIFIER_VERSION bump" in message


def test_guard_does_not_hide_swapped_decisions_under_one_five_field_key(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [replace(BASE, description="first"), replace(BASE, description="second")]
    first = _unit_guard_answer()
    second = replace(first, orientation="reversed")
    monkeypatch.setitem(
        globals(), "classify", lambda record: first if record.description == "first" else second
    )
    reference = [
        {"record": asdict(records[0]), "decision": list(_UNIT_GUARD_DECISION)},
        {
            "record": asdict(records[1]),
            "decision": [*_UNIT_GUARD_DECISION[:1], "reversed", *_UNIT_GUARD_DECISION[2:]],
        },
    ]
    pinned = _guard_digests(records)
    monkeypatch.setitem(
        globals(), "classify", lambda record: second if record.description == "first" else first
    )
    message = _guard_failure_message(records, reference, pinned)
    assert "classification decision changes" in message
    assert "variant 1" in message
    assert "variant 2" in message
    assert "orientation" in message


def test_guard_rejects_a_corrupt_reference_before_diagnosing_current_rules(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records, reference, pinned = _unit_guard_data(monkeypatch)
    reference[0]["decision"][0] = "corrupted"
    message = _guard_failure_message(records, reference, pinned)
    assert "frozen reference" in message
    assert "pinned digests" in message


def test_guard_rejects_a_reference_from_a_different_classifier_version(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records, reference, pinned = _unit_guard_data(monkeypatch)
    message = _guard_failure_message(
        records, reference, pinned, reference_version="contrast-orientation-v0"
    )
    assert "reference version" in message
