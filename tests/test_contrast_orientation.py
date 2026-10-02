"""Contrast classes and orientation for BugSigDB experiments (#188)."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from bsdbng.contrast_orientation import (
    ARM_GROUP,
    DISEASE_VS_DISEASE,
    DISEASE_VS_HEALTHY,
    INTERVENTION_EVIDENCE,
    LONGITUDINAL,
    NATIVE,
    NOT_APPLICABLE,
    OTHER,
    REVERSED,
    TREATMENT,
    UNRESOLVED,
    Classification,
    ExperimentRecord,
    classify,
    oriented_profiles,
    read_experiment_records,
    write_manifest,
)
from bsdbng.trait_guilds import Profile

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
