"""Contrast classes and orientation for BugSigDB experiments (#188)."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from bsdbng.contrast_orientation import (
    DISEASE_VS_DISEASE,
    DISEASE_VS_HEALTHY,
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
    assert _classify(group_0_name="HPV-/HIV-", group_1_name="HPV+/HIV-") == (
        DISEASE_VS_HEALTHY,
        NATIVE,
    )
    # A hyphen used as punctuation is not a negative test.
    assert _classify(group_0_name="NAFLD- [fatty liver]", group_1_name="NASH") != (
        DISEASE_VS_HEALTHY,
        NATIVE,
    )


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
