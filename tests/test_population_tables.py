"""Explicit population selection preserves supplied tables and direction identity (#1344)."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from mtp import contrast_orientation as orientation
from mtp import population_tables as population

TIER12 = "disease_vs_healthy_tier1_2_subject_controls"
NATIVE = "group_1_vs_group_0"
ORIENTED = "disease_vs_comparison"
KEYS = ["study_id", "experiment_id", "inc_dec"]


@dataclass
class Inputs:
    profiles: Path
    experiments: Path
    raw: Path
    frame: pd.DataFrame
    raw_rows: list[dict[str, str]]

    def selection(self, name: str = TIER12, basis: str = NATIVE) -> population.PopulationSelection:
        return population.get_population_selection(
            self.profiles,
            self.experiments,
            contrast_population=name,
            input_direction_basis=basis,
            raw_dump_path=self.raw,
        )

    def write_raw(self) -> None:
        with self.raw.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(self.raw_rows[0]))
            writer.writeheader()
            writer.writerows(self.raw_rows)


@pytest.fixture
def inputs(tmp_path: Path) -> Inputs:
    examples = [
        (1, "Asthma", "Healthy controls", "Asthma patients", "Homo sapiens", "16S", True),
        (2, "Asthma", "Asthma patients", "Healthy controls", "Homo sapiens", "WMS", True),
        (3, "HIV infection", "HIV negative", "HIV positive", "Homo sapiens", "16S", True),
        (
            4,
            "HIV infection",
            "Controls",
            "HIV Patients treated with Maraviroc",
            "Homo sapiens",
            "16S",
            True,
        ),
        (5, "Asthma", "Healthy controls", "Asthma patients", "Mus musculus", "16S", True),
        (6, "Asthma", "Healthy controls", "Asthma patients", "Homo sapiens", "ITS / ITS2", True),
        (7, "Asthma", "Healthy controls", "Asthma patients", "Homo sapiens", "16S", False),
    ]
    rows, metadata, raw = [], [], []
    for number, condition, group0, group1, host, sequencing, paired in examples:
        key = ("bsdb:1", f"bsdb:1-{number}")
        metadata.append(
            dict(
                zip(
                    ["study_id", "experiment_id", "host_species", "sequencing_type"],
                    [*key, host, sequencing],
                    strict=True,
                )
            )
        )
        for direction in ["inc", "dec"] if paired else ["inc"]:
            rows.append([*key, direction, number * 10 + (1 if direction == "inc" else 2)])
        raw.append(
            {
                "BSDB ID": f"bsdb:1/{number}/1",
                "PMID": "",
                "DOI": "",
                "Study design": "case-control",
                "Condition": condition,
                "EFO ID": "",
                "Body site": "Feces",
                "UBERON ID": "",
                "Host species": host,
                "Group 0 name": group0,
                "Group 1 name": group1,
                "Group 1 definition": group1,
                "Group 0 sample size": "12",
                "Group 1 sample size": "15",
                "Description": "",
            }
        )
    frame = pd.DataFrame(rows, columns=[*KEYS, "trait"])
    result = Inputs(
        tmp_path / "profiles.tsv", tmp_path / "experiments.tsv", tmp_path / "raw.csv", frame, raw
    )
    frame.to_csv(result.profiles, sep="\t", index=False)
    pd.DataFrame(metadata).to_csv(result.experiments, sep="\t", index=False)
    result.write_raw()
    return result


def test_subset_uses_host_method_pairing_control_tier_and_native_direction(inputs: Inputs) -> None:
    selection = inputs.selection()
    assert selection.selected_keys == frozenset(("bsdb:1", f"bsdb:1-{i}") for i in [1, 2, 3])
    assert selection.keys_to_flip == frozenset({("bsdb:1", "bsdb:1-2")})
    assert selection.input_direction_basis == NATIVE
    assert selection.output_direction_basis == ORIENTED
    before = inputs.frame.copy(deep=True)
    actual = population.select_profiles_pandas(inputs.frame, selection)
    assert actual["inc_dec"].tolist() == ["inc", "dec", "dec", "inc", "inc", "dec"]
    assert actual["trait"].tolist() == [11, 12, 21, 22, 31, 32]
    assert actual.index.tolist() == [0, 1, 2, 3, 4, 5]
    assert actual.dtypes.equals(inputs.frame.dtypes)
    assert_frame_equal(inputs.frame, before)


def test_other_declared_subsets_keep_tier_three(inputs: Inputs) -> None:
    for name in ("disease_vs_healthy_oriented", "case_vs_control_subject"):
        assert inputs.selection(name).selected_keys == frozenset(
            ("bsdb:1", f"bsdb:1-{i}") for i in [1, 2, 3, 4]
        )


def test_supplied_disease_basis_never_flips_again(inputs: Inputs) -> None:
    oriented = population.select_profiles_pandas(inputs.frame, inputs.selection())
    selection = population.selection_from_verified_core(
        oriented, contrast_population=TIER12, direction_basis=ORIENTED
    )
    assert not selection.keys_to_flip
    assert_frame_equal(population.select_profiles_pandas(oriented, selection), oriented)
    assert not inputs.selection(basis=ORIENTED).keys_to_flip


def test_verified_core_constructor_does_not_read_or_classify(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame([("s1", "e1", "inc", 2), ("s2", "e1", "dec", 5)], columns=[*KEYS, "trait"])

    def unexpected(*_args: object, **_kwargs: object) -> None:
        pytest.fail("Verified-core selection must not reclassify or read external files")

    monkeypatch.setattr(orientation, "read_experiment_records", unexpected)
    monkeypatch.setattr(orientation, "classify", unexpected)
    selection = population.selection_from_verified_core(
        frame, contrast_population=TIER12, direction_basis=ORIENTED
    )
    assert selection.selected_keys == frozenset({("s1", "e1"), ("s2", "e1")})
    assert_frame_equal(population.select_profiles_pandas(frame, selection), frame)


def test_all_paired_is_identity_including_unpaired_nonhuman_and_extra_taxa(inputs: Inputs) -> None:
    selection = population.get_population_selection(
        inputs.profiles,
        Path("unused-experiments.tsv"),
        contrast_population="all_paired",
        input_direction_basis=NATIVE,
    )
    assert not selection.keys_to_flip
    assert selection.output_direction_basis == NATIVE
    assert_frame_equal(population.select_profiles_pandas(inputs.frame, selection), inputs.frame)
    taxa = pd.DataFrame(
        [("other", "e", "increased", "t1"), ("bsdb:1", "bsdb:1-7", "decreased", "t2")],
        columns=["study_id", "experiment_id", "direction", "taxon_id"],
    )
    assert_frame_equal(population.select_taxa_pandas(taxa, selection), taxa)


def test_taxon_occurrences_preserved_and_both_direction_columns_oriented(inputs: Inputs) -> None:
    taxa = pd.DataFrame(
        [
            ("bsdb:1", "bsdb:1-2", "increased", "inc", "t1"),
            ("bsdb:1", "bsdb:1-2", "increased", "inc", "t1"),
            ("bsdb:1", "bsdb:1-1", "decreased", "dec", "t2"),
            ("bsdb:1", "bsdb:1-5", "increased", "inc", "t3"),
        ],
        columns=["study_id", "experiment_id", "direction", "inc_dec", "taxon_id"],
        index=[9, 3, 8, 6],
    )
    before = taxa.copy(deep=True)
    actual = population.select_taxa_pandas(taxa, inputs.selection())
    assert actual.index.tolist() == [9, 3, 8]
    assert actual["direction"].tolist() == ["decreased", "decreased", "decreased"]
    assert actual["inc_dec"].tolist() == ["dec", "dec", "dec"]
    assert actual["taxon_id"].tolist() == ["t1", "t1", "t2"]
    assert actual.dtypes.equals(taxa.dtypes)
    assert_frame_equal(taxa, before)


def test_replaced_files_at_same_path_are_not_cached(inputs: Inputs) -> None:
    assert ("bsdb:1", "bsdb:1-2") in inputs.selection().keys_to_flip
    inputs.raw_rows[1]["Group 0 name"] = "Healthy controls"
    inputs.raw_rows[1]["Group 1 name"] = "Asthma patients"
    inputs.raw_rows[1]["Group 1 definition"] = "Asthma patients"
    inputs.write_raw()
    assert not inputs.selection().keys_to_flip
    inputs.frame.loc[inputs.frame["experiment_id"] != "bsdb:1-3"].to_csv(
        inputs.profiles, sep="\t", index=False
    )
    assert ("bsdb:1", "bsdb:1-3") not in inputs.selection().selected_keys


@pytest.mark.parametrize("missing", ["experiment", "classification"])
def test_profile_orphans_are_refused(inputs: Inputs, missing: str) -> None:
    if missing == "experiment":
        metadata = pd.read_csv(inputs.experiments, sep="\t", dtype=str)
        metadata.loc[metadata["experiment_id"] != "bsdb:1-1"].to_csv(
            inputs.experiments, sep="\t", index=False
        )
    else:
        inputs.raw_rows = inputs.raw_rows[1:]
        inputs.write_raw()
    with pytest.raises(ValueError, match=r"Missing.*(metadata|classification)"):
        inputs.selection()


@pytest.mark.parametrize(
    ("population_name", "basis"),
    [("unknown", NATIVE), ("all_paired", ORIENTED), (TIER12, NATIVE), (TIER12, "as_reported")],
)
def test_verified_core_rejects_invalid_population_basis(
    inputs: Inputs, population_name: str, basis: str
) -> None:
    with pytest.raises(ValueError, match=r"population|basis"):
        population.selection_from_verified_core(
            inputs.frame, contrast_population=population_name, direction_basis=basis
        )


@pytest.mark.parametrize("mutation", ["blank", "null", "duplicate", "direction"])
def test_explicit_profile_selection_rejects_invalid_keys(inputs: Inputs, mutation: str) -> None:
    frame = inputs.frame.copy()
    if mutation == "blank":
        frame.loc[0, "study_id"] = ""
    elif mutation == "null":
        frame.loc[0, "study_id"] = None
    elif mutation == "duplicate":
        frame = pd.concat([frame, frame.iloc[:1]])
    else:
        frame.loc[0, "inc_dec"] = "up"
    with pytest.raises(ValueError, match=r"identity|direction|Duplicate"):
        population.selection_from_verified_core(
            frame, contrast_population=TIER12, direction_basis=ORIENTED
        )


def test_taxa_conflicting_direction_columns_are_refused(inputs: Inputs) -> None:
    frame = pd.DataFrame(
        [("bsdb:1", "bsdb:1-2", "increased", "dec")],
        columns=["study_id", "experiment_id", "direction", "inc_dec"],
    )
    with pytest.raises(ValueError, match=r"direction.*disagree"):
        population.select_taxa_pandas(frame, inputs.selection())


def test_explicit_subset_requires_raw_source(inputs: Inputs) -> None:
    with pytest.raises(ValueError, match="raw"):
        population.get_population_selection(
            inputs.profiles,
            inputs.experiments,
            contrast_population=TIER12,
            input_direction_basis=NATIVE,
        )


def test_composite_selection_does_not_admit_another_study() -> None:
    frame = pd.DataFrame(
        [
            ("s1", "e1", "inc", 2),
            ("s1", "e1", "dec", 1),
            ("s2", "e1", "inc", 8),
            ("s2", "e1", "dec", 5),
        ],
        columns=[*KEYS, "trait"],
    )
    selection = population.selection_from_verified_core(
        frame.iloc[:2],
        contrast_population=TIER12,
        direction_basis=ORIENTED,
    )
    assert_frame_equal(population.select_profiles_pandas(frame, selection), frame.iloc[:2])


def test_explicit_selection_preserves_string_dtype_and_ignores_count_values(inputs: Inputs) -> None:
    frame = inputs.frame.astype(
        {"study_id": "string", "experiment_id": "string", "inc_dec": "string", "trait": "Int64"}
    )
    actual = population.select_profiles_pandas(frame, inputs.selection())
    assert actual.dtypes.equals(frame.dtypes)
    assert actual["trait"].tolist() == frame.iloc[:6]["trait"].tolist()
