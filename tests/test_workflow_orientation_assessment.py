"""Selected-core orientation evidence stays bound and is never flipped twice (#658)."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pytest

from mtp.contrast_orientation import CLASSIFIER_VERSION, REVERSED
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_guilds import TABLE1_TRAITS
from mtp.workflow.arms.noel import build_commands
from mtp.workflow.arms.types import ArmContext
from mtp.workflow.operations import OPERATIONS, validate_execution_product
from mtp.workflow.orientation_assessment import main, write_assessment
from mtp.workflow.runtime import artifact_path, execute_workflow
from mtp.workflow.verification import verify_core_run

POPULATION = "disease_vs_healthy_tier1_2_subject_controls"
TRAITS = [terms[0] for terms in TABLE1_TRAITS.values()]


def _table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _write_raw(path: Path) -> None:
    """Keep native, reversed, excluded, and zero-annotation experiments distinct."""
    fields = [
        "BSDB ID",
        "Study",
        "Title",
        "Year",
        "PMID",
        "DOI",
        "Study design",
        "Host species",
        "Sequencing type",
        "Body site",
        "UBERON ID",
        "Condition",
        "EFO ID",
        "Group 0 name",
        "Group 1 name",
        "Group 1 definition",
        "Group 0 sample size",
        "Group 1 sample size",
        "Description",
        "Abundance in Group 1",
        "NCBI Taxonomy IDs",
        "MetaPhlAn taxon names",
    ]
    with path.open("w", encoding="utf-8", newline="") as stream:
        stream.write("# BugSigDB synthetic-orientation-assessment, License: CC0\n")
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for study in range(1, 5):
            group_0, group_1 = "Healthy controls", "CRC patients"
            if study == 2:
                group_0, group_1 = group_1, group_0
            elif study == 3:
                group_0, group_1 = "Placebo", "Probiotic"
            for signature, direction in enumerate(("increased", "decreased"), start=1):
                taxon = 100 + signature
                if study == 2:
                    taxon = 103 - signature
                elif study == 4:
                    taxon = 300
                name = {
                    101: "s__Escherichia_coli",
                    102: "s__Bacillus_subtilis",
                    300: "s__Bacteroides_fragilis",
                }[taxon]
                writer.writerow(
                    {
                        "BSDB ID": f"bsdb:{study}/1/{signature}",
                        "Study": f"Study {study}",
                        "Title": "Synthetic orientation assessment",
                        "Year": "2026",
                        "PMID": str(study),
                        "DOI": "",
                        "Study design": "case-control",
                        "Host species": "Homo sapiens",
                        "Sequencing type": "16S",
                        "Body site": "Feces",
                        "UBERON ID": "UBERON:0001988",
                        "Condition": "Colorectal cancer",
                        "EFO ID": "EFO:0005842",
                        "Group 0 name": group_0,
                        "Group 1 name": group_1,
                        "Group 1 definition": group_1,
                        "Group 0 sample size": "10",
                        "Group 1 sample size": "10",
                        "Description": "",
                        "Abundance in Group 1": direction,
                        "NCBI Taxonomy IDs": str(taxon),
                        "MetaPhlAn taxon names": name,
                    }
                )


def _write_kg(directory: Path) -> tuple[Path, Path]:
    """All six traits occur in the KG; the third reported taxon has no annotations."""
    edges = directory / "merged-kg_edges.tsv"
    nodes = directory / "merged-kg_nodes.tsv"
    with edges.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["id", "subject", "predicate", "object", "primary_knowledge_source"])
        for index, trait in enumerate(TRAITS):
            writer.writerow(
                [
                    f"annotation{index}",
                    "NCBITaxon:102" if index == 0 else "NCBITaxon:101",
                    "biolink:has_phenotype" if trait.startswith("METPO:") else "METPO:2000103",
                    trait,
                    "fixture:annotation",
                ]
            )
    with nodes.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["id", "category", "name"])
        for identifier, name in (
            ("NCBITaxon:101", "Escherichia coli"),
            ("NCBITaxon:102", "Bacillus subtilis"),
            ("NCBITaxon:300", "Bacteroides fragilis"),
        ):
            writer.writerow([identifier, "biolink:OrganismTaxon", name])
        for index, trait in enumerate(TRAITS):
            writer.writerow([trait, "biolink:PhenotypicQuality", f"Fixture trait {index}"])
        writer.writerow(["METPO:2000103", "biolink:OntologyClass", "Fixture trait relationship"])
    return edges, nodes


@pytest.fixture
def inputs(tmp_path: Path) -> dict[str, Path]:
    raw = tmp_path / "full_dump.csv"
    _write_raw(raw)
    edges, nodes = _write_kg(tmp_path)
    run = execute_workflow(
        REPO_ROOT / "workflows/bugsigdb_core_subset.yaml",
        {"raw_dump": raw, "kg_edges": edges, "kg_nodes": nodes},
        tmp_path / "core",
        OPERATIONS,
        parameter_overrides={"source_mode": "fixture", "contrast_population": POPULATION},
        validator=validate_execution_product,
        repository_root=REPO_ROOT,
    )
    review = tmp_path / "review.tsv"
    review.write_text(
        "study_id\texperiment_id\tsampled_class\tsampled_orientation\t"
        "verdict\tnote\treviewed_class\treviewed_orientation\n"
        "bsdb:1\tbsdb:1-1\tdisease_vs_healthy\tnative\tagree\t\tdisease_vs_healthy\tnative\n"
        "bsdb:3\tbsdb:3-1\ttreatment_exposure\tnot_applicable\t"
        "disagree\tfixture disagreement\tdisease_vs_healthy\tnative\n",
        encoding="utf-8",
    )
    return {
        "core_run": Path(run.manifest_path),
        "profile_input": artifact_path(run, "count_profiles"),
        "experiments_input": artifact_path(run, "experiments"),
        "full_dump_input": artifact_path(run, "raw_dump"),
        "studies_input": artifact_path(run, "studies"),
        "source_studies_input": artifact_path(run, "source_studies"),
        "selection_input": artifact_path(run, "contrast_selection"),
        "output_dir": tmp_path / "assessment",
        "review_input": review,
    }


def _run(inputs: dict[str, Path], population: str = POPULATION) -> Path:
    return write_assessment(
        core_run=inputs["core_run"],
        profile_input=inputs["profile_input"],
        experiments_input=inputs["experiments_input"],
        full_dump_input=inputs["full_dump_input"],
        studies_input=inputs["studies_input"],
        source_studies_input=inputs["source_studies_input"],
        selection_input=inputs["selection_input"],
        output_dir=inputs["output_dir"],
        review_input=inputs["review_input"],
        population=population,
        code_commit="a" * 40,
    )


def test_original_reversal_is_reported_but_selected_profiles_are_not_flipped(
    inputs: dict[str, Path],
) -> None:
    path = _run(inputs)
    provenance = json.loads(path.read_text(encoding="utf-8"))
    manifest = _table(path.parent / "contrast_manifest.tsv")
    reversed_row = next(row for row in manifest if row["study_id"] == "bsdb:2")
    assert reversed_row["orientation"] == REVERSED
    assert reversed_row["group_0_name"] == "CRC patients"
    assert reversed_row["group_1_name"] == "Healthy controls"
    assert reversed_row["in_analysed_corpus"] == "true"
    comparison = _table(path.parent / "disease_control_comparison.tsv")
    fraction = next(
        row
        for row in comparison
        if row["comparison"] == "fraction" and row["numerator"] == "facultatively_aerobic"
    )
    # Both nonzero studies are negative in already-oriented input. A second flip
    # of the reversed study would instead report one positive and one negative.
    assert fraction["studies_up"] == "0"
    assert fraction["studies_down"] == "2"
    assert fraction["studies_scored"] == "2"
    assert fraction["population_experiments"] == "3"
    assert {row["population"] for row in comparison} == {POPULATION}
    assert provenance["classifier_version"] == CLASSIFIER_VERSION
    assert provenance["profile_direction_basis"] == "disease_vs_comparison"
    assert "disease group" in provenance["direction_convention"]
    assert "Group 1" in provenance["raw_direction_convention"]
    assert provenance["profile_directions_relabelled_by_assessment"] is False
    assert provenance["analysed_experiments"] == 3
    assert provenance["selection_evidence"]["reversed_experiments"] == [
        {"study_id": "bsdb:2", "experiment_id": "bsdb:2-1"}
    ]


def test_actual_input_identity_and_only_in_population_review_are_retained(
    inputs: dict[str, Path],
) -> None:
    path = _run(inputs)
    record = json.loads(path.read_text(encoding="utf-8"))
    for name, field in (
        ("profile_input", "profile_sha256"),
        ("experiments_input", "experiments_sha256"),
        ("full_dump_input", "full_dump_sha256"),
    ):
        assert record[field] == file_sha256(inputs[name])
    assert record["selection_evidence"]["sha256"] == file_sha256(inputs["selection_input"])
    assert record["core_execution"]["manifest_sha256"] == file_sha256(inputs["core_run"])
    assert record["core_execution"]["contrast_population"] == POPULATION
    assert record["core_execution"]["direction_basis"] == "disease_vs_comparison"
    assert set(record["core_execution"]["bound_products"]) == {
        "count_profiles",
        "experiments",
        "raw_dump",
        "studies",
        "source_studies",
        "contrast_selection",
    }
    review = record["review"]
    assert review["sha256"] == file_sha256(inputs["review_input"])
    assert review["reviewed_experiments"] == review["agreed"] == 1
    assert review["disagreed"] == 0
    assert review["reviewed_experiments_outside_population"] == 1
    assert review["not_reviewed_experiments"] == 2
    assert "checksum agreement do not establish biological validity" in " ".join(
        record["limitations"]
    )


def test_missing_hand_review_stays_explicitly_unreviewed(inputs: dict[str, Path]) -> None:
    inputs["review_input"].unlink()
    record = json.loads(_run(inputs).read_text(encoding="utf-8"))
    assert record["review"]["sha256"] is None
    assert record["review"]["reviewed_experiments"] == 0
    assert record["review"]["not_reviewed_experiments"] == 3
    manifest = _table(inputs["output_dir"] / "contrast_manifest.tsv")
    assert {row["validation_status"] for row in manifest} == {"not_reviewed"}


def test_missing_profile_metadata_is_not_silently_dropped(inputs: dict[str, Path]) -> None:
    path = inputs["experiments_input"]
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text(
        "\n".join(line for line in lines if not line.startswith("bsdb:2\t")) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="content or inventory changed"):
        _run(inputs)
    assert not inputs["output_dir"].exists()


def test_outside_experiment_cannot_be_assessed_as_selected(inputs: dict[str, Path]) -> None:
    metadata = inputs["experiments_input"]
    original = metadata.read_text(encoding="utf-8")
    row = next(line for line in original.splitlines() if line.startswith("bsdb:1\t"))
    metadata.write_text(original + row.replace("bsdb:1", "bsdb:3") + "\n", encoding="utf-8")
    profiles = inputs["profile_input"]
    original = profiles.read_text(encoding="utf-8")
    added = [
        line.replace("bsdb:1", "bsdb:3")
        for line in original.splitlines()
        if line.startswith("bsdb:1\t")
    ]
    profiles.write_text(original + "\n".join(added) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="content or inventory changed"):
        _run(inputs)
    assert not inputs["output_dir"].exists()


def test_missing_selected_pair_is_not_silently_omitted(inputs: dict[str, Path]) -> None:
    path = inputs["profile_input"]
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text(
        "\n".join(line for line in lines if not line.startswith("bsdb:2\t")) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="content or inventory changed"):
        _run(inputs)
    assert not inputs["output_dir"].exists()


def test_duplicate_profile_cannot_double_weight_a_selected_direction(
    inputs: dict[str, Path],
) -> None:
    path = inputs["profile_input"]
    original = path.read_text(encoding="utf-8")
    repeated = next(line for line in original.splitlines() if line.startswith("bsdb:2\t"))
    path.write_text(original + repeated + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="content or inventory changed"):
        _run(inputs)
    assert not inputs["output_dir"].exists()


@pytest.mark.parametrize("field", ["raw_dump_sha256", "output_direction_basis", "population"])
def test_tampered_selection_evidence_fails_before_publication(
    inputs: dict[str, Path], field: str
) -> None:
    path = inputs["selection_input"]
    selection = json.loads(path.read_text(encoding="utf-8"))
    selection[field] = "incorrect"
    path.write_text(json.dumps(selection), encoding="utf-8")
    with pytest.raises(ValueError, match="content or inventory changed"):
        _run(inputs)
    assert not inputs["output_dir"].exists()


@pytest.mark.parametrize("population", ["all_paired", "unknown"])
def test_adapter_refuses_native_or_unknown_direction_basis(
    inputs: dict[str, Path], population: str
) -> None:
    with pytest.raises(ValueError, match="requires a disease-oriented population"):
        _run(inputs, population)
    assert not inputs["output_dir"].exists()


def test_prior_assessment_is_preserved(inputs: dict[str, Path]) -> None:
    output = inputs["output_dir"]
    output.mkdir()
    sentinel = output / "previous.txt"
    sentinel.write_text("previous assessment\n", encoding="utf-8")
    with pytest.raises(FileExistsError):
        _run(inputs)
    assert sentinel.read_text(encoding="utf-8") == "previous assessment\n"


def test_empty_precreated_output_directory_is_accepted(inputs: dict[str, Path]) -> None:
    output = inputs["output_dir"]
    output.mkdir()
    provenance = _run(inputs)
    assert provenance.is_file()
    assert (output / "contrast_manifest.tsv").is_file()


@pytest.mark.parametrize("overwrite_core", [False, True])
def test_native_profiles_with_identical_keys_cannot_claim_core_orientation(
    inputs: dict[str, Path], tmp_path: Path, overwrite_core: bool
) -> None:
    source_rows = _table(inputs["profile_input"])
    native_rows = [
        {
            **row,
            "inc_dec": (
                {"inc": "dec", "dec": "inc"}[row["inc_dec"]]
                if row["study_id"] == "bsdb:2"
                else row["inc_dec"]
            ),
        }
        for row in source_rows
    ]
    key = ("study_id", "experiment_id", "inc_dec")
    assert {tuple(row[field] for field in key) for row in native_rows} == {
        tuple(row[field] for field in key) for row in source_rows
    }
    assert native_rows != source_rows
    native = inputs["profile_input"] if overwrite_core else tmp_path / "native_profiles.tsv"
    with native.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(source_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(native_rows)
    changed = {**inputs, "profile_input": native}
    diagnostic = (
        "content or inventory changed"
        if overwrite_core
        else "core-run conflicts with explicit --profile-input"
    )
    with pytest.raises(ValueError, match=diagnostic):
        _run(changed)
    assert not inputs["output_dir"].exists()


@pytest.mark.parametrize(
    "argument",
    [
        "profile_input",
        "experiments_input",
        "full_dump_input",
        "studies_input",
        "source_studies_input",
        "selection_input",
    ],
)
def test_each_explicit_input_must_be_the_declared_core_product(
    inputs: dict[str, Path], tmp_path: Path, argument: str
) -> None:
    source = inputs[argument]
    other = tmp_path / f"unbound-{source.name}"
    if source.is_dir():
        shutil.copytree(source, other)
    else:
        shutil.copyfile(source, other)
    changed = {**inputs, argument: other}
    with pytest.raises(
        ValueError, match=f"core-run conflicts with explicit --{argument.replace('_', '-')}"
    ):
        _run(changed)
    assert not inputs["output_dir"].exists()


def test_population_must_match_the_verified_execution(inputs: dict[str, Path]) -> None:
    with pytest.raises(
        ValueError, match="population or direction basis differs from verified core"
    ):
        _run(inputs, "case_vs_control_subject")
    assert not inputs["output_dir"].exists()


def test_output_cannot_be_written_inside_the_verified_core(inputs: dict[str, Path]) -> None:
    changed = {**inputs, "output_dir": inputs["core_run"].parent / "assessment"}
    with pytest.raises(ValueError, match="overlaps verified core inputs or run bundle"):
        _run(changed)
    assert not changed["output_dir"].exists()


def test_noel_plan_binds_the_adapter_to_its_manifest(
    inputs: dict[str, Path], tmp_path: Path
) -> None:
    run = verify_core_run(inputs["core_run"])
    context = ArmContext(
        repository=REPO_ROOT,
        manifest=inputs["core_run"],
        output=tmp_path / "arm",
        inputs=tmp_path / "arm" / "inputs",
        products={
            artifact.product_id: artifact_path(run, artifact.product_id)
            for artifact in run.artifacts or []
        },
        code_commit="a" * 40,
        replicates=None,
        seed=None,
        figures=False,
        contrast_population=POPULATION,
    )
    command = next(
        command for command in build_commands(context) if command.id == "contrast_orientation"
    )
    assert command.argv[2] == "mtp.workflow.orientation_assessment"
    assert command.argv[command.argv.index("--core-run") + 1] == str(inputs["core_run"])
    assert inputs["core_run"] in command.inputs


def test_real_cli_accepts_explicit_selection_inputs(inputs: dict[str, Path]) -> None:
    argv = [
        part for key, path in inputs.items() for part in (f"--{key.replace('_', '-')}", str(path))
    ]
    assert main([*argv, "--population", POPULATION, "--code-commit", "a" * 40]) == 0
    assert (inputs["output_dir"] / "contrast_orientation_provenance.json").is_file()
