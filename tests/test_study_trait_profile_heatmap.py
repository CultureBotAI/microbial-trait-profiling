import importlib.util
import runpy
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import QuadMesh

from mtp.paths import file_sha256

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "create_study_trait_profile_heatmap.py"
)
SPEC = importlib.util.spec_from_file_location("create_study_trait_profile_heatmap", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
heatmap = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(heatmap)


def test_main_preserves_direction_basis_in_figures_and_drawn_values(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exercise both real main() plotting calls and compare their drawn cells (#955)."""
    fixture = runpy.run_path(
        str(MODULE_PATH.parent / "adhoc/issue842_render_direction_fixtures.py")
    )
    inputs = fixture["write_inputs"](tmp_path / "inputs")
    products = {
        "--output": "profiles.png",
        "--matrix-output": "matrix.tsv",
        "--study-metadata-output": "profiles.tsv",
        "--trait-metadata-output": "traits.tsv",
        "--paired-output": "paired.png",
        "--paired-matrix-output": "paired_matrix.tsv",
        "--paired-metadata-output": "paired_metadata.tsv",
        "--paired-required-trait-check-output": "required_traits.tsv",
    }
    drawn = {}
    for basis, label in (
        ("group_1_vs_group_0", "Group 1 vs Group 0"),
        ("disease_vs_comparison", "disease vs comparison"),
    ):
        target = tmp_path / basis
        target.mkdir()
        options = {**inputs, **{flag: target / name for flag, name in products.items()}}
        argv = [str(MODULE_PATH), *(str(part) for item in options.items() for part in item)]
        monkeypatch.setattr(
            sys, "argv", [*argv, "--direction-basis", basis, "--min-profile-support", "1"]
        )
        with (
            patch.object(heatmap.plt.Figure, "savefig", autospec=True) as savefig,
            patch.object(heatmap.plt, "close"),
        ):
            heatmap.main()
        figures = {Path(call.args[1]).name: call.args[0] for call in savefig.call_args_list}
        try:
            metadata = pd.read_csv(target / "paired_metadata.tsv", sep="\t")
            assert len(metadata) == 3
            assert "heatmap_parameter_direction_basis" in metadata, (
                "CLI basis is missing from paired metadata"
            )
            assert metadata["heatmap_parameter_direction_basis"].eq(basis).all(), (
                "CLI basis did not reach paired metadata"
            )
            assert set(figures) == {"profiles.png", "paired.png"}
            titles = {
                axis.get_legend().get_title().get_text()
                for axis in figures["profiles.png"].axes
                if axis.get_legend()
            }
            assert f"BugSigDB direction\n({label})" in titles, (
                "CLI basis did not reach profile legend"
            )
            labels = {axis.get_ylabel() for axis in figures["paired.png"].axes}
            assert any(f"BugSigDB direction: {label}" in text for text in labels), (
                "CLI basis did not reach paired colorbar"
            )
            for name, figure in figures.items():
                arrays = [
                    np.asarray(mesh.get_array()).copy()
                    for axis in figure.axes
                    if axis.get_xlabel() == "Trait"
                    for mesh in axis.collections
                    if isinstance(mesh, QuadMesh)
                ]
                assert len(arrays) == 1
                assert arrays[0].size == (18 if name == "profiles.png" else 9)
                assert np.ptp(arrays[0]) > 0
                drawn[basis, name] = arrays[0]
        finally:
            for figure in figures.values():
                heatmap.plt.close(figure)
    for name in ("profiles.png", "paired.png"):
        np.testing.assert_array_equal(
            drawn["group_1_vs_group_0", name], drawn["disease_vs_comparison", name]
        )


@pytest.mark.parametrize(
    ("direction_basis", "label"),
    [
        ("group_1_vs_group_0", "Group 1 vs Group 0"),
        ("disease_vs_comparison", "disease vs comparison"),
    ],
)
@pytest.mark.parametrize("paired", [False, True])
def test_rendered_legends_name_the_profile_direction_basis(
    tmp_path: Path, direction_basis: str, label: str, paired: bool
) -> None:
    """Both real plotting paths name the supplied basis without changing values (#842)."""
    matrix = pd.DataFrame(
        [[0.7, -0.3, 0.1], [-0.2, 0.8, -0.4], [0.3, -0.6, 0.9]],
        index=["a", "b", "c"],
        columns=["Trait A", "Trait B", "Trait C"],
    )
    metadata = pd.DataFrame(
        {
            "inc_dec": ["inc", "dec", "dec"],
            "body_site_group": ["Feces"] * 3,
            "condition_group": ["Condition A", "Condition B", "Condition C"],
        },
        index=matrix.index,
    )
    traits = pd.DataFrame(
        {"display_category": ["Shape", "Shape", "Metabolism"]},
        index=["METPO:1", "METPO:2", "METPO:3"],
    )
    before = matrix.copy(deep=True)
    with (
        patch.object(heatmap.plt.Figure, "savefig", autospec=True) as savefig,
        patch.object(heatmap.plt, "close"),
    ):
        draw = heatmap.save_paired_enrichment_heatmap if paired else heatmap.save_heatmap
        draw(matrix, metadata, traits, tmp_path / "heatmap.png", 3, 3, direction_basis)
    figure = savefig.call_args.args[0]
    try:
        figure.canvas.draw()
        if paired:
            labels = [axis.get_ylabel() for axis in figure.axes]
            assert any(f"BugSigDB direction: {label}" in text for text in labels)
            assert any("increased minus decreased" in text for text in labels)
        else:
            legends = [axis.get_legend() for axis in figure.axes if axis.get_legend()]
            assert f"BugSigDB direction\n({label})" in {
                legend.get_title().get_text() for legend in legends
            }
        pd.testing.assert_frame_equal(matrix, before)
    finally:
        heatmap.plt.close(figure)


@pytest.mark.parametrize("arguments", [[], ["--direction-basis", "unknown"]])
def test_cli_requires_a_known_direction_basis_before_reading_inputs(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], arguments: list[str]
) -> None:
    def reject_input_read(*_args: object) -> None:
        raise AssertionError("CLI must reject an absent/invalid basis before reading inputs")

    monkeypatch.setattr(heatmap, "load_profile_matrix", reject_input_read)
    monkeypatch.setattr(sys, "argv", ["heatmap", *arguments])
    with pytest.raises(SystemExit) as error:
        heatmap.main()
    assert error.value.code == 2
    assert "--direction-basis" in capsys.readouterr().err


def test_paired_export_preserves_the_count_input_identity(tmp_path: Path) -> None:
    """The actual pair export carries identities checked by site-stat regeneration (#824)."""
    profiles = tmp_path / "profiles.tsv"
    profiles.write_text("synthetic counts\n")
    source_hashes = {"profile": file_sha256(profiles)}
    matrix = pd.DataFrame({"METPO:1": [1.0]}, index=["pair1"])
    metadata = pd.DataFrame({"study_id": ["study1"]}, index=matrix.index)
    output = tmp_path / "pairs.tsv"
    heatmap.write_paired_outputs(
        matrix,
        matrix,
        metadata,
        pd.DataFrame({"present": [True]}),
        1,
        tmp_path / "figure.png",
        tmp_path / "matrix.tsv",
        output,
        tmp_path / "required.tsv",
        source_hashes=source_hashes,
    )
    observed = pd.read_csv(output, sep="\t", index_col=0)
    assert observed["profile_sha256"].eq(source_hashes["profile"]).all()
    assert "profile_sha256" not in metadata


def test_filter_human_host_profiles_keeps_only_homo_sapiens() -> None:
    profile_matrix = pd.DataFrame(
        {"METPO:1": [1, 2, 3, 4]},
        index=[
            "study1|exp1|inc",
            "study1|exp1|dec",
            "study2|exp1|inc",
            "study3|exp1|inc",
        ],
    )
    profile_metadata = pd.DataFrame(
        {
            "study_id": ["study1", "study1", "study2", "study3"],
            "experiment_id": ["exp1", "exp1", "exp1", "exp1"],
            "inc_dec": ["inc", "dec", "inc", "inc"],
            "host_species": [
                "Homo sapiens",
                " Homo sapiens ",
                "Mus musculus",
                "Canis lupus familiaris",
            ],
            "sequencing_type": ["16S", "16S", "16S", "WMS"],
        },
        index=profile_matrix.index,
    )

    filtered_matrix, filtered_metadata = heatmap.filter_human_host_profiles(
        profile_matrix,
        profile_metadata,
    )

    assert list(filtered_matrix.index) == ["study1|exp1|inc", "study1|exp1|dec"]
    assert list(filtered_metadata.index) == ["study1|exp1|inc", "study1|exp1|dec"]
    assert set(filtered_metadata["host_species"]) == {"Homo sapiens"}


def test_filter_human_host_profiles_drops_other_sequencing_types() -> None:
    """A human PCR or untyped profile in a non-default input is left out (#754)."""
    profile_matrix = pd.DataFrame(
        {"METPO:1": [1, 2, 3]}, index=["s1|e1|inc", "s2|e1|inc", "s3|e1|inc"]
    )
    profile_metadata = pd.DataFrame(
        {"host_species": "Homo sapiens", "sequencing_type": ["WMS", "PCR", None]},
        index=profile_matrix.index,
    )
    filtered_matrix, _ = heatmap.filter_human_host_profiles(profile_matrix, profile_metadata)
    assert list(filtered_matrix.index) == ["s1|e1|inc"]


def test_filter_human_host_profiles_fails_without_host_species() -> None:
    profile_matrix = pd.DataFrame({"METPO:1": [1]}, index=["study1|exp1|inc"])
    profile_metadata = pd.DataFrame({"study_id": ["study1"]}, index=profile_matrix.index)

    with pytest.raises(ValueError, match="missing host_species"):
        heatmap.filter_human_host_profiles(profile_matrix, profile_metadata)


def test_load_raw_experiment_metadata_maps_bsdb_id_to_experiment_id(tmp_path: Path) -> None:
    full_dump_path = tmp_path / "full_dump.csv"
    full_dump_path.write_text(
        "\n".join(
            [
                "# comment",
                "BSDB ID,Host species",
                "bsdb:study1/1/1,Homo sapiens",
                "bsdb:study1/2/1,Mus musculus",
                "bsdb:study2/3/1, Homo sapiens ",
            ]
        ),
        encoding="utf-8",
    )

    hosts = heatmap.load_raw_experiment_metadata(full_dump_path)

    assert hosts.loc[("bsdb:study1", "bsdb:study1-1"), "host_species"] == "Homo sapiens"
    assert hosts.loc[("bsdb:study1", "bsdb:study1-2"), "host_species"] == "Mus musculus"
    assert hosts.loc[("bsdb:study2", "bsdb:study2-3"), "host_species"] == "Homo sapiens"
    assert hosts.loc[("bsdb:study2", "bsdb:study2-3"), "condition"] == "unknown"


def write_two_experiment_study(tmp_path: Path) -> tuple[Path, Path]:
    """One study whose two experiments differ in condition and body site (#555)."""
    full_dump_path = tmp_path / "full_dump.csv"
    full_dump_path.write_text(
        "\n".join(
            [
                "# BugSigDB 2026-04-24_00:41_UTC",
                "BSDB ID,Condition,Body site,Location of subjects,Host species",
                "bsdb:study1/1/1,Crohn's disease,Feces,Sweden,Homo sapiens",
                "bsdb:study1/1/2,Crohn's disease,Feces,Sweden,Homo sapiens",
                "bsdb:study1/2/1,Ulcerative colitis,Colonic mucosa,Norway,Homo sapiens",
            ]
        ),
        encoding="utf-8",
    )
    studies_path = tmp_path / "studies.tsv"
    studies_path.write_text(
        "study_id\ttitle\tpublication_year\tsignature_count\ttaxon_count\n"
        "bsdb:study1\tInflammatory bowel disease cohort\t2020\t3\t12\n",
        encoding="utf-8",
    )
    return full_dump_path, studies_path


def test_join_profile_labels_takes_each_profiles_own_experiment_labels(tmp_path: Path) -> None:
    full_dump_path, studies_path = write_two_experiment_study(tmp_path)
    profile_ids = [
        "bsdb:study1|bsdb:study1-1|inc",
        "bsdb:study1|bsdb:study1-1|dec",
        "bsdb:study1|bsdb:study1-2|inc",
        "bsdb:study1|bsdb:study1-2|dec",
    ]
    profile_metadata = pd.DataFrame(
        {
            "study_id": ["bsdb:study1"] * 4,
            "experiment_id": ["bsdb:study1-1", "bsdb:study1-1", "bsdb:study1-2", "bsdb:study1-2"],
            "inc_dec": ["inc", "dec", "inc", "dec"],
        },
        index=profile_ids,
    )

    labelled = heatmap.join_profile_labels(profile_metadata, studies_path, full_dump_path)

    assert labelled["condition"].tolist() == [
        "Crohn's disease",
        "Crohn's disease",
        "Ulcerative colitis",
        "Ulcerative colitis",
    ]
    assert labelled["body_site"].tolist() == ["Feces", "Feces", "Colonic mucosa", "Colonic mucosa"]
    assert labelled["location"].tolist() == ["Sweden", "Sweden", "Norway", "Norway"]
    study_columns = [
        "study_majority_condition",
        "study_majority_body_site",
        "study_majority_location",
        "study_condition_count",
        "study_body_site_count",
    ]
    assert labelled[study_columns].drop_duplicates().to_dict("records") == [
        {
            "study_majority_condition": "Crohn's disease",
            "study_majority_body_site": "Feces",
            "study_majority_location": "Sweden",
            "study_condition_count": 2,
            "study_body_site_count": 2,
        }
    ]
    assert set(labelled["title"]) == {"Inflammatory bowel disease cohort"}


LABEL_VALUES = {
    "condition": ("Crohn's disease", "Ulcerative colitis"),
    "body_site": ("Feces", "Colonic mucosa"),
    "location": ("Sweden", "Norway"),
    "host_species": ("Homo sapiens", "Mus musculus"),
}


def write_one_experiment(tmp_path: Path, column: str, second_value: str) -> Path:
    """Two signatures of one experiment that differ only in `column`."""
    first = {field: values[0] for field, values in LABEL_VALUES.items()}
    second = {**first, column: second_value}
    header = ",".join(["BSDB ID", *heatmap.EXPERIMENT_LABEL_COLUMNS.values()])
    rows = [
        ",".join([f"bsdb:study1/1/{number}", *(labels[field] for field in LABEL_VALUES)])
        for number, labels in [(1, first), (2, second)]
    ]
    full_dump_path = tmp_path / "full_dump.csv"
    full_dump_path.write_text("\n".join([header, *rows]), encoding="utf-8")
    return full_dump_path


@pytest.mark.parametrize("column", list(LABEL_VALUES))
def test_load_raw_experiment_metadata_rejects_conflicting_experiment_labels(
    tmp_path: Path, column: str
) -> None:
    full_dump_path = write_one_experiment(tmp_path, column, LABEL_VALUES[column][1])

    with pytest.raises(ValueError, match=f"more than one {column} label"):
        heatmap.load_raw_experiment_metadata(full_dump_path)


@pytest.mark.parametrize("blank", ["", "NA"])
@pytest.mark.parametrize("column", list(LABEL_VALUES))
def test_load_raw_experiment_metadata_takes_the_value_beside_a_blank(
    tmp_path: Path, column: str, blank: str
) -> None:
    full_dump_path = write_one_experiment(tmp_path, column, blank)

    labels = heatmap.load_raw_experiment_metadata(full_dump_path)

    assert labels.loc[("bsdb:study1", "bsdb:study1-1"), column] == LABEL_VALUES[column][0]


def test_select_traits_retains_required_paired_traits_above_max_support_cap() -> None:
    profiles = [f"study{idx}|exp1|inc" for idx in range(20)]
    matrix = pd.DataFrame(
        {
            "METPO:1000615": [10] * 19 + [1],
            "METPO:1000698": [8, 2] * 9 + [0, 0],
            "METPO:1000703": [5, 5] * 9 + [0, 0],
            "GO:0030245": [4, 1] * 8 + [0] * 4,
        },
        index=profiles,
    )
    labels = pd.DataFrame(
        {
            "display_label": [
                "mesophilic (METPO:1000615)",
                "gram positive (METPO:1000698)",
                "non motile (METPO:1000703)",
                "cellulose catabolic process (GO:0030245)",
            ],
            "column_label": [
                "METPO:1000615",
                "METPO:1000698",
                "METPO:1000703",
                "GO:0030245",
            ],
            "object_name": [
                "mesophilic",
                "gram positive",
                "non motile",
                "cellulose catabolic process",
            ],
            "object_category": ["METPO", "METPO", "METPO", "GO"],
            "display_category": ["Temperature", "Gram stain", "Motility", "Metabolism"],
        },
        index=["METPO:1000615", "METPO:1000698", "METPO:1000703", "GO:0030245"],
    )

    selected_ids, selected_meta = heatmap.select_traits(
        selected_profile_matrix=matrix,
        trait_labels=labels,
        max_traits=4,
        min_profile_support=5,
        max_profile_support_fraction=0.90,
        include_media_traits=False,
    )

    assert set(heatmap.REQUIRED_PAIRED_TRAITS).issubset(selected_ids)
    assert "METPO:1000615" in selected_meta.index
