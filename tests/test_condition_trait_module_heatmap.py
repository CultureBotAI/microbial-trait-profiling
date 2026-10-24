import importlib.util
import sys
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "create_condition_trait_module_heatmap.py"
)
SPEC = importlib.util.spec_from_file_location("create_condition_trait_module_heatmap", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
heatmap = importlib.util.module_from_spec(SPEC)
with patch.object(sys, "path", [str(MODULE_PATH.parent), *sys.path]):
    SPEC.loader.exec_module(heatmap)


def test_load_aggregated_study_profiles_filters_to_human_experiments(tmp_path: Path) -> None:
    profile_path = tmp_path / "profiles.tsv"
    profile_path.write_text(
        "\n".join(
            [
                "study_id\texperiment_id\tinc_dec\tMETPO:1\tMETPO:2",
                "study1\tstudy1-1\tinc\t1\t0",
                "study1\tstudy1-1\tdec\t0\t2",
                "study2\tstudy2-1\tinc\t8\t8",
                "study3\tstudy3-1\tinc\t4\t0",
            ]
        ),
        encoding="utf-8",
    )

    study_matrix, profile_metadata = heatmap.load_aggregated_study_profiles(
        profile_path,
        human_experiment_ids={"study1-1", "study3-1"},
    )

    assert list(study_matrix.index) == ["study1", "study3"]
    assert study_matrix.loc["study1", "METPO:1"] == 1
    assert study_matrix.loc["study1", "METPO:2"] == 2
    assert profile_metadata.loc["study1", "profile_count"] == 2
    assert profile_metadata.loc["study3", "profile_count"] == 1

    # study_id is a column only, so the main() join and reset_index() work (#602).
    assert profile_metadata.index.name != "study_id"
    assert profile_metadata.reset_index()["study_id"].tolist() == ["study1", "study3"]
    titles = pd.DataFrame(
        {"title": ["first", "third"]}, index=pd.Index(["study1", "study3"], name="study_id")
    )
    assert profile_metadata.join(titles, on="study_id")["title"].tolist() == ["first", "third"]


PROFILE_HEADER = "study_id\texperiment_id\tinc_dec\tMETPO:1\tMETPO:2"
TWO_CONDITION_PROFILES = [
    "study1\tstudy1-1\tinc\t3\t0",
    "study1\tstudy1-1\tdec\t2\t1",
    "study1\tstudy1-2\tinc\t1\t2",
    "study1\tstudy1-3\tinc\t0\t5",
    "study2\tstudy2-1\tinc\t4\t4",
]


def two_condition_experiments() -> pd.DataFrame:
    """Study 1 has two Crohn's disease experiments and one ulcerative colitis experiment."""
    return pd.DataFrame(
        {
            "study_id": ["study1", "study1", "study1", "study2", "study3"],
            "experiment_id": ["study1-1", "study1-2", "study1-3", "study2-1", "study3-1"],
            "condition": [
                "Crohn's disease",
                "Crohn's disease",
                "Ulcerative colitis",
                "Ulcerative colitis",
                "Crohn's disease",
            ],
            "body_site": ["Feces", "Feces", "Colonic mucosa", "Colonic mucosa", "Feces"],
            "location": ["Sweden", "Sweden", "Sweden", "Norway", "Norway"],
            "host_species": ["Homo sapiens"] * 4 + ["Mus musculus"],
        }
    ).set_index(["study_id", "experiment_id"])


def assert_two_condition_modules(study_matrix: pd.DataFrame, study_metadata: pd.DataFrame) -> None:
    groups, group_metadata = heatmap.build_condition_modules(
        study_matrix, study_metadata, ["condition", "body_site"], 1, 1, 1, 1, False, 10
    )
    crohns = "Crohn's disease | Feces"
    colitis = "Ulcerative colitis | Colonic mucosa"
    counts = group_metadata[["study_count", "profile_count", "experiment_count"]]
    assert groups.loc[crohns].tolist() == [6, 3]
    assert groups.loc[colitis].tolist() == [4, 9]
    assert counts.loc[crohns].tolist() == [1, 3, 2]
    assert counts.loc[colitis].tolist() == [2, 2, 2]


def test_modules_take_each_experiments_own_condition_and_body_site(tmp_path: Path) -> None:
    """A two-condition study adds each experiment to its own module, not its majority's (#555)."""
    profile_path = tmp_path / "profiles.tsv"
    profile_path.write_text("\n".join([PROFILE_HEADER, *TWO_CONDITION_PROFILES]), encoding="utf-8")

    study_matrix, study_metadata = heatmap.load_aggregated_study_profiles(
        profile_path,
        experiment_metadata=two_condition_experiments(),
        group_by=["condition", "body_site"],
    )

    assert_two_condition_modules(study_matrix, study_metadata)


def test_modules_stay_aligned_when_the_human_host_filter_drops_rows(tmp_path: Path) -> None:
    """Dropped non-human rows must not shift labels onto other profiles' counts (#646)."""
    profile_path = tmp_path / "profiles.tsv"
    mouse_rows = ["study3\tstudy3-1\tinc\t50\t0", "study3\tstudy3-1\tdec\t0\t50"]
    rows = [*TWO_CONDITION_PROFILES[:2], *mouse_rows, *TWO_CONDITION_PROFILES[2:]]
    profile_path.write_text("\n".join([PROFILE_HEADER, *rows]), encoding="utf-8")

    study_matrix, study_metadata = heatmap.load_aggregated_study_profiles(
        profile_path,
        human_experiment_ids={"study1-1", "study1-2", "study1-3", "study2-1"},
        experiment_metadata=two_condition_experiments(),
        group_by=["condition", "body_site"],
    )

    assert "study3" not in set(study_metadata["study_id"])
    assert_two_condition_modules(study_matrix, study_metadata)


def test_experiment_metadata_without_group_by_raises(tmp_path: Path) -> None:
    """Without grouping fields every experiment would get its study's majority label (#646)."""
    profile_path = tmp_path / "profiles.tsv"
    profile_path.write_text("\n".join([PROFILE_HEADER, *TWO_CONDITION_PROFILES]), encoding="utf-8")

    no_group_by: list[list[str] | None] = [None, []]
    for group_by in no_group_by:
        with pytest.raises(ValueError, match="needs group_by"):
            heatmap.load_aggregated_study_profiles(
                profile_path,
                experiment_metadata=two_condition_experiments(),
                group_by=group_by,
            )


def test_main_builds_modules_from_each_experiments_own_labels(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """main() passes the grouping fields, so modules follow experiment labels (#646)."""
    heatmap.plt.switch_backend("Agg")
    profiles = tmp_path / "profiles.tsv"
    profiles.write_text(
        "\n".join(
            [
                "study_id\texperiment_id\tinc_dec\tMETPO:1\tMETPO:2\tMETPO:3",
                "bsdb:1\tbsdb:1-1\tinc\t30\t1\t5",
                "bsdb:1\tbsdb:1-1\tdec\t20\t2\t5",
                "bsdb:1\tbsdb:1-2\tinc\t25\t3\t5",
                "bsdb:1\tbsdb:1-3\tinc\t1\t30\t5",
                "bsdb:2\tbsdb:2-1\tinc\t2\t40\t5",
                "bsdb:3\tbsdb:3-1\tinc\t9\t9\t30",
            ]
        ),
        encoding="utf-8",
    )
    full_dump = tmp_path / "full_dump.csv"
    full_dump.write_text(
        "\n".join(
            [
                "# BugSigDB 2026-04-24_00:41_UTC",
                "BSDB ID,Condition,Body site,Location of subjects,Host species,Sequencing type",
                "bsdb:1/1/1,Crohn's disease,Feces,Sweden,Homo sapiens,16S",
                "bsdb:1/1/2,Crohn's disease,Feces,Sweden,Homo sapiens,16S",
                "bsdb:1/2/1,Crohn's disease,Feces,Sweden,Homo sapiens,16S",
                "bsdb:1/3/1,Ulcerative colitis,Colonic mucosa,Sweden,Homo sapiens,16S",
                "bsdb:2/1/1,Ulcerative colitis,Colonic mucosa,Norway,Homo sapiens,WMS",
                "bsdb:3/1/1,Obesity,Feces,Norway,Homo sapiens,16S",
            ]
        ),
        encoding="utf-8",
    )
    studies = tmp_path / "studies.tsv"
    studies.write_text(
        "study_id\ttitle\tpublication_year\ttaxon_count\tsignature_count\n"
        "bsdb:1\tIBD cohort\t2020\t10\t4\n"
        "bsdb:2\tColitis cohort\t2021\t8\t1\n"
        "bsdb:3\tObesity cohort\t2022\t6\t1\n",
        encoding="utf-8",
    )
    columns = tmp_path / "columns.tsv"
    columns.write_text(
        "column_id\tobject_name\tcolumn_label\tobject_category\n"
        "METPO:1\ttrait one\ttrait one\tbiolink:PhenotypicQuality\n"
        "METPO:2\ttrait two\ttrait two\tbiolink:PhenotypicQuality\n"
        "METPO:3\ttrait three\ttrait three\tbiolink:PhenotypicQuality\n",
        encoding="utf-8",
    )
    group_metadata = tmp_path / "group_metadata.tsv"
    arguments = {
        "--study-trait-profiles-input": profiles,
        "--trait-columns-input": columns,
        "--studies-input": studies,
        "--full-dump-input": full_dump,
        "--output": tmp_path / "heatmap.png",
        "--matrix-output": tmp_path / "matrix.tsv",
        "--group-metadata-output": group_metadata,
        "--trait-metadata-output": tmp_path / "traits.tsv",
        "--guild-definition-table-output": tmp_path / "guild.png",
        "--guild-definition-table-svg-output": tmp_path / "guild.svg",
        "--min-studies": 1,
        "--min-profiles": 1,
        "--min-traits": 1,
        "--min-total-signal": 1,
        "--min-group-support": 1,
        "--max-group-support-fraction": 1.0,
    }
    argv = [str(MODULE_PATH)]
    for flag, value in arguments.items():
        argv.extend([flag, str(value)])
    monkeypatch.setattr(sys, "argv", argv)
    # Match the import path of the standalone CLI; main owns execution-only imports.
    monkeypatch.syspath_prepend(str(MODULE_PATH.parent))

    heatmap.main()

    modules = pd.read_csv(group_metadata, sep="\t", index_col="module_id")
    assert modules.loc["Crohn's disease | Feces", "profile_count"] == 3
    assert modules.loc["Crohn's disease | Feces", "study_count"] == 1
    assert modules.loc["Ulcerative colitis | Colonic mucosa", "profile_count"] == 2
    assert modules.loc["Ulcerative colitis | Colonic mucosa", "study_count"] == 2
