import importlib.util
import json
import math
import sys
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.patches import Circle

from mtp.integrated_evidence import native_composition_records
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_guilds import TABLE1_TRAITS
from mtp.workflow.consumers import CoreBundle, open_core_bundle, verify_consumer_record
from mtp.workflow.operations import OPERATIONS, validate_core_product
from mtp.workflow.runtime import execute_workflow

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "create_trait_profile_infographic.py"
)
SPEC = importlib.util.spec_from_file_location("create_trait_profile_infographic", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
infographic = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(infographic)


@pytest.fixture
def stratified_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Synthetic rebuilt counts and heatmap-selected pairs for #824."""
    profiles = tmp_path / "profiles.tsv"
    labels = tmp_path / "labels.tsv"
    pairs = tmp_path / "pairs.tsv"
    rows = [
        {
            "study_id": f"s{i}",
            "experiment_id": "e1",
            "inc_dec": direction,
            "METPO:1": i + 1 if direction == "inc" else 0,
            "METPO:2": i + 1 if direction == "dec" else 0,
        }
        for i in range(16)
        for direction in ("inc", "dec")
    ]
    pd.DataFrame(rows).to_csv(profiles, sep="\t", index=False)
    pd.DataFrame(
        {
            "column_id": ["METPO:1", "METPO:2"],
            "column_label": ["oxygen", "motility"],
            "object_name": ["oxygen", "motility"],
            "object_category": ["biolink:OntologyClass"] * 2,
        }
    ).to_csv(labels, sep="\t", index=False)
    pd.DataFrame(
        {
            "study_id": [f"s{i}" for i in range(16)],
            "experiment_id": ["e1"] * 16,
            "body_site_group": ["Feces"] * 8 + ["Mouth"] * 8,
            "condition": ["Disease"] * 16,
            "increased_profile_ids": [f"s{i}|e1|inc" for i in range(16)],
            "decreased_profile_ids": [f"s{i}|e1|dec" for i in range(16)],
            "profile_sha256": file_sha256(profiles),
            "trait_labels_sha256": file_sha256(labels),
        }
    ).to_csv(pairs, sep="\t")
    return profiles, pairs, labels


def test_stratified_stats_only_records_outputs_and_stale_files(
    stratified_inputs: tuple[Path, Path, Path], tmp_path: Path
) -> None:
    profiles, pairs, labels = stratified_inputs
    output = tmp_path / "analysis"
    output.mkdir()
    stale = output / "trait_profile_paired_contrast_trait_stats_saliva.tsv"
    stale_contents = {
        stale: b"historical table\n",
        output / infographic.CONDITION_SUMMARY_TSV.name: b"historical condition summary\n",
        output / infographic.BODY_SITE_SUMMARY_TSV.name: b"historical body-site summary\n",
        output / infographic.OUTPUT_PNG.name: b"historical infographic bytes\n",
    }
    for path, content in stale_contents.items():
        path.write_bytes(content)
    unrelated = output / "unrelated.png"
    unrelated.write_bytes(b"unrelated image\n")
    pooled = infographic.trait_stats_path(output, "pooled")
    pooled_alias = output / "trait_profile_paired_contrast_trait_stats.tsv"
    for path in (pooled, pooled_alias):
        path.write_text("historical pooled table\n")
    infographic.main(
        [
            "--stratified-stats-only",
            "--profiles",
            str(profiles),
            "--pair-metadata",
            str(pairs),
            "--trait-labels",
            str(labels),
            "--output-dir",
            str(output),
            "--min-stratified-pairs",
            "2",
        ]
    )
    manifest = pd.read_csv(infographic.manifest_path(output), sep="\t")
    assert set(manifest["filename_token"]) == {"feces", "feces_excluded", "mouth"}
    assert manifest["pair_count"].eq(8).all()
    record = json.loads((output / "stratified_stats_provenance.json").read_text())
    assert record["verification_status"] == "legacy_inputs_hashed_not_core_verified"
    assert record["inputs"]["profiles"]["sha256"] == file_sha256(profiles)
    assert record["inputs"]["pair_metadata"]["sha256"] == file_sha256(pairs)
    assert record["inputs"]["trait_labels"]["sha256"] == file_sha256(labels)
    assert record["parameters"]["heatmap_selection_status"] == "unavailable"
    assert len(record["outputs"]) == 6
    for item in record["outputs"]:
        assert item["sha256"] == file_sha256(output / Path(item["path"]).name)
        assert not Path(item["path"]).is_absolute()
    assert not Path(record["producer"]["path"]).is_absolute()
    assert all(not Path(item["path"]).is_absolute() for item in record["inputs"].values())
    assert {Path(item["path"]).name for item in record["not_regenerated"]} == {
        path.name for path in stale_contents
    }
    for item in record["not_regenerated"]:
        path = output / Path(item["path"]).name
        assert path.read_bytes() == stale_contents[path]
        assert item["sha256"] == file_sha256(path)
        assert item["size_bytes"] == len(stale_contents[path])
    assert set(output.glob("*.png")) == {output / infographic.OUTPUT_PNG.name, unrelated}
    assert pooled.read_bytes() == pooled_alias.read_bytes()
    native_output = tmp_path / "native-pooled"
    infographic.main(
        [
            "--stats-only",
            "--profiles",
            str(profiles),
            "--pair-metadata",
            str(pairs),
            "--trait-labels",
            str(labels),
            "--output-dir",
            str(native_output),
        ]
    )
    assert pooled.read_bytes() == infographic.trait_stats_path(native_output, "pooled").read_bytes()


@pytest.mark.parametrize(
    "problem",
    [
        "changed_counts",
        "changed_labels",
        "unbound_pairs",
        "wrong_pair",
        "duplicate_pair",
        "mixed_parameters",
    ],
)
def test_stratified_stats_rejects_stale_or_mismatched_pair_metadata(
    stratified_inputs: tuple[Path, Path, Path], tmp_path: Path, problem: str
) -> None:
    profiles, pairs, labels = stratified_inputs
    metadata = pd.read_csv(pairs, sep="\t", index_col=0)
    if problem == "changed_counts":
        profiles.write_text(profiles.read_text().replace("\t1\t0", "\t2\t0"))
    elif problem == "changed_labels":
        labels.write_text(labels.read_text().replace("oxygen", "other label"))
    elif problem == "unbound_pairs":
        metadata = metadata.drop(columns="profile_sha256")
    elif problem == "duplicate_pair":
        metadata.loc[1] = metadata.loc[0]
    elif problem == "mixed_parameters":
        metadata["heatmap_parameter_max_profiles"] = 0
        metadata.loc[0, "heatmap_parameter_max_profiles"] = 24
    else:
        metadata.loc[0, "increased_profile_ids"] = "s1|e1|inc"
    metadata.to_csv(pairs, sep="\t")
    output = tmp_path / "analysis"
    with pytest.raises(ValueError, match=r"Rebuild|mismatched"):
        infographic.main(
            [
                "--stratified-stats-only",
                "--profiles",
                str(profiles),
                "--pair-metadata",
                str(pairs),
                "--trait-labels",
                str(labels),
                "--output-dir",
                str(output),
                "--min-stratified-pairs",
                "2",
            ]
        )
    assert not output.exists()


def test_stratified_stats_no_eligible_sites_has_readable_empty_manifest(
    stratified_inputs: tuple[Path, Path, Path], tmp_path: Path
) -> None:
    profiles, pairs, labels = stratified_inputs
    output = tmp_path / "analysis"
    infographic.write_legacy_stratified_stats(profiles, pairs, labels, output)
    assert pd.read_csv(infographic.manifest_path(output), sep="\t").empty
    record = json.loads((output / "stratified_stats_provenance.json").read_text())
    assert len(record["outputs"]) == 3


@pytest.mark.parametrize("mode", [[], ["--stats-only"]])
def test_min_stratified_pairs_is_rejected_outside_site_stats(
    stratified_inputs: tuple[Path, Path, Path],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    mode: list[str],
) -> None:
    profiles, pairs, labels = stratified_inputs
    output = tmp_path / "isolated-legacy-output"
    monkeypatch.setattr(infographic, "INVESTIGATION_DIR", output)
    for name in (
        "TRAIT_STATS_TSV",
        "POOLED_TRAIT_STATS_TSV",
        "CONDITION_SUMMARY_TSV",
        "BODY_SITE_SUMMARY_TSV",
        "STRATIFIED_TRAIT_STATS_MANIFEST_TSV",
        "OUTPUT_PNG",
    ):
        monkeypatch.setattr(infographic, name, output / getattr(infographic, name).name)
    monkeypatch.setattr(infographic, "create_infographic", lambda *_args, **_kwargs: None)
    with pytest.raises(SystemExit, match="2"):
        infographic.main(
            [
                *mode,
                "--profiles",
                str(profiles),
                "--pair-metadata",
                str(pairs),
                "--trait-labels",
                str(labels),
                "--output-dir",
                str(output),
                "--min-stratified-pairs",
                "2",
            ]
        )
    assert "--min-stratified-pairs requires --stratified-stats-only" in capsys.readouterr().err
    assert not output.exists()


@pytest.mark.parametrize("mode", [[], ["--stats-only"]])
def test_other_modes_invalidate_old_site_provenance(
    stratified_inputs: tuple[Path, Path, Path],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mode: list[str],
) -> None:
    profiles, pairs, labels = stratified_inputs
    output = tmp_path / "full"
    output.mkdir()
    provenance = output / "stratified_stats_provenance.json"
    provenance.write_text('{"obsolete": true}\n')
    monkeypatch.setattr(infographic, "INVESTIGATION_DIR", output)
    for name in (
        "TRAIT_STATS_TSV",
        "POOLED_TRAIT_STATS_TSV",
        "CONDITION_SUMMARY_TSV",
        "BODY_SITE_SUMMARY_TSV",
        "STRATIFIED_TRAIT_STATS_MANIFEST_TSV",
        "OUTPUT_PNG",
    ):
        monkeypatch.setattr(infographic, name, output / getattr(infographic, name).name)
    monkeypatch.setattr(infographic, "create_infographic", lambda *_args, **_kwargs: None)
    infographic.main(
        [
            *mode,
            "--profiles",
            str(profiles),
            "--pair-metadata",
            str(pairs),
            "--trait-labels",
            str(labels),
        ]
    )
    assert not provenance.exists()
    assert infographic.trait_stats_path(output, "pooled").exists()


@pytest.mark.parametrize("direction_basis", ["group_1_vs_group_0", "disease_vs_comparison"])
def test_heatmap_main_to_site_stats_binds_hashes_and_selection_parameters(
    stratified_inputs: tuple[Path, Path, Path],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    direction_basis: str,
) -> None:
    """Exercise the producer mapping instead of hand-writing hash columns (#875)."""
    heatmap_path = MODULE_PATH.with_name("create_study_trait_profile_heatmap.py")
    heatmap_spec = importlib.util.spec_from_file_location("site_stats_heatmap", heatmap_path)
    assert heatmap_spec is not None
    assert heatmap_spec.loader is not None
    heatmap = importlib.util.module_from_spec(heatmap_spec)
    heatmap_spec.loader.exec_module(heatmap)
    profiles, pairs, labels = stratified_inputs
    frame = pd.read_csv(profiles, sep="\t")
    frame["study_id"] = "bsdb:" + frame["study_id"]
    frame["experiment_id"] = frame["study_id"] + "-1"
    frame["METPO:1000615"] = frame["METPO:1"]
    frame["METPO:1000698"] = frame["METPO:2"]
    frame["METPO:1000703"] = frame["METPO:1"] * 2 + frame["METPO:2"]
    frame.to_csv(profiles, sep="\t", index=False)
    label_frame = pd.read_csv(labels, sep="\t")
    pd.concat(
        [
            label_frame,
            pd.DataFrame(
                {
                    "column_id": list(heatmap.REQUIRED_PAIRED_TRAITS),
                    "column_label": list(heatmap.REQUIRED_PAIRED_TRAITS.values()),
                    "object_name": list(heatmap.REQUIRED_PAIRED_TRAITS.values()),
                    "object_category": "biolink:OntologyClass",
                }
            ),
        ],
        ignore_index=True,
    ).to_csv(labels, sep="\t", index=False)
    studies = tmp_path / "studies.tsv"
    pd.DataFrame(
        {
            "study_id": [f"bsdb:s{i}" for i in range(16)],
            "title": "Synthetic study",
            "publication_year": 2026,
            "taxon_count": 2,
            "signature_count": 2,
        }
    ).to_csv(studies, sep="\t", index=False)
    dump = tmp_path / "dump.csv"
    pd.DataFrame(
        {
            "BSDB ID": [f"bsdb:s{i}/1/1" for i in range(16)],
            "Condition": "Disease",
            "Body site": ["Feces"] * 8 + ["Mouth"] * 8,
            "Host species": "Homo sapiens",
            "Sequencing type": "16S",
        }
    ).to_csv(dump, index=False)
    argv = [
        "heatmap",
        "--direction-basis",
        direction_basis,
        "--study-trait-profiles-input",
        str(profiles),
        "--trait-columns-input",
        str(labels),
        "--studies-input",
        str(studies),
        "--full-dump-input",
        str(dump),
        "--paired-metadata-output",
        str(pairs),
        "--max-profiles",
        "24",
        "--max-body-site-categories",
        "1",
    ]
    for flag in (
        "output",
        "matrix-output",
        "study-metadata-output",
        "trait-metadata-output",
        "paired-output",
        "paired-matrix-output",
        "paired-required-trait-check-output",
    ):
        argv.extend([f"--{flag}", str(tmp_path / f"{flag}.tsv")])
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(heatmap, "save_heatmap", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(heatmap, "save_paired_enrichment_heatmap", lambda *_args, **_kwargs: None)
    heatmap.main()
    output = tmp_path / "site-stats"
    infographic.write_legacy_stratified_stats(profiles, pairs, labels, output, min_pairs=6)
    record = json.loads((output / "stratified_stats_provenance.json").read_text())
    assert record["inputs"]["profiles"]["sha256"] == file_sha256(profiles)
    assert record["parameters"]["heatmap_selection"]["max_profiles"] == 24
    assert record["parameters"]["heatmap_selection"]["max_body_site_categories"] == 1
    assert "direction_basis" in record["parameters"]["heatmap_selection"], (
        "Paired metadata did not preserve the basis in site-stat provenance"
    )
    assert record["parameters"]["heatmap_selection"]["direction_basis"] == direction_basis
    assert 0 < record["pair_denominator_audit"]["n_experiments"] < 16


@pytest.mark.parametrize(
    ("values", "studies", "expected_state"),
    [
        ([1.0, -1.0, 1.0, -1.0], ["s1", "s1", "s2", "s2"], "all_tied"),
        ([2.0, 0.0, 2.0, 0.0], ["s1", "s1", "s2", "s2"], "insufficient_variance"),
        ([1.0, -1.0], ["s1", "s1"], "insufficient_support"),
        ([0.0, 0.0], ["s1", "s2"], "all_tied"),
    ],
)
@pytest.mark.parametrize("legacy", [False, True])
def test_native_producer_untestable_study_results_integrate_without_invented_tests(
    tmp_path: Path, values: list[float], studies: list[str], expected_state: str, legacy: bool
) -> None:
    """Actual producer output preserves study versus pair test eligibility (#706)."""
    columns = {
        feature: [
            (terms[0], terms[0] if terms[0].startswith("METPO:") else f"METPO:2000103|{terms[0]}")
        ]
        for feature, terms in TABLE1_TRAITS.items()
    }
    matrix = pd.DataFrame({column[0][1]: values for column in columns.values()})
    stats = infographic.compute_trait_stats(
        matrix,
        study_ids=studies,
        selected_profile_count=2 * len(values),
        denominator_summary={
            "n_experiments": len(values),
            "n_nonzero_experiments": len(values),
            "n_zero_total_experiments": 0,
            "n_zero_total_excluded": 0,
        },
    )
    source = tmp_path / "native.tsv"
    if legacy:
        stats = stats.drop(columns=["pair_test_status", "study_test_status"])
    stats.to_csv(source, sep="\t", index=False)
    records = native_composition_records(
        source,
        columns=columns,
        q_threshold=0.05,
        strict_q_threshold=0.001,
        profile_sha256="fixture",
    )
    assert len(records) == 6
    for record in records:
        assert record.test_status == expected_state
        assert record.decision == record.strict_decision == "undecided"
        assert record.p_value is record.q_value is None
        assert record.interval_low is record.interval_high is None
        assert record.study_effect_mean == float(
            stats.iloc[0]["study_mean_increased_minus_decreased"]
        )
        assert record.family_size == 6
        if len(set(values)) > 1:
            assert record.native_pair_p_value is not None
            assert record.native_pair_q_value is not None
        else:
            assert record.native_pair_p_value is record.native_pair_q_value is None


def test_load_paired_inputs_maps_display_columns_to_stable_trait_ids(tmp_path: Path) -> None:
    matrix_path = tmp_path / "paired_matrix.tsv"
    metadata_path = tmp_path / "paired_metadata.tsv"
    trait_metadata_path = tmp_path / "trait_metadata.tsv"

    pd.DataFrame(
        {
            "oxygen\nMETPO:1": [1.0, 2.0],
            "nitrate reduction\npathways:nitrate_reduction": [0.0, 1.0],
        },
        index=["pair1", "pair2"],
    ).to_csv(matrix_path, sep="\t")
    pd.DataFrame(
        {
            "study_id": ["study1", "study2"],
            "increased_profile_ids": ["study1|exp1|inc", "study2|exp1|inc"],
            "decreased_profile_ids": ["study1|exp1|dec", "study2|exp1|dec"],
        },
        index=["pair1", "pair2"],
    ).to_csv(metadata_path, sep="\t")
    pd.DataFrame(
        {
            "trait_id": ["METPO:1", "METPO:2000103|pathways:nitrate_reduction"],
            "display_label": [
                "oxygen\nMETPO:1",
                "nitrate reduction\npathways:nitrate_reduction",
            ],
        },
        index=["METPO:1", "METPO:2000103|pathways:nitrate_reduction"],
    ).to_csv(trait_metadata_path, sep="\t")

    matrix, _metadata, trait_metadata = infographic.load_paired_inputs(
        matrix_path,
        metadata_path,
        trait_metadata_path,
    )

    assert list(matrix.columns) == ["METPO:1", "METPO:2000103|pathways:nitrate_reduction"]
    assert trait_metadata.loc["METPO:2000103|pathways:nitrate_reduction", "trait_label"] == (
        "nitrate reduction pathways:nitrate_reduction"
    )


def test_compute_trait_stats_writes_stable_ids_and_annotation_depth() -> None:
    matrix = pd.DataFrame(
        {
            "METPO:1": [1.0, 2.0, 3.0, 4.0],
            "METPO:2000103|pathways:nitrate_reduction": [0.0, 0.0, 0.0, 0.0],
        }
    )
    trait_metadata = pd.DataFrame(
        {
            "trait_id": ["METPO:1", "METPO:2000103|pathways:nitrate_reduction"],
            "trait_label": ["oxygen METPO:1", "nitrate reduction pathways:nitrate_reduction"],
            "selected_trait_rank": [1, 2],
            "profile_support": [4, 0],
            "total_signal": [10.0, 0.0],
            "composition_variance": [0.2, 0.0],
            "selection_score": [0.5, 0.0],
        },
        index=["METPO:1", "METPO:2000103|pathways:nitrate_reduction"],
    )
    annotation_counts = pd.Series({"METPO:1": 4, "METPO:2000103|pathways:nitrate_reduction": 0})

    stats = infographic.compute_trait_stats(
        matrix,
        trait_metadata,
        annotation_counts,
        stratum_kind="pooled",
        stratum="pooled",
        selected_profile_count=8,
    ).set_index("trait")

    assert "nitrate reduction\npathways:nitrate_reduction" not in stats.index
    assert stats.loc["METPO:1", "trait_label"] == "oxygen METPO:1"
    assert stats.loc["METPO:1", "n_profiles_where_annotated"] == 4
    assert stats.loc["METPO:1", "annotated_profile_fraction"] == 0.5
    assert stats.loc["METPO:2000103|pathways:nitrate_reduction", "q_value"] == 1.0
    assert stats.loc["METPO:1", "pair_count"] == 4
    assert stats.loc["METPO:1", "study_count"] == 4
    assert "ci95_lower" in stats.columns
    assert "ci95_upper" in stats.columns
    assert "study_ci95_lower" in stats.columns
    assert "study_ci95_upper" in stats.columns
    assert "study_p_value" in stats.columns
    assert "study_q_value" in stats.columns


def test_write_stratified_trait_stats_emits_feces_excluded_and_body_site_files(
    tmp_path: Path,
) -> None:
    profile_matrix = pd.DataFrame(
        {
            "METPO:1": [10, 1, 9, 1, 0, 0, 1, 0],
            "METPO:2": [1, 10, 1, 9, 5, 1, 5, 1],
            "METPO:3": [0, 0, 5, 0, 10, 1, 9, 1],
        },
        index=[
            "study1|exp1|inc",
            "study1|exp1|dec",
            "study2|exp1|inc",
            "study2|exp1|dec",
            "study3|exp1|inc",
            "study3|exp1|dec",
            "study4|exp1|inc",
            "study4|exp1|dec",
        ],
    )
    pair_metadata = pd.DataFrame(
        {
            "study_id": ["study1", "study2", "study3", "study4"],
            "body_site_group": ["Feces", "Feces", "Oral", "Oral"],
            "increased_profile_ids": [
                "study1|exp1|inc",
                "study2|exp1|inc",
                "study3|exp1|inc",
                "study4|exp1|inc",
            ],
            "decreased_profile_ids": [
                "study1|exp1|dec",
                "study2|exp1|dec",
                "study3|exp1|dec",
                "study4|exp1|dec",
            ],
        },
        index=["pair1", "pair2", "pair3", "pair4"],
    )
    trait_labels = pd.DataFrame(
        {
            "trait_label": ["oxygen METPO:1", "nitrate METPO:2", "motility METPO:3"],
            "display_label": ["oxygen METPO:1", "nitrate METPO:2", "motility METPO:3"],
            "column_label": ["oxygen METPO:1", "nitrate METPO:2", "motility METPO:3"],
            "object_name": ["oxygen", "nitrate", "motility"],
            "object_category": ["biolink:OntologyClass"] * 3,
            "display_category": ["OntologyClass"] * 3,
        },
        index=["METPO:1", "METPO:2", "METPO:3"],
    )

    manifest = infographic.write_stratified_trait_stats(
        profile_matrix,
        pair_metadata,
        trait_labels,
        output_dir=tmp_path,
        min_pairs=2,
        max_traits=2,
        min_profile_support=1,
        max_profile_support_fraction=1.0,
    )

    assert set(manifest["filename_token"]) == {"feces", "feces_excluded", "oral"}
    for token in ["feces", "feces_excluded", "oral"]:
        output = tmp_path / f"trait_profile_paired_contrast_trait_stats_{token}.tsv"
        assert output.exists()
        stats = pd.read_csv(output, sep="\t")
        assert len(stats) == 2
        assert stats["trait"].str.startswith("METPO:").all()
        assert "selected_trait_rank" in stats.columns
        assert "n_profiles_where_annotated" in stats.columns
        assert "study_count" in stats.columns
        assert "study_mean_increased_minus_decreased" in stats.columns
        assert "ci95_lower" in stats.columns
        assert "ci95_upper" in stats.columns
        assert "study_ci95_lower" in stats.columns
        assert "study_ci95_upper" in stats.columns
        assert "study_p_value" in stats.columns
        assert "study_q_value" in stats.columns

    assert (tmp_path / "trait_profile_paired_contrast_trait_stats_manifest.tsv").exists()

    uncapped_dir = tmp_path / "uncapped"
    uncapped_manifest = infographic.write_stratified_trait_stats(
        profile_matrix,
        pair_metadata,
        trait_labels,
        output_dir=uncapped_dir,
        min_pairs=2,
        min_profile_support=1,
        max_profile_support_fraction=1.0,
    )

    assert set(uncapped_manifest["selected_trait_count"]) == {3}
    for token in ["feces", "feces_excluded", "oral"]:
        output = uncapped_dir / f"trait_profile_paired_contrast_trait_stats_{token}.tsv"
        stats = pd.read_csv(output, sep="\t")
        assert len(stats) == 3


def test_compute_trait_stats_study_aggregation_and_confidence_intervals() -> None:
    # 4 pairs across 3 studies: study1 has 2 pairs, study2 has 1 pair, study3 has 1 pair
    matrix = pd.DataFrame(
        {
            "METPO:1": [1.0, 3.0, 2.0, 4.0],
            "METPO:2": [-2.0, -4.0, -1.0, -3.0],
        },
        index=["p1", "p2", "p3", "p4"],
    )
    study_ids = pd.Series(["s1", "s1", "s2", "s3"], index=["p1", "p2", "p3", "p4"])
    trait_metadata = pd.DataFrame(
        {
            "trait_id": ["METPO:1", "METPO:2"],
            "trait_label": ["oxygen METPO:1", "anaerobic METPO:2"],
            "selected_trait_rank": [1, 2],
            "profile_support": [4, 4],
            "total_signal": [10.0, 10.0],
            "composition_variance": [0.2, 0.2],
            "selection_score": [0.5, 0.5],
        },
        index=["METPO:1", "METPO:2"],
    )
    annotation_counts = pd.Series({"METPO:1": 4, "METPO:2": 4})

    stats = infographic.compute_trait_stats(
        matrix,
        trait_metadata,
        annotation_counts,
        stratum_kind="pooled",
        stratum="pooled",
        selected_profile_count=8,
        study_ids=study_ids,
    ).set_index("trait")

    # Pair-level:
    assert stats.loc["METPO:1", "pair_count"] == 4
    assert np.isclose(stats.loc["METPO:1", "mean_increased_minus_decreased"], 2.5)
    # Pair std = sqrt(5/3) = 1.290994, SE = 0.645497, df = 3, t_crit = 3.182446
    # 95% CI: [0.44575, 4.55425]
    assert np.isclose(stats.loc["METPO:1", "ci95_lower"], 2.5 - 3.182446 * 0.645497, atol=1e-3)
    assert np.isclose(stats.loc["METPO:1", "ci95_upper"], 2.5 + 3.182446 * 0.645497, atol=1e-3)

    # Study-level:
    # s1 mean: (1+3)/2 = 2.0; s2: 2.0; s3: 4.0 -> values [2.0, 2.0, 4.0]
    # study count = 3
    assert stats.loc["METPO:1", "study_count"] == 3
    assert np.isclose(stats.loc["METPO:1", "study_mean_increased_minus_decreased"], 8.0 / 3.0)
    # study variance = 4/3, std = sqrt(4/3), SE = 2/3, t = 4.0, df = 2
    assert np.isclose(stats.loc["METPO:1", "study_t_statistic"], 4.0)
    # t_crit(df=2) = 4.30265
    assert np.isclose(
        stats.loc["METPO:1", "study_ci95_lower"],
        8.0 / 3.0 - 4.30265 * (2.0 / 3.0),
        atol=1e-3,
    )
    assert np.isclose(
        stats.loc["METPO:1", "study_ci95_upper"],
        8.0 / 3.0 + 4.30265 * (2.0 / 3.0),
        atol=1e-3,
    )
    assert stats.loc["METPO:1", "study_p_value"] < 0.06
    assert np.isfinite(stats.loc["METPO:1", "study_q_value"])
    assert np.isclose(
        stats.loc["METPO:1", "study_minus_log10_p"],
        -np.log10(stats.loc["METPO:1", "study_p_value"]),
    )


def test_pair_denominators_distinguish_zero_total_from_zero_contrast() -> None:
    # p1 has equal, positive profiles: a zero contrast is still an eligible pair.
    # p2 is MediaDive-only on one side; s1 survives through p1. s2 has no
    # nonzero pairs because p3 has a genuinely empty side.
    profiles = pd.DataFrame(
        {
            "METPO:1": [2, 2, 0, 3, 1, 0],
            "METPO:2": [0, 0, 0, 1, 4, 0],
            "mediadive.medium:1": [0, 0, 100, 0, 0, 0],
        },
        index=["s1|e1|inc", "s1|e1|dec", "s1|e2|inc", "s1|e2|dec", "s2|e1|inc", "s2|e1|dec"],
    )
    metadata = pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s2"],
            "increased_profile_ids": ["s1|e1|inc", "s1|e2|inc", "s2|e1|inc"],
            "decreased_profile_ids": ["s1|e1|dec", "s1|e2|dec", "s2|e1|dec"],
        },
        index=["p1", "p2", "p3"],
    )
    audit = infographic.pair_denominator_summary(profiles, metadata)
    assert audit["n_experiments"] == 3
    assert audit["n_nonzero_experiments"] == 1
    assert audit["n_zero_total_experiments"] == 2
    assert audit["n_zero_total_excluded"] == 0
    assert audit["n_nonzero_studies"] == 1
    assert audit["n_studies_without_nonzero_experiments"] == 1
    assert audit["n_zero_total_profiles"] == 2
    assert audit["retained_trait_count"] == 2
    assert audit["zero_total_policy"] == "included_in_historical_native_composition"

    matrix = infographic.build_pair_matrix_from_profiles(
        infographic.build_zscore_matrix(profiles, pd.Index(["METPO:1", "METPO:2"])),
        metadata,
    )
    assert matrix.loc["p1"].eq(0).all()
    stats = infographic.compute_trait_stats(
        matrix, study_ids=metadata["study_id"], denominator_summary=audit
    )
    assert stats["pair_count"].eq(3).all()
    assert stats["n_nonzero_experiments"].eq(1).all()
    assert stats["study_count"].eq(2).all()
    assert stats["correction_family_size"].eq(2).all()
    assert stats["test_method"].eq("two_sided_one_sample_t_test").all()
    with pytest.raises(ValueError, match="same pairs"):
        infographic.compute_trait_stats(matrix.iloc[:2], denominator_summary=audit)
    with pytest.raises(ValueError, match="missing raw profiles"):
        infographic.pair_denominator_summary(profiles.iloc[:-1], metadata)


def test_stats_only_writes_pooled_table_with_native_denominators(tmp_path: Path) -> None:
    profiles_path = tmp_path / "profiles.tsv"
    metadata_path = tmp_path / "pairs.tsv"
    labels_path = tmp_path / "labels.tsv"
    output_dir = tmp_path / "output"
    pd.DataFrame(
        {
            "study_id": ["s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4"],
            "experiment_id": ["e1"] * 8,
            "inc_dec": ["inc", "dec"] * 4,
            "METPO:1": [10, 1, 8, 1, 7, 3, 0, 0],
            "METPO:2": [1, 10, 2, 9, 2, 8, 0, 0],
        }
    ).to_csv(profiles_path, sep="\t", index=False)
    pd.DataFrame(
        {
            "study_id": ["s1", "s2", "s3", "s4"],
            "increased_profile_ids": [f"s{i}|e1|inc" for i in range(1, 5)],
            "decreased_profile_ids": [f"s{i}|e1|dec" for i in range(1, 5)],
        }
    ).to_csv(metadata_path, sep="\t")
    pd.DataFrame(
        {
            "column_id": ["METPO:1", "METPO:2"],
            "column_label": ["oxygen", "motility"],
            "object_name": ["oxygen", "motility"],
            "object_category": ["biolink:OntologyClass"] * 2,
        }
    ).to_csv(labels_path, sep="\t", index=False)
    infographic.main(
        [
            "--stats-only",
            "--profiles",
            str(profiles_path),
            "--pair-metadata",
            str(metadata_path),
            "--trait-labels",
            str(labels_path),
            "--output-dir",
            str(output_dir),
        ]
    )
    output_path = output_dir / "trait_profile_paired_contrast_trait_stats_pooled.tsv"
    assert list(output_dir.iterdir()) == [output_path]
    stats = pd.read_csv(output_path, sep="\t")
    assert len(stats) == 2
    assert stats["pair_count"].eq(4).all()
    assert stats["profile_count"].eq(8).all()
    assert stats["n_nonzero_experiments"].eq(3).all()
    assert stats["n_zero_total_experiments"].eq(1).all()
    assert stats["n_zero_total_excluded"].eq(0).all()
    assert stats["n_nonzero_studies"].eq(3).all()
    assert stats["correction_family_size"].eq(len(stats)).all()
    assert stats["profile_sha256"].eq(file_sha256(profiles_path)).all()
    assert stats["pair_metadata_sha256"].eq(file_sha256(metadata_path)).all()


def test_plot_top_traits_renders_with_study_level_stats() -> None:
    trait_stats = pd.DataFrame(
        {
            "trait": ["METPO:1", "METPO:2"],
            "trait_label": ["aerobic METPO:1", "anaerobic METPO:2"],
            "study_minus_log10_p": [8.5, 5.2],
            "study_mean_increased_minus_decreased": [0.21, -0.15],
            "study_ci95_lower": [0.14, -0.22],
            "study_ci95_upper": [0.28, -0.08],
            "minus_log10_p": [12.0, 8.0],
            "mean_increased_minus_decreased": [0.25, -0.18],
            "ci95_lower": [0.18, -0.25],
            "ci95_upper": [0.32, -0.11],
        }
    )
    fig, ax = plt.subplots()
    try:
        infographic.plot_top_traits(ax, trait_stats)
        assert "study-level p" in ax.get_xlabel()
        texts = [t.get_text() for t in ax.texts]
        assert any("[+0.14, +0.28]" in t for t in texts)
        assert any("[-0.22, -0.08]" in t for t in texts)
        assert all(t.get_bbox_patch() is not None for t in ax.texts)
    finally:
        plt.close(fig)


def test_body_site_panel_keeps_labels_and_bubbles_apart() -> None:
    """Panel C at its size in the infographic: nothing overlaps (#612, #620)."""
    summary = pd.DataFrame(
        {
            "body_site_group": [
                "Throat",
                "Uterine cervix",
                "Feces",
                "Subgingival dental plaque",
                "Mouth",
                "Saliva",
                "Vagina",
                "Oral cavity",
                "Nasopharynx",
            ],
            "pair_count": [28, 35, 1460, 35, 34, 158, 48, 27, 25],
            "absolute_contrast_load": [0.65, 0.62, 0.60, 0.58, 0.56, 0.55, 0.48, 0.44, 0.40],
            "net_increased_minus_decreased": [
                0.08,
                0.04,
                0.02,
                0.10,
                -0.13,
                -0.07,
                -0.04,
                -0.18,
                0.01,
            ],
        }
    )
    # Panel C is the bottom-left cell of create_infographic's grid.
    fig = plt.figure(figsize=infographic.FIGURE_SIZE, dpi=infographic.FIGURE_DPI)
    grid = fig.add_gridspec(2, 2, **infographic.PANEL_GRID)
    ax = fig.add_subplot(grid[1, 0])
    try:
        infographic.plot_body_site_summary(ax, summary)
        fig.canvas.draw()
        boxes = [text.get_window_extent() for text in ax.texts]
        assert len(boxes) == 2 * len(summary) + 1
        touching_text = [
            (ax.texts[a].get_text(), ax.texts[b].get_text())
            for a, b in combinations(range(len(boxes)), 2)
            if boxes[a].overlaps(boxes[b])
        ]
        assert touching_text == []
        bubbles = [
            patch for patch in ax.patches if isinstance(patch, Circle) and patch.get_zorder() == 4
        ]
        assert len(bubbles) == len(summary)
        touching_bubbles = [
            (a.center, b.center)
            for a, b in combinations(bubbles, 2)
            if math.dist(a.center, b.center) < a.radius + b.radius
        ]
        assert touching_bubbles == []
        # Leader lines run above the bubbles, so a line crossing one stays visible (#712).
        leaders = [line for line in ax.lines if line.get_zorder() > bubbles[0].get_zorder()]
        assert len(leaders) == len(summary)
        # Throat sits above Subgingival dental plaque under leader_order_key (#718).
        label_y = {text.get_text(): text.get_position()[1] for text in ax.texts}
        assert label_y["Throat"] > label_y["Subgingival dental plaque"]
    finally:
        plt.close(fig)


@pytest.fixture
def pooled_core(tmp_path: Path) -> CoreBundle:
    """A real core bundle with enough supported traits plus excluded/zero cohorts."""
    source = REPO_ROOT / "tests" / "fixtures" / "workflow"
    raw = pd.read_csv(source / "full_dump.csv", comment="#", dtype=str).head(7)
    raw.loc[raw["BSDB ID"] == "bsdb:1/1/2", "NCBI Taxonomy IDs"] = "101"
    raw.loc[raw["BSDB ID"] == "bsdb:1/1/2", "MetaPhlAn taxon names"] = "s__Examplegenus_alpha"
    frames = []
    for offset in range(0, 8, 2):
        frame = raw.copy()
        frame["BSDB ID"] = frame["BSDB ID"].map(
            lambda value, offset=offset: (
                (f"bsdb:{int(value.split('/')[0].removeprefix('bsdb:')) + offset}/")
                + "/".join(value.split("/")[1:])
            )
        )
        frame["Study"] = frame["Study"].map(
            lambda value, offset=offset: f"Study {int(value.removeprefix('Study ')) + offset}"
        )
        frames.append(frame)
    dump = tmp_path / "dump.csv"
    dump.write_text("# BugSigDB fixture-551, License: CC0 synthetic test data\n")
    pd.concat(frames, ignore_index=True).to_csv(dump, mode="a", index=False)
    edges = tmp_path / "edges.tsv"
    edges.write_text(
        (source / "merged-kg_edges.tsv").read_text()
        + "motile\tNCBITaxon:101\tbiolink:has_phenotype\tMETPO:1000700\tfixture:annotation\n"
    )
    nodes = tmp_path / "nodes.tsv"
    nodes.write_text(
        (source / "merged-kg_nodes.tsv").read_text()
        + "METPO:1000700\tbiolink:PhenotypicQuality\tmotile\n"
    )
    run = execute_workflow(
        REPO_ROOT / "workflows" / "bugsigdb_core.yaml",
        {"raw_dump": dump, "kg_edges": edges, "kg_nodes": nodes},
        tmp_path / "runs",
        OPERATIONS,
        validator=validate_core_product,
        repository_root=REPO_ROOT,
    )
    return open_core_bundle(Path(run.manifest_path))


def test_stats_only_core_derives_cohort_and_records_exact_outputs(
    pooled_core: CoreBundle, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Core mode must derive human pairs, retain zero totals and ignore heatmap defaults."""
    for name in ("PROFILE_MATRIX_PATH", "PAIRED_METADATA_PATH", "TRAIT_COLUMNS_PATH"):
        monkeypatch.setattr(infographic, name, tmp_path / "no-legacy-input.tsv")
    output = tmp_path / "analysis"
    infographic.main(
        [
            "--stats-only",
            "--core-run",
            str(pooled_core.manifest),
            "--profiles",
            str(pooled_core.path("count_profiles")),
            "--trait-labels",
            str(pooled_core.path("feature_metadata")),
            "--output-dir",
            str(output),
        ]
    )
    pair_path = output / "trait_profile_paired_contrast_pair_metadata.tsv"
    pairs = pd.read_csv(pair_path, sep="\t", index_col=0)
    assert len(pairs) == 8
    assert pairs["study_id"].nunique() == 4
    assert pairs["host_species"].eq("Homo sapiens").all()
    assert pairs["condition"].eq("fixture condition").all()
    assert pairs["input_profile_rows"].eq(28).all()
    assert pairs["nonhuman_profile_rows_excluded"].eq(8).all()
    assert pairs["unpaired_profile_rows_excluded"].eq(4).all()
    assert pairs["missing_metadata_profile_rows_excluded"].eq(0).all()
    assert pairs["profile_sha256"].eq(file_sha256(pooled_core.path("count_profiles"))).all()
    assert pairs["experiments_sha256"].eq(file_sha256(pooled_core.path("experiments"))).all()
    stats_path = output / "trait_profile_paired_contrast_trait_stats_pooled.tsv"
    result = pd.read_csv(stats_path, sep="\t")
    assert len(result) == 2
    assert result["trait_label"].str.contains("motile").any()
    assert result["pair_count"].eq(8).all()
    assert result["study_count"].eq(4).all()
    assert result["n_nonzero_experiments"].eq(4).all()
    assert result["n_zero_total_experiments"].eq(4).all()
    assert result["n_zero_total_excluded"].eq(0).all()
    assert result["profile_count"].eq(16).all()
    assert result["pair_metadata_sha256"].eq(file_sha256(pair_path)).all()
    assert result["experiments_sha256"].eq(file_sha256(pooled_core.path("experiments"))).all()
    record = verify_consumer_record(output / "core_consumption.json")
    assert record.consumer == "native-pooled-composition-stats"
    assert record.feature_policy == "exclude_mediadive"
    # Coverage describes all source profiles; pair metadata documents the
    # narrower human analysis cohort without mislabelling source coverage.
    assert record.cohort == "all"
    assert record.coverage_scope == "input_cohort_before_analysis_specific_exclusions"
    assert {Path(item.path).name for item in record.outputs} == {
        pair_path.name,
        stats_path.name,
        "core_consumption.coverage.tsv",
    }
    pair_path.write_text(pair_path.read_text().replace("fixture condition", "modified label"))
    with pytest.raises(ValueError, match="Consumer output changed"):
        verify_consumer_record(output / "core_consumption.json")


@pytest.mark.parametrize("argument", ["--profiles", "--trait-labels"])
def test_stats_only_core_rejects_conflicting_inputs(
    pooled_core: CoreBundle, tmp_path: Path, argument: str
) -> None:
    output = tmp_path / "analysis"
    with pytest.raises(ValueError, match="conflicts with explicit"):
        infographic.main(
            [
                "--stats-only",
                "--core-run",
                str(pooled_core.manifest),
                argument,
                str(tmp_path / "legacy.tsv"),
                "--output-dir",
                str(output),
            ]
        )
    assert not output.exists()


@pytest.mark.parametrize("arguments", [[], ["--stats-only", "--pair-metadata", "old.tsv"]])
def test_core_mode_rejects_figures_and_external_pair_metadata(arguments: list[str]) -> None:
    with pytest.raises(SystemExit):
        infographic.main(["--core-run", "unused-run.json", *arguments])


@pytest.mark.parametrize("collision", ["core_directory", "symlink", "hardlink", "sidecar"])
def test_stats_only_core_preflight_protects_inputs_and_sidecars(
    pooled_core: CoreBundle, tmp_path: Path, collision: str
) -> None:
    protected = pooled_core.path("count_profiles")
    original = protected.read_bytes()
    output = tmp_path / "analysis"
    output.mkdir()
    if collision == "core_directory":
        output = protected.parent
    else:
        target = output / "trait_profile_paired_contrast_pair_metadata.tsv"
        if collision == "symlink":
            target.symlink_to(protected)
        elif collision == "hardlink":
            target.hardlink_to(protected)
        else:
            sidecar = output / "core_consumption.json"
            sidecar.write_text("previous consumer record\n")
            target.hardlink_to(sidecar)
    with pytest.raises(ValueError, match=r"overlap|collid"):
        infographic.main(
            [
                "--stats-only",
                "--core-run",
                str(pooled_core.manifest),
                "--output-dir",
                str(output),
            ]
        )
    assert protected.read_bytes() == original
    assert not (output / "trait_profile_paired_contrast_trait_stats_pooled.tsv").exists()
    if collision == "sidecar":
        assert (output / "core_consumption.json").read_text() == "previous consumer record\n"
