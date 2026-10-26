"""Signed direction contrasts for a declared trait list (#1054).

Synthetic profiles stand in for a core run; nothing here reads data/derived or
renders a figure.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from mtp.trait_lists import load_trait_list

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
with patch.object(sys, "path", [str(SCRIPTS_DIR), *sys.path]):
    sankeys = importlib.import_module("create_six_shared_trait_sankeys")
    heatmaps = importlib.import_module("create_six_shared_trait_cytoscape_heatmaps")

SIX = [trait.trait_id for trait in sankeys.SHARED_TRAITS]
EXTRA = {
    "METPO:1000605": "facultatively anaerobic",
    "METPO:1000603": "anaerobic",
    "METPO:2000103|GO:0045493": "xylan catabolic process",
}
ABSENT_FROM_PROFILES = "METPO:1000447"
MEDIA = "METPO:2000517|mediadive.medium:1"
CONDITIONS = ["Parkinson's disease", "Obesity", "Colorectal cancer", "HIV infection"]
SITES = [("Feces", "Feces"), ("Saliva", "Oral cavity / dental")]
# Columns compared between the six-trait path and a trait-list run: every count and
# test that does not depend on which other traits share the BH family.
LIST_INDEPENDENT = [
    "view_key",
    "condition",
    "body_site_group",
    "trait_id",
    "study_count",
    "increased_positive_profile_count",
    "increased_profile_count",
    "decreased_positive_profile_count",
    "decreased_profile_count",
    "increased_positive_fraction",
    "decreased_positive_fraction",
    "increased_minus_decreased_positive_fraction",
    "fisher_p_value",
    "passes_support",
]


@pytest.fixture
def run_inputs(tmp_path: Path) -> dict[str, Path]:
    """Profiles, matching profile metadata and a trait columns table."""
    rng = np.random.default_rng(1054)
    traits = [*SIX, *EXTRA]
    profiles: list[dict[str, object]] = []
    metadata: list[dict[str, object]] = []
    study = 0
    for condition in CONDITIONS:
        for site, group in SITES:
            for _ in range(3):
                study += 1
                for experiment in (1, 2):
                    for inc_dec in ("inc", "dec"):
                        counts = rng.poisson(rng.uniform(1, 30, size=len(traits)))
                        keys = {
                            "study_id": f"bsdb:{study}",
                            "experiment_id": f"bsdb:{study}-{experiment}",
                            "inc_dec": inc_dec,
                        }
                        profiles.append(
                            {
                                **keys,
                                **dict(zip(traits, counts.tolist(), strict=True)),
                                MEDIA: int(rng.poisson(40)),
                            }
                        )
                        metadata.append(
                            {
                                "profile_id": "|".join(keys.values()),
                                **keys,
                                "total_trait_signal": int(counts.sum()),
                                "condition": condition,
                                "body_site": site,
                                "body_site_group": group,
                            }
                        )
    paths = {
        "profiles": tmp_path / "study_trait_profiles_wide.tsv",
        "metadata": tmp_path / "profile_metadata.tsv",
        "columns": tmp_path / "taxa_trait_matrix_columns.tsv",
    }
    pd.DataFrame(profiles).to_csv(paths["profiles"], sep="\t", index=False)
    pd.DataFrame(metadata).set_index("profile_id").to_csv(paths["metadata"], sep="\t")
    labels = {trait.trait_id: trait.label for trait in sankeys.SHARED_TRAITS} | EXTRA
    pd.DataFrame(
        {
            "column_id": [*traits, ABSENT_FROM_PROFILES, MEDIA],
            "object_name": [*(labels[trait] for trait in traits), "thermophilic", "medium"],
            "object_category": ["biolink:PhenotypicQuality"] * (len(traits) + 1)
            + ["METPO:1004005"],
        }
    ).to_csv(paths["columns"], sep="\t", index=False)
    return paths


def declare(path: Path, trait_ids: Sequence[str]) -> Path:
    path.write_text("\n".join(trait_ids) + "\n", encoding="utf-8")
    return path


def run_trait_list(run_inputs: dict[str, Path], output: Path, list_path: Path) -> pd.DataFrame:
    contrasts = sankeys.write_trait_list_contrasts(
        load_trait_list(list_path),
        output,
        profiles_path=run_inputs["profiles"],
        metadata_path=run_inputs["metadata"],
        columns_path=run_inputs["columns"],
        code_commit="abc1234",
    )
    assert isinstance(contrasts, pd.DataFrame)
    return contrasts


def six_trait_path_contrasts(run_inputs: dict[str, Path], output: Path) -> pd.DataFrame:
    """The default six-trait computation: no trait argument anywhere, TSV round trip."""
    output.mkdir()
    metadata = sankeys.load_profile_metadata(run_inputs["metadata"])
    raw = sankeys.load_raw_trait_counts(sankeys.SHARED_TRAITS, metadata, run_inputs["profiles"])
    zscores = sankeys.build_zscore_matrix(raw, metadata)
    edges = {
        "grouped": sankeys.build_trait_direction_tables(metadata, zscores)[0],
        "individual_conditions": sankeys.build_individual_condition_direction_tables(
            metadata, zscores
        )[0],
    }
    trait_order = pd.DataFrame(
        {
            "trait_id": SIX,
            "trait_label": [trait.label for trait in sankeys.SHARED_TRAITS],
            "trait_rank": range(len(SIX)),
        }
    )
    tables: list[pd.DataFrame] = []
    for mode in heatmaps.CONTRAST_MODES:
        frames = []
        for direction in ("increased", "decreased"):
            path = output / f"{mode.key}_{direction}.tsv"
            mode_edges = edges[mode.key]
            mode_edges[mode_edges["focus_direction"] == direction].to_csv(
                path, sep="\t", index=False
            )
            frames.append(heatmaps.read_tsv(path))
        tables.append(heatmaps.prepare_contrast_cells(*frames, mode, trait_order))
    return pd.concat(tables, ignore_index=True)


def test_six_trait_list_reproduces_the_six_trait_contrast_cells(
    tmp_path: Path, run_inputs: dict[str, Path]
) -> None:
    """Declaring the default six in default order gives the default cells, BH included."""
    declared = run_trait_list(run_inputs, tmp_path / "declared", declare(tmp_path / "six.txt", SIX))
    expected = six_trait_path_contrasts(run_inputs, tmp_path / "default")

    assert len(declared) == len(expected) > 0
    key = ["view_key", "condition_body_site_pair", "trait_id"]
    pd.testing.assert_frame_equal(
        declared.sort_values(key).reset_index(drop=True),
        expected.sort_values(key).reset_index(drop=True),
    )
    written = pd.read_csv(
        tmp_path / "declared" / heatmaps.TRAIT_LIST_CONTRAST_TABLE_FILENAME, sep="\t"
    )
    assert len(written) == len(declared)


def test_longer_trait_list_keeps_shared_trait_cells_and_adds_declared_traits(
    tmp_path: Path, run_inputs: dict[str, Path]
) -> None:
    declared_ids = ["METPO:1000605", *SIX, "METPO:1000603", "METPO:2000103|GO:0045493"]
    six = run_trait_list(run_inputs, tmp_path / "six", declare(tmp_path / "six.txt", SIX))
    longer = run_trait_list(
        run_inputs, tmp_path / "longer", declare(tmp_path / "longer.txt", declared_ids)
    )

    assert set(longer["trait_id"]) == set(declared_ids)
    ranks = longer.drop_duplicates("trait_id").set_index("trait_id")["trait_rank"]
    assert ranks.loc[declared_ids].tolist() == list(range(len(declared_ids)))
    labels = longer.drop_duplicates("trait_id").set_index("trait_id")["trait_label"]
    assert labels["METPO:1000605"] == "facultatively anaerobic"
    key = ["view_key", "condition", "body_site_group", "trait_id"]
    shared = longer[longer["trait_id"].isin(SIX)][LIST_INDEPENDENT].sort_values(key)
    pd.testing.assert_frame_equal(
        shared.reset_index(drop=True),
        six[LIST_INDEPENDENT].sort_values(key).reset_index(drop=True),
    )
    grouped_feces = longer[
        longer["view_key"].eq("grouped_direction_contrast") & longer["body_site_group"].eq("Feces")
    ]
    assert set(grouped_feces["trait_id"]) == set(declared_ids)


def test_trait_list_provenance_records_the_list_inputs_and_table(
    tmp_path: Path, run_inputs: dict[str, Path]
) -> None:
    declared_ids = ["METPO:1000605", SIX[1]]
    list_path = declare(tmp_path / "headline.txt", declared_ids)
    output = tmp_path / "output"
    contrasts = run_trait_list(run_inputs, output, list_path)

    provenance = json.loads((output / heatmaps.TRAIT_LIST_PROVENANCE_FILENAME).read_text())

    def sha256(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    assert provenance["trait_list_file"] == "headline.txt"
    assert provenance["trait_list_sha256"] == sha256(list_path)
    assert provenance["trait_ids"] == declared_ids
    assert provenance["profiles_sha256"] == sha256(run_inputs["profiles"])
    assert provenance["profile_metadata_sha256"] == sha256(run_inputs["metadata"])
    assert provenance["trait_columns_sha256"] == sha256(run_inputs["columns"])
    assert provenance["code_commit"] == "abc1234"
    table = output / heatmaps.TRAIT_LIST_CONTRAST_TABLE_FILENAME
    assert provenance["contrast_table_sha256"] == sha256(table)
    assert list(provenance["cells_per_trait"]) == declared_ids
    assert sum(provenance["cells_per_trait"].values()) == len(contrasts)
    assert sum(provenance["cells_per_mode"].values()) == len(contrasts)
    declared_ranks = {trait_id: rank for rank, trait_id in enumerate(declared_ids)}
    for mode in heatmaps.CONTRAST_MODES:
        for direction, view_key in (
            ("increased", mode.increased_view_key),
            ("decreased", mode.decreased_view_key),
        ):
            edges = pd.read_csv(output / heatmaps.trait_list_edges_filename(view_key), sep="\t")
            assert set(edges["trait_id"]) == set(declared_ids)
            assert edges["focus_direction"].eq(direction).all()
            assert edges["trait_rank"].tolist() == edges["trait_id"].map(declared_ranks).tolist()
    assert not any(output.glob("*.png"))


@pytest.mark.parametrize(
    ("trait_id", "message"),
    [
        ("METPO:9999999", "absent from the trait columns table"),
        (ABSENT_FROM_PROFILES, "missing from the full raw matrix"),
        (MEDIA, "missing from the full raw matrix"),
    ],
)
def test_trait_list_rejects_traits_the_run_cannot_score(
    tmp_path: Path, run_inputs: dict[str, Path], trait_id: str, message: str
) -> None:
    list_path = declare(tmp_path / "traits.txt", [SIX[0], trait_id])
    with pytest.raises(ValueError, match=message):
        run_trait_list(run_inputs, tmp_path / "output", list_path)
    assert not (tmp_path / "output" / heatmaps.TRAIT_LIST_CONTRAST_TABLE_FILENAME).exists()


def test_contrast_writer_rejects_order_drift_and_traits_without_cells(tmp_path: Path) -> None:
    list_path = declare(tmp_path / "traits.txt", ["METPO:1", "METPO:2"])
    trait_list = load_trait_list(list_path)
    edges = pd.DataFrame(
        {
            "condition": ["Condition A"],
            "body_site_group": ["Feces"],
            "trait_id": ["METPO:1"],
            "trait_label": ["one"],
            "focus_direction": ["increased"],
            "study_count": [3],
            "focus_profile_count": [10],
            "focus_enriched_profile_count": [6],
            "comparator_profile_count": [10],
            "comparator_enriched_profile_count": [2],
        }
    )
    for mode in heatmaps.CONTRAST_MODES:
        for view_key in (mode.increased_view_key, mode.decreased_view_key):
            edges.to_csv(
                tmp_path / heatmaps.trait_list_edges_filename(view_key), sep="\t", index=False
            )
    order = pd.DataFrame(
        {"trait_id": ["METPO:2", "METPO:1"], "trait_label": ["two", "one"], "trait_rank": [0, 1]}
    )
    with pytest.raises(ValueError, match="declared order"):
        heatmaps.write_trait_list_direction_contrasts(tmp_path, order, trait_list, {})
    with pytest.raises(ValueError, match=r"no contrast cell: \['METPO:2'\]"):
        heatmaps.write_trait_list_direction_contrasts(
            tmp_path, order.iloc[::-1].reset_index(drop=True), trait_list, {}
        )
    assert not (tmp_path / heatmaps.TRAIT_LIST_CONTRAST_TABLE_FILENAME).exists()


def test_cli_writes_trait_list_contrasts_and_requires_an_output_dir(
    tmp_path: Path, run_inputs: dict[str, Path], capsys: pytest.CaptureFixture[str]
) -> None:
    list_path = declare(tmp_path / "traits.txt", list(EXTRA))
    argv = [
        "--trait-list",
        str(list_path),
        "--profiles",
        str(run_inputs["profiles"]),
        "--profile-metadata",
        str(run_inputs["metadata"]),
        "--columns",
        str(run_inputs["columns"]),
        "--code-commit",
        "abc1234",
    ]
    with patch.object(sankeys, "main") as six_trait_main, pytest.raises(SystemExit):
        sankeys.cli(argv)
    assert "--trait-list requires --output-dir" in capsys.readouterr().err

    output = tmp_path / "output"
    with patch.object(sankeys, "main") as six_trait_main:
        sankeys.cli([*argv, "--output-dir", str(output), "--allow-temporary-output"])
    six_trait_main.assert_not_called()
    contrasts = pd.read_csv(output / heatmaps.TRAIT_LIST_CONTRAST_TABLE_FILENAME, sep="\t")
    assert set(contrasts["trait_id"]) == set(EXTRA)
    provenance = json.loads((output / heatmaps.TRAIT_LIST_PROVENANCE_FILENAME).read_text())
    assert provenance["code_commit"] == "abc1234"


@pytest.mark.parametrize(
    "argv",
    [
        ["--profiles", "profiles.tsv"],
        ["--output-dir", "output"],
        ["--code-commit", ""],
    ],
)
def test_trait_list_inputs_never_reach_the_six_trait_run(
    argv: list[str], capsys: pytest.CaptureFixture[str]
) -> None:
    """Trait-list inputs are rejected without --trait-list, an empty commit included (#1083)."""
    with patch.object(sankeys, "main") as six_trait_main, pytest.raises(SystemExit):
        sankeys.cli(argv)
    six_trait_main.assert_not_called()
    assert "apply only with --trait-list" in capsys.readouterr().err


@pytest.mark.parametrize("relative", [".", "nested", "nested/deeper"])
def test_trait_list_refuses_the_six_trait_directory_and_anything_inside_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, relative: str
) -> None:
    """The six-trait run clears its directory and stops on nested ones (#1083)."""
    six_trait_dir = tmp_path / "six-shared-traits"
    six_trait_dir.mkdir()
    monkeypatch.setattr(sankeys, "OUTPUT_DIR", six_trait_dir)
    list_path = declare(tmp_path / "traits.txt", SIX)
    with pytest.raises(ValueError, match="outside the six-trait one"):
        sankeys.write_trait_list_contrasts(load_trait_list(list_path), six_trait_dir / relative)
    assert list(six_trait_dir.iterdir()) == []
