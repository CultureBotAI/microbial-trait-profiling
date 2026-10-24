"""Canonical population, data-derived labels and provenance regressions (#820/#821)."""

import importlib.util
import inspect
import json
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.text import Annotation

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/create_taxa_coverage_histogram.py"
SPEC = importlib.util.spec_from_file_location("taxa_coverage_histogram", SCRIPT)
assert SPEC is not None
assert SPEC.loader is not None
coverage = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(coverage)


def readme_regeneration_argv(text: str) -> list[str]:
    """Read the runnable shell block rather than searching for option substrings."""
    if "```sh\n" in text:
        command = text.split("```sh\n", 1)[1].split("\n```", 1)[0]
    else:
        # Keep the historical inline command readable for failing-first replay proof.
        command = "just taxa-coverage" + text.split("`just taxa-coverage", 1)[1].split("`", 1)[0]
    return shlex.split(command)


@pytest.mark.parametrize("scope", ["off", "phenotype"])
def test_readme_command_replays_explicit_inputs_and_options(
    inputs: dict[str, Path], tmp_path: Path, scope: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = {}
    for name, path in inputs.items():
        paths[name] = path.rename(tmp_path / f"input's {name}.tsv")
    raw = tmp_path / "raw dump's data.csv"
    raw.write_text("# BugSigDB fixture-970, License: CC0\nStudy\nfixture\n")
    sidecar = tmp_path / "KG's provenance.json"
    sidecar.write_text(json.dumps({"release": "fixture-release", "files": {}}))
    output = tmp_path / f"{scope}'s coverage export"
    arguments = [
        str(SCRIPT),
        "--strain-rollup",
        scope,
        "--tables-only",
        "--output-dir",
        str(output),
        "--raw-dump",
        str(raw),
        "--kg-provenance",
        str(sidecar),
        "--bugsigdb-vintage",
        "declared fixture vintage",
    ]
    for name, path in paths.items():
        arguments.extend([f"--{name}", str(path)])
    monkeypatch.setattr(sys, "argv", arguments)
    coverage.main()
    replay = readme_regeneration_argv((output / "README.md").read_text())
    assert replay[:2] == ["just", "taxa-coverage"]
    captured = {}
    monkeypatch.setattr(coverage, "run_coverage", lambda **kwargs: captured.update(kwargs))
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), *replay[2:]])
    coverage.main()
    for name, path in paths.items():
        assert Path(captured[name]).resolve() == path.resolve(), (
            f"Regeneration changed input: {name}"
        )
    assert captured["strain_rollup"] == scope
    assert Path(captured["output_dir"]).resolve() == output.resolve()
    assert captured["tables_only"] is True
    assert Path(captured["raw_dump"]).resolve() == raw.resolve()
    assert Path(captured["kg_provenance"]).resolve() == sidecar.resolve()
    assert captured["bugsigdb_vintage"] == "declared fixture vintage"


@pytest.mark.parametrize("scope", ["off", "phenotype"])
def test_scope_is_explicitly_declared_in_readme_and_provenance(
    inputs: dict[str, Path], tmp_path: Path, scope: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The same one-feature inputs cannot identify upstream strain rollup (#959)."""
    output = tmp_path / scope
    arguments = [
        str(SCRIPT),
        "--output-dir",
        str(output),
        "--tables-only",
        "--strain-rollup",
        scope,
    ]
    for name, path in inputs.items():
        arguments.extend([f"--{name}", str(path)])
    monkeypatch.setattr(sys, "argv", arguments)
    coverage.main()
    record = json.loads((output / "provenance.json").read_text())
    assert record["strain_rollup"] == scope
    assert record["strain_rollup_evidence"]["status"] == "declared_unverified"
    assert record["inputs"]["matrix"]["sha256"] == coverage.file_sha256(inputs["matrix"])
    assert record["inputs"]["profiles"]["sha256"] == coverage.file_sha256(inputs["profiles"])
    readme = (output / "README.md").read_text()
    assert f"Declared upstream strain-rollup scope: `{scope}`" in readme
    command = readme_regeneration_argv(readme)
    assert command[command.index("--strain-rollup") + 1] == scope
    assert "do not verify the upstream rollup scope" in readme


@pytest.mark.parametrize("arguments", [[], ["--strain-rollup", "unknown"]])
def test_cli_requires_valid_scope_before_reading_inputs(
    monkeypatch: pytest.MonkeyPatch, arguments: list[str]
) -> None:
    def reject_execution(**_kwargs: object) -> None:
        raise AssertionError("Scope validation must precede input consumption")

    monkeypatch.setattr(coverage, "run_coverage", reject_execution)
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), *arguments])
    with pytest.raises(SystemExit) as error:
        coverage.main()
    assert error.value.code == 2


def test_regeneration_paths_preserve_repository_local_symlink_spelling(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    assert (
        coverage.regeneration_path(coverage.TAXA_TRAIT_MATRIX_PATH)
        == "data/derived/taxa_trait_matrix.tsv"
    )
    root = tmp_path / "checkout"
    root.mkdir()
    external = tmp_path / "external source's data"
    external.mkdir()
    (root / "data").symlink_to(external, target_is_directory=True)
    monkeypatch.setattr(coverage, "REPO_ROOT", root)
    assert coverage.regeneration_path(root / "data/matrix.tsv") == "data/matrix.tsv"
    assert coverage.regeneration_path(external / "matrix.tsv") == str(external / "matrix.tsv")
    (tmp_path / "matrix.tsv").write_text("intended symlink-parent input")
    (root / "matrix.tsv").write_text("wrong lexically collapsed input")
    parent_path = root / "data/../matrix.tsv"
    assert parent_path.read_text() == "intended symlink-parent input"
    spelling = coverage.regeneration_path(parent_path)
    assert spelling == "data/../matrix.tsv"
    assert (root / spelling).read_text() == "intended symlink-parent input"


def test_api_rejects_missing_or_invalid_scope_before_writing(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def reject_read(_path: Path) -> str:
        raise AssertionError("Scope validation must precede input hashing")

    monkeypatch.setattr(coverage, "file_sha256", reject_read)
    output = tmp_path / "invalid-scope"
    with pytest.raises(TypeError, match="strain_rollup"):
        coverage.run_coverage(output_dir=output, tables_only=True, **inputs)
    with pytest.raises(ValueError, match="strain_rollup"):
        coverage.run_coverage(
            output_dir=output, tables_only=True, strain_rollup="unknown", **inputs
        )
    assert not output.exists()


@pytest.fixture
def inputs(tmp_path: Path) -> dict[str, Path]:
    """Three eligible profiles, plus a mouse and an ITS profile with unique taxa."""
    frames = {
        "experiments": pd.DataFrame(
            [
                ("s1", "e1", "Homo sapiens", "16S"),
                ("s2", "e2", "Homo sapiens", "WMS"),
                ("s3", "e3", "Mus musculus", "16S"),
                ("s4", "e4", "Homo sapiens", "ITS"),
            ],
            columns=["study_id", "experiment_id", "host_species", "sequencing_type"],
        ),
        "profiles": pd.DataFrame(
            [
                ("s1", "e1", "inc", 2, 1),
                ("s1", "e1", "dec", 1, 0),
                ("s2", "e2", "inc", 1, 1),
                ("s3", "e3", "inc", 1, 7),
                ("s4", "e4", "inc", 1, 9),
            ],
            columns=["study_id", "experiment_id", "inc_dec", "taxa_count", "METPO:1"],
        ),
        "matrix": pd.DataFrame(
            [
                ("t1", "Mapped species", "species", 1),
                ("t2", "No annotation", "genus", 0),
                ("t3", "Mouse only", "species", 7),
                ("t4", "ITS only", "species", 9),
                ("t5", "Not reported", "species", 0),
            ],
            columns=["taxon_id", "taxon_name", "taxonomic_rank", "METPO:1"],
        ),
        "taxa": pd.DataFrame(
            [
                ("s1", "e1", "sig1", "increased", "t1", "Mapped species", "species"),
                ("s1", "e1", "sig1", "increased", "t2", "No annotation", "genus"),
                ("s1", "e1", "sig2", "increased", "t1", "Mapped species", "species"),
                ("s1", "e1", "sig3", "decreased", "t2", "No annotation", "genus"),
                ("s2", "e2", "sig4", "increased", "t1", "Mapped species", "species"),
                ("s3", "e3", "sig5", "increased", "t3", "Mouse only", "species"),
                ("s4", "e4", "sig6", "increased", "t4", "ITS only", "species"),
            ],
            columns=[
                "study_id",
                "experiment_id",
                "signature_id",
                "direction",
                "taxon_id",
                "taxon_name",
                "taxonomic_rank",
            ],
        ),
    }
    paths = {}
    for name, frame in frames.items():
        paths[name] = tmp_path / f"{name}.tsv"
        frame.to_csv(paths[name], sep="\t", index=False)
    return paths


def test_population_and_table_contract(inputs: dict[str, Path]) -> None:
    taxa, ranks, profile_ranks, profiles, _, zero = coverage.load_and_process_data(**inputs)
    assert set(taxa.taxon_id) == {"t1", "t2"}
    assert ranks.total_taxa.sum() == 2
    assert ranks.total_occurrences.sum() == 5
    assert len(profiles) == 3
    assert profiles.mapped_trait_count.tolist() == [1, 0, 1]
    assert profiles.total_unique_taxa.tolist() == [2, 1, 1]
    assert zero.iloc[0].profile_count == 2
    assert "pct_of_mapped_profiles" in profile_ranks
    assert "pct_of_7102_mapped_profiles" not in profile_ranks
    # Existing deck consumers' two table contracts are retained.
    assert {"taxonomic_rank", "total_taxa", "zero_trait_taxa", "trait_bearing_pct"} <= set(ranks)
    assert {"taxon_id", "taxon_name", "taxonomic_rank", "trait_count"} <= set(taxa)


def test_core_counts_exclude_qualified_mediadive_columns(
    inputs: dict[str, Path], tmp_path: Path
) -> None:
    for name in ("matrix", "profiles"):
        frame = pd.read_csv(inputs[name], sep="\t")
        media = [0, 1, 0, 0, 0] if name == "matrix" else [1, 1, 0, 0, 0]
        for term in ("mediadive.medium:1", "MediaDive.medium:2"):
            frame[f"METPO:2000517|{term}"] = media
        frame.to_csv(inputs[name], sep="\t", index=False)
    output = tmp_path / "core-coverage"
    coverage.run_coverage(output_dir=output, tables_only=True, strain_rollup="off", **inputs)
    taxa = pd.read_csv(output / "taxa_kgmicrobe_coverage_details.tsv", sep="\t").set_index(
        "taxon_id"
    )
    assert taxa.loc["t2", "trait_count"] == 2
    assert taxa.loc["t2", "core_trait_count"] == 0
    assert taxa.loc["t1", "core_trait_count"] == 1
    profiles = pd.read_csv(output / "profile_trait_mapping_summary.tsv", sep="\t")
    assert profiles.mapped_trait_count.tolist() == [3, 2, 1]
    assert profiles.mapped_core_trait_count.tolist() == [1, 0, 1]
    text = (output / "README.md").read_text()
    assert "| Taxa with MediaDive annotations only | 1 |" in text
    assert "| Profiles with MediaDive annotations only | 1 |" in text


def test_vintage_is_read_and_hashed_from_the_explicit_raw_dump(
    inputs: dict[str, Path], tmp_path: Path
) -> None:
    raw = tmp_path / "full_dump.csv"
    raw.write_text("# BugSigDB fixture-912, License: CC0\nStudy\nfixture\n")
    output = tmp_path / "vintage"
    coverage.run_coverage(
        output_dir=output, tables_only=True, raw_dump=raw, strain_rollup="off", **inputs
    )
    record = json.loads((output / "provenance.json").read_text())
    assert record["corpus_identity"]["bugsigdb_vintage"] == "fixture-912"
    assert record["inputs"]["raw_dump"]["sha256"] == coverage.file_sha256(raw)
    assert record["bugsigdb_vintage_source"] == "raw_dump_header"
    command = readme_regeneration_argv((output / "README.md").read_text())
    assert Path(command[command.index("--raw-dump") + 1]).resolve() == raw.resolve()
    assert "--bugsigdb-vintage" not in command


def test_empty_ranks_and_profile_status_colors_match_across_figures(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    matrix = pd.read_csv(inputs["matrix"], sep="\t")
    matrix.loc[matrix.taxon_id.eq("t2"), "taxonomic_rank"] = "species"
    matrix.to_csv(inputs["matrix"], sep="\t", index=False)
    taxa, ranks, profile_ranks, profiles, _, _ = coverage.load_and_process_data(**inputs)
    captured = []
    monkeypatch.setattr(coverage, "save_figure", lambda fig, _path: captured.append(fig))
    coverage.plot_rank_stratified_histogram(ranks, tmp_path / "rank.png", profiles)
    coverage.plot_occurrences_vs_unique(ranks, taxa, profiles, tmp_path / "summary.png")
    coverage.plot_profile_distribution_by_rank(profile_ranks, profiles, tmp_path / "profile.png")
    try:
        assert [tick.get_text() for tick in captured[0].axes[0].get_xticklabels()] == ["Species"]
        assert [tick.get_text() for tick in captured[2].axes[0].get_xticklabels()] == ["Species"]
        mapped = coverage.matplotlib.colors.to_rgba(coverage.COLOR_ALL_MAPPED_PROFILE)
        assert captured[0].axes[2].containers[2][0].get_facecolor() == mapped
        assert captured[1].axes[0].containers[2][2].get_facecolor() == mapped
        assert captured[2].axes[0].containers[2][0].get_facecolor() == mapped
        mixed = next(
            text for text in captured[1].axes[0].texts if text.get_text().startswith("Mixed:")
        )
        assert coverage.matplotlib.colors.to_rgba(mixed.get_color()) == (1, 1, 1, 1)
        middle = captured[2].axes[0].containers[1][0]
        label = next(
            text
            for text in captured[2].axes[0].texts
            if text.get_position() == (0, middle.get_y() + middle.get_height() / 2)
        )
        assert coverage.matplotlib.colors.to_rgba(label.get_color()) == (1, 1, 1, 1)
    finally:
        for fig in captured:
            coverage.plt.close(fig)


def test_sparse_figure_annotations_stay_visible_without_covering_bars(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    matrix = pd.read_csv(inputs["matrix"], sep="\t")
    matrix.loc[matrix.taxon_id.eq("t2"), "taxonomic_rank"] = "species"
    matrix.to_csv(inputs["matrix"], sep="\t", index=False)
    taxa, ranks, _, profiles, _, _ = coverage.load_and_process_data(**inputs)
    captured = []
    monkeypatch.setattr(coverage, "save_figure", lambda fig, _path: captured.append(fig))
    coverage.plot_rank_stratified_histogram(ranks, tmp_path / "rank.png", profiles)
    coverage.plot_occurrences_vs_unique(ranks, taxa, profiles, tmp_path / "summary.png")
    try:
        for fig in captured:
            fig.canvas.draw()
        proportion_axis = captured[0].axes[1]
        parity_label = next(t for t in proportion_axis.texts if "parity threshold" in t.get_text())
        assert parity_label.get_window_extent().x1 <= proportion_axis.get_window_extent().x1
        for ax in captured[0].axes:
            bars = [patch for group in ax.containers for patch in group if patch.get_height() > 0]
            for text in ax.texts:
                box = text.get_bbox_patch()
                if box is not None:
                    assert not any(
                        box.get_window_extent().overlaps(bar.get_window_extent()) for bar in bars
                    )
        ax = captured[1].axes[0]
        bounds = ax.get_window_extent()
        for text in ax.texts:
            if coverage.matplotlib.colors.to_rgba(text.get_color()) == (1, 1, 1, 1):
                extent = text.get_window_extent()
                assert bounds.x0 <= extent.x0 < extent.x1 <= bounds.x1
    finally:
        for fig in captured:
            coverage.plt.close(fig)


def narrow_segment_frames(zero_count: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Synthetic summaries with narrow taxon segments and two tiny profile segments."""
    taxa = pd.DataFrame(
        {"coverage_category": ["zero_traits"] * zero_count + ["rich_gte_3"] * (1000 - zero_count)}
    )
    ranks = pd.DataFrame(
        {
            "total_occurrences": [1000],
            "zero_trait_occurrences": [zero_count],
            "sparse_occurrences": [0],
            "rich_occurrences": [1000 - zero_count],
        }
    )
    profiles = pd.DataFrame(
        {"mapped_trait_count": [0] * 998 + [1, 1], "zero_trait_taxa": [1] * 999 + [0]}
    )
    return ranks, taxa, profiles


def assert_fitted_segment_label(axis: Axes, expected: str, group: int, row: int) -> None:
    """Require complete text on its own row and inside its actual segment (#970)."""
    labels = [
        text
        for text in axis.texts
        if text.get_text() == expected and text.get_position()[1] == pytest.approx(row)
    ]
    assert len(labels) == 1, f"Missing fitted segment label: {expected!r} on row {row}"
    box = labels[0].get_window_extent()
    segment = axis.containers[group][row].get_window_extent()
    assert segment.contains(box.x0, box.y0), expected
    assert segment.contains(box.x1, box.y1), expected


@pytest.mark.parametrize("zero_count", [78, 2, 998])
def test_narrow_segments_keep_labels_inside_their_colored_regions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, zero_count: int
) -> None:
    """Exercise the reviewed 7.8% width, tiny zeros and a tiny mapped segment (#947)."""
    ranks, taxa, profiles = narrow_segment_frames(zero_count)
    captured = []
    monkeypatch.setattr(coverage, "save_figure", lambda fig, _path: captured.append(fig))
    coverage.plot_occurrences_vs_unique(ranks, taxa, profiles, tmp_path / "narrow.png")
    figure = captured[0]
    try:
        figure.canvas.draw()
        axis = figure.axes[0]
        bounds = axis.get_window_extent()
        inside = {
            78: [("7.8%\n(78)", 0), ("Trait-bearing: 92.2%\n(922)", 2)],
            2: [("Trait-bearing: 99.8%\n(998)", 2)],
            998: [("0 traits: 99.8%\n(998)", 0)],
        }
        for expected, group in inside[zero_count]:
            for row in (0, 1):
                assert_fitted_segment_label(axis, expected, group, row)
        for label in axis.texts:
            if not label.get_text():
                continue
            box = label.get_window_extent()
            assert bounds.x0 <= box.x0 < box.x1 <= bounds.x1, label.get_text()
            assert bounds.y0 <= box.y0 < box.y1 <= bounds.y1, label.get_text()
            if coverage.matplotlib.colors.to_rgba(label.get_color()) == (1, 1, 1, 1):
                # These fixtures have no sparse/rich boundary, so one segment must
                # actually contain all white text; axis inclusion alone is too weak.
                assert any(
                    bar.get_window_extent().contains(box.x0, box.y0)
                    and bar.get_window_extent().contains(box.x1, box.y1)
                    for group in axis.containers
                    for bar in group
                    if bar.get_width() > 0
                ), label.get_text()
        expected_callouts = {
            "Mixed: 0.1% (1)": [axis.containers[1][2]],
            "All mapped: 0.1% (1)": [axis.containers[2][2]],
        }
        if zero_count == 2:
            expected_callouts["0 traits: 0.2% (2)"] = [axis.containers[0][i] for i in (0, 1)]
        elif zero_count == 998:
            expected_callouts["Trait-bearing: 0.2% (2)"] = [axis.containers[2][i] for i in (0, 1)]
        for description, segments in expected_callouts.items():
            labels = sorted(
                (text for text in axis.texts if text.get_text() == description),
                key=lambda text: text.get_position()[1],
            )
            assert len(labels) == len(segments), f"Missing segment callout: {description}"
            for label, segment in zip(labels, segments, strict=True):
                label_box = label.get_window_extent()
                bar_box = segment.get_window_extent()
                assert label_box.y0 > bar_box.y1, f"Callout is on its bar: {description}"
                assert all(
                    not label_box.overlaps(bar.get_window_extent())
                    for group in axis.containers
                    for bar in group
                    if bar.get_width() > 0
                ), f"Callout overlaps another bar: {description}"
                endpoint = (
                    segment.get_x() + segment.get_width() / 2,
                    segment.get_y() + segment.get_height(),
                )
                leaders = [
                    text
                    for text in axis.texts
                    if isinstance(text, Annotation)
                    and text.arrow_patch is not None
                    and text.xy == pytest.approx(endpoint)
                ]
                assert len(leaders) == 1, f"Missing segment leader: {description}"
                leader = leaders[0]
                start = axis.transData.transform(leader.get_position())
                assert start[0] == pytest.approx((label_box.x0 + label_box.x1) / 2, abs=1)
                assert start[1] == pytest.approx(label_box.y0, abs=1)
                end = axis.transData.transform(leader.xy)
                assert end[0] == pytest.approx((bar_box.x0 + bar_box.x1) / 2, abs=1)
                assert end[1] == pytest.approx(bar_box.y1, abs=1)
                arrow = leader.arrow_patch
                assert arrow is not None
                vertices = np.asarray(
                    arrow.get_path().transformed(arrow.get_transform()).vertices, dtype=float
                )
                assert np.isfinite(vertices).all()
                assert np.linalg.norm(vertices[-1] - vertices[0]) > 0
                assert vertices[0, 1] > vertices[-1, 1]
        boxes = [text.get_window_extent() for text in axis.texts if text.get_text()]
        assert all(
            not left.overlaps(right) for i, left in enumerate(boxes) for right in boxes[i + 1 :]
        )
    finally:
        coverage.plt.close(figure)


def test_fitted_profile_labels_keep_their_percentage_count_and_segment(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    ranks, taxa, _ = narrow_segment_frames(78)
    profiles = pd.DataFrame(
        {"mapped_trait_count": [1] * 1000, "zero_trait_taxa": [1] * 500 + [0] * 500}
    )
    captured = []
    monkeypatch.setattr(coverage, "save_figure", lambda fig, _path: captured.append(fig))
    coverage.plot_occurrences_vs_unique(ranks, taxa, profiles, tmp_path / "inside-profiles.png")
    figure = captured[0]
    try:
        figure.canvas.draw()
        axis = figure.axes[0]
        assert_fitted_segment_label(axis, "Mixed:\n50.0%\n(500)", 1, 2)
        assert_fitted_segment_label(axis, "All mapped:\n50.0%\n(500)", 2, 2)
    finally:
        coverage.plt.close(figure)


def dense_rank_frame(original: pd.DataFrame) -> pd.DataFrame:
    """Seven fully mapped ranks expose parity boxes over the 100% label at y=50."""
    ranks = pd.concat([original.iloc[:1]] * 7, ignore_index=True)
    ranks["taxonomic_rank"] = coverage.RANK_ORDER[:7]
    ranks["zero_trait_pct"] = 0
    ranks["sparse_1_2_pct"] = 0
    ranks["rich_gte_3_pct"] = 100
    ranks["total_taxa"] = 100
    ranks["zero_trait_taxa"] = 0
    ranks["sparse_1_2_taxa"] = 0
    ranks["rich_gte_3_taxa"] = 100
    ranks["trait_bearing_pct"] = 100
    return ranks


def four_digit_rank_frame(original: pd.DataFrame) -> pd.DataFrame:
    """Synthetic seven-rank frame with the reviewed four-digit label widths (#959)."""
    ranks = dense_rank_frame(original).astype(
        dict.fromkeys(
            ("zero_trait_pct", "sparse_1_2_pct", "rich_gte_3_pct", "trait_bearing_pct"), float
        )
    )
    species = ranks.taxonomic_rank.eq("species")
    ranks.loc[species, ["total_taxa", "zero_trait_taxa", "sparse_1_2_taxa", "rich_gte_3_taxa"]] = [
        3165,
        1636,
        204,
        1325,
    ]
    ranks.loc[species, ["zero_trait_pct", "sparse_1_2_pct", "rich_gte_3_pct"]] = (
        np.array([1636, 204, 1325]) / 3165 * 100
    )
    ranks.loc[species, "trait_bearing_pct"] = (204 + 1325) / 3165 * 100
    return ranks


@pytest.mark.parametrize("dpi", [100, 300])
def test_four_digit_species_labels_fit_their_actual_bars(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, dpi: int
) -> None:
    _, original, _, profiles, _, _ = coverage.load_and_process_data(**inputs)
    captured = []
    monkeypatch.setattr(coverage, "save_figure", lambda fig, _path: captured.append(fig))
    coverage.plot_rank_stratified_histogram(
        four_digit_rank_frame(original), tmp_path / "four-digits.png", profiles
    )
    figure = captured[0]
    try:
        figure.set_dpi(dpi)
        figure.canvas.draw()
        axis = figure.axes[1]
        for value, percentage, group in (("1,636", "51.7%", 0), ("1,325", "41.9%", 2)):
            labels = [text for text in axis.texts if value in text.get_text()]
            assert len(labels) == 1
            label = labels[0]
            assert percentage in label.get_text()
            box = label.get_window_extent()
            species_index = [tick.get_text() for tick in axis.get_xticklabels()].index("Species")
            bar = axis.containers[group][species_index].get_window_extent()
            assert bar.contains(box.x0, box.y0)
            assert bar.contains(box.x1, box.y1)
            assert label.get_fontsize() >= 6
    finally:
        coverage.plt.close(figure)


def test_parity_annotation_cannot_cover_dense_rank_bars_or_their_labels(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, original, _, profiles, _, _ = coverage.load_and_process_data(**inputs)
    ranks = dense_rank_frame(original)
    captured = []
    monkeypatch.setattr(coverage, "save_figure", lambda fig, _path: captured.append(fig))
    coverage.plot_rank_stratified_histogram(ranks, tmp_path / "dense.png", profiles)
    figure = captured[0]
    try:
        figure.canvas.draw()
        axis = figure.axes[1]
        parity = next(text for text in axis.texts if "parity threshold" in text.get_text())
        patch = parity.get_bbox_patch()
        box = patch.get_window_extent() if patch is not None else parity.get_window_extent()
        bars = [bar for group in axis.containers for bar in group if bar.get_height() > 0]
        assert not any(box.overlaps(bar.get_window_extent()) for bar in bars)
        labels = [text for text in axis.texts if text is not parity]
        assert len(labels) == 7
        assert not any(box.overlaps(text.get_window_extent()) for text in labels)
        for label in labels:
            bounds = label.get_window_extent()
            assert any(
                bar.get_window_extent().contains(bounds.x0, bounds.y0)
                and bar.get_window_extent().contains(bounds.x1, bounds.y1)
                for bar in bars
            ), label.get_text()
    finally:
        coverage.plt.close(figure)


def test_readme_reconciles_with_written_tables_and_provenance(
    inputs: dict[str, Path],
    tmp_path: Path,
) -> None:
    output = tmp_path / "coverage"
    coverage.run_coverage(output_dir=output, tables_only=True, strain_rollup="off", **inputs)
    profiles = pd.read_csv(output / "profile_trait_mapping_summary.tsv", sep="\t")
    taxa = pd.read_csv(output / "taxa_kgmicrobe_coverage_details.tsv", sep="\t")
    readme = (output / "README.md").read_text()
    assert f"| Total profiles | {len(profiles)} |" in readme
    assert f"| Unique taxa | {len(taxa)} |" in readme
    assert "| Profiles with mapped traits | 2 |" in readme
    assert "| Zero-trait profiles | 1 |" in readme
    assert "| Taxon signature mentions | 5 |" in readme
    assert "Homo sapiens" in readme
    assert "16S" in readme
    assert "WMS" in readme
    provenance = json.loads((output / "provenance.json").read_text())
    assert provenance["population"]["contrast_filter"] == "none"
    assert provenance["counts"]["excluded_profiles"] == 2
    assert provenance["counts"]["selected_profiles"] == 3
    assert len(provenance["inputs"]["profiles"]["sha256"]) == 64
    assert provenance["outputs"]["README.md"]["sha256"] == coverage.file_sha256(
        output / "README.md"
    )


@pytest.mark.parametrize(
    "problem",
    [
        "missing_metadata",
        "missing_profile",
        "missing_taxon",
        "duplicate_profile",
        "empty_population",
        "wrong_counts",
    ],
)
def test_mismatched_inputs_fail_before_writing(
    inputs: dict[str, Path],
    tmp_path: Path,
    problem: str,
) -> None:
    name = {"missing_metadata": "experiments", "missing_taxon": "matrix"}.get(problem, "profiles")
    frame = pd.read_csv(inputs[name], sep="\t")
    if problem in {"missing_metadata", "missing_profile", "missing_taxon"}:
        frame = frame.iloc[1:]
    elif problem == "duplicate_profile":
        frame = pd.concat([frame, frame.iloc[:1]])
    elif problem == "empty_population":
        frame = frame.iloc[3:]
    else:
        frame.loc[0, "METPO:1"] = 99
    frame.to_csv(inputs[name], sep="\t", index=False)
    output = tmp_path / "invalid"
    with pytest.raises(ValueError, match=r"metadata|missing|Duplicate|No profiles|reconcile"):
        coverage.run_coverage(output_dir=output, tables_only=True, strain_rollup="off", **inputs)
    assert not output.exists()


def test_zero_mapped_profiles_still_produce_tables(inputs: dict[str, Path], tmp_path: Path) -> None:
    for name in ("matrix", "profiles"):
        frame = pd.read_csv(inputs[name], sep="\t")
        frame["METPO:1"] = 0
        frame.to_csv(inputs[name], sep="\t", index=False)
    output = tmp_path / "zeros"
    coverage.run_coverage(output_dir=output, tables_only=True, strain_rollup="off", **inputs)
    ranks = pd.read_csv(output / "profile_coverage_by_rank.tsv", sep="\t")
    assert ranks.pct_of_mapped_profiles.eq(0).all()
    assert "| Profiles with mapped traits | 0 |" in (output / "README.md").read_text()


def test_cli_accepts_canonical_profiles_without_taxa_count(
    inputs: dict[str, Path],
    tmp_path: Path,
) -> None:
    profiles = pd.read_csv(inputs["profiles"], sep="\t").drop(columns="taxa_count")
    profiles.to_csv(inputs["profiles"], sep="\t", index=False)
    output = tmp_path / "cli"
    args = [
        sys.executable,
        str(SCRIPT),
        "--tables-only",
        "--output-dir",
        str(output),
        "--strain-rollup",
        "off",
    ]
    for name, path in inputs.items():
        args.extend([f"--{name}", str(path)])
    result = subprocess.run(args, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert "| Total profiles | 3 |" in (output / "README.md").read_text()


@pytest.mark.parametrize("all_zero", [False, True])
def test_all_figure_paths_handle_sparse_ranks_and_empty_mapped_sets(
    inputs: dict[str, Path],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    all_zero: bool,
) -> None:
    if all_zero:
        for name in ("profiles", "matrix"):
            frame = pd.read_csv(inputs[name], sep="\t")
            frame["METPO:1"] = 0
            frame.to_csv(inputs[name], sep="\t", index=False)
    taxa, ranks, profile_ranks, profiles, _, _ = coverage.load_and_process_data(**inputs)
    written = []

    def inspect_figure(fig: object, path: Path) -> None:
        written.append(path.stem)
        coverage.plt.close(fig)

    monkeypatch.setattr(coverage, "save_figure", inspect_figure)
    coverage.render_figures(taxa, ranks, profile_ranks, profiles, tmp_path)
    assert set(written) == set(coverage.FIGURE_STEMS)


def test_unrecognized_rank_is_counted_in_other(inputs: dict[str, Path]) -> None:
    frame = pd.read_csv(inputs["matrix"], sep="\t")
    frame.loc[0, "taxonomic_rank"] = "subspecies"
    frame.to_csv(inputs["matrix"], sep="\t", index=False)
    _, ranks, _, profiles, _, _ = coverage.load_and_process_data(**inputs)
    assert ranks.loc[ranks.taxonomic_rank == "other", "total_taxa"].item() == 1
    assert profiles.mapped_taxa_other.sum() == 2


def test_tables_only_refuses_to_leave_old_figures(
    inputs: dict[str, Path],
    tmp_path: Path,
) -> None:
    output = tmp_path / "old"
    output.mkdir()
    (output / "old.svg").write_text("old figure")
    with pytest.raises(ValueError, match="existing figures"):
        coverage.run_coverage(output_dir=output, tables_only=True, strain_rollup="off", **inputs)
    assert not (output / "provenance.json").exists()


def test_dynamic_figure_labels_and_median_positions(
    inputs: dict[str, Path],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _, ranks, _, profiles, _, _ = coverage.load_and_process_data(**inputs)
    captured = []
    monkeypatch.setattr(coverage, "save_figure", lambda fig, _path: captured.append(fig))
    coverage.plot_rank_stratified_histogram(ranks, tmp_path / "rank.png", profiles)
    coverage.plot_histogram1_mapped_profiles(profiles, tmp_path / "mapped.png")
    text = "\n".join(t.get_text() for fig in captured for ax in fig.axes for t in ax.texts)
    assert captured[0].axes[0].get_title(loc="left") == "A. Unique Taxa by Rank (N=2)"
    assert captured[0].axes[2].get_title(loc="left") == "C. Profiles by Rank and Mapping (N=3)"
    assert "7,329" not in text
    assert "7,102" not in text
    assert captured[1].axes[0].lines[0].get_xdata()[0] == 0
    assert captured[1].axes[1].lines[0].get_xdata()[0] == 0
    for fig in captured:
        coverage.plt.close(fig)


def test_svg_export_strips_date_metadata(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fig = coverage.plt.figure()
    calls = []
    monkeypatch.setattr(fig, "savefig", lambda path, **kwargs: calls.append((path, kwargs)))
    coverage.save_figure(fig, tmp_path / "figure.png")
    svg_calls = [kwargs for path, kwargs in calls if path.suffix == ".svg"]
    assert svg_calls == [{"bbox_inches": "tight", "facecolor": "white", "metadata": {"Date": None}}]


def test_histogram2_family_and_mixed_legend_colors_are_distinct(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    taxa, _, _, profiles, _, _ = coverage.load_and_process_data(**inputs)
    captured = []
    monkeypatch.setattr(coverage, "save_figure", lambda fig, _path: captured.append(fig))
    coverage.plot_histogram2_taxa_traits_and_zero_frequency(taxa, profiles, tmp_path / "zero.png")
    legend = captured[0].legends[0]
    colors = {
        label.get_text(): tuple(handle.get_facecolor())
        for label, handle in zip(legend.get_texts(), legend.legend_handles, strict=True)
    }
    assert colors["Family"] != colors["Mixed Profile (≥1 Trait)"]
    assert colors["All-mapped Profile"] != colors["Species"]
    assert coverage.RANK_COLORS["family"] != coverage.COLOR_SPARSE
    # Rank and status appear in the same legend: avoid two near-identical slates.
    distance = (
        sum(
            (first - second) ** 2
            for first, second in zip(
                colors["Superkingdom"][:3], colors["Mixed Profile (≥1 Trait)"][:3], strict=True
            )
        )
        ** 0.5
    )
    assert distance > 0.3
    coverage.plt.close(captured[0])


def test_export_omits_empty_ranks_and_records_corpus_identity(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "identity"
    sidecar = tmp_path / "kg.json"
    sidecar.write_text(json.dumps({"release": "fixture-release", "files": {}}))
    monkeypatch.setattr(coverage, "repository_state", lambda _root: ("a" * 40, True))
    coverage.run_coverage(
        output_dir=output,
        strain_rollup="off",
        tables_only=True,
        kg_provenance=sidecar,
        bugsigdb_vintage="fixture-vintage",
        **inputs,
    )
    ranks = pd.read_csv(output / "taxa_coverage_by_rank.tsv", sep="\t")
    assert ranks.total_taxa.gt(0).all()
    profile_ranks = pd.read_csv(output / "profile_coverage_by_rank.tsv", sep="\t")
    assert profile_ranks.profiles_containing_rank.gt(0).all()
    record = json.loads((output / "provenance.json").read_text())
    assert record["producer"]["code_commit"].endswith("+dirty")
    assert record["corpus_identity"]["bugsigdb_vintage"] == "fixture-vintage"
    assert record["corpus_identity"]["kg_microbe"]["release_tag"] == "fixture-release"
    assert record["corpus_identity"]["matrix"]["sha256"] == coverage.file_sha256(inputs["profiles"])
    assert record["inputs"]["kg_provenance"]["sha256"] == coverage.file_sha256(sidecar)
    assert record["corpus_identity"]["matrix"]["path"] == "external/profiles.tsv"
    assert str(tmp_path) not in json.dumps(record)


def test_plate6_reads_compact_remediation_and_nonempty_rank_rows(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Check the unchanged consumer's rendered artists, not just TSV column names."""
    for name in ("taxa", "matrix"):
        frame = pd.read_csv(inputs[name], sep="\t")
        frame.loc[frame.taxon_id.eq("t2"), "taxon_id"] = "NCBITaxon:541000"
        frame.to_csv(inputs[name], sep="\t", index=False)
    output = tmp_path / "trait-visualization/trait-profile-investigation/taxa-coverage"
    coverage.run_coverage(output_dir=output, tables_only=True, strain_rollup="off", **inputs)
    spec = importlib.util.spec_from_file_location(
        "coverage_plate_contract", SCRIPT.with_name("generate_tufte_plates.py")
    )
    assert spec is not None
    assert spec.loader is not None
    plate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(plate)
    captured = []
    monkeypatch.setattr(
        plate.plt.Figure, "savefig", lambda fig, *_args, **_kwargs: captured.append(fig)
    )
    plate.build_plate6(tmp_path, tmp_path)
    fig = captured[0]
    fig.canvas.draw()
    text = "\n".join(t.get_text() for ax in fig.axes for t in ax.texts)
    assert "0/0" not in text
    assert "Other / unrecorded (0 taxa)" not in text
    alias = [t for t in fig.axes[2].texts if "NCBITaxon:216572" in t.get_text()]
    assert len(alias) == 1
    assert "Review possible" not in alias[0].get_text()
    assert alias[0].get_window_extent().x1 <= fig.axes[2].get_window_extent().x1
    # Downward propagation is deliberately not restored; #864 owns that consumer label.
    _, remediation = coverage.diagnose_zero_trait_taxon("t6", "uncultured bacterium", "species")
    assert "do not establish leaf-taxon traits" in remediation
    coverage.plt.close(fig)

    # Also verify that the three new #957 / #1009 diagnosis categories render with distinct
    # Panel B short names, distinct bar colors, and Panel C remediation labels (#1009):
    zero_path = output / "top_zero_trait_bugsigdb_taxa.tsv"
    zero_df = pd.read_csv(zero_path, sep="\t")
    extra_rows = pd.DataFrame(
        [
            {
                **zero_df.iloc[0].to_dict(),
                "zero_trait_rank": len(zero_df) + 1,
                "taxon_id": "NCBITaxon:155900",
                "taxon_name": "uncultured organism",
                "taxonomic_rank": "species",
                "signature_occurrences": 17,
                "profile_count": 17,
                "study_count": 8,
                "diagnosis_category": "Unresolved Taxonomic Label",
                "suggested_remediation": (
                    "Review the reported taxonomic label and source assignment before trait "
                    "matching; an unresolved label does not establish a multi-taxon composite"
                ),
            },
            {
                **zero_df.iloc[0].to_dict(),
                "zero_trait_rank": len(zero_df) + 2,
                "taxon_id": "NCBITaxon:1505663",
                "taxon_name": "Thomasclavelia",
                "taxonomic_rank": "genus",
                "signature_occurrences": 14,
                "profile_count": 14,
                "study_count": 11,
                "diagnosis_category": "Taxonomy Unresolved",
                "suggested_remediation": (
                    "Review placement in the supplied taxonomy; missing lineage does not "
                    "establish rarity, novelty or a bacterial isolate"
                ),
            },
            {
                **zero_df.iloc[0].to_dict(),
                "zero_trait_rank": len(zero_df) + 3,
                "taxon_id": "NCBITaxon:10239",
                "taxon_name": "Caudoviricetes sp.",
                "taxonomic_rank": "species",
                "signature_occurrences": 11,
                "profile_count": 11,
                "study_count": 5,
                "diagnosis_category": "Non-Bacterial / Viral Taxon",
                "suggested_remediation": (
                    "Review viral annotation resources; missing retained annotations do not "
                    "establish biological absence"
                ),
            },
        ]
    )
    pd.concat([zero_df, extra_rows], ignore_index=True).to_csv(zero_path, sep="\t", index=False)
    captured.clear()
    plate.build_plate6(tmp_path, tmp_path)
    fig_new = captured[0]
    fig_new.canvas.draw()
    panel_b_yticks = [t.get_text() for t in fig_new.axes[1].get_yticklabels()]
    assert "Unresolved Taxon Label" in panel_b_yticks
    assert "Taxonomy Unresolved" in panel_b_yticks
    assert "Viral Taxon (Non-Bacterial)" in panel_b_yticks
    panel_c_texts = "\n".join(t.get_text() for t in fig_new.axes[2].texts)
    assert "Review taxonomic label / source" in panel_c_texts
    assert "Review taxonomy placement" in panel_c_texts
    assert "Review viral annotation resources" in panel_c_texts
    coverage.plt.close(fig_new)


@pytest.mark.parametrize("all_zero", [False, True])
def test_every_readme_table_cell_matches_the_exported_tables(
    inputs: dict[str, Path], tmp_path: Path, all_zero: bool
) -> None:
    if all_zero:
        for name in ("profiles", "matrix"):
            frame = pd.read_csv(inputs[name], sep="\t")
            frame["METPO:1"] = 0
            frame.to_csv(inputs[name], sep="\t", index=False)
    output = tmp_path / "readme"
    coverage.run_coverage(output_dir=output, tables_only=True, strain_rollup="off", **inputs)
    profiles = pd.read_csv(output / "profile_trait_mapping_summary.tsv", sep="\t")
    taxa = pd.read_csv(output / "taxa_kgmicrobe_coverage_details.tsv", sep="\t")
    ranks = pd.read_csv(output / "taxa_coverage_by_rank.tsv", sep="\t")
    mapped = profiles.loc[profiles.has_mapped_traits]
    zero = taxa.loc[taxa.trait_count.eq(0)]
    expected = {
        "Total profiles": [str(len(profiles))],
        "Profiles with mapped traits": [str(len(mapped))],
        "Zero-trait profiles": [str(len(profiles) - len(mapped))],
        "Unique taxa": [str(len(taxa))],
        "Trait-bearing taxa": [str(taxa.has_traits.sum())],
        "Taxa with non-MediaDive traits": [str(taxa.core_trait_count.gt(0).sum())],
        "Taxa with MediaDive annotations only": [
            str((taxa.has_traits & taxa.core_trait_count.eq(0)).sum())
        ],
        "Profiles with non-MediaDive traits": [str(profiles.mapped_core_trait_count.gt(0).sum())],
        "Profiles with MediaDive annotations only": [
            str((profiles.has_mapped_traits & profiles.mapped_core_trait_count.eq(0)).sum())
        ],
        "Zero-trait taxa": [str(len(zero))],
        "Taxon signature mentions": [str(taxa.signature_occurrences.sum())],
        "Zero-trait taxa recurring in ≥20 profiles": [str(zero.profile_count.ge(20).sum())],
    }
    for row in ranks.itertuples():
        expected[str(row.display_rank)] = [
            str(row.total_taxa),
            str(row.zero_trait_taxa),
            f"{row.trait_bearing_pct:.2f}",
            str(row.total_occurrences),
            str(row.profiles_with_rank),
        ]
    for label, values in {
        "Traits per mapped profile": mapped.mapped_trait_count,
        "Mapped taxa per mapped profile": mapped.mapped_trait_taxa,
        "Traits per reported taxon": taxa.trait_count,
        "Traits per trait-bearing taxon": taxa.loc[taxa.has_traits, "trait_count"],
    }.items():
        expected[label] = (
            [
                f"{value:g}"
                for value in [
                    values.median(),
                    values.quantile(0.25),
                    values.quantile(0.75),
                    values.mean(),
                ]
            ]
            if len(values)
            else ["not estimable"] * 4
        )
    actual = {}
    for line in (output / "README.md").read_text().splitlines():
        if not line.startswith("| "):
            continue
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        if cells[0] in {"Metric", "Rank", "Distribution", "---"}:
            continue
        actual[cells[0]] = cells[1:]
    assert actual == expected


def test_two_fixture_runs_have_identical_tables_readme_and_svgs(
    inputs: dict[str, Path], tmp_path: Path
) -> None:
    first, second = tmp_path / "first", tmp_path / "second"
    # Serialize real SVGs, including their IDs and metadata; do not mock savefig.
    with coverage.plt.rc_context({"font.family": "DejaVu Sans"}):
        coverage.run_coverage(output_dir=first, strain_rollup="off", **inputs)
        coverage.run_coverage(output_dir=second, strain_rollup="off", **inputs)
    compared = []
    for path in first.iterdir():
        if path.name == "README.md":
            normalized = []
            for directory in (first, second):
                text = (directory / path.name).read_text()
                arguments = readme_regeneration_argv(text)
                command = shlex.join(arguments)
                index = arguments.index("--output-dir") + 1
                assert Path(arguments[index]).resolve() == directory.resolve()
                arguments[index] = "OUTPUT_DIRECTORY"
                normalized.append(text.replace(command, shlex.join(arguments)))
            # Only the explicitly recorded destination differs between these runs.
            assert normalized[0] == normalized[1]
            compared.append(path.name)
        elif path.suffix in {".tsv", ".svg"}:
            assert path.read_bytes() == (second / path.name).read_bytes(), path.name
            compared.append(path.name)
    assert len(compared) == 12  # Five TSVs, README, six SVGs.


def test_pre_742_experiment_metadata_has_actionable_error(
    inputs: dict[str, Path], tmp_path: Path
) -> None:
    frame = pd.read_csv(inputs["experiments"], sep="\t").drop(columns="sequencing_type")
    frame.to_csv(inputs["experiments"], sep="\t", index=False)
    with pytest.raises(ValueError, match=r"#742.*sequencing_type"):
        coverage.run_coverage(
            output_dir=tmp_path / "old", tables_only=True, strain_rollup="off", **inputs
        )
    assert not (tmp_path / "old").exists()


def test_coverage_and_builder_use_the_same_population_predicate(
    inputs: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    from mtp import study_trait_profiles
    from mtp.human_host_corpus import experiment_in_profile_population

    assert coverage.experiment_in_profile_population is experiment_in_profile_population
    builder_globals = inspect.getclosurevars(study_trait_profiles.collect_profile_taxa).globals
    assert builder_globals["experiment_in_profile_population"] is experiment_in_profile_population
    assert experiment_in_profile_population(
        " Homo sapiens ", " WMS ", host_species="Homo sapiens", sequencing_types=("16S", "WMS")
    )
    assert not experiment_in_profile_population(
        "Homo sapiens", None, host_species="Homo sapiens", sequencing_types=("16S", "WMS")
    )
    monkeypatch.setattr(
        coverage, "experiment_in_profile_population", lambda *_args, **_kwargs: False
    )
    with pytest.raises(ValueError, match="No profiles"):
        coverage.load_and_process_data(**inputs)


def test_unrecognized_directions_and_empty_taxon_ids_are_skipped_like_builder(
    inputs: dict[str, Path],
) -> None:
    taxa = pd.read_csv(inputs["taxa"], sep="\t")
    invalid = taxa.iloc[:2].copy()
    invalid.loc[invalid.index[0], "direction"] = "no difference"
    invalid.loc[invalid.index[1], "taxon_id"] = " "
    pd.concat([taxa, invalid]).to_csv(inputs["taxa"], sep="\t", index=False)
    selected, _, _, profiles, _, _ = coverage.load_and_process_data(**inputs)
    assert selected.signature_occurrences.sum() == 5
    assert len(profiles) == 3


def test_histogram2_summary_boxes_do_not_cover_bar_labels(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    taxa, _, _, profiles, _, _ = coverage.load_and_process_data(**inputs)
    captured = []
    monkeypatch.setattr(coverage, "save_figure", lambda fig, _path: captured.append(fig))
    coverage.plot_histogram2_taxa_traits_and_zero_frequency(taxa, profiles, tmp_path / "labels.png")
    fig = captured[0]
    fig.canvas.draw()
    for axis in fig.axes[1:]:
        summary = next(text for text in axis.texts if text.get_bbox_patch() is not None)
        bounds = summary.get_bbox_patch().get_window_extent()
        assert all(
            not bounds.overlaps(text.get_window_extent())
            for text in axis.texts
            if text is not summary
        )
    coverage.plt.close(fig)
