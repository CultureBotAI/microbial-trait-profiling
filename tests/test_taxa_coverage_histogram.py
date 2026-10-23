"""Canonical population, data-derived labels and provenance regressions (#820/#821)."""

import importlib.util
import inspect
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/create_taxa_coverage_histogram.py"
SPEC = importlib.util.spec_from_file_location("taxa_coverage_histogram", SCRIPT)
assert SPEC is not None
assert SPEC.loader is not None
coverage = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(coverage)


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


def test_readme_reconciles_with_written_tables_and_provenance(
    inputs: dict[str, Path],
    tmp_path: Path,
) -> None:
    output = tmp_path / "coverage"
    coverage.run_coverage(output_dir=output, tables_only=True, **inputs)
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
        coverage.run_coverage(output_dir=output, tables_only=True, **inputs)
    assert not output.exists()


def test_zero_mapped_profiles_still_produce_tables(inputs: dict[str, Path], tmp_path: Path) -> None:
    for name in ("matrix", "profiles"):
        frame = pd.read_csv(inputs[name], sep="\t")
        frame["METPO:1"] = 0
        frame.to_csv(inputs[name], sep="\t", index=False)
    output = tmp_path / "zeros"
    coverage.run_coverage(output_dir=output, tables_only=True, **inputs)
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
    args = [sys.executable, str(SCRIPT), "--tables-only", "--output-dir", str(output)]
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
        coverage.run_coverage(output_dir=output, tables_only=True, **inputs)
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
    assert "N = 3" in text
    assert "N = 2 unique taxa" in text
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
        tables_only=True,
        kg_provenance=sidecar,
        bugsigdb_vintage="fixture-vintage",
        **inputs,
    )
    ranks = pd.read_csv(output / "taxa_coverage_by_rank.tsv", sep="\t")
    assert ranks.total_taxa.gt(0).all()
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
    coverage.run_coverage(output_dir=output, tables_only=True, **inputs)
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
    coverage.run_coverage(output_dir=output, tables_only=True, **inputs)
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
    with coverage.plt.rc_context(
        {"svg.hashsalt": "coverage-fixture", "font.family": "DejaVu Sans"}
    ):
        coverage.run_coverage(output_dir=first, **inputs)
        coverage.run_coverage(output_dir=second, **inputs)
    compared = []
    for path in first.iterdir():
        if path.suffix in {".tsv", ".svg"} or path.name == "README.md":
            assert path.read_bytes() == (second / path.name).read_bytes(), path.name
            compared.append(path.name)
    assert len(compared) == 12  # Five TSVs, README, six SVGs.


def test_pre_742_experiment_metadata_has_actionable_error(
    inputs: dict[str, Path], tmp_path: Path
) -> None:
    frame = pd.read_csv(inputs["experiments"], sep="\t").drop(columns="sequencing_type")
    frame.to_csv(inputs["experiments"], sep="\t", index=False)
    with pytest.raises(ValueError, match=r"#742.*sequencing_type"):
        coverage.run_coverage(output_dir=tmp_path / "old", tables_only=True, **inputs)
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
