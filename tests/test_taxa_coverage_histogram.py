"""Canonical population, data-derived labels and provenance regressions (#820/#821)."""

import importlib.util
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
    coverage.plt.close(captured[0])
