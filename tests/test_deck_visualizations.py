"""Tests for publication slide deck visualization suite (#347, #354)."""

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from mtp.deck_visualizations import (
    TABLE_1_TRAITS,
    generate_all_deck_figures,
    generate_figure1_condition_anatomy_ora_heatmap,
    generate_figure2_calibrated_ora_infographic,
    generate_figure3_lead_findings_ora_evidence,
    generate_figure4_annotation_density_pc1_landscape,
    generate_figure5_pd_crc_taxon_attribution_bubble,
    main,
)


@pytest.fixture
def mock_pipeline_data(tmp_path: Path) -> dict[str, Path]:
    """Create isolated, minimal fixture tables modeling the pipeline outputs."""
    derived_dir = tmp_path / "derived"
    derived_dir.mkdir(parents=True, exist_ok=True)

    # 1. Experiments
    exp_rows = [
        ("study_1", "exp_1", "Parkinson's disease", "Feces", "Homo sapiens"),
        ("study_1", "exp_2", "Parkinson's disease", "Feces", "Homo sapiens"),
        ("study_2", "exp_3", "Parkinson's disease", "Feces", "Homo sapiens"),
        ("study_3", "exp_4", "HIV infection", "Feces", "Homo sapiens"),
        ("study_4", "exp_5", "HIV infection", "Feces", "Homo sapiens"),
        ("study_5", "exp_6", "Colorectal cancer", "Feces", "Homo sapiens"),
        ("study_6", "exp_7", "Colorectal cancer", "Feces", "Homo sapiens"),
        ("study_7", "exp_8", "Crohn's disease", "Ileum", "Homo sapiens"),
    ]
    exp_df = pd.DataFrame(
        exp_rows,
        columns=["study_id", "experiment_id", "condition", "body_site", "host_species"],
    )
    exp_df["sequencing_type"] = "16S"
    exp_path = derived_dir / "experiments.tsv"
    exp_df.to_csv(exp_path, sep="\t", index=False)

    # 2. Profiles
    trait_cols = list(TABLE_1_TRAITS.keys())
    profile_rows = []
    # For each experiment, add 'inc' and 'dec' rows
    for s_id, e_id, cond, _, _ in exp_rows:
        # inc row
        inc_vals = dict.fromkeys(trait_cols, 5.0)
        if "Parkinson" in cond:
            inc_vals["METPO:1000623"] = 12.0  # Halophilic enriched in inc
            inc_vals["METPO:2000103|GO:0030245"] = 1.0  # Cellulose lower in inc
        elif "HIV" in cond:
            inc_vals["METPO:2000103|GO:0030245"] = 0.5
        elif "Colorectal" in cond:
            inc_vals["METPO:2000103|GO:0030245"] = 1.0

        profile_rows.append({"study_id": s_id, "experiment_id": e_id, "inc_dec": "inc", **inc_vals})

        # dec row
        dec_vals = dict.fromkeys(trait_cols, 5.0)
        if "Parkinson" in cond:
            dec_vals["METPO:1000623"] = 2.0
            dec_vals["METPO:2000103|GO:0030245"] = 10.0  # Cellulose higher in dec
        elif "HIV" in cond:
            dec_vals["METPO:2000103|GO:0030245"] = 8.0
        elif "Colorectal" in cond:
            dec_vals["METPO:2000103|GO:0030245"] = 9.0

        profile_rows.append({"study_id": s_id, "experiment_id": e_id, "inc_dec": "dec", **dec_vals})

    prof_df = pd.DataFrame(profile_rows)
    prof_path = derived_dir / "study_trait_profiles_ora_wide.tsv"
    prof_df.to_csv(prof_path, sep="\t", index=False)

    # 3. Taxa trait matrix & columns
    taxa_rows = []
    for i in range(50):
        t_row: dict[str, Any] = {
            "taxon_id": f"NCBI:{1000 + i}",
            "taxon_name": f"Taxon_species_{i}",
        }
        for idx, t in enumerate(trait_cols):
            # assign non-zero to create variation
            t_row[t] = 1 if (i + idx) % 3 == 0 else 0
        taxa_rows.append(t_row)

    taxa_df = pd.DataFrame(taxa_rows)
    taxa_path = derived_dir / "taxa_trait_matrix.tsv"
    taxa_df.to_csv(taxa_path, sep="\t", index=False)

    cols_rows = []
    for t_id, t_label in TABLE_1_TRAITS.items():
        cols_rows.append(
            {
                "column_id": t_id,
                "column_label": t_label,
                "predicate_id": "",
                "predicate_name": "",
                "object_id": t_id,
                "object_name": t_label,
                "object_category": "biolink:PhenotypicQuality",
            }
        )
    cols_df = pd.DataFrame(cols_rows)
    cols_path = derived_dir / "taxa_trait_matrix_columns.tsv"
    cols_df.to_csv(cols_path, sep="\t", index=False)

    # 4. BugSigDB Taxa
    bst_rows = [
        ("study_1", "exp_1", "increased", "NCBI:1000", "Taxon_species_0", "species"),
        ("study_1", "exp_1", "decreased", "NCBI:1001", "Taxon_species_1", "species"),
        ("study_2", "exp_3", "increased", "NCBI:1000", "Taxon_species_0", "species"),
        ("study_5", "exp_6", "increased", "NCBI:1002", "Taxon_species_2", "species"),
        ("study_6", "exp_7", "decreased", "NCBI:1002", "Taxon_species_2", "species"),
    ]
    bst_df = pd.DataFrame(
        bst_rows,
        columns=[
            "study_id",
            "experiment_id",
            "direction",
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
        ],
    )
    bst_path = derived_dir / "bugsigdb_taxa.tsv"
    bst_df.to_csv(bst_path, sep="\t", index=False)

    # 5. Density adjustment outputs
    density_dir = derived_dir / "density_adjustment"
    density_dir.mkdir(parents=True, exist_ok=True)
    t1_effects_rows = [
        {
            "trait_id": t_id,
            "trait_label": t_label,
            "raw_mean_diff": 10.0,
            "raw_wilcoxon_p": 0.001,
            "adj_mean_diff": 8.0,
            "adj_wilcoxon_p": 0.005,
            "relative_change_pct": -20.0,
            "status": "survived (significant post-adjustment)",
            # Six equal p-values have BH q equal to p over the six-trait family (#1321).
            "raw_wilcoxon_bh_q": 0.001,
            "adj_wilcoxon_bh_q": 0.005,
            "bh_family_size": len(TABLE_1_TRAITS),
        }
        for t_id, t_label in TABLE_1_TRAITS.items()
    ]
    t1_effects_df = pd.DataFrame(t1_effects_rows)
    t1_effects_path = density_dir / "table1_density_adjusted_effects.tsv"
    t1_effects_df.to_csv(t1_effects_path, sep="\t", index=False)

    # 6. Ranked enrichment & depletion results
    enc_rows = [
        {"trait_id": t, "trait_label": TABLE_1_TRAITS[t], "q_value": 0.01} for t in trait_cols[:4]
    ]
    enc_path = derived_dir / "trait_enrichment_results.tsv"
    with enc_path.open("w", encoding="utf-8") as stream:
        stream.write("# Stat-Type: fdr_q\n")
        pd.DataFrame(enc_rows).to_csv(stream, sep="\t", index=False)

    dep_rows = [
        {"trait_id": t, "trait_label": TABLE_1_TRAITS[t], "q_value": 0.01} for t in trait_cols[2:]
    ]
    dep_path = derived_dir / "trait_depletion_results.tsv"
    with dep_path.open("w", encoding="utf-8") as stream:
        stream.write("# Stat-Type: fdr_q\n")
        pd.DataFrame(dep_rows).to_csv(stream, sep="\t", index=False)

    return {
        "derived_dir": derived_dir,
        "experiments": exp_path,
        "profiles": prof_path,
        "taxa_matrix": taxa_path,
        "columns": cols_path,
        "bugsigdb_taxa": bst_path,
        "table1_effects": t1_effects_path,
        "enrichment_results": enc_path,
        "depletion_results": dep_path,
    }


def test_figure1_generation(mock_pipeline_data: dict[str, Path], tmp_path: Path) -> None:
    """Verify Figure 1 generates valid PNG with custom min_studies threshold."""
    out_file = tmp_path / "fig1.png"
    result = generate_figure1_condition_anatomy_ora_heatmap(
        profiles_path=mock_pipeline_data["profiles"],
        experiments_path=mock_pipeline_data["experiments"],
        output_path=out_file,
        min_studies=1,
    )
    assert result == out_file
    assert out_file.exists()
    assert out_file.stat().st_size > 1000


def test_figure2_generation(mock_pipeline_data: dict[str, Path], tmp_path: Path) -> None:
    """Verify Figure 2 generates valid PNG with dynamic M and metadata."""
    out_file = tmp_path / "fig2.png"
    result = generate_figure2_calibrated_ora_infographic(
        taxa_matrix_path=mock_pipeline_data["taxa_matrix"],
        experiments_path=mock_pipeline_data["experiments"],
        output_path=out_file,
        enrichment_results_path=mock_pipeline_data["enrichment_results"],
        depletion_results_path=mock_pipeline_data["depletion_results"],
    )
    assert result == out_file
    assert out_file.exists()
    assert out_file.stat().st_size > 1000


def test_figure3_generation(mock_pipeline_data: dict[str, Path], tmp_path: Path) -> None:
    """Verify Figure 3 generates valid PNG with exact non-tied statistics."""
    out_file = tmp_path / "fig3.png"
    result = generate_figure3_lead_findings_ora_evidence(
        profiles_path=mock_pipeline_data["profiles"],
        experiments_path=mock_pipeline_data["experiments"],
        output_path=out_file,
    )
    assert result == out_file
    assert out_file.exists()
    assert out_file.stat().st_size > 1000


def test_figure4_generation(mock_pipeline_data: dict[str, Path], tmp_path: Path) -> None:
    """Verify Figure 4 generates valid PNG with real SVD PC1 and dynamic Jaccard."""
    out_file = tmp_path / "fig4.png"
    result = generate_figure4_annotation_density_pc1_landscape(
        taxa_matrix_path=mock_pipeline_data["taxa_matrix"],
        columns_path=mock_pipeline_data["columns"],
        profiles_path=mock_pipeline_data["profiles"],
        table1_effects_path=mock_pipeline_data["table1_effects"],
        output_path=out_file,
    )
    assert result == out_file
    assert out_file.exists()
    assert out_file.stat().st_size > 1000


def test_figure5_generation(mock_pipeline_data: dict[str, Path], tmp_path: Path) -> None:
    """Verify Figure 5 generates valid bubble map PNG."""
    out_file = tmp_path / "fig5.png"
    result = generate_figure5_pd_crc_taxon_attribution_bubble(
        bugsigdb_taxa_path=mock_pipeline_data["bugsigdb_taxa"],
        experiments_path=mock_pipeline_data["experiments"],
        taxa_matrix_path=mock_pipeline_data["taxa_matrix"],
        output_path=out_file,
    )
    assert result == out_file
    assert out_file.exists()
    assert out_file.stat().st_size > 1000


def _figure5_population_inputs(
    tmp_path: Path,
    *,
    excluded_host: str | None = "Mus musculus",
    excluded_method: str | None = "16S",
    excluded_site: str = "Feces",
    excluded_condition: str | None = None,
) -> dict[str, Path]:
    """Six eligible studies; excluded taxa dominate ranking and contaminate shared counts."""
    experiments: list[tuple[str, str, str, str, str | None, str | None]] = []
    occurrences: list[tuple[str, str, str, str, str, str]] = []
    for panel, condition in enumerate(("Parkinson's disease", "Colorectal cancer")):
        directions = (
            ("increased", "increased", "decreased")
            if panel == 0
            else ("decreased", "decreased", "increased")
        )
        valid_methods = ("16S", "WMS", " WMS " if panel == 0 else " 16S ")
        for index, (direction, method) in enumerate(zip(directions, valid_methods, strict=True)):
            study = f"valid-{panel}-{index}"
            host = " Homo sapiens " if index == 2 else "Homo sapiens"
            experiments.append((study, "e1", condition, "Feces", host, method))
            occurrences.append((study, "e1", direction, "NCBITaxon:1", "Eligible taxon", "species"))
        contaminating_direction = "decreased" if panel == 0 else "increased"
        for index in range(8):
            study = f"excluded-{panel}-{index}"
            experiments.append(
                (
                    study,
                    "e1",
                    excluded_condition or condition,
                    excluded_site,
                    excluded_host,
                    excluded_method,
                )
            )
            occurrences.append(
                (
                    study,
                    "e1",
                    contaminating_direction,
                    "NCBITaxon:2",
                    "Excluded-only taxon",
                    "species",
                )
            )
            if index < 3:
                occurrences.append(
                    (
                        study,
                        "e1",
                        contaminating_direction,
                        "NCBITaxon:1",
                        "Eligible taxon",
                        "species",
                    )
                )
    paths = {
        "experiments": tmp_path / "population_experiments.tsv",
        "bugsigdb_taxa": tmp_path / "population_taxa.tsv",
        "taxa_matrix": tmp_path / "population_matrix.tsv",
    }
    pd.DataFrame(
        experiments,
        columns=[
            "study_id",
            "experiment_id",
            "condition",
            "body_site",
            "host_species",
            "sequencing_type",
        ],
    ).to_csv(paths["experiments"], sep="\t", index=False)
    pd.DataFrame(
        occurrences,
        columns=[
            "study_id",
            "experiment_id",
            "direction",
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
        ],
    ).to_csv(paths["bugsigdb_taxa"], sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "taxon_id": taxon_id,
                "taxon_name": taxon_name,
                **dict.fromkeys(TABLE_1_TRAITS, 0),
                "METPO:2000103|GO:0030245": 1,
            }
            for taxon_id, taxon_name in (
                ("NCBITaxon:1", "Eligible taxon"),
                ("NCBITaxon:2", "Excluded-only taxon"),
            )
        ]
    ).to_csv(paths["taxa_matrix"], sep="\t", index=False)
    return paths


def _capture_figure5_population(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> dict[str, Any]:
    """Persist the real figure JSON in place of raster export; all calculations remain real."""
    import plotly.graph_objects as go

    output = tmp_path / "figure5.capture.json"
    before = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in inputs.items()}
    (tmp_path / "input_hashes.json").write_text(
        json.dumps(before, sort_keys=True), encoding="utf-8"
    )
    written: list[Path] = []

    def capture(figure: Any, path: str, **_kwargs: Any) -> None:
        assert Path(path) == output
        Path(path).write_text(figure.to_json(), encoding="utf-8")
        written.append(Path(path))

    monkeypatch.setattr(go.Figure, "write_image", capture)
    result = generate_figure5_pd_crc_taxon_attribution_bubble(
        bugsigdb_taxa_path=inputs["bugsigdb_taxa"],
        experiments_path=inputs["experiments"],
        taxa_matrix_path=inputs["taxa_matrix"],
        output_path=output,
        top_n_taxa=1,
    )
    assert result == output
    assert written == [output]
    assert output.stat().st_size > 0
    assert {
        name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in inputs.items()
    } == before
    figure = json.loads(output.read_text(encoding="utf-8"))
    assert isinstance(figure, dict)
    return figure


@pytest.mark.parametrize(
    ("host", "method", "site", "condition"),
    [
        pytest.param("Mus musculus", "16S", "Feces", None, id="mouse-16s"),
        pytest.param("Mus musculus", "WMS", "Feces", None, id="mouse-wms"),
        pytest.param(None, "16S", "Feces", None, id="missing-host"),
        pytest.param("unknown", "16S", "Feces", None, id="unknown-host"),
        pytest.param("Homo sapiens", "ITS", "Feces", None, id="human-its"),
        pytest.param("Homo sapiens", "ITS2", "Feces", None, id="human-its2"),
        pytest.param("Homo sapiens", "PCR", "Feces", None, id="human-pcr"),
        pytest.param("Homo sapiens", "18S", "Feces", None, id="human-18s"),
        pytest.param("Homo sapiens", None, "Feces", None, id="missing-method"),
        pytest.param("Homo sapiens", "untyped", "Feces", None, id="untyped-method"),
        pytest.param("Homo sapiens", "16s", "Feces", None, id="method-is-case-sensitive"),
        pytest.param("Homo sapiens", "16S", "Blood", None, id="blood"),
        pytest.param("Homo sapiens", "WMS", "Feces,Blood", None, id="mixed-body-site"),
        pytest.param("Homo sapiens", "16S", "Feces", "Asthma", id="other-condition"),
    ],
)
def test_figure5_filters_before_ranking_and_direction_counts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    host: str | None,
    method: str | None,
    site: str,
    condition: str | None,
) -> None:
    inputs = _figure5_population_inputs(
        tmp_path,
        excluded_host=host,
        excluded_method=method,
        excluded_site=site,
        excluded_condition=condition,
    )
    figure = _capture_figure5_population(inputs, tmp_path, monkeypatch)
    assert len(figure["data"]) == 2
    for trace, expected_direction in zip(figure["data"], (1 / 3, -1 / 3), strict=True):
        assert trace["y"] == ["Eligible taxon"]
        assert trace["x"] == [0]
        assert trace["marker"]["size"] == [13.5]
        assert trace["marker"]["color"] == pytest.approx([expected_direction])
        assert len(trace["text"]) == 1
        assert "Studies: 3<br>" in trace["text"][0]
        assert f"Net Direction: {expected_direction:+.2f}" in trace["text"][0]
    assert [annotation["text"] for annotation in figure["layout"]["annotations"]] == [
        "Parkinson's disease x Feces",
        "Colorectal cancer x Feces",
    ]


@pytest.mark.parametrize("missing_column", ["host_species", "sequencing_type"])
def test_figure5_requires_population_metadata(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, missing_column: str
) -> None:
    import plotly.graph_objects as go

    inputs = _figure5_population_inputs(tmp_path)
    experiments = pd.read_csv(inputs["experiments"], sep="\t")
    experiments.drop(columns=[missing_column]).to_csv(inputs["experiments"], sep="\t", index=False)
    before = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in inputs.items()}
    output = tmp_path / "figure5.capture.json"
    written: list[Path] = []

    def capture(figure: Any, path: str, **_kwargs: Any) -> None:
        Path(path).write_text(figure.to_json(), encoding="utf-8")
        written.append(Path(path))

    monkeypatch.setattr(go.Figure, "write_image", capture)
    with pytest.raises((KeyError, ValueError), match=missing_column):
        generate_figure5_pd_crc_taxon_attribution_bubble(
            bugsigdb_taxa_path=inputs["bugsigdb_taxa"],
            experiments_path=inputs["experiments"],
            taxa_matrix_path=inputs["taxa_matrix"],
            output_path=output,
            top_n_taxa=1,
        )
    assert written == []
    assert not output.exists()
    assert {
        name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in inputs.items()
    } == before


def test_generate_all_deck_figures(mock_pipeline_data: dict[str, Path], tmp_path: Path) -> None:
    """Verify end-to-end figure suite orchestration."""
    out_dir = tmp_path / "deck_out"
    paths = generate_all_deck_figures(
        derived_dir=mock_pipeline_data["derived_dir"],
        output_dir=out_dir,
    )
    assert len(paths) == 5
    for _key, p in paths.items():
        assert p.exists()
        assert p.stat().st_size > 1000


def test_cli_main(mock_pipeline_data: dict[str, Path], tmp_path: Path) -> None:
    """Verify CLI exit code and parameter handling."""
    out_dir = tmp_path / "cli_deck_out"
    exit_code = main(
        [
            "--derived-dir",
            str(mock_pipeline_data["derived_dir"]),
            "--output-dir",
            str(out_dir),
        ]
    )
    assert exit_code == 0
    assert (out_dir / "figure1_condition_anatomy_ora_heatmap.png").exists()
    assert (out_dir / "figure4_annotation_density_pc1_landscape.png").exists()


def test_figure4_draws_only_significant_adjusted_effects_as_surviving(
    mock_pipeline_data: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Only rows significant both raw and adjusted in the same direction are green (#1084-#1098)."""
    import plotly.graph_objects as go

    from mtp.density_adjustment import SIGNIFICANT_POST_ADJUSTMENT_STATUSES, adjustment_status

    # The classifier reads BH q-values (#1321): raw nonsignificant, adjusted significant.
    nonsig_raw_status = adjustment_status(-2.181818, raw_q=0.42, adj_mean=-5.834487, adj_q=0.0004)
    assert nonsig_raw_status == "significant only post-adjustment"
    assert nonsig_raw_status not in SIGNIFICANT_POST_ADJUSTMENT_STATUSES

    effects_path = mock_pipeline_data["table1_effects"]
    effects = pd.read_csv(effects_path, sep="\t")
    effects = pd.concat([effects, effects.iloc[[0]]], ignore_index=True)
    effects["status"] = [
        "survived (significant post-adjustment)",
        "enhanced (significant, larger magnitude)",
        nonsig_raw_status,
        "enhanced magnitude (nonsignificant post-adjustment)",
        "collapsed (nonsignificant post-adjustment, larger magnitude)",
        "reversed (significant post-adjustment)",
        "undefined magnitude change (raw mean zero, significant post-adjustment)",
    ]
    effects.to_csv(effects_path, sep="\t", index=False)
    figures: list[Any] = []

    def capture(figure: Any, path: str, **_kwargs: Any) -> None:
        figures.append(figure)
        Path(path).write_bytes(b"")

    monkeypatch.setattr(go.Figure, "write_image", capture)
    generate_figure4_annotation_density_pc1_landscape(
        taxa_matrix_path=mock_pipeline_data["taxa_matrix"],
        columns_path=mock_pipeline_data["columns"],
        profiles_path=mock_pipeline_data["profiles"],
        table1_effects_path=effects_path,
        output_path=tmp_path / "fig4.png",
    )

    (figure,) = figures
    green, red = "#27ae60", "#e74c3c"
    expected = [green, green, red, red, red, red, red]
    adjusted = next(trace for trace in figure.data if trace.name == "PC1-Adjusted")
    assert list(adjusted.marker.color) == expected
    dumbbells = [
        trace.line.color
        for trace in figure.data
        if trace.type == "scatter" and trace.mode == "lines" and len(trace.x) == 2
    ]
    assert dumbbells == expected


def test_facultative_union_label_in_deck_visualizations(
    mock_pipeline_data: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Deck visualizations render union label only from member metadata (#1138, #1232, #1235)."""
    import json

    import plotly.graph_objects as go

    from mtp.deck_visualizations import TABLE_1_TRAITS, TRAIT_SHORT_NAMES, clean_trait_label

    assert TABLE_1_TRAITS["METPO:1000608"] == "facultatively aerobic"
    assert TRAIT_SHORT_NAMES["METPO:1000608"] == "Facultatively aerobic"
    assert clean_trait_label("facultatively aerobic") == "Facultatively aerobic"
    assert (
        clean_trait_label("facultative (aerobe or anaerobe)") == "Facultative (aerobe or anaerobe)"
    )
    import inspect

    assert '"facultative (aerobe or anaerobe)"' not in inspect.getsource(clean_trait_label)

    figures: list[Any] = []

    def capture(figure: Any, path: str, **_kwargs: Any) -> None:
        figures.append(figure)
        Path(path).write_bytes(b"")

    monkeypatch.setattr(go.Figure, "write_image", capture)

    # Raw/legacy columns.tsv (no extra union members) preserves 'Facultatively aerobic' on Figure 4.
    generate_figure4_annotation_density_pc1_landscape(
        taxa_matrix_path=mock_pipeline_data["taxa_matrix"],
        columns_path=mock_pipeline_data["columns"],
        profiles_path=mock_pipeline_data["profiles"],
        table1_effects_path=mock_pipeline_data["table1_effects"],
        output_path=tmp_path / "fig4_raw.png",
    )
    raw_fig = figures.pop()
    raw_scatter_texts = [
        str(lbl)
        for trace in raw_fig.data
        if getattr(trace, "mode", None) == "markers+text"
        for lbl in (trace.text or [])
    ]
    assert "Facultatively aerobic" in raw_scatter_texts
    assert "Facultative (aerobe or anaerobe)" not in raw_scatter_texts

    # Union columns.tsv (with METPO:1000605 and METPO:1000612) renders the union label on Figure 4.
    cols_df = pd.read_csv(mock_pipeline_data["columns"], sep="\t")
    cols_df["trait_member_ids_json"] = [
        json.dumps(["METPO:1000605", "METPO:1000608", "METPO:1000612"])
        if cid == "METPO:1000608"
        else json.dumps([cid])
        for cid in cols_df["column_id"]
    ]
    union_cols_path = tmp_path / "union_columns.tsv"
    cols_df.to_csv(union_cols_path, sep="\t", index=False)

    generate_figure4_annotation_density_pc1_landscape(
        taxa_matrix_path=mock_pipeline_data["taxa_matrix"],
        columns_path=union_cols_path,
        profiles_path=mock_pipeline_data["profiles"],
        table1_effects_path=mock_pipeline_data["table1_effects"],
        output_path=tmp_path / "fig4_union.png",
    )
    union_fig = figures.pop()
    union_scatter_texts = [
        str(lbl)
        for trace in union_fig.data
        if getattr(trace, "mode", None) == "markers+text"
        for lbl in (trace.text or [])
    ]
    assert "Facultative (aerobe or anaerobe)" in union_scatter_texts
    assert "Facultatively aerobic" not in union_scatter_texts


def test_deck_figures_label_uncorrected_p_values_as_nominal_with_test_counts(
    mock_pipeline_data: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Figures 1, 3, and 4 label binomial, Wilcoxon, and Spearman p as nominal (#1335)."""
    import plotly.graph_objects as go

    figures: list[Any] = []

    def capture(figure: Any, path: str, **_kwargs: Any) -> None:
        figures.append(figure)
        Path(path).write_bytes(b"")

    monkeypatch.setattr(go.Figure, "write_image", capture)

    import re

    def _text_width_px(text: str, size: float) -> float:
        from matplotlib import font_manager
        from PIL import ImageFont

        fallback_px = len(text) * 0.58 * size
        measured_px = 0.0
        for font_path in (
            font_manager.findfont("DejaVu Sans", fallback_to_default=True),
            "/System/Library/Fonts/Supplemental/Arial.ttf",
        ):
            try:
                font = ImageFont.truetype(font_path, round(size))
                measured_px = max(measured_px, float(font.getlength(text)))
            except OSError:
                continue
        return max(measured_px, fallback_px)

    def assert_plotly_figure_text_fits(fig: Any) -> None:
        width = float(fig.layout.width)
        margin_l = float(fig.layout.margin.l or 80)
        margin_r = float(fig.layout.margin.r or 60)
        title_x = float(fig.layout.title.x or 0.02)
        title_avail_px = width * (1.0 - title_x)
        font_base = (
            float(fig.layout.title.font.size or 17.0)
            if getattr(fig.layout.title, "font", None)
            else 17.0
        )
        title_lines = str(fig.layout.title.text).split("<br>")
        for idx, raw_line in enumerate(title_lines):
            clean = re.sub(r"<[^>]+>", "", raw_line)
            font_px = font_base if idx == 0 else font_base * 0.7
            est_px = _text_width_px(clean, font_px)
            assert est_px <= title_avail_px, (
                f"Title line {idx} overflows canvas ({est_px:.1f}px > {title_avail_px:.1f}px): "
                f"{clean!r}"
            )

        plot_w = width - margin_l - margin_r
        for ann in fig.layout.annotations:
            ann_lines = str(ann.text).split("<br>")
            for raw_line in ann_lines:
                clean = re.sub(r"<[^>]+>", "", raw_line)
                est_px = _text_width_px(clean, 16.0)
                avail_px = plot_w * 0.36 if len(fig.layout.annotations) >= 3 else plot_w * 0.58
                assert est_px <= avail_px, (
                    f"Annotation line overflows subplot ({est_px:.1f}px > {avail_px:.1f}px): "
                    f"{clean!r}"
                )

        for trace in fig.data:
            if getattr(trace, "type", None) != "table":
                continue
            dom_x = trace.domain.x if trace.domain and trace.domain.x else (0.0, 1.0)
            table_px = (float(dom_x[1]) - float(dom_x[0])) * plot_w
            col_weights = [float(w) for w in trace.columnwidth]
            total_w = sum(col_weights)
            header_font_px = float(trace.header.font.size or 10.0)
            header_height = float(trace.header.height or 48.0)
            cell_font_px = float(trace.cells.font.size or 9.5)
            cell_height = float(trace.cells.height or 40.0)
            for col_idx, weight in enumerate(col_weights):
                col_px = (weight / total_w) * table_px
                h_lines = str(trace.header.values[col_idx]).split("<br>")
                assert len(h_lines) * header_font_px * 1.35 + 10.0 <= header_height, (
                    f"Table header col {col_idx} overflows height {header_height}px: {h_lines!r}"
                )
                for raw_line in h_lines:
                    clean = re.sub(r"<[^>]+>", "", raw_line)
                    est_px = _text_width_px(clean, header_font_px) + 10.0
                    assert est_px <= col_px, (
                        f"Table header col {col_idx} line overflows width "
                        f"({est_px:.1f}px > {col_px:.1f}px): {clean!r}"
                    )
                for row_val in trace.cells.values[col_idx]:
                    c_lines = str(row_val).split("<br>")
                    assert len(c_lines) * cell_font_px * 1.35 + 8.0 <= cell_height, (
                        f"Table cell col {col_idx} overflows height {cell_height}px: {c_lines!r}"
                    )
                    for raw_line in c_lines:
                        clean = re.sub(r"<[^>]+>", "", raw_line)
                        est_px = _text_width_px(clean, cell_font_px) + 10.0
                        assert est_px <= col_px, (
                            f"Table cell col {col_idx} line overflows width "
                            f"({est_px:.1f}px > {col_px:.1f}px): {clean!r}"
                        )

    generate_figure1_condition_anatomy_ora_heatmap(
        profiles_path=mock_pipeline_data["profiles"],
        experiments_path=mock_pipeline_data["experiments"],
        output_path=tmp_path / "fig1.png",
        min_studies=1,
    )
    fig1 = figures.pop()
    fig1_title = str(fig1.layout.title.text)
    assert "(* nominal p < 0.05, ** nominal p < 0.01, " in fig1_title
    assert (
        "*** nominal p < 0.001;<br>4 nominal tests across 4 x 6 condition-site x trait cells)"
        in fig1_title
    )
    assert "24 nominal tests" not in fig1_title
    assert_plotly_figure_text_fits(fig1)

    generate_figure3_lead_findings_ora_evidence(
        profiles_path=mock_pipeline_data["profiles"],
        experiments_path=mock_pipeline_data["experiments"],
        output_path=tmp_path / "fig3.png",
    )
    fig3 = figures.pop()
    fig3_title = str(fig3.layout.title.text)
    assert "8 nominal tests: 4 Wilcoxon + 4 binomial, uncorrected" in fig3_title
    table_trace = next(trace for trace in fig3.data if trace.type == "table")
    header_vals = [str(v) for v in table_trace.header.values]
    assert "Nominal Wilcoxon p" in header_vals[2]
    assert "Nominal Binomial p" in header_vals[3]
    for cell in table_trace.cells.values[2]:
        assert "nominal p=" in str(cell)
        assert str(cell).count("<br>") == 2
    for cell in table_trace.cells.values[3]:
        assert "nominal p=" in str(cell)
        assert str(cell).count("<br>") == 1
    assert_plotly_figure_text_fits(fig3)

    # Also verify Figure 3 when one lead condition is absent (untested) and one is significant
    exp_df = pd.read_csv(mock_pipeline_data["experiments"], sep="\t")
    prof_df = pd.read_csv(mock_pipeline_data["profiles"], sep="\t")
    exp_df = exp_df[exp_df["condition"] != "HIV infection"].copy()
    trait_cols = list(TABLE_1_TRAITS.keys())
    extra_exp_rows = []
    extra_prof_rows = []
    for idx in range(3, 8):
        s_id = f"study_pd_{idx}"
        e_id = f"exp_pd_{idx}"
        extra_exp_rows.append(
            {
                "study_id": s_id,
                "experiment_id": e_id,
                "condition": "Parkinson's disease",
                "body_site": "Feces",
                "host_species": "Homo sapiens",
            }
        )
        inc_vals = dict.fromkeys(trait_cols, 5.0)
        inc_vals["METPO:2000103|GO:0030245"] = 1.0
        dec_vals = dict.fromkeys(trait_cols, 5.0)
        dec_vals["METPO:2000103|GO:0030245"] = 10.0 + idx * 0.1
        extra_prof_rows.append(
            {"study_id": s_id, "experiment_id": e_id, "inc_dec": "inc", **inc_vals}
        )
        extra_prof_rows.append(
            {"study_id": s_id, "experiment_id": e_id, "inc_dec": "dec", **dec_vals}
        )
    exp_partial_path = tmp_path / "exp_partial.tsv"
    prof_partial_path = tmp_path / "prof_partial.tsv"
    pd.concat([exp_df, pd.DataFrame(extra_exp_rows)], ignore_index=True).to_csv(
        exp_partial_path, sep="\t", index=False
    )
    pd.concat([prof_df, pd.DataFrame(extra_prof_rows)], ignore_index=True).to_csv(
        prof_partial_path, sep="\t", index=False
    )
    generate_figure3_lead_findings_ora_evidence(
        profiles_path=prof_partial_path,
        experiments_path=exp_partial_path,
        output_path=tmp_path / "fig3_partial.png",
    )
    fig3_partial = figures.pop()
    fig3_partial_title = str(fig3_partial.layout.title.text)
    assert "6 nominal tests: 3 Wilcoxon + 3 binomial, uncorrected" in fig3_partial_title
    partial_scatter = next(trace for trace in fig3_partial.data if trace.type == "scatter")
    assert pd.isna(partial_scatter.x[1])
    assert partial_scatter.text[1] == ""
    assert "50.0%" not in [str(t) for t in partial_scatter.text]
    partial_table = next(trace for trace in fig3_partial.data if trace.type == "table")
    assert "(not tested)" in str(partial_table.cells.values[2][1])
    assert "Mean: NA" in str(partial_table.cells.values[2][1])
    assert "Mean: +0.00" not in str(partial_table.cells.values[2][1])
    assert "(not tested)" in str(partial_table.cells.values[3][1])
    assert "1.0000" not in str(partial_table.cells.values[2][1])
    assert "1.0000" not in str(partial_table.cells.values[3][1])
    assert partial_table.cells.values[4][0] == "Concordant (nominal p < 0.05)"
    assert partial_table.cells.values[4][1] == "Not tested"
    assert partial_table.cells.values[4][2] == "Directionally Aligned"
    assert_plotly_figure_text_fits(fig3_partial)

    generate_figure4_annotation_density_pc1_landscape(
        taxa_matrix_path=mock_pipeline_data["taxa_matrix"],
        columns_path=mock_pipeline_data["columns"],
        profiles_path=mock_pipeline_data["profiles"],
        table1_effects_path=mock_pipeline_data["table1_effects"],
        output_path=tmp_path / "fig4.png",
    )
    fig4 = figures.pop()
    fig4_title = str(fig4.layout.title.text)
    assert "1 nominal Spearman test across 6 traits" in fig4_title
    panel_a_title = str(fig4.layout.annotations[0].text)
    assert "nominal p =" in panel_a_title
    assert "<br>1 nominal test, n = 6 traits)" in panel_a_title
    assert_plotly_figure_text_fits(fig4)

    # Simulate Linux CI where macOS Arial and TrueType loading fail (#1512)
    from PIL import ImageFont

    real_truetype = ImageFont.truetype

    def _no_macos_fonts(font: Any = None, *args: Any, **kwargs: Any) -> Any:
        if "/System/Library/Fonts" in str(font):
            raise OSError(f"simulated Linux: no {font}")
        return real_truetype(font, *args, **kwargs)

    monkeypatch.setattr(ImageFont, "truetype", _no_macos_fonts)
    for rendered_fig in (fig1, fig3, fig3_partial, fig4):
        assert_plotly_figure_text_fits(rendered_fig)

    def _no_fonts_at_all(font: Any = None, *_args: Any, **_kwargs: Any) -> Any:
        raise OSError(f"simulated fallback: no {font}")

    monkeypatch.setattr(ImageFont, "truetype", _no_fonts_at_all)
    for rendered_fig in (fig1, fig3, fig3_partial, fig4):
        assert_plotly_figure_text_fits(rendered_fig)


def test_analyze_pd_trait_profiles_labels_nominal_p_values_and_test_counts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PD trait-profile script labels Mann-Whitney and Wilcoxon p as nominal (#1335)."""
    import importlib.util
    import sys

    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "mplconfig"))
    script_path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "analyze_parkinsons_disease_trait_profiles.py"
    )
    spec = importlib.util.spec_from_file_location("analyze_pd_trait_profiles", script_path)
    assert spec is not None
    assert spec.loader is not None
    pd_profiles = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = pd_profiles
    spec.loader.exec_module(pd_profiles)

    scores = pd.DataFrame(
        {
            "study_id": ["S1", "S1", "S2", "S2", "S3", "S3"],
            "experiment_id": ["E1", "E1", "E2", "E2", "E3", "E3"],
            "inc_dec": ["inc", "dec", "inc", "dec", "inc", "dec"],
            "aerobe_v2_score": [0.8, -0.2, 0.6, -0.4, 0.7, -0.1],
            "anaerobe_v2_score": [-0.5, 0.7, -0.3, 0.9, -0.4, 0.6],
            "aerobe_v2_positive_traits": [3, 1, 3, 1, 2, 1],
            "anaerobe_v2_positive_traits": [1, 3, 1, 4, 1, 3],
        }
    )
    paired = pd.DataFrame(
        {
            "study_id": ["S1", "S2", "S3"],
            "experiment_id": ["E1", "E2", "E3"],
            "aerobe_v2_increased_minus_decreased": [1.0, 1.0, 0.8],
            "anaerobe_v2_decreased_minus_increased": [1.2, 1.2, 1.0],
        }
    )
    tests = pd_profiles.summarize_paired_tests(scores, paired)
    assert tests["p_value_type"].tolist() == ["nominal", "nominal", "nominal", "nominal"]
    assert tests["n_tests"].tolist() == [4, 4, 4, 4]

    monkeypatch.setattr(pd_profiles, "OUTPUT_DIR", tmp_path)
    direction_summary = pd.DataFrame([{"inc_dec": "inc", "profile_count": 3}])
    taxa = pd.DataFrame(
        {
            "analysis": [pd_profiles.AEROBE_CONFIG.key, pd_profiles.ANAEROBE_CONFIG.key],
            "target_taxon_name": ["Streptococcus", "Roseburia"],
            "decreased_profile_count": [1, 6],
            "increased_profile_count": [5, 1],
            "decreased_study_count": [1, 3],
            "increased_study_count": [3, 1],
        }
    )
    taxa_summary = pd.DataFrame([{"guild_analysis": pd_profiles.AEROBE_CONFIG.key}])
    pd_profiles.write_readme(scores, direction_summary, tests, paired, taxa, taxa_summary)
    readme = (tmp_path / "README.md").read_text(encoding="utf-8")
    assert readme.count("nominal Mann-Whitney p=") == 2
    assert readme.count("(uncorrected; 2 profile-level tests, 4 tests in table).") == 2
    assert readme.count("nominal Wilcoxon p=") == 2
    assert readme.count("(uncorrected; 2 paired-contrast tests, 4 tests in table).") == 2
    assert (
        "(2 profile-level tests). Paired contrasts use two-sided Wilcoxon signed-rank "
        "tests against zero (2 paired-contrast tests; 4 nominal tests in total)."
    ) in readme
    assert (
        "- `parkinsons_feces_trait_score_tests.tsv`: profile-level and paired tests "
        "(4 nominal tests, uncorrected)."
    ) in readme

    # Verify unpaired / zero-contrast edge cases do not overcount untested rows
    zero_paired = pd.DataFrame(
        {
            "study_id": ["S1"],
            "experiment_id": ["E1"],
            "aerobe_v2_increased_minus_decreased": [0.0],
            "anaerobe_v2_decreased_minus_increased": [0.0],
        }
    )
    tests_partial = pd_profiles.summarize_paired_tests(scores, zero_paired)
    assert tests_partial["p_value_type"].tolist() == [
        "nominal",
        "not_tested",
        "nominal",
        "not_tested",
    ]
    assert tests_partial["n_tests"].tolist() == [2, 2, 2, 2]
    pd_profiles.write_readme(
        scores, direction_summary, tests_partial, zero_paired, taxa, taxa_summary
    )
    readme_partial = (tmp_path / "README.md").read_text(encoding="utf-8")
    assert readme_partial.count("(uncorrected; 2 profile-level tests, 2 tests in table).") == 2
    assert "(2 profile-level tests; 2 nominal tests in total)." in readme_partial
    assert "(2 nominal tests, uncorrected)." in readme_partial
    assert "p=NA" not in readme_partial

    # Verify mixed tested/untested Mann-Whitney and Wilcoxon rows print 'not tested' (#1518)
    mixed_scores = scores.copy()
    mixed_scores.loc[mixed_scores["inc_dec"] == "inc", "aerobe_v2_score"] = float("nan")
    mixed_paired = pd.DataFrame(
        {
            "study_id": ["S1", "S2", "S3"],
            "experiment_id": ["E1", "E2", "E3"],
            "aerobe_v2_increased_minus_decreased": [0.0, 0.0, 0.0],
            "anaerobe_v2_decreased_minus_increased": [1.2, 1.2, 1.0],
        }
    )
    tests_mixed = pd_profiles.summarize_paired_tests(mixed_scores, mixed_paired)
    assert tests_mixed["p_value_type"].tolist() == [
        "not_tested",
        "not_tested",
        "nominal",
        "nominal",
    ]
    assert tests_mixed["n_tests"].tolist() == [2, 2, 2, 2]
    pd_profiles.write_readme(
        mixed_scores, direction_summary, tests_mixed, mixed_paired, taxa, taxa_summary
    )
    readme_mixed = (tmp_path / "README.md").read_text(encoding="utf-8")
    assert "Mann-Whitney not tested (1 profile-level tests, 2 tests in table)." in readme_mixed
    assert "Wilcoxon not tested (1 paired-contrast tests, 2 tests in table)." in readme_mixed
    assert "p=NA" not in readme_mixed


def test_figure3_lead_concordance_depends_on_shift_sign(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """#1413, #1615, #1616: Figure 3 lead status uses Wilcoxon direction, Mixed, and < 0.05 pins."""
    import math

    import plotly.graph_objects as go
    from scipy import stats

    from mtp.deck_visualizations import lead_consistency_status, wilcoxon_signed_rank_diff

    assert (
        lead_consistency_status("Decreased", -2.5, 0.008, 0.9, 0.015)
        == "Concordant (nominal p < 0.05)"
    )
    assert (
        lead_consistency_status("Decreased", 2.5, 0.008, 0.1, 0.015)
        == "Discordant (nominal p < 0.05)"
    )
    assert (
        lead_consistency_status("Increased", -2.5, 0.008, 0.6, 0.5)
        == "Discordant (nominal p < 0.05)"
    )
    assert (
        lead_consistency_status("Increased", 0.5, 0.4, 0.1, 0.015)
        == "Discordant (nominal p < 0.05)"
    )
    assert lead_consistency_status("Decreased", -0.5, 0.4, 3 / 7, 1.0) == "Mixed"
    assert lead_consistency_status("Decreased", 0.5, 0.4, 4 / 7, 1.0) == "Mixed"
    assert lead_consistency_status("Decreased", 0.5, 0.4, 3 / 7, 1.0) == "Opposite Direction"
    assert lead_consistency_status("Decreased", 0.5, 0.4, 0.5, 1.0) == "Opposite Direction"
    assert lead_consistency_status("Decreased", 0.0, 1.0, 3 / 7, 1.0) == "Opposite Direction"
    assert lead_consistency_status("Decreased", -0.5, 0.4, 4 / 7, 1.0) == "Directionally Aligned"
    assert lead_consistency_status("Decreased", -0.5, 0.4, 0.5, 1.0) == "Directionally Aligned"
    assert lead_consistency_status("Decreased", 0.0, 1.0, 0.5, 1.0) == "No Net Shift"
    assert (
        lead_consistency_status("Decreased", float("nan"), float("nan"), float("nan"), float("nan"))
        == "Not tested"
    )
    with pytest.raises(ValueError, match="Unexpected target_dir"):
        lead_consistency_status("Flat", 0.0, 1.0, 0.5, 1.0)

    # #1616: boundary pins at nominal p = 0.049, 0.050, 0.051 for both tests and both directions,
    # plus precedence when the two tests are significant in opposite directions.
    assert (
        lead_consistency_status("Increased", 1.0, 0.049, 0.6, 0.50)
        == "Concordant (nominal p < 0.05)"
    )
    assert lead_consistency_status("Increased", 1.0, 0.050, 0.6, 0.50) == "Directionally Aligned"
    assert lead_consistency_status("Increased", 1.0, 0.051, 0.6, 0.50) == "Directionally Aligned"
    assert (
        lead_consistency_status("Increased", 1.0, 0.50, 0.8, 0.049)
        == "Concordant (nominal p < 0.05)"
    )
    assert lead_consistency_status("Increased", 1.0, 0.50, 0.8, 0.050) == "Directionally Aligned"
    assert lead_consistency_status("Increased", 1.0, 0.50, 0.8, 0.051) == "Directionally Aligned"
    assert (
        lead_consistency_status("Increased", -1.0, 0.049, 0.4, 0.50)
        == "Discordant (nominal p < 0.05)"
    )
    assert lead_consistency_status("Increased", -1.0, 0.050, 0.4, 0.50) == "Opposite Direction"
    assert lead_consistency_status("Increased", -1.0, 0.051, 0.4, 0.50) == "Opposite Direction"
    assert (
        lead_consistency_status("Increased", -1.0, 0.50, 0.2, 0.049)
        == "Discordant (nominal p < 0.05)"
    )
    assert lead_consistency_status("Increased", -1.0, 0.50, 0.2, 0.050) == "Opposite Direction"
    assert lead_consistency_status("Increased", -1.0, 0.50, 0.2, 0.051) == "Opposite Direction"
    assert (
        lead_consistency_status("Increased", 2.0, 0.001, 0.1, 0.001)
        == "Discordant (nominal p < 0.05)"
    )
    assert (
        lead_consistency_status("Increased", -2.0, 0.001, 0.9, 0.001)
        == "Discordant (nominal p < 0.05)"
    )

    # #1615: reviewer's outlier and mixed study-shift vectors.
    outlier_inc = [0.1] * 20 + [-5.0]
    assert sum(outlier_inc) / len(outlier_inc) < 0
    assert wilcoxon_signed_rank_diff(outlier_inc) == pytest.approx(189.0)
    assert (
        lead_consistency_status(
            "Increased",
            wilcoxon_signed_rank_diff(outlier_inc),
            float(stats.wilcoxon(outlier_inc).pvalue),
            20 / 21,
            float(stats.binomtest(20, 21, p=0.5).pvalue),
        )
        == "Concordant (nominal p < 0.05)"
    )
    outlier_dec = [-0.1] * 20 + [5.0]
    assert sum(outlier_dec) / len(outlier_dec) > 0
    assert wilcoxon_signed_rank_diff(outlier_dec) == pytest.approx(-189.0)
    assert (
        lead_consistency_status(
            "Decreased",
            wilcoxon_signed_rank_diff(outlier_dec),
            float(stats.wilcoxon(outlier_dec).pvalue),
            20 / 21,
            float(stats.binomtest(20, 21, p=0.5).pvalue),
        )
        == "Concordant (nominal p < 0.05)"
    )
    mixed_inc = [3.0] * 9 + [-0.1] * 11
    assert wilcoxon_signed_rank_diff(mixed_inc) == pytest.approx(78.0)
    assert (
        lead_consistency_status(
            "Increased",
            wilcoxon_signed_rank_diff(mixed_inc),
            float(stats.wilcoxon(mixed_inc).pvalue),
            9 / 20,
            float(stats.binomtest(9, 20, p=0.5).pvalue),
        )
        == "Mixed"
    )
    assert wilcoxon_signed_rank_diff([0.0, 0.0]) == 0.0
    assert math.isnan(wilcoxon_signed_rank_diff([]))

    # #1643: pin zero-shift exclusion, average tie ranks, and Concordant-before-Mixed precedence.
    assert wilcoxon_signed_rank_diff([5, 4, -0.1, -0.2, -0.3, 0, 0, 0, 0]) == pytest.approx(3.0)
    assert wilcoxon_signed_rank_diff([1, -1]) == pytest.approx(0.0)
    assert (
        lead_consistency_status("Increased", 1.0, 0.01, 0.4, 0.6) == "Concordant (nominal p < 0.05)"
    )
    assert (
        lead_consistency_status("Decreased", -1.0, 0.01, 0.4, 0.6)
        == "Concordant (nominal p < 0.05)"
    )
    assert (
        lead_consistency_status("Increased", -1.0, 0.6, 0.8, 0.01)
        == "Concordant (nominal p < 0.05)"
    )
    assert (
        lead_consistency_status("Decreased", 1.0, 0.6, 0.8, 0.01) == "Concordant (nominal p < 0.05)"
    )

    figures: list[Any] = []

    def _capture(fig: Any, path: str, **_kw: Any) -> None:
        figures.append(fig)
        Path(path).write_bytes(b"")

    monkeypatch.setattr(go.Figure, "write_image", _capture)

    trait_cols = list(TABLE_1_TRAITS.keys())
    exp_rows = []
    prof_rows = []
    # 21 PD studies:
    # - Cellulose (expected Decreased) shifted UP in all 21 (opposite, p < 0.05 -> Discordant);
    # - Moderately halophilic (expected Increased) has 20 studies at +0.1 and 1 outlier at -5.0
    #   (negative mean -0.14, but W+ > W- and 20/21 sign test -> Concordant, #1615).
    for idx in range(21):
        s_id = f"study_pd_{idx}"
        e_id = f"exp_pd_{idx}"
        exp_rows.append((s_id, e_id, "Parkinson's disease", "Feces", "Homo sapiens"))
        inc_vals = dict.fromkeys(trait_cols, 5.0)
        dec_vals = dict.fromkeys(trait_cols, 5.0)
        inc_vals["METPO:2000103|GO:0030245"] = 10.0 + idx * 0.1
        dec_vals["METPO:2000103|GO:0030245"] = 1.0
        halo_shift = 0.1 if idx < 20 else -5.0
        inc_vals["METPO:1000623"] = 5.0 + halo_shift
        dec_vals["METPO:1000623"] = 5.0
        prof_rows.append({"study_id": s_id, "experiment_id": e_id, "inc_dec": "inc", **inc_vals})
        prof_rows.append({"study_id": s_id, "experiment_id": e_id, "inc_dec": "dec", **dec_vals})

    # 2 HIV studies: Cellulose (expected Decreased) shifted UP
    # (opposite, p >= 0.05 -> Opposite Direction).
    for idx in range(2):
        s_id = f"study_hiv_{idx}"
        e_id = f"exp_hiv_{idx}"
        exp_rows.append((s_id, e_id, "HIV infection", "Feces", "Homo sapiens"))
        inc_vals = dict.fromkeys(trait_cols, 5.0)
        dec_vals = dict.fromkeys(trait_cols, 5.0)
        inc_vals["METPO:2000103|GO:0030245"] = 8.0 + idx * 0.1
        dec_vals["METPO:2000103|GO:0030245"] = 1.0
        prof_rows.append({"study_id": s_id, "experiment_id": e_id, "inc_dec": "inc", **inc_vals})
        prof_rows.append({"study_id": s_id, "experiment_id": e_id, "inc_dec": "dec", **dec_vals})

    # 20 CRC studies: Cellulose (expected Decreased) has 9 studies at -3.0 and 11 at +0.1
    # (Wilcoxon aligned W- > W+ p = 0.13, sign test 9/20 opposite p = 0.82 -> Mixed, #1615).
    for idx in range(20):
        s_id = f"study_crc_{idx}"
        e_id = f"exp_crc_{idx}"
        exp_rows.append((s_id, e_id, "Colorectal cancer", "Feces", "Homo sapiens"))
        inc_vals = dict.fromkeys(trait_cols, 5.0)
        dec_vals = dict.fromkeys(trait_cols, 5.0)
        crc_shift = -3.0 if idx < 9 else 0.1
        inc_vals["METPO:2000103|GO:0030245"] = 5.0 + crc_shift
        dec_vals["METPO:2000103|GO:0030245"] = 5.0
        prof_rows.append({"study_id": s_id, "experiment_id": e_id, "inc_dec": "inc", **inc_vals})
        prof_rows.append({"study_id": s_id, "experiment_id": e_id, "inc_dec": "dec", **dec_vals})

    exp_path = tmp_path / "exp_sign.tsv"
    prof_path = tmp_path / "prof_sign.tsv"
    pd.DataFrame(
        exp_rows,
        columns=["study_id", "experiment_id", "condition", "body_site", "host_species"],
    ).to_csv(exp_path, sep="\t", index=False)
    pd.DataFrame(prof_rows).to_csv(prof_path, sep="\t", index=False)

    generate_figure3_lead_findings_ora_evidence(
        profiles_path=prof_path,
        experiments_path=exp_path,
        output_path=tmp_path / "fig3_sign.png",
    )
    fig = figures.pop()
    table_trace = next(trace for trace in fig.data if trace.type == "table")
    mean_cells = list(table_trace.cells.values[2])
    statuses = list(table_trace.cells.values[4])
    assert "Mean: -0.14" in mean_cells[3]
    assert statuses == [
        "Discordant (nominal p < 0.05)",
        "Opposite Direction",
        "Mixed",
        "Concordant (nominal p < 0.05)",
    ]

    import mtp.deck_visualizations as dv

    monkeypatch.setattr(dv, "EXPLORATORY_LEADS", dv.EXPLORATORY_LEADS[:3])
    generate_figure3_lead_findings_ora_evidence(
        profiles_path=prof_path,
        experiments_path=exp_path,
        output_path=tmp_path / "fig3_sign_3.png",
    )
    fig3 = figures.pop()
    table3 = next(trace for trace in fig3.data if trace.type == "table")
    assert list(table3.cells.fill.color[0]) == ["#ffffff", "#f8f9fa", "#ffffff"]
    assert list(table3.cells.fill.color[4]) == ["#e8f5e9", "#e8f5e9", "#e8f5e9"]
