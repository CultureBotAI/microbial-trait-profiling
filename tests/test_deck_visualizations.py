"""Tests for publication slide deck visualization suite (#347, #354)."""

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
    pd.DataFrame(enc_rows).to_csv(enc_path, sep="\t", index=False)

    dep_rows = [
        {"trait_id": t, "trait_label": TABLE_1_TRAITS[t], "q_value": 0.01} for t in trait_cols[2:]
    ]
    dep_path = derived_dir / "trait_depletion_results.tsv"
    pd.DataFrame(dep_rows).to_csv(dep_path, sep="\t", index=False)

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
    """Only significant effects that keep the raw sign are drawn as surviving (#1084, #1097)."""
    import plotly.graph_objects as go

    effects_path = mock_pipeline_data["table1_effects"]
    effects = pd.read_csv(effects_path, sep="\t")
    effects["status"] = [
        "survived (significant post-adjustment)",
        "enhanced (significant, larger magnitude)",
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
    expected = [green, green, red, red, red, red]
    adjusted = next(trace for trace in figure.data if trace.name == "PC1-Adjusted")
    assert list(adjusted.marker.color) == expected
    dumbbells = [
        trace.line.color
        for trace in figure.data
        if trace.type == "scatter" and trace.mode == "lines" and len(trace.x) == 2
    ]
    assert dumbbells == expected
