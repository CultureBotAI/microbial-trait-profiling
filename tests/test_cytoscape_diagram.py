import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd

# Load module dynamically from scripts/
spec = importlib.util.spec_from_file_location(
    "generate_cytoscape_diagram", "scripts/generate_cytoscape_diagram.py"
)
assert spec is not None
assert spec.loader is not None
generate_cytoscape_diagram = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generate_cytoscape_diagram)

DEFAULT_TRAITS = generate_cytoscape_diagram.DEFAULT_TRAITS
map_body_site = generate_cytoscape_diagram.map_body_site
plot_cytoscape_diagram = generate_cytoscape_diagram.plot_cytoscape_diagram


def test_map_body_site() -> None:
    assert map_body_site("feces") == "Feces"
    assert map_body_site("Stool sample") == "Feces"
    assert "Oral cavity" in map_body_site("saliva")
    assert "Vaginal" in map_body_site("cervicovaginal")
    assert "Upper airway" in map_body_site("nasopharynx")
    assert map_body_site(None) == "Other"


def test_plot_cytoscape_diagram(tmp_path: Path) -> None:
    enrich_df = pd.DataFrame(
        [
            {
                "study_id": "s1",
                "experiment_id": "s1-1",
                "profile_id": "s1|e1|dec",
                "condition": "Parkinson's disease",
                "body_site": "Feces",
                "trait": DEFAULT_TRAITS[0],
                "inc_dec": "dec",
                "p_value": 0.0001,
                "fold": 2.5,
            },
            {
                "study_id": "s2",
                "experiment_id": "s2-1",
                "profile_id": "s2|e1|dec",
                "condition": "Parkinson's disease",
                "body_site": "Feces",
                "trait": DEFAULT_TRAITS[0],
                "inc_dec": "dec",
                "p_value": 0.0005,
                "fold": 2.0,
            },
            {
                "study_id": "s1",
                "experiment_id": "s1-1",
                "profile_id": "s1|e1|dec",
                "condition": "Colorectal cancer",
                "body_site": "Feces",
                "trait": DEFAULT_TRAITS[1],
                "inc_dec": "dec",
                "p_value": 0.0002,
                "fold": 3.0,
            },
            {
                "study_id": "s3",
                "experiment_id": "s3-1",
                "profile_id": "s3|e1|dec",
                "condition": "Colorectal cancer",
                "body_site": "Feces",
                "trait": DEFAULT_TRAITS[1],
                "inc_dec": "dec",
                "p_value": 0.0003,
                "fold": 2.8,
            },
            {
                "study_id": "s1",
                "experiment_id": "s1-1",
                "profile_id": "s1|e1|inc",
                "condition": "Obesity",
                "body_site": "Saliva",
                "trait": DEFAULT_TRAITS[2],
                "inc_dec": "inc",
                "p_value": 0.001,
                "fold": 1.8,
            },
            {
                "study_id": "s4",
                "experiment_id": "s4-1",
                "profile_id": "s4|e1|inc",
                "condition": "Obesity",
                "body_site": "Saliva",
                "trait": DEFAULT_TRAITS[2],
                "inc_dec": "inc",
                "p_value": 0.0015,
                "fold": 2.1,
            },
        ]
    )

    metadata_df = pd.DataFrame(
        [
            {
                "study_id": "s1",
                "experiment_id": "s1-1",
                "condition": "Parkinson's disease",
                "body_site": "Feces",
            },
            {
                "study_id": "s2",
                "experiment_id": "s2-1",
                "condition": "Parkinson's disease",
                "body_site": "Feces",
            },
            {
                "study_id": "s3",
                "experiment_id": "s3-1",
                "condition": "Colorectal cancer",
                "body_site": "Feces",
            },
            {
                "study_id": "s4",
                "experiment_id": "s4-1",
                "condition": "Obesity",
                "body_site": "Saliva",
            },
        ]
    )

    trait_labels = {
        DEFAULT_TRAITS[0]: "facultatively\naerobic",
        DEFAULT_TRAITS[1]: "cellulose\ncatabolic",
        DEFAULT_TRAITS[2]: "coccus shaped",
    }

    out_dec_png = tmp_path / "cytoscape_six_traits_decreased.png"
    out_dec_svg = tmp_path / "cytoscape_six_traits_decreased.svg"

    plot_cytoscape_diagram(
        enrich_df=enrich_df,
        metadata_df=metadata_df,
        selected_traits=DEFAULT_TRAITS[:3],
        trait_labels=trait_labels,
        direction="dec",
        output_png=out_dec_png,
        output_svg=out_dec_svg,
        min_profiles=1,
        min_studies=1,
    )

    assert out_dec_png.exists()
    assert out_dec_png.stat().st_size > 0
    assert out_dec_svg.exists()
    assert out_dec_svg.stat().st_size > 0

    out_inc_png = tmp_path / "cytoscape_six_traits_increased.png"
    plot_cytoscape_diagram(
        enrich_df=enrich_df,
        metadata_df=metadata_df,
        selected_traits=DEFAULT_TRAITS[:3],
        trait_labels=trait_labels,
        direction="inc",
        output_png=out_inc_png,
        min_profiles=1,
        min_studies=1,
    )

    assert out_inc_png.exists()
    assert out_inc_png.stat().st_size > 0


def test_cli_execution(tmp_path: Path) -> None:
    dummy_enrich = tmp_path / "trait_enrichment_results.tsv"
    pd.DataFrame(
        [
            {
                "study_id": "s1",
                "experiment_id": "s1-1",
                "profile_id": "s1|e1|dec",
                "condition": "Parkinson's disease",
                "body_site": "Feces",
                "trait": DEFAULT_TRAITS[0],
                "inc_dec": "dec",
                "p_value": 0.0001,
                "fold": 2.5,
            },
            {
                "study_id": "s2",
                "experiment_id": "s2-1",
                "profile_id": "s2|e1|dec",
                "condition": "Parkinson's disease",
                "body_site": "Feces",
                "trait": DEFAULT_TRAITS[0],
                "inc_dec": "dec",
                "p_value": 0.0005,
                "fold": 2.0,
            },
            {
                "study_id": "s1",
                "experiment_id": "s1-1",
                "profile_id": "s1|e1|inc",
                "condition": "Obesity",
                "body_site": "Saliva",
                "trait": DEFAULT_TRAITS[1],
                "inc_dec": "inc",
                "p_value": 0.001,
                "fold": 1.8,
            },
            {
                "study_id": "s4",
                "experiment_id": "s4-1",
                "profile_id": "s4|e1|inc",
                "condition": "Obesity",
                "body_site": "Saliva",
                "trait": DEFAULT_TRAITS[1],
                "inc_dec": "inc",
                "p_value": 0.0015,
                "fold": 2.1,
            },
        ]
    ).to_csv(dummy_enrich, sep="\t", index=False)

    out_dir = tmp_path / "diagrams"
    cmd = [
        sys.executable,
        "scripts/generate_cytoscape_diagram.py",
        "--enrichment",
        str(dummy_enrich),
        "--output-dir",
        str(out_dir),
        "--direction",
        "both",
        "--min-profiles",
        "1",
        "--min-studies",
        "1",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    assert (out_dir / "cytoscape_six_traits_decreased.png").exists()
    assert (out_dir / "cytoscape_six_traits_increased.png").exists()
