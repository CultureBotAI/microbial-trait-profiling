import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd

# Load module dynamically from scripts/
spec = importlib.util.spec_from_file_location(
    "generate_six_trait_sankey", "scripts/generate_six_trait_sankey.py"
)
assert spec is not None
assert spec.loader is not None
generate_six_trait_sankey = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generate_six_trait_sankey)

DEFAULT_TRAITS = generate_six_trait_sankey.DEFAULT_TRAITS
map_body_site = generate_six_trait_sankey.map_body_site
plot_grouped_six_trait_sankey = generate_six_trait_sankey.plot_grouped_six_trait_sankey


def test_map_body_site() -> None:
    assert map_body_site("feces") == "Feces"
    assert map_body_site("Stool sample") == "Feces"
    assert "Oral cavity" in map_body_site("saliva")
    assert "Vaginal" in map_body_site("cervicovaginal")
    assert "Upper airway" in map_body_site("nasopharynx")
    assert map_body_site(None) == "Other"


def test_plot_grouped_six_trait_sankey(tmp_path: Path) -> None:
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

    meta_df = pd.DataFrame(
        [
            {"experiment_id": "s1-1", "condition": "Parkinson's disease", "body_site": "Feces"},
            {"experiment_id": "s2-1", "condition": "Parkinson's disease", "body_site": "Feces"},
            {"experiment_id": "s3-1", "condition": "Colorectal cancer", "body_site": "Feces"},
            {"experiment_id": "s4-1", "condition": "Obesity", "body_site": "Saliva"},
        ]
    )

    trait_labels = {t: f"Label {i}" for i, t in enumerate(DEFAULT_TRAITS)}

    out_png_dec = tmp_path / "sankey_dec.png"
    out_svg_dec = tmp_path / "sankey_dec.svg"

    plot_grouped_six_trait_sankey(
        enrich_df=enrich_df,
        metadata_df=meta_df,
        selected_traits=DEFAULT_TRAITS,
        trait_labels=trait_labels,
        direction="dec",
        output_png=out_png_dec,
        output_svg=out_svg_dec,
        min_profiles=1,
        min_studies=1,
    )

    assert out_png_dec.exists()
    assert out_png_dec.stat().st_size > 0
    assert out_svg_dec.exists()
    assert out_svg_dec.stat().st_size > 0

    out_png_inc = tmp_path / "sankey_inc.png"
    plot_grouped_six_trait_sankey(
        enrich_df=enrich_df,
        metadata_df=meta_df,
        selected_traits=DEFAULT_TRAITS,
        trait_labels=trait_labels,
        direction="inc",
        output_png=out_png_inc,
        output_svg=None,
        min_profiles=1,
        min_studies=1,
    )

    assert out_png_inc.exists()
    assert out_png_inc.stat().st_size > 0


def test_cli_execution(tmp_path: Path) -> None:
    tsv_path = tmp_path / "test_enrichment.tsv"
    tsv_path.write_text(
        "study_id\texperiment_id\tprofile_id\tcondition\tbody_site\ttrait\tinc_dec\tp_value\tfold\n"
        f"s1\ts1-1\ts1|e1|dec\tParkinson's disease\tFeces\t{DEFAULT_TRAITS[0]}\tdec\t0.0001\t2.5\n"
        f"s2\ts2-1\ts2|e1|dec\tParkinson's disease\tFeces\t{DEFAULT_TRAITS[0]}\tdec\t0.0002\t2.5\n",
        encoding="utf-8",
    )

    out_dir = tmp_path / "sankey_out"

    subprocess.run(
        [
            sys.executable,
            "scripts/generate_six_trait_sankey.py",
            "--enrichment",
            str(tsv_path),
            "--output-dir",
            str(out_dir),
            "--direction",
            "dec",
            "--min-profiles",
            "1",
            "--min-studies",
            "1",
        ],
        capture_output=True,
        text=True,
        check=True,
    )

    assert (out_dir / "six_trait_sankey_decreased.png").exists()
    assert (out_dir / "six_trait_sankey_decreased.svg").exists()
