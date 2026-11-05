"""The PD cellulose analysis refuses stale six-shared-trait Sankey edge tables (#687)."""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mtp-mplconfig")
MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "create_pd_cellulose_analysis.py"
sys.path.insert(0, str(MODULE_PATH.parent))
SPEC = importlib.util.spec_from_file_location("pd_cellulose_analysis", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
cellulose = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = cellulose
SPEC.loader.exec_module(cellulose)


def _display(increased: tuple[int, int], decreased: tuple[int, int]) -> pd.DataFrame:
    rows = [
        (direction, cellulose.CELLULOSE_LABEL, cellulose.FECES, profiles, studies)
        for direction, (profiles, studies) in [("increased", increased), ("decreased", decreased)]
    ]
    rows.append(("decreased", "coccus shaped", cellulose.FECES, 99, 99))
    return pd.DataFrame(
        rows,
        columns=["direction", "trait", "body_site", "focus_enriched_profiles", "focus_studies"],
    )


def _recounted(increased: tuple[int, int], decreased: tuple[int, int]) -> pd.DataFrame:
    return pd.DataFrame(
        [("increased", *increased), ("decreased", *decreased)],
        columns=["direction", "focus_enriched_profiles", "focus_studies"],
    )


def test_matching_edges_pass() -> None:
    cellulose.check_display_edges_current(
        _display((34, 22), (69, 34)), _recounted((34, 22), (69, 34))
    )


@pytest.mark.parametrize(
    "recounted",
    [_recounted((34, 22), (70, 34)), _recounted((34, 21), (69, 34))],
)
def test_stale_edges_are_refused_with_the_recipe_to_run(recounted: pd.DataFrame) -> None:
    with pytest.raises(ValueError, match=r"just six-shared-trait-sankeys"):
        cellulose.check_display_edges_current(_display((34, 22), (69, 34)), recounted)


def test_a_missing_edge_is_refused() -> None:
    display = _display((34, 22), (69, 34))
    display = display[~display["direction"].eq("increased")]
    with pytest.raises(ValueError, match="stale"):
        cellulose.check_display_edges_current(display, _recounted((34, 22), (69, 34)))


def test_cellulose_tables_plots_and_interpretations_label_nominal_p_values(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Every uncorrected Fisher, Mann-Whitney, and OLS p-value is labelled nominal (#1335)."""
    rows = []
    z_rows = []
    for study_idx in range(6):
        for cond, base in [(cellulose.PD_CONDITION, 0.8), ("other condition", -0.4)]:
            for inc_dec, delta in [("inc", -0.3), ("dec", 0.9)]:
                rows.append(
                    {
                        "study_id": f"{cond}-{study_idx}",
                        "condition": cond,
                        "body_site_group": cellulose.FECES,
                        "inc_dec": inc_dec,
                    }
                )
                z_rows.append(
                    {
                        cellulose.CELLULOSE_ID: base + delta + 0.05 * study_idx,
                        "METPO:1000668": -base,
                    }
                )
    metadata = pd.DataFrame(rows)
    zscores = pd.DataFrame(z_rows)

    pd_vs_other = cellulose.build_pd_vs_other_cellulose_summary(metadata, zscores)
    assert pd_vs_other["p_value_type"].eq("nominal").all()
    assert pd_vs_other["n_tests"].eq(4).all()

    model_table, group_summary = cellulose.build_cellulose_model_tables(metadata, zscores)
    assert model_table["p_value_type"].eq("nominal").all()
    assert model_table["n_model_terms"].eq(4).all()
    assert model_table["n_tests"].eq(8).all()

    import matplotlib.figure
    import matplotlib.pyplot as plt

    captured_texts: dict[str, list[str]] = {}

    def capture_figure(fig: matplotlib.figure.Figure, output_path: Path) -> None:
        axes = fig.axes
        texts = [t.get_text() for ax in axes for t in ax.texts]
        texts.extend(ax.get_title() for ax in axes)
        texts.extend(ax.get_xlabel() for ax in axes)
        captured_texts[output_path.name] = texts
        plt.close(fig)

    monkeypatch.setattr(cellulose, "save_figure", capture_figure)

    pd_display = pd.DataFrame(
        [
            {
                "direction": "increased",
                "trait": cellulose.CELLULOSE_LABEL,
                "body_site": cellulose.FECES,
                "focus_enriched_profiles": 10,
                "comparator_enriched_profiles": 16,
                "focus_total_profiles": 20,
                "focus_studies": 6,
                "focus_enriched_mean_z": 0.45,
                "focus_minus_comparator_z": -0.25,
            },
            {
                "direction": "decreased",
                "trait": cellulose.CELLULOSE_LABEL,
                "body_site": cellulose.FECES,
                "focus_enriched_profiles": 16,
                "comparator_enriched_profiles": 10,
                "focus_total_profiles": 20,
                "focus_studies": 6,
                "focus_enriched_mean_z": 0.85,
                "focus_minus_comparator_z": 0.75,
            },
        ]
    )
    cellulose.plot_direction_specificity(pd_display, pd_vs_other, tmp_path / "dir.png")
    dir_joined = "\n".join(captured_texts["dir.png"])
    assert "nominal Fisher p=" in dir_joined
    assert "2 directions tested" in dir_joined

    cellulose.plot_enrichment_test(model_table, group_summary, tmp_path / "enrich.png")
    enrich_joined = "\n".join(captured_texts["enrich.png"])
    assert "nominal p" in enrich_joined
    assert "3 model terms tested" in enrich_joined
    assert "3 nominal tests" in enrich_joined

    trait_summary = pd.DataFrame(
        [
            {
                "trait": cellulose.CELLULOSE_LABEL,
                "increased_minus_decreased_mean_z": -0.6,
            },
            {
                "trait": "coccus shaped",
                "increased_minus_decreased_mean_z": 0.4,
            },
        ]
    )
    drivers = pd.DataFrame(
        [
            {
                "taxon_name": "Roseburia intestinalis",
                "taxonomic_rank": "species",
                "decreased_profile_count": 12,
                "increased_profile_count": 2,
                "decreased_study_count": 6,
                "increased_study_count": 2,
                "decreased_minus_increased_profiles": 10,
                "total_profile_count": 14,
            }
        ]
    )
    readme_path = tmp_path / "README.md"
    cellulose.write_interpretations(
        readme_path,
        pd_display,
        pd_vs_other,
        trait_summary,
        model_table,
        drivers,
    )
    readme = readme_path.read_text(encoding="utf-8")
    assert "nominal p=" in readme
    assert "4 nominal tests" in readme
    assert "3 non-intercept model terms" in readme
    assert "8 terms across 2 outcomes" in readme
