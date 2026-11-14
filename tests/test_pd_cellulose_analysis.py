"""The PD cellulose analysis refuses stale six-shared-trait Sankey edge tables (#687)."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pandas as pd
import pytest

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
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "mplconfig"))
    import numpy as np

    cell_z = {
        ("other condition", "inc"): [-0.8, -0.5, -0.3, -0.1, 0.2, 0.4],
        ("other condition", "dec"): [-0.6, -0.2, 0.1, -0.3, 0.3, 0.6],
        (cellulose.PD_CONDITION, "inc"): [-0.4, -0.1, 0.2, -0.2, 0.5, -0.3],
        (cellulose.PD_CONDITION, "dec"): [0.5, 0.8, 1.1, 0.9, -0.1, 1.3],
    }
    rows = []
    z_rows = []
    for study_idx in range(6):
        for cond, base in [(cellulose.PD_CONDITION, 0.8), ("other condition", -0.4)]:
            for inc_dec in ("inc", "dec"):
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
                        cellulose.CELLULOSE_ID: cell_z[(cond, inc_dec)][study_idx],
                        "METPO:1000668": -base,
                    }
                )
    metadata = pd.DataFrame(rows)
    zscores = pd.DataFrame(z_rows)

    pd_vs_other = cellulose.build_pd_vs_other_cellulose_summary(metadata, zscores)
    assert pd_vs_other["p_value_type"].eq("nominal").all()
    assert pd_vs_other["n_tests"].eq(4).all()

    # When one direction has no PD feces profiles, its 2 tests are not run (n_tests == 2)
    dec_only_mask = ~(
        metadata["condition"].eq(cellulose.PD_CONDITION) & metadata["inc_dec"].eq("inc")
    )
    pd_vs_other_partial = cellulose.build_pd_vs_other_cellulose_summary(
        metadata.loc[dec_only_mask].reset_index(drop=True),
        zscores.loc[dec_only_mask].reset_index(drop=True),
    )
    assert pd_vs_other_partial["n_tests"].tolist() == [2, 2, 2, 2]
    assert pd_vs_other_partial.loc[
        pd_vs_other_partial["direction"].eq("increased"), "p_value_type"
    ].tolist() == ["not_tested", "not_tested"]
    assert (
        pd_vs_other_partial.loc[
            pd_vs_other_partial["direction"].eq("increased"),
            ["pd_vs_other_fisher_p", "pd_vs_other_mannwhitney_p"],
        ]
        .isna()
        .all()
        .all()
    )
    assert pd_vs_other_partial.loc[
        pd_vs_other_partial["direction"].eq("decreased"), "p_value_type"
    ].tolist() == ["nominal", "nominal"]

    # When one PD increased row has a NaN continuous z-score,
    # Fisher runs while Mann-Whitney propagates NaN (n_tests == 3)
    zscores_nan_inc = zscores.copy()
    first_pd_inc = metadata.index[
        metadata["condition"].eq(cellulose.PD_CONDITION) & metadata["inc_dec"].eq("inc")
    ][0]
    zscores_nan_inc.loc[first_pd_inc, cellulose.CELLULOSE_ID] = float("nan")
    pd_vs_other_three = cellulose.build_pd_vs_other_cellulose_summary(metadata, zscores_nan_inc)
    assert pd_vs_other_three["n_tests"].tolist() == [3, 3, 3, 3]
    assert pd_vs_other_three["p_value_type"].eq("nominal").all()
    assert (
        pd_vs_other_three.loc[
            pd_vs_other_three["direction"].eq("increased"), "pd_vs_other_mannwhitney_p"
        ]
        .isna()
        .all()
    )

    model_table, group_summary = cellulose.build_cellulose_model_tables(metadata, zscores)
    assert model_table["p_value_type"].eq("nominal").all()
    assert model_table["n_model_terms"].eq(3).all()
    assert model_table["n_tests"].eq(8).all()

    # Row-order permutation stability across 25 shuffles (#1515)
    rng = np.random.default_rng(1515)
    for _ in range(25):
        perm = rng.permutation(len(metadata))
        perm_model, _ = cellulose.build_cellulose_model_tables(
            metadata.iloc[perm].reset_index(drop=True),
            zscores.iloc[perm].reset_index(drop=True),
        )
        assert perm_model["p_value_type"].eq("nominal").all()
        assert perm_model["n_tests"].eq(8).all()

    # Degenerate zero-variance non-zero outcome produces not_tested and NaN p-values (#1515)
    const_zscores = zscores.assign(**{cellulose.CELLULOSE_ID: 2.0})
    deg_model, _ = cellulose.build_cellulose_model_tables(metadata, const_zscores)
    assert deg_model["p_value_type"].eq("not_tested").all()
    assert deg_model["p_value"].isna().all()
    assert deg_model["n_tests"].eq(0).all()

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

    # Exercise both p < 0.001 ("nominal p<0.001") and p >= 0.001 ("nominal p=") branches
    model_for_plot = model_table.copy()
    model_for_plot.loc[
        model_for_plot["outcome"].eq("cellulose_z") & model_for_plot["term"].eq("is_pd"),
        "p_value",
    ] = 0.042
    model_for_plot.loc[
        model_for_plot["outcome"].eq("cellulose_z") & model_for_plot["term"].eq("pd_x_decreased"),
        "p_value",
    ] = 0.0001
    cellulose.plot_enrichment_test(model_for_plot, group_summary, tmp_path / "enrich.png")
    enrich_joined = "\n".join(captured_texts["enrich.png"])
    assert "nominal p<0.001" in enrich_joined
    assert "nominal p=0.042" in enrich_joined
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
    assert readme.count("nominal p=") == 4
    assert (
        "This is the interaction's nominal p, one of 3 non-intercept model terms tested "
        "in the figure (8 terms across 2 outcomes in `pd_cellulose_interaction_model.tsv`), "
        "not corrected for the other conditions and traits examined"
    ) in readme
    assert (
        "uncorrected across 2 directions x 2 tests = 4 nominal tests in "
        "`pd_cellulose_vs_other_feces_summary.tsv`"
    ) in readme
    assert (
        "uncorrected across 3 non-intercept model terms plotted in "
        "`pd_cellulose_enrichment_test.png` and 8 terms across 2 outcomes in "
        "`pd_cellulose_interaction_model.tsv`"
    ) in readme
    assert "It tests whether the decreased direction has" in readme


def _degenerate_consumer_tables(
    *,
    all_untested: bool,
    rank_deficient: bool = False,
    reverse_directions: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Use a real clustered model with different test availability across outcomes."""
    metadata_rows = []
    zscore_rows = []
    pd_inc_values = (-0.5, 0.2, -0.3, 0.4, -0.2, 0.3)
    other_inc_values = (-0.6, 0.3, -0.4, 0.2, -0.3, 0.4)
    for study in range(6):
        jitter = (1.0 if study % 2 == 0 else -1.0) + 0.05 * study
        for condition, base in [(cellulose.PD_CONDITION, 0.8), ("other condition", -0.4)]:
            for direction, delta in [("inc", -0.3), ("dec", 0.9)]:
                if rank_deficient and condition != cellulose.PD_CONDITION and direction == "dec":
                    continue
                metadata_rows.append(
                    {
                        "study_id": f"{condition}-{study}",
                        "condition": condition,
                        "body_site_group": cellulose.FECES,
                        "inc_dec": direction,
                    }
                )
                value = 0.0 if all_untested else base + delta + jitter
                if rank_deficient:
                    if direction == "dec":
                        value = 0.5 + 0.15 * study
                    elif condition == cellulose.PD_CONDITION:
                        value = pd_inc_values[study]
                    else:
                        value = other_inc_values[study]
                zscore_rows.append({cellulose.CELLULOSE_ID: value})
    metadata = pd.DataFrame(metadata_rows)
    if reverse_directions:
        metadata["inc_dec"] = metadata["inc_dec"].map({"inc": "dec", "dec": "inc"})
    zscores = pd.DataFrame(zscore_rows)
    model, groups = cellulose.build_cellulose_model_tables(metadata, zscores)
    comparisons = cellulose.build_pd_vs_other_cellulose_summary(metadata, zscores)
    return model, groups, comparisons


@pytest.mark.parametrize(
    ("all_untested", "rank_deficient"),
    [(False, False), (True, False), (False, True)],
    ids=["mixed-outcomes", "zero-tests", "rank-deficient"],
)
def test_degenerate_model_figure_labels_only_tested_cellulose_terms(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, all_untested: bool, rank_deficient: bool
) -> None:
    """An untested coefficient stays visible without becoming a nominal test (#1539)."""
    import matplotlib.figure
    import matplotlib.pyplot as plt

    model, groups, _ = _degenerate_consumer_tables(
        all_untested=all_untested, rank_deficient=rank_deficient
    )
    expected_tested = 0 if all_untested else 1
    plotted = model[model["outcome"].eq("cellulose_z") & ~model["term"].eq("intercept")]
    assert int(plotted["p_value_type"].eq("nominal").sum()) == expected_tested
    expected_total = 0 if all_untested else (4 if rank_deficient else 5)
    assert model["n_tests"].eq(expected_total).all()
    captured: list[str] = []
    unavailable_positions: list[tuple[float, float]] = []

    def capture(fig: matplotlib.figure.Figure, _output: Path) -> None:
        captured.extend(text.get_text() for ax in fig.axes for text in ax.texts)
        for ax in fig.axes:
            for label in ax.texts:
                if "not estimable" in label.get_text():
                    background = label.get_bbox_patch()
                    assert background is not None, "Keep the zero-reference line behind the label"
                    assert background.get_facecolor() == (1.0, 1.0, 1.0, 1.0)
        unavailable_positions.extend(
            text.get_position()
            for ax in fig.axes
            for text in ax.texts
            if "not estimable" in text.get_text()
        )
        captured.extend(ax.get_title() for ax in fig.axes)
        captured.extend(ax.get_xlabel() for ax in fig.axes)
        captured.extend(text.get_text() for text in fig.texts)
        plt.close(fig)

    monkeypatch.setattr(cellulose, "save_figure", capture)
    cellulose.plot_enrichment_test(model, groups, tmp_path / "degenerate.png")
    text = "\n".join(captured)
    assert "nominal p=nan" not in text
    assert "nominal p=inf" not in text
    assert text.count("not tested") == 3 - expected_tested
    assert text.count("nominal p=") + text.count("nominal p<") == expected_tested
    term_noun = "term" if expected_tested == 1 else "terms"
    test_noun = "test" if expected_tested == 1 else "tests"
    assert f"{expected_tested} model {term_noun} tested" in text
    assert f"({expected_tested} nominal {test_noun})" in text
    assert "1 model terms" not in text
    assert "1 nominal tests" not in text
    if rank_deficient:
        assert text.count("not estimable") == 2
        assert "nan" not in text.lower()
        assert len(unavailable_positions) == 2
        assert all(math.isfinite(value) for position in unavailable_positions for value in position)


@pytest.mark.parametrize(
    ("all_untested", "rank_deficient"),
    [(False, False), (True, False), (False, True)],
    ids=["mixed-outcomes", "zero-tests", "rank-deficient"],
)
def test_degenerate_model_interpretations_distinguish_displayed_and_total_tests(
    tmp_path: Path, all_untested: bool, rank_deficient: bool
) -> None:
    """The other outcome cannot make a cellulose-z term count as tested (#1539)."""
    model, _groups, comparisons = _degenerate_consumer_tables(
        all_untested=all_untested, rank_deficient=rank_deficient
    )
    expected_tested = 0 if all_untested else 1
    expected_total = 0 if all_untested else (4 if rank_deficient else 5)
    if not all_untested and not rank_deficient:
        nominal_non_intercepts = model[
            model["p_value_type"].eq("nominal") & ~model["term"].eq("intercept")
        ]
        assert nominal_non_intercepts["term"].nunique() == 2
    text = _write_synthetic_interpretations(tmp_path, model, comparisons)
    assert "nominal p=nan" not in text
    assert "nominal p=inf" not in text
    assert "nan%" not in text
    assert "odds ratio of `nan`" not in text
    assert "This is the interaction's nominal p" not in text
    assert "95% CI" not in text
    assert "interval includes zero" not in text
    assert "interval excludes zero" not in text
    availability_verb = "has" if expected_tested == 1 else "have"
    plotted_scope = (
        f"{expected_tested} of the 3 plotted cellulose-z terms {availability_verb} a nominal test"
    )
    assert f"The interaction was not tested; {plotted_scope}" in text
    assert f"{expected_total} of the 8 terms across 2 outcomes" in text
    model_paragraph = next(
        paragraph
        for paragraph in text.split("\n\n")
        if paragraph.startswith("The feces-only linear model")
    )
    summary_paragraph = next(
        paragraph
        for paragraph in text.split("\n\n")
        if paragraph.startswith("In PD fecal profiles the mean cellulose annotation z-score")
    )
    assert (summary_paragraph + model_paragraph).count("`not tested`") == 2
    assert f"{plotted_scope} in `pd_cellulose_enrichment_test.png`" in model_paragraph
    assert (
        f"{expected_total} of the 8 table terms have a nominal test across 2 outcomes in "
        "`pd_cellulose_interaction_model.tsv`"
    ) in model_paragraph
    assert "uncorrected across" not in model_paragraph
    assert "It tests whether" not in text
    comparison_section = text.split("## PD vs Other Fecal Profiles", 1)[1].split(
        "The feces-only linear model", 1
    )[0]
    if all_untested:
        assert comparisons["p_value_type"].eq("nominal").all()
        assert comparisons["pd_vs_other_fisher_p"].eq(1.0).all()
        assert "an odds ratio that is not estimable" in comparison_section
        assert comparison_section.count("nominal p=1") == 2
        assert "not tested" not in comparison_section
        assert "2 directions x 2 tests = 4 nominal tests" in comparison_section
    if rank_deficient:
        assert text.count("not estimable") == 2
        assert "nan" not in text.lower()
        decreased = comparisons[
            comparisons["direction"].eq("decreased") & comparisons["group"].eq("PD feces")
        ].iloc[0]
        assert decreased.p_value_type == "not_tested"
        assert decreased.positive_z_rate == 1.0
        increased = comparisons[
            comparisons["direction"].eq("increased") & comparisons["group"].eq("PD feces")
        ].iloc[0]
        assert increased.p_value_type == "nominal"
        assert (
            comparison_section.count(
                "The decreased-direction PD-vs-other comparison was not tested."
            )
            == 1
        )
        assert comparison_section.count("was not tested") == 1
        assert "The increased-direction PD-vs-other Fisher test" in comparison_section
        assert f"`nominal p={increased.pd_vs_other_fisher_p:.3g}`" in comparison_section
        assert (
            f"Mann-Whitney `nominal p={increased.pd_vs_other_mannwhitney_p:.3g}`"
        ) in comparison_section
        assert comparison_section.count("nominal p=") == 2
        assert "`100.0%`" not in comparison_section
        assert "1 direction x 2 tests = 2 nominal tests" in comparison_section
        assert "1 directions" not in comparison_section


def test_interpretations_omit_rates_for_untested_increased_direction(tmp_path: Path) -> None:
    """The other unavailable direction uses its own model and comparison inputs."""
    model, _groups, comparisons = _degenerate_consumer_tables(
        all_untested=False, rank_deficient=True, reverse_directions=True
    )
    increased = comparisons[
        comparisons["direction"].eq("increased") & comparisons["group"].eq("PD feces")
    ].iloc[0]
    assert increased.p_value_type == "not_tested"
    assert increased.positive_z_rate == 1.0
    text = _write_synthetic_interpretations(tmp_path, model, comparisons)
    comparison_section = text.split("## PD vs Other Fecal Profiles", 1)[1].split(
        "The feces-only linear model", 1
    )[0]
    assert (
        comparison_section.count("increased-direction PD-vs-other comparison was not tested") == 1
    )
    assert "Among decreased fecal profiles" in comparison_section
    assert "`100.0%`" not in comparison_section
    assert "nan" not in comparison_section.lower()
    assert comparison_section.count("nominal p=") == 2
    assert "1 direction x 2 tests = 2 nominal tests" in comparison_section
    assert "1 directions" not in comparison_section


def test_interpretations_do_not_invent_a_missing_comparison_test(tmp_path: Path) -> None:
    """A direction-level nominal flag does not make its missing second p available."""
    model, _groups, comparisons = _degenerate_consumer_tables(
        all_untested=False, rank_deficient=True
    )
    increased = comparisons["direction"].eq("increased")
    comparisons.loc[increased, "pd_vs_other_mannwhitney_p"] = math.nan
    comparisons.loc[:, "n_tests"] = 1
    original = comparisons.copy(deep=True)

    text = _write_synthetic_interpretations(tmp_path, model, comparisons)
    section = text.split("## PD vs Other Fecal Profiles", 1)[1].split(
        "The feces-only linear model", 1
    )[0]
    assert "The increased-direction PD-vs-other Fisher test" in section
    assert "`nominal p=1`" in section
    assert "Mann-Whitney p unavailable" in section
    assert "1 nominal test" in section
    assert "1 nominal tests" not in section
    assert "1 direction x 2 tests = 1" not in section
    assert "nominal p=nan" not in section
    assert section.count("was not tested") == 1
    pd.testing.assert_frame_equal(comparisons, original)


@pytest.mark.parametrize(
    "case",
    [
        "untested-decreased",
        "untested-increased",
        "tested-infinite-odds",
        "tested-undefined-odds",
    ],
)
def test_direction_figure_marks_availability_without_changing_bar_values(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, case: str
) -> None:
    """Consumer labels distinguish absent tests and boundary odds from plotted rates."""
    import matplotlib.figure
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.text import Annotation

    untested = case.startswith("untested-")
    _model, _groups, comparisons = _degenerate_consumer_tables(
        all_untested=case == "tested-undefined-odds",
        rank_deficient=untested,
        reverse_directions=case == "untested-increased",
    )
    if case == "tested-infinite-odds":
        decreased = comparisons["direction"].eq("decreased")
        comparisons.loc[decreased, "pd_vs_other_positive_z_odds_ratio"] = math.inf
        comparisons.loc[decreased, "pd_vs_other_fisher_p"] = 0.04
    display = pd.DataFrame(
        [
            {
                "direction": direction,
                "trait": cellulose.CELLULOSE_LABEL,
                "body_site": cellulose.FECES,
                "focus_enriched_profiles": focus,
                "comparator_enriched_profiles": comparator,
                "focus_studies": 3,
                "focus_enriched_mean_z": 0.2,
                "focus_minus_comparator_z": 0.1,
            }
            for direction, focus, comparator in (("increased", 4, 6), ("decreased", 6, 4))
        ]
    )
    original = comparisons.copy(deep=True)
    texts: list[str] = []
    bars: list[tuple[float, float]] = []
    labels: dict[float, str] = {}

    def capture(fig: matplotlib.figure.Figure, _output: Path) -> None:
        assert len(fig.axes) == 2
        right = fig.axes[1]
        texts.extend(text.get_text() for ax in fig.axes for text in ax.texts)
        texts.extend(ax.get_title() for ax in fig.axes)
        texts.extend(text.get_text() for text in fig.texts)
        bars.extend(
            (round(patch.get_x() + patch.get_width() / 2, 2), patch.get_height())
            for patch in right.patches
            if isinstance(patch, Rectangle)
        )
        labels.update(
            (round(float(label.xy[0]), 2), label.get_text())
            for label in right.texts
            if isinstance(label, Annotation)
        )
        plt.close(fig)

    monkeypatch.setattr(cellulose, "save_figure", capture)
    cellulose.plot_direction_specificity(display, comparisons, tmp_path / "direction.png")
    assert len(bars) == 4
    expected: dict[float, float] = {}
    for row in comparisons.to_dict("records"):
        position = (0 if row["direction"] == "increased" else 1) + (
            -0.18 if row["group"] == "Other feces" else 0.18
        )
        expected[round(position, 2)] = float(row["positive_z_rate"]) * 100
    assert {position for position, _height in bars} == set(expected)
    for position, height in bars:
        assert height == pytest.approx(expected[position], nan_ok=True)
    pd.testing.assert_frame_equal(comparisons, original)

    text = "\n".join(texts)
    assert "nan" not in text.lower()
    assert "1 directions" not in text
    if untested:
        direction = "increased" if case == "untested-increased" else "decreased"
        position = 0.18 if direction == "increased" else 1.18
        assert "1 direction tested" in text
        assert "not tested" in labels[position].lower()
        assert labels[position] != "100.0%"
        if direction == "decreased":
            assert "decreased comparison not tested" in text.lower()
            assert "nominal Fisher p=" not in text
    elif case == "tested-infinite-odds":
        assert "OR=inf" in text
        assert "nominal Fisher p=0.040" in text
        assert "not tested" not in text
    else:
        assert "OR not estimable" in text
        assert "nominal Fisher p=1.000" in text
        assert "not tested" not in text


def test_interpretations_preserve_tested_infinite_odds_ratio(tmp_path: Path) -> None:
    """A tested boundary odds estimate is different from an undefined odds ratio."""
    model, _groups, comparisons = _degenerate_consumer_tables(all_untested=False)
    decreased = comparisons["direction"].eq("decreased")
    comparisons.loc[decreased, "pd_vs_other_positive_z_odds_ratio"] = math.inf
    comparisons.loc[decreased, "pd_vs_other_fisher_p"] = 0.04
    assert comparisons.loc[decreased, "p_value_type"].eq("nominal").all()
    text = _write_synthetic_interpretations(tmp_path, model, comparisons)
    comparison_section = text.split("## PD vs Other Fecal Profiles", 1)[1].split(
        "The feces-only linear model", 1
    )[0]
    assert "an odds ratio of `inf`" in comparison_section
    assert "nominal p=0.04" in comparison_section
    assert "not estimable" not in comparison_section
    assert "not tested" not in comparison_section


def _write_synthetic_interpretations(
    directory: Path, model: pd.DataFrame, comparisons: pd.DataFrame
) -> str:
    """Exercise the real writer with small saved-table-shaped context."""
    display = pd.DataFrame(
        [
            {
                "direction": direction,
                "trait": cellulose.CELLULOSE_LABEL,
                "body_site": cellulose.FECES,
                "focus_enriched_profiles": 4,
                "focus_studies": 3,
                "focus_minus_comparator_z": 0.2,
            }
            for direction in ("increased", "decreased")
        ]
    )
    traits = pd.DataFrame(
        {"trait": [cellulose.CELLULOSE_LABEL], "increased_minus_decreased_mean_z": [-0.4]}
    )
    drivers = pd.DataFrame(
        {
            "taxon_name": ["Roseburia intestinalis"],
            "taxonomic_rank": ["species"],
            "decreased_profile_count": [4],
            "increased_profile_count": [2],
        }
    )
    output = directory / "interpretations.md"
    cellulose.write_interpretations(output, display, comparisons, traits, model, drivers)
    return output.read_text(encoding="utf-8")


@pytest.mark.parametrize("interaction_tested", [True, False], ids=["tested", "untested"])
def test_interpretations_use_singular_grammar_for_one_tablewide_nominal_term(
    tmp_path: Path, interaction_tested: bool
) -> None:
    """Writer availability fixtures cover a table-wide count the real cases do not."""
    model, _groups, comparisons = _degenerate_consumer_tables(all_untested=False)
    model.loc[:, "p_value_type"] = "not_tested"
    model.loc[:, "p_value"] = math.nan
    chosen = (
        model["outcome"].eq("cellulose_z") & model["term"].eq("pd_x_decreased")
        if interaction_tested
        else model["outcome"].eq("cellulose_positive_z") & model["term"].eq("intercept")
    )
    assert int(chosen.sum()) == 1
    model.loc[chosen, "p_value_type"] = "nominal"
    model.loc[chosen, ["p_value", "coefficient", "ci_low", "ci_high"]] = [0.04, 0.4, 0.1, 0.7]
    model.loc[:, "n_tests"] = 1
    assert int(model["p_value_type"].eq("nominal").sum()) == 1
    text = _write_synthetic_interpretations(tmp_path, model, comparisons)
    assert "1 terms across" not in text
    if interaction_tested:
        assert "one of 1" not in text
        assert "1 term across 2 outcomes" not in text
        assert "the only non-intercept model term tested in the figure" in text
        assert (
            "the only nominal term in `pd_cellulose_interaction_model.tsv`, which has 2 outcomes"
        ) in text
        assert "It tests whether" in text
    else:
        assert "1 of the 8 table terms has a nominal test across 2 outcomes" in text
        assert "1 of the 8 table terms have" not in text
        assert "It tests whether" not in text


def test_rank_deficient_and_cluster_degenerate_designs_mark_unestimable_terms_not_tested() -> None:
    """Rank-deficient and cluster-degenerate designs mark unestimable terms not_tested (#1515)."""
    import numpy as np

    # 1. Rank-deficient design X: only PD has decreased profiles, so pd_x_decreased == is_decreased
    collinear_rows = []
    collinear_z = []
    pd_inc_vals = (-0.5, 0.2, -0.3, 0.4, -0.2, 0.3)
    other_inc_vals = (-0.6, 0.3, -0.4, 0.2, -0.3, 0.4)
    for study_idx in range(6):
        for cond, z_inc, z_dec in [
            (cellulose.PD_CONDITION, pd_inc_vals[study_idx], 0.5 + 0.15 * study_idx),
            ("other condition", other_inc_vals[study_idx], None),
        ]:
            collinear_rows.append(
                {
                    "study_id": f"{cond}-{study_idx}",
                    "condition": cond,
                    "body_site_group": cellulose.FECES,
                    "inc_dec": "inc",
                }
            )
            collinear_z.append({cellulose.CELLULOSE_ID: z_inc, "METPO:1000668": -z_inc})
            if z_dec is not None:
                collinear_rows.append(
                    {
                        "study_id": f"{cond}-{study_idx}",
                        "condition": cond,
                        "body_site_group": cellulose.FECES,
                        "inc_dec": "dec",
                    }
                )
                collinear_z.append({cellulose.CELLULOSE_ID: z_dec, "METPO:1000668": -z_dec})

    collinear_meta = pd.DataFrame(collinear_rows)
    collinear_zscores = pd.DataFrame(collinear_z)
    rng = np.random.default_rng(151501)
    for _ in range(40):
        perm = rng.permutation(len(collinear_meta))
        model_df, _ = cellulose.build_cellulose_model_tables(
            collinear_meta.iloc[perm].reset_index(drop=True),
            collinear_zscores.iloc[perm].reset_index(drop=True),
        )
        assert model_df["n_tests"].eq(4).all()
        for outcome in ("cellulose_z", "cellulose_positive_z"):
            sub = model_df[model_df["outcome"].eq(outcome)].set_index("term")
            assert sub.loc[["intercept", "is_pd"], "p_value_type"].tolist() == [
                "nominal",
                "nominal",
            ]
            assert sub.loc[["intercept", "is_pd"], "p_value"].notna().all()
            assert sub.loc[["is_decreased", "pd_x_decreased"], "p_value_type"].tolist() == [
                "not_tested",
                "not_tested",
            ]
            assert (
                sub.loc[
                    ["is_decreased", "pd_x_decreased"],
                    ["coefficient", "std_error_clustered_by_study", "t_stat", "p_value"],
                ]
                .isna()
                .all()
                .all()
            )

    # 2. Cluster-degenerate design: within-study jitter is shared between inc and dec,
    # so is_decreased and pd_x_decreased in cellulose_z and is_decreased in cellulose_positive_z
    # have analytically zero cluster-robust SE and are marked not_tested under every row order.
    deg_rows = []
    deg_z = []
    for study_idx in range(6):
        jitter = (1.0 if study_idx % 2 == 0 else -1.0) + 0.05 * study_idx
        for cond, base in [(cellulose.PD_CONDITION, 0.8), ("other condition", -0.4)]:
            for inc_dec, delta in [("inc", -0.3), ("dec", 0.9)]:
                deg_rows.append(
                    {
                        "study_id": f"{cond}-{study_idx}",
                        "condition": cond,
                        "body_site_group": cellulose.FECES,
                        "inc_dec": inc_dec,
                    }
                )
                deg_z.append(
                    {
                        cellulose.CELLULOSE_ID: base + delta + jitter,
                        "METPO:1000668": -base,
                    }
                )
    deg_meta = pd.DataFrame(deg_rows)
    deg_zscores = pd.DataFrame(deg_z)
    expected_types = [
        ("cellulose_z", "intercept", "nominal"),
        ("cellulose_z", "is_pd", "nominal"),
        ("cellulose_z", "is_decreased", "not_tested"),
        ("cellulose_z", "pd_x_decreased", "not_tested"),
        ("cellulose_positive_z", "intercept", "nominal"),
        ("cellulose_positive_z", "is_pd", "nominal"),
        ("cellulose_positive_z", "is_decreased", "not_tested"),
        ("cellulose_positive_z", "pd_x_decreased", "nominal"),
    ]
    for scale_factor in (1.0, 1e6):
        scaled_z = deg_zscores.assign(
            **{cellulose.CELLULOSE_ID: deg_zscores[cellulose.CELLULOSE_ID] * scale_factor}
        )
        for _ in range(25):
            perm = rng.permutation(len(deg_meta))
            model_df, _ = cellulose.build_cellulose_model_tables(
                deg_meta.iloc[perm].reset_index(drop=True),
                scaled_z.iloc[perm].reset_index(drop=True),
            )
            assert model_df["n_tests"].eq(5).all()
            got = list(
                model_df[["outcome", "term", "p_value_type"]].itertuples(index=False, name=None)
            )
            assert got == expected_types
            not_tested_rows = model_df[model_df["p_value_type"].eq("not_tested")]
            assert not_tested_rows[["t_stat", "p_value"]].isna().all().all()


@pytest.mark.parametrize("profiles", [400, 20000])
def test_estimability_is_invariant_to_outcome_shift_and_holds_for_large_designs(
    profiles: int,
) -> None:
    """Shifting the outcome or growing the design never marks real slopes untested (#1538)."""
    import numpy as np

    rng = np.random.default_rng(1538 + profiles)
    is_pd = rng.integers(0, 2, profiles)
    is_decreased = rng.integers(0, 2, profiles)
    frame = pd.DataFrame(
        {
            "study_id": rng.integers(0, max(profiles // 6, 12), profiles),
            "is_pd": is_pd,
            "is_decreased": is_decreased,
            "pd_x_decreased": is_pd * is_decreased,
        }
    )
    frame["outcome"] = (
        0.3 * frame["is_pd"] - 0.2 * frame["is_decreased"] + rng.normal(0.0, 1.0, profiles)
    )
    frame["shifted"] = frame["outcome"] + 1000.0
    terms = ["is_pd", "is_decreased", "pd_x_decreased"]
    base = cellulose.fit_cluster_ols(frame, "outcome", terms, "study_id")
    shifted = cellulose.fit_cluster_ols(frame, "shifted", terms, "study_id")
    slopes = base.iloc[1:]
    assert slopes["p_value"].notna().all()
    assert shifted.iloc[1:]["p_value"].to_numpy() == pytest.approx(slopes["p_value"].to_numpy())
