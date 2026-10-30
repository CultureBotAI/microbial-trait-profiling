"""Saved study-weight sensitivity preserves full-family estimates and ties (#1046)."""

from __future__ import annotations

import re
from types import SimpleNamespace
from typing import TYPE_CHECKING, cast

import pandas as pd
import pytest
from matplotlib.text import Text

from mtp.publication_figures import rank_sensitivity
from mtp.publication_figures.plots import validate_authored_wording
from mtp.publication_figures.rank_sensitivity import (
    rank_sensitivity_pages,
    rank_sensitivity_table,
)

if TYPE_CHECKING:
    from mtp.publication_figures.inputs import FigureInputs
    from mtp.publication_figures.plots import Page


MEANS = ("experiment_mean_delta_z", "study_mean_delta_z")
RANKS = tuple(
    f"{mean}_rank_{direction}" for mean in MEANS for direction in ("ascending", "descending")
)


def source() -> FigureInputs:
    # Independently specified persisted values/ranks; no rank computation in the fixture.
    values = [
        ("A", 9.0, -1.0, 9, 1, 3, 7),
        ("B", 8.0, 8.0, 8, 2, 9, 1),
        ("C", 1.0, 7.0, 7, 3, 7, 2),
        ("D", 0.5, 7.0, 6, 4, 7, 2),
        ("E", -0.5, -8.0, 4, 6, 1, 9),
        ("F", -1.0, -7.0, 3, 7, 2, 8),
        ("G", -8.0, 0.0, 2, 8, 4, 5),
        ("H", -9.0, 1.0, 1, 9, 6, 4),
        ("I", 0.0, 0.0, 5, 5, 4, 5),
    ]
    frame = pd.DataFrame(values, columns=["trait", *MEANS, *RANKS])
    frame["trait_label"] = [f"Native label {trait}" for trait in frame.trait]
    return cast(
        "FigureInputs",
        SimpleNamespace(
            parent=SimpleNamespace(
                evidence=frame, record={"retained_traits": frame.trait.tolist()}
            ),
            record={"strain_rollup": "phenotype", "contrast_population": "fixture"},
            labels=dict(zip(frame.trait, frame.trait_label, strict=True)),
            tables={},
        ),
    )


def rendered(pages: list[Page]) -> str:
    return "\n".join(text.get_text() for page in pages for text in page.figure.findobj(Text))


def test_table_preserves_every_saved_value_and_unselected_trait() -> None:
    bundle = source()
    before = bundle.parent.evidence.copy(deep=True)
    result = rank_sensitivity_table(bundle, top_n=2)
    pd.testing.assert_frame_equal(
        result[["trait", "trait_label", *MEANS, *RANKS]],
        before[["trait", "trait_label", *MEANS, *RANKS]],
    )
    pd.testing.assert_frame_equal(bundle.parent.evidence, before)
    assert result.trait.tolist() == list("ABCDEFGHI")
    by_id = result.set_index("trait")
    assert set(result.loc[result.selected_by_experiment_top_n, "trait"]) == set("ABGH")
    assert set(result.loc[result.selected_by_study_top_n, "trait"]) == set("BCDEF")
    assert set(result.loc[result.study_only_selection, "trait"]) == set("CDEF")
    assert not bool(by_id.loc["I", "selected_for_display"])
    # A's sign reverses, but each selection uses its own saved mean/ranks.
    assert bool(by_id.loc["A", "selected_by_experiment_top_n"])
    assert not bool(by_id.loc["A", "selected_by_study_top_n"])


def test_renderer_uses_complete_tied_union_in_experiment_directional_order() -> None:
    pages = list(rank_sensitivity_pages(source(), top_n=2))
    assert [trait for page in pages for trait in page.row_ids] == list("ABCDHGFE")
    assert sum(page.mark_count for page in pages) == 16
    assert all(page.details["retained_family_count"] == 9 for page in pages)
    assert all(page.details["selected_trait_ids"] == list("ABCDHGFE") for page in pages)
    assert "Study only" in rendered(pages)
    assert "equal study weights" in rendered(pages)
    assert "ties retained" in rendered(pages)
    assert "Native label I" not in rendered(pages)


def test_saved_rank_is_not_recomputed_from_selected_rows(monkeypatch: pytest.MonkeyPatch) -> None:
    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("Renderer must not recompute ranks or means")

    monkeypatch.setattr(pd.Series, "rank", forbidden)
    monkeypatch.setattr(pd.DataFrame, "rank", forbidden)
    monkeypatch.setattr(pd.Series, "mean", forbidden)
    result = rank_sensitivity_table(source(), top_n=2)
    assert result.loc[result.trait.eq("C"), "study_mean_delta_z_rank_descending"].item() == 2
    assert list(rank_sensitivity_pages(source(), top_n=2))


@pytest.mark.parametrize("column", ["study_mean_delta_z", "experiment_mean_delta_z", *RANKS])
def test_missing_saved_column_refuses_without_ledger_fallback(column: str) -> None:
    bundle = source()
    bundle.parent.evidence.drop(columns=column, inplace=True)
    with pytest.raises(ValueError, match="Missing saved"):
        rank_sensitivity_table(bundle)


@pytest.mark.parametrize("top_n", [0, -1, True, 1.5])
def test_invalid_selection_cutoff_refuses(top_n: object) -> None:
    with pytest.raises(ValueError, match="positive integer"):
        rank_sensitivity_table(source(), top_n=cast("int", top_n))


def test_decimal_saved_values_are_never_rounded_in_audit() -> None:
    bundle = source()
    bundle.parent.evidence["study_mean_delta_z"] = bundle.parent.evidence.study_mean_delta_z.map(
        str
    )
    bundle.parent.evidence.loc[0, "study_mean_delta_z"] = "-1.0000000000000002"
    result = rank_sensitivity_table(bundle)
    assert result.loc[0, "study_mean_delta_z"] == "-1.0000000000000002"


def test_missing_saved_values_remain_unavailable_and_unselected() -> None:
    bundle = source()
    bundle.parent.evidence.loc[8, [*MEANS, *RANKS]] = float("nan")
    result = rank_sensitivity_table(bundle, top_n=2)
    assert result.loc[8, [*MEANS, *RANKS]].isna().all()
    assert not result.loc[8, "selected_for_display"]


def test_all_rendered_authored_text_is_guarded_and_source_label_is_exempt(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bundle = source()
    bundle.labels["A"] = "Native redox label"
    seen: list[str] = []
    original = validate_authored_wording

    def record(value: str) -> str:
        seen.append(value)
        return original(value)

    monkeypatch.setattr(rank_sensitivity, "validate_authored_wording", record)
    pages = list(rank_sensitivity_pages(bundle, top_n=2))
    text = rendered(pages)
    assert "Native redox label" in text
    assert all("Native redox label" not in value for value in seen)
    assert any("Secondary weighting sensitivity" in value for value in seen)
    assert any("equal study weights" in value for value in seen)
    # Independent lexical check after removing only the known literal source label.
    assert not re.search(
        r"\b(redox|coupled|enriched|depleted|fecal|faecal|abundance)\b",
        text.replace("Native redox label", ""),
        re.IGNORECASE,
    )
    assert all(
        float(text.get_fontsize()) >= 8
        for page in pages
        for text in page.figure.findobj(Text)
        if text.get_text()
    )
    assert all(page.figure.get_figwidth() == pytest.approx(7.2) for page in pages)


def test_large_tie_family_paginates_without_omitting_rows() -> None:
    bundle = source()
    frame = pd.concat([bundle.parent.evidence.iloc[[1]]] * 35, ignore_index=True)
    frame["trait"] = [f"trait-{number:03d}" for number in range(35)]
    frame["trait_label"] = [f"Verified trait {number}" for number in range(35)]
    cast("SimpleNamespace", bundle.parent).evidence = frame
    bundle.parent.record["retained_traits"] = frame.trait.tolist()
    bundle.labels.update(dict(zip(frame.trait, frame.trait_label, strict=True)))
    pages = list(rank_sensitivity_pages(bundle, top_n=2))
    assert len(pages) > 1
    assert [trait for page in pages for trait in page.row_ids] == frame.trait.tolist()
    assert sum(page.mark_count for page in pages) == 70
    assert all(page.figure.get_figheight() <= 10.5 for page in pages)


def test_zero_only_family_keeps_audit_and_renders_explicit_empty_state() -> None:
    bundle = source()
    bundle.parent.evidence.loc[:, list(MEANS)] = 0.0
    result = rank_sensitivity_table(bundle)
    assert len(result) == 9
    assert not result.selected_for_display.any()
    pages = list(rank_sensitivity_pages(bundle))
    assert len(pages) == 1
    assert pages[0].row_ids == []
    assert pages[0].mark_count == 0
    assert "No traits meet" in rendered(pages)
