"""Display persisted weighting sensitivity without recalculating ranks (#1046).

The full-family audit preserves saved comparison and equal-study means/ranks.
Pages show the complete union of both signed selections, including rank ties.
This secondary display does not alter the comparison-ranked headline forest.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import TYPE_CHECKING, Any

import matplotlib
import pandas as pd
from matplotlib.figure import Figure

from mtp.publication_figures.plots import (
    ACCENT,
    FONT,
    INK,
    STYLE,
    Page,
    _context,
    _figure,
    _format,
    _label,
    _number,
    _wrap,
    validate_authored_wording,
)
from mtp.trait_key_policy import parse_trait_key_policy

if TYPE_CHECKING:
    from mtp.publication_figures.inputs import FigureInputs

MEANS = ("experiment_mean_delta_z", "study_mean_delta_z")
RANKS = tuple(
    f"{mean}_rank_{direction}" for mean in MEANS for direction in ("ascending", "descending")
)
MAX_BODY_HEIGHT = 7.8


def _cutoff(top_n: int) -> None:
    if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n <= 0:
        raise ValueError("top_n must be a positive integer")


def _selected(row: dict[str, Any], mean: str, top_n: int) -> bool:
    effect = _number(row[mean])
    ranks = [_number(row[f"{mean}_rank_{direction}"]) for direction in ("ascending", "descending")]
    if any(rank is not None and (rank < 1 or not rank.is_integer()) for rank in ranks):
        raise ValueError("Saved ranks must be positive integers or missing")
    if effect is None:
        if any(rank is not None for rank in ranks):
            raise ValueError("Missing saved mean has a defined rank")
        return False
    if any(rank is None for rank in ranks):
        raise ValueError("Defined saved mean lacks a rank")
    rank = ranks[1 if effect > 0 else 0]
    return effect != 0 and rank is not None and rank <= top_n


def rank_sensitivity_table(bundle: FigureInputs, *, top_n: int = 5) -> pd.DataFrame:
    """Preserve the complete saved family; append only explicit selection flags."""
    _cutoff(top_n)
    required = ["trait", "trait_label", *MEANS, *RANKS]
    evidence = bundle.parent.evidence
    missing = set(required) - set(evidence.columns)
    if missing:
        raise ValueError(f"Missing saved rank-sensitivity columns: {sorted(missing)}")
    if evidence.trait.duplicated().any() or any(
        not isinstance(trait, str) or not trait.strip() for trait in evidence.trait
    ):
        raise ValueError("Rank-sensitivity trait identities must be unique and nonblank")
    result = evidence[required].copy(deep=True)
    labels = [bundle.labels.get(trait) for trait in evidence.trait]
    if any(not isinstance(label, str) or not label.strip() for label in labels):
        raise ValueError("Rank-sensitivity trait lacks a verified label")
    result["trait_label"] = labels
    rows = [{str(key): value for key, value in row.items()} for row in result.to_dict("records")]
    result["selected_by_experiment_top_n"] = [_selected(row, MEANS[0], top_n) for row in rows]
    result["selected_by_study_top_n"] = [_selected(row, MEANS[1], top_n) for row in rows]
    result["selected_for_display"] = (
        result.selected_by_experiment_top_n | result.selected_by_study_top_n
    )
    result["study_only_selection"] = (
        result.selected_by_study_top_n & ~result.selected_by_experiment_top_n
    )
    return result


def _order(row: dict[str, Any]) -> tuple[int, float, str]:
    effect = _number(row[MEANS[0]])
    if effect is not None and effect != 0:
        suffix = "descending" if effect > 0 else "ascending"
        rank = _number(row[f"{MEANS[0]}_rank_{suffix}"])
        return (
            0 if effect > 0 else 1,
            rank if rank is not None else float("inf"),
            str(row["trait"]),
        )
    # A study-selected trait with zero/unavailable comparison mean remains visible.
    return (2, float("inf"), str(row["trait"]))


def _tail_rank(row: dict[str, Any], mean: str) -> str:
    value = _number(row[mean])
    if value is None:
        return "unavailable"
    if value == 0:
        return "zero\nnot selected"
    direction = "increase" if value > 0 else "decrease"
    suffix = "descending" if value > 0 else "ascending"
    return f"{direction}\n{_format(row[f'{mean}_rank_{suffix}'])}"


def _membership(row: dict[str, Any]) -> str:
    if row["study_only_selection"]:
        return "Study only"
    return "Both" if row["selected_by_study_top_n"] else "Comparison\nonly"


def _authored(figure: Figure, x: float, y: float, value: str, **kwargs: Any) -> None:
    figure.text(x, y, validate_authored_wording(value), fontsize=FONT, **kwargs)


def rank_sensitivity_pages(
    bundle: FigureInputs, *, trait_key_policy: str, top_n: int = 5
) -> Iterator[Page]:
    """Display saved selections and the policy supplied from the verified core."""
    policy = parse_trait_key_policy(trait_key_policy)
    audit = rank_sensitivity_table(bundle, top_n=top_n)
    rows = [
        {str(key): value for key, value in row.items()}
        for row in audit.loc[audit.selected_for_display].to_dict("records")
    ]
    rows.sort(key=_order)
    selected_ids = [str(row["trait"]) for row in rows]
    prepared = [(row, _label(bundle, str(row["trait"]), width=32)) for row in rows]
    groups: list[list[tuple[dict[str, Any], str, float]]] = [[]]
    used = 0.0
    for row, label in prepared:
        row_height = max(0.48, len(label.splitlines()) * 0.145 + 0.12)
        if row_height > MAX_BODY_HEIGHT:
            raise ValueError(
                "Verified trait label exceeds a complete page; no label truncation is permitted"
            )
        if groups[-1] and used + row_height > MAX_BODY_HEIGHT:
            groups.append([])
            used = 0.0
        groups[-1].append((row, label, row_height))
        used += row_height
    for page_number, group in enumerate(groups, start=1):
        body_height = sum(height for _, _, height in group)
        height = max(3.5, body_height + 2.35)
        with matplotlib.rc_context(STYLE):
            figure = _figure(height)

            _authored(figure, 0.04, 1 - 0.20 / height, "Study-weighted rank sensitivity", va="top")
            _authored(
                figure,
                0.04,
                1 - 0.48 / height,
                _wrap(
                    f"Top {top_n} saved ranks per direction; union of both weightings. "
                    + _context(bundle),
                    105,
                )
                + f"\ntrait keys: {policy}",
                va="top",
            )
            columns = (0.04, 0.45, 0.56, 0.67, 0.78, 0.91)
            headings = (
                "Trait annotation",
                "Comparison\nmean Δz",
                "Study\nmean Δz",
                "Comparison\nrank (tail)",
                "Study\nrank (tail)",
                "Selected by",
            )
            for x, heading in zip(columns, headings, strict=True):
                _authored(
                    figure,
                    x,
                    1 - 0.99 / height,
                    heading,
                    ha="left" if x == columns[0] else "center",
                    va="top",
                )
            offset = 1.48
            for row, label, row_height in group:
                y = 1 - (offset + row_height / 2) / height
                # Verified metadata labels and IDs are literal source text.
                figure.text(columns[0], y, label, fontsize=FONT, color=INK, va="center")
                cells = (
                    _format(row[MEANS[0]]),
                    _format(row[MEANS[1]]),
                    _tail_rank(row, MEANS[0]),
                    _tail_rank(row, MEANS[1]),
                    _membership(row),
                )
                for x, value in zip(columns[1:], cells, strict=True):
                    _authored(
                        figure,
                        x,
                        y,
                        value,
                        ha="center",
                        va="center",
                        color=ACCENT if row["study_only_selection"] else INK,
                    )
                offset += row_height
            if not group:
                _authored(
                    figure,
                    0.04,
                    0.5,
                    "No traits meet the declared rank and direction rule",
                    va="center",
                )
            _authored(
                figure,
                0.04,
                0.48 / height,
                "Secondary sensitivity: equal study weights after averaging comparisons "
                "within each study.\nSaved full-family ranks; ties retained. "
                "Headline forest remains comparison-ranked.",
                va="bottom",
            )
            _authored(
                figure,
                0.96,
                0.14 / height,
                f"Page {page_number} of {len(groups)}",
                ha="right",
                va="bottom",
            )
            marks = sum(_number(row[mean]) is not None for row, _, _ in group for mean in MEANS)
            page = Page(
                figure,
                [str(row["trait"]) for row, _, _ in group],
                marks,
                {
                    "trait_key_policy": policy,
                    "top_n": top_n,
                    "retained_family_count": len(audit),
                    "selected_trait_ids": selected_ids,
                    "selection_rule": (
                        "union of saved signed experiment/study ranks <= top_n; "
                        "ties retained; zero excluded"
                    ),
                    "ordering": "saved experiment directional ranks; trait ID breaks display ties",
                    "headline_ranking": "experiment_mean_delta_z",
                    "sensitivity_weighting": "equal_study_weights",
                    "study_only_trait_ids": [
                        str(row["trait"]) for row in rows if row["study_only_selection"]
                    ],
                    "point_estimands": list(MEANS),
                    "inference_recomputed": False,
                },
            )
        yield page
