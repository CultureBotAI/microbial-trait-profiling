"""Regression tests for #1033 Item 2: final rerun record citation in figure_methods.md."""

from __future__ import annotations

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURE_METHODS_MD = REPO_ROOT / "trait-visualization" / "figure_methods.md"
RERUN_RECORD_REL = "docs/reviews/2026-10-23-final-rerun-record.md"
RERUN_RECORD_LINK = "../docs/reviews/2026-10-23-final-rerun-record.md"
RERUN_RECORD_MD = REPO_ROOT / RERUN_RECORD_REL


def _clustered_heatmaps_intro(text: str) -> str:
    heading = "## Clustered Heatmaps"
    next_heading = "### `taxa_guild_discovery_clustered_heatmap.png`"
    assert heading in text, "Missing '## Clustered Heatmaps' section in figure_methods.md"
    after_heading = text.split(heading, 1)[1]
    assert next_heading in after_heading, (
        f"Missing subsection {next_heading!r} after {heading!r} in figure_methods.md"
    )
    intro = after_heading.split(next_heading, 1)[0].strip()
    assert intro, "Empty intro under '## Clustered Heatmaps' in figure_methods.md"
    return " ".join(intro.split())


def _assert_clustered_heatmaps_rerun_record_citation(
    figure_methods_text: str, rerun_record_text: str
) -> None:
    assert rerun_record_text.startswith("# Final rerun record: commit 511c891 (2026-10-23)"), (
        f"Unexpected header in {RERUN_RECORD_REL}"
    )
    resolved_link = (FIGURE_METHODS_MD.parent / RERUN_RECORD_LINK).resolve()
    assert resolved_link == RERUN_RECORD_MD.resolve(), (
        f"Relative link {RERUN_RECORD_LINK!r} does not resolve to {RERUN_RECORD_REL}"
    )
    assert resolved_link.is_file(), f"Missing canonical final rerun record: {RERUN_RECORD_REL}"

    flat_intro = _clustered_heatmaps_intro(figure_methods_text)
    assert "from the final rerun record before reusing them" not in flat_intro, (
        "Unnamed 'from the final rerun record before reusing them' must be replaced with "
        f"an explicit citation to {RERUN_RECORD_REL}"
    )
    assert "frozen code commit `511c891` (#1004)" in flat_intro, (
        "Clustered Heatmaps intro must preserve frozen code commit `511c891` (#1004)"
    )
    assert f"`{RERUN_RECORD_REL}`" in flat_intro, (
        f"Clustered Heatmaps intro must explicitly cite `{RERUN_RECORD_REL}`"
    )
    assert f"({RERUN_RECORD_LINK})" in flat_intro, (
        f"Clustered Heatmaps intro must include relative link ({RERUN_RECORD_LINK})"
    )
    assert "canonical final rerun record" in flat_intro, (
        "Clustered Heatmaps intro must identify the document as the canonical final rerun record"
    )
    assert (
        "reconciling the historical pre-fix tree and the corrected 19-table rerun" in flat_intro
    ), (
        "Clustered Heatmaps intro must state that the final rerun record reconciles "
        "the historical pre-fix tree and the corrected 19-table rerun"
    )
    assert "Read their population and strain-rollup scope from" in flat_intro, (
        "Clustered Heatmaps intro must preserve the population and strain-rollup scope instruction"
    )
    assert "before reusing them." in flat_intro, (
        "Clustered Heatmaps intro must preserve the 'before reusing them.' clause"
    )
    assert (
        "The shared-label implementation and new paired-metadata basis field added in #962 "
        "apply to subsequent runs." in flat_intro
    ), "Clustered Heatmaps intro must preserve the #962 subsequent-runs scope sentence"
    assert (
        "This follow-up does not regenerate those figures or establish their legibility "
        "at slide size." in flat_intro
    ), "Clustered Heatmaps intro must preserve the non-regeneration / slide-legibility caveat"


def test_figure_methods_clustered_heatmaps_cites_canonical_final_rerun_record() -> None:
    """#1033 Item 2: figure_methods.md:21 explicitly cites the canonical final rerun record."""
    _assert_clustered_heatmaps_rerun_record_citation(
        FIGURE_METHODS_MD.read_text(encoding="utf-8"),
        RERUN_RECORD_MD.read_text(encoding="utf-8"),
    )


@pytest.mark.parametrize(
    ("case_id", "target", "old", "new"),
    [
        pytest.param(
            "revert-unnamed-final-rerun-record",
            "figure_methods",
            (
                "from the canonical final rerun record in\n"
                "[`docs/reviews/2026-10-23-final-rerun-record.md`]"
                "(../docs/reviews/2026-10-23-final-rerun-record.md)\n"
                "(#1004, reconciling the historical pre-fix tree and the corrected 19-table rerun)"
            ),
            "from the final rerun record",
            id="revert-unnamed-final-rerun-record",
        ),
        pytest.param(
            "drop-canonical-record-path",
            "figure_methods",
            f"[`{RERUN_RECORD_REL}`]({RERUN_RECORD_LINK})",
            f"[final rerun record]({RERUN_RECORD_LINK})",
            id="drop-canonical-record-path",
        ),
        pytest.param(
            "drop-relative-link",
            "figure_methods",
            f"({RERUN_RECORD_LINK})",
            "(docs/reviews/2026-10-23-final-rerun-record.md)",
            id="drop-relative-link",
        ),
        pytest.param(
            "drop-canonical-label",
            "figure_methods",
            "canonical final rerun record",
            "final rerun note",
            id="drop-canonical-label",
        ),
        pytest.param(
            "drop-reconciliation-clause",
            "figure_methods",
            "reconciling the historical pre-fix tree and the corrected 19-table rerun",
            "recording the rerun",
            id="drop-reconciliation-clause",
        ),
        pytest.param(
            "drop-commit-511c891",
            "figure_methods",
            "frozen code commit `511c891` (#1004)",
            "frozen code commit (#1004)",
            id="drop-commit-511c891",
        ),
        pytest.param(
            "drop-962-sentence",
            "figure_methods",
            "added in #962",
            "added earlier",
            id="drop-962-sentence",
        ),
        pytest.param(
            "corrupt-rerun-record-header",
            "rerun_record",
            "# Final rerun record: commit 511c891 (2026-10-23)",
            "# Draft rerun notes",
            id="corrupt-rerun-record-header",
        ),
    ],
)
def test_figure_methods_1033_guard_rejects_unnamed_or_incomplete_citation(
    case_id: str, target: str, old: str, new: str
) -> None:
    """#1033 Item 2: guard rejects unnamed record, missing path/link, or missing reconciliation."""
    figure_methods_text = FIGURE_METHODS_MD.read_text(encoding="utf-8")
    rerun_record_text = RERUN_RECORD_MD.read_text(encoding="utf-8")
    if target == "figure_methods":
        assert figure_methods_text.count(old) == 1, f"Expected unique anchor for {case_id}"
        figure_methods_text = figure_methods_text.replace(old, new, 1)
    else:
        assert rerun_record_text.count(old) == 1, f"Expected unique anchor for {case_id}"
        rerun_record_text = rerun_record_text.replace(old, new, 1)

    with pytest.raises(AssertionError):
        _assert_clustered_heatmaps_rerun_record_citation(figure_methods_text, rerun_record_text)
