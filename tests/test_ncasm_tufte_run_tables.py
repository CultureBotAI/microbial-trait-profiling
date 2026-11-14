"""Reject internally contradictory NCASM study-screen rows (#2380)."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path

import pytest

from mtp.ncasm_tufte_plates import run_tables

FIXTURE_RUN = Path(__file__).parent / "fixtures/ncasm_tufte_plates/run_tier12_off_511c891"
PLOTTED_TRAIT = "METPO:2000103|GO:0030245"
SCREEN_ONLY_TRAIT = "METPO:1000603"


def _rows(root: Path, relative: Path) -> list[dict[str, str]]:
    with (root / relative).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _edited_screen(tmp_path: Path, trait: str, changes: dict[str, str]) -> Path:
    """Change exactly one real fixture row, including rows not drawn in the cloud."""
    root = tmp_path / "run"
    shutil.copytree(FIXTURE_RUN, root)
    rows = _rows(root, run_tables.DIRECTION_SCREEN)
    matches = [row for row in rows if row["trait"] == trait]
    assert len(matches) == 1
    plotted = {row["trait"] for row in _rows(root, run_tables.ARM1_POOLED_STATS)}
    assert (trait in plotted) == (trait == PLOTTED_TRAIT)
    matches[0].update(changes)
    with (root / run_tables.DIRECTION_SCREEN).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return root


def test_committed_screen_has_producer_rounded_shares_and_excludes_ties() -> None:
    cloud = run_tables.load_trait_cloud(FIXTURE_RUN)
    point = next(point for point in cloud.points if point.trait == PLOTTED_TRAIT)
    row = next(
        row
        for row in _rows(FIXTURE_RUN, run_tables.DIRECTION_SCREEN)
        if row["trait"] == PLOTTED_TRAIT
    )
    assert (point.studies_up, point.studies_down, int(row["studies_scored"])) == (181, 265, 446)
    assert int(row["studies_tied"]) == 94
    assert point.share_up == 0.40583
    assert point.share_up != 181 / 446
    assert point.share_up != pytest.approx(181 / (446 + 94))
    assert (len(cloud.points), cloud.screen_traits) == (311, 317)


@pytest.mark.parametrize(
    "trait", [PLOTTED_TRAIT, SCREEN_ONLY_TRAIT], ids=["plotted", "screen-only"]
)
@pytest.mark.parametrize("field", ["share_up", "by_q_nominal"])
@pytest.mark.parametrize("value", ["-0.01", "1.01", "nan", "inf", "-inf"])
def test_loader_refuses_invalid_screen_probabilities(
    tmp_path: Path, trait: str, field: str, value: str
) -> None:
    root = _edited_screen(tmp_path, trait, {field: value})
    with pytest.raises(run_tables.RunTableError, match=field):
        run_tables.load_trait_cloud(root)


@pytest.mark.parametrize(
    "trait", [PLOTTED_TRAIT, SCREEN_ONLY_TRAIT], ids=["plotted", "screen-only"]
)
@pytest.mark.parametrize("field", ["studies_up", "studies_down", "studies_scored"])
@pytest.mark.parametrize("value", ["-1", "1.5"])
def test_loader_refuses_invalid_counts_even_for_screen_only_traits(
    tmp_path: Path, trait: str, field: str, value: str
) -> None:
    root = _edited_screen(tmp_path, trait, {field: value})
    with pytest.raises(run_tables.RunTableError, match=field):
        run_tables.load_trait_cloud(root)


@pytest.mark.parametrize(
    "trait", [PLOTTED_TRAIT, SCREEN_ONLY_TRAIT], ids=["plotted", "screen-only"]
)
@pytest.mark.parametrize(
    ("changes", "message"),
    [
        pytest.param({"studies_up": "9999"}, "studies_scored", id="up-disagrees"),
        pytest.param({"studies_down": "9999"}, "studies_scored", id="down-disagrees"),
        pytest.param({"studies_scored": "9999"}, "studies_scored", id="scored-disagrees"),
        pytest.param(
            {"studies_up": "0", "studies_down": "0", "studies_scored": "0", "share_up": "0"},
            "studies_scored",
            id="tested-with-no-votes",
        ),
        pytest.param(
            {"studies_up": "1", "studies_down": "2", "studies_scored": "3", "share_up": "0.6"},
            "share_up",
            id="share-disagrees",
        ),
        pytest.param(
            {
                "studies_up": "1",
                "studies_down": "2",
                "studies_tied": "7",
                "studies_scored": "10",
                "share_up": "0.1",
            },
            "studies_scored",
            id="ties-counted-as-scored",
        ),
    ],
)
def test_loader_refuses_contradictory_tested_votes(
    tmp_path: Path, trait: str, changes: dict[str, str], message: str
) -> None:
    root = _edited_screen(tmp_path, trait, changes)
    with pytest.raises(run_tables.RunTableError, match=message):
        run_tables.load_trait_cloud(root)


@pytest.mark.parametrize(
    "trait", [PLOTTED_TRAIT, SCREEN_ONLY_TRAIT], ids=["plotted", "screen-only"]
)
@pytest.mark.parametrize(
    ("up", "down", "tied", "share", "q"),
    [
        pytest.param(0, 5, 7, "0", "0", id="zero-probabilities"),
        pytest.param(5, 0, 7, "1", "1", id="one-probabilities"),
        pytest.param(1, 2, 7, "0.333333", "0.25", id="rounded-down"),
        pytest.param(2, 1, 7, "0.666667", "0.25", id="rounded-up"),
        pytest.param(181, 265, 94, "0.40583", "0.0407707", id="fixture-rounding"),
        pytest.param(1, 538, 1, "0.00185529", "0.75", id="small-six-significant-digits"),
    ],
)
def test_loader_accepts_boundary_probabilities_and_serialized_vote_shares(
    tmp_path: Path, trait: str, up: int, down: int, tied: int, share: str, q: str
) -> None:
    root = _edited_screen(
        tmp_path,
        trait,
        {
            "studies_up": str(up),
            "studies_down": str(down),
            "studies_tied": str(tied),
            "studies_scored": str(up + down),
            "share_up": share,
            "by_q_nominal": q,
        },
    )
    cloud = run_tables.load_trait_cloud(root)
    assert (len(cloud.points), cloud.screen_traits) == (311, 317)
    if trait == PLOTTED_TRAIT:
        point = next(point for point in cloud.points if point.trait == trait)
        assert (point.studies_up, point.studies_down, point.share_up, point.by_q) == (
            up,
            down,
            float(share),
            float(q),
        )
    else:
        assert trait not in {point.trait for point in cloud.points}
