import math
from itertools import combinations

import matplotlib.pyplot as plt
import pytest
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from mtp.body_site_diagram import (
    MAX_BUBBLE_RADIUS,
    MIN_BUBBLE_AREA_FRACTION,
    OTHER_ANCHOR,
    SITE_ANCHORS,
    UNKNOWN_ANCHOR,
    add_colorbar_beside,
    anatomical_site_base,
    bubble_radius,
    leader_line,
)

# The 14 largest body-site labels of study_trait_profile_study_metadata.tsv with
# per-experiment labels (#555), with their profile counts. The default top-8 grouping
# draws the first eight; a larger --max-body-site-categories adds the rest (#620).
SITE_PROFILE_COUNTS = {
    "Feces": 3604,
    "Saliva": 374,
    "Vagina": 133,
    "Uterine cervix": 81,
    "Mouth": 80,
    "Subgingival dental plaque": 77,
    "Oral cavity": 75,
    "Throat": 66,
    "Nasopharynx": 60,
    "Rectum": 46,
    "Stomach": 44,
    "Skin of body": 44,
    "Intestine": 41,
    "Lung": 37,
}
# A tenth of the largest group's count: every common site but Feces stays below it.
MINOR_SHARE = 0.1
DEFAULT_TOP_GROUPS = list(SITE_PROFILE_COUNTS)[:8]
# Label column (top and bottom label y) and leader-line elbow x of the v2 site
# diagrams, the anaerobe enrichment-filtered diagram and the infographic's panel C.
LABEL_LAYOUTS = {
    "v2 sites": (0.84, 0.14, 0.52),
    "v2 enriched profiles": (0.84, 0.24, 0.52),
    "infographic": (0.90, 0.10, 0.55),
}


def distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.dist(a, b)


def segment_distance(
    point: tuple[float, float], start: tuple[float, float], end: tuple[float, float]
) -> float:
    dx, dy = end[0] - start[0], end[1] - start[1]
    along = ((point[0] - start[0]) * dx + (point[1] - start[1]) * dy) / (dx * dx + dy * dy)
    along = min(max(along, 0.0), 1.0)
    return distance(point, (start[0] + along * dx, start[1] + along * dy))


def test_mouth_and_oral_cavity_have_separate_anchors() -> None:
    mouth = anatomical_site_base("Mouth")
    oral_cavity = anatomical_site_base("Oral cavity")
    assert mouth != oral_cavity
    assert anatomical_site_base("Mouth mucosa") == mouth
    assert anatomical_site_base("Mucosa of oral region") == oral_cavity


def test_throat_and_pharynx_sit_on_the_neck() -> None:
    neck = anatomical_site_base("Throat")
    assert anatomical_site_base("Oropharynx") == neck
    assert anatomical_site_base("Nasopharynx") != neck
    assert neck != OTHER_ANCHOR


def test_cervix_and_uterine_sites_sit_in_the_pelvis() -> None:
    pelvis = anatomical_site_base("Uterine cervix")
    assert anatomical_site_base("Uterus") == pelvis
    assert pelvis != anatomical_site_base("Vagina")
    assert pelvis != OTHER_ANCHOR


def test_every_common_site_has_its_own_anchor() -> None:
    anchors = {site: anatomical_site_base(site) for site in SITE_PROFILE_COUNTS}
    assert len(set(anchors.values())) == len(anchors)
    assert not {OTHER_ANCHOR, UNKNOWN_ANCHOR} & set(anchors.values())


def test_gut_lung_and_skin_sites_follow_the_body() -> None:
    x, y = 0, 1
    throat = anatomical_site_base("Throat")
    lung = anatomical_site_base("Lung")
    stomach = anatomical_site_base("Stomach")
    intestine = anatomical_site_base("Intestine")
    rectum = anatomical_site_base("Rectum")
    feces = anatomical_site_base("Feces")
    # Lung in the chest, then the gut down the abdomen to the rectum, above feces.
    assert throat[y] > lung[y] > stomach[y] > intestine[y] > rectum[y] > feces[y]
    # The outline faces the viewer: the subject's left (stomach, colon) has x > 0.30
    # and the caecum, on the subject's right, has x < 0.30.
    assert stomach[x] > 0.30
    assert anatomical_site_base("Colon")[x] > 0.30
    assert anatomical_site_base("Caecum")[x] < 0.30
    assert rectum[x] == pytest.approx(0.30)
    # Skin sits on an arm, off the torso, which spans x = 0.18-0.42 at that height.
    assert abs(anatomical_site_base("Skin of body")[x] - 0.30) > 0.12
    assert anatomical_site_base("Gastric juice") == stomach
    assert anatomical_site_base("Small intestine") == intestine
    assert anatomical_site_base("Colorectum") == anatomical_site_base("Colon")


def test_bubbles_of_the_common_sites_do_not_overlap() -> None:
    largest = max(SITE_PROFILE_COUNTS.values())
    bubbles = [
        (site, anatomical_site_base(site), bubble_radius(count, largest))
        for site, count in SITE_PROFILE_COUNTS.items()
    ]
    overlaps = [
        (site_a, site_b)
        for (site_a, anchor_a, radius_a), (site_b, anchor_b, radius_b) in combinations(bubbles, 2)
        if distance(anchor_a, anchor_b) < radius_a + radius_b
    ]
    assert overlaps == []


def test_anchors_leave_room_for_minor_sites_and_the_feces_bubble() -> None:
    minor = bubble_radius(MINOR_SHARE, 1.0)
    anchors = [anchor for _, anchor in SITE_ANCHORS] + [OTHER_ANCHOR, UNKNOWN_ANCHOR]
    crowded = [(a, b) for a, b in combinations(anchors, 2) if distance(a, b) < 2 * minor]
    assert crowded == []
    feces = anatomical_site_base("Feces")
    near_feces = [
        anchor
        for anchor in anchors
        if anchor != feces and distance(anchor, feces) < MAX_BUBBLE_RADIUS + minor
    ]
    assert near_feces == []


def test_leader_line_starts_on_the_bubble_edge() -> None:
    xs, ys = leader_line((0.30, 0.80), 0.02, 0.50, 0.52, 0.57)
    start = (xs[0], ys[0])
    assert distance(start, (0.30, 0.80)) == pytest.approx(0.02)
    assert segment_distance(start, (0.30, 0.80), (0.52, 0.50)) == pytest.approx(0.0)
    assert list(zip(xs[1:], ys[1:], strict=True)) == [(0.52, 0.50), (0.57, 0.50)]


@pytest.mark.parametrize("layout", LABEL_LAYOUTS)
@pytest.mark.parametrize(
    "groups",
    [DEFAULT_TOP_GROUPS, [*DEFAULT_TOP_GROUPS[:7], "Nasopharynx"]],
    ids=["top 8", "nasopharynx for throat"],
)
def test_leader_lines_clear_the_other_bubbles(layout: str, groups: list[str]) -> None:
    """Labels follow the anchors top to bottom, as the diagrams sort them (#612)."""
    top, bottom, elbow = LABEL_LAYOUTS[layout]
    largest = max(SITE_PROFILE_COUNTS[group] for group in groups)
    sites = sorted(
        groups,
        key=lambda group: (-anatomical_site_base(group)[1], -SITE_PROFILE_COUNTS[group]),
    )
    step = (top - bottom) / (len(sites) - 1)
    crossed = [
        (site, other)
        for index, site in enumerate(sites)
        for other in sites
        if other != site
        and segment_distance(
            anatomical_site_base(other), anatomical_site_base(site), (elbow, top - index * step)
        )
        < bubble_radius(SITE_PROFILE_COUNTS[other], largest)
    ]
    assert crossed == []


def test_bubble_radius_grows_with_the_count_share() -> None:
    smallest = MAX_BUBBLE_RADIUS * math.sqrt(MIN_BUBBLE_AREA_FRACTION)
    assert bubble_radius(0, 100) == pytest.approx(smallest)
    assert bubble_radius(100, 100) == pytest.approx(MAX_BUBBLE_RADIUS)
    assert bubble_radius(10, 100) < bubble_radius(40, 100) < bubble_radius(90, 100)
    assert bubble_radius(150, 100) == pytest.approx(MAX_BUBBLE_RADIUS)
    assert bubble_radius(5, 0) == pytest.approx(smallest)


def test_other_and_unknown_stay_off_the_body() -> None:
    assert anatomical_site_base("Other") == OTHER_ANCHOR
    for value in ["unknown", None, " NA ", "nan"]:
        assert anatomical_site_base(value) == UNKNOWN_ANCHOR


@pytest.mark.parametrize("group", ["Thyroid gland", "Urine"])
def test_unmapped_body_site_group_raises(group: str) -> None:
    with pytest.raises(ValueError, match=repr(group)):
        anatomical_site_base(group)


def test_colorbar_sits_right_of_the_labels() -> None:
    fig, ax = plt.subplots(figsize=(11, 7.5))
    try:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal", adjustable="box")
        ax.axis("off")
        texts = [
            ax.text(0.58, 0.8, "Subgingival dental plaque", fontsize=10, fontweight="bold"),
            ax.text(0.58, 0.7, "3,604 profiles; increased-profile trait z +0.39", fontsize=8.5),
        ]
        # A wide title makes tight_layout hand the axes a box wider than it is tall,
        # which the equal aspect narrows only at draw time (#612).
        title = "Aerobe v2 Trait Enrichment in Increased Profiles: Anatomical Sites"
        ax.text(0.02, 1.005, title, fontsize=17, fontweight="bold", va="bottom")
        fig.tight_layout()
        mappable = ScalarMappable(norm=Normalize(-0.5, 0.5), cmap="RdBu_r")
        colorbar = add_colorbar_beside(fig, ax, mappable, texts, "trait z")
        fig.canvas.draw()
        colorbar_left = colorbar.ax.get_window_extent().x0
        assert all(text.get_window_extent().x1 < colorbar_left for text in texts)
        assert colorbar.ax.get_ylabel() == "trait z"
    finally:
        plt.close(fig)
