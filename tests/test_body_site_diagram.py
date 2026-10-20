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
    leader_order_key,
)

Point = tuple[float, float]

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
CORE_GROUPS = list(SITE_PROFILE_COUNTS)[:8]
NEXT_GROUPS = list(SITE_PROFILE_COUNTS)[8:]
# Label column (top and bottom label y), leader-line elbow x and the most groups drawn
# by the v2 site diagrams (12), the anaerobe enrichment-filtered diagram (12) and the
# infographic's panel C (9).
LABEL_LAYOUTS = {
    "v2 sites": (0.84, 0.14, 0.52, 12),
    "v2 enriched profiles": (0.84, 0.24, 0.52, 12),
    "infographic": (0.90, 0.10, 0.55, 9),
}
GROUP_SETS = [
    (layout, [*CORE_GROUPS, *extra])
    for layout, (_, _, _, most) in LABEL_LAYOUTS.items()
    for size in range(most - len(CORE_GROUPS) + 1)
    for extra in combinations(NEXT_GROUPS, size)
]


def distance(a: Point, b: Point) -> float:
    return math.dist(a, b)


def segment_distance(point: Point, start: Point, end: Point) -> float:
    dx, dy = end[0] - start[0], end[1] - start[1]
    along = ((point[0] - start[0]) * dx + (point[1] - start[1]) * dy) / (dx * dx + dy * dy)
    along = min(max(along, 0.0), 1.0)
    return distance(point, (start[0] + along * dx, start[1] + along * dy))


def turn(origin: Point, a: Point, b: Point) -> float:
    return (a[0] - origin[0]) * (b[1] - origin[1]) - (a[1] - origin[1]) * (b[0] - origin[0])


def segments_cross(first: tuple[Point, Point], second: tuple[Point, Point]) -> bool:
    (p1, p2), (q1, q2) = first, second
    return turn(q1, q2, p1) * turn(q1, q2, p2) < 0 and turn(p1, p2, q1) * turn(p1, p2, q2) < 0


def leader_segments(groups: list[str], layout: str) -> dict[str, tuple[Point, Point]]:
    """Each group's leader from its bubble edge to the elbow, in the diagrams' order."""
    top, bottom, elbow, _ = LABEL_LAYOUTS[layout]
    largest = max(SITE_PROFILE_COUNTS[group] for group in groups)
    order = sorted(
        groups,
        key=lambda group: (
            -leader_order_key(anatomical_site_base(group), elbow),
            -SITE_PROFILE_COUNTS[group],
        ),
    )
    step = (top - bottom) / (len(order) - 1)
    segments = {}
    for index, group in enumerate(order):
        radius = bubble_radius(SITE_PROFILE_COUNTS[group], largest)
        xs, ys = leader_line(anatomical_site_base(group), radius, top - index * step, elbow, 1.0)
        segments[group] = ((xs[0], ys[0]), (xs[1], ys[1]))
    return segments


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
    y = 1
    throat = anatomical_site_base("Throat")
    lung = anatomical_site_base("Lung")
    stomach = anatomical_site_base("Stomach")
    intestine = anatomical_site_base("Intestine")
    rectum = anatomical_site_base("Rectum")
    feces = anatomical_site_base("Feces")
    # Lung in the chest, then the gut down the abdomen to the rectum, above feces.
    assert throat[y] > lung[y] > stomach[y] > intestine[y] > rectum[y] > feces[y]
    assert anatomical_site_base("Gastric juice") == stomach
    assert anatomical_site_base("Small intestine") == intestine


@pytest.mark.parametrize(
    ("site", "subject_side"),
    [
        ("Saliva", "left"),
        ("Lung", "left"),
        ("Stomach", "left"),
        ("Colon", "left"),
        ("Skin of body", "left"),
        ("Subgingival dental plaque", "right"),
        ("Caecum", "right"),
    ],
)
def test_lateral_anchors_sit_on_their_side(site: str, subject_side: str) -> None:
    """The outline faces the viewer, so the subject's left is at x > 0.30 (#712)."""
    x = anatomical_site_base(site)[0]
    assert x > 0.30 if subject_side == "left" else x < 0.30


def test_skin_sits_on_the_forearm_and_other_off_the_body() -> None:
    # The subject's left forearm spans about x = 0.41-0.46 at this height; the right
    # arm reaches out to x = 0.105, so Other at x = 0.07 is clear of it (#714).
    skin_x, skin_y = anatomical_site_base("Skin of body")
    assert 0.41 < skin_x < 0.46
    assert 0.40 < skin_y < 0.60
    assert OTHER_ANCHOR[0] < 0.10


@pytest.mark.parametrize(
    ("label", "site"),
    [
        ("Nasopharynx", "Nasopharynx"),
        ("Colorectal mucosa", "Colon"),
        ("Colorectum", "Colon"),
        ("Large intestine", "Colon"),
        ("Mucosa of small intestine", "Intestine"),
        ("Oropharynx", "Throat"),
    ],
)
def test_rule_order_settles_overlapping_keywords(label: str, site: str) -> None:
    """Pins the rule order: e.g. "nasopharynx" before "pharynx", colon before rectum (#712)."""
    assert anatomical_site_base(label) == anatomical_site_base(site)


@pytest.mark.parametrize(
    ("label", "site"),
    [
        ("Vagina,Uterine cervix", "Vagina"),
        ("Feces,Oral cavity", "Feces"),
        ("Feces,Caecum", "Feces"),
        ("Bronchus,Mouth", "Lung"),
        ("Ileum,Right colon,Rectum", "Intestine"),
    ],
)
def test_multi_site_labels_take_their_first_listed_site(label: str, site: str) -> None:
    assert anatomical_site_base(label) == anatomical_site_base(site)


def test_multi_site_label_with_unplaced_first_site_raises() -> None:
    with pytest.raises(ValueError, match="'meconium'"):
        anatomical_site_base("Meconium,Feces")


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


def test_leader_order_puts_throat_before_dental_plaque() -> None:
    """Ordering by height alone crossed the dental-plaque and throat leaders (#711)."""
    assert leader_order_key((0.30, 0.80), 0.52) == pytest.approx(0.58)
    throat = leader_order_key(anatomical_site_base("Throat"), 0.52)
    dental = leader_order_key(anatomical_site_base("Subgingival dental plaque"), 0.52)
    assert throat > dental


@pytest.mark.parametrize(
    ("layout", "groups"),
    GROUP_SETS,
    ids=[f"{layout}: {', '.join(groups[8:]) or 'core 8'}" for layout, groups in GROUP_SETS],
)
def test_leader_lines_clear_each_other_and_the_bubbles(layout: str, groups: list[str]) -> None:
    """The current 8 groups plus any of the next six, up to the most each figure draws."""
    segments = leader_segments(groups, layout)
    largest = max(SITE_PROFILE_COUNTS[group] for group in groups)
    # Keep at least 0.002 outline units (~1.5 pt, > the 0.85 pt combined half-strokes)
    # between each leader segment and every other bubble's edge (#718).
    through_bubbles = [
        (group, other)
        for group, (start, end) in segments.items()
        for other in groups
        if other != group
        and segment_distance(anatomical_site_base(other), start, end)
        < bubble_radius(SITE_PROFILE_COUNTS[other], largest) + 0.002
    ]
    assert through_bubbles == []
    crossing = [
        (first, second)
        for first, second in combinations(segments, 2)
        if segments_cross(segments[first], segments[second])
    ]
    assert crossing == []


def test_bubble_radius_grows_with_the_count_share() -> None:
    smallest = MAX_BUBBLE_RADIUS * math.sqrt(MIN_BUBBLE_AREA_FRACTION)
    assert bubble_radius(0, 100) == pytest.approx(smallest)
    assert bubble_radius(100, 100) == pytest.approx(MAX_BUBBLE_RADIUS)
    # The area grows with the square root of the share, not linearly (#712).
    quarter_area = MIN_BUBBLE_AREA_FRACTION + (1 - MIN_BUBBLE_AREA_FRACTION) * 0.5
    assert bubble_radius(25, 100) == pytest.approx(MAX_BUBBLE_RADIUS * math.sqrt(quarter_area))
    assert bubble_radius(10, 100) < bubble_radius(40, 100) < bubble_radius(90, 100)
    assert bubble_radius(150, 100) == pytest.approx(MAX_BUBBLE_RADIUS)
    assert bubble_radius(5, 0) == pytest.approx(smallest)


def test_other_and_unknown_labels_map_to_their_anchors() -> None:
    assert anatomical_site_base("Other") == OTHER_ANCHOR
    for value in ["unknown", None, " NA ", "nan"]:
        assert anatomical_site_base(value) == UNKNOWN_ANCHOR


@pytest.mark.parametrize("group", ["Thyroid gland", "Urine"])
def test_unmapped_body_site_group_raises(group: str) -> None:
    with pytest.raises(ValueError, match=repr(group)):
        anatomical_site_base(group)


def test_colorbar_sits_right_of_the_labels_and_bars() -> None:
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
        # A count bar can end right of every label once the axes is large (#711).
        (bar,) = ax.plot([0.58, 1.05], [0.6, 0.6], lw=7, solid_capstyle="round")
        # A wide title makes tight_layout hand the axes a box wider than it is tall,
        # which the equal aspect narrows only at draw time (#612).
        title = "Aerobe v2 Trait Enrichment in Increased Profiles: Anatomical Sites"
        ax.text(0.02, 1.005, title, fontsize=17, fontweight="bold", va="bottom")
        fig.tight_layout()
        mappable = ScalarMappable(norm=Normalize(-0.5, 0.5), cmap="RdBu_r")
        colorbar = add_colorbar_beside(fig, ax, mappable, texts, "trait z", [bar])
        fig.canvas.draw()
        colorbar_left = colorbar.ax.get_window_extent().x0
        assert all(text.get_window_extent().x1 < colorbar_left for text in texts)
        cap = bar.get_linewidth() * fig.dpi / 72 / 2
        assert bar.get_window_extent().x1 + cap < colorbar_left
        assert colorbar.ax.get_ylabel() == "trait z"
    finally:
        plt.close(fig)
