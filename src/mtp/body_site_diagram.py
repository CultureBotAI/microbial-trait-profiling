"""Positions of body-site groups on the front-view body outline.

The aerobe and anaerobe v2 anatomical-site diagrams and the paired-contrast
infographic draw the same body outline in a unit square, centred at x = 0.30,
and place each body-site group at an approximate anatomical anchor. The outline
faces the viewer, so the subject's left side is at x > 0.30.

Bubbles are sized in the same unit-square coordinates (`bubble_radius`), so the
spacing between anchors holds in the full-page diagrams and in the infographic's
smaller panel alike (#620).
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # matplotlib is the optional "viz" extra, so it is imported for annotations only;
    # the figure scripts pass in their own figure, axes and artists.
    from matplotlib.axes import Axes
    from matplotlib.cm import ScalarMappable
    from matplotlib.colorbar import Colorbar
    from matplotlib.figure import Figure
    from matplotlib.text import Text

__all__ = [
    "MAX_BUBBLE_RADIUS",
    "MIN_BUBBLE_AREA_FRACTION",
    "OTHER_ANCHOR",
    "SITE_ANCHORS",
    "UNKNOWN_ANCHOR",
    "add_colorbar_beside",
    "anatomical_site_base",
    "bubble_radius",
    "leader_line",
]

OTHER_ANCHOR = (0.16, 0.47)
UNKNOWN_ANCHOR = (0.16, 0.16)

# The first rule with a keyword in the lower-cased label wins, so "nasopharynx"
# precedes the neck rule's "pharynx" and "colorect" keeps "Colorectum" with the colon.
# Leader lines run down and to the right, so the face sites rise from left to right:
# dental plaque at the jaw, mouth and oral cavity side by side (#612), and saliva at
# the parotid gland beside the ear. Each line then passes below the bubbles to its
# right. Stomach and the descending colon are on the subject's left, the caecum on
# the right; the small intestine is central and the rectum sits low on the midline,
# above feces (#620).
SITE_ANCHORS: tuple[tuple[tuple[str, ...], tuple[float, float]], ...] = (
    (("nasopharynx", "airway"), (0.300, 0.888)),
    (("subgingival", "dental"), (0.243, 0.822)),
    (("saliva",), (0.364, 0.890)),
    (("mouth",), (0.281, 0.846)),
    (("oral",), (0.321, 0.850)),
    (("throat", "pharynx"), (0.300, 0.800)),
    (("lung", "bronch"), (0.258, 0.685)),
    (("stomach", "gastric"), (0.348, 0.590)),
    (("colon", "colorect"), (0.350, 0.490)),
    (("caecum", "cecum"), (0.252, 0.460)),
    (("rectum", "rectal"), (0.300, 0.432)),
    (("intestin", "ileum", "jejun", "duoden"), (0.300, 0.520)),
    (("skin",), (0.435, 0.500)),
    (("feces", "faeces", "stool"), (0.300, 0.372)),
    (("cervix", "uterus", "uterine"), (0.300, 0.316)),
    (("vagina", "urogenital"), (0.300, 0.268)),
)

# Radius of the largest group's bubble, in outline units. The smallest bubble keeps
# MIN_BUBBLE_AREA_FRACTION of its area, the 95/995 floor of the point-sized markers
# that these bubbles replace.
MAX_BUBBLE_RADIUS = 0.031
MIN_BUBBLE_AREA_FRACTION = 95 / 995
# Gap, in axes fractions, between the rightmost label and a colorbar.
COLORBAR_GAP = 0.04


def anatomical_site_base(value: object) -> tuple[float, float]:
    """Return the outline anchor for a body-site group label.

    `Other` and `unknown` sit off the body. Any other label without a rule raises,
    so a new top body-site group cannot land on the arm unnoticed (#601).
    """
    label = str(value or "").strip().lower()
    if label in {"", "nan", "na", "unknown"}:
        return UNKNOWN_ANCHOR
    if label == "other":
        return OTHER_ANCHOR
    for keywords, anchor in SITE_ANCHORS:
        if any(keyword in label for keyword in keywords):
            return anchor
    msg = f"No body-outline position for body-site group {value!r}; add a rule to SITE_ANCHORS."
    raise ValueError(msg)


def bubble_radius(count: float, max_count: float) -> float:
    """Return the radius, in outline units, of a group's bubble.

    The area grows with the square root of `count / max_count`, from
    MIN_BUBBLE_AREA_FRACTION of the largest bubble's area up to the whole of it.
    """
    share = min(max(count / max_count, 0.0), 1.0) if max_count > 0 else 0.0
    area_fraction = MIN_BUBBLE_AREA_FRACTION + (1 - MIN_BUBBLE_AREA_FRACTION) * math.sqrt(share)
    return MAX_BUBBLE_RADIUS * math.sqrt(area_fraction)


def leader_line(
    anchor: tuple[float, float],
    radius: float,
    label_y: float,
    elbow_x: float,
    end_x: float,
) -> tuple[list[float], list[float]]:
    """Return the x and y points of the line from a bubble to its label.

    The line runs toward (`elbow_x`, `label_y`), then level to `end_x`. It starts on
    the bubble's edge, so the diagrams draw it above the bubbles: a line that crosses
    another bubble stays visible and still ends at its own.
    """
    x, y = anchor
    length = math.dist(anchor, (elbow_x, label_y))
    start_x = x + radius * (elbow_x - x) / length
    start_y = y + radius * (label_y - y) / length
    return [start_x, elbow_x, end_x], [start_y, label_y, label_y]


def add_colorbar_beside(
    fig: Figure,
    ax: Axes,
    mappable: ScalarMappable,
    texts: Iterable[Text],
    label: str,
) -> Colorbar:
    """Add a colorbar as an inset of `ax`, right of where `texts` end (#612).

    Call it once the layout is final, after `tight_layout`: the inset is placed from
    the drawn extent of the labels, so the colorbar cannot cover their values.
    """
    # tight_layout leaves an equal-aspect axes at its full box until the next draw;
    # narrowing it first measures the labels where they will be drawn.
    ax.apply_aspect()
    right_edge = max(text.get_window_extent().x1 for text in texts)
    left = float(ax.transAxes.inverted().transform((right_edge, 0.0))[0]) + COLORBAR_GAP
    colorbar = fig.colorbar(mappable, cax=ax.inset_axes((left, 0.22, 0.025, 0.56)))
    colorbar.set_label(label)
    return colorbar
