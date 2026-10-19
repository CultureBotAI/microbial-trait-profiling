"""Positions of body-site groups on the front-view body outline.

The aerobe and anaerobe v2 anatomical-site diagrams and the paired-contrast
infographic draw the same body outline in a unit square, centred at x = 0.30,
and place each body-site group at an approximate anatomical anchor.
"""

from __future__ import annotations

__all__ = [
    "OTHER_ANCHOR",
    "SITE_ANCHORS",
    "UNKNOWN_ANCHOR",
    "anatomical_site_base",
]

OTHER_ANCHOR = (0.16, 0.47)
UNKNOWN_ANCHOR = (0.16, 0.16)

# The first rule with a keyword in the lower-cased label wins, so "nasopharynx"
# precedes the neck rule's "pharynx".
SITE_ANCHORS: tuple[tuple[tuple[str, ...], tuple[float, float]], ...] = (
    (("nasopharynx", "airway"), (0.30, 0.86)),
    (("subgingival", "dental"), (0.275, 0.825)),
    (("saliva",), (0.325, 0.825)),
    (("mouth", "oral"), (0.30, 0.825)),
    (("throat", "pharynx"), (0.30, 0.79)),
    (("colon",), (0.285, 0.47)),
    (("caecum", "cecum"), (0.33, 0.46)),
    (("feces", "faeces", "stool"), (0.31, 0.37)),
    (("cervix", "uterus", "uterine"), (0.31, 0.31)),
    (("vagina", "urogenital"), (0.31, 0.25)),
)


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
