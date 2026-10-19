"""Table 1 study sign tests as the deck and the Tufte plates quote them (#579).

``just contrast-orientation`` writes ``disease_control_comparison.tsv``. For every paired
human study (``historical_human_paired``) and for each disease-versus-healthy subset, it
gives the study sign test of the six Table 1 trait fractions (Benjamini-Hochberg over the
six) and of their 15 pairwise log-ratios (Benjamini-Yekutieli over the 15). The deck and
the plates used to type these numbers for one contrast classifier (569 oriented studies
under v2), so a rerun on another classifier left them stale. Reading them here keeps the
text in step with the table.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

COMPARISON = Path("data/derived/contrast_orientation/disease_control_comparison.tsv")
CONTRAST_PROVENANCE = Path("data/derived/contrast_orientation/contrast_orientation_provenance.json")
ALL_STUDIES = "historical_human_paired"
ORIENTED = "disease_vs_healthy_oriented"
CELLULOSE = "cellulose_catabolic_process"
ALPHA = 0.05
#: The six Table 1 traits as the table names them, with the codes the plates print.
CODES = {
    "facultatively_aerobic": "FA",
    "cellulose_catabolic_process": "CEL",
    "coccus_shaped": "COC",
    "sphere_shaped": "SPH",
    "arsenate_reduction": "ARS",
    "moderately_halophilic": "HAL",
}
PAIRS = frozenset(frozenset(pair) for pair in itertools.combinations(CODES, 2))


@dataclass(frozen=True)
class SignTest:
    """The study sign test of one trait fraction, or of the log-ratio of two traits."""

    numerator: str
    reference: str
    studies_scored: int
    studies_up: int
    studies_down: int
    p_value: float
    q_value: float

    @property
    def label(self) -> str:
        """'CEL' for a fraction, 'CEL/HAL' for a log-ratio, in the table's order."""
        code = CODES[self.numerator]
        return f"{code}/{CODES[self.reference]}" if self.reference else code

    @property
    def share_up(self) -> float:
        return self.studies_up / self.studies_scored

    @property
    def share_down(self) -> float:
        return self.studies_down / self.studies_scored

    def passes(self, alpha: float = ALPHA) -> bool:
        return self.q_value <= alpha


@dataclass(frozen=True)
class Population:
    """One population's six fraction tests and 15 log-ratio tests."""

    name: str
    studies: int
    fractions: Mapping[str, SignTest]
    log_ratios: Mapping[frozenset[str], SignTest]

    def ratio(self, a: str, b: str) -> SignTest:
        """The log-ratio of two traits, in whichever order the table holds it."""
        return self.log_ratios[frozenset((a, b))]

    def ratio_survivors(self, alpha: float = ALPHA) -> list[SignTest]:
        """The log-ratios at BY q ≤ alpha, in table order."""
        return [test for test in self.log_ratios.values() if test.passes(alpha)]


def sign_test(row: Mapping[str, str]) -> SignTest:
    try:
        return SignTest(
            numerator=row["numerator"],
            reference=row["reference"],
            studies_scored=int(row["studies_scored"]),
            studies_up=int(row["studies_up"]),
            studies_down=int(row["studies_down"]),
            p_value=float(row["p_value"]),
            q_value=float(row["q_value"]),
        )
    except (KeyError, ValueError) as error:
        raise ValueError(
            f"{row.get('population')} {row.get('numerator')}/{row.get('reference')}: "
            f"untested or malformed row ({error})"
        ) from error


def populations(rows: Iterable[Mapping[str, str]]) -> dict[str, Population]:
    """The comparison table's populations, which must include all studies and the oriented
    subset. Each needs all six fractions and all 15 log-ratios, so no text quotes part of
    a table."""
    grouped: dict[str, list[Mapping[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["population"], []).append(row)
    found = {}
    for name, members in grouped.items():
        fractions: dict[str, SignTest] = {}
        ratios: dict[frozenset[str], SignTest] = {}
        for row in members:
            test = sign_test(row)
            if row["comparison"] == "fraction":
                fractions[test.numerator] = test
            elif row["comparison"] == "log_ratio":
                ratios[frozenset((test.numerator, test.reference))] = test
            else:
                raise ValueError(f"{name}: unknown comparison {row['comparison']!r}")
        if set(fractions) != set(CODES) or set(ratios) != PAIRS:
            raise ValueError(f"{name}: expected the six Table 1 fractions and their 15 log-ratios")
        sizes = {int(row["population_studies"]) for row in members}
        if len(sizes) != 1:
            raise ValueError(f"{name}: rows disagree on the number of studies {sorted(sizes)}")
        found[name] = Population(name, sizes.pop(), fractions, ratios)
    missing = [name for name in (ALL_STUDIES, ORIENTED) if name not in found]
    if missing:
        raise ValueError(f"the comparison table has no {' or '.join(missing)} rows")
    return found
