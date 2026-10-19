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
    studies_tied: int
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
    """One population's six fraction tests and 15 log-ratio tests.

    ``studies`` counts the population's studies. ``tested_studies`` counts those left after
    experiments with a zero annotation total are dropped, which every fraction test scores
    or counts as tied, so the fraction sign tests are quoted with it (#683). Log-ratio
    eligibility is decided pair by pair, so they are quoted with ``studies``.
    """

    name: str
    studies: int
    tested_studies: int
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
            studies_tied=int(row["studies_tied"]),
            p_value=float(row["p_value"]),
            q_value=float(row["q_value"]),
        )
    except (KeyError, ValueError) as error:
        raise ValueError(
            f"{row.get('population')} {row.get('numerator')}/{row.get('reference')}: "
            f"untested or malformed row ({error})"
        ) from error


#: The populations the deck and the plates quote. The table's other subsets, such as the
#: control-tier ones, are not read, so a row there that no text quotes cannot stop a build.
QUOTED = (ALL_STUDIES, ORIENTED)


def populations(rows: Iterable[Mapping[str, str]]) -> dict[str, Population]:
    """All studies and the oriented subset, from the comparison table. Each needs all six
    fractions and all 15 log-ratios, each once, so no text quotes part of a table or one
    of two copies of a row."""
    grouped: dict[str, list[Mapping[str, str]]] = {name: [] for name in QUOTED}
    for row in rows:
        if row["population"] in grouped:
            grouped[row["population"]].append(row)
    missing = [name for name, members in grouped.items() if not members]
    if missing:
        raise ValueError(f"the comparison table has no {' or '.join(missing)} rows")
    found = {}
    for name, members in grouped.items():
        fractions: dict[str, SignTest] = {}
        ratios: dict[frozenset[str], SignTest] = {}
        for row in members:
            test = sign_test(row)
            if row["comparison"] == "fraction":
                if test.numerator in fractions:
                    raise ValueError(f"{name}: two fraction rows for {test.label}")
                fractions[test.numerator] = test
            elif row["comparison"] == "log_ratio":
                pair = frozenset((test.numerator, test.reference))
                if pair in ratios:
                    raise ValueError(f"{name}: two log_ratio rows for {test.label}")
                ratios[pair] = test
            else:
                raise ValueError(f"{name}: unknown comparison {row['comparison']!r}")
        if set(fractions) != set(CODES) or set(ratios) != PAIRS:
            raise ValueError(f"{name}: expected the six Table 1 fractions and their 15 log-ratios")
        sizes = {int(row["population_studies"]) for row in members}
        if len(sizes) != 1:
            raise ValueError(f"{name}: rows disagree on the number of studies {sorted(sizes)}")
        tested = {test.studies_scored + test.studies_tied for test in fractions.values()}
        if len(tested) != 1:
            raise ValueError(
                f"{name}: the fraction tests disagree on their studies {sorted(tested)}"
            )
        found[name] = Population(name, sizes.pop(), tested.pop(), fractions, ratios)
    return found
