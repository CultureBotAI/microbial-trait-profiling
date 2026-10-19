"""Results the deck and the Tufte plates quote from tables other than the Table 1 sign tests.

The cover, Methods Plate II and Plate V typed Arm 2's test and call counts, the oxygen-guild
McNemar test and the post-selection rejection rate of the screened sign test (#579). A rerun
whose analysis code changed any of them left the text stale. Reading them here keeps the
text in step with the tables:

- ``trait_enrichment_results.tsv`` and ``trait_depletion_results.tsv`` (``just analyze`` and
  ``just analyze-depletion``) head their rows with ``# Key: value`` lines that give the tests
  and calls;
- ``guild_aerobe_anaerobe_correlation.tsv`` (``just guild-maps``) gives the McNemar test;
- ``calibration_criteria.tsv`` (``just calibrate``) gives the rejection rate.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

ORA_ENRICHMENT = Path("data/derived/trait_enrichment_results.tsv")
ORA_DEPLETION = Path("data/derived/trait_depletion_results.tsv")
GUILD_CORRELATION = Path("data/derived/guild_maps/guild_aerobe_anaerobe_correlation.tsv")
CALIBRATION_CRITERIA = Path("data/derived/calibration/calibration_criteria.tsv")


@dataclass(frozen=True)
class OraCalls:
    """One ORA table's tests and calls, with the thresholds that define a call."""

    tests: int
    calls: int
    p_threshold: str
    fold: str

    @classmethod
    def from_header(cls, lines: Iterable[str]) -> OraCalls:
        """Read the ``# Key: value`` lines heading the table; the rows are not read."""
        header: dict[str, str] = {}
        for line in lines:
            if not line.startswith("#"):
                break
            key, _, value = line[1:].partition(":")
            header[key.strip()] = value.strip()
        try:
            return cls(
                tests=int(header["Total-Tests-Performed"]),
                calls=int(header["Significant-Hits"]),
                p_threshold=header["P-Threshold"],
                fold=header.get("Min-Fold") or header["Max-Fold"],
            )
        except (KeyError, ValueError) as error:
            raise ValueError(f"an ORA table header without its counts ({error})") from error


@dataclass(frozen=True)
class McNemar:
    """The paired oxygen-guild test: studies whose aerobe and anaerobe shares moved in
    opposite directions, each way, and the exact McNemar test on them."""

    aerobe_up_anaerobe_down: int
    anaerobe_up_aerobe_down: int
    odds_ratio: float
    p_value: float

    @property
    def discordant(self) -> int:
        return self.aerobe_up_anaerobe_down + self.anaerobe_up_aerobe_down

    @classmethod
    def from_rows(cls, rows: Iterable[Mapping[str, str]]) -> McNemar:
        discordant: dict[str, int] = {}
        test: Mapping[str, str] | None = None
        for row in rows:
            if row["statistic"] == "discordant_studies":
                discordant[row["comparison"]] = int(row["value"])
            elif row["statistic"] == "mcnemar_odds_ratio":
                test = row
        if test is None or not test["value"] or not test["p_value"]:
            raise ValueError("the guild table has no McNemar odds ratio and p-value")
        try:
            return cls(
                aerobe_up_anaerobe_down=discordant["aerobe_up_anaerobe_down"],
                anaerobe_up_aerobe_down=discordant["anaerobe_up_aerobe_down"],
                odds_ratio=float(test["value"]),
                p_value=float(test["p_value"]),
            )
        except KeyError as error:
            raise ValueError(f"the guild table has no {error} discordant count") from error


def screened_rejection_rate(rows: Iterable[Mapping[str, str]]) -> float:
    """The pooled per-test rejection rate of the screened fraction sign test under the
    complete null: what screening the traits on the same data does to its size."""
    for row in rows:
        if (row["criterion"], row["procedure"]) == (
            "per_test_level",
            "screened_fraction_sign_test",
        ):
            return float(row["estimate"])
    raise ValueError("the calibration criteria have no per-test level for the screened sign test")
