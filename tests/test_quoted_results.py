"""Reading the Arm 2, McNemar and calibration results the deck and plates quote (#579)."""

from __future__ import annotations

import pytest

from mtp.quoted_results import McNemar, OraCalls, screened_rejection_rate

ENRICHMENT_HEADER = [
    "# Tail: upper\n",
    "# FDR: True\n",
    "# P-Threshold: 0.001\n",
    "# Min-Fold: 1.5\n",
    "# Total-Tests-Performed: 1480390\n",
    "# Significant-Hits: 66240\n",
    "study_id\texperiment_id\n",
    "# Significant-Hits: 1\n",
]


def test_ora_counts_come_from_the_table_header() -> None:
    calls = OraCalls.from_header(ENRICHMENT_HEADER)
    assert calls == OraCalls(tests=1480390, calls=66240, p_threshold="0.001", fold="1.5")
    depletion = [line.replace("Min-Fold: 1.5", "Max-Fold: 0.5") for line in ENRICHMENT_HEADER]
    assert OraCalls.from_header(depletion).fold == "0.5"


def test_an_ora_header_without_counts_is_refused() -> None:
    with pytest.raises(ValueError, match="without its counts"):
        OraCalls.from_header(["# FDR: True\n", "study_id\n"])


def test_calls_without_fdr_control_are_refused() -> None:
    """The deck quotes them as BH q (#704)."""
    raw = [line.replace("FDR: True", "FDR: False") for line in ENRICHMENT_HEADER]
    with pytest.raises(ValueError, match="not FDR-controlled"):
        OraCalls.from_header(raw)


def _guild_rows(odds_ratio: str = "1.448276") -> list[dict[str, str]]:
    return [
        {"statistic": "studies_shifted", "comparison": "aerobe_inc", "value": "536", "p_value": ""},
        {
            "statistic": "discordant_studies",
            "comparison": "aerobe_up_anaerobe_down",
            "value": "420",
            "p_value": "",
        },
        {
            "statistic": "discordant_studies",
            "comparison": "anaerobe_up_aerobe_down",
            "value": "290",
            "p_value": "",
        },
        {
            "statistic": "mcnemar_odds_ratio",
            "comparison": "aerobe_vs_anaerobe_paired_signs (tested)",
            "value": odds_ratio,
            "p_value": "1.20881e-06",
        },
    ]


def test_mcnemar_reads_the_guild_table() -> None:
    test = McNemar.from_rows(_guild_rows())
    assert (test.aerobe_up_anaerobe_down, test.anaerobe_up_aerobe_down) == (420, 290)
    assert test.discordant == 710
    assert (test.odds_ratio, test.p_value) == (1.448276, 1.20881e-06)


def test_an_undefined_mcnemar_ratio_is_refused() -> None:
    """guild_maps writes no ratio when no study is reverse-discordant."""
    with pytest.raises(ValueError, match="no McNemar odds ratio"):
        McNemar.from_rows(_guild_rows(odds_ratio=""))


def test_the_screened_rejection_rate_is_read_by_criterion_and_procedure() -> None:
    rows = [
        {"criterion": "per_test_level", "procedure": "fraction_sign_test", "estimate": "0.0457"},
        {
            "criterion": "family_error",
            "procedure": "screened_fraction_sign_test",
            "estimate": "0.0893",
        },
        {
            "criterion": "per_test_level",
            "procedure": "screened_fraction_sign_test",
            "estimate": "0.0910",
            "replicates_recorded": "2800",
            "replicates_required": "2800",
        },
    ]
    assert screened_rejection_rate(rows) == 0.0910
    # The documented quick pass certifies nothing (#704)
    quick = [*rows[:2], {**rows[2], "replicates_recorded": "5"}]
    with pytest.raises(ValueError, match="ran 5 of the 2800 replicates"):
        screened_rejection_rate(quick)
    with pytest.raises(ValueError, match="no per-test level"):
        screened_rejection_rate(rows[:2])
