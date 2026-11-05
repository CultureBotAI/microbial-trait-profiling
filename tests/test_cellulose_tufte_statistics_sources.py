"""Bind historical plate statistics to independent tracked sources (#1551, #1555).

These tests construct artists, not rendered files. Display rounding is applied to
source values; they do not recompute a statistical result. Historical fixture
extracts are verbatim tracked documents, with their source commits recorded beside them.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure

import test_cellulose_tufte_historical_labels as historical_labels
from test_cellulose_tufte_historical_labels import label_plates as label_plates

ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "methods/table1_log_ratio_validation.md"
ORA_AUDIT = ROOT / "docs/reviews/2026-10-19-ora-provenance-refresh.md"
HISTORICAL = ROOT / "tests/fixtures/cellulose_tufte_historical_sources"


def _row(path: Path, key: str, *, section: str | None = None) -> list[str]:
    source = path.read_text(encoding="utf-8")
    if section is not None:
        assert section in source, (path, section)
        source = source.split(section, 1)[1].split("\n## ", 1)[0]
    rows = [
        [cell.strip().strip("*") for cell in line.strip().strip("|").split("|")]
        for line in source.splitlines()
        if line.startswith("|")
    ]
    matching = [row for row in rows if row[0] == key]
    assert len(matching) == 1, (path, key, matching)
    return matching[0]


def _text(ax: Axes) -> list[str]:
    return [historical_labels._plain(artist.get_text()) for artist in ax.texts]


def _label_and_summary(ax: Axes, prefix: str) -> tuple[str, str]:
    labels = [
        artist
        for artist in ax.texts
        if artist.get_text().split("(", 1)[0].strip() == prefix.strip()
    ]
    assert len(labels) == 1, prefix
    label = labels[0]
    summaries = [
        artist.get_text()
        for artist in ax.texts
        if artist is not label and artist.get_position()[1] == label.get_position()[1]
    ]
    assert len(summaries) == 1, (prefix, summaries)
    return historical_labels._plain(label.get_text()), historical_labels._plain(summaries[0])


def _scientific(value: float) -> str:
    coefficient, exponent = f"{value:.1e}".split("e")
    superscripts = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"{coefficient}\u00d710{str(int(exponent)).translate(superscripts)}"


def _family_sizes(source: str) -> tuple[int, int]:
    fractions = source.split("## Original fraction directions", 1)[1].split("## All ", 1)[0]
    ratios = source.split("paired log-ratios (matched-corpus analog)", 1)[1].split(
        "## Implications", 1
    )[0]
    fraction_rows = re.findall(r"^\| [A-Z]{2,3} \| \d+/\d+/\d+ \|", fractions, re.MULTILINE)
    ratio_rows = re.findall(r"^\| [A-Z]{2,3}/[A-Z]{2,3} \|", ratios, re.MULTILINE)
    assert fraction_rows
    assert ratio_rows
    return len(fraction_rows), len(ratio_rows)


def test_plate_five_shared_corpus_and_cascade_match_the_pinned_methods(
    label_plates: dict[str, Figure],
) -> None:
    source = historical_labels._plain(VALIDATION.read_text(encoding="utf-8"))
    corpus = re.search(
        r"files: ([\d,]+) profile rows read, ([\d,]+) dropped as non-human and "
        r"([\d,]+) experiments dropped as unpaired, leaving ([\d,]+) retained profiles, "
        r"([\d,]+) matched experiments and ([\d,]+) unique studies\.",
        source,
    )
    assert corpus is not None
    raw, nonhuman, unpaired, profiles, pairs, studies = corpus.groups()
    retained = re.search(r"header columns \(3 identifier, [\d,]+ MediaDive, (\d+) retained", source)
    assert retained is not None
    figure = label_plates["V"]
    subtitles = [
        historical_labels._plain(artist.get_text())
        for artist in figure.texts
        if artist.get_text().startswith("All three arms start")
    ]
    assert subtitles == [
        "All three arms start from the same paired human trait matrix "
        f"({pairs} experiments across {studies} studies and {retained[1]} non-MediaDive "
        "traits), but test different sampling units, transformations, and null hypotheses."
    ]
    cascade = _text(figure.axes[0])
    assert f"{raw} profiles (all hosts, paired & unpaired)" in cascade
    assert f"2. Drop {nonhuman} non-human rows" in cascade
    assert f"3. Drop {unpaired} unpaired rows" in cascade
    assert f"{profiles} paired profiles = {pairs} pairs in {studies} studies" in cascade


@pytest.mark.parametrize("trait", ["FA", "HAL", "ARS", "SPH", "COC", "CEL"])
def test_plate_five_fraction_labels_match_the_tracked_validation_rows(
    label_plates: dict[str, Figure], trait: str
) -> None:
    # methods/table1_log_ratio_validation.md, "Original fraction directions" table.
    row = _row(VALIDATION, trait, section="## Original fraction directions")
    share = float(row[3].split()[0])
    q_value = float(row[5])
    q_text = _scientific(q_value) if trait in {"FA", "HAL", "CEL"} else f"{q_value:.4f}"
    stars = "**" if q_value <= 0.001 else "*"
    _, summary = _label_and_summary(label_plates["V"].axes[3], f"{trait} ")
    if trait == "CEL":
        expected = f"{share:.1%} up / {1 - share:.1%} down (q = {q_text} {stars})"
    else:
        expected = f"{share:.1%} (n={row[2]}, BH q = {q_text} {stars})"
    assert summary == expected


@pytest.mark.parametrize("ratio", ["FA/CEL", "FA/COC", "CEL/SPH", "CEL/ARS"])
def test_plate_five_ratio_labels_match_the_tracked_validation_rows(
    label_plates: dict[str, Figure], ratio: str
) -> None:
    # Same document, "All 15 paired log-ratios" table; no CI geometry assertion.
    row = _row(VALIDATION, ratio)
    share = float(row[5].split()[0])
    q_value = float(row[7])
    label, summary = _label_and_summary(label_plates["V"].axes[3], ratio.replace("/", " / "))
    if ratio == "FA/COC":
        expected = f"{share:.1%} up (n={row[1].split()[0]}, BY q = {q_value:.2f} n.s.)"
    elif ratio == "FA/CEL":
        expected = f"{share:.1%} up (n={row[1].split()[0]}, BY q = {q_value:.4f} *)"
    else:
        expected = f"{share:.1%} up / {1 - share:.1%} down (BY q = {q_value:.4f} *)"
    assert summary == expected
    if ratio != "FA/COC":
        effect = f"{float(row[4]):+.2f}".replace("-", "\u2212")
        assert f"({effect} log₂)" in label


def test_plate_five_retained_trait_family_matches_the_matrix_source(
    label_plates: dict[str, Figure],
) -> None:
    source = VALIDATION.read_text(encoding="utf-8")
    match = re.search(r"header columns \(3 identifier, ([\d,]+) MediaDive, (\d+) retained", source)
    assert match is not None
    media, retained = match.groups()
    columns = re.search(r"`: ([\d,]+)\s+header columns", source)
    assert columns is not None
    fractions, ratios = _family_sizes(source)
    assert (
        f"Column filter: {columns[1]} matrix columns = 3 ID columns + {media} MediaDive "
        f"growth-medium columns (excluded) + {retained} retained trait columns."
    ) in _text(label_plates["V"].axes[0])
    assert any(
        f"BH ({fractions} traits) & BY ({ratios} ratios, {retained} traits)" in text
        for text in _text(label_plates["V"].axes[2])
    )
    panel = _text(label_plates["V"].axes[3])
    assert f"Six Table 1 Single-Trait Fractions (BH across {fractions} traits):" in panel
    assert f"Key Table 1 Paired Log-Ratios (BY across {ratios} pairs; cancels T):" in panel


def test_plate_five_ora_numbers_match_the_historical_export_audit(
    label_plates: dict[str, Figure],
) -> None:
    enrichment = _row(ORA_AUDIT, "Pooled enrichment")
    depletion = _row(ORA_AUDIT, "Pooled depletion")
    assert enrichment[1:3] == depletion[1:3]
    profiles, tests, hits = enrichment[1:]
    retained = re.search(
        r"header columns \(3 identifier, [\d,]+ MediaDive, (\d+) retained",
        VALIDATION.read_text(encoding="utf-8"),
    )
    assert retained is not None
    assert (
        f"{profiles} profiles \u00d7 {retained[1]} traits "
        f"({int(tests.replace(',', '')) / 1_000_000:.2f}M tests; #131)"
    ) in _text(label_plates["V"].axes[1])
    cards = _text(label_plates["V"].axes[2])
    assert any(f"({profiles} \u00d7 {retained[1]} = {tests} tests)" in text for text in cards)
    assert any(
        f"BH across {tests} tests; {hits} enrichments" in text
        and f"& {depletion[3]} depletions" in text
        for text in cards
    )


@pytest.mark.parametrize(
    ("display", "state", "prefix"),
    [
        ("drawn", "untested", "drawn"),
        ("drawn", "zero_annotation_signal", "sampled-without-signal"),
        ("masked_below_min_studies", "insufficient_support", "masked_below_min_studies"),
        ("not_drawn", "not_sampled", "not_sampled"),
    ],
)
def test_plate_six_cell_counts_match_the_historical_inventory(
    label_plates: dict[str, Figure], display: str, state: str, prefix: str
) -> None:
    source = (HISTORICAL / "integrated_evidence_a3c32cf.md").read_text(encoding="utf-8")
    matches = re.findall(rf"^\| {display} \| {state} \| (\d+) \|$", source, re.MULTILINE)
    assert len(matches) == 1
    count = int(matches[0])
    total = int(_row(HISTORICAL / "integrated_evidence_a3c32cf.md", "oxygen_guild_cell_support")[1])
    assert (f"{prefix} ({count:,} cells, {count / total:.1%}) [test_status = {state}]") in _text(
        label_plates["VI"].axes[2]
    )


def test_plate_six_total_rows_match_the_historical_methods_snapshot(
    label_plates: dict[str, Figure],
) -> None:
    source = (HISTORICAL / "methods_integrated_evidence_a3c32cf.md").read_text(encoding="utf-8")
    totals = re.search(r"\*\*([\d,]+) rows: (\d+) method results and ([\d,]+) figure", source)
    assert totals is not None
    total, methods, cells = totals.groups()
    subtitles = [
        historical_labels._plain(artist.get_text())
        for artist in label_plates["VI"].texts
        if artist.get_text().startswith("Rather than pooling")
    ]
    assert subtitles == [
        "Rather than pooling incompatible p-values across sampling units, the pipeline aligns "
        "corpus denominators and records each method, threshold family, and grid-cell state "
        f"in a single {total}-row table."
    ]
    panel = _text(label_plates["VI"].axes[2])
    assert any(f"integrated_evidence.tsv ({total} Rows)" in text for text in panel)
    assert any(f"Method Rows ({methods} Rows" in text for text in panel)
    assert any(f"Grid Cells ({cells} Rows" in text for text in panel)


@pytest.mark.parametrize(
    ("source_key", "plate_key"),
    [
        ("directional consistency, fraction sign test", "table1_Burst_direction"),
        ("directional consistency, log-ratio sign test", "table1_paired_log_ratio"),
    ],
)
def test_plate_six_method_counts_match_the_historical_methods_snapshot(
    label_plates: dict[str, Figure], source_key: str, plate_key: str
) -> None:
    source = _row(HISTORICAL / "methods_integrated_evidence_a3c32cf.md", source_key)
    panel = label_plates["VI"].axes[2]
    keys = [artist for artist in panel.texts if artist.get_text() == plate_key]
    assert len(keys) == 1, (plate_key, keys)
    key = keys[0]
    row = [artist for artist in panel.texts if artist.get_position()[1] == key.get_position()[1]]
    for heading, expected in (
        ("Rows", f"{source[1]} rows"),
        ("q ≤ 0.05", f"{source[4]} / {source[1]}"),
        ("q ≤ 0.001", f"{source[5]} / {source[1]}"),
    ):
        headers = [artist for artist in panel.texts if artist.get_text() == heading]
        assert len(headers) == 1, heading
        column_x = headers[0].get_position()[0]
        cells = [
            historical_labels._plain(artist.get_text())
            for artist in row
            if artist.get_position()[0] == column_x
        ]
        assert cells == [expected], (plate_key, heading, cells)
