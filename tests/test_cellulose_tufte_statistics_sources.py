"""Bind selected historical plate statistics to independent tracked sources (#1551, #1555).

These tests construct artists, not rendered files. Display rounding is applied to
source values; they do not recompute a statistical result. Historical fixture
extracts are verbatim tracked documents, with their source commits recorded beside them.
"""

from __future__ import annotations

import hashlib
import re
import runpy
from collections.abc import Callable
from pathlib import Path

import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure

import test_cellulose_tufte_historical_labels as historical_labels
from mtp.local_path_patterns import compile_absolute_local_paths
from test_cellulose_tufte_historical_labels import label_plates as label_plates

ROOT = Path(__file__).resolve().parents[1]
HISTORICAL = ROOT / "tests/fixtures/cellulose_tufte_historical_sources"
ORA_AUDIT = HISTORICAL / "ora_provenance_refresh_781e92e4.md"
VALIDATION = HISTORICAL / "table1_log_ratio_validation_1e0d0eef.md"


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


def _fixture_manifest() -> dict[str, tuple[str, str, str]]:
    readme = (HISTORICAL / "README.md").read_text(encoding="utf-8")
    assert "| Fixture |" in readme
    manifest = readme.split("| Fixture |", 1)[1].split("\n\n", 1)[0]
    rows = [line.strip("| ").split("|") for line in manifest.splitlines() if line.startswith("| `")]
    assert rows, "Document every historical fixture and its original Git blob"
    records = {}
    for cells in rows:
        assert len(cells) == 4, (
            "Fixture table must record file, full commit, source path and blob ID"
        )
        filename, commit, path, blob = (cell.strip().strip("`") for cell in cells)
        assert re.fullmatch(r"[0-9a-f]{40}", commit), commit
        assert re.fullmatch(r"[0-9a-f]{40}", blob), blob
        assert filename not in records
        records[filename] = (commit, path, blob)
    return records


def _historical_source_names(directory: Path) -> list[str]:
    """Discover non-hidden sources recursively; only the root README is the manifest."""
    sources = []
    for path in directory.rglob("*"):
        relative = path.relative_to(directory)
        if not path.is_file():
            continue
        if relative == Path("README.md") or any(part.startswith(".") for part in relative.parts):
            continue
        sources.append(relative.as_posix())
    return sorted(sources)


def test_historical_fixture_manifest_covers_every_preserved_source() -> None:
    records = _fixture_manifest()
    assert set(records) == set(_historical_source_names(HISTORICAL))


@pytest.mark.parametrize(
    "fixture",
    _historical_source_names(HISTORICAL),
)
def test_historical_fixture_bytes_equal_the_documented_git_blob(fixture: str) -> None:
    expected = _fixture_manifest()[fixture][2]
    content = (HISTORICAL / fixture).read_bytes()
    header = f"blob {len(content)}\0".encode("ascii")
    actual = hashlib.sha1(header + content, usedforsecurity=False).hexdigest()
    assert actual == expected, f"Historical fixture changed: {fixture}"


def test_plate_five_validation_is_a_documented_historical_fixture() -> None:
    assert VALIDATION.parent == HISTORICAL, (
        "Plate V's methods-validation numbers must read the historical snapshot, not HEAD"
    )
    assert ORA_AUDIT.parent == HISTORICAL, (
        "Plate V's ORA numbers must read the historical snapshot, not HEAD"
    )
    assert _fixture_manifest()[VALIDATION.name][1] == "methods/table1_log_ratio_validation.md"


def test_plate_five_ora_audit_is_a_documented_historical_fixture() -> None:
    """The ORA source is frozen independently of subsequent live-document edits (#1993)."""
    assert ORA_AUDIT.parent == HISTORICAL, "Plate V's ORA source must be a historical snapshot"
    assert ORA_AUDIT.name == "ora_provenance_refresh_781e92e4.md"
    commit, source, _ = _fixture_manifest()[ORA_AUDIT.name]
    assert commit.startswith("781e92e4")
    assert source == "docs/reviews/2026-10-19-ora-provenance-refresh.md"


def test_historical_ora_snapshot_has_no_machine_local_paths() -> None:
    pattern = compile_absolute_local_paths(flags=0, classified=True)
    assert pattern.search(ORA_AUDIT.read_text(encoding="utf-8")) is None


def _require_guard_failure(check: Callable[[], None], message: str) -> None:
    failure = None
    try:
        check()
    except AssertionError as exc:
        failure = exc
    assert failure is not None, message


def test_historical_source_parent_guard_refuses_live_ora_audit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setitem(
        globals(), "ORA_AUDIT", ROOT / "docs/reviews/2026-10-19-ora-provenance-refresh.md"
    )
    _require_guard_failure(
        test_plate_five_validation_is_a_documented_historical_fixture,
        "The parent guard accepted the live ORA audit instead of its historical snapshot",
    )


def _write_historical_manifest(directory: Path, sources: dict[str, bytes]) -> None:
    """Make tiny source fixtures and their independent Git-blob manifest for guard tests."""
    rows = [
        "| Fixture | Source commit | Original repository path | Git blob ID |",
        "|---|---|---|---|",
    ]
    for filename, content in sources.items():
        path = directory / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        header = f"blob {len(content)}\0".encode("ascii")
        blob = hashlib.sha1(header + content, usedforsecurity=False).hexdigest()
        rows.append(f"| `{filename}` | `{'a' * 40}` | `docs/source.md` | `{blob}` |")
    (directory / "README.md").write_text("\n".join(rows) + "\n\n", encoding="utf-8")


@pytest.mark.parametrize(
    "sources",
    [
        {"source.md": b"root\n", "nested/source.md": b"nested\n"},
        {"first/source.md": b"first\n", "second/source.md": b"second\n"},
        {"nested/README.md": b"historical source, not the root manifest\n"},
    ],
    ids=["nested-source", "same-basename", "nested-readme"],
)
def test_historical_manifest_accepts_documented_nested_sources(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, sources: dict[str, bytes]
) -> None:
    _write_historical_manifest(tmp_path, sources)
    monkeypatch.setitem(globals(), "HISTORICAL", tmp_path)
    test_historical_fixture_manifest_covers_every_preserved_source()
    for filename in sources:
        test_historical_fixture_bytes_equal_the_documented_git_blob(filename)


def test_historical_manifest_refuses_undocumented_nested_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _write_historical_manifest(tmp_path, {"source.md": b"documented\n"})
    nested = tmp_path / "nested" / "unlisted.md"
    nested.parent.mkdir()
    nested.write_bytes(b"unlisted\n")
    monkeypatch.setitem(globals(), "HISTORICAL", tmp_path)
    _require_guard_failure(
        test_historical_fixture_manifest_covers_every_preserved_source,
        "The manifest guard ignored an undocumented nested historical source",
    )


@pytest.mark.parametrize("hidden", [".DS_Store", ".cache/source.md", "nested/.cache/source.md"])
def test_historical_manifest_ignores_hidden_entries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, hidden: str
) -> None:
    _write_historical_manifest(tmp_path, {"source.md": b"documented\n"})
    path = tmp_path / hidden
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"operating-system or tool metadata\n")
    monkeypatch.setitem(globals(), "HISTORICAL", tmp_path)
    test_historical_fixture_manifest_covers_every_preserved_source()


def test_historical_source_discovery_uses_relative_names_for_every_byte_case(
    tmp_path: Path,
) -> None:
    sources = {"source.md": b"root\n", "nested/source.md": b"nested\n"}
    _write_historical_manifest(tmp_path, sources)
    (tmp_path / ".DS_Store").write_bytes(b"ignored\n")
    discover = globals().get("_historical_source_names")
    assert callable(discover), "Manifest and byte tests must share recursive source discovery"
    assert discover(tmp_path) == sorted(sources)


def test_historical_byte_parametrization_includes_nested_sources(tmp_path: Path) -> None:
    """Evaluate the actual decorator in a tiny checkout, independently of helper assertions."""
    tests_dir = tmp_path / "checkout" / "tests"
    historical = tests_dir / "fixtures" / "cellulose_tufte_historical_sources"
    sources = {
        "first/source.md": b"first\n",
        "second/source.md": b"second\n",
        "nested/README.md": b"nested source\n",
    }
    _write_historical_manifest(historical, sources)
    (historical / ".DS_Store").write_bytes(b"metadata\n")
    hidden = historical / ".cache" / "source.md"
    hidden.parent.mkdir()
    hidden.write_bytes(b"metadata\n")
    copied_module = tests_dir / Path(__file__).name
    copied_module.write_bytes(Path(__file__).read_bytes())

    namespace = runpy.run_path(str(copied_module))
    byte_test = namespace.get("test_historical_fixture_bytes_equal_the_documented_git_blob")
    assert callable(byte_test)
    parametrizations = [
        mark
        for mark in getattr(byte_test, "pytestmark", ())
        if mark.name == "parametrize" and mark.args[0] == "fixture"
    ]
    assert len(parametrizations) == 1
    assert list(parametrizations[0].args[1]) == sorted(sources), (
        "The byte-test decorator must enumerate every nested source using its relative path"
    )


def test_historical_byte_guard_refuses_changed_nested_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    filename = "nested/source.md"
    _write_historical_manifest(tmp_path, {filename: b"original\n"})
    monkeypatch.setitem(globals(), "HISTORICAL", tmp_path)
    test_historical_fixture_bytes_equal_the_documented_git_blob(filename)
    (tmp_path / filename).write_bytes(b"modified\n")
    _require_guard_failure(
        lambda: test_historical_fixture_bytes_equal_the_documented_git_blob(filename),
        "The byte guard accepted a changed nested historical source",
    )
