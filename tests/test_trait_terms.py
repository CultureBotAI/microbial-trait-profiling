"""The grouped trait identifiers and their label check (#165)."""

from pathlib import Path

from bsdbng.trait_terms import (
    EXPECTED_LABELS,
    check_term_labels,
    grouped_identifiers,
    main,
)


def _nodes(path: Path, rows: list[tuple[str, str]]) -> Path:
    lines = ["id\tcategory\tname"] + [
        f"{term}\tbiolink:PhenotypicQuality\t{name}" for term, name in rows
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_every_grouped_identifier_has_an_expected_label() -> None:
    assert grouped_identifiers() <= set(EXPECTED_LABELS)


def test_a_wrong_or_missing_label_is_reported(tmp_path: Path) -> None:
    nodes = _nodes(tmp_path / "nodes.tsv", [("METPO:1", "aerobic"), ("METPO:2", "anaerobe")])
    checks = {
        c.identifier: c
        for c in check_term_labels(
            nodes, {"METPO:1": "aerobic", "METPO:2": "anaerobic", "METPO:3": "x"}
        )
    }
    assert checks["METPO:1"].ok
    assert not checks["METPO:2"].ok
    assert checks["METPO:3"].found is None


def test_the_cli_fails_on_any_mismatch(tmp_path: Path) -> None:
    rows = list(EXPECTED_LABELS.items())
    good = _nodes(tmp_path / "good.tsv", rows)
    assert main(["--nodes", str(good)]) == 0
    bad = _nodes(tmp_path / "bad.tsv", [(term, "renamed") for term, _ in rows])
    assert main(["--nodes", str(bad)]) == 1
