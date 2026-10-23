"""The grouped trait identifiers and their label check (#165, #289)."""

import json
from pathlib import Path

import pytest

from mtp.taxa_trait_matrix import NEGATED_PREDICATES, PHENOTYPE_EXCLUDED_TRAITS
from mtp.trait_terms import (
    EXPECTED_LABELS,
    check_term_labels,
    file_md5,
    grouped_identifiers,
    identify_node_table,
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


def _provenance(path: Path, nodes: Path) -> Path:
    """A sidecar pinning ``nodes`` by its real MD5, as `just download-kg` writes."""
    entry = {"files": {nodes.name: {"md5": file_md5(nodes), "size": nodes.stat().st_size}}}
    path.write_text(json.dumps(entry), encoding="utf-8")
    return path


def test_the_cli_fails_on_any_mismatch(tmp_path: Path) -> None:
    rows = list(EXPECTED_LABELS.items())
    good = _nodes(tmp_path / "good.tsv", rows)
    provenance = _provenance(tmp_path / "provenance.json", good)
    assert main(["--nodes", str(good), "--provenance", str(provenance)]) == 0

    bad = _nodes(tmp_path / "bad.tsv", [(term, "renamed") for term, _ in rows])
    # Pin the renamed table too, so the run fails on the labels rather than on
    # the hash - otherwise this would pass for the wrong reason.
    bad_provenance = _provenance(tmp_path / "bad_provenance.json", bad)
    assert main(["--nodes", str(bad), "--provenance", str(bad_provenance)]) == 1


def test_an_unpinned_node_table_stops_before_any_label_is_read(tmp_path: Path) -> None:
    """A different KG build must not be reported as a missing term (#289).

    Without the hash gate the only signal is `MISMATCH ... found None`, which
    reads as "the pinned build lost this term" when it means "this is not the
    pinned build".
    """
    rows = list(EXPECTED_LABELS.items())
    nodes = _nodes(tmp_path / "nodes.tsv", rows)
    other = _nodes(tmp_path / "other.tsv", [("METPO:9999", "something else")])
    provenance = _provenance(tmp_path / "provenance.json", other)

    # Every label is correct, so a label-only check would return 0.
    assert main(["--nodes", str(nodes), "--provenance", str(provenance)]) == 1
    # The override reports the labels anyway, and they pass.
    assert (
        main(["--nodes", str(nodes), "--provenance", str(provenance), "--allow-unpinned-nodes"])
        == 0
    )


def test_a_missing_provenance_sidecar_is_unknown_not_a_match(tmp_path: Path) -> None:
    """ "No md5 recorded" is not the same as "md5 agrees", and must not pass silently."""
    nodes = _nodes(tmp_path / "nodes.tsv", list(EXPECTED_LABELS.items()))
    identity = identify_node_table(nodes, tmp_path / "absent.json")

    assert identity.declared_md5 is None
    assert not identity.pinned
    assert identity.undetermined is not None
    assert main(["--nodes", str(nodes), "--provenance", str(tmp_path / "absent.json")]) == 1


def test_the_pinned_md5_is_read_from_the_sidecar_by_filename(tmp_path: Path) -> None:
    """The sidecar keys by file name, so the entry for another file must not match."""
    nodes = _nodes(tmp_path / "nodes.tsv", list(EXPECTED_LABELS.items()))
    sidecar = tmp_path / "provenance.json"
    sidecar.write_text(
        json.dumps({"files": {"merged-kg_edges.tsv": {"md5": file_md5(nodes)}}}),
        encoding="utf-8",
    )

    identity = identify_node_table(nodes, sidecar)

    assert identity.declared_md5 is None
    assert identity.actual_md5 == file_md5(nodes)


def test_the_strain_rollup_rule_identifiers_are_label_checked() -> None:
    assert set(PHENOTYPE_EXCLUDED_TRAITS) | set(NEGATED_PREDICATES) <= set(EXPECTED_LABELS)


def test_a_negated_relation_missing_from_the_rule_fails(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A KG that adds a "does not" relation must not slip past the full scope (#788)."""
    nodes = _nodes(
        tmp_path / "nodes.tsv", [*EXPECTED_LABELS.items(), ("METPO:2999999", "does not swim")]
    )
    provenance = _provenance(tmp_path / "provenance.json", nodes)

    assert main(["--nodes", str(nodes), "--provenance", str(provenance)]) == 1
    assert "METPO:2999999" in capsys.readouterr().out
