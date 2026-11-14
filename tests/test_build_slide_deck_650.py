"""Exercise deck verification orchestration and review artifacts (#650)."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
from types import ModuleType
from typing import Any
from unittest.mock import Mock

import pytest


def load_deck(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location("deck_650", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_build_with_entries(
    module: ModuleType,
    root: Path,
    patch: pytest.MonkeyPatch,
    entries: list[dict[str, Any]],
) -> tuple[list[dict[str, str]], dict[str, Mock], Path]:
    patch.setattr(module, "REPO_ROOT", root)
    table = root / "quoted.tsv"
    table.write_text("value\n2\n")
    metrics = {
        "contrast_classifier": "test",
        "matrix_taxa": 1,
        "trait_count": 1,
        "contrasts_by_key": {},
        "pooled_pairs": 1,
        "pooled_studies": 1,
        "ind_total_tests": 0,
    }

    def load_metrics(_root: Path) -> dict[str, object]:
        module.TABLES_READ.add(table)
        return metrics

    patch.setattr(module, "load_rerun_metrics", load_metrics)
    stubs: dict[str, object] = {
        "load_cellulose_captions": {"statements": {"overlap": "test", "fragile": ""}},
        "cellulose_revision": "test",
        "plate_revision": "test",
        "record_problems": [],
        "extract_committed_figures": None,
        "cellulose_decr_slides": [],
        "get_remaining_appendix_slides": [],
        "check_coverage_table_sync": None,
        "drawn_provenance": None,
        "arm_two_text": "test",
        "cover_table1_text": "test",
        "mcnemar_text": "test",
        "cover_wilcoxon_text": "test",
        "link_problem": None,
    }
    for name, value in stubs.items():
        patch.setattr(module, name, Mock(return_value=value))
    patch.setattr(
        module,
        "revision_problem",
        lambda _t, _r, rel, _p: "changed table" if rel == "quoted.tsv" else None,
    )
    for entry in entries:
        (root / str(entry["figure_file_path"])).write_bytes(b"image fixture")
    patch.setattr(module, "get_new_slides_metadata", Mock(return_value=entries))
    for relative, content in {
        "Justfile": "",
        f"{module.PLATES_DIR}plates_provenance.json": "{}",
        f"{module.CELLULOSE_DIR}/figures/figures_provenance.json": "{}",
        f"{module.CELLULOSE_DIR}/cellulose_decr_provenance.json": "{}",
        f"{module.DECK_DIR}/slides/slide_manifest_v1.tsv": "figure_file_path\n",
        f"{module.DECK_DIR}/slides/slide_manifest_v2.tsv": "figure_file_path\n",
        f"{module.DECK_DIR}/provenance/{module.PROVENANCE_NAME}": "figures: []\n",
    }.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    guards = {}
    for name in (
        "refuse_or_review",
        "check_results_cite_the_deck",
        "verified_links",
        "placed_slides",
        "check_slides",
        "plate_problems",
        "cellulose_problems",
    ):
        guards[name] = Mock(wraps=getattr(module, name))
        patch.setattr(module, name, guards[name])
    output = root / "output"
    module.build_deck(
        root,
        output,
        root / "provenance",
        trees={module.COMMIT: {"quoted.tsv": "old"}},
        allow_unverified=True,
    )
    with (output / "slide_manifest_v6.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows, guards, output


def exercise_build(module: ModuleType, root: Path, patch: pytest.MonkeyPatch) -> None:
    entries = [
        {
            "figure_file_path": f"{name}.png",
            "owner": "test",
            "title": name,
            "tag": "test",
            "takeaway": "test",
            "caveat": "test",
            "code_path": "",
            "quoted_tables": json.dumps(["quoted.tsv"] if name == "dependent" else []),
        }
        for name in ("dependent", "independent")
    ]
    rows, guards, output = run_build_with_entries(module, root, patch, entries)
    for name, guard in guards.items():
        assert guard.called, name
    assert guards["check_slides"].call_args.args[3] == module.TABLES_READ
    assert guards["refuse_or_review"].call_args.args[0]
    assert len(rows) == 2
    assert "quoted.tsv" in rows[0]["unverified"]
    assert rows[1]["unverified"] == ""
    assert json.loads(rows[0]["quoted_tables"]) == ["quoted.tsv"]
    assert "NOT VERIFIED" in module.revision_stamp(rows[0])
    validation = json.loads((output / "deck_validation_v6.json").read_text())
    assert validation["review_build"] is True
    assert {
        "html_v6",
        "html_default",
        "manifest_v6",
        "manifest_default",
        "provenance_yaml",
    } <= validation["outputs"].keys()
    review = Path(validation["outputs"]["review_html_v6"])
    assert review.name == "selected_figures_v6_review.html"
    assert review.read_text() == (output / "selected_figures_v6.html").read_text()
    assert "REVIEW BUILD" in review.read_text()


def test_build_guards_and_quoted_table_review(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = load_deck(Path(__file__).resolve().parents[1] / "scripts/build_slide_deck.py")
    exercise_build(module, tmp_path, monkeypatch)


def test_check_slides_normalizes_quoted_tables(tmp_path: Path) -> None:
    module = load_deck(Path(__file__).resolve().parents[1] / "scripts/build_slide_deck.py")
    table = tmp_path / "quoted.tsv"
    table.write_text("value\n2\n")
    tree: dict[str, str] = {"quoted.tsv": "0" * 40}
    cases: list[tuple[str, Any, bool]] = [
        ("empty_str.png", "", False),
        ("none_val.png", None, False),
        ("empty_seq.png", [], False),
        ("json_str.png", json.dumps(["quoted.tsv"]), True),
        ("list_seq.png", ["quoted.tsv"], True),
        ("tuple_seq.png", ("quoted.tsv",), True),
    ]
    entries: list[dict[str, Any]] = []
    for fig_rel, raw_quoted, _ in cases:
        fig_path = tmp_path / fig_rel
        fig_path.write_bytes(b"image fixture")
        tree[fig_rel] = module.git_blob_id(fig_path)
        entries.append({"figure_file_path": fig_rel, "quoted_tables": raw_quoted})
    omitted_fig = tmp_path / "omitted.png"
    omitted_fig.write_bytes(b"image fixture")
    tree["omitted.png"] = module.git_blob_id(omitted_fig)
    entries.append({"figure_file_path": "omitted.png"})

    unverified, quoted = module.check_slides(tmp_path, entries, {module.COMMIT: tree}, [table])
    expected_reason = f"differs from {module.SHORT_COMMIT}"
    assert quoted == [f"quoted.tsv: {expected_reason}"]
    assert "omitted.png" not in unverified
    for fig_rel, _, expect_flagged in cases:
        if expect_flagged:
            assert unverified[fig_rel] == [("quoted.tsv", expected_reason)]
        else:
            assert fig_rel not in unverified


def test_build_deck_normalizes_quoted_tables(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = load_deck(Path(__file__).resolve().parents[1] / "scripts/build_slide_deck.py")
    shared = [
        f"{module.DECK_DIR}/slides/slide_manifest_v1.tsv",
        f"{module.DECK_DIR}/slides/slide_manifest_v2.tsv",
        "quoted.tsv",
    ]
    cases: list[tuple[str, Any, list[str]]] = [
        ("blank_str.png", "", shared),
        ("none_val.png", None, shared),
        ("seq_quoted.png", ["quoted.tsv"], ["quoted.tsv"]),
        ("seq_empty.png", [], []),
        ("json_empty.png", "[]", []),
    ]
    entries: list[dict[str, Any]] = [
        {
            "figure_file_path": fig_rel,
            "owner": "test",
            "title": fig_rel,
            "tag": "test",
            "takeaway": "test",
            "caveat": "test",
            "code_path": "",
            "quoted_tables": raw_quoted,
        }
        for fig_rel, raw_quoted, _ in cases
    ]
    entries.append(
        {
            "figure_file_path": "omitted.png",
            "owner": "test",
            "title": "omitted.png",
            "tag": "test",
            "takeaway": "test",
            "caveat": "test",
            "code_path": "",
        }
    )
    rows, _, _ = run_build_with_entries(module, tmp_path, monkeypatch, entries)
    by_fig = {row["figure_file_path"]: row for row in rows}
    for fig_rel, _, expected_tables in cases:
        assert json.loads(by_fig[fig_rel]["quoted_tables"]) == expected_tables
        if expected_tables:
            assert "quoted.tsv" in by_fig[fig_rel]["unverified"]
        else:
            assert by_fig[fig_rel]["unverified"] == ""
    assert json.loads(by_fig["omitted.png"]["quoted_tables"]) == shared
    assert "quoted.tsv" in by_fig["omitted.png"]["unverified"]
