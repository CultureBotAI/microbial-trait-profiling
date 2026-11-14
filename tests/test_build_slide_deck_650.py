"""Exercise deck verification orchestration and review artifacts (#650)."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
from types import ModuleType
from unittest.mock import Mock

import pytest


def load_deck(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location("deck_650", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exercise_build(module: ModuleType, root: Path, patch: pytest.MonkeyPatch) -> None:
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
    entries = []
    for name in ("dependent", "independent"):
        (root / f"{name}.png").write_bytes(b"image fixture")
        entries.append(
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
        )
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
    for name, guard in guards.items():
        assert guard.called, name
    assert guards["check_slides"].call_args.args[3] == module.TABLES_READ
    assert guards["refuse_or_review"].call_args.args[0]
    with (output / "slide_manifest_v6.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
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
