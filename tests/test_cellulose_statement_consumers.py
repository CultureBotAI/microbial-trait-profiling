"""Keep exported cellulose sentences usable by the current deck and Plate IV (#1867)."""

from __future__ import annotations

import importlib.util
import json
import shutil
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from matplotlib.figure import Figure

from mtp.cellulose_decr import figures
from test_cellulose_nominal_labels import synthetic_results

REPO = Path(__file__).resolve().parents[1]
ANALYSIS = Path("data/analysis/cellulose_decr")


def _load_script(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / f"{name}.py")
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def consumers() -> tuple[ModuleType, ModuleType]:
    return _load_script("build_slide_deck"), _load_script("generate_tufte_plates")


@pytest.fixture(params=["committed", "synthetic-sensitive", "synthetic-stable"])
def exported(
    request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> tuple[Path, figures.Results, dict[str, Any], dict[str, int]]:
    """Run the real JSON writer; replace only drawing and file-save boundaries."""
    analysis = tmp_path / ANALYSIS
    analysis.mkdir(parents=True)
    if request.param == "committed":
        for name in figures.FIGURE_INPUTS:
            shutil.copyfile(REPO / ANALYSIS / name, analysis / name)
        counts = {"overlap": 2, "crc_compartments": 4, "fragile": 7}
    else:
        results = synthetic_results()
        sensitive = request.param == "synthetic-sensitive"
        if sensitive:
            selected = results.scenarios["short_condition"].eq("HIV") & results.scenarios[
                "scenario"
            ].eq("named_treated_cases_excluded")
            results.scenarios.loc[selected, "binom_p_two_sided"] = 0.2
        (analysis / "cellulose_decr_summary.json").write_text(
            json.dumps(results.summary), encoding="utf-8"
        )
        tables = {
            "pairwise_jaccard_significance.tsv": results.pairwise,
            "sensitivity_definitions_jaccard.tsv": results.sensitivity,
            "compartment_jaccard_and_direction_stats.tsv": results.compartments,
            "compartment_direction_interaction.tsv": results.interaction,
            "cellulose_decreased_taxa_catalog.tsv": results.catalog,
            "condition_definition_sensitivity.tsv": results.scenarios,
            "condition_definition_interaction.tsv": results.scenario_interaction,
        }
        for name, frame in tables.items():
            frame.to_csv(analysis / name, sep="\t", index=False)
        # The fixture has no descriptive benchmark; its caption is outside this test.
        (analysis / "trait_convergence_benchmark.tsv").write_text("trait\n", encoding="utf-8")
        counts = {"overlap": 2, "crc_compartments": 4, "fragile": 3 if sensitive else 0}

    results = figures.load(analysis)
    for name in (
        "graphical_convergence",
        "graphical_core",
        "graphical_compartments",
        "graphical_benchmark",
        "detailed_nulls",
        "detailed_compartments",
        "detailed_core_and_benchmark",
    ):
        monkeypatch.setattr(figures, name, lambda _: Figure(figsize=(1, 1)))
    for name in ("caption_shared_taxa", "caption_ranking", "caption_shared_by_two"):
        monkeypatch.setattr(
            figures,
            name,
            lambda _: figures.Caption(
                "Descriptive", "Fixture", "No test", figures.FIGURE_INPUTS[0]
            ),
        )

    @contextmanager
    def saved_book(path: Path, **_kwargs: Any) -> Iterator[None]:
        yield None
        path.write_bytes(b"book save boundary\n")

    def saved_figure(_figure: Figure, directory: Path, stem: str, _book: Any) -> list[Path]:
        directory.mkdir(parents=True, exist_ok=True)
        paths = [directory / f"{stem}{suffix}" for suffix in (".svg", ".pdf", ".png")]
        for path in paths:
            path.write_bytes(b"figure save boundary\n")
        return paths

    monkeypatch.setattr(figures, "PdfPages", saved_book)
    monkeypatch.setattr(figures, "save", saved_figure)
    assert figures.main(["--analysis-dir", str(analysis)]) == 0
    saved = json.loads((analysis / "figures/figure_captions.json").read_text(encoding="utf-8"))
    return tmp_path, results, saved, counts


def test_main_statements_remain_usable_by_real_deck_and_plate_consumers(
    exported: tuple[Path, figures.Results, dict[str, Any], dict[str, int]],
    consumers: tuple[ModuleType, ModuleType],
) -> None:
    root, _results, saved, counts = exported
    deck, plates = consumers
    statements = saved["statements"]
    fragile = statements["fragile"]
    sentence = deck.fragile_sentence({"cellulose_statements": statements})
    # The deck's HTML note uses this exact expression and appends a full stop.
    note = statements["overlap"] + (f"; {fragile}" if fragile else "") + "."
    callout = plates.cellulose_decr_notes(root)["plate4_callout"]
    for text in (sentence, note, callout):
        assert ".." not in text
        assert ".;" not in text
        assert "Nominal p-values:" not in text
        assert "This display scope" not in text
    if counts["fragile"]:
        assert fragile
        assert sentence == f" In the cellulose deep-dive, {fragile}."
        bullet = callout.split("• Cellulose deep-dive: ", 1)[1]
        assert len(bullet.splitlines()) <= 2
        assert " ".join(line.strip() for line in bullet.splitlines()) == fragile
    else:
        assert fragile == ""
        assert sentence == ""
        assert "Cellulose deep-dive:" not in callout


def test_main_retains_each_statement_scope_separately(
    exported: tuple[Path, figures.Results, dict[str, Any], dict[str, int]],
) -> None:
    _root, results, saved, counts = exported
    assert "statement_scopes" in saved, "Separating a scope must not discard it"
    statements = saved["statements"]
    scopes = saved["statement_scopes"]
    writers = {
        "overlap": figures.overlap_statement,
        "crc_compartments": figures.crc_compartment_statement,
        "fragile": figures.fragile_statement,
    }
    assert set(statements) == set(writers)
    assert set(scopes) == set(writers)
    for name, writer in writers.items():
        assert statements[name] == writer(results)
        assert "Nominal p-values:" not in statements[name]
        if statements[name]:
            assert not statements[name].endswith((".", ";"))
            expected = (
                f"Nominal p-values: {counts[name]} finite endpoints; 0 unavailable. "
                "This display scope is not a correction family."
            )
            assert scopes[name] == expected
            assert writer(results, include_scope=True) == f"{statements[name]}. {scopes[name]}"
        else:
            assert scopes[name] == ""
