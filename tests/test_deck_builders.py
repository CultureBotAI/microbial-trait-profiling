"""The selected-figure deck builders in scripts/ (#556, #573)."""

import importlib.util
import json
import shutil
from pathlib import Path
from types import ModuleType

import pytest

REPO = Path(__file__).resolve().parents[1]
CELLULOSE = REPO / "data" / "analysis" / "cellulose_decr"


def _load(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / f"{name}.py")
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


deck = _load("build_slide_deck")
plates = _load("generate_tufte_plates")

COMMIT = "0123456789" * 4


# --------------------------------------------------------------------------- deck


@pytest.mark.parametrize(
    ("recorded", "expected"),
    [(COMMIT, COMMIT), (f"{COMMIT}+dirty", None), ("unrecorded", None), ("", None), (None, None)],
)
def test_cellulose_slides_link_only_a_clean_drawing_commit(
    tmp_path: Path, recorded: str | None, expected: str | None
) -> None:
    """A branch name moves and a modified tree has no commit to link (#556)."""
    figures = tmp_path / deck.CELLULOSE_DIR / "figures"
    figures.mkdir(parents=True)
    record = {"code_commit": recorded} if recorded != "" else {}
    (figures / "figures_provenance.json").write_text(json.dumps(record))
    assert deck.cellulose_revision(tmp_path) == expected


def test_cellulose_slides_carry_their_revision() -> None:
    caption = {"title": "t", "takeaway": "a", "caveat": "c", "evidence": "x.tsv"}
    captions = {"figures": {"figures/fig1_x.svg": caption}, "statements": {}}
    (slide,) = deck.cellulose_decr_slides(captions, "figures", "Appendix", "Marcin", COMMIT)
    assert slide["link_ref"] == COMMIT
    assert COMMIT[:12] in slide["source_note"]
    (slide,) = deck.cellulose_decr_slides(captions, "figures", "Appendix", "Marcin", None)
    assert slide["link_ref"] is None
    assert "no revision is linked" in slide["source_note"]


def test_stamp_names_the_drawing_commit_of_cellulose_figures() -> None:
    cellulose = "data/analysis/cellulose_decr/figures/fig1_x.svg"
    rerun = {"figure_file_path": "trait-visualization/x.png", "commit": deck.COMMIT}
    assert deck.revision_stamp(rerun) == f"Rerun revision {deck.SHORT_COMMIT}"
    drawn = {"figure_file_path": cellulose, "commit": COMMIT}
    assert deck.revision_stamp(drawn) == f"Figures drawn at {COMMIT[:12]}"
    unpinned = {"figure_file_path": cellulose, "commit": ""}
    assert deck.revision_stamp(unpinned) == "Figures drawn without a clean commit"


def test_figure_folders_match_whole_names() -> None:
    base = "data/analysis/cellulose_decr"
    assert f"{base}/figures/a.svg".startswith(deck.CELLULOSE_FIGURES)
    assert f"{base}/figures_graphical/a.svg".startswith(deck.CELLULOSE_FIGURES)
    assert not f"{base}/figures_tufte/a.svg".startswith(deck.CELLULOSE_FIGURES)
    assert f"{base}/figures_tufte/a.svg".startswith(deck.UNLINKED_FIGURES)
    assert not f"{base}/figures_extra/a.svg".startswith(deck.UNLINKED_FIGURES)


def test_cellulose_provenance_replaces_a_carried_over_record() -> None:
    """Records from the source YAML named a generator and tables that do not exist (#556)."""
    record = {
        "figure_id": "fig_slide_146",
        "generator_script_path": "data/analysis/cellulose_decr/plot_cellulose_decr_figures.py",
        "just_recipe": "uv run python data/analysis/cellulose_decr/plot_cellulose_decr_figures.py",
        "command_series": ["uv run python data/analysis/cellulose_decr/plot_cellulose_decr.py"],
        "pipeline_step_ids": [],
        "upstream_raw_inputs": ["data/raw/merged-kg_nodes.tsv"],
        "input_data_files": ["data/analysis/cellulose_decr/cellulose_decr_core17_taxa.tsv"],
        "input_data_sha256s": [],
        "methods_doc_path": "data/analysis/cellulose_decr/EXECUTIVE_SUMMARY.md",
    }
    deck.cellulose_figure_provenance(record, {"pairwise_jaccard_significance.tsv": "f" * 64})
    assert record["figure_id"] == "fig_slide_146"
    assert record["generator_script_path"] == "src/mtp/cellulose_decr/figures.py"
    assert record["just_recipe"] == "just cellulose-decr-figures"
    assert record["command_series"][-3:] == [
        "just contrast-orientation",
        "just cellulose-decr",
        "just cellulose-decr-figures",
    ]
    assert record["input_data_files"] == [
        "data/analysis/cellulose_decr/pairwise_jaccard_significance.tsv"
    ]
    assert record["input_data_sha256s"] == ["f" * 64]
    assert record["methods_doc_path"] == "methods/cellulose_decr.md"
    text = json.dumps(record)
    for stale in ("plot_cellulose_decr", "core17", "EXECUTIVE_SUMMARY", "merged-kg_nodes"):
        assert stale not in text


def test_plate_two_text_reads_the_stratum_tables() -> None:
    def row(label: str, q: float, mean: float) -> dict[str, str]:
        return {
            "trait_label": label,
            "study_q_value": str(q),
            "study_mean_increased_minus_decreased": str(mean),
        }

    pooled = [
        row("cellulose catabolic process GO:0030245", 0.01, -0.2),
        row("coccus shaped METPO:1000668", 0.02, 0.3),
        row("anaerobic METPO:1000603", 0.5, -0.4),
        row("euryhaline METPO:1000627", 0.04, 0.1),
        row("arsenate reduction arsenate_reduction", 0.03, 0.2),
        row("aerobic METPO:1000602", 0.01, 0.25),
    ]
    feces = [row("x X:1", 0.01, 0.1), row("y Y:1", 0.2, 0.1)]
    text = deck.paired_contrast_text({"paired_contrasts": {"pooled": pooled, "feces": feces}})
    assert text == (
        "Of 6 traits tested pooled, 5 pass study-level BH q ≤ 0.05 (1 in Feces): 1 with mean "
        "Δz < 0, higher in decreased signatures (cellulose catabolic process), and 4 with "
        "Δz > 0, led by coccus shaped, aerobic and arsenate reduction."
    )


@pytest.mark.parametrize(("p_value", "verdict"), [("0.0800", "is not shown"), ("0.01", "is shown")])
def test_pd_direction_caveat_quotes_the_clustered_model(p_value: str, verdict: str) -> None:
    """The figure's title says direction-specific; the study-clustered model decides (#590)."""
    model = {"coefficient": "0.4098", "ci_low": "-0.0492", "ci_high": "0.8687", "p_value": p_value}
    caveat = deck.pd_direction_caveat(model)
    assert "+0.41 cellulose z (95% CI -0.05 to +0.87" in caveat
    assert f"a PD-specific direction {verdict} at p < 0.05" in caveat


def test_pd_direction_caveat_keeps_two_significant_digits() -> None:
    model = {"coefficient": "0.41", "ci_low": "-0.05", "ci_high": "0.87", "p_value": "0.08004"}
    assert "p = 0.080)" in deck.pd_direction_caveat(model)


# --------------------------------------------------------------------------- plates


def _cellulose_repo(tmp_path: Path, fragile: str | None) -> Path:
    target = tmp_path / plates.CELLULOSE_DIR
    (target / "figures").mkdir(parents=True)
    for name in ("compartment_direction_interaction.tsv", "cellulose_decr_summary.json"):
        shutil.copy(CELLULOSE / name, target / name)
    if fragile is not None:
        captions = {"figures": {}, "statements": {"fragile": fragile}}
        (target / "figures" / "figure_captions.json").write_text(json.dumps(captions))
    return tmp_path


def test_plate_four_callout_names_definition_dependence_and_fits(tmp_path: Path) -> None:
    """Plate IV stars HIV cellulose, so its callout says when that depends on definitions."""
    fragile = (
        "the HIV result depends on the condition definition (study-level p 0.0020\N{EN DASH}0.18)"
    )
    callout = plates.cellulose_decr_notes(_cellulose_repo(tmp_path, fragile))["plate4_callout"]
    lines = callout.splitlines()
    assert max(len(line) for line in lines) <= 88
    assert sum(line.startswith("• ") for line in lines) == 3
    assert f"• Cellulose deep-dive: {fragile}" in " ".join(line.strip() for line in lines)


def test_plate_four_callout_without_definition_dependence(tmp_path: Path) -> None:
    callout = plates.cellulose_decr_notes(_cellulose_repo(tmp_path, ""))["plate4_callout"]
    assert "deep-dive" not in callout
    assert sum(line.startswith("• ") for line in callout.splitlines()) == 2


def test_plates_need_the_cellulose_captions(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="just cellulose-decr-figures"):
        plates.cellulose_decr_notes(_cellulose_repo(tmp_path, None))
