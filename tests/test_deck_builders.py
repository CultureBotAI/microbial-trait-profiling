"""The selected-figure deck builders in scripts/ (#556, #573)."""

import hashlib
import importlib.util
import itertools
import json
import shutil
from dataclasses import replace
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from mtp.table1_directions import ALL_STUDIES, CODES, ORIENTED, populations
from mtp.table1_directions import CELLULOSE as CELLULOSE_TRAIT

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


def test_plate_slides_stamp_the_drawing_commit() -> None:
    plate = {
        "figure_file_path": f"{deck.PLATES_DIR}tufte_plate2.png",
        "commit": deck.COMMIT,
        "code_commit": COMMIT,
    }
    assert deck.revision_stamp(plate) == (
        f"Rerun revision {deck.SHORT_COMMIT} · plates drawn at {COMMIT[:12]}"
    )
    assert deck.revision_stamp({**plate, "code_commit": ""}).endswith(
        "plates drawn without a clean commit"
    )


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
    # `just taxa-trait-matrix` reads both KG files from data/kg-microbe (#609).
    assert record["upstream_raw_inputs"] == [
        "data/raw/full_dump.csv",
        "data/kg-microbe/merged-kg_nodes.tsv",
        "data/kg-microbe/merged-kg_edges.tsv",
    ]
    text = json.dumps(record)
    for stale in ("plot_cellulose_decr", "core17", "EXECUTIVE_SUMMARY", "data/raw/merged-kg"):
        assert stale not in text


def _paired(label: str, q: object, mean: float, **extra: str) -> dict[str, str]:
    return {
        "trait_label": label,
        "study_q_value": str(q),
        "study_mean_increased_minus_decreased": str(mean),
        **extra,
    }


def test_plate_two_text_reads_the_stratum_tables() -> None:
    pooled = [
        _paired("cellulose catabolic process GO:0030245", 0.01, -0.2),
        _paired("coccus shaped METPO:1000668", 0.02, 0.3),
        _paired("anaerobic METPO:1000603", 0.5, -0.4),
        _paired("euryhaline METPO:1000627", 0.04, 0.1),
        _paired("arsenate reduction arsenate_reduction", 0.03, 0.2),
        _paired("aerobic METPO:1000602", 0.01, 0.25),
    ]
    feces = [_paired("x X:1", 0.01, 0.1), _paired("y Y:1", 0.2, 0.1)]
    metrics = {"paired_contrasts": {"pooled": pooled, "feces": feces}}
    assert deck.paired_contrast_text(metrics) == (
        "Of 6 traits tested pooled, 5 pass study-level BH q ≤ 0.05 (1 in Feces): in 1 the "
        "annotation z-score is higher in decreased signatures (mean Δz < 0: cellulose catabolic "
        "process); in 4 it is higher in increased signatures (Δz > 0), led by coccus shaped, "
        "aerobic and arsenate reduction."
    )
    assert deck.paired_contrast_title(metrics) == (
        "Tufte Plate II: 5 traits pass study-level BH pooled, by direction and body site"
    )


def test_plate_two_text_edge_cases() -> None:
    """One side empty, one trait, a missing q and a zero mean Δz (#608)."""
    pooled = [_paired("aerobic METPO:1000602", 0.01, 0.25), _paired("x X:1", "", -0.3)]
    metrics = {"paired_contrasts": {"pooled": pooled, "feces": []}}
    assert deck.paired_contrast_text(metrics) == (
        "Of 1 trait tested pooled, 1 passes study-level BH q ≤ 0.05 (0 in Feces): in 1 the "
        "annotation z-score is higher in increased signatures (Δz > 0), led by aerobic."
    )
    level = {"paired_contrasts": {"pooled": [_paired("z Z:1", 0.01, 0.0)], "feces": []}}
    assert deck.paired_contrast_text(level).endswith(": in 1 the mean Δz is 0.")
    none = {"paired_contrasts": {"pooled": [_paired("z Z:1", 0.5, 0.2)], "feces": []}}
    assert deck.paired_contrast_text(none) == (
        "Of 1 trait tested pooled, 0 pass study-level BH q ≤ 0.05 (0 in Feces)."
    )
    # A literal NaN q is not a test, and one passing trait "passes" (#617).
    nan = {"paired_contrasts": {"pooled": [*pooled, _paired("n N:1", "NaN", 0.2)], "feces": []}}
    assert deck.paired_contrast_text(nan).startswith("Of 1 trait tested pooled, 1 passes")
    assert deck.paired_contrast_title(metrics).startswith("Tufte Plate II: 1 trait passes")
    # Slide 5 names a compound as Plate II does.
    raffinose = _paired("raffinose CHEBI:16634", 0.01, -0.2, trait="METPO:2000006|CHEBI:16634")
    compound = {"paired_contrasts": {"pooled": [raffinose], "feces": []}}
    assert "(mean Δz < 0: raffinose (carbon source))" in deck.paired_contrast_text(compound)


NITRATE = "METPO:2000103|pathways:nitrate_reduction"


def _nitrate(t: float, p: float, q: float, studies: int) -> list[dict[str, str]]:
    return [
        _paired(
            "nitrate reduction nitrate_reduction",
            q,
            t / 10,
            trait=NITRATE,
            study_t_statistic=str(t),
            study_p_value=str(p),
            study_count=str(studies),
        )
    ]


@pytest.mark.parametrize(
    ("saliva", "expected"),
    [
        (
            (-2.30, 0.0252, 0.979),
            "Saliva's nitrate reduction has the opposite sign to Feces' at a nominal p only "
            "(t = -2.30, p = 0.025, q = 0.98; Feces t = +3.97, q = 0.0013)",
        ),
        ((2.5, 0.001, 0.01), "the same sign as Feces' at BH q ≤ 0.05"),
        ((-1.0, 0.32, 0.99), "the opposite sign to Feces' without reaching p ≤ 0.05"),
    ],
)
def test_plate_two_caveat_reads_saliva_and_feces(
    saliva: tuple[float, float, float], expected: str
) -> None:
    """A nominal saliva sign is not headlined as a site difference (#607)."""
    metrics = {
        "paired_contrasts": {
            "feces": _nitrate(3.97, 8.04e-5, 0.00126, 594),
            "feces_excluded": _nitrate(-1.14, 0.256, 0.994, 365),
            "saliva": _nitrate(*saliva, 55),
        }
    }
    caveat = deck.paired_contrast_caveat(metrics)
    assert caveat.startswith("Without Feces (365 studies), 0 of 1 trait pass study BH q ≤ 0.05")
    assert expected in caveat
    assert caveat.endswith("and no test compares the two sites.")


def test_plate_two_caveat_needs_the_nitrate_rows() -> None:
    metrics = {
        "paired_contrasts": {
            "feces": _nitrate(3.97, 8.04e-5, 0.00126, 594),
            "feces_excluded": _nitrate(-1.14, 0.256, 0.994, 365),
            "saliva": [],
        }
    }
    with pytest.raises(ValueError, match="saliva paired-contrast table needs one"):
        deck.paired_contrast_caveat(metrics)


def test_plate_five_text_names_the_exceptions() -> None:
    """Fecal strata are not all aerobe-up; the slide says which are not (#610)."""

    def stratum(condition: str, site: str, score: float, rank: float) -> dict[str, str]:
        return {
            "condition": condition,
            "body_site_display": site,
            "aerobe_up_anaerobe_down_score": str(score),
            "rank_score": str(rank),
        }

    rows = [
        stratum("C", "Lower airway / lung", 0.97, 2.0),
        stratum("D", "Oral cavity / dental", -1.01, 1.5),
        stratum("A", "Feces", 1.18, 3.0),
        stratum("B", "Feces", -0.51, 1.0),
        stratum("E", "Feces", 0.1, 0.5),
    ]
    assert deck.oxygen_switch_text({"oxygen_switch": rows}) == (
        "Of the 4 strata shown, 1 of 2 fecal strata score aerobe-up / anaerobe-down and 1 of 2 "
        "other strata score anaerobe-up / aerobe-down. The exceptions: B · Feces -0.51 and "
        "C · Lower airway / lung +0.97."
    )
    clean = {"oxygen_switch": [rows[2], rows[1]]}
    assert "exceptions" not in deck.oxygen_switch_text(clean)
    # The selection follows rank_score, not the table's row order (#617).
    up = [stratum(f"U{i}", "Feces", 0.5, float(i)) for i in range(12)]
    chosen = deck.oxygen_switch_strata(up)
    assert [r["condition"] for r in chosen] == [f"U{i}" for i in range(11, 1, -1)]
    only_other = {"oxygen_switch": [rows[1]]}
    assert deck.oxygen_switch_text(only_other) == (
        "Of the 1 stratum shown, 1 of 1 other stratum scores anaerobe-up / aerobe-down."
    )


def _model(coefficient: float, p: float) -> dict[str, str]:
    return {
        "coefficient": str(coefficient),
        "ci_low": str(coefficient - 0.46),
        "ci_high": str(coefficient + 0.46),
        "p_value": str(p),
    }


@pytest.mark.parametrize(
    ("model", "positive", "verdict"),
    [
        ((0.4098, 0.08004), (0.1859, 0.0553), "a PD-specific direction is not shown at p < 0.05"),
        ((0.41, 0.01), (0.19, 0.02), "a PD-specific direction is shown at p < 0.05."),
        ((0.41, 0.01), (0.19, 0.2), "is shown at p < 0.05 for one measure only"),
        ((-0.41, 0.04), (0.19, 0.2), "points the other way"),
    ],
)
def test_pd_direction_caveat_quotes_the_clustered_model(
    model: tuple[float, float], positive: tuple[float, float], verdict: str
) -> None:
    """The figure's title says direction-specific; the clustered model decides (#590, #608)."""
    caveat = deck.pd_direction_caveat(_model(*model), _model(*positive))
    assert f"{model[0]:+.2f} cellulose z" in caveat
    assert "the share of profiles with cellulose z > 0" in caveat
    assert verdict in caveat


def test_pd_direction_caveat_gives_the_share_in_percentage_points() -> None:
    """The positive-z outcome is 0/1, so its term is a share, with its interval (#617)."""
    model = {"coefficient": "0.4098", "ci_low": "-0.0492", "ci_high": "0.8687", "p_value": "0.08"}
    positive = {
        "coefficient": "0.1859",
        "ci_low": "-0.0042",
        "ci_high": "0.3761",
        "p_value": "0.0553",
    }
    assert "+18.6 percentage points (95% CI -0.4 to +37.6, p = 0.055)" in (
        deck.pd_direction_caveat(model, positive)
    )


def test_p_values_do_not_round_across_the_threshold() -> None:
    assert deck.p_text(0.08004) == "0.080"
    assert deck.p_text(0.0495) == "0.0495"
    assert deck.p_text(0.04996) == "0.04996"
    assert deck.p_text(0.049999) == "0.049999"
    assert deck.p_text(0.00126) == "0.0013"
    assert deck.p_text(0.979) == "0.98"
    assert deck.signed(-0.004) == "0.0"
    assert deck.signed(18.59) == "+18.6"
    assert deck.signed(-0.004, 2) == "0.00"


def test_plate_records_name_the_plate_script() -> None:
    """The source YAML named a plate script outside the repository (#621)."""
    old_script = "data/derived/figure_review_2026-09-21/slides/generate_tufte_plates.py"
    record = {
        "generator_script_path": old_script,
        "just_recipe": f"python {old_script}",
        "command_series": [*deck.UPSTREAM_COMMANDS, "just analyze", f"uv run python {old_script}"],
    }
    inputs = {"trait-visualization/a.tsv": "aa", "data/derived/b.tsv": "bb"}
    deck.plate_provenance(record, "methods/table1_paired_log_ratios.md", inputs)
    assert record["generator_script_path"] == "scripts/generate_tufte_plates.py"
    assert record["methods_doc_path"] == "methods/table1_paired_log_ratios.md"
    assert record["just_recipe"] == "just figure-review-plates"
    assert record["command_series"][-5:] == [
        "just analyze",
        "just contrast-orientation",
        "just cellulose-decr",
        "just cellulose-decr-figures",
        "just figure-review-plates",
    ]
    assert "data/derived/figure_review" not in json.dumps(record)
    assert record["input_data_files"] == list(inputs)
    assert record["input_data_sha256s"] == ["aa", "bb"]
    # A rebuild reads the previous build's record and must not append the recipes again (#639)
    before = json.dumps(record)
    deck.plate_provenance(record, "methods/table1_paired_log_ratios.md", inputs)
    assert json.dumps(record) == before
    # trait-profile-investigation reads what trait-heatmaps writes (#681)
    investigation = {"command_series": ["just trait-profile-investigation", "just analyze"]}
    deck.plate_provenance(investigation, "", {})
    assert investigation["command_series"][:2] == [
        "just trait-heatmaps",
        "just trait-profile-investigation",
    ]


def test_carried_over_records_get_the_current_commands() -> None:
    """The source YAML described the upstream recipes with a nonexistent mtp.cli (#621)."""
    series = [
        "just download  # uv run python -m mtp.cli download",
        "just ingest  # uv run python -m mtp.cli ingest",
        "just study-tsvs  # uv run python -m mtp.cli study-tsvs",
        "just taxa-trait-matrix  # uv run python -m mtp.cli taxa-trait-matrix",
        "just study-trait-profiles  # uv run python -m mtp.cli study-trait-profiles",
        "just analyze  # uv run python -m mtp.cli analyze --mode enrichment",
        "uv run python scripts/x.py",
    ]
    assert deck.current_commands(series) == [
        *deck.UPSTREAM_COMMANDS,
        "just analyze",
        "uv run python scripts/x.py",
    ]
    assert not any("mtp.cli" in command for command in deck.current_commands(series))


def test_plate_slides_link_their_code_at_the_drawing_commit(tmp_path: Path) -> None:
    """The plate script is not in the rerun revision the plates' tables come from (#613)."""
    plates_dir = tmp_path / deck.PLATES_DIR
    plates_dir.mkdir(parents=True)
    assert deck.plate_revision(tmp_path) is None
    (plates_dir / "plates_provenance.json").write_text(json.dumps({"code_commit": COMMIT}))
    assert deck.plate_revision(tmp_path) == COMMIT

    plate = {
        "code_path": "scripts/generate_tufte_plates.py",
        "methods_path": "methods/table1_paired_log_ratios.md",
        "evidence_path": "trait-visualization/x.tsv",
        "code_ref": COMMIT,
    }
    links = deck.slide_links(plate, f"{deck.PLATES_DIR}tufte_plate3.png")
    assert links["code_url"].endswith(f"/blob/{COMMIT}/scripts/generate_tufte_plates.py")
    assert links["methods_url"].endswith(f"/blob/{deck.COMMIT}/methods/table1_paired_log_ratios.md")
    assert links["figure_url"] == "NA"
    assert links["commit"] == deck.COMMIT
    unpinned = deck.slide_links({**plate, "code_ref": None}, "x.png")
    assert unpinned["code_url"] == ""
    assert deck.slide_links({"code_path": "a.py"}, "x.png")["code_url"].endswith(
        f"/blob/{deck.COMMIT}/a.py"
    )


@pytest.mark.parametrize("argv", [["--over"], ["--out", "x"], ["--prov", "x"]])
def test_deck_options_are_never_abbreviated(tmp_path: Path, argv: list[str]) -> None:
    """The recipe's guards match whole option names (#617)."""
    # Every folder is in tmp_path, so a regression cannot rebuild the real deck (#623).
    safe = ["--input-root", str(tmp_path), "--output-dir", str(tmp_path / "slides")]
    safe += ["--provenance-dir", str(tmp_path / "provenance")]
    with pytest.raises(SystemExit) as stopped:
        deck.main([*safe, *argv])
    # argparse's usage error, not a later guard
    assert stopped.value.code == 2


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


def test_plate_labels_name_carbon_sources(tmp_path: Path) -> None:
    """A compound's object_name alone names no trait (#607)."""
    columns = tmp_path / plates.TRAIT_COLUMNS
    columns.parent.mkdir(parents=True)
    columns.write_text(
        "column_id\tobject_name\tpredicate_name\n"
        "METPO:2000006|CHEBI:16634\traffinose\tuses as carbon source\n"
        "METPO:2000103|pathways:nitrate_reduction\tnitrate_reduction\tcapable of\n"
        "METPO:1000606\tobligately aerobic\t\n"
    )
    assert plates.trait_labels(tmp_path) == {
        "METPO:2000006|CHEBI:16634": "raffinose (carbon source)",
        "METPO:2000103|pathways:nitrate_reduction": "nitrate reduction",
        "METPO:1000606": "obligately aerobic",
    }


def test_plate_two_rows_must_be_in_the_pooled_table() -> None:
    plates.check_plate_two_rows(["a", "b"], {"a": {}, "b": {}})
    with pytest.raises(ValueError, match=r"missing from the pooled table: \['b'\]"):
        plates.check_plate_two_rows(["a", "b"], {"a": {}})


def test_plates_record_their_drawing_commit_and_inputs(tmp_path: Path) -> None:
    inputs = {"b.png": {"t/y.tsv": "2", "t/x.tsv": "1"}, "a.png": {}}
    outputs = {"b.png": "4", "a.png": "3"}
    path = plates.write_plates_provenance(tmp_path, inputs, COMMIT, outputs, {"n": 17})
    assert json.loads(path.read_text()) == {
        "script": "scripts/generate_tufte_plates.py",
        "code_commit": COMMIT,
        "rerun_commit": plates.RERUN_COMMIT,
        "plates": ["a.png", "b.png"],
        "inputs": {"a.png": {}, "b.png": {"t/x.tsv": "1", "t/y.tsv": "2"}},
        "outputs": {"a.png": "3", "b.png": "4"},
        "counts": {"n": 17},
    }


def test_plates_and_deck_cite_one_rerun() -> None:
    """The plate footers and the deck stamps name the same revision (#637)."""
    assert plates.RERUN_COMMIT == deck.COMMIT


def test_plate_inputs_are_recorded_by_root(tmp_path: Path) -> None:
    table = tmp_path / "trait-visualization" / "t.tsv"
    table.parent.mkdir()
    table.write_text("a\tb\n1\t2\n")
    plates.INPUTS_READ.clear()
    plates.read_tsv(table)
    assert plates.input_record(tmp_path) == {
        "trait-visualization/t.tsv": hashlib.sha256(table.read_bytes()).hexdigest()
    }


def test_plate_and_figure_options_are_never_abbreviated(tmp_path: Path) -> None:
    """The recipes' guards match whole option names (#617)."""
    from mtp.cellulose_decr import figures

    with pytest.raises(SystemExit) as stopped:
        figures.parse_args(["--analysis", "/tmp/x"])
    assert stopped.value.code == 2
    with pytest.raises(SystemExit) as stopped:
        plates.main(["--input-root", str(tmp_path), "--output-dir", str(tmp_path), "--code", "x"])
    assert stopped.value.code == 2


def test_plates_need_the_cellulose_captions(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="just cellulose-decr-figures"):
        plates.cellulose_decr_notes(_cellulose_repo(tmp_path, None))


# --------------------------------------------------------------------------- final review


def test_a_missing_table_stops_the_text_that_quotes_it(tmp_path: Path) -> None:
    """An empty table gave slide text such as "the 0-test family" (#639)."""
    with pytest.raises(FileNotFoundError, match="missing"):
        deck.read_tsv(tmp_path / "absent.tsv")


def test_git_blob_id_of_known_content(tmp_path: Path) -> None:
    hello = tmp_path / "hello.txt"
    hello.write_text("hello\n")
    assert deck.git_blob_id(hello) == "ce013625030ba8dba906f756967f9e9ca394464a"


def _trees(tmp_path: Path, entries: dict[str, str], ref: str = COMMIT) -> dict[str, dict[str, str]]:
    folder = tmp_path / "trees"
    folder.mkdir(exist_ok=True)
    listing = "".join(f"100644 blob {blob}\t{name}\0" for name, blob in entries.items())
    (folder / f"{ref}.tree").write_text(listing)
    trees: dict[str, dict[str, str]] = deck.load_git_trees(folder)
    return trees


def test_revision_problems_name_what_differs(tmp_path: Path) -> None:
    """Slides were stamped and linked at a revision whatever the file was (#637)."""
    figure = tmp_path / "figure.png"
    figure.write_text("drawn\n")
    trees = _trees(tmp_path, {"t/figure.png": deck.git_blob_id(figure), "t/old.png": "0" * 40})
    assert deck.revision_problem(trees, COMMIT, "t/figure.png", figure) is None
    assert deck.revision_problem(trees, COMMIT, "t/old.png", figure) == (
        f"differs from {COMMIT[:12]}"
    )
    assert deck.revision_problem(trees, COMMIT, "t/new.png", figure) == (
        f"not committed at {COMMIT[:12]}"
    )
    assert deck.revision_problem(trees, "f" * 40, "t/figure.png", figure) == (
        f"no listing of {'f' * 12}"
    )
    assert deck.revision_problem(trees, None, "t/figure.png", figure) == (
        "drawn without a clean commit"
    )
    assert deck.revision_problem(trees, COMMIT, "t/figure.png", tmp_path / "gone") == "missing"


def test_plate_problems_need_the_rerun_and_unchanged_inputs(tmp_path: Path) -> None:
    table = tmp_path / "trait-visualization" / "t.tsv"
    table.parent.mkdir()
    table.write_text("a\n1\n")
    sha = hashlib.sha256(table.read_bytes()).hexdigest()
    ignored = tmp_path / "data" / "derived" / "d.tsv"
    ignored.parent.mkdir(parents=True)
    ignored.write_text("b\n2\n")
    inputs = {
        "trait-visualization/t.tsv": sha,
        "data/derived/d.tsv": hashlib.sha256(ignored.read_bytes()).hexdigest(),
    }
    plates_dir = tmp_path / "plates"
    plates_dir.mkdir()
    for name in ("p.png", "p.svg"):
        (plates_dir / name).write_text(name)
    outputs = {
        name: hashlib.sha256((plates_dir / name).read_bytes()).hexdigest()
        for name in ("p.png", "p.svg")
    }
    record = {
        "rerun_commit": deck.COMMIT,
        "plates": ["p.png"],
        "inputs": {"p.png": inputs},
        "outputs": outputs,
    }
    trees = _trees(tmp_path, {"trait-visualization/t.tsv": deck.git_blob_id(table)}, deck.COMMIT)
    assert deck.plate_problems(tmp_path, record, trees, None, plates_dir) == {}
    # A gitignored input is checked against the recorded hash only.
    ignored.write_text("changed\n")
    assert deck.plate_problems(tmp_path, record, trees, None, plates_dir) == {
        "p.png": ["data/derived/d.tsv: changed since the plate was drawn"]
    }
    other = {**record, "rerun_commit": COMMIT}
    assert deck.plate_problems(tmp_path, other, trees, None, plates_dir) == {
        "p.png": [f"plates drawn for rerun {COMMIT[:12]}, not {deck.SHORT_COMMIT}"]
    }


def test_plates_must_record_their_inputs_and_be_the_plates_drawn(tmp_path: Path) -> None:
    """A plate record written before #637, or a plate overwritten after drawing, is a
    problem before anything is written, not a crash after (#637)."""
    plates_dir = tmp_path / "plates"
    plates_dir.mkdir()
    for name in ("p.png", "p.svg"):
        (plates_dir / name).write_text(name)
    old = {"rerun_commit": deck.COMMIT, "plates": ["p.png"]}
    assert deck.plate_problems(tmp_path, old, {}, None, plates_dir) == {
        "p.png": ["no inputs recorded", "PNG: no hash recorded", "SVG: no hash recorded"]
    }
    drawn = {"p.png": "0" * 64, "p.svg": hashlib.sha256(b"p.svg").hexdigest()}
    record = {**old, "inputs": {"p.png": {}}, "outputs": drawn}
    assert deck.plate_problems(tmp_path, record, {}, None, plates_dir) == {
        "p.png": ["PNG: not the plate that was drawn"]
    }


def test_gitignored_outputs_are_not_checked_against_a_revision() -> None:
    assert deck.committed("trait-visualization/x.tsv")
    assert deck.committed("data/analysis/cellulose_decr/pairwise_jaccard_significance.tsv")
    assert not deck.committed("data/derived/x.tsv")
    assert not deck.committed("data/analysis/cellulose_decr/figures/figure_captions.json")
    # A file the revision does commit under data/derived is checked (#637).
    assert deck.committed("data/derived/x.tsv", {"data/derived/x.tsv": "0" * 40})


def test_results_must_cite_figures_in_the_deck() -> None:
    """A dropped figure's results were re-attached to the next figure's slide (#638)."""
    doc = {
        "figures": [
            {"figure_id": "fig_slide_004", "figure_file_path": "a.png"},
            {"figure_id": "fig_slide_005", "figure_file_path": "b.png"},
        ],
        "statistical_results": [{"figure_ids": ["fig_slide_004"]}, {"figure_ids": []}],
    }
    assert deck.stale_figure_ids(doc, {"a.png", "b.png"}) == []
    assert deck.stale_figure_ids(doc, {"b.png"}) == ["fig_slide_004"]


def test_unverified_slides_say_so() -> None:
    row = {"figure_file_path": "trait-visualization/x.png", "commit": deck.COMMIT}
    assert deck.revision_stamp({**row, "unverified": "x.png: differs"}).endswith("NOT VERIFIED")
    assert "NOT VERIFIED" not in deck.revision_stamp({**row, "unverified": ""})


def _root(tmp_path: Path) -> Path:
    root = tmp_path / "root"
    root.mkdir()
    return root


def _same_classifier(root: Path) -> None:
    """The contrast classifier this checkout's cellulose_decr tables were built on."""
    provenance = REPO / "data/analysis/cellulose_decr/cellulose_decr_provenance.json"
    version = json.loads(provenance.read_text())["contrast_classifier_version"]
    contrasts = root / plates.CONTRAST_PROVENANCE
    contrasts.parent.mkdir(parents=True)
    contrasts.write_text(json.dumps({"classifier_version": version}))


def _deck_args(tmp_path: Path, root: Path) -> list[str]:
    return [
        "--input-root",
        str(root),
        "--output-dir",
        str(tmp_path / "slides"),
        "--provenance-dir",
        str(tmp_path / "provenance"),
    ]


def test_deck_refuses_to_replace_a_deck(tmp_path: Path) -> None:
    root = _root(tmp_path)
    existing = tmp_path / "slides" / "selected_figures_v6.html"
    existing.parent.mkdir()
    existing.write_text("old deck")
    with pytest.raises(SystemExit, match="--overwrite"):
        deck.main(_deck_args(tmp_path, root))
    assert existing.read_text() == "old deck"


def test_deck_needs_git_listings_or_a_review_build(tmp_path: Path) -> None:
    root = _root(tmp_path)
    with pytest.raises(SystemExit, match="--allow-unverified"):
        deck.main(_deck_args(tmp_path, root))
    assert not (tmp_path / "slides").exists()


def test_plates_refuse_to_replace_plates(tmp_path: Path) -> None:
    root = _root(tmp_path)
    out = tmp_path / "plates"
    out.mkdir()
    (out / "old.png").write_bytes(b"old")
    with pytest.raises(SystemExit, match="--overwrite"):
        plates.main(["--input-root", str(root), "--output-dir", str(out)])
    assert (out / "old.png").read_bytes() == b"old"


def test_a_failed_plate_run_drops_the_old_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """New plates must not sit beside the previous run's commit and inputs (#621)."""
    root = _root(tmp_path)
    _same_classifier(root)
    out = tmp_path / "plates"
    out.mkdir()
    (out / "plates_provenance.json").write_text("{}")

    def fail(_repo: Path, _out: Path) -> Path:
        raise RuntimeError("plate failed")

    monkeypatch.setattr(plates, "build_methods_plate1", fail)
    with pytest.raises(RuntimeError, match="plate failed"):
        plates.main(["--input-root", str(root), "--output-dir", str(out), "--overwrite"])
    assert not (out / "plates_provenance.json").exists()


def _oxygen_rows() -> list[dict[str, str]]:
    rows = []
    for index, (score, rank) in enumerate(
        [(0.9, 1.0), (0.5, 9.0), (0.25, 5.0), (0.21, 5.0), (0.19, 50.0), (-0.16, 2.0)]
        + [(-0.9 + i / 100, 1.0) for i in range(12)]
        + [(0.3 + i / 100, float(i)) for i in range(12)]
    ):
        rows.append(
            {
                "condition": f"c{index}",
                "body_site_display": "Feces" if index % 2 else "Saliva",
                "aerobe_up_anaerobe_down_score": str(score),
                "rank_score": str(rank),
            }
        )
    return rows


def test_slide_eight_and_plate_five_select_the_same_strata() -> None:
    """Slide 8's text is worded from the strata Plate V draws (#639)."""
    rows = _oxygen_rows()
    up, down = plates.plate5_strata(rows)
    assert deck.oxygen_switch_strata(rows) == up + down
    assert len(up) == 10
    assert len(down) == 10


def test_plate_p_values_do_not_round_across_the_threshold() -> None:
    assert plates.fmt_p(0.0496) == "0.0496"
    assert plates.fmt_p(0.04996) == "0.04996"
    assert plates.fmt_p(0.0561) == "0.056"
    assert plates.fmt_p(0.000167) == "1.7e-04"


def _plate_one_tables(halophily_selected: bool) -> tuple[Any, ...]:
    pooled = {
        tid: {
            "q_value": "0.001",
            "study_q_value": "0.002",
            "pair_count": "2335",
            "study_count": "959",
        }
        for tid, *_ in plates.TABLE1_TRAITS
    }
    feces = {
        tid: {"study_q_value": "0.01", "study_count": "594"} for tid, *_ in plates.TABLE1_TRAITS
    }
    sign = {
        label: {
            "bh_q": "0.001",
            "family_size": "6",
            "studies_up": "500",
            "studies_down": "400",
            "studies_tied": "39",
        }
        for _, label, *_ in plates.TABLE1_TRAITS
    }
    cell = {
        "view_key": "individual_conditions_direction_contrast",
        "passes_support": "True",
        "increased_profile_count": "40",
        "decreased_profile_count": "50",
        "body_site_group": "Feces",
    }
    contrasts = [
        {
            **cell,
            "trait_label": "cellulose catabolic process",
            "condition": condition,
            "fisher_q_value": "0.008",
            "selected_after_bh": "True",
        }
        for condition in ("Parkinson's disease", "HIV infection")
    ] + [
        {
            **cell,
            "trait_label": "moderately halophilic",
            "condition": "Parkinson's disease",
            "fisher_q_value": "0.026",
            "selected_after_bh": str(halophily_selected),
        },
        {
            **cell,
            "trait_label": "coccus shaped",
            "condition": "Aging",
            "passes_support": "False",
            "fisher_q_value": "0.007",
            "selected_after_bh": "False",
        },
    ]
    wilcoxon = [
        {
            "trait": label,
            "condition": "Parkinson's disease",
            "q_value": q,
            "paired_study_count": n,
            "selected_after_bh": "False",
        }
        for (_, label, *_), q, n in zip(
            plates.TABLE1_TRAITS,
            ["0.6", "0.5", "0.4", "0.3", "0.2", "0.0561"],
            ["3", "45", "9", "9", "9", "45"],
            strict=True,
        )
    ]
    return pooled, feces, sign, contrasts, wilcoxon


def test_plate_one_counts_survivors_from_its_tables() -> None:
    """Plate I claimed survivors its own tables do not have (#636)."""
    summary = plates.plate1_summary(*_plate_one_tables(halophily_selected=True))
    assert summary["survivors"] == [6, 6, 6, 6, 2, 0]
    assert summary["m"] == [6, 6, 6, 6, 4, 6]
    assert summary["wilcoxon_studies"] == (3, 45)
    # An unsupported cell's q is not plotted as the trait's best Scale 5 result.
    assert summary["points"]["coccus shaped"][4] == 1.0
    assert plates.plate1_subtitle(summary) == (
        "All 6 pass Scales 1\N{EN DASH}4; 2 pass Scale 5 (cellulose: HIV, PD; halophily: PD); "
        "none passes Scale 6 (lowest q = 0.056, halophily in PD)"
    )
    unselected = plates.plate1_summary(*_plate_one_tables(halophily_selected=False))
    assert unselected["survivors"][4] == 1


def test_plate_one_slide_text_reads_the_tables() -> None:
    pooled = [
        {"trait": tid, "study_q_value": "0.01" if tid != "METPO:1000683" else "0.2"}
        for tid in deck.TABLE1_TRAIT_IDS
    ]
    metrics = {
        "paired_contrasts": {"pooled": pooled},
        "pd_halo_wilcoxon": {"p_value": "0.000167", "q_value": "0.0561"},
        "ind_sig_count": 1,
        "pooled_pairs": 2335,
        "pooled_studies": 959,
        "feces_studies": 602,
        "ind_total_tests": 3144,
        "pd_net_total": 408,
        "ind_conditions": 356,
        "ind_sites": 20,
        "ind_median_studies": 1,
    }
    takeaway, caveat = deck.plate_one_text(metrics)
    assert "602 fecal studies" in takeaway
    assert "408 fecal study Wilcoxon tests" in takeaway
    assert caveat.startswith("5 of 6 traits pass study-level BH")
    assert "(median 1 study per cell), 1 individual profile cell is selected after BH" in caveat
    assert "p = 0.00017 (BH q = 0.056)" in caveat


def test_prior_link_says_when_the_study_is_in_the_corpus() -> None:
    studies = {"28536926": {"publication_year": "2017", "doi": "10.1007/s11427-016-9001-4"}}
    url = "https://pubmed.ncbi.nlm.nih.gov/28536926/"
    assert deck.prior_label(url, studies) == (
        "PubMed 28536926 (2017) · DOI 10.1007/s11427-016-9001-4 · a study in the analysed "
        "corpus (bsdb:28536926)"
    )
    assert deck.prior_label("https://pubmed.ncbi.nlm.nih.gov/1/", studies) == "PubMed 1"


def test_crc_oxygen_split_reads_the_table() -> None:
    rows = [
        {
            "condition": "Colorectal cancer",
            "body_site_display": "Feces",
            "aerobe_up_anaerobe_down_score": "0.683",
        },
        {
            "condition": "Colorectal cancer",
            "body_site_display": "Lower GI tract",
            "aerobe_up_anaerobe_down_score": "-0.793",
        },
    ]
    assert "(Feces +0.68 vs Lower GI tract -0.79)" in deck.crc_oxygen_text(rows)
    assert plates.crc_oxygen_split(rows) == (
        "Feces +0.68 (aerobe↑) vs Lower GI tract -0.79 (anaerobe↑)"
    )
    assert deck.crc_oxygen_text(rows[:1]) == ""
    assert plates.crc_oxygen_split(rows[:1]) == ""


def test_lollipop_title_names_its_profiles() -> None:
    summary = [
        {"increased_profiles": "130", "decreased_profiles": "123"},
        {"increased_profiles": "130", "decreased_profiles": "123"},
    ]
    assert deck.pd_profiles_text({"pd_six_trait_summary": summary}) == (
        "130 increased vs 123 decreased fecal profiles"
    )


def test_slides_are_checked_link_by_link(tmp_path: Path) -> None:
    """Figure, evidence, code and methods links are each checked at their own revision;
    quoted tables at COMMIT (#637)."""
    good = tmp_path / "trait-visualization" / "good.png"
    good.parent.mkdir()
    good.write_text("good")
    stale = tmp_path / "trait-visualization" / "stale.tsv"
    stale.write_text("now")
    tree = {
        "trait-visualization/good.png": deck.git_blob_id(good),
        "trait-visualization/stale.tsv": "0" * 40,
    }
    trees = _trees(tmp_path, tree, deck.COMMIT)
    entry = {
        "figure_file_path": "trait-visualization/good.png",
        "evidence_path": "trait-visualization/stale.tsv",
        "code_path": "scripts/missing.py",
        "methods_path": "methods/missing.md",
    }
    unverified, quoted = deck.check_slides(tmp_path, [entry], trees, [stale])
    assert unverified == {
        "trait-visualization/good.png": [
            ("trait-visualization/stale.tsv", f"differs from {deck.SHORT_COMMIT}"),
            ("scripts/missing.py", f"not committed at {deck.SHORT_COMMIT}"),
            ("methods/missing.md", f"not committed at {deck.SHORT_COMMIT}"),
        ]
    }
    # A code or methods link only has to open: no local copy is compared.
    linked_only = {**entry, "code_path": "trait-visualization/stale.tsv"}
    found, _ = deck.check_slides(tmp_path, [{**linked_only, "evidence_path": ""}], trees, [])
    assert ("trait-visualization/stale.tsv", f"differs from {deck.SHORT_COMMIT}") not in found.get(
        "trait-visualization/good.png", []
    )
    assert quoted == [f"trait-visualization/stale.tsv: differs from {deck.SHORT_COMMIT}"]
    links = deck.slide_links(entry, entry["figure_file_path"])
    checked = deck.verified_links(links, unverified["trait-visualization/good.png"], entry)
    assert checked["figure_url"] == links["figure_url"]
    assert checked["figure_url"] != "NA"
    assert checked["evidence_url"] == "NA"
    assert checked["code_url"] == ""
    assert checked["methods_url"] == ""


def test_a_slide_not_rerun_is_checked_by_its_figure(tmp_path: Path) -> None:
    """Its figure must be the pull request's committed file; the methods it quotes are
    only linked there, since this checkout holds main's version (#614)."""
    figure = tmp_path / "methods" / "figures" / "table1_traits_focus_oral.png"
    figure.parent.mkdir(parents=True)
    figure.write_text("pr173")
    notes = tmp_path / "methods" / "notes.md"
    notes.write_text("main's version")
    ref = deck.UNMERGED_TABLE1_MAPS[1]
    trees = _trees(
        tmp_path,
        {
            "methods/figures/table1_traits_focus_oral.png": deck.git_blob_id(figure),
            "methods/notes.md": "0" * 40,
            "src/bsdbng/table1_trait_maps.py": "1" * 40,
        },
        ref,
    )
    entry: dict[str, Any] = {
        "figure_file_path": "methods/figures/table1_traits_focus_oral.png",
        "code_path": "src/bsdbng/table1_trait_maps.py",
        "methods_path": "methods/notes.md",
        "evidence_path": "methods/notes.md",
    }
    deck.mark_unmerged(entry)
    unverified, _ = deck.check_slides(tmp_path, [entry], trees, [])
    assert unverified == {}
    figure.write_text("main's map")
    unverified, _ = deck.check_slides(tmp_path, [entry], trees, [])
    assert unverified == {
        "methods/figures/table1_traits_focus_oral.png": [
            ("methods/figures/table1_traits_focus_oral.png", f"differs from {ref[:12]}")
        ]
    }


def test_files_that_differ_stop_the_build_unless_it_is_a_review() -> None:
    assert deck.refuse_or_review([], allow_unverified=False) == []
    with pytest.raises(SystemExit, match="1 files differ"):
        deck.refuse_or_review(["x.png: differs from 402ef5c25b71"], allow_unverified=False)
    assert deck.refuse_or_review(["x.png: differs"], allow_unverified=True) == [
        "1 files differ from the revision their slides cite"
    ]


def test_stale_results_are_named(tmp_path: Path) -> None:
    doc = {
        "figures": [{"figure_id": "fig_slide_004", "figure_file_path": "a.png"}],
        "statistical_results": [{"result_id": "res_1", "figure_ids": ["fig_slide_004"]}],
    }
    deck.check_results_cite_the_deck(doc, {"a.png"}, tmp_path / "p.yaml")
    with pytest.raises(ValueError, match=r"res_1.*fig_slide_004.*p\.yaml"):
        deck.check_results_cite_the_deck(doc, set(), tmp_path / "p.yaml")


def test_a_missing_curated_figure_stops_the_build(tmp_path: Path) -> None:
    """A missing plate dropped its slide and moved its results to the next one (#638)."""
    (tmp_path / "a.png").write_text("a")
    items = [{"figure_file_path": "a.png"}, {"figure_file_path": "a.png"}]
    included: set[str] = set()
    placed = deck.placed_slides(items, tmp_path, included, "plate")
    assert [item["figure_file_path"] for item in placed] == ["a.png"]
    assert included == {"a.png"}
    with pytest.raises(FileNotFoundError, match="missing plate"):
        deck.placed_slides([{"figure_file_path": "b.png"}], tmp_path, set(), "plate")


def test_cellulose_figures_must_be_the_ones_drawn(tmp_path: Path) -> None:
    figure = f"{deck.CELLULOSE_DIR}/figures/fig1.svg"
    table = tmp_path / deck.CELLULOSE_DIR / "t.tsv"
    (tmp_path / deck.CELLULOSE_DIR / "figures").mkdir(parents=True)
    (tmp_path / figure).write_text("drawn")
    captions = tmp_path / deck.CELLULOSE_DIR / "figures" / "figure_captions.json"
    captions.write_text("{}")
    table.write_text("t")
    record = {
        "inputs": {"t.tsv": hashlib.sha256(b"t").hexdigest()},
        "outputs": {
            "figures/fig1.svg": hashlib.sha256(b"drawn").hexdigest(),
            "figures/figure_captions.json": hashlib.sha256(b"{}").hexdigest(),
        },
    }
    trees = _trees(tmp_path, {f"{deck.CELLULOSE_DIR}/t.tsv": deck.git_blob_id(table)}, COMMIT)
    assert deck.cellulose_problems(record, [figure], trees, COMMIT, tmp_path) == []
    (tmp_path / figure).write_text("redrawn")
    captions.write_text('{"edited": true}')
    assert deck.cellulose_problems(record, [figure], trees, COMMIT, tmp_path) == [
        f"{figure}: not the file that was drawn",
        f"{deck.CELLULOSE_DIR}/figures/figure_captions.json: not the file that was drawn",
    ]
    old = {**record, "outputs": ["figures/fig1.svg"]}
    assert deck.cellulose_problems(old, [figure], trees, COMMIT, tmp_path) == [
        "the figures' provenance records no output hashes"
    ]


def test_carried_over_hashes_are_refreshed_in_one_format(tmp_path: Path) -> None:
    """Records mixed "path:sha256" with bare hashes, and nothing checked them (#637)."""
    first = tmp_path / "trait-visualization" / "a.tsv"
    first.parent.mkdir()
    first.write_text("a")
    (tmp_path / "trait-visualization" / "b.tsv").write_text("b")
    sha_a = hashlib.sha256(b"a").hexdigest()
    files = ["trait-visualization/a.tsv", "trait-visualization/b.tsv", "trait-visualization/c.tsv"]
    recorded = [f"trait-visualization/a.tsv:{sha_a}", "0" * 64, "1" * 64]
    current, changed = deck.refreshed_hashes(tmp_path, files, recorded)
    assert current == [sha_a, hashlib.sha256(b"b").hexdigest(), "1" * 64]
    assert changed == ["trait-visualization/b.tsv", "trait-visualization/c.tsv"]


def test_slide_three_quotes_the_counts_methods_plate_two_draws() -> None:
    selected = [
        {"all3_concordant": True, "a1_st_q": 0.01, "a3_bh_q": 0.01, "a3_by_q": 0.01},
        {"all3_concordant": True, "a1_st_q": 0.01, "a3_bh_q": 0.01, "a3_by_q": 0.2},
        {"all3_concordant": False, "a1_st_q": 0.01, "a3_bh_q": 0.01, "a3_by_q": 0.01},
        # A Table 1 trait drawn on Plate VII without passing both arms
        {"all3_concordant": True, "a1_st_q": 0.3, "a3_bh_q": 0.01, "a3_by_q": 0.01},
    ]
    assert plates.three_arm_counts(selected) == (2, 1)


def test_slide_eight_says_whether_the_mucosa_stratum_is_shown() -> None:
    rows = _oxygen_rows()
    assert deck.crc_mucosa_shown({"oxygen_switch": rows}) == (
        "the second is not among the strata shown"
    )
    rows.append(
        {
            "condition": "Colorectal cancer",
            "body_site_display": "Intestinal / GI mucosa",
            "aerobe_up_anaerobe_down_score": "-2.0",
            "rank_score": "1",
        }
    )
    assert deck.crc_mucosa_shown({"oxygen_switch": rows}) == "the second is among the strata shown"


def test_plate_one_subtitle_grammar_for_one_survivor() -> None:
    summary = plates.plate1_summary(*_plate_one_tables(halophily_selected=False))
    assert "1 passes Scale 5 (cellulose: HIV, PD)" in plates.plate1_subtitle(summary)


# ------------------------------------------------------------------------ Table 1 text (#579)

FA_CEL = ("facultatively_aerobic", "cellulose_catabolic_process")
CEL_SPH = ("cellulose_catabolic_process", "sphere_shaped")
CEL_ARS = ("cellulose_catabolic_process", "arsenate_reduction")
CEL_HAL = ("cellulose_catabolic_process", "moderately_halophilic")
ALL_PASSING = {FA_CEL, CEL_SPH, CEL_ARS}


def _comparison(
    population: str, studies: int, passing: set[tuple[str, str]], cellulose: tuple[int, int]
) -> list[dict[str, str]]:
    """disease_control_comparison.tsv rows: cellulose shifts down, the other five up, and
    the ``passing`` log-ratios at BY q = 0.01 (the rest 0.3). 20 studies have a zero
    annotation total; the rest are scored or tied in every fraction test."""
    up, down = cellulose
    tested = studies - 20
    rows = []
    for trait in CODES:
        shifted = (up, down) if trait == CELLULOSE_TRAIT else (60, 40)
        rows.append(
            {
                "population": population,
                "comparison": "fraction",
                "numerator": trait,
                "reference": "",
                "population_studies": str(studies),
                "studies_scored": str(sum(shifted)),
                "studies_up": str(shifted[0]),
                "studies_down": str(shifted[1]),
                "studies_tied": str(tested - sum(shifted)),
                "p_value": "4.1302e-06" if trait == CELLULOSE_TRAIT else "0.001",
                "q_value": "0.002",
            }
        )
    for numerator, reference in itertools.combinations(CODES, 2):
        rows.append(
            {
                "population": population,
                "comparison": "log_ratio",
                "numerator": numerator,
                "reference": reference,
                "population_studies": str(studies),
                "studies_scored": "237",
                "studies_up": "96",
                "studies_down": "141",
                "studies_tied": "0",
                "p_value": "0.004",
                "q_value": "0.01" if (numerator, reference) in passing else "0.3",
            }
        )
    return rows


def _signs(oriented: set[tuple[str, str]]) -> dict[str, Any]:
    return populations(
        _comparison(ALL_STUDIES, 959, ALL_PASSING, (333, 441))
        + _comparison(ORIENTED, 592, oriented, (187, 288))
    )


def test_ratio_text_names_what_the_oriented_subset_changes() -> None:
    same = deck.ratio_survivors_text({"table1": _signs(ALL_PASSING)})
    assert same == (
        "3 pairs pass 15-test BY q ≤ 0.05 across all 959 studies (FA/CEL, CEL/SPH and "
        "CEL/ARS) and 3 in the 592 oriented studies"
    )
    added = deck.ratio_survivors_text({"table1": _signs(ALL_PASSING | {CEL_HAL})})
    assert added.endswith(", adding CEL/HAL (96↑:141↓ [40.5%], BY q = 0.010)")
    lost = deck.ratio_survivors_text({"table1": _signs({FA_CEL, CEL_SPH})})
    assert lost.endswith("and 2 in the 592 oriented studies; CEL/ARS does not pass there")


def test_increased_side_caveat_counts_the_pairs_among_them() -> None:
    caveat = deck.increased_side_caveat({"table1": _signs(ALL_PASSING)})
    assert caveat.startswith("None of the 10 pairs among the 5 increased-side traits passes")
    coccus = ("facultatively_aerobic", "coccus_shaped")
    caveat = deck.increased_side_caveat({"table1": _signs(ALL_PASSING | {coccus})})
    assert caveat.startswith("Among the 5 increased-side traits, FA/COC passes BY q ≤ 0.05")


def test_cover_quotes_the_comparison_table() -> None:
    text = deck.cover_table1_text({"table1": _signs(ALL_PASSING)})
    assert text == (
        "Arm 3 study-level Table 1: all 6 traits pass BH q &lt; 0.01 across the 939 studies with "
        "a nonzero annotation total (in the 572 such disease-vs-healthy oriented studies, the "
        "cellulose annotation share is larger in the decreased signature in 288 of 475 scored "
        "studies (60.6%; 97 tied), p = 4.1 \N{MULTIPLICATION SIGN} 10⁻⁶, and 3 log-ratios pass "
        "BY q ≤ 0.05: FA/CEL, CEL/SPH and CEL/ARS)"
    )


def test_oriented_cellulose_text_reads_both_populations() -> None:
    assert deck.oriented_cellulose_text({"table1": _signs(ALL_PASSING)}) == (
        "the cellulose catabolism annotation share is larger in the decreased signature in "
        "57.0% of scored studies (BH q = 0.0020) and in 60.6% of the oriented ones (p = "
        "4.1e-06, BH q = 0.0020), BH over the six Table 1 traits"
    )


def test_direction_shift_slide_quotes_its_map_and_the_oriented_subset() -> None:
    """The typed 283/465 matched no version of the table (#606)."""
    shift = {
        trait: {
            "trait": trait,
            "studies_up": "532" if trait == "facultatively_aerobic" else "333",
            "studies_down": "383" if trait == "facultatively_aerobic" else "441",
            "studies_scored": "915" if trait == "facultatively_aerobic" else "774",
            "share_up": "0.581421" if trait == "facultatively_aerobic" else "0.430233",
            "bh_q": "0.03" if trait == "coccus_shaped" else "5.66164e-06",
        }
        for trait in CODES
    }
    metrics = {
        "table1": _signs(ALL_PASSING),
        "table1_shift": shift,
        "contrast_classifier": "contrast-orientation-v4",
    }
    title, takeaway, caveat = deck.direction_shift_text(metrics)
    assert title.startswith("6 of the 6 Table 1 traits shift in one direction")
    assert "532/915 scored studies (58.1%, BH q = 5.7e-06)" in takeaway
    assert "decreased signature in 441/774 (43.0%↑ / 57.0%↓" in takeaway
    oriented = (
        "572 disease-vs-healthy oriented studies with a nonzero annotation total, not drawn "
        "here, the cellulose share is larger in the decreased signature in 288 of 475 scored "
        "studies (60.6%; 97 tied)"
    )
    assert f"{oriented}, p = 4.1e-06, BH q = 0.0020." in takeaway
    assert "contrast-orientation-v4" in caveat
    assert "317" not in takeaway + caveat


def test_sci_writes_out_powers_of_ten() -> None:
    assert deck.sci(4.1302e-06) == "4.1 \N{MULTIPLICATION SIGN} 10⁻⁶"
    assert deck.sci(0.0518) == "0.052"


def test_plate_seven_slide_reads_the_plate_counts() -> None:
    counts = {
        "table1_concordant": 6,
        "screened_traits": 317,
        "three_arm_additional": 11,
        "three_arm_concordant": 17,
        "three_arm_strict_by": 8,
        "ora_inverted": 2,
        "not_concordant": 0,
    }
    assert deck.three_arm_text({"plate_counts": counts}) == (
        "6 of the 6 original Table 1 traits (★) pass all 3 arms with concordant direction; "
        "screening all 317 traits finds 11 more 3-arm-concordant traits (17 in all at study BH "
        "q ≤ 0.05; 8 at strict BY q ≤ 0.05), 2 traits whose paired and sign-test directions "
        "agree while ORA points the other way."
    )
    counts["not_concordant"] = 1
    assert deck.three_arm_text({"plate_counts": counts}).endswith(
        "the other way and 1 trait whose arms disagree otherwise."
    )
    # Plates drawn before the counts were recorded say so rather than print a number
    assert deck.three_arm_text({"plate_counts": {}}).startswith("? of the 6 original")


def test_plates_count_plate_seven_groups() -> None:
    def trait(in_t1: bool, concordant: bool, group: str) -> dict[str, Any]:
        return {
            "in_t1": in_t1,
            "all3_concordant": concordant,
            "a1_st_q": 0.01,
            "a3_bh_q": 0.01,
            "group": group,
        }

    selected = [trait(True, True, "table1")] * 5 + [trait(True, False, "table1")]
    selected += [trait(False, True, "concordant")] * 3
    selected += [trait(False, False, "ora_inverted")] * 2
    selected += [trait(False, False, "not_concordant")]
    assert plates.plate_seven_groups(selected) == {
        "table1_traits": 6,
        "table1_concordant": 5,
        "three_arm_additional": 3,
        "ora_inverted": 2,
        "not_concordant": 1,
        "forest_traits": 12,
    }


def test_plate_three_labels_follow_the_table() -> None:
    signs = _signs(ALL_PASSING)
    every, oriented = signs[ALL_STUDIES], signs[ORIENTED]
    cellulose = [(CELLULOSE_TRAIT, trait) for trait in plates.increased_side(every)]
    increased = list(itertools.combinations(plates.increased_side(every), 2))
    assert plates.oriented_cellulose_pairs_text(every, oriented, cellulose) == (
        "Oriented disease-vs-healthy filter: the same 3 of 5 Cellulose pairs pass BY; nearest "
        "other CEL/COC BY q = 0.30"
    )
    assert plates.increased_pairs_text(every, oriented, increased) == (
        "all BY q ≥ 0.30 in both all & oriented"
    )
    signs = _signs(ALL_PASSING | {CEL_HAL})
    assert plates.oriented_cellulose_pairs_text(signs[ALL_STUDIES], signs[ORIENTED], cellulose) == (
        "Oriented disease-vs-healthy filter adds CEL/HAL (4 of 5 Cellulose pairs pass BY)"
    )
    assert plates.survivor_codes(signs[ORIENTED]) == "FA/CEL, CEL/SPH, CEL/ARS, CEL/HAL"


# ----------------------------------------------------------------------- links (#614, #661)


def test_only_the_direction_shift_map_moves_to_main() -> None:
    """Its text is rebuilt from main's table; main draws the other maps without PR #173's
    outlines and triangles, which their text describes (#614)."""
    assert list(deck.MOVED_FIGURES) == ["methods/figures/table1_traits_direction_shift.png"]
    row = {"figure_file_path": "methods/figures/table1_traits_direction_shift.png"}
    deck.moved_figure(row)
    assert row["figure_file_path"] == "data/derived/table1_traits/table1_traits_direction_shift.png"
    assert row["evidence_path"].endswith("table1_traits_direction_shift.tsv")
    for linked in (row["code_path"], row["methods_path"]):
        assert (REPO / linked).exists(), linked


def test_pr173_maps_are_shown_and_checked_at_its_head() -> None:
    entry: dict[str, Any] = {
        "figure_file_path": "methods/figures/table1_traits_focus_oral.png",
        "code_path": "src/bsdbng/table1_trait_maps.py",
        "methods_path": "methods/table1_trait_condition_body_site_maps.md",
        "evidence_path": "methods/table1_trait_condition_body_site_maps.md",
    }
    deck.mark_unmerged(entry)
    name, ref = deck.UNMERGED_TABLE1_MAPS
    assert entry["link_ref"] == entry["code_ref"] == ref
    links = deck.slide_links(entry, entry["figure_file_path"])
    for key in ("figure_url", "methods_url", "evidence_url", "code_url"):
        assert f"/blob/{ref}/" in links[key], key
    row = {"commit": links["commit"], **entry}
    assert deck.revision_stamp(row) == f"Not rerun on main · code in {name} at {ref[:12]}"
    assert ref in deck.cited_refs()
    assert deck.drawing_command(entry["code_path"], {}, entry["figure_file_path"]) == (
        f"# not rerun on main: drawn by src/bsdbng/table1_trait_maps.py as {name} has it, at {ref}"
    )
    # Main's own map from the shim of the same name is rerun
    main_map: dict[str, Any] = {
        "figure_file_path": "data/derived/table1_traits/table1_traits_focus_oral.png",
        "code_path": "src/bsdbng/table1_trait_maps.py",
    }
    deck.mark_unmerged(main_map)
    assert "not_rerun" not in main_map
    # The feces-vs-other slide shows PR #173's map too, and its record follows it
    assert deck.RENAMED_FIGURES[
        "data/derived/table1_traits/table1_traits_focus_feces_other.png"
    ] == (deck.PR173_FECES_OTHER)


def test_relinked_code_and_evidence_are_on_main() -> None:
    """The never-committed create_trait_Heatmaps.py and create_oral_six_trait_heatmap.py
    are replaced by the scripts main runs, with the tables they write (#614)."""
    for linked in (*deck.TRAIT_HEATMAP_CODE.values(), deck.ORAL_HEATMAP_CODE):
        assert (REPO / linked).exists(), linked
    assert (REPO / deck.ORAL_HEATMAP_EVIDENCE).exists()
    for figure in deck.TRAIT_HEATMAP_CODE:
        if "/bsdb-taxa-only/" in figure:
            assert (REPO / figure.replace(".png", "_matrix.tsv")).exists(), figure


def test_unmerged_code_is_linked_at_its_pull_request() -> None:
    entry: dict[str, Any] = {
        "figure_file_path": "data/derived/x/gram_positive_inc_sankey.png",
        "code_path": "src/bsdbng/gram_sankey.py",
    }
    deck.mark_unmerged(entry)
    name, ref = deck.UNMERGED_GRAM_SANKEY
    assert entry["code_ref"] == ref
    assert entry["not_rerun"] == f"{name} at {ref[:12]}"
    assert deck.slide_links(entry, entry["figure_file_path"])["code_url"].endswith(
        f"/blob/{ref}/src/bsdbng/gram_sankey.py"
    )
    row = {"figure_file_path": entry["figure_file_path"], "commit": deck.COMMIT, **entry}
    assert deck.revision_stamp(row) == f"Not rerun on main · code in {name} at {ref[:12]}"
    assert {ref, deck.UNMERGED_ORA_FIGURES[1]} <= set(deck.cited_refs())
    merged: dict[str, Any] = {"code_path": deck.ORAL_HEATMAP_CODE}
    deck.mark_unmerged(merged)
    assert "code_ref" not in merged


def test_unmerged_code_is_not_on_main() -> None:
    """Code merged into main is linked and rerun there, not at a pull request."""
    for code in deck.UNMERGED_CODE:
        assert not (REPO / code).exists(), code


def test_a_slide_shows_the_figure_it_records(tmp_path: Path) -> None:
    """A static svg_sources/ copy would show an older image than the PNG (#661)."""
    (tmp_path / "assets").mkdir()
    (tmp_path / "assets_svg").mkdir()
    png = tmp_path / "figure.png"
    png.write_bytes(b"png")
    assert deck.copy_asset(png, tmp_path, 2) == "assets/v6_slide_002_figure.png"
    png.with_suffix(".svg").write_text("<svg/>")
    assert deck.copy_asset(png, tmp_path, 2) == "assets_svg/v6_slide_002_figure.svg"
    assert (tmp_path / "assets/v6_slide_002_figure.png").read_bytes() == b"png"


# ----------------------------------------------------------------- provenance (#615, #622)


def _recipes() -> dict[str, str]:
    recipes: dict[str, str] = deck.justfile_recipes((REPO / "Justfile").read_text("utf-8"))
    return recipes


def test_pipeline_steps_name_real_recipes_and_code() -> None:
    recipes = _recipes()
    seen: set[str] = set()
    for step in deck.PIPELINE_STEPS:
        assert set(step["recipes"]) <= set(recipes), step["step_id"]
        assert (REPO / step["script_or_module"]).exists(), step["script_or_module"]
        assert set(step["prerequisite_step_ids"]) <= seen, step["step_id"]
        seen.add(step["step_id"])
    steps = deck.pipeline_steps()
    assert [step["step_order"] for step in steps] == list(range(1, len(steps) + 1))
    assert steps[3]["command_line"] == "just download-kg && just taxa-trait-matrix"
    assert "recipes" not in steps[0]


def test_recipes_are_found_by_the_code_they_run() -> None:
    by_code = deck.recipes_by_code(_recipes())
    assert by_code["scripts/create_table1_trait_maps.py"] == ["table1-trait-maps"]
    assert by_code["src/mtp/trait_screen.py"] == ["trait-direction-screen"]
    assert by_code[deck.ORAL_HEATMAP_CODE] == [
        "six-shared-trait-heatmaps",
        "oral-six-trait-heatmap",
    ]


def test_command_series_drop_what_no_longer_exists() -> None:
    recipes = _recipes()
    by_code = deck.recipes_by_code(recipes)
    series = [
        *deck.UPSTREAM_COMMANDS,
        "just trait-profile-investigation  # uv run python scripts/create_trait_Heatmaps.py",
        "uv run python scripts/create_trait_Heatmaps.py",
        "uv run python -m mtp.cli analyze --mode enrichment",
        "just no-such-recipe",
    ]
    assert deck.figure_commands(
        series, "scripts/create_taxa_guild_heatmap.py", recipes, by_code
    ) == [*deck.UPSTREAM_COMMANDS, "just trait-heatmaps"]
    name, ref = deck.UNMERGED_ORA_FIGURES
    assert deck.figure_commands(
        series, "scripts/plot_enrichment_histograms.py", recipes, by_code
    ) == [f"# not rerun on main: scripts/plot_enrichment_histograms.py is only in {name}, at {ref}"]


def test_command_series_do_not_repeat_the_drawing_step() -> None:
    recipes = _recipes()
    series = [*deck.UPSTREAM_COMMANDS, "just oral-six-trait-heatmap  # old note"]
    assert deck.figure_commands(
        series, deck.ORAL_HEATMAP_CODE, recipes, deck.recipes_by_code(recipes)
    ) == [
        *deck.UPSTREAM_COMMANDS,
        "just six-shared-trait-heatmaps && just oral-six-trait-heatmap",
    ]


def test_figure_records_name_their_code_and_existing_steps() -> None:
    recipes = _recipes()
    by_code = deck.recipes_by_code(recipes)
    row = {
        "code_path": "scripts/create_taxa_guild_heatmap.py",
        "figure_file_path": "trait-visualization/taxa_guild_discovery_clustered_heatmap.png",
        "evidence_path": "",
        "methods_path": "",
    }
    rec: dict[str, Any] = {"pipeline_step_ids": ["step_01_download", "step_99_gone"]}
    deck.drawn_provenance(rec, row, recipes, by_code)
    assert rec["generator_script_path"] == row["code_path"]
    assert rec["just_recipe"] == "just trait-heatmaps"
    assert rec["pipeline_step_ids"] == ["step_01_download"]
    unmerged = {**row, "code_path": "src/bsdbng/gram_sankey.py"}
    rec = {"pipeline_step_ids": ["step_01_download"]}
    deck.drawn_provenance(rec, unmerged, recipes, by_code)
    assert rec["just_recipe"].startswith("# not rerun on main")
    assert rec["pipeline_step_ids"] == []


def test_only_recipes_that_run_without_arguments_are_steps() -> None:
    """A command series names steps that run as written (#622): `core-screen run output`
    needs arguments, so it is not one, while defaults and *ARGS are fine."""
    text = (
        'deck_slides := "data/derived/x"\n'
        'set shell := ["bash", "-cu"]\n'
        "\n"
        "# Screen\n"
        "trait-direction-screen *ARGS:\n"
        "\tuv run python -m mtp.trait_screen --code-commit x $*\n"
        "\n"
        "core-screen run output:\n"
        '\tuv run python -m mtp.trait_screen --core-run "{{run}}"\n'
        "\n"
        "[positional-arguments]\n"
        'workflow-check spec="workflows/core.yaml":\n'
        '\tuv run python -m mtp.workflow check "{{spec}}"\n'
        "\n"
        "@quiet:\n"
        "\tuv run python scripts/quiet.py\n"
    )
    recipes = deck.justfile_recipes(text)
    assert set(recipes) == {"trait-direction-screen", "workflow-check", "quiet"}
    assert deck.recipes_by_code(recipes) == {
        "src/mtp/trait_screen.py": ["trait-direction-screen"],
        "src/mtp/workflow.py": ["workflow-check"],
        "scripts/quiet.py": ["quiet"],
    }


def test_exploratory_heatmap_slide_reads_its_manifest() -> None:
    """The September text typed 23 literature matches; after #584 there are 25 (#579)."""
    manifest = {
        "displayed_rows": 30,
        "taxonomic_ranks": ["species", "genus", "family"],
        "displayed_rows_with_literature_reference": 25,
        "literature_bug_sig_direction_agreement_counts": {"agrees": 22, "differs": 3},
        "cells": 180,
        "status_counts": {"mixed": 161, "low support": 12, "dec": 6, "inc": 1},
    }
    taxa = {
        "Escherichia coli": {"increased_profiles": "1", "decreased_profiles": "2"},
        "Escherichia": {"increased_profiles": "9", "decreased_profiles": "2"},
    }
    title, takeaway, caveat = deck.exploratory_heatmap_text(
        {"pd_exploratory": {"manifest": manifest, "taxa": taxa}}
    )
    assert title.startswith("30 PD taxa")
    assert takeaway == (
        "30 taxa across 3 ranks; 25 have a PD literature direction, which agrees with the "
        "BugSigDB direction for 22 and differs for 3. Recurrence: E. coli 1 increased / 2 "
        "decreased; Escherichia 9 increased / 2 decreased."
    )
    assert caveat.startswith("Influence is exploratory, non-additive and non-causal: 7 of 180")
    assert "161 are mixed and 12 low support" in caveat


def test_plate_four_arrows_follow_the_means() -> None:
    assert plates.lead_arrow({"study_weighted_mean_contrast": "-0.677"}) == "↓"
    assert plates.lead_arrow({"study_weighted_mean_contrast": "1.136"}) == "↑"
    assert plates.lead_arrow({"study_weighted_mean_contrast": "NA"}) == ""
    assert plates.lead_arrow({}) == ""


def test_appendix_skips_byte_identical_copies(tmp_path: Path) -> None:
    """Copies of one heatmap under several names became separate slides (#586)."""
    folder = tmp_path / "trait-visualization" / "trait-profile-investigation" / "taxa-coverage"
    folder.mkdir(parents=True)
    for name in ("a_heatmap.png", "b_heatmap.png"):
        (folder / name).write_bytes(b"same image")
    (folder / "c_heatmap.png").write_bytes(b"another image")
    slides = deck.get_remaining_appendix_slides(set(), tmp_path, {})
    assert [Path(slide["figure_file_path"]).name for slide in slides] == [
        "a_heatmap.png",
        "c_heatmap.png",
    ]


# ------------------------------------------------------------- review of #665 (#677-#686)


def test_kept_commands_run_in_pipeline_order() -> None:
    """trait-profile-investigation reads what trait-heatmaps writes (#681)."""
    recipes = _recipes()
    series = ["just trait-profile-investigation", "just analyze", "just trait-heatmaps"]
    by_code = deck.recipes_by_code(recipes)
    assert deck.figure_commands(series, "scripts/create_guild_maps.py", recipes, by_code)[
        -4:-1
    ] == [
        "just trait-heatmaps",
        "just trait-profile-investigation",
        "just analyze",
    ]
    assert deck.step_order("just no-such-recipe") == len(deck.PIPELINE_STEPS)


def test_plate_two_rows_keep_their_group_sign() -> None:
    pooled = {
        "a": {"study_mean_increased_minus_decreased": "0.2"},
        "b": {"study_mean_increased_minus_decreased": "-0.1"},
    }
    plates.check_plate_two_directions([("a", "A", "Redox"), ("b", "B", "Plant glycan ↓")], pooled)
    with pytest.raises(ValueError, match=r"other sign: \['b in Redox at -0.100'\]"):
        plates.check_plate_two_directions([("b", "B", "Redox")], pooled)


def test_pd_split_reads_the_switch_table() -> None:
    """The typed 23 and 32 matched only 402ef5c (#680)."""

    def genus(name: str, inc: str, dec: str) -> dict[str, str]:
        return {
            "taxon_name": name,
            "top_increased_conditions": inc,
            "top_decreased_conditions": dec,
        }

    pd = "Parkinson's disease ({} profiles, 15 studies); Obesity (3 profiles, 2 studies)"
    rows = [
        genus("Bifidobacterium", pd.format(26), "Eczema (8 profiles, 2 studies)"),
        genus("Roseburia", "Glaucoma (5 profiles, 2 studies)", pd.format(33)),
        genus("Lactobacillus", pd.format(16), ""),
        genus("Blautia", "", "COVID-19 (9 profiles, 7 studies); " + pd.format(17)),
        genus("Prevotella", pd.format(4), pd.format(2)),
    ]
    assert plates.pd_split(rows) == (
        "Bifidobacterium (26↑) & Lactobacillus (16↑) vs Roseburia (33↓) & Blautia (17↓)"
    )
    assert plates.pd_split(rows[:1]) == "Bifidobacterium (26↑) vs none"


def test_plate_three_pairs_follow_what_panel_b_draws() -> None:
    """Cellulose on the increased side made ratio(CEL, CEL) raise KeyError (#686)."""
    rows = _comparison(ALL_STUDIES, 959, ALL_PASSING, (441, 333)) + _comparison(
        ORIENTED, 592, ALL_PASSING, (288, 187)
    )
    signs = populations(rows)
    every, oriented = signs[ALL_STUDIES], signs[ORIENTED]
    assert CELLULOSE_TRAIT in plates.increased_side(every)
    others = [trait for trait in every.fractions if trait != CELLULOSE_TRAIT]
    pairs = [(CELLULOSE_TRAIT, trait) for trait in others]
    assert plates.oriented_cellulose_pairs_text(every, oriented, pairs).startswith(
        "Oriented disease-vs-healthy filter: the same 3 of 5 Cellulose pairs pass BY"
    )
    assert plates.increased_pairs_text(every, oriented, []) == "no pairs"
    coccus = ("facultatively_aerobic", "coccus_shaped")
    one = populations(
        _comparison(ALL_STUDIES, 959, ALL_PASSING | {coccus}, (333, 441))
        + _comparison(ORIENTED, 592, ALL_PASSING, (187, 288))
    )
    assert plates.increased_pairs_text(one[ALL_STUDIES], one[ORIENTED], [coccus]) == (
        "FA/COC passes BY in all or oriented"
    )


def test_deck_text_counts_one_trait_and_one_pair() -> None:
    signs = _signs(ALL_PASSING)
    every = signs[ALL_STUDIES]
    fractions = dict(every.fractions)
    for trait in list(fractions)[1:]:
        fractions[trait] = replace(fractions[trait], q_value=0.2)
    single = {**signs, ALL_STUDIES: replace(every, fractions=fractions)}
    assert "1 of the 6 traits passes BH" in deck.cover_table1_text({"table1": single})
    # Two traits on the increased side (FA, COC) leave one pair; none leaves no pair to test
    for trait in list(fractions)[3:]:
        fractions[trait] = replace(fractions[trait], studies_up=10, studies_down=90)
    two = {**signs, ALL_STUDIES: replace(every, fractions=fractions)}
    assert deck.increased_side_caveat({"table1": two}).startswith(
        "The one pair of increased-side traits does not pass BY q ≤ 0.05"
    )
    for trait in list(fractions)[:3]:
        fractions[trait] = replace(fractions[trait], studies_up=10, studies_down=90)
    none = {**signs, ALL_STUDIES: replace(every, fractions=fractions)}
    assert deck.increased_side_caveat({"table1": none}) == (
        "0 Table 1 traits are on the increased side, so no pair among them is tested."
    )


def test_plates_refuse_mixed_contrast_classifiers(tmp_path: Path) -> None:
    """The cellulose_decr notes and the Table 1 text must define disease vs healthy alike
    (#679)."""
    root = _root(tmp_path)
    _same_classifier(root)
    plates.check_classifiers(REPO, root)
    (root / plates.CONTRAST_PROVENANCE).write_text(json.dumps({"classifier_version": "v0"}))
    with pytest.raises(SystemExit, match="Table 1 tables v0"):
        plates.check_classifiers(REPO, root)


def test_appendix_refuses_a_figure_it_cannot_credit(tmp_path: Path) -> None:
    """Any data/derived PNG used to be credited to PR #315's histogram script (#686)."""
    histogram = tmp_path / "data/derived/trait_enrichment_x_histogram.png"
    histogram.parent.mkdir(parents=True)
    histogram.write_bytes(b"histogram")
    (slide,) = deck.get_remaining_appendix_slides(set(), tmp_path, {})
    assert slide["code_path"] == "scripts/plot_enrichment_histograms.py"
    assert "September rerun by PR #315" in slide["takeaway"]
    unknown = tmp_path / "data/derived/new_figure.png"
    unknown.write_bytes(b"new")
    with pytest.raises(ValueError, match=r"no generator is known for data/derived/new_figure\.png"):
        deck.get_remaining_appendix_slides(set(), tmp_path, {})


def test_the_svg_a_slide_shows_is_checked_like_its_png(tmp_path: Path) -> None:
    """The deck certified the PNG but showed whatever SVG sat beside it (#684)."""
    png = tmp_path / "trait-visualization" / "map.png"
    png.parent.mkdir()
    png.write_text("png")
    svg = png.with_suffix(".svg")
    svg.write_text("committed svg")
    tree = {
        "trait-visualization/map.png": deck.git_blob_id(png),
        "trait-visualization/map.svg": deck.git_blob_id(svg),
    }
    trees = _trees(tmp_path, tree, deck.COMMIT)
    entry = {"figure_file_path": "trait-visualization/map.png"}
    assert deck.check_slides(tmp_path, [entry], trees, [])[0] == {}
    svg.write_text("an older run's svg")
    unverified, _ = deck.check_slides(tmp_path, [entry], trees, [])
    assert unverified == {
        "trait-visualization/map.png": [
            ("trait-visualization/map.svg", f"differs from {deck.SHORT_COMMIT}")
        ]
    }
