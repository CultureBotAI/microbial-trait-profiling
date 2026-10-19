"""The selected-figure deck builders in scripts/ (#556, #573)."""

import hashlib
import importlib.util
import json
import shutil
from pathlib import Path
from types import ModuleType
from typing import Any

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
        ((-0.41, 0.04), (0.19, 0.2), "points the other way, against the figure's title"),
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
    assert record["command_series"][-3:] == [
        "just analyze",
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


def _deck_root(tmp_path: Path, version: str) -> Path:
    root = tmp_path / "root"
    provenance = root / deck.CONTRAST_PROVENANCE
    provenance.parent.mkdir(parents=True)
    provenance.write_text(json.dumps({"classifier_version": version}))
    return root


def _deck_args(tmp_path: Path, root: Path) -> list[str]:
    return [
        "--input-root",
        str(root),
        "--output-dir",
        str(tmp_path / "slides"),
        "--provenance-dir",
        str(tmp_path / "provenance"),
    ]


def test_deck_refuses_another_contrast_version(tmp_path: Path) -> None:
    root = _deck_root(tmp_path, "contrast-orientation-v3")
    with pytest.raises(SystemExit, match="#579"):
        deck.main(_deck_args(tmp_path, root))
    assert not (tmp_path / "slides").exists()


def test_deck_refuses_to_replace_a_deck(tmp_path: Path) -> None:
    root = _deck_root(tmp_path, deck.TYPED_CONTRAST_VERSION)
    existing = tmp_path / "slides" / "selected_figures_v6.html"
    existing.parent.mkdir()
    existing.write_text("old deck")
    with pytest.raises(SystemExit, match="--overwrite"):
        deck.main(_deck_args(tmp_path, root))
    assert existing.read_text() == "old deck"


def test_deck_needs_git_listings_or_a_review_build(tmp_path: Path) -> None:
    root = _deck_root(tmp_path, deck.TYPED_CONTRAST_VERSION)
    with pytest.raises(SystemExit, match="--allow-unverified"):
        deck.main(_deck_args(tmp_path, root))
    assert not (tmp_path / "slides").exists()


def _plates_root(tmp_path: Path, version: str) -> Path:
    root = tmp_path / "root"
    provenance = root / plates.CONTRAST_PROVENANCE
    provenance.parent.mkdir(parents=True)
    provenance.write_text(json.dumps({"classifier_version": version}))
    return root


def test_plates_refuse_another_contrast_version(tmp_path: Path) -> None:
    root = _plates_root(tmp_path, "contrast-orientation-v3")
    out = tmp_path / "plates"
    with pytest.raises(SystemExit, match="#579"):
        plates.main(["--input-root", str(root), "--output-dir", str(out)])
    assert not out.exists()


def test_plates_refuse_to_replace_plates(tmp_path: Path) -> None:
    root = _plates_root(tmp_path, plates.TYPED_CONTRAST_VERSION)
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
    root = _plates_root(tmp_path, plates.TYPED_CONTRAST_VERSION)
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
