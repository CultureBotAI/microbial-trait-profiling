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
    deck.plate_provenance(record)
    assert record["generator_script_path"] == "scripts/generate_tufte_plates.py"
    assert record["just_recipe"] == "just figure-review-plates"
    assert record["command_series"][-3:] == [
        "just analyze",
        "just cellulose-decr-figures",
        "just figure-review-plates",
    ]
    assert "data/derived/figure_review" not in json.dumps(record)


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


@pytest.mark.parametrize("argv", [["--over"], ["--out", "/tmp/x"], ["--prov", "/tmp/x"]])
def test_deck_options_are_never_abbreviated(argv: list[str]) -> None:
    """The recipe's guards match whole option names (#617)."""
    with pytest.raises(SystemExit) as stopped:
        deck.main(argv)
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


def test_plates_record_their_drawing_commit(tmp_path: Path) -> None:
    path = plates.write_plates_provenance(
        tmp_path, [tmp_path / "b.png", tmp_path / "a.png"], COMMIT
    )
    assert json.loads(path.read_text()) == {
        "script": "scripts/generate_tufte_plates.py",
        "code_commit": COMMIT,
        "plates": ["a.png", "b.png"],
    }


def test_plate_and_figure_options_are_never_abbreviated() -> None:
    """The recipes' guards match whole option names (#617)."""
    from mtp.cellulose_decr import figures

    with pytest.raises(SystemExit) as stopped:
        figures.parse_args(["--analysis", "/tmp/x"])
    assert stopped.value.code == 2
    with pytest.raises(SystemExit) as stopped:
        plates.main(["--code", "abc"])
    assert stopped.value.code == 2


def test_plates_need_the_cellulose_captions(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="just cellulose-decr-figures"):
        plates.cellulose_decr_notes(_cellulose_repo(tmp_path, None))
