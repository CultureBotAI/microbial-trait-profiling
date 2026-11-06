"""Historical CELLULOSE assets cannot re-enter through old deck manifests (#1572)."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[1]
CELLULOSE = "data/analysis/cellulose_decr"
HISTORICAL_STEMS = (
    "fig1_tufte_convergence_slope_and_overlap",
    "fig2_tufte_core_guild_dot_table",
    "fig3_tufte_why_only_stool",
    "fig4_tufte_trait_landscape_and_four_pairs",
    "fig5_tufte_three_statistical_arms",
    "fig6_tufte_combining_results_integrated_evidence",
)
HISTORICAL_OUTPUTS = (
    *(
        f"{CELLULOSE}/figures_tufte/{stem}.{ext}"
        for stem in HISTORICAL_STEMS
        for ext in ("svg", "png", "pdf")
    ),
    f"{CELLULOSE}/figures_tufte/fig7_tufte_convergence_network_and_venn.svg",
    f"{CELLULOSE}/figures_tufte/fig8_literature_support_four_pairs.svg",
    f"{CELLULOSE}/figures_tufte/fig9_set_theory_17_core_vs_pairwise.svg",
    f"{CELLULOSE}/figures_tufte/cellulose_decr_tufte_figures.pdf",
    f"{CELLULOSE}/cellulose_decr_slide_deck.html",
)
PRESERVED = (
    f"{CELLULOSE}/figures/fig1_pairwise_and_3way_jaccard_nulls.svg",
    f"{CELLULOSE}/figures_graphical/fig1_overlap_and_null_models.svg",
    "data/derived/tufte_plates/tufte_plate_1.png",
    "data/derived/tufte_plates/fig1_tufte_convergence_slope_and_overlap.svg",
    "data/analysis/other/figures_tufte/fig1_tufte_convergence_slope_and_overlap.svg",
    "trait-visualization/trait-profile-investigation/PD-cellulose/pd_cellulose_enrichment_test.png",
)


@pytest.fixture
def deck() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "deck_retirement_fixture", REPO / "scripts" / "build_slide_deck.py"
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SelectionCaptured(Exception):
    """Stop after real entry selection, before provenance checks or asset rendering."""


def selected_inherited_entries(
    deck: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    version: str,
    section: str,
    paths: tuple[str, ...],
    checkout: Path | None = None,
) -> list[dict[str, str]]:
    """Inject actual TSV inputs while stubbing work outside manifest selection."""
    checkout = checkout or tmp_path
    monkeypatch.setattr(deck, "REPO_ROOT", checkout)
    monkeypatch.setattr(
        deck,
        "load_rerun_metrics",
        lambda _repo: {"contrast_classifier": "fixture", "matrix_taxa": 1, "trait_count": 1},
    )
    monkeypatch.setattr(deck, "load_cellulose_captions", lambda _repo: {"statements": {}})
    monkeypatch.setattr(deck, "cellulose_revision", lambda _repo: "fixture")
    monkeypatch.setattr(deck, "plate_revision", lambda _repo: "fixture")
    monkeypatch.setattr(deck, "record_problems", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(deck, "extract_committed_figures", lambda _path: None)
    monkeypatch.setattr(deck, "get_new_slides_metadata", lambda _metrics: [])
    monkeypatch.setattr(deck, "cellulose_decr_slides", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(deck, "get_remaining_appendix_slides", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(
        deck,
        "update_previous_slide_text",
        lambda row, _metrics: (row["title"], row["takeaway"], row["caveat"]),
    )
    for relative in (
        "data/derived/tufte_plates/plates_provenance.json",
        f"{CELLULOSE}/figures/figures_provenance.json",
        f"{CELLULOSE}/cellulose_decr_provenance.json",
    ):
        path = checkout / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("{}", encoding="utf-8")

    # Existing user artifacts must be excluded too; this does not render them.
    for relative in HISTORICAL_OUTPUTS:
        path = checkout / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("legacy artifact fixture", encoding="utf-8")

    previous = tmp_path / deck.DECK_DIR / "slides"
    previous.mkdir(parents=True)
    fields = (
        "slide",
        "figure_file_path",
        "deck_section",
        "code_path",
        "title",
        "takeaway",
        "caveat",
    )
    for manifest_version in ("v1", "v2"):
        with (previous / f"slide_manifest_{manifest_version}.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            if manifest_version == version:
                for index, figure_path in enumerate(paths, start=1):
                    writer.writerow(
                        {
                            "slide": str(index),
                            "figure_file_path": figure_path,
                            "deck_section": section,
                            "code_path": "scripts/source_fixture.py",
                            "title": figure_path,
                            "takeaway": "fixture takeaway",
                            "caveat": "fixture caveat",
                        }
                    )

    selected: list[dict[str, str]] = []

    def capture(_metrics: dict[str, Any], entries: list[dict[str, str]], **_kwargs: Any) -> None:
        selected.extend(dict(entry) for entry in entries)
        raise SelectionCaptured

    monkeypatch.setattr(deck, "check_coverage_table_sync", capture)
    with pytest.raises(SelectionCaptured):
        deck.build_deck(tmp_path, tmp_path / "out", tmp_path / "provenance", trees={})
    return selected


@pytest.mark.parametrize("version", ["v1", "v2"])
@pytest.mark.parametrize("section", ["Main", "Appendix"])
def test_inherited_manifests_exclude_retired_cellulose_outputs(
    deck: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    version: str,
    section: str,
) -> None:
    aliases = tuple(
        alias
        for path in HISTORICAL_OUTPUTS
        for alias in (
            f"./{path}",
            path.replace("data/analysis/", "data//analysis/./"),
            path.replace("cellulose_decr/", "cellulose_decr/figures/../", 1),
        )
    )
    paths = (*PRESERVED[:2], *HISTORICAL_OUTPUTS, *aliases, *PRESERVED[2:])
    entries = selected_inherited_entries(deck, tmp_path, monkeypatch, version, section, paths)

    assert [entry["figure_file_path"] for entry in entries] == list(PRESERVED)
    assert [entry["deck_section"] for entry in entries] == [section] * len(PRESERVED)
    assert [entry["title"] for entry in entries] == list(PRESERVED)
    for entry in entries:
        assert entry[f"{version}_slide"] != "NA"
        assert entry[f"{'v2' if version == 'v1' else 'v1'}_slide"] == "NA"


@pytest.mark.parametrize("version", ["v1", "v2"])
@pytest.mark.parametrize("section", ["Main", "Appendix"])
def test_absolute_retired_assets_are_excluded_from_both_roots(
    deck: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    version: str,
    section: str,
) -> None:
    input_root = tmp_path / "input"
    checkout = tmp_path / "checkout"
    retired = tuple(
        str(root / path) for root in (input_root, checkout) for path in HISTORICAL_OUTPUTS
    )
    preserved = tuple(str(input_root / path) for path in PRESERVED)
    entries = selected_inherited_entries(
        deck,
        input_root,
        monkeypatch,
        version,
        section,
        (*preserved, *retired),
        checkout=checkout,
    )
    assert [entry["figure_file_path"] for entry in entries] == list(preserved)


@pytest.mark.parametrize("version", ["v1", "v2"])
@pytest.mark.parametrize("section", ["Main", "Appendix"])
def test_symlink_aliases_to_retired_assets_are_excluded_from_both_roots(
    deck: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    version: str,
    section: str,
) -> None:
    """Resolve file and directory links while retaining links to current figures."""
    input_root = tmp_path / "input"
    checkout = tmp_path / "checkout"
    retired_aliases: list[str] = []
    preserved_aliases: list[str] = []
    for root in (input_root, checkout):
        for relative in HISTORICAL_OUTPUTS:
            target = root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text("legacy artifact fixture", encoding="utf-8")
        directory_alias = root / "retired-directory-alias"
        directory_alias.symlink_to(root / CELLULOSE, target_is_directory=True)
        for relative in HISTORICAL_OUTPUTS:
            file_alias = root / "file-aliases" / Path(relative).name
            file_alias.parent.mkdir(parents=True, exist_ok=True)
            file_alias.symlink_to(root / relative)
            retired_aliases.extend(
                (
                    str(file_alias),
                    str(directory_alias / Path(relative).relative_to(CELLULOSE)),
                )
            )
        for index, relative in enumerate((PRESERVED[0], PRESERVED[2])):
            target = root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text("current figure fixture", encoding="utf-8")
            alias = root / f"current-figure-alias-{index}{target.suffix}"
            alias.symlink_to(target)
            preserved_aliases.append(str(alias))

    preserved = (*PRESERVED, *preserved_aliases)
    entries = selected_inherited_entries(
        deck,
        input_root,
        monkeypatch,
        version,
        section,
        (*retired_aliases, *preserved),
        checkout=checkout,
    )
    assert [entry["figure_file_path"] for entry in entries] == list(preserved)


@pytest.mark.parametrize("version", ["v1", "v2"])
@pytest.mark.parametrize("section", ["Main", "Appendix"])
def test_retirement_filter_preserves_paths_outside_exact_root_and_name_scope(
    deck: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    version: str,
    section: str,
) -> None:
    """Pin the existing two-root, 23-name, case-sensitive policy without widening it."""
    input_root = tmp_path / "input"
    checkout = tmp_path / "checkout"
    third_root = tmp_path / "other-checkout"
    preserved = (
        *(str(third_root / relative) for relative in HISTORICAL_OUTPUTS),
        *(
            str(root / relative)
            for root in (input_root, checkout)
            for relative in (
                HISTORICAL_OUTPUTS[0].replace("cellulose_decr/", "Cellulose_decr/"),
                HISTORICAL_OUTPUTS[0].replace(".svg", ".SVG"),
                f"{CELLULOSE}/figures_tufte/unlisted_figure.svg",
            )
        ),
    )
    assert len(HISTORICAL_OUTPUTS) == 23
    assert frozenset(HISTORICAL_OUTPUTS) == deck.RETIRED_CELLULOSE_OUTPUTS
    entries = selected_inherited_entries(
        deck,
        input_root,
        monkeypatch,
        version,
        section,
        preserved,
        checkout=checkout,
    )
    assert [entry["figure_file_path"] for entry in entries] == list(preserved)
