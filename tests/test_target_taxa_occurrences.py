"""The v2 scripts keep human-host target-taxon occurrences and note the rest (#587, #712, #742)."""

import csv
import importlib.util
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))


def load_script(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(f"occurrences_{name}", SCRIPTS_DIR / f"{name}.py")
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


aerobe = load_script("investigate_facultative_aerobe_expansion")
anaerobe = load_script("investigate_gut_anaerobe_depletion")
# Each script's focus direction and the score column its profile metadata carries.
SCRIPTS = {
    "aerobe": (aerobe, "increased", "oxygen_tolerance_z_mean"),
    "anaerobe": (anaerobe, "decreased", "anaerobe_z_mean"),
}
HUMAN_PROFILES = ["bsdb:1|bsdb:1-1|inc", "bsdb:1|bsdb:1-1|dec"]


def write_inputs(root: Path, taxon_ids: list[str]) -> None:
    """A human experiment with one occurrence per direction; two mouse experiments,
    one with an increased and one with two decreased occurrences; and a human PCR
    experiment with an increased occurrence, outside the 16S and shotgun corpus (#742)."""
    first, second = taxon_ids[:2]
    derived = root / "data" / "derived"
    raw = root / "data" / "raw"
    derived.mkdir(parents=True)
    raw.mkdir(parents=True)
    taxa = pd.DataFrame(
        [
            ("bsdb:1", "bsdb:1-1", "bsdb:1-1-1", "increased", first),
            ("bsdb:1", "bsdb:1-1", "bsdb:1-1-2", "decreased", second),
            ("bsdb:2", "bsdb:2-1", "bsdb:2-1-1", "increased", first),
            ("bsdb:3", "bsdb:3-1", "bsdb:3-1-1", "decreased", first),
            ("bsdb:3", "bsdb:3-1", "bsdb:3-1-1", "decreased", second),
            ("bsdb:4", "bsdb:4-1", "bsdb:4-1-1", "increased", first),
        ],
        columns=["study_id", "experiment_id", "signature_id", "direction", "taxon_id"],
    )
    taxa["taxon_name"] = "taxon"
    taxa["taxonomic_rank"] = "species"
    taxa.to_csv(derived / "bugsigdb_taxa.tsv", sep="\t", index=False)
    with (raw / "full_dump.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "BSDB ID",
                "Host species",
                "Condition",
                "Body site",
                "Location of subjects",
                "Sequencing type",
            ]
        )
        writer.writerow(
            ["bsdb:1/1/1", "Homo sapiens", "Colorectal cancer", "Feces", "China", "16S"]
        )
        writer.writerow(["bsdb:2/1/1", "Mus musculus", "Obesity", "Feces", "China", "16S"])
        writer.writerow(["bsdb:3/1/1", "Mus musculus", "Obesity", "Feces", "China", "WMS"])
        writer.writerow(["bsdb:4/1/1", "Homo sapiens", "Obesity", "Feces", "China", "PCR"])


def scored_profiles(score_column: str) -> pd.DataFrame:
    """Human-host profile metadata for the one human experiment."""
    return pd.DataFrame(
        {
            "study_id": "bsdb:1",
            "experiment_id": "bsdb:1-1",
            "inc_dec": ["inc", "dec"],
            "direction_label": ["increased", "decreased"],
            "title": "A study",
            "publication_year": 2020,
            "condition": "Colorectal cancer",
            "body_site": "Feces",
            "body_site_group": "Feces",
            "condition_group": "Colorectal cancer",
            "location": "China",
            "host_species": "Homo sapiens",
            "guild_z_mean": [0.5, -0.5],
            score_column: [0.5, -0.5],
        },
        index=pd.Index(HUMAN_PROFILES, name="profile_id"),
    )


def load_occurrences(
    script: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[pd.DataFrame, pd.DataFrame]:
    module, _, score_column = SCRIPTS[script]
    write_inputs(tmp_path, list(module.TARGET_TAXA))
    monkeypatch.setattr(module, "REPO_ROOT", tmp_path)
    occurrences, excluded = module.load_target_taxa_occurrences(scored_profiles(score_column))
    assert isinstance(occurrences, pd.DataFrame)
    assert isinstance(excluded, pd.DataFrame)
    return occurrences, excluded


@pytest.mark.parametrize("script", SCRIPTS)
def test_loader_keeps_human_host_occurrences_and_sets_aside_the_rest(
    script: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    occurrences, excluded = load_occurrences(script, tmp_path, monkeypatch)
    assert sorted(occurrences["profile_id"]) == sorted(HUMAN_PROFILES)
    assert occurrences["condition"].eq("Colorectal cancer").all()
    reasons = excluded.set_index("study_id")["exclusion_reason"]
    assert reasons.to_dict() == {
        "bsdb:2": "host",
        "bsdb:3": "host",
        "bsdb:4": "sequencing_type",
    }
    assert excluded.loc[excluded["study_id"] == "bsdb:4", "sequencing_type"].tolist() == ["PCR"]
    assert excluded["direction"].value_counts().to_dict() == {"decreased": 2, "increased": 2}


@pytest.mark.parametrize("script", SCRIPTS)
def test_condition_charts_note_the_set_aside_focus_direction(
    script: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module, direction, _ = SCRIPTS[script]
    occurrences, excluded = load_occurrences(script, tmp_path, monkeypatch)
    calls: dict[str, dict[str, Any]] = {}

    def record(name: str) -> Any:
        def save(*args: Any, **kwargs: Any) -> None:
            calls[name] = {"args": args, **kwargs}

        return save

    monkeypatch.setattr(module, "save_condition_breakdown_bar_chart", record("breakdown"))
    monkeypatch.setattr(module, "save_condition_target_taxa_stacked_bar_chart", record("stacked"))
    module.save_target_taxa_condition_charts(occurrences, excluded, pd.DataFrame())

    expected = {
        "increased": "2 increased-profile occurrences (2 profiles) from non-human hosts and "
        "other sequencing methods are",
        "decreased": "2 decreased-profile occurrences (1 profile) from non-human hosts are",
    }[direction]
    note = f"Human-host 16S and shotgun profiles only: {expected} excluded."
    assert calls["breakdown"]["note"] == note
    assert calls["stacked"]["note"] == note
    counts = calls["breakdown"]["args"][0]
    assert counts.to_dict() == {"Colorectal cancer": 1}


@pytest.mark.parametrize("script", SCRIPTS)
def test_anatomical_sites_diagram_orders_labels_and_clears_count_bars(
    script: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Caller diagrams order labels by leader_order_key and clear the count bars (#718)."""
    module, _, _ = SCRIPTS[script]
    summary = pd.DataFrame(
        {
            "body_site_group": ["Dental", "Throat"],
            "profile_count": [10, 10],
            "signal": [0.1, -0.1],
        }
    )
    real_close = module.plt.close
    monkeypatch.setattr(module.plt, "close", lambda *_a, **_k: None)
    try:
        module.save_anatomical_sites_diagram(
            summary,
            output_path=tmp_path / f"{script}_sites.png",
            title="Anatomical Sites",
            signal_column="signal",
            signal_label="z",
        )
        fig = module.plt.gcf()
        fig.canvas.draw()
        ax = fig.axes[0]
        label_y = {text.get_text(): text.get_position()[1] for text in ax.texts}
        assert label_y["Throat"] > label_y["Dental"]
        bars = [line for line in ax.lines if line.get_zorder() == 3]
        assert len(bars) == 2
        bar_right = max(
            ax.transData.transform((max(line.get_xdata()), line.get_ydata()[0]))[0]
            + 0.5 * line.get_linewidth() * fig.dpi / 72.0
            for line in bars
        )
        cbar_ax = ax.child_axes[0]
        assert cbar_ax.get_window_extent().x0 > bar_right
    finally:
        real_close("all")


def test_anatomical_enriched_diagram_orders_labels_by_leader_key(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The anaerobe enriched-profile diagram also orders labels by leader_order_key (#718)."""
    summary = pd.DataFrame(
        {
            "body_site_group": ["Subgingival dental plaque", "Throat"],
            "enriched_profile_count": [10, 10],
            "increased_enriched_profile_count": [4, 4],
            "decreased_enriched_profile_count": [6, 6],
        }
    )
    real_close = anaerobe.plt.close
    monkeypatch.setattr(anaerobe.plt, "close", lambda *_a, **_k: None)
    try:
        anaerobe.save_anatomical_enriched_profile_diagram(
            summary,
            {"inc": 1, "dec": 2},
            output_path=tmp_path / "enriched_sites.png",
        )
        fig = anaerobe.plt.gcf()
        fig.canvas.draw()
        ax = fig.axes[0]
        label_y = {text.get_text(): text.get_position()[1] for text in ax.texts}
        assert label_y["Throat"] > label_y["Subgingival dental plaque"]
    finally:
        real_close("all")
