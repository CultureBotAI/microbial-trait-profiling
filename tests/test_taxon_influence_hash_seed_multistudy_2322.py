"""Seeded reproducibility over every table written by run_influence (#2322)."""

from __future__ import annotations

import csv
import gzip
import io
import os
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

import mtp.taxon_influence as taxon_influence_mod
from mtp.contribution_inputs import Assertion, Base, Graph
from mtp.taxon_contributions import reconcile
from mtp.taxon_influence import run_influence

TOP_LEVEL_FILES = {
    "scenarios.tsv",
    "predictions.tsv",
    "contributing_universe.tsv",
    "concentration.tsv",
    "cross_trait_overlap.tsv",
    "effects.tsv.gz",
    "pooled_fraction_effects.tsv",
    "influence_summary.md",
}

SCENARIO_FILES = {
    "pair_eligibility.tsv.gz",
    "removed_source_counts.tsv.gz",
    "ambiguous_retained_sources.tsv",
    "profile_counts.tsv.gz",
    "ora_background_counts.tsv",
}

PRIMARY_SOURCES = (
    "infores:bacdive",
    "infores:ijsem",
    "infores:faprotax",
    "infores:madin",
)


def _read_table_bytes(path: Path) -> bytes:
    raw = path.read_bytes()
    return gzip.decompress(raw) if path.suffix == ".gz" else raw


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    text = _read_table_bytes(path).decode("utf-8")
    return list(csv.DictReader(io.StringIO(text, newline=""), delimiter="\t"))


def _build_multistudy_fixture() -> tuple[Base, Graph, dict[str, object]]:
    """Build a 5-study, 2-experiment-per-study fixture with unequal counts and assertions."""
    traits = [f"METPO:fixture{index}" for index in range(6)]
    families = ("NCBITaxon:20001", "NCBITaxon:20002", "NCBITaxon:20003")
    parents: dict[str, set[str]] = {}
    annotations: dict[str, set[str]] = {}
    assertions: dict[tuple[str, str], list[Assertion]] = {}
    ranks: dict[str, str] = {}
    names: dict[str, str] = {}
    genera: list[str] = []
    species_by_genus: dict[str, list[str]] = {}

    for fam_idx, family in enumerate(families):
        ranks[family] = "family"
        names[family] = f"Family_{fam_idx}"
        # Give family direct traits so ambiguous_retained_sources.tsv is non-empty on removal.
        fam_traits = {traits[(fam_idx + offset) % len(traits)] for offset in range(2)}
        annotations[family] = fam_traits
        for trait in sorted(fam_traits):
            assertions[(family, trait)] = [
                Assertion(
                    f"edge_{family}_{trait}",
                    "biolink:has_phenotype",
                    trait,
                    PRIMARY_SOURCES[fam_idx % len(PRIMARY_SOURCES)],
                    "infores:kg-microbe",
                    "knowledge_assertion",
                    10 + fam_idx,
                )
            ]
        for g_offset in range(4):
            genus_num = 3000 + fam_idx * 10 + g_offset
            genus = f"NCBITaxon:{genus_num}"
            genera.append(genus)
            parents[genus] = {family}
            ranks[genus] = "genus"
            names[genus] = (
                "GenusAlpha"
                if (fam_idx == 0 and g_offset == 0)
                else "GenusBeta"
                if (fam_idx == 1 and g_offset == 0)
                else f"Genus_{genus_num}"
            )
            species_by_genus[genus] = []
            # Unequal species counts across genera (1..5 species) for non-trivial HHI/mean/stdev.
            species_count = 1 + (fam_idx * 2 + g_offset) % 5
            for sp_idx in range(species_count):
                species = f"NCBITaxon:{genus_num * 100 + sp_idx}"
                species_by_genus[genus].append(species)
                parents[species] = {genus}
                ranks[species] = "species"
                names[species] = f"Species_{genus_num}_{sp_idx}"
                sp_traits = {
                    trait
                    for t_idx, trait in enumerate(traits)
                    if (sp_idx + g_offset + fam_idx + t_idx) % (t_idx % 3 + 2) != 0
                }
                if not sp_traits:
                    sp_traits = {traits[(sp_idx + g_offset) % len(traits)]}
                annotations[species] = sp_traits
                for t_idx, trait in enumerate(sorted(sp_traits)):
                    src_items: list[Assertion] = []
                    # Exercise multi-assertion, single-assertion, and missing-primary-source cases.
                    if (sp_idx + t_idx) % 4 != 0:
                        primary = PRIMARY_SOURCES[(sp_idx + t_idx + fam_idx) % len(PRIMARY_SOURCES)]
                        src_items.append(
                            Assertion(
                                f"edge_{species}_{t_idx}_a",
                                "biolink:has_phenotype",
                                trait,
                                primary,
                                "infores:kg-microbe",
                                "knowledge_assertion",
                                100 + sp_idx * 10 + t_idx,
                            )
                        )
                    if (sp_idx + g_offset) % 3 == 0:
                        secondary = PRIMARY_SOURCES[(t_idx + 1) % len(PRIMARY_SOURCES)]
                        src_items.append(
                            Assertion(
                                f"edge_{species}_{t_idx}_b",
                                "biolink:has_phenotype",
                                trait,
                                secondary,
                                "infores:kg-microbe",
                                "knowledge_assertion",
                                500 + sp_idx * 10 + t_idx,
                            )
                        )
                    if src_items:
                        assertions[(species, trait)] = src_items

    graph = Graph(parents, annotations, assertions, names, False)

    # 5 studies x 2 experiments = 10 pairs (20 profiles) with distinct unequal taxon mixes.
    studies = [f"study_{i}" for i in range(1, 6)]
    reported: dict[tuple[str, str, str], dict[str, set[str]]] = {}
    metadata: dict[tuple[str, str], dict[str, str]] = {}
    conditions = ("cond_a", "cond_b")
    body_sites = ("gut", "oral")

    all_species = [sp for genus in genera for sp in species_by_genus[genus]]
    for s_idx, study_id in enumerate(studies):
        for exp_idx in range(2):
            exp_id = f"exp_{exp_idx + 1}"
            pair = (study_id, exp_id)
            metadata[pair] = {
                "condition": conditions[(s_idx + exp_idx) % len(conditions)],
                "body_site": body_sites[(s_idx + exp_idx) % len(body_sites)],
            }
            for dir_idx, direction in enumerate(("inc", "dec")):
                key = (study_id, exp_id, direction)
                chosen: dict[str, set[str]] = {}
                # Ensure all 5 studies share lineages (minimum_lineage_studies=3)
                # while varying per-experiment and per-direction taxon sets and totals.
                step = 2 if direction == "inc" else 3
                start = (s_idx + exp_idx + dir_idx) % 3
                span = max(1, len(all_species) - 4)
                count_pick = 4 + ((s_idx * 3 + exp_idx * 2 + dir_idx) % span)
                for offset in range(count_pick):
                    taxon = all_species[(start + offset * step) % len(all_species)]
                    chosen[taxon] = {f"sig_{study_id}_{exp_id}_{direction}"}
                # Include a family directly in alternating runs to trigger ambiguous retention.
                if (s_idx + exp_idx) % 2 == 0:
                    fam = families[(s_idx + dir_idx) % len(families)]
                    chosen[fam] = {f"sig_{study_id}_{exp_id}_{direction}"}
                reported[key] = chosen

    targets = {taxon for taxa in reported.values() for taxon in taxa}
    taxa_counts: dict[str, dict[str, int]] = {target: {} for target in targets}
    reconstructed: dict[str, Counter[str]] = {target: Counter() for target in targets}
    for source, source_traits in graph.traits.items():
        for target in graph.ancestors(source) & targets:
            reconstructed[target].update(source_traits)
    for target, counter in reconstructed.items():
        taxa_counts[target] = {trait: counter[trait] for trait in traits}

    profiles: dict[tuple[str, str, str], dict[str, int]] = {}
    for key, taxon_map in reported.items():
        prof_counter: Counter[str] = Counter()
        for target in taxon_map:
            prof_counter.update(reconstructed[target])
        profiles[key] = {trait: prof_counter[trait] for trait in traits}

    base = Base(
        traits,
        dict(zip(traits, traits, strict=True)),
        profiles,
        metadata,
        reported,
        taxa_counts,
        names,
        ranks,
        [],
    )
    protocol: dict[str, object] = {
        "minimum_lineage_studies": 3,
        "top_lineages_per_rank": 1,
        "include_trait_leading_genus": True,
        "predicted_genera": ["GenusAlpha", "GenusBeta"],
        "scopes": [
            {"id": "pooled", "condition": "*", "body_site": "*"},
            {"id": "cond_a", "condition": "cond_a", "body_site": "*"},
        ],
    }
    return base, graph, protocol


def write_multistudy_fixture(output: Path) -> None:
    """Run ledger reconciliation and full run_influence on the multi-study fixture."""
    base, graph, protocol = _build_multistudy_fixture()
    ledger = reconcile(base, graph)
    assert not ledger.mismatches
    output.mkdir(parents=True, exist_ok=True)
    run_influence(base, graph, ledger, protocol, output)


@pytest.fixture(scope="module")
def seeded_outputs(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    root = tmp_path_factory.mktemp("hash-seeds-multistudy")
    outputs = (root / "seed-1", root / "seed-2")
    for seed, output in zip((1, 2), outputs, strict=True):
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), str(output)],
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
            capture_output=True,
            text=True,
            check=False,
            timeout=60,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        if result.stdout:
            print(f"seed={seed}: {result.stdout.strip()}", flush=True)
    return outputs


def test_every_run_influence_output_is_hash_seed_independent(
    seeded_outputs: tuple[Path, Path],
) -> None:
    first_root, second_root = seeded_outputs
    first_rel = {
        path.relative_to(first_root).as_posix() for path in first_root.rglob("*") if path.is_file()
    }
    second_rel = {
        path.relative_to(second_root).as_posix()
        for path in second_root.rglob("*")
        if path.is_file()
    }
    assert first_rel == second_rel

    scenario_rows = _read_tsv_rows(first_root / "scenarios.tsv")
    assert scenario_rows
    active_scenarios = [
        row["scenario"].replace(":", "_") for row in scenario_rows if row["status"] == "run"
    ]
    assert any(name.startswith("exclude_genus_") for name in active_scenarios)
    assert any(name.startswith("exclude_family_") for name in active_scenarios)
    scenario_dirs = sorted({"baseline", *active_scenarios})
    expected_rel = set(TOP_LEVEL_FILES) | {
        f"{scenario_dir}/{filename}"
        for scenario_dir in scenario_dirs
        for filename in SCENARIO_FILES
    }
    assert first_rel == expected_rel
    assert second_rel == expected_rel

    # Confirm primary_source_incidence assertions and ambiguous retention were genuinely exercised.
    universe_rows = _read_tsv_rows(first_root / "contributing_universe.tsv")
    primary_rows = [row for row in universe_rows if row["rank"] == "primary_source_incidence"]
    assert {row["taxon"] for row in primary_rows} & set(PRIMARY_SOURCES)
    assert any(row["taxon"] == "unresolved_source" for row in primary_rows)
    concentration_rows = _read_tsv_rows(first_root / "concentration.tsv")
    assert any(
        row["rank"] == "primary_source_incidence" and row["hhi_resolved"] != ""
        for row in concentration_rows
    )
    prediction_rows = _read_tsv_rows(first_root / "predictions.tsv")
    assert [row["prediction"] for row in prediction_rows] == ["GenusAlpha", "GenusBeta"]
    assert all(row["matched_lineages"] for row in prediction_rows)
    assert _read_tsv_rows(first_root / "cross_trait_overlap.tsv")
    assert _read_tsv_rows(first_root / "effects.tsv.gz")
    assert _read_tsv_rows(first_root / "pooled_fraction_effects.tsv")
    assert _read_table_bytes(first_root / "influence_summary.md")
    assert any(
        _read_tsv_rows(first_root / scenario_dir / "removed_source_counts.tsv.gz")
        for scenario_dir in scenario_dirs
        if scenario_dir != "baseline"
    )
    assert any(
        _read_tsv_rows(first_root / scenario_dir / "ambiguous_retained_sources.tsv")
        for scenario_dir in scenario_dirs
        if scenario_dir != "baseline"
    )

    byte_mismatches: list[str] = []
    for rel in sorted(first_rel):
        left_bytes = _read_table_bytes(first_root / rel)
        right_bytes = _read_table_bytes(second_root / rel)
        if left_bytes != right_bytes:
            if rel.endswith((".tsv", ".tsv.gz")):
                left_rows = _read_tsv_rows(first_root / rel)
                right_rows = _read_tsv_rows(second_root / rel)
                diffs = [
                    {"row": idx, "column": key, "seed_1": a[key], "seed_2": b[key]}
                    for idx, (a, b) in enumerate(zip(left_rows, right_rows, strict=True), start=2)
                    for key in a
                    if a[key] != b[key]
                ]
                byte_mismatches.append(f"{rel}: {diffs[:10]}")
            else:
                byte_mismatches.append(
                    f"{rel}: bytes differ ({len(left_bytes)} vs {len(right_bytes)})"
                )

    assert not byte_mismatches, "\n".join(byte_mismatches)


def test_omitted_active_scenario_directory_fails(
    seeded_outputs: tuple[Path, Path], tmp_path: Path
) -> None:
    """Omitting a complete status=run scenario directory in both seeds must fail."""
    first_src, second_src = seeded_outputs
    first_copy = tmp_path / "seed-1"
    second_copy = tmp_path / "seed-2"
    shutil.copytree(first_src, first_copy)
    shutil.copytree(second_src, second_copy)
    scenario_rows = _read_tsv_rows(first_copy / "scenarios.tsv")
    target_dir = next(
        row["scenario"].replace(":", "_") for row in scenario_rows if row["status"] == "run"
    )
    shutil.rmtree(first_copy / target_dir)
    shutil.rmtree(second_copy / target_dir)
    with pytest.raises(AssertionError):
        test_every_run_influence_output_is_hash_seed_independent((first_copy, second_copy))


def test_header_only_predictions_fails(seeded_outputs: tuple[Path, Path], tmp_path: Path) -> None:
    """Header-only predictions.tsv in both seeds must fail the population guard."""
    first_src, second_src = seeded_outputs
    first_copy = tmp_path / "seed-1"
    second_copy = tmp_path / "seed-2"
    shutil.copytree(first_src, first_copy)
    shutil.copytree(second_src, second_copy)
    for root in (first_copy, second_copy):
        pred_path = root / "predictions.tsv"
        header = pred_path.read_text(encoding="utf-8").splitlines()[0]
        pred_path.write_text(f"{header}\n", encoding="utf-8")
    with pytest.raises(AssertionError):
        test_every_run_influence_output_is_hash_seed_independent((first_copy, second_copy))


if __name__ == "__main__":
    print(f"TAXON_INFLUENCE_FILE={taxon_influence_mod.__file__}", flush=True)
    write_multistudy_fixture(Path(sys.argv[1]))
