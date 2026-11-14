"""Influence table bytes must not depend on the interpreter hash seed (#2320)."""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

from mtp.contribution_inputs import Base, Graph
from mtp.taxon_contributions import reconcile
from mtp.taxon_influence import run_influence


def write_fixture(output: Path) -> None:
    """Exercise real set-valued source rollup with unequal resolved genus counts."""
    traits = [f"METPO:fixture{index}" for index in range(6)]
    family = "NCBITaxon:10000"
    parents: dict[str, set[str]] = {}
    annotations: dict[str, set[str]] = {}
    ranks = {family: "family"}
    for genus_index in range(1, 12):
        genus = f"NCBITaxon:{1000 + genus_index}"
        parents[genus] = {family}
        ranks[genus] = "genus"
        for species_index in range(genus_index):
            species = f"NCBITaxon:{100000 + genus_index * 100 + species_index}"
            parents[species] = {genus}
            ranks[species] = "species"
            annotations[species] = {
                trait
                for trait_index, trait in enumerate(traits)
                if (species_index + genus_index + trait_index) % (trait_index + 2) != 0
            }
    totals: Counter[str] = Counter()
    for values in annotations.values():
        totals.update(values)
    profiles = {("study", "experiment", direction): dict(totals) for direction in ("inc", "dec")}
    base = Base(
        traits,
        dict(zip(traits, traits, strict=True)),
        profiles,
        {("study", "experiment"): {"condition": "condition", "body_site": "site"}},
        {key: {family: {"signature"}} for key in profiles},
        {family: dict(totals)},
        {},
        ranks,
        [],
    )
    graph = Graph(parents, annotations, {}, {}, False)
    ledger = reconcile(base, graph)
    assert not ledger.mismatches
    output.mkdir()
    run_influence(
        base,
        graph,
        ledger,
        {
            "minimum_lineage_studies": 1,
            "top_lineages_per_rank": 0,
            "include_trait_leading_genus": False,
            "predicted_genera": [],
            "scopes": [{"id": "pooled", "condition": "*", "body_site": "*"}],
        },
        output,
    )


@pytest.fixture(scope="module")
def seeded_outputs(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    root = tmp_path_factory.mktemp("hash-seeds")
    outputs = (root / "seed-1", root / "seed-2")
    for seed, output in zip((1, 2), outputs, strict=True):
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), str(output)],
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
            capture_output=True,
            text=True,
            check=False,
            timeout=30,
        )
        assert result.returncode == 0, result.stdout + result.stderr
    return outputs


@pytest.mark.parametrize(
    ("table", "column"),
    [("concentration.tsv", "hhi_resolved"), ("cross_trait_overlap.tsv", "weighted_overlap")],
)
def test_influence_tables_are_hash_seed_independent(
    seeded_outputs: tuple[Path, Path], table: str, column: str
) -> None:
    first, second = (root / table for root in seeded_outputs)
    with first.open(encoding="utf-8") as handle:
        left = list(csv.DictReader(handle, delimiter="\t"))
    with second.open(encoding="utf-8") as handle:
        right = list(csv.DictReader(handle, delimiter="\t"))
    assert left
    assert len(left) == len(right)
    assert any(row[column] for row in left)
    differences = [
        {"row": index, "column": key, "seed_1": a[key], "seed_2": b[key]}
        for index, (a, b) in enumerate(zip(left, right, strict=True), start=2)
        for key in a
        if a[key] != b[key]
    ]
    assert first.read_bytes() == second.read_bytes(), f"{table} differing cells: {differences}"


if __name__ == "__main__":
    write_fixture(Path(sys.argv[1]))
