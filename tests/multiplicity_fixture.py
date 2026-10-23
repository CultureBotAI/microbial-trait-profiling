"""Small synthetic native exports for #853; no production data or fitted claims."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pytest

from mtp.datamodel.workflow_pydantic import (
    AnalysisBasis,
    ArmCommandExecution,
    ArmRun,
    ComparisonRun,
    ConsumerRecord,
    ContrastPopulation,
    DirectionBasis,
    FileInventoryEntry,
    ParameterValue,
    RecordedPath,
    StorageKind,
    WorkflowRun,
)
from mtp.multiplicity import native, sources, stage
from mtp.multiplicity.sources import VerifiedSources
from mtp.paths import file_sha256
from mtp.trait_composition import compute_log_ratio_shifts, write_log_ratio_outputs
from mtp.trait_guilds import TABLE1_TRAITS, load_profiles
from mtp.trait_screen import candidate_slots, run_screen, write_screen_outputs
from mtp.workflow.coverage import feature_filter_id
from mtp.workflow.runtime import inventory


def make_basis(catalog: Path, features: tuple[str, ...], product: str, scope: str) -> AnalysisBasis:
    return AnalysisBasis(
        id=f"basis:{product}",
        feature_catalog_sha256=file_sha256(catalog),
        feature_filter_id=feature_filter_id(features),
        retained_feature_count=len(features),
        profile_unit=product,
        aggregation_policy="fixture",
        denominator_policy="fixture",
        background_inventory_id="fixture",
        weighting_policy="fixture",
        rollup_policy=f"ncbi_taxonomy_closure_rollup_{scope}",
    )


def write_table(
    path: Path, records: list[dict[str, Any]], comments: dict[str, str] | None = None
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        for key, value in (comments or {}).items():
            handle.write(f"# {key}: {value}\n")
        writer = csv.DictWriter(
            handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(records)


def file_entry(path: Path, *, relative_to: Path | None = None) -> FileInventoryEntry:
    return FileInventoryEntry(
        path=str(path.relative_to(relative_to)) if relative_to else str(path),
        size_bytes=path.stat().st_size,
        sha256=file_sha256(path),
    )


def file_snapshot(path: Path) -> RecordedPath:
    digest, size, members = inventory(path, "file")
    return RecordedPath(
        path=str(path),
        storage_kind=StorageKind.file,
        sha256=digest,
        size_bytes=size,
        inventory=members,
    )


def fixture_sources(
    root: Path, monkeypatch: pytest.MonkeyPatch, *, subset: bool = False
) -> VerifiedSources:
    """Stub only execution verification; exercise real importers, inventories and replay.

    Runtime verification has its own execution tests. These exports deliberately
    exercise an oriented subset and differing catalog sizes without manufacturing
    completed production manifests. Noel/Jed tables use actual native exporters.
    """
    core_dir, comparison_dir = root / "core", root / "comparison"
    core_dir.mkdir(parents=True)
    comparison_dir.mkdir()
    products = {
        name: core_dir / f"{name}.tsv"
        for name in (
            "count_profiles",
            "ora_profiles",
            "experiments",
            "feature_catalog",
            "feature_metadata",
            "taxa_matrix",
        )
    }
    features = tuple(f"qualified|{terms[0]}" for terms in TABLE1_TRAITS.values())
    if not subset:
        features += ("qualified|extra",)
    records: list[dict[str, Any]] = []
    for index in range(8):
        for direction in ("inc", "dec"):
            # The oriented subset already swaps the reversed source experiment.
            oriented = (
                "dec"
                if subset and index == 0 and direction == "inc"
                else "inc"
                if subset and index == 0
                else direction
            )
            counts = {
                feature: (
                    1 + (index + column) % 3 if (index + column + (direction == "inc")) % 5 else 0
                )
                for column, feature in enumerate(features)
            }
            if "qualified|extra" in counts:
                counts["qualified|extra"] = 1  # Ubiquitous: excluded by native support filter.
            records.append(
                {"study_id": f"s{index}", "experiment_id": "e1", "inc_dec": oriented, **counts}
            )
    write_table(products["count_profiles"], records)
    write_table(
        products["ora_profiles"],
        [
            {
                **{key: row[key] for key in ("study_id", "experiment_id", "inc_dec")},
                "taxa_count": 2,
                **{feature: int(row[feature] > 0) for feature in features},
            }
            for row in records
        ],
    )
    write_table(
        products["experiments"],
        [
            {
                "study_id": f"s{index}",
                "experiment_id": "e1",
                "condition": "fixture",
                "body_site": "Feces",
                "host_species": "Homo sapiens",
            }
            for index in range(8)
        ],
    )
    labels = [
        {
            "column_id": feature,
            "object_name": feature,
            "column_label": feature,
            "object_category": "biolink:PhenotypicQuality",
        }
        for feature in reversed(features)
    ]
    write_table(products["feature_metadata"], labels)
    write_table(
        products["feature_catalog"],
        [{**row, "feature_class": "trait_annotation"} for row in labels],
    )
    write_table(
        products["taxa_matrix"],
        [
            {
                "taxon_id": f"NCBITaxon:{index}",
                "taxon_name": f"Fixture {index}",
                "taxonomic_rank": "species",
                **dict.fromkeys(features, int(index == 1)),
            }
            for index in (1, 2)
        ],
    )
    full_dump = core_dir / "raw.csv"
    full_dump.write_text("fixture\n", encoding="utf-8")
    profile, experiments = products["count_profiles"], products["experiments"]
    noel_dir, jed_dir, ilyas_dir = root / "noel", root / "jed", root / "ilyas"
    ilyas_dir.mkdir()
    noel_profiles, corpus = load_profiles(profile, experiments, candidate_slots(features))
    write_screen_outputs(
        run_screen(noel_profiles, features),
        noel_dir,
        profile_path=profile,
        experiments_path=experiments,
        full_dump_path=full_dump,
        corpus_stats=corpus,
    )
    jed_profiles, corpus = load_profiles(profile, experiments, TABLE1_TRAITS)
    write_log_ratio_outputs(
        compute_log_ratio_shifts(jed_profiles),
        jed_dir,
        profile_path=profile,
        experiments_path=experiments,
        full_dump_path=full_dump,
        corpus_stats=corpus,
    )
    for tail, prefix in (("upper", "enrichment"), ("lower", "depletion")):
        ora_records = [
            {
                "study_id": row["study_id"],
                "experiment_id": "e1",
                "inc_dec": row["inc_dec"],
                "trait": feature,
                "p_value": 0.1 if tail == "upper" else 0.9,
                "fdr_q": 0.1 if tail == "upper" else 0.9,
                "fold": 1,
                "count": 1,
                "expected": 1,
                "testable_profiles": 16,
                "background_m": 2,
                "background_k": 1,
                "is_significant": False,
                "correction_family": "global",
                "eligibility_status": "testable",
            }
            for row in records
            for feature in features
        ]
        write_table(
            ilyas_dir / f"trait_{prefix}_all_results.tsv",
            ora_records,
            {
                "Profiles-SHA256": file_sha256(products["ora_profiles"]),
                "Taxa-Matrix-SHA256": file_sha256(products["taxa_matrix"]),
                "Tail": tail,
                "Testable-Set": "K>0 & N>0",
                "FDR": "True",
                "FDR-Scope": "global",
                "P-Threshold": "0.001",
                "Testable-Profiles": "16",
                "Background-Taxa-Universe-M": "2",
                "Total-Tests-Performed": str(len(ora_records)),
                "Background-Mode": "TaxaMatrix(fixture)",
                "Min-Fold" if tail == "upper" else "Max-Fold": "1.5" if tail == "upper" else "0.5",
            },
        )
    native_dir = comparison_dir / "native"
    pair_path = native_dir / "trait_profile_paired_contrast_pair_metadata.tsv"
    hashes = {
        f"{name}_sha256": file_sha256(products[product])
        for name, product in (
            ("profile", "count_profiles"),
            ("experiments", "experiments"),
            ("trait_labels", "feature_metadata"),
        )
    }
    write_table(
        pair_path,
        [
            {
                "study_id": f"s{index}",
                "experiment_id": "e1",
                "zero_total_policy": "included_in_historical_native_composition",
                **hashes,
            }
            for index in range(8)
        ],
    )
    selected, _ = native.native_candidates(
        profile, experiments, features, products["feature_metadata"]
    )
    stats = native_dir / "trait_profile_paired_contrast_trait_stats_pooled.tsv"
    write_table(
        stats,
        [
            {
                "trait": feature,
                "study_p_value": ".01",
                "study_q_value": ".01",
                "study_test_status": "tested",
                "stratum_kind": "pooled",
                "stratum": "pooled",
                "test_method": "two_sided_one_sample_t_test",
                "correction_method": "BH",
                "correction_family_size": len(selected),
                "retained_trait_count": len(features),
                "pair_metadata_sha256": file_sha256(pair_path),
                **hashes,
            }
            for feature in selected
        ],
    )
    consumer_path = native_dir / "core_consumption.json"
    consumer_path.write_text("{}\n")
    population = (
        ContrastPopulation.disease_vs_healthy_tier1_2_subject_controls
        if subset
        else ContrastPopulation.all_paired
    )
    direction_basis = (
        DirectionBasis.disease_vs_comparison if subset else DirectionBasis.group_1_vs_group_0
    )
    core_manifest = core_dir / "run.json"
    core_manifest.write_text(json.dumps({"synthetic": True, "population": population}))
    core = WorkflowRun.model_construct(
        run_id="fixture-core",
        manifest_path=str(core_manifest),
        contrast_population=population,
        direction_basis=direction_basis,
        parameters=[
            ParameterValue(
                name="strain_rollup", value="phenotype" if subset else "off", source="override"
            )
        ],
    )
    arms = {}
    for name, directory, steps in (
        ("noel", noel_dir, ["directional_screen"]),
        ("jed", jed_dir, ["composition"]),
        ("ilyas", ilyas_dir, ["enrichment", "depletion"]),
    ):
        manifest = directory / "arm_run.json"
        manifest.write_text("{}\n")
        snapshot = RecordedPath(
            path=str(directory),
            storage_kind=StorageKind.collection,
            sha256="a" * 64,
            size_bytes=0,
            inventory=[file_entry(path, relative_to=directory) for path in directory.iterdir()],
        )
        arms[name] = ArmRun.model_construct(
            arm=name,
            manifest_path=str(manifest),
            steps=[
                ArmCommandExecution.model_construct(
                    id=step,
                    outputs=[snapshot, *(file_snapshot(path) for path in directory.iterdir())],
                )
                for step in steps
            ],
            core_manifest_sha256=file_sha256(core_manifest),
        )
    bases = {
        product: make_basis(
            products["feature_catalog"], features, product, "phenotype" if subset else "off"
        )
        for product in ("count_profiles", "ora_profiles")
    }
    comparison_manifest = comparison_dir / "comparison_run.json"
    comparison_manifest.write_text("{}\n")
    comparison = ComparisonRun.model_construct(
        manifest_path=str(comparison_manifest),
        core_manifest_path=str(core_manifest),
        core_manifest_sha256=file_sha256(core_manifest),
        core_run_id=core.run_id,
        contrast_population=population,
        direction_basis=direction_basis,
        analysis_bases=list(bases.values()),
        outputs=[file_entry(path, relative_to=comparison_dir) for path in native_dir.iterdir()],
        noel_manifest_path=arms["noel"].manifest_path,
        jed_manifest_path=arms["jed"].manifest_path,
        ilyas_manifest_path=arms["ilyas"].manifest_path,
    )
    verified = VerifiedSources(
        comparison,
        core,
        arms,
        bases,
        dict.fromkeys(bases, features),
        "phenotype" if subset else "off",
    )
    monkeypatch.setattr(sources, "artifact_path", lambda _run, name: products[name])
    monkeypatch.setattr(stage, "open_verified_sources", lambda _path: verified)
    monkeypatch.setattr(
        native,
        "verify_consumer_record",
        lambda _path: ConsumerRecord.model_construct(
            consumer="native-pooled-composition-stats",
            run_id=core.run_id,
            run_manifest_sha256=file_sha256(core_manifest),
            outputs=[file_entry(stats), file_entry(pair_path)],
        ),
    )
    return verified
