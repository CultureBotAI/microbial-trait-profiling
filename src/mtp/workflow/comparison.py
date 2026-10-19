"""Compare a verified core rebuild with legacy tables without inventing provenance (#629)."""

from __future__ import annotations

import json
from pathlib import Path

from mtp.count_input import profile_count_reader, taxa_trait_count_reader
from mtp.paths import file_sha256
from mtp.workflow.runtime import artifact_path
from mtp.workflow.verification import verify_core_run

BASELINE_FILES = {
    "taxa_matrix": "taxa_trait_matrix.tsv",
    "count_profiles": "study_trait_profiles_wide.tsv",
    "ora_profiles": "study_trait_profiles_ora_wide.tsv",
}


def compare_count_tables(
    baseline: Path, current: Path, *, matrix: bool = False
) -> dict[str, object]:
    """Compare typed counts by row and qualified feature ID, independent of column order."""
    reader_factory = taxa_trait_count_reader if matrix else profile_count_reader
    before = reader_factory(baseline)
    after = reader_factory(current)
    baseline_rows = {row.key: row.counts for row in before}
    columns_before = {column: index for index, column in enumerate(before.count_columns)}
    columns_after = {column: index for index, column in enumerate(after.count_columns)}
    shared = sorted(columns_before.keys() & columns_after.keys())
    selected = [(column, columns_before[column], columns_after[column]) for column in shared]
    added_rows: list[list[str]] = []
    changed_cells = 0
    changed_nonmedia_cells = 0
    compared_cells = 0
    examples: list[dict[str, object]] = []
    rows = 0
    for row in after:
        rows += 1
        old = baseline_rows.pop(row.key, None)
        if old is None:
            added_rows.append(list(row.key))
            continue
        compared_cells += len(selected)
        for feature, left, right in selected:
            if old[left] != row.counts[right]:
                changed_cells += 1
                if "mediadive" not in feature.lower():
                    changed_nonmedia_cells += 1
                if len(examples) < 20:
                    examples.append(
                        {
                            "row": list(row.key),
                            "feature": feature,
                            "baseline": old[left],
                            "current": row.counts[right],
                        }
                    )
    removed_rows = [list(key) for key in sorted(baseline_rows)]
    added_features = sorted(columns_after.keys() - columns_before.keys())
    removed_features = sorted(columns_before.keys() - columns_after.keys())
    return {
        "baseline_path": str(baseline.resolve()),
        "baseline_sha256": file_sha256(baseline),
        "current_path": str(current.resolve()),
        "current_sha256": file_sha256(current),
        "current_rows": rows,
        "shared_features": len(shared),
        "compared_cells": compared_cells,
        "changed_cells": changed_cells,
        "changed_nonmedia_cells": changed_nonmedia_cells,
        "added_rows": added_rows,
        "removed_rows": removed_rows,
        "added_features": added_features,
        "removed_features": removed_features,
        "changed_examples": examples,
        "identical_counts": not any(
            (changed_cells, added_rows, removed_rows, added_features, removed_features)
        ),
    }


def compare_legacy_bundle(manifest: Path, baseline_dir: Path, output: Path) -> Path:
    """Write factual numerical migration evidence; baseline execution remains unverified."""
    if output.exists():
        raise FileExistsError(output)
    run = verify_core_run(manifest)
    comparisons = {
        product: compare_count_tables(
            baseline_dir / filename, artifact_path(run, product), matrix=product == "taxa_matrix"
        )
        for product, filename in BASELINE_FILES.items()
    }
    report = {
        "run_id": run.run_id,
        "run_manifest": str(manifest.resolve()),
        "run_manifest_sha256": file_sha256(manifest),
        "execution_fingerprint": run.execution_fingerprint,
        "baseline_verification_status": "imported_unverified",
        "interpretation": "This compares actual table values. Matching values do not establish "
        "the historical baseline's producing inputs, execution, or biological validity.",
        "comparisons": comparisons,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return output
