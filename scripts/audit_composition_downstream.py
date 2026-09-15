"""Replay the legacy downstream workflows with only the #174 normalization change."""

from __future__ import annotations

import argparse
import importlib
import importlib.metadata
import json
import shutil
import sys
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from bsdbng.contribution_inputs import digest
from bsdbng.trait_composition import DENOMINATOR_DEFINITION

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / ".tmp/issue174"
BASE = WORK / "baseline"
REPORT = ROOT / "trait-visualization/denominator-harmonization/downstream"
AFTER = ROOT / "trait-visualization/trait-profile-investigation"
BEFORE = WORK / "downstream-before"
MATRIX = WORK / "study_trait_profiles_wide.tsv"


@contextmanager
def arguments(args: list[str]):
    original = sys.argv
    sys.argv = [original[0], *args]
    try:
        yield
    finally:
        sys.argv = original


def module(name: str):
    loaded = importlib.import_module(name)
    # Imported renderers also use this root when writing relative artifact links.
    for imported in list(sys.modules.values()):
        path = getattr(imported, "__file__", None)
        if (
            path
            and Path(path).parent in {ROOT / "scripts", BASE / "scripts"}
            and hasattr(imported, "REPO_ROOT")
        ):
            imported.REPO_ROOT = ROOT
    return loaded


def run(stage: str, data: Path, resume_sankey_exports: bool = False) -> None:
    before = stage == "before"
    code = BASE if before else ROOT
    vis = code / "trait-visualization"
    output = BEFORE if before else AFTER
    output.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(code / "scripts"))
    auxiliary = [
        data / "raw/full_dump.csv",
        data / "raw/experiments.csv",
        data / "derived/studies.tsv",
        data / "derived/taxa.tsv",
        data / "derived/taxa_trait_matrix.tsv",
        ROOT / "trait-visualization/denominator-harmonization/inputs/taxa_trait_matrix_columns.tsv",
    ]
    labels = auxiliary[-1]
    inputs = {str(path): digest(path) for path in [MATRIX, *auxiliary]}
    expected_matrix = "c07d589f6cc88724e43c683ebd132bcd5c7f316f36b95df3ff9acad24ebbf62f"
    if inputs[str(MATRIX)] != expected_matrix:
        raise ValueError("The downstream audit requires the exact issue #174 matrix")

    sankey = module("create_six_shared_trait_sankeys")
    sankey.INVESTIGATION_DIR = output
    sankey.OUTPUT_DIR = output / "six-shared-traits"
    for name in [
        "create_six_shared_trait_graphs",
        "create_six_shared_trait_static_graphs",
        "create_six_shared_trait_cytoscape_heatmaps",
    ]:
        module(name).OUTPUT_DIR = sankey.OUTPUT_DIR
    sankey.PROFILE_MATRIX_PATH = MATRIX
    sankey.PROFILE_METADATA_PATH = vis / "study_trait_profile_study_metadata.tsv"
    sankey.TRAIT_STATS_PATH = (
        vis / "trait-profile-investigation/trait_profile_paired_contrast_trait_stats.tsv"
    )
    if resume_sankey_exports:
        sankey.write_six_shared_trait_graphs(sankey.OUTPUT_DIR, ROOT)
        sankey.write_six_shared_trait_static_graphs(sankey.OUTPUT_DIR)
        sankey.write_six_shared_trait_cytoscape_heatmaps(sankey.OUTPUT_DIR, ROOT)
    else:
        sankey.main()

    pd_analysis = module("create_pd_cellulose_analysis")
    pd_analysis.PROFILE_MATRIX_PATH = MATRIX
    pd_analysis.PROFILE_METADATA_PATH = sankey.PROFILE_METADATA_PATH
    pd_analysis.TAXA_PATH = data / "derived/taxa.tsv"
    pd_analysis.TAXA_TRAIT_MATRIX_PATH = data / "derived/taxa_trait_matrix.tsv"
    pd_analysis.SIX_SHARED_DIR = sankey.OUTPUT_DIR
    pd_analysis.OUTPUT_DIR = output / "PD-cellulose"
    pd_analysis.main()

    guild = module("create_condition_body_site_guild_maps")
    guild.MATRIX_PATH = vis / "study_trait_profile_clustered_heatmap_matrix.tsv"
    guild.METADATA_PATH = sankey.PROFILE_METADATA_PATH
    guild.TRAIT_METADATA_PATH = vis / "study_trait_profile_trait_metadata.tsv"
    guild.PROFILE_MATRIX_PATH = MATRIX
    for setting in ["FACULTATIVE_CONFIG", "ANAEROBE_CONFIG"]:
        config = getattr(guild, setting)
        directory = output / config.output_dir.name
        directory.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(
            BASE
            / "trait-visualization/trait-profile-investigation"
            / directory.name
            / "guild_trait_selection.tsv",
            directory / "guild_trait_selection.tsv",
        )
        setattr(guild, setting, replace(config, output_dir=directory))
    guild.main()

    condition = module("create_condition_trait_module_heatmap")
    destination = output / "condition-modules"
    destination.mkdir(parents=True, exist_ok=True)
    options = {
        "study-trait-profiles-input": MATRIX,
        "trait-columns-input": labels,
        "studies-input": data / "derived/studies.tsv",
        "full-dump-input": data / "raw/full_dump.csv",
        "output": destination / "condition_trait_module_clustered_heatmap.png",
        "matrix-output": destination / "condition_trait_module_clustered_heatmap_matrix.tsv",
        "group-metadata-output": destination / "condition_trait_module_group_metadata.tsv",
        "trait-metadata-output": destination / "condition_trait_module_trait_metadata.tsv",
        "guild-definition-table-output": destination / "guild_definition_table.png",
        "guild-definition-table-svg-output": destination / "guild_definition_table.svg",
    }
    with arguments([part for key, value in options.items() for part in [f"--{key}", str(value)]]):
        condition.main()

    followup = module("analyze_enrichment_infographic_followup")
    options = {
        "profiles-input": MATRIX,
        "raw-experiments-input": data / "raw/experiments.csv",
        "full-dump-input": data / "raw/full_dump.csv",
        "studies-input": data / "derived/studies.tsv",
        "taxa-trait-matrix-input": data / "derived/taxa_trait_matrix.tsv",
        "taxa-trait-columns-input": labels,
        "output-dir": output / "infographic-followup",
    }
    with arguments([part for key, value in options.items() for part in [f"--{key}", str(value)]]):
        followup.main()
    sources = {
        str(path.relative_to(code)): digest(path) for path in (code / "scripts").glob("*.py")
    }
    inputs.update(
        {
            str(path): digest(path)
            for path in [
                sankey.PROFILE_METADATA_PATH,
                sankey.TRAIT_STATS_PATH,
                guild.MATRIX_PATH,
                guild.TRAIT_METADATA_PATH,
            ]
        }
    )
    REPORT.mkdir(parents=True, exist_ok=True)
    manifest = {
        "stage": stage,
        "code_commit": (WORK / "code_commit.txt").read_text().strip(),
        "baseline_commit": "a2b2df98e9beec98b772c9cb5af7aa090d033a98",
        "denominator": "all loaded trait columns" if before else DENOMINATOR_DEFINITION,
        "inputs_sha256": inputs,
        "source_sha256": sources,
        "package_versions": {
            name: importlib.metadata.version(name)
            for name in ["pandas", "numpy", "scipy", "matplotlib", "seaborn"]
        },
        "settings": "Baseline defaults; frozen guilds, paired metadata and Table 1 selection. "
        "Condition modules use the human-study aggregate; follow-up uses the full matrix. "
        "Core pairing, z-score ddof, significance families and support thresholds are unchanged.",
    }
    (REPORT / f"{stage}_provenance.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


def compare() -> None:
    before_manifest = json.loads((REPORT / "before_provenance.json").read_text())
    after_manifest = json.loads((REPORT / "after_provenance.json").read_text())
    for path, checksum in before_manifest["inputs_sha256"].items():
        if (
            path in after_manifest["inputs_sha256"]
            and after_manifest["inputs_sha256"][path] != checksum
        ):
            raise ValueError(f"Input changed between stages: {path}")
    rows = []
    changes = []
    for old_path in sorted(BEFORE.rglob("*.tsv")):
        relative = old_path.relative_to(BEFORE)
        new_path = AFTER / relative
        if not new_path.exists():
            raise ValueError(f"Missing regenerated output: {relative}")
        old = pd.read_csv(old_path, sep="\t")
        new = pd.read_csv(new_path, sep="\t")
        snapshot = REPORT / "before" / relative
        snapshot.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(old_path, snapshot)
        keys = [
            column
            for column in [
                "analysis",
                "view_key",
                "trait_id",
                "focus_direction",
                "condition",
                "body_site_group",
                "body_site",
                "condition_body_site_pair",
                "guild",
                "outcome",
                "term",
                "group",
                "taxon_id",
                "profile_id",
                "neighbor_profile_id",
                "neighbor_rank",
            ]
            if column in old and column in new
        ]
        if not keys or old.duplicated(keys).any() or new.duplicated(keys).any():
            first = old.columns[0]
            keys = [first] if first in new and old[first].is_unique and new[first].is_unique else []
        record = {
            "artifact": relative.as_posix(),
            "rows_before": len(old),
            "rows_after": len(new),
            "exactly_equal": old.equals(new),
            "comparison_keys": ";".join(keys),
        }
        if keys:
            left, right = old.set_index(keys), new.set_index(keys)
            record.update(
                {
                    "rows_removed": len(left.index.difference(right.index)),
                    "rows_added": len(right.index.difference(left.index)),
                }
            )
            common = left.index.intersection(right.index)
            for column in left.columns.intersection(right.columns):
                a, b = left.loc[common, column], right.loc[common, column]
                if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
                    equal = np.isclose(
                        a.to_numpy(dtype=float),
                        b.to_numpy(dtype=float),
                        rtol=1e-9,
                        atol=1e-12,
                        equal_nan=True,
                    )
                else:
                    equal = (a.fillna("") == b.fillna("")).to_numpy()
                changed = int((~equal).sum())
                if changed:
                    changes.append(
                        {
                            "artifact": relative.as_posix(),
                            "column": column,
                            "matched_rows": len(common),
                            "changed_rows": changed,
                        }
                    )
        rows.append(record)
    pd.DataFrame(rows).to_csv(REPORT / "artifact_comparison.tsv", sep="\t", index=False)
    pd.DataFrame(changes).to_csv(REPORT / "changed_columns.tsv", sep="\t", index=False)
    compare_heatmap_families()
    sys.path.insert(0, str(ROOT / "scripts"))
    guild = module("create_condition_body_site_guild_maps")
    guild.PROFILE_MATRIX_PATH = MATRIX
    guild.load_study_heatmap()
    provenance = {
        "implementation_base_commit": (WORK / "code_commit.txt").read_text().strip(),
        "denominator": DENOMINATOR_DEFINITION,
        "comparison_source_sha256": digest(Path(__file__)),
        "shared_source_sha256": {
            path.name: digest(path)
            for path in [
                ROOT / "src/bsdbng/trait_composition.py",
                ROOT / "src/bsdbng/contribution_inputs.py",
            ]
        },
        "package_versions": {
            name: importlib.metadata.version(name)
            for name in ["pandas", "numpy", "scipy", "matplotlib", "seaborn"]
        },
        "stage_provenance_sha256": {
            stage: digest(REPORT / f"{stage}_provenance.json") for stage in ["before", "after"]
        },
        "compared_tables": len(rows),
        "exactly_equal_tables": sum(row["exactly_equal"] for row in rows),
        "transformed_cache_verified_against_raw_counts": True,
        "source_commit_note": "Base commit and source hashes identify the rerun implementation.",
    }
    (REPORT / "comparison_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"Compared {len(rows)} downstream tables; "
        f"{sum(not row['exactly_equal'] for row in rows)} differ"
    )


def compare_heatmap_families() -> None:
    """Keep each original view's correction family, including unselected cells."""
    frames = {}
    manifests = {}
    for stage, root in [("before", BEFORE), ("after", AFTER)]:
        directory = root / "six-shared-traits"
        frames[stage] = pd.concat(
            [
                pd.read_csv(
                    directory / f"six_shared_trait_condition_body_site_{suffix}.tsv", sep="\t"
                )
                for suffix in ["trait_pattern_matrix", "direction_contrasts"]
            ],
            ignore_index=True,
        )
        manifests[stage] = pd.read_csv(
            directory / "six_shared_trait_cytoscape_heatmap_manifest.tsv", sep="\t"
        )
    keys = ["view_key", "condition_body_site_pair", "trait_id"]
    cells = frames["before"].merge(
        frames["after"],
        on=keys,
        how="outer",
        suffixes=("_before", "_after"),
        indicator=True,
        validate="one_to_one",
    )
    cells.to_csv(REPORT / "heatmap_cells_before_after.tsv", sep="\t", index=False)
    columns = ["view_key", "tested_cell_count", "displayed_cell_count", "bh_selected_cell_count"]
    families = manifests["before"][columns].merge(
        manifests["after"][columns],
        on="view_key",
        how="outer",
        suffixes=("_before", "_after"),
        validate="one_to_one",
    )
    families["family_definition"] = (
        "all emitted Fisher cells within this view, before p/q selection"
    )
    families["q_cutoff"] = 0.05
    families.to_csv(REPORT / "heatmap_family_comparison.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--stage", choices=["before", "after", "compare"], required=True)
    parser.add_argument(
        "--resume-sankey-exports",
        action="store_true",
        help="Resume an interrupted renderer export with existing Sankey tables.",
    )
    args = parser.parse_args()
    pins = json.loads((REPORT.parent / "inputs/downstream_auxiliary_sha256.json").read_text())
    for relative, expected in pins.items():
        if digest(args.data_root / relative) != expected:
            raise ValueError(f"Pinned downstream auxiliary input mismatch: {relative}")
    if args.stage == "compare":
        compare()
    else:
        run(args.stage, args.data_root.resolve(), args.resume_sankey_exports)


if __name__ == "__main__":
    main()
