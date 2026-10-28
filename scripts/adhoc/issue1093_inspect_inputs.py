"""Inspect #1093 inputs without regenerating profiles or figures (2026-10-27).

Inputs: current checkout, explicit column labels, F archive root, and the
historical reference JSON. Optional historical working-tree inventory requires
an explicit argument; current development uses only the tracked canonical files.
Outputs: diagnostic JSON, native stool cells,
group ranking TSV, and a short Markdown findings file in the requested directory.
All output paths identify source roots symbolically; no input is modified.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib
import json
import sys
from pathlib import Path
from typing import Any

import git
import pandas as pd

from mtp.composition_denominators import retained_trait_counts
from mtp.trait_lists import load_trait_list

KEYS = ["study_id", "experiment_id", "inc_dec"]
OXYGEN_IDS = ["METPO:1000605", "METPO:1000608", "METPO:1000612"]
PROFILE = Path("data/derived/study_trait_profiles_wide.tsv")
METADATA = Path("trait-visualization/study_trait_profile_study_metadata.tsv")
COLUMNS = Path("data/derived/taxa_trait_matrix_columns.tsv")


def digest(path: Path) -> str:
    """Hash existing bytes without interpreting or changing them."""
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def profile_inventory(path: Path) -> dict[str, Any]:
    """Count actual profile rows, studies, features, and oxygen annotations."""
    frame = pd.read_csv(path, sep="\t")
    traits = frame.drop(columns=KEYS)
    retained = retained_trait_counts(traits)
    return {
        "profile_count": len(frame),
        "study_count": int(frame["study_id"].nunique()),
        "experiment_count": len(frame[["study_id", "experiment_id"]].drop_duplicates()),
        "column_count": len(frame.columns),
        "all_trait_columns": len(traits.columns),
        "retained_non_mediadive_traits": len(retained.columns),
        "oxygen_columns": {
            key: {
                "present": key in retained,
                "positive_profiles": int(retained[key].gt(0).sum()) if key in retained else 0,
                "total_count": int(retained[key].sum()) if key in retained else 0,
            }
            for key in OXYGEN_IDS
        },
    }


def records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    """Convert pandas scalars and missing values to strict JSON values."""
    result: list[dict[str, Any]] = json.loads(frame.to_json(orient="records"))
    return result


def compare_baseline(baseline_path: Path, output: Path, result: dict[str, Any]) -> dict[str, Any]:
    """Check actual baseline bytes, input binding, and all four regenerated tables."""
    payload = baseline_path.read_bytes()
    baseline = json.loads(payload)
    old_tables = baseline["default_six_baseline"]
    new_tables = result["default_six_baseline"]
    if set(old_tables) != set(new_tables) or len(old_tables) != 4:
        raise ValueError("Expected the same four default-six baseline tables")
    for key in ("$CHECKOUT/" + str(PROFILE) + ".gz", "$CHECKOUT/" + str(METADATA)):
        if baseline["inputs"][key] != result["inputs"][key]:
            raise ValueError("Baseline comparison requires identical canonical profiles/metadata")
    comparison: dict[str, Any] = {}
    for name, expected in old_tables.items():
        previous = (baseline_path.parent / name).read_bytes()
        current = (output / name).read_bytes()
        if hashlib.sha256(previous).hexdigest() != expected["sha256"]:
            raise ValueError("Recorded baseline table hash does not match its current bytes")
        comparison[name] = {
            "baseline_sha256": expected["sha256"],
            "current_sha256": new_tables[name]["sha256"],
            "baseline_rows": expected["rows"],
            "current_rows": new_tables[name]["rows"],
            "bytes_equal": previous == current,
        }
    report = {
        "baseline_record_sha256": hashlib.sha256(payload).hexdigest(),
        "baseline_code_commit": baseline["code_commit"],
        "current_code_commit": result["code_commit"],
        "canonical_profile_and_metadata_hashes_equal": True,
        "implementation_changes": {
            name: {
                "baseline_sha256": baseline["source_hashes"].get(name),
                "current_sha256": checksum,
            }
            for name, checksum in result["source_hashes"].items()
            if baseline["source_hashes"].get(name) != checksum
        },
        "tables": comparison,
        "all_four_tables_byte_identical": all(row["bytes_equal"] for row in comparison.values()),
    }
    (output / "default_six_comparison.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    if not report["all_four_tables_byte_identical"]:
        raise ValueError("Default-six output bytes changed; see default_six_comparison.json")
    return report


def metadata_causes(args: argparse.Namespace) -> None:
    """Attribute the reference count differences using raw canonical metadata."""
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "scripts"))
    guild = importlib.import_module("create_condition_body_site_guild_maps")
    metadata_path = repo_root / METADATA
    code_path = repo_root / "scripts/create_condition_body_site_guild_maps.py"
    hashes = {
        "$CHECKOUT/" + str(METADATA): digest(metadata_path),
        "$CHECKOUT/scripts/create_condition_body_site_guild_maps.py": digest(code_path),
        "$REFERENCE/ten_traits_stats_phenotype.json": digest(args.reference),
    }
    frame = pd.read_csv(metadata_path, sep="\t", index_col=0)
    frame["native_group"] = frame["condition"].map(
        lambda value: guild.v3_fully_grouped_condition_display_label(value)[0]
    )
    frame["native_body_group"] = [
        guild.v3_body_site_display_label(
            row.body_site, row.body_site_group, guild.FACULTATIVE_CONFIG
        )[0]
        for row in frame.itertuples()
    ]
    masks = {
        "exact_raw_body_site": frame["body_site"].eq("Feces"),
        "cached_body_site_group": frame["body_site_group"].eq("Feces"),
        "native_grouped_body_site": frame["native_body_group"].eq("Feces"),
    }
    reference = json.loads(args.reference.read_bytes())["s14b_canonical"]
    counts = []
    for row in reference:
        for name, mask in masks.items():
            group = frame.loc[mask & frame["native_group"].eq(row["group"])]
            counts.append(
                {
                    "group": row["group"],
                    "rule": name,
                    "studies": int(group["study_id"].nunique()),
                    "increased_profiles": int(group["inc_dec"].eq("inc").sum()),
                    "decreased_profiles": int(group["inc_dec"].eq("dec").sum()),
                }
            )
    excluded = frame.loc[
        ~masks["exact_raw_body_site"]
        & (masks["cached_body_site_group"] | masks["native_grouped_body_site"])
    ].copy()
    columns = [
        "study_id",
        "experiment_id",
        "inc_dec",
        "condition",
        "body_site",
        "body_site_group",
        "native_group",
        "native_body_group",
    ]
    result = {
        "input_hashes": hashes,
        "counts_by_stool_rule": counts,
        "excluded_by_exact_rule": records(excluded[columns]),
        "note": "Raw body_site strings are compared literally, without clean_text.",
    }
    if hashes != {
        "$CHECKOUT/" + str(METADATA): digest(metadata_path),
        "$CHECKOUT/scripts/create_condition_body_site_guild_maps.py": digest(code_path),
        "$REFERENCE/ten_traits_stats_phenotype.json": digest(args.reference),
    }:
        raise RuntimeError("Metadata attribution inputs changed while read")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    (args.output_dir / "metadata_causes.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, allow_nan=False))


def run(args: argparse.Namespace) -> None:
    """Inspect explicit inputs using the unchanged native analysis functions."""
    if args.compare_baseline is not None and not args.six_baseline:
        raise ValueError("--compare-baseline requires --six-baseline")
    repo_root = Path(__file__).resolve().parents[2]
    repo = git.Repo(repo_root)
    sys.path.insert(0, str(repo_root / "scripts"))
    sankey = importlib.import_module("create_six_shared_trait_sankeys")
    heatmap = importlib.import_module("create_six_shared_trait_cytoscape_heatmaps")
    guild = importlib.import_module("create_condition_body_site_guild_maps")
    source_names = [
        "scripts/create_six_shared_trait_sankeys.py",
        "scripts/create_six_shared_trait_cytoscape_heatmaps.py",
        "scripts/create_condition_body_site_guild_maps.py",
        "src/mtp/composition_denominators.py",
        "src/mtp/trait_lists.py",
        "methods/trait_lists/headline_ten_phenotype.txt",
    ]
    added_source = "src/mtp/stool_group_contrasts.py"
    if (repo_root / added_source).is_file():
        source_names.append(added_source)
    source_hashes = {name: digest(repo_root / name) for name in source_names}
    full_candidates = list(
        args.f_root.glob(
            "runs/full-phenotype/core/run-*/steps/count_profiles/derived/"
            "study_trait_profiles_wide.tsv"
        )
    )
    if len(full_candidates) != 1:
        raise ValueError("Expected exactly one archived full-phenotype core profile table")
    paths = {
        "$CHECKOUT/" + str(PROFILE) + ".gz": repo_root / (str(PROFILE) + ".gz"),
        "$CHECKOUT/" + str(METADATA): repo_root / METADATA,
        "$F/wt/" + str(PROFILE): args.f_root / "wt" / PROFILE,
        "$F/" + str(full_candidates[0].relative_to(args.f_root)): full_candidates[0],
        "$REFERENCE/ten_traits_stats_phenotype.json": args.reference,
    }
    if args.compare_baseline is not None:
        paths["$BASELINE/input_diagnostics.json"] = args.compare_baseline
    if args.main_input_root is not None:
        paths.update(
            {
                "$MAIN_INPUTS/" + str(relative): args.main_input_root / relative
                for relative in (PROFILE, METADATA, COLUMNS)
            }
        )
    columns_path = args.columns_file
    if columns_path is None:
        if args.main_input_root is None:
            raise ValueError("Supply an existing archived --columns-file")
        columns_path = args.main_input_root / COLUMNS
    paths["$COLUMN_METADATA/" + columns_path.name] = columns_path
    before = {
        name: {"sha256": digest(path), "size_bytes": path.stat().st_size}
        for name, path in paths.items()
    }
    inventories = {
        name: {**before[name], **profile_inventory(path)}
        for name, path in paths.items()
        if path.name in {PROFILE.name, PROFILE.name + ".gz"}
    }
    compressed_path = repo_root / (str(PROFILE) + ".gz")
    with gzip.open(compressed_path, "rb") as handle:
        decompressed_digest = hashlib.file_digest(handle, "sha256").hexdigest()
    native_paths = {
        "profiles": compressed_path,
        "metadata": repo_root / METADATA,
        "columns": columns_path,
    }
    if args.input_basis == "main-defaults":
        if args.main_input_root is None:
            raise ValueError("Historical main-defaults mode requires its explicit input root")
        native_paths["profiles"] = args.main_input_root / PROFILE
        native_paths["metadata"] = args.main_input_root / METADATA
    native_choice = (
        "$CHECKOUT tracked gzip + tracked metadata; explicit $COLUMN_METADATA labels"
        if args.input_basis == "tracked-canonical"
        else "$MAIN_INPUTS actual default files; no automatic gzip fallback"
    )
    trait_list = load_trait_list(repo_root / "methods/trait_lists/headline_ten_phenotype.txt")
    traits = sankey.declared_shared_traits(trait_list, native_paths["columns"])
    metadata = sankey.load_profile_metadata(native_paths["metadata"])
    raw_counts = sankey.load_raw_trait_counts(traits, metadata, native_paths["profiles"])
    zscores = sankey.build_zscore_matrix(raw_counts, metadata, traits)
    stool = metadata.loc[metadata["body_site"].eq("Feces")].copy()
    edges = sankey.build_trait_direction_tables(stool, zscores.loc[stool.index], traits)[0]
    trait_order = pd.DataFrame(
        {
            "trait_id": [trait.trait_id for trait in traits],
            "trait_label": [trait.label for trait in traits],
            "slug": [trait.slug for trait in traits],
            "base_color": [trait.base_color for trait in traits],
            "trait_rank": range(len(traits)),
        }
    )
    cells = heatmap.prepare_contrast_cells(
        edges.loc[edges["focus_direction"].eq("increased")],
        edges.loc[edges["focus_direction"].eq("decreased")],
        next(mode for mode in heatmap.CONTRAST_MODES if mode.key == "grouped"),
        trait_order,
    )
    scored = sankey.score_profiles(stool, zscores.loc[stool.index], traits[0])
    grouped = guild.add_v3_display_rollups(
        scored,
        sankey.single_trait_map_config(traits[0], "inc"),
        display_mode=guild.V3_GROUPED_MODE,
    )[0]
    groups = grouped.groupby("condition", as_index=False).agg(
        stool_studies=("study_id", "nunique"),
        stool_profiles=("profile_id", "count"),
        increased_profiles=("inc_dec", lambda values: int(values.eq("inc").sum())),
        decreased_profiles=("inc_dec", lambda values: int(values.eq("dec").sum())),
    )
    support = cells.loc[cells["passes_support"].astype(bool)].groupby("condition").size()
    groups["supported_cells"] = groups["condition"].map(support).fillna(0).astype(int)
    groups["passes_support"] = groups["supported_cells"].gt(0)
    groups = groups.sort_values(["stool_studies", "condition"], ascending=[False, True])
    selected = groups.loc[groups["passes_support"]].head(6)
    selected_names = selected["condition"].tolist()
    groups["selected"] = groups["condition"].isin(selected_names)
    reference = json.loads(args.reference.read_bytes())["s14b_canonical"]
    reference_traits = list(reference[0]["traits"])
    comparisons: list[dict[str, Any]] = []
    cell_comparisons: list[dict[str, Any]] = []
    for old in reference:
        current_group = groups.loc[groups["condition"].eq(old["group"])]
        if len(current_group) != 1:
            raise ValueError("Historical group does not have one native current group")
        new = current_group.iloc[0]
        comparisons.append(
            {
                "group": old["group"],
                "old_studies": old["n_stud"],
                "new_studies": int(new["stool_studies"]),
                "old_increased_profiles": old["inc_tot"],
                "new_increased_profiles": int(new["increased_profiles"]),
                "old_decreased_profiles": old["dec_tot"],
                "new_decreased_profiles": int(new["decreased_profiles"]),
                "selected_by_current_rule": old["group"] in selected_names,
            }
        )
        for key in trait_list.trait_ids:
            match = cells.loc[cells["condition"].eq(old["group"]) & cells["trait_id"].eq(key)]
            if len(match) != 1:
                raise ValueError("Expected exactly one historical-group/current-trait cell")
            current = match.iloc[0]
            old_cell = old["traits"][key]
            cell_comparisons.append(
                {
                    "group": old["group"],
                    "trait_id": key,
                    "old_delta_pp": old_cell["delta_pp"],
                    "new_delta_pp": 100
                    * float(current["increased_minus_decreased_positive_fraction"]),
                    "old_p": old_cell["p"],
                    "new_p": float(current["fisher_p_value"]),
                    "old_inc_pos": old_cell["inc_pos"],
                    "new_inc_pos": int(current["increased_positive_profile_count"]),
                    "old_dec_pos": old_cell["dec_pos"],
                    "new_dec_pos": int(current["decreased_positive_profile_count"]),
                }
            )
    six_tables: dict[str, pd.DataFrame] = {}
    if args.six_baseline:
        six_traits = sankey.SHARED_TRAITS
        six_zscores = sankey.build_zscore_matrix(raw_counts, metadata, six_traits)
        six_order = pd.DataFrame(
            {
                "trait_id": [trait.trait_id for trait in six_traits],
                "trait_label": [trait.label for trait in six_traits],
                "slug": [trait.slug for trait in six_traits],
                "base_color": [trait.base_color for trait in six_traits],
                "trait_rank": range(len(six_traits)),
            }
        )
        six_edges = {
            "grouped": sankey.build_trait_direction_tables(metadata, six_zscores)[0],
            "individual_conditions": sankey.build_individual_condition_direction_tables(
                metadata, six_zscores
            )[0],
        }
        for mode in heatmap.CONTRAST_MODES:
            mode_edges = six_edges[mode.key]
            six_tables[f"six_{mode.key}_edges.tsv"] = mode_edges
            six_tables[f"six_{mode.key}_contrasts.tsv"] = heatmap.prepare_contrast_cells(
                mode_edges.loc[mode_edges["focus_direction"].eq("increased")],
                mode_edges.loc[mode_edges["focus_direction"].eq("decreased")],
                mode,
                six_order,
            )
    metadata_raw = pd.read_csv(native_paths["metadata"], sep="\t", index_col=0)
    columns = pd.read_csv(native_paths["columns"], sep="\t")
    after = {
        name: {"sha256": digest(path), "size_bytes": path.stat().st_size}
        for name, path in paths.items()
    }
    if before != after or source_hashes != {
        name: digest(repo_root / name) for name in source_names
    }:
        raise RuntimeError("An input or native implementation changed during inspection")
    result: dict[str, Any] = {
        "issue": 1093,
        "code_commit": repo.head.commit.hexsha,
        "git_dirty": repo.is_dirty(untracked_files=True),
        "verification_scope": "read-only input inspection, not a core or arm rerun",
        "source_hashes": source_hashes,
        "input_hashes_before_and_after_equal": True,
        "inputs": before,
        "profile_inventories": inventories,
        "tracked_gzip_decompressed_sha256": decompressed_digest,
        "tracked_gzip_last_commit": next(
            repo.iter_commits(paths=str(PROFILE) + ".gz", max_count=1)
        ).hexsha,
        "tracked_gzip_matches_main_uncompressed": (
            decompressed_digest == digest(args.main_input_root / PROFILE)
            if args.main_input_root is not None
            else None
        ),
        "checkout_default_uncompressed_exists": (repo_root / PROFILE).is_file(),
        "checkout_default_columns_exists": (repo_root / COLUMNS).is_file(),
        "native_input_choice": native_choice,
        "input_annotation_policy": "pre-union; no read-time facultative canonicalization",
        "metadata": {
            "rows": len(metadata_raw),
            "valid_direction_rows": len(metadata),
            "studies": int(metadata["study_id"].nunique()),
            "exact_raw_feces_profiles": int(metadata_raw["body_site"].eq("Feces").sum()),
            "native_cleaned_exact_feces_profiles": len(stool),
            "stool_studies": int(stool["study_id"].nunique()),
            "host_species_counts": metadata_raw["host_species"].value_counts().to_dict(),
            "sequencing_type_counts": (
                metadata_raw["sequencing_type"].value_counts().to_dict()
                if "sequencing_type" in metadata_raw
                else "not present in this metadata file"
            ),
            "multisite_feces_profiles_excluded": int(
                metadata["body_site"].str.contains("Feces", regex=False).sum() - len(stool)
            ),
        },
        "column_metadata": {
            "rows": len(columns),
            "oxygen_rows": records(columns.loc[columns["column_id"].isin(OXYGEN_IDS)]),
        },
        "trait_list": trait_list.provenance(),
        "reference_traits": reference_traits,
        "reference_only_traits": sorted(set(reference_traits) - set(trait_list.trait_ids)),
        "reference_anaerobic_column": "METPO:1000603",
        "standardization": {
            "population": "all metadata-aligned profiles before exact raw Feces filter",
            "profile_count": len(metadata),
            "retained_traits_in_denominator": len(raw_counts.columns),
            "formula": (
                "z(log1p(1000 * trait_count / sum(all retained non-MediaDive counts))), ddof=0"
            ),
            "z_mean_by_trait": zscores.mean().to_dict(),
            "z_population_sd_by_trait": zscores.std(ddof=0).to_dict(),
        },
        "group_rule": (
            "any existing support-passing cell; descending distinct stool study count, "
            "then ascending group name; first six"
        ),
        "support_rule": {
            "minimum_positive_profiles_in_either_direction": heatmap.MIN_POSITIVE_PROFILES,
            "minimum_studies": heatmap.MIN_STUDIES,
        },
        "groups": records(groups),
        "selected_groups": records(selected),
        "reference_group_comparisons": comparisons,
        "reference_cell_comparisons": cell_comparisons,
    }
    output: Path = args.output_dir
    output.mkdir(parents=True, exist_ok=False)
    for name, table in six_tables.items():
        if name.endswith("_contrasts.tsv"):
            heatmap.write_tsv(table, output / name)
        else:
            table.to_csv(output / name, sep="\t", index=False)
    result["default_six_baseline"] = {
        name: {"sha256": digest(output / name), "rows": len(table)}
        for name, table in six_tables.items()
    }
    if args.compare_baseline is not None:
        result["default_six_comparison"] = compare_baseline(args.compare_baseline, output, result)
    (output / "input_diagnostics.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    cells.to_csv(output / "native_exact_feces_cells.tsv", sep="\t", index=False)
    groups.to_csv(output / "native_group_ranking.tsv", sep="\t", index=False)
    lines = [
        "# Issue 1093 input review",
        "",
        f"Native code: `{repo.head.commit.hexsha}`. All inputs were unchanged after inspection.",
        "Native input choice: " + native_choice + ". No regeneration or figures.",
        "",
        "## Input population",
        "",
    ]
    for name, inventory in inventories.items():
        lines.append(
            f"- `{name}`: {inventory['profile_count']} profiles, "
            f"{inventory['study_count']} studies, "
            f"{inventory['retained_non_mediadive_traits']} retained non-MediaDive traits."
        )
    lines.extend(
        [
            "",
            "## Reference group counts",
            "",
            "| Group | Studies old → current | Increased lists old → current "
            "| Decreased lists old → current |",
            "|---|---:|---:|---:|",
        ]
    )
    for item in comparisons:
        lines.append(
            f"| {item['group']} | {item['old_studies']} → {item['new_studies']} "
            f"| {item['old_increased_profiles']} → {item['new_increased_profiles']} "
            f"| {item['old_decreased_profiles']} → {item['new_decreased_profiles']} |"
        )
    lines.extend(
        [
            "",
            "The current declaration contains "
            + str(len(trait_list.trait_ids))
            + " traits. The historical reference has "
            + str(len(reference_traits))
            + "; reference-only terms: "
            + ", ".join(result["reference_only_traits"])
            + ".",
            "The reference's anaerobic column is `METPO:1000603`.",
            "The diagnostic JSON records exact hashes, oxygen columns, support and ranking, "
            "and every comparable cell.",
            "Z-scores use all metadata-aligned profiles before filtering exact raw "
            "`body_site == Feces`.",
            "",
        ]
    )
    (output / "FINDINGS.md").write_text("\n".join(lines), encoding="utf-8")
    print(
        json.dumps(
            {
                "profiles": len(metadata),
                "stool_profiles": len(stool),
                "selected_groups": records(selected),
            },
            indent=2,
        )
    )


def main() -> None:
    """Accept explicit read-only source roots and one fresh output directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main-input-root", type=Path)
    parser.add_argument("--columns-file", type=Path)
    parser.add_argument("--f-root", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--input-basis",
        choices=("main-defaults", "tracked-canonical"),
        default="tracked-canonical",
    )
    parser.add_argument("--six-baseline", action="store_true")
    parser.add_argument("--compare-baseline", type=Path)
    parser.add_argument("--metadata-causes-only", action="store_true")
    args = parser.parse_args()
    if args.metadata_causes_only:
        metadata_causes(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
