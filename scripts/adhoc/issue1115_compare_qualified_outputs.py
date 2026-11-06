"""Compare completed trial outputs with pinned F outputs for #1115 (2026-10-28).

Inputs: pinned trial/parent/baseline bundles, bound experiments/list, three Git SHAs.
Outputs: a fresh durable JSON parity report. The actual trial command performs
full source verification; this helper only checks completed outputs and bindings.
It never reconstructs or reruns Delta-z or the derived statistical tables.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
from pathlib import Path
from typing import Any

import pandas as pd
from git import Blob, Repo

from mtp.multiple_testing import canonical_sha256
from mtp.workflow.output_roots import temporary_root

PAIR_LEDGER = "composition_pair_delta_z.tsv"
STUDY_LEDGER = "composition_study_delta_z.tsv"
EVIDENCE = "per_trait_evidence.tsv"
TABLE_KEYS = {
    "trait_summary.tsv": ["trait"],
    "quadrant_pairs.tsv": ["study_id", "experiment_id"],
    "quadrant_summary.tsv": ["population"],
    "stool_disease_shifts.tsv": ["body_site", "condition", "trait"],
    "stool_study_sensitivity.tsv": ["body_site", "condition", "trait"],
    "site_disease_shifts.tsv": ["body_site", "condition", "trait"],
}
LABEL_FIELDS = {"trait_label", "fa_label", "oa_label"}
REMOVED_TRAIT = "METPO:1000605"
COMPLETE_COVERAGE = "complete"
HISTORICAL_OFF_SUMMARY_COVERAGE = "historical-off-summary"
HISTORICAL_PHENOTYPE_PRE_1115_COVERAGE = "historical-phenotype-pre-1115"
COVERAGE_MODES = (
    COMPLETE_COVERAGE,
    HISTORICAL_OFF_SUMMARY_COVERAGE,
    HISTORICAL_PHENOTYPE_PRE_1115_COVERAGE,
)
HISTORICAL_PHENOTYPE_PRE_1115_TABLES = tuple(
    name for name in TABLE_KEYS if name != "stool_study_sensitivity.tsv"
)
COMPLETE_BASELINE_OUTPUTS = frozenset({PAIR_LEDGER, STUDY_LEDGER, *TABLE_KEYS})
HISTORICAL_OFF_SUMMARY_OUTPUTS = frozenset({PAIR_LEDGER, STUDY_LEDGER, "trait_summary.tsv"})
HISTORICAL_PHENOTYPE_PRE_1115_OUTPUTS = frozenset(
    {PAIR_LEDGER, STUDY_LEDGER, *HISTORICAL_PHENOTYPE_PRE_1115_TABLES}
)


def validate_baseline_coverage(
    coverage_mode: Any,
    strain_rollup: str,
    baseline: dict[str, Any],
    baseline_outputs: dict[str, bytes],
) -> str:
    if not isinstance(coverage_mode, str) or coverage_mode not in COVERAGE_MODES:
        raise ValueError(
            "Explicit baseline coverage mode must be 'complete', "
            "'historical-off-summary' or 'historical-phenotype-pre-1115'"
        )
    observed = set(baseline_outputs)
    if coverage_mode == COMPLETE_COVERAGE:
        if baseline.get("summary_only") is not False:
            raise ValueError(
                "Complete baseline coverage requires baseline summary_only to be False"
            )
        missing = sorted(COMPLETE_BASELINE_OUTPUTS - observed)
        unexpected = sorted(observed - COMPLETE_BASELINE_OUTPUTS)
        if missing or unexpected:
            raise ValueError(
                "Complete baseline coverage missing or unexpected tables: "
                f"missing={missing}, unexpected={unexpected}"
            )
        return COMPLETE_COVERAGE
    if coverage_mode == HISTORICAL_PHENOTYPE_PRE_1115_COVERAGE:
        if strain_rollup != "phenotype":
            raise ValueError(
                "historical-phenotype-pre-1115 baseline coverage is restricted to phenotype scope"
            )
        if baseline.get("summary_only") is not False:
            raise ValueError(
                "historical-phenotype-pre-1115 baseline coverage requires "
                "baseline summary_only to be False"
            )
        missing = sorted(HISTORICAL_PHENOTYPE_PRE_1115_OUTPUTS - observed)
        unexpected = sorted(observed - HISTORICAL_PHENOTYPE_PRE_1115_OUTPUTS)
        if missing or unexpected:
            raise ValueError(
                "historical-phenotype-pre-1115 baseline coverage requires "
                "the five pre-#1115 derived tables and both ledgers only: "
                f"missing={missing}, unexpected={unexpected}"
            )
        return HISTORICAL_PHENOTYPE_PRE_1115_COVERAGE
    if strain_rollup != "off":
        raise ValueError("historical-off-summary baseline coverage is restricted to off scope")
    if baseline.get("summary_only") is not True:
        raise ValueError(
            "historical-off-summary baseline coverage requires baseline summary_only to be True"
        )
    missing = sorted(HISTORICAL_OFF_SUMMARY_OUTPUTS - observed)
    unexpected = sorted(observed - HISTORICAL_OFF_SUMMARY_OUTPUTS)
    if missing or unexpected:
        raise ValueError(
            "historical-off-summary baseline coverage requires trait_summary.tsv "
            f"and both ledgers only: missing={missing}, unexpected={unexpected}"
        )
    return HISTORICAL_OFF_SUMMARY_COVERAGE


def identity(raw: bytes) -> dict[str, str | int]:
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def snapshot(path: Path, snapshots: dict[Path, bytes]) -> bytes:
    raw = path.read_bytes()
    if path in snapshots and snapshots[path] != raw:
        raise ValueError(f"Input changed during qualification: {path}")
    snapshots[path] = raw
    return raw


def bundle(
    path: Path, pin: str, snapshots: dict[Path, bytes]
) -> tuple[dict[str, Any], dict[str, bytes]]:
    raw = snapshot(path, snapshots)
    if hashlib.sha256(raw).hexdigest() != pin:
        raise ValueError(f"Manifest differs from explicit pin: {path}")
    record: dict[str, Any] = json.loads(raw)
    outputs = {}
    for name, expected in record["outputs"].items():
        if Path(name).name != name:
            raise ValueError("Expected flat completed output inventory")
        content = snapshot(path.parent / name, snapshots)
        if identity(content) != expected:
            raise ValueError(f"Output differs from pinned inventory: {name}")
        outputs[name] = content
    actual = {item.name for item in path.parent.iterdir()}
    if actual != {path.name, *outputs}:
        raise ValueError(f"Completed bundle has unexpected/missing files: {path.parent}")
    return record, outputs


def bind_input(name: str, raw: bytes, record: dict[str, Any]) -> None:
    expected = identity(raw)
    named = [value for key, value in record["inputs"].items() if key.startswith(f"{name}@")]
    if named != [expected]:
        raise ValueError(f"Input bytes differ from the pinned inventory: {name}")


def frame(raw: bytes, keys: list[str], *, text: bool = False) -> pd.DataFrame:
    header = next(csv.reader(io.StringIO(raw.decode("utf-8")), delimiter="\t"))
    if not header or len(set(header)) != len(header) or not all(header):
        raise ValueError("Empty or duplicate TSV columns")
    result: pd.DataFrame = pd.read_csv(
        io.BytesIO(raw),
        sep="\t",
        dtype=str if text else dict.fromkeys(keys, str),
        keep_default_na=False,
        float_precision="round_trip",
    )
    if result.empty or result.duplicated(keys).any():
        raise ValueError(f"Empty or duplicate TSV identities: {keys}")
    return result


def git_blob(repository: Repo, commit_sha: str, name: str) -> bytes:
    item = repository.commit(commit_sha).tree / name
    if not isinstance(item, Blob):
        raise ValueError(f"Expected a Git file blob: {name}")
    return bytes(item.data_stream.read())


def check_trial(current: dict[str, Any], args: argparse.Namespace, project: Path) -> dict[str, str]:
    commits = {
        "issue1115_commit": args.issue1115_sha,
        "pr1130_commit": args.pr1130_sha,
        "local_trial_merge_commit": args.trial_sha,
    }
    if any(re.fullmatch(r"[0-9a-f]{40}", value) is None for value in commits.values()):
        raise ValueError("All three commits must be explicit full Git SHAs")
    if current["code_commit"] != args.trial_sha or current["code_dirty"] is not False:
        raise ValueError("Completed trial producer must be the declared clean merge")
    required = {"src/mtp/delta_z_derived/tables.py", "src/mtp/delta_z_derived/stage.py", "uv.lock"}
    if not required.issubset(current["code_inventory"]):
        raise ValueError("Trial code inventory omits required producer files")
    with Repo(project) as repository:
        parents = [parent.hexsha for parent in repository.commit(args.trial_sha).parents]
        if len(parents) != 2 or set(parents) != {args.issue1115_sha, args.pr1130_sha}:
            raise ValueError("Trial merge must have precisely the declared #1115 and #1130 parents")
        expected_inventory = {
            str(item.path): hashlib.sha256(item.data_stream.read()).hexdigest()
            for item in repository.commit(args.trial_sha).tree.traverse()
            if isinstance(item, Blob)
            and str(item.path).startswith("src/mtp/")
            and str(item.path).endswith(".py")
        }
        for name in (
            "methods/three_test_table.md",
            "methods/delta_z_derived.md",
            "uv.lock",
            "pyproject.toml",
            "Justfile",
        ):
            expected_inventory[name] = hashlib.sha256(
                git_blob(repository, args.trial_sha, name)
            ).hexdigest()
    if current["code_inventory"] != expected_inventory:
        raise ValueError("Trial producer inventory differs from the complete Git tree inventory")
    if current["implementation_sha256"] != canonical_sha256(expected_inventory):
        raise ValueError("Trial implementation digest differs from the complete producer inventory")
    if current["methods_sha256"] != expected_inventory["methods/delta_z_derived.md"]:
        raise ValueError("Trial method digest differs from the derived method blob")
    if current["lock_sha256"] != expected_inventory["uv.lock"]:
        raise ValueError("Trial lock identity differs from code inventory")
    return commits


def compare(
    name: str,
    old_bytes: bytes,
    new_bytes: bytes,
    traits: list[str],
    labels: dict[str, str],
) -> dict[str, Any]:
    keys = TABLE_KEYS[name]
    old, new = frame(old_bytes, keys, text=True), frame(new_bytes, keys, text=True)
    removed_rows = []
    if name == "stool_disease_shifts.tsv":
        excluded = old.loc[~old.trait.isin(traits)]
        if set(excluded.trait) - {REMOVED_TRAIT}:
            raise ValueError("Historical stool output contains an undeclared removed trait")
        removed_rows = excluded[keys].to_dict("records")
        old = old.loc[old.trait.isin(traits)].reset_index(drop=True)
    if list(old[keys].itertuples(index=False, name=None)) != list(
        new[keys].itertuples(index=False, name=None)
    ):
        raise ValueError(f"Historical output identities/order differ: {name}")
    if not set(old.columns).issubset(new.columns):
        raise ValueError(f"Historical fields were removed: {name}")
    changes, unexpected = [], []
    numeric_cells = 0
    for index in range(len(old)):
        old_row, new_row = old.iloc[index], new.iloc[index]
        key = {field: str(old_row[field]) for field in keys}
        for column in old.columns:
            before, after = str(old_row[column]), str(new_row[column])
            if before == after:
                continue
            change = {"row": key, "field": column, "old": before, "new": after}
            if (
                column == "wilcoxon_method"
                and before == "auto"
                and after in {"exact", "permutation", "asymptotic", "not_tested"}
            ):
                changes.append({**change, "kind": "resolved_auto_method_metadata"})
                continue
            if column in LABEL_FIELDS:
                trait_field = {
                    "trait_label": "trait",
                    "fa_label": "fa_trait",
                    "oa_label": "oa_trait",
                }[column]
                if after == labels[str(new_row[trait_field])]:
                    changes.append({**change, "kind": "label_from_pinned_parent_evidence"})
                    continue
            try:
                old_value, new_value = float(before), float(after)
            except ValueError:
                unexpected.append(change)
                continue
            numeric_cells += 1
            if old_value != new_value and not (math.isnan(old_value) and math.isnan(new_value)):
                unexpected.append({**change, "kind": "numeric_difference"})
    return {
        "old_rows_compared": len(old),
        "new_rows": len(new),
        "old_fields_compared": list(old.columns),
        "added_fields": [column for column in new.columns if column not in old.columns],
        "removed_605_stool_rows": removed_rows,
        "intentional_changes": changes,
        "numeric_format_only_differences": numeric_cells
        - sum(change.get("kind") == "numeric_difference" for change in unexpected),
        "unexpected_changes": unexpected,
        "numeric_comparison_rule": "exact_parsed_float_equality_no_tolerance",
    }


def check_additive_counts(
    outputs: dict[str, bytes],
    pairs: pd.DataFrame,
    experiments: pd.DataFrame,
    traits: list[str],
    top_diseases: int,
) -> None:
    """Bind new-output identities/counts to source metadata without retesting effects."""
    pair_keys = TABLE_KEYS["quadrant_pairs.tsv"]
    quadrants = frame(outputs["quadrant_pairs.tsv"], pair_keys, text=True).set_index(pair_keys)
    if set(quadrants.index) != set(pairs.index):
        raise ValueError("Quadrant identities differ from parent comparison keys")
    metadata = experiments.set_index(pair_keys).loc[pairs.index]
    for field in ("body_site", "condition"):
        if not quadrants.reindex(metadata.index)[field].equals(metadata[field]):
            raise ValueError(f"Quadrant {field} differs from bound experiment metadata")
    summary = frame(outputs["trait_summary.tsv"], ["trait"], text=True)
    if list(summary.trait) != list(pairs.columns):
        raise ValueError("Summary identities differ from parent retained traits")
    if isinstance(top_diseases, bool) or not isinstance(top_diseases, int) or top_diseases < 1:
        raise ValueError("Invalid declared number of stool conditions")
    stool = metadata.loc[metadata.body_site == "Feces"]
    groups = {str(condition): group for condition, group in stool.groupby("condition", sort=True)}
    selected = sorted(groups, key=lambda condition: (-len(groups[condition]), condition))[
        :top_diseases
    ]
    expected_cells = [("Feces", condition, trait) for condition in selected for trait in traits]
    for name in ("stool_disease_shifts.tsv", "stool_study_sensitivity.tsv"):
        cells = frame(outputs[name], TABLE_KEYS[name], text=True)
        observed_cells = list(cells[TABLE_KEYS[name]].itertuples(index=False, name=None))
        if observed_cells != expected_cells:
            raise ValueError(f"Stool table omits, adds or reorders declared cells: {name}")
        for _, row in cells.iterrows():
            group = groups[row["condition"]]
            if int(row["n_comparisons"]) != len(group) or int(row["n_studies"]) != len(
                set(group.index.get_level_values("study_id"))
            ):
                raise ValueError(f"Stool cell counts differ from bound experiment metadata: {name}")


def qualify(args: argparse.Namespace) -> dict[str, Any]:
    output = args.output.resolve()
    sources = [
        args.current_manifest,
        args.parent_manifest,
        args.baseline_manifest,
        args.experiments,
        args.trait_list,
    ]
    if output.exists() or temporary_root(output) is not None:
        raise ValueError("Output must be a fresh durable directory")
    for path in sources:
        root = path.resolve().parent
        if output == root or root in output.parents or output in root.parents:
            raise ValueError("Output must be disjoint from input directories")
    snapshots: dict[Path, bytes] = {}
    parent, parent_outputs = bundle(args.parent_manifest, args.parent_sha256, snapshots)
    baseline, baseline_outputs = bundle(args.baseline_manifest, args.baseline_sha256, snapshots)
    current, current_outputs = bundle(args.current_manifest, args.current_sha256, snapshots)
    if set(parent_outputs) != {PAIR_LEDGER, STUDY_LEDGER, EVIDENCE}:
        raise ValueError("Parent must contain exactly the three completed tables")
    for field in (
        "strain_rollup",
        "contrast_population",
        "direction_basis",
        "retained_traits",
        "source_commits",
        "source_run_ids",
        "core_run_id",
        "comparison_run_id",
        "analysis_bases",
    ):
        if parent[field] != baseline[field] or parent[field] != current[field]:
            raise ValueError(f"Parent, baseline and current bundle differ: {field}")
    if current["parent_manifest_sha256"] != args.parent_sha256:
        raise ValueError("Trial bundle is not bound to the explicitly pinned parent")
    if (
        current["parent_run_id"] != parent["run_id"]
        or current["parent_input_inventory"] != parent["inputs"]
    ):
        raise ValueError("Trial parent identity or input inventory differs")
    if current["summary_only"] is not False:
        raise ValueError("Trial qualification requires all derived tables")
    if set(current_outputs) != {PAIR_LEDGER, STUDY_LEDGER, *TABLE_KEYS}:
        raise ValueError("Trial output inventory differs from expected full bundle")
    if parent["strain_rollup"] not in {"off", "phenotype"}:
        raise ValueError("Expected archived off or phenotype scope")
    coverage_mode = validate_baseline_coverage(
        getattr(args, "baseline_coverage", None),
        parent["strain_rollup"],
        baseline,
        baseline_outputs,
    )
    for name in (PAIR_LEDGER, STUDY_LEDGER):
        if (
            baseline_outputs[name] != parent_outputs[name]
            or current_outputs[name] != parent_outputs[name]
        ):
            raise ValueError(f"Parent, current and historical baseline ledger differ: {name}")
    bind_input(args.parent_manifest.name, snapshots[args.parent_manifest], current)
    for name, raw in parent_outputs.items():
        bind_input(name, raw, current)
    experiments_raw = snapshot(args.experiments, snapshots)
    for record in (parent, baseline, current):
        bind_input("experiments.tsv", experiments_raw, record)
    list_raw = snapshot(args.trait_list, snapshots)
    traits = [
        value
        for line in list_raw.decode("utf-8").splitlines()
        if (value := line.strip()) and not value.startswith("#")
    ]
    project = Path(__file__).resolve().parents[2]
    commits = check_trial(current, args, project)
    helper = Path(__file__).resolve()
    helper_identity = identity(snapshot(helper, snapshots))
    with Repo(project) as repository:
        if (
            git_blob(repository, args.issue1115_sha, str(helper.relative_to(project)))
            != snapshots[helper]
        ):
            raise ValueError("Comparison helper differs from the declared #1115 commit")
    if not traits or len(traits) != len(set(traits)):
        raise ValueError("Current trait list must be nonempty and unique")
    bind_input(args.trait_list.name, list_raw, current)
    if (
        current["trait_list"]["trait_ids"] != traits
        or current["trait_list"]["trait_list_sha256"] != identity(list_raw)["sha256"]
    ):
        raise ValueError("Trial declared trait list differs from explicitly supplied bytes")
    list_name = "methods/trait_lists/headline_ten_phenotype.txt"
    with Repo(project) as repository:
        if git_blob(repository, args.issue1115_sha, list_name) != list_raw:
            raise ValueError("Trait list differs from the declared #1115 commit")
    pairs = frame(parent_outputs[PAIR_LEDGER], ["study_id", "experiment_id"]).set_index(
        ["study_id", "experiment_id"]
    )
    studies = frame(parent_outputs[STUDY_LEDGER], ["study_id"]).set_index("study_id")
    if (
        list(pairs.columns) != parent["retained_traits"]
        or list(studies.columns) != parent["retained_traits"]
    ):
        raise ValueError("Ledger headers differ from pinned retained catalog")
    evidence = frame(parent_outputs[EVIDENCE], ["trait"])
    experiments = frame(experiments_raw, ["study_id", "experiment_id"], text=True)
    generated = {name: current_outputs[name] for name in TABLE_KEYS}
    check_additive_counts(
        generated, pairs, experiments, traits, current["policies"]["row3_top_diseases"]
    )
    labels = dict(zip(evidence.trait, evidence.trait_label, strict=True))
    comparisons = {
        name: compare(name, baseline_outputs[name], raw, traits, labels)
        for name, raw in generated.items()
        if name in baseline_outputs
    }
    unhandled = set(baseline_outputs) - {PAIR_LEDGER, STUDY_LEDGER, *comparisons}
    if unhandled:
        raise ValueError(f"Uncompared historical outputs: {sorted(unhandled)}")
    new_tables_without_old_reference = [name for name in generated if name not in comparisons]
    is_complete = coverage_mode == COMPLETE_COVERAGE
    if is_complete and (set(comparisons) != set(TABLE_KEYS) or new_tables_without_old_reference):
        raise ValueError("Complete baseline coverage did not compare all derived tables")
    if coverage_mode == HISTORICAL_PHENOTYPE_PRE_1115_COVERAGE and (
        set(comparisons) != set(HISTORICAL_PHENOTYPE_PRE_1115_TABLES)
        or new_tables_without_old_reference != ["stool_study_sensitivity.tsv"]
    ):
        raise ValueError(
            "historical-phenotype-pre-1115 baseline coverage must compare "
            "the five pre-#1115 derived tables only"
        )
    if coverage_mode == HISTORICAL_OFF_SUMMARY_COVERAGE and (
        set(comparisons) != {"trait_summary.tsv"}
        or new_tables_without_old_reference
        != [name for name in TABLE_KEYS if name != "trait_summary.tsv"]
    ):
        raise ValueError(
            "historical-off-summary baseline coverage must compare trait_summary.tsv only"
        )
    summary = frame(generated["trait_summary.tsv"], ["trait"], text=True)
    flags = list(summary.sdz_active_at_least_50)
    if set(flags) - {"True", "False"}:
        raise ValueError("Exported active-study flags must be Boolean text")
    counts = {
        "comparisons": len(pairs),
        "studies": len(studies),
        "traits": len(studies.columns),
        "traits_with_at_least_50_active_studies": flags.count("True"),
    }
    if counts != current["counts"]:
        raise ValueError("Current manifest counts differ from supplied tables")
    for path, raw in snapshots.items():
        if path.read_bytes() != raw:
            raise ValueError(f"Input or implementation changed: {path}")
    status = (
        "FAIL" if any(value["unexpected_changes"] for value in comparisons.values()) else "PASS"
    )
    result = {
        "issue": 1115,
        "description": (
            "complete parity comparison of completed verified trial outputs "
            "against pinned historical outputs"
            if is_complete
            else (
                "partial-reference comparison of completed verified trial outputs "
                "against pinned historical pre-#1115 phenotype baseline"
                if coverage_mode == HISTORICAL_PHENOTYPE_PRE_1115_COVERAGE
                else (
                    "partial-reference comparison of completed verified trial outputs "
                    "against pinned historical off summary baseline"
                )
            )
        ),
        "baseline_coverage_mode": coverage_mode,
        "comparison_scope": ("complete_parity" if is_complete else "partial_reference_comparison"),
        "complete_parity": is_complete and status == "PASS",
        "partial_reference_comparison": not is_complete,
        "full_upstream_verification_performed_by_this_helper": False,
        "qualification_source": "actual delta-z-derived command on local trial merge",
        "merge_commits": commits,
        "current_run_id": current["run_id"],
        "current_manifest_sha256": args.current_sha256,
        "delta_z_recomputed": False,
        "status": status,
        "strain_rollup": parent["strain_rollup"],
        "parent_run_id": parent["run_id"],
        "parent_manifest_sha256": args.parent_sha256,
        "baseline_run_id": baseline["run_id"],
        "baseline_manifest_sha256": args.baseline_sha256,
        "baseline_parent_manifest_sha256": baseline["parent_manifest_sha256"],
        "parent_current_and_baseline_pair_study_bytes_identical": True,
        "inputs_and_code_unchanged": True,
        "additive_output_identities_and_counts_match_parent_metadata": True,
        "input_inventory": {str(path.resolve()): identity(raw) for path, raw in snapshots.items()},
        "helper_identity": helper_identity,
        "trial_code_inventory": current["code_inventory"],
        "trial_environment": current["environment"],
        "declared_traits": traits,
        "counts": counts,
        "comparisons": comparisons,
        "new_tables_without_old_reference": new_tables_without_old_reference,
        "outputs": {name: identity(raw) for name, raw in generated.items()},
        "output_rows": {
            name: len(frame(raw, TABLE_KEYS[name], text=True)) for name, raw in generated.items()
        },
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "table_parity.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for option in (
        "current-manifest",
        "parent-manifest",
        "baseline-manifest",
        "experiments",
        "trait-list",
        "output",
    ):
        parser.add_argument(f"--{option}", type=Path, required=True)
    parser.add_argument("--parent-sha256", required=True)
    parser.add_argument("--baseline-sha256", required=True)
    parser.add_argument("--current-sha256", required=True)
    parser.add_argument("--issue1115-sha", required=True)
    parser.add_argument("--pr1130-sha", required=True)
    parser.add_argument("--trial-sha", required=True)
    parser.add_argument("--baseline-coverage", choices=COVERAGE_MODES, required=True)
    result = qualify(parser.parse_args())
    print(json.dumps(result, indent=2))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
