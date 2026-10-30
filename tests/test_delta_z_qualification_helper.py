"""Reject incomplete producer evidence and unexpected baseline coverage (#1156, #1187, #1213)."""

from __future__ import annotations

import argparse
import hashlib
import json
import runpy
import sys
import tempfile
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

SCRIPT_PATH = Path(__file__).parents[1] / "scripts/adhoc/issue1115_compare_qualified_outputs.py"
PAIR_LEDGER = "composition_pair_delta_z.tsv"
STUDY_LEDGER = "composition_study_delta_z.tsv"
EVIDENCE = "per_trait_evidence.tsv"
TABLE_KEYS: dict[str, list[str]] = {
    "trait_summary.tsv": ["trait"],
    "quadrant_pairs.tsv": ["study_id", "experiment_id"],
    "quadrant_summary.tsv": ["population"],
    "stool_disease_shifts.tsv": ["body_site", "condition", "trait"],
    "stool_study_sensitivity.tsv": ["body_site", "condition", "trait"],
    "site_disease_shifts.tsv": ["body_site", "condition", "trait"],
}
PRE_1115_PHENOTYPE_TABLES: tuple[str, ...] = tuple(
    name for name in TABLE_KEYS if name != "stool_study_sensitivity.tsv"
)
EXPLICIT_COVERAGE_MODE_MATCH = (
    r"^Explicit baseline coverage mode must be 'complete', "
    r"'historical-off-summary' or 'historical-phenotype-pre-1115'$"
)
TRAIT_LIST_REL = "methods/trait_lists/headline_ten_phenotype.txt"
TRAIT_LIST_BYTES = b"METPO:1000608\nMETPO:1000607\n"
PRODUCER_FILES = {
    "src/mtp/delta_z_derived/tables.py": b"tables = 1\n",
    "src/mtp/delta_z_derived/stage.py": b"stage = 1\n",
    "src/mtp/extra.py": b"extra = 1\n",
    "src/mtp/nested/dependency.py": b"dependency = 1\n",
    "methods/three_test_table.md": b"Parent method\n",
    "methods/delta_z_derived.md": b"Derived method\n",
    "uv.lock": b"version = 1\n",
    "pyproject.toml": b"[project]\n",
    "Justfile": b"default:\n",
}


def _identity(raw: bytes) -> dict[str, str | int]:
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def digest_inventory(inventory: dict[str, str]) -> str:
    raw = json.dumps(inventory, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


@pytest.fixture(scope="module")
def _shared_trial_repo() -> Iterator[tuple[dict[str, Any], argparse.Namespace, Path]]:
    scratch_Base = Path(__file__).resolve().parents[1] / ".tmp"
    scratch_root = scratch_Base / "delta-z-qualification-tests"
    scratch_root.mkdir(parents=True, exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(prefix="repo-", dir=scratch_root) as directory:
            repo_dir = Path(directory)
            files = {
                **PRODUCER_FILES,
                "src/other/non_producer.py": b"outside = 1\n",
                "scripts/adhoc/issue1115_compare_qualified_outputs.py": SCRIPT_PATH.read_bytes(),
                TRAIT_LIST_REL: TRAIT_LIST_BYTES,
            }
            for name, raw in files.items():
                path = repo_dir / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(raw)
            with pytest.MonkeyPatch.context() as mp:
                mp.setenv("GIT_CONFIG_GLOBAL", "/dev/null")
                mp.setenv("GIT_CONFIG_NOSYSTEM", "1")
                from git import Actor, Repo

                author = Actor("Qualification fixture", "qualification@example.invalid")
                with Repo.init(repo_dir) as repository:
                    repository.index.add(list(files))
                    base = repository.index.commit("base", author=author, committer=author)
                    issue = repository.index.commit(
                        "1115", parent_commits=[base], head=False, author=author, committer=author
                    )
                    policy = repository.index.commit(
                        "1130", parent_commits=[base], head=False, author=author, committer=author
                    )
                    merged = repository.index.commit(
                        "trial", parent_commits=[issue, policy], author=author, committer=author
                    )
            inventory = {
                name: hashlib.sha256(raw).hexdigest() for name, raw in PRODUCER_FILES.items()
            }
            record = {
                "code_commit": merged.hexsha,
                "code_dirty": False,
                "code_inventory": inventory,
                "implementation_sha256": digest_inventory(inventory),
                "methods_sha256": inventory["methods/delta_z_derived.md"],
                "lock_sha256": inventory["uv.lock"],
            }
            args = argparse.Namespace(
                issue1115_sha=issue.hexsha, pr1130_sha=policy.hexsha, trial_sha=merged.hexsha
            )
            yield record, args, repo_dir
    finally:
        if scratch_root.exists() and not any(scratch_root.iterdir()):
            scratch_root.rmdir()
        if scratch_Base.exists() and not any(scratch_Base.iterdir()):
            scratch_Base.rmdir()


@pytest.fixture
def trial(
    _shared_trial_repo: tuple[dict[str, Any], argparse.Namespace, Path],
    monkeypatch: pytest.MonkeyPatch,
) -> tuple[dict[str, Any], argparse.Namespace, Path]:
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", "/dev/null")
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    record, args, repo_dir = _shared_trial_repo
    return (
        {
            **record,
            "code_inventory": dict(record["code_inventory"]),
        },
        argparse.Namespace(**vars(args)),
        repo_dir,
    )


def _write_bundle(
    bundle_dir: Path,
    manifest_name: str,
    record: dict[str, Any],
    files: dict[str, bytes],
) -> tuple[Path, str]:
    if bundle_dir.exists():
        for existing in bundle_dir.iterdir():
            existing.unlink()
    else:
        bundle_dir.mkdir(parents=True, exist_ok=True)
    for name, raw in files.items():
        (bundle_dir / name).write_bytes(raw)
    manifest_Record = {
        **record,
        "outputs": {name: _identity(raw) for name, raw in files.items()},
    }
    manifest_bytes = (json.dumps(manifest_Record, indent=2, sort_keys=True) + "\n").encode("utf-8")
    manifest_path = bundle_dir / manifest_name
    manifest_path.write_bytes(manifest_bytes)
    return manifest_path, hashlib.sha256(manifest_bytes).hexdigest()


def _load_repo_helper(repo_root: Path) -> dict[str, Any]:
    return runpy.run_path(str(repo_root / "scripts/adhoc/issue1115_compare_qualified_outputs.py"))


def _build_qualification_scenario(
    trial_fixture: tuple[dict[str, Any], argparse.Namespace, Path],
    *,
    strain_rollup: str = "off",
    baseline_ledger_names: tuple[str, ...] = (PAIR_LEDGER, STUDY_LEDGER),
    baseline_table_names: tuple[str, ...] = tuple(TABLE_KEYS),
    baseline_table_overrides: dict[str, bytes] | None = None,
    baseline_summary_only: bool = False,
    baseline_coverage: str | None = "complete",
) -> tuple[argparse.Namespace, Path, Path]:
    producer_record, commit_args, repo_root = trial_fixture
    cases_root = repo_root / "cases"
    cases_root.mkdir(parents=True, exist_ok=True)
    case_dir = Path(tempfile.mkdtemp(prefix="case-", dir=cases_root))
    work_root = case_dir / "qualification_inputs"
    work_root.mkdir(parents=True, exist_ok=True)

    pair_bytes = (
        b"study_id\texperiment_id\tMETPO:1000608\tMETPO:1000607\nstudy_1\texp_1\t0.25\t-0.25\n"
    )
    study_bytes = b"study_id\tMETPO:1000608\tMETPO:1000607\nstudy_1\t0.25\t-0.25\n"
    evidence_bytes = (
        b"trait\ttrait_label\n"
        b"METPO:1000608\tfacultatively aerobic\n"
        b"METPO:1000607\tobligately anaerobic\n"
    )
    experiments_bytes = (
        b"study_id\texperiment_id\tbody_site\tcondition\nstudy_1\texp_1\tFeces\tCRC\n"
    )
    experiments_dir = work_root / "inputs"
    experiments_dir.mkdir(parents=True, exist_ok=True)
    experiments_path = experiments_dir / "experiments.tsv"
    experiments_path.write_bytes(experiments_bytes)
    trait_list_path = repo_root / TRAIT_LIST_REL

    derived_tables: dict[str, bytes] = {
        "trait_summary.tsv": (
            b"trait\ttrait_label\tsdz_active_at_least_50\tmean_sdz\n"
            b"METPO:1000608\tfacultatively aerobic\tTrue\t0.25\n"
            b"METPO:1000607\tobligately anaerobic\tFalse\t-0.25\n"
        ),
        "quadrant_pairs.tsv": (
            b"study_id\texperiment_id\tbody_site\tcondition\tfa_trait\tfa_label\t"
            b"oa_trait\toa_label\tfa_dz\toa_dz\n"
            b"study_1\texp_1\tFeces\tCRC\tMETPO:1000608\tfacultatively aerobic\t"
            b"METPO:1000607\tobligately anaerobic\t0.25\t-0.25\n"
        ),
        "quadrant_summary.tsv": b"population\tn_pairs\tpearson_r\nall_pairs\t1\t-1.0\n",
        "stool_disease_shifts.tsv": (
            b"body_site\tcondition\ttrait\ttrait_label\tn_comparisons\tn_studies\t"
            b"wilcoxon_method\tmean_dz\n"
            b"Feces\tCRC\tMETPO:1000608\tfacultatively aerobic\t1\t1\texact\t0.25\n"
            b"Feces\tCRC\tMETPO:1000607\tobligately anaerobic\t1\t1\texact\t-0.25\n"
        ),
        "stool_study_sensitivity.tsv": (
            b"body_site\tcondition\ttrait\ttrait_label\tn_comparisons\tn_studies\t"
            b"wilcoxon_method\tmean_sdz\n"
            b"Feces\tCRC\tMETPO:1000608\tfacultatively aerobic\t1\t1\texact\t0.25\n"
            b"Feces\tCRC\tMETPO:1000607\tobligately anaerobic\t1\t1\texact\t-0.25\n"
        ),
        "site_disease_shifts.tsv": (
            b"body_site\tcondition\ttrait\ttrait_label\tn_comparisons\tn_studies\tmean_dz\n"
            b"Feces\tCRC\tMETPO:1000608\tfacultatively aerobic\t1\t1\t0.25\n"
            b"Feces\tCRC\tMETPO:1000607\tobligately anaerobic\t1\t1\t-0.25\n"
        ),
    }

    shared_meta = {
        "strain_rollup": strain_rollup,
        "contrast_population": "all_retained",
        "direction_basis": "case_minus_control",
        "retained_traits": ["METPO:1000608", "METPO:1000607"],
        "source_commits": {"core": "a" * 40},
        "source_run_ids": {"core": "core-1"},
        "core_run_id": "core-1",
        "comparison_run_id": "comp-1",
        "analysis_bases": {"jed": "log1p"},
    }
    parent_inputs = {
        f"experiments.tsv@sha256:{_identity(experiments_bytes)['sha256']}": _identity(
            experiments_bytes
        )
    }
    parent_manifest, parent_sha256 = _write_bundle(
        work_root / "parent_bundle",
        "three_test_run.json",
        {
            **shared_meta,
            "run_id": "three-test-parent-1",
            "inputs": parent_inputs,
        },
        {
            PAIR_LEDGER: pair_bytes,
            STUDY_LEDGER: study_bytes,
            EVIDENCE: evidence_bytes,
        },
    )
    parent_manifest_bytes = parent_manifest.read_bytes()

    derived_inputs = {
        **parent_inputs,
        f"three_test_run.json@sha256:{parent_sha256}": _identity(parent_manifest_bytes),
        f"{PAIR_LEDGER}@sha256:{_identity(pair_bytes)['sha256']}": _identity(pair_bytes),
        f"{STUDY_LEDGER}@sha256:{_identity(study_bytes)['sha256']}": _identity(study_bytes),
        f"{EVIDENCE}@sha256:{_identity(evidence_bytes)['sha256']}": _identity(evidence_bytes),
        f"{trait_list_path.name}@sha256:{_identity(TRAIT_LIST_BYTES)['sha256']}": _identity(
            TRAIT_LIST_BYTES
        ),
    }
    counts = {
        "comparisons": 1,
        "studies": 1,
        "traits": 2,
        "traits_with_at_least_50_active_studies": 1,
    }
    ledger_bytes_by_name = {
        PAIR_LEDGER: pair_bytes,
        STUDY_LEDGER: study_bytes,
    }
    baseline_tables = {
        **derived_tables,
        **(baseline_table_overrides or {}),
    }
    baseline_files = {
        **{name: ledger_bytes_by_name[name] for name in baseline_ledger_names},
        **{name: baseline_tables[name] for name in baseline_table_names},
    }
    baseline_manifest, baseline_sha256 = _write_bundle(
        work_root / "baseline_bundle",
        "delta_z_derived_run.json",
        {
            **shared_meta,
            "run_id": "derived-baseline-1",
            "parent_run_id": "three-test-parent-1",
            "parent_manifest_sha256": parent_sha256,
            "parent_input_inventory": parent_inputs,
            "summary_only": baseline_summary_only,
            "inputs": derived_inputs,
            "counts": counts,
        },
        baseline_files,
    )

    current_manifest, current_sha256 = _write_bundle(
        work_root / "current_bundle",
        "delta_z_derived_run.json",
        {
            **shared_meta,
            **producer_record,
            "run_id": "derived-current-1",
            "parent_run_id": "three-test-parent-1",
            "parent_manifest_sha256": parent_sha256,
            "parent_input_inventory": parent_inputs,
            "summary_only": False,
            "inputs": derived_inputs,
            "trait_list": {
                "trait_ids": ["METPO:1000608", "METPO:1000607"],
                "trait_list_sha256": _identity(TRAIT_LIST_BYTES)["sha256"],
            },
            "policies": {"row3_top_diseases": 1},
            "counts": counts,
            "environment": {"python": "3.11"},
        },
        {
            PAIR_LEDGER: pair_bytes,
            STUDY_LEDGER: study_bytes,
            **derived_tables,
        },
    )

    output_dir = case_dir / "qualification_output"
    namespace_kwargs: dict[str, Any] = {
        "current_manifest": current_manifest,
        "parent_manifest": parent_manifest,
        "baseline_manifest": baseline_manifest,
        "experiments": experiments_path,
        "trait_list": trait_list_path,
        "output": output_dir,
        "parent_sha256": parent_sha256,
        "baseline_sha256": baseline_sha256,
        "current_sha256": current_sha256,
        "issue1115_sha": commit_args.issue1115_sha,
        "pr1130_sha": commit_args.pr1130_sha,
        "trial_sha": commit_args.trial_sha,
    }
    if baseline_coverage is not None:
        namespace_kwargs["baseline_coverage"] = baseline_coverage
    return argparse.Namespace(**namespace_kwargs), repo_root, output_dir


def test_complete_trial_inventory_passes(
    trial: tuple[dict[str, Any], argparse.Namespace, Path],
) -> None:
    record, args, root = trial
    module = _load_repo_helper(root)
    result = module["check_trial"](record, args, root)
    assert result["local_trial_merge_commit"] == args.trial_sha


@pytest.mark.parametrize("corruption", ["omitted_module", "aggregate_digest", "method_digest"])
def test_incomplete_or_stale_producer_identity_rejected(
    trial: tuple[dict[str, Any], argparse.Namespace, Path], corruption: str
) -> None:
    record, args, root = trial
    module = _load_repo_helper(root)
    if corruption == "omitted_module":
        del record["code_inventory"]["src/mtp/extra.py"]
        record["implementation_sha256"] = digest_inventory(record["code_inventory"])
    else:
        field = "implementation_sha256" if corruption == "aggregate_digest" else "methods_sha256"
        record[field] = "0" * 64
    with pytest.raises(ValueError, match=r"inventory|implementation|method"):
        module["check_trial"](record, args, root)


@pytest.mark.parametrize("corruption", ["omitted_module", "extra_module", "changed_hash"])
def test_inventory_mismatch_with_full_implementation_digest_rejected(
    trial: tuple[dict[str, Any], argparse.Namespace, Path],
    corruption: str,
) -> None:
    record, args, root = trial
    module = _load_repo_helper(root)
    assert record["implementation_sha256"] == digest_inventory(record["code_inventory"])
    if corruption == "omitted_module":
        del record["code_inventory"]["src/mtp/extra.py"]
    elif corruption == "extra_module":
        record["code_inventory"]["src/mtp/not_in_tree.py"] = "0" * 64
    else:
        record["code_inventory"]["src/mtp/extra.py"] = "0" * 64
    # Keep the full digest valid so only exact inventory equality can reject this.
    assert record["implementation_sha256"] != digest_inventory(record["code_inventory"])
    with pytest.raises(
        ValueError,
        match=r"^Trial producer inventory differs from the complete Git tree inventory$",
    ):
        module["check_trial"](record, args, root)


@pytest.mark.parametrize("invalid_mode", [None, "partial", "auto"])
def test_qualify_rejects_missing_or_invalid_baseline_coverage_mode(
    trial: tuple[dict[str, Any], argparse.Namespace, Path],
    invalid_mode: str | None,
) -> None:
    args, repo_root, output_dir = _build_qualification_scenario(
        trial,
        baseline_coverage=invalid_mode,
    )
    module = _load_repo_helper(repo_root)
    with pytest.raises(ValueError, match=EXPLICIT_COVERAGE_MODE_MATCH):
        module["qualify"](args)
    assert not output_dir.exists()


def test_main_requires_explicit_baseline_coverage_argument(
    trial: tuple[dict[str, Any], argparse.Namespace, Path],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    args, repo_root, output_dir = _build_qualification_scenario(trial, baseline_coverage=None)
    module = _load_repo_helper(repo_root)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(repo_root / "scripts/adhoc/issue1115_compare_qualified_outputs.py"),
            "--current-manifest",
            str(args.current_manifest),
            "--parent-manifest",
            str(args.parent_manifest),
            "--baseline-manifest",
            str(args.baseline_manifest),
            "--experiments",
            str(args.experiments),
            "--trait-list",
            str(args.trait_list),
            "--output",
            str(args.output),
            "--parent-sha256",
            args.parent_sha256,
            "--baseline-sha256",
            args.baseline_sha256,
            "--current-sha256",
            args.current_sha256,
            "--issue1115-sha",
            args.issue1115_sha,
            "--pr1130-sha",
            args.pr1130_sha,
            "--trial-sha",
            args.trial_sha,
        ],
    )
    with pytest.raises(SystemExit) as exc_info:
        module["main"]()
    assert exc_info.value.code == 2
    assert not output_dir.exists()


@pytest.mark.parametrize(
    ("baseline_ledgers", "baseline_tables", "summary_only"),
    [
        ((PAIR_LEDGER, STUDY_LEDGER), (), False),
        ((PAIR_LEDGER, STUDY_LEDGER), ("trait_summary.tsv",), True),
        ((PAIR_LEDGER, STUDY_LEDGER), ("trait_summary.tsv",), False),
        (
            (PAIR_LEDGER, STUDY_LEDGER),
            tuple(name for name in TABLE_KEYS if name != "site_disease_shifts.tsv"),
            False,
        ),
        ((PAIR_LEDGER, STUDY_LEDGER), tuple(TABLE_KEYS), True),
        ((PAIR_LEDGER,), tuple(TABLE_KEYS), False),
        ((STUDY_LEDGER,), tuple(TABLE_KEYS), False),
    ],
)
def test_complete_coverage_mode_rejects_missing_derived_tables_or_summary_only_baseline(
    trial: tuple[dict[str, Any], argparse.Namespace, Path],
    baseline_ledgers: tuple[str, ...],
    baseline_tables: tuple[str, ...],
    summary_only: bool,
) -> None:
    args, repo_root, output_dir = _build_qualification_scenario(
        trial,
        strain_rollup="off",
        baseline_ledger_names=baseline_ledgers,
        baseline_table_names=baseline_tables,
        baseline_summary_only=summary_only,
        baseline_coverage="complete",
    )
    module = _load_repo_helper(repo_root)
    with pytest.raises(
        ValueError,
        match=(
            r"^Complete baseline coverage "
            r"(requires baseline summary_only to be False|"
            r"missing or unexpected tables: )"
        ),
    ):
        module["qualify"](args)
    assert not output_dir.exists()


@pytest.mark.parametrize(
    ("strain_rollup", "baseline_ledgers", "baseline_tables", "summary_only"),
    [
        ("phenotype", (PAIR_LEDGER, STUDY_LEDGER), ("trait_summary.tsv",), True),
        ("off", (PAIR_LEDGER, STUDY_LEDGER), (), True),
        ("off", (PAIR_LEDGER, STUDY_LEDGER), ("trait_summary.tsv", "quadrant_summary.tsv"), True),
        ("off", (PAIR_LEDGER, STUDY_LEDGER), tuple(TABLE_KEYS), True),
        ("off", (PAIR_LEDGER, STUDY_LEDGER), ("trait_summary.tsv",), False),
        ("off", (PAIR_LEDGER,), ("trait_summary.tsv",), True),
        ("off", (STUDY_LEDGER,), ("trait_summary.tsv",), True),
    ],
)
def test_historical_off_summary_mode_rejects_non_off_scope_or_unexpected_baseline_inventory(
    trial: tuple[dict[str, Any], argparse.Namespace, Path],
    strain_rollup: str,
    baseline_ledgers: tuple[str, ...],
    baseline_tables: tuple[str, ...],
    summary_only: bool,
) -> None:
    args, repo_root, output_dir = _build_qualification_scenario(
        trial,
        strain_rollup=strain_rollup,
        baseline_ledger_names=baseline_ledgers,
        baseline_table_names=baseline_tables,
        baseline_summary_only=summary_only,
        baseline_coverage="historical-off-summary",
    )
    module = _load_repo_helper(repo_root)
    with pytest.raises(
        ValueError,
        match=(
            r"^historical-off-summary baseline coverage "
            r"(is restricted to off scope|"
            r"requires baseline summary_only to be True|"
            r"requires trait_summary\.tsv and both ledgers only: )"
        ),
    ):
        module["qualify"](args)
    assert not output_dir.exists()


def test_complete_coverage_mode_labels_report_as_complete_parity(
    trial: tuple[dict[str, Any], argparse.Namespace, Path],
) -> None:
    args, repo_root, output_dir = _build_qualification_scenario(
        trial,
        strain_rollup="phenotype",
        baseline_table_names=tuple(TABLE_KEYS),
        baseline_summary_only=False,
        baseline_coverage="complete",
    )
    module = _load_repo_helper(repo_root)
    result = module["qualify"](args)
    assert result["status"] == "PASS"
    assert result["baseline_coverage_mode"] == "complete"
    assert result["comparison_scope"] == "complete_parity"
    assert result["complete_parity"] is True
    assert result["partial_reference_comparison"] is False
    assert "complete parity" in result["description"]
    assert set(result["comparisons"]) == set(TABLE_KEYS)
    assert result["new_tables_without_old_reference"] == []
    persisted = json.loads((output_dir / "table_parity.json").read_text(encoding="utf-8"))
    assert persisted == result


def test_historical_off_summary_mode_labels_report_as_partial_reference_comparison(
    trial: tuple[dict[str, Any], argparse.Namespace, Path],
) -> None:
    args, repo_root, output_dir = _build_qualification_scenario(
        trial,
        strain_rollup="off",
        baseline_table_names=("trait_summary.tsv",),
        baseline_summary_only=True,
        baseline_coverage="historical-off-summary",
    )
    module = _load_repo_helper(repo_root)
    result = module["qualify"](args)
    assert result["status"] == "PASS"
    assert result["baseline_coverage_mode"] == "historical-off-summary"
    assert result["comparison_scope"] == "partial_reference_comparison"
    assert result["complete_parity"] is False
    assert result["partial_reference_comparison"] is True
    assert "partial-reference comparison" in result["description"]
    assert set(result["comparisons"]) == {"trait_summary.tsv"}
    assert result["new_tables_without_old_reference"] == [
        name for name in TABLE_KEYS if name != "trait_summary.tsv"
    ]
    persisted = json.loads((output_dir / "table_parity.json").read_text(encoding="utf-8"))
    assert persisted == result


@pytest.mark.parametrize(
    ("strain_rollup", "baseline_ledgers", "baseline_tables", "summary_only"),
    [
        ("off", (PAIR_LEDGER, STUDY_LEDGER), PRE_1115_PHENOTYPE_TABLES, False),
        ("phenotype", (PAIR_LEDGER, STUDY_LEDGER), PRE_1115_PHENOTYPE_TABLES, True),
        ("phenotype", (PAIR_LEDGER, STUDY_LEDGER), tuple(TABLE_KEYS), False),
        (
            "phenotype",
            (PAIR_LEDGER, STUDY_LEDGER),
            tuple(name for name in PRE_1115_PHENOTYPE_TABLES if name != "site_disease_shifts.tsv"),
            False,
        ),
        ("phenotype", (PAIR_LEDGER, STUDY_LEDGER), ("trait_summary.tsv",), False),
        ("phenotype", (PAIR_LEDGER,), PRE_1115_PHENOTYPE_TABLES, False),
        ("phenotype", (STUDY_LEDGER,), PRE_1115_PHENOTYPE_TABLES, False),
    ],
)
def test_historical_phenotype_pre_1115_rejects_non_phenotype_or_unexpected_baseline(
    trial: tuple[dict[str, Any], argparse.Namespace, Path],
    strain_rollup: str,
    baseline_ledgers: tuple[str, ...],
    baseline_tables: tuple[str, ...],
    summary_only: bool,
) -> None:
    args, repo_root, output_dir = _build_qualification_scenario(
        trial,
        strain_rollup=strain_rollup,
        baseline_ledger_names=baseline_ledgers,
        baseline_table_names=baseline_tables,
        baseline_summary_only=summary_only,
        baseline_coverage="historical-phenotype-pre-1115",
    )
    module = _load_repo_helper(repo_root)
    with pytest.raises(
        ValueError,
        match=(
            r"^historical-phenotype-pre-1115 baseline coverage "
            r"(is restricted to phenotype scope|"
            r"requires baseline summary_only to be False|"
            r"requires the five pre-#1115 derived tables and both ledgers only: )"
        ),
    ):
        module["qualify"](args)
    assert not output_dir.exists()


def test_historical_phenotype_pre_1115_mode_labels_report_as_partial_reference_comparison(
    trial: tuple[dict[str, Any], argparse.Namespace, Path],
) -> None:
    args, repo_root, output_dir = _build_qualification_scenario(
        trial,
        strain_rollup="phenotype",
        baseline_table_names=PRE_1115_PHENOTYPE_TABLES,
        baseline_summary_only=False,
        baseline_coverage="historical-phenotype-pre-1115",
    )
    module = _load_repo_helper(repo_root)
    result = module["qualify"](args)
    assert result["status"] == "PASS"
    assert result["baseline_coverage_mode"] == "historical-phenotype-pre-1115"
    assert result["comparison_scope"] == "partial_reference_comparison"
    assert result["complete_parity"] is False
    assert result["partial_reference_comparison"] is True
    assert "partial-reference comparison" in result["description"]
    assert set(result["comparisons"]) == set(PRE_1115_PHENOTYPE_TABLES)
    assert result["new_tables_without_old_reference"] == ["stool_study_sensitivity.tsv"]
    persisted = json.loads((output_dir / "table_parity.json").read_text(encoding="utf-8"))
    assert persisted == result


def test_complete_coverage_mode_sets_complete_parity_false_and_exits_nonzero_on_fail(
    trial: tuple[dict[str, Any], argparse.Namespace, Path],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    args, repo_root, output_dir = _build_qualification_scenario(
        trial,
        strain_rollup="phenotype",
        baseline_table_names=tuple(TABLE_KEYS),
        baseline_table_overrides={
            "trait_summary.tsv": (
                b"trait\ttrait_label\tsdz_active_at_least_50\tmean_sdz\n"
                b"METPO:1000608\tfacultatively aerobic\tTrue\t0.50\n"
                b"METPO:1000607\tobligately anaerobic\tFalse\t-0.25\n"
            ),
        },
        baseline_summary_only=False,
        baseline_coverage="complete",
    )
    module = _load_repo_helper(repo_root)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(repo_root / "scripts/adhoc/issue1115_compare_qualified_outputs.py"),
            "--current-manifest",
            str(args.current_manifest),
            "--parent-manifest",
            str(args.parent_manifest),
            "--baseline-manifest",
            str(args.baseline_manifest),
            "--experiments",
            str(args.experiments),
            "--trait-list",
            str(args.trait_list),
            "--output",
            str(args.output),
            "--parent-sha256",
            args.parent_sha256,
            "--baseline-sha256",
            args.baseline_sha256,
            "--current-sha256",
            args.current_sha256,
            "--issue1115-sha",
            args.issue1115_sha,
            "--pr1130-sha",
            args.pr1130_sha,
            "--trial-sha",
            args.trial_sha,
            "--baseline-coverage",
            "complete",
        ],
    )
    with pytest.raises(SystemExit) as exc_info:
        module["main"]()
    assert exc_info.value.code == 1
    persisted = json.loads((output_dir / "table_parity.json").read_text(encoding="utf-8"))
    assert persisted["status"] == "FAIL"
    assert persisted["baseline_coverage_mode"] == "complete"
    assert persisted["comparison_scope"] == "complete_parity"
    assert persisted["complete_parity"] is False
    assert persisted["partial_reference_comparison"] is False
    assert persisted["comparisons"]["trait_summary.tsv"]["unexpected_changes"] != []
