"""Qualify archived parity using bound synthetic bytes and real Git parents (#1168)."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
from collections.abc import Iterator
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any

import pytest
from git import Repo

from mtp.workflow import output_roots

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/adhoc/issue1117_compare_archived_output.py"
HELPER_PATH = "scripts/adhoc/issue1117_compare_archived_output.py"
SPEC = importlib.util.spec_from_file_location("issue1117_parity_under_test", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Cannot load the archived parity helper")
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)

SCRIPTS = (
    "src/mtp/matrix_overlap_coverage.py",
    "src/mtp/trait_lists.py",
    "src/mtp/matrix_overlap_inputs.py",
    "src/mtp/human_host_corpus.py",
    "src/mtp/workflow/coverage.py",
    "src/mtp/workflow/output_roots.py",
)
METHOD = "methods/matrix_coverage_and_pd_overlap.md"
OLD_EXPORT_LABEL = "All-hosts BugSigDB export table (bugsigdb_experiments.tsv)"
CURRENT_EXPORT_LABEL = "All-hosts BugSigDB export table (experiments.tsv)"
LABEL_ADAPTER = {
    "file": "filter_cascade.tsv",
    "row": {"dimension": "studies", "step_id": "all_hosts_bugsigdb_export"},
    "column": "step_label",
    "old": OLD_EXPORT_LABEL,
    "new": CURRENT_EXPORT_LABEL,
    "cause_commit": "8832c67ee6022a1eb181cbf95c87f516bc9f20e3",
    "integration_commit": "c1e374e99517caf473df01cec4dcb2db5d29ff4a",
}
CURRENT_TRAITS = tuple(f"METPO:2000{number:03d}" for number in range(9))
REMOVED_TRAIT = "METPO:1000605"
OLD_TRAITS = (*CURRENT_TRAITS, REMOVED_TRAIT)
POPULATION = "disease_vs_healthy_tier1_2_subject_controls"
TABLES = (
    "matrix_trait_overlap.tsv",
    "tier12_coverage_by_rank.tsv",
    "filter_cascade.tsv",
    "parkinsons_stool_study_overlap_summary.tsv",
    "parkinsons_stool_study_pairs.tsv",
)


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def write_json(path: Path, record: object) -> None:
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def read_json(path: Path) -> Any:
    return json.loads(path.read_bytes())


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    path.write_text(stream.getvalue(), encoding="utf-8")


def rows(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(path.read_text()), delimiter="\t"))


def declared(path: Path, trait_ids: tuple[str, ...]) -> dict[str, Any]:
    checksum = sha256(path.read_bytes())
    return {
        "trait_list_file": path.name,
        "trait_list_sha256": checksum,
        "trait_count": len(trait_ids),
        "trait_ids": list(trait_ids),
        "source_name": path.name,
        "sha256": checksum,
    }


@dataclass
class ArchiveCase:
    args: argparse.Namespace
    repo: Repo
    reference: Path
    current: Path
    core: Path

    def reseal(self, side: str = "current") -> None:
        root = self.current if side == "current" else self.reference
        manifest = root / "core_consumption.json"
        record = read_json(manifest)
        for entry in record["outputs"]:
            raw = (root / Path(entry["path"]).name).read_bytes()
            entry.update(sha256=sha256(raw), size_bytes=len(raw))
        write_json(manifest, record)
        setattr(self.args, f"{side}_sha256", sha256(manifest.read_bytes()))

    def change_json(self, name: str, section: str, key: str, value: object) -> None:
        path = self.current / name
        record = read_json(path)
        record[section][key] = value
        write_json(path, record)
        self.reseal()

    def compare(self) -> dict[str, Any]:
        result: dict[str, Any] = HELPER.compare(self.args)
        return result


@pytest.fixture
def archive_case(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[ArchiveCase]:
    """Only repository location and OS temp roots are controlled; proofs remain real."""
    producer = tmp_path / "producer"
    producer.mkdir()
    repo = Repo.init(producer)
    for relative in (*SCRIPTS, METHOD):
        path = producer / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"Synthetic public producer bytes for {relative}\n", encoding="utf-8")
    helper_path = producer / HELPER_PATH
    helper_path.parent.mkdir(parents=True)
    helper_path.write_bytes(SCRIPT.read_bytes())
    repo.index.add([*SCRIPTS, METHOD, HELPER_PATH])
    base = repo.index.commit("Base producer")
    (producer / SCRIPTS[0]).write_text("Synthetic issue implementation\n", encoding="utf-8")
    repo.index.add([SCRIPTS[0]])
    issue = repo.index.commit("Issue implementation")
    compatibility = repo.index.commit("Compatibility branch", parent_commits=[base], head=False)
    trial = repo.index.commit("Trial merge", parent_commits=[issue, compatibility])
    assert [parent.hexsha for parent in trial.parents] == [issue.hexsha, compatibility.hexsha]
    assert not repo.is_dirty(untracked_files=True)

    monkeypatch.setattr(HELPER, "Repo", lambda _path: repo)
    monkeypatch.setattr(HELPER, "__file__", str(helper_path))
    monkeypatch.setattr(
        output_roots, "system_temporary_directories", lambda: [tmp_path / "volatile"]
    )
    core = tmp_path / "core" / "run.json"
    core.parent.mkdir()
    write_json(
        core,
        {
            "run_id": "core-fixture",
            "git_commit": base.hexsha,
            "git_dirty": False,
            "status": "completed",
            "definition": {"id": "bugsigdb-core-subset"},
            "contrast_population": POPULATION,
            "direction_basis": "disease_vs_comparison",
            "parameters": [
                {"name": "strain_rollup", "value": "phenotype"},
                {"name": "contrast_population", "value": POPULATION},
            ],
        },
    )
    current_list, old_list = tmp_path / "nine.txt", tmp_path / "ten.txt"
    current_list.write_text("\n".join(CURRENT_TRAITS) + "\n", encoding="utf-8")
    old_list.write_text("\n".join(OLD_TRAITS) + "\n", encoding="utf-8")
    reference, current = tmp_path / "reference", tmp_path / "current"
    common = {
        "consumer": "matrix-overlap-coverage",
        "core_run_id": "core-fixture",
        "workflow_id": "bugsigdb-core-subset",
        "strain_rollup": "phenotype",
        "contrast_population": POPULATION,
    }
    overlap = [
        {"trait_a_id": first, "trait_b_id": second, "intersection_taxa": index + 1}
        for index, (first, second) in enumerate(combinations(OLD_TRAITS, 2))
    ]
    pd_summary: dict[str, Any] = {"scope": "all_pd", "study_count": 3, "shared_taxa": 4}
    coverage = {"source_occurrence_rows": 17, "retained_occurrence_rows": 11}
    for root, trait_path, traits in (
        (reference, old_list, OLD_TRAITS),
        (current, current_list, CURRENT_TRAITS),
    ):
        root.mkdir()
        chosen = [
            row for row in overlap if row["trait_a_id"] in traits and row["trait_b_id"] in traits
        ]
        write_rows(root / TABLES[0], chosen)
        write_rows(root / TABLES[1], [{"rank": "genus", "taxa": 11, "matched_taxa": 7}])
        write_rows(
            root / TABLES[2],
            [
                {
                    "dimension": "studies",
                    "step_id": "all_hosts_bugsigdb_export",
                    "step_label": CURRENT_EXPORT_LABEL,
                    "studies": 3,
                    "experiments": 4,
                }
            ],
        )
        pd_row = dict(pd_summary)
        cover = dict(coverage)
        if root == current:
            pd_row.update(
                input_occurrence_rows=17,
                excluded_outside_active_profile_rows=6,
                excluded_outside_active_profile_keys=2,
                input_experiments=4,
                excluded_experiments_without_active_profiles=0,
                source_host_case_variant_experiments=0,
                source_host_case_variant_profile_keys=0,
            )
            cover.update(
                tier12_input_occurrence_rows=17,
                tier12_excluded_occurrence_rows=6,
                tier12_excluded_profile_keys=2,
                tier12_input_experiments=4,
                tier12_excluded_experiments_without_active_profiles=0,
                source_host_case_variant_experiments=0,
                source_host_case_variant_profile_keys=0,
            )
        write_rows(root / TABLES[3], [pd_row])
        write_rows(
            root / TABLES[4],
            [{"study_a_id": "study-a", "study_b_id": "study-b", "jaccard_fraction": "0.25"}],
        )
        declared_summary: dict[str, Any] = {
            **declared(trait_path, traits),
            "pair_count": len(chosen),
        }
        write_json(
            root / "summary.json",
            {
                **common,
                "declared_traits": declared_summary,
                "coverage_and_cascade": cover,
                "parkinsons_stool_overlap": [pd_row],
            },
        )
        outputs = [*TABLES, "summary.json"]
        if root == current:
            write_json(
                root / "provenance.json",
                {
                    **common,
                    "code_commit": trial.hexsha,
                    "core_git_commit_sha": base.hexsha,
                    "scripts": {
                        name: {"sha256": sha256((producer / name).read_bytes())} for name in SCRIPTS
                    },
                    "methods": {METHOD: {"sha256": sha256((producer / METHOD).read_bytes())}},
                    "script_sha256": sha256((producer / SCRIPTS[0]).read_bytes()),
                    "methods_sha256": sha256((producer / METHOD).read_bytes()),
                    "declared_traits": declared(trait_path, traits),
                },
            )
            outputs.append("provenance.json")
        write_json(
            root / "core_consumption.json",
            {
                "consumer": "matrix-overlap-coverage",
                "run_id": "core-fixture",
                "run_manifest_sha256": sha256(core.read_bytes()),
                "outputs": [
                    {"path": f"/archived/{root.name}/{name}", "sha256": "", "size_bytes": 0}
                    for name in outputs
                ],
            },
        )
    args = argparse.Namespace(
        reference_consumer=reference / "core_consumption.json",
        current_consumer=current / "core_consumption.json",
        trait_list=current_list,
        output_dir=tmp_path / "uncreated" / "parity",
        trial_commit=trial.hexsha,
        issue_commit=issue.hexsha,
        compatibility_commit=compatibility.hexsha,
        helper_commit=trial.hexsha,
        core_manifest=core,
        core_sha256=sha256(core.read_bytes()),
    )
    case = ArchiveCase(args, repo, reference, current, core)
    case.reseal("reference")
    case.reseal()
    try:
        yield case
    finally:
        repo.close()


def test_only_nine_of_ten_trait_projection_and_declared_counters_are_allowed(
    archive_case: ArchiveCase,
) -> None:
    result = archive_case.compare()
    report_path = archive_case.args.output_dir / "parity.json"
    report = read_json(report_path)
    assert result["sha256"] == sha256(report_path.read_bytes())
    assert report["removed_traits"] == [REMOVED_TRAIT]
    assert report["producer"]["parents"] == [
        archive_case.args.issue_commit,
        archive_case.args.compatibility_commit,
    ]
    assert report["historical_clean_producer_verified"] is False
    assert report["unexpected_changes"] == []
    assert report["new_exclusion_counters"] == {
        "tier12_input_occurrence_rows": 17,
        "tier12_excluded_occurrence_rows": 6,
        "tier12_excluded_profile_keys": 2,
        "tier12_input_experiments": 4,
        "tier12_excluded_experiments_without_active_profiles": 0,
        "source_host_case_variant_experiments": 0,
        "source_host_case_variant_profile_keys": 0,
    }
    matrix = next(row for row in result["tables"] if row["table"] == TABLES[0])
    assert (matrix["reference_rows"], matrix["current_rows"]) == (45, 36)
    assert matrix["exact_projected_row_equality"] is True
    assert matrix["byte_identical"] is False


@pytest.mark.parametrize("name", TABLES)
def test_bound_but_undeclared_table_change_is_refused(archive_case: ArchiveCase, name: str) -> None:
    path = archive_case.current / name
    changed = rows(path)
    column = next(
        key
        for key in changed[0]
        if key
        not in {
            "trait_a_id",
            "trait_b_id",
            "rank",
            "stage",
            "scope",
            "study_a_id",
            "study_b_id",
            "dimension",
            "step_id",
            "step_label",
        }
    )
    changed[0][column] = "999"
    write_rows(path, changed)
    archive_case.reseal()
    with pytest.raises(ValueError, match="Unexpected numerical"):
        archive_case.compare()
    assert not archive_case.args.output_dir.parent.exists()


def test_changed_bound_table_bytes_are_refused(archive_case: ArchiveCase) -> None:
    path = archive_case.current / TABLES[1]
    path.write_bytes(path.read_bytes().replace(b"\t11\t", b"\t12\t"))
    with pytest.raises(ValueError, match="Output differs from bound bytes"):
        archive_case.compare()
    assert not archive_case.args.output_dir.parent.exists()


@pytest.mark.parametrize(
    "fault", ["reversed_parents", "not_head", "dirty", "wrong_producer", "wrong_script_hash"]
)
def test_real_git_producer_guards_fail_closed(archive_case: ArchiveCase, fault: str) -> None:
    args = archive_case.args
    if fault == "reversed_parents":
        args.issue_commit, args.compatibility_commit = args.compatibility_commit, args.issue_commit
    elif fault == "not_head":
        archive_case.repo.git.checkout(args.issue_commit)
    elif fault == "dirty":
        assert archive_case.repo.working_tree_dir is not None
        (Path(archive_case.repo.working_tree_dir) / "untracked.txt").write_text("uncommitted\n")
    else:
        path = archive_case.current / "provenance.json"
        record = read_json(path)
        if fault == "wrong_producer":
            record["code_commit"] = args.issue_commit
        else:
            record["scripts"][SCRIPTS[0]]["sha256"] = "0" * 64
        write_json(path, record)
        archive_case.reseal()
    with pytest.raises(ValueError, match=r"parents|clean|Producer code"):
        archive_case.compare()
    assert not args.output_dir.parent.exists()


@pytest.mark.parametrize(
    "fault",
    [
        "empty_scripts",
        "missing_script",
        "extra_script",
        "script_alias",
        "method_alias",
        "method_hash",
        "missing_methods",
    ],
)
def test_complete_producer_inventory_is_required(archive_case: ArchiveCase, fault: str) -> None:
    path = archive_case.current / "provenance.json"
    record = read_json(path)
    if fault == "empty_scripts":
        record["scripts"] = {}
    elif fault == "missing_script":
        del record["scripts"][SCRIPTS[-1]]
    elif fault == "extra_script":
        record["scripts"][METHOD] = record["methods"][METHOD]
    elif fault == "script_alias":
        record["script_sha256"] = "0" * 64
    elif fault == "method_alias":
        record["methods_sha256"] = "0" * 64
    elif fault == "method_hash":
        record["methods"][METHOD]["sha256"] = "0" * 64
    else:
        record["methods"] = {}
    write_json(path, record)
    archive_case.reseal()
    with pytest.raises(ValueError, match=r"(?i)script|method|inventory"):
        archive_case.compare()
    assert not archive_case.args.output_dir.parent.exists()


@pytest.mark.parametrize(
    "field", ["trait_count", "pair_count", "source_name", "trait_list_file", "sha256"]
)
def test_declared_metadata_cannot_be_hidden_by_summary_projection(
    archive_case: ArchiveCase, field: str
) -> None:
    value: object = 999 if field.endswith("count") else "incorrect"
    archive_case.change_json("summary.json", "declared_traits", field, value)
    with pytest.raises(ValueError, match=r"(?i)trait|declared|pair"):
        archive_case.compare()
    assert not archive_case.args.output_dir.parent.exists()


@pytest.mark.parametrize(
    ("field", "value"),
    [
        (None, None),
        ("trait_ids", [REMOVED_TRAIT]),
        ("trait_count", 999),
        ("trait_list_sha256", "0" * 64),
    ],
)
def test_provenance_declared_traits_are_checked_independently(
    archive_case: ArchiveCase, field: str | None, value: object
) -> None:
    summary_before = (archive_case.current / "summary.json").read_bytes()
    path = archive_case.current / "provenance.json"
    provenance = read_json(path)
    if field is None:
        provenance.pop("declared_traits")
    else:
        provenance["declared_traits"][field] = value
    write_json(path, provenance)
    archive_case.reseal()
    assert (archive_case.current / "summary.json").read_bytes() == summary_before
    with pytest.raises(ValueError, match="Current provenance declared trait metadata differs"):
        archive_case.compare()
    assert not archive_case.args.output_dir.parent.exists()


def test_new_exclusion_diagnostics_have_explicit_parity_allowances(
    archive_case: ArchiveCase,
) -> None:
    coverage = {
        "tier12_input_experiments": 7,
        "tier12_excluded_experiments_without_active_profiles": 3,
        "source_host_case_variant_experiments": 2,
        "source_host_case_variant_profile_keys": 4,
    }
    pd_diagnostics = {
        "input_experiments": 7,
        "excluded_experiments_without_active_profiles": 3,
        "source_host_case_variant_experiments": 2,
        "source_host_case_variant_profile_keys": 4,
    }
    summary_path = archive_case.current / "summary.json"
    summary = read_json(summary_path)
    summary["coverage_and_cascade"].update(coverage)
    summary["parkinsons_stool_overlap"][0].update(pd_diagnostics)
    write_json(summary_path, summary)
    table_path = archive_case.current / "parkinsons_stool_study_overlap_summary.tsv"
    pd_rows = rows(table_path)
    pd_rows[0].update({key: str(value) for key, value in pd_diagnostics.items()})
    write_rows(table_path, pd_rows)
    archive_case.reseal()
    archive_case.compare()
    report = read_json(archive_case.args.output_dir / "parity.json")
    assert {key: report["new_exclusion_counters"][key] for key in coverage} == coverage
    assert report["unexpected_changes"] == []


def test_unreviewed_summary_counter_is_not_discarded(archive_case: ArchiveCase) -> None:
    archive_case.change_json("summary.json", "coverage_and_cascade", "unreviewed_counter", 42)
    with pytest.raises(ValueError, match="Unexpected summary difference"):
        archive_case.compare()
    assert not archive_case.args.output_dir.parent.exists()


@pytest.mark.parametrize("side", ["reference", "current"])
@pytest.mark.parametrize("field", ["run_id", "run_manifest_sha256"])
def test_consumers_must_bind_the_same_pinned_core(
    archive_case: ArchiveCase, side: str, field: str
) -> None:
    root = archive_case.reference if side == "reference" else archive_case.current
    path = root / "core_consumption.json"
    record = read_json(path)
    record[field] = "wrong-run" if field == "run_id" else "0" * 64
    write_json(path, record)
    archive_case.reseal(side)
    with pytest.raises(ValueError, match=r"(?i)core|consumer"):
        archive_case.compare()
    assert not archive_case.args.output_dir.parent.exists()


@pytest.mark.parametrize(
    "field", ["core_run_id", "core_git_commit_sha", "strain_rollup", "contrast_population"]
)
def test_current_provenance_must_match_the_pinned_core(
    archive_case: ArchiveCase, field: str
) -> None:
    path = archive_case.current / "provenance.json"
    record = read_json(path)
    record[field] = "wrong-core-value"
    write_json(path, record)
    archive_case.reseal()
    with pytest.raises(ValueError, match=r"(?i)core|provenance|scope|population"):
        archive_case.compare()
    assert not archive_case.args.output_dir.parent.exists()


def test_core_pin_is_checked_against_actual_bytes(archive_case: ArchiveCase) -> None:
    archive_case.core.write_bytes(archive_case.core.read_bytes() + b"\n")
    with pytest.raises(ValueError, match=r"(?i)core|pin"):
        archive_case.compare()
    assert not archive_case.args.output_dir.parent.exists()


def historical_cascade(case: ArchiveCase) -> None:
    path = case.reference / "filter_cascade.tsv"
    entries = rows(path)
    entries[0]["step_label"] = OLD_EXPORT_LABEL
    write_rows(path, entries)
    case.reseal("reference")


def test_exact_historical_filename_label_adapter_is_declared(archive_case: ArchiveCase) -> None:
    historical_cascade(archive_case)
    archive_case.compare()
    report = read_json(archive_case.args.output_dir / "parity.json")
    assert report["reference_adapters"] == [LABEL_ADAPTER]
    assert report["unexpected_changes"] == []
    assert rows(archive_case.reference / "filter_cascade.tsv")[0]["step_label"] == OLD_EXPORT_LABEL
    assert (
        rows(archive_case.current / "filter_cascade.tsv")[0]["step_label"] == CURRENT_EXPORT_LABEL
    )


@pytest.mark.parametrize("fault", ["new_label", "old_label", "dimension", "step_id", "numerical"])
def test_reference_adapter_cannot_hide_other_changes(archive_case: ArchiveCase, fault: str) -> None:
    historical_cascade(archive_case)
    old_path = archive_case.reference / "filter_cascade.tsv"
    new_path = archive_case.current / "filter_cascade.tsv"
    before, after = rows(old_path), rows(new_path)
    if fault == "new_label":
        after[0]["step_label"] = "A different current label"
    elif fault == "old_label":
        before[0]["step_label"] = "A different historical label"
    elif fault in {"dimension", "step_id"}:
        before[0][fault] = after[0][fault] = "unrelated-row"
    else:
        after[0]["studies"] = "999"
    write_rows(old_path, before)
    write_rows(new_path, after)
    archive_case.reseal("reference")
    archive_case.reseal()
    with pytest.raises(ValueError, match="Unexpected numerical"):
        archive_case.compare()
    assert not archive_case.args.output_dir.parent.exists()


def test_helper_commit_is_explicit_and_distinct_from_output_producer(
    archive_case: ArchiveCase,
) -> None:
    archive_case.compare()
    report = read_json(archive_case.args.output_dir / "parity.json")
    assert HELPER.__file__ is not None
    assert report["helper"] == {
        "commit": archive_case.args.helper_commit,
        "tree": archive_case.repo.head.commit.tree.hexsha,
        "sha256": sha256(Path(HELPER.__file__).read_bytes()),
        "size_bytes": len(Path(HELPER.__file__).read_bytes()),
    }
    assert report["producer"]["trial"] == archive_case.args.trial_commit
    assert report["reference_adapters"] == []


def test_later_clean_helper_only_commit_can_check_unchanged_trial_outputs(
    archive_case: ArchiveCase,
) -> None:
    assert HELPER.__file__ is not None
    helper = Path(HELPER.__file__)
    helper.write_bytes(helper.read_bytes() + b"\n# Helper-only qualification revision.\n")
    archive_case.repo.index.add([HELPER_PATH])
    helper_commit = archive_case.repo.index.commit("Qualification helper only")
    archive_case.args.helper_commit = helper_commit.hexsha
    assert helper_commit.parents[0].hexsha == archive_case.args.trial_commit
    assert not archive_case.repo.is_dirty(untracked_files=True)
    archive_case.compare()
    report = read_json(archive_case.args.output_dir / "parity.json")
    assert report["helper"]["commit"] == helper_commit.hexsha
    assert report["helper"]["sha256"] == sha256(helper.read_bytes())
    assert report["producer"]["trial"] == archive_case.args.trial_commit
    assert report["helper"]["commit"] != report["producer"]["trial"]
    assert report["unexpected_changes"] == []


@pytest.mark.parametrize("fault", ["wrong_pin", "dirty_helper", "uncommitted_helper", "bytes"])
def test_helper_commit_and_actual_execution_bytes_fail_closed(
    archive_case: ArchiveCase, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, fault: str
) -> None:
    assert HELPER.__file__ is not None
    helper = Path(HELPER.__file__)
    if fault == "wrong_pin":
        archive_case.args.helper_commit = archive_case.args.issue_commit
    elif fault == "dirty_helper":
        helper.write_bytes(helper.read_bytes() + b"\n# Uncommitted revision.\n")
    elif fault == "uncommitted_helper":
        archive_case.repo.index.remove([HELPER_PATH])
        commit = archive_case.repo.index.commit("Remove tracked qualification helper")
        archive_case.args.helper_commit = commit.hexsha
        assert helper.exists()
    else:
        executable = tmp_path / "different-execution" / HELPER_PATH
        executable.parent.mkdir(parents=True)
        executable.write_bytes(helper.read_bytes() + b"\n# Different loaded helper bytes.\n")
        monkeypatch.setattr(HELPER, "__file__", str(executable))
        assert not archive_case.repo.is_dirty(untracked_files=True)
    with pytest.raises(ValueError, match=r"(?i)helper|clean|blob|commit"):
        archive_case.compare()
    assert not archive_case.args.output_dir.parent.exists()
