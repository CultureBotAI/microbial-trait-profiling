"""Which corpus a run read, recorded so a number can be traced to it (#222)."""

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from mtp.corpus_identity import (
    DEFAULT_KG_PROVENANCE_PATH,
    PINNED_KG_MICROBE_RELEASE,
    UNKNOWN_RELEASE,
    CorpusIdentity,
    corpus_identity,
    measure_matrix,
    read_kg_inputs,
    summary_lines,
    verify_provenance_source_hashes,
    write_identity,
)
from mtp.download import BUGSIGDB_EXPORTS_VERSION, PINNED_FULL_DUMP_MD5
from mtp.human_host_corpus import CORPUS_SEQUENCING_TYPES
from mtp.paths import file_sha256
from mtp.strain_rollup_record import StrainRollupRecordError, write_strain_rollup_record
from mtp.trait_composition import (
    LogRatioResult,
    LogRatioSummary,
    StudyLogRatio,
    write_log_ratio_outputs,
)


def _matrix(path: Path, *, mediadive: int, retained: int, rows: int) -> Path:
    header = ["study_id", "experiment_id", "inc_dec"]
    header += [f"mediadive.medium:{i}" for i in range(mediadive)]
    header += [f"METPO:{i}" for i in range(retained)]
    lines = ["\t".join(header)]
    counts = ["0"] * (mediadive + retained)
    for row in range(rows):
        lines.append("\t".join([f"bsdb:{row}", f"bsdb:{row}-1", "inc", *counts]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_the_identity_states_the_scope_the_table_records(tmp_path: Path) -> None:
    """A provenance record names the profile table's strain rollup scope, is null when
    the table records none, and refuses a table replaced without its record (#847, #939)."""
    table = _matrix(tmp_path / "m.tsv", mediadive=1, retained=2, rows=1)
    assert measure_matrix(table).strain_rollup is None
    write_strain_rollup_record(table, "phenotype")
    identity = measure_matrix(table)
    assert identity.strain_rollup == "phenotype"
    assert "  strain rollup    phenotype" in summary_lines(
        CorpusIdentity(identity, None, None, "provenance.json")
    )
    _matrix(table, mediadive=1, retained=2, rows=2)
    with pytest.raises(StrainRollupRecordError, match="not the table its record"):
        measure_matrix(table)


def test_the_three_column_groups_account_for_the_whole_header(tmp_path: Path) -> None:
    identity = measure_matrix(_matrix(tmp_path / "m.tsv", mediadive=4, retained=7, rows=3))

    assert identity.identifier_columns == 3
    assert identity.mediadive_columns == 4
    assert identity.retained_trait_columns == 7
    assert identity.header_columns == 14
    assert identity.profile_rows == 3


def test_two_matrices_differing_only_in_mediadive_columns_are_distinguishable(
    tmp_path: Path,
) -> None:
    """The 2,290-vs-2,333 case: same retained count, different matrix (#222)."""
    first = measure_matrix(_matrix(tmp_path / "a.tsv", mediadive=4, retained=7, rows=2))
    second = measure_matrix(_matrix(tmp_path / "b.tsv", mediadive=6, retained=7, rows=2))

    assert first.retained_trait_columns == second.retained_trait_columns
    # Identical retained counts, so only the full record tells them apart.
    assert first.mediadive_columns != second.mediadive_columns
    assert first.header_columns != second.header_columns
    assert first.sha256 != second.sha256


def test_a_header_that_does_not_add_up_raises(tmp_path: Path) -> None:
    identity = measure_matrix(_matrix(tmp_path / "m.tsv", mediadive=2, retained=2, rows=1))
    broken = type(identity)(**{**identity.__dict__, "retained_trait_columns": 99})

    with pytest.raises(ValueError, match="not the"):
        broken.check()


def test_measure_matrix_rejects_duplicate_or_missing_identifier_headers(tmp_path: Path) -> None:
    """measure_matrix must fail on real malformed input (#334)."""
    dup_path = tmp_path / "dup.tsv"
    dup_path.write_text(
        "study_id\tstudy_id\tfoo\tbar\tmediadive.medium:1\nbsdb:1\tbsdb:1\t0\t0\t0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="not the 5 columns in the header"):
        measure_matrix(dup_path)

    missing_id_path = tmp_path / "missing_id.tsv"
    missing_id_path.write_text(
        "study_id\texperiment_id\tMETPO:1\nbsdb:1\tbsdb:1-1\t0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="expected 3 identifier columns"):
        measure_matrix(missing_id_path)


def test_measure_matrix_ignores_blank_lines(tmp_path: Path) -> None:
    """Blank lines are skipped just like load_profiles (#338)."""
    path = tmp_path / "blank_lines.tsv"
    path.write_text(
        "study_id\texperiment_id\tinc_dec\tMETPO:1\n\nbsdb:1\tbsdb:1-1\tinc\t0\n\n",
        encoding="utf-8",
    )
    identity = measure_matrix(path)
    assert identity.profile_rows == 1


def test_a_missing_release_tag_is_recorded_rather_than_omitted(tmp_path: Path) -> None:
    """ "Not known" and "not applicable" differ, so the record says which (#222)."""
    sidecar = tmp_path / "provenance.json"
    sidecar.write_text(
        json.dumps({"files": {"merged-kg_nodes.tsv": {"md5": "abc", "size": 1}}}),
        encoding="utf-8",
    )

    release, inputs = read_kg_inputs(sidecar)

    assert release == UNKNOWN_RELEASE
    assert [item.name for item in inputs] == ["merged-kg_nodes.tsv"]
    assert inputs[0].declared_md5 == "abc"
    # The sidecar names a file that is not beside it, and the record says so.
    assert inputs[0].present is False


def test_a_declared_release_tag_is_used_when_present(tmp_path: Path) -> None:
    sidecar = tmp_path / "provenance.json"
    sidecar.write_text(json.dumps({"release": "kg-microbe-2026-04", "files": {}}), encoding="utf-8")

    release, _ = read_kg_inputs(sidecar)

    assert release == "kg-microbe-2026-04"


def test_pinned_20260120_release_is_resolved_from_tracked_or_legacy_sidecar(
    tmp_path: Path,
) -> None:
    """Tracked provenance.json and legacy sidecars matching 20260120 resolve to 20260120 (#336)."""
    tracked_release, _ = read_kg_inputs(DEFAULT_KG_PROVENANCE_PATH)
    assert tracked_release == PINNED_KG_MICROBE_RELEASE

    legacy_sidecar = tmp_path / "legacy_provenance.json"
    legacy_data = json.loads(DEFAULT_KG_PROVENANCE_PATH.read_text(encoding="utf-8"))
    legacy_data.pop("release", None)
    legacy_sidecar.write_text(json.dumps(legacy_data), encoding="utf-8")
    inferred_release, _ = read_kg_inputs(legacy_sidecar)
    assert inferred_release == PINNED_KG_MICROBE_RELEASE


def test_an_absent_or_malformed_sidecar_does_not_fail_the_run(tmp_path: Path) -> None:
    """A recipe must still record the matrix when the KG sidecar is unusable."""
    absent, inputs = read_kg_inputs(tmp_path / "nope.json")
    assert UNKNOWN_RELEASE in absent
    assert inputs == []

    broken = tmp_path / "broken.json"
    broken.write_text("not json", encoding="utf-8")
    release, inputs = read_kg_inputs(broken)
    assert UNKNOWN_RELEASE in release
    assert inputs == []


def test_the_written_record_names_both_sides_of_the_join(tmp_path: Path) -> None:
    matrix = _matrix(tmp_path / "m.tsv", mediadive=1, retained=2, rows=1)
    sidecar = tmp_path / "provenance.json"
    sidecar.write_text(
        json.dumps({"files": {"merged-kg_nodes.tsv": {"md5": "abc"}}}), encoding="utf-8"
    )

    identity = corpus_identity(matrix, bugsigdb_vintage="v", kg_provenance_path=sidecar)
    path = write_identity(identity, tmp_path / "out", code_commit="deadbeef")
    record = json.loads(path.read_text(encoding="utf-8"))

    assert record["matrix"]["sha256"] == identity.matrix.sha256
    assert record["bugsigdb_vintage"] == "v"
    assert record["kg_microbe"]["inputs"][0]["declared_md5"] == "abc"
    assert "sidecar" in record["kg_microbe"]["provenance_source"]
    assert record["code_commit"] == "deadbeef"
    # A vintage alone does not identify a matrix, and the record says so.
    assert "join" in record["note"]


def test_the_printed_summary_states_the_matrix_and_both_inputs(tmp_path: Path) -> None:
    matrix = _matrix(tmp_path / "m.tsv", mediadive=1, retained=2, rows=1)
    sidecar = tmp_path / "provenance.json"
    sidecar.write_text(
        json.dumps({"files": {"merged-kg_nodes.tsv": {"md5": "abc"}}}), encoding="utf-8"
    )

    text = "\n".join(
        summary_lines(corpus_identity(matrix, bugsigdb_vintage="v", kg_provenance_path=sidecar))
    )

    assert "sha256" in text
    assert "1 MediaDive" in text
    assert "2 retained" in text
    assert "merged-kg_nodes.tsv" in text
    assert "declared md5 (sidecar)" in text


def test_verify_provenance_source_hashes_detects_stale_and_fresh_records(
    tmp_path: Path,
) -> None:
    """Stale Justfile or source hashes in provenance records are flagged (#499)."""
    justfile = tmp_path / "Justfile"
    justfile.write_text("check:\n    pytest\n", encoding="utf-8")
    current_sha = file_sha256(justfile)

    fresh_record: dict[str, object] = {
        "source_hashes": {"Justfile": current_sha},
        "justfile_sha256": current_sha,
    }
    assert verify_provenance_source_hashes(fresh_record, repo_root=tmp_path) == {}

    stale_record: dict[str, object] = {
        "source_hashes": {"Justfile": "0" * 64, "missing.py": "f" * 64},
    }
    mismatches = verify_provenance_source_hashes(stale_record, repo_root=tmp_path)
    assert set(mismatches) == {"Justfile", "missing.py"}
    assert mismatches["Justfile"]["recorded"] == "0" * 64
    assert mismatches["Justfile"]["actual"] == current_sha
    assert mismatches["missing.py"]["actual"] == "(missing)"


def test_log_ratio_outputs_include_corpus_identity(tmp_path: Path) -> None:
    """write_log_ratio_outputs embeds corpus_identity in trait_log_ratio_provenance.json (#337)."""
    matrix = _matrix(tmp_path / "m.tsv", mediadive=1, retained=2, rows=2)
    experiments = tmp_path / "experiments.tsv"
    experiments.write_text("study_id\texperiment_id\nhost\n", encoding="utf-8")
    full_dump = tmp_path / "full_dump.csv"
    full_dump.write_text("# v1.3.1\n", encoding="utf-8")
    identity = corpus_identity(matrix, bugsigdb_vintage="v1.3.1")
    summary = LogRatioSummary(
        numerator="METPO:0",
        reference="METPO:1",
        studies_total=1,
        studies_eligible=1,
        studies_excluded=0,
        studies_partially_observed=0,
        studies_up=1,
        studies_down=0,
        studies_tied=0,
        studies_scored=1,
        paired_experiments=1,
        eligible_experiments=1,
        excluded_zero_experiments=0,
        mean_log2_ratio_change=0.5,
        median_log2_ratio_change=0.5,
        share_ratio_up=1.0,
        share_ci_low=0.025,
        share_ci_high=1.0,
        sign_p=1.0,
        by_q=1.0,
        status="tested",
    )
    study_row = StudyLogRatio(
        numerator="METPO:0",
        reference="METPO:1",
        study_id="bsdb:0",
        paired_experiments=1,
        eligible_experiments=1,
        excluded_zero_experiments=0,
        mean_log2_ratio_change=0.5,
    )
    result = LogRatioResult(summaries=[summary], studies=[study_row], family_size=1)
    write_log_ratio_outputs(
        result,
        tmp_path / "lr",
        profile_path=matrix,
        experiments_path=experiments,
        full_dump_path=full_dump,
        code_commit="abc1234",
        identity=identity,
    )
    prov = json.loads((tmp_path / "lr" / "trait_log_ratio_provenance.json").read_text("utf-8"))
    assert "corpus_identity" in prov
    assert prov["corpus_identity"]["matrix"]["sha256"] == identity.matrix.sha256
    assert prov["corpus_identity"]["bugsigdb_vintage"] == "v1.3.1"


def _load_script_module(script_name: str) -> Any:
    script_path = Path(__file__).resolve().parents[1] / "scripts" / script_name
    spec = importlib.util.spec_from_file_location(script_path.stem, script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_export_trait_enrichment_rankings_uses_n_profiles_tested_denominator() -> None:
    """build_overall_table divides by per-trait n_profiles_tested and validates fdr_scope (#119)."""
    mod = _load_script_module("export_trait_enrichment_rankings.py")
    parsed = mod.parse_args(["--all-results", "--fdr-scope", "global"])
    assert parsed.fdr_scope == "global"
    assert parsed.unfiltered is True

    df = pd.DataFrame(
        [
            {
                "study_id": "s1",
                "experiment_id": "e1",
                "inc_dec": "inc",
                "profile_id": "s1:e1:inc",
                "trait": "METPO:1",
                "p_value": 0.0001,
                "fdr_q": 0.001,
                "fold": 2.5,
                "obs": 5.0,
                "exp": 2.0,
                "testable_profiles": 50,
                "trait_total": 20,
                "is_significant": True,
            },
            {
                "study_id": "s2",
                "experiment_id": "e1",
                "inc_dec": "inc",
                "profile_id": "s2:e1:inc",
                "trait": "METPO:1",
                "p_value": 0.0005,
                "fdr_q": 0.005,
                "fold": 2.0,
                "obs": 4.0,
                "exp": 2.0,
                "testable_profiles": 50,
                "trait_total": 20,
                "is_significant": True,
            },
            {
                "study_id": "s3",
                "experiment_id": "e1",
                "inc_dec": "inc",
                "profile_id": "s3:e1:inc",
                "trait": "METPO:1",
                "p_value": 0.8,
                "fdr_q": 0.9,
                "fold": 0.5,
                "obs": 1.0,
                "exp": 2.0,
                "testable_profiles": 50,
                "trait_total": 20,
                "is_significant": False,
            },
            {
                "study_id": "s4",
                "experiment_id": "e1",
                "inc_dec": "dec",
                "profile_id": "s4:e1:dec",
                "trait": "METPO:1",
                "p_value": 0.9,
                "fdr_q": 0.95,
                "fold": 0.5,
                "obs": 1.0,
                "exp": 2.0,
                "testable_profiles": 50,
                "trait_total": 20,
                "is_significant": False,
            },
        ]
    )
    out = mod.build_overall_table(
        df,
        "fdr_q",
        total_profiles=50,
        fdr_scope="global",
        complete_results=True,
    )
    row = out.iloc[0]
    assert int(row["testable_profiles"]) == 50
    assert int(row["n_profiles_tested"]) == 4
    assert int(row["hit_rows"]) == 2
    assert float(row["prevalence"]) == pytest.approx(0.5)
    assert row["prevalence_summary"] == "2 / 4"

    # Significant-only inputs carry a flag too; it is not proof of completeness.
    selected = mod.build_overall_table(df[df["is_significant"]], "fdr_q")
    assert selected.iloc[0]["prevalence_summary"] == "2 / 50"

    with pytest.raises(ValueError, match="Unsupported fdr_scope"):
        mod.build_overall_table(
            df,
            "fdr_q",
            total_profiles=50,
            fdr_scope="invalid",
        )


def test_rerun_recipe_hygiene_invariants() -> None:
    """Justfile and arm rerun scripts preserve pinned inputs and scope cleanliness checks."""
    repo_root = Path(__file__).resolve().parents[1]
    justfile_text = (repo_root / "Justfile").read_text(encoding="utf-8")
    assert 'code_paths := "src scripts methods schema Justfile pyproject.toml uv.lock"' in (
        justfile_text
    )
    assert "git status --porcelain" not in justfile_text
    assert "bsdbng.integrated_evidence" not in justfile_text
    assert "python -m mtp.integrated_evidence" in justfile_text
    assert (
        'calibrate replicates="10" seed="42" scenario="all": study-trait-profiles-ora'
        in justfile_text
    )

    noel_script = (repo_root / "scripts" / "rerun_directional_consistency_arm_noel.sh").read_text(
        encoding="utf-8"
    )
    # The standalone canonical recipe keeps its human-only contract (#656), on 16S and
    # shotgun experiments only (#742); reproducing the committed table uses the same
    # flags (#847).
    assert (
        'profile_flags := \'--host-species "Homo sapiens" '
        f"--sequencing-types {' '.join(CORPUS_SEQUENCING_TYPES)} --drop-empty-traits'\n"
        in justfile_text
    )
    assert (
        "study-trait-profiles: taxa-trait-matrix\n"
        "\tuv run python scripts/create_study_trait_profiles.py {{profile_flags}}\n"
        in justfile_text
    )
    assert justfile_text.count("{{profile_flags}}") == 2
    assert (
        "study-trait-profiles-ora:\n"
        "\tuv run python scripts/create_study_trait_profiles.py --binary --include-taxa-count "
        f'--host-species "Homo sapiens" --sequencing-types {" ".join(CORPUS_SEQUENCING_TYPES)}\n'
        in justfile_text
    )
    jed_script = (repo_root / "scripts" / "rerun_composition_arm_jed.sh").read_text(
        encoding="utf-8"
    )
    assert 'exec just rerun-noel "$@"' in noel_script
    assert 'exec just rerun-jed "$@"' in jed_script
    assert "\njust table1-taxon-influence\n" not in noel_script
    # Figures are drawn headless, whatever the session's display (#719)
    assert '\nexport MPLBACKEND := "Agg"\n' in justfile_text
    ilyas_script = (repo_root / "scripts" / "rerun_enrichment_ora_arm_ilyas.sh").read_text(
        encoding="utf-8"
    )
    for script in (noel_script, jed_script, ilyas_script):
        assert "\nexport MPLBACKEND=Agg\n" in script

    prov_path = repo_root / "data" / "raw" / "download_provenance.json"
    assert prov_path.exists()
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    assert prov["bugsigdb_release"]["version"] == BUGSIGDB_EXPORTS_VERSION
    assert prov["full_dump.csv"]["md5"] == PINNED_FULL_DUMP_MD5
