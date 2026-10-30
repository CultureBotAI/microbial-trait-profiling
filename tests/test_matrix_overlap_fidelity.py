"""Independent input-fidelity and population regressions for #1117 (2026-10-28)."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd
import pytest

from mtp import matrix_overlap_coverage as analysis
from mtp.workflow import consumers, runtime
from test_matrix_overlap_coverage import FAC_AER, FAC_ANA, _synthetic_matrix_and_columns

POPULATION = "disease_vs_healthy_tier1_2_subject_controls"


@pytest.fixture
def corpus() -> dict[str, pd.DataFrame]:
    matrix, columns = _synthetic_matrix_and_columns()
    experiments = pd.DataFrame(
        [
            {
                "study_id": f"s{i}",
                "experiment_id": f"e{i}",
                "body_site": "Feces",
                "condition": "Parkinson's disease",
                "host_species": "Homo sapiens",
                "sequencing_type": "16S",
            }
            for i in range(1, 4)
        ]
    )
    taxa = pd.DataFrame(
        [
            {
                "study_id": study,
                "experiment_id": experiment,
                "direction": direction,
                "taxon_id": taxon,
                "taxonomic_rank": "genus",
            }
            for study, experiment, direction, taxon in [
                ("s1", "e1", "increased", "NCBITaxon:10"),
                ("s1", "e1", "increased", "NCBITaxon:10"),
                ("s1", "e1", "decreased", "NCBITaxon:11"),
                ("s2", "e2", "increased", "NCBITaxon:10"),
                ("s2", "e2", "decreased", "NCBITaxon:12"),
                ("s3", "e3", "increased", "NCBITaxon:13"),
            ]
        ]
    )
    keys = pd.DataFrame(
        [
            {"study_id": f"s{i}", "experiment_id": f"e{i}", "inc_dec": direction}
            for i in (1, 2)
            for direction in ("inc", "dec")
        ]
    )
    return {
        "matrix": matrix,
        "columns": columns,
        "tier12_experiments": experiments.copy(),
        "tier12_taxa": taxa.copy(),
        "tier12_profile_keys": keys,
        "source_experiments": experiments,
        "source_taxa": taxa,
    }


def coverage(corpus: dict[str, pd.DataFrame]) -> tuple[Any, Any, Any]:
    return analysis.compute_tier12_coverage_and_cascade(
        matrix=corpus["matrix"],
        columns=corpus["columns"],
        tier12_experiments=corpus["tier12_experiments"],
        tier12_taxa=corpus["tier12_taxa"],
        tier12_profile_keys=corpus["tier12_profile_keys"],
        source_experiments=corpus["source_experiments"],
        source_taxa=corpus["source_taxa"],
    )


def overlap(corpus: dict[str, pd.DataFrame]) -> tuple[Any, Any]:
    return analysis.compute_parkinsons_stool_overlap(
        **{name: frame for name, frame in corpus.items() if name not in {"matrix", "columns"}}
    )


@pytest.mark.parametrize("bad", [None, "", " ", " s1", 1, float("nan")])
@pytest.mark.parametrize(
    ("table", "key"),
    [
        ("matrix", "taxon_id"),
        ("columns", "column_id"),
        ("tier12_experiments", "study_id"),
        ("tier12_taxa", "taxon_id"),
        ("tier12_profile_keys", "experiment_id"),
        ("source_taxa", "study_id"),
    ],
)
def test_invalid_identifiers_refused(
    corpus: dict[str, pd.DataFrame], table: str, key: str, bad: str | int | float | None
) -> None:
    corpus[table] = corpus[table].astype({key: object})
    corpus[table].loc[0, key] = bad
    with pytest.raises(ValueError, match=r"identity|identifier|ID"):
        coverage(corpus)


@pytest.mark.parametrize(
    "table", ["matrix", "columns", "tier12_experiments", "tier12_profile_keys"]
)
def test_duplicate_unique_identities_refused(corpus: dict[str, pd.DataFrame], table: str) -> None:
    corpus[table] = pd.concat([corpus[table], corpus[table].iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match=r"Duplicate|duplicate"):
        coverage(corpus)


@pytest.mark.parametrize("bad", ["sideways", "", None])
@pytest.mark.parametrize(
    ("table", "key"),
    [
        ("source_taxa", "direction"),
        ("tier12_taxa", "direction"),
        ("tier12_profile_keys", "inc_dec"),
    ],
)
def test_unknown_direction_refused(
    corpus: dict[str, pd.DataFrame], table: str, key: str, bad: str | None
) -> None:
    corpus[table].loc[0, key] = bad
    with pytest.raises(ValueError, match="direction"):
        coverage(corpus)


@pytest.mark.parametrize("bad", ["nonsense", None, float("inf"), float("-inf"), -1, True])
def test_invalid_matrix_cells_refused_even_outside_focal_traits(
    corpus: dict[str, pd.DataFrame], bad: str | int | float | None
) -> None:
    corpus["matrix"] = corpus["matrix"].astype({FAC_ANA: object})
    corpus["matrix"].loc[0, FAC_ANA] = bad
    with pytest.raises(ValueError, match=r"numeric|finite|negative"):
        coverage(corpus)
    with pytest.raises(ValueError, match=r"numeric|finite|negative"):
        analysis.compute_matrix_trait_overlap(
            corpus["matrix"], corpus["columns"], [FAC_AER, FAC_ANA]
        )


@pytest.mark.parametrize("table", ["source_taxa", "tier12_taxa", "tier12_profile_keys"])
def test_orphan_experiment_keys_refused(corpus: dict[str, pd.DataFrame], table: str) -> None:
    corpus[table].loc[0, "study_id"] = "wrong-study"
    with pytest.raises(ValueError, match=r"metadata|experiment|orphan"):
        coverage(corpus)
    with pytest.raises(ValueError, match=r"metadata|experiment|orphan"):
        overlap(corpus)


@pytest.mark.parametrize("consumer", [coverage, overlap])
def test_active_profile_needs_tier12_metadata_even_when_source_metadata_exists(
    corpus: dict[str, pd.DataFrame], consumer: Any
) -> None:
    corpus["tier12_experiments"] = corpus["tier12_experiments"].iloc[1:].copy()
    with pytest.raises(ValueError, match="Active profile contains an orphan experiment"):
        consumer(corpus)


def test_inactive_experiments_are_counted_even_without_reported_taxa(
    corpus: dict[str, pd.DataFrame],
) -> None:
    extra = corpus["tier12_experiments"].iloc[[0]].copy()
    extra.loc[:, "study_id"] = "s4"
    extra.loc[:, "experiment_id"] = "e4"
    extra.loc[:, "body_site"] = "Blood"
    extra.loc[:, "condition"] = "Other condition"
    for name in ("tier12_experiments", "source_experiments"):
        corpus[name] = pd.concat([corpus[name], extra], ignore_index=True)
    _, _, summary = coverage(corpus)
    assert summary["tier12_input_experiments"] == 4
    assert summary["tier12_excluded_experiments_without_active_profiles"] == 2
    assert summary["tier12_experiments"] == 2
    assert summary["tier12_excluded_occurrence_rows"] == 1
    summaries, _ = overlap(corpus)
    tier12 = next(row for row in summaries if row["scope"] == "tier12_paired")
    assert tier12["input_experiments"] == 4
    assert tier12["excluded_experiments_without_active_profiles"] == 2
    assert tier12["n_experiments"] == 2
    assert next(row for row in summaries if row["scope"] == "full_export")["n_studies"] == 3


@pytest.mark.parametrize("positive", [False, True])
def test_matrix_feature_without_metadata_is_refused(
    corpus: dict[str, pd.DataFrame], positive: bool
) -> None:
    matrix = corpus["matrix"]
    matrix["METPO:unmapped"] = 0
    if positive:
        matrix.loc[matrix.index[-1], "METPO:unmapped"] = 1
    with pytest.raises(ValueError, match=r"Matrix feature.*metadata"):
        coverage(corpus)
    with pytest.raises(ValueError, match=r"Matrix feature.*metadata"):
        analysis.compute_matrix_trait_overlap(matrix, corpus["columns"], [FAC_AER, FAC_ANA])


@pytest.mark.parametrize("host", ["homo sapiens", "HOMO SAPIENS", " homo Sapiens "])
def test_case_variant_hosts_are_counted_without_admission(
    corpus: dict[str, pd.DataFrame], host: str
) -> None:
    corpus["source_experiments"].loc[0, "host_species"] = host
    # This variant has no reported taxa: the experiment diagnostic must still count it.
    extra = corpus["source_experiments"].iloc[[0]].copy()
    extra.loc[:, "study_id"] = "s4"
    extra.loc[:, "experiment_id"] = "e4"
    extra.loc[:, "body_site"] = "Blood"
    extra.loc[:, "condition"] = "Other condition"
    extra.loc[:, "sequencing_type"] = "ITS"
    corpus["source_experiments"] = pd.concat(
        [corpus["source_experiments"], extra], ignore_index=True
    )
    _, cascade, summary = coverage(corpus)
    assert summary["source_host_case_variant_experiments"] == 2
    assert summary["source_host_case_variant_profile_keys"] == 2
    human = next(row for row in cascade if row["step_id"] == "full_human_16s_wms")
    assert human["retained_count"] == 2
    assert human["signatures_count"] == 3
    summaries, pairs = overlap(corpus)
    human_pd = next(row for row in summaries if row["scope"] == "export_human_hosts")
    assert human_pd["source_host_case_variant_experiments"] == 2
    assert human_pd["source_host_case_variant_profile_keys"] == 2
    assert human_pd["n_studies"] == 2
    assert all(
        "s1" not in (row["study_a_id"], row["study_b_id"])
        for row in pairs
        if row["scope"] == "export_human_hosts"
    )
    assert next(row for row in summaries if row["scope"] == "full_export")["n_studies"] == 3


@pytest.mark.parametrize(
    ("host", "method", "human_studies", "human_pd_studies"),
    [(" Homo sapiens ", "16S", 3, 3), ("Homo sapiens", "ITS", 2, 3), ("Mus musculus", "16S", 2, 2)],
)
def test_host_case_diagnostics_do_not_reclassify_other_exclusions(
    corpus: dict[str, pd.DataFrame],
    host: str,
    method: str,
    human_studies: int,
    human_pd_studies: int,
) -> None:
    corpus["source_experiments"].loc[0, ["host_species", "sequencing_type"]] = [host, method]
    _, cascade, summary = coverage(corpus)
    assert summary["source_host_case_variant_experiments"] == 0
    assert summary["source_host_case_variant_profile_keys"] == 0
    human = next(row for row in cascade if row["step_id"] == "full_human_16s_wms")
    assert human["retained_count"] == human_studies
    summaries, _ = overlap(corpus)
    human_pd = next(row for row in summaries if row["scope"] == "export_human_hosts")
    assert human_pd["source_host_case_variant_experiments"] == 0
    assert human_pd["source_host_case_variant_profile_keys"] == 0
    assert human_pd["n_studies"] == human_pd_studies


def test_legitimate_inactive_rows_are_counted_and_repeated_occurrences_deduplicated(
    corpus: dict[str, pd.DataFrame],
) -> None:
    rows, _, summary = coverage(corpus)
    assert summary["tier12_input_occurrence_rows"] == 6
    assert summary["tier12_excluded_occurrence_rows"] == 1
    assert summary["tier12_excluded_profile_keys"] == 1
    assert summary["tier12_raw_occurrence_rows"] == 5
    assert rows[0]["total_mentions"] == 4
    assert rows[0]["total_taxa"] == 3


@pytest.fixture
def asymmetric_population(corpus: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """Keep e1 increased active while its decreased-only taxon remains in the input."""
    keys = corpus["tier12_profile_keys"]
    inactive_direction = (
        keys.study_id.eq("s1") & keys.experiment_id.eq("e1") & keys.inc_dec.eq("dec")
    )
    corpus["tier12_profile_keys"] = keys.loc[~inactive_direction].copy()
    return corpus


def test_inactive_direction_is_excluded_from_coverage(
    asymmetric_population: dict[str, pd.DataFrame],
) -> None:
    keys = asymmetric_population["tier12_profile_keys"]
    assert set(keys.itertuples(index=False, name=None)) == {
        ("s1", "e1", "inc"),
        ("s2", "e2", "inc"),
        ("s2", "e2", "dec"),
    }
    rows, _, summary = coverage(asymmetric_population)
    assert summary["tier12_input_occurrence_rows"] == 6
    assert summary["tier12_excluded_occurrence_rows"] == 2
    assert summary["tier12_excluded_profile_keys"] == 2
    assert summary["tier12_raw_occurrence_rows"] == 4
    assert summary["tier12_signatures"] == 3
    assert summary["tier12_studies"] == 2
    assert summary["tier12_experiments"] == 2
    # Only taxa 10 and 12 contribute: dec-only taxon 11 and inactive-study 13 do not.
    assert rows[0]["total_taxa"] == rows[0]["labeled_taxa"] == 2
    assert rows[0]["total_mentions"] == rows[0]["labeled_mentions"] == 3


def test_inactive_direction_is_excluded_from_tier12_pd_overlap(
    asymmetric_population: dict[str, pd.DataFrame],
) -> None:
    summaries, pairs = overlap(asymmetric_population)
    tier12 = next(row for row in summaries if row["scope"] == "tier12_paired")
    assert tier12["mean_jaccard"] == tier12["median_jaccard"] == 0.5
    assert tier12["excluded_outside_active_profile_rows"] == 2
    assert tier12["excluded_outside_active_profile_keys"] == 2
    assert (tier12["n_studies"], tier12["n_experiments"], tier12["n_study_pairs"]) == (2, 2, 1)
    tier12_pairs = [row for row in pairs if row["scope"] == "tier12_paired"]
    assert [(row["study_a_id"], row["study_b_id"]) for row in tier12_pairs] == [("s1", "s2")]
    [pair] = tier12_pairs
    assert (pair["study_a_taxa"], pair["study_b_taxa"]) == (1, 2)
    assert (pair["intersection_taxa"], pair["union_taxa"]) == (1, 2)
    # Full-export evidence retains the decreased row; only the active population excludes it.
    full = next(row for row in summaries if row["scope"] == "full_export")
    assert full["n_studies"] == 3
    full_pair = next(
        row
        for row in pairs
        if row["scope"] == "full_export" and (row["study_a_id"], row["study_b_id"]) == ("s1", "s2")
    )
    assert (full_pair["study_a_taxa"], full_pair["study_b_taxa"]) == (2, 2)
    assert (full_pair["intersection_taxa"], full_pair["union_taxa"]) == (1, 3)


def test_pd_profile_keys_exclude_a_decoy_taxon_that_changes_the_overlap(
    corpus: dict[str, pd.DataFrame],
) -> None:
    # Extra e4 belongs to the same PD study, but has no active profile key.
    experiment = {**corpus["tier12_experiments"].iloc[0].to_dict(), "experiment_id": "e4"}
    taxon = {
        **corpus["tier12_taxa"].iloc[0].to_dict(),
        "experiment_id": "e4",
        "taxon_id": "NCBITaxon:99",
    }
    for name in ("tier12_experiments", "source_experiments"):
        corpus[name] = pd.concat([corpus[name], pd.DataFrame([experiment])], ignore_index=True)
    for name in ("tier12_taxa", "source_taxa"):
        corpus[name] = pd.concat([corpus[name], pd.DataFrame([taxon])], ignore_index=True)
    summaries, pairs = overlap(corpus)
    tier12 = summaries[0]
    assert tier12["mean_jaccard"] == pytest.approx(1 / 3)
    assert tier12["excluded_outside_active_profile_rows"] == 2
    assert tier12["excluded_outside_active_profile_keys"] == 2
    row = next(row for row in pairs if row["scope"] == "tier12_paired")
    assert (row["intersection_taxa"], row["union_taxa"]) == (1, 3)


@pytest.mark.parametrize("site", ["feces", " Feces ", "FECES", "Feces,Blood", "Oral cavity,Feces"])
def test_literal_stool_rule_controls_cascade_and_every_pd_scope(
    corpus: dict[str, pd.DataFrame],
    site: str,
) -> None:
    for name in ("tier12_experiments", "source_experiments"):
        corpus[name].loc[1, "body_site"] = site
    _, _, summary = coverage(corpus)
    assert summary["tier12_stool_studies"] == 1
    assert summary["tier12_stool_experiments"] == 1
    summaries, pairs = overlap(corpus)
    assert summaries[0]["n_studies"] == 1
    assert summaries[0]["mean_jaccard"] is None
    assert all("s2" not in (row["study_a_id"], row["study_b_id"]) for row in pairs)


def test_shared_host_method_and_medium_rules(corpus: dict[str, pd.DataFrame]) -> None:
    corpus["source_experiments"].loc[1, ["host_species", "sequencing_type"]] = [
        "Homo sapiens",
        "ITS",
    ]
    corpus["source_experiments"].loc[2, "host_species"] = "homo sapiens"
    _, cascade, _ = coverage(corpus)
    human = next(row for row in cascade if row["step_id"] == "full_human_16s_wms")
    assert human["retained_count"] == 1
    summaries, _ = overlap(corpus)
    assert next(row for row in summaries if row["scope"] == "export_human_hosts")["n_studies"] == 2
    mixed_case_medium = "METPO:999999|MediaDive.Medium:42"
    matrix = corpus["matrix"].rename(columns={FAC_ANA: mixed_case_medium})
    columns = corpus["columns"].replace({FAC_ANA: mixed_case_medium})
    with pytest.raises(ValueError, match="MediaDive"):
        analysis.compute_matrix_trait_overlap(matrix, columns, [FAC_AER, mixed_case_medium])


def test_empty_population_has_defined_schemas_and_missing_inference(
    corpus: dict[str, pd.DataFrame],
) -> None:
    for name in corpus:
        if name not in {"matrix", "columns"}:
            corpus[name] = corpus[name].iloc[:0]
    rows, _, summary = coverage(corpus)
    assert summary["tier12_studies"] == 0
    assert rows[0]["total_taxa"] == 0
    assert rows[0]["labeled_taxa_pct"] is None
    summaries, pairs = overlap(corpus)
    assert len(summaries) == 4
    assert pairs == []
    assert all(row["mean_jaccard"] is None for row in summaries)


def test_compound_experiment_keys_do_not_mix_studies(corpus: dict[str, pd.DataFrame]) -> None:
    for name in corpus:
        if name not in {"matrix", "columns"}:
            corpus[name] = corpus[name].copy()
            corpus[name]["experiment_id"] = "shared-local-id"
    _, _, summary = coverage(corpus)
    assert summary["tier12_experiments"] == 2
    summaries, _ = overlap(corpus)
    assert summaries[0]["n_studies"] == 2
    assert summaries[0]["n_experiments"] == 2
    assert summaries[0]["mean_jaccard"] == pytest.approx(1 / 3)


def test_undefined_metrics_are_null_but_real_zero_remains_zero(
    corpus: dict[str, pd.DataFrame],
) -> None:
    matrix = corpus["matrix"].copy()
    matrix[[FAC_AER, FAC_ANA]] = 0
    empty = analysis.compute_matrix_trait_overlap(matrix, corpus["columns"], [FAC_AER, FAC_ANA])[0]
    assert empty["jaccard_fraction"] is None
    assert empty["jaccard_pct"] is None
    matrix.loc[0, FAC_AER], matrix.loc[1, FAC_ANA] = 1, 1
    zero = analysis.compute_matrix_trait_overlap(matrix, corpus["columns"], [FAC_AER, FAC_ANA])[0]
    assert zero["jaccard_fraction"] == 0.0
    assert json.loads(json.dumps(empty))["jaccard_fraction"] is None
    assert pd.DataFrame([empty]).to_csv(sep="\t", index=False).splitlines()[1].endswith("\t\t")
    corpus["tier12_profile_keys"] = corpus["tier12_profile_keys"].iloc[:0]
    rows, _, summary = coverage(corpus)
    assert rows[0]["labeled_taxa_pct"] is None
    assert rows[0]["labeled_median_studies"] is None
    assert summary["tier12_studies"] == 0


@pytest.fixture
def bundle(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, corpus: dict[str, pd.DataFrame]) -> Any:
    root = tmp_path / "core"
    root.mkdir()
    names = {
        "taxa_matrix": "matrix",
        "feature_metadata": "columns",
        "count_profiles": "tier12_profile_keys",
        "experiments": "tier12_experiments",
        "reported_taxa": "tier12_taxa",
        "source_experiments": "source_experiments",
        "source_reported_taxa": "source_taxa",
        "feature_catalog": "columns",
    }
    paths = {name: root / f"{name}.tsv" for name in names}
    for name, key in names.items():
        corpus[key].to_csv(paths[name], sep="\t", index=False)
    paths["contrast_selection"] = root / "contrast_selection.tsv"
    paths["contrast_selection"].write_text("unread\n")
    manifest = root / "run.json"
    manifest.write_text("{}\n")
    run = SimpleNamespace(
        run_id="core-fixture",
        git_commit="core-parent",
        contrast_population=POPULATION,
        direction_basis="disease_vs_comparison",
        parameters=[SimpleNamespace(name="strain_rollup", value="phenotype")],
        definition=SimpleNamespace(
            id="fixture", products=[SimpleNamespace(id=name) for name in paths]
        ),
    )
    value = SimpleNamespace(run=run, manifest=manifest, path=lambda name: paths[name])
    monkeypatch.setattr(consumers, "open_core_bundle", lambda _path: value)
    monkeypatch.setattr(consumers, "preflight_consumer_outputs", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(consumers, "write_consumer_record", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(runtime, "repository_state", lambda _root: ("consumer-current", False))
    traits = tmp_path / "traits.txt"
    traits.write_text(f"{FAC_AER}\n{FAC_ANA}\n")
    return SimpleNamespace(value=value, traits=traits, output=tmp_path / "output")


def run_bundle(bundle: Any) -> dict[str, Any]:
    return analysis.run_matrix_overlap_coverage(
        bundle.value.manifest,
        bundle.output,
        trait_list_path=bundle.traits,
        allow_temporary_output=True,
    )


def test_matrix_api_and_cli_guard_share_the_expanded_home_destination(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from mtp.workflow import output_roots

    supplied = Path("~/matrix-output-test")
    expanded = tmp_path / "synthetic-home" / "matrix-output-test"
    real_expanduser = Path.expanduser

    def expanduser(path: Path) -> Path:
        return expanded if path == supplied else real_expanduser(path)

    monkeypatch.setattr(Path, "expanduser", expanduser)
    monkeypatch.setattr(
        output_roots, "system_temporary_directories", lambda: [tmp_path / "scratch"]
    )
    assert analysis.validate_durable_output_dir(supplied) == expanded
    guarded: list[Path] = []
    consumed: list[Path] = []
    real_guard = output_roots.require_durable_output

    def guard(parser: Any, output: Path, *, allowed: bool) -> None:
        guarded.append(output)
        real_guard(parser, output, allowed=allowed)

    def consume(_core: Path, output: Path, **_kwargs: Any) -> dict[str, Any]:
        consumed.append(output)
        return {}

    monkeypatch.setattr(analysis, "require_durable_output", guard)
    monkeypatch.setattr(analysis, "run_matrix_overlap_coverage", consume)
    assert (
        analysis.main(
            [
                "--core-run",
                str(tmp_path / "core.json"),
                "--output-dir",
                str(supplied),
                "--trait-list",
                str(tmp_path / "traits.txt"),
            ]
        )
        == 0
    )
    assert guarded == [expanded]
    assert consumed == [expanded]
    assert not expanded.parent.exists()


def test_wrong_direction_basis_is_refused_before_publication(bundle: Any) -> None:
    bundle.value.run.direction_basis = "as_reported"
    with pytest.raises(ValueError, match="direction_basis"):
        run_bundle(bundle)
    assert not bundle.output.exists()


def test_contradictory_direction_parameter_is_refused_independently(bundle: Any) -> None:
    assert bundle.value.run.direction_basis == "disease_vs_comparison"
    bundle.value.run.parameters.append(SimpleNamespace(name="direction_basis", value="as_reported"))
    with pytest.raises(ValueError, match="direction_basis parameter contradicts"):
        run_bundle(bundle)
    assert not bundle.output.exists()


def test_contradictory_population_parameter_is_refused_independently(bundle: Any) -> None:
    assert bundle.value.run.contrast_population == POPULATION
    bundle.value.run.parameters.append(
        SimpleNamespace(name="contrast_population", value="all_paired")
    )
    with pytest.raises(ValueError, match="contrast_population parameter contradicts"):
        run_bundle(bundle)
    assert not bundle.output.exists()


@pytest.mark.parametrize(
    ("field", "bad"), [("contrast_population", "all_paired"), ("direction_basis", "as_reported")]
)
def test_authoritative_run_fields_cannot_be_overridden_by_parameters(
    bundle: Any, field: str, bad: str
) -> None:
    original = getattr(bundle.value.run, field)
    setattr(bundle.value.run, field, bad)
    bundle.value.run.parameters.append(SimpleNamespace(name=field, value=original))
    with pytest.raises(ValueError, match=r"contradict|differs|requires"):
        run_bundle(bundle)
    assert not bundle.output.exists()


@pytest.mark.parametrize("missing", ["source_experiments", "source_reported_taxa"])
def test_missing_full_source_product_never_falls_back(bundle: Any, missing: str) -> None:
    bundle.value.run.definition.products = [
        item for item in bundle.value.run.definition.products if item.id != missing
    ]
    with pytest.raises(ValueError, match=r"source_experiments|source_reported_taxa"):
        run_bundle(bundle)
    assert not bundle.output.exists()


@pytest.mark.parametrize("dirty", [False, True])
def test_automatic_commit_is_consumer_and_unread_product_is_not_declared_consumed(
    bundle: Any,
    monkeypatch: pytest.MonkeyPatch,
    dirty: bool,
) -> None:
    monkeypatch.setattr(runtime, "repository_state", lambda _root: ("consumer-current", dirty))
    run_bundle(bundle)
    record = json.loads((bundle.output / "provenance.json").read_text())
    assert record["code_commit"] == "consumer-current" + ("+dirty" if dirty else "")
    assert record["core_git_commit_sha"] == "core-parent"
    assert "contrast_selection" not in record["consumed_products"]
    assert {
        "src/mtp/matrix_overlap_inputs.py",
        "src/mtp/human_host_corpus.py",
        "src/mtp/workflow/coverage.py",
        "src/mtp/workflow/output_roots.py",
    } <= set(record["scripts"])
