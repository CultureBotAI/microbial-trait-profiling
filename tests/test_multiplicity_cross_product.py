"""Producer-ledger consistency beyond individually plausible rows (#976)."""

import runpy
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from scipy.stats import false_discovery_control

from mtp.multiplicity.families import jed_family, noel_family, rows
from mtp.multiplicity.native import native_family
from mtp.multiplicity.ora import native_bh_rows, ora_comments, ora_family
from mtp.multiplicity.sources import (
    VerifiedSources,
    recorded_arm_output,
    recorded_comparison_output,
)
from mtp.paths import file_sha256
from mtp.trait_composition import compute_log_ratio_shifts, write_log_ratio_outputs
from mtp.trait_guilds import TABLE1_TRAITS, load_profiles
from mtp.trait_screen import candidate_slots, run_screen, write_screen_outputs
from multiplicity_fixture import fixture_sources, write_table


@pytest.mark.parametrize("analysis", ["noel", "jed", "native", "ora"])
def test_cross_product_status_and_correction_tampering_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, analysis: str
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    adapter = {"noel": noel_family, "jed": jed_family, "native": native_family, "ora": ora_family}[
        analysis
    ]
    family = adapter(sources)
    control = list(family.observations())
    assert control
    assert any(item.raw.status == "tested" for item in control)
    comments = None
    if analysis == "native":
        path = recorded_comparison_output(
            sources.comparison, "trait_profile_paired_contrast_trait_stats_pooled.tsv"
        )
        records = list(rows(path))
        row = next(row for row in records if row["study_test_status"] == "tested")
        row.update(
            study_test_status="all_tied",
            study_mean_increased_minus_decreased="0",
            study_std_increased_minus_decreased="0",
            study_p_value="",
            study_q_value="1",
        )
    elif analysis == "ora":
        path = recorded_arm_output(
            sources.arms["ilyas"], "enrichment", "trait_enrichment_all_results.tsv"
        )
        records = list(rows(path))
        comments = ora_comments(path)
        records[0]["p_value"] = "0.0000001"
    else:
        path = recorded_arm_output(
            sources.arms[analysis],
            "directional_screen" if analysis == "noel" else "composition",
            "trait_direction_screen.tsv" if analysis == "noel" else "trait_log_ratio_stats.tsv",
        )
        records = list(rows(path))
        status = "test_status" if analysis == "noel" else "status"
        row = next(row for row in records if row[status] == "tested")
        total = sum(int(row[f"studies_{kind}"]) for kind in ("up", "down", "tied"))
        row.update(studies_up="0", studies_down="0", studies_scored="0")
        if analysis == "noel":
            row.update(
                test_status="all_tied", studies_tied=str(total), binomial_p="", by_q_nominal=""
            )
        else:
            row.update(
                status="no_positive_pairs",
                studies_tied="0",
                studies_eligible="0",
                studies_excluded=row["studies_total"],
                studies_partially_observed="0",
                eligible_experiments="0",
                excluded_zero_experiments=row["paired_experiments"],
                sign_p="",
                by_q="",
            )
    write_table(path, records, comments)
    # Byte-binding is deliberately absent from these standalone fixture adapters:
    # the independently consumed core/ledger must reject the plausible export.
    with pytest.raises(ValueError, match=r"core|ledger|moments|native.*BH"):
        list(family.observations())


def rewrite_native_export(sources: VerifiedSources) -> None:
    """Use the actual producer after a synthetic core change, without guessed statistics."""
    producer = runpy.run_path(
        str(Path(__file__).parents[1] / "scripts/create_trait_profile_infographic.py")
    )
    profile, experiments = sources.product("count_profiles"), sources.product("experiments")
    labels = sources.product("feature_metadata")
    pair_path = recorded_comparison_output(
        sources.comparison, "trait_profile_paired_contrast_pair_metadata.tsv"
    )
    table = recorded_comparison_output(
        sources.comparison, "trait_profile_paired_contrast_trait_stats_pooled.tsv"
    )
    hashes = {
        "profile_sha256": file_sha256(profile),
        "experiments_sha256": file_sha256(experiments),
        "trait_labels_sha256": file_sha256(labels),
    }
    pairs = producer["core_pair_metadata"](profile, experiments)
    for key, value in hashes.items():
        pairs[key] = value
    pairs.to_csv(pair_path, sep="\t")
    stats = producer["compute_stratum_trait_stats"](
        producer["load_profile_matrix"](profile),
        pairs,
        producer["load_trait_labels"](labels),
        stratum_kind="pooled",
        stratum="pooled",
    )
    for key, value in {**hashes, "pair_metadata_sha256": file_sha256(pair_path)}.items():
        stats[key] = value
    stats.to_csv(table, sep="\t", index=False)


@pytest.mark.parametrize(
    ("analysis", "reason"),
    [
        ("noel", "all_tied"),
        ("noel", "no_eligible_studies"),
        ("jed", "all_tied"),
        ("jed", "no_positive_pairs"),
        ("native", "all_tied"),
        ("native", "insufficient_support"),
        ("native", "insufficient_variance"),
    ],
)
def test_actual_exporters_preserve_valid_untestable_controls(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, analysis: str, reason: str
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    profile, experiments = sources.product("count_profiles"), sources.product("experiments")
    features = sources.features["count_profiles"]
    records = list(rows(profile))
    if analysis == "native":
        if reason == "all_tied":
            decreased = {row["study_id"]: row for row in records if row["inc_dec"] == "dec"}
            for row in records:
                row.update({feature: decreased[row["study_id"]][feature] for feature in features})
        elif reason == "insufficient_variance":
            for row in records:
                row.update(dict.fromkeys(features, "0"))
                row[features[0]] = str(int(row["inc_dec"] == "inc"))
                row[features[1]] = str(int(row["inc_dec"] == "dec"))
        else:
            metadata = list(rows(experiments))
            for row in records:
                row["experiment_id"], row["study_id"] = row["study_id"], "one-study"
            for row in metadata:
                row["experiment_id"], row["study_id"] = row["study_id"], "one-study"
            write_table(experiments, metadata)
        write_table(profile, records)
        rewrite_native_export(sources)
        observations = list(native_family(sources).observations())
    else:
        for row in records:
            row.update(dict.fromkeys(features, "1" if reason == "all_tied" else "0"))
        write_table(profile, records)
        arm = sources.arms[analysis]
        if analysis == "noel":
            paired, corpus = load_profiles(profile, experiments, candidate_slots(features))
            write_screen_outputs(
                run_screen(paired, features),
                Path(arm.manifest_path).parent,
                profile_path=profile,
                experiments_path=experiments,
                full_dump_path=profile.parent / "raw.csv",
                corpus_stats=corpus,
            )
            observations = list(noel_family(sources).observations())
        else:
            paired, corpus = load_profiles(profile, experiments, TABLE1_TRAITS)
            write_log_ratio_outputs(
                compute_log_ratio_shifts(paired),
                Path(arm.manifest_path).parent,
                profile_path=profile,
                experiments_path=experiments,
                full_dump_path=profile.parent / "raw.csv",
                corpus_stats=corpus,
            )
            observations = list(jed_family(sources).observations())
    assert observations
    assert {item.raw.reason for item in observations} == {reason}
    assert all(item.raw.p_value is None for item in observations)
    assert {item.raw.reason_vocabulary for item in observations} == {"analysis_reason_v1"}


def test_heterogeneous_actual_ora_exports_use_complete_separate_bh_families(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = fixture_sources(tmp_path / "inputs", monkeypatch)
    features = sources.features["ora_profiles"]
    matrix, profiles = sources.product("taxa_matrix"), sources.product("ora_profiles")
    background: list[dict[str, Any]] = [
        {
            "taxon_id": f"NCBITaxon:{taxon + 1}",
            "taxon_name": f"Fixture species {taxon + 1}",
            "taxonomic_rank": "species",
            **{
                feature: int(taxon < column + 1) if column < len(features) - 1 else 0
                for column, feature in enumerate(features)
            },
        }
        for taxon in range(8)
    ]
    write_table(matrix, background)
    records = list(rows(profiles))
    for index, row in enumerate(records):
        selected = [(index + offset) % len(background) for offset in range(1 + index % 4)]
        row["taxa_count"] = str(len(selected))
        row.update(
            {
                feature: str(sum(int(background[taxon][feature]) for taxon in selected))
                for feature in features
            }
        )
    write_table(profiles, records)
    producer = runpy.run_path(str(Path(__file__).parents[1] / "scripts/trait_enrichment.py"))[
        "main"
    ]
    for tail, step in (("upper", "enrichment"), ("lower", "depletion")):
        path = recorded_arm_output(sources.arms["ilyas"], step, f"trait_{step}_all_results.tsv")
        producer(
            [
                "--input",
                str(profiles),
                "--taxa-matrix",
                str(matrix),
                "--output",
                str(path.with_suffix(".txt")),
                "--all-results-output",
                str(path),
                "--tail",
                tail,
                "--max-fold",
                "0.5",
            ]
        )
    family = ora_family(sources)
    imported = list(family.observations())
    assert len(imported) == 2 * len(records) * (len(features) - 1)
    by_tail = {}
    for tail in ("upper", "lower"):
        values = [item for item in imported if item.native_family == f"ilyas-{tail}-global"]
        probabilities = np.array([item.raw.p_value for item in values], dtype=float)
        assert len(set(probabilities)) > 3
        assert [item.native_q for item in values] == pytest.approx(
            false_discovery_control(probabilities, method="bh"), rel=1e-14, abs=0
        )
        by_tail[tail] = probabilities
        assert all(float(item.row["p_value"]) == item.raw.p_value for item in values)
    assert not np.array_equal(by_tail["upper"], by_tail["lower"])
    # Display sorting cannot define or partition a correction family.
    upper = recorded_arm_output(
        sources.arms["ilyas"], "enrichment", "trait_enrichment_all_results.tsv"
    )
    original = list(rows(upper))
    write_table(upper, original[::-1], ora_comments(upper))
    assert len(list(ora_family(sources).observations())) == len(imported)


@pytest.mark.parametrize("mutation", ["duplicate", "missing", "unknown_profile", "missing_q"])
def test_native_bh_validates_exact_membership_and_missing_conventions(
    tmp_path: Path, mutation: str
) -> None:
    path = tmp_path / "tail.tsv"
    records = [
        {
            "study_id": "s",
            "experiment_id": "e",
            "inc_dec": direction,
            "trait": "t",
            "p_value": "0.2",
            "fdr_q": "0.2",
        }
        for direction in ("inc", "dec")
    ]
    if mutation == "duplicate":
        records[1] = records[0].copy()
    elif mutation == "missing":
        records.pop()
    elif mutation == "unknown_profile":
        records[0]["study_id"] = "other"
    else:
        records[0]["fdr_q"] = ""
    write_table(path, records)
    with pytest.raises(ValueError, match=r"duplicate|truncated|undeclared|lacks"):
        list(
            native_bh_rows(
                path,
                expected_sha256=None,
                profile_keys={("s", "e", "inc"), ("s", "e", "dec")},
                features={"t"},
            )
        )
