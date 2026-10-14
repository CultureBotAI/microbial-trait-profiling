"""Reproduce #174 against frozen inputs, changing only the composition denominator."""

from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

import create_study_trait_profile_heatmap as heatmap
import create_trait_profile_infographic as corrected
from bsdbng.table1_trait_maps import compute_direction_shifts, write_shift_table
from bsdbng.trait_composition import DENOMINATOR_DEFINITION, retained_trait_counts
from bsdbng.trait_guilds import TABLE1_TRAITS, load_profiles

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / ".tmp/issue174"
OUTPUT = ROOT / "trait-visualization/denominator-harmonization"
BASELINE_COMMIT = "a2b2df98e9beec98b772c9cb5af7aa090d033a98"
INPUT_BLOB = "d837f449e7e967c478d10f92938923e57eec9a19"
SIX_IDS = [
    "METPO:1000608",
    "METPO:2000103|GO:0030245",
    "METPO:1000668",
    "METPO:1000683",
    "METPO:2000103|pathways:arsenate_reduction",
    "METPO:1000623",
]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValueError(f"Cannot import pinned source: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_table(frame: pd.DataFrame, name: str, *, index: bool = False) -> None:
    path = OUTPUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = {"method": "gzip", "mtime": 0} if name.endswith(".gz") else None
    frame.to_csv(path, sep="\t", index=index, na_rep="NA", compression=compression)


def sha256(path: Path) -> str:
    content = path.read_bytes()
    if b"\r\n" in content:
        content = content.replace(b"\r\n", b"\n")
    return hashlib.sha256(content).hexdigest()


def strata(pairs: pd.DataFrame):
    yield "pooled", "pooled", "pooled", pairs
    yield "feces", "body_site_group", "Feces", pairs[pairs.body_site_group.eq("Feces")]
    yield (
        "feces_excluded",
        "body_site_group_excluded",
        "Feces excluded",
        pairs[~pairs.body_site_group.eq("Feces")],
    )
    for stratum_id, label, group in corrected.body_site_strata(
        pairs, corrected.MIN_STRATIFIED_PAIRS
    ):
        if stratum_id != "feces":
            yield stratum_id, "body_site_group", label, group


def compare_stratum(raw, labels, pairs, baseline, token, kind, label):
    ids = corrected.profile_ids_for_pairs(pairs)
    subset = raw.loc[ids]
    records = []
    stages = {}
    for stage, module in [("before", baseline), ("after", corrected)]:
        selected = module.select_stratum_traits(subset, labels)
        zscores = module.build_zscore_matrix(subset, selected.index)
        contrasts = module.build_pair_matrix_from_profiles(zscores, pairs)
        stats = module.compute_trait_stats(
            contrasts,
            selected,
            module.annotation_counts_for_traits(subset, selected.index),
            stratum_kind=kind,
            stratum=label,
            selected_profile_count=len(subset),
        ).set_index("trait")
        write_table(stats, f"{stage}/trait_stats_{token}.tsv", index=True)
        write_table(selected, f"{stage}/selection_{token}.tsv")
        write_table(contrasts, f"{stage}/paired_contrasts_{token}.tsv.gz", index=True)
        stages[stage] = (stats, selected, contrasts)
        records.append(
            {
                "stratum": token,
                "version": stage,
                "pairs": len(pairs),
                "profiles": len(subset),
                "studies": pairs.study_id.nunique(),
                "correction_family": f"composition:{token}",
                "tested_traits": len(stats),
                "q_le_0_001": int(stats.q_value.le(0.001).sum()),
                "q_le_0_05": int(stats.q_value.le(0.05).sum()),
            }
        )
    all_traits = retained_trait_counts(raw).columns
    comparison = pd.DataFrame(index=all_traits)
    comparison.index.name = "trait"
    comparison["trait_label"] = labels.trait_label.reindex(all_traits)
    comparison["stratum"] = token
    comparison["correction_family"] = f"composition:{token}"
    for stage, (stats, _, _) in stages.items():
        comparison[f"{stage}_eligible"] = comparison.index.isin(stats.index)
        for field in [
            "mean_increased_minus_decreased",
            "direction",
            "p_value",
            "q_value",
            "composition_variance",
            "selection_score",
            "selected_trait_rank",
            "profile_support",
        ]:
            comparison[f"{stage}_{field}"] = stats[field].reindex(all_traits)
        comparison[f"{stage}_family_size"] = len(stats)
        for cutoff in [0.001, 0.05]:
            comparison[f"{stage}_q_le_{cutoff}"] = comparison[f"{stage}_q_value"].le(cutoff)
    comparison["direction_changed"] = (
        comparison.before_eligible
        & comparison.after_eligible
        & comparison.before_direction.ne(comparison.after_direction)
    )
    # Compare the rerun baseline to the published stratum, rather than assuming reproducibility.
    published = pd.read_csv(
        WORK
        / "baseline/trait-visualization/trait-profile-investigation"
        / f"trait_profile_paired_contrast_trait_stats_{token}.tsv",
        sep="\t",
    ).set_index("trait")
    actual = stages["before"][0].sort_index()
    published = published.sort_index()
    if not actual.index.equals(published.index):
        raise ValueError(f"Published eligibility not reproduced for {token}")
    for field in ["mean_increased_minus_decreased", "p_value", "q_value"]:
        np.testing.assert_allclose(actual[field], published[field], rtol=1e-9, atol=1e-12)
    return comparison.reset_index(), records, stages


def denominator_audit(raw, pairs, metadata):
    retained = retained_trait_counts(raw)
    totals = retained.sum(axis=1)
    media = raw.sum(axis=1) - totals
    paired_ids = corrected.profile_ids_for_pairs(pairs)
    scope_rows = []
    for scope, ids in [("full_matrix", raw.index), ("paired", pd.Index(paired_ids))]:
        scope_rows.append(
            {
                "scope": scope,
                "profiles": len(ids),
                "retained_trait_columns": len(retained.columns),
                "mediadive_columns": len(raw.columns) - len(retained.columns),
                "profiles_with_mediadive": int(media.loc[ids].gt(0).sum()),
                "zero_retained_signal": int(totals.loc[ids].eq(0).sum()),
                "mediadive_only_profiles": int(
                    (totals.loc[ids].eq(0) & media.loc[ids].gt(0)).sum()
                ),
                "retained_total": int(totals.loc[ids].sum()),
                "media_total": int(media.loc[ids].sum()),
            }
        )
    scope_frame = pd.DataFrame(scope_rows)
    write_table(scope_frame, "input_scope.tsv")

    # Run Noel's actual loader on the same frozen keys; labels do not enter its denominator.
    experiment_metadata = metadata[
        ["study_id", "experiment_id", "condition", "body_site", "host_species"]
    ].drop_duplicates(["study_id", "experiment_id"])
    experiment_path = WORK / "matched_experiments.tsv"
    experiment_metadata.to_csv(experiment_path, sep="\t", index=False)
    noel, corpus = load_profiles(
        WORK / "study_trait_profiles_wide.tsv", experiment_path, TABLE1_TRAITS
    )
    noel_totals = pd.Series({f"{p.study_id}|{p.experiment_id}|{p.inc_dec}": p.total for p in noel})
    if set(noel_totals.index) != set(paired_ids):
        raise ValueError("Noel and composition paired profile keys differ")
    # Apply Ilyas's pinned column selection and positive-count accumulation literally.
    ora_columns = [c for c in raw.columns if "mediadive.medium:" not in c]
    ora_totals = raw[ora_columns].where(raw[ora_columns].gt(0), 0).sum(axis=1)
    audit = pd.DataFrame(
        {
            "legacy_total": raw.sum(axis=1),
            "harmonized_total": totals,
            "mediadive_total": media,
            "ilyas_total": ora_totals,
            "noel_total": noel_totals,
            "paired": raw.index.isin(paired_ids),
        }
    )
    audit.index.name = "profile_id"
    np.testing.assert_array_equal(audit.harmonized_total, audit.ilyas_total)
    np.testing.assert_array_equal(
        audit.loc[paired_ids].harmonized_total, noel_totals.reindex(paired_ids)
    )
    fraction_checks = []
    for profile in noel:
        key = f"{profile.study_id}|{profile.experiment_id}|{profile.inc_dec}"
        for slug, terms in TABLE1_TRAITS.items():
            trait = next(c for c in SIX_IDS if c.rsplit("|", 1)[-1] in terms)
            value = float(retained.at[key, trait] / totals[key]) if totals[key] else 0.0
            difference = value - profile.fraction(slug)
            if difference != 0:
                raise ValueError(f"Noel fraction mismatch: {key}, {trait}")
            fraction_checks.append(
                {
                    "profile_id": key,
                    "trait": trait,
                    "composition_fraction": value,
                    "noel_fraction": profile.fraction(slug),
                    "absolute_difference": abs(difference),
                }
            )
    write_table(pd.DataFrame(fraction_checks), "table1_fraction_parity.tsv.gz")
    noel_dir = OUTPUT / "noel_reference"
    noel_dir.mkdir(exist_ok=True)
    write_shift_table(compute_direction_shifts(noel), noel_dir, "BugSigDBExports v1.3.1")
    write_table(audit, "denominator_parity.tsv", index=True)
    return scope_frame, corpus


def integration(comparison, raw, pairs):
    ora = pd.read_csv(WORK / "ora-q05.tsv", sep="\t", comment="#")
    legacy_ora = pd.read_csv(
        WORK / "ora/data/derived/trait_enrichment_results.tsv", sep="\t", comment="#"
    )
    keys = ["study_id", "experiment_id", "inc_dec", "trait"]
    now = ora[ora.fdr_q.le(0.001)].set_index(keys).sort_index()
    old = legacy_ora.set_index(keys).sort_index()
    if not now.index.equals(old.index):
        raise ValueError("Pinned ORA q<=0.001 membership was not reproduced")
    np.testing.assert_allclose(now.fdr_q, old.fdr_q, rtol=1e-9, atol=1e-12)
    ora["profile_id"] = ora.study_id + "|" + ora.experiment_id + "|" + ora.inc_dec
    records = []
    detail = []
    for token, _, _, group in strata(pairs):
        comp = comparison[comparison.stratum.eq(token)]
        ora_group = ora[ora.profile_id.isin(corrected.profile_ids_for_pairs(group))]
        for cutoff in [0.001, 0.05]:
            ora_set = set(ora_group.loc[ora_group.fdr_q.le(cutoff), "trait"])
            before = set(comp.loc[comp[f"before_q_le_{cutoff}"], "trait"])
            after = set(comp.loc[comp[f"after_q_le_{cutoff}"], "trait"])
            records.append(
                {
                    "stratum": token,
                    "q_cutoff": cutoff,
                    "ora_family": "global observed-trait-by-paired-profile tests; fold>=1.5",
                    "composition_family": f"composition:{token}",
                    "ora_traits": len(ora_set),
                    "before_composition_traits": len(before),
                    "after_composition_traits": len(after),
                    "before_overlap": len(before & ora_set),
                    "after_overlap": len(after & ora_set),
                    "lost_overlap": ";".join(sorted((before & ora_set) - after)),
                    "gained_overlap": ";".join(sorted((after & ora_set) - before)),
                }
            )
            for trait in retained_trait_counts(raw).columns:
                detail.append(
                    {
                        "stratum": token,
                        "q_cutoff": cutoff,
                        "trait": trait,
                        "ora_any_significant_profile": trait in ora_set,
                        "before_composition_significant": trait in before,
                        "after_composition_significant": trait in after,
                    }
                )
    write_table(pd.DataFrame(records), "cross_framework_overlap.tsv")
    write_table(pd.DataFrame(detail), "cross_framework_membership.tsv")
    return pd.DataFrame(records)


def refresh_caches(raw, metadata, labels, pairs, baseline):
    """Rebuild producer caches, including totals consumed by six-trait and PD scripts."""
    retained = retained_trait_counts(raw)
    current = metadata.loc[retained.index].copy()
    current["total_trait_signal"] = retained.sum(axis=1)
    current["trait_count"] = retained.gt(0).sum(axis=1)
    current["denominator_definition"] = DENOMINATOR_DEFINITION
    selected, selected_meta = heatmap.select_traits(raw, labels, 60, 6, 0.95, False)
    zscores = heatmap.build_plot_matrix(raw, selected)
    pair_matrix = corrected.build_pair_matrix_from_profiles(zscores, pairs)
    pair_metadata = pairs.copy()
    pair_metadata["mean_paired_trait_enrichment"] = pair_matrix.mean(axis=1)
    pair_metadata["max_abs_paired_trait_enrichment"] = pair_matrix.abs().max(axis=1)
    write_table(current, "after/cache/study_trait_profile_study_metadata.tsv", index=True)
    write_table(selected_meta, "after/cache/study_trait_profile_trait_metadata.tsv", index=True)
    display = heatmap.build_display_matrix(zscores, current, selected_meta)
    write_table(display, "after/cache/study_trait_profile_clustered_heatmap_matrix.tsv", index=True)
    write_table(zscores, "after/cache/profile_zscores_stable_ids.tsv.gz", index=True)
    paired_display = heatmap.build_paired_display_matrix(pair_matrix, pair_metadata, selected_meta)
    write_table(
        paired_display, "after/cache/study_trait_profile_paired_enrichment_matrix.tsv", index=True
    )
    write_table(
        pair_metadata, "after/cache/study_trait_profile_paired_enrichment_metadata.tsv", index=True
    )
    write_table(
        heatmap.summarize_required_paired_traits(pair_matrix, selected_meta),
        "after/cache/study_trait_profile_paired_required_trait_check.tsv",
    )
    # Six-trait consumers have a full-human-corpus z-score reference, not the pooled-pair reference.
    six_zscores = corrected.build_zscore_matrix(raw, pd.Index(SIX_IDS))
    write_table(six_zscores, "after/cache/six_trait_profile_zscores.tsv", index=True)
    legacy_six_zscores = baseline.build_zscore_matrix(raw, pd.Index(SIX_IDS))
    write_table(legacy_six_zscores, "before/cache/six_trait_profile_zscores.tsv", index=True)
    cached_pair_before = corrected.build_pair_matrix_from_profiles(legacy_six_zscores, pairs)
    cached_pair_after = corrected.build_pair_matrix_from_profiles(six_zscores, pairs)
    cache_impact = []
    for trait in SIX_IDS:
        cache_impact.append(
            {
                "trait": trait,
                "trait_label": labels.at[trait, "trait_label"],
                "reference_profiles": len(raw),
                "before_positive_profiles": int(legacy_six_zscores[trait].gt(0).sum()),
                "after_positive_profiles": int(six_zscores[trait].gt(0).sum()),
                "profile_zscore_sign_changes": int(
                    np.sign(legacy_six_zscores[trait]).ne(np.sign(six_zscores[trait])).sum()
                ),
                "paired_contrast_sign_changes": int(
                    np.sign(cached_pair_before[trait]).ne(np.sign(cached_pair_after[trait])).sum()
                ),
                "maximum_absolute_zscore_change": float(
                    (legacy_six_zscores[trait] - six_zscores[trait]).abs().max()
                ),
            }
        )
    write_table(pd.DataFrame(cache_impact), "cached_six_trait_impact.tsv")
    ordered = heatmap.order_profiles_by_direction(zscores, current)
    rendered, rendered_meta = heatmap.build_rendered_display_matrix(
        zscores.loc[ordered], current.loc[ordered], selected_meta, 320
    )
    heatmap.save_heatmap(
        rendered,
        rendered_meta,
        selected_meta,
        OUTPUT / "after/cache/study_trait_profile_clustered_heatmap.png",
        full_profile_count=len(current),
        full_study_count=current.study_id.nunique(),
    )
    ordered_pairs = heatmap.cluster_index(pair_matrix)
    rendered_pairs, rendered_pair_meta = heatmap.build_rendered_paired_enrichment_matrix(
        pair_matrix.loc[ordered_pairs], pair_metadata.loc[ordered_pairs], selected_meta, 320
    )
    heatmap.save_paired_enrichment_heatmap(
        rendered_pairs,
        rendered_pair_meta,
        selected_meta,
        OUTPUT / "after/cache/study_trait_profile_paired_enrichment_heatmap.png",
        full_pair_count=len(pairs),
        full_study_count=pairs.study_id.nunique(),
    )
    return current, selected_meta, zscores, pair_matrix


def publish_current_outputs(families: pd.DataFrame, pairs: pd.DataFrame) -> None:
    """Install corrected producer artifacts at the paths read by downstream figures."""
    target = ROOT / "trait-visualization"
    investigation = target / "trait-profile-investigation"
    investigation.mkdir(exist_ok=True)
    for path in (OUTPUT / "after/cache").glob("study_trait_profile_*"):
        shutil.copyfile(path, target / path.name)
    manifest = []
    definitions = {name: (kind, label) for name, kind, label, _ in strata(pairs)}
    roles = {
        "pooled": "pooled",
        "feces": "primary-analysis",
        "feces_excluded": "robustness-feces-excluded",
    }
    for row in families[families.version.eq("after")].itertuples(index=False):
        filename = f"trait_profile_paired_contrast_trait_stats_{row.stratum}.tsv"
        shutil.copyfile(
            OUTPUT / "after" / f"trait_stats_{row.stratum}.tsv", investigation / filename
        )
        kind, label = definitions[row.stratum]
        manifest.append(
            {
                "stratum_kind": kind,
                "stratum": label,
                "filename_token": row.stratum,
                "analysis_role": roles.get(row.stratum, "body-site"),
                "pair_count": row.pairs,
                "study_count": row.studies,
                "selected_trait_count": row.tested_traits,
                "correction_family": row.correction_family,
                "path": (investigation / filename).relative_to(ROOT).as_posix(),
            }
        )
    shutil.copyfile(
        investigation / "trait_profile_paired_contrast_trait_stats_pooled.tsv",
        investigation / "trait_profile_paired_contrast_trait_stats.tsv",
    )
    for name in ["condition_summary", "body_site_summary"]:
        shutil.copyfile(
            OUTPUT / "after" / f"{name}.tsv",
            investigation / f"trait_profile_paired_contrast_{name}.tsv",
        )
    shutil.copyfile(
        OUTPUT / "after/paired_contrast_infographic.png",
        investigation / "trait_profile_paired_contrast_infographic.png",
    )
    pd.DataFrame(manifest).to_csv(
        investigation / "trait_profile_paired_contrast_trait_stats_manifest.tsv",
        sep="\t",
        index=False,
    )


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    matrix_path = WORK / "study_trait_profiles_wide.tsv"
    data = matrix_path.read_bytes()
    actual_blob = hashlib.sha1(
        f"blob {len(data)}\0".encode() + data, usedforsecurity=False
    ).hexdigest()
    if actual_blob != INPUT_BLOB:
        raise ValueError(f"Input blob changed: {actual_blob}")
    baseline = load_module(
        WORK / "baseline/scripts/create_trait_profile_infographic.py", "legacy_composition"
    )
    raw = baseline.load_profile_matrix(matrix_path)
    meta_root = WORK / "baseline/trait-visualization"
    pairs = pd.read_csv(
        meta_root / "study_trait_profile_paired_enrichment_metadata.tsv", sep="\t", index_col=0
    )
    metadata = pd.read_csv(
        meta_root / "study_trait_profile_study_metadata.tsv", sep="\t", index_col=0
    )
    labels_path = OUTPUT / "inputs/taxa_trait_matrix_columns.tsv"
    expected_labels_hash = "70b8bc1c2889d0e451bcaf604e0a986d4c29ffb802efeae5bdc64b67fc446ed3"
    if sha256(labels_path) != expected_labels_hash:
        raise ValueError("Frozen trait labels changed")
    labels = corrected.load_trait_labels(labels_path)
    if not raw.index.is_unique or set(raw.index) != set(metadata.index):
        raise ValueError("Pinned matrix and metadata profile keys differ")
    if labels.reindex(raw.columns).trait_label.isna().any():
        raise ValueError("Trait labels missing from local pinned annotation metadata")
    comparisons = []
    families = []
    pooled = None
    for stratum_id, kind, label, group in strata(pairs):
        print(f"Comparing {stratum_id}: {len(group)} pairs", flush=True)
        comparison, records, stages = compare_stratum(
            raw, labels, group, baseline, stratum_id, kind, label
        )
        comparisons.append(comparison)
        families.extend(records)
        if stratum_id == "pooled":
            pooled = stages
    comparison = pd.concat(comparisons, ignore_index=True)
    write_table(comparison, "trait_before_after.tsv")
    changed = comparison.before_eligible.ne(comparison.after_eligible)
    for cutoff in [0.001, 0.05]:
        changed |= comparison[f"before_q_le_{cutoff}"].ne(comparison[f"after_q_le_{cutoff}"])
    write_table(comparison[changed], "changed_significance_membership.tsv")
    write_table(comparison[comparison.direction_changed], "changed_directions.tsv")
    write_table(comparison[comparison.trait.isin(SIX_IDS)], "table1_six_traits_before_after.tsv")
    family_frame = pd.DataFrame(families)
    write_table(family_frame, "correction_families.tsv")
    scope, noel_corpus = denominator_audit(raw, pairs, metadata)
    overlap = integration(comparison, raw, pairs)
    refresh_caches(raw, metadata, labels, pairs, baseline)
    if pooled is None:
        raise ValueError("Pooled statistics missing")
    for stage in ["before", "after"]:
        scores = corrected.add_pair_scores(pooled[stage][2], pairs)
        conditions = corrected.summarize_groups(scores, "condition", min_pairs=10)
        sites = corrected.summarize_groups(scores, "body_site_group", min_pairs=3)
        write_table(conditions, f"{stage}/condition_summary.tsv")
        write_table(sites, f"{stage}/body_site_summary.tsv")
        renderer = baseline if stage == "before" else corrected
        renderer.OUTPUT_PNG = OUTPUT / stage / "paired_contrast_infographic.png"
        renderer.create_infographic(
            pooled[stage][0].reset_index(),
            conditions,
            sites,
            pair_count=len(pairs),
            study_count=pairs.study_id.nunique(),
            trait_count=len(pooled[stage][0]),
        )
    sources = [
        matrix_path,
        meta_root / "study_trait_profile_paired_enrichment_metadata.tsv",
        meta_root / "study_trait_profile_study_metadata.tsv",
        labels_path,
        OUTPUT / "inputs/download_provenance.json",
        ROOT / "data/kg-microbe/provenance.json",
        WORK / "baseline/scripts/create_trait_profile_infographic.py",
        WORK / "baseline/scripts/create_study_trait_profile_heatmap.py",
        WORK / "ora/scripts/trait_enrichment.py",
        WORK / "ora-q05.tsv",
    ]
    provenance = {
        "baseline_code_commit": BASELINE_COMMIT,
        "implementation_base_commit": (WORK / "code_commit.txt").read_text().strip(),
        "implementation_identification": "base commit plus SHA256 of implementation files below",
        "ora_code_commit": (WORK / "ora_commit.txt").read_text().strip(),
        "input_git_blob": INPUT_BLOB,
        "denominator": DENOMINATOR_DEFINITION,
        "inputs": {p.relative_to(ROOT).as_posix(): sha256(p) for p in sources},
        "original_full_dump_sha256": (
            "af84bfca7bfe76913a20c1fedbaadbfad4e4bf78605a2cfd7d00eb0c0fdd5788"
        ),
        "original_full_dump_usage": "recorded source lineage; frozen metadata used in rerun",
        "ora_test_count": int(
            retained_trait_counts(raw.loc[corrected.profile_ids_for_pairs(pairs)]).gt(0).sum().sum()
        ),
        "ora_nonzero_profiles": int(
            scope.loc[1, "profiles"] - scope.loc[1, "zero_retained_signal"]
        ),
        "implementation_sha256": {
            p.relative_to(ROOT).as_posix(): sha256(p)
            for p in [
                ROOT / "src/bsdbng/trait_composition.py",
                ROOT / "src/bsdbng/trait_guilds.py",
                ROOT / "src/bsdbng/table1_trait_maps.py",
                Path(__file__),
                ROOT / "scripts/create_trait_profile_infographic.py",
                ROOT / "scripts/create_study_trait_profile_heatmap.py",
                ROOT / "uv.lock",
                ROOT / "Justfile",
            ]
        },
        "runtime": {
            pkg: importlib.metadata.version(pkg)
            for pkg in ["numpy", "pandas", "scipy", "matplotlib"]
        },
        "pairing": (
            "frozen human-host study/experiment inc-dec pairs; no study reweighting; "
            "original study-level display labels"
        ),
        "zero_signal_policy": "retain rows; zero fractions/log1p inputs; standardize with all rows",
        "noel_matched_profiles": noel_corpus.kept_profiles,
        "baseline_validation": (
            "all emitted stratum eligibility, effects, p and q match published tables "
            "at rtol=1e-9, atol=1e-12"
        ),
        "ora_validation": "q<=0.001 membership and q match frozen global-BH/fold>=1.5 output",
    }
    (OUTPUT / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    publish_current_outputs(family_frame, pairs)
    print(scope.to_string(index=False))
    print(family_frame.to_string(index=False))
    print(overlap.to_string(index=False))
    print(f"Direction changes: {comparison.direction_changed.sum()}")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
