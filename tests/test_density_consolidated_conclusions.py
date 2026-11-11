"""Consolidated section 6 describes computed diagnostics without fixed claims (#1328).

Synthetic report inputs cover correlations and dimensionality that contradict the
historical narrative. These are rendering tests; no audit pipeline or data is read.
"""

from __future__ import annotations

import importlib.util
import re
from dataclasses import replace
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from mtp.density_adjustment import (
    ConditionalTwinSummary,
    CorrelationMetrics,
    Table1TraitEffect,
    adjustment_status,
)
from mtp.trait_lists import DeclaredTraitList

_SECTION = "6. CONCLUSION & METHODOLOGICAL IMPLICATIONS\n"


@pytest.fixture(scope="module")
def audit_script() -> ModuleType:
    path = Path(__file__).resolve().parents[1] / "scripts/audit_trait_density_correlation.py"
    spec = importlib.util.spec_from_file_location("density_conclusions_under_test", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _metrics(median: float, pc1: float) -> CorrelationMetrics:
    return CorrelationMetrics(
        n_samples=100,
        n_traits=40,
        constant_traits_excluded=0,
        median_pairwise_r=median,
        mean_pairwise_r=median,
        pct_pairs_gt_05=0.0,
        pct_pairs_gt_09=0.0,
        pc1_variance_pct=pc1,
        pc1_5_variance_pct=70.0,
        m_eff_li_ji=30.0,
        m_eff_galwey=28.0,
        m_eff_nyholt=32.0,
        test_pair_r=-0.23,
        near_twin_r=-0.8071,
    )


def _effect(
    trait_id: str, label: str, raw_q: float, adj_q: float, adj_mean: float = 0.01
) -> Table1TraitEffect:
    return Table1TraitEffect(
        trait_id=trait_id,
        trait_label=label,
        raw_mean_diff=0.05,
        raw_wilcoxon_p=0.001,
        adj_mean_diff=adj_mean,
        adj_wilcoxon_p=0.01,
        relative_change_pct=-80.0,
        status=adjustment_status(0.05, raw_q, adj_mean, adj_q),
        magnitude_change_pct=-80.0,
        sign_flipped=adj_mean < 0,
        raw_wilcoxon_bh_q=raw_q,
        adj_wilcoxon_bh_q=adj_q,
        bh_family_size=4,
    )


def synthetic_results(case: str = "negative_correlation") -> tuple[dict[str, Any], dict[str, Any]]:
    """Reusable inputs for report tests and the independent before/after proof."""
    density: dict[str, Any] = {
        "raw_pearson": _metrics(-0.0234, 7.12),
        "raw_spearman": _metrics(-0.0312, 6.54),
        "adj_pc1_pearson": _metrics(-0.0456, 6.13),
        "table1_effects": [
            _effect("METPO:1000608", "facultative union (fixture)", 0.0496, 0.05004),
            _effect("T:survivor", "survivor fixture", 0.04, 0.05),
            _effect("T:adjusted_only", "adjusted-only fixture", 0.05004, 0.0496),
            _effect("T:reversed", "reversed fixture", 0.01, 0.02, -0.01),
        ],
        "twin_summary": ConditionalTwinSummary(
            trait_a="METPO:1000668",
            trait_b="METPO:1000683",
            trait_a_label="coccus shaped",
            trait_b_label="sphere shaped",
            profile_correlation=-0.8071,
            profile_partial_r_depth=-0.4,
            vif_profile=2.87,
            profile_slope_a_on_b=-0.8,
            profile_intercept_a_on_b=0.01,
            n_paired_experiments=50,
            paired_correlation=-0.3,
            vif_paired=1.10,
            paired_slope_a_on_b=-0.3,
            paired_slope_pvalue=0.1,
        ),
    }
    ontology: dict[str, Any] = {
        "n_redundant_pairs": 8,
        "n_clusters": 3,
        "n_collapsed_features": 6,
        "bounded_trait_space": 34,
        "native_metrics": _metrics(-0.0123, 4.23),
        "max_metrics": _metrics(-0.0098, 3.34),
        "adj_collapsed_metrics": _metrics(-0.0045, 2.45),
    }
    if case == "low_vif":
        density["twin_summary"] = replace(
            density["twin_summary"],
            profile_correlation=0.0,
            vif_profile=1.0,
            paired_correlation=0.2,
            vif_paired=1.04,
        )
    elif case == "undefined":
        density["raw_pearson"] = _metrics(float("nan"), float("nan"))
        density["adj_pc1_pearson"] = _metrics(float("nan"), float("nan"))
        density["twin_summary"] = replace(
            density["twin_summary"],
            profile_correlation=float("nan"),
            vif_profile=float("nan"),
            paired_correlation=-1.0,
            vif_paired=float("inf"),
        )
        for key in ("native_metrics", "max_metrics", "adj_collapsed_metrics"):
            ontology[key] = _metrics(float("nan"), float("nan"))
    elif case != "negative_correlation":
        raise ValueError(f"Unknown synthetic case: {case}")
    return density, ontology


def _report_and_conclusions(
    audit_script: ModuleType,
    output_dir: Path,
    density: dict[str, Any],
    ontology: dict[str, Any],
    trait_list: DeclaredTraitList | None = None,
) -> tuple[str, str]:
    report: Path = audit_script.write_consolidated_report(output_dir, density, ontology, trait_list)
    assert report.is_file()
    content = report.read_text(encoding="utf-8")
    assert content.count(_SECTION) == 1
    if trait_list is None:
        _assert_no_fixed_claims(content)
    return content, content.split(_SECTION, 1)[1]


def _conclusions(
    audit_script: ModuleType,
    output_dir: Path,
    density: dict[str, Any],
    ontology: dict[str, Any],
    trait_list: DeclaredTraitList | None = None,
) -> str:
    _, conclusions = _report_and_conclusions(
        audit_script, output_dir, density, ontology, trait_list
    )
    return conclusions


def _point(conclusions: str, number: int) -> str:
    match = re.search(rf"(?ms)^{number}\. (.*?)(?=^\d+\. |^={{5,}}|\Z)", conclusions)
    assert match is not None, f"Missing conclusion {number}: {conclusions}"
    return " ".join(match.group(1).split())


def _assert_no_fixed_claims(conclusions: str) -> None:
    text = " ".join(conclusions.lower().split())
    for claim in (
        "vast majority",
        "driven by curation density",
        "driven by pc1",
        "separating biology from curation effort",
        "confirms pc1 is curation density",
        "severe collinearity",
        "indicating strong collinearity",
        "conditional testing confirms",
        "unstable coefficients",
        "collapsing alone does not collapse pc1",
        "does not attenuate",
        "necessary and complementary",
        "strictly necessary",
    ):
        assert claim not in text
    # A VIF and a pairwise regression do not evaluate a joint outcome model.
    assert re.search(r"(?:stability|joint.{0,25}model).{0,80}not (?:evaluated|assessed)", text)


@pytest.mark.parametrize(
    ("profile_r", "profile_vif", "paired_r", "paired_vif"),
    [
        (-0.8071, 2.87, -0.3, 1.10),
        (0.0, 1.0, 0.2, 1.04),
        (0.9594, 12.56, -0.8660, 4.0),
    ],
)
def test_conclusions_report_profile_and_paired_diagnostics_without_severity_claims(
    tmp_path: Path,
    audit_script: ModuleType,
    profile_r: float,
    profile_vif: float,
    paired_r: float,
    paired_vif: float,
) -> None:
    density, ontology = synthetic_results()
    density["twin_summary"] = replace(
        density["twin_summary"],
        profile_correlation=profile_r,
        vif_profile=profile_vif,
        paired_correlation=paired_r,
        vif_paired=paired_vif,
    )
    conclusions = _conclusions(audit_script, tmp_path, density, ontology)
    twin = _point(conclusions, 4)
    for label, correlation, vif in (
        ("profile", profile_r, profile_vif),
        ("paired", paired_r, paired_vif),
    ):
        r_text = re.escape(f"{correlation:+.4f}")
        vif_text = re.escape(f"{vif:.2f}")
        assert re.search(rf"{label}.{{0,60}}r\s*=\s*{r_text}", twin, re.I)
        assert re.search(rf"{label}.{{0,100}}VIF\s*=\s*{vif_text}", twin, re.I)
    _assert_no_fixed_claims(conclusions)


@pytest.mark.parametrize(
    ("raw_pc1", "native_pc1", "max_pc1", "adjusted_pc1"),
    [
        (7.12, 4.23, 3.34, 2.45),
        (82.12, 43.23, 14.34, 9.45),
        (7.12, 40.23, 51.34, 62.45),
    ],
)
def test_conclusions_report_observed_correlations_and_pc1_changes(
    tmp_path: Path,
    audit_script: ModuleType,
    raw_pc1: float,
    native_pc1: float,
    max_pc1: float,
    adjusted_pc1: float,
) -> None:
    density, ontology = synthetic_results()
    density["raw_pearson"] = replace(density["raw_pearson"], pc1_variance_pct=raw_pc1)
    for key, value in (
        ("native_metrics", native_pc1),
        ("max_metrics", max_pc1),
        ("adj_collapsed_metrics", adjusted_pc1),
    ):
        ontology[key] = replace(ontology[key], pc1_variance_pct=value)
    content, conclusions = _report_and_conclusions(audit_script, tmp_path, density, ontology)
    assert (f"collapsed + PC1 regression PC1 variance explained = {adjusted_pc1:.2f}%).") in content
    assert "collapsed + PC1 regression variance explained" not in content
    correlation = _point(conclusions, 1)
    assert "raw = -0.0234; post-PC1 regression = -0.0456." in correlation
    assert f"Raw PC1 variance explained (%) = {raw_pc1:.2f}." in correlation
    collapse = _point(conclusions, 5)
    assert "8 redundant pairs (6 features collapsed into 3 clusters)" in collapse
    assert "bounding distinct trait space to 34 features" in collapse
    assert (
        f"PC1 variance explained (%): raw = {raw_pc1:.2f}; "
        f"collapsed native = {native_pc1:.2f}; "
        f"collapsed max-member = {max_pc1:.2f}; "
        f"collapsed + PC1 regression = {adjusted_pc1:.2f}."
    ) in collapse
    assert (
        "Median pairwise Pearson r: collapsed native = -0.0123; "
        "collapsed max-member = -0.0098; "
        "collapsed + PC1 regression = -0.0045."
    ) in collapse
    _assert_no_fixed_claims(conclusions)


@pytest.mark.parametrize("unavailable", [float("nan"), float("-inf")])
def test_conclusions_distinguish_undefined_from_infinite_diagnostics(
    tmp_path: Path, audit_script: ModuleType, unavailable: float
) -> None:
    density, ontology = synthetic_results("undefined")
    density["raw_pearson"] = replace(_metrics(unavailable, unavailable), test_pair_r=unavailable)
    density["raw_spearman"] = _metrics(unavailable, unavailable)
    density["adj_pc1_pearson"] = replace(
        _metrics(unavailable, unavailable), test_pair_r=unavailable
    )
    density["twin_summary"] = replace(
        density["twin_summary"],
        profile_correlation=unavailable,
        vif_profile=unavailable,
        profile_slope_a_on_b=unavailable,
        paired_correlation=-1.0,
        vif_paired=float("inf"),
    )
    for key in ("native_metrics", "max_metrics", "adj_collapsed_metrics"):
        ontology[key] = _metrics(unavailable, unavailable)
    content, conclusions = _report_and_conclusions(audit_script, tmp_path, density, ontology)
    twin = _point(conclusions, 4)
    assert re.search(r"profile.{0,60}r\s*=\s*NA", twin, re.I)
    assert re.search(r"profile.{0,100}VIF\s*=\s*NA", twin, re.I)
    assert re.search(r"paired.{0,60}r\s*=\s*-1\.0000", twin, re.I)
    assert re.search(r"paired.{0,100}VIF\s*=\s*inf", twin, re.I)
    assert "Variance Inflation Factor (VIF): NA" in content
    assert "Paired VIF: inf" in content
    assert "Post-PC1 regression r: NA (change: NA)" in content
    assert "Profile OLS regression slope (coccus ~ sphere): NA" in content
    assert "NA" in _point(conclusions, 1)
    assert "NA" in _point(conclusions, 5)
    assert not re.search(r"\bnan\b", content, re.I)
    _assert_no_fixed_claims(content)


def test_conclusions_report_no_paired_experiments_as_unavailable(
    tmp_path: Path, audit_script: ModuleType
) -> None:
    density, ontology = synthetic_results()
    density["twin_summary"] = replace(
        density["twin_summary"],
        n_paired_experiments=0,
        paired_correlation=float("nan"),
        vif_paired=float("nan"),
    )
    content, conclusions = _report_and_conclusions(audit_script, tmp_path, density, ontology)
    twin = _point(conclusions, 4)
    assert re.search(r"paired.{0,60}N=0", twin, re.I)
    assert re.search(r"paired.{0,60}r\s*=\s*NA", twin, re.I)
    assert re.search(r"paired.{0,100}VIF\s*=\s*NA", twin, re.I)
    assert not re.search(r"\bnan\b", content, re.I)
    _assert_no_fixed_claims(content)


@pytest.mark.parametrize("n_constant", [38, 39, 40])
def test_conclusions_do_not_treat_unavailable_zero_sentinels_as_observations(
    tmp_path: Path, audit_script: ModuleType, n_constant: int
) -> None:
    """Two or more nonconstant traits report computed values; <2 report NA (#1406, #1514)."""
    density, ontology = synthetic_results()
    base_metrics = _metrics(0.1234, 45.67) if n_constant == 38 else _metrics(0.0, 0.0)
    sentinel = replace(base_metrics, constant_traits_excluded=n_constant)
    density["raw_pearson"] = sentinel
    density["raw_spearman"] = sentinel
    density["adj_pc1_pearson"] = sentinel
    for key in ("native_metrics", "max_metrics", "adj_collapsed_metrics"):
        ontology[key] = sentinel
    content, conclusions = _report_and_conclusions(audit_script, tmp_path, density, ontology)
    correlation = _point(conclusions, 1)
    collapse = _point(conclusions, 5)
    if n_constant == 38:
        assert "raw = +0.1234; post-PC1 regression = +0.1234." in correlation
        assert "Raw PC1 variance explained (%) = 45.67." in correlation
        assert (
            "PC1 variance explained (%): raw = 45.67; collapsed native = 45.67; "
            "collapsed max-member = 45.67; collapsed + PC1 regression = 45.67."
        ) in collapse
        assert (
            "Median pairwise Pearson r: collapsed native = +0.1234; "
            "collapsed max-member = +0.1234; collapsed + PC1 regression = +0.1234."
        ) in collapse
        assert "= NA" not in correlation + collapse
        return
    assert re.search(r"raw\s*=\s*NA", correlation)
    assert re.search(r"post-PC1 regression\s*=\s*NA", correlation)
    assert re.search(r"PC1 variance explained \(%\)\s*=\s*NA", correlation)
    for label in (
        "raw",
        "collapsed native",
        "collapsed max-member",
        "collapsed + PC1 regression",
    ):
        assert re.search(rf"{re.escape(label)}\s*=\s*NA", collapse)
    assert "0.00" not in correlation + collapse
    assert "Pearson pairwise correlation: median NA, mean NA" in content
    assert "PC1 variance explained: NA | PC1-PC5 variance explained: NA" in content
    assert "Median pairwise rho: NA, Mean pairwise rho: NA" in content
    assert "Baseline Pearson r: NA" in content
    assert "Post-PC1 regression r: NA (change: NA)" in content
    assert "PC1 var = NA | Med r = NA" in content


def test_conclusions_keep_bh_statuses_labels_and_survivor_membership(
    tmp_path: Path, audit_script: ModuleType
) -> None:
    density, ontology = synthetic_results()
    conclusions = _conclusions(audit_script, tmp_path, density, ontology)
    contrast = _point(conclusions, 2)
    assert "facultative union (fixture)" in contrast
    assert "'collapsed (attenuated to nonsignificant)'" in contrast
    assert "raw BH q=4.96e-02, adjusted BH q=5.004e-02" in contrast
    survivors = _point(conclusions, 3)
    assert "BH q <= 0.05" in survivors
    assert "survivor fixture" in survivors
    assert "adjusted-only fixture" not in survivors
    assert "reversed fixture" not in survivors
    assert "facultative union (fixture)" not in survivors


def test_declared_trait_list_keeps_the_existing_section6_guard(
    tmp_path: Path, audit_script: ModuleType
) -> None:
    density, ontology = synthetic_results()
    declared = DeclaredTraitList(
        tuple(effect.trait_id for effect in density["table1_effects"]), "fixture.txt", "a" * 64
    )
    conclusions = _conclusions(audit_script, tmp_path, density, ontology, declared)
    assert conclusions == (
        "-" * 80
        + "\nNot written for a declared trait list; section 3 reports each listed\n"
        + "trait's computed effects.\n"
        + "=" * 80
    )


def test_consolidated_report_pluralizes_singular_and_plural_counts(
    tmp_path: Path, audit_script: ModuleType
) -> None:
    """#1407: singular counts use singular nouns ('1 redundant pair', '1 feature', '1 cluster')."""
    density, ontology = synthetic_results()
    ontology = {
        **ontology,
        "n_redundant_pairs": 1,
        "n_clusters": 1,
        "n_collapsed_features": 1,
        "bounded_trait_space": 1,
    }
    content, conclusions = _report_and_conclusions(audit_script, tmp_path, density, ontology)
    assert "Redundant pairs identified: 1 pair across 1 concept cluster" in content
    assert "Bounded distinct trait space: 1 feature" in content
    collapse = _point(conclusions, 5)
    assert "1 redundant pair (1 feature collapsed into 1 cluster)" in collapse
    assert "bounding distinct trait space to 1 feature." in collapse

    # #1826, #1857: when a collapsed member (e.g. sphere) is present without its cluster
    # survivor (coccus), n_redundant_pairs == 0 and n_clusters == 0 while
    # n_collapsed_features == 1; item 5 must state that 1 feature was removed because its
    # cluster survivor is absent, and must count orphan removals separately even when
    # n_clusters > 0.
    orphan_ontology = {
        **ontology,
        "n_redundant_pairs": 0,
        "n_clusters": 0,
        "n_collapsed_features": 1,
        "n_orphan_features": 1,
        "bounded_trait_space": 2,
    }
    orphan_content, orphan_conclusions = _report_and_conclusions(
        audit_script, tmp_path, density, orphan_ontology
    )
    orphan_collapse = _point(orphan_conclusions, 5)
    assert "collapsed into 0 clusters" not in orphan_content
    assert "across 0 clusters" not in orphan_collapse
    assert (
        "0 redundant pairs (1 feature removed because its cluster survivor is absent)"
        in orphan_collapse
    )

    multi_orphan_ontology = {
        **ontology,
        "n_redundant_pairs": 0,
        "n_clusters": 0,
        "n_collapsed_features": 2,
        "n_orphan_features": 2,
        "bounded_trait_space": 2,
    }
    _, multi_orphan_conclusions = _report_and_conclusions(
        audit_script, tmp_path, density, multi_orphan_ontology
    )
    assert (
        "0 redundant pairs (2 features removed because their cluster survivors are absent)"
        in _point(multi_orphan_conclusions, 5)
    )

    mixed_ontology = {
        **ontology,
        "n_redundant_pairs": 1,
        "n_clusters": 1,
        "n_collapsed_features": 2,
        "n_orphan_features": 1,
        "bounded_trait_space": 2,
    }
    _, mixed_conclusions = _report_and_conclusions(audit_script, tmp_path, density, mixed_ontology)
    mixed_collapse = _point(mixed_conclusions, 5)
    assert "2 features collapsed into 1 cluster" not in mixed_collapse
    assert (
        "1 redundant pair (1 feature collapsed into 1 cluster; "
        "1 feature removed because its cluster survivor is absent)"
    ) in mixed_collapse

    # When n_orphan_features is recorded on a multi-pair corpus where n_pairs > n_collapsed
    # (so max(0, n_collapsed - n_pairs) would give 0), the recorded orphan count is used.
    motility_orphan_ontology = {
        **ontology,
        "n_redundant_pairs": 21,
        "n_clusters": 9,
        "n_collapsed_features": 16,
        "n_orphan_features": 3,
        "bounded_trait_space": 301,
    }
    _, motility_orphan_conclusions = _report_and_conclusions(
        audit_script, tmp_path, density, motility_orphan_ontology
    )
    motility_orphan_collapse = _point(motility_orphan_conclusions, 5)
    assert "16 features collapsed into 9 clusters" not in motility_orphan_collapse
    assert (
        "21 redundant pairs (13 features collapsed into 9 clusters; "
        "3 features removed because their cluster survivors are absent)"
    ) in motility_orphan_collapse

    # When n_orphan_features is absent from ontology_res and n_clusters > 0, keep base's
    # combined wording rather than estimating orphans from n_collapsed - n_pairs; when
    # n_clusters == 0 (#1906), every present collapsed member is an orphan so item 5 must
    # not fall back to "collapsed into 0 clusters".
    legacy_ontology = {k: v for k, v in mixed_ontology.items() if k != "n_orphan_features"}
    _, legacy_conclusions = _report_and_conclusions(
        audit_script, tmp_path, density, legacy_ontology
    )
    legacy_collapse = _point(legacy_conclusions, 5)
    assert "removed because" not in legacy_collapse
    assert "1 redundant pair (2 features collapsed into 1 cluster)" in legacy_collapse

    legacy_zero_cluster_ontology = {
        k: v for k, v in orphan_ontology.items() if k != "n_orphan_features"
    }
    legacy_zero_content, legacy_zero_conclusions = _report_and_conclusions(
        audit_script, tmp_path, density, legacy_zero_cluster_ontology
    )
    legacy_zero_collapse = _point(legacy_zero_conclusions, 5)
    assert "collapsed into 0 clusters" not in legacy_zero_content
    assert (
        "0 redundant pairs (1 feature removed because its cluster survivor is absent)"
        in legacy_zero_collapse
    )

    legacy_plural_zero_cluster_ontology = {
        **legacy_zero_cluster_ontology,
        "n_collapsed_features": 3,
    }
    legacy_plural_zero_content, legacy_plural_zero_conclusions = _report_and_conclusions(
        audit_script, tmp_path, density, legacy_plural_zero_cluster_ontology
    )
    legacy_plural_zero_collapse = _point(legacy_plural_zero_conclusions, 5)
    assert "collapsed into 0 clusters" not in legacy_plural_zero_content
    assert "1 feature" not in legacy_plural_zero_collapse
    assert (
        "0 redundant pairs (3 features removed because their cluster survivors are absent)"
        in legacy_plural_zero_collapse
    )


def test_density_interpretation_and_audit_docs_qualify_historical_vif_and_list_provenance() -> None:
    """#1391, #1392, #1400, #1827, #1857: docs qualify historical binarized VIF 5.88, omit
    severity and curation-mechanism claims, state largely overlapping taxa, and list
    density_adjustment_provenance.json."""
    repo_root = Path(__file__).resolve().parents[1]
    bio_audit = (repo_root / "methods/biological_interpretation_audit.md").read_text(
        encoding="utf-8"
    )
    assert "show substantial collinearity" not in bio_audit
    g16_line = next(line for line in bio_audit.splitlines() if line.startswith("| G16 "))
    assert "historical binarized-profile VIF = 5.88" in g16_line
    assert "2.74" in g16_line

    coupling_doc = (repo_root / "methods/trait_coupling_sensitivity.md").read_text(encoding="utf-8")
    assert "historical binarized-profile r = +0.911" in coupling_doc
    assert "Historical binarized-profile VIF = 5.88" in coupling_doc
    assert "Curators assign" not in coupling_doc
    assert "annotated to the same taxa" not in coupling_doc
    assert "Coccus and sphere are annotated to largely overlapping taxa" in coupling_doc
    assert (
        "Both labels are annotated to largely overlapping taxa in the trait matrix." in coupling_doc
    )

    density_doc = (repo_root / "docs/density_adjustment_audit.md").read_text(encoding="utf-8")
    assert "data/derived/density_adjustment/density_adjustment_provenance.json" in density_doc
