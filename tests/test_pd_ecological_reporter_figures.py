import importlib.util
from pathlib import Path

import pandas as pd
import pytest

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts/create_pd_ecological_reporter_figures.py"
)
SPEC = importlib.util.spec_from_file_location("pd_ecological_reporter_figures", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
figures = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(figures)


def _traits() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "trait_id": ["fac", "cell", "coccus", "sphere", "ars", "halo"],
            "trait_label": [
                "facultatively aerobic",
                "cellulose catabolic process",
                "coccus shaped",
                "sphere shaped",
                "arsenate reduction",
                "moderately halophilic",
            ],
        }
    )


def test_condition_summary_uses_paired_study_contrasts_and_oriented_composite() -> None:
    metadata = pd.DataFrame(
        [
            ["s1", "e1", "inc", "Parkinson's disease", "Neurologic", "Feces"],
            ["s1", "e1", "dec", "Parkinson's disease", "Neurologic", "Feces"],
            ["s2", "e2", "inc", "Parkinson's disease", "Neurologic", "Feces"],
            ["s2", "e2", "dec", "Parkinson's disease", "Neurologic", "Feces"],
            ["s3", "e3", "inc", "Parkinson's disease", "Neurologic", "Oral"],
        ],
        columns=[
            "study_id",
            "experiment_id",
            "inc_dec",
            "condition",
            "condition_group",
            "body_site_group",
        ],
        index=["p1", "p2", "p3", "p4", "p5"],
    )
    zscores = pd.DataFrame(
        [
            [2, 0, 0, 0, 1, 3],
            [0, 2, 0, 0, 0, 1],
            [1, 1, 0, 0, 3, 2],
            [0, 2, 0, 0, 1, 0],
            [100, 100, 100, 100, 100, 100],
        ],
        columns=["fac", "cell", "coccus", "sphere", "ars", "halo"],
        index=metadata.index,
    )
    summary, study = figures.build_condition_trait_summaries(metadata, zscores, _traits())
    pd_rows = study[study["condition"].eq("Parkinson's disease")]
    assert pd_rows["study_id"].nunique() == 2
    s1 = pd_rows[pd_rows["study_id"].eq("s1")].set_index("trait")
    assert s1.loc["cellulose catabolic process", "contrast"] == -2
    assert s1.loc[figures.POSITIVE_COMPOSITE, "contrast"] == pytest.approx(5 / 3)
    assert s1.loc[figures.REPORTER_COMPOSITE, "contrast"] == pytest.approx(11 / 3)
    reporter = summary[summary["trait"].eq(figures.REPORTER_COMPOSITE)].iloc[0]
    assert reporter["paired_study_count"] == 2
    assert reporter["increased_profile_count"] == 2
    assert reporter["decreased_profile_count"] == 2


def test_network_keeps_exact_genus_species_and_excludes_overlapping_family() -> None:
    metadata = pd.DataFrame(
        [
            ["s1", "e1", "inc", "PD", "Feces"],
            ["s2", "e2", "inc", "PD", "Feces"],
            ["s1", "e1", "dec", "PD", "Feces"],
            ["s2", "e2", "dec", "PD", "Feces"],
            ["s3", "e3", "inc", "Diet", "Feces"],
            ["s4", "e4", "inc", "Diet", "Feces"],
        ],
        columns=["study_id", "experiment_id", "inc_dec", "condition", "body_site_group"],
        index=["p1", "p2", "p3", "p4", "p5", "p6"],
    )
    metadata["profile_id"] = metadata.index
    taxa = pd.DataFrame(
        [
            ["s1", "e1", "increased", "g1"],
            ["s2", "e2", "increased", "g1"],
            ["s1", "e1", "decreased", "g1"],
            ["s2", "e2", "decreased", "g1"],
            ["s3", "e3", "increased", "sp1"],
            ["s4", "e4", "increased", "sp1"],
            ["s1", "e1", "increased", "f1"],
            ["s2", "e2", "increased", "f1"],
        ],
        columns=["study_id", "experiment_id", "direction", "taxon_id"],
    )
    matrix = pd.DataFrame(
        [
            ["g1", "Genus one", "genus", 1, 1, 0, 0, 0, 0],
            ["sp1", "Species one", "species", 1, 0, 0, 0, 0, 0],
            ["f1", "Family one", "family", 1, 1, 1, 1, 1, 1],
        ],
        columns=["taxon_id", "taxon_name", "taxonomic_rank", *(_traits()["trait_id"])],
    )
    nodes, edges, selection = figures.build_condition_taxon_trait_network(
        metadata,
        taxa,
        matrix,
        _traits(),
        ["PD", "Diet"],
        min_context_studies=2,
        top_taxa_per_context=1,
        max_taxa=4,
    )
    taxon_nodes = nodes[nodes["layer"].eq("taxon")]
    assert set(taxon_nodes["node_id"]) == {"taxon::g1", "taxon::sp1"}
    assert "taxon::f1" not in set(edges["target"])
    assert set(selection["taxonomic_rank"]) == {"genus", "species"}
    annotation_edges = edges[edges["evidence_type"].eq("KG-Microbe taxon-trait annotation")]
    assert set(annotation_edges["source"]) == {"taxon::g1", "taxon::sp1"}
