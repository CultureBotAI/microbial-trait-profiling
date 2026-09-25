import importlib.util
import sys
from pathlib import Path

import pandas as pd

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts/create_pd_exploratory_taxa_influence_heatmap.py"
)
sys.path.insert(0, str(MODULE_PATH.parent))
SPEC = importlib.util.spec_from_file_location("pd_exploratory_taxa_heatmap", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
heatmap = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(heatmap)


def test_exploratory_row_selection_prioritizes_literature_matches() -> None:
    rows = pd.DataFrame(
        [
            ["taxon-lit-family", "Known family", "family", 1, 1, 2, 2, 0.1],
            ["taxon-lit-species", "Known species", "species", 1, 0, 1, 1, 0.2],
            ["taxon-top-genus", "Top genus", "genus", 8, 3, 5, 1, 0.3],
            ["taxon-other", "Other species", "species", 5, 0, 4, 1, 0.4],
        ],
        columns=[
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
            "increased_profiles",
            "decreased_profiles",
            "paired_support_studies",
            "annotated_trait_count",
            "max_absolute_influence_z",
        ],
    )
    reference_table = pd.DataFrame(
        [
            [
                "taxon-lit-family",
                "Known family",
                "decreased",
                "Known family",
                0,
                "similar",
                "agrees",
            ],
            [
                "taxon-lit-species",
                "Known species",
                "increased",
                "Known species",
                0,
                "similar",
                "agrees",
            ],
            [
                "taxon-top-genus",
                "not listed",
                "not listed",
                "",
                float("nan"),
                "not applicable",
                "not in literature reference",
            ],
            [
                "taxon-other",
                "not listed",
                "not listed",
                "",
                float("nan"),
                "not applicable",
                "not in literature reference",
            ],
        ],
        columns=[
            "taxon_id",
            "literature_taxon",
            "literature_direction",
            "literature_match_term",
            "direction_hamming_distance",
            "hamming_interpretation",
            "direction_agreement",
        ],
    )

    selected = heatmap.select_exploratory_rows(rows, reference_table, max_rows=3)

    assert selected["taxon_id"].tolist() == [
        "taxon-lit-family",
        "taxon-lit-species",
        "taxon-top-genus",
    ]
    assert selected["selection_reason"].tolist() == [
        "PD literature seed",
        "PD literature seed",
        "top supported exploratory roll-up",
    ]


def test_rows_are_grouped_into_species_genus_family_blocks() -> None:
    rows = pd.DataFrame(
        [
            ["family-a", "family"],
            ["species-a", "species"],
            ["genus-a", "genus"],
            ["species-b", "species"],
            ["family-b", "family"],
        ],
        columns=["taxon_id", "taxonomic_rank"],
    )

    grouped = heatmap.group_rows_by_taxonomic_rank(rows)

    assert grouped["taxon_id"].tolist() == [
        "species-a",
        "species-b",
        "genus-a",
        "family-a",
        "family-b",
    ]
