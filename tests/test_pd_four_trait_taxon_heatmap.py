import importlib.util
import math
from pathlib import Path

import pandas as pd

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts/create_pd_four_trait_taxon_heatmap.py"
SPEC = importlib.util.spec_from_file_location("pd_four_trait_taxon_heatmap", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
heatmap = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(heatmap)


def test_species_counts_keep_zero_denominators_and_missing_annotations_distinct() -> None:
    source = pd.DataFrame(
        [
            ["s1", "Species one", "species", "coccus shaped", 1, 5, 0],
            ["s1", "Species one", "species", "facultatively aerobic", 1, 5, 0],
            ["s2", "Species two", "species", "coccus shaped", 1, 1, 2],
            ["s3", "Decreased only", "species", "coccus shaped", 1, 0, 4],
            ["g1", "Genus one", "genus", "coccus shaped", 10, 12, 1],
            ["s4", "Other trait", "species", "sphere shaped", 1, 2, 1],
        ],
        columns=[
            "taxon_id",
            "taxon_name",
            "taxonomic_rank",
            "trait",
            "annotation_count",
            "increased_profiles",
            "decreased_profiles",
        ],
    )
    cells, rows = heatmap.prepare_cells(source)
    assert rows["taxon_id"].tolist() == ["s1", "s2"]
    assert len(cells) == 8
    indexed = cells.set_index(["taxon_id", "trait"])
    zero = indexed.loc[("s1", "coccus shaped")]
    assert zero["inc_dec_counts"] == "5/0"
    assert math.isinf(zero["inc_dec_ratio"])
    assert zero["color_log2_ratio"] == heatmap.COLOR_LIMIT
    decreased = indexed.loc[("s2", "coccus shaped")]
    assert decreased["inc_dec_counts"] == "1/2"
    assert decreased["inc_dec_ratio"] == 0.5
    assert decreased["color_log2_ratio"] == -1
    missing = indexed.loc[("s1", "arsenate reduction")]
    assert not missing["trait_annotated"]
    assert pd.isna(missing["inc_dec_ratio"])
    assert missing["inc_dec_counts"] == ""
