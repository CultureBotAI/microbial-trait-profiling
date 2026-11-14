"""Verify that methods/table1_trait_condition_body_site_maps.md distinguishes the
canonical human-only profile matrix (0 non-human rows dropped) from the earlier
all-host wide matrices where 1,500 or 1,482 non-human rows were dropped (#333)."""

from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
METHODS_PATH = REPO_ROOT / "methods/table1_trait_condition_body_site_maps.md"


def _input_corpus_section(text: str) -> str:
    after_heading = text.split("## Input corpus\n", 1)[1]
    return after_heading.split("\n## ", 1)[0]


def test_table1_maps_methods_distinguishes_human_only_and_historical_all_host_matrices() -> None:
    """#333: the Input corpus section must bind the 0-row and 1,500-row host-filter drops
    to their respective matrices rather than leaving contradictory statements side by side."""
    text = METHODS_PATH.read_text(encoding="utf-8")
    section = " ".join(_input_corpus_section(text).split())

    # Canonical human-only count_input_16s_wms_v1 matrix: 5,629 rows, 0 non-human rows dropped.
    assert "5,629 profile rows read under the `count_input_16s_wms_v1` profile policy" in section
    assert (
        "none was dropped as non-human, because the matrix holds human profiles only (#656)"
        in section
    )
    assert (
        "The host filter drops 0 rows on the canonical `count_input_16s_wms_v1` profile "
        "matrix identified below (`519a095b…` / `43b2e230…`, 5,629 profile rows)"
    ) in section

    # Historical all-host wide matrices: 7,329 rows -> 1,500 dropped (#333) and
    # 7,111 rows -> 1,482 dropped (#222).
    assert (
        "`6e70724e…` over 7,329 profile rows, where 1,500 rows were dropped as non-human (#333)"
        in section
    )
    assert (
        "`a0c28978…` over 7,111 profile rows, where 1,482 rows were dropped as non-human (#222)"
        in section
    )

    # Guard against the old unqualified parenthetical or an unattributed 1,500-row drop.
    assert (
        "the v1.3.1 wide matrix is already human-only upstream, so that filter drops 0 rows"
        not in section
    )
    sentences = [s.strip() for s in re.split(r"(?<=[.!?])\s+", section) if s.strip()]
    drop_1500_sentences = [s for s in sentences if "1,500" in s]
    assert len(drop_1500_sentences) == 1
    assert "7,329" in drop_1500_sentences[0]
    assert "6e70724e" in drop_1500_sentences[0]


def test_table1_maps_methods_names_canonical_producer_module() -> None:
    """#2193: producer references must name src/mtp/table1_trait_maps.py, not the shim."""
    text = METHODS_PATH.read_text(encoding="utf-8")
    normalized = " ".join(text.split())
    assert (REPO_ROOT / "src/mtp/table1_trait_maps.py").is_file()
    assert (
        "Produced by `src/mtp/table1_trait_maps.py`, driven by `just table1-trait-maps`"
        in normalized
    )
    assert (
        "are produced by `src/mtp/table1_trait_maps.py`, driven deterministically by "
        "`just table1-trait-maps`" in normalized
    )
    assert "src/bsdbng/table1_trait_maps.py" not in text
    assert "src/bsdbng/guild_maps.py" not in text
