"""Pin active producer links in the six zero-overlap contributor docs (Refs #2026)."""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CANONICAL_REFERENCES: dict[str, dict[str, int]] = {
    ".claude/skills/review-open-issues/SKILL.md": {
        "src/mtp/download.py": 1,
        "src/mtp/ingest.py": 2,
        "src/mtp/study_tsv.py": 1,
        "src/mtp/taxa_trait_matrix.py": 1,
        "src/mtp/study_trait_profiles.py": 1,
        "src/mtp/trait_guilds.py": 1,
        "src/mtp/trait_composition.py": 1,
        "src/mtp/stats.py": 1,
        "src/mtp/guild_maps.py": 1,
        "src/mtp/table1_trait_maps.py": 1,
        "src/mtp/trait_screen.py": 1,
        "src/mtp/enrichment_stats.py": 1,
        "src/mtp/calibration.py": 1,
        "src/mtp/corpus_identity.py": 1,
    },
    "methods/README.md": {
        "src/mtp/guild_maps.py": 1,
        "src/mtp/table1_trait_maps.py": 3,
        "src/mtp/trait_screen.py": 1,
    },
    "methods/integrated_evidence.md": {"src/mtp/integrated_evidence.py": 1},
    "methods/table1_log_ratio_sensitivity.md": {"src/mtp/log_ratio_sensitivity.py": 1},
    "methods/table1_paired_log_ratios.md": {"src/mtp/trait_composition.py": 1},
    "docs/rosalind-six-trait-interpretation.md": {
        "src/mtp/study_trait_profiles.py": 1,
        "src/mtp/trait_guilds.py": 2,
    },
}


@pytest.mark.parametrize("document", CANONICAL_REFERENCES)
def test_active_producer_references_use_canonical_mtp(document: str) -> None:
    """Repeated references must all migrate, and every destination must exist."""
    text = (ROOT / document).read_text(encoding="utf-8")
    for canonical, count in CANONICAL_REFERENCES[document].items():
        assert text.count(canonical) == count, (
            f"canonical producer reference count: {document}: {canonical}: "
            f"expected {count}, found {text.count(canonical)}"
        )
        assert (ROOT / canonical).is_file(), f"Missing canonical producer: {canonical}"
    assert "src/bsdbng/" not in text, f"Active producer still points to a shim: {document}"
