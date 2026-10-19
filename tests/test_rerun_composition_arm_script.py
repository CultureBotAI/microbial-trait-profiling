"""The composition arm checks the pinned KG before step 1, not at step 8 (#695)."""

from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "rerun_composition_arm_jed.sh"


def test_prerequisites_include_the_pinned_kg_and_its_md5_check() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    prerequisites, _, steps = text.partition('echo "=== [1/9]')

    loop = prerequisites.split("for required_file in", 1)[1].split("; do", 1)[0]
    assert '"data/kg-microbe/merged-kg_edges.tsv"' in loop
    assert '"data/kg-microbe/merged-kg_nodes.tsv"' in loop
    assert "\njust verify-pinned-kg\n" in prerequisites
    assert "just pd-trait-contributors" in steps
    assert "methods/table1_taxon_influence_protocol.json" in prerequisites
