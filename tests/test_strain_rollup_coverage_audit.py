"""Account for the traits the BacDive strain rollup adds (#774)."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

from mtp.paths import REPO_ROOT

SCRIPT = REPO_ROOT / "scripts" / "adhoc" / "issue774_strain_rollup_coverage.py"


def write_inputs(root: Path) -> dict[str, Path]:
    paths = {
        "reported_taxa": root / "bugsigdb_taxa.tsv",
        "edges": root / "edges.tsv",
        "nodes": root / "nodes.tsv",
        "study_dir": root / "studies",
    }
    paths["reported_taxa"].write_text(
        "direction\ttaxon_id\ttaxon_name\ttaxonomic_rank\n"
        "increased\tNCBITaxon:100\tGenus\tgenus\n"
        "increased\tNCBITaxon:101\tSpecies a\tspecies\n"
        "decreased\tNCBITaxon:102\tSpecies b\tspecies\n",
        encoding="utf-8",
    )
    paths["edges"].write_text(
        "id\tsubject\tpredicate\tobject\n"
        "e1\tNCBITaxon:101\tbiolink:subclass_of\tNCBITaxon:100\n"
        "e2\tNCBITaxon:102\tbiolink:subclass_of\tNCBITaxon:100\n"
        "e3\tNCBITaxon:101\tbiolink:has_phenotype\tMETPO:1000603\n"
        "e4\tNCBITaxon:102\tbiolink:has_phenotype\tMETPO:1000605\n"
        "s1\tkgmicrobe.strain:1\tbiolink:subclass_of\tNCBITaxon:101\n"
        "s2\tkgmicrobe.strain:1\tbiolink:has_phenotype\tMETPO:1000605\n"
        "s3\tkgmicrobe.strain:2\tbiolink:subclass_of\tNCBITaxon:102\n"
        "s4\tkgmicrobe.strain:2\tMETPO:2000303\tEC:3.2.1.21\n"
        "s5\tkgmicrobe.strain:3\tbiolink:has_phenotype\tMETPO:1000608\n",
        encoding="utf-8",
    )
    paths["nodes"].write_text(
        "id\tcategory\tname\n"
        "METPO:1000603\tbiolink:PhenotypicQuality\tanaerobic\n"
        "METPO:1000605\tbiolink:PhenotypicQuality\tfacultatively anaerobic\n"
        "METPO:1000608\tbiolink:PhenotypicQuality\tfacultatively aerobic\n"
        "METPO:2000303\tbiolink:OntologyClass\tdoes not show activity of\n",
        encoding="utf-8",
    )
    paths["study_dir"].mkdir()
    (paths["study_dir"] / "bsdb_1.yaml").write_text(
        "id: bsdb:1\n"
        "experiments:\n"
        "- id: bsdb:1/1/1\n"
        "  host_species: Homo sapiens\n"
        "  sequencing_type: 16S\n"
        "  signatures:\n"
        "  - direction: increased\n"
        "    taxa:\n"
        "    - id: NCBITaxon:101\n"
        "  - direction: decreased\n"
        "    taxa:\n"
        "    - id: NCBITaxon:102\n",
        encoding="utf-8",
    )
    return paths


def test_audit_reports_what_the_strain_rollup_adds(tmp_path: Path) -> None:
    paths = write_inputs(tmp_path)
    output = tmp_path / "out"
    arguments = [f"--{name.replace('_', '-')}={path}" for name, path in paths.items()]
    result = subprocess.run(
        [sys.executable, str(SCRIPT), *arguments, f"--output-dir={output}"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert summary["strains"]["by_ncbitaxon_link_count"] == {"0": 1, "1": 2}
    assert summary["strains"]["strain_trait_pairs_without_ncbitaxon_link"] == 1
    assert summary["vocabulary"] == {
        "matrix_columns_off": 2,
        "matrix_columns_on": 3,
        "retained_off": 2,
        "retained_on": 3,
        "entering": 1,
        "leaving": 0,
        "entering_with_negated_predicate": 1,
        "entering_by_predicate": {"does not show activity of": 1},
    }
    anaerobic = summary["focal_traits"]["METPO:1000605"]
    assert (anaerobic["ncbitaxon_nodes_off"], anaerobic["ncbitaxon_nodes_on"]) == (1, 2)
    assert (anaerobic["reported_taxa_off"], anaerobic["reported_taxa_on"]) == (2, 3)
    human = summary["signature_taxa"]["cohorts"][1]
    assert human["host_species"] == "Homo sapiens"
    assert (human["signature_taxa"], human["gaining_a_trait_label"]) == (2, 2)
    assert human["focal_traits"]["METPO:1000605"] == {
        "label": "facultatively anaerobic",
        "taxa_off": 1,
        "taxa_on": 2,
    }

    with (output / "trait_coverage.tsv").open(encoding="utf-8", newline="") as handle:
        rows = {row["trait_key"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["METPO:2000303|EC:3.2.1.21"]["vocabulary_change"] == "enters"
    unlinked = rows["METPO:1000608"]
    assert (unlinked["ncbitaxon_nodes_on"], unlinked["unlinked_strains"]) == ("0", "1")
    assert unlinked["vocabulary_change"] == "outside"
