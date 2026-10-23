"""Account for the traits each BacDive strain rollup scope adds (#774, #788)."""

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
    # Strain 1 adds a phenotype its species lacks and a clinical class; strain 2 adds a
    # phenotype and both an enzyme activity and its absence; strain 3 has no taxon link;
    # strain 4 links to a node outside the taxonomy.
    paths["edges"].write_text(
        "id\tsubject\tpredicate\tobject\n"
        "e1\tNCBITaxon:101\tbiolink:subclass_of\tNCBITaxon:100\n"
        "e2\tNCBITaxon:102\tbiolink:subclass_of\tNCBITaxon:100\n"
        "e3\tNCBITaxon:101\tbiolink:has_phenotype\tMETPO:1000603\n"
        "e4\tNCBITaxon:102\tbiolink:has_phenotype\tMETPO:1000605\n"
        "s1\tkgmicrobe.strain:1\tbiolink:subclass_of\tNCBITaxon:101\n"
        "s2\tkgmicrobe.strain:1\tbiolink:has_phenotype\tMETPO:1000605\n"
        "s3\tkgmicrobe.strain:1\tbiolink:has_phenotype\tMETPO:1004004\n"
        "s4\tkgmicrobe.strain:2\tbiolink:subclass_of\tNCBITaxon:102\n"
        "s5\tkgmicrobe.strain:2\tbiolink:has_phenotype\tMETPO:1000678\n"
        "s6\tkgmicrobe.strain:2\tMETPO:2000302\tEC:3.2.1.21\n"
        "s7\tkgmicrobe.strain:2\tMETPO:2000303\tEC:3.2.1.21\n"
        "s8\tkgmicrobe.strain:3\tbiolink:has_phenotype\tMETPO:1000608\n"
        "s9\tkgmicrobe.strain:4\tbiolink:subclass_of\tNCBITaxon:999\n"
        "s10\tkgmicrobe.strain:4\tbiolink:has_phenotype\tMETPO:1000668\n",
        encoding="utf-8",
    )
    paths["nodes"].write_text(
        "id\tcategory\tname\n"
        "METPO:1000603\tbiolink:PhenotypicQuality\tanaerobic\n"
        "METPO:1000605\tbiolink:PhenotypicQuality\tfacultatively anaerobic\n"
        "METPO:1000608\tbiolink:PhenotypicQuality\tfacultatively aerobic\n"
        "METPO:1000668\tbiolink:PhenotypicQuality\tcoccus shaped\n"
        "METPO:1000678\tbiolink:PhenotypicQuality\toval shaped\n"
        "METPO:1004004\tbiolink:PhenotypicQuality\thuman pathogen\n"
        "METPO:2000302\tbiolink:OntologyClass\tshows activity of\n"
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


def test_audit_reports_what_each_scope_adds(tmp_path: Path) -> None:
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
    strains = summary["strains"]
    assert strains["by_ncbitaxon_link_count"] == {"0": 1, "1": 3}
    assert strains["strain_trait_pairs_without_ncbitaxon_link"] == 1
    assert strains["orphan_links"] == {
        "linked_nodes": 3,
        "orphan_nodes": 1,
        "orphan_strains": 1,
        "orphan_nodes_reported": 0,
    }
    scopes = summary["scopes"]
    assert {scope: scopes[scope]["retained_vocabulary"] for scope in scopes} == {
        "off": 2,
        "phenotype": 3,
        "full": 5,
    }
    assert scopes["phenotype"]["entering_by_predicate"] == {"phenotype (METPO object)": 1}
    assert scopes["full"]["entering_by_predicate"] == {
        "phenotype (METPO object)": 2,
        "shows activity of": 1,
    }
    anaerobic = {scope: scopes[scope]["focal_traits"]["METPO:1000605"] for scope in scopes}
    assert [anaerobic[s]["ncbitaxon_nodes"] for s in ("off", "phenotype", "full")] == [1, 2, 2]
    assert [anaerobic[s]["reported_taxa"] for s in ("off", "phenotype", "full")] == [2, 3, 3]
    human = summary["signature_taxa"]["cohorts"][1]
    assert human["phenotype"]["host_species"] == "Homo sapiens"
    assert [human[s]["gaining_a_trait_label"] for s in ("off", "phenotype", "full")] == [0, 2, 2]

    with (output / "trait_coverage.tsv").open(encoding="utf-8", newline="") as handle:
        rows = {row["trait_key"]: row for row in csv.DictReader(handle, delimiter="\t")}
    absent = rows["METPO:2000303|EC:3.2.1.21"]
    assert (absent["negated_predicate"], absent["vocabulary_full"]) == ("True", "False")
    assert absent["ncbitaxon_nodes_full"] == "0"
    clinical = rows["METPO:1004004"]
    assert (clinical["excluded_phenotype"], clinical["vocabulary_phenotype"]) == ("True", "False")
    assert clinical["vocabulary_full"] == "True"
    unlinked = rows["METPO:1000608"]
    assert (unlinked["unlinked_strains"], unlinked["vocabulary_full"]) == ("1", "False")
