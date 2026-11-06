"""Literal baseline outputs for paper readers intentionally unchanged by #1344."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from mtp.composition_denominators import load_retained_profiles
from mtp.contribution_inputs import load_graph
from mtp.native_composition import load_profile_matrix

PROFILES = (
    b"study_id\texperiment_id\tinc_dec\tMETPO:1000608\tMETPO:1000605\tMETPO:1000612\t"
    b"METPO:2000517|mediadive.medium:1\n"
    b"zero\tz1\tinc\t0\t0\t0\t0\n"
    b"zero\tz1\tdec\t0\t0\t0\t0\n"
    b"media_only\tm1\tinc\t0\t0\t0\t11\n"
    b"media_only\tm1\tdec\t0\t0\t0\t5\n"
    b"unpaired\tu1\tinc\t2\t3\t4\t1\n"
    b"mouse_study\tm2\tinc\t1\t0\t0\t0\n"
    b"mouse_study\tm2\tdec\t0\t0\t2\t0\n"
    b"reversed\tr1\tinc\t1\t2\t0\t0\n"
    b"reversed\tr1\tdec\t7\t0\t4\t0\n"
    b"normal\tn1\tinc\t4\t1\t0\t0\n"
    b"normal\tn1\tdec\t0\t3\t2\t0\n"
)

PATHS = ("data/derived/study_trait_profiles_wide.tsv", "isolated/profiles.tsv")


def _profiles(tmp_path: Path, relative: str) -> Path:
    path = tmp_path / relative
    path.parent.mkdir(parents=True)
    path.write_bytes(PROFILES)
    # These standalone readers do not infer a new population or orientation from
    # neighboring files or canonical-looking names. Selection belongs to callers.
    (path.parent / "experiments.tsv").write_text(
        "study_id\texperiment_id\thost_species\torientation_action\n"
        "mouse_study\tm2\tMus musculus\tkeep\n"
        "reversed\tr1\tHomo sapiens\tswap\n",
        encoding="utf-8",
    )
    return path


@pytest.mark.parametrize("relative", PATHS, ids=("canonical_name", "alternate_name"))
def test_native_composition_default_preserves_literal_profile_bytes(
    tmp_path: Path, relative: str
) -> None:
    path = _profiles(tmp_path, relative)

    matrix = load_profile_matrix(path)

    assert matrix.to_csv(sep="\t", lineterminator="\n").encode() == (
        b"\tMETPO:1000608\tMETPO:1000605\tMETPO:1000612\n"
        b"zero|z1|inc\t0\t0\t0\n"
        b"zero|z1|dec\t0\t0\t0\n"
        b"media_only|m1|inc\t0\t0\t0\n"
        b"media_only|m1|dec\t0\t0\t0\n"
        b"unpaired|u1|inc\t2\t3\t4\n"
        b"mouse_study|m2|inc\t1\t0\t0\n"
        b"mouse_study|m2|dec\t0\t0\t2\n"
        b"reversed|r1|inc\t1\t2\t0\n"
        b"reversed|r1|dec\t7\t0\t4\n"
        b"normal|n1|inc\t4\t1\t0\n"
        b"normal|n1|dec\t0\t3\t2\n"
    )
    assert path.read_bytes() == PROFILES


@pytest.mark.parametrize("relative", PATHS, ids=("canonical_name", "alternate_name"))
def test_retained_profiles_default_preserves_rows_and_metadata_order(
    tmp_path: Path, relative: str
) -> None:
    path = _profiles(tmp_path, relative)
    metadata = pd.DataFrame(
        {"total_trait_signal": [11, 0, 2, 9, 0, 5, 0, 1, 0, 5, 3]},
        index=[
            "reversed|r1|dec",
            "zero|z1|inc",
            "mouse_study|m2|dec",
            "unpaired|u1|inc",
            "media_only|m1|inc",
            "normal|n1|dec",
            "zero|z1|dec",
            "mouse_study|m2|inc",
            "media_only|m1|dec",
            "normal|n1|inc",
            "reversed|r1|inc",
        ],
    )

    matrix = load_retained_profiles(path, metadata)

    assert matrix.to_csv(sep="\t", lineterminator="\n").encode() == (
        b"\tMETPO:1000608\tMETPO:1000605\tMETPO:1000612\n"
        b"reversed|r1|dec\t7\t0\t4\n"
        b"zero|z1|inc\t0\t0\t0\n"
        b"mouse_study|m2|dec\t0\t0\t2\n"
        b"unpaired|u1|inc\t2\t3\t4\n"
        b"media_only|m1|inc\t0\t0\t0\n"
        b"normal|n1|dec\t0\t3\t2\n"
        b"zero|z1|dec\t0\t0\t0\n"
        b"mouse_study|m2|inc\t1\t0\t0\n"
        b"media_only|m1|dec\t0\t0\t0\n"
        b"normal|n1|inc\t4\t1\t0\n"
        b"reversed|r1|inc\t1\t2\t0\n"
    )
    assert path.read_bytes() == PROFILES


def test_contribution_graph_default_keeps_native_union_and_phenotype_rollup(
    tmp_path: Path,
) -> None:
    edges = tmp_path / "edges.tsv"
    edges.write_text(
        "id\tsubject\tpredicate\tobject\tprimary_knowledge_source\n"
        "d608\tNCBITaxon:1\tbiolink:has_phenotype\tMETPO:1000608\tinfores:gold\n"
        "d605\tNCBITaxon:2\tbiolink:has_phenotype\tMETPO:1000605\tinfores:gold\n"
        "d612\tNCBITaxon:3\tbiolink:has_phenotype\tMETPO:1000612\tinfores:bactotraits\n"
        "parent\tNCBITaxon:1\tbiolink:subclass_of\tNCBITaxon:10\tinfores:ncbitaxon\n"
        "strain\tkgmicrobe.strain:4\tbiolink:subclass_of\tNCBITaxon:4\tinfores:bacdive\n"
        "s605\tkgmicrobe.strain:4\tbiolink:has_phenotype\tMETPO:1000605\tinfores:bacdive\n"
        "clinical\tkgmicrobe.strain:4\tbiolink:has_phenotype\tMETPO:1004004\tinfores:bacdive\n"
        "spath\tkgmicrobe.strain:4\tMETPO:2000103\tGO:0030245\tinfores:bacdive\n"
        "dpath\tNCBITaxon:1\tMETPO:2000103\tGO:0030245\tinfores:gold\n",
        encoding="utf-8",
    )
    nodes = tmp_path / "nodes.tsv"
    nodes.write_text(
        "id\tname\nNCBITaxon:1\tSpecies one\nNCBITaxon:2\tSpecies two\n"
        "NCBITaxon:3\tSpecies three\nNCBITaxon:4\tSpecies four\nNCBITaxon:10\tGenus\n",
        encoding="utf-8",
    )
    columns = {
        "METPO:1000605",
        "METPO:1000608",
        "METPO:1000612",
        "METPO:1004004",
        "METPO:2000103|GO:0030245",
    }

    graph = load_graph(edges, nodes, columns, columns)
    serialized = {
        "parents": {key: sorted(value) for key, value in graph.parents.items()},
        "traits": {key: sorted(value) for key, value in graph.traits.items()},
        "assertions": [
            [
                taxon,
                trait,
                assertion.edge_id,
                assertion.predicate,
                assertion.object_id,
                assertion.primary_source,
                assertion.knowledge_source,
                assertion.knowledge_level,
                assertion.source_row,
            ]
            for (taxon, trait), assertions in sorted(graph.assertions.items())
            for assertion in assertions
        ],
        "names": graph.names,
        "header_repaired": graph.header_repaired,
        "continued_records": graph.continued_records,
    }

    assert json.dumps(serialized, sort_keys=True, separators=(",", ":")).encode() == (
        b'{"assertions":['
        b'["NCBITaxon:1","METPO:1000608","d608","biolink:has_phenotype",'
        b'"METPO:1000608","infores:gold","","",2],'
        b'["NCBITaxon:1","METPO:2000103|GO:0030245","dpath","METPO:2000103",'
        b'"GO:0030245","infores:gold","","",10],'
        b'["NCBITaxon:2","METPO:1000608","d605","biolink:has_phenotype",'
        b'"METPO:1000605","infores:gold","","",3],'
        b'["NCBITaxon:3","METPO:1000608","d612","biolink:has_phenotype",'
        b'"METPO:1000612","infores:bactotraits","","",4],'
        b'["NCBITaxon:4","METPO:1000608","s605","biolink:has_phenotype",'
        b'"METPO:1000605","infores:bacdive","","",7]],'
        b'"continued_records":0,"header_repaired":false,'
        b'"names":{"NCBITaxon:1":"Species one","NCBITaxon:10":"Genus",'
        b'"NCBITaxon:2":"Species two","NCBITaxon:3":"Species three",'
        b'"NCBITaxon:4":"Species four"},'
        b'"parents":{"NCBITaxon:1":["NCBITaxon:10"]},'
        b'"traits":{"NCBITaxon:1":["METPO:1000608","METPO:2000103|GO:0030245"],'
        b'"NCBITaxon:2":["METPO:1000608"],"NCBITaxon:3":["METPO:1000608"],'
        b'"NCBITaxon:4":["METPO:1000608"]}}'
    )
