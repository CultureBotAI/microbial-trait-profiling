"""Explicit archived-run regression for #1044; no owner arm is executed."""

from __future__ import annotations

import os
import re
from pathlib import Path

import pandas as pd
import pytest

from mtp.three_test_table.stage import verify_table


@pytest.mark.integration
def test_verified_off_cellulose_target() -> None:
    required = ("MTP_THREE_TEST_COMPARISON", "MTP_THREE_TEST_OUTPUT", "MTP_THREE_TEST_SCRATCH")
    optional = ("MTP_THREE_TEST_MULTIPLICITY", "MTP_THREE_TEST_LABEL_SUPPORT")
    if not any(name in os.environ for name in (*required, *optional, "MTP_THREE_TEST_BASIS")):
        pytest.skip("Requires a verified archived off comparison and completed sidecar")
    missing = [name for name in required if not os.environ.get(name, "").strip()]
    if missing:
        raise ValueError(f"Missing reproduction inputs: {', '.join(missing)}")
    basis = os.environ.get("MTP_THREE_TEST_BASIS", "")
    if re.fullmatch(r"[0-9a-f]{7,40}", basis, flags=re.IGNORECASE) is None:
        raise ValueError(
            "Reproduction basis MTP_THREE_TEST_BASIS must be 7-40 hexadecimal characters"
        )
    comparison, output, scratch = (Path(os.environ[name]) for name in required)
    multiplicity, label_support = (
        Path(os.environ[name]) if os.environ.get(name) else None for name in optional
    )
    record = verify_table(
        output,
        comparison,
        scratch=scratch,
        multiplicity=multiplicity,
        label_support=label_support,
    )
    assert record["strain_rollup"] == "off"
    assert record["contrast_population"] == "disease_vs_healthy_tier1_2_subject_controls"
    assert record["family_size"] == 317
    assert record["source_commits"]["core"].startswith(basis.lower())
    table = pd.read_csv(output / "per_trait_evidence.tsv", sep="\t").set_index("trait")
    assert len(table) == 317
    cellulose = table.loc["METPO:2000103|GO:0030245"]
    assert cellulose["composition_studies"] == 548
    assert cellulose["composition_pairs"] == 935
    assert cellulose["study_t_status"] == "tested"
    # The historical reference was published at these precisions; do not tune to it.
    assert cellulose["study_t_statistic"] == pytest.approx(-5.24, abs=0.005)
    assert cellulose["study_t_q_by"] == pytest.approx(0.00047, abs=0.000005)
