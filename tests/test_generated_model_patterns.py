"""Characterize generated regex boundaries before the LinkML upgrade (#1043)."""

from __future__ import annotations

import runpy
from pathlib import Path
from typing import Any

import pytest
from pydantic import BaseModel, ValidationError

from mtp.datamodel import (
    core_products_pydantic as core,
)
from mtp.datamodel import (
    mtp_pydantic as current,
)
from mtp.datamodel import (
    multiple_testing_pydantic as multiple,
)
from mtp.datamodel import (
    workflow_pydantic as workflow,
)

# The public bsdbng.datamodel compatibility alias points to mtp. Load the
# separately generated legacy file explicitly so its validators are covered.
LEGACY = runpy.run_path(str(Path(__file__).parents[1] / "src/bsdbng/datamodel/bsdbng_pydantic.py"))

DIGEST = "a" * 64
STUDY = {"id": "bsdb:1", "source_record_id": "1", "experiments": []}
BAD_DIGESTS = ("a" * 63, "a" * 65, "A" * 64, "x" + DIGEST, DIGEST + "x", DIGEST + "\nextra")


@pytest.mark.parametrize("model", [LEGACY["Study"], current.Study], ids=["legacy", "current"])
@pytest.mark.parametrize(
    ("field", "valid", "invalid"),
    [
        (
            "doi",
            ("10.1234/example", "10.any nonempty suffix", "10.1234/example\n"),
            ("x10.1234/example", "10.", "\n10.1234/example", "10.example\nextra"),
        ),
        (
            "id",
            ("bsdb:1", "bsdb:PMC1", "bsdb:any nonempty suffix", "bsdb:1\n"),
            ("xbsdb:1", "bsdb:", "\nbsdb:1", "bsdb:1\nextra"),
        ),
    ],
)
def test_study_patterns_keep_existing_permissive_suffixes(
    model: type[BaseModel], field: str, valid: tuple[str, ...], invalid: tuple[str, ...]
) -> None:
    for value in valid:
        assert model.model_validate({**STUDY, field: value}).model_dump()[field] == value
    for value in invalid:
        with pytest.raises(ValidationError) as failure:
            model.model_validate({**STUDY, field: value})
        assert {error["loc"] for error in failure.value.errors()} == {(field,)}


@pytest.mark.parametrize("model", [LEGACY["Taxon"], current.Taxon], ids=["legacy", "current"])
def test_taxon_patterns_reject_embedded_curie_and_nondigit_suffix(model: type[BaseModel]) -> None:
    for value in ("NCBITaxon:1", "NCBITaxon:9606", "NCBITaxon:1\n"):
        assert (
            model.model_validate({"id": value, "taxon_name": "Example taxon"}).model_dump()["id"]
            == value
        )
    for value in ("xNCBITaxon:1", "NCBITaxon:1x", "NCBITaxon:", "\nNCBITaxon:1"):
        with pytest.raises(ValidationError) as failure:
            model.model_validate({"id": value, "taxon_name": "Example taxon"})
        assert {error["loc"] for error in failure.value.errors()} == {("id",)}


@pytest.mark.parametrize(
    ("model", "required"),
    [
        pytest.param(
            workflow.FileInventoryEntry, {"path": "input.tsv", "size_bytes": 0}, id="workflow"
        ),
        pytest.param(
            core.SourceFileIdentity,
            {"product_id": "raw_dump", "path": "input.tsv", "size_bytes": 0},
            id="core-products",
        ),
        pytest.param(multiple.SourceArtifact, {"identity": "fixture:input"}, id="multiple-testing"),
    ],
)
def test_generated_digest_patterns_preserve_length_case_and_boundaries(
    model: type[BaseModel], required: dict[str, Any]
) -> None:
    # Python's existing ^...$ patterns allow a final newline; this upgrade does
    # not silently replace the schema policy with fullmatch or stricter anchors.
    for value in (DIGEST, DIGEST + "\n"):
        assert model.model_validate({**required, "sha256": value}).model_dump()["sha256"] == value
    for value in BAD_DIGESTS:
        with pytest.raises(ValidationError) as failure:
            model.model_validate({**required, "sha256": value})
        assert {error["loc"] for error in failure.value.errors()} == {("sha256",)}


def test_generated_profile_filter_id_requires_exact_prefix_and_digest() -> None:
    required = {
        "study_id": "bsdb:1",
        "experiment_id": "bsdb:1/1",
        "inc_dec": "inc",
        "feature_counts": [],
    }
    for value in ("sha256:" + DIGEST, "sha256:" + DIGEST + "\n"):
        result = core.CountProfileRow.model_validate({**required, "trait_filter_id": value})
        assert result.trait_filter_id == value
    invalid = ("xsha256:" + DIGEST, "SHA256:" + DIGEST, DIGEST)
    for value in (*invalid, *("sha256:" + digest for digest in BAD_DIGESTS)):
        with pytest.raises(ValidationError) as failure:
            core.CountProfileRow.model_validate({**required, "trait_filter_id": value})
        assert {error["loc"] for error in failure.value.errors()} == {("trait_filter_id",)}
