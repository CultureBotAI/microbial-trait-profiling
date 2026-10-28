"""Workflow validation must require run-bound policy and metadata format (#1135)."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import NoReturn

import pytest

from mtp.workflow.contracts import validate_product
from mtp.workflow.coverage import write_feature_catalog
from mtp.workflow.definition import read_definition
from mtp.workflow.matrix_contracts import (
    validate_feature_node_metadata,
    validate_matrix_sources,
)
from mtp.workflow.operations import trace_profile_feature, validate_core_product

ROOT = Path(__file__).resolve().parents[1]


def _unexpected_read(*_args: object, **_kwargs: object) -> NoReturn:
    raise AssertionError("Policy omission reached product I/O")


@pytest.fixture
def calls(tmp_path: Path) -> dict[str, tuple[Callable[..., object], tuple[object, ...]]]:
    product = next(
        product
        for product in read_definition(ROOT / "workflows" / "bugsigdb_core.yaml").products
        if product.id == "kg_nodes"
    )
    missing = tmp_path / "must-not-read.tsv"
    return {
        "core_product": (validate_core_product, (product, missing, {})),
        "product": (validate_product, (missing, "feature_metadata")),
        "matrix_sources": (validate_matrix_sources, (missing,) * 5),
        "feature_nodes": (validate_feature_node_metadata, (missing, missing)),
        "feature_catalog": (write_feature_catalog, (missing, missing)),
        "trace": (
            trace_profile_feature,
            (missing, missing, missing, missing, ("study", "experiment", "inc"), "METPO:1000605"),
        ),
    }


@pytest.mark.parametrize(
    "api", ["core_product", "product", "matrix_sources", "feature_nodes", "feature_catalog"]
)
@pytest.mark.parametrize("omitted", ["trait_key_policy", "include_policy_metadata", "both"])
def test_validation_requires_explicit_policy_and_format_before_io(
    calls: dict[str, tuple[Callable[..., object], tuple[object, ...]]],
    api: str,
    omitted: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    function, arguments = calls[api]
    keywords: dict[str, object] = {
        "trait_key_policy": "raw_object_v1",
        "include_policy_metadata": False,
    }
    if omitted == "both":
        keywords.clear()
    else:
        keywords.pop(omitted)
    if api == "core_product":
        # A parameters mapping is not a substitute for resolving the producing run.
        keywords["parameters"] = {"trait_key_policy": "raw_object_v1"}
    monkeypatch.setattr(Path, "open", _unexpected_read)
    match = "trait_key_policy" if omitted == "both" else omitted
    with pytest.raises(TypeError, match=match):
        function(*arguments, **keywords)


def test_trace_requires_explicit_policy_before_io(
    calls: dict[str, tuple[Callable[..., object], tuple[object, ...]]],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    function, arguments = calls["trace"]
    monkeypatch.setattr(Path, "open", _unexpected_read)
    with pytest.raises(TypeError, match="trait_key_policy"):
        function(*arguments)


def test_core_validation_cannot_replace_explicit_none_with_parameter_policy(
    calls: dict[str, tuple[Callable[..., object], tuple[object, ...]]],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    function, arguments = calls["core_product"]
    monkeypatch.setattr(Path, "open", _unexpected_read)
    with pytest.raises(ValueError, match="trait_key_policy"):
        function(
            *arguments,
            trait_key_policy=None,
            include_policy_metadata=False,
            parameters={"trait_key_policy": "raw_object_v1"},
        )


@pytest.mark.parametrize("policy", ["raw_object_v1", "facultative_union_v1"])
def test_new_execution_adapter_forwards_resolved_policy_and_modern_format(
    calls: dict[str, tuple[Callable[..., object], tuple[object, ...]]],
    monkeypatch: pytest.MonkeyPatch,
    policy: str,
) -> None:
    from mtp.workflow import operations

    _, arguments = calls["core_product"]
    forwarded: list[dict[str, object]] = []

    def capture(*_args: object, **keywords: object) -> list[object]:
        forwarded.append(keywords)
        return []

    monkeypatch.setattr(operations, "validate_core_product", capture)
    adapter: Callable[..., object] = operations.validate_execution_product
    parameters = {"trait_key_policy": policy, "strain_rollup": "off"}
    adapter(*arguments, parameters=parameters)
    assert forwarded == [
        {
            "parameters": parameters,
            "trait_key_policy": policy,
            "include_policy_metadata": True,
        }
    ]


@pytest.mark.parametrize("parameters", [{}, {"trait_key_policy": ""}])
def test_new_execution_adapter_has_no_policy_fallback(
    calls: dict[str, tuple[Callable[..., object], tuple[object, ...]]],
    monkeypatch: pytest.MonkeyPatch,
    parameters: dict[str, str],
) -> None:
    from mtp.workflow import operations

    _, arguments = calls["core_product"]
    monkeypatch.setattr(Path, "open", _unexpected_read)
    adapter: Callable[..., object] = operations.validate_execution_product
    with pytest.raises((KeyError, ValueError), match="trait_key_policy"):
        adapter(*arguments, parameters=parameters)
