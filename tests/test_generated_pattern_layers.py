"""Exercise #1316 through LinkML, its JSON Schema output, and generated model serialization."""

from __future__ import annotations

import copy
import importlib
import json
import runpy
import shutil
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol, cast

import pytest
from pydantic import BaseModel, ValidationError

from mtp.datamodel import (
    core_products_pydantic,
    mtp_pydantic,
    multiple_testing_pydantic,
    workflow_pydantic,
)

ROOT = Path(__file__).resolve().parents[1]
DIGEST = "a" * 64


class PatternError(Protocol):
    @property
    def validator(self) -> str: ...

    @property
    def absolute_path(self) -> Sequence[str | int]: ...


class ValidationFinding(Protocol):
    @property
    def severity(self) -> str: ...

    @property
    def source(self) -> PatternError | None: ...


class ValidationReport(Protocol):
    @property
    def results(self) -> list[ValidationFinding]: ...


class LinkMLValidator(Protocol):
    def validate(self, instance: object, target_class: str | None = None) -> ValidationReport: ...


class ValidatorFactory(Protocol):
    def __call__(self, schema: Path, *, validation_plugins: list[object]) -> LinkMLValidator: ...


class PluginFactory(Protocol):
    def __call__(self, *, closed: bool = False, json_schema_path: Path | None = None) -> object: ...


class SlotPattern(Protocol):
    @property
    def pattern(self) -> str | None: ...


class SchemaView(Protocol):
    def induced_slot(self, slot_name: str, class_name: str) -> SlotPattern: ...


class JSONSchemaGenerator(Protocol):
    @property
    def schemaview(self) -> SchemaView: ...

    def serialize(self) -> str: ...


class GeneratorFactory(Protocol):
    def __call__(
        self, schema: str, *, mergeimports: bool, top_class: str
    ) -> JSONSchemaGenerator: ...


def _factory(module: str, name: str) -> object:
    """Bridge the used API of the installed, untyped LinkML dependency explicitly."""
    value: object = getattr(importlib.import_module(module), name)
    assert callable(value), f"Missing callable LinkML API: {module}.{name}"
    return value


@dataclass(frozen=True)
class Case:
    schema: str
    class_name: str
    field: str
    instance: dict[str, Any]

    @property
    def key(self) -> str:
        return f"{self.schema}/{self.class_name}/{self.field}"

    def values(self) -> list[str]:
        value = self.instance[self.field]
        return value if isinstance(value, list) else [value]

    def with_final_newline(self) -> dict[str, Any]:
        altered = copy.deepcopy(self.instance)
        value = altered[self.field]
        if isinstance(value, list):
            value[-1] += "\n"
        else:
            altered[self.field] = value + "\n"
        return altered


TAXON = {"id": "NCBITaxon:9606", "taxon_name": "Fixture taxon"}
SIGNATURE = {"id": "bsdb:PMC1-2-3", "direction": "increased", "taxa": [TAXON]}
EXPERIMENT = {
    "id": "bsdb:PMC1-2",
    "signatures": [SIGNATURE],
    "body_site_ontology_id": ["UBERON:1", "UBERON:0001988"],
}
STUDY = {
    "id": "bsdb:nonempty suffix with spaces",
    "source_record_id": "fixture",
    "doi": "10.any nonempty suffix",
    "experiments": [EXPERIMENT],
}
COUNT = {
    "study_id": "bsdb:1",
    "experiment_id": "bsdb:1/1",
    "inc_dec": "inc",
    "feature_counts": [{"feature_id": "METPO:1000605", "count": 2}],
    "trait_filter_id": "sha256:" + DIGEST,
    "unit": "rolled_up_source_taxon_count",
}
CASES = (
    *(
        Case(schema, class_name, field, instance)
        for schema in ("bsdbng", "mtp")
        for class_name, field, instance in (
            ("Study", "id", STUDY),
            ("Study", "doi", STUDY),
            ("Experiment", "id", EXPERIMENT),
            ("Experiment", "body_site_ontology_id", EXPERIMENT),
            ("Signature", "id", SIGNATURE),
            ("Taxon", "id", TAXON),
            (
                "TaxonTraitMatrixRow",
                "taxon_id",
                {"taxon_id": "NCBITaxon:1", "taxon_name": "Fixture taxon"},
            ),
            ("Phenotype", "phenotype_id", {"phenotype_id": "METPO:1000605"}),
        )
    ),
    Case(
        "core_products",
        "SourceFileIdentity",
        "sha256",
        {"product_id": "raw_dump", "path": "input.tsv", "sha256": DIGEST, "size_bytes": 0},
    ),
    *(
        Case("core_products", "CountProfileRow", field, COUNT)
        for field in ("trait_filter_id", "inc_dec", "unit")
    ),
    Case(
        "multiple_testing",
        "SourceArtifact",
        "sha256",
        {"identity": "fixture:input", "sha256": DIGEST},
    ),
    Case(
        "workflow",
        "FileInventoryEntry",
        "sha256",
        {"path": "input.tsv", "sha256": DIGEST, "size_bytes": 0},
    ),
    Case("workflow", "ParameterDefinition", "name", {"name": "Alpha_1"}),
)


def _patterns(value: Any) -> list[str]:
    if isinstance(value, dict):
        result = [value["pattern"]] if "pattern" in value else []
        for child in value.values():
            result.extend(_patterns(child))
        return result
    if isinstance(value, list):
        return [pattern for child in value for pattern in _patterns(child)]
    return []


@dataclass
class Layers:
    linkml: dict[str, LinkMLValidator]
    exported: dict[tuple[str, str], LinkMLValidator]
    patterns: dict[str, str]


@pytest.fixture(scope="module")
def layers(tmp_path_factory: pytest.TempPathFactory) -> Layers:
    """Generate each complete source schema once, then select real class definitions."""
    output = tmp_path_factory.mktemp("pattern-json-schema")
    validator_factory = cast(ValidatorFactory, _factory("linkml.validator", "Validator"))
    plugin_factory = cast(
        PluginFactory, _factory("linkml.validator.plugins", "JsonschemaValidationPlugin")
    )
    generator_factory = cast(
        GeneratorFactory, _factory("linkml.generators.jsonschemagen", "JsonSchemaGenerator")
    )
    live: dict[str, LinkMLValidator] = {}
    exported: dict[tuple[str, str], LinkMLValidator] = {}
    patterns: dict[str, str] = {}
    for name in sorted({case.schema for case in CASES}):
        source = ROOT / "schema" / f"{name}.yaml"
        cases = [case for case in CASES if case.schema == name]
        live[name] = validator_factory(source, validation_plugins=[plugin_factory(closed=True)])
        generator = generator_factory(str(source), mergeimports=True, top_class=cases[0].class_name)
        document = json.loads(generator.serialize())
        assert "$defs" in document, f"LinkML omitted class definitions for {name}"
        for case in cases:
            assert case.class_name in document["$defs"], case.key
            definition = document["$defs"][case.class_name]
            assert case.field in definition["properties"], case.key
            actual = _patterns(definition["properties"][case.field])
            source_pattern = generator.schemaview.induced_slot(case.field, case.class_name).pattern
            assert isinstance(source_pattern, str), case.key
            assert source_pattern, case.key
            assert actual == [source_pattern], f"JSON Schema lost the source constraint: {case.key}"
            patterns[case.key] = source_pattern
            pair = name, case.class_name
            if pair in exported:
                continue
            selected = {
                "$schema": document["$schema"],
                "$defs": document["$defs"],
                "$ref": f"#/$defs/{case.class_name}",
            }
            schema_path = output / f"{name}-{case.class_name}.json"
            with schema_path.open("x", encoding="utf-8") as handle:
                json.dump(selected, handle)
            assert schema_path.stat().st_size > 0
            exported[pair] = validator_factory(
                source,
                validation_plugins=[plugin_factory(json_schema_path=schema_path)],
            )
    return Layers(live, exported, patterns)


@pytest.mark.parametrize("case", CASES, ids=lambda case: case.key)
@pytest.mark.parametrize("route", ["linkml", "exported_json_schema"])
def test_linkml_and_exported_json_schema_refuse_final_newline(
    layers: Layers, case: Case, route: str
) -> None:
    validator = (
        layers.linkml[case.schema]
        if route == "linkml"
        else layers.exported[(case.schema, case.class_name)]
    )
    valid = validator.validate(copy.deepcopy(case.instance), target_class=case.class_name)
    assert valid.results == [], f"Invalid positive control for {case.key}: {valid.results}"
    invalid = validator.validate(case.with_final_newline(), target_class=case.class_name)
    failures = [result for result in invalid.results if result.severity == "ERROR"]
    assert len(failures) == 1, f"{route} accepted final LF in {case.key}: {invalid.results}"
    # A different required field or format failure must not stand in for the pattern witness.
    error = failures[0].source
    assert error is not None
    assert error.validator == "pattern"
    expected_path: list[str | int] = [case.field]
    if isinstance(case.instance[case.field], list):
        expected_path.append(len(case.instance[case.field]) - 1)
    assert list(error.absolute_path) == expected_path


@pytest.fixture(scope="module")
def model_classes() -> dict[tuple[str, str], type[BaseModel]]:
    # Load the tracked legacy product itself; its public package alias points to mtp.
    legacy = runpy.run_path(str(ROOT / "src/bsdbng/datamodel/bsdbng_pydantic.py"))
    modules = {
        "mtp": mtp_pydantic,
        "core_products": core_products_pydantic,
        "multiple_testing": multiple_testing_pydantic,
        "workflow": workflow_pydantic,
    }
    return {
        (case.schema, case.class_name): cast(
            type[BaseModel],
            legacy[case.class_name]
            if case.schema == "bsdbng"
            else getattr(modules[case.schema], case.class_name),
        )
        for case in CASES
    }


@pytest.mark.parametrize("case", CASES, ids=lambda case: case.key)
def test_valid_synthetic_values_keep_exact_serialization_and_invalid_values_are_not_trimmed(
    model_classes: dict[tuple[str, str], type[BaseModel]], case: Case
) -> None:
    model = model_classes[(case.schema, case.class_name)]
    instance = copy.deepcopy(case.instance)
    result = model.model_validate(instance)
    dumped = result.model_dump(mode="json", by_alias=True, exclude_unset=True)
    assert dumped == instance
    assert json.loads(result.model_dump_json(by_alias=True, exclude_unset=True)) == instance
    with pytest.raises(ValidationError) as failure:
        model.model_validate(case.with_final_newline())
    expected_location: tuple[str | int, ...] = (case.field,)
    # Generated list validators deliberately report their slot, not the element index.
    assert {error["loc"] for error in failure.value.errors()} == {expected_location}


def test_exported_patterns_have_the_same_lf_boundary_in_an_available_ecmascript_engine(
    layers: Layers,
) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip(
            "Optional ECMAScript probe: no existing Node executable; no dependency installed"
        )
    rows = [
        {"key": case.key, "pattern": layers.patterns[case.key], "value": value}
        for case in CASES
        for value in case.values()
    ]
    javascript = """
const fs = require('node:fs');
const rows = JSON.parse(fs.readFileSync(0, 'utf8'));
const result = rows.map(row => {
  const pattern = new RegExp(row.pattern);
  return {key: row.key, valid: pattern.test(row.value), invalid: pattern.test(row.value + '\\n')};
});
process.stdout.write(JSON.stringify(result));
"""
    completed = subprocess.run(
        [node, "-e", javascript],
        input=json.dumps(rows),
        text=True,
        capture_output=True,
        timeout=15,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    observed = json.loads(completed.stdout)
    assert observed == [{"key": row["key"], "valid": True, "invalid": False} for row in rows]
