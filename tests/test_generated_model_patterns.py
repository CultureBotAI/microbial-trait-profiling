"""Guard LinkML schema regex anchors and actual generated-model boundaries (#1316)."""

from __future__ import annotations

import ast
import importlib.util
import re
import runpy
import shutil
from pathlib import Path
from typing import Any

import pytest
import yaml
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

# Public bsdbng.datamodel model imports forward to canonical mtp modules. Load the
# separately generated legacy file explicitly so its validators are covered.
REPO_ROOT = Path(__file__).parents[1]
_BSDBNG_SPEC = importlib.util.find_spec("bsdbng")
assert _BSDBNG_SPEC is not None
assert _BSDBNG_SPEC.origin is not None
_LEGACY_PATH = Path(_BSDBNG_SPEC.origin).parent / "datamodel/bsdbng_pydantic.py"
LEGACY = runpy.run_path(str(_LEGACY_PATH))
assert core.__file__ is not None
assert current.__file__ is not None
assert multiple.__file__ is not None
assert workflow.__file__ is not None
_MODEL_ROOT = Path(multiple.__file__).resolve().parents[3]
SCHEMA_DIR = _MODEL_ROOT / "schema" if (_MODEL_ROOT / "schema").is_dir() else REPO_ROOT / "schema"
GENERATED_MODEL_PATHS: tuple[Path, ...] = (
    _LEGACY_PATH,
    Path(core.__file__),
    Path(current.__file__),
    Path(multiple.__file__),
    Path(workflow.__file__),
)

DIGEST = "a" * 64
STUDY = {"id": "bsdb:1", "source_record_id": "1", "experiments": []}
BAD_DIGESTS = ("a" * 63, "a" * 65, "A" * 64, "x" + DIGEST, DIGEST + "x", DIGEST + "\nextra")
ACTUAL_END = r"(?![\s\S])"
# Explicit witnesses pin the accepted body language, including every alternative.
# A new body must acquire a witness rather than silently escaping the inventory guard.
PATTERN_WITNESSES: dict[str, tuple[str, ...]] = {
    "bsdb:.+": ("bsdb:1", "bsdb:any nonempty suffix"),
    "bsdb:.+-[0-9]+": ("bsdb:PMC1-42",),
    "bsdb:.+-[0-9]+-[0-9]+": ("bsdb:PMC1-42-1",),
    "NCBITaxon:[0-9]+": ("NCBITaxon:9606",),
    r"10\..+": ("10.1234/example", "10.any nonempty suffix"),
    "UBERON:[0-9]+": ("UBERON:0001988",),
    "METPO:[0-9]+": ("METPO:1000447",),
    "[A-Za-z][A-Za-z0-9_.-]*": ("Product_1.2-3",),
    "[A-Za-z][A-Za-z0-9_]*": ("Parameter_1",),
    "[A-Za-z][A-Za-z0-9_:.-]*": ("Analysis_1:2.3-4",),
    "[0-9a-f]{64}": (DIGEST,),
    "sha256:[0-9a-f]{64}": ("sha256:" + DIGEST,),
    "(raw_object_v1|facultative_union_v1)": ("raw_object_v1", "facultative_union_v1"),
    "input_cohort_before_analysis_specific_exclusions": (
        "input_cohort_before_analysis_specific_exclusions",
    ),
    "imported_unverified": ("imported_unverified",),
    "1": ("1",),
    "(ilyas|noel|jed)": ("ilyas", "noel", "jed"),
    "standard_analysis_calibration_and_optional_figures": (
        "standard_analysis_calibration_and_optional_figures",
    ),
    "(count_profiles|ora_profiles)": ("count_profiles", "ora_profiles"),
    "(off|phenotype)": ("off", "phenotype"),
    "(tested|untestable)": ("tested", "untestable"),
    r"1\.2\.0": ("1.2.0",),
    "(inc|dec)": ("inc", "dec"),
    "rolled_up_source_taxon_count": ("rolled_up_source_taxon_count",),
    "distinct_matched_reported_taxon_count": ("distinct_matched_reported_taxon_count",),
    "(unknown|fixture|live|pinned)": ("unknown", "fixture", "live", "pinned"),
    "(pinned_bytes_verified|origin_unverified)": ("pinned_bytes_verified", "origin_unverified"),
    "(increased|decreased)": ("increased", "decreased"),
    "(direct_or_descendant_unresolved|direct_descendant_or_strain_unresolved)": (
        "direct_or_descendant_unresolved",
        "direct_descendant_or_strain_unresolved",
    ),
    "unresolved": ("unresolved",),
}


@pytest.mark.parametrize("model", [LEGACY["Study"], current.Study], ids=["legacy", "current"])
@pytest.mark.parametrize(
    ("field", "valid", "invalid"),
    [
        (
            "doi",
            ("10.1234/example", "10.any nonempty suffix", "10.suffix ", "10.suffix\t"),
            (
                "x10.1234/example",
                "10.",
                "\n10.1234/example",
                "10.example\nextra",
                "10.1234/example\n",
            ),
        ),
        (
            "id",
            ("bsdb:1", "bsdb:PMC1", "bsdb:any nonempty suffix", "bsdb:suffix ", "bsdb:suffix\t"),
            ("xbsdb:1", "bsdb:", "\nbsdb:1", "bsdb:1\nextra", "bsdb:1\n"),
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
    for value in ("NCBITaxon:1", "NCBITaxon:9606"):
        assert (
            model.model_validate({"id": value, "taxon_name": "Example taxon"}).model_dump()["id"]
            == value
        )
    for value in ("xNCBITaxon:1", "NCBITaxon:1x", "NCBITaxon:", "\nNCBITaxon:1", "NCBITaxon:1\n"):
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
    assert model.model_validate({**required, "sha256": DIGEST}).model_dump()["sha256"] == DIGEST
    for value in (
        *BAD_DIGESTS,
        *(DIGEST + suffix for suffix in ("\n", "\r\n", "\r", " ", "\t", "\u2028", "\u2029")),
    ):
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
    value = "sha256:" + DIGEST
    result = core.CountProfileRow.model_validate({**required, "trait_filter_id": value})
    assert result.trait_filter_id == value
    invalid = ("xsha256:" + DIGEST, "SHA256:" + DIGEST, DIGEST, value + "\n")
    for value in (*invalid, *("sha256:" + digest for digest in BAD_DIGESTS)):
        with pytest.raises(ValidationError) as failure:
            core.CountProfileRow.model_validate({**required, "trait_filter_id": value})
        assert {error["loc"] for error in failure.value.errors()} == {("trait_filter_id",)}


@pytest.mark.parametrize(
    "model", [LEGACY["Experiment"], current.Experiment], ids=["legacy", "current"]
)
@pytest.mark.parametrize("bad_index", [0, 1])
@pytest.mark.parametrize("suffix", ["\n", "\r\n", "\r", " ", "\t", "\u2028", "\u2029"])
def test_multivalued_ontology_ids_reject_suffix_without_normalizing(
    model: type[BaseModel], bad_index: int, suffix: str
) -> None:
    values = ["UBERON:0001988", "UBERON:0000178"]
    required = {"id": "bsdb:1-1", "signatures": []}
    validated = model.model_validate({**required, "body_site_ontology_id": values})
    assert validated.model_dump()["body_site_ontology_id"] == values
    malformed = values.copy()
    malformed[bad_index] += suffix
    with pytest.raises(ValidationError) as failure:
        model.model_validate({**required, "body_site_ontology_id": malformed})
    assert {error["loc"] for error in failure.value.errors()} == {("body_site_ontology_id",)}
    assert malformed[bad_index] == values[bad_index] + suffix


@pytest.mark.parametrize("models", [LEGACY, vars(current)], ids=["legacy", "current"])
@pytest.mark.parametrize(
    ("model_name", "field", "value", "required"),
    [
        ("Signature", "id", "bsdb:PMC1-42-1", {"direction": "increased", "taxa": []}),
        ("Phenotype", "phenotype_id", "METPO:1000447", {}),
        ("TaxonTraitMatrixRow", "taxon_id", "NCBITaxon:9606", {"taxon_name": "Example taxon"}),
    ],
)
def test_remaining_biological_identifiers_reject_final_lf(
    models: dict[str, Any], model_name: str, field: str, value: str, required: dict[str, Any]
) -> None:
    model = models[model_name]
    assert model.model_validate({**required, field: value}).model_dump()[field] == value
    with pytest.raises(ValidationError) as failure:
        model.model_validate({**required, field: value + "\n"})
    assert {error["loc"] for error in failure.value.errors()} == {(field,)}


@pytest.mark.parametrize(
    ("model", "field", "value", "required"),
    [
        pytest.param(
            workflow.WorkflowDefinition,
            "id",
            "Workflow_1.2-3",
            {"version": "1", "products": [], "steps": []},
            id="workflow-id",
        ),
        pytest.param(
            workflow.AnalysisSpec,
            "id",
            "Analysis_1:2.3-4",
            {
                "estimand": "annotation",
                "kernel": "fixture",
                "transformation": "none",
                "sampling_unit": "study",
                "weighting": "equal",
                "eligibility_rule": "all",
                "multiple_testing_policy": "fixture",
            },
            id="analysis-id",
        ),
        pytest.param(
            core.CountProfileRow,
            "inc_dec",
            "inc",
            {
                "study_id": "bsdb:1",
                "experiment_id": "bsdb:1/1",
                "feature_counts": [],
                "trait_filter_id": "sha256:" + DIGEST,
            },
            id="count-direction",
        ),
        pytest.param(
            core.CountProfileRow,
            "unit",
            "rolled_up_source_taxon_count",
            {
                "study_id": "bsdb:1",
                "experiment_id": "bsdb:1/1",
                "inc_dec": "inc",
                "feature_counts": [],
                "trait_filter_id": "sha256:" + DIGEST,
            },
            id="count-unit",
        ),
        pytest.param(
            core.ORAProfileRow,
            "unit",
            "distinct_matched_reported_taxon_count",
            {
                "study_id": "bsdb:1",
                "experiment_id": "bsdb:1/1",
                "inc_dec": "dec",
                "sample_size": 0,
                "feature_counts": [],
                "trait_filter_id": "sha256:" + DIGEST,
            },
            id="ora-unit",
        ),
    ],
)
@pytest.mark.parametrize("suffix", ["\n", "\r\n", "\r", " ", "\t", "\u2028", "\u2029"])
def test_generated_identifier_and_fixed_token_patterns_reject_trailing_characters(
    model: type[BaseModel], field: str, value: str, required: dict[str, Any], suffix: str
) -> None:
    assert model.model_validate({**required, field: value}).model_dump()[field] == value
    with pytest.raises(ValidationError) as failure:
        model.model_validate({**required, field: value + suffix})
    assert {error["loc"] for error in failure.value.errors()} == {(field,)}


def _format_pattern_target(schema_name: str, trail: tuple[str, ...]) -> str:
    if len(trail) == 2 and trail[0] == "slots":
        return f"{schema_name}: slot '{trail[1]}' (slots.{trail[1]})"
    if len(trail) == 2 and trail[0] == "types":
        return f"{schema_name}: type '{trail[1]}' (types.{trail[1]})"
    if len(trail) == 4 and trail[0] == "classes" and trail[2] in {"attributes", "slot_usage"}:
        _, class_name, section, slot_name = trail
        return (
            f"{schema_name}: slot '{class_name}.{slot_name}' "
            f"(classes.{class_name}.{section}.{slot_name})"
        )
    return f"{schema_name}: {'.'.join(trail)}"


def _collect_schema_patterns(
    schema_dir: Path = SCHEMA_DIR,
) -> list[tuple[str, str, object]]:
    entries: list[tuple[str, str, object]] = []
    for path in sorted(schema_dir.glob("*.yaml")):
        document = yaml.safe_load(path.read_text(encoding="utf-8"))

        def walk(node: object, trail: tuple[str, ...], *, schema_name: str = path.name) -> None:
            if isinstance(node, dict):
                if "pattern" in node and not isinstance(node["pattern"], dict):
                    entries.append(
                        (
                            schema_name,
                            _format_pattern_target(schema_name, trail),
                            node["pattern"],
                        )
                    )
                for key, value in node.items():
                    if key == "pattern" and not isinstance(value, dict):
                        continue
                    walk(value, (*trail, str(key)), schema_name=schema_name)
            elif isinstance(node, list):
                for index, item in enumerate(node):
                    walk(item, (*trail, str(index)), schema_name=schema_name)

        walk(document, ())
    return entries


def _pattern_anchor_violation(pattern: object) -> str | None:
    """Check outer-anchor syntax independently of actual-end behavior."""
    if not isinstance(pattern, str) or not pattern:
        return "pattern must be a non-empty string"
    try:
        compiled = re.compile(pattern)
    except re.error as exc:
        return f"invalid regex ({exc})"
    if compiled.flags & re.MULTILINE:
        return "MULTILINE flag makes ^ and $ match internal line breaks under re.search"

    escaped = False
    in_char_class = False
    char_class_content_start = 0
    group_depth = 0
    has_top_level_alternation = False
    trailing_unescaped_dollar = False

    for index, char in enumerate(pattern):
        if escaped:
            escaped = False
            continue
        if char == "\\":
            escaped = True
            continue
        if in_char_class:
            if char == "]" and index > char_class_content_start:
                in_char_class = False
            continue
        if char == "[":
            in_char_class = True
            char_class_content_start = (
                index + 2 if pattern[index + 1 : index + 2] == "^" else index + 1
            )
        elif char == "(":
            group_depth += 1
        elif char == ")":
            group_depth -= 1
        elif char == "|" and group_depth == 0:
            has_top_level_alternation = True
        elif char == "$" and index == len(pattern) - 1 and group_depth == 0:
            trailing_unescaped_dollar = True

    if not pattern.startswith("^"):
        return "missing leading '^' anchor"
    if not trailing_unescaped_dollar:
        return "missing trailing '$' anchor"
    if pattern == "^$":
        return "empty anchored pattern '^$'"
    if has_top_level_alternation:
        return "ungrouped top-level alternation; wrap alternatives as ^(a|b)$"
    return None


def _pattern_actual_end_violation(pattern: object) -> str | None:
    violation = _pattern_anchor_violation(pattern)
    if violation is not None:
        return violation
    assert isinstance(pattern, str)
    body = pattern[1:-1].removesuffix(ACTUAL_END)
    if body not in PATTERN_WITNESSES:
        return f"no valid-body witnesses declared for {body!r}"
    compiled = re.compile(pattern)
    for value in PATTERN_WITNESSES[body]:
        if compiled.search(value) is None:
            return f"rejects valid body witness {value!r}"
        if compiled.search(value + "\n") is not None:
            return f"accepts final LF after valid body witness {value!r}"
    return None


def find_unanchored_schema_patterns(schema_dir: Path = SCHEMA_DIR) -> list[str]:
    entries = _collect_schema_patterns(schema_dir)
    assert entries, f"No pattern declarations found under {schema_dir}"
    violations: list[str] = []
    for _, target, pattern in entries:
        violation = _pattern_actual_end_violation(pattern)
        if violation is not None:
            violations.append(f"{target}: {pattern!r} ({violation})")
    return violations


def assert_schema_patterns_anchored(schema_dir: Path = SCHEMA_DIR) -> None:
    violations = find_unanchored_schema_patterns(schema_dir)
    assert not violations, "Unanchored LinkML schema patterns:\n" + "\n".join(violations)


def _collect_generated_model_patterns(
    model_paths: tuple[Path, ...] = GENERATED_MODEL_PATHS,
) -> list[tuple[str, str, object]]:
    entries: list[tuple[str, str, object]] = []
    for path in model_paths:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in tree.body:
            if not isinstance(node, ast.ClassDef):
                continue
            for item in node.body:
                if not isinstance(item, ast.FunctionDef) or not item.name.startswith("pattern_"):
                    continue
                slot_name = item.name.removeprefix("pattern_")
                compiled_arg: object = None
                for subnode in ast.walk(item):
                    if (
                        isinstance(subnode, ast.Call)
                        and isinstance(subnode.func, ast.Attribute)
                        and isinstance(subnode.func.value, ast.Name)
                        and subnode.func.value.id == "re"
                        and subnode.func.attr == "compile"
                        and subnode.args
                        and isinstance(subnode.args[0], ast.Constant)
                    ):
                        compiled_arg = subnode.args[0].value
                        break
                target = f"{path.name}:{item.lineno}: slot '{node.name}.{slot_name}'"
                entries.append((path.name, target, compiled_arg))
    return entries


def find_unanchored_generated_model_patterns(
    model_paths: tuple[Path, ...] = GENERATED_MODEL_PATHS,
) -> list[str]:
    entries = _collect_generated_model_patterns(model_paths)
    assert entries, "No generated pattern validators found"
    violations: list[str] = []
    for _, target, pattern in entries:
        violation = _pattern_actual_end_violation(pattern)
        if violation is not None:
            violations.append(f"{target}: {pattern!r} ({violation})")
    return violations


def test_all_linkml_schema_patterns_are_anchored_at_both_ends() -> None:
    """Every source pattern accepts its witnesses and rejects each witness plus final LF."""
    assert_schema_patterns_anchored(SCHEMA_DIR)
    entries = _collect_schema_patterns(SCHEMA_DIR)
    assert {name for name, _, _ in entries} == {
        "bsdbng.yaml",
        "core_products.yaml",
        "mtp.yaml",
        "multiple_testing.yaml",
        "workflow.yaml",
    }
    assert len(entries) >= 86


def test_all_generated_model_pattern_validators_are_anchored_at_both_ends() -> None:
    """All 99 generated patterns must enforce actual end under their re.search API."""
    violations = find_unanchored_generated_model_patterns(GENERATED_MODEL_PATHS)
    assert not violations, "Unanchored generated pattern validators:\n" + "\n".join(violations)
    entries = _collect_generated_model_patterns(GENERATED_MODEL_PATHS)
    assert len(entries) >= 99


def test_schema_pattern_anchor_guard_reports_offending_slot_or_type(tmp_path: Path) -> None:
    """#1315: unanchored slots, types, and ungrouped alternations name each offending target."""
    scratch_schema_dir = tmp_path / "schema"
    shutil.copytree(REPO_ROOT / "schema", scratch_schema_dir)

    multiple_testing_path = scratch_schema_dir / "multiple_testing.yaml"
    multiple_testing_schema = yaml.safe_load(multiple_testing_path.read_text(encoding="utf-8"))
    family = multiple_testing_schema["classes"]["HypothesisFamilyMetadata"]["attributes"]
    run = multiple_testing_schema["classes"]["PostArmMultiplicityRun"]["attributes"]
    assert family["strain_rollup"]["pattern"].replace(ACTUAL_END, "") == "^(off|phenotype)$"
    assert run["strain_rollup"]["pattern"].replace(ACTUAL_END, "") == "^(off|phenotype)$"
    assert run["version"]["pattern"].replace(ACTUAL_END, "") == r"^1\.2\.0$"
    family["strain_rollup"]["pattern"] = "(off|phenotype)$"
    run["strain_rollup"]["pattern"] = "^off|phenotype$"
    run["version"]["pattern"] = r"^1\.2\.0"
    multiple_testing_path.write_text(
        yaml.safe_dump(multiple_testing_schema),
        encoding="utf-8",
    )

    mtp_path = scratch_schema_dir / "mtp.yaml"
    mtp_schema = yaml.safe_load(mtp_path.read_text(encoding="utf-8"))
    body_site_slot = mtp_schema["slots"]["body_site_ontology_id"]
    assert body_site_slot["pattern"].replace(ACTUAL_END, "") == "^UBERON:[0-9]+$"
    body_site_slot["pattern"] = "UBERON:[0-9]+$"
    mtp_schema.setdefault("types", {})["curie_code"] = {
        "typeof": "string",
        "pattern": "[A-Z]+:[0-9]+",
    }
    mtp_path.write_text(
        yaml.safe_dump(mtp_schema),
        encoding="utf-8",
    )

    with pytest.raises(AssertionError, match="Unanchored LinkML schema patterns:") as failure:
        assert_schema_patterns_anchored(scratch_schema_dir)
    message = str(failure.value)
    assert (
        "multiple_testing.yaml: slot 'HypothesisFamilyMetadata.strain_rollup' "
        "(classes.HypothesisFamilyMetadata.attributes.strain_rollup): "
        "'(off|phenotype)$' (missing leading '^' anchor)"
    ) in message
    assert (
        "multiple_testing.yaml: slot 'PostArmMultiplicityRun.strain_rollup' "
        "(classes.PostArmMultiplicityRun.attributes.strain_rollup): "
        "'^off|phenotype$' (ungrouped top-level alternation; wrap alternatives as ^(a|b)$)"
    ) in message
    assert (
        "multiple_testing.yaml: slot 'PostArmMultiplicityRun.version' "
        "(classes.PostArmMultiplicityRun.attributes.version): "
        "'^1\\\\.2\\\\.0' (missing trailing '$' anchor)"
    ) in message
    assert (
        "mtp.yaml: slot 'body_site_ontology_id' (slots.body_site_ontology_id): "
        "'UBERON:[0-9]+$' (missing leading '^' anchor)"
    ) in message
    assert (
        "mtp.yaml: type 'curie_code' (types.curie_code): "
        "'[A-Z]+:[0-9]+' (missing leading '^' anchor)"
    ) in message


def test_generated_alternation_parameter_and_experiment_patterns_reject_affixes() -> None:
    """#1315: guard the generated validators where unanchored mutants previously survived."""
    family_required = {
        "family_id": "fam-1",
        "analysis_id": "arm-1",
        "analysis_basis_sha256": DIGEST,
        "source_membership_sha256": DIGEST,
        "retained_feature_count": 6,
        "feature_filter_id": "sha256:" + DIGEST,
        "expected_family_size": 6,
        "population": "human",
        "direction_basis": "curated",
        "profile_unit": "count_profiles",
        "inference_unit": "study",
        "hypothesis_unit": "trait",
        "estimand": "annotation",
        "selection_rule": "all",
        "selection_history": [],
        "sources": [],
    }
    for value in ("off", "phenotype"):
        validated = multiple.HypothesisFamilyMetadata.model_validate(
            {**family_required, "strain_rollup": value}
        )
        assert validated.strain_rollup == value
    for value in (
        "xoff",
        "offx",
        "xphenotype",
        "phenotypex",
        "off|phenotype",
        "\noff",
        "off\n",
        "phenotype\n",
    ):
        with pytest.raises(ValidationError) as failure:
            multiple.HypothesisFamilyMetadata.model_validate(
                {**family_required, "strain_rollup": value}
            )
        assert {error["loc"] for error in failure.value.errors()} == {("strain_rollup",)}

    for value in ("alpha", "Alpha_1"):
        param = workflow.ParameterDefinition.model_validate({"name": value})
        assert param.name == value
    for value in ("1alpha", "_alpha", "alpha-beta", "alpha beta", "\nalpha", "alpha\n"):
        with pytest.raises(ValidationError) as failure:
            workflow.ParameterDefinition.model_validate({"name": value})
        assert {error["loc"] for error in failure.value.errors()} == {("name",)}

    for model in (LEGACY["Experiment"], current.Experiment):
        for value in ("bsdb:10202341-1", "bsdb:PMC1-42"):
            exp = model.model_validate({"id": value, "signatures": []})
            assert exp.model_dump()["id"] == value
        for value in ("xbsdb:1-1", "bsdb:1-1x", "bsdb:1", "bsdb:-1", "\nbsdb:1-1", "bsdb:1-1\n"):
            with pytest.raises(ValidationError) as failure:
                model.model_validate({"id": value, "signatures": []})
            assert {error["loc"] for error in failure.value.errors()} == {("id",)}


@pytest.mark.parametrize(
    ("pattern", "expected_violation"),
    [
        ("^(off|phenotype)$", None),
        ("^(?:off|phenotype)$", None),
        ("^[a|b]+$", None),
        ("^[]|a]+$", None),
        ("^[^]|a]+$", None),
        (r"^a\|b$", None),
        (r"^a\\$", None),
        ("^off|phenotype$", "ungrouped top-level alternation; wrap alternatives as ^(a|b)$"),
        ("^(off)|(phenotype)$", "ungrouped top-level alternation; wrap alternatives as ^(a|b)$"),
        ("^[a]|b$", "ungrouped top-level alternation; wrap alternatives as ^(a|b)$"),
        (r"^\[a]|b$", "ungrouped top-level alternation; wrap alternatives as ^(a|b)$"),
        ("(off|phenotype)$", "missing leading '^' anchor"),
        ("(^a)$", "missing leading '^' anchor"),
        ("^1\\.2\\.0", "missing trailing '$' anchor"),
        (r"^a\$", "missing trailing '$' anchor"),
        (r"^a\\\$", "missing trailing '$' anchor"),
        ("^[a$]", "missing trailing '$' anchor"),
        ("^(a$)", "missing trailing '$' anchor"),
        ("^$", "empty anchored pattern '^$'"),
        ("(?m)^a$", "MULTILINE flag makes ^ and $ match internal line breaks under re.search"),
    ],
)
def test_pattern_anchor_violation_scanner_tracks_escapes_classes_and_groups(
    pattern: str, expected_violation: str | None
) -> None:
    """#1360: scanner tracks escapes, character classes, and group depth without re._parser."""
    assert _pattern_anchor_violation(pattern) == expected_violation


@pytest.mark.parametrize("body", PATTERN_WITNESSES)
def test_actual_end_guard_rejects_bare_dollar_and_accepts_hardened_body(body: str) -> None:
    hardened = "^" + body + ACTUAL_END + "$"
    assert _pattern_actual_end_violation(hardened) is None
    violation = _pattern_actual_end_violation("^" + body + "$")
    assert violation is not None
    assert violation.startswith("accepts final LF after valid body witness ")


def test_actual_end_guard_refuses_unknown_bodies() -> None:
    violation = _pattern_actual_end_violation("^new-unwitnessed-token" + ACTUAL_END + "$")
    assert violation == "no valid-body witnesses declared for 'new-unwitnessed-token'"


def test_actual_end_guard_checks_valid_body_before_rejection(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setitem(PATTERN_WITNESSES, "[0-9a-f]{64}", ("not-a-digest",))
    violation = _pattern_actual_end_violation("^[0-9a-f]{64}" + ACTUAL_END + "$")
    assert violation == "rejects valid body witness 'not-a-digest'"


def test_schema_actual_end_guard_names_bare_dollar_slot(tmp_path: Path) -> None:
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    (schema_dir / "bare.yaml").write_text(
        "slots:\n  digest:\n    pattern: '^[0-9a-f]{64}$'\n", encoding="utf-8"
    )
    violations = find_unanchored_schema_patterns(schema_dir)
    assert len(violations) == 1
    assert "bare.yaml: slot 'digest' (slots.digest)" in violations[0]
    assert "accepts final LF after valid body witness" in violations[0]


def test_generated_actual_end_guard_names_bare_dollar_validator(tmp_path: Path) -> None:
    model_path = tmp_path / "bare.py"
    model_path.write_text(
        "class Source:\n    def pattern_sha256(cls, v):\n"
        "        pattern = re.compile(r'^[0-9a-f]{64}$')\n"
        "        return v\n",
        encoding="utf-8",
    )
    violations = find_unanchored_generated_model_patterns((model_path,))
    assert len(violations) == 1
    assert "bare.py:2: slot 'Source.sha256'" in violations[0]
    assert "accepts final LF after valid body witness" in violations[0]
