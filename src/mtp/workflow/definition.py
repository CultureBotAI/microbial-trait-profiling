"""Schema and graph validation for tracked workflow definitions (#629)."""

from __future__ import annotations

from pathlib import Path, PurePosixPath

import yaml

from mtp.workflow.models import WorkflowDefinition


def unique(values: list[str], label: str) -> None:
    """Reject duplicates rather than silently collapsing identifiers."""
    if len(values) != len(set(values)):
        raise ValueError(f"Duplicate {label}")


def safe_relative_path(value: str) -> None:
    """Only normalized portable relative paths can name published products."""
    path = PurePosixPath(value)
    if (
        not value
        or path.is_absolute()
        or "\\" in value
        or any(part in {"", ".", ".."} for part in value.split("/"))
        or value != path.as_posix()
    ):
        raise ValueError(f"Unsafe relative output path: {value!r}")


def validate_definition(definition: WorkflowDefinition) -> list[str]:
    """Validate references and return a deterministic topological step order."""
    unique([item.id for item in definition.products], "product identifiers")
    unique([item.id for item in definition.steps], "step identifiers")
    unique([item.name for item in definition.parameters or []], "parameter names")
    products = {item.id: item for item in definition.products}
    steps = {item.id: item for item in definition.steps}
    parameters = {item.name for item in definition.parameters or []}
    if not products or not steps:
        raise ValueError("A workflow needs products and steps")
    producers: dict[str, str] = {}
    for product in products.values():
        safe_relative_path(product.filename)
        if PurePosixPath(product.filename).parts[0] == "_step.json":
            raise ValueError("Output path is reserved for the step publication record")
        if product.cohort is not None and product.cohort != product.cohort.strip():
            raise ValueError(f"Whitespace-padded cohort contract for {product.id}")
        if product.cohort_parameter and product.cohort_parameter not in parameters:
            raise ValueError(f"Unknown cohort parameter for {product.id}")
    for step in steps.values():
        for label, references in (
            ("inputs", step.inputs),
            ("outputs", step.outputs),
            ("dependencies", step.depends_on or []),
            ("parameters", step.parameters or []),
        ):
            unique(references, f"{step.id} {label}")
        if not step.outputs:
            raise ValueError(f"Step {step.id} declares no outputs")
        if set(step.inputs) & set(step.outputs):
            raise ValueError(f"Step {step.id} overwrites an input")
        for product_id in [*step.inputs, *step.outputs]:
            if product_id not in products:
                raise ValueError(f"Unknown product {product_id} in step {step.id}")
        if not set(step.parameters or []).issubset(parameters):
            raise ValueError(f"Unknown parameter in step {step.id}")
        if not set(step.depends_on or []).issubset(steps):
            raise ValueError(f"Unknown dependency in step {step.id}")
        filenames: list[PurePosixPath] = []
        for product_id in step.outputs:
            product = products[product_id]
            if product.external or product_id in producers:
                raise ValueError(f"Duplicate or external output product {product_id}")
            if product.cohort_parameter and product.cohort_parameter not in (step.parameters or []):
                raise ValueError(
                    f"Producer {step.id} must resolve cohort parameter for {product_id}"
                )
            filename = PurePosixPath(product.filename)
            if any(
                filename == prior or prior in filename.parents or filename in prior.parents
                for prior in filenames
            ):
                raise ValueError(f"Overlapping output paths in step {step.id}")
            filenames.append(filename)
            producers[product_id] = step.id
    unused = (
        set(products)
        - set(producers)
        - {key for step in steps.values() for key in step.inputs if products[key].external}
    )
    if unused:
        raise ValueError(f"Products have no producer or use: {sorted(unused)}")
    used_parameters = {name for step in steps.values() for name in step.parameters or []}
    if parameters - used_parameters:
        raise ValueError(f"Unused declared parameters: {sorted(parameters - used_parameters)}")
    dependencies: dict[str, set[str]] = {}
    for step in steps.values():
        dependencies[step.id] = set(step.depends_on or [])
        for product_id in step.inputs:
            if not products[product_id].external:
                if product_id not in producers:
                    raise ValueError(f"No producer for input {product_id}")
                dependencies[step.id].add(producers[product_id])
    ordered: list[str] = []
    while len(ordered) < len(steps):
        ready = sorted(
            key
            for key, deps in dependencies.items()
            if key not in ordered and deps.issubset(ordered)
        )
        if not ready:
            raise ValueError("Workflow dependency cycle")
        ordered.extend(ready)
    return ordered


def read_definition(path: str | Path) -> WorkflowDefinition:
    """Read a schema-valid definition and reject structurally valid bad graphs."""
    value = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    definition = WorkflowDefinition.model_validate(value)
    validate_definition(definition)
    return definition
