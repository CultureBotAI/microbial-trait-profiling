"""Reject supplied-root layout errors before workflow writes (#1817).

These bounded routing tests stub scientific verification, but use real portable
path checks. The comparison suite separately verifies complete runs and retries.
"""

from __future__ import annotations

import importlib
import importlib.util
import re
from dataclasses import dataclass, replace
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

from mtp import integrated_evidence, provenance_paths
from mtp.paths import file_sha256
from mtp.workflow import comparison
from mtp.workflow.arms import runtime
from mtp.workflow.arms.types import ArmContext, Command


@dataclass(frozen=True)
class Inputs:
    root: Path
    manifests: dict[str, Path]
    pooled: Path
    products: dict[str, Path]
    evidence: dict[tuple[str, str, str], Path]

    def known_sources(self) -> tuple[Path, ...]:
        return (
            *self.manifests.values(),
            self.pooled,
            *self.products.values(),
            *self.evidence.values(),
        )


@pytest.fixture
def inputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Inputs:
    root = tmp_path / "source-bundle"
    manifests = {name: root / f"{name}.json" for name in ("core", "noel", "jed", "ilyas")}
    pooled = root / "pooled.tsv"
    products = {
        "count_profiles": root / "core-products/count.tsv",
        "ora_profiles": root / "core-products/ora.tsv",
    }
    keys = {
        *((arm, step, name) for _, arm, step, name in comparison.ARM_EVIDENCE_INPUTS),
        *comparison.CALIBRATION_SOURCES,
    }
    evidence = {(arm, step, name): root / arm / step / name for arm, step, name in keys}
    value = Inputs(root, manifests, pooled, products, evidence)
    for path in value.known_sources():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("synthetic routing input\n", encoding="utf-8")
    monkeypatch.setenv("TMPDIR", str(tmp_path))
    core = SimpleNamespace(
        run_id="synthetic-core",
        contrast_population="all_paired",
        direction_basis="group_1_vs_group_0",
        parameters=[],
        artifacts=[SimpleNamespace(product_id=name) for name in products],
    )
    arms = {
        path: SimpleNamespace(
            arm=name,
            core_run_id=core.run_id,
            core_manifest_sha256=file_sha256(manifests["core"]),
            contrast_population=core.contrast_population,
            direction_basis=core.direction_basis,
        )
        for name, path in manifests.items()
        if name != "core"
    }
    monkeypatch.setattr(comparison, "verify_core_run", lambda _path: core)
    monkeypatch.setattr(comparison, "verify_arm", lambda path: arms[path])
    monkeypatch.setattr(comparison, "artifact_path", lambda _run, product: products[product])
    monkeypatch.setattr(
        comparison, "recorded_output", lambda arm, step, name: evidence[(arm.arm, step, name)]
    )
    for name in (
        "noel_calibration_verdict",
        "jed_calibration_verdict",
        "shared_basis_ledger",
        "core_analysis_bases",
    ):
        monkeypatch.setattr(comparison, name, lambda *_args: None)
    bundle = SimpleNamespace(manifest=manifests["core"], run=core, path=products.__getitem__)
    monkeypatch.setattr(runtime, "open_core_bundle", lambda *_args, **_kwargs: bundle)
    monkeypatch.setattr(runtime, "preflight_consumer_outputs", lambda *_args, **_kwargs: None)
    for module in (comparison, runtime):
        monkeypatch.setattr(module, "repository_state", lambda _root: ("a" * 40, False))
    return value


def bad_root(kind: str, inputs: Inputs, boundary: Path) -> Path:
    if kind == "missing":
        return boundary / "missing-root"
    if kind == "file":
        return inputs.pooled
    if kind == "home":
        return Path.home()
    if kind == "temporary-boundary":
        return boundary
    if kind == "filesystem-root":
        return Path(boundary.anchor)
    assert kind == "excludes-sources"
    excluded = inputs.root / "unrelated-directory"
    excluded.mkdir()
    return excluded


def producer_refusal(record: Path, sources: tuple[Path, ...], root: Path) -> str:
    with pytest.raises(
        ValueError,
        match=r"^(?:--provenance-root|Provenance root|The filesystem root|Source is outside)",
    ) as stopped:
        provenance_paths.PortablePaths.create(record, sources, root=root)
    return str(stopped.value)


def invoke_comparison(inputs: Inputs, output: Path, root: Path, *, generated: bool = False) -> None:
    comparison.execute_comparison_run(
        inputs.manifests["core"],
        inputs.manifests["noel"],
        inputs.manifests["jed"],
        inputs.manifests["ilyas"],
        output,
        arm1_pooled_stats=None if generated else inputs.pooled,
        provenance_root=root,
    )


class ReachedNextStage(Exception):
    """A valid retry reached the real wrapper's next stage, without scientific work."""


@pytest.mark.parametrize("workflow", ["comparison", "noel"])
@pytest.mark.parametrize(
    "kind",
    ["missing", "file", "home", "temporary-boundary", "filesystem-root", "excludes-sources"],
)
def test_invalid_supplied_root_leaves_no_output_and_allows_a_retry(
    inputs: Inputs,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    workflow: str,
    kind: str,
) -> None:
    output = (
        tmp_path / "comparison-output" if workflow == "comparison" else inputs.root / "noel-output"
    )
    supplied = bad_root(kind, inputs, tmp_path)
    record = output / "integrated_evidence/integrated_evidence_provenance.json"
    expected = producer_refusal(record, inputs.known_sources(), supplied)
    allowed = False
    reached: list[Path | None] = []
    processes: list[str] = []

    def planner(context: ArmContext) -> list[Command]:
        assert allowed, "An invalid supplied root reached the Noel planner"
        reached.append(context.provenance_root)
        raise ReachedNextStage

    def integrated(arguments: list[str]) -> int:
        assert allowed, "An invalid supplied root reached the integrated producer"
        index = arguments.index("--provenance-root")
        reached.append(Path(arguments[index + 1]))
        raise ReachedNextStage

    async def forbidden_process(*_args: object) -> int:
        processes.append("process")
        raise AssertionError("The root-layout test must not launch any subprocess")

    monkeypatch.setattr(integrated_evidence, "main", integrated)
    monkeypatch.setattr(comparison, "_exec_command", forbidden_process)
    monkeypatch.setattr(runtime, "_run_process", forbidden_process)

    def invoke(root: Path) -> None:
        if workflow == "comparison":
            invoke_comparison(inputs, output, root)
        else:
            runtime.execute_arm(
                "noel", inputs.manifests["core"], output, planner, provenance_root=root
            )

    with pytest.raises(ValueError, match=f"^{re.escape(expected)}$") as stopped:
        invoke(supplied)
    assert str(stopped.value) == expected
    assert not output.exists()
    assert reached == []
    assert processes == []

    allowed = True
    try:
        with pytest.raises(ReachedNextStage):
            invoke(inputs.root)
    except (ValueError, OSError) as error:
        raise AssertionError("The corrected root must allow retrying the same output") from error
    assert reached == [inputs.root.resolve()]
    assert processes == []
    if workflow == "comparison":
        assert list(output.iterdir()) == []
    else:
        assert not output.exists()


@pytest.mark.parametrize("kind", ["missing", "excludes-generated-source"])
def test_generated_pooled_source_is_checked_before_creating_or_running_anything(
    inputs: Inputs, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kind: str
) -> None:
    output = tmp_path / "outside-comparison"
    supplied = tmp_path / "missing-root" if kind == "missing" else inputs.root
    planned = output / "arm1_pooled_stats/trait_profile_paired_contrast_trait_stats_pooled.tsv"
    record = output / "integrated_evidence/integrated_evidence_provenance.json"
    expected = producer_refusal(record, (*inputs.known_sources(), planned), supplied)
    calls: list[str] = []

    async def forbidden_generator(*_args: object) -> int:
        calls.append("generator")
        raise AssertionError("Invalid root must be refused before generating pooled statistics")

    monkeypatch.setattr(comparison, "_exec_command", forbidden_generator)
    with pytest.raises(ValueError, match=f"^{re.escape(expected)}$") as stopped:
        invoke_comparison(inputs, output, supplied, generated=True)
    assert str(stopped.value) == expected
    assert calls == []
    assert not output.exists()
    assert not planned.exists()


@pytest.mark.parametrize("category", ["core-product", "arm-evidence", "manual-pooled"])
def test_comparison_checks_each_source_category_before_writing(
    inputs: Inputs, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, category: str
) -> None:
    outside = tmp_path / "external" / "single-source.tsv"
    outside.parent.mkdir()
    outside.write_text("synthetic external source\n", encoding="utf-8")
    if category == "core-product":
        inputs.products["count_profiles"] = outside
    elif category == "arm-evidence":
        _, arm, step, name = comparison.ARM_EVIDENCE_INPUTS[0]
        inputs.evidence[(arm, step, name)] = outside
    else:
        inputs = replace(inputs, pooled=outside)
    output = tmp_path / "comparison-output"
    record = output / "integrated_evidence/integrated_evidence_provenance.json"
    expected = producer_refusal(record, inputs.known_sources(), inputs.root)
    reached: list[str] = []

    def forbidden_integrated(_arguments: list[str]) -> int:
        reached.append("integrated")
        raise AssertionError("An outside source must be rejected before the producer")

    async def forbidden_generator(*_args: object) -> int:
        reached.append("generator")
        raise AssertionError("An outside source must be rejected before a subprocess")

    monkeypatch.setattr(integrated_evidence, "main", forbidden_integrated)
    monkeypatch.setattr(comparison, "_exec_command", forbidden_generator)
    with pytest.raises(ValueError, match=f"^{re.escape(expected)}$") as stopped:
        invoke_comparison(inputs, output, inputs.root)
    assert str(stopped.value) == expected
    assert reached == []
    assert not output.exists()


@pytest.mark.parametrize("category", ["output", "core-product"])
def test_noel_checks_output_and_core_products_before_planning(
    inputs: Inputs, tmp_path: Path, category: str
) -> None:
    output = inputs.root / "noel-output"
    if category == "output":
        output = tmp_path / "outside-noel-output"
    else:
        outside = tmp_path / "external" / "single-core-product.tsv"
        outside.parent.mkdir()
        outside.write_text("synthetic external core product\n", encoding="utf-8")
        inputs.products["count_profiles"] = outside
    record = output / "integrated_evidence/integrated_evidence_provenance.json"
    expected = producer_refusal(record, (*inputs.known_sources(), output), inputs.root)
    reached: list[str] = []

    def forbidden_planner(_context: ArmContext) -> list[Command]:
        reached.append("planner")
        raise AssertionError("An outside Noel source or output must be rejected before planning")

    with pytest.raises(ValueError, match=f"^{re.escape(expected)}$") as stopped:
        runtime.execute_arm(
            "noel",
            inputs.manifests["core"],
            output,
            forbidden_planner,
            provenance_root=inputs.root,
        )
    assert str(stopped.value) == expected
    assert reached == []
    assert not output.exists()


def preflight_module() -> ModuleType:
    assert importlib.util.find_spec("mtp.workflow.provenance") is not None, (
        "The workflows need one shared provenance-root preflight"
    )
    module = importlib.import_module("mtp.workflow.provenance")
    assert callable(getattr(module, "preflight_provenance_root", None))
    assert module.__file__ is not None
    assert Path(module.__file__).resolve() == (
        Path(__file__).resolve().parents[1] / "src/mtp/workflow/provenance.py"
    )
    return module


@pytest.mark.parametrize("account", ["tsv", "json"])
def test_future_file_extension_is_not_mistaken_for_the_account_name(
    inputs: Inputs, monkeypatch: pytest.MonkeyPatch, account: str
) -> None:
    module = preflight_module()
    monkeypatch.setattr(provenance_paths, "_account_names", lambda: {account})
    planned = inputs.root / "future" / ("statistics." + account)
    record = inputs.root / "report/integrated_evidence_provenance.json"
    try:
        module.preflight_provenance_root(
            record, (inputs.manifests["core"],), root=inputs.root, planned_files=(planned,)
        )
    except (ValueError, OSError) as error:
        raise AssertionError(
            "A declared future file must keep its format-suffix exemption"
        ) from error
    assert not planned.parent.exists()
    assert not record.parent.exists()


@pytest.mark.parametrize("account", ["tsv", "json"])
@pytest.mark.parametrize("location", ["directory", "stem"])
def test_future_file_context_does_not_exempt_account_tokens_elsewhere(
    inputs: Inputs, monkeypatch: pytest.MonkeyPatch, account: str, location: str
) -> None:
    module = preflight_module()
    monkeypatch.setattr(provenance_paths, "_account_names", lambda: {account})
    planned = (
        inputs.root / ("chunk." + account) / "statistics.tsv"
        if location == "directory"
        else inputs.root / "future" / (account + ".tsv")
    )
    record = inputs.root / "report/integrated_evidence_provenance.json"
    with pytest.raises(ValueError, match="Provenance contains a local account name"):
        module.preflight_provenance_root(
            record, (inputs.manifests["core"],), root=inputs.root, planned_files=(planned,)
        )
    assert not planned.parent.exists()
    assert not record.parent.exists()
