"""Scratch-copy mutant verifier for #1148 and #1146 (2026-11-08).

Purpose:
    Verify that tests/test_historical_trait_policy_core.py passes on a scratch
    copy of scripts/adhoc/issue1134_generate_historical_core_fixture.py and
    kills mutants that reintroduce cast(Validator, ...), pass wrong historical
    policy arguments in adapt_historical_validator, or fail to adapt the
    pre-union 511c891 validator signature.
Inputs:
    - scripts/adhoc/issue1134_generate_historical_core_fixture.py
    - tests/test_historical_trait_policy_core.py
Outputs:
    - Exit 0 only when the unmutated scratch copy passes and every mutant fails
      at the expected test assertion with PYTEST_ORIGIN_OK.
Issues:
    - #1148, #1146
"""

from __future__ import annotations

import importlib.util
import multiprocessing
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
HELPER_REL = Path("scripts/adhoc/issue1134_generate_historical_core_fixture.py")
TEST_REL = Path("tests/test_historical_trait_policy_core.py")


@dataclass(frozen=True)
class Mutant:
    name: str
    target: str
    replacement: str
    expected_substring: str


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        name="M1_reintroduce_cast",
        target="validator=adapt_historical_validator(validate_core_product),",
        replacement="validator=cast(Validator, validate_core_product),",
        expected_substring="cast(",
    ),
    Mutant(
        name="M2_wrong_trait_key_policy",
        target='trait_key_policy="raw_object_v1",',
        replacement='trait_key_policy="facultative_union_v1",',
        expected_substring="Feature metadata",
    ),
    Mutant(
        name="M3_wrong_include_policy_metadata",
        target="include_policy_metadata=False,",
        replacement="include_policy_metadata=True,",
        expected_substring="missing columns",
    ),
    Mutant(
        name="M4_unconditional_policy_kwargs_on_legacy",
        target="if accepts_policy:",
        replacement="if True:",
        expected_substring="unexpected keyword argument",
    ),
    Mutant(
        name="M5_bypass_validator_adapter",
        target="validator=adapt_historical_validator(validate_core_product),",
        replacement="validator=validate_core_product,",
        expected_substring="validator=adapt_historical_validator(validate_core_product)",
    ),
)


def _stage_scratch(scratch_root: Path) -> Path:
    helper_dst = scratch_root / HELPER_REL
    helper_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(REPO_ROOT / HELPER_REL, helper_dst)

    test_dst = scratch_root / TEST_REL
    test_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(REPO_ROOT / TEST_REL, test_dst)

    fixtures_dst = scratch_root / "tests/fixtures/trait_key_policy"
    fixtures_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(REPO_ROOT / "tests/fixtures/trait_key_policy", fixtures_dst)
    return helper_dst


def _child_worker(scratch_root_str: str, queue: Any) -> None:
    scratch_root = Path(scratch_root_str).resolve()
    sys.path.insert(0, str(REPO_ROOT / "src"))
    test_file = scratch_root / TEST_REL
    spec = importlib.util.spec_from_file_location("scratch_test_historical_1146", test_file)
    if spec is None or spec.loader is None:
        queue.put(("ERR", "Could not load scratch test module"))
        return
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    origin = Path(getattr(mod, "__file__", "")).resolve()
    if origin != test_file or not str(origin).startswith(str(scratch_root)):
        queue.put(("ERR", f"Origin mismatch: {origin} != {test_file}"))
        return

    tmp_dir = scratch_root / "tmp_case"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    try:
        mod.test_real_verifier_accepts_relocated_policyless_raw_core_without_changing_bytes(
            tmp_dir / "t1"
        )
        mod.test_historical_fixture_generator_adapter_and_decoupled_helper_receipt(tmp_dir / "t2")
    except Exception as exc:
        queue.put(("FAIL", f"{type(exc).__name__}: {exc}"))
        return
    queue.put(("PASS", f"PYTEST_ORIGIN_OK:{origin}"))


def _run_in_spawn(scratch_root: Path) -> tuple[str, str]:
    ctx = multiprocessing.get_context("spawn")
    queue = ctx.Queue()
    proc = ctx.Process(target=_child_worker, args=(str(scratch_root), queue))
    proc.start()
    proc.join(timeout=60)
    if proc.is_alive():
        proc.terminate()
        proc.join()
        return ("TIMEOUT", "Child process timed out")
    if queue.empty():
        return ("CRASH", f"Child exited with code {proc.exitcode}")
    status, detail = queue.get()
    return (str(status), str(detail))


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="verify_1146_mutants_") as tmp_str:
        scratch_root = Path(tmp_str).resolve()
        helper_dst = _stage_scratch(scratch_root)
        original = helper_dst.read_text(encoding="utf-8")

        status, detail = _run_in_spawn(scratch_root)
        if status != "PASS" or not detail.startswith("PYTEST_ORIGIN_OK:"):
            print(f"Baseline unmutated scratch run failed: {status} {detail}", file=sys.stderr)
            return 1
        print(f"Baseline unmutated scratch run passed ({detail})")

        killed = 0
        for mutant in MUTANTS:
            if mutant.target not in original:
                print(
                    f"Mutant {mutant.name} target not found in {HELPER_REL}",
                    file=sys.stderr,
                )
                return 1
            mutated = original.replace(mutant.target, mutant.replacement, 1)
            helper_dst.write_text(mutated, encoding="utf-8")
            try:
                m_status, m_detail = _run_in_spawn(scratch_root)
            finally:
                helper_dst.write_text(original, encoding="utf-8")

            if m_status != "FAIL" or mutant.expected_substring not in m_detail:
                print(
                    f"Mutant {mutant.name} survived or failed unexpectedly: "
                    f"status={m_status}, detail={m_detail!r}, "
                    f"expected={mutant.expected_substring!r}",
                    file=sys.stderr,
                )
                return 1
            killed += 1
            print(f"Killed {mutant.name}: {m_detail}")

        print(f"Summary: {killed}/{len(MUTANTS)} mutants killed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
