"""Compatibility shim forwarding bsdbng.coupling to mtp.coupling."""

import sys
from typing import Any

from mtp.coupling import (
    CANONICAL_LCO_CLUSTER,
    COCCUS_ID,
    DEFAULT_OUTPUT_DIR,
    SPHERE_ID,
    ArmResult,
    DecisionRuleResult,
    LCOSummaryRow,
    LOTORow,
    build_coupling_report_text,
    compute_composition_arm,
    compute_directional_arm,
    compute_ora_arm,
    evaluate_decision_rule,
    evaluate_lco,
    evaluate_loto,
    run_coupling_audit,
)

__all__ = [
    "CANONICAL_LCO_CLUSTER",
    "COCCUS_ID",
    "DEFAULT_OUTPUT_DIR",
    "SPHERE_ID",
    "ArmResult",
    "DecisionRuleResult",
    "LCOSummaryRow",
    "LOTORow",
    "build_coupling_report_text",
    "compute_composition_arm",
    "compute_directional_arm",
    "compute_ora_arm",
    "evaluate_decision_rule",
    "evaluate_lco",
    "evaluate_loto",
    "run_coupling_audit",
]

import mtp.coupling as _target

sys.modules[__name__] = _target


def __getattr__(name: str) -> Any:
    return getattr(_target, name)
