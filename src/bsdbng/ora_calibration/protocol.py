"""Statistical calibration protocol and interfaces for multi-method trait evaluation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Protocol

import pandas as pd


@dataclass
class MethodInferenceResult:
    """Standardized output from a method adapter."""

    method_id: str
    owner: str
    target_null: str
    results_df: pd.DataFrame  # profile_id, trait, p_value, q_value, is_significant, effect_size
    metadata: dict[str, Any]


class MethodAdapter(Protocol):
    """Protocol that each framework adapter must satisfy."""

    method_id: str
    owner: str
    target_null: str

    def run_inference(
        self,
        profiles_df: pd.DataFrame,
        taxa_matrix_df: pd.DataFrame,
        config: dict[str, Any] | None = None,
    ) -> MethodInferenceResult:
        """Run the end-to-end inference procedure on provided inputs."""
        ...
