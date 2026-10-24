"""Display labels for the workflow's declared profile direction basis (#962).

These describe the supplied data. They do not infer, verify, or change its orientation.
Keep this hand-maintained display map outside the generated workflow datamodel.
"""

from collections.abc import Mapping
from types import MappingProxyType

from mtp.datamodel.workflow_pydantic import DirectionBasis

DIRECTION_BASIS_LABELS: Mapping[str, str] = MappingProxyType(
    {
        DirectionBasis.group_1_vs_group_0.value: "Group 1 vs Group 0",
        DirectionBasis.disease_vs_comparison.value: "disease vs comparison",
    }
)
