"""Reject non-finite step-up inputs while preserving finite BH/BY results (#1032)."""

import numpy as np
import pytest

from mtp.trait_screen import step_up, step_up_array


@pytest.mark.parametrize(
    ("api", "method", "bad", "position"),
    [
        pytest.param(api, method, bad, position, id=f"{api}-{method}-{label}-{position}")
        for api in ("array", "list")
        for method in ("bh", "by")
        for label, bad in (("nan", np.nan), ("inf", np.inf), ("negative-inf", -np.inf))
        for position in ("start", "middle", "end", "single")
    ],
)
def test_nonfinite(api: str, method: str, bad: float, position: str) -> None:
    values = {
        "start": [bad, 0.01, 0.03],
        "middle": [0.01, bad, 0.03],
        "end": [0.01, 0.03, bad],
        "single": [bad],
    }[position]
    if api == "array":
        with pytest.raises(ValueError, match=r"^p-values must be finite$"):
            step_up_array(np.asarray(values), method=method, family_size=3)
    else:
        with pytest.raises(ValueError, match=r"^p-values must be finite$"):
            step_up(values, method=method, family_size=3)


@pytest.mark.parametrize("api", ["array", "list"])
@pytest.mark.parametrize("method", ["bh", "by"])
@pytest.mark.parametrize("values", [[], [0.01, 1.0, 0.03], [0.0, 0.5, 1.0]])
def test_finite_values(api: str, method: str, values: list[float]) -> None:
    if not values:
        expected = []
    elif values[0] == 0.01:
        expected = [0.03, 1.0, 0.045] if method == "bh" else [0.055, 1.0, 0.0825]
    else:
        expected = [0.0, 0.75, 1.0] if method == "bh" else [0.0, 1.0, 1.0]
    result = (
        step_up_array(np.asarray(values), method=method, family_size=3)
        if api == "array"
        else step_up(values, method=method, family_size=3)
    )
    assert np.all(np.isfinite(result))
    np.testing.assert_allclose(result, expected, rtol=1e-14, atol=0.0)


def test_finite_contract_documented() -> None:
    doc = " ".join((step_up_array.__doc__ or "").split())
    assert "finite 1D p-values" in doc
    assert "non-finite values are rejected" in doc
