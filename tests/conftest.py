"""Suite-wide guards."""

from collections.abc import Iterator
from pathlib import Path

import pytest

from bsdbng.paths import DERIVED_DATA_DIR


def _snapshot(root: Path) -> dict[str, tuple[int, int]]:
    if not root.is_dir():
        return {}
    return {
        str(path.relative_to(root)): (path.stat().st_mtime_ns, path.stat().st_size)
        for path in root.rglob("*")
        if path.is_file()
    }


@pytest.fixture(autouse=True, scope="session")
def _tests_never_write_real_derived_data() -> Iterator[None]:
    """Fail the run if any test creates or modifies a file under data/derived (#417).

    A test that falls back to a default output path overwrites real pipeline
    outputs: ingest() once wrote its fixture log over data/derived/ingest_log.json.
    """
    before = _snapshot(DERIVED_DATA_DIR)
    yield
    after = _snapshot(DERIVED_DATA_DIR)
    changed = sorted(name for name in after if before.get(name) != after[name])
    assert not changed, f"tests wrote to real {DERIVED_DATA_DIR}: {', '.join(changed)}"
