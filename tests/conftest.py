"""Suite-wide guards."""

import importlib.util
from collections.abc import Iterator
from pathlib import Path

import pytest

from bsdbng.paths import REPO_ROOT

# The cellulose_decr figure tests import matplotlib, from the viz extra. Without it they
# are left out rather than failing at collection (#631).
collect_ignore = (
    [] if importlib.util.find_spec("matplotlib") else ["test_cellulose_decr_figures.py"]
)

# Every real pipeline input and output: raw downloads and their provenance,
# study YAML, derived tables, and the KG-Microbe snapshot.
GUARDED_DATA_DIR = REPO_ROOT / "data"


def _snapshot(root: Path) -> dict[str, tuple[int, int]]:
    if not root.is_dir():
        return {}
    return {
        str(path.relative_to(root)): (path.stat().st_mtime_ns, path.stat().st_size)
        for path in root.rglob("*")
        if path.is_file()
    }


@pytest.fixture(autouse=True, scope="session")
def _tests_never_write_real_data() -> Iterator[None]:
    """Fail the run if any test creates, modifies or deletes a file under data/ (#417, #421).

    A test that falls back to a default path writes over real pipeline data:
    ingest() once wrote its fixture log over data/derived/ingest_log.json, and
    ingest() and download_exports() also default to data/studies and data/raw.
    """
    before = _snapshot(GUARDED_DATA_DIR)
    yield
    after = _snapshot(GUARDED_DATA_DIR)
    changed = sorted(
        name for name in before.keys() | after.keys() if before.get(name) != after.get(name)
    )
    assert not changed, f"tests changed real files under {GUARDED_DATA_DIR}: {', '.join(changed)}"
