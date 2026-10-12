"""Validate all per-study YAML files against the Pydantic model."""

from __future__ import annotations

import sys
from pathlib import Path

import yaml

from mtp.datamodel import Study
from mtp.paths import STUDY_DATA_DIR


def validate_all_studies(study_dir: Path | None = None) -> None:
    """Validate every YAML file in study_dir through the Study Pydantic model.

    Prints pass/fail counts to stderr. Exits non-zero if any fail.
    """
    study_dir = study_dir or STUDY_DATA_DIR
    files = sorted(study_dir.glob("*.yaml"))

    if not files:
        print("no study YAML files found", file=sys.stderr)
        sys.exit(1)

    passed = 0
    failed = 0
    for path in files:
        try:
            data = yaml.safe_load(path.read_text())
            Study.model_validate(data)
            passed += 1
        except Exception as exc:
            print(f"FAIL: {path.name}: {exc}", file=sys.stderr)
            failed += 1

    print(f"validated {passed + failed} files: {passed} passed, {failed} failed", file=sys.stderr)
    if failed:
        sys.exit(1)
