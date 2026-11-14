"""Probe real wrapper diagnostics with synthetic queries (#1185/#1163, 2026-11-08).

Inputs: wrapper, pinned synthetic requirements, cache/trace paths, and scenario.
Outputs: real entrypoint streams/exit plus a JSON attempt and sleep trace.
Only PyPIService.query and time.sleep are replaced; no vulnerability feed is read.
"""

from __future__ import annotations

import argparse
import importlib
import json
import logging
import runpy
import sys
import time
from pathlib import Path
from typing import Any
from unittest.mock import patch


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wrapper", required=True, type=Path)
    parser.add_argument("--requirements", required=True, type=Path)
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--trace", required=True, type=Path)
    parser.add_argument(
        "--scenario", required=True, choices=("recover", "exhaust", "exit-message", "usage")
    )
    args = parser.parse_args()
    if "pip_audit._cli" in sys.modules:
        raise RuntimeError("Probe requires a fresh process with pip_audit._cli unimported")

    requests: Any = importlib.import_module("requests")
    service: Any = importlib.import_module("pip_audit._service.pypi")
    attempts: list[int] = []
    sleeps: list[float] = []
    original_argv = sys.argv

    def query(_self: Any, spec: Any) -> tuple[Any, list[Any]]:
        attempt = len(attempts) + 1
        attempts.append(attempt)
        logging.getLogger().warning("PROBE root attempt=%s", attempt)
        print(f"PROBE stdout-before attempt={attempt}", flush=True)
        logging.getLogger("urllib3.issue1185").warning("PROBE urllib3 attempt=%s", attempt)
        print(f"PROBE stderr attempt={attempt}", file=sys.stderr, flush=True)
        logging.getLogger("pip_audit.issue1185").warning("PROBE pip-audit attempt=%s", attempt)
        print(f"PROBE stdout-after attempt={attempt}", flush=True)
        if args.scenario == "exit-message":
            raise SystemExit("PROBE string-valued exit message")
        if args.scenario == "exhaust" or attempt == 1:
            raise requests.exceptions.ConnectionError("PROBE synthetic disconnect")
        return spec, []

    def sleep(delay: float) -> None:
        sleeps.append(delay)
        print(f"PROBE sleep={delay:g}", file=sys.stderr, flush=True)

    sys.argv = [str(args.wrapper)]
    if args.scenario == "usage":
        sys.argv.append("--issue1185-invalid-option")
    else:
        sys.argv.extend(
            [
                "-r",
                str(args.requirements),
                "--disable-pip",
                "--no-deps",
                "--progress-spinner",
                "off",
                "--cache-dir",
                str(args.cache_dir),
            ]
        )
    try:
        with patch.object(service.PyPIService, "query", query), patch.object(time, "sleep", sleep):
            runpy.run_path(str(args.wrapper), run_name="__main__")
    finally:
        sys.argv = original_argv
        with args.trace.open("x", encoding="utf-8") as handle:
            json.dump(
                {"cli_preloaded": False, "attempts": attempts, "sleeps": sleeps},
                handle,
                indent=2,
            )
            handle.write("\n")


if __name__ == "__main__":
    main()
