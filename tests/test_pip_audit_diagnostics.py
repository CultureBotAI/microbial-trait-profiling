"""Fresh-process, real-entrypoint diagnostic witnesses for #1163 and #1185."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PROBE = REPO_ROOT / "scripts" / "adhoc" / "issue1185_pip_audit_probe.py"
WRAPPER = REPO_ROOT / "scripts" / "run_pip_audit.py"


def _run_probe(
    tmp_path: Path, scenario: str
) -> tuple[subprocess.CompletedProcess[str], dict[str, Any]]:
    requirements = tmp_path / "requirements.txt"
    requirements.write_text("demo-pkg==1.0.0\n", encoding="utf-8")
    trace_path = tmp_path / "trace.json"
    completed = subprocess.run(
        [
            sys.executable,
            "-I",
            "-B",
            "-u",
            str(PROBE),
            "--wrapper",
            str(WRAPPER),
            "--requirements",
            str(requirements),
            "--cache-dir",
            str(tmp_path / "cache"),
            "--trace",
            str(trace_path),
            "--scenario",
            scenario,
        ],
        cwd=tmp_path,
        env={**os.environ, "PIP_AUDIT_LOGLEVEL": "INFO"},
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=30,
    )
    assert trace_path.is_file(), completed.stdout
    trace = json.loads(trace_path.read_text(encoding="utf-8"))
    assert trace["cli_preloaded"] is False
    return completed, trace


def _attempt_markers(attempt: int) -> list[str]:
    return [
        f"WARNING:root:PROBE root attempt={attempt}",
        f"PROBE stdout-before attempt={attempt}",
        f"WARNING:urllib3.issue1185:PROBE urllib3 attempt={attempt}",
        f"PROBE stderr attempt={attempt}",
        f"WARNING:pip_audit.issue1185:PROBE pip-audit attempt={attempt}",
        f"PROBE stdout-after attempt={attempt}",
    ]


def _assert_ordered_once(output: str, markers: list[str]) -> None:
    positions = []
    for marker in markers:
        assert output.count(marker) == 1, (marker, output)
        positions.append(output.index(marker))
    assert positions == sorted(positions), output


@pytest.mark.parametrize(
    ("scenario", "attempt_count", "exit_code", "sleeps"),
    [("recover", 2, 0, [2.0]), ("exhaust", 4, 3, [2.0, 4.0, 8.0])],
)
def test_real_entrypoint_preserves_every_attempt_stream_and_logger_in_order(
    tmp_path: Path, scenario: str, attempt_count: int, exit_code: int, sleeps: list[float]
) -> None:
    completed, trace = _run_probe(tmp_path, scenario)
    assert completed.returncode == exit_code, completed.stdout
    assert trace["attempts"] == list(range(1, attempt_count + 1))
    assert trace["sleeps"] == sleeps
    markers = []
    for attempt in range(1, attempt_count + 1):
        markers.extend(_attempt_markers(attempt))
        if attempt < attempt_count:
            markers.extend(
                [
                    f"[pip-audit] NETWORK ERROR (attempt {attempt}/4)",
                    f"PROBE sleep={sleeps[attempt - 1]:g}",
                ]
            )
    if scenario == "recover":
        markers.append("No known vulnerabilities found")
        assert "NETWORK FAILURE" not in completed.stdout
    else:
        markers.append("[pip-audit] NETWORK FAILURE")
        assert "No known vulnerabilities found" not in completed.stdout
    _assert_ordered_once(completed.stdout, markers)


def test_real_entrypoint_prints_string_system_exit_message_once(tmp_path: Path) -> None:
    completed, trace = _run_probe(tmp_path, "exit-message")
    assert completed.returncode == 1, completed.stdout
    assert trace["attempts"] == [1]
    assert trace["sleeps"] == []
    _assert_ordered_once(
        completed.stdout, [*_attempt_markers(1), "PROBE string-valued exit message"]
    )
    assert "NETWORK ERROR" not in completed.stdout
    assert "NETWORK FAILURE" not in completed.stdout


def test_real_entrypoint_usage_exit_is_distinct_from_network_exhaustion(tmp_path: Path) -> None:
    completed, trace = _run_probe(tmp_path, "usage")
    assert completed.returncode == 2, completed.stdout
    assert completed.returncode != 3
    assert trace["attempts"] == []
    assert trace["sleeps"] == []
    assert "unrecognized arguments: --issue1185-invalid-option" in completed.stdout
    assert "NETWORK ERROR" not in completed.stdout
    assert "NETWORK FAILURE" not in completed.stdout
