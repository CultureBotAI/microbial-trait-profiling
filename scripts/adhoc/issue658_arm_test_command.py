"""Exercise arm recording with real subprocesses for issue #658 (2026-10-19).

Inputs: explicit fixture file/directory and optional mutation target.
Outputs: explicit small text artifact; selected modes deliberately fail or corrupt
test-owned fixture files to check that the runtime rejects the attempt.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--payload", default="")
    parser.add_argument("--mode", default="copy")
    parser.add_argument("--mutate-file", type=Path)
    args = parser.parse_args()
    if args.mode == "missing":
        return 0
    if args.mode == "symlink":
        args.output.symlink_to(args.input)
        return 0
    if args.mode == "hardlink":
        args.output.hardlink_to(args.input)
        return 0
    payload = args.input.read_bytes() if args.input.is_file() else b"collection"
    args.output.write_bytes(payload + args.payload.encode("utf-8"))
    if args.mode == "fail":
        return 7
    if args.mode == "mutate_input":
        args.input.write_text("modified input", encoding="utf-8")
    if args.mode == "add_input":
        (args.input / "undeclared.txt").write_text("added", encoding="utf-8")
    if args.mutate_file is not None:
        args.mutate_file.write_text("modified by subprocess", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
