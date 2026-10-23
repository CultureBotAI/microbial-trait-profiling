"""Compatibility entry point for the label-support sidecar."""

from mtp.label_support_audit import main

__all__ = ["main"]

if __name__ == "__main__":
    raise SystemExit(main())
