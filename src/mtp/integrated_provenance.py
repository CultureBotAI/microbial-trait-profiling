"""Verify every recorded file reference in integrated-evidence provenance (#1755).

Verification resolves and hashes existing file references, including references
inside assessment metadata. Scalar fingerprints without a locator cannot identify
a file to open. Explicitly unavailable and producer-scoped locators are reported
as skips, never counted as verified. The record itself is parsed and hashed from
one buffer. This command reads evidence; it never modifies or regenerates it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path, PureWindowsPath

from mtp.provenance_paths import (
    assert_portable_text,
    resolve_provenance_root,
    resolve_reference,
    verify_reference,
)


def _object(value: object, location: str) -> Mapping[str, object]:
    if not isinstance(value, dict):
        raise ValueError(f"Integrated provenance {location} must be an object")
    return value


def _relative_locator(value: object) -> bool:
    return (
        isinstance(value, str)
        and bool(value)
        and not Path(value).is_absolute()
        and not PureWindowsPath(value).drive
    )


def _input_origin(value: Mapping[str, object], location: str) -> str | None:
    """Validate selection metadata when present; old records have no such field."""
    if "input_origin" not in value:
        return None
    origin = value["input_origin"]
    supplied = value.get("supplied")
    if (
        not isinstance(origin, str)
        or origin not in {"supplied", "checkout_default", "absent"}
        or not isinstance(supplied, bool)
        or supplied != (origin == "supplied")
    ):
        raise ValueError(f"Invalid assessment input origin or supplied flag: {location}")
    if origin == "absent" and value.get("path") is not None:
        raise ValueError(f"Absent input origin has a present source: {location}")
    if value.get("path") is None:
        expected_binding = "missing" if supplied else "not_supplied"
        if value.get("binding") != expected_binding:
            raise ValueError(f"Input origin contradicts unavailable source binding: {location}")
    return origin


def _assessment_source_state(assessment: Mapping[str, object], location: str) -> None:
    """A skipped assessment cannot retain claims that require reading its input.

    This is an implication from an unavailable source, not from an unreadable
    verdict: malformed but readable JSON remains hash-bound and verifiable.
    Legacy records without an explicit availability binding keep their contract.
    """
    source = assessment.get("source")
    if not isinstance(source, Mapping):
        return  # The normal reference walk reports the malformed source object.
    binding = source.get("binding")
    if binding not in ("unreadable", "missing", "not_supplied"):
        return
    if set(source) - {"path", "sha256", "binding", "supplied", "input_origin"}:
        raise ValueError(f"Unavailable assessment source retains read metadata: {location}")
    if binding == "unreadable":
        if not isinstance(source.get("path"), str) or source.get("sha256") is not None:
            raise ValueError(f"Unreadable source has contradictory reference fields: {location}")
    elif source.get("path") is not None or source.get("sha256") is not None:
        raise ValueError(f"Unavailable source has contradictory reference fields: {location}")
    verdicts = assessment.get("verdicts_by_procedure")
    default = assessment.get("status_for_a_procedure_it_is_silent_about")
    state = "unreadable" if binding == "unreadable" else "pending"
    if (
        not isinstance(verdicts, dict)
        or bool(verdicts)
        or not isinstance(default, str)
        or not (default == state or default.startswith(state + " "))
    ):
        raise ValueError(f"Assessment source binding contradicts its status: {location}")


def verify_integrated_provenance(
    record_path: Path,
    *,
    provenance_root: Path | None = None,
) -> dict[str, object]:
    """Verify recorded references and enumerate references that cannot be resolved.

    Legacy records retain their absolute or record-relative source semantics.
    Their nested relative review/protocol locators belonged to another producer,
    whose root was not recorded, so those remain explicitly unverified. Declared
    root records must represent such locators with ``producer_relative_path``.
    """
    try:
        raw = record_path.read_bytes()
        content = json.loads(raw)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError("Cannot read integrated provenance JSON record") from exc
    record = _object(content, "record")
    if record.get("version") != "integrated-evidence-v2":
        raise ValueError("Unsupported integrated provenance record version")
    native = _object(record.get("native_sources"), "native_sources")
    assessments = _object(record.get("assessment_status"), "assessment_status")
    if "ora_eligibility_source" not in record:
        raise ValueError("Integrated provenance record lacks ora_eligibility_source")
    basis = _object(record["path_basis"], "path_basis") if "path_basis" in record else None
    resolve_provenance_root(record_path, path_basis=basis, provenance_root=provenance_root)
    if basis is not None:
        assert_portable_text(raw.decode("utf-8"))
    forced = {f"native_sources.{key}" for key in native}
    for key, value in assessments.items():
        assessment = _object(value, f"assessment_status.{key}")
        if "source" not in assessment:
            raise ValueError(f"Integrated assessment source is absent: {key}")
        _assessment_source_state(assessment, f"assessment_status.{key}.source")
        forced.add(f"assessment_status.{key}.source")
    if record["ora_eligibility_source"] is not None:
        forced.add("ora_eligibility_source")

    verified: list[str] = []
    skipped: list[dict[str, str]] = []

    def skip(location: str, reason: str) -> None:
        skipped.append({"location": location, "reason": reason})

    def reference(value: Mapping[str, object], location: str) -> None:
        origin = _input_origin(value, location)
        if "producer_relative_path" in value:
            if "path" in value or not _relative_locator(value["producer_relative_path"]):
                raise ValueError(f"Ambiguous producer-scoped reference: {location}")
            skip(location, "producer_relative_path")
            return
        if "path" not in value:
            if not value:
                skip(location, "no_reference_recorded")
                return
            raise ValueError(f"Recorded source reference lacks a path: {location}")
        path = value["path"]
        digest = value.get("sha256")
        if value.get("binding") == "unreadable":
            if (
                not location.startswith("assessment_status.")
                or not location.endswith(".source")
                or not isinstance(path, str)
                or digest is not None
            ):
                raise ValueError(
                    f"Unreadable source has contradictory reference fields: {location}"
                )
            # Validate the locator even though there were no readable bytes to
            # attest. It remains explicitly unverified, including after relocation.
            resolve_reference(record_path, path, path_basis=basis, provenance_root=provenance_root)
            skip(location, "unreadable")
            return
        if path is None:
            binding = value.get("binding")
            if digest is not None or binding not in (None, "not_supplied", "missing"):
                raise ValueError(
                    f"Unavailable source has contradictory reference fields: {location}"
                )
            if binding == "missing":
                if value.get("supplied") is not True:
                    raise ValueError(f"Missing source does not declare supplied=true: {location}")
                skip(location, "supplied_but_missing")
            elif binding == "not_supplied":
                if value.get("supplied") is True:
                    raise ValueError(f"Source is both supplied and not supplied: {location}")
                skip(
                    location,
                    "defaulted_but_missing" if origin == "checkout_default" else "not_supplied",
                )
            else:
                # Older native-source records omitted a presence declaration.
                skip(location, "unavailable_at_recording")
            return
        if basis is None and location.endswith(".source.review") and _relative_locator(path):
            skip(location, "legacy_producer_relative_review")
            return
        try:
            verify_reference(record_path, value, path_basis=basis, provenance_root=provenance_root)
        except ValueError as exc:
            raise ValueError(f"Invalid source reference {location}: {exc}") from exc
        verified.append(location)

    def walk(value: object, location: str) -> None:
        if location in forced and not isinstance(value, Mapping):
            raise ValueError(f"Recorded source reference must be an object: {location}")
        if isinstance(value, Mapping):
            if location in forced or any(
                key in value for key in ("path", "sha256", "producer_relative_path")
            ):
                reference(value, location)
            for key, item in value.items():
                child = f"{location}.{key}" if location else str(key)
                if (
                    location in forced
                    and location.startswith("assessment_status.")
                    and key in ("review", "protocol")
                    and item is not None
                    and not (basis is None and key == "protocol" and isinstance(item, str))
                    and (
                        not isinstance(item, Mapping)
                        or not any(
                            locator in item for locator in ("path", "producer_relative_path")
                        )
                    )
                ):
                    raise ValueError(
                        f"Optional assessment locator lacks a reference scope: {child}"
                    )
                if key == "protocol" and isinstance(item, str):
                    if basis is not None:
                        raise ValueError(
                            f"Protocol locator lacks an explicit reference scope: {child}"
                        )
                    if _relative_locator(item):
                        skip(child, "legacy_producer_relative_protocol")
                    else:
                        reference({"path": item, "sha256": value.get("protocol_sha256")}, child)
                else:
                    walk(item, child)
        elif isinstance(value, list):
            for index, item in enumerate(value):
                walk(item, f"{location}[{index}]")
        elif isinstance(value, str) and value.startswith(("{", "[")):
            try:
                decoded = json.loads(value)
            except json.JSONDecodeError:
                return
            walk(decoded, location)

    walk(record, "")
    return {
        "record_sha256": hashlib.sha256(raw).hexdigest(),
        "verified_reference_count": len(verified),
        "verified_references": verified,
        "skipped_reference_count": len(skipped),
        "skipped_references": skipped,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record", type=Path, help="Integrated-evidence provenance JSON")
    parser.add_argument(
        "--provenance-root",
        type=Path,
        help="Safe source root required by a supplied-root record",
    )
    args = parser.parse_args(argv)
    try:
        result = verify_integrated_provenance(args.record, provenance_root=args.provenance_root)
    except ValueError as exc:
        print(f"Integrated provenance verification failed: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
