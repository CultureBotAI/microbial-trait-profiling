"""How unusual is the PD, HIV and CRC overlap among all condition triples? (#835)

Parkinson's disease, HIV infection and colorectal cancer, each with cellulose catabolic
process decreased, were picked as exploratory leads from this corpus. This scan repeats
cellulose_decr's three-way core and mean pairwise Jaccard distance, against its
whole-study-signature null (U6), for every triple of eligible conditions and for the top
decreasing traits of the Noel arm's directional screen. It reports where PD, HIV and CRC
fall among the triples, and where cellulose falls among the traits.

The rule is declared in ``methods/cellulose_decr_selection_scan.md`` and fixed in
``methods/cellulose_decr_selection_scan_protocol.json`` (v1), both committed before this
code. Protocol v2 (``methods/cellulose_decr_selection_scan_protocol_v2.json``) amends v1
once, on 2026-10-25, after the results of both scopes had been seen: phenotype becomes the
primary strain-rollup scope and off the declared sensitivity. As a post-results choice of
primary scope, it reports the two side by side, always (#1050;
``mtp.cellulose_decr.selection_scan_scopes``). A run follows v2. v1 stays byte-identical
and readable, so the records written under it still verify.
The population is a verified Tier 1+2 core run rather than cellulose_decr's own
``data/derived`` inputs, which give different set sizes (#818).

The sets, the U6 units and the U6 draws come from ``mtp.cellulose_decr.analysis``. One
sequence of U6 draws serves every triple and trait: each eligible condition takes its own
block of each draw's studies, so the draws a triple sees are U6's draws for that triple.
A triple whose blocks would overlap gets its own U6 draws.

The selection-aware lane of #835 (#854) splits the studies. Load a core once with
``load_core_inputs``, apply the declared rule to the training studies with ``select``, and
test the chosen hypotheses on the held-out studies with ``evaluate``. Each reads only the
studies it is given. The CLI is ``select`` on every study.

The counts are annotation counts of curated signature taxa, not measurements of patients,
cells or pathway activity (methods/biological_interpretation.md). Every p-value is nominal.
"""

from __future__ import annotations

import argparse
import dataclasses
import itertools
import json
import sys
import time
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, date, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy

from mtp.cellulose_decr.analysis import (
    CONDITIONS,
    FECES,
    FLOAT_FORMAT,
    TRAIT_ID,
    clean,
    clean_label,
    compartment_of,
    read_tsv,
    signature_draws,
    trait_taxa,
    u6_units,
    write_json,
    write_tsv,
)
from mtp.paths import REPO_ROOT, file_sha256
from mtp.trait_guilds import Profile, load_profiles
from mtp.trait_screen import (
    TESTED,
    ColumnLabel,
    candidate_slots,
    load_column_labels,
    run_screen,
    screen_candidates,
)
from mtp.workflow.output_roots import add_temporary_output_option, require_durable_output

if TYPE_CHECKING:
    from mtp.datamodel.workflow_pydantic import WorkflowRun
    from mtp.workflow.consumers import CoreBundle

PROTOCOL_V1_PATH = REPO_ROOT / "methods" / "cellulose_decr_selection_scan_protocol.json"
PROTOCOL_PATH = REPO_ROOT / "methods" / "cellulose_decr_selection_scan_protocol_v2.json"
METHODS_PATH = REPO_ROOT / "methods" / "cellulose_decr_selection_scan.md"
VERSION_V1 = "cellulose-decr-selection-scan-v1"
VERSION = "cellulose-decr-selection-scan-v2"
# Every protocol version this code reads, with its committed file. A run follows VERSION;
# v1 stays readable so that the records written under it still verify (#1050).
PROTOCOLS = {VERSION_V1: PROTOCOL_V1_PATH, VERSION: PROTOCOL_PATH}
# The version a later version amends. It names that version's file by path and SHA-256,
# and adds amendments to the ones in force there.
AMENDS = {VERSION: VERSION_V1}
# The fields of a file that amends an earlier version: nothing else, so no rule field can
# be set there except through an amendment.
AMENDING_FILE_KEYS = {
    "protocol",
    "issue",
    "methods",
    "amends",
    "declared",
    "drafted_by",
    "amendments",
}
TIER1_2 = "disease_vs_healthy_tier1_2_subject_controls"
CORE_PRODUCTS = (
    "experiments",
    "reported_taxa",
    "count_profiles",
    "taxa_matrix",
    "feature_metadata",
)

MEAN_DISTANCE = "mean_pairwise_jaccard_distance"
CORE = "three_way_core"
ALL_TRIPLES = "all_tested_triples"
DISJOINT = "study_disjoint_triples"
REFERENCE_SETS = (ALL_TRIPLES, DISJOINT)

# A triple's test status for one trait.
NO_TAXA = "all_sets_empty"
POOL_TOO_SMALL = "pool_too_small"
# A held-out hypothesis with a condition that no held-out study reports (``evaluate``).
NO_UNITS = "no_units"
# Why a selection's family has triples or none (``Selection.status``, #901).
SELECTED = "selected"
NO_TRAITS = "no_traits"
TOO_FEW_CONDITIONS = "fewer_than_three_eligible_conditions"
# Where a triple's U6 draws come from.
SHARED = "shared"
OWN = "own"
NOT_DRAWN = "none"

# Draws vectorized at once. It bounds memory and leaves every number unchanged.
CHUNK = 1000
QUANTILES = (0.01, 0.05, 0.1, 0.25, 0.5)

OUTPUTS = {
    "conditions": "selection_scan_conditions.tsv",
    "traits": "selection_scan_traits.tsv",
    "triples": "selection_scan_triples.tsv.gz",
    "distribution": "selection_scan_distribution.tsv",
    "summary": "selection_scan_summary.json",
    "provenance": "selection_scan_provenance.json",
}
CONSUMER_RECORD = "core_consumption.json"
CONSUMER = "mtp.cellulose_decr.selection_scan"
# The core run's external sources, recorded with their hashes so that two scans can be
# shown to read the same corpus and KG (#1050).
CORE_SOURCES = ("raw_dump", "kg_edges", "kg_nodes")

BoolArray = npt.NDArray[np.bool_]
FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]
WordArray = npt.NDArray[np.uint64]
# A seed for numpy's SeedSequence: an integer, or integers such as (seed, split, stage).
SeedLike = int | Sequence[int]


def seed_entropy(seed: SeedLike) -> list[int]:
    """A seed as a list of integers; ``default_rng(s)`` and ``default_rng([s])`` agree."""
    return [int(seed)] if isinstance(seed, int | np.integer) else [int(part) for part in seed]


@dataclass(frozen=True)
class Protocol:
    """The declared rule's parameters, read from a protocol file with its amendments.

    ``results_seen`` lists the scopes whose scan results had been seen when the
    strain-rollup roles were declared; ``side_by_side`` is whether the scopes must be
    reported together (both set by v2, #1050).
    """

    version: str
    min_studies: int
    top_n: int
    always_included: str
    permutations: int
    seed: int
    tolerance: float
    lead_triple: tuple[str, ...]
    lead_trait: str
    bins: int
    scope_primary: str
    scope_sensitivity: tuple[str, ...]
    scope_absent: str
    results_seen: tuple[str, ...]
    side_by_side: bool


STRAIN_ROLLUP = "strain_rollup"
PRIMARY = "primary"
SENSITIVITY = "declared_sensitivity"
SIDE_BY_SIDE = "side_by_side"
# Where a run's scope came from: the core's recorded parameter, or, for a core built
# before #828 that records none, the protocol's declared ``absent_means`` (#927).
RECORDED = "recorded"
ABSENT_MEANS = "absent_means"
# What this code implements for each field the protocol fixes. load_protocol refuses a
# protocol that says otherwise, so an edited field cannot pass as followed (#872).
FIXED_FIELDS: dict[tuple[str, ...], object] = {
    ("population", "core_contrast_population"): TIER1_2,
    ("population", "stool_body_site"): FECES,
    ("population", "core_products"): list(CORE_PRODUCTS),
    ("statistics", "primary"): MEAN_DISTANCE,
    ("statistics", "secondary"): CORE,
    ("reference_sets",): list(REFERENCE_SETS),
    ("lead_triple",): list(CONDITIONS.values()),
    ("lead_trait",): TRAIT_ID,
}
# The strain-rollup scope each version declares, as this code implements it: v1's
# amendment of 2026-10-22 (#872) made off primary; v2's of 2026-10-25, made after the
# results of both scopes had been seen, makes phenotype primary and both scopes
# reported (#1050).
SCOPE_DECLARATIONS: dict[str, dict[str, object]] = {
    VERSION_V1: {
        "parameter": STRAIN_ROLLUP,
        "primary": "off",
        "sensitivity": ["phenotype"],
        "absent_means": "off",
    },
    VERSION: {
        "parameter": STRAIN_ROLLUP,
        "primary": "phenotype",
        "sensitivity": ["off"],
        "absent_means": "off",
        "reported": SIDE_BY_SIDE,
        "declared_after_results": ["off", "phenotype"],
    },
}
# The secondary rank the amendment of 2026-10-22 declares (#872), in force in both versions.
SECONDARY_RANK: dict[str, object] = {
    "by": "u6_z",
    MEAN_DISTANCE: "ascending",
    CORE: "descending",
    "ties": "range",
}
AMENDABLE_FIELDS = (STRAIN_ROLLUP, "secondary_rank")
AMENDMENT_KEYS = {"date", "change", "reason", "known_when_amended", "sets"}


def repository_path(path: Path) -> str:
    """A committed file's path within the repository, as the records give it."""
    return str(path.relative_to(REPO_ROOT))


def declaration(path: Path) -> tuple[str, dict[str, Any], list[dict[str, Any]]]:
    """A protocol file's version, the fields it declares, and every amendment in force.

    A version that amends an earlier one (``AMENDS``) must name that version's committed
    file by path and SHA-256 and hold nothing but amendments; the earlier version's fields
    and amendments come first, unchanged. A file of an unknown version is refused.
    """
    data = json.loads(path.read_text(encoding="utf-8"))
    version = data["protocol"]
    if version not in PROTOCOLS:
        raise ValueError(f"{path} declares {version!r}; this code reads {sorted(PROTOCOLS)}")
    base = AMENDS.get(version)
    if base is None:
        if "amends" in data:
            raise ValueError(f"{path}: {version} amends no earlier version")
        return version, data, list(data["amendments"])
    base_path = PROTOCOLS[base]
    reference = data.get("amends", {})
    if set(reference) != {"protocol", "path", "sha256"} or (
        reference["protocol"],
        reference["path"],
    ) != (base, repository_path(base_path)):
        raise ValueError(f"{path}: {version} must amend {base} at {repository_path(base_path)}")
    if file_sha256(base_path) != reference["sha256"]:
        raise ValueError(f"{repository_path(base_path)} has changed since {version} amended it")
    extra = sorted(set(data) - AMENDING_FILE_KEYS)
    if extra:
        raise ValueError(f"{path}: {version} may set fields only by amendment, not {extra}")
    _, declared, amendments = declaration(base_path)
    return version, declared, [*amendments, *data["amendments"]]


def load_protocol(path: Path = PROTOCOL_PATH) -> Protocol:
    """The declared parameters of the protocol in ``path``, with its amendments applied in
    date order: by default v2, which a run follows; ``PROTOCOL_V1_PATH`` reads v1.

    A protocol of an unknown version is refused, and so is one whose fields say something
    this code does not do for its version (``FIXED_FIELDS``, ``SCOPE_DECLARATIONS``,
    ``SECONDARY_RANK``), whose amendments are malformed, out of date order or set an
    unknown field, or that does not amend the earlier version as declared (``declaration``).
    """
    version, data, amendments = declaration(path)
    effective = dict(data)
    dates = []
    for number, amendment in enumerate(amendments):
        if set(amendment) != AMENDMENT_KEYS:
            raise ValueError(f"{path}: amendment {number} needs exactly {sorted(AMENDMENT_KEYS)}")
        dates.append(date.fromisoformat(amendment["date"]))
        unknown = sorted(set(amendment["sets"]) - set(AMENDABLE_FIELDS))
        if unknown:
            raise ValueError(f"{path}: amendment {number} sets unknown fields {unknown}")
        effective.update(amendment["sets"])
    if dates != sorted(dates):
        raise ValueError(f"{path}: amendments must be in date order")
    for keys, expected in FIXED_FIELDS.items():
        value: Any = effective
        for key in keys:
            value = value[key]
        if value != expected:
            raise ValueError(f"{path}: {'.'.join(keys)} is {value!r}; this code does {expected!r}")
    amended = {STRAIN_ROLLUP: SCOPE_DECLARATIONS[version], "secondary_rank": SECONDARY_RANK}
    for field, expected in amended.items():
        if effective.get(field) != expected:
            raise ValueError(
                f"{path}: {field} is {effective.get(field)!r}; this code does {expected!r}"
            )
    rollup = effective[STRAIN_ROLLUP]
    return Protocol(
        version=version,
        min_studies=int(data["eligibility"]["min_studies"]),
        top_n=int(data["traits"]["top_n"]),
        always_included=str(data["traits"]["always_included"]),
        permutations=int(data["null"]["permutations"]),
        seed=int(data["null"]["seed"]),
        tolerance=float(data["null"]["distance_tolerance"]),
        lead_triple=tuple(data["lead_triple"]),
        lead_trait=str(data["lead_trait"]),
        bins=int(data["distribution_bins"]),
        scope_primary=str(rollup["primary"]),
        scope_sensitivity=tuple(rollup["sensitivity"]),
        scope_absent=str(rollup["absent_means"]),
        results_seen=tuple(rollup.get("declared_after_results", ())),
        side_by_side=rollup.get("reported") == SIDE_BY_SIDE,
    )


def protocol_record(version: str = VERSION) -> dict[str, Any]:
    """A version's committed file, its hash, and the file and hash of the version it
    amends, as a scan's provenance records them."""
    path = PROTOCOLS[version]
    record: dict[str, Any] = {
        "path": repository_path(path),
        "sha256": file_sha256(path),
        "version": version,
    }
    base = AMENDS.get(version)
    if base is not None:
        record["amends"] = {
            "protocol": base,
            "path": repository_path(PROTOCOLS[base]),
            "sha256": file_sha256(PROTOCOLS[base]),
        }
    return record


def scope_role(scope: str, protocol: Protocol) -> str:
    """A strain-rollup scope's declared role, primary or sensitivity; another is refused.

    The scope is the core's (``core_scope``).
    """
    if scope == protocol.scope_primary:
        return PRIMARY
    if scope in protocol.scope_sensitivity:
        return SENSITIVITY
    raise ValueError(
        f"strain_rollup {scope!r} is neither the declared primary scope "
        f"{protocol.scope_primary!r} nor a declared sensitivity {list(protocol.scope_sensitivity)}"
    )


def core_scope(run: WorkflowRun, protocol: Protocol) -> tuple[str, str]:
    """A core run's strain-rollup scope, and where it came from.

    A run that records the parameter has that scope (``recorded``), read with
    ``mtp.workflow.operations.run_strain_rollup``. A core built before #828 records none:
    its scope is the protocol's declared ``absent_means`` (``absent_means``), off, which
    v2 makes the sensitivity. No default outside the protocol decides it (#927).
    """
    from mtp.workflow.operations import run_strain_rollup

    if STRAIN_ROLLUP in {item.name for item in run.parameters or []}:
        return str(run_strain_rollup(run)), RECORDED
    return protocol.scope_absent, ABSENT_MEANS


def scope_record(scope: str, source: str, protocol: Protocol) -> dict[str, Any]:
    """What a run records about its scope, in its summary, provenance and ``--plan``.

    The scope and where it came from (``core_scope``); the role the protocol declares for
    it, and the version that declared it after which results had been seen; and, when the
    scopes are reported side by side (v2, #1050), the scopes it must be reported beside.
    """
    declared = [protocol.scope_primary, *protocol.scope_sensitivity]
    return {
        "scope": scope,
        "scope_source": source,
        "role": scope_role(scope, protocol),
        "declared_primary": protocol.scope_primary,
        "declared_sensitivity": list(protocol.scope_sensitivity),
        "declared_by": protocol.version,
        "declared_after_results": list(protocol.results_seen),
        "reported": SIDE_BY_SIDE if protocol.side_by_side else None,
        "report_beside": [other for other in declared if other != scope]
        if protocol.side_by_side
        else [],
    }


def recorded_protocol(summary: Mapping[str, Any], provenance: Mapping[str, Any]) -> Protocol:
    """The protocol a scan's summary and provenance name, once they agree with it.

    Both must name the same known version; the provenance must give that version's
    committed file, its current hash and, for v2, the v1 file and hash it amends; and the
    recorded scope's role and the declared scopes must be the ones that version declares.
    A record written under v1, such as those of 2026-10-23, verifies against the v1 file,
    which v2 leaves unchanged (#1050).
    """
    version = summary["protocol"]
    recorded = provenance["protocol"]
    if recorded["version"] != version:
        raise ValueError(
            f"the summary names {version!r} but the provenance names {recorded['version']!r}"
        )
    if version not in PROTOCOLS:
        raise ValueError(f"the records name an unknown protocol version {version!r}")
    committed = protocol_record(version)
    if recorded["path"] != committed["path"]:
        raise ValueError(f"the records give {version} the path {recorded['path']!r}")
    if recorded["sha256"] != committed["sha256"]:
        raise ValueError(f"the recorded hash of {committed['path']} is not the committed file's")
    if recorded != committed:
        raise ValueError(f"the records of {version} do not match {committed}")
    protocol = load_protocol(PROTOCOLS[version])
    rollup = summary[STRAIN_ROLLUP]
    if provenance[STRAIN_ROLLUP] != rollup:
        raise ValueError("the summary and the provenance differ in the strain-rollup scope")
    expected = {
        "role": scope_role(rollup["scope"], protocol),
        "declared_primary": protocol.scope_primary,
        "declared_sensitivity": list(protocol.scope_sensitivity),
    }
    actual = {key: rollup.get(key) for key in expected}
    if actual != expected:
        raise ValueError(
            f"the records give strain_rollup {rollup['scope']!r} the role and declaration "
            f"{actual}; {version} declares {expected}"
        )
    return protocol


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument(
        "--core-run",
        type=Path,
        required=True,
        help=f"Verified core run manifest (run.json) of the {TIER1_2} population.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="New or empty directory for the scan's outputs.",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        help="U6 draws (default: the protocol's). Another value is recorded as off-protocol.",
    )
    parser.add_argument(
        "--seed", type=int, help="U6 seed (default: the protocol's); recorded like --permutations."
    )
    parser.add_argument(
        "--plan",
        action="store_true",
        help="Verify the core and print the eligible conditions, traits, triples and a timed "
        "estimate, without running the null or writing outputs.",
    )
    parser.add_argument(
        "--code-commit",
        default="unrecorded",
        help="Revision of this script, e.g. $(git rev-parse HEAD) with +dirty when modified.",
    )
    add_temporary_output_option(parser)
    args = parser.parse_args(argv)
    # A plan writes nothing; a scan's outputs are cited, so they must outlive a restart (#1087)
    if not args.plan:
        require_durable_output(parser, args.output_dir, allowed=args.allow_temporary_output)
    return args


# --------------------------------------------------------------------------- population


def load_population(
    experiments_path: Path,
    taxa_path: Path,
    profile_path: Path,
    matrix_taxa: pd.Index,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read a core run's experiments, profile keys and signature taxa (``build_population``)."""
    keys = ["study_id", "experiment_id"]
    keyed = pd.read_csv(
        profile_path, sep="\t", usecols=[*keys, "inc_dec"], dtype=str, keep_default_na=False
    )
    return build_population(read_tsv(experiments_path), keyed, read_tsv(taxa_path), matrix_taxa)


def build_population(
    experiments: pd.DataFrame,
    keyed: pd.DataFrame,
    reported: pd.DataFrame,
    matrix_taxa: pd.Index,
    studies: Iterable[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Profiles and signature taxa of a core's population, as cellulose_decr frames them.

    ``experiments``, ``keyed`` (the profiles' study, experiment and ``inc_dec``) and
    ``reported`` (the signature taxa) are the core's tables in memory; ``studies``, when
    given, keeps only those studies, whole. The core has already selected and oriented
    the experiments, so a profile's own direction is the used one: ``dec`` means lower
    in the disease group. The population is human and paired, which is checked.
    """
    keys = ["study_id", "experiment_id"]
    if studies is not None:
        keep = frozenset(studies)
        unknown = sorted(keep - set(experiments["study_id"]))
        if unknown:
            raise ValueError(f"{len(unknown)} studies are not in the population: {unknown[:3]}")
        experiments = experiments[experiments["study_id"].isin(keep)]
        keyed = keyed[keyed["study_id"].isin(keep)]
        reported = reported[reported["study_id"].isin(keep)]
    experiments = experiments[[*keys, "condition", "body_site", "host_species"]]
    if experiments.duplicated(keys).any():
        raise ValueError("the experiments table repeats an experiment")
    hosts = set(experiments["host_species"])
    if hosts - {"Homo sapiens"}:
        raise ValueError(f"expected human experiments only, found hosts {sorted(hosts)}")
    if not set(keyed["inc_dec"]) <= {"inc", "dec"}:
        raise ValueError(f"unexpected profile directions {sorted(set(keyed['inc_dec']))}")
    profiled = set(keyed[keys].itertuples(index=False, name=None))
    listed = set(experiments[keys].itertuples(index=False, name=None))
    if profiled != listed:
        raise ValueError(
            f"{len(listed - profiled)} experiments have no profile and {len(profiled - listed)} "
            "profiles have no experiment"
        )
    if keyed.duplicated([*keys, "inc_dec"]).any():
        raise ValueError("the profiles repeat a profile")
    if not keyed.groupby(keys)["inc_dec"].nunique().eq(2).all():
        raise ValueError("the population must be paired: every experiment needs both directions")
    merged = keyed.merge(experiments, on=keys, how="left", validate="many_to_one")
    merged.index = pd.Index(
        merged["study_id"] + "|" + merged["experiment_id"] + "|" + merged["inc_dec"],
        name="profile_id",
    )
    profiles = merged.assign(
        compartment=merged["body_site"].map(compartment_of), dir=merged["inc_dec"]
    )[[*keys, "condition", "body_site", "compartment", "dir"]]

    taxa = reported.copy()
    direction = taxa["direction"].map({"increased": "inc", "decreased": "dec"})
    if direction.isna().any():
        raise ValueError(
            f"unexpected directions {sorted(set(taxa.loc[direction.isna(), 'direction']))}"
        )
    taxa["profile_id"] = taxa["study_id"] + "|" + taxa["experiment_id"] + "|" + direction
    taxa = taxa[taxa["profile_id"].isin(profiles.index)]
    taxa = taxa.drop_duplicates(["profile_id", "taxon_id"])[["profile_id", "taxon_id"]]
    taxa = taxa[taxa["taxon_id"].isin(matrix_taxa)].reset_index(drop=True)
    labels = profiles.loc[taxa["profile_id"], [*keys, "condition", "compartment", "dir"]]
    taxa = pd.concat([taxa, labels.reset_index(drop=True)], axis=1)
    return profiles, taxa


def condition_table(profiles: pd.DataFrame, min_studies: int) -> pd.DataFrame:
    """Every named condition's stool comparisons, studies and U6 units, and its eligibility.

    A condition is eligible when its stool comparisons come from at least ``min_studies``
    studies; three comparisons from one study are one cohort and one U6 unit.
    """
    stool = profiles[profiles["compartment"].eq(FECES)]
    stool = stool[stool["condition"].map(clean_label).ne("unknown")]
    comparisons = stool[["study_id", "experiment_id", "condition"]].drop_duplicates()
    decreased = stool[stool["dir"].eq("dec")]
    frame = pd.DataFrame(
        {
            "stool_comparisons": comparisons.groupby("condition").size(),
            "stool_studies": comparisons.groupby("condition")["study_id"].nunique(),
            "u6_units": decreased.groupby("condition")["study_id"].nunique(),
        }
    )
    frame = frame.fillna(0).astype(np.int64).rename_axis("condition").reset_index()
    frame["eligible"] = frame["stool_studies"] >= min_studies
    return frame.sort_values("condition").reset_index(drop=True)


# --------------------------------------------------------------------------- traits


def screen_table(
    profiles: Sequence[Profile],
    columns: Sequence[str],
    labels: Mapping[str, ColumnLabel],
    studies: Iterable[str] | None = None,
) -> pd.DataFrame:
    """The Noel arm's directional screen on these profiles, one row per candidate trait.

    ``studies``, when given, keeps only the profiles of those studies. The screen votes
    over every body site; the scan's sets are stool-only.
    """
    if studies is not None:
        keep = frozenset(studies)
        profiles = [profile for profile in profiles if profile.study_id in keep]
    result = run_screen(list(profiles), columns, labels=dict(labels))
    return pd.DataFrame(
        [
            {
                "trait_id": row.trait,
                "trait_label": row.trait_label,
                "studies_up": row.studies_up,
                "studies_down": row.studies_down,
                "studies_tied": row.studies_tied,
                "studies_scored": row.studies_scored,
                "share_up": row.share_up,
                "screen_p": row.p_value,
                "screen_test_status": row.test_status,
            }
            for row in result.rows
        ]
    )


def choose_traits(screen: pd.DataFrame, top_n: int, always: str | None) -> pd.DataFrame:
    """The top decreasing traits in the screen report's order, and ``always`` if missing.

    A trait is decreasing when the screen tested it and more studies shift down than up.
    The order is the screen report's display order: p ascending, then more scored
    studies; trait ID breaks any remaining tie.
    """
    decreasing = screen[
        screen["screen_test_status"].eq(TESTED) & screen["studies_down"].gt(screen["studies_up"])
    ]
    decreasing = decreasing.sort_values(
        ["screen_p", "studies_scored", "trait_id"], ascending=[True, False, True]
    ).reset_index(drop=True)
    decreasing.insert(0, "decreasing_rank", np.arange(1, len(decreasing) + 1))
    chosen = decreasing.head(top_n).assign(in_top_n=True)
    if always is not None and always not in set(chosen["trait_id"]):
        row = screen[screen["trait_id"].eq(always)]
        if row.empty:
            raise ValueError(f"{always} is not a screen candidate")
        ranks = decreasing[["trait_id", "decreasing_rank"]]
        row = row.merge(ranks, on="trait_id", how="left").assign(in_top_n=False)
        chosen = pd.concat([chosen, row[chosen.columns]], ignore_index=True)
    chosen["decreasing_rank"] = chosen["decreasing_rank"].astype("Int64")
    return chosen.reset_index(drop=True)


# --------------------------------------------------------------------------- packed sets


def pack(bits: BoolArray) -> WordArray:
    """Each row's booleans as 64-bit words; an empty row still has one word."""
    rows, width = bits.shape
    words = max(1, -(-width // 64))
    padded = np.zeros((rows, words * 64), dtype=bool)
    padded[:, :width] = bits
    packed: WordArray = np.packbits(padded, axis=1).view(np.uint64)
    return packed


def word_counts(words: WordArray, offsets: IntArray) -> IntArray:
    """Set sizes per trait: the bits set in each trait's run of words (last axis)."""
    counts: IntArray = np.add.reduceat(np.bitwise_count(words), offsets, axis=-1, dtype=np.int64)
    return counts


def triple_blocks(
    masks: WordArray, offsets: IntArray
) -> Iterator[tuple[slice, FloatArray, IntArray]]:
    """Each triple's mean pairwise Jaccard distance and three-way core, for every trait.

    ``masks`` holds each condition's packed sets, ``(draws, conditions, words)``, with
    each trait's words starting at its ``offsets`` entry. Triples come in
    ``itertools.combinations`` order, one block for each pair ``a < b`` with every
    ``c > b``. A pair whose union is empty has no distance and the mean is over the
    pairs that have one. The arithmetic is ``NullDraws.add``'s, term for term, so a
    triple laid out as (PD, HIV, CRC) gets cellulose_decr's values exactly.
    """
    k = masks.shape[1]
    sizes = word_counts(masks, offsets)
    distance: list[FloatArray] = []
    for a in range(k - 1):
        inter = word_counts(masks[:, a : a + 1] & masks[:, a + 1 :], offsets)
        union = sizes[:, a : a + 1] + sizes[:, a + 1 :] - inter
        with np.errstate(invalid="ignore", divide="ignore"):
            distance.append(np.where(union > 0, 1 - inter / np.maximum(union, 1), np.nan))
    start = 0
    for a in range(k - 2):
        for b in range(a + 1, k - 1):
            width = k - b - 1
            pairs = np.broadcast_arrays(
                distance[a][:, b - a - 1 : b - a], distance[a][:, b - a :], distance[b]
            )
            defined = sum((~np.isnan(pair)).astype(np.int64) for pair in pairs)
            total = (np.nan_to_num(pairs[0]) + np.nan_to_num(pairs[1])) + np.nan_to_num(pairs[2])
            with np.errstate(invalid="ignore", divide="ignore"):
                mean = np.where(defined > 0, total / np.maximum(defined, 1), np.nan)
            core = word_counts((masks[:, a] & masks[:, b])[:, None] & masks[:, b + 1 :], offsets)
            yield slice(start, start + width), mean, core
            start += width


def triple_values(masks: WordArray, offsets: IntArray) -> tuple[FloatArray, IntArray]:
    """``triple_blocks`` for one set of masks, as ``(triples, traits)`` arrays."""
    blocks = list(triple_blocks(masks, offsets))
    if not blocks:
        empty = (0, len(offsets))
        return np.empty(empty), np.empty(empty, dtype=np.int64)
    return (
        np.concatenate([mean[0] for _, mean, _ in blocks]),
        np.concatenate([core[0] for _, _, core in blocks]),
    )


def block_positions(counts: Sequence[int], n_studies: int) -> list[IntArray]:
    """Consecutive blocks of a draw's study order, one per condition, wrapping at the end."""
    offsets = np.concatenate([[0], np.cumsum(counts)[:-1]])[: len(counts)].astype(np.int64)
    return [
        (offset + np.arange(count, dtype=np.int64)) % n_studies
        for offset, count in zip(offsets, counts, strict=True)
    ]


def draw_masks(
    unit_bits: WordArray, unit_of_study: IntArray, order: IntArray, positions: Sequence[IntArray]
) -> WordArray:
    """Each block's synthetic set in each draw: the union of its studies' chosen units.

    This is ``signature_null``'s step from a draw to a set, on packed words.
    """
    masks = np.empty((order.shape[0], len(positions), unit_bits.shape[1]), dtype=np.uint64)
    for column, block in enumerate(positions):
        chosen = np.take_along_axis(unit_of_study, order[:, block], axis=1)
        masks[:, column, :] = np.bitwise_or.reduce(unit_bits[chosen], axis=1)
    return masks


# --------------------------------------------------------------------------- the null


@dataclass
class Totals:
    """Running U6 counts and sums for every triple and trait."""

    distance_at_most: IntArray
    distance_defined: IntArray
    distance_sum: FloatArray
    distance_squares: FloatArray
    core_at_least: IntArray
    core_sum: IntArray
    core_squares: IntArray
    draws: IntArray

    @classmethod
    def zeros(cls, shape: tuple[int, int]) -> Totals:
        def ints() -> IntArray:
            return np.zeros(shape, dtype=np.int64)

        def floats() -> FloatArray:
            return np.zeros(shape, dtype=np.float64)

        return cls(ints(), ints(), floats(), floats(), ints(), ints(), ints(), ints())

    def add(
        self,
        rows: slice,
        mean: FloatArray,
        core: IntArray,
        observed_mean: FloatArray,
        observed_core: IntArray,
        tolerance: float,
    ) -> None:
        """Count a chunk of draws, ``(draws, triples, traits)``, against the observed values."""
        defined = ~np.isnan(mean)
        values = np.where(defined, mean, 0.0)
        self.distance_at_most[rows] += (mean <= observed_mean + tolerance).sum(axis=0)
        self.distance_defined[rows] += defined.sum(axis=0)
        self.distance_sum[rows] += values.sum(axis=0)
        self.distance_squares[rows] += (values * values).sum(axis=0)
        self.core_at_least[rows] += (core >= observed_core).sum(axis=0)
        self.core_sum[rows] += core.sum(axis=0)
        self.core_squares[rows] += (core * core).sum(axis=0)
        self.draws[rows] += mean.shape[0]


@dataclass(frozen=True)
class ScanResult:
    """Observed statistics and U6 totals for every triple of eligible conditions and trait."""

    conditions: tuple[str, ...]
    traits: tuple[str, ...]
    triples: IntArray
    units: IntArray
    sizes: IntArray
    union: IntArray
    observed_mean: FloatArray
    observed_core: IntArray
    totals: Totals
    source: npt.NDArray[np.str_]
    shared_studies: IntArray
    pool_studies: int
    u3_taxa: IntArray


@dataclass(frozen=True)
class Prepared:
    """Everything the null needs, built once from the population.

    ``triples`` index ``conditions``; ``complete`` says they are every combination, in
    ``itertools.combinations`` order, which ``triple_blocks`` scores a block at a time.
    """

    conditions: tuple[str, ...]
    traits: tuple[str, ...]
    triples: IntArray
    complete: bool
    units: IntArray
    unit_bits: WordArray
    offsets: IntArray
    unit_study: IntArray
    observed: WordArray
    incidence: BoolArray
    sets: list[BoolArray]
    u3_taxa: IntArray
    positions: list[IntArray]
    source: npt.NDArray[np.str_]


def prepare(
    profiles: pd.DataFrame,
    taxa: pd.DataFrame,
    matrix: pd.DataFrame,
    traits: Sequence[str],
    conditions: Sequence[str],
    triples: IntArray | None = None,
) -> Prepared:
    """U6 units and observed sets for each trait, and the block layout of the conditions.

    ``conditions`` are laid out in the given order; every one needs a unit. ``triples``
    (rows of increasing indices into ``conditions``) default to every combination. Each
    trait's units come from cellulose_decr's ``u6_units``; the units, and so the draws,
    are the same for every trait, and only the taxa each unit carries differ.
    """
    _, unit_study, _, frame = u6_units(profiles, taxa, set())
    built = [u6_units(profiles, taxa, trait_taxa(matrix, trait)) for trait in traits]
    for _, study, _, other in built:
        if not (other.equals(frame) and np.array_equal(study, unit_study)):
            raise ValueError("U6 units differ between traits")
    rows = frame.reset_index(drop=True)
    members = [np.flatnonzero(rows["condition"].eq(name).to_numpy()) for name in conditions]
    units = np.array([len(member) for member in members], dtype=np.int64)
    if (units == 0).any():
        missing = [name for name, count in zip(conditions, units, strict=True) if not count]
        raise ValueError(f"conditions without a stool decreased signature: {missing}")
    packed = [pack(matrix_bits) for matrix_bits, *_ in built]
    widths = np.array([words.shape[1] for words in packed], dtype=np.int64)
    offsets = np.concatenate([[0], np.cumsum(widths)[:-1]])[: len(widths)].astype(np.int64)
    unit_bits = np.hstack(packed) if packed else np.zeros((len(unit_study), 0), dtype=np.uint64)
    observed = np.zeros((1, len(conditions), unit_bits.shape[1]), dtype=np.uint64)
    for column, member in enumerate(members):
        observed[0, column] = np.bitwise_or.reduce(unit_bits[member], axis=0)
    # Filled row by row, so no eligible condition gives empty sets, not an error (#901).
    sets = [
        np.zeros((len(conditions), matrix_bits.shape[1]), dtype=bool) for matrix_bits, *_ in built
    ]
    for (matrix_bits, *_), condition_sets in zip(built, sets, strict=True):
        for column, member in enumerate(members):
            condition_sets[column] = matrix_bits[member].any(axis=0)
    n_studies = int(unit_study.max()) + 1 if unit_study.size else 0
    incidence = np.zeros((n_studies, len(conditions)), dtype=bool)
    for column, member in enumerate(members):
        incidence[unit_study[member], column] = True
    every = np.array(list(itertools.combinations(range(len(conditions)), 3)), dtype=np.int64)
    every = every.reshape(-1, 3)
    if triples is None:
        triples = every
    elif (
        triples.ndim != 2
        or triples.shape[1] != 3
        or (
            len(triples)
            and (
                (np.diff(triples, axis=1) <= 0).any()
                or triples.min() < 0
                or triples.max() >= len(conditions)
            )
        )
    ):
        raise ValueError("triples must be rows of three increasing indices into conditions")
    positions = block_positions(units.tolist(), n_studies)
    occupied = np.zeros((n_studies, len(conditions)), dtype=np.int64)
    for column, block in enumerate(positions):
        occupied[block, column] = 1
    overlap = (occupied.T @ occupied) > 0
    np.fill_diagonal(overlap, False)
    a, b, c = triples.T
    clash = overlap[a, b] | overlap[a, c] | overlap[b, c]
    fits = units[a] + units[b] + units[c] <= n_studies
    source = np.where(~fits, NOT_DRAWN, np.where(clash, OWN, SHARED))
    return Prepared(
        conditions=tuple(conditions),
        traits=tuple(traits),
        triples=triples,
        complete=bool(np.array_equal(triples, every)),
        units=units,
        unit_bits=unit_bits,
        offsets=offsets,
        unit_study=unit_study,
        observed=observed,
        incidence=incidence,
        sets=sets,
        u3_taxa=np.array([matrix_bits.shape[1] for matrix_bits, *_ in built], dtype=np.int64),
        positions=positions,
        source=source,
    )


def observed_values(prepared: Prepared) -> tuple[FloatArray, IntArray]:
    """Each triple's observed mean pairwise distance and core, ``(triples, traits)``."""
    if not prepared.traits or not len(prepared.triples):
        empty = (len(prepared.triples), len(prepared.traits))
        return np.full(empty, np.nan), np.zeros(empty, dtype=np.int64)
    if prepared.complete:
        return triple_values(prepared.observed, prepared.offsets)
    values = [
        triple_values(prepared.observed[:, row], prepared.offsets) for row in prepared.triples
    ]
    return (
        np.concatenate([mean for mean, _ in values]),
        np.concatenate([core for _, core in values]),
    )


def score_draws(
    totals: Totals,
    masks: WordArray,
    prepared: Prepared,
    observed_mean: FloatArray,
    observed_core: IntArray,
    tolerance: float,
) -> None:
    """Count one chunk of shared draws for every triple of the layout.

    Every combination is scored a block at a time; listed triples one at a time, with
    the same arithmetic.
    """
    if prepared.complete:
        for block, mean, core in triple_blocks(masks, prepared.offsets):
            totals.add(block, mean, core, observed_mean[block], observed_core[block], tolerance)
        return
    for index, row in enumerate(prepared.triples):
        target = slice(index, index + 1)
        for _, mean, core in triple_blocks(masks[:, row], prepared.offsets):
            totals.add(target, mean, core, observed_mean[target], observed_core[target], tolerance)


def run_null(
    prepared: Prepared, permutations: int, seed: SeedLike, tolerance: float
) -> tuple[FloatArray, IntArray, Totals]:
    """Observed statistics and their U6 totals, from shared draws and each own-draw triple's.

    The shared draws use ``seed``. A triple whose blocks overlap is drawn on its own with
    the seed ``[*seed, triple index]`` (``[seed, triple index]`` for an integer seed),
    laid out as (a, b, c).
    """
    observed_mean, observed_core = observed_values(prepared)
    totals = Totals.zeros((len(prepared.triples), len(prepared.traits)))
    if not prepared.traits or not len(prepared.triples):
        return observed_mean, observed_core, totals
    entropy = seed_entropy(seed)
    rng = np.random.default_rng(entropy)
    for unit_of_study, order in signature_draws(rng, prepared.unit_study, permutations):
        for start in range(0, order.shape[0], CHUNK):
            rows = slice(start, start + CHUNK)
            masks = draw_masks(
                prepared.unit_bits, unit_of_study[rows], order[rows], prepared.positions
            )
            score_draws(totals, masks, prepared, observed_mean, observed_core, tolerance)
    n_studies = int(prepared.unit_study.max()) + 1
    for index in np.flatnonzero(prepared.source != SHARED):
        for field in totals.__dataclass_fields__:
            getattr(totals, field)[index] = 0
        if prepared.source[index] == NOT_DRAWN:
            continue
        own = block_positions(prepared.units[prepared.triples[index]].tolist(), n_studies)
        rng = np.random.default_rng([*entropy, int(index)])
        row = slice(int(index), int(index) + 1)
        for unit_of_study, order in signature_draws(rng, prepared.unit_study, permutations):
            for start in range(0, order.shape[0], CHUNK):
                rows = slice(start, start + CHUNK)
                masks = draw_masks(prepared.unit_bits, unit_of_study[rows], order[rows], own)
                for _, mean, core in triple_blocks(masks, prepared.offsets):
                    totals.add(row, mean, core, observed_mean[row], observed_core[row], tolerance)
    return observed_mean, observed_core, totals


def scan(
    profiles: pd.DataFrame,
    taxa: pd.DataFrame,
    matrix: pd.DataFrame,
    traits: Sequence[str],
    conditions: Sequence[str],
    *,
    permutations: int,
    seed: SeedLike,
    tolerance: float,
    triples: IntArray | None = None,
) -> ScanResult:
    """Every triple of ``conditions`` (in that layout order) and trait, against U6.

    ``triples`` restricts the scan to listed triples (see ``prepare``); the layout, and
    so each triple's draws, still comes from all of ``conditions``.
    """
    prepared = prepare(profiles, taxa, matrix, traits, conditions, triples)
    observed_mean, observed_core, totals = run_null(prepared, permutations, seed, tolerance)
    a, b, c = prepared.triples.T
    union = np.zeros((len(prepared.triples), len(prepared.traits)), dtype=np.int64)
    sizes = np.zeros((len(prepared.conditions), len(prepared.traits)), dtype=np.int64)
    for column, sets in enumerate(prepared.sets):
        union[:, column] = (sets[a] | sets[b] | sets[c]).sum(axis=1)
        sizes[:, column] = sets.sum(axis=1)
    incidence = prepared.incidence
    pair_shared = (
        (incidence[:, a] & incidence[:, b])
        | (incidence[:, a] & incidence[:, c])
        | (incidence[:, b] & incidence[:, c])
    )
    return ScanResult(
        conditions=prepared.conditions,
        traits=prepared.traits,
        triples=prepared.triples,
        units=prepared.units,
        sizes=sizes,
        union=union,
        observed_mean=observed_mean,
        observed_core=observed_core,
        totals=totals,
        source=prepared.source,
        shared_studies=pair_shared.sum(axis=0).astype(np.int64),
        pool_studies=int(prepared.unit_study.max()) + 1 if prepared.unit_study.size else 0,
        u3_taxa=prepared.u3_taxa,
    )


# --------------------------------------------------------------------------- tables


def null_summary(result: ScanResult) -> dict[str, FloatArray]:
    """U6 expected value, SD, z and p for both statistics, ``(triples, traits)``."""
    totals = result.totals
    with np.errstate(invalid="ignore", divide="ignore"):
        defined = totals.distance_defined.astype(np.float64)
        expected_mean = totals.distance_sum / defined
        sd_mean = np.sqrt(
            np.maximum(totals.distance_squares - totals.distance_sum**2 / defined, 0)
            / (defined - 1)
        )
        draws = totals.draws.astype(np.float64)
        expected_core = totals.core_sum / draws
        sd_core = np.sqrt(
            np.maximum(totals.core_squares - totals.core_sum.astype(np.float64) ** 2 / draws, 0)
            / (draws - 1)
        )
        z_mean = np.where(sd_mean > 0, (result.observed_mean - expected_mean) / sd_mean, np.nan)
        z_core = np.where(sd_core > 0, (result.observed_core - expected_core) / sd_core, np.nan)
        # No draws (permutations=0) means no null, not p = 1.
        p_mean = np.where(draws > 0, (1 + totals.distance_at_most) / (1 + draws), np.nan)
        p_core = np.where(draws > 0, (1 + totals.core_at_least) / (1 + draws), np.nan)
    return {
        "u6_expected_mean_distance": expected_mean,
        "u6_sd_mean_distance": sd_mean,
        "u6_z_mean_distance": z_mean,
        "u6_p_mean_distance": p_mean,
        "u6_expected_core": expected_core,
        "u6_sd_core": sd_core,
        "u6_z_core": z_core,
        "u6_p_core": p_core,
    }


P_COLUMNS = {MEAN_DISTANCE: "u6_p_mean_distance", CORE: "u6_p_core"}
Z_COLUMNS = {MEAN_DISTANCE: "u6_z_mean_distance", CORE: "u6_z_core"}
RANK_COLUMNS = {MEAN_DISTANCE: "rank_mean_distance", CORE: "rank_core"}
# The secondary rank by U6 z (amendment of 2026-10-22): more overlap first, which is
# ascending z for the mean distance and descending z for the core.
Z_RANK_COLUMNS = {MEAN_DISTANCE: "z_rank_mean_distance", CORE: "z_rank_core"}
Z_ASCENDING = {MEAN_DISTANCE: True, CORE: False}
# The columns of ``triple_table``, in order; a test keeps them in step with it.
TRIPLE_COLUMNS = (
    "trait_id",
    "condition_a",
    "condition_b",
    "condition_c",
    "n_a",
    "n_b",
    "n_c",
    "u6_units_a",
    "u6_units_b",
    "u6_units_c",
    "shared_studies",
    "three_way_core",
    "three_way_union",
    MEAN_DISTANCE,
    "u6_expected_mean_distance",
    "u6_sd_mean_distance",
    "u6_z_mean_distance",
    "u6_p_mean_distance",
    "u6_expected_core",
    "u6_sd_core",
    "u6_z_core",
    "u6_p_core",
    "u6_draws",
    "u6_draws_without_distance",
    "u6_draws_source",
    "test_status",
    "is_lead_triple",
    *RANK_COLUMNS.values(),
    *Z_RANK_COLUMNS.values(),
)


def triple_table(result: ScanResult, lead: Sequence[str]) -> pd.DataFrame:
    """One row per triple and trait, in trait order and then triple order."""
    n_triples = len(result.triples)
    labels = np.array(result.conditions, dtype=object)
    a, b, c = result.triples.T
    summary = null_summary(result)
    lead_set = set(lead)
    is_lead = np.array(
        [set(labels[row]) == lead_set for row in result.triples], dtype=bool
    ).reshape(n_triples)
    frames = []
    for column, trait in enumerate(result.traits):
        status = np.where(
            result.source == NOT_DRAWN,
            POOL_TOO_SMALL,
            np.where(np.isnan(result.observed_mean[:, column]), NO_TAXA, TESTED),
        )
        tested = status == TESTED
        frame = pd.DataFrame(
            {
                "trait_id": trait,
                "condition_a": labels[a],
                "condition_b": labels[b],
                "condition_c": labels[c],
                "n_a": result.sizes[a, column],
                "n_b": result.sizes[b, column],
                "n_c": result.sizes[c, column],
                "u6_units_a": result.units[a],
                "u6_units_b": result.units[b],
                "u6_units_c": result.units[c],
                "shared_studies": result.shared_studies,
                "three_way_core": result.observed_core[:, column],
                "three_way_union": result.union[:, column],
                MEAN_DISTANCE: result.observed_mean[:, column],
                **{
                    name: np.where(tested, values[:, column], np.nan)
                    for name, values in summary.items()
                },
                "u6_draws": np.where(tested, result.totals.draws[:, column], 0),
                "u6_draws_without_distance": np.where(
                    tested,
                    result.totals.draws[:, column] - result.totals.distance_defined[:, column],
                    0,
                ),
                "u6_draws_source": result.source,
                "test_status": status,
                "is_lead_triple": is_lead,
            }
        )
        for statistic, p_column in P_COLUMNS.items():
            ranks = frame.loc[tested, p_column].rank(method="min")
            frame[RANK_COLUMNS[statistic]] = ranks.reindex(frame.index).astype("Int64")
        for statistic, z_column in Z_COLUMNS.items():
            z = frame.loc[tested, z_column].dropna()
            z_ranks = z.rank(method="min", ascending=Z_ASCENDING[statistic])
            frame[Z_RANK_COLUMNS[statistic]] = z_ranks.reindex(frame.index).astype("Int64")
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=TRIPLE_COLUMNS)
    return pd.concat(frames, ignore_index=True)


def reference(frame: pd.DataFrame, name: str) -> pd.DataFrame:
    """The tested triples of one reference set."""
    tested = frame[frame["test_status"].eq(TESTED)]
    return tested if name == ALL_TRIPLES else tested[tested["shared_studies"].eq(0)]


def lead_rank(p_values: pd.Series, lead: float, *, counted: bool = True) -> dict[str, Any]:
    """Where the lead triple's p falls among the reference triples, smallest first.

    Ties make the rank a range: ``rank_min`` is one more than the triples with a smaller
    p, and ``rank_max`` counts those whose p is at most the lead's. When the lead is not
    among ``p_values`` (``counted`` false: it shares a study and the reference is the
    study-disjoint triples), it is placed among them as one more triple.
    """
    extra = 0 if counted else 1
    at_most = int((p_values <= lead).sum()) + extra
    triples = len(p_values) + extra
    return {
        "triples": triples,
        "rank_min": int((p_values < lead).sum()) + 1,
        "rank_max": at_most,
        "share_at_most": at_most / triples,
    }


def bin_counts(p_values: pd.Series, bins: int) -> list[int]:
    """Counts of p in (0, 1/bins], (1/bins, 2/bins], ... (closed above, so 0.05 is bin 1)."""
    upper = np.linspace(0, 1, bins + 1)[1:]
    index = np.searchsorted(upper, p_values.to_numpy(dtype=float), side="left")
    counts: list[int] = np.bincount(np.minimum(index, bins - 1), minlength=bins).tolist()
    return counts


Z_QUANTILES = (0.05, 0.25, 0.5, 0.75, 0.95)


def describe(p_values: pd.Series, z_values: pd.Series) -> dict[str, Any]:
    """Quantiles of the triples' U6 p-values and z-scores, and the shares at 0.05 and 0.01."""
    if p_values.empty:
        return {"triples": 0}
    z = z_values.dropna()
    return {
        "triples": len(p_values),
        "p_quantiles": dict(
            zip(map(str, QUANTILES), np.quantile(p_values, QUANTILES), strict=True)
        ),
        "share_p_at_most_0.05": float((p_values <= 0.05).mean()),
        "share_p_at_most_0.01": float((p_values <= 0.01).mean()),
        "z_quantiles": (
            dict(zip(map(str, Z_QUANTILES), np.quantile(z, Z_QUANTILES), strict=True))
            if len(z)
            else None
        ),
    }


LEAD_COLUMNS = (
    "condition_a",
    "condition_b",
    "condition_c",
    "n_a",
    "n_b",
    "n_c",
    "shared_studies",
    "three_way_core",
    "three_way_union",
    MEAN_DISTANCE,
    "u6_expected_mean_distance",
    "u6_z_mean_distance",
    "u6_p_mean_distance",
    "u6_expected_core",
    "u6_z_core",
    "u6_p_core",
    "test_status",
)


def lead_and_distribution(
    triples: pd.DataFrame, traits: Sequence[str], bins: int
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Per trait: the lead triple's results and ranks, the binned p-values, and quantiles.

    Each statistic is ranked over two reference sets: every tested triple, and the
    tested triples whose conditions share no study (the declared sensitivity). The lead
    is ranked by U6 p and, by the amendment of 2026-10-22, by U6 z (``*_z_*`` fields),
    among the triples whose z is defined.
    """
    rows: list[dict[str, Any]] = []
    histogram: list[dict[str, Any]] = []
    quantiles: dict[str, Any] = {}
    edges = np.linspace(0, 1, bins + 1)
    for trait in traits:
        mine = triples[triples["trait_id"].eq(trait)]
        lead = mine[mine["is_lead_triple"]]
        row: dict[str, Any] = {"trait_id": trait, "lead_in_scan": not lead.empty}
        if not lead.empty:
            row.update({f"lead_{column}": lead[column].iloc[0] for column in LEAD_COLUMNS})
        tested = not lead.empty and lead["test_status"].iloc[0] == TESTED
        quantiles[trait] = {}
        for statistic, p_column in P_COLUMNS.items():
            quantiles[trait][statistic] = {}
            for name in REFERENCE_SETS:
                chosen = reference(mine, name)
                if tested:
                    ranked = lead_rank(
                        chosen[p_column],
                        float(lead[p_column].iloc[0]),
                        counted=bool(chosen["is_lead_triple"].any()),
                    )
                    row.update(
                        {f"{statistic}_{name}_{key}": value for key, value in ranked.items()}
                    )
                    z_column = Z_COLUMNS[statistic]
                    lead_z = float(lead[z_column].iloc[0])
                    if not np.isnan(lead_z):
                        # lead_rank ranks small values first; flip z where large is first.
                        sign = 1.0 if Z_ASCENDING[statistic] else -1.0
                        defined = chosen[chosen[z_column].notna()]
                        z_ranked = lead_rank(
                            sign * defined[z_column],
                            sign * lead_z,
                            counted=bool(defined["is_lead_triple"].any()),
                        )
                        row.update(
                            {
                                f"{statistic}_{name}_z_{key}": value
                                for key, value in z_ranked.items()
                            }
                        )
                for index, count in enumerate(bin_counts(chosen[p_column], bins)):
                    histogram.append(
                        {
                            "trait_id": trait,
                            "statistic": statistic,
                            "reference_set": name,
                            "p_above": edges[index],
                            "p_at_most": edges[index + 1],
                            "triples": count,
                        }
                    )
                quantiles[trait][statistic][name] = describe(
                    chosen[p_column], chosen[Z_COLUMNS[statistic]]
                )
        rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(histogram), quantiles


def trait_table(
    chosen: pd.DataFrame, lead: pd.DataFrame, result: ScanResult, lead_trait: str
) -> pd.DataFrame:
    """The trait set with its screen statistics, U3 size and the lead triple's results.

    ``rank_among_traits_*`` ranks the traits by the lead triple's U6 p (smallest first),
    by the share of tested triples whose p is at most the lead's, and by the share whose
    z is at least as extreme, so cellulose can be read against the other top decreasing
    traits.
    """
    u3 = dict(zip(result.traits, result.u3_taxa.tolist(), strict=True))
    frame = chosen.assign(
        is_lead_trait=chosen["trait_id"].eq(lead_trait), u3_taxa=chosen["trait_id"].map(u3)
    ).merge(lead, on="trait_id", how="left", validate="one_to_one")
    for statistic, p_column in P_COLUMNS.items():
        for basis, column in (
            ("lead_p", f"lead_{p_column}"),
            ("lead_share", f"{statistic}_{ALL_TRIPLES}_share_at_most"),
            ("lead_z_share", f"{statistic}_{ALL_TRIPLES}_z_share_at_most"),
        ):
            values = frame[column] if column in frame else pd.Series(np.nan, index=frame.index)
            frame[f"rank_among_traits_by_{basis}_{statistic}"] = values.rank(method="min").astype(
                "Int64"
            )
    return frame


# --------------------------------------------------------------------------- run


def read_matrix(path: Path, traits: Sequence[str]) -> pd.DataFrame:
    """The matrix's taxon rows with the chosen trait columns only."""
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    missing = sorted(set(traits) - set(header))
    if missing:
        raise ValueError(f"traits missing from {path}: {missing}")
    matrix: pd.DataFrame = pd.read_csv(path, sep="\t", index_col=0, usecols=[header[0], *traits])
    return matrix


def time_draws(prepared: Prepared, tolerance: float) -> dict[str, float]:
    """Seconds per draw for the shared draws and for one own-draw triple, from a burst."""
    if not prepared.traits or not len(prepared.triples):
        return {SHARED: 0.0}
    observed_mean, observed_core = observed_values(prepared)
    totals = Totals.zeros((len(prepared.triples), len(prepared.traits)))
    n_studies = int(prepared.unit_study.max()) + 1
    first = block_positions(prepared.units[prepared.triples[0]].tolist(), n_studies)
    seconds: dict[str, float] = {}
    for name, positions in ((SHARED, prepared.positions), (OWN, first)):
        started = time.perf_counter()
        for unit_of_study, order in signature_draws(
            np.random.default_rng(0), prepared.unit_study, CHUNK
        ):
            masks = draw_masks(prepared.unit_bits, unit_of_study, order, positions)
            if name == SHARED:
                score_draws(totals, masks, prepared, observed_mean, observed_core, tolerance)
                continue
            for _, mean, core in triple_blocks(masks, prepared.offsets):
                row = slice(0, 1)
                totals.add(row, mean, core, observed_mean[row], observed_core[row], tolerance)
        seconds[name] = (time.perf_counter() - started) / CHUNK
    return seconds


def plan(
    conditions: pd.DataFrame,
    chosen: pd.DataFrame,
    prepared: Prepared,
    protocol: Protocol,
    permutations: int,
) -> dict[str, Any]:
    """What a full run would draw, and a timed estimate of how long it takes."""
    seconds = time_draws(prepared, protocol.tolerance)
    own = int((prepared.source == OWN).sum())
    estimate = permutations * (seconds[SHARED] + own * seconds.get(OWN, 0.0))
    return {
        "eligible_conditions": len(prepared.conditions),
        "conditions_with_stool_comparisons": len(conditions),
        "triples": len(prepared.triples),
        "traits": chosen[["trait_id", "trait_label", "decreasing_rank", "in_top_n"]].to_dict(
            "records"
        ),
        "pool_studies": int(prepared.unit_study.max()) + 1 if prepared.unit_study.size else 0,
        "eligible_u6_units": int(prepared.units.sum()),
        "own_draw_triples": own,
        "untestable_triples": int((prepared.source == NOT_DRAWN).sum()),
        "lead_conditions_eligible": {
            label: label in prepared.conditions for label in protocol.lead_triple
        },
        "permutations": permutations,
        "seconds_per_shared_draw": seconds[SHARED],
        "seconds_per_own_triple_draw": seconds.get(OWN),
        "estimated_null_minutes": estimate / 60,
    }


# --------------------------------------------------------------------------- study splits (#854)


@dataclass(frozen=True)
class ScanInputs:
    """One core population, read once and restricted to any set of its studies.

    ``profiles`` and ``taxa`` are every selected experiment's profiles and signature taxa
    (``build_population``). ``matrix`` holds every screen candidate trait. ``screen``
    holds the paired human profiles the directional screen votes on, over ``columns``,
    and ``labels`` name the traits. ``contrast_population`` must be the Tier 1+2
    population; ``strain_rollup`` is the core's scope and ``strain_rollup_role`` its
    declared role (``scope_role``), which ``select`` and ``evaluate`` check against the
    protocol (#907).
    """

    profiles: pd.DataFrame
    taxa: pd.DataFrame
    matrix: pd.DataFrame
    screen: tuple[Profile, ...]
    columns: tuple[str, ...]
    labels: Mapping[str, ColumnLabel]
    contrast_population: str
    strain_rollup: str
    strain_rollup_role: str

    def __post_init__(self) -> None:
        if self.contrast_population != TIER1_2:
            raise ValueError(
                f"the scan's population is {TIER1_2}, not {self.contrast_population!r}"
            )
        if self.strain_rollup_role not in (PRIMARY, SENSITIVITY):
            raise ValueError(f"unknown strain-rollup role {self.strain_rollup_role!r}")

    @property
    def studies(self) -> frozenset[str]:
        """Every study of the population; a split divides these."""
        return frozenset(self.profiles["study_id"])

    def restricted(self, studies: Iterable[str]) -> ScanInputs:
        """The population without the studies outside ``studies``.

        Whole studies are kept or dropped, with every experiment, profile and signature
        taxon they hold. A study that is not in the population is an error.
        """
        keep = frozenset(studies)
        unknown = sorted(keep - self.studies)
        if unknown:
            raise ValueError(f"{len(unknown)} studies are not in the population: {unknown[:3]}")
        return dataclasses.replace(
            self,
            profiles=self.profiles[self.profiles["study_id"].isin(keep)],
            taxa=self.taxa[self.taxa["study_id"].isin(keep)],
            screen=tuple(profile for profile in self.screen if profile.study_id in keep),
        )


def load_inputs(
    paths: Mapping[str, Path],
    *,
    contrast_population: str,
    strain_rollup: str,
    protocol: Protocol | None = None,
) -> ScanInputs:
    """Read a core run's ``CORE_PRODUCTS`` once, for the scan and for any study split.

    The core's population and strain-rollup scope must be stated; a population other
    than Tier 1+2 or an undeclared scope is refused. ``load_core_inputs`` reads both
    from the verified run.
    """
    role = scope_role(strain_rollup, protocol or load_protocol())
    if contrast_population != TIER1_2:
        raise ValueError(f"the scan's population is {TIER1_2}, not {contrast_population!r}")
    columns = screen_candidates(paths["count_profiles"])
    screen, _ = load_profiles(
        paths["count_profiles"],
        paths["experiments"],
        candidate_slots(columns),
        reject_missing_counts=True,
    )
    matrix = read_matrix(paths["taxa_matrix"], columns)
    profiles, taxa = load_population(
        paths["experiments"], paths["reported_taxa"], paths["count_profiles"], matrix.index
    )
    labels = load_column_labels(paths["feature_metadata"])
    return ScanInputs(
        profiles,
        taxa,
        matrix,
        tuple(screen),
        tuple(columns),
        labels,
        contrast_population=contrast_population,
        strain_rollup=strain_rollup,
        strain_rollup_role=role,
    )


def require_scope(inputs: ScanInputs, protocol: Protocol) -> None:
    """Refuse inputs of an undeclared scope, or whose role is not the protocol's (#907)."""
    role = scope_role(inputs.strain_rollup, protocol)
    if role != inputs.strain_rollup_role:
        raise ValueError(
            f"strain_rollup {inputs.strain_rollup!r} is the protocol's {role}, "
            f"not {inputs.strain_rollup_role}"
        )


def open_tier1_2(core_run: Path, protocol: Protocol | None = None) -> tuple[CoreBundle, str, str]:
    """Verify a core run; return it with its strain-rollup scope and where the scope came
    from (``core_scope``).

    A core of another population, direction basis or undeclared scope is refused.
    """
    from mtp.workflow.consumers import open_core_bundle
    from mtp.workflow.population import ORIENTED_DIRECTION_BASIS

    bundle = open_core_bundle(core_run)
    if (
        bundle.run.contrast_population != TIER1_2
        or bundle.run.direction_basis != ORIENTED_DIRECTION_BASIS
    ):
        raise ValueError(
            f"the scan reads a {TIER1_2} core oriented {ORIENTED_DIRECTION_BASIS}; "
            f"{core_run} is {bundle.run.contrast_population}, {bundle.run.direction_basis}"
        )
    protocol = protocol or load_protocol()
    scope, source = core_scope(bundle.run, protocol)
    scope_role(scope, protocol)  # refuses an undeclared scope
    return bundle, scope, source


def load_core_inputs(
    core_run: Path, protocol: Protocol | None = None
) -> tuple[CoreBundle, ScanInputs]:
    """Verify a Tier 1+2 core run (its ``run.json``) of a declared scope; read it once."""
    bundle, scope, _ = open_tier1_2(core_run, protocol)
    paths = {product: bundle.path(product) for product in CORE_PRODUCTS}
    inputs = load_inputs(
        paths,
        contrast_population=str(bundle.run.contrast_population),
        strain_rollup=scope,
        protocol=protocol,
    )
    return bundle, inputs


def declared_sets(
    inputs: ScanInputs, protocol: Protocol, always_included: bool = True
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """The declared rule's trait set, condition table and eligible conditions.

    ``always_included`` false leaves the protocol's always-included trait out unless
    the screen of these studies selects it.
    """
    screen = screen_table(inputs.screen, inputs.columns, inputs.labels)
    always = protocol.always_included if always_included else None
    traits = choose_traits(screen, protocol.top_n, always)
    conditions = condition_table(inputs.profiles, protocol.min_studies)
    eligible = sorted(conditions.loc[conditions["eligible"], "condition"])
    return traits, conditions, eligible


@dataclass(frozen=True)
class Selection:
    """The declared rule applied to one set of studies.

    ``conditions`` is ``condition_table`` on those studies; ``traits`` the trait set,
    chosen by the directional screen of those studies; ``triples`` one row per triple of
    eligible conditions and trait (``triple_table``), with U6 drawn from those studies.
    ``result`` holds the arrays behind ``triples``.
    """

    studies: frozenset[str]
    conditions: pd.DataFrame
    traits: pd.DataFrame
    triples: pd.DataFrame
    result: ScanResult

    @property
    def status(self) -> str:
        """Whether the family has triples: ``selected``, or why not.

        ``no_traits`` when the screen of these studies selects no trait, else
        ``fewer_than_three_eligible_conditions`` when no triple can be formed. An empty
        family is a result of the declared rule, not a failed run (#901).
        """
        if self.traits.empty:
            return NO_TRAITS
        if len(self.result.conditions) < 3:
            return TOO_FEW_CONDITIONS
        return SELECTED


def select(
    inputs: ScanInputs,
    studies: Iterable[str] | None = None,
    *,
    protocol: Protocol | None = None,
    permutations: int | None = None,
    seed: SeedLike | None = None,
    always_included: bool = True,
) -> Selection:
    """Apply the declared rule (eligibility, trait set, U6 scan) to ``studies`` only.

    ``studies`` defaults to every study of the population, which is the scan the CLI
    writes. For a split, pass the training studies: nothing outside them is read, for
    eligibility, the trait screen, the observed sets or the U6 pool. ``always_included``
    false lets the training screen alone choose the traits, so the historical lead trait
    is not forced into a split's selection. ``permutations`` and ``seed`` default to the
    protocol's.
    """
    protocol = protocol or load_protocol()
    require_scope(inputs, protocol)
    chosen = inputs if studies is None else inputs.restricted(studies)
    traits, conditions, eligible = declared_sets(chosen, protocol, always_included)
    result = scan(
        chosen.profiles,
        chosen.taxa,
        chosen.matrix,
        traits["trait_id"].tolist(),
        eligible,
        permutations=protocol.permutations if permutations is None else permutations,
        seed=protocol.seed if seed is None else seed,
        tolerance=protocol.tolerance,
    )
    triples = triple_table(result, protocol.lead_triple)
    return Selection(chosen.studies, conditions, traits, triples, result)


KEY_COLUMNS = ["trait_id", "condition_a", "condition_b", "condition_c"]
# Explicit held-out outcomes (#854), beside ``test_status``.
OUTCOME_COLUMNS = ("missing_condition", "shares_studies", "empty_set", "null_undefined")


def evaluate(
    inputs: ScanInputs,
    studies: Iterable[str],
    hypotheses: Iterable[Sequence[str]],
    *,
    protocol: Protocol | None = None,
    permutations: int | None = None,
    seed: SeedLike | None = None,
) -> pd.DataFrame:
    """U6 tests of hypotheses fixed beforehand, on ``studies`` only (a held-out half).

    Each hypothesis is (condition, condition, condition, trait ID), for example a row
    ``select`` chose on the training studies. Only ``studies`` are read: the sets, the U6
    pool and its draws. No eligibility rule is applied and no hypothesis is replaced. One
    whose condition has no stool decreased signature among ``studies`` is returned with
    ``test_status`` ``no_units``.

    Returns one row per hypothesis in the given order, numbered by ``hypothesis``, with
    the conditions sorted and the columns of ``triple_table`` except the within-trait
    ranks, then four outcome flags (``OUTCOME_COLUMNS``):

    - ``missing_condition``: a condition has no held-out stool decreased signature
      (``test_status`` ``no_units``);
    - ``shares_studies``: a held-out study reports two of the conditions, which U6's
      distinct studies cannot reproduce;
    - ``empty_set``: a condition's held-out set has no taxon bearing the trait;
    - ``null_undefined``: no U6 draw gave a distance (not drawn, or every draw empty).

    Set sizes and shared studies are given for every hypothesis, drawn or not. The
    involved conditions, sorted, are the layout of the shared draws.
    """
    protocol = protocol or load_protocol()
    require_scope(inputs, protocol)
    held_out = inputs.restricted(studies)
    asked: list[dict[str, Any]] = []
    for number, hypothesis in enumerate(hypotheses):
        if len(hypothesis) != 4 or len(set(hypothesis[:3])) != 3:
            raise ValueError(f"hypothesis {number} is not three conditions and a trait")
        if hypothesis[3] not in inputs.columns:
            raise ValueError(f"hypothesis {number} names a trait that is not a screen candidate")
        labels = sorted(hypothesis[:3])
        asked.append(
            {"hypothesis": number, **dict(zip(KEY_COLUMNS, [hypothesis[3], *labels], strict=True))}
        )
    frame = pd.DataFrame(asked, columns=["hypothesis", *KEY_COLUMNS])
    units = condition_table(held_out.profiles, 1).set_index("condition")["u6_units"]
    drawable = sorted(
        {label for label in frame[KEY_COLUMNS[1:]].to_numpy().ravel() if units.get(label, 0)}
    )
    position = {label: index for index, label in enumerate(drawable)}
    triples = sorted(
        {
            tuple(position[label] for label in row)
            for row in frame[KEY_COLUMNS[1:]].itertuples(index=False, name=None)
            if all(label in position for label in row)
        }
    )
    result = scan(
        held_out.profiles,
        held_out.taxa,
        held_out.matrix,
        list(dict.fromkeys(frame["trait_id"])),
        drawable,
        permutations=protocol.permutations if permutations is None else permutations,
        seed=protocol.seed if seed is None else seed,
        tolerance=protocol.tolerance,
        triples=np.array(triples, dtype=np.int64).reshape(-1, 3),
    )
    table = triple_table(result, protocol.lead_triple)
    ranks = [*RANK_COLUMNS.values(), *Z_RANK_COLUMNS.values()]
    table = table.drop(columns=ranks, errors="ignore")
    merged = frame.merge(table, on=KEY_COLUMNS, how="left", validate="many_to_one")
    # A hypothesis the scan did not draw has a condition without held-out units.
    drawn = merged["test_status"].notna()
    merged["test_status"] = merged["test_status"].where(drawn, NO_UNITS)
    merged["u6_draws_source"] = merged["u6_draws_source"].where(drawn, NOT_DRAWN)
    for letter in "abc":
        held = merged[f"condition_{letter}"].map(units).fillna(0)
        merged[f"u6_units_{letter}"] = merged[f"u6_units_{letter}"].where(drawn, held)
    lead = set(protocol.lead_triple)
    merged["is_lead_triple"] = [
        set(labels) == lead for labels in merged[KEY_COLUMNS[1:]].itertuples(index=False, name=None)
    ]
    # Sizes and shared studies straight from the held-out studies, for every hypothesis;
    # for a drawn one they are the scan's own.
    stool = held_out.taxa[held_out.taxa["compartment"].eq(FECES) & held_out.taxa["dir"].eq("dec")]
    reported = stool.groupby("condition")["taxon_id"].agg(set).to_dict()
    decreased = held_out.profiles[
        held_out.profiles["compartment"].eq(FECES) & held_out.profiles["dir"].eq("dec")
    ]
    studies_of = decreased.groupby("condition")["study_id"].agg(set).to_dict()
    bearing = {
        trait: trait_taxa(held_out.matrix, trait) for trait in dict.fromkeys(frame["trait_id"])
    }
    for letter in "abc":
        merged[f"n_{letter}"] = [
            len(reported.get(label, set()) & bearing[trait])
            for label, trait in zip(merged[f"condition_{letter}"], merged["trait_id"], strict=True)
        ]
    merged["shared_studies"] = [
        sum(
            sum(study in studies_of.get(label, set()) for label in labels) >= 2
            for study in set().union(*(studies_of.get(label, set()) for label in labels))
        )
        for labels in merged[KEY_COLUMNS[1:]].itertuples(index=False, name=None)
    ]
    counts = ["n_a", "n_b", "n_c", "u6_units_a", "u6_units_b", "u6_units_c", "shared_studies"]
    counts += ["three_way_core", "three_way_union", "u6_draws", "u6_draws_without_distance"]
    merged[counts] = merged[counts].astype("Int64")
    merged["missing_condition"] = merged["test_status"].eq(NO_UNITS)
    merged["shares_studies"] = merged["shared_studies"].gt(0)
    merged["empty_set"] = merged[["n_a", "n_b", "n_c"]].eq(0).any(axis=1)
    with_distance = merged["u6_draws"].fillna(0) - merged["u6_draws_without_distance"].fillna(0)
    merged["null_undefined"] = with_distance.le(0)
    merged[list(OUTCOME_COLUMNS)] = merged[list(OUTCOME_COLUMNS)].astype(bool)
    return merged[["hypothesis", *table.columns, *OUTCOME_COLUMNS]]


def closest_disjoint_triples(triples: pd.DataFrame, tolerance: float = 0.0) -> pd.DataFrame:
    """#854's selection rule: per trait, the closest study-disjoint triple.

    Among the tested triples whose conditions share no study (``shared_studies`` 0),
    the one with the smallest observed mean pairwise Jaccard distance. A tie goes to the
    first in lexical order of the three (sorted) conditions. Ties are exact by default,
    as #854 declares; a positive ``tolerance`` also ties distances that close. A trait
    without such a triple selects nothing. U6 is not used, so
    ``select(..., permutations=0)`` is enough for it.
    """
    usable = triples[triples["test_status"].eq(TESTED) & triples["shared_studies"].eq(0)]
    rows = []
    for _, group in usable.groupby("trait_id", sort=False):
        tied = group[group[MEAN_DISTANCE] <= group[MEAN_DISTANCE].min() + tolerance]
        rows.append(tied.sort_values(KEY_COLUMNS[1:], kind="mergesort").iloc[:1])
    chosen: pd.DataFrame = pd.concat(rows, ignore_index=True) if rows else triples.iloc[:0]
    return chosen


def rank_candidates(triples: pd.DataFrame, statistic: str = MEAN_DISTANCE) -> pd.DataFrame:
    """Tested triple-and-trait rows in a fixed order of evidence, strongest first.

    The order is U6 p ascending, then z (lower for the mean distance, higher for the
    core: more overlap first), then trait ID and the three conditions. It is an order,
    not a selection: how many rows to take is the caller's declared rule.
    """
    p_column, z_column = P_COLUMNS[statistic], Z_COLUMNS[statistic]
    tested = triples[triples["test_status"].eq(TESTED)]
    sign = 1 if statistic == MEAN_DISTANCE else -1
    ordered = tested.assign(_z=sign * tested[z_column]).sort_values(
        [p_column, "_z", *KEY_COLUMNS], na_position="last", kind="mergesort"
    )
    ordered = ordered.drop(columns="_z").reset_index(drop=True)
    ordered.insert(0, "candidate_rank", np.arange(1, len(ordered) + 1))
    return ordered


def main(argv: Sequence[str] | None = None) -> int:
    from mtp.workflow.consumers import preflight_consumer_outputs, write_consumer_record
    from mtp.workflow.population import ORIENTED_DIRECTION_BASIS

    args = parse_args(argv)
    protocol = load_protocol()
    permutations = protocol.permutations if args.permutations is None else args.permutations
    seed = protocol.seed if args.seed is None else args.seed
    out = args.output_dir
    if not args.plan and out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(f"{out} must be a new or empty directory")
    bundle, scope, source = open_tier1_2(args.core_run, protocol)
    rollup = scope_record(scope, source, protocol)
    outputs = [out / name for name in OUTPUTS.values()]
    if not args.plan:
        out.mkdir(parents=True, exist_ok=True)
        preflight_consumer_outputs(bundle, out / CONSUMER_RECORD, outputs, output_directory=out)
    paths = {product: bundle.path(product) for product in CORE_PRODUCTS}
    # Paths relative to the core run, so the record names no local directory.
    hashes = {
        product: {
            "path": str(path.resolve().relative_to(bundle.manifest.parent)),
            "sha256": file_sha256(path),
            "bytes": path.stat().st_size,
        }
        for product, path in paths.items()
    }

    inputs = load_inputs(
        paths,
        contrast_population=str(bundle.run.contrast_population),
        strain_rollup=scope,
        protocol=protocol,
    )
    if args.plan:
        chosen, conditions, eligible = declared_sets(inputs, protocol)
        traits = chosen["trait_id"].tolist()
        prepared = prepare(inputs.profiles, inputs.taxa, inputs.matrix, traits, eligible)
        report = {
            "protocol": protocol.version,
            STRAIN_ROLLUP: rollup,
            **plan(conditions, chosen, prepared, protocol, permutations),
        }
        print(json.dumps(clean(report), indent=1, ensure_ascii=False))
        return 0

    selection = select(inputs, protocol=protocol, permutations=permutations, seed=seed)
    result, triples, chosen = selection.result, selection.triples, selection.traits
    conditions = selection.conditions.copy()
    traits = chosen["trait_id"].tolist()
    eligible = list(result.conditions)
    lead, histogram, quantiles = lead_and_distribution(triples, traits, protocol.bins)
    trait_frame = trait_table(chosen, lead, result, protocol.lead_trait)
    layout = {label: position for position, label in enumerate(eligible, start=1)}
    conditions["layout_position"] = conditions["condition"].map(layout).astype("Int64")
    conditions["is_lead_condition"] = conditions["condition"].isin(protocol.lead_triple)
    follows = permutations == protocol.permutations and seed == protocol.seed
    notes = [
        "Every p-value is nominal; no multiplicity correction is applied here.",
        "One sequence of U6 draws serves every triple and trait, so p-values share draws.",
        "Counts describe curated signature taxa and their annotations, not physiology.",
    ]
    if protocol.side_by_side:
        seen = " and ".join(protocol.results_seen)
        notes.append(
            f"The strain-rollup roles were declared after the {seen} results had been seen "
            f"({protocol.version}), a post-results choice of primary scope. Report this "
            f"{scope} scan only beside the {', '.join(rollup['report_beside'])} scan of the "
            "same protocol (just cellulose-decr-selection-scan-scopes)."
        )
    summary = {
        "protocol": protocol.version,
        "follows_protocol": follows,
        "rule": {
            "min_studies": protocol.min_studies,
            "top_n": protocol.top_n,
            "always_included": protocol.always_included,
            "permutations": permutations,
            "seed": seed,
            "protocol_permutations": protocol.permutations,
            "protocol_seed": protocol.seed,
            "distance_tolerance": protocol.tolerance,
            "primary_statistic": MEAN_DISTANCE,
            "secondary_statistic": CORE,
            "secondary_rank": "u6_z",
        },
        "core_run": {
            "run_id": bundle.run.run_id,
            "contrast_population": TIER1_2,
            "direction_basis": ORIENTED_DIRECTION_BASIS,
        },
        STRAIN_ROLLUP: rollup,
        "counts": {
            "conditions_with_stool_comparisons": len(conditions),
            "eligible_conditions": len(eligible),
            "triples": len(result.triples),
            "traits": len(traits),
            "pool_studies": result.pool_studies,
            "eligible_u6_units": int(result.units.sum()),
            "shared_draw_triples": int((result.source == SHARED).sum()),
            "own_draw_triples": int((result.source == OWN).sum()),
            "untestable_triples": int((result.source == NOT_DRAWN).sum()),
            "study_disjoint_triples": int((result.shared_studies == 0).sum()),
        },
        "lead_triple": {
            "conditions": list(protocol.lead_triple),
            "eligible": {label: label in layout for label in protocol.lead_triple},
        },
        "lead_trait": protocol.lead_trait,
        "traits": trait_frame.to_dict("records"),
        "p_value_distribution": quantiles,
        "notes": notes,
    }
    write_tsv(conditions, out / OUTPUTS["conditions"])
    write_tsv(trait_frame, out / OUTPUTS["traits"])
    triples.to_csv(
        out / OUTPUTS["triples"],
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
        lineterminator="\n",
        compression={"method": "gzip", "mtime": 0},
    )
    write_tsv(histogram, out / OUTPUTS["distribution"])
    write_json(summary, out / OUTPUTS["summary"])
    modules = sorted(Path(__file__).resolve().parent.glob("*.py"))
    provenance = {
        "generated_at_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "code_commit": args.code_commit,
        "script": "src/mtp/cellulose_decr/selection_scan.py",
        "scripts": {f"src/mtp/cellulose_decr/{path.name}": file_sha256(path) for path in modules},
        "protocol": protocol_record(protocol.version),
        "methods": {
            "path": repository_path(METHODS_PATH),
            "sha256": file_sha256(METHODS_PATH),
        },
        "core_run": {
            "manifest_sha256": file_sha256(bundle.manifest),
            "run_id": bundle.run.run_id,
            "git_commit": bundle.run.git_commit,
            "git_dirty": bundle.run.git_dirty,
            "contrast_population": TIER1_2,
            "direction_basis": ORIENTED_DIRECTION_BASIS,
            "parameters": {item.name: item.value for item in bundle.run.parameters or []},
            "sources": {
                artifact.product_id: artifact.sha256
                for artifact in bundle.run.artifacts or []
                if artifact.product_id in CORE_SOURCES
            },
        },
        STRAIN_ROLLUP: rollup,
        "inputs": hashes,
        "permutations": permutations,
        "seed": seed,
        "follows_protocol": follows,
        "chunk": CHUNK,
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
        "outputs": sorted([*OUTPUTS.values(), CONSUMER_RECORD]),
    }
    write_json(provenance, out / OUTPUTS["provenance"])
    write_consumer_record(bundle, out / CONSUMER_RECORD, consumer=CONSUMER, outputs=outputs)
    print(f"wrote {len(outputs) + 1} files to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
