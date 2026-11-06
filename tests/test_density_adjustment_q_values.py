"""Density adjustment labels follow BH q-values over the reported traits (#1321).

Raw and adjusted Wilcoxon p-values are corrected separately, each over every trait
the run reports: the six Table 1 traits by default, or the declared list. A trait
with ten or fewer non-zero pairs stays in the family at p = 1 and exports no q.
"""

from __future__ import annotations

import json
import math
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control

import mtp.density_adjustment as density_module
from mtp.density_adjustment import (
    TABLE_1_TRAITS,
    TRAIT_LIST_PROVENANCE_FILENAME,
    adjustment_status,
    evaluate_table1_paired_effects,
    run_density_audit,
)
from mtp.trait_lists import load_trait_list

_STEPS = np.arange(1.0, 16.0)
# Ranks 1-30 with ranks 1-14 and 22 negative: a two-sided Wilcoxon p of about 0.029.
_BORDERLINE = np.arange(1.0, 31.0) * np.where(
    np.isin(np.arange(1, 31), [*range(1, 15), 22]), -1.0, 1.0
)


def _strong(mean: float) -> np.ndarray:
    """Thirty same-signed differences averaging ``mean``; Wilcoxon p is about 1.9e-9."""
    return mean + np.concatenate([_STEPS, -_STEPS]) * abs(mean) / 100.0


def _borderline(mean: float) -> np.ndarray:
    """Thirty differences averaging ``mean`` with a nominal p of about 0.029 at any scale."""
    return _BORDERLINE * (mean / float(_BORDERLINE.mean()))


def _null(mean: float = 1.0) -> np.ndarray:
    """Thirty differences straddling zero symmetrically; Wilcoxon p is about 0.77."""
    return mean + np.concatenate([_STEPS, -_STEPS]) * 2.0


def _untestable() -> np.ndarray:
    """Ten non-zero differences among thirty: too few for a Wilcoxon test."""
    return np.array([1.0] * 10 + [0.0] * 20)


def _family_effects(diffs: Mapping[str, tuple[np.ndarray, np.ndarray]]) -> dict[str, Any]:
    """Evaluate every trait of ``diffs`` (raw, adjusted paired differences) as one family."""
    trait_ids = list(diffs)
    n_pairs = len(next(iter(diffs.values()))[0])
    rows: list[dict[str, object]] = []
    adjusted: list[list[float]] = []
    for index in range(n_pairs):
        for inc_dec in ("inc", "dec"):
            row: dict[str, object] = {
                "study_id": f"s{index}",
                "experiment_id": f"e{index}",
                "inc_dec": inc_dec,
            }
            adjusted_row = []
            for trait_id, (raw, adj) in diffs.items():
                row[trait_id] = 1000.0 + raw[index] if inc_dec == "inc" else 1000.0
                adjusted_row.append(float(adj[index]) if inc_dec == "inc" else 0.0)
            rows.append(row)
            adjusted.append(adjusted_row)
    effects = evaluate_table1_paired_effects(
        pd.DataFrame(rows), np.array(adjusted), trait_ids, {t: t for t in trait_ids}
    )
    return {effect.trait_id: effect for effect in effects}


def _bh(p_values: Sequence[float] | np.ndarray) -> np.ndarray:
    """Independent oracle: SciPy's BH over the given vector."""
    return np.asarray(false_discovery_control(np.asarray(p_values, dtype=float), method="bh"))


def test_q_values_are_bh_over_the_reported_family_and_relabel_nominal_traits() -> None:
    """A raw or adjusted p <= 0.05 whose q exceeds 0.05 no longer counts as significant."""
    effects = _family_effects(
        {
            "T:strong": (_strong(12.0), _strong(3.0)),
            "T:nominal_raw": (_borderline(12.0), _strong(3.0)),
            "T:nominal_adj": (_strong(12.0), _borderline(3.0)),
            "T:null_1": (_null(), _null()),
            "T:null_2": (_null(), _null()),
            "T:null_3": (_null(), _null()),
        }
    )
    rows = list(effects.values())
    raw_p = [row.raw_wilcoxon_p for row in rows]
    adj_p = [row.adj_wilcoxon_p for row in rows]
    np.testing.assert_allclose([row.raw_wilcoxon_bh_q for row in rows], _bh(raw_p))
    np.testing.assert_allclose([row.adj_wilcoxon_bh_q for row in rows], _bh(adj_p))
    assert {row.bh_family_size for row in rows} == {6}

    nominal_raw = effects["T:nominal_raw"]
    assert nominal_raw.raw_wilcoxon_p <= 0.05 < nominal_raw.raw_wilcoxon_bh_q
    assert nominal_raw.adj_wilcoxon_bh_q <= 0.05
    assert (
        adjustment_status(
            nominal_raw.raw_mean_diff,
            nominal_raw.raw_wilcoxon_p,
            nominal_raw.adj_mean_diff,
            nominal_raw.adj_wilcoxon_p,
        )
        == "survived (significant post-adjustment)"
    )
    assert nominal_raw.status == "significant only post-adjustment"

    nominal_adj = effects["T:nominal_adj"]
    assert nominal_adj.adj_wilcoxon_p <= 0.05 < nominal_adj.adj_wilcoxon_bh_q
    assert nominal_adj.raw_wilcoxon_bh_q <= 0.05
    assert nominal_adj.status == "collapsed (attenuated to nonsignificant)"

    assert effects["T:strong"].status == "survived (significant post-adjustment)"


def test_untestable_traits_stay_in_the_family_at_p_one_and_export_nan_q() -> None:
    effects = _family_effects(
        {
            "T:strong": (_strong(12.0), _strong(3.0)),
            "T:nominal": (_borderline(12.0), _borderline(3.0)),
            "T:untestable": (_untestable(), _untestable()),
        }
    )
    untestable = effects["T:untestable"]
    assert math.isnan(untestable.raw_wilcoxon_p)
    assert math.isnan(untestable.raw_wilcoxon_bh_q)
    assert math.isnan(untestable.adj_wilcoxon_bh_q)
    assert untestable.status == "insufficient non-zero pairs"
    assert {effect.bh_family_size for effect in effects.values()} == {3}

    nominal = effects["T:nominal"]
    padded = _bh([effects["T:strong"].raw_wilcoxon_p, nominal.raw_wilcoxon_p, 1.0])
    dropped = _bh([effects["T:strong"].raw_wilcoxon_p, nominal.raw_wilcoxon_p])
    assert nominal.raw_wilcoxon_bh_q == pytest.approx(padded[1])
    assert nominal.raw_wilcoxon_bh_q != pytest.approx(dropped[1])


def test_family_bh_q_values_holds_undefined_p_at_one_and_needs_every_member() -> None:
    family_bh_q_values = density_module.family_bh_q_values
    q = family_bh_q_values([0.01, float("nan"), 0.04], family_size=3)
    np.testing.assert_allclose([q[0], q[2]], _bh([0.01, 1.0, 0.04])[[0, 2]])
    assert math.isnan(q[1])
    assert family_bh_q_values([], family_size=0) == []
    # A declared member without a p (NaN for untestable) is an error, not padding (#1324).
    for p_values, family_size in (([0.01], 6), ([], 6), ([0.01, 0.02], 1)):
        with pytest.raises(ValueError, match="family_size"):
            family_bh_q_values(p_values, family_size=family_size)


def test_family_bh_q_values_accepts_a_numpy_array() -> None:
    """A length check, not truthiness: an array of two or more values is valid (#1325)."""
    p_values = np.array([0.01, np.nan, 0.04])
    q = density_module.family_bh_q_values(p_values, family_size=3)
    np.testing.assert_allclose([q[0], q[2]], _bh([0.01, 1.0, 0.04])[[0, 2]])
    assert math.isnan(q[1])
    assert density_module.family_bh_q_values(np.array([]), family_size=0) == []


@pytest.mark.parametrize(
    ("q", "text"),
    [
        (0.0504, "5.04e-02"),
        (0.0496, "4.96e-02"),
        (0.05, "5.00e-02"),
        # Three significant digits would round these onto 0.05, across the threshold.
        (0.05004, "5.004e-02"),
        (0.050000001, "5.0000001e-02"),
        (0.000345, "3.45e-04"),
        (float("nan"), "NA"),
    ],
)
def test_q_text_never_rounds_across_the_threshold(q: float, text: str) -> None:
    """A q just above 0.05 must not display as 0.05 next to a nonsignificant label (#1325)."""
    assert density_module.format_q_value(q) == text
    assert density_module.format_p_and_q(0.01, q) == f"nominal p=1.00e-02, q={text}"
    if not math.isnan(q):
        assert (float(text) <= 0.05) is (q <= 0.05)


@pytest.mark.parametrize("scalar", [np.float64, np.float32])
@pytest.mark.parametrize(
    ("q", "text"),
    [
        (0.0498, "4.98e-02"),
        (0.0496, "4.96e-02"),
        (0.0504, "5.04e-02"),
        (0.05004, "5.004e-02"),
        (float("nan"), "NA"),
    ],
)
def test_q_text_formats_numpy_scalars_as_numbers(scalar: Any, q: float, text: str) -> None:
    """A q read back from the TSV is a NumPy scalar; it printed as its repr (#1330)."""
    value = scalar(q)
    assert density_module.format_q_value(value) == text
    assert density_module.format_p_and_q(scalar(0.01), value) == f"nominal p=1.00e-02, q={text}"


def test_float32_q_at_the_threshold_displays_on_the_side_of_its_label() -> None:
    """float32(0.05) is 0.0500000007 as a float: the label and the display both read it as
    above 0.05, although NumPy compares it equal to 0.05 in float32 (#1330)."""
    value = np.float32(0.05)
    text = density_module.format_q_value(value)
    assert text == "5.0000001e-02"
    assert float(text) > 0.05
    assert adjustment_status(12.0, 0.001, 3.0, value) == "collapsed (attenuated to nonsignificant)"


def test_q_value_records_accept_numpy_integers() -> None:
    """A family size taken from a NumPy count is still written as a JSON integer (#1330)."""
    record = density_module.q_value_provenance(np.int64(6))
    assert json.loads(json.dumps(record))["family_size"] == 6
    assert type(record["family_size"]) is int
    assert "over all 6 reported traits" in density_module.q_value_statement(np.int64(6))[0]
    q = density_module.family_bh_q_values(
        np.array([0.01, np.nan], dtype=np.float32), family_size=np.int64(2)
    )
    assert all(type(value) is float for value in q)
    assert q[0] == pytest.approx(0.02)
    assert math.isnan(q[1])


def test_default_run_rejects_a_table1_trait_absent_from_the_columns(tmp_path: Path) -> None:
    """An absent Table 1 trait fails the run, as in list mode, rather than joining the
    family without a row (#1324)."""
    profiles, columns = _write_fixture(tmp_path)
    table = pd.read_csv(columns, sep="\t")
    table[table["column_id"] != CONSTANT_TRAIT].to_csv(columns, sep="\t", index=False)
    with pytest.raises(ValueError, match=r"Table 1 traits .*METPO:1000683"):
        run_density_audit(profiles, columns, tmp_path / "output")
    assert not (tmp_path / "output" / "table1_density_adjusted_effects.tsv").exists()


def test_effects_reject_reported_traits_absent_from_the_columns() -> None:
    """The family is the reported rows: a trait that cannot have a row is an error (#1324)."""
    family = {"T:strong": (_strong(12.0), _strong(3.0)), "T:null": (_null(), _null())}
    effects = _family_effects(family)
    assert {effect.bh_family_size for effect in effects.values()} == {len(effects)}
    profiles = pd.DataFrame(
        {"study_id": ["s0", "s0"], "experiment_id": ["e0", "e0"], "inc_dec": ["inc", "dec"]}
        | {trait: [1.0, 0.0] for trait in family}
    )
    with pytest.raises(ValueError, match=r"not non-MediaDive trait columns: \['T:absent'\]"):
        evaluate_table1_paired_effects(
            profiles,
            np.zeros((2, len(family))),
            list(family),
            {**{t: t for t in family}, "T:absent": "absent"},
        )


_STATUS_CASES: list[tuple[Callable[[], np.ndarray], Callable[[], np.ndarray], str, str]] = [
    (
        lambda: np.zeros(30),
        lambda: _borderline(5.0),
        "undefined magnitude change (raw mean zero, significant post-adjustment)",
        "undefined magnitude change (raw mean zero, nonsignificant post-adjustment)",
    ),
    (
        lambda: _strong(4.0),
        lambda: _borderline(-6.0),
        "reversed (significant post-adjustment)",
        "reversed (nonsignificant post-adjustment)",
    ),
    (
        lambda: _borderline(12.0),
        lambda: _strong(3.0),
        "survived (significant post-adjustment)",
        "significant only post-adjustment",
    ),
    (
        lambda: _borderline(-5.2),
        lambda: _strong(-9.7),
        "enhanced (significant, larger magnitude)",
        "significant only post-adjustment",
    ),
    (
        lambda: _strong(12.0),
        lambda: _borderline(3.0),
        "survived (significant post-adjustment)",
        "collapsed (attenuated to nonsignificant)",
    ),
    (
        lambda: _strong(-5.0),
        lambda: _borderline(-9.0),
        "enhanced (significant, larger magnitude)",
        "collapsed (nonsignificant post-adjustment, larger magnitude)",
    ),
    (
        lambda: _borderline(12.0),
        lambda: _borderline(3.0),
        "survived (significant post-adjustment)",
        "attenuated (nonsignificant post-adjustment)",
    ),
    (
        lambda: _borderline(6.95),
        lambda: _borderline(187.8),
        "enhanced (significant, larger magnitude)",
        "enhanced magnitude (nonsignificant post-adjustment)",
    ),
    (
        lambda: _strong(12.0),
        lambda: _strong(3.0),
        "survived (significant post-adjustment)",
        "survived (significant post-adjustment)",
    ),
    (
        lambda: _strong(-5.2),
        lambda: _strong(-9.7),
        "enhanced (significant, larger magnitude)",
        "enhanced (significant, larger magnitude)",
    ),
]


@pytest.mark.parametrize(("raw", "adj", "p_status", "q_status"), _STATUS_CASES)
def test_every_status_branch_follows_q_not_p(
    raw: Callable[[], np.ndarray],
    adj: Callable[[], np.ndarray],
    p_status: str,
    q_status: str,
) -> None:
    """One focal trait among five null traits: a nominal p of 0.029 has q of about 0.18."""
    family = {"T:focal": (raw(), adj())}
    family.update({f"T:null_{i}": (_null(), _null()) for i in range(5)})
    focal = _family_effects(family)["T:focal"]

    nominal = adjustment_status(
        focal.raw_mean_diff, focal.raw_wilcoxon_p, focal.adj_mean_diff, focal.adj_wilcoxon_p
    )
    assert nominal == p_status
    assert focal.status == q_status
    assert focal.bh_family_size == 6
    assert focal.status == adjustment_status(
        focal.raw_mean_diff,
        focal.raw_wilcoxon_bh_q,
        focal.adj_mean_diff,
        focal.adj_wilcoxon_bh_q,
    )


EXTRA_TRAITS = ("METPO:1000605", "METPO:1000607")
CONSTANT_TRAIT = "METPO:1000683"  # sphere shaped: no paired difference, so untestable
MEDIA_TRAIT = "METPO:2000517|mediadive.medium:1"


def _write_fixture(tmp_path: Path) -> tuple[Path, Path]:
    """Paired profiles over the six Table 1 traits (one constant), two more and a medium."""
    rng = np.random.default_rng(1321)
    traits = [*TABLE_1_TRAITS, *EXTRA_TRAITS]
    rows = []
    for experiment in range(40):
        for inc_dec in ("inc", "dec"):
            counts = dict(zip(traits, rng.poisson(rng.uniform(2, 30, len(traits))), strict=True))
            counts[CONSTANT_TRAIT] = 5
            rows.append(
                {
                    "study_id": f"s{experiment // 2}",
                    "experiment_id": f"e{experiment}",
                    "inc_dec": inc_dec,
                    **{trait: int(value) for trait, value in counts.items()},
                    MEDIA_TRAIT: int(rng.poisson(40)),
                }
            )
    profiles = tmp_path / "profiles.tsv"
    pd.DataFrame(rows).to_csv(profiles, sep="\t", index=False)
    names = {**TABLE_1_TRAITS, **{trait: f"extra {trait}" for trait in EXTRA_TRAITS}}
    columns = tmp_path / "columns.tsv"
    pd.DataFrame(
        {
            "column_id": [*traits, MEDIA_TRAIT],
            "column_label": [*(f"{names[t]} {t}" for t in traits), "medium 1"],
            "object_name": [*(names[t] for t in traits), "medium 1"],
            "object_category": ["biolink:PhenotypicQuality"] * len(traits) + ["METPO:1004005"],
        }
    ).to_csv(columns, sep="\t", index=False)
    return profiles, columns


def _assert_table_q_values(table: pd.DataFrame, family_size: int) -> None:
    """Each table's q-values are BH over its own rows, NaN p held at 1, NaN q exported."""
    assert set(table["bh_family_size"]) == {family_size}
    assert len(table) == family_size
    for p_column, q_column in (
        ("raw_wilcoxon_p", "raw_wilcoxon_bh_q"),
        ("adj_wilcoxon_p", "adj_wilcoxon_bh_q"),
    ):
        p = table[p_column].to_numpy(dtype=float)
        expected = _bh(np.where(np.isnan(p), 1.0, p))
        expected[np.isnan(p)] = np.nan
        np.testing.assert_allclose(table[q_column].to_numpy(dtype=float), expected)


def test_tsv_reports_q_values_over_the_six_traits_and_the_declared_list(tmp_path: Path) -> None:
    profiles, columns = _write_fixture(tmp_path)
    declared_ids = [*EXTRA_TRAITS, *TABLE_1_TRAITS]
    list_path = tmp_path / "traits.txt"
    list_path.write_text("\n".join(declared_ids) + "\n", encoding="utf-8")

    run_density_audit(profiles, columns, tmp_path / "default")
    run_density_audit(
        profiles, columns, tmp_path / "declared", trait_list=load_trait_list(list_path)
    )

    default_path = tmp_path / "default/table1_density_adjusted_effects.tsv"
    default = pd.read_csv(default_path, sep="\t")
    declared = pd.read_csv(tmp_path / "declared/table1_density_adjusted_effects.tsv", sep="\t")
    _assert_table_q_values(default, 6)
    _assert_table_q_values(declared, 8)

    header = default_path.read_text(encoding="utf-8").split("\n")[0].split("\t")
    assert header[-3:] == ["raw_wilcoxon_bh_q", "adj_wilcoxon_bh_q", "bh_family_size"]
    # Untestable rows export an empty q, not 1.0.
    constant_line = next(
        line
        for line in default_path.read_text(encoding="utf-8").splitlines()
        if line.startswith(CONSTANT_TRAIT)
    ).split("\t")
    assert constant_line[header.index("raw_wilcoxon_bh_q")] == ""
    assert constant_line[header.index("adj_wilcoxon_bh_q")] == ""

    # The p-values of a shared trait do not depend on the list; its q-values may.
    shared = declared.set_index("trait_id").loc[list(TABLE_1_TRAITS)].reset_index()
    for column in ("raw_mean_diff", "raw_wilcoxon_p", "adj_mean_diff", "adj_wilcoxon_p"):
        np.testing.assert_array_equal(shared[column].to_numpy(), default[column].to_numpy())

    for name, family_size in (("default", 6), ("declared", 8)):
        report = (tmp_path / name / "density_adjustment_report.txt").read_text(encoding="utf-8")
        assert "Benjamini-Hochberg q <= 0.05" in report
        assert f"over all {family_size} reported traits" in report
        assert "q=NA" in report

    provenance = json.loads(
        (tmp_path / "declared" / TRAIT_LIST_PROVENANCE_FILENAME).read_text(encoding="utf-8")
    )
    assert provenance["q_values"] == {
        "procedure": "benjamini_hochberg",
        "family": "every reported trait",
        "family_size": 8,
        "corrected_separately": ["raw_wilcoxon_p", "adj_wilcoxon_p"],
        "untestable_p_held_at": 1.0,
        "label_threshold_q": 0.05,
    }


def test_default_run_records_machine_readable_provenance_matching_trait_list_schema(
    tmp_path: Path,
) -> None:
    """Default (arm) density runs write machine-readable BH provenance reusing trait-list
    fields and names, while leaving trait-list provenance byte-for-byte unchanged (#1327)."""
    import gzip

    import bsdbng.density_adjustment as bsdbng_density
    from mtp.trait_lists import file_sha256
    from mtp.workflow.arms import ilyas, jed
    from mtp.workflow.arms.types import ArmContext

    profiles, columns = _write_fixture(tmp_path)
    declared_ids = [*EXTRA_TRAITS, *TABLE_1_TRAITS]
    list_path = tmp_path / "traits.txt"
    list_path.write_text("\n".join(declared_ids) + "\n", encoding="utf-8")
    trait_list = load_trait_list(list_path)

    default_dir = tmp_path / "default"
    declared_dir = tmp_path / "declared"
    default_res = run_density_audit(profiles, columns, default_dir, code_commit="abc123")
    declared_res = run_density_audit(
        profiles, columns, declared_dir, code_commit="abc123", trait_list=trait_list
    )

    assert (default_dir / "density_adjustment_provenance.json").is_file()
    filename = density_module.DENSITY_PROVENANCE_FILENAME
    assert filename == "density_adjustment_provenance.json"
    assert filename == bsdbng_density.DENSITY_PROVENANCE_FILENAME
    assert bsdbng_density.density_provenance is density_module.density_provenance
    assert bsdbng_density.write_density_provenance is density_module.write_density_provenance

    default_prov_path = default_dir / filename
    declared_prov_path = declared_dir / TRAIT_LIST_PROVENANCE_FILENAME
    assert default_prov_path.is_file()
    assert not (default_dir / TRAIT_LIST_PROVENANCE_FILENAME).exists()
    assert default_res["density_provenance"] == default_prov_path
    assert default_res["trait_list_provenance"] is None
    assert default_res["provenance"] == default_prov_path

    assert declared_prov_path.is_file()
    assert not (declared_dir / filename).exists()
    assert declared_res["density_provenance"] is None
    assert declared_res["trait_list_provenance"] == declared_prov_path
    assert declared_res["provenance"] == declared_prov_path

    default_prov = json.loads(default_prov_path.read_text(encoding="utf-8"))
    declared_prov = json.loads(declared_prov_path.read_text(encoding="utf-8"))

    # Default provenance reuses the exact field names and order from trait-list mode,
    # omitting only the two declared-list file fields.
    assert list(default_prov) == [
        key for key in declared_prov if key not in {"trait_list_file", "trait_list_sha256"}
    ]
    expected_default_q = {
        "procedure": "benjamini_hochberg",
        "family": "every reported trait",
        "family_size": 6,
        "corrected_separately": ["raw_wilcoxon_p", "adj_wilcoxon_p"],
        "untestable_p_held_at": 1.0,
        "label_threshold_q": 0.05,
    }
    assert default_prov == {
        "trait_count": 6,
        "trait_ids": list(TABLE_1_TRAITS),
        "reported_table": "table1_density_adjusted_effects.tsv",
        "pc1_fit_traits": "every non-MediaDive trait column, independent of the list",
        "q_values": expected_default_q,
        "profiles_file": "profiles.tsv",
        "profiles_sha256": file_sha256(profiles),
        "columns_file": "columns.tsv",
        "columns_sha256": file_sha256(columns),
        "code_commit": "abc123",
    }
    assert default_res["q_values"] == expected_default_q
    assert declared_res["q_values"] == declared_prov["q_values"]

    # Trait-list provenance is byte-for-byte unchanged.
    expected_declared_text = (
        json.dumps(
            {
                **trait_list.provenance(),
                "reported_table": "table1_density_adjusted_effects.tsv",
                "pc1_fit_traits": "every non-MediaDive trait column, independent of the list",
                "q_values": density_module.q_value_provenance(len(declared_ids)),
                "profiles_file": "profiles.tsv",
                "profiles_sha256": file_sha256(profiles),
                "columns_file": "columns.tsv",
                "columns_sha256": file_sha256(columns),
                "code_commit": "abc123",
            },
            indent=2,
        )
        + "\n"
    )
    assert declared_prov_path.read_text(encoding="utf-8") == expected_declared_text

    # Gzip fallback hashes the actual .gz files read in default mode too.
    gz_dir = tmp_path / "gz_inputs"
    gz_dir.mkdir()
    profiles_gz = gz_dir / "profiles.tsv.gz"
    columns_gz = gz_dir / "columns.tsv.gz"
    profiles_gz.write_bytes(gzip.compress(profiles.read_bytes()))
    columns_gz.write_bytes(gzip.compress(columns.read_bytes()))
    gz_out = tmp_path / "gz_default"
    run_density_audit(gz_dir / "profiles.tsv", gz_dir / "columns.tsv", gz_out, code_commit="gz123")
    gz_prov = json.loads((gz_out / filename).read_text(encoding="utf-8"))
    assert gz_prov["profiles_file"] == "profiles.tsv.gz"
    assert gz_prov["profiles_sha256"] == file_sha256(profiles_gz)
    assert gz_prov["columns_file"] == "columns.tsv.gz"
    assert gz_prov["columns_sha256"] == file_sha256(columns_gz)

    # Both arm plans (jed and ilyas) declare the default density provenance file.
    products = dict.fromkeys(
        (
            "raw_dump",
            "count_profiles",
            "experiments",
            "feature_metadata",
            "reported_taxa",
            "taxa_matrix",
            "study_table",
            "kg_edges",
            "kg_nodes",
            "source_identity",
            "ora_profiles",
        ),
        tmp_path / "core" / "in.tsv",
    )
    arm_ctx = ArmContext(
        repository=Path(__file__).resolve().parents[1],
        manifest=tmp_path / "core" / "run.json",
        output=tmp_path / "arm",
        inputs=tmp_path / "arm" / "inputs",
        products=products,
        code_commit="abc123",
        replicates=None,
        seed=None,
        figures=True,
    )
    jed_density = next(c for c in jed.build_commands(arm_ctx) if c.id == "density_adjustment")
    ilyas_density = next(
        c for c in ilyas.build_commands(arm_ctx) if c.id == "figure_density_adjustment"
    )
    assert (arm_ctx.output / "density_adjustment" / filename) in jed_density.outputs
    assert (arm_ctx.output / "density_adjustment" / filename) in ilyas_density.outputs
