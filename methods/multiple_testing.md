# Additive full-family multiple testing

This is the contract for the additive sidecar. The project-wide rules, and the registry of
every correction family, are in [multiple_testing_policy.md](multiple_testing_policy.md).

Policy version 1 is defined in [multiple_testing_policy_v1.yaml](multiple_testing_policy_v1.yaml)
and validated by [the LinkML schema](../schema/multiple_testing.yaml). The engine
in `mtp.multiple_testing` reports BY as primary **within this additive sidecar**,
and BH for comparison, at inclusive adjusted-p thresholds 0.05 and 0.001. Native
arm calls and arm execution plans remain unchanged. Changing which calls are
primary in the scientific report requires a separate decision.

These thresholds apply separately to each executed adapter family, not to a pooled
set of all sidecar rows. The registered families and source-output routes are listed
in [the adapter notes](../docs/multiple_testing_sidecars.md). Quote the family ID,
its executed `expected_family_size`, the correction procedure and the threshold
beside any reported q; source p-values remain nominal. Recomputed sidecar q-values
must remain distinguishable from native q-values, especially for the joint-tail
ORA family, whose membership differs from the separate native tail families.

BY controls false discovery rate under arbitrary dependence when the individual
null p-values are valid. BH requires independence or suitable positive dependence.
Neither method repairs invalid tests, study pseudoreplication, data-driven lead
selection, or multiple unreported family choices. Correction is within each
declared family, not a claim of global control across every arm or family. The
implementation delegates arithmetic to
[SciPy `false_discovery_control`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html),
which documents these assumptions and the original BH and BY papers.

An adapter must verify upstream execution, derive the complete declared
hypothesis inventory from its AnalysisBasis and selection rule, and supply that
inventory independently of the observed p-value rows. `retained_feature_count`
comes from that basis; it is not a fixed corpus constant and is distinct from
`expected_family_size`. A family may contain one test per trait, fifteen unordered
focal pairs with a fixed ratio orientation, or profile-by-trait-by-tail tests. An adapter declaring a
joint upper/lower-tail family must independently enumerate both tails for every
eligible coordinate and reject a missing counterpart. The generic engine checks
the supplied inventory exactly; it permits one-sided families and cannot infer
that a caller intended a joint-tail inventory. These scientific choices belong
in the adapter, not an inference from whichever rows survive filtering.

Each member declares qualified features, their order (numerator then denominator
for ratios), test, tail and optional eligibility-ledger identity. Schema v1.1
also requires explicit coordinates: profile, stratum, condition set and split.
Empty coordinates denote a pooled test. Condition identifiers are exact source
identifiers or labels, not necessarily ontology mappings; their sorted order is
canonical. Adapters validate required coordinates and their cardinality. The
same ordered features, test, tail and coordinates cannot be declared twice under
different row IDs or eligibility hashes. Families record profile and inference
units, population and direction basis, estimand, selection rule and history,
upstream byte identities and membership hashes. `profile_unit` and
`feature_filter_id` are copied verbatim from the verified AnalysisBasis;
`strain_rollup` explicitly distinguishes `off` and `phenotype`. Full strain
rollup is unsupported.
These represent annotation-based analyses; they do not measure physiological
activity or establish independent cross-arm replication.

Every declared member must have exactly one raw record. An untestable member
must explicitly state its reason and have undefined p. It reserves one member
using conservative p=1 padding internally; its exported p, adjusted p and
decisions remain undefined.
Omitted, duplicated and unexpected members fail validation. Tested members must
have finite p in [0,1] with no undefined reason; negative zero is rejected.
Schema and engine version 1.2 add a per-record `reason_vocabulary`. Its default,
`legacy_free_text_v1`, preserves the foundation API's nonblank free-text reasons,
including existing callers that omit the field. Explicit null and unknown
vocabularies are rejected. `analysis_reason_v1` requires exact membership in the
authoritative LinkML `AnalysisReasonV1` enum; spelling, whitespace and case are
not normalized. Tested rows require no reason under either vocabulary.
The shared semantic validator applies this rule in both correction paths and
when reading corrected JSONL records. Adjusted records inherit the discriminator;
raw/output hashes and replay bind it. Per-row versions may coexist in a family
without changing its declared membership or correction arithmetic.

This preserves caller compatibility, not the hashes of historical records.
Newly serialized legacy rows explicitly include their vocabulary, changing
canonical raw/output identities. Persisted version-1.1 records must be replayed
with their original implementation and environment; they are not silently
upgraded. The statistical policy remains version 1 and its immutable bytes are
unchanged. The post-arm adapters explicitly emit `analysis_reason_v1`; other
foundation callers retain the legacy contract until they deliberately opt in.

An empty family is permitted only
with explicit size zero and empty raw membership. Significance, effect-size and
presentation filters never reduce the correction family.

Canonical hashes ignore hypothesis/source row order while preserving ordered
ratio features and chronological selection history. The execution captures
policy, scientific family, complete raw and adjusted results, engine/model code,
NumPy and SciPy versions, UTC time and output hashes. The policy values are read
from the versioned YAML, whose immutable bytes are pinned by the engine;
thresholds and correction methods come from that validated policy.
`verify_correction` recomputes every
adjusted value and decision against an externally supplied trusted family and
raw results, rejects policy overrides and checks the complete execution record.
It does not authenticate a caller's invented upstream identities: verified
source adapters are required before using this foundation on scientific data.
Archival replay requires the recorded code and numerical environment.

`CorrectionRun` is the complete in-memory contract. Large adapters should write
the inventory, raw and adjusted tables once as separate hashed artifacts, and
keep their provenance summary compact. The foundation does not run any arm or
rewrite any existing result. Run its fixtures with `just test-multiple-testing`;
`just workflow-check` checks the schema and policy; `just check-generated` checks
the generated model against its authoritative schema.

## Pre-testing selection rules and BH family padding in producer scripts

Outside the additive sidecar engine, producer modules that compute native
Benjamini-Hochberg (BH) q-values declare their pre-testing family membership
and untestable-member padding explicitly:

- **Native taxon influence (`src/mtp/native_taxon_influence.py`)**:
  New runs use protocol `190-native-v3`, whose machine-carried composition
  declaration names the selected-trait correction family and p=1 padding.
  Already-written `190-native-v2` manifests retain their finite-test family
  meaning and remain readable; the current producer refuses to execute v2.
  The publication handoff reports `native_family_size` separately from
  `native_tested_count`, using the recorded protocol to choose the denominator.
  Historical v2 tables use `family_tests` even if they also have
  `family_selected_traits`; legacy manifests without a version use the finite-test
  count with an explicit legacy basis. V3 requires the selected-trait count.
  - **Composition t-test family (`transformed`)**: Pre-testing eligibility
    (`row["selected"]`) requires `profile_support >= minimum_profile_support`,
    `profile_support <= max(2, int(len(keys) * maximum_profile_support_fraction))`,
    and positive composition variance (`composition_variance > 0`). Every
    selected trait (`selected_rows`, counted in `family_selected_traits`) remains
    in the BH family passed to `adjust(...)`: selected traits that are untestable
    in paired t-testing (`row["selected"] and row["p"] == ""`, e.g. fewer than
    `minimum_pairs` or contrast SD $\le 10^{-8}$) are padded at `p = 1.0` inside
    the correction vector and receive `q = 1.0` (`status = "untestable"`), while
    `family_tests` records the count of selected traits with finite t-test
    p-values and unselected traits (`status = "not_selected"`) stay outside the
    family with `p = ""` and `q = ""` (#1334).
  - **ORA hypergeometric tail families (`ora_cells`)**: Separate global BH
    families are evaluated per tail (`enrichment` and `depletion`) across all
    active profiles (`row_total > 0`) and retained traits (`base.columns`) before
    filtering to the six Table 1 display traits. The `enrichment` family
    (`enrichment_tests`) tests all positive profile-by-trait cells
    (`profiles[key][trait] > 0`); the `depletion` family (`depletion_tests`)
    tests all cells in active profiles, including zero-count cells.
  - **Shared directional sign-test families (`correct_shared_comparisons`)**:
    Grouped by `("scenario", "scope", "reference")` (composition comparisons) or
    `("scenario", "scope", "background")` (ORA fold comparisons) and adjusted
    separately per population (`before`, `after`, `before_matched`,
    `after_matched`, `matched_change`) over the full declared panel
    (`BH{len(family)}_directional_sign_test`), holding untestable members (`None`)
    in the `len(family)` denominator while exporting `""` for their `q`.
- **PD significant condition-taxon-trait network (`scripts/create_pd_significant_condition_taxon_trait_network.py`)**:
  - **Condition-taxon McNemar/binomial family (`build_condition_taxon_tests`)**:
    Pre-testing support eligibility (`passes_min_taxon_studies`, requiring
    `reported_study_count >= MIN_TAXON_STUDIES`) defines the declared BH family
    of size `bh_family_size = int(support.sum())` (distinct from `n_tests`, the
    count of rows with finite p-values). Supported rows with zero discordant
    studies (`p_value` is NaN) are held in the BH family via conservative
    `fillna(1.0)` padding (`benjamini_hochberg(tests.loc[support, "p_value"].fillna(1.0))`)
    while their exported `q_value` remains NaN (#1321, #1334).
  - **Condition-trait Wilcoxon family (`build_condition_trait_tests`)**:
    Pre-testing support eligibility (`passes_min_paired_studies`, requiring
    `paired_study_count >= MIN_PAIRED_STUDIES`) similarly defines the declared BH
    family of size `bh_family_size = int(support.sum())` and pads untestable
    supported rows at `1.0` via `fillna(1.0)` while keeping their exported
    `q_value` as NaN.
