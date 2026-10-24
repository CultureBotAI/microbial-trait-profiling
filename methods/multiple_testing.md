# Additive full-family multiple testing

Policy version 1 is defined in [multiple_testing_policy_v1.yaml](multiple_testing_policy_v1.yaml)
and validated by [the LinkML schema](../schema/multiple_testing.yaml). The engine
in `mtp.multiple_testing` reports BY as primary **within this additive sidecar**,
and BH for comparison, at inclusive adjusted-p thresholds 0.05 and 0.001. Native
arm calls and arm execution plans remain unchanged. Changing which calls are
primary in the scientific report requires a separate decision.

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
`expected_family_size`. A family may contain one test per trait, fifteen ordered
ratio comparisons, or profile-by-trait-by-tail tests. An adapter declaring a
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
