# Comparing off and phenotype evidence (#1041)

This additive consumer compares completed `three-test-per-trait-v1` bundles
on a common archived population. Phenotype is primary and off is a sensitivity
analysis. This primary-scope choice was made after both off and phenotype results were seen;
it is not the off-primary choice in the historical selection-scan protocol.
The comparison preserves that decision in its manifest. It does not rewrite
historical protocols or their execution records.

## Comparable inputs

Supply both parent manifests, explicit SHA-256 pins, and their comparison
manifests. The shared completed-table reader validates each bundle, its clean
historical producer, artifacts, upstream core and arms, trait catalog, and
persisted experiment and study contrasts. This consumer additionally requires
the same parent producer commit, complete implementation inventory, methods,
dependency lock and numerical environment. Corresponding upstream code and
parameters must agree except for the declared rollup scope. Populations,
direction, exact experiment and study identities, source input identities,
selection audit and experiment metadata must agree. Scope-dependent profile,
catalog and matrix products may differ. Source-identity timestamps and paths
are not biological differences; the corresponding semantic input identities
must match. Analysis-basis differences are limited to catalog, feature count,
feature filtering and rollup policy.

The two scopes need not have the same retained traits. Join on qualified trait
IDs, never on display labels. Export the intersection and the complete union,
including off-only and phenotype-only traits. Membership and counts are derived
from the supplied catalogs. For the development F units, `METPO:1000608`
(facultatively aerobic) is shared; `METPO:1000605` (facultatively anaerobic) is
phenotype-only. Labels remain separately recorded for each scope.

## Preserved quantities and families

The shared-trait table places every parent evidence column side by side with
`off_` and `phenotype_` prefixes. Original parent evidence files are also copied
byte-for-byte. Preserve native Noel study votes and nominal BY q values, the
full-family study-mean t-test and its BY q values, support, explicit statuses,
and all six signed/magnitude ranks. Each test retains its own full family
within its original scope. Never correct again on the shared intersection or
across different tests, and never rank the intersection as though it were the
original candidate family.

`study_mean_delta_z` gives every study equal weight. Its
`study_mean_delta_z_ci95_low` and `study_mean_delta_z_ci95_high` columns are the
parent's Student-t 95% interval across study means. The separately named
`comparison_delta_z_mean` is an explicit alias of the parent's
`experiment_mean_delta_z`, giving every comparison equal weight. No interval
is attached to this comparison mean. A study with several comparisons can
make these two means differ. Untestable values remain missing with their
reasons; they are not replaced by zero or a nonsignificant q value.

Delta-z describes standardized annotation shares among taxa reported increased
versus decreased. Even a shared trait with unchanged direct annotations can
change when the retained catalog changes: the denominator sums all retained
annotation columns, and each column is standardized over its scope's profiles.
A scope difference is therefore not necessarily new biological evidence or
merely a change in multiple-testing correction.

## Declared selection views

There is no cross-test combined ranking. Export ten separate views:

- Native Noel nominal BY q <= 0.05, using its original eligible study votes.
- Full-family study-mean t-test BY q <= 0.05.
- Top 5 and top 10 positive comparison means, and top 5 and top 10 negative
  comparison means, each for all traits and for traits with >=50 active studies.

An active study has `abs(study_mean_delta_z) > 1e-12`. Apply the activity and
positive/negative sign filters before ranking the corresponding view within
each scope. Positive values rank descending; negative values rank ascending.
Exact ties share the minimum rank, so a top-N view may include more than N
traits. Zero and opposite-sign values never enter a directional view. These
`view_rank` values are separate from the unchanged parent full-family ranks.
No significance requirement is added to an effect-ranked view.

For every union trait and view, record availability, eligibility, selection,
and entered/left/retained/neither membership moving from off to phenotype.
An unavailable trait is distinct from an available but unselected or untestable
trait. Selection counts come from the qualified inputs, not historical prose.
Do not compare p-values across tests. Agreement using shared data is consistency,
not independent replication. ORA remains a labelled directional evidence column
and never selects a headline or enters these selection views.

## Execution and provenance

```bash
just rollup-comparison --off-three-test-run OFF_MANIFEST --off-sha256 OFF_SHA256 \
  --off-comparison-run OFF_COMPARISON \
  --phenotype-three-test-run PHENOTYPE_MANIFEST --phenotype-sha256 PHENOTYPE_SHA256 \
  --phenotype-comparison-run PHENOTYPE_COMPARISON --output FRESH_DURABLE_DIRECTORY
```

Products are `shared_traits.tsv`, `trait_catalog.tsv`, `rank_transitions.tsv`,
exact parent copies `off_per_trait_evidence.tsv` and
`phenotype_per_trait_evidence.tsv`, and `rollup_comparison_run.json`.
The manifest records both parent pins and identities, upstream run IDs and
commits, scopes and families, input and output hashes, selection policies,
and the comparison producer's Git and numerical environment identity.

Both source trees and their aliases are protected before publication probes.
Inputs and producer are checked again before an exclusive atomic directory
publication; an existing destination is never overwritten. Development outputs
use durable storage outside version control. Tracked results wait for the
selected new main basis P. No core/arm run or Delta-z transformation is executed
by this consumer.

`just test-rollup-reproduction OFF_PARENT PHENOTYPE_PARENT RESULT_DIRECTORY
CONSUMER_COMMIT SOURCE_COMMIT OFF_SHA256 PHENOTYPE_SHA256` checks completed
development outputs read-only. It requires both full commit pins and independent
SHA-256 pins for both parent manifests, verifies copied parent bytes and every
imported shared-trait field, and independently derives each selection view from
the full parent tables. Its study-t selection counts are derived rather than
fixed to a historical claim. It does not replay either parent or any arm.
