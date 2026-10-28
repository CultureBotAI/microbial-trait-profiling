# Analyses of completed Delta-z tables

Issue #1052 covers rows 1, 2, 3 and 5 of #1045. The input is a completed
`three-test-per-trait-v1` bundle, its explicitly supplied manifest SHA-256,
and the verified upstream comparison. This consumer reads the persisted
per-comparison and per-study contrasts. It does not reconstruct fractions,
standardize profiles, recompute Delta-z, or replay the parent analysis.

The parent definition and the retained-feature denominator are described in
[`three_test_table.md`](three_test_table.md). Delta-z compares trait annotation
shares in taxa reported increased and decreased. It is a standardized
annotation contrast, not a change in measured physiological activity.
Phenotype is primary; off is sensitivity. Full remains exploratory and is not
supported by this consumer or its shared source verifier.

## Inputs and verification

The reader accepts an older clean parent producer. It checks the caller's
manifest hash, the exact bundle inventory, file sizes and hashes, the recorded
producer's Git inventory, and the parent contracts. It verifies the comparison,
core and arms and checks their identities, source commits, retained catalog,
scope, analysis bases, population and direction against the parent.

Pair keys are `(study_id, experiment_id)`; study keys are `study_id`. Keys and
traits must be unique, complete and nonblank, and contrast values must be
finite. The supplied study ledger must agree with equal-weight means of its
supplied pair ledger. These consistency checks do not reproduce the original
profile transformation. TSV floats are parsed with round-trip precision;
mean consistency uses relative tolerance `1e-12` and absolute tolerance `1e-14`.
Noel votes are imported from the parent's native
screen, with its different denominator and exclusion policy preserved.

The SHA-256 pin identifies the completed evidence being consumed. Inventory
and contract checks are distinct from a replay of the original computation;
the existing parent replay command retains its stricter producer requirement.

## Row 1: study activity and two labelled votes

For each retained trait, `sdz_up` counts supplied study-mean contrasts greater
than `1e-12`; `sdz_down` counts values less than `-1e-12`. All remaining study
means are tied for this summary. Active studies equal up plus down. The support
reference counts traits with at least 50 active studies, including exactly 50.

These counts appear beside `noel_studies_up` and `noel_studies_down`. The former
sign the mean transformed contrast; Noel signs the mean raw-fraction
difference under its native eligibility rules. Multi-comparison studies can
therefore receive different votes. Neither vote replaces the other. Existing
experiment-weighted and study-weighted means and ranks remain labelled.

## Row 2: oxygen-label quadrants

The x trait is `METPO:1000608`, facultatively aerobic. The y trait is
`METPO:1000607`, obligately anaerobic. `METPO:1000603` is anaerobic and is never
substituted for the y trait.

Quadrants use exact positive and negative signs, preserving the historical
calculation. A pair is on an axis if either contrast equals zero. This differs
deliberately from row 1's activity threshold. The two discordant quadrants
are x-positive/y-negative and x-negative/y-positive. Their counts enter a
two-sided exact binomial sign test with null probability 0.5. Concordant and
axis pairs do not enter that test. Pearson correlation uses all available
pairs, including axis pairs, without rounding or clipping.

Report these quantities for all supplied comparisons and for stool. Stool is
the exact condition `experiments.body_site == "Feces"`; it is not inferred
from `body_site_category`, a case-insensitive alias, or a condition name.
The output includes both distinct study and comparison counts. A sign test
with no discordant pairs and a correlation with a constant or insufficient
vector have explicit untestable status and missing inferential values.

## Row 3: declared traits in the most frequent stool conditions

The trait list is an explicit hashed input. The standard phenotype list is
[`trait_lists/headline_ten_phenotype.txt`](trait_lists/headline_ten_phenotype.txt).
Preserve its order; do not select traits using these results. Missing declared
traits cause an error. Despite its historical filename, the current list contains
nine traits: its facultative entry is the `METPO:1000608` union under
[`facultative_union_v1`](trait_key_policy.md). The historical F list had ten entries,
including a separate `METPO:1000605` entry absent under off. Historical reproduction
must use that producer's hashed list and raw-object policy, not the current list.

Among the supplied stool comparisons, rank exact `condition` strings by
comparison count descending, then condition name ascending for ties. Select
the first seven by default. There is no special disease exclusion or
hard-coded minimum count. The selection and its counts are recorded.

For each selected condition and declared trait, report the arithmetic mean
of per-comparison Delta-z, distinct study count, comparison count, and a
two-sided one-sample Wilcoxon signed-rank test. Use `zero_method="wilcox"`,
`method="auto"`, and `correction=False`, with the numerical environment
recorded. Exact zeros are excluded from its signed ranks. An all-zero vector
has an explicit untestable status and missing inference. These are nominal
p-values; no corrected family is claimed.

## Row 5: body site by condition

Export every exact `(body_site, condition, trait)` cell represented in the
supplied comparisons. The effect is the arithmetic mean of comparison
Delta-z, alongside comparison and distinct study counts. A study can
contribute to several cells. No additional body-site pooling, disease
normalization, minimum cell size, or inferential test is introduced. Blank
metadata is retained explicitly rather than silently dropped by grouping.

## Dependence and reporting

The row 2 and row 3 tests use comparisons, which can share studies. Their
nominal p-values describe the specified historical analyses; they do not
establish independence of comparisons, provide study-cluster calibration, or
determine a headline call. Report effect estimates, denominators and the
named representation. Do not compare p-values across different tests.

## Execution and products

Run `just delta-z-derived --three-test-run <manifest> --three-test-sha256
<sha256> --comparison-run <comparison> --trait-list <declared-list> --output
<fresh-durable-directory>`. Use `--summary-only` in place of `--trait-list`
for row 1 alone, including the off support reference. `--top-diseases` changes
the declared maximum number of stool conditions, not their ranking rule.

The output includes exact byte copies of the parent pair and study ledgers,
the selected derived tables, and `delta_z_derived_run.json`. Provenance records
parent and upstream identities, hashes, scope, population, direction, policies,
trait-list identity when applicable, producer identity, numerical environment,
and output hashes. Inputs and producer identity are checked again before
publication. Outputs must be durable and disjoint from source bundles.
Exclusive atomic publication refuses an existing destination and cleans only
staging created by this run.

Development qualification uses the verified F units read-only. Generated
results remain outside version control until the selected main basis P.
`just test-delta-z-reproduction` reads completed off and phenotype consumer
bundles, the historical numerical reference JSON and its explicit SHA-256,
the expected consumer commit, and the independently bound core experiments
TSV. It checks the ordered trait list against that consumer's committed list
and derives expected site/condition groups from the core metadata, rather
than from a second output of the consumer under test. It never replays Delta-z.
