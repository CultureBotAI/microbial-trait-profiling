# Table 1 taxon contributions and lineage influence

Issue [#190](https://github.com/CultureBotAI/microbial-trait-profiling/issues/190).
The [measured run summary](table1_taxon_influence_results.md) records the reconciled
corpus, observed contribution patterns, conditional uncertainty and remaining work.
The [native rerun extension](table1_native_influence.md) now adds composition and
both ORA background treatments. Statements below about missing native results
describe the original contribution/shared-effect stage; consult the completed
`native_manifest.json` and native summary for the extension and its review limits.
Run `just table1-taxon-influence`, or pass a directory containing the **same pinned**
KG files: `just table1-taxon-influence C:/path/to/pinned-kg`.
The recipe writes TSVs, compressed TSVs, `manifest.json` and
`influence_summary.md` to `data/derived/table1_taxon_influence/`. It neither
downloads a replacement KG nor changes input profiles. `just check` validates
the implementation; the full recipe additionally validates the local corpus.

This is a provisional annotation influence diagnostic. The archived Rosalind
discussion and separate review in [PR #192](https://github.com/CultureBotAI/microbial-trait-profiling/pull/192)
motivate tests, not established contributors or mechanisms. Contrast eligibility
and orientation remain subject to #188; feature semantics remain subject to #189.
Increased/decreased here are curated signature labels, not verified disease/control
labels. Neither taxon removal nor robustness establishes causation, oxygen exposure,
pathway activity, SCFA production, or independent confirmation. Apply the shared
interpretation guidance in [#179 / PR #193](https://github.com/CultureBotAI/microbial-trait-profiling/pull/193).

## Evidence and count unit

The chain is source study → curated signature → distinct reported NCBI taxon →
KG source taxon and trait key → taxonomy propagation → profile count → fraction
or ratio → experiment contrast → study average → across-study summary.

The implementation mirrors `taxa_trait_matrix.make_trait_key` and its aggregation:
one **distinct source-taxon/trait-key combination** contributes one count to each
reported ancestor, including itself. Multiple KG edges supporting that combination
do not each add a count. A source can nevertheless contribute repeatedly when both
its parent and child occur among reported taxa. Profile aggregation deduplicates
reported taxon IDs within an experiment/direction, including across signatures.
These are not independent organisms or patient-level abundance measurements.

The audit reconstructs **every non-MediaDive column**, first for every taxon in the
matrix and then for every retained human paired profile. Any mismatch stops all
influence calculations and is recorded in `reconciliation.tsv`. Empty mismatch
output means exact equality, including zeros. It does not validate the upstream
source assertions or source count contract (#184).

The protocol pins the 20260120 KG by file MD5; the run manifest additionally records
SHA-256 for KG files, profile/matrix/taxon/experiment inputs, protocol and implementation.
No current taxonomy snapshot is substituted. The pinned edge export has a split header
and unquoted line breaks following `manual_validation_of_automated_agent`. The reader
joins only those recognized continuations after checking their field structure and
records their number. Other malformed records are errors. Original edge IDs and
starting physical line numbers remain available. This file-format repair does not
change the source-count algorithm.

KG edge IDs are not guaranteed to identify original source assertions. Original
assertion counts and source release versions remain unresolved when unavailable.
`primary_source` is the original KG provenance cell, kept without translating it
into an inferred gene, phenotype measurement, or independent evidence source.

## Relational contribution export

The following joins are intentional: expanding every signature × source × supporting
edge into a flat table would incorrectly invite counting supporting edges repeatedly.
All tables in a run share the KG and input revisions in `manifest.json`.
The source-contribution and supporting-edge reference tables include the full
matrix/KG supersets. Restrict them through `profile_reported_taxa.tsv.gz` when
interpreting the retained human corpus; summing the reference tables directly
does not recover a human profile count.

| File | Key and interpretation |
| --- | --- |
| `input_eligibility.tsv` | Every original profile; retained, unpaired, nonhuman or missing metadata. Orientation explicitly unreviewed. |
| `profile_reported_taxa.tsv.gz` | Study, experiment, direction, exact condition/site, reported taxon/name/rank and JSON list of source signature IDs. |
| `source_contributions.tsv.gz` | Join on `reported_taxon`; one row per reported taxon, contributing source taxon and six-trait key. Count contribution is one; includes deterministic rollup path and alternative-parent-path flag. |
| `source_assertions.tsv.gz` | Join on `source_taxon, trait_key`; retains every supporting KG edge, source provenance cell and physical line. This is an evidence join, **not an additional count multiplier**. |
| `profile_counts.tsv.gz` | All six traits, both directions, including zeros; existing count, distinct reported/source taxa, descendant sources, KG edges, repeated propagation and full denominator. Original-assertion count is missing when unresolved. |
| `contributing_universe.tsv` | Full positive contribution universe for each declared scope, trait, direction and resolution level; no display cutoff. |
| `concentration.tsv` | All six traits and both directions; resolved/unresolved mass, resolved taxon count, HHI and largest share conditional on resolved mass. |
| `cross_trait_overlap.tsv` | All 15 trait pairs at species/genus/family, separately by curated direction and scope. |
| `scenarios.tsv`, `predictions.tsv` | Full lineage inventory, support, selection/status and the four predicted genera, whether or not resolved. |
| `effects.tsv.gz` | Scenario/scope/representation/feature effects, uncertainty, directional support, eligibility changes and separate correction families. |

`distinct_source_taxa` includes directly reported source taxa; `descendant_source_taxa`
counts contributing taxa not themselves reported in that profile. A reported source
can also propagate into a reported ancestor. Inspect the joined rollup paths for this
overlap; repeated propagation is the total count minus distinct source-taxon/trait
units in that profile. A source's directly asserted rank is not inferred from its name.

Species/genus/family resolution uses explicit ranks in the reported-taxon/matrix
inputs together with pinned KG ancestry. A missing rank or multiple candidate ancestors
stays unresolved/ambiguous. The largest-share and HHI columns describe the **resolved
subset** and must travel with coverage counts. Jaccard uses presence sets; weighted
overlap is the sum of minimum normalized contribution shares over shared taxa.
Neither statistic tests independence. Source concentration uses provenance-cell
incidence: multiple provenance cells can support one deduplicated count unit, so its
denominator differs from the taxonomic count denominator.

## Rule-declared exclusions and scope

[`table1_taxon_influence_protocol.json`](table1_taxon_influence_protocol.json) is an
input, hashed before effects are calculated. It declares pooled human paired profiles,
exact `Feces`, exact `Saliva`, and exact `Parkinson's disease`/`Feces`; it does not merge
site or condition labels. All source records retain their original exact condition/site
even outside those three strata. An empty stratum produces untestable results.

The full genus/family inventory is retained. The exclusion panel comprises the ten
largest genera and ten largest families by **full non-MediaDive propagated contribution
count**, plus Streptococcus, Enterococcus, Staphylococcus, Veillonella and their resolved
families. The rule-declared extension adds the leading genus for each of the six Table 1
trait keys by propagated contribution count summed across both directions. A lineage
needs contributions in at least three retained studies. Ties are broken by taxon ID;
signs, effects and p-values do not select scenarios. The extension was recorded on
[#190](https://github.com/CultureBotAI/microbial-trait-profiling/issues/190#issuecomment-5703892447)
before the extended run. This bounded panel tests the named predictions alongside large
overall and per-trait contributors. It does not claim exhaustive influence testing of
every lineage or freedom from data selection. Unselected and unsupported lineages
remain explicitly recorded.

For each lineage, removal deletes all known source-taxon/trait units whose source
falls within that lineage, wherever they propagate, including into higher-rank
reported records. Both arms are reconstructed. Every one of the 317 retained columns
contributes to the new denominator; removing only the six selected columns would
answer a different and incorrect question. Direct annotations to ancestors of the
excluded lineage cannot be allocated to a descendant confidently: retain them and
flag the ambiguous source/target combinations. The result is a diagnostic conditional
on those retained ambiguous counts, not complete biological removal of a lineage.

Each scenario directory contains:

| File | Purpose |
| --- | --- |
| `removed_source_counts.tsv.gz` | Sparse removal vector for **all** non-MediaDive traits by reported taxon; join through the profile-reported-taxon table to reconstruct full remaining profiles. |
| `ambiguous_retained_sources.tsv` | Explicit source/target records retained because higher-rank direct assertions cannot be apportioned. |
| `profile_counts.tsv.gz` | Original/remaining six-trait counts and full totals, both arms, with ambiguity flag. |
| `pair_eligibility.tsv.gz` | Original/remaining totals and zero-valued six-trait keys for each arm; establishes the precise membership and exclusion reason for every fraction and all 15 ratios. |
| `ora_background_counts.tsv` | Full count backgrounds separately labeled `fixed_baseline` and `recomputed_remaining`, per scope, using both arms of all retained profiles including zero totals. |

The ORA exports are count contracts for Ilyas, not native ORA reruns or odds ratios.
A fixed baseline asks about residual annotations against the original corpus; a
recomputed background asks about residual annotations against the residual corpus.
Their native tests, support rules, sampling units and correction families must be
declared by the ORA implementation. Do not mix these backgrounds in a shared ranking.

## Effects, uncertainty and integration

Shared fraction comparisons require positive full totals in both arms. Individual
trait zeros are valid. The contrast is `inc count/inc total − dec count/dec total`,
in fraction units (multiply by 100 for percentage points). Native Noel comparisons
retain his zero-total-as-zero convention and six-trait directional test; an undefined
share/interval stays blank with a status rather than being displayed as zero.

All 15 ratios follow the manifest's six-trait order. For each ordered pair A/B,
the contrast is `log2(A_inc/B_inc) − log2(A_dec/B_dec)`. All four counts must be
positive; no pseudocount is introduced. In each representation, average eligible
experiments within study, then weight eligible studies equally.

The output distinguishes `before`, `after`, `before_matched`, `after_matched` and
`matched_change`. The matched versions re-average the **same surviving experiments**
in each study; matching study IDs alone is insufficient. `matched_change` averages
experiment-level differences on this common set. Own-population change decomposes as:

`after − before = before_eligibility_component + matched_change + after_eligibility_component`.

Components are blank when not identifiable, including no common eligible pairs.
New zero totals and trait zeros remain visible in the pair manifest. Counts removed
are not additive contributions to a nonlinear fraction effect.

Every population provides pair/study counts, positive/negative/tied study counts,
mean effect with conditional 95% Student t interval across study means, and positive
share with exact 95% binomial interval excluding ties (`1e-12` tolerance). The mean
interval concerns effect magnitude; the binomial interval concerns directional share.
They are not interchangeable. The two-sided exact sign-test p-value tests directional
balance, not the mean magnitude. BH6 is separate for each six-trait family and BY15
for each ratio family, separately by scenario, scope, representation and population.
All declared family members count toward multiplicity; unavailable values remain
missing. There is no simultaneous guarantee across scenarios or scopes.

These intervals are conditional on the observed data, selected features and protocol;
they are not full-pipeline calibrated guarantees (#182). Constant/tied, opposite,
null, unsupported and untestable outputs remain visible. Loss of significance alone
does not establish that a lineage explains a signal.

For #178, join on the **full qualified trait key**, scenario, exact scope, representation
and population; retain the effect unit, all uncertainty/support fields, family identity,
input hashes, missingness and ambiguity flags. Native Jed transformed statistics
require a pinned transformation/reference-population policy and full native correction
family; they remain explicitly **not rerun** here. Native ORA remains with Ilyas.
Their absence must remain a missing native result, not be replaced by the shared
sign-test p-value. These owner integrations, #188/#189 prerequisites and review of
resolved source semantics remain required before closing #190 or assigning mechanisms.
