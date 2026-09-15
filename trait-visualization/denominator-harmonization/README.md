# Issue #174: denominator-only comparison

## Findings

Removing MediaDive from the denominator changes numerical effects and some
thresholded sets, but **does not change the six pooled Table 1 traits or
their directions at q <= 0.001**. The baseline run is reproduced before
applying the correction. No eligible/tested set changes in the five analysis
strata. Thirteen trait-by-stratum effects change sign; every such effect is
nonsignificant at q <= 0.05 both before and after.

### Input scope

| Scope | Profiles | With MediaDive | Zero retained signal | MediaDive only |
| --- | ---: | ---: | ---: | ---: |
| Full matrix | 5,829 | 5,590 | 206 | 11 |
| Paired subset | 4,670 | 4,495 | 153 | 10 |

Both scopes have 2,290 MediaDive columns and 317 retained columns. The paired
subset contains 2,335 pairs from 959 studies. Zero-retained-signal profiles
remain in the comparison with zero fractions, including MediaDive-only rows.
See [input_scope.tsv](input_scope.tsv) for aggregate count totals.

### Composition significance

Each row is its own BH family; family sizes are unchanged. Counts use <=.

| Stratum | Pairs | Studies | Tested traits | q <= .001 before / after | q <= .05 before / after |
| --- | ---: | ---: | ---: | ---: | ---: |
| Pooled | 2,335 | 959 | 315 | 6 / 6 | 48 / 50 |
| Feces | 1,485 | 594 | 314 | 28 / 27 | 93 / 90 |
| Feces excluded | 850 | 365 | 298 | 0 / 0 | 0 / 0 |
| Other | 521 | 223 | 299 | 0 / 0 | 0 / 0 |
| Saliva | 129 | 55 | 285 | 0 / 0 | 0 / 0 |

At q <= .001, fecal mannitol (`METPO:2000006|CHEBI:29864`) leaves the set
(q .0009846 to .0021453). At q <= .05, the pooled set adds pH delta very low
(`METPO:1000473`) and pleomorphic shaped (`METPO:1000679`). The fecal set
loses NaCl optimum mid1 (`METPO:1000466`), temperature range mid1
(`METPO:1000450`), fructose (`CHEBI:28757`), and ribose (`CHEBI:33942`),
and gains pH delta low (`METPO:1000474`). No direction reversal occurs among
traits significant at q <= .05 in either version.

### Six Table 1 traits: pooled

Effects are mean increased-minus-decreased standardized log-fractions, not
raw fraction differences. All six remain eligible and significant at both
cutoffs. Their unrounded p-values, q-values, effects, eligibility, and
membership across all five strata are in
[table1_six_traits_before_after.tsv](table1_six_traits_before_after.tsv).

| Trait | Effect before | Effect after | q before | q after |
| --- | ---: | ---: | ---: | ---: |
| Facultatively aerobic | +0.193772 | +0.186310 | 2.780e-10 | 1.686e-9 |
| Cellulose catabolic process | -0.161938 | -0.164876 | 6.901e-7 | 3.315e-7 |
| Coccus shaped | +0.161071 | +0.156906 | 1.129e-6 | 2.469e-6 |
| Sphere shaped | +0.143819 | +0.139652 | 1.947e-5 | 3.338e-5 |
| Arsenate reduction | +0.151411 | +0.147797 | 1.947e-5 | 3.338e-5 |
| Moderately halophilic | +0.138125 | +0.131709 | 3.107e-5 | 9.411e-5 |

### Cross-framework membership

The ORA set is the union of traits with at least one significant profile
within the matched stratum, subject to fold >= 1.5. Increased and decreased
profiles both contribute. This is a trait-presence overlap, not a claim of
directional concordance. ORA retains its original global BH family of
556,516 observed-trait/profile tests; its q-values are not recomputed per
stratum. Composition retains the stratum-specific families above.

| Stratum | Cutoff | ORA traits (unchanged) | Overlap before | Overlap after |
| --- | --- | ---: | ---: | ---: |
| Pooled | q <= .001 | 278 | 6 | 6 |
| Pooled | q <= .05 | 310 | 45 | 47 |
| Feces | q <= .001 | 269 | 20 | 19 |
| Feces | q <= .05 | 306 | 89 | 86 |

Other emitted strata have zero composition overlap at either cutoff.
The manuscript's pooled six-trait overlap is retained, but any .05 integration
counts/figures and fecal .001 counts need these updated memberships. This does
not resolve the separate study-dependence or ORA sampling-unit concerns.

### Denominator parity and remaining differences

[denominator_parity.tsv](denominator_parity.tsv) verifies exact retained-total
agreement with Ilyas's pinned column/count filtering for all 5,829 keys and
with Noel's actual loader on all 4,670 paired keys. Composition and ORA use
106,993,585 retained counts across that paired corpus. Noel's six-trait
fractions match exactly in all 28,020 profile/trait comparisons, recorded in
`table1_fraction_parity.tsv.gz`. The original ORA model omits the 153
zero-total profiles from testing, leaving 4,517 testable profiles; those rows
remain in composition and Noel's paired calculations with zero fractions.

Parity is checked with frozen experiment keys and the original labels, not
with live experiment metadata. This isolates arithmetic from metadata drift.
The three intended tests still differ:

- Composition: log1p-transformed fractions, column z-scores, experiment-pair
  t-tests, no study weighting, BH separately by stratum.
- Ilyas: count-based hypergeometric upper-tail tests against the pooled paired
  background, including each profile's own counts; global BH, fold >= 1.5.
  The earlier ORA code is frozen to avoid confounding this fix with #133.
- Noel: raw fractions averaged by study and direction, sign of the study
  contrast, exact two-sided binomial test, ties excluded, BH over six traits.

The matched Noel reference in `noel_reference/` has the same six direction
signs and all six q <= .05, but only facultative aerobiosis, cellulose
catabolism, and moderate halophily pass q <= .001. This difference remains
despite identical fractions and is not a denominator mismatch. Do not replace
Noel's six-trait q-values with composition's 315-trait q-values, or vice versa.

## Artifacts and downstream audit

- [trait_before_after.tsv](trait_before_after.tsv): all 317 retained traits x
  five strata, eligibility, selection variance/score/rank, effect, direction,
  p, q, family size and .001/.05 membership; ineligible tests are NA, not p=1.
- `before/` and `after/`: all stratum statistics, selection tables, compressed
  pair-contrast matrices, condition/body-site summaries and infographics.
- [changed_significance_membership.tsv](changed_significance_membership.tsv)
  and [changed_directions.tsv](changed_directions.tsv): focused change lists.
- [cross_framework_overlap.tsv](cross_framework_overlap.tsv) and
  [cross_framework_membership.tsv](cross_framework_membership.tsv): original
  family-preserving overlap counts and trait-level membership.
- `after/cache/`: corrected full-profile totals, 60-trait metadata, full
  transformed matrix, paired matrix/metadata, required-trait checks, six-trait
  reference z-scores, profile heatmap and paired heatmap. Conventional producer
  filenames are also installed in the parent visualization directory.
- [cached_six_trait_impact.tsv](cached_six_trait_impact.tsv): changes in positive
  profile z-scores and pair-contrast signs using the full-profile reference.
  These are cache-impact diagnostics, not p-values or graph-edge counts.

The follow-up on this issue branch ports the baseline downstream renderers,
applies the retained-column definition, and replays their outputs with the
same frozen matrix and auxiliary inputs. See the
[downstream comparison](downstream/README.md) for measured changes.

| Consumer at baseline commit | Dependency | Status / required action |
| --- | --- | --- |
| `create_six_shared_trait_sankeys.py` | Cached totals and all raw retained columns | Rebuilt Sankeys, graph exports, and p/q heatmaps; mismatched totals now fail |
| `create_pd_cellulose_analysis.py` | Cached totals and all raw retained columns | Rebuilt PD summaries/figures; one exploratory Fisher test crosses raw p=.05 |
| `create_condition_body_site_guild_maps.py` | Transformed profile cache and trait metadata | Rebuilt both guilds; row keys, totals, and transformed values are checked against raw data |
| `create_condition_trait_module_heatmap.py` | Independent study/module aggregation | Shared retained-column/fraction helpers applied before aggregation and selection; rerun complete |
| `analyze_enrichment_infographic_followup.py` | Independent totals, guild fractions, cosine vectors | MediaDive excluded before those calculations; rerun complete |

The corrected downstream PNGs, SVGs, graph HTML, and backing tables are now
under `trait-visualization/trait-profile-investigation/` on this branch.
Across the full-profile six-trait cache, 168 profile/trait z-score signs and
94 experiment-pair/trait contrast signs change (strict sign comparison,
including zero-to-nonzero changes). Counts are summed across traits, not unique
profiles or graph edges. Thus unchanged pooled membership does not imply
unchanged inputs for direction-based downstream displays.
Raw taxon-to-trait annotations and unweighted presence counts are unchanged by
the denominator correction alone; any transformed thresholds or selected
sets used to build edges can nevertheless change. Do not assume an existing
PD or Sankey figure is unchanged because the pooled six-trait set is stable.
Old orphan stratum tables from the baseline Git tree are not newly tested
families; consume the corrected manifest, not every matching filename.

## Reproduction and limits

Run `just sync`, `just composition-denominator-audit`,
`just composition-input-replay /path/to/data /path/to/kg`,
`just composition-downstream /path/to/data`, then `just check`.
See [methods and captions](../../methods/composition_denominator_harmonization.md)
and [provenance.json](provenance.json) for the exact matrix, metadata, source
hashes, package versions, and fixed modeling choices. The baseline numerical
validation is enforced by the recipe, not just asserted in this report.

The matrix is the exact Git blob cited in #174, not the newer local active
matrix. The raw-input replay now verifies the original BugSigDBExports
v1.3.1 and KG-Microbe hashes and reconstructs all 1,847,793 retained
profile/trait counts exactly (5,829 x 317), including all 4,670 paired
profiles. See [raw_input_replay_provenance.json](raw_input_replay_provenance.json).
The two producers and downstream renderers are ported from the inspected
baseline because they were absent from this branch's starting `main`.

Remaining work is review and adoption: inspect the linked PR/CI, replace
manuscript and integration figure references with the corrected artifacts,
and apply any separately agreed interpretation (#179), study-dependence
(#133), taxon-influence (#190), or Sankey-layout (#195) changes independently.
Those modeling and presentation questions are not resolved by denominator
parity. Do not combine their numerical changes into the denominator-only
comparison above.
