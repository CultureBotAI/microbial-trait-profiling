# The integrated evidence table

Produced by `src/bsdbng/integrated_evidence.py`, driven by `just integrated-evidence`,
which writes `data/derived/integrated_evidence/`: the machine-readable table
(`integrated_evidence.tsv`), three tidy side-cars for the feature-to-column map,
the correction families and the assessments, a generated view with the data
dictionary (`integrated_evidence.md`) and a provenance record. It implements
[#178](https://github.com/CultureBotAI/microbial-trait-profiling/issues/178) and follows
the shared comparison basis in [`three_arm_statistics.md`](three_arm_statistics.md).

That basis requires excluding a whole experiment pair whose retained-annotation
total is zero on either side, and this table now does (#294). `Profile.fraction`
returns 0.0 for a zero total, so such a pair would otherwise enter the contrast
as a measured share of nothing, produce a contrast of exactly zero, and be
counted by the sign test as a tie - merging zero-annotation studies into the
all-tied count, which is the collapse the basis forbids. On the pinned corpus
this removes **94 of 2,335 experiments** and takes the study count from **959 to
939**. It changes the tie counts rather than the directions: facultatively
aerobic 42 to 24, arsenate reduction 343 to 325, moderately halophilic 114 to 96.
All six fraction traits remain significant at q <= 0.05.

## What problem it solves

Three frameworks report on one matrix and mean different things by "significant", and
a blank cell in any of their outputs can mean at least five different things:
untested, all tied, no annotation signal, masked from a figure, or tested and not
significant. Joining them into one "results" table is where those distinctions get
lost, and losing them is how a figure ends up implying that a method found nothing
when it was never run.

This table keeps every result in its own row, with its own family, units and states.
It is a **contract, not a verdict**: it does not pool p-values, average effects across
methods, rank methods against each other, or claim a false discovery rate for a set
defined as significant under any method. Separate corrections do not combine into one
guarantee, and the table says so on every row and in its own summary.

## The contract

One row per (method, feature, inference level). The fields, grouped:

- **Identity**: `method`, `method_version`, `inference_level`, `feature_type`,
  `trait` and `reference`, each with two identifier columns, `stratum_kind` and
  `stratum`, and `source` (the input file and a hash prefix). `trait_qualifier`
  and `reference_qualifier` carry the term ids the analysis declared;
  `trait_columns` and `reference_columns` carry the raw matrix columns those
  terms resolved to. The two differ: `GO:0030245` is carried in the matrix as
  `METPO:2000103|GO:0030245`, so only the raw column joins back to it, and a
  ratio row that named the numerator alone described half of its own hypothesis
  (#294). `integrated_evidence_features.tsv` is the same mapping in tidy form,
  one row per (feature, term, column), because a pooled guild lists several
  columns in one cell of the main table.
- **Population**: `population`, `eligibility_rule`, `denominator`, `background`, and
  the counts that are never interchangeable - `unique_studies`,
  `paired_experiments`, `profiles`, `testable_units`, `ties`. A ratio population
  is conditional on the pair, so it also carries `studies_total`,
  `studies_excluded`, `studies_partially_observed`, `experiments_eligible` and
  `experiments_excluded_zero_count`. These are blank, not zero, for a method that
  declares no such rule.
- **Effect**: `effect_value` with its `effect_units`, and `interval_low`,
  `interval_high` and `interval_type`. A standardized mean contrast, an ORA fold and
  a directional proportion are different quantities and are never put in one column.
  Where no interval was computed, `interval_type` is `unavailable` rather than blank.
  `study_effect_mean`, `study_effect_median` and `study_effect_units` describe the
  per-study distribution the test ran over rather than the tested estimate, which
  is why they survive on a row whose `effect_value` is blanked. The ratio mean
  previously existed only as a sentence inside `notes`.
- **Inference**: `p_value`, `q_value`, `correction_method`, `family_id`,
  `family_size`, `effect_threshold`, and two thresholds with their decisions -
  `q_threshold`/`decision` and `strict_q_threshold`/`strict_decision`, each with a
  `*_source` naming where the cutoff comes from (#294). A decision holds only
  inside its own family, and neither threshold is native to all three families.
- **States**, kept apart because they overlap: `test_status` (`tested`, `untested`,
  `insufficient_support`, `all_tied`, `zero_annotation_signal`, `not_sampled`,
  `source_not_supplied`), `decision` (`significant`, `not_significant`,
  `undecided`), `display_status` (`drawn`, `masked_below_min_studies`, `not_drawn`,
  `not_applicable`) and `zero_annotation_state`.
- **Standing**: `selection_note`, and one status field per open assessment -
  `calibration_status` (#182), `sensitivity_status` (#183), `coverage_status` (#184),
  `orientation_status` (#188) and `taxon_influence_status` (#190). Each reads from
  that analysis's own output when it exists and is `pending` otherwise, and each
  carries the state for **this row's own procedure**. An assessment does not have
  a single verdict: #182 reports one per procedure and separately names the arms
  it does not cover, which read `not_calibrated_here`, and #183 reports counts per
  scenario, which is not a verdict at all and says so. A concatenation of every
  verdict copied onto every row attached a calibration result to procedures that
  assessment had explicitly excluded (#294). Per-scenario and per-contrast-class
  detail stays in `integrated_evidence_assessments.tsv`, attached to the scope it
  describes. A status file that is not a JSON object, or that holds something
  other than a string where a verdict belongs, reads `unreadable` with the
  reason rather than crashing the run or being written through verbatim.

Two rules are enforced in code rather than left to the reader. A row that was not
tested has its effect, interval, p and q blanked, because a source may carry a
placeholder - the fraction analysis reports share 0 and a q derived from p = 1 for an
all-tied trait before #165 - and a placeholder read as an effect is the confusion this
table exists to prevent. And a method whose output was not supplied gets a row per
feature with `source_not_supplied`, so that its absence is visible rather than being
read as a negative result.

## What is in it now

On the pinned v1.3.1 corpus, **1,912 rows: 40 method results and 1,872 figure
cells**.

| Method | Records | Family | Tested | q <= 0.05 | q <= 0.001 |
| --- | --- | --- | --- | --- | --- |
| directional consistency, fraction sign test | 6 | BH over six traits | 6 | 6 | **3** |
| directional consistency, log-ratio sign test | 15 | BY over 15 pairs | 15 | 3 | **0** |
| oxygen guild McNemar | 1 | single uncorrected test | 1 | no q | no q |
| composition arm, native transformed contrast | 6 | - | 0 | not supplied | not supplied |
| enrichment arm, profile hypergeometric | 6 | - | 0 | not supplied | not supplied |
| aggregated profile support | 6 | - | 0 | not supplied | not supplied |
| oxygen guild cell support | 1,872 | - | 0 | not a test | not a test |

The guild McNemar is a single uncorrected test. It carries p = 2.2e-06 and leaves
q blank, so it has no threshold decision in either column: a decision in this
table is a decision on a corrected value and there is no corrected value here.
An earlier version wrote the uncorrected p into `q_value` as well, which made it
look like a one-member FDR family (#294).

The cell rows are the oxygen guild figure's own inventory, one per cell of the
condition-by-body-site grid at the figure's parameters: **316 drawn, 7 drawn as
sampled-without-signal, 209 masked for too little support and 1,340 never
sampled**. They report support and display, never a test, which is why they are
counted in the generated view rather than listed row by row. Before they
existed, nothing in the table set `drawn` or `masked_below_min_studies` and
`not_sampled` was never emitted at all, so three of the four display states named
something the table could not produce.

Both columns are reported because **no single cutoff is native to all three
families** (#294). An earlier version applied one `--q-threshold` to all of them
and reported only the q <= 0.05 count, which read as a property of the results
rather than of the threshold chosen for display. Every row now carries a decision
at each, and each threshold names its source: q <= 0.05 is chosen by this
integration for display, matching the overlay of PR #173, while q <= 0.001 is the
manuscript's rule as the enrichment arm's committed run applies it. The
composition arm publishes q with no cutoff at all, so neither is that arm's own
rule. At the stricter threshold coccus shaped (q = 0.0073), sphere shaped
(0.0029) and arsenate reduction (0.0038) drop out of the fraction family, and the
whole ratio family does.

The three unsupplied methods are the point of the design rather than a gap in it. The
composition arm's native output moves under #174/#197, the enrichment arm's rests on
the urn of #131, and the aggregated support rule lives on the #173 branch. Each is
recorded as not supplied with the reason.

**There is no way to supply them yet.** The command line takes the profile
matrix, the experiments table, the full dump, the two thresholds and the five
assessment status files, and nothing else; there is no `--composition-stats`,
`--ora` or `--support` input. Importing an arm means adding its reader and
deciding how its effect, family and eligibility map onto these columns, which is
the work those issues describe - not pointing an existing recipe at a file. An
earlier version of this note said otherwise (#294).

The placeholder rows also cover the six Table 1 traits only. That is their
declared scope, not an inventory of what each method would report: the enrichment
arm covers every retained trait in the matrix, and the corpus-wide directional
screen of #276 is not integrated here at all.

At q <= 0.05, six of the six fraction results and three of the fifteen ratio
results are significant in their own families; at q <= 0.001 it is three and
none. Those are not nine findings, or four: the ratios are
measured on their own eligible populations, the families are corrected separately,
and the two ask different questions. The table reports them side by side precisely so
that a reader cannot silently add them up.

## Standing of the results in it

Every row carries the same caveats as fields rather than as prose. **On `main`,
every one of the five status columns reads `pending`**, because none of the five
assessments has landed. The findings below come from open pull requests and are
reported here as what those branches currently say, not as settled properties of
these results; an earlier version of this note stated them as established
(#294).

- **Selection.** The six Table 1 traits were chosen on this corpus by the composition
  arm, which every row records in `selection_note` regardless of any assessment.
  These results are exploratory on that ground alone.
- **Calibration (#182, open as PR #240).** Its branch measures what selecting on
  the data and then testing the selection costs - a per-test error of roughly 0.08
  against a nominal 0.05, a lower bound - and reports the procedures this
  repository owns as holding their nominal rates under the paired null. The same
  output names the enrichment arm's urn, the composition arm's native test and the
  aggregated support rule as `not_calibrated_here`, and those rows say so rather
  than carrying any verdict. Until #240 merges, every `calibration_status` on
  `main` is `pending #182`.
- **Population (#183 as PR #237, #188 as PR #238).** #237's branch reports that
  most fraction directions do not survive restriction to a ratio's eligible
  population, and #238's gives a disease-versus-healthy subset in which the
  directions do hold. Both are populations rather than corrections, and #237
  reports per scenario without reaching a verdict at all. Neither is on `main`.
- **Coverage and influence.** #184 and #190 remain pending; until they land, no row
  should be read as accounting for annotation coverage or taxonomic concentration.
- **Guild grouping.** The guild rows previously claimed the oxygen guilds were
  fixed a priori. Nothing in this repository records where or when they were
  declared, and the aerobe guild contains METPO:1000608, which is also one of the
  six traits selected on this corpus - so that grouping is not independent of the
  selection either. The claim is removed rather than restated (#294).

## Reproducing

`just integrated-evidence` on the pinned inputs writes six files into
`data/derived/integrated_evidence/`: the table (`integrated_evidence.tsv`), the
feature-to-column map (`integrated_evidence_features.tsv`), the correction
families with their declared sizes (`integrated_evidence_families.tsv`), each
assessment against the procedure and scenario it assessed
(`integrated_evidence_assessments.tsv`), the generated view
(`integrated_evidence.md`) and the provenance record. The family
export is the #236/#289 requirement that no exported q appears without the
family it was corrected in and that family's size; the run fails if a declared
size and the rows exported for that family disagree, because a q corrected over
a different set of tests than the table shows cannot be read. The provenance
record carries the zero-total pair exclusion as counts and a rule, both decision
thresholds with their sources, the resolved feature columns, and any declared
term that matched no column.

It reads
`data/derived/calibration/calibration_status.json`,
`data/derived/table1_composition_sensitivity/sensitivity_status.json` and
`data/derived/contrast_orientation/contrast_orientation_provenance.json` when they
exist, and `--coverage-status` and `--influence-status` take the #184 and #190
outputs once those exist. The provenance record carries the code revision, input
hashes, the corpus cascade, the assessment statuses and the per-method summary.

`tests/test_integrated_evidence.py` covers the states the contract exists to
preserve: an all-tied trait, a pair with no eligible experiment, a method that was
not supplied, two methods keeping separate families and units, the four figure-cell
support states taken from the real matrix, a masked cell leaving a tested result
alone, each assessment reaching only the procedure it assessed, ten malformed
status files, both thresholds deciding independently, the ratio populations and
effect summaries reconciling with their source run, the McNemar keeping its p with
q blank, an undefined odds ratio coming out missing, and the CLI end to end with
its exports and provenance record. Each test was checked against a deliberately
broken version of the code it covers.

The table, the contract and this description were drafted with the assistance of a
coding agent (Claude Code) and require confirmation by the method owners before any
part of it goes into the paper.
