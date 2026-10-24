# The integrated evidence table

Produced by `src/bsdbng/integrated_evidence.py`, driven by `just integrated-evidence`,
which writes `data/derived/integrated_evidence/`: the machine-readable table
(`integrated_evidence.tsv`), three tidy side-cars for the feature-to-column map,
the correction families and the assessments, a generated view with the data
dictionary (`integrated_evidence.md`) and a provenance record. It implements
[#178](https://github.com/CultureBotAI/microbial-trait-profiling/issues/178) and follows
the shared comparison basis in [`three_arm_statistics.md`](three_arm_statistics.md).

That basis requires excluding a whole experiment pair whose retained-annotation
total is zero on either side, and this table's own tests now do (#294). `Profile.fraction`
returns 0.0 for a zero total, so such a pair would otherwise enter the contrast
as a measured share of nothing, produce a contrast of exactly zero, and be
counted by the sign test as a tie - merging zero-annotation studies into the
all-tied count, which is the collapse the basis forbids. On the committed corpus
(identified under "What is in it now") this removes **55 of 2,257 experiments**
and takes the study count from **932 to 922**. It changes the tie counts rather
than the directions: facultatively aerobic 32 to 24, arsenate reduction 324 to
316, moderately halophilic 101 to 93. All six fraction traits remain significant
at q <= 0.05. Since #783 the oxygen guild figure, and with it the figure-cell
inventory, leaves those pairs out too. Every test and figure cell computed here
shares one population, chosen by the repository owner on 2026-10-22 for the
reasons given in [`three_arm_statistics.md`](three_arm_statistics.md). The rows
say so: their `population` reads "human, paired experiments; pairs with a zero
retained-annotation total excluded". The imported native rows keep their
sources' populations, and the provenance record lists both sets. Before, the
fraction and McNemar rows read "human, paired experiments" while running on the
excluded set.

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

On the committed corpus, **4,200 rows: 40 method results and 4,160 figure
cells**. That is the run that wrote `data/derived/integrated_evidence/`. Its
profiles are sha256
`519a095b09e179dc92fb51f9d892984a52dd100dece4dadf6c94b036a5323e6c` (the committed
`study_trait_profiles_wide.tsv.gz`), its experiments
`16271237a70f61c8195d9f2acdd16f71e88f65859160f95bb734626bc1328e06`, and its
BugSigDB export `2026-04-24_00:41_UTC` (full dump
`af84bfca7bfe76913a20c1fedbaadbfad4e4bf78605a2cfd7d00eb0c0fdd5788`). The
provenance record carries every other input's hash. Main's code at 35ece98
(#808) writes every number in these files. #810 regenerated them on main's code
at 68d4e27 with its own change, which reworded one sentence of the view, added
the provenance fields #813 introduced and changed no number. The record's
`code_commit` names the branch commit the run was made from, which a squash
merge does not keep on main (#643).

The final rerun at 511c891 did not regenerate these files: their assessment and
arm inputs were not rebuilt there, so they still describe the `off`-scope profiles
above, with 317 retained traits rather than the 344 of the primary `phenotype` scope.
Regenerating them with every input is #1003.

| Method | Records | Family | Tested | q <= 0.05 | q <= 0.001 |
| --- | --- | --- | --- | --- | --- |
| directional consistency, fraction sign test | 6 | BH over six traits | 6 | 6 | **3** |
| directional consistency, log-ratio sign test | 15 | BY over 15 pairs | 15 | 3 | **0** |
| oxygen guild McNemar | 1 | single uncorrected test | 1 | no q | no q |
| composition arm, native transformed contrast | 6 | the arm's BH over 313 support/variance-eligible traits | 6 | 6 | **4** |
| enrichment arm, profile hypergeometric | 6 | global BH over 1,430,938 profile tests | 0 (6 profile summaries) | no trait-level q | no trait-level q |
| aggregated profile support | 6 | - | 0 | not a test | not a test |
| oxygen guild cell support | 4,160 | - | 0 | not a test | not a test |

The guild McNemar is a single uncorrected test. It carries p = 3.6e-06 over 922
studies and leaves q blank, so it has no threshold decision in either column: a
decision in this table is a decision on a corrected value and there is no
corrected value here. Since #783 the figure's own McNemar is the same test on the
same population: 410 and 287 discordant studies, p = 3.6e-06. Before, the figure
kept zero-total pairs and reported p = 2.0e-06 over 932 studies.
An earlier version wrote the uncorrected p into `q_value` as well, which made it
look like a one-member FDR family (#294).

The cell rows are the oxygen guild figure's own inventory, one per cell of its
condition-by-body-site grid: 52 condition rows (51 conditions and the pooled
Other) by 20 body-site blocks (19 sites and the pooled Other), each block four
guild-by-direction columns wide, of which
**340 are drawn, 10 drawn as sampled-without-signal, 302 masked for too little
support and 3,508 never sampled**. The figure's own builder,
`guild_maps.build_figure_matrix`, makes them at the figure's parameters and over
its population. Since #783 that population is the tests' shared basis, so the 55
zero-total pairs are left out of the cells as they are of the tests. With those
pairs kept, the grid was a different one: it had rows for Hypertension and
Non-alcoholic fatty liver disease where it now has Attention deficit-hyperactivity
disorder and Periodontitis, and it counted 20 sampled-without-signal, 296 masked
and 3,504 never sampled. An earlier version rebuilt the grid here, at an 8-site cut
the figure had left for 20, so its 1,872 cells described a grid the figure never
drew (#770). They report support and
display, never a test, which is why they are
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
rule. At the stricter threshold coccus shaped (q = 0.0083), sphere shaped
(0.0021) and arsenate reduction (0.0027) drop out of the fraction family, and the
whole ratio family does.

The composition and enrichment arms are read from their own exports (#551).
`--arm1-pooled-stats` takes the composition arm's native study-level t-tests,
which keep their own BH family over the arm's support/variance-eligible traits
(313 here). This table exports six of them (`source_subset`), and the pair-level
tests stay secondary fields. `--arm2-enrichment-results` and
`--arm2-depletion-results` take the enrichment arm's complete profile-level
exports, summarized per trait (`tested_profiles_summarized`). They come from a
global BH family over every profile test, and no trait-level q or decision is
defined from them. The aggregated profile support rows are descriptive and stay
untested; the support rule itself lives on the #173 branch. A method whose
export is not supplied still gets `source_not_supplied` rows. An earlier
version of this note said none of the three could be supplied.

These rows cover the six Table 1 traits only. That is their declared scope, not
an inventory of what each arm reports: the composition arm's family covers 313
traits, and the enrichment arm covers every retained trait in the matrix. The
corpus-wide directional screen (#255, PR #276) is attached to the six
fraction rows only, as their `full_screen_*` columns from its own 317-trait
family (the `off` scope of these files; 344 traits under `phenotype`), with BY as
its primary correction.

At q <= 0.05, six of the six fraction results and three of the fifteen ratio
results are significant in their own families; at q <= 0.001 it is three and
none. Those are not nine findings, or four: the ratios are
measured on their own eligible populations, the families are corrected separately,
and the two ask different questions. The table reports them side by side precisely so
that a reader cannot silently add them up.

## Standing of the results in it

Every row carries the same caveats as fields rather than as prose. The statuses
below are those of the committed rows (`integrated_evidence_method_results.tsv`
and `integrated_evidence_assessments.tsv`), written by the run named under "What is
in it now" on profiles `519a095b…`, experiments `16271237…` and BugSigDB export
`2026-04-24_00:41_UTC`. Each status reports one assessment's state for the row's own
procedure. None of them is a verdict on a trait, and none makes a selected result
confirmatory ([`biological_interpretation.md`](biological_interpretation.md)). An
earlier version of this note said every status column read `pending` (#800), and
one before that stated the assessments' branch findings as established (#294).

- **Selection.** The six Table 1 traits were chosen on this corpus by the composition
  arm, which every row records in `selection_note` regardless of any assessment.
  These results are exploratory on that ground alone.
- **Calibration (#182).** The status file is `calibration-v3` at code commit
  8b280a3, sha256 `3007e1f1…`, bound to this table's corpus (`matched`).
  - It reads `pass #182` on 22 rows: the six fraction sign tests, the 15 log-ratio
    sign tests and the McNemar test. The criteria are the per-test level (a sanity
    check), the family error and, except for the McNemar, the reported power. The
    side-car keeps each criterion and the completed replicates: 1,400 for each
    complete-null stream and 400 for each injected effect.
  - It names the enrichment arm's urn, the composition arm's native test and the
    aggregated support rule as `not_calibrated_here`, so those 18 rows carry no
    verdict.
  - What it does not settle: a pass says these procedures hold their error rates
    on traits fixed in advance. The same file reports a screened fraction sign test,
    which selects traits on the data it then tests, as `fail` on both criteria, with
    the selection inflation measured. No row here is that procedure, but the six
    traits were selected on this corpus.
  - A status file that gives one verdict string per procedure is refused as
    unreadable. That was the contract under which a five-replicate run with a failed
    family criterion reached this table as a pass (#240 review, 2026-10-03).
- **Sensitivity (#183).** The status file, sha256 `043a0e95…`, declares no input
  hashes, so its binding reads `inputs_not_declared` and is accepted pending #816.
  - It gives a verdict string for the log-ratio sign test only: `scenarios reported,
    no verdict #183`, on the 15 ratio rows. The other 25 rows read `pending #183`.
  - The side-car keeps one record per correction family, with its scenario, size
    and correction. These are read from #237's `family_status` records because its
    `scenario_status` pools the two S4 families (6 fractions under BH, 15 ratios under
    BY) into tested=21. For example, all six fractions reach q <= 0.05 in the
    baseline scenarios (S0, S1) but two do on the six traits' common support (S4).
    Three of the 15 ratios do on complete cases (S2) and four on common support (S4).
  - What it does not settle: it reports counts per scenario and reaches no verdict,
    so no row is assessed as robust or fragile to these population choices.
- **Orientation (#188).** The status file is `contrast-orientation-v5` at code commit
  8b280a3, sha256 `12de1b52…`, bound to this table's corpus (`matched`).
  - All 40 rows read `contrast population described, hand_review_agreed, no verdict
    #188`: the 60 hand-reviewed experiments all agree with the classifier.
  - The side-car keeps the contrast classes with their orientation and counts. For
    example, 1,005 experiments in 574 studies are natively oriented disease vs
    healthy, and one is reversed.
  - What it does not settle: it describes the population and gives no verdict. These
    rows run on the historical paired population, not on the disease-vs-healthy
    subset, so none of them is a disease-versus-control result. The verified-subset
    comparison is `disease_control_comparison.tsv` in
    `data/derived/contrast_orientation/`, which is not imported here.
- **Coverage (#184).** Nothing produces this assessment yet, so all 40 rows read
  `pending #184`. No row should be read as accounting for annotation coverage.
- **Influence (#190).** All 40 rows read `pending #190: influence evidence not
  supplied`, because no influence output was given to this run. No row should be read
  as accounting for taxonomic concentration.
- **Validation data (#191).** Whether measured inflammation or SCFA data could
  validate the trait directions independently is assessed in
  [`validation_feasibility.md`](validation_feasibility.md). It identifies evidence
  pointers, not validation: an independent check would need sample-level data, so
  no row here carries a validation status from it.
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
family it was corrected in and that family's size. The Table 1 direction-shift TSV
written by `table1_trait_maps.write_shift_table` now carries the same
`family_id`, `family_size` and `correction_method` columns, which closes the
source-export half of that requirement that #236 deferred here. The run fails if a declared
size and the rows exported for that family disagree, because a q corrected over
a different set of tests than the table shows cannot be read. The provenance
record carries the zero-total pair exclusion as counts, the kept study IDs and a rule, both decision
thresholds with their sources, the resolved feature columns, and any declared
term that matched no column.

It reads
`data/derived/calibration/calibration_status.json`,
`data/derived/table1_composition_sensitivity/sensitivity_status.json` and
`data/derived/contrast_orientation/contrast_orientation_provenance.json` when they
exist, and `--coverage-status` and `--influence-status` take the #184 and #190
outputs once those exist. Under `--require-complete-sources`, which the cross-arm
comparison stage passes, a native table or a calibration, sensitivity, orientation or
influence status that is not supplied or does not exist is an error rather than a
pending row (#814). So is a supplied status that reads as unreadable or whose binding is
`mismatch_or_missing_identity` or `missing_producer_identity` (#844). A status that
declares no input hashes is still accepted, pending #816. Nothing produces the #184
coverage status yet, so it may be left out and then reads as not supplied; a coverage
path that is supplied must exist.
The provenance record carries the code revision, input
hashes, the corpus cascade, the assessment statuses and the per-method summary.

`tests/test_integrated_evidence.py` covers the states the contract exists to
preserve: an all-tied trait, a pair with no eligible experiment, a method that was
not supplied, two methods keeping separate families and units, the four figure-cell
support states taken from the real matrix, the cell inventory, in all four
states, and the McNemar equal to the figure's own cell-status table and McNemar
on one corpus, no
zero-total pair reaching the figure's matrix or either McNemar whatever the caller
passes, a masked cell leaving a tested result
alone, each assessment reaching only the procedure it assessed, ten malformed
status files, both thresholds deciding independently, the ratio populations and
effect summaries reconciling with their source run, the McNemar keeping its p with
q blank, an undefined odds ratio coming out missing, and the CLI end to end with
its exports and provenance record. Each test was checked against a deliberately
broken version of the code it covers.

The table, the contract and this description were drafted with the assistance of a
coding agent (Claude Code) and require confirmation by the method owners before any
part of it goes into the paper.
