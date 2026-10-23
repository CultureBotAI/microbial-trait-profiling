# Selection scan: how unusual is the PD, HIV and CRC overlap? (#835)

Parkinson's disease (PD), HIV infection (HIV) and colorectal cancer (CRC), each with
cellulose catabolic process decreased, were picked as exploratory leads from this
BugSigDB corpus (`EXPLORATORY_LEADS` in `src/mtp/deck_visualizations.py`). Cellulose
was chosen afterwards, as the top decreasing trait. The criterion that picked the three
conditions is not recorded in the repository, so every p-value in
[`cellulose_decr.md`](cellulose_decr.md) is nominal.

This scan puts that selection in context. It repeats cellulose_decr's overlap
statistics and its whole-study-signature null (U6) for every triple of eligible
conditions and for the top decreasing traits. It then reports where PD, HIV and CRC
fall among all triples, and where cellulose falls among the traits.

Like cellulose_decr, the scan counts curated signature taxa and their KG-Microbe
annotations. It does not measure cellulose degradation or any other physiology
([`biological_interpretation.md`](biological_interpretation.md)).

## Declared rule

This section was written on 2026-10-22 and committed on its own, before the scan's
code and before any run on real data. Its parameters are in
[`cellulose_decr_selection_scan_protocol.json`](cellulose_decr_selection_scan_protocol.json),
which the code reads. A later change to the rule is added to that file's `amendments`,
with its date and what was known when it was made, rather than edited in here.

### Population

- **Source.** The scan reads one verified core run (`--core-run`) whose contrast
  population is `disease_vs_healthy_tier1_2_subject_controls`, and refuses a run of any
  other population. That population keeps paired human disease-versus-healthy
  experiments with subject controls (`contrast_type` `case_vs_control_subject`), Tier 1
  healthy/normal or Tier 2 condition-negative controls, and 16S or shotgun sequencing.
  Reversed experiments are relabelled, so *decreased* means lower in the disease group
  ([`docs/workflow.md`](../docs/workflow.md), [`contrast_orientation.md`](contrast_orientation.md)).
- **Not cellulose_decr's own inputs.** cellulose_decr reads the checkout's
  `data/derived` tables, including its own contrast manifest. Its `control_tiers_1_2`
  definition keeps Tier 1–2 controls without requiring `case_vs_control_subject` or a
  sequencing type. #818 found that the two populations give different set sizes. The scan uses the core's population, so that
  its numbers agree with the core run and the three arms. Its PD, HIV and CRC sets
  can therefore differ from those in `cellulose_decreased_taxa_catalog_tier1_2.tsv`.
- **Products.** From the core run: `experiments` (each experiment's condition and body
  site), `reported_taxa` (signature taxa, directions already oriented),
  `count_profiles` (the profiles, and the trait screen), and `taxa_matrix` with
  `feature_metadata` (annotations and trait labels).
- **Stool** is the body site exactly `Feces`, as in cellulose_decr.
- **A condition** is the experiment's exact BugSigDB condition label, as in
  cellulose_decr's primary analysis. Synonym labels are not merged. A multi-condition
  label such as `Constipation,Parkinson's disease` is a condition of its own. A blank
  label is not a condition.

### Eligible conditions

- A condition is eligible when its stool disease-versus-healthy comparisons come from
  at least **k = 3 studies**. It then has at least three such comparisons.
- **Why 3.** It is the repository's existing minimum support for a condition in one
  body site. The guild maps, the Table 1 maps and the integrated-evidence inventory
  mask a condition × body-site cell backed by fewer than three studies
  (`FIGURE_MIN_STUDIES` in `src/mtp/guild_maps.py`, `--min-studies` in
  `src/mtp/table1_trait_maps.py`). An eligible condition is such a cell, with stool as
  the body site.
- **Why studies.** That minimum counts studies, and U6 draws studies. Three
  comparisons from one study are one cohort and one U6 unit.
- Every eligible condition enters every triple. There is no other filter.

### Traits

- **The original ranking is not recoverable.** No code, methods note, output or commit
  message in the repository records how cellulose was ranked first among decreasing
  traits (searched 2026-10-22). As #835 asks for that case, the scan uses the Noel
  arm's declared directional screen instead
  ([`trait_direction_screen.md`](trait_direction_screen.md),
  [`trait_screen_manifest.json`](trait_screen_manifest.json)).
- **The screen** runs on the same core's count profiles, with the screen's own code
  (`mtp.trait_screen`): its candidate rule, one vote per study, and a two-sided exact
  binomial test. It votes over every body site of the population, while the sets
  below are stool-only. (Clarified on 2026-10-22, #872; nothing computed changes.)
- **A decreasing trait** is one the screen tested in which more studies shift down
  than up (share up below 0.5). Its annotations make up a larger share of the
  decreased signatures than of the increased ones.
- **Order.** Decreasing traits are ranked in the screen report's display order:
  p-value ascending, then more scored studies, then trait ID.
- **The trait set** is the top **N = 10**, plus cellulose catabolic process
  (`METPO:2000103|GO:0030245`) when it is not among them. A trait added this way is
  flagged.
- **Why 10.** cellulose_decr already reports a ten-trait list, its benchmark's top ten
  (`benchmark_top_10_excluding_growth_media` in its summary JSON; README section 7).

### Statistics

For each trait and each triple of eligible conditions:

- **A condition's set** is cellulose_decr's: the matrix taxa that bear the trait (a
  positive matrix count, so the taxon or a descendant has the KG edge) and that at
  least one of the condition's stool decreased signatures reports, at any rank.
- **The three-way core** is the number of taxa in all three sets.
- **The mean pairwise Jaccard distance** is the mean of d_J = 1 − |X ∩ Y| / |X ∪ Y|
  over the three pairs. A pair whose union is empty has no distance, and the mean is
  taken over the pairs that have one, as in U6's draws. A triple whose three sets are
  all empty is not tested.
- **Primary and secondary.** The mean pairwise distance is the primary statistic: its
  U6 p-value is cellulose_decr's specificity result (p = 0.19 over all contrasts). The
  core is secondary.

### Null: U6

- **U6 is cellulose_decr's whole-study-signature null**, and the scan calls its code
  (`u6_units`, and the random draws of `signature_null`):
  - The pool is every study in the population with a stool decreased signature, in any
    condition, eligible or not.
  - Each study has one unit per condition. A unit's taxa are the trait-annotated taxa
    that its stool decreased signatures report.
  - For a triple, each iteration draws as many distinct studies as the three
    conditions have units. It takes one unit of each study at random and splits them
    among the three synthetic conditions, so no study contributes to two of them.
- **Permutations and seed:** 100,000 iterations, the number cellulose_decr's U6 uses,
  with seed 20261022.
- **p-values:** p = (1 + draws at least as extreme) / (1 + draws). For the mean
  distance, a draw is as extreme when its distance is at most the observed one (within
  1e-12). For the core, it is as extreme when its core is at least the observed one. A
  draw whose three synthetic sets are all empty has no distance and counts in the
  denominator only, as in cellulose_decr.
- **Shared draws.** One sequence of U6 draws serves every triple and every trait. Each
  eligible condition takes its own block of each draw's studies, so every triple's
  null is exactly U6. A triple whose blocks would overlap gets its own U6 draws. Each
  p-value is valid on its own, but the p-values of different triples and traits share
  draws and are not independent.
- **Untestable triples.** A triple whose conditions have more units between them than
  the pool has studies cannot be drawn. It is not tested.

### What is reported

- **One row per triple and trait.** Each row gives the three sets' sizes, their union,
  the core and the mean pairwise distance; each condition's U6 units; and U6's
  expected value, SD, z and p for both statistics. It also gives the number of studies
  that report stool decreased signatures in two or more of the triple's conditions, and
  the test status.
- **PD/HIV/CRC's rank among all tested triples,** for each trait and statistic, by U6
  p, smallest first. The lead triple uses the exact labels `Parkinson's disease`,
  `HIV infection` and `Colorectal cancer`. Ties are reported as a range of ranks,
  together with the share of triples whose p is at most PD/HIV/CRC's. If any of the
  three conditions is not eligible, the scan still runs and says so, with no rank.
- **The null-calibrated distribution.** For each trait and statistic, the distribution
  of U6 p-values over all tested triples is given as quantiles and as counts in twenty
  bins of width 0.05. If overlap were no larger than U6 expects between any three
  conditions, the mean-distance p-values would be roughly uniform. The core's would not:
  a triple with no taxon in all three sets has p = 1, so its p-values pile up at 1.
  Mass near 0 shows how common overlap beyond U6 is between arbitrary conditions.
  (Corrected on 2026-10-22, #872: this said "roughly uniform" for both statistics.)
- **Cellulose in context.** For each trait, the table gives PD/HIV/CRC's U6 p-values and
  ranks among triples, and cellulose's rank among the traits by those p-values.
- **Study-disjoint triples (declared sensitivity).** U6 draws distinct studies, so it
  cannot reproduce two conditions that share a study. The observed sets can share one:
  for example, a study comparing Crohn's disease and ulcerative colitis against the
  same controls. In cellulose_decr's data, no study contributes to two of PD, HIV and
  CRC. The ranks and the distributions are therefore also given over the triples whose
  conditions share no study.
- **No multiplicity correction.** Every p-value is nominal. Selection-aware inference
  for the leads (a split-sample replication, or Benjamini–Bogomolov) belongs to the
  unified multiple-testing layer, which reads the per-triple rows.

## Amendments

Both amendments were made on 2026-10-22, after the review of PR #861 (#872) and before
any run on real data. They are recorded in the protocol file's `amendments`, and the
declared rule above is left as it was committed. When they were made, no scan result on
real data existed. The eligible counts on the pinned export and timings on synthetic
fixtures were known.

1. **The strain-rollup scope is pinned.**
   - The rule did not name the scope of the core run. After #828 the default core is
     `phenotype`, and choosing between `off` and `phenotype` after seeing both would
     reopen the forking path this scan exists to close.
   - The primary scope is `off`: the annotation the leads were chosen from, and the
     scope of the talk's numbers. `phenotype` is a declared sensitivity.
   - A core of any other scope is refused. A core that records no `strain_rollup`
     parameter, built before #828, is read as `off`.
   - The summary, the provenance and `--plan` report the scope and its role.
2. **Triples are also ranked by U6 z.**
   - PD is the largest condition, so its triples have tighter nulls, and smaller p for
     the same excess overlap. A rank by z separates large from unusual.
   - The z rank sits beside the p rank, over the same two reference sets, with ties
     as a range. It ranks the mean pairwise distance by ascending z and the core by
     descending z, so more overlap comes first.
   - A triple whose null has no spread has no z and no z rank.
   - The trait table also ranks the traits by the share of triples whose z is at least
     as extreme as the lead's.

## How the scan runs

This section and the next describe the implementation. They were written with the code,
after the rule above was committed, and do not change it.

- **Command.** `just cellulose-decr-selection-scan RUN OUTPUT` runs
  `mtp.cellulose_decr.selection_scan` (`src/mtp/cellulose_decr/selection_scan.py`) on a
  core run's `run.json`. It writes to a new or empty directory and refuses any other.
- **Dry run.** `--plan` verifies the core and prints the eligible conditions, the trait
  set, the number of triples and a timed estimate of the run. It writes nothing and runs
  no null.
- **Off-protocol runs.** `--permutations` and `--seed` exist for dry runs and tests. A
  run with values other than the protocol's records `follows_protocol: false`.
- **Verification.** The core run is verified before the scan reads it, and again before
  the scan writes its consumer record (`core_consumption.json`), as for the other core
  consumers. `just core-verify-consumer` checks that record later. A core of another
  population or direction basis is refused.
- **Shared draws.** Each U6 draw is one call to cellulose_decr's draw function
  (`signature_draws`, which `signature_null` also uses). It gives a random order of the
  pool's studies and one random unit for each study.
  - The eligible conditions, in alphabetical order, take consecutive blocks of that
    order, one block per condition, as large as its number of units.
  - A triple's three blocks hold distinct studies, so its draws are U6 draws. The
    synthetic sets of every condition are built once per draw, and every triple and
    trait is scored from them.
  - When the eligible conditions together have more units than the pool has studies,
    the blocks wrap around the order. A triple whose blocks then overlap is drawn on its
    own, with the same code and the seed (seed, triple index).
- **Agreement with cellulose_decr.** Laid out as (PD, HIV, CRC), the scan gives
  cellulose_decr's sets, statistics and U6 p-values exactly, for the same seed and
  number of draws (`tests/test_cellulose_decr_selection_scan.py`).
- **Expected values.** The scan's U6 expected distance and SD are taken over the draws
  that have a distance. cellulose_decr takes them over every draw, so its expected value
  is missing when a draw has none. The number of draws without a distance is reported.
- **Runtime.** The null takes time in proportion to the number of triples. On a
  synthetic corpus of 34 eligible conditions (5,984 triples), 11 traits and 100,000
  draws, it took about 4.5 minutes on one core of an Apple-silicon laptop. Verifying a
  Tier 1+2 core run takes a few more minutes, and happens twice.

## Outputs

| File | Contents |
| --- | --- |
| `selection_scan_triples.tsv.gz` | One row per trait and triple: the three conditions, set sizes, union, core and mean pairwise distance; U6 units; the studies the conditions share; U6 expected value, SD, z and p for both statistics; the number of draws and where they came from; the test status; whether it is the lead triple; and the triple's rank within its trait |
| `selection_scan_traits.tsv` | The trait set with its screen statistics and U3 size; PD/HIV/CRC's results and ranks over all tested triples and over study-disjoint triples; and each trait's rank among the traits by PD/HIV/CRC's p and by the share of triples at most as large |
| `selection_scan_conditions.tsv` | Every named condition with a stool comparison: comparisons, studies, U6 units, eligibility and layout position |
| `selection_scan_distribution.tsv` | U6 p-values over triples in twenty bins, per trait, statistic and reference set |
| `selection_scan_summary.json` | The rule as run, the counts, the trait table, and quantiles of p and z per trait, statistic and reference set |
| `selection_scan_provenance.json` | Code commit, script hashes, protocol and methods hashes; the core run's ID, manifest hash and parameters; each input's path within the core run, with its hash; permutations, seed and package versions |
| `core_consumption.json`, `core_consumption.coverage.tsv` | The binding to the verified core run |

The triple table is gzip-compressed with a fixed timestamp, so the same core, rule and
seed give the same bytes.

## Study splits (#854)

The selection-aware lane of #835 (#854) selects on some studies and tests on the others,
over repeated study splits. The scan exposes both steps as functions in
`mtp.cellulose_decr.selection_scan`. The CLI runs `select` on every study, so a split's
selection applies the declared rule above, unchanged, to its training studies. The split
protocol itself (the number of splits, which candidates to carry forward, the
multiplicity family and the reporting) belongs to #854.

| Function | Signature | Returns |
| --- | --- | --- |
| `load_core_inputs` | `(core_run: Path) -> tuple[CoreBundle, ScanInputs]` | Verifies a Tier 1+2 core run and reads its population once. A core of another population is refused. `ScanInputs.studies` lists the population's studies; `ScanInputs.restricted(studies)` keeps whole studies. |
| `select` | `(inputs, studies=None, *, protocol=None, permutations=None, seed=None, always_included=True) -> Selection` | The declared rule on `studies` only. `Selection.conditions` is the condition table, `Selection.traits` the trait set from those studies' screen, and `Selection.triples` one row per triple and trait, with the columns `TRIPLE_COLUMNS`. `Selection.result` holds the arrays and `Selection.studies` the studies used. |
| `rank_candidates` | `(triples, statistic="mean_pairwise_jaccard_distance") -> DataFrame` | The tested rows in a fixed order, numbered `candidate_rank`: U6 p ascending, then z (more overlap first), then trait ID and the three conditions. It is an order, not a selection. |
| `evaluate` | `(inputs, studies, hypotheses, *, protocol=None, permutations=None, seed=None) -> DataFrame` | One row per hypothesis, `(condition, condition, condition, trait ID)`, in the given order and numbered `hypothesis`, with the columns `TRIPLE_COLUMNS` except the two within-trait ranks. |

- **Study subsets.** Each function reads only the studies it is given: for eligibility
  (k = 3 within those studies), the trait screen, the observed sets and the U6 pool.
  Studies are kept or dropped whole. A study outside the population is an error.
- **Held-out tests.**
  - `evaluate` applies no eligibility rule and never replaces a hypothesis.
  - A hypothesis with a condition that no held-out study reports in a stool decreased
    signature is returned with `test_status` `no_units`. The other statuses are
    `tested`, `all_sets_empty` and `pool_too_small`.
  - The involved conditions that have held-out units, in alphabetical order, are the
    layout of the shared draws.
- **Historical leads.**
  - With `always_included=False`, the training screen alone chooses the traits, so
    cellulose is not forced into a split's selection.
  - PD, HIV and CRC are never forced in: they enter a split's triples only when they are
    eligible in its training studies.
- **Seeds.**
  - `seed` is an integer, or a sequence of integers for numpy's `SeedSequence`. An
    integer `s` gives the same draws as `[s]`.
  - The shared draws use `default_rng(seed)`. A triple drawn on its own uses
    `default_rng([*seed, index])`, where `index` is its row among the call's triples; for
    `evaluate` these are its distinct triples, sorted.
  - For split `i`, the suggested seeds are `(20261022, i, 0)` for `select` and
    `(20261022, i, 1)` for `evaluate`. The caller assigns the splits, with its own seed.
- **Draws.** `permutations` defaults to the protocol's 100,000; a split run may pass
  fewer and say so.
