# The integrated evidence table

Produced by `src/bsdbng/integrated_evidence.py`, driven by `just integrated-evidence`,
which writes `data/derived/integrated_evidence/`: the machine-readable table
(`integrated_evidence.tsv`), a generated view with the data dictionary
(`integrated_evidence.md`) and a provenance record. It implements
[#178](https://github.com/CultureBotAI/microbial-trait-profiling/issues/178) and follows
the shared comparison basis in [`three_arm_statistics.md`](three_arm_statistics.md).

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
  `trait`, `trait_qualifier` (the raw identifiers), `reference` for a comparison that
  names one, `stratum_kind` and `stratum`, and `source` (the input file and a hash
  prefix).
- **Population**: `population`, `eligibility_rule`, `denominator`, `background`, and
  the counts that are never interchangeable - `unique_studies`,
  `paired_experiments`, `profiles`, `testable_units`, `ties`.
- **Effect**: `effect_value` with its `effect_units`, and `interval_low`,
  `interval_high` and `interval_type`. A standardized mean contrast, an ORA fold and
  a directional proportion are different quantities and are never put in one column.
  Where no interval was computed, `interval_type` is `unavailable` rather than blank.
- **Inference**: `p_value`, `q_value`, `correction_method`, `family_id`,
  `family_size`, `q_threshold`, `effect_threshold`. A decision holds only inside its
  own family.
- **States**, kept apart because they overlap: `test_status` (`tested`, `untested`,
  `insufficient_support`, `all_tied`, `zero_annotation_signal`, `not_sampled`,
  `source_not_supplied`), `decision` (`significant`, `not_significant`,
  `undecided`), `display_status` (`drawn`, `masked_below_min_studies`, `not_drawn`,
  `not_applicable`) and `zero_annotation_state`.
- **Standing**: `selection_note`, and one status field per open assessment -
  `calibration_status` (#182), `sensitivity_status` (#183), `coverage_status` (#184),
  `orientation_status` (#188) and `taxon_influence_status` (#190). Each reads from
  that analysis's own output when it exists and is `pending` otherwise.

Two rules are enforced in code rather than left to the reader. A row that was not
tested has its effect, interval, p and q blanked, because a source may carry a
placeholder - the fraction analysis reports share 0 and a q derived from p = 1 for an
all-tied trait before #165 - and a placeholder read as an effect is the confusion this
table exists to prevent. And a method whose output was not supplied gets a row per
feature with `source_not_supplied`, so that its absence is visible rather than being
read as a negative result.

## What is in it now

On the pinned v1.3.1 corpus, 40 records:

| Method | Records | Family | Tested | Significant at q <= 0.05 |
| --- | --- | --- | --- | --- |
| directional consistency, fraction sign test | 6 | BH over six traits | 6 | 6 |
| directional consistency, log-ratio sign test | 15 | BY over 15 pairs | 15 | 3 |
| oxygen guild McNemar | 1 | single test | 1 | 1 |
| composition arm, native transformed contrast | 6 | - | 0 | not supplied |
| enrichment arm, profile hypergeometric | 6 | - | 0 | not supplied |
| aggregated profile support | 6 | - | 0 | not supplied |

The three unsupplied methods are the point of the design rather than a gap in it. The
composition arm's native output moves under #174/#197, the enrichment arm's rests on
the urn of #131, and the aggregated support rule lives on the #173 branch. Each is
recorded as not supplied with the reason, and pointing the recipe at those outputs
adds them without changing any other row.

Six of the six fraction results and three of the fifteen ratio results are
significant in their own families. Those are not nine findings: the ratios are
measured on their own eligible populations, the families are corrected separately,
and the two ask different questions. The table reports them side by side precisely so
that a reader cannot silently add them up.

## Standing of the results in it

Every row carries the same caveats as fields rather than as prose:

- **Selection.** The six Table 1 traits were chosen on this corpus by the composition
  arm. #182 measures what selecting on the data and testing the selection costs:
  per-test error roughly 0.08 against a nominal 0.05, a lower bound. These results
  are exploratory.
- **Calibration.** Where #182 has run, `calibration_status` carries its verdict per
  procedure. The procedures this repository owns hold their nominal rates under the
  paired null; the other arms' procedures are not assessed.
- **Population.** #183 shows most fraction directions do not survive restriction to a
  ratio's eligible population, and #188 gives the disease-versus-healthy subset in
  which the directions do hold. Both are populations, not corrections.
- **Coverage and influence.** #184 and #190 remain pending; until they land, no row
  should be read as accounting for annotation coverage or taxonomic concentration.

## Reproducing

`just integrated-evidence` on the pinned inputs. It reads
`data/derived/calibration/calibration_status.json`,
`data/derived/table1_composition_sensitivity/sensitivity_status.json` and
`data/derived/contrast_orientation/contrast_orientation_provenance.json` when they
exist, and `--coverage-status` and `--influence-status` take the #184 and #190
outputs once those exist. The provenance record carries the code revision, input
hashes, the corpus cascade, the assessment statuses and the per-method summary.

`tests/test_integrated_evidence.py` covers the states the contract exists to
preserve: an all-tied trait, a pair with no positive counts, a method that was not
supplied, two methods keeping separate families and units, display state independent
of decision, statuses attached without changing a result, an unreadable status file,
and a full export in which counts reconcile with the source runs.

The table, the contract and this description were drafted with the assistance of a
coding agent (Claude Code) and require confirmation by the method owners before any
part of it goes into the paper.
