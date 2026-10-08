# Biological interpretation of annotation results

This guide applies to every contributor and reviewer, across analyses, tables,
figures, captions, methods, manuscripts, generated narratives, talks, and public
descriptions. It implements [#179](https://github.com/CultureBotAI/microbial-trait-profiling/issues/179).
Use it with the [claim audit](biological_interpretation_audit.md).

The observations are **annotations linked through KG-Microbe to taxa in curated
BugSigDB differential-abundance signatures**. They describe annotation
representation under specified mappings and analytical choices. They do not
directly measure patient-level trait abundance, pathway activity, oxygen
exposure, or causation. Fixing a denominator, sampling unit, or statistical test
does not supply those biological measurements.

The code-level account below and
[Jed's paired contrast interpretation notes at `a2b2df98`](https://github.com/CultureBotAI/microbial-trait-profiling/blob/a2b2df98e9beec98b772c9cb5af7aa090d033a98/trait-visualization/trait-profile-investigation/paired_contrast_interpretation_notes.md)
provide the starting evidence: standardized annotation signals differ from raw
abundance, curated direction differs from z-score sign, and the arms share data.
This guide makes background enrichment distinct from an above-average z-score
and applies these limits to all outputs, beyond the paired-contrast analysis.

## Evidence chain and measurement units

Source study -> curated signature -> taxon mapping -> annotation/rollup -> statistic.

| Stage | What enters the analysis | Required traceability and limits |
| --- | --- | --- |
| Source study | Published comparison and reported differentially abundant taxa | Publication identifier, study ID, actual groups, body site, host, intervention/time point, and source measurement methods. A condition label alone does not establish disease versus healthy control. |
| Curated signature | Taxa labeled increased or decreased for a study/experiment | BugSigDB study, experiment, and original signature IDs; group definitions and orientation; export version, timestamp and checksum. `inc`/`dec` preserve the curated direction, not a verified disease/control contrast. Multiple signatures can aggregate into one profile. |
| Taxon mapping | Reported taxon IDs matched to NCBITaxon/KG nodes | Original ID and rank, resolved ID, unmatched/ambiguous mappings, and taxonomy snapshot. A higher-rank taxon is not a species-level measurement. |
| Annotation/rollup | Trait keys attached to mapped taxa and their descendants | KG release/checksums, source assertion and predicate/object IDs, full qualified column key, ontology version, direct versus rolled-up origin, deduplication rule, and source evidence where available. Display labels are not source evidence. |
| Statistic | Counts, fractions, paired contrasts, enrichment or study-level direction summaries | Input hashes, code revision/recipe, filters, denominator, reference/background, transformation, inference unit, support, effect, uncertainty, test and correction family. Retain machine-readable backing tables and links to upstream IDs. |

In [`taxa_trait_matrix.py`](../src/bsdbng/taxa_trait_matrix.py), `read_kg_edges`
stores a set of trait keys per KG taxon; repeated identical keys on that taxon
are deduplicated. `compute_trait_counts` counts trait-bearing taxa at or below
each reported taxon. In
[`study_trait_profiles.py`](../src/bsdbng/study_trait_profiles.py),
`collect_profile_taxa` deduplicates reported taxon IDs within a
study/experiment/direction by default, and `aggregate_profile_counts` sums their
vectors. This does not deduplicate overlapping descendant contributions from
different reported taxa, make traits mutually exclusive, or weight counts by
patient-level abundance. A count is therefore a **rolled-up annotation count**,
not a count of patients, cells, unique signature taxa, or independent assertions.

| Reported quantity | Units and interpretation |
| --- | --- |
| Profile trait count | Sum of mapped annotation counts for one study/experiment/direction; record rollup and source multiplicity. |
| Trait fraction `c[t] / sum(c[retained traits])` | Dimensionless annotation share; multiplying by 100 gives percentage of retained annotation counts, not percentage of taxa or patients. Declare all retained columns, including how MediaDive and grouped traits are handled. A zero denominator needs an explicit policy, not an inferred biological zero. |
| Increased-minus-decreased fraction | Difference in annotation shares between the two curated directions; percentage points only if the fractions were multiplied by 100. |
| Transformed composition contrast | Difference on the specified scale, for example `log1p(1000 * fraction)` followed by per-trait z-scoring. A value of 0.3 on that scale is not a 30% biological increase. |
| Enrichment/depletion versus background | Over-/under-representation of annotation counts against the declared urn/background and its inclusion rules; report observed/expected counts and fold, as well as p/q. It is distinct from the curated direction. |
| Study sign share | Fraction of non-tied studies with a positive annotation contrast; report positive, negative, tied and excluded studies and an interval. It does not measure the magnitude of patient-level change. |
| Correlation or agreement across analyses | Association or robustness to named analytical choices on shared evidence. Shared taxa, denominators, studies and trait selection create dependence. It is not independent replication. |

The aggregate profile file does not retain every original assertion or signature
ID. Keep upstream artifacts and a mapping/contribution sidecar when those are
needed for a claim. If provenance cannot be recovered, mark it unresolved and
limit the claim; do not manufacture gene-level evidence or replace the pinned KG
with a convenient newer snapshot. Source semantics and contribution reconciliation
are tracked in #189/#190; input coverage is tracked in #184.

## Four different meanings of direction

1. **Signature direction:** a source comparison labels a taxon increased or
   decreased. Human host plus paired signatures does not verify disease/control
   eligibility or orientation; use the manifest required by #188 before assigning
   disease-enriched/depleted labels.
2. **Relative annotation contrast:** a trait has a larger annotation share in one
   signature direction than the other, using an explicit denominator and scale.
3. **Enrichment against a background:** a profile's trait annotations are more
   represented than expected under a specified reference and test.
4. **Direct biological change:** quantitative observations in samples establish a
   change in an organism, trait, function or exposure. This requires additional
   measurements and a suitable design.

For example, a decreased signature with annotation share 0.20 can be enriched
against a background share of 0.05, even if its paired increased signature has
share 0.30. The decreased signature is enriched; the paired relative contrast is
positive. Neither establishes that a trait became more abundant in patients.
These are illustrative values, not a repository result or significance test.

## Wording and the evidence needed to go further

| Observed result | Defensible wording | Additional evidence for the stronger claim |
| --- | --- | --- |
| Higher annotation share in an increased signature | Greater representation of trait annotations among taxa in that curated increased signature | Quantitative sample-level trait measurements, with eligible and oriented comparisons, to claim increased trait abundance in patients |
| Taxa map to a pathway term | The signature contains taxa with pathway-associated annotations | Source-specific genomic evidence for capacity; expression, metabolites, flux or other functional measurements for pathway activity |
| Shift in oxygen-tolerance annotations | A shift in representation of annotated oxygen-tolerance traits | Anatomically matched environmental/physiological measurements for oxygen exposure, and an appropriate design to identify oxygen as a driver |
| No mapped trait count | No recorded annotation signal under the selected sources and mappings | Coverage assessment, detection limits and explicit negative evidence for biological absence |
| Multiple analyses agree | Robustness to the specified choices on shared annotation evidence | Non-overlapping data/cohorts and appropriate validation for replication; a suitable causal design for causation |
| Significant p/q-value | Evidence against the specified statistical null under its assumptions | Effect magnitude, uncertainty, support, calibration, population scope and independent evidence for biological importance or mechanism |

Do not mechanically ban words such as abundance or respiration: source studies
may measure them. Attribute each such claim to the actual source, population,
measurement and design, and distinguish it from this repository's result.

## Minimum result and caption contract

Every result-bearing statement must carry these details in its text/caption or
an explicitly linked table. A short public description must retain the annotation
unit and scoped claim and link the full evidence record.

- **Identity and scope:** qualified trait IDs, study/experiment/signature links,
  exact condition and anatomical stratum, host, contrast eligibility/orientation,
  vintage and relevant corpus restrictions. Report pooled `Other` constituents,
  exclusions, pairing, zero-signal/zero-denominator policies and coverage gaps.
- **Effect and uncertainty:** effect value, units, interval and interval method;
  p/q, null, thresholds and correction family where inference is made. If no
  interval was computed, say so and keep the summary descriptive. Never invent
  uncertainty or treat a significance marker as an effect size.
- **Support and inference level:** unique studies, matched experiments/profiles,
  mapped/unmapped taxa or annotation support as applicable, ties/exclusions, and
  the unit the test treats as an observation. Assertion counts are not independent
  biological replicates; multiple experiments or cells may share studies.
- **Comparison basis:** full retained denominator, transformation, background,
  population and selection rule for each arm. Matching export dates alone does
  not establish comparability. Keep the shared comparison basis and pairwise ratio
  sensitivity results across all 15 unique pairs among the six traits (as evaluated
  in the repository review for PR #192); report native effects and error-rate
  families separately. A selected six-trait set does not gain a new error-rate
  guarantee when displayed in a combined figure.
- **Evidence limit:** distinguish association from mechanism, and record unresolved
  calibration, semantic and coverage limitations. Significance alone does not
  establish biological importance, generality, replication or a mechanism.

Template for a result or generated narrative (replace every bracket; use
`not estimated` or `unresolved` where appropriate):

> In [eligible corpus, condition/site, vintage], [trait ID and label] had
> [effect, units, interval/method] for [precisely defined annotation contrast or
> background comparison], supported by [studies, pairs/profiles, annotation or
> mapping support; ties/exclusions]. The analysis used [denominator, transform,
> background, inference unit, p/q and correction family]. This describes
> annotation representation in curated signatures; [specific unresolved limit].
> Backing result and provenance: [artifact/row and source IDs].

Heatmap captions must identify fills as annotation shares and provide support and
cell-status tables. A zero signal, an untested/low-support cell and an unsampled
cell are different states even if all render blank. Machine-readable status tables
use `no_annotation_signal` for sampled zero-signal cells across guild and individual
trait outputs. State when the heatmap has no cell-level uncertainty. In PR #173,
the proposed enrichment outline summarizes
profile-level calls across supporting studies; this aggregation has no calibrated
cell-level p-value. The implemented mark kinds and verdicts must match the legend
(the audit records a discrepancy, including composition-only marks).
A body-site triangle summarizes a pooled site-level composition contrast and
does not establish significance for each condition in that column. The
[audit](biological_interpretation_audit.md#pr-173-proposed-overlay-and-narrative)
provides replacement wording and remaining limitations.

## Worked example: the Rosalind transcript and repository review

The [PR #192 archive](https://github.com/CultureBotAI/microbial-trait-profiling/pull/192)
preserves a user-supplied, hypothesis-generating conversation separately from a
repository review. Its original date, model version and conversation URL were
not supplied. Preserve the archival transcript; use the review, not the chat's
conclusion, when writing results. This guide adopts the following distinctions
from that supplied review; it does not claim a new literature review or reanalysis.

| Archived interpretation | Correction and next evidence |
| --- | --- |
| "probably causal ecological readouts" | Unsupported by these annotation data. The contrasts motivate redox, substrate and taxonomic-turnover hypotheses requiring source-specific and independent validation. |
| A unifying gut-redox mechanism for the six directions | Candidate mechanism in selected intestinal settings. Require eligible disease/control comparisons, exact anatomical scope, retained counts, effects and uncertainty. A pooled multi-site result does not establish the mechanism in oral, skin, airway or other sites. |
| Veillonella as "facultative-adjacent" support for oxygen use | The cited Rojas-Tapias example concerns anaerobic nitrate respiration, not evidence of facultative aerobiosis or oxygen use for this repository's annotated taxa. See the primary paper linked in the PR review. |
| Merge coccus and sphere into an opportunist signal | Preserve METPO:1000668 and METPO:1000683 separately until ontology relations, source assertions and empirical overlap justify grouping. Shared labels or correlations do not prove synonyms. Any union needs deduplication and a newly declared denominator/testing family (#189). |
| Arsenate reduction means an arsC detoxification/plasmid passenger | `METPO:2000103\|pathways:arsenate_reduction` does not identify arsC or arrAB, plasmids or co-selection. Trace the assertion; classify detoxification, respiratory, mixed or unresolved support (#189). |
| Moderately halophilic means halotolerant/osmotolerant | For METPO:1000623 recover source salt requirement, optimum/range, units and culture conditions before translating the term (#189). Neither local salt exposure nor contributing genera follow from the label. |
| Cellulose signal indicates SCFA loss or impaired cross-feeding | `METPO:2000103\|GO:0030245` does not measure complete cellulose degradation, butyrate, expression or flux. Inventory actual gene, expression and metabolite observations (#191). |
| Named genera should drive the signal; exclusions should weaken it | Falsifiable predictions, not expected findings. Reconcile all trait contributions, allow null/opposite outcomes and recompute denominator/background/eligibility under declared exclusions. Exclusions diagnose influence, not a causal intervention (#190). |

A defensible mechanistic paragraph is:

> These annotation contrasts motivate a hypothesis that redox or substrate
> availability contributes to taxonomic turnover in selected intestinal settings.
> Testing it requires verified contrast orientation, source-resolved trait and
> taxon attribution, and anatomically matched measurements of respiratory
> resources, inflammation and substrate/metabolite availability in suitable
> independent data. The present annotations do not establish exposure, pathway
> activity or causation; patterns in other anatomical sites require their own
> evaluation.

## Review responsibilities and open evidence work

Jed, Ilyas and Noel coordinate adoption and identify the authors of every affected
manuscript section, report, caption, slide deck and public description. They are
not the only reviewers responsible for this guide. Each output owner updates the
claim audit with an exact location, old claim, correction or evidence, revision,
owner and disposition. Correct generators as well as their rendered output;
label historical artifacts and keep their provenance. An unreviewed external
output is not implicitly approved by this guide.

| Follow-up | Deliverable and lead | Boundary |
| --- | --- | --- |
| #188 | Contrast eligibility/orientation manifest; Noel with Jed | Required before disease/control language; host filtering alone is insufficient |
| #189 | Trait semantics, source evidence and justified grouping; Ilyas with Jed | Required before label translations, gene claims or merging features |
| #190 | Reconciled contributions and genus/family influence; Jed with Noel and Ilyas | Draft PR #194 supplies provisional contribution/influence artifacts; orientation, semantics, native Jed/ORA integration and calibration remain unresolved, and sensitivity is not causal validation |
| #191 | Inventory of measured inflammation/SCFA evidence and feasible independent validation; Noel with Ilyas | Starts with feasibility and protocol; does not assume measured markers or successful validation |

These deliverables feed #178 and #179; none marks mechanisms validated. Preserve
the existing scope of #118/#120/#176 (scope/strata), #182/#183 (calibration and
sensitivity), #184 (input coverage), #131/#132 (annotation sampling/background
limitations), #165 (dependence), and #172/#173 (integrated display). Statistical
implementation changes and biological evidence are separate review questions.

