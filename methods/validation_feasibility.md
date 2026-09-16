# Measured evidence for independent validation: what exists, and what a test would need

Produced by `src/bsdbng/validation_evidence.py`, driven by `just validation-evidence`,
which writes `data/derived/validation_evidence/`: the per-study evidence-availability
matrix, its summary, the marker-stratified experiments, the inferred-function traits,
the external cohort overlap, and a provenance record. It implements
[#191](https://github.com/CultureBotAI/microbial-trait-profiling/issues/191) and
acquires no new data.

## The question, and the short answer

The Rosalind discussion proposed testing the six Table 1 trait signals against
measured inflammation and SCFA data. The profile matrix cannot supply that: its
cells are counts of trait annotations inferred from taxa in curated signatures.
This note asks what measured evidence is reachable, and concludes that **no
correlation between a trait signal and a measured marker is feasible on the current
data**, for a reason that is structural rather than a matter of effort: BugSigDB
records, per experiment, a list of taxa that differed between two groups. It holds
no per-sample abundances. A measured calprotectin or butyrate value is per sample,
so there is nothing in this corpus to correlate it against, whatever a source paper
reports. Any real test needs sample-level abundances and sample-level measurements
from the same individuals, which means reprocessing raw sequence data and obtaining
matched clinical measurements - a separately scoped project.

## What the curated metadata mentions

Searching the export's curated text - condition, both group names, the Group 1
definition, description, matching and confounder fields, title and keywords - over
1,080 human studies, 959 of them in the discovery corpus:

| Evidence | Discovery studies mentioning | Most frequent words |
| --- | --- | --- |
| Clinical covariates (age, sex, BMI, smoking, medication) | 480 | age, sex, bmi |
| SCFA or metabolomics | 94 | short-chain fatty acids, metabolome, butyrate |
| Diet, fibre or salt | 78 | diet, dietary, fibre |
| Disease activity or severity | 41 | severity, disease activity, UPDRS |
| Actual antibiotic exposure | 33 | antibiotic exposure, antibiotic treatment |
| Inflammation markers | 28 | cytokines, calprotectin (4 studies), C-reactive protein (3) |
| Transit or stool consistency | 11 | Bristol stool, transit |
| Nitrate, oxygen or redox | 6 | oxidative stress, oxygen |

A mention is a pointer to a source paper, not data. The export holds no values, no
units, no time points and no sample linkage for any of them, and "not mentioned"
does not mean "not measured": the fields are curated summaries, not data
dictionaries. Two fields are often mistaken for exposure and are reported
separately: `Antibiotics exclusion` is an eligibility rule, not a record of who took
antibiotics, and a disease label is not an inflammation measurement.

147 of the 959 discovery studies used whole-metagenome sequencing, so gene or
pathway content could in principle be measured from their raw reads. That is a
different quantity from a metabolite concentration, and it needs the raw data.

## Evidence that is inside the corpus, and therefore not independent

**Marker-stratified experiments.** 26 experiments across the export have a condition
that is itself a measurement; 11 are human and 4 are in the discovery corpus: a
cytokine comparison in type 2 diabetes, an interleukin-4 measurement comparison, and
two oxidative-stress comparisons (baseline against placebo and against probiotic).
These do link a curated signature to a measured marker, and they are the only such
link available without new data - but they are part of the same corpus the signals
were derived from, they are four experiments, and two of them are intervention
arms. They cannot serve as independent validation, and they are far too few for a
correlation.

**Inferred-function traits.** Six retained trait columns concern SCFA or
fermentation, and what they assert matters more than that they exist:

| Column | Trait | What the annotation asserts |
| --- | --- | --- |
| `METPO:2000006\|CHEBI:30089` | acetate | uses as carbon source |
| `METPO:2000006\|CHEBI:17272` | propionate | uses as carbon source |
| `METPO:2000006\|CHEBI:24996` | lactate | uses as carbon source |
| `METPO:2000006\|CHEBI:26806` | succinate | uses as carbon source |
| `METPO:2000103\|pathways:propionate_oxidation` | propionate oxidation | capable of |
| `METPO:1002005` | Fermentation | (unqualified) |

Every SCFA-specific column in the matrix is **consumption**, not production: taxa
annotated as using acetate or propionate as a carbon source, or as oxidising
propionate. A hypothesis about SCFA *production* cannot be tested even indirectly
with these columns, and reading a higher "propionate" annotation share as more
propionate produced would invert the assertion. The one production-adjacent column
is the unqualified Fermentation trait. This also bears on the trait-semantics work
in #189.

## External cohorts

The issue names [Sun et al. (2024)](https://doi.org/10.1038/s41522-024-00567-9), a
reanalysis of 6,314 fecal metagenomes from 36 Chinese case-control studies plus two
cohorts the authors recruited themselves. Its sample table and Supplementary Table 1
are recorded in [`external_cohorts_sun2024.tsv`](external_cohorts_sun2024.tsv),
regenerated from pinned sources by
`just validation-evidence -- --snapshot-from <sample.group> --snapshot-supplement <table.txt>`
(the sample table is pinned by commit and SHA-256 in the module; Supplementary
Table 1 is the PDF's text layer).

Matching the source studies' publication titles against BugSigDB titles:

- **7 of the 36 overlap the discovery corpus** (pulmonary tuberculosis, Parkinson's
  disease, end-stage renal disease, young-onset colorectal cancer, COVID-19,
  schizophrenia, breast cancer).
- 1 more is curated in BugSigDB but outside the analysed corpus.
- 28 have no title match, which is not proof of independence: a study may be
  curated under a different title, and later papers can reuse the same participants.
- 2 are the authors' own cohorts, which have no prior publication title and cannot
  be assessed this way.

More decisive than the overlap: the collection's sample metadata carries age, sex,
BMI, province and sequencing platform, and **no inflammation, SCFA, nitrate or
oxygen measurement at all**. It could support an independent test of taxon-level or
annotation-level contrasts on its non-overlapping studies, after those are
confirmed, but it cannot test a redox or SCFA mechanism.

## The feasibility gap, stated as a data requirement

A test of the redox or SCFA interpretation needs, for one cohort:

1. Per-sample microbiome abundances, not a differential-abundance taxon list.
2. Per-sample measurements from the same individuals: fecal calprotectin or another
   inflammation marker; SCFA concentrations by chromatography; ideally nitrate or a
   redox proxy, transit, and diet.
3. Enough samples for the estimand below, with disease and control eligibility and
   orientation verified as in #188.
4. No overlap with the discovery studies, or an explicit statement that the test is
   internal rather than independent.

None of the sources examined here supplies 1 and 2 together. This is the deliverable
the issue asks for: an explicit gap, not an invented correlation. Public cohorts
that plausibly hold both (inflammatory-bowel-disease cohorts with paired
metagenomes and faecal calprotectin, or intervention trials with SCFA panels) exist
in the literature, but no such cohort has been assessed here, and none should be
named as available until its data have been checked.

## Protocol, to be frozen before any such data are examined

If a cohort meeting the requirement is obtained, the analysis is declared first:

- **Estimand.** The association between a sample's measured marker (calprotectin, or
  an SCFA concentration) and that sample's annotation share for the trait in
  question, computed with the same denominator as the Table 1 analyses (the summed
  non-MediaDive trait counts of that sample's taxa).
- **Observation unit.** One sample. Repeated samples from a subject, and subjects
  within a cohort, are nested: a mixed model with subject and cohort random effects,
  or one value per subject with the aggregation stated.
- **Compositional treatment.** Annotation shares are compositional, so the primary
  analysis uses a log-ratio against a declared reference trait, with the
  fraction-scale result reported beside it, as in the log-ratio sensitivity work.
- **Uncertainty and multiplicity.** Effect estimates with intervals; the family is
  the six traits by the declared markers, corrected once, declared before testing.
- **Independence.** Cohort overlap with the discovery studies checked and reported
  first. If any overlap exists, the result is labelled internal consistency, not
  replication. Hypotheses, trait definitions and the reference trait are frozen
  before outcomes are examined; a retrospective split of the discovery corpus is not
  independent validation.
- **Outcomes that would weaken the hypothesis**, declared in advance: no association
  between the oxygen-tolerance annotation share and measured inflammation; an
  association of the opposite sign; an association that disappears when transit,
  diet or antibiotic exposure is adjusted for; or a cellulose-annotation association
  with measured SCFA that does not follow the annotation's direction of assertion
  (consumption, per the table above). An inconclusive estimate - a wide interval
  including zero - is reported as inconclusive, not as evidence against.
- **What a positive result would not show.** An observational association does not
  establish that oxygen, inflammation or SCFA production drives the annotation
  contrasts, and annotation shares remain annotations, not patient measurements
  (#179).

## Reproducing

`just validation-evidence` on the pinned BugSigDBExports v1.3.1 inputs. The
provenance record gives the export, profile and cohort-table SHA-256 values, the
pinned external source and its checksum, and states that mentions are pointers
rather than data. The KG-Microbe node table is needed only for the trait labels; the
run records whether it was available.

The inventory, the protocol and this description were drafted with the assistance of
a coding agent (Claude Code) and require review by the author before manuscript use.
