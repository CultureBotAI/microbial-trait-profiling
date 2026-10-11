# Measured evidence for independent validation: what exists, and what a test would need

Produced by `src/bsdbng/validation_evidence.py`, driven by `just validation-evidence`,
which writes `data/derived/validation_evidence/`: the per-study evidence-availability
matrix, its summary, the marker-split candidates, the inferred-function traits,
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
1,080 human studies, 959 of them in the discovery corpus. Every figure below is
read back from the artifacts by
`just validation-evidence -- --report`; none is maintained by hand.

| Evidence | Discovery studies mentioning | Human studies mentioning | Most frequent words |
| --- | --- | --- | --- |
| Clinical covariates (age, sex, BMI, smoking, medication) | 481 | 542 | age (418), sex (276), body mass index (163) |
| Whole-metagenome sequencing | 147 | 167 | - |
| Diet, fibre or salt | 79 | 87 | diet (62), dietary (26), fiber (15) |
| SCFA named as a metabolite | 47 | 56 | short-chain fatty acids (22), butyrate (14) |
| Untargeted metabolomics | 45 | 48 | metabolome (21), metabolomics (20) |
| Disease activity or severity | 41 | 45 | severity (29), disease activity (5), UPDRS (4) |
| Actual antibiotic exposure | 33 | 34 | antibiotic exposure (21), antibiotic treatment (5) |
| Inflammation markers | 28 | 32 | cytokines (9), cytokine (6), calprotectin (4) |
| Expression, proteome or flux | 16 | 17 | gene expression (5), metatranscriptome (3) |
| Transit or stool consistency | 11 | 11 | transit (6), Bristol stool (4) |
| Nitrate, oxygen or redox | 6 | 7 | oxidative stress (3), oxygen (3) |
| SCFA named inside a gene or enzyme name | 1 | 1 | BCoAT, Butyryl-CoA:acetate |

A mention is a pointer to a source paper, not data. The export holds no values, no
units, no time points and no sample linkage for any of them, and "not mentioned"
does not mean "not measured": the fields are curated summaries, not data
dictionaries. Each row records which fields were searched for that study, so an
absent mention can be read against the text that existed to search.

**Three things that look like one row and are not.** An SCFA named as a metabolite
(47 discovery studies), untargeted metabolomics (45) and an SCFA inside an enzyme
name (1) are counted separately, because they are different evidence. A single
lumped pattern previously reported 94 discovery studies, of which 44 matched only
`metabolom*` - a metabolome survey is not an SCFA panel. A further three studies
matched on a name that is not a metabolite at all, and are excluded by span and
recorded in `scfa_metabolite_excluded_terms`: the condition
"Medroxyprogesterone acetate" (`bsdb:30866773`) and the title phrase "fluticasone
propionate cream" (`bsdb:27543211`) are steroid esters, and
"BCoAT: Butyryl-CoA:acetate CoAtransferase" (`bsdb:24296500`) is an enzyme name,
which is why genes and metabolites are kept distinct (#191 AC2).

Two further fields are often mistaken for exposure and are reported separately:
`Antibiotics exclusion` is an eligibility rule, not a record of who took
antibiotics, and a disease label is not an inflammation measurement.

147 of the 959 discovery studies used whole-metagenome sequencing, so gene or
pathway content could in principle be measured from their raw reads; 16 mention
transcriptomics, proteomics or flux, which is the expression and flux inventory
#191 AC2 asks for. Both are pointers to source papers, not datasets, and both need
the raw data. Neither is a metabolite concentration.

Whether the source articles' own supplementary data are obtainable has not been
assessed, and no per-study access value is emitted: a column saying "not assessed"
in every row would be a sentinel standing in for work nobody has done.

## Evidence that is inside the corpus, and therefore not independent

**Marker-split candidates.** 75 experiments across the export name a marker of one
of these evidence types in their condition or group fields. Whether the groups were
actually cut on that marker's measured value is decided from the group names and
the Group 1 definition, not from the condition label, which is unreliable in both
directions: `bsdb:37180439-9` is labelled "C-reactive protein measurement" but
compares healthy controls with liver-disease patients, while `bsdb:35888992-6` is
labelled "Increased inflammatory response" and names no marker in its condition at
all, though its groups are inflammatory bowel disease patients either side of fecal
calprotectin 1,000 mg/kg. Candidates that fail are kept in the table with the
reason rather than dropped.

| Classification | Experiments | Human | In discovery corpus |
| --- | --- | --- | --- |
| `marker_stratified` | 31 | 24 | 11 |
| `marker_named_groups_not_measured` | 29 | 18 | 15 |
| `intervention_contrast` | 15 | 4 | 2 |

The 24 human measured splits come from nine studies: a multiplex cytokine and
growth-factor panel in schizophrenia (`bsdb:39103909`, 13 experiments, five in the
discovery corpus), fecal calprotectin in Parkinson's disease (`bsdb:31815177`,
three), the inflammatory-bowel-disease calprotectin threshold (`bsdb:35888992-6`),
genital inflammation defined by measured cytokine levels (`bsdb:33220709-1`),
systemic inflammation in children (`bsdb:38890302-6`), C-reactive protein in
obesity (`bsdb:23526699-4`) and in ageing (`bsdb:39237540-9`), interleukin-6 and
malondialdehyde in bipolar disorder (`bsdb:30051546`, two), and the corpus's only
measured SCFA split, low against high fecal butyrate before haematopoietic stem
cell transplantation (`bsdb:31065565-2`).

Excluded, with the reason recorded: `bsdb:879-4` ("Oxygen") compares patients on
oxygen therapy with patients not on it, which is an intervention;
`bsdb:39748236-1/2/4` ("Oxidative Stress") compare a baseline against a placebo or
probiotic arm; `bsdb:36687728-2` is conditioned on "Cytokine" but its groups are
prediabetes against type 2 diabetes.

These do link curated signatures to measured markers, and they are the only such
link available without new data - but they are part of the same corpus the signals
were derived from, 11 of them are in the discovery corpus at all, and 22 of the 24
measure inflammation rather than SCFAs or redox. They cannot serve as independent
validation, and they are far too few for a correlation. This is a feasibility
inventory of what exists, not a patient-level correlation analysis; no correlation
is computed anywhere in this work, because none can be.

**Inferred-function traits.** Six retained trait columns concern SCFA or
fermentation, and what they assert matters more than that they exist. The labels
come from the KG-Microbe node table, whose MD5 is verified against
`data/kg-microbe/provenance.json` before they are read; the run records that
verdict, and refuses to label anything from an unidentified node file.

| Column | Trait | What the annotation asserts |
| --- | --- | --- |
| `METPO:2000006\|CHEBI:30089` | acetate | uses as carbon source |
| `METPO:2000006\|CHEBI:17272` | propionate | uses as carbon source |
| `METPO:2000006\|CHEBI:24996` | lactate | uses as carbon source |
| `METPO:2000006\|CHEBI:26806` | succinate | uses as carbon source |
| `METPO:2000103\|pathways:propionate_oxidation` | propionate_oxidation | capable of |
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
reanalysis of Chinese fecal metagenomes. Three sample counts are in play and it is
worth keeping them apart: the paper's abstract says it reanalysed **6,314**
metagenomes from 36 case-control studies, while the sample table the authors
published lists **6,616** samples in all - **6,463** across those 36 studies plus
**153** in the two cohorts they recruited themselves. The 6,314 is the paper's own
figure and 149 short of its own table; the difference is unexplained in the
supplement and is presumably post-QC. The counts in
[`external_cohorts_sun2024.tsv`](external_cohorts_sun2024.tsv) are the table's.

That snapshot is regenerated from two pinned inputs by

```
just validation-evidence -- --snapshot-from <sample.group> --snapshot-supplement <table.txt>
```

- `sample.group`, pinned by commit and SHA-256 in the module;
- Supplementary Table 1, which is a PDF. A PDF text layer is produced by the
  extractor rather than stored in the file, so the input's identity is the PDF plus
  the command: `41522_2024_567_MOESM1_ESM.pdf` (1,865,194 bytes, SHA-256
  `67fa73e0...`) from the Europe PMC supplementary-files zip, extracted with
  `pdftotext -layout` under poppler 26.02.0, whose output SHA-256 is pinned too.

Both are checked before anything is read, and the regenerated table is checked
against the committed snapshot's SHA-256 afterwards, so these exact bytes are known
to produce exactly the table in version control.

Overlap with the discovery corpus is decided by **normalised publication-title
similarity at 0.9 and nothing else** - case, spacing and punctuation removed, which
absorbs PDF artefacts such as "Ametagenome" for "A metagenome". Author surnames,
year and location are *not* used: a title either names the same paper or it does
not, and guessing from a surname would manufacture overlap. A source study with no
publication title cannot be assessed this way and is reported as such.

| Overlap status | Source studies |
| --- | --- |
| `no_title_match_found` | 28 |
| `overlaps_discovery_corpus` | 7 |
| `not_assessable_without_title` | 2 |
| `curated_in_bugsigdb_outside_discovery_corpus` | 1 |

The 7 overlapping studies are pulmonary tuberculosis, Parkinson's disease,
end-stage renal disease, young-onset colorectal cancer, COVID-19, schizophrenia and
breast cancer. The 28 without a title match are not thereby independent: a study
may be curated under a different title, and later papers can reuse the same
participants. The 2 unassessable ones are the authors' own cohorts, which have no
prior publication.

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

This is the pre-registration #191 AC7 asks for, and it is the protocol the
annotation-versus-measurement distinction in
[#179](https://github.com/CultureBotAI/microbial-trait-profiling/issues/179) and the
eligibility and orientation checks in
[#178](https://github.com/CultureBotAI/microbial-trait-profiling/issues/178) feed
into: #179 is why the estimand below is an annotation share rather than a patient
measurement, and #178 is why requirement 3 names eligibility and orientation
verification. If a cohort meeting the requirement is obtained, the analysis is
declared first:

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
  diet or antibiotic exposure is adjusted for; or an association between measured
  SCFA and one of the acetate or propionate "uses as carbon source" columns in the
  table above that runs opposite to what consumption implies. An earlier draft cited
  a cellulose row here; there is no cellulose column in that table - cellulose
  catabolism is one of the six Table 1 traits, not one of the SCFA or fermentation
  columns. An inconclusive estimate - a wide interval including zero - is reported
  as inconclusive, not as evidence against.
- **What a positive result would not show.** An observational association does not
  establish that oxygen, inflammation or SCFA production drives the annotation
  contrasts, and annotation shares remain annotations, not patient measurements
  (#179).

## Reproducing

```
just validation-evidence            # the five tables and the provenance record
just validation-evidence -- --report  # the tables above, as markdown, from those files
```

on the pinned BugSigDBExports v1.3.1 inputs. Every count in this note comes from
`--report` reading the artifacts back, so a stale number here is a diff rather than
a discrepancy nobody notices.

The provenance record gives the export's vintage and SHA-256, the profile SHA-256,
the cohort-table SHA-256, the pinned external sources with their checksums and the
supplement's extraction command, the per-classification marker-split counts, and
the node table's verification verdict. It also states that mentions are pointers
rather than data, and that source-article access was not assessed.

The KG-Microbe node table is needed only for the inferred-function labels, and its
MD5 is verified against `data/kg-microbe/provenance.json` before they are read.
Without it the run warns on stderr and writes that table with its header and no
rows - never a headerless empty file, which cannot be told from a truncated write -
and `--require-nodes` turns the warning into an error. Fetch it with
`just download-kg`.

**Scope.** This is a feasibility inventory: it counts what measured evidence exists
and concludes that a patient-level correlation between a trait signal and a measured
marker is not possible on this corpus. No such correlation is estimated here, and
none of the counts above should be read as one.

The inventory, the protocol and this description were drafted with the assistance of
a coding agent (Claude Code) and require review by the author before manuscript use.
