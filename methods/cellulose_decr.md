# Cross-condition overlap of stool-decreased cellulose-annotated taxa

This note describes how the tables and figures in
[`data/analysis/cellulose_decr/`](../data/analysis/cellulose_decr/README.md) are
produced (#552). The analysis asks two questions about taxa annotated
`METPO:2000103|GO:0030245` (cellulose catabolic process) in KG-Microbe:

1. Within Parkinson's disease (PD), HIV infection and colorectal cancer (CRC), do
   curated decreased signatures carry a larger cellulose annotation share than the
   paired increased signatures?
2. Do the three conditions report the same cellulose-annotated taxa decreased in
   feces more often than chance, and against which chance model?

Both are questions about curated signature annotations. Neither measures
cellulose degradation, fiber intake, short-chain fatty acids or any other
physiology (see [`biological_interpretation.md`](biological_interpretation.md)).

## Inputs and population

| Input | Use |
| --- | --- |
| `trait-visualization/study_trait_profile_study_metadata.tsv` | Defines the 5,829 retained human profiles (one per study, experiment and direction) and their non-MediaDive annotation totals. Its condition and body-site columns are **not** used by default (see below). |
| `data/derived/experiments.tsv` | Condition, body site and host of each experiment. It matches `data/raw/full_dump.csv` for all 4,435 experiments. |
| `data/derived/study_trait_profiles_wide.tsv` | Annotation counts per profile, for the cellulose annotation share. |
| `data/derived/bugsigdb_taxa.tsv` | Taxa reported in each signature. |
| `data/derived/taxa_trait_matrix.tsv` and `taxa_trait_matrix_columns.tsv` | Which of the 5,829 matrix taxa carry each trait, and the trait labels. |
| `data/derived/contrast_orientation/contrast_manifest.tsv` | Contrast class and orientation of every experiment ([`contrast_orientation.md`](contrast_orientation.md)). |
| `data/derived/studies.tsv`, `data/raw/full_dump.csv` | Study identifiers, authors, journals and locations for the evidence tables. |

The run's `cellulose_decr_provenance.json` records the SHA-256 of every input, the
code commit, the seeds and the package versions. Numbers in this note come from
BugSigDB vintage `2026-04-24_00:41_UTC`.

**Condition and body site come from each experiment's own record.** The profile
metadata table holds one modal condition and body site per *study*
(`scripts/create_study_trait_profile_heatmap.py`, `load_raw_study_metadata`), so
1,112 of its 5,829 profiles carry a label that belongs to another experiment of
the same study. The October 16 outputs used those labels; 61 experiments in the
three conditions changed when the analysis switched to `experiments.tsv`, among
them saliva and nasal experiments that had been counted as PD feces.
`label_source_corrections.tsv` lists them, and `--label-source study_modal`
reproduces the October 16 numbers.

A condition matches only its exact BugSigDB label, so multi-condition labels such
as `Constipation,Parkinson's disease` are excluded (listed in the summary JSON).
Feces is the exact body site `Feces`; other sites are grouped by keyword into
non-fecal GI (mucosa or lumen), oral or oropharyngeal, and other non-GI, and the
grouping actually applied to the three conditions is written to the summary JSON.

## Direction tests

For each profile the cellulose annotation share is the cellulose count divided
by the sum of all retained non-MediaDive counts. For each experiment with both
signatures, the paired difference is the share in the increased signature minus
the share in the decreased signature; *down* means the decreased signature had
the larger share. Ties are excluded.

- **Experiment level:** exact two-sided sign test and Wilcoxon signed-rank test
  over paired experiments. Experiments from one study are not independent, so
  these p-values are optimistic.
- **Study level:** one vote per study, the sign of the mean paired difference over
  its experiments; exact two-sided sign test, and a Wilcoxon test on the study
  means. No study contributes to more than one of the three conditions.
- **Populations:** all contrasts; disease-versus-healthy contrasts with reversed
  experiments relabelled so *decreased* means lower in disease
  (`disease_vs_healthy_oriented`); the same without experiments flagged
  `weak_comparison_evidence`; and, for HIV, the oriented population without
  experiments whose group names or descriptions mention MSM. That last filter is
  lexical; it does not adjust for sexual behaviour, which was not recorded for most
  studies.
- **Compartments:** for each compartment and condition the table reports the
  smallest two-sided p a sign test can reach with the observed untied pairs and the
  power to detect the pooled fecal share at alpha 0.05. Fisher's exact test of
  direction by compartment compares feces with each other compartment, at the
  experiment and the study level, with the conditional maximum-likelihood odds
  ratio and its exact 95% interval. Units are treated as independent, and a few
  studies sampled two compartments.

## Overlap tests

A condition's set holds the matrix taxa annotated with the trait that at least
one of its feces decreased signatures reports, at any rank. Overlap is the
Jaccard similarity J and distance d_J = 1 - J for each pair, their mean over the
three pairs, and the number of taxa in all three sets.

| Null | What is drawn | What it holds fixed |
| --- | --- | --- |
| U1 | Uniform draws from the 189 matrix taxa with the annotation | set sizes |
| U2 | Uniform draws from the annotated taxa ever reported in a feces signature | set sizes |
| U3 | Uniform draws from the annotated taxa ever reported decreased in feces | set sizes |
| U4 | A random trait placed on as many of the conditions' decreased feces taxa as the observed union | each condition's decreased feces taxa |
| U5 | Draws from the U2 universe weighted by the number of other feces decreased studies reporting the taxon, plus 0.5 | set sizes, reporting prevalence |
| U6 | Whole study signatures: 49, 19 and 31 distinct feces-decreased studies drawn from all 646, each giving one of its condition-level signatures | number of studies, taxa reported together |
| U6 size-matched | U6 draws whose three set sizes are each within max(2, 10%) of the observed | as U6, plus set sizes |

U1 to U3 have exact hypergeometric tails for each pair and an exact convolution
for the 3-way core. Every Monte Carlo null uses 100,000 draws with fixed seeds,
and p = (1 + draws at least as extreme) / (1 + draws). U1 to U4 treat taxa as
independent draws, which they are not: genera, their families and their phyla
are reported together, and some taxa are reported in most feces studies whatever
the condition. U6 keeps both features and is the reference null for claims about
overlap.

A sensitivity table repeats the overlap under other set definitions: profiles
with a positive cellulose z-score only, genus and species ranks only, taxa with
more decreased than increased profiles or studies, all body sites,
disease-versus-healthy contrasts, and every decreased feces taxon regardless of
trait. The trait benchmark ranks every matrix trait with at least 10 decreased
feces taxa in each condition by mean pairwise d_J. Its ranks describe this corpus
and are not a test; traits share taxa.

## Outputs

`cellulose_decr_analysis.py` writes the TSV and JSON tables listed in the
directory README; `cellulose_decr_figures.py` draws the figures from those tables
only. To regenerate from a checkout whose `data/` inputs are present:

```bash
COMMIT="$(git rev-parse HEAD)$(git diff --quiet HEAD -- data/analysis/cellulose_decr/scripts || echo +dirty)"
uv run python data/analysis/cellulose_decr/scripts/cellulose_decr_analysis.py --code-commit "$COMMIT"
uv run python data/analysis/cellulose_decr/scripts/cellulose_decr_figures.py
```

`--input-root` points the analysis at another checkout's `data/` when the
scripts run from a worktree. The scripts live beside their outputs for now; moving
them into `src/mtp` behind a `just` recipe is tracked separately.

## Limits

- Trait annotations are rolled up: a family or phylum carries the annotation if a
  descendant does. Eight of the 17 taxa in all three sets are such nodes.
- A taxon enters a set once any study reports it decreased, so the sets measure
  overlap among reported taxonomic labels, not consistent loss. Several core taxa
  are reported increased in more studies than decreased in some condition.
- Sampling effort differs greatly between compartments, so stronger results in
  feces do not show that the pattern is specific to feces.
- The three conditions and the trait are exploratory leads picked from this corpus
  (`EXPLORATORY_LEADS` in `src/mtp/deck_visualizations.py`), so none of these tests
  is confirmatory.
