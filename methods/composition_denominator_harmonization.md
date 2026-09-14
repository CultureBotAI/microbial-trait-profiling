# Harmonized composition denominator (#174)

## Definition and analysis

For profile i and retained trait j, define f_ij = c_ij / sum_k(c_ik),
where k ranges over **all non-MediaDive columns in the input matrix**.
Columns whose identifiers contain `mediadive.medium`, case-insensitively,
are removed before computing totals, support summaries, fraction variances,
selection scores, or transformed values. The denominator is not restricted
to eligible, displayed, or tested traits, nor to the six Table 1 traits.
These fractions describe shares of trait-annotation counts, not organism
abundances or percentages of taxa.

Profiles with zero retained counts, including MediaDive-only profiles, remain
in the fixed comparison corpus. Their retained fractions and log1p inputs are
zero. They participate in column standardization and pairing; their final
z-scores need not be zero. No profile or pair is dropped solely to remove a
zero denominator in this rerun. A constant column receives zero z-scores.
The general heatmap CLI still supports explicit profile-support filters;
those filters now use non-MediaDive counts. The audit freezes the original
profile keys instead of reapplying filters that could change the cohort.

Within each analysis stratum, transform fractions as log1p(1000 * f_ij),
then center and scale each trait using the population standard deviation
(ddof=0). Subtract decreased-profile z-scores from the increased-profile
z-scores within each experiment and apply the original two-sided one-sample
t-test across pairs. Retain the original minimum support of 6 profiles,
maximum support of 95% of profiles, positive fraction variance, uncapped
statistical trait selection, and separate Benjamini-Hochberg correction
within each stratum. Ranking uses fraction variance times log1p(total
retained-trait signal for that column). Experiment pairs are not reweighted
by study; issue #133 remains a distinct modeling question.

The diagnostic profile and pair heatmaps display 60 ranked traits with the
same three required display traits as the original implementation. They
standardize against the full 5,829-profile reference, unlike the statistical
tables, which standardize within each paired analysis stratum. Their
display selection does not define a statistical correction family.

## Reproduction and provenance

Run `just sync`, `just composition-denominator-audit`, and `just check`.
The audit requires the repository history containing composition commit
`a2b2df98e9beec98b772c9cb5af7aa090d033a98` and ORA commit
`3ae7ac847d5d1ebc3bdd388bb2793d5e60bfdc81`; fetch the analysis branches if
using a shallow or single-branch clone. It extracts snapshots into
`.tmp/issue174/` and does not replace `data/derived/study_trait_profiles_wide.tsv`.

The matrix is verified as Git blob
`d837f449e7e967c478d10f92938923e57eec9a19`, SHA256
`c07d589f6cc88724e43c683ebd132bcd5c7f316f36b95df3ff9acad24ebbf62f`.
Frozen paired metadata holds human-host, study/experiment inc/dec pairing
and the original study-level condition/body-site labels fixed. This yields
2,335 pairs from 959 studies. Only the five strata emitted by the original
run are compared, not older orphaned statistics files in its Git tree.

Input hashes, implementation base commit and source-file hashes, package
versions, denominator definition, and pairing filters are recorded in
`trait-visualization/denominator-harmonization/provenance.json`; tested family
sizes are in `correction_families.tsv`. The original BugSigDB provenance
records v1.3.1, tag commit `5815c2638d6b56134ac8db19edc2bb3c4b17f1ca`.
The KG-Microbe provenance records edges MD5
`555ecf5ba5d89288cdc44e368375374f` and nodes MD5
`dd76b4656f9341db25a7c4afa96cd58e`. This is a rerun from the exact derived
matrix, not a new raw-KG ingestion: the original large KG TSVs were not
available locally. Source lineage is recorded, not independently rederived.

Before accepting the corrected results, the audit reproduces baseline
eligibility, effects, p-values, and q-values in all five strata within
rtol=1e-9 and atol=1e-12. It also reproduces the committed ORA q <= 0.001
membership and q-values. The ORA rerun extends the reporting cutoff to 0.05
without changing its tests, background, fold floor, or global BH family.

## Figure captions

**Profile heatmap.** Column z-scores of log1p(1000 times the trait fraction),
with each fraction normalized by all non-MediaDive trait-annotation counts
in its profile before trait selection. Rows are clustered display bins of
the fixed 5,829-profile corpus; colors are capped at +/-1 for display only.
A colored cell represents a transformed trait fraction, not a p-value or
a declaration of significance.

**Paired heatmap.** Increased-minus-decreased profile z-score contrasts for
2,335 matched study/experiment pairs from 959 studies. Profile fractions
use all non-MediaDive counts before selection; z-scores use the full-profile
reference. Display bins summarize adjacent clustered pairs and colors are
capped at +/-1. Positive values indicate greater transformed trait fraction
on the increased side. Color alone does not establish statistical significance.

**Composition infographic.** Paired transformed-fraction contrasts using
the harmonized non-MediaDive denominator. Statistical tests use the paired
stratum reference, with BH correction across all eligible tested traits
within that stratum. Displaying top-ranked traits does not make the figure
a thresholded significant set.

## Integration and interpretation

The [measured report](../trait-visualization/denominator-harmonization/README.md)
documents the unchanged pooled six-trait q <= 0.001 overlap and the changed
q <= 0.05 sets. Do not recompute BH after choosing the overlap or selecting
the six traits. Compare membership at the declared cutoff while retaining
each result's original correction family. P-values, q-values, and effect
units are not interchangeable across the three frameworks.

Regenerated producer caches replace their conventional paths under
`trait-visualization/`. Old Sankey, PD, guild, and integration figures on
other branches must not be represented as harmonized merely because their
six trait identifiers are unchanged. See the report's downstream audit for
the caches rebuilt here and the separate renderers still requiring reruns.
