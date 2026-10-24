# Paired Contrast Interpretation Notes

This note consolidates the interpretation guidance discussed for the paired
trait-profile contrast outputs in this folder. It is intended as a companion to
`../figure_methods.md`, the paired contrast TSVs, and the increased/decreased
profile guild figures.

## Key Output Files

- `trait_profile_paired_contrast_trait_stats.tsv`: pooled paired contrast trait
  statistics.
- `trait_profile_paired_contrast_trait_stats_pooled.tsv`: explicit pooled copy
  of the same analysis.
- `trait_profile_paired_contrast_trait_stats_<body_site_group>.tsv`: body-site
  or robustness-stratified paired contrast trait statistics (`feces`,
  `feces_excluded`, `other`, and `saliva` are the four current manifest-listed
  strata with `>= 50` pairs; `caecum`, `colon`, `mouth`,
  `subgingival_dental_plaque`, and `vagina` are legacy `< 50`-pair TSVs retained
  on disk from before the human-host-only and experiment-level body-site
  restriction in #555).
- `trait_profile_paired_contrast_trait_stats_manifest.tsv`: emitted strata
  meeting the `>= 50`-pair threshold (`feces`, `feces_excluded`, `other`,
  `saliva`), pair counts, study counts, trait counts, and paths.
- `trait_profile_paired_contrast_condition_summary.tsv`: condition-level paired
  contrast load summary.
- `trait_profile_paired_contrast_body_site_summary.tsv`: body-site-level paired
  contrast load summary.
- `aerobe-enrichment-v2/condition_increased_guild_signal.png`: condition
  figure for the increased-profile aerobe v2 story.
- `aerobe-enrichment-v2/anatomical_sites_diagram.png`: body-site figure for
  the increased-profile aerobe v2 story.
- `anaerobe-depletion-v2/condition_decreased_guild_signal.png`: condition
  figure for the decreased-profile anaerobe story.
- `anaerobe-depletion-v2/anatomical_sites_diagram.png`: body-site figure for
  the decreased-profile anaerobe story.

## Methods Paragraph

The paired contrast data were generated from
`data/derived/study_trait_profiles_wide.tsv` using matched BugSigDB
study/experiment profiles that contained both an `increased` and a `decreased`
direction. For each analysis set, including pooled data, Feces, Feces-excluded,
and each body-site group with at least 50 matched pairs (`feces`,
`feces_excluded`, `other`, and `saliva` in
`trait_profile_paired_contrast_trait_stats_manifest.tsv`), traits were selected
within that set by requiring annotation support in at least 6 profiles, presence
in no more than 95% of profiles, nonzero composition variance, and excluding
MediaDive medium traits; all eligible traits were retained with no top-60 cap.
Raw trait counts were converted to per-profile fractions of total trait signal,
scaled per 1,000 counts, log1p-transformed, and standardized to trait-wise
z-scores. For each matched study/experiment pair, the contrast value was
calculated as the mean increased-profile z-score minus the mean
decreased-profile z-score for each trait. Trait-level summary tables report
descriptive paired contrast statistics, including mean, median, standard
deviation, directional pair fractions, one-sample t-test against zero,
Benjamini-Hochberg q-value, annotation support, and direction label. These
p-values are screening summaries of curated profile-level trait contrasts, not
new differential-abundance tests from raw sequencing data.

## What the Mouth and Legacy Body-Site Labels Mean

In `trait_profile_paired_contrast_trait_stats_mouth.tsv`, `mouth` means the
`body_site_group` stratum is exactly `Mouth`. It is not a broad synonym search
for all oral labels. Separate body-site groups such as `Saliva` and
`Subgingival dental plaque` have their own TSVs. The `Mouth` file uses matched
increased/decreased study-experiment pairs whose upstream paired metadata has
`body_site_group = Mouth`. In the current human-host experiment-level corpus
(#555), only `Feces` (`1,418` pairs), `Feces excluded` (`839` pairs), `Other`
(`480` pairs), and `Saliva` (`156` pairs) meet the `>= 50`-pair manifest
threshold; the `caecum`, `colon`, `mouth`, `subgingival_dental_plaque`, and
`vagina` TSVs on disk are legacy outputs from before that restriction and are
not listed in `trait_profile_paired_contrast_trait_stats_manifest.tsv`.

The filename token is produced by lowercasing and sanitizing the stratum label,
so `Mouth` becomes `mouth`. The body-site grouping itself is inherited from the
upstream paired metadata, where common body-site labels are kept as groups and
rarer labels are collapsed to `Other`.

## Interpreting Trait Rows

Each trait row is a profile-level directional signal, not proof that the trait
was directly measured or biologically increased in raw samples. A trait appears
in a profile because taxa in the curated BugSigDB signature map to that trait
through the derived KG-Microbe trait annotations.

Use these columns together:

- `mean_increased_minus_decreased`: average direction and size of the paired
  contrast.
- `positive_pair_fraction`: fraction of matched pairs where the trait is higher
  in the increased profile.
- `negative_pair_fraction`: fraction of matched pairs where the trait is higher
  in the decreased profile.
- `p_value`: descriptive one-sample t-test of paired contrasts against zero.
- `q_value`: Benjamini-Hochberg adjusted p-value across traits in that TSV.
- `n_profiles_where_annotated` and `annotated_profile_fraction`: annotation
  coverage for the trait.
- `profile_support`, `total_signal`, and `selection_score`: support and
  selection context for the trait in that analysis set.

A positive `mean_increased_minus_decreased` means the trait is more represented
in increased profiles than in matched decreased profiles. A negative value means
the trait is more represented in decreased profiles than in matched increased
profiles. Near-zero values indicate that the trait may be common but does not
strongly separate the two profile directions in matched contrasts.

## Interpreting Z-Scores

The z-scores are standardized profile-level trait signals. They are computed
after raw trait counts are converted to fractions of each profile's total trait
signal, scaled, log-transformed, and standardized by trait column.

- `z = 0`: average trait or trait-set signal relative to the analyzed profile
  universe.
- `z > 0`: above-average trait or trait-set signal.
- `z < 0`: below-average trait or trait-set signal.

A negative z-score in an increased-profile diagram does not mean negative
abundance, and it does not mean the profile is decreased. The word `increased`
refers to the BugSigDB signature direction. The z-score refers to whether the
selected trait or trait set is above or below average among profile-level trait
signals.

For the increased-profile aerobe v2 diagram, positive mean trait-set
z-scores support the idea that increased-signature taxa are relatively enriched
for aerobic, facultative, or microaerophilic traits. Negative values mean the
increased-signature taxa in that condition or body site are not especially
aerobe-like relative to the broader dataset.

For the decreased-profile anaerobe diagram, positive mean trait-set z-scores
support the idea that decreased-signature taxa are relatively enriched for anaerobic or
obligately anaerobic traits. Negative values mean the decreased-signature taxa
are not especially anaerobe-like in that context.

## Conditions and Body Sites

Conditions and body sites are grouping labels for profiles. They are not used to
create the original profile-level z-scores. The z-scores are calculated first,
then profiles are grouped by BugSigDB `condition` or by `body_site_group`, and
the relevant z-scores are averaged within each group.

For conditions, a condition-level mean trait-set z-score answers: in profiles
associated with this condition, is the selected trait or trait-set signal above
or below the broader profile-level average?

For body sites, a body-site mean trait-set z-score answers the same question after
grouping profiles by anatomical context such as `Feces`, `Mouth`, `Saliva`,
`Colon`, or `Vagina`.

The main directional contrasts are:

- Increased-profile analysis: `increased mean trait-set z - decreased mean trait-set z`.
- Decreased-profile analysis: `decreased mean trait-set z - increased mean trait-set z`.

Interpretation grid:

| Pattern | Interpretation |
| --- | --- |
| Focus mean z > 0 and focus-minus-comparator > 0 | Strongest support: above-average trait signal in the expected direction. |
| Focus mean z > 0 and contrast <= 0 | Trait is enriched in the condition or body site, but not specifically in the focus direction. |
| Focus mean z < 0 and contrast > 0 | Focus direction has the larger mean than the comparator, but still below average overall. |
| Focus mean z < 0 and contrast < 0 | Does not support the proposed trait-direction story for that condition or body site. |

## Relationship to Single-Trait Enrichment

The single-trait enrichment analysis and paired contrast analysis answer related
but different questions.

The single-trait enrichment analysis is mainly a prevalence or enrichment view:
does this trait show up often or strongly among taxa in the profiles?

The paired contrast analysis is a directionality view: when a study/experiment
has both increased and decreased signatures, is this trait higher on one side?

Integration grid:

| Single-trait enrichment | Paired contrast | Interpretation |
| --- | --- | --- |
| High enrichment | Positive paired contrast | Strong candidate trait for increased-profile biology. |
| High enrichment | Negative paired contrast | Common trait, but more associated with decreased profiles in matched contrasts. |
| High enrichment | Near-zero paired contrast | Broadly common trait, but not direction-specific. |
| Low or moderate enrichment | Strong paired contrast | Direction-specific but less prevalent; possible niche or body-site-specific signal. |
| Low enrichment | Near-zero contrast | Weak support for the profile story. |

For the increased-profile aerobe v2 story, the strongest support comes
from traits that are enriched in the single-trait analysis and also have
positive `mean_increased_minus_decreased`, positive-pair fraction above 0.5,
adequate annotation coverage, and supporting body-site or condition context.

For the decreased-profile anaerobe story, the strongest support comes from
traits that are enriched in the single-trait analysis and have negative
`mean_increased_minus_decreased`, high negative-pair fraction, adequate
annotation coverage, and agreement in relevant strata such as feces, colon, or
caecum.

Avoid simply combining p-values from the two analyses because they are not
independent and they answer different questions. A better integration is an
evidence table with one row per trait, containing enrichment rank, enrichment
q-value, paired contrast mean, paired q-value, directional pair fractions,
annotation coverage, supporting strata, and an interpretation class such as
`concordant increased`, `concordant decreased`, `common not directional`,
`directional niche`, or `weak/unsupported`.

## Practical Interpretation Rules

- Treat all p-values and q-values as descriptive screening statistics.
- Require annotation coverage before interpreting absence of signal.
- Prefer traits that show both directionality and consistency across pairs.
- Use body-site TSVs to identify context-specific signals that are hidden in the
  pooled analysis.
- Do not interpret a positive z-score as raw abundance or direct measurement.
- Do not interpret a negative z-score in an increased or decreased profile as
  the opposite BugSigDB direction; it only means below-average trait signal for
  that profile group.
- Conditions may mix body sites, studies, hosts, and designs, so condition
  z-scores should be treated as screening summaries rather than disease effect
  sizes.
