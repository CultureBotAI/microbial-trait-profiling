# Pinned-corpus validation of the 15 Table 1-trait log-ratios (matched-corpus analog)

> [!IMPORTANT]
> **Matched-Corpus Table 1-Trait Analog (Issue #210)**
> The pinned-corpus analyses reported here are a matched-corpus Table 1-trait analog sensitivity analysis, not a reproduction of the published Table 1 baseline. They use the pinned human paired KG-Microbe profile corpus (959 studies, 2,335 pairs) and do not reproduce published Table 1 baseline values (e.g. facultatively aerobic $q = 2.8 \times 10^{-10}$).

The 2026-09-13 rerun reports all 15 unordered pairs among the six Table 1-named traits.
Three ratios have BY q <= 0.05 and none have BY q <= 0.001. These are exploratory,
relative annotation comparisons conditional on positive counts. Both cutoffs
are reported for transparency; no cutoff was selected to define a preferred set.

This snapshot compares Noel's original fraction-direction analysis with the
paired log-ratio extension in [PR #181](https://github.com/CultureBotAI/microbial-trait-profiling/pull/181).
The [three-framework comparison](three_arm_statistics.md) describes Jed's and
Ilyas's different hypotheses. The [log-ratio methods](table1_paired_log_ratios.md)
define eligibility, study weighting, ties and the correction family.

## Corpus and reproduction

Both recipes used the same pinned profile and experiment files: 7,329 profile rows
read, 1,500 dropped as non-human and 1,159 experiments dropped as unpaired, leaving
4,670 retained profiles, 2,335 matched experiments and 959 unique studies. The
retained corpus is human and paired. The fraction denominator sums the 317
non-MediaDive columns; 2,333 MediaDive columns are excluded. The log-ratio
contrast uses only its named numerator and reference counts.

The matrix these counts were measured on is
`data/derived/study_trait_profiles_wide.tsv`, sha256
`6e70724ee807cf24e9c1f371a51421fa5b09d6283fa7f3abcd6620e5492c6a73`: 2,653
header columns (3 identifier, 2,333 MediaDive, 317 retained) over 7,329 profile
rows, joined from BugSigDB `2026-04-24_00:41_UTC` and a KG-Microbe build whose
nodes are md5 `dd76b4656f9341db25a7c4afa96cd58e` and which records no release
tag. Naming the dump does not identify the matrix, so the run writes this
record to `corpus_identity.json` beside its outputs (#222).


| Provenance field | Value |
| --- | --- |
| BugSigDB vintage | `2026-04-24_00:41_UTC` (BugSigDBExports v1.3.1) |
| Analysis source revision | [`b1f66ea`](https://github.com/CultureBotAI/microbial-trait-profiling/commit/b1f66ea4efcad88c82af06ce1561674f5dda5f5d) |
| Log-ratio version | `paired-log-ratio-sign-v1` |
| Profile matrix SHA256 | `c07d589f6cc88724e43c683ebd132bcd5c7f316f36b95df3ff9acad24ebbf62f` |
| Experiment metadata SHA256 | `403afea98d9b0a5b56939f68bcc9aace56a8eda350232a586806fe7bf30ad230` |
| Log-ratio module SHA256 | `df664a9dade43092d7d5b03df70554b6b29340ef8af335ecd7547736bf3515cf` |
| Shared loader SHA256 | `70ecbb78e03f1b34649fd8dc84cd486fb9e0817c59d31fadcd2d00cdae794683` |
| Tested level | All body sites and conditions pooled; one contrast per study |

The module hash moved after `b1f66ea` because the provenance block gained the
shared-loader hash; the statistics are untouched. Rerunning on the same pinned
inputs reproduces `trait_log_ratio_stats.tsv`, `trait_log_ratio_studies.tsv` and
`trait_log_ratio_report.md` byte-identically, so the table below still holds.

Regenerate with the existing recipes on these input files:

```bash
just table1-trait-maps
just table1-compositional-stats
```

The first recipe also renders figures using the project's installed Chrome.
Its pooled direction statistics use the retained profiles before heatmap display
cuts. Neither table below establishes significance at a condition or body site.

The source tables are
`data/derived/table1_traits/table1_traits_direction_shift.tsv` and
`data/derived/table1_composition/trait_log_ratio_stats.tsv`.
The latter directory also contains per-study coverage/effects, provenance JSON
and a generated Markdown report. These outputs are gitignored; this checked-in
snapshot makes the results reviewable without regenerating figures. Values below
are rounded; retain the source precision for threshold comparisons.

## Trait identities

Abbreviations keep the comparison table readable; they are not new trait groups.

| Label | Trait | Identifier resolved by the loader |
| --- | --- | --- |
| FA | facultatively aerobic | `METPO:1000608` |
| CEL | cellulose catabolic process | `GO:0030245` |
| COC | coccus shaped | `METPO:1000668` |
| SPH | sphere shaped | `METPO:1000683` |
| ARS | arsenate reduction | `pathways:arsenate_reduction` |
| HAL | moderately halophilic | `METPO:1000623` |

Qualified source columns are resolved by the shared loader. In particular,
cellulose and arsenate columns carry a `METPO:2000103|` prefix in this matrix.

## Original fraction directions (matched-corpus analog)

Each study contributes the sign of its mean increased-minus-decreased trait
fraction. The null concerns the share of non-tied studies shifting up, and BH
adjustment spans these six traits. Ties are retained in the counts and excluded
from the directional test. The exact 95% interval is for the share up.

| Trait | Up/down/tied studies | Scored studies | Share up (95% CI) | p | BH q, six traits |
| --- | --- | --- | --- | --- | --- |
| FA | 534/383/42 | 917 | 0.582 (0.550, 0.614) | 6.901e-7 | 4.140e-6 |
| CEL | 334/440/185 | 774 | 0.432 (0.396, 0.467) | 0.0001571 | 0.0003142 |
| COC | 499/416/44 | 915 | 0.545 (0.512, 0.578) | 0.00668 | 0.00668 |
| SPH | 467/377/115 | 844 | 0.553 (0.519, 0.587) | 0.002169 | 0.002982 |
| ARS | 346/270/343 | 616 | 0.562 (0.521, 0.601) | 0.002485 | 0.002982 |
| HAL | 485/360/114 | 845 | 0.574 (0.540, 0.608) | 1.928e-5 | 5.785e-5 |

All six have BH q <= 0.05; FA, CEL and HAL also have BH q <= 0.001. This describes
the matched-corpus Table 1-trait analog fraction hypotheses and their correction
family, not a reproduction of the published Table 1 baseline.

## All 15 paired log-ratios (matched-corpus analog)

For A/B, the experiment contrast is
`log2(A_inc/B_inc) - log2(A_dec/B_dec)`; eligible contrasts are averaged once per
study. Both traits must be positive on both sides of a complete experiment pair.
Positive contrasts mean A/B is higher in the increased signature.

Every row has 959 total studies and 2,335 total paired experiments. The number in
parentheses is the subset of eligible studies that lost some experiments to zero
exclusions. Entirely excluded studies equal 959 minus eligible studies. Up,
down and tied counts partition eligible studies; only up and down enter the test.
All 15 comparisons have status `tested` in this corpus.

| A/B | Eligible studies (partial) | Eligible pairs | Up/down/tied studies | Mean log2 change | Share up (95% CI) | p | BY q, 15 pairs |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FA/CEL | 467 (211) | 780 | 278/189/0 | 0.791 | 0.595 (0.549, 0.640) | 4.434e-5 | 0.001103 |
| FA/COC | 676 (247) | 1305 | 343/332/1 | 0.115 | 0.508 (0.470, 0.546) | 0.7003 | 1 |
| FA/SPH | 579 (240) | 1050 | 319/258/2 | 0.316 | 0.553 (0.511, 0.594) | 0.01243 | 0.1237 |
| FA/ARS | 255 (130) | 382 | 138/117/0 | 0.313 | 0.541 (0.478, 0.604) | 0.2103 | 1 |
| FA/HAL | 582 (242) | 998 | 309/272/1 | 0.304 | 0.532 (0.490, 0.573) | 0.1352 | 0.7479 |
| CEL/COC | 470 (210) | 792 | 216/252/2 | -0.575 | 0.462 (0.416, 0.508) | 0.1056 | 0.6569 |
| CEL/SPH | 429 (193) | 705 | 174/253/2 | -0.538 | 0.407 (0.360, 0.456) | 0.000154 | 0.002555 |
| CEL/ARS | 243 (124) | 358 | 88/153/2 | -0.848 | 0.365 (0.304, 0.429) | 3.383e-5 | 0.001103 |
| CEL/HAL | 402 (187) | 623 | 174/226/2 | -0.535 | 0.435 (0.386, 0.485) | 0.01068 | 0.1237 |
| COC/SPH | 609 (239) | 1133 | 324/281/4 | 0.161 | 0.536 (0.495, 0.576) | 0.08764 | 0.6231 |
| COC/ARS | 256 (131) | 385 | 122/133/1 | -0.096 | 0.478 (0.416, 0.542) | 0.5313 | 1 |
| COC/HAL | 554 (239) | 942 | 303/250/1 | 0.242 | 0.548 (0.505, 0.590) | 0.02693 | 0.2234 |
| SPH/ARS | 253 (128) | 381 | 121/130/2 | -0.218 | 0.482 (0.419, 0.546) | 0.6137 | 1 |
| SPH/HAL | 502 (231) | 830 | 249/248/5 | -0.034 | 0.501 (0.456, 0.546) | 1 | 1 |
| ARS/HAL | 255 (130) | 384 | 127/127/1 | -0.050 | 0.500 (0.437, 0.563) | 1 | 1 |

The mean includes tied eligible studies. The exact 95% interval concerns the
directional share among non-tied studies, not uncertainty in that mean; it is
also not a simultaneous interval adjusted across the 15 pairs. A mean change and
a directional share can differ in sign around their nulls: SPH/HAL has a mean
of -0.034 but a share up of 0.501. Its sign test does not test the mean effect.

## Implications for integration

FA/CEL shifts upward, while CEL/SPH and CEL/ARS shift downward, at BY q <= 0.05.
These results describe cellulose annotation representation relative to those
specific references among eligible studies. The other cellulose ratios also
have shares below 0.5, but their adjusted tests do not meet that cutoff.

FA and COC both have upward original fraction directions, yet FA/COC has a
directional share of 0.508 and BY q = 1. Such results can coexist: an annotation's
share of the total and its ratio to another annotation answer different questions.
The latter result does not establish equivalence or absence of a relative effect.

Positive-count eligibility retains 243-676 studies per ratio (25.3%-70.5% of the
959 studies) and can retain only part of a study's experiments. Consequently,
differences from the fraction results mix changes in the tested quantity,
eligible population and multiplicity procedure. They cannot quantify the
denominator mismatch in Jed's code or isolate the effect of normalization.

Before interpreting stability across references, Noel's
[#183](https://github.com/CultureBotAI/microbial-trait-profiling/issues/183) will compare the full
fraction sample with fractions and ratios on identical eligible pairs, and
assess common support and study influence. Jed's denominator harmonization is
tracked separately in [#174](https://github.com/CultureBotAI/microbial-trait-profiling/issues/174).
Calibration and selection remain follow-up work in
[#182](https://github.com/CultureBotAI/microbial-trait-profiling/issues/182).

The integrated evidence table in
[#178](https://github.com/CultureBotAI/microbial-trait-profiling/issues/178) must preserve these
results as separate hypotheses, with their references, denominators, populations,
uncertainty types and families. BY adjustment within the 15 ratios does not
establish FDR control for a union with the fraction or enrichment discoveries.
Agreement on shared annotations is not independent replication; these results
do not measure patient-level abundance, pathway activity or causation.
