# Autism and Autoimmune Disease Follow-Up: Aerobe v2

This follow-up starts from aerobe v2 trait enrichment in increased profiles and inspects the BugSigDB taxa reported in the autism/autism spectrum disorder and autoimmune disease condition labels.

## Plot Guild Specificity

| plot | guild-specific? | basis |
| --- | --- | --- |
| `autism_autoimmune_taxa_signed_heatmap.png` | no | taxon direction-count heatmap; selected traits are annotation context |

## Included Condition Labels

| condition group | exact condition label | profiles | studies | increased / decreased profiles | increased - decreased trait z | body-site groups |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Autism / autism spectrum disorder | Autism | 14 | 4 | 8 / 6 | 0.755 | Feces |
| Autism / autism spectrum disorder | Autism spectrum disorder | 78 | 15 | 42 / 36 | -0.125 | Feces; Saliva |
| Autoimmune disease | Autoimmune disease | 8 | 2 | 4 / 4 | 0.423 | Feces |

## Top Taxa

The tables below count unique taxon/profile occurrences after deduplicating the same taxon within the same study/experiment/direction profile.

### Autism / autism spectrum disorder

Top increased taxa:

| taxon | rank | profiles | studies | selected traits |
| --- | --- | ---: | ---: | --- |
| taxon_186803 | genus | 7 | 6 | facultatively aerobic; microaerophilic |
| taxon_1224 | genus | 6 | 4 | facultatively aerobic; aerobic; microaerophilic |
| taxon_543 | genus | 6 | 4 | facultatively aerobic; aerobic; microaerophilic |
| taxon_1678 | genus | 6 | 3 | facultatively aerobic; aerobic; microaerophilic |
| taxon_1578 | genus | 5 | 4 | facultatively aerobic; aerobic; microaerophilic |
| taxon_375288 | genus | 5 | 4 |  |
| taxon_31953 | genus | 5 | 3 | facultatively aerobic; aerobic; microaerophilic |
| taxon_128827 | genus | 4 | 4 | facultatively aerobic; aerobic; microaerophilic |
| taxon_186802 | genus | 4 | 4 | facultatively aerobic; aerobic; microaerophilic |
| taxon_216572 | genus | 4 | 4 | facultatively aerobic; microaerophilic |
| taxon_1506553 | genus | 4 | 4 |  |
| taxon_216851 | genus | 4 | 3 | facultatively aerobic |

Top decreased taxa:

| taxon | rank | profiles | studies | selected traits |
| --- | --- | ---: | ---: | --- |
| taxon_186803 | genus | 9 | 6 | facultatively aerobic; microaerophilic |
| taxon_171552 | genus | 8 | 4 | aerobic |
| taxon_838 | genus | 7 | 5 |  |
| taxon_31977 | genus | 6 | 4 | aerobic |
| taxon_577309 | genus | 6 | 3 |  |
| taxon_39948 | genus | 5 | 5 | aerobic |
| taxon_201174 | genus | 5 | 4 | facultatively aerobic; aerobic; microaerophilic |
| taxon_171549 | genus | 5 | 3 | facultatively aerobic; aerobic; microaerophilic |
| taxon_158846 | genus | 5 | 3 |  |
| taxon_39485 | genus | 5 | 3 |  |
| taxon_228604 | genus | 5 | 2 |  |
| taxon_33042 | genus | 5 | 2 |  |

### Autoimmune disease

Top increased taxa:

| taxon | rank | profiles | studies | selected traits |
| --- | --- | ---: | ---: | --- |
| taxon_816 | genus | 3 | 2 | facultatively aerobic; aerobic; microaerophilic |
| taxon_28116 | genus | 3 | 2 |  |
| taxon_1224 | genus | 3 | 1 | facultatively aerobic; aerobic; microaerophilic |
| taxon_1236 | genus | 3 | 1 | facultatively aerobic; aerobic; microaerophilic |
| taxon_31977 | genus | 3 | 1 | aerobic |
| taxon_1239 | genus | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |
| taxon_1300 | genus | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |
| taxon_1301 | genus | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |
| taxon_171549 | genus | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |
| taxon_186801 | genus | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |
| taxon_186802 | genus | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |
| taxon_186826 | genus | 2 | 1 | facultatively aerobic; aerobic; microaerophilic |

Top decreased taxa:

| taxon | rank | profiles | studies | selected traits |
| --- | --- | ---: | ---: | --- |
| taxon_186806 | genus | 3 | 1 |  |
| taxon_817 | genus | 3 | 1 |  |
| taxon_818 | genus | 2 | 2 |  |
| taxon_186803 | genus | 2 | 1 | facultatively aerobic; microaerophilic |
| taxon_216572 | genus | 2 | 1 | facultatively aerobic; microaerophilic |
| taxon_729 | genus | 2 | 1 | facultatively aerobic; microaerophilic |
| taxon_171551 | genus | 2 | 1 | aerobic |
| taxon_165137 | genus | 2 | 1 |  |
| taxon_165261 | genus | 2 | 1 |  |
| taxon_1730 | genus | 2 | 1 |  |
| taxon_33038 | genus | 2 | 1 |  |
| taxon_39485 | genus | 2 | 1 |  |

## Red/Blue Pattern Labels

These lists make the profile-similarity interpretation auditable. Shared red labels have positive signed counts in both exact `Autism` and `Autoimmune disease`, supporting the increased-profile similarity readout. Shared blue labels would have negative signed counts in both exact labels. The ASD blue-contrast list marks taxa where the broader `Autism spectrum disorder` label is decreased-skewed while at least one exact comparison label is neutral or increased-skewed.

### Shared red labels: exact Autism and Autoimmune disease

| heatmap label | Autism | ASD | Autoimmune disease | selected traits | in heatmap |
| --- | ---: | ---: | ---: | --- | --- |
| NCBI:1224 genus | 3 | 3 | 3 | aerobic; facultatively aerobic; microaerophilic | yes |
| NCBI:543 genus | 1 | 1 | 2 | aerobic; facultatively aerobic; microaerophilic | yes |
| NCBI:200643 genus | 1 | -2 | 1 | aerobic; facultatively aerobic; microaerophilic | yes |
| NCBI:816 genus | 1 | 1 | 2 | aerobic; facultatively aerobic; microaerophilic | no |
| NCBI:976 genus | 1 | 0 | 1 | aerobic; facultatively aerobic; microaerophilic | no |
| NCBI:203492 genus | 1 | 0 | 1 |  | no |

### Shared blue labels: exact Autism and Autoimmune disease

| heatmap label | Autism | ASD | Autoimmune disease | selected traits | in heatmap |
| --- | ---: | ---: | ---: | --- | --- |
| NCBI:186803 genus | -1 | -1 | -1 | facultatively aerobic; microaerophilic | no |

### ASD blue-contrast labels

| heatmap label | Autism | ASD | Autoimmune disease | selected traits | in heatmap |
| --- | ---: | ---: | ---: | --- | --- |
| NCBI:171552 genus | 0 | -7 | -1 | aerobic | yes |
| NCBI:838 genus | 0 | -6 | 0 |  | yes |
| NCBI:577309 genus | 0 | -5 | 0 |  | yes |
| NCBI:171549 genus | -1 | -4 | 1 | aerobic; facultatively aerobic; microaerophilic | yes |
| NCBI:33039 genus | 0 | -4 | 0 |  | yes |
| NCBI:39488 genus | 0 | -4 | 0 |  | yes |
| NCBI:31977 genus | -1 | -3 | 3 | aerobic | yes |
| NCBI:1300 genus | 0 | -3 | 1 | aerobic; facultatively aerobic; microaerophilic | yes |
| NCBI:39485 genus | 0 | -3 | -2 |  | yes |
| NCBI:201174 genus | 0 | -3 | 0 | aerobic; facultatively aerobic; microaerophilic | yes |
| NCBI:158846 genus | 0 | -3 | 0 |  | yes |
| NCBI:283168 genus | 0 | -3 | 0 |  | yes |
| NCBI:853 genus | 0 | -3 | 0 | facultatively aerobic | yes |
| NCBI:39491 genus | 0 | -3 | -1 |  | yes |
| NCBI:47678 genus | 0 | -3 | -1 |  | yes |

## Inferred Guild Evidence

This section uses a data-derived trait-enrichment view: for each exact condition label, it identifies v2 traits enriched in that condition's increased profiles, then lists taxa observed in those same profiles that also belong to the KG-Microbe-wide candidate universe and carry the enriched trait annotations. The universe contains 2,667 taxa with at least one core oxygen/respiration criterion annotation, not the fixed focal taxa list.

| condition | increased profiles | universe taxa / observed taxa | enriched traits | supporting taxa |
| --- | ---: | ---: | --- | --- |
| Autism | 8 | 24 / 36 | aerobic (+1.05); facultatively aerobic (+0.80); microaerophilic (+0.42) | taxon_1224 (3 traits, 3 profiles); taxon_1350 (3 traits, 3 profiles); taxon_31953 (3 traits, 3 profiles); taxon_85023 (3 traits, 3 profiles); taxon_1653 (3 traits, 2 profiles); taxon_976 (3 traits, 2 profiles); taxon_1678 (3 traits, 2 profiles); taxon_128827 (3 traits, 1 profile); taxon_1578 (3 traits, 1 profile); taxon_1716 (3 traits, 1 profile); taxon_1760 (3 traits, 1 profile); taxon_186817 (3 traits, 1 profile) |
| Autoimmune disease | 4 | 45 / 84 | facultatively aerobic (+0.61); aerobic (+0.45); microaerophilic (+0.21) | taxon_816 (3 traits, 3 profiles); taxon_1224 (3 traits, 3 profiles); taxon_1236 (3 traits, 3 profiles); taxon_1239 (3 traits, 2 profiles); taxon_1300 (3 traits, 2 profiles); taxon_1301 (3 traits, 2 profiles); taxon_171549 (3 traits, 2 profiles); taxon_186801 (3 traits, 2 profiles); taxon_186802 (3 traits, 2 profiles); taxon_186826 (3 traits, 2 profiles); taxon_200643 (3 traits, 2 profiles); taxon_543 (3 traits, 2 profiles) |

## Interpretation

This v2 follow-up uses only the narrowed aerobe v2 selected traits: facultatively aerobic, aerobic, microaerophilic. It supports a label-dependent signal check rather than a broad causal claim.

Included condition labels with positive increased-minus-decreased v2 trait-enrichment shifts are `Autism` (0.755), `Autoimmune disease` (0.423).

At the taxa-list level, overlap should still be interpreted cautiously. The signed heatmap counts direction balance by taxon/profile occurrence, while the guild annotation columns only indicate whether those taxa carry one or more of the v2 oxygen-tolerance traits.

`autism_autoimmune_taxa_signed_pattern_labels.tsv` extracts the heatmap labels that directly support the red/blue readout above so the visual pattern can be audited without manually reading the plot.

The most defensible interpretation is that the exact autism and autoimmune disease labels share a candidate gut dysbiosis axis involving increased-profile facultative or aerotolerant trait enrichment. This analysis does not show that autism and autoimmune disease are causally linked, driven by the same taxa, or even consistently connected across broader ASD metadata. The autoimmune subset is small, and many taxa are locally unresolved as `taxon_<NCBI id>`, so the taxa-level results should be treated as candidate NCBI IDs for follow-up rather than polished biological names.

## Heatmap

`autism_autoimmune_taxa_signed_heatmap.png` shows the top 45 genus/species taxa selected from 923 observed taxa. Cell values are `increased unique taxon-profile count - decreased unique taxon-profile count`; red is increased-skewed and blue is decreased-skewed. Colored x-axis labels and the marker strip above the heatmap identify taxa included in the red/blue pattern-label extraction.

## Files

- `autism_autoimmune_profile_summary.tsv`
- `autism_autoimmune_taxa_direction_summary.tsv`
- `autism_autoimmune_taxa_direction_summary_grouped.tsv`
- `autism_autoimmune_taxa_signed_balance.tsv`
- `autism_autoimmune_taxa_signed_heatmap_matrix.tsv`
- `autism_autoimmune_taxa_signed_pattern_labels.tsv`
- `autism_autoimmune_taxa_signed_heatmap.png`
- `autism_autoimmune_inferred_guild_summary.tsv`
- `autism_autoimmune_inferred_guild_trait_evidence.tsv`
- `autism_autoimmune_inferred_guild_taxa_evidence.tsv`
- `kg_microbe_candidate_taxon_universe.tsv`
