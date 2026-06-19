# Autism and Autoimmune Disease Follow-Up: Facultative-Aerobe v2

This follow-up starts from the facultative-aerobe increased-profile guild signal and inspects the BugSigDB taxa reported in the autism/autism spectrum disorder and autoimmune disease condition labels.

## Plot Guild Specificity

| plot | guild-specific? | basis |
| --- | --- | --- |
| `autism_autoimmune_taxa_signed_heatmap.png` | no | taxon direction-count heatmap; guild traits are annotation context |

## Included Condition Labels

| condition group | exact condition label | profiles | studies | increased / decreased profiles | increased - decreased guild z | body-site groups |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Autism / autism spectrum disorder | Autism | 14 | 4 | 8 / 6 | 0.755 | Feces |
| Autism / autism spectrum disorder | Autism spectrum disorder | 78 | 15 | 42 / 36 | -0.125 | Feces; Saliva |
| Autoimmune disease | Autoimmune disease | 8 | 2 | 4 / 4 | 0.423 | Feces |

## Top Taxa

The tables below count unique taxon/profile occurrences after deduplicating the same taxon within the same study/experiment/direction profile.

### Autism / autism spectrum disorder

Top increased taxa:

| taxon | rank | profiles | studies | guild traits |
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

| taxon | rank | profiles | studies | guild traits |
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

| taxon | rank | profiles | studies | guild traits |
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

| taxon | rank | profiles | studies | guild traits |
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

## Inferred Guild Evidence

This section uses a data-derived guild view: for each exact condition label, it identifies v2 guild traits enriched in that condition's increased profiles, then lists taxa observed in those same profiles that also belong to the KG-Microbe-wide candidate universe and carry the enriched trait annotations. The universe contains 2,667 taxa with at least one core oxygen/respiration criterion annotation, not the fixed focal taxa list.

| condition | increased profiles | universe taxa / observed taxa | enriched traits | supporting taxa |
| --- | ---: | ---: | --- | --- |
| Autism | 8 | 24 / 36 | aerobic (+1.05); facultatively aerobic (+0.80); microaerophilic (+0.42) | taxon_1224 (3 traits, 3 profiles); taxon_1350 (3 traits, 3 profiles); taxon_31953 (3 traits, 3 profiles); taxon_85023 (3 traits, 3 profiles); taxon_1653 (3 traits, 2 profiles); taxon_976 (3 traits, 2 profiles); taxon_1678 (3 traits, 2 profiles); taxon_128827 (3 traits, 1 profile); taxon_1578 (3 traits, 1 profile); taxon_1716 (3 traits, 1 profile); taxon_1760 (3 traits, 1 profile); taxon_186817 (3 traits, 1 profile) |
| Autoimmune disease | 4 | 45 / 84 | facultatively aerobic (+0.61); aerobic (+0.45); microaerophilic (+0.21) | taxon_816 (3 traits, 3 profiles); taxon_1224 (3 traits, 3 profiles); taxon_1236 (3 traits, 3 profiles); taxon_1239 (3 traits, 2 profiles); taxon_1300 (3 traits, 2 profiles); taxon_1301 (3 traits, 2 profiles); taxon_171549 (3 traits, 2 profiles); taxon_186801 (3 traits, 2 profiles); taxon_186802 (3 traits, 2 profiles); taxon_186826 (3 traits, 2 profiles); taxon_200643 (3 traits, 2 profiles); taxon_543 (3 traits, 2 profiles) |

## Interpretation

This v2 follow-up uses only the narrowed facultative-aerobe guild traits: facultatively aerobic, aerobic, microaerophilic. It supports a label-dependent signal check rather than a broad causal claim.

Included condition labels with positive increased-minus-decreased v2 guild shifts are `Autism` (0.755), `Autoimmune disease` (0.423).

At the taxa-list level, overlap should still be interpreted cautiously. The signed heatmap counts direction balance by taxon/profile occurrence, while the guild annotation columns only indicate whether those taxa carry one or more of the v2 oxygen-tolerance traits.

The most defensible interpretation is that the exact autism and autoimmune disease labels share a candidate gut dysbiosis axis involving increased facultative or aerotolerant guild signal. This analysis does not show that autism and autoimmune disease are causally linked, driven by the same taxa, or even consistently connected across broader ASD metadata. The autoimmune subset is small, and many taxa are locally unresolved as `taxon_<NCBI id>`, so the taxa-level results should be treated as candidate NCBI IDs for follow-up rather than polished biological names.

## Heatmap

`autism_autoimmune_taxa_signed_heatmap.png` shows the top 45 genus/species taxa selected from 923 observed taxa. Cell values are `increased unique taxon-profile count - decreased unique taxon-profile count`; red is increased-skewed and blue is decreased-skewed.

## Files

- `autism_autoimmune_profile_summary.tsv`
- `autism_autoimmune_taxa_direction_summary.tsv`
- `autism_autoimmune_taxa_direction_summary_grouped.tsv`
- `autism_autoimmune_taxa_signed_balance.tsv`
- `autism_autoimmune_taxa_signed_heatmap_matrix.tsv`
- `autism_autoimmune_taxa_signed_heatmap.png`
- `autism_autoimmune_inferred_guild_summary.tsv`
- `autism_autoimmune_inferred_guild_trait_evidence.tsv`
- `autism_autoimmune_inferred_guild_taxa_evidence.tsv`
- `kg_microbe_candidate_taxon_universe.tsv`
