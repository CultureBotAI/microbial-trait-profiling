# Six Shared Traits Interpretations

## Overall Interpretation

The grouped Cytoscape diagrams show that the six shared traits recur across broad disease domains, especially infectious disease, cancer/neoplasia, metabolic/endocrine disease, neurologic/cognitive disease, colorectal neoplasia, respiratory/airway disease, and reproductive/hormonal health. The individual-condition diagrams identify the condition labels driving those grouped patterns, including Parkinson's disease, colorectal cancer, COVID-19, obesity, human papillomavirus infection, gastric cancer, and dry eye syndrome.

These diagrams are descriptive maps of positive-z trait memberships. Their density shows cross-condition recurrence, not disease specificity or statistical significance. Conditions, studies, experiments, and profiles contribute unequal numbers of observations, and repeated observations are not independent. The figures therefore support interpreting the six traits as candidate cross-condition microbial phenotype signals rather than established shared disease biomarkers.

## Draft-paper figure description: grouped condition-body-site signed direction heatmap

**Figure X. Grouped condition-by-body-site signed direction contrasts across the six shared microbial traits** (`six_shared_traits_grouped_condition_body_site_signed_direction_contrast_clustered_heatmap.png`, per-experiment condition and body-site labels, #555). Rows are the 91 grouped condition-by-anatomy pairs, spanning 36 condition labels and 15 anatomy groups, that have at least one supported and testable trait contrast above the absolute-effect cutoff. Columns are the six traits shared by the composition and enrichment analyses. Each cell is the fraction of increased profiles with a positive trait z-score minus the corresponding fraction of decreased profiles. Red cells therefore favor increased profiles, blue cells favor decreased profiles, and the printed value gives the signed fraction difference on a fixed `-1` to `+1` scale. Grey cells are unavailable, under-supported, or have an absolute effect of `0.01` or less. Displayed cells require at least three positive-z profiles in either direction and at least two contributing studies. Dark borders mark Fisher exact tests passing Benjamini-Hochberg correction within this grouped heatmap view (`q<=0.05`); unbordered colored cells are descriptive effects.

Of 872 candidate cells, 811 have finite Fisher p-values, 425 are displayed, 30 have nominal `p<0.01`, and eight retain dark borders after heatmap-wide correction. Five bordered cellulose contrasts favor decreased profiles: neurologic/cognitive disease in feces (`-0.235`, `q=0.0019`), infectious disease in feces (`-0.226`, `q=0.0346`), stroke in feces (`-0.550`, `q=0.0326`), cancer/neoplasia in feces (`-0.264`, `q=0.0326`), and colorectal neoplasia in feces (`-0.288`, `q=0.0346`). The three increased-profile contrasts are infectious disease with coccus shaped in feces (`+0.230`, `q=0.0346`), cancer/neoplasia with facultatively aerobic in feces (`+0.278`, `q=0.0326`), and neurologic/cognitive disease with moderately halophilic in feces (`+0.179`, `q=0.0346`). These cells reinforce the broader descriptive opposition between decreased-profile cellulose signal and increased-profile aerotolerance, morphology and halophily traits.

The left dendrogram orders the 80 clusterable rows by average-linkage clustering on pairwise-complete Spearman distance (`1 - rho`) across the six signed trait contrasts, requiring at least three shared supported traits per row pair. Blank cells are never filled with zero: the 11 rows with no defined Spearman distance to any other row are shown in a separate insufficient-overlap block, and remaining undefined row-pair distances among clusterable rows use neutral distance `1.0`. Because linkage distances have high tie excess across six traits, no dendrogram cut, cluster ID, or branch label is interpreted; the dendrogram orders rows by shared trait-direction similarity rather than defining discrete disease or body-site clusters.

The heatmap-wide q-values use profile counts and are specific to this grouped visualization. They are not equivalent to the manuscript's stricter global condition-trait screen, which collapses profiles to paired, equally weighted study contrasts and corrects across all eligible condition-trait tests. Because profiles and experiments can repeat within studies, the heatmap supports cross-condition pattern discovery rather than patient-level prevalence, disease specificity, or independent biomarker claims.

## Global Condition-Trait Evidence

A separate inferential screen evaluates the six traits across all eligible human fecal conditions. It includes `66` conditions with at least three studies containing both increased and decreased signatures, producing `396` condition-trait comparisons (`392` with finite p-values). Each comparison uses paired study-level increased-minus-decreased z scores and a two-sided Wilcoxon signed-rank test. Benjamini-Hochberg correction is applied across the complete family of `396` comparisons.

No condition-trait result passes the global `q<=0.05` threshold.

Six results fall in the prespecified exploratory tier (`0.05<q<=0.20`) and must not be called significant:

| condition | trait | direction and mean contrast | raw p | global q | paired studies |
| --- | --- | ---: | ---: | ---: | ---: |
| Parkinson's disease | moderately halophilic | increased, `+0.515` | `0.000207` | `0.0790` | `45` |
| Parkinson's disease | arsenate reduction | increased, `+0.435` | `0.00120` | `0.0790` | `45` |
| colorectal cancer | cellulose catabolic process | decreased, `-0.847` | `0.000876` | `0.0790` | `28` |
| COVID-19 | cellulose catabolic process | decreased, `-1.214` | `0.000977` | `0.0790` | `12` |
| colorectal cancer | moderately halophilic | increased, `+0.457` | `0.00115` | `0.0790` | `28` |
| COVID-19 | facultatively aerobic | increased, `+1.001` | `0.000977` | `0.0790` | `12` |

All other condition-trait combinations, including Parkinson's disease with cellulose catabolic process (`q=0.480`), remain descriptive after global correction. The earlier diagrams and unadjusted summaries help locate patterns for follow-up, but they do not supersede these global q-values.

See the [evidence-landscape methods](../PD-cellulose/pd_condition_taxon_trait_network_methods.md) and [complete condition-trait test table](../PD-cellulose/pd_condition_taxon_trait_network_condition_trait_tests.tsv).

## Increased vs Decreased

The increased diagrams are generally denser: grouped increased profiles have `n=4,261` positive-z memberships versus `n=3,623` for decreased; individual increased has `n=1,628` versus `n=1,352` decreased. These are membership counts, not effect-size estimates or independent-sample tests.

Five traits are descriptively stronger in increased profiles: facultatively aerobic, coccus shaped, sphere shaped, arsenate reduction, and moderately halophilic. Cellulose catabolic process is the main exception and is comparatively more prominent in decreased profiles. The global screen narrows this broad pattern substantially: no condition-trait contrast passes global `q<=0.05`, and only six contrasts (Parkinson's disease with moderately halophilic at `q=0.0790`, Parkinson's disease with arsenate reduction at `q=0.0790`, colorectal cancer with cellulose catabolic process at `q=0.0790`, COVID-19 with cellulose catabolic process at `q=0.0790`, colorectal cancer with moderately halophilic at `q=0.0790`, and COVID-19 with facultatively aerobic at `q=0.0790`) reach the exploratory `0.05<q<=0.20` tier.

Feces dominates both directions in the grouped diagrams, contributing `2,659` increased memberships versus `2,084` decreased. The largest increased-dominant fecal patterns involve facultatively aerobic, moderately halophilic, arsenate reduction, coccus shaped, and sphere shaped, whereas fecal cellulose catabolic process is decreased-dominant. Oral cavity/dental is one of the clearest descriptive decreased skews (`415` increased versus `463` decreased); ocular/eye, intestinal/GI mucosa, and stomach/esophagus lean more increased.

## Grouped vs Individual

Grouped diagrams are best for seeing broad disease-domain structure. They show the same six traits spanning infectious disease, cancer/neoplasia, metabolic/endocrine, neurologic/cognitive, and colorectal neoplasia categories.

Individual diagrams identify the disease anchors behind those groups: for example, Parkinson's disease within neurologic/cognitive disease, obesity within metabolic/endocrine disease, COVID-19 and HPV infection within infectious disease, and several cancer labels within cancer/neoplasia. The individual diagrams are sparser because they retain raw disease labels and fewer body-site nodes.

Grouping can make a signal look broadly shared even when one or a few conditions contribute much of its weight. Neither view provides an inferential ranking on its own, so condition-level claims should be checked against the globally corrected test table.

## Cross-Disease Patterns

The clearest descriptive cross-disease pattern is a fecal multi-trait axis spanning neurologic, metabolic, colorectal/GI, and infectious conditions. Other site-specific recurrences include COVID-19 across upper-airway/nasal and fecal samples, HPV infection across vaginal/reproductive and oral/dental samples, gastric cancer in stomach/esophagus, and dry eye syndrome and glaucoma in ocular/eye samples.

Taken together, the maps motivate a shared microbial phenotype hypothesis involving aerotolerance/oxygen exposure, coccoid or spherical morphology, halophily/osmotic tolerance, arsenate reduction, and cellulose metabolism. Increased profiles descriptively emphasize the aerobic/coccus/sphere/halophilic/arsenate pattern, while decreased profiles emphasize cellulose catabolic process and oral/dental-associated signals.

The evidence does not support treating that whole axis as a validated cross-disease signature. After global correction across the fecal condition-trait screen, no individual condition-trait contrast passes `q<=0.05` (six exploratory results tie at the smallest global `q=0.0790`; increased moderately halophilic signal in Parkinson's disease has the smallest nominal `p=0.000207`); the recurring patterns remain descriptive or exploratory hypotheses requiring independent validation.
