# Six Shared Traits Interpretations

## Overall Interpretation

The grouped Cytoscape diagrams show that the six shared traits recur across broad disease domains, especially infectious disease, cancer/neoplasia, metabolic/endocrine disease, neurologic/cognitive disease, colorectal neoplasia, respiratory/airway disease, and reproductive/hormonal health. The individual-condition diagrams identify the condition labels driving those grouped patterns, including Parkinson's disease, colorectal cancer, COVID-19, obesity, human papillomavirus infection, gastric cancer, and dry eye syndrome.

These diagrams are descriptive maps of positive-z trait memberships. Their density shows cross-condition recurrence, not disease specificity or statistical significance. Conditions, studies, experiments, and profiles contribute unequal numbers of observations, and repeated observations are not independent. The figures therefore support interpreting the six traits as candidate cross-condition microbial phenotype signals rather than established shared disease biomarkers.

## Draft-paper figure description: grouped condition-body-site signed direction heatmap

**Figure X. Grouped condition-by-body-site signed direction contrasts across the six shared microbial traits** (`clustered-heatmap-generation` at `4a418ab`). Rows are the 78 grouped condition-by-anatomy pairs, spanning 34 condition labels and 15 anatomy groups, that have at least one supported and testable trait contrast above the absolute-effect cutoff. Columns are the six traits shared by the composition and enrichment analyses. Each cell is the fraction of increased profiles with a positive trait z-score minus the corresponding fraction of decreased profiles. Red cells therefore favor increased profiles, blue cells favor decreased profiles, and the printed value gives the signed fraction difference on a fixed `-1` to `+1` scale. Grey cells are unavailable, under-supported, or have an absolute effect of `0.01` or less. Displayed cells require at least three positive-z profiles in either direction and at least two contributing studies. Dark borders mark Fisher exact tests passing Benjamini-Hochberg correction within this grouped heatmap view (`q<=0.05`); unbordered colored cells are descriptive effects.

Of 747 tested cells, 384 are displayed, 29 have nominal `p<0.01`, and eight retain dark borders after heatmap-wide correction. Four bordered cellulose contrasts favor decreased profiles: neurologic/cognitive disease in feces (`-0.231`, `q=0.00096`), infectious disease in feces (`-0.246`, `q=0.0093`), stroke in feces (`-0.556`, `q=0.0148`), and cancer/neoplasia in feces (`-0.247`, `q=0.0350`). The four increased-profile contrasts are infectious disease with coccus shaped in feces (`+0.257`, `q=0.0093`), neurologic/cognitive disease with moderately halophilic in feces (`+0.183`, `q=0.0209`), age/development/sex with facultatively aerobic in feces (`+0.315`, `q=0.0350`), and pregnancy/neonatal with sphere shaped in feces (`+0.475`, `q=0.0484`). These cells reinforce the broader descriptive opposition between decreased-profile cellulose signal and increased-profile aerotolerance, morphology, and stress-tolerance traits.

The left dendrogram orders rows by average-linkage clustering of Spearman-correlation distance across the six signed trait contrasts; missing or unsupported cells are treated as neutral zero only for clustering. Branch labels are shown only when a recurrent trait-direction pattern passes a separate branch-pattern FDR test, has at least three matching displayed cells, reaches at least two-thirds directional agreement, and contains at least one heatmap-wide BH-selected cell. Five branches qualify, producing seven nested labels for increased-associated facultatively aerobic, coccus shaped, sphere shaped, or moderately halophilic patterns. No branch qualifies for a condition- or anatomy-enrichment label, so the dendrogram should be interpreted as shared trait-direction structure rather than a taxonomy of diseases or body sites.

The heatmap-wide q-values use profile counts and are specific to this grouped visualization. They are not equivalent to the manuscript's stricter global condition-trait screen, which collapses profiles to paired, equally weighted study contrasts and corrects across all eligible condition-trait tests. Because profiles and experiments can repeat within studies, the heatmap supports cross-condition pattern discovery rather than patient-level prevalence, disease specificity, or independent biomarker claims.

## Global Condition-Trait Evidence

A separate inferential screen evaluates the six traits across all eligible human fecal conditions. It includes `56` conditions with at least three studies containing both increased and decreased signatures, producing `336` condition-trait comparisons. Each comparison uses paired study-level increased-minus-decreased z scores and a two-sided Wilcoxon signed-rank test. Benjamini-Hochberg correction is applied across the complete family of `336` tests.

Only one condition-trait result passes the global `q<=0.05` threshold: Parkinson's disease with moderately halophilic, favoring increased profiles (`mean contrast=+0.530`, raw `p=0.000143`, `q=0.048`, `n=45` paired studies).

Three additional results fall in the prespecified exploratory tier (`0.05<q<=0.20`) and must not be called significant:

| condition | trait | direction and mean contrast | raw p | global q | paired studies |
| --- | --- | ---: | ---: | ---: | ---: |
| Parkinson's disease | arsenate reduction | increased, `+0.506` | `0.00159` | `0.164` | `45` |
| COVID-19 | cellulose catabolic process | decreased, `-1.307` | `0.00195` | `0.164` | `11` |
| colorectal cancer | cellulose catabolic process | decreased, `-0.821` | `0.00194` | `0.164` | `24` |

All other condition-trait combinations, including Parkinson's disease with cellulose catabolic process (`q=0.597`), remain descriptive after global correction. The earlier diagrams and unadjusted summaries help locate patterns for follow-up, but they do not supersede these global q-values.

See the [evidence-landscape methods](../PD-cellulose/pd_condition_taxon_trait_network_methods.md) and [complete condition-trait test table](../PD-cellulose/pd_condition_taxon_trait_network_condition_trait_tests.tsv).

## Increased vs Decreased

The increased diagrams are generally denser: grouped increased profiles have `n=4,495` positive-z memberships versus `n=3,795` for decreased; individual increased has `n=1,645` versus `n=1,434` decreased. These are membership counts, not effect-size estimates or independent-sample tests.

Five traits are descriptively stronger in increased profiles: facultatively aerobic, coccus shaped, sphere shaped, arsenate reduction, and moderately halophilic. Cellulose catabolic process is the main exception and is comparatively more prominent in decreased profiles. The global screen narrows this broad pattern substantially: only the Parkinson's disease-moderately-halophilic contrast is FDR-significant, while the COVID-19-cellulose and colorectal cancer-cellulose contrasts are exploratory.

Feces dominates both directions in the grouped diagrams, contributing `2,857` increased memberships versus `2,260` decreased. The largest increased-dominant fecal patterns involve facultatively aerobic, moderately halophilic, arsenate reduction, coccus shaped, and sphere shaped, whereas fecal cellulose catabolic process is decreased-dominant. Oral cavity/dental is one of the clearest descriptive decreased skews (`396` increased versus `478` decreased); ocular/eye, intestinal/GI mucosa, and stomach/esophagus lean more increased.

## Grouped vs Individual

Grouped diagrams are best for seeing broad disease-domain structure. They show the same six traits spanning infectious disease, cancer/neoplasia, metabolic/endocrine, neurologic/cognitive, and colorectal neoplasia categories.

Individual diagrams identify the disease anchors behind those groups: for example, Parkinson's disease within neurologic/cognitive disease, obesity within metabolic/endocrine disease, COVID-19 and HPV infection within infectious disease, and several cancer labels within cancer/neoplasia. The individual diagrams are sparser because they retain raw disease labels and fewer body-site nodes.

Grouping can make a signal look broadly shared even when one or a few conditions contribute much of its weight. Neither view provides an inferential ranking on its own, so condition-level claims should be checked against the globally corrected test table.

## Cross-Disease Patterns

The clearest descriptive cross-disease pattern is a fecal multi-trait axis spanning neurologic, metabolic, colorectal/GI, infectious, and diet-related conditions. Other site-specific recurrences include COVID-19 across upper-airway/nasal and fecal samples, HPV infection across vaginal/reproductive and oral/dental samples, gastric cancer in stomach/esophagus, and dry eye syndrome and glaucoma in ocular/eye samples.

Taken together, the maps motivate a shared microbial phenotype hypothesis involving aerotolerance/oxygen exposure, coccoid or spherical morphology, halophily/osmotic tolerance, arsenate reduction, and cellulose metabolism. Increased profiles descriptively emphasize the aerobic/coccus/sphere/halophilic/arsenate pattern, while decreased profiles emphasize cellulose catabolic process and oral/dental-associated signals.

The evidence does not support treating that whole axis as a validated cross-disease signature. After global correction, the only confirmed condition-trait direction in the fecal screen is increased moderately halophilic signal in Parkinson's disease; the remaining recurring patterns are descriptive or exploratory hypotheses requiring independent validation.
