# Table 1 source contributions and lineage influence

Reconciled 4,670 human paired profiles from 959 studies across all 317 non-MediaDive traits. Enumerated 1233 resolved genus/family lineages; 22 entered the predeclared panel with at least three supporting studies. See manifest.json for hashes, exact scopes and protocol.

The full contributing universe, concentration and all 15 cross-trait overlaps precede any display cutoff. Effects include all six raw annotation keys, both curated directions, shared fractions, all 15 log ratios and Noel's native zero-total convention. Each scenario preserves matched-pair changes and eligibility changes separately.

These are provisional annotation influence diagnostics. Contrast orientation (#188) and source semantics (#189) remain unreviewed. Explicit reported ranks only resolve part of the pinned taxonomy; unresolved rank and source identity remain visible. Original source assertions cannot be counted from KG edge IDs alone. Higher-rank direct annotations are retained and flagged when their allocation to an excluded lineage is ambiguous.

Native Jed transformed statistics and native ORA reruns remain for their owners; no substitute native q-values are supplied. Fixed and recomputed ORA backgrounds are distinct export contracts. Conditional intervals and separate BH6/BY15 families do not validate the full selection pipeline (#182). Loss of significance is not proof of explanation. These outputs feed #178 and the evidence limits in #179; they do not validate the Rosalind mechanisms, oxygen exposure, SCFA production, or causal ecological readouts.

Largest contributing resolved genus, with share of **all** propagated counts:

| Trait key | Direction | Genus | Count share | Genus resolution coverage |
| --- | --- | --- | --- | --- |
| METPO:1000608 | inc | Streptococcus | 22.7% | 95.0% |
| METPO:1000608 | dec | Streptococcus | 23.6% | 94.4% |
| METPO:2000103\|GO:0030245 | inc | Ruminococcus | 11.2% | 80.3% |
| METPO:2000103\|GO:0030245 | dec | Ruminococcus | 13.1% | 82.5% |
| METPO:1000668 | inc | Staphylococcus | 54.8% | 97.7% |
| METPO:1000668 | dec | Staphylococcus | 58.2% | 97.7% |
| METPO:1000683 | inc | Streptococcus | 20.9% | 76.7% |
| METPO:1000683 | dec | Streptococcus | 17.8% | 77.9% |
| METPO:2000103\|pathways:arsenate_reduction | inc | Escherichia | 58.7% | 88.1% |
| METPO:2000103\|pathways:arsenate_reduction | dec | Escherichia | 41.7% | 82.4% |
| METPO:1000623 | inc | Vibrio | 8.9% | 76.6% |
| METPO:1000623 | dec | Vibrio | 7.4% | 74.9% |

Measured run: implementation commit `8c78a5a`, protocol `190-v1`, pinned KG build
20260120. Input and implementation hashes match the completed local manifest.
Regenerate the full backing tables using the recipe in
[the methods](table1_taxon_influence.md); generated data remain under
`data/derived/table1_taxon_influence/`.

The contribution ledger contains 73,758 reported/source/trait rows and 12,154
supporting KG edge records; neither is a patient or organism count. Shared fractions
retain 2,241 pairs from 939 studies. The native Noel convention starts from all
2,335 pairs / 959 studies: 94 pairs have a zero full total in at least one arm.
The complete output contains 2,484 scenario/scope/representation/feature rows.

Pooled shared fraction mean contrasts below are **percentage points**, with
conditional 95% Student t intervals across equally weighted study averages.
Directional-share intervals and sign-test families are separate fields in the TSV.
In particular, a directional significance statement does not imply that the
mean-effect interval excludes zero.

| Trait key | Baseline mean [95% CI] | Pairs / studies |
| --- | --- | --- |
| METPO:1000608 | 0.8901 [0.5992, 1.1811] | 2241 / 939 |
| METPO:2000103\|GO:0030245 | -0.0636 [-0.1209, -0.0064] | 2241 / 939 |
| METPO:1000668 | 0.7249 [0.3415, 1.1082] | 2241 / 939 |
| METPO:1000683 | 0.0213 [-0.0174, 0.0601] | 2241 / 939 |
| METPO:2000103\|pathways:arsenate_reduction | 0.0064 [0.0033, 0.0096] | 2241 / 939 |
| METPO:1000623 | 0.0186 [0.0044, 0.0328] | 2241 / 939 |

Largest absolute matched changes **within the declared 22-scenario panel** follow.
This is a descriptive selection, not a calibrated ranking or an exhaustive lineage
screen. The panel uses global denominator contribution counts plus the four named
predictions; it does not necessarily contain each trait's leading genus (for example,
Ruminococcus and Vibrio). The full universe is available for separately declared
extensions. These conditional intervals do not account for selecting this display.
All opposite, null and untestable results remain in the backing TSVs.

| Trait key | Excluded lineage | Matched change [95% CI], pp | Remaining pairs / studies | Newly lost pairs |
| --- | --- | --- | --- | --- |
| METPO:1000608 | Enterobacteriaceae | -0.2494 [-0.3948, -0.1040] | 2203 / 931 | 38 |
| METPO:2000103\|GO:0030245 | Clostridiaceae | -0.0098 [-0.0195, -0.0001] | 2221 / 937 | 20 |
| METPO:1000668 | Streptococcaceae | -0.3717 [-0.5990, -0.1445] | 2215 / 935 | 26 |
| METPO:1000683 | Lactobacillaceae | 0.0139 [0.0001, 0.0277] | 2193 / 919 | 48 |
| METPO:2000103\|pathways:arsenate_reduction | Escherichia | -0.0065 [-0.0092, -0.0037] | 2237 / 937 | 4 |
| METPO:1000623 | Enterococcaceae | -0.0082 [-0.0121, -0.0043] | 2233 / 938 | 8 |

These tables illustrate why count concentration and fraction influence are distinct:
Staphylococcus contributes the largest coccus annotation mass, while removing
Streptococcaceae produces the largest matched coccus fraction change in this panel.
Every exclusion changes the full denominator, and newly excluded pairs also change
study averages. For example, the facultative-aerobic mean changes by -0.2494 pp on
matched pairs after Enterobacteriaceae removal, with a further -0.0770 pp from
eligibility changes. This does not establish that oxygen exposure or a biological
mechanism caused the original associations.

Native Jed transformed results and native ORA reruns remain missing pending owner
integration. The run does not complete #188, #189, #182, or all acceptance criteria
of #190. #178 must preserve these missing native results; #179 must preserve the
annotation-level inference, anatomical/corpus restrictions and unresolved provenance.
