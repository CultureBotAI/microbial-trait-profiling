# Matched-corpus Table 1-trait analog source contributions and lineage influence

> [!IMPORTANT]
> **Matched-Corpus Table 1-Trait Analog (Issue #210)**
> The source contributions and lineage influence results reported here are a matched-corpus Table 1-trait analog sensitivity analysis, not influence on the published manuscript Table 1. They use the pinned human paired KG-Microbe profile corpus and are not a reproduction of the published Table 1 baseline (they do not reproduce published Table 1 baseline values such as facultatively aerobic $q = 2.8 \times 10^{-10}$).

This summary reports the extended 190-v2 contribution/shared-effect run and its
[native extension](table1_native_influence.md). The source ledger and baseline effects
are unchanged from 190-v1; the rule-declared exclusion panel and all scenario-dependent
shared and native outputs were rerun.

Reconciled 4,514 human paired profiles from 932 studies across all 317 non-MediaDive traits of the `off` scope (`count_input_16s_wms_v1`). Enumerated 1,218 resolved genus/family lineages; 24 entered the rule-declared panel with at least three supporting studies. The per-trait leading-genus extension was [recorded on #190](https://github.com/CultureBotAI/microbial-trait-profiling/issues/190#issuecomment-5703892447) before this run. See manifest.json for hashes, exact scopes and protocol.

The full contributing universe, concentration and all 15 cross-trait overlaps precede any display cutoff. Effects include all six raw annotation keys, both curated directions, shared fractions, all 15 log ratios and Noel's native zero-total convention. Each scenario preserves matched-pair changes and eligibility changes separately.

These are provisional annotation influence diagnostics. Contrast orientation (#188) and source semantics (#189) remain unreviewed. Explicit reported ranks only resolve part of the pinned taxonomy; unresolved rank and source identity remain visible. Original source assertions cannot be counted from KG edge IDs alone. Higher-rank direct annotations are retained and flagged when their allocation to an excluded lineage is ambiguous.

Native transformed statistics and native ORA reruns are published in the
[compact review bundle](table1-native-results/README.md). Fixed and recomputed ORA
backgrounds remain distinct analyses. Conditional intervals and separate native,
BH6 and BY15 families do not validate the full selection pipeline (#182). Loss of
significance is not proof of explanation. These outputs feed #178 and the evidence
limits in #179; they do not validate the Rosalind mechanisms, oxygen exposure, SCFA
production, or causal ecological readouts.

Largest contributing resolved genus, with share of **all** propagated counts:

| Trait key | Direction | Genus | Count share | Genus resolution coverage |
| --- | --- | --- | --- | --- |
| METPO:1000608 | inc | Streptococcus | 22.8% | 94.9% |
| METPO:1000608 | dec | Streptococcus | 23.6% | 94.4% |
| METPO:2000103\|GO:0030245 | inc | Ruminococcus | 11.4% | 80.2% |
| METPO:2000103\|GO:0030245 | dec | Ruminococcus | 13.2% | 82.6% |
| METPO:1000668 | inc | Staphylococcus | 54.6% | 97.6% |
| METPO:1000668 | dec | Staphylococcus | 58.0% | 97.7% |
| METPO:1000683 | inc | Streptococcus | 21.0% | 76.5% |
| METPO:1000683 | dec | Streptococcus | 17.8% | 77.9% |
| METPO:2000103\|pathways:arsenate_reduction | inc | Escherichia | 58.9% | 88.2% |
| METPO:2000103\|pathways:arsenate_reduction | dec | Escherichia | 41.9% | 82.5% |
| METPO:1000623 | inc | Vibrio | 8.9% | 76.7% |
| METPO:1000623 | dec | Vibrio | 7.3% | 74.9% |

Measured run: implementation commit `8b280a3`, protocol `190-v2`, pinned KG build
20260120. Input and implementation hashes match the completed local manifest.
Regenerate the full backing tables using the recipe in
[the methods](table1_taxon_influence.md); generated data remain under
`data/derived/table1_taxon_influence/`.

The contribution ledger contains 73,081 reported/source/trait rows and 12,154
supporting KG edge records; neither is a patient or organism count. Shared fractions
retain 2,202 pairs from 922 studies. The native Noel convention starts from all
2,257 pairs / 932 studies: 55 pairs have a zero full total in at least one arm.
The complete output contains 2,700 scenario/scope/representation/feature rows.

Pooled shared fraction mean contrasts below are **percentage points**, with
conditional 95% Student t intervals across equally weighted study averages.
Directional-share intervals and sign-test families are separate fields in the TSV.
In particular, a directional significance statement does not imply that the
mean-effect interval excludes zero.

| Trait key | Baseline mean [95% CI] | Pairs / studies |
| --- | --- | --- |
| METPO:1000608 | 0.8853 [0.5931, 1.1775] | 2202 / 922 |
| METPO:2000103\|GO:0030245 | -0.0631 [-0.1213, -0.0048] | 2202 / 922 |
| METPO:1000668 | 0.6544 [0.2715, 1.0373] | 2202 / 922 |
| METPO:1000683 | 0.0197 [-0.0195, 0.0590] | 2202 / 922 |
| METPO:2000103\|pathways:arsenate_reduction | 0.0064 [0.0032, 0.0096] | 2202 / 922 |
| METPO:1000623 | 0.0174 [0.0038, 0.0311] | 2202 / 922 |

Largest absolute matched changes **within the rule-declared 24-scenario panel** follow.
This is a descriptive selection, not a calibrated ranking or an exhaustive lineage
screen. The panel uses global denominator contribution counts, the four named
predictions and each trait's leading genus. The last rule added Ruminococcus for
cellulose catabolism and Vibrio for moderate halophily. These conditional intervals do
not account for selecting this display. All opposite, null and untestable results
remain in the backing TSVs.

| Trait key | Excluded lineage | Matched change [95% CI], pp | Remaining pairs / studies | Newly lost pairs |
| --- | --- | --- | --- | --- |
| METPO:1000608 | Enterobacteriaceae | -0.2393 [-0.3863, -0.0924] | 2166 / 915 | 36 |
| METPO:2000103\|GO:0030245 | Clostridiaceae | -0.0100 [-0.0199, -0.0001] | 2182 / 920 | 20 |
| METPO:1000668 | Streptococcaceae | -0.3704 [-0.6010, -0.1398] | 2177 / 918 | 25 |
| METPO:1000683 | Lactobacillaceae | 0.0141 [0.0001, 0.0281] | 2155 / 902 | 47 |
| METPO:2000103\|pathways:arsenate_reduction | Escherichia | -0.0064 [-0.0091, -0.0037] | 2199 / 920 | 3 |
| METPO:1000623 | Enterococcaceae | -0.0085 [-0.0124, -0.0046] | 2195 / 921 | 7 |

The two newly added leading-genus scenarios produced small pooled matched changes in
both the shared fraction and native fixed-baseline-SD representations; every listed
conditional interval includes zero. The native q-value is the post-exclusion refitted
native model's full-family q-value, not a test of the matched change.

| Trait key | Excluded lineage | Shared matched change [95% CI], pp | Remaining pairs / studies; lost | Native fixed-SD matched change [95% CI] | Native after q |
| --- | --- | --- | --- | --- | --- |
| METPO:2000103\|GO:0030245 | Ruminococcus | -0.0006 [-0.0273, 0.0261] | 2194 / 921; 8 | 0.00412 [-0.02116, 0.02941] | 7.67e-7 |
| METPO:1000623 | Vibrio | -0.0007 [-0.0022, 0.0009] | 2202 / 922; 0 | -0.00487 [-0.01088, 0.00114] | 2.66e-4 |

Study-weighted native ORA matched-fold intervals also include zero under both
backgrounds: Ruminococcus is 0.00253 [-0.27579, 0.28085] with the fixed baseline and
-0.08074 [-0.36045, 0.19898] after recomputing the remaining background; Vibrio is
-0.00443 [-0.01417, 0.00531] and 0.00442 [-0.00709, 0.01594], respectively. These
diagnostics do not show a large matched change from either leading-genus removal, but
they do not establish lineage independence or a causal mechanism.

## Cross-scale interpretation flags

Two rows require explicit qualification before these results enter an integrated table
or manuscript claim. Their exact source fields are published in
[`table1-native-results/interpretation_handoff.tsv`](table1-native-results/interpretation_handoff.tsv).

**Arsenate reduction is Escherichia-sensitive on this corpus.** Escherichia supplies
58.9% of the increased-arm propagated counts. Its exclusion changes the matched
study-weighted raw-fraction mean by -0.0064 percentage points, leaving -0.0007 pp
[-0.0019, 0.0005]. The matched-corpus native transformed q-value likewise changes
from 4.24e-05 to 0.728 under this genus exclusion. This supports describing an
Escherichia-sensitive annotation signal, not a lineage-independent arsenate mechanism.
Taxon removal is an influence diagnostic rather than causal attribution, and unresolved
higher-rank assertions remain in the denominator.

**Sphere shaped does not reproduce as a nonzero study-weighted raw-fraction mean.**
Its baseline mean is 0.0197 pp [-0.0195, 0.0590], so that interval crosses zero before
any exclusion. This does not make every sphere-shaped statistic null: the separate
study-level directional sign test has q = 0.00213, and the experiment-weighted native
transformed analysis has q = 6.91e-05. The tests differ in transformation, weighting,
estimand and uncertainty. The raw-fraction interval is therefore cross-scale
non-reproduction, not refutation or independent replication of the native result.

These tables illustrate why count concentration and fraction influence are distinct:
Staphylococcus contributes the largest coccus annotation mass, while removing
Streptococcaceae produces the largest matched coccus fraction change in this panel.
Every exclusion changes the full denominator, and newly excluded pairs also change
study averages. For example, the facultative-aerobic mean changes by -0.2393 pp on
matched pairs after Enterobacteriaceae removal, with a further -0.0757 pp from
eligibility changes. This does not establish that oxygen exposure or a biological
mechanism caused the original associations.

The run does not complete #188, #189, #182, #184 or native-owner review. #178 must
preserve the distinction among shared, native composition and the two ORA backgrounds;
#179 must preserve the annotation-level inference, anatomical/corpus restrictions and
unresolved provenance.
