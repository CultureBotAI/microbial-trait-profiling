# Paired Composition Arm Calibration against Permutation Null (#263)

**Issue Reference:** [#263](https://github.com/CultureBotAI/microbial-trait-profiling/issues/263)  
**Lead (Composition Arm Calibration):** Jed Dongjin Kim-Ozaeta  
**Related Issues:** [#182](https://github.com/CultureBotAI/microbial-trait-profiling/issues/182), [#264](https://github.com/CultureBotAI/microbial-trait-profiling/issues/264), [#118](https://github.com/CultureBotAI/microbial-trait-profiling/issues/118), [#179](https://github.com/CultureBotAI/microbial-trait-profiling/issues/179)  
**Machine-Readable Summary Artifacts:**  
- Summary JSON: [`data/derived/calibration/composition/composition_calibration_summary.json`](../data/derived/calibration/composition/composition_calibration_summary.json)  
- Summary TSV: [`data/derived/calibration/composition/composition_calibration_summary.tsv`](../data/derived/calibration/composition/composition_calibration_summary.tsv)  
- Per-Replicate TSV: [`data/derived/calibration/composition/composition_calibration_replicates.tsv`](../data/derived/calibration/composition/composition_calibration_replicates.tsv)  
- Null p-values: `data/derived/calibration/composition/composition_null_pvalues.tsv.gz` (683 KB, not committed; regenerate with `just calibrate-composition`)  
- Table 1 Calibration: [`data/derived/calibration/composition/table1_composition_calibration.tsv`](../data/derived/calibration/composition/table1_composition_calibration.tsv)  

---

## 1. Executive Summary & Problem Formulation

Issue [#263](https://github.com/CultureBotAI/microbial-trait-profiling/issues/263) posed the central calibration question for the paired composition arm:
> *"Does the paired composition arm control the empirical false discovery rate under a complete null, or does it suffer from permutation runaway? If permuted runs produce dozens of $q \le 0.001$ traits, calibration is settled."*

To resolve this question conclusively, we benchmarked the complete direction-label permutation null (`complete_null_study_flip` stream, seed 182, 100 Monte Carlo replicates) across the full corpus of 959 paired BugSigDB studies, 2,335 paired experiments, and 317 non-MediaDive traits under the `off` scope (31,700 total hypothesis tests under the null):

1. **Zero Runaway at $q \le 0.001$:**  
   Across all 100 null replicates, the mean number of false discoveries at Benjamini-Hochberg $q \le 0.001$ was **0.00** (max: 0; 0/100 replicates produced any call; empirical family-wise error = **0.0%**, 95% Clopper-Pearson CI: 0.0%–3.6%). Within the 100-replicate empirical null benchmark, no permutation runaway was observed at stringent discovery thresholds.

2. **Conservative Control at $q \le 0.05$:**  
   At the standard Benjamini-Hochberg $q \le 0.05$ threshold, the mean number of false discoveries was **0.09** out of 317 `off`-scope tests (max: 7 in a single replicate; 97 out of 100 replicates produced exactly zero false calls; empirical family-wise error = **3.0%**, 95% CI: 0.6%–8.5%). With Benjamini-Yekutieli (BY) correction, the mean false call count fell to **0.02** (family-wise error = **1.0%**). None of the 100 null replicates reached or approached the 132 calls observed in unpermuted data.

3. **Nominal Type I Error Rate:**  
   Across all 31,700 null tests, the empirical Type I error rate matched or stayed strictly below nominal rates:
   - Nominal $\alpha = 0.05 \implies \text{empirical } 0.0460$
   - Nominal $\alpha = 0.01 \implies \text{empirical } 0.0081$
   - Nominal $\alpha = 0.001 \implies \text{empirical } 0.0007$

4. **Table 1 Canonical Trait Survival:**  
   Five of the six canonical Table 1 traits survive calibration with zero null replicates reaching the observed test statistic ($0 / 100$ null replicates reach $|t_{\text{obs}}|$):
   - *Cellulose catabolic process* ($t = -3.770, q = 1.06 \times 10^{-3}$)
   - *Coccus shaped* ($t = +3.764, q = 1.06 \times 10^{-3}$)
   - *Sphere shaped* ($t = +4.663, q = 5.37 \times 10^{-5}$)
   - *Arsenate reduction* ($t = +4.341, q = 1.66 \times 10^{-4}$)
   - *Moderately halophilic* ($t = +4.798, q = 3.54 \times 10^{-5}$)  
   *Facultatively aerobic* ($t = +1.663, p = 0.0967, q = 0.1836$) was not significant in the unpermuted paired composition arm (11/100 null replicates reached $|t_{\text{obs}}|$), consistent with the density adjustment findings in [#264](https://github.com/CultureBotAI/microbial-trait-profiling/issues/264) where this trait attenuated under curation-density controls.

---

## 2. Paired Composition Arm Mathematical Formulation

The paired composition arm evaluates systematic directional shifts in trait composition across paired BugSigDB studies. The inference pipeline proceeds as follows:

1. **Filtering to Paired Experiments:**  
   A profile dataset is restricted to experiment IDs that contain both an `increased` and a `decreased` signature. Unpaired single-direction experiments are excluded to preserve within-experiment matching.

2. **Log-Scale Transformation:**  
   Apply pseudocount scaling ($C = 1000.0$) and logarithmic transformation:
   $$y_{ij} = \log\left(1 + C \cdot f_{ij}\right)$$

3. **Standardization Across Profiles:**  
   Standardize transformed fractions across all paired profiles to zero mean and unit variance:
   $$z_{ij} = \frac{y_{ij} - \bar{y}_j}{s_j} \quad (\text{masked to } 0.0 \text{ if } s_j \le 10^{-12})$$

4. **Matched Experiment Contrast:**  
   For each experiment $p = (s, e)$ with both increased and decreased profiles:
   $$\Delta_{pj} = z_{\text{inc}, j}^{(p)} - z_{\text{dec}, j}^{(p)}$$

5. **Study-Level Aggregation:**  
   To prevent experiments within multi-experiment studies from violating statistical independence, average experiment contrasts within each study $s$:
   $$\bar{\Delta}_{sj} = \frac{1}{|P_s|} \sum_{p \in P_s} \Delta_{pj}$$

6. **Hypothesis Testing:**  
   Under the null hypothesis $H_0: \mathbb{E}[\bar{\Delta}_{sj}] = 0$, evaluate the one-sample Student's t-statistic across the $S = 959$ paired studies:
   $$t_j = \frac{\frac{1}{S} \sum_{s=1}^S \bar{\Delta}_{sj}}{\hat{s}_{\bar{\Delta}, j} / \sqrt{S}}, \quad p_j = 2 \cdot \left(1 - F_t(|t_j|, S - 1)\right)$$

7. **Multiple Testing Correction:**  
   Control FDR across the $P$ trait tests ($P = 317$ under `off`, 344 under `phenotype`) using Benjamini-Hochberg (BH) and Benjamini-Yekutieli (BY).

---

## 3. Direction-Label Permutation Null Design

The complete direction-label permutation null (`complete_null_study_flip`) tests the specific hypothesis that trait composition shifts correlate with curated direction labels (`increased` vs. `decreased`):

- **Study-Level Direction Flips:**  
  For each study $s \in \{1, \dots, S\}$ and replicate $r \in \{1, \dots, R\}$, sample an independent Rademacher random variable:
  $$\sigma_s^{(r)} \in \{-1, +1\}, \quad \mathbb{P}(\sigma_s = -1) = \mathbb{P}(\sigma_s = +1) = 0.5$$
  The null study mean contrast vector becomes:
  $$\bar{\Delta}_{sj}^{(r)} = \sigma_s^{(r)} \cdot \bar{\Delta}_{sj}$$
- **Exact Vectorization & Equivalence:**  
  Because negating all experiment contrasts within a study $\Delta'_{pj} = z_{\text{dec}, j}^{(p)} - z_{\text{inc}, j}^{(p)} = -\Delta_{pj}$ is mathematically equivalent to flipping $\bar{\Delta}_{sj} \to -\bar{\Delta}_{sj}$, the null t-statistics can be computed over the entire $(R, S, P)$ tensor instantaneously while guaranteeing mathematical identity to individual profile-level inversion down to machine epsilon ($\max |\Delta t| < 4.44 \times 10^{-16}$).
- **Preserved Properties under the Null:**  
  - Profile-level trait correlation structure (PC1 curation density and ontology overlaps) is fully preserved within each profile.
  - Taxon counts, signature sizes, and study sampling distributions are unaltered.
  - Directional alignment with condition labels is completely decoupled.

---

## 4. Quantitative Calibration Benchmark Results

Across 100 Monte Carlo replicates (`stream_key: complete_null_study_flip`, `seed: 182`), 959 studies, and 317 traits (`off`):

| Correction | Threshold | Observed Calls | Null Mean Calls | Null Max Calls | Replicates With $\ge 1$ Call | Family-Wise Error (95% CI) | Replicates Reaching Observed |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **BH** | $q \le 0.05$ | **132** | $0.09$ | 7 | 3 / 100 | **3.0%** (0.6%–8.5%) | 0 / 100 |
| **BH** | $q \le 0.01$ | **95** | $0.05$ | 5 | 1 / 100 | **1.0%** (0.0%–5.4%) | 0 / 100 |
| **BH** | $q \le 0.001$ | **49** | **0.00** | **0** | **0 / 100** | **0.0%** (0.0%–3.6%) | 0 / 100 |
| **BY** | $q \le 0.05$ | **90** | $0.02$ | 2 | 1 / 100 | **1.0%** (0.0%–5.4%) | 0 / 100 |
| **BY** | $q \le 0.01$ | **58** | **0.00** | **0** | **0 / 100** | **0.0%** (0.0%–3.6%) | 0 / 100 |
| **BY** | $q \le 0.001$ | **29** | **0.00** | **0** | **0 / 100** | **0.0%** (0.0%–3.6%) | 0 / 100 |

### Type I Error Rate across Null Hypotheses

Across the 31,700 null hypothesis tests:
- $\alpha = 0.05$: 1,457 rejections $\implies \hat{\alpha} = 0.0460$ (Nominal 0.050)
- $\alpha = 0.01$: 258 rejections $\implies \hat{\alpha} = 0.0081$ (Nominal 0.010)
- $\alpha = 0.001$: 21 rejections $\implies \hat{\alpha} = 0.0007$ (Nominal 0.001)

All empirical rejection rates are strictly bounded by nominal significance levels, demonstrating that one-sample study-level t-testing is well-calibrated and slightly conservative under empirical direction exchangeability.

---

## 5. Table 1 Canonical Traits Calibration Behavior

We evaluated the observed effect size, t-statistic, and survival under null calibration for the canonical Table 1 traits:

| Trait Identifier | Trait Label | Observed Mean Contrast | Observed $t$ | Observed BH $q$ | Null $t$ Range | Null Type I ($\alpha=0.05$) | Null Reaching $|t_{\text{obs}}|$ | Calibrated Survival |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `METPO:2000103\|GO:0030245` | Cellulose catabolic process | $-0.1287$ | **$-3.770$** | **$1.06 \times 10^{-3}$** | $[-2.018, +3.289]$ | $0.06$ | **0 / 100** | **Survives** |
| `METPO:1000668` | Coccus shaped | $+0.1296$ | **$+3.764$** | **$1.06 \times 10^{-3}$** | $[-2.421, +2.348]$ | $0.05$ | **0 / 100** | **Survives** |
| `METPO:1000683` | Sphere shaped | $+0.1503$ | **$+4.663$** | **$5.37 \times 10^{-5}$** | $[-3.822, +2.075]$ | $0.03$ | **0 / 100** | **Survives** |
| `METPO:2000103\|pathways:arsenate_reduction` | Arsenate reduction | $+0.1594$ | **$+4.341$** | **$1.66 \times 10^{-4}$** | $[-2.868, +2.597]$ | $0.04$ | **0 / 100** | **Survives** |
| `METPO:1000623` | Moderately halophilic | $+0.1594$ | **$+4.798$** | **$3.54 \times 10^{-5}$** | $[-2.239, +3.233]$ | $0.06$ | **0 / 100** | **Survives** |
| `METPO:1000608` | Facultatively aerobic | $+0.0555$ | $+1.663$ | $0.1836$ | $[-2.896, +3.320]$ | $0.09$ | 11 / 100 | **Not Significant** |

### Key Observations:
- **Consistent Convergence Across 5 Traits:**  
  Cellulose depletion, coccus shape, spherical morphology, arsenate reduction, and moderate halophily each show clear separation between the observed contrast and the 100-replicate empirical null distribution. In all five cases, no permuted replicate among the 100 reached the observed test statistic magnitude ($p_{\text{perm}} < 0.01$).
- **Facultatively Aerobic Concordance with #264:**  
  In the paired composition arm, *facultatively aerobic* achieves an observed $t = +1.663$ with $p = 0.0967$ and BH $q = 0.1836$, failing to reach significance even before permutation testing. Under null study flips, 11% of replicates exceed the observed test statistic. This aligns with the findings of Issue [#264](https://github.com/CultureBotAI/microbial-trait-profiling/issues/264): the apparent unadjusted association of *facultatively aerobic* is sensitive to annotation density and profile representation, and paired within-study composition contrast analysis does not detect a significant shift.

---

## 6. Methodological Takeaways & Recommendations

1. **Paired Composition Arm Empirical Null Calibration:**  
   Across 100 study-flip replicates, the paired composition adapter (`JedCompositionAdapter`, `jed_composition_ttest`) maintains empirical false discovery rates at or below nominal thresholds under the study-level direction-flip null. In this benchmark (#263), zero false discoveries are generated across 100 replicates at $q \le 0.001$, and 3% of null replicates generate at least one call at $q \le 0.05$.

2. **Study Aggregation is Statistically Essential:**  
   Aggregating matched experiment contrasts to study-level means ($\bar{\Delta}_{sj}$) before performing the Student's t-test successfully accounts for within-study correlation across multi-experiment datasets, preventing inflated degrees of freedom and ensuring valid empirical Type I error rates ($\hat{\alpha} = 0.0460 \le 0.05$).

3. **Reproducibility & CLI Tooling:**  
   The calibration battery is integrated into the build pipeline via `just calibrate-composition` and `mtp.composition_calibration.composition`. All random number streams adhere strictly to seed 182 with deterministic stream isolation (`complete_null_study_flip`).
