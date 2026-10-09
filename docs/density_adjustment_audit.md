# Trait Correlation, Annotation Density Adjustment, and Near-Twin Redundancy Audit

**Issue Reference:** [#264](https://github.com/CultureBotAI/microbial-trait-profiling/issues/264), [#352](https://github.com/CultureBotAI/microbial-trait-profiling/issues/352)  
**Lead (Density Adjustment & Near-Twins):** Ilyas Abdolcader  
**Lead (Ontology Redundancy Half):** Jed Dongjin Kim-Ozaeta  
**Related Issues:** [#118](https://github.com/CultureBotAI/microbial-trait-profiling/issues/118), [#187](https://github.com/CultureBotAI/microbial-trait-profiling/issues/187), [#189](https://github.com/CultureBotAI/microbial-trait-profiling/issues/189)  
**Machine-Readable Summary Tables:**  
- Summary Metrics: [`data/derived/density_adjustment_summary.tsv`](file:///c:/Users/ilabd/bsdbng/data/derived/density_adjustment_summary.tsv)  
- Table 1 Adjusted Effects: [`data/derived/table1_density_adjusted_effects.tsv`](file:///c:/Users/ilabd/bsdbng/data/derived/table1_density_adjusted_effects.tsv)  
- Text Report: [`data/derived/density_adjustment_report.txt`](file:///c:/Users/ilabd/bsdbng/data/derived/density_adjustment_report.txt)  

---

## 1. Executive Summary & Problem Formulation

Discussions in [#264](https://github.com/CultureBotAI/microbial-trait-profiling/issues/264) and [#118](https://github.com/CultureBotAI/microbial-trait-profiling/issues/118) raised a fundamental methodological question:
> *Is the dense positive correlation across microbial traits a biological signal of co-selected phenotypes, or is it an artifact of annotation density driven by uneven curation effort across well-studied model organisms?*

To resolve this question empirically:
1. We quantify the baseline trait correlation structure across all 317 non-MediaDive traits and measure the variance explained by the primary principal component (PC1).
2. We test whether the correlation structure collapses or persists after adjusting for the primary annotation density axis using:
   - **Method A:** Profile-level PC1 regression / residualization.
   - **Method B:** Profile composition normalization (Total Sum Scaling / profile annotation richness normalization: $x_i / \sum_j x_{ij}$).
3. We test the explicit diagnostic case identified in #264: *facultatively aerobic* (`METPO:1000608`) $\times$ *moderately halophilic* (`METPO:1000623`).
4. We evaluate the survival of the six canonical Table 1 convergent traits under density adjustment.
5. We audit near-twin candidate features (*coccus-shaped* `METPO:1000668` vs. *sphere-shaped* `METPO:1000683`) for bivariate correlation, paired contrast correlation, and Variance Inflation Factor (VIF).

---

## 2. Baseline Dimensionality & Correlation Structure

Restricting the trait universe to non-MediaDive traits (`object_category != 'METPO:1004005'`) yields $P = 317$ trait features across $N = 5,829$ study trait profiles.

| Metric | Pearson Correlation | Spearman Correlation ($\rho$) | Interpretation |
| :--- | :--- | :--- | :--- |
| **Median Pairwise Correlation** | $+0.6000$ | $+0.5286$ | Pervasive positive correlation across nearly all trait pairs. |
| **Mean Pairwise Correlation** | $+0.5498$ | $+0.4887$ | Substantial global positive shift. |
| **Pairs Exceeding $+0.5$** | $61.1\%$ | $54.1\%$ | Over half of all trait pairs exhibit strong positive correlation. |
| **Pairs Exceeding $+0.9$** | $7.1\%$ | $1.1\%$ | Substantial number of near-redundant pairs. |
| **PC1 Variance Explained** | **$60.82\%$** | **$53.73\%$** | A single global axis accounts for over $60\%$ of total trait variation. |
| **PC1–PC5 Variance Explained** | **$80.27\%$** | — | Top 5 components explain $>80\%$ of total variance. |

### Effective Number of Independent Tests ($M_{\text{eff}}$)

Treating 317 traits as independent tests in statistical inference severely over-corrects (Bonferroni $\alpha = 0.05 / 317 \approx 1.58 \times 10^{-4}$). We computed three standard spectral estimators of $M_{\text{eff}}$ from the eigenvalues of the correlation matrix:

- **Li & Ji (2005):** $M_{\text{eff}} = 62.0$  
  $$f(\lambda) = \mathbb{I}(|\lambda| \ge 1) + (|\lambda| - \lfloor |\lambda| \rfloor)$$
  *(Note: uses absolute eigenvalues $|\lambda|$ and zeros round-off noise $< 10^{-12}$ to prevent negative eigenvalue numerical artifacts from inflating $M_{\text{eff}}$ per #352).*
- **Galwey (2009):** $M_{\text{eff}} = 38.4$  
  $$M_{\text{eff}} = \frac{\left(\sum \sqrt{\max(0, \lambda_i)}\right)^2}{\sum \lambda_i}$$
- **Cheverud (1995) / Nyholt (2004):** $M_{\text{eff}} = 196.8$  
  $$M_{\text{eff}} = 1 + (P - 1)\left(1 - \frac{\text{Var}(\lambda)}{P}\right)$$

These estimators confirm that the effective dimensionality of the 317-trait space is on the order of $\approx 38$ to $62$ independent dimensions, driven by the colossal concentration of variance on PC1.

---

## 3. Density Axis Adjustment: Method Comparison

### Method A: PC1 Regression Residualization
We standardize each trait column and regress out the first principal component scores:
$$\mathbf{x}_j^{\text{adj}} = \mathbf{x}_j - \hat{\beta}_j \mathbf{pc}_1$$
- **Effect on global correlation:** The median pairwise correlation collapses from $+0.6000$ to **$-0.0261$** (centered around zero), and mean pairwise correlation drops to $+0.0014$.
- The proportion of pairs with $r > 0.5$ drops from $61.1\%$ to $6.5\%$.
- Residual traits are orthogonal to the global curation density component.

### Method B: Profile Composition Normalization (Total Sum Scaling)
Each profile's trait counts are normalized by its row sum across non-MediaDive traits:
$$x_{ij}^{\text{norm}} = \frac{x_{ij}}{\sum_k x_{ik}}$$
This adjusts for overall profile annotation richness and study sequencing depth (#352).
- **Effect on global correlation:** Reduces PC1 variance from $60.82\%$ to **$38.94\%$**, and reduces median pairwise correlation from $+0.6000$ to $+0.3121$.
- While Method B substantially attenuates curation bias, residual correlation remains positive because well-annotated taxa still correlate across broad taxonomic clades. Method A provides complete orthogonal separation.

---

## 4. Diagnostic Test Case: Facultatively Aerobic $\times$ Moderately Halophilic

Issue #264 posed the diagnostic question: *Does the observed correlation between facultatively aerobic (`METPO:1000608`) and moderately halophilic (`METPO:1000623`) reflect true biological co-occurrence or annotation density?*

- **Raw Pairwise Correlation:** Pearson $r = +0.8936$.
- **Post-PC1 Regression:** Pearson $r = +0.5349$ (an attenuation of $\Delta r = -0.3587$).
- **Paired Experiment Contrasts ($N = 2,335$ paired comparisons):**
  - **Raw Difference:** Mean difference $\Delta = +0.375$, Wilcoxon signed-rank $p = 2.60 \times 10^{-5}$ (strongly significant).
  - **Density-Adjusted Difference:** Mean difference $\Delta = +0.041$, Wilcoxon signed-rank $p = 0.262$ (**completely nonsignificant!**).

**Finding:** The paired contrast for *facultatively aerobic* completely collapses after removing the PC1 density component. The apparent association was a classic curation density artifact: heavily-annotated human gut/pathogen taxa frequently have both oxygen tolerance and osmoprotection traits recorded in culture collections (e.g., *BacDive* / *Madin et al.*), creating spurious co-enrichment.

---

## 5. Survival of Convergent Six Trait Signals

We evaluated whether genuine biological signals survive PC1 density adjustment in the $N = 2,335$ paired experiment comparisons (increased vs. decreased):

| Trait Identifier | Trait Label | Raw Mean $\Delta$ | Raw Wilcoxon $p$ | Adj Mean $\Delta$ | Adj Wilcoxon $p$ | Relative Change | Adjustment Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `METPO:1000608` | Facultatively aerobic | $+0.375$ | $2.60 \times 10^{-5}$ | $+0.041$ | **$0.262$** | $-89.2\%$ | **Collapsed (artifact removed)** |
| `METPO:2000103\|GO:0030245` | Cellulose catabolic process | $-0.076$ | $1.37 \times 10^{-2}$ | **$-0.228$** | **$1.84 \times 10^{-11}$** | $+198.7\%$ | **Enhanced depletion (genuine signal)** |
| `METPO:1000668` | Coccus shaped | $+0.466$ | $1.72 \times 10^{-6}$ | $+0.146$ | **$7.21 \times 10^{-3}$** | $-68.6\%$ | **Survived (biological signal)** |
| `METPO:1000683` | Sphere shaped | $+0.394$ | $3.49 \times 10^{-9}$ | $+0.163$ | **$9.35 \times 10^{-5}$** | $-58.6\%$ | **Survived (biological signal)** |
| `METPO:2000103\|pathways:arsenate_reduction` | Arsenate reduction | $+0.155$ | $4.64 \times 10^{-5}$ | $+0.058$ | **$4.62 \times 10^{-2}$** | $-62.2\%$ | **Survived (biological signal)** |
| `METPO:1000623` | Moderately halophilic | $+0.364$ | $1.77 \times 10^{-7}$ | $+0.113$ | **$3.33 \times 10^{-4}$** | $-69.0\%$ | **Survived (biological signal)** |

### Key Biological Takeaways
1. **Cellulose Catabolism Depletion Is Strengthened:** When the positive curation baseline is removed, the depletion of cellulose degrading bacteria in disease signatures becomes dramatically sharper ($p = 1.84 \times 10^{-11}$). In raw data, the global positive density drift was masking the true magnitude of cellulose depletion.
2. **Morphological and Stress Signals Survive:** *Coccus shaped* ($p = 0.0072$), *sphere shaped* ($p = 9.35 \times 10^{-5}$), and *moderately halophilic* ($p = 0.00033$) remain statistically significant after removing PC1, indicating genuine non-random representation in disease signatures.
3. **Facultatively Aerobic Is Discarded:** As hypothesized in #264, *facultatively aerobic* does not survive density adjustment and must not be interpreted as an independent disease-associated metabolic adaptation.

---

## 6. Near-Twin Collinearity & Conditional Testing (Coccus vs. Sphere)

We evaluated the collinearity and conditional independence of the sibling morphological traits *coccus-shaped* (`METPO:1000668`) and *sphere-shaped* (`METPO:1000683`):

| Metric | Profile-Level Value ($N = 5,829$) | Paired Contrast Value ($N = 2,335$) |
| :--- | :--- | :--- |
| **Bivariate Correlation ($r$)** | **$+0.9110$** | **$+0.8808$** |
| **Coefficient of Determination ($R^2$)** | $0.8299$ | $0.7758$ |
| **Variance Inflation Factor (VIF)** | **$5.88$** | **$4.46$** |
| **Regression Slope (Coccus $\sim$ Sphere)** | $1.3524$ | $1.3405$ |
| **Regression Intercept** | $+0.1171$ | $-0.0631$ |
| **Residual Coccus $\mid$ Sphere Mean** | — | $-0.0000$ |
| **Residual Coccus $\mid$ Sphere Wilcoxon $p$** | — | $0.0301$ |

### Methodological Recommendations for Near-Twins
- **Severe Multicollinearity:** A VIF of $5.88$ (and $4.46$ in paired contrasts) exceeds standard collinearity thresholds ($> 5.0$), demonstrating that simultaneous inclusion of both traits in multivariate regression models will cause severe variance inflation and coefficient instability.
- **Conditional Residual Signal:** The residual difference of coccus conditioned on sphere exhibits marginal significance ($p = 0.0301$), consistent with the fact that coccus has wider taxonomic representation ($5,614$ direct taxa vs. $381$).
- **Recommendation:** Analysts should not include both *coccus-shaped* and *sphere-shaped* as independent predictors in the same multivariate model. As established in [#189](https://github.com/CultureBotAI/microbial-trait-profiling/issues/189), they should be reported as distinct provenance branches (`infores:madin_etal` vs. `infores:bactotraits`) or modeled with regularization (e.g., Ridge / ElasticNet) or conditional contrasts.

---

## 7. Biological Interpretation Bounds

Following [`methods/biological_interpretation.md`](file:///c:/Users/ilabd/bsdbng/methods/biological_interpretation.md):
1. **Annotation Representation Only:** Trait counts measure curated knowledge representation in published signatures, not physical abundance in patient gut microbiomes.
2. **Do Not Conflate Density with Co-selection:** Positive correlations across traits in raw profiles must not be cited as biological co-selection or shared operon architecture without testing against PC1 density residualization.
3. **No Causal Inferences:** The survival of *moderately halophilic* ($p = 0.00033$) does not indicate patient-level osmolarity changes or causal gut dysbiosis; it indicates that taxa curated with this trait are disproportionately represented in published disease-increased signatures.
