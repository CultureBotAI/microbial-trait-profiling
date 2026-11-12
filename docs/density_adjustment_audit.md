# Trait Correlation, Annotation Density Adjustment, and Near-Twin Redundancy Audit

**Issue Reference:** [#264](https://github.com/CultureBotAI/microbial-trait-profiling/issues/264), [#352](https://github.com/CultureBotAI/microbial-trait-profiling/issues/352)  
**Lead (Density Adjustment & Near-Twins):** Ilyas Abdolcader  
**Lead (Ontology Redundancy Half):** Jed Dongjin Kim-Ozaeta  
**Related Issues:** [#118](https://github.com/CultureBotAI/microbial-trait-profiling/issues/118), [#179](https://github.com/CultureBotAI/microbial-trait-profiling/issues/179), [#187](https://github.com/CultureBotAI/microbial-trait-profiling/issues/187), [#189](https://github.com/CultureBotAI/microbial-trait-profiling/issues/189), [#325](https://github.com/CultureBotAI/microbial-trait-profiling/issues/325), [#357](https://github.com/CultureBotAI/microbial-trait-profiling/issues/357), [#375](https://github.com/CultureBotAI/microbial-trait-profiling/issues/375), [#381](https://github.com/CultureBotAI/microbial-trait-profiling/issues/381), [#382](https://github.com/CultureBotAI/microbial-trait-profiling/issues/382)  
**Machine-Readable Summary Tables:**  
- Summary Metrics: `data/derived/density_adjustment/density_adjustment_summary.tsv`  
- Table 1 Adjusted Effects: `data/derived/density_adjustment/table1_density_adjusted_effects.tsv`  
- Text Report: `data/derived/density_adjustment/density_adjustment_report.txt`  
- Provenance JSON: `data/derived/density_adjustment/density_adjustment_provenance.json`  
- Ontology Redundancy Pairs: `data/derived/density_adjustment/ontology_redundancy_pairs.tsv`  
- Ontology Redundancy Summary: `data/derived/density_adjustment/ontology_redundancy_summary.tsv`  
- Ontology Collapsing Comparison: `data/derived/density_adjustment/ontology_collapsing_comparison.tsv`  

Regenerate all of the above with `just audit-trait-density`. Only the three ontology tables are committed (CONTRIBUTING.md, "Generated outputs on `main`").

**Reporting scope (#1449).** The numerical examples below retain the historical
`off`-scope audit; they are not current-run family sizes or a phenotype-scope
validation. The rounded nominal p-values printed in section 5 are the source
for the documentation-only BH reconstruction by
[`issue1321_doc_table_bh_q.py`](../scripts/adhoc/issue1321_doc_table_bh_q.py).
Its raw and adjusted columns are separate families,
`mt-script-doc-table-bh-raw` and `mt-script-doc-table-bh-adjusted`: every one of
the six printed trait rows belongs to each column's family, with BH and
q ≤ 0.05. These are the families of the reconstructed q-values quoted in
sections 4, 5 and 8, not full-precision q-values from a new producer execution.
The earlier density run's complete producer binding is not supplied on this
page; recover and verify it before reusing its numerical examples as run-bound
results. The [MT policy](../methods/multiple_testing_policy.md) and
[registry](../methods/multiple_testing_families.yaml) govern new reporting.

`just density-adjustment-trait-list` reports the paired effects for a declared list of trait identifiers instead of the six Table 1 traits (#1054). PC1 is still fit on every non-MediaDive trait, so a trait's means and Wilcoxon p-values are the same whichever list it appears in, while its BH q-values and status are taken over the list it is reported in (#1321); each listed trait must be a non-MediaDive column, and labels come from the columns table's `object_name`. The six-trait default keeps its fixed labels, so a shared trait can carry a different label in the two tables while its means and p-values are identical: arsenate reduction is `arsenate_reduction` in the current columns table (#1084). The run writes `trait_list_provenance.json` with the list file's SHA-256, the input hashes and the q-value family (`q_values`). `just audit-trait-density --trait-list FILE --output-dir DIR` passes the same list through the jed arm's script; its consolidated report computes section 6 from the same reported effects and statuses, including declared lists (#1069). A list run refuses the default `data/derived/density_adjustment/` and any directory with a `density_adjustment` path component, including arm directories, nested destinations and symlink aliases (#1093). Component comparisons use `casefold()`, so spellings such as `Density_Adjustment` are also protected (#1144). A destination already containing `table1_density_adjusted_effects.tsv` is reserved regardless of the entry contents or type; use a fresh destination for each declared-list run. The audit API checks the destination before profile-input reads or directory creation; the module CLI also checks before loading the trait list. An empty `--output-dir` is rejected rather than read as the current directory (#1084). Without a list, default runs write `density_adjustment_provenance.json` with the same input-hash and q-value family fields (without the list-file fields, #1327).

### Consolidated report contract (#1069)

The report section 3 corresponds to this document's section 5: both contain
paired trait effects. The report section 6 contains conclusions; this document's
section 6 remains the separate near-twin discussion. Their section numbers are
not interchangeable.

For every reported trait set, report section 6 repeats each reported label,
status and raw/adjusted BH q-value from section 3. Its survivor list includes only
`survived (significant post-adjustment)` and
`enhanced (significant, larger magnitude)`, in report order. It prints `none`
when no reported row has either status. Adjusted-only, reversed, undefined
magnitude and nonsignificant rows do not enter that survivor list. Declared lists
use the same rule; the paired-effect/status summary in report section 6, point 2,
introduces no trait absent from section 3. Global near-twin diagnostics remain
separate in point 4 and can describe traits outside the reported effects list.
The survivor list describes annotation contrasts, not measured trait activity
or causal effects.

Applied to the historical statuses retained in section 5, the summary is:

Historical surviving labels: Facultatively aerobic; Cellulose catabolic process; Coccus shaped; Sphere shaped; Arsenate reduction; Moderately halophilic.

This is a consistency summary of the historical table, not a new audit run.
Its numerical cells and the near-twin diagnostics are unchanged. The reporting
scope above explains the documentation-only BH reconstruction and missing
historical producer binding.

The CLI's `--code-commit` is passed explicitly to the consolidated report writer
and printed as `code_commit: VALUE`. This records the caller-supplied identifier;
it does not prove a clean checkout or verify the identifier against Git. When
omitted, the report prints:

`Provenance note: code_commit not supplied; producer commit is unverified.`

No Git discovery or dirty-working-tree check is required. This unbound note
also applies when rendering historical results whose producer commit was not
recorded; it does not turn this page's examples into verified current outputs.

> [!NOTE]
> **Step 5 Decomposition Boundary (#355, #382):** This document and PR #370 cover trait correlation, annotation density adjustment, and near-twin collinearity. The trait-semantics audit (#311) is tracked separately in the ontology redundancy workstream.

---

## 1. Executive Summary & Problem Formulation

Discussions in [#264](https://github.com/CultureBotAI/microbial-trait-profiling/issues/264) and [#118](https://github.com/CultureBotAI/microbial-trait-profiling/issues/118) raised a fundamental methodological question:
> *Is the dense positive correlation across microbial traits a biological signal of co-selected phenotypes, or is it an artifact of annotation density driven by uneven curation effort across well-studied model organisms?*

To characterize the associations without identifying their cause:
1. We quantify the baseline trait correlation structure across all non-MediaDive traits and measure the variance explained by the primary principal component (PC1).
2. We test whether the correlation structure collapses or persists after adjusting for the primary annotation density axis using:
   - **Method A:** Profile-level PC1 regression / residualization.
   - **Method B:** Profile composition normalization (Total Sum Scaling / profile annotation richness normalization: $x_i / \sum_j x_{ij}$).
3. We test the explicit diagnostic case identified in #264: *facultatively aerobic* (`METPO:1000608`) $\times$ *moderately halophilic* (`METPO:1000623`).
4. We evaluate the paired contrasts of the six selected Table 1 traits under density adjustment.
5. We audit near-twin candidate features (*coccus-shaped* `METPO:1000668` vs. *sphere-shaped* `METPO:1000683`) for bivariate correlation, partial correlation controlling for study annotation depth, paired contrast correlation, and Variance Inflation Factor (VIF).

---

## 2. Baseline Dimensionality & Correlation Structure

Restricting the trait universe to non-MediaDive traits (`object_category != 'METPO:1004005'`) yields $P = 317$ trait features under the `off` strain-rollup scope this audit used (344 under the primary `phenotype` scope, not audited here) across $N = 5,629$ study trait profiles ($0$ constant columns excluded). Note that exact correlation figures depend on the profile representation (raw count sums vs. binarized presence/absence per #365); the audit evaluates the wide profile representations generated by the pipeline.

| Metric | Pearson Correlation | Spearman Correlation ($\rho$) | Interpretation |
| :--- | :--- | :--- | :--- |
| **Median Pairwise Correlation** | $+0.7049$ | $+0.6105$ | Pervasive positive correlation across trait pairs. |
| **Mean Pairwise Correlation** | $+0.5946$ | $+0.5471$ | Substantial global positive shift. |
| **Pairs Exceeding $+0.5$** | 67.7% | 63.4% | Majority of trait pairs exhibit strong positive correlation. |
| **Pairs Exceeding $+0.9$** | 17.8% | 3.5% | Substantial number of near-redundant pairs in raw counts. |
| **PC1 Variance Explained** | **66.50%** | **60.42%** | A single global axis accounts for over 66% of total trait variation. |
| **PC1–PC5 Variance Explained** | **89.35%** | **78.81%** | Top 5 components explain $\approx 89$% of total variance. |

### Effective-Test-Count Sensitivity Diagnostics ($M_{\text{eff}}$)

Trait features are correlated. Bonferroni ($\alpha = 0.05 / P$) does not require
independence and can be conservative under dependence; the correlations alone
do not quantify that conservatism. We evaluate spectral estimators of
$M_{\text{eff}}$ as sensitivity diagnostics, not replacement family sizes:

- **Li & Ji (2005):** Canonical estimator in `mtp.trait_screen.effective_family_size`  
  $$f(\lambda) = \mathbb{I}(|\lambda| \ge 1) + (|\lambda| - \lfloor |\lambda| \rfloor)$$
  *(Uses absolute eigenvalues $|\lambda|$ and excludes constant columns to prevent numerical roundoff artifacts from inflating $M_{\text{eff}}$ per #352, #377, #378).*
  **$M_{\text{eff}} = 30.0$**
- **Galwey (2009):** Alternative spectral estimator  
  $$M_{\text{eff}} = \frac{\left(\sum \sqrt{\max(0, \lambda_i)}\right)^2}{\sum \lambda_i}$$
  **$M_{\text{eff}} = 16.0$**
- **Cheverud (1995) / Nyholt (2004):** Excluded per #325 and #357 ($M_{\text{eff}} = 172.6$); the documented concern is sensitivity to collinear eigenvalue distributions and failure of the audited boundary constraints.

The historical estimates for the 317-trait `off` matrix range from
$\approx 16$ to $30$. They summarize its eigenvalue spectrum and do not establish
that this many hypotheses are independent or that an effective-M correction
controls FDR. Under policy R9, `mt-core-noel-sign-effective` is a sensitivity
bound with no FDR guarantee; the full-catalog BY family
`mt-core-noel-sign-nominal` remains primary. Its size comes from the producing
screen's declared candidate inventory, not either spectral estimate.

---

## 3. Density Axis Adjustment: Method Comparison

### Method A: PC1 Regression Residualization
We standardize each trait column and regress out the first principal component scores:
$$\mathbf{x}_j^{\text{adj}} = \mathbf{x}_j - \hat{\beta}_j \mathbf{pc}_1$$
- **Effect on global correlation:** The median pairwise correlation collapses from $+0.7049$ to **$-0.0268$** (centered around zero), and mean pairwise correlation drops to $+0.0024$.
- The proportion of pairs with $r > 0.5$ drops from $67.7\%$ to $16.9\%$.
- Residual traits are orthogonal to the fitted PC1 scores. Identifying PC1 with curation effort is an interpretation, not a causal result of PCA.

### Method B: Profile Composition Normalization (Total Sum Scaling)
Each profile's trait counts are normalized by its row sum across non-MediaDive traits:
$$x_{ij}^{\text{norm}} = \frac{x_{ij}}{\sum_k x_{ik}}$$
This rescales overall profile annotation totals (#352, #379); it does not identify the biological or curation causes of those totals.
- **Effect on diagnostic pair correlation:** Reduces Pearson correlation between *facultatively aerobic* and *moderately halophilic* from $+0.8957$ to **$+0.1651$** (median pairwise $r = +0.0165$, pairs $> 0.5 = 2.4\%$).
- Method B reduces the displayed correlation under a compositional representation. Method A makes the fitted residuals orthogonal to PC1. Neither calculation alone shows that curation bias has been removed.

---

## 4. Diagnostic Test Case: Facultatively Aerobic $\times$ Moderately Halophilic

Issue #264 posed the diagnostic question: *Does the observed correlation between facultatively aerobic (`METPO:1000608`) and moderately halophilic (`METPO:1000623`) reflect true biological co-occurrence or annotation density?*

- **Raw Pairwise Correlation:** Pearson $r = +0.8957$.
- **Post-PC1 Regression:** Pearson $r = +0.6256$ (an attenuation of $\Delta r = -0.2701$).
- **Post-Profile-Composition:** Pearson $r = +0.1651$ (an attenuation of $\Delta r = -0.7306$).
- **Paired Experiment Contrasts ($N = 2,257$):**
  - **Raw Difference:** Nominal Wilcoxon signed-rank $p = 8.05 \times 10^{-13}$, BH $q = 4.83 \times 10^{-12}$ in `mt-script-doc-table-bh-raw`, the six printed raw-p rows (significant at $q \le 0.05$, #1321; documentation reconstruction described above).
  - **Density-Adjusted Difference:** Nominal Wilcoxon signed-rank $p = 2.97 \times 10^{-12}$, BH $q = 1.78 \times 10^{-11}$ in `mt-script-doc-table-bh-adjusted`, the six printed adjusted-p rows (passes $q \le 0.05$; status `survived (significant post-adjustment)`).

**Finding:** The correlation between *facultatively aerobic* and *moderately halophilic* decreases by $\Delta r = -0.27$ under PC1 regression and drops to $+0.17$ under profile composition normalization. This describes sensitivity to the transformation, not a causal attribution to curation. The raw-count paired contrast passes the documentation-family BH threshold above. The earlier binarized-profile claim in #375 was read from nominal p-values; it is not a corrected-family replication of this result.

---

## 5. Paired Contrasts of the Six Selected Traits

The historical table records paired experiment comparisons before and after PC1
adjustment (increased vs. decreased, $N = 2,257$ pairs). Its twelve displayed
Wilcoxon p-values are nominal. The statuses are retained from that historical
table; their documentation-only BH comparison uses the two separate six-row
families defined above, not a joint family of all displayed tests:

| Trait Identifier | Trait Label | Raw Mean $\Delta$ | Raw nominal Wilcoxon $p$ | Adj Mean $\Delta$ | Adj nominal Wilcoxon $p$ | Relative Change | Magnitude Change | Adjustment Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `METPO:1000608` | Facultatively aerobic | $+356.65$ | $8.05 \times 10^{-13}$ | $+204.56$ | **$2.97 \times 10^{-12}$** | -42.6% | -42.6% | **Survived (significant post-adjustment)** |
| `METPO:2000103\|GO:0030245` | Cellulose catabolic process | $-1.78$ | $4.98 \times 10^{-2}$ | **$-4.85$** | **$5.61 \times 10^{-6}$** | -172.5% | +172.5% | **Enhanced (significant, larger magnitude)** |
| `METPO:1000668` | Coccus shaped | $+210.05$ | $2.54 \times 10^{-8}$ | $+89.64$ | **$1.32 \times 10^{-3}$** | -57.3% | -57.3% | **Survived (significant post-adjustment)** |
| `METPO:1000683` | Sphere shaped | $+9.18$ | $3.80 \times 10^{-7}$ | $+2.26$ | **$1.11 \times 10^{-3}$** | -75.4% | -75.4% | **Survived (significant post-adjustment)** |
| `METPO:2000103\|pathways:arsenate_reduction` | Arsenate reduction | $+2.42$ | $2.55 \times 10^{-8}$ | $+1.32$ | **$1.01 \times 10^{-4}$** | -45.5% | -45.5% | **Survived (significant post-adjustment)** |
| `METPO:1000623` | Moderately halophilic | $+12.38$ | $1.29 \times 10^{-6}$ | $+3.93$ | **$6.47 \times 10^{-9}$** | -68.3% | -68.3% | **Survived (significant post-adjustment)** |

*Relative Change* is `relative_change_pct`, the signed change $(\text{adj} - \text{raw}) / |\text{raw}|$. It reads as attenuation only for a positive raw mean: cellulose's raw mean is negative, so its $-172.5\%$ means the effect grew. *Magnitude Change* is `magnitude_change_pct`, $(|\text{adj}| - |\text{raw}|) / |\text{raw}|$, negative exactly when the effect shrank whatever its sign; `sign_flipped` marks an adjusted mean of the opposite sign (#1084). A mean within $10^{-9}$ of zero has no sign, and a raw one leaves the magnitude change undefined (NaN). The Magnitude Change values above are computed from the same raw and adjusted means. The status combines the Benjamini-Hochberg $q$-values of the Wilcoxon $p$-values (significant at $q \le 0.05$, #1321) with the magnitude change:
- A raw mean of zero gets *undefined magnitude change (raw mean zero, …)*, significant or nonsignificant post-adjustment.
- A sign reversal gets *reversed (significant post-adjustment)* or *reversed (nonsignificant post-adjustment)*, whatever its magnitude (#1097).
- Otherwise a significant adjusted contrast whose raw contrast was also significant ($q \le 0.05$) is *survived*, or *enhanced* if $|\text{adj}| > |\text{raw}|$, while one whose raw contrast was nonsignificant ($q > 0.05$) or had ten or fewer non-zero pairs (no raw $p$ or $q$) is *significant only post-adjustment* (#1098); one that lost significance is *collapsed*, attenuated or with a larger magnitude; one nonsignificant before and after is *attenuated* or *enhanced magnitude*.

**Runtime correction families (#1321).** `mt-supp-density-raw-wilcoxon` and
`mt-supp-density-adjusted-wilcoxon` each use BH across every reported trait, raw
and adjusted columns separately, at $q \le 0.05$. They are distinct from the
documentation-reconstruction families above. The reported list is the six
Table 1 traits by default, or the declared list of a
`density-adjustment-trait-list` run. A trait with ten or fewer non-zero pairs
stays in the family, held at $p = 1$ internally, and has no exported $q$.
The family size is the number of reported traits, including these untestable
members. Every requested trait must have a non-MediaDive column; an absent
Table 1 or list member fails the run (#1324).

For a quoted runtime q, read `bh_family_size` beside `raw_wilcoxon_bh_q` or
`adj_wilcoxon_bh_q` in `table1_density_adjusted_effects.tsv` and bind that table
to `density_adjustment_provenance.json` or `trait_list_provenance.json`, whose
`q_values` records the family and procedure. Preserve its input hashes, list
and threshold. Do not substitute the historical six-row size for a different
list run. Because the family is the reported set, a trait can retain its means
and nominal p-values across lists while its q and status change. BH is the
registry's deliberate choice for these density sensitivity families; it does
not confer a full-catalog BY result or measure trait activity.

The table above predates #1321: its statuses were read from nominal $p \le 0.05$. Recomputed separately with BH over each column's six documented $p$-values (`scripts/adhoc/issue1321_doc_table_bh_q.py`), every status is unchanged. The closest is cellulose catabolic process, whose raw $p = 4.98 \times 10^{-2}$ is the largest of the six and so equals its $q$; every other raw and adjusted $q$ is at most $1.32 \times 10^{-3}$.

A raw contrast with ten or fewer non-zero pairs has no raw $p$ and is treated as not raw-significant: it becomes *significant only post-adjustment* when the adjusted contrast is significant, and *attenuated* or *enhanced magnitude* when the adjusted contrast is nonsignificant. Deck figure 4 draws only *survived* and *enhanced* rows as surviving (green); *significant only post-adjustment*, reversed, and raw-mean-zero rows are drawn red with the rest.

### Key Methodological Takeaways
The $q$-values below are BH over the six traits of the table above, recomputed from its $p$-values by `scripts/adhoc/issue1321_doc_table_bh_q.py` (#1326): raw q belongs to `mt-script-doc-table-bh-raw`, adjusted q to `mt-script-doc-table-bh-adjusted`, each at $q \le 0.05$. Each takeaway describes the paired contrast (increased minus decreased profile of the same experiment) of rolled-up trait annotation counts in raw count profiles, before and after PC1 adjustment; none measures trait activity in the sampled communities.
1. **The cellulose catabolic process contrast grows under PC1 adjustment:** The contrast is negative: profiles of taxa reported decreased carry more *cellulose catabolic process* annotations than profiles of taxa reported increased. Its mean moves from $-1.78$ (raw $q = 4.98 \times 10^{-2}$, just inside $q \le 0.05$) to $-4.85$ after adjustment (adjusted $q = 1.12 \times 10^{-5}$; magnitude change $+172.5\%$). The fitted PC1 contribution partly offsets this negative contrast; that decomposition does not identify the source of PC1.
2. **The coccus shaped and sphere shaped contrasts shrink and stay at $q \le 0.05$:** *Coccus shaped* moves from raw $q = 5.10 \times 10^{-8}$ to adjusted $q = 1.32 \times 10^{-3}$ (magnitude change $-57.3\%$), and *sphere shaped* from raw $q = 5.70 \times 10^{-7}$ to adjusted $q = 1.32 \times 10^{-3}$ ($-75.4\%$). #375 reported that neither passed its nominal p-value threshold under binarized presence/absence profiles; that statement does not report a corrected family.
3. **The facultatively aerobic, moderately halophilic and arsenate reduction contrasts shrink and stay at $q \le 0.05$:** *Facultatively aerobic* moves from raw $q = 4.83 \times 10^{-12}$ to adjusted $q = 1.78 \times 10^{-11}$ ($-42.6\%$), *moderately halophilic* from $1.55 \times 10^{-6}$ to $1.94 \times 10^{-8}$ ($-68.3\%$), and *arsenate reduction* from $5.10 \times 10^{-8}$ to $1.51 \times 10^{-4}$ ($-45.5\%$). These describe how reported taxa align with oxygen, salt and arsenate annotations, not measured oxygen tolerance, salinity or arsenate reduction in the sampled communities.
4. **Representation Sensitivity (#375):** Whether a contrast keeps its label (e.g., *facultatively aerobic* and *arsenate reduction*) depends on profile representation (raw count sums vs. binarized presence/absence per #365 opt-in flags). Pipelines should report both representations.

---

## 6. Near-Twin Collinearity & Conditional Testing (Coccus vs. Sphere)

We evaluated correlation and regression diagnostics for the sibling morphological traits *coccus-shaped* (`METPO:1000668`) and *sphere-shaped* (`METPO:1000683`); these diagnostics do not establish conditional independence:

- **Profile Correlation:** $r = +0.7970$
- **Profile Composition Normalized Correlation:** $r = +0.2309$ (PC1 residualized $r = +0.7124$)
- **Variance Inflation Factor (VIF):** $2.74$ (paired contrast VIF = $2.34$, paired $r = +0.7573$)
- **Profile Regression Slope (Coccus ~ Sphere):** $\beta = 20.41$
- **Paired Regression Slope:** $\beta = 22.69$ (nominal diagnostic $p < 10^{-15}$; one displayed regression test, no family correction claimed)

### Methodological Recommendations for Near-Twins
- **Collinearity Diagnostics:** The reported $r \approx 0.80$ and VIF of $2.74$ describe association and variance inflation in the fitted design. No severity cutoff or model-instability result is established here.
- **Conditioning:** The traits retain positive association ($r = +0.2309$ under composition normalization, $+0.7124$ after PC1 regression). Persistence after these transformations does not identify an ontological or biological cause.
- **Annotation Disparity:** The regression slope ($\beta \approx 22.7$) describes a scale difference; the separate source audit (#189) identifies `infores:madin_etal` and `infores:bactotraits`. The slope alone does not attribute that difference to source curation.
- **Recommendation:** Evaluate the two terms' source overlap and the intended model before choosing separate, regularized or conditional predictors. A correlation or VIF alone does not justify merging their biological definitions.

---

## 7. METPO Ontology Redundancy Audit & Branch Collapsing Rule

**Issue Reference:** [#264](https://github.com/CultureBotAI/microbial-trait-profiling/issues/264) (Task 2)  
**Lead:** Jed Dongjin Kim-Ozaeta  
**Related Semantic Audits:** [#189](https://github.com/CultureBotAI/microbial-trait-profiling/issues/189) (trait semantics & provenance branches), [#131](https://github.com/CultureBotAI/microbial-trait-profiling/issues/131), [#179](https://github.com/CultureBotAI/microbial-trait-profiling/issues/179)  

### 7.1 Problem Formulation: Ontology Redundancy vs. Annotation Density

Issue #264 raised whether the pervasive correlation structure across traits is an artifact of ontology hierarchy (ancestor-descendant subsumption and synonym duplication) or curation effort (well-studied taxa carrying numerous curated labels at once).

To separate these two hypotheses:
1. **Audit Ontology Redundancy:** Identify ancestor-descendant hierarchies, chemical/lexical synonyms, and provenance near-twin branches across all 317 non-MediaDive traits of the `off` scope.
2. **Formulate Explicit Branch Collapsing Rules:** Define which term survives and the mathematical semantics of its collapsed count.
3. **Recompute Trait Space Dimensionality:** Evaluate whether collapsing redundant branches attenuates the dominant correlation axis and PC1 variance, or if PC1 density residualization remains necessary.

---

### 7.2 Redundancy Audit Results: 24 Redundant Pairs Across 10 Canonical Clusters

Auditing the $P = 317$ non-MediaDive traits of the `off` scope across $N = 5,629$ study trait profiles identified **24 redundant pairs** structured into **10 canonical concept clusters**:

| Cluster ID | Concept Label | Redundancy Tier | Specific Members & Subsumption | Pairwise $r$ Range | Surviving Representative Trait |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `motility` | Bacterial Motility | Ancestor-Descendant | `METPO:1000702` (*motile*) subsumes `METPO:1000704` (*flagellated*), `METPO:1000705` (*axially filamented*), `METPO:1000706` (*gliding*) | $+0.180$ to $+0.964$ | `METPO:1000702` (*motile*, support $= 3,369,467$) |
| `obligate_anaerobe` | Obligate Anaerobiosis | Ancestor-Descendant | `METPO:1000607` (*obligately anaerobic*) subsumes `METPO:1000611` (*strictly anaerobic*) | $+0.543$ | `METPO:1000607` (*obligately anaerobic*, support $= 229,367$) |
| `dextrin_utilization` | Dextrin Utilization | Ancestor-Descendant | `CHEBI:28675` (*dextrin*) subsumes `CHEBI:40585` (*alpha-cyclodextrin*) | $+0.910$ | `CHEBI:28675` (*dextrin*, support $= 157,198$) |
| `butan2ol_utilization` | Butan-2-ol Utilization | Synonym / Subclass | `CHEBI:35687` (*butan-2-ol*) vs. `CHEBI:45475` (*(2S)-butan-2-ol*) | $+1.000$ | `CHEBI:35687` (*butan-2-ol*, support $= 93$) |
| `citrate_utilization` | Citrate Utilization | Chemical Synonyms | `CHEBI:16947` (*citrate*) vs. `CHEBI:133748` (*citrate anion*) vs. `CHEBI:50744` (*citrate salt*) | $+1.000$ | `CHEBI:16947` (*citrate*, support $= 3,219$) |
| `glutamate_utilization` | Glutamate Utilization | Chemical Synonyms | `CHEBI:18237` (*glutamic acid*) vs. `CHEBI:16015` (*L-glutamic acid*) vs. `CHEBI:29987` (*glutamate*) vs. `CHEBI:14321` (*glutamate(1-)*) | $+1.000$ | `CHEBI:18237` (*glutamic acid*, support $= 1,004$) |
| `malonate_utilization` | Malonate Utilization | Chemical Synonyms | `CHEBI:15792` (*malonate*) vs. `CHEBI:38083` (*malonate ester*) | $+1.000$ | `CHEBI:15792` (*malonate*, support $= 780$) |
| `co_oxidation_pathway` | CO Oxidation Pathway | Lexical Duplicates | `pathways:carbonmonoxide_oxidation` vs. `pathways:carbon_monoxide_oxidation` | $+0.683$ | `METPO:2000103\|pathways:carbon_monoxide_oxidation` (*carbon_monoxide_oxidation*, support $= 2,585$) |
| `coccus_sphere_shape` | Coccus / Sphere Morphology | Sibling Near-Twins | `METPO:1000668` (*coccus shaped*) vs. `METPO:1000683` (*sphere shaped*) | $+0.797$ | `METPO:1000668` (*coccus shaped*, support $= 4,430,794$) |
| `bacillus_rod_shape` | Bacillus / Rod Morphology | Sibling Near-Twins | `METPO:1000667` (*bacillus shaped*) vs. `METPO:1000681` (*rod shaped*) vs. `METPO:1000688` (*coccobacillus shaped*) | $+0.908$ to $+0.931$ | `METPO:1000667` (*bacillus shaped*, support $= 6,884,150$) |

*Note:* In an expanded unconstrained pairwise scan (e.g., matching across discretized temperature, pH, and salinity bins), pairwise redundancy can reach 111 pairs. The 10 canonical clusters above capture all structurally distinct non-MediaDive redundancies within the core metabolic and morphological hierarchy.

---

### 7.3 Formal Stated Branch Collapsing Rules

In accordance with #264's directive that branch collapse requires a stated rule rather than assuming terms are identical measurements, we formalize:

#### Rule 1: Survivor Selection (Max-Support Evidence Rule)
For each redundant cluster $\mathcal{C}_k = \{m_1, m_2, \dots, m_K\}$, the surviving representative feature $s \in \mathcal{C}_k$ is chosen as the feature maximizing total annotation support across all profiles:
$$s = \arg\max_{m \in \mathcal{C}_k} \sum_{i=1}^N \mathbf{X}_{im}$$
*Rationale:* As shown in #189, terms like *coccus* ($4,430,794$ total profile count support) and *sphere* ($200,030$ total profile count support) differ by more than twenty-fold in total evidence. The max-support rule chooses the term with the larger observed total; it does not establish improved stability, power or biological equivalence.

#### Rule 2: Collapsed Count Semantics
We evaluated two distinct count collapse semantics:
1. **Max Member Count (Lower Bound on the Union, Recommended):**
   $$C_{ik} = \max_{m \in \mathcal{C}_k} \mathbf{X}_{im}$$
   *Interpretation:* The count of the best-annotated member. The number of taxa annotated with any trait in $\mathcal{C}_k$ is at least this (the union contains its largest member) and at most $\sum_m \mathbf{X}_{im}$; the max equals the union only when members' taxa are nested. Summing double-counts taxa co-annotated with several members, so the max is the conservative choice.
2. **Native Surviving Representative Count:**
   $$C_{ik} = \mathbf{X}_{is}$$
   *Interpretation:* Directly takes the raw count of the surviving highest-evidence feature $s$, ignoring the subordinate features.

---

### 7.4 Empirical Impact on Dimensionality & Correlation Structure

We recomputed the complete correlation matrix, principal components, and effective test numbers ($M_{\text{eff}}$) across all $N = 5,629$ profiles under raw and collapsed configurations:

| Configuration | Trait Space ($P$) | Median Pairwise $r$ | Mean Pairwise $r$ | Pairs $>0.5$ | Pairs $>0.9$ | PC1 Variance Explained | PC1–PC5 Variance | $M_{\text{eff}}^{\text{Galwey}}$ | $M_{\text{eff}}^{\text{Li&Ji}}$ | Test Pair $r$ (*facultative* $\times$ *halophilic*) | Near-Twin $r$ (*coccus* $\times$ *sphere*) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Raw Uncollapsed** | 317 (`off`) | $+0.7049$ | $+0.5946$ | $67.7\%$ | $17.8\%$ | **$66.50\%$** | $89.35\%$ | $16.0$ | $30.0$ | $+0.8957$ | $+0.7970$ |
| **Collapsed (Native Survivor)** | 301 | $+0.7087$ | $+0.5973$ | $68.1\%$ | $18.5\%$ | **$66.80\%$** | $89.58\%$ | $15.5$ | $28.0$ | $+0.8957$ | — (collapsed) |
| **Collapsed (Deduplicated Max)** | 301 | $+0.7092$ | $+0.5975$ | $68.1\%$ | $18.5\%$ | **$66.81\%$** | $89.59\%$ | $15.5$ | $28.0$ | $+0.8957$ | — (collapsed) |
| **Collapsed + PC1 Adjusted** | 301 | **$-0.0283$** | $+0.0026$ | $17.4\%$ | $1.8\%$ | **$38.66\%$** | $79.13\%$ | **$28.3$** | **$43.0$** | $+0.6181$ | — (collapsed) |

---

### 7.5 Critical Scientific Finding: Redundancy vs. Density

Comparing the collapsed configurations against the raw baseline shows a consistent empirical pattern:
1. **Correlation Persists Under This Collapse Rule:** Collapsing the 24 redundant pairs shrinks the `off` trait universe from 317 to 301 features, while the median pairwise correlation ($+0.7049 \to +0.7092$) and PC1 variance explained ($66.50\% \to 66.81\%$) change little. This does not establish that the removed variation was noise.
2. **Curation Density Remains a Candidate Explanation:** Persistence after this ontology collapse is consistent with shared annotation density, including the *facultatively aerobic* $\times$ *moderately halophilic* correlation at $r = +0.8957$. It does not exclude biological co-occurrence, shared taxa, other ontology effects or sampling differences, and does not assign a causal share to curation.
3. **Complementary Roles of Collapsing and Density Adjustment:**
   - **Ontology collapsing** changes feature definitions and requires a separately justified semantic and deduplication rule; VIF $= 2.74$ alone does not require it.
   - **PC1 adjustment** removes the fitted first-component contribution and can be reported as a sensitivity analysis. It does not establish that PC1 contains only curation effort.
   These operations address different representations and must not be treated as validated causal corrections.

---

## 8. Biological Interpretation Bounds

Following the project interpretation guidelines in `AGENTS.md` and [#179](https://github.com/CultureBotAI/microbial-trait-profiling/issues/179):
1. **Annotation Representation Only:** Trait counts measure curated knowledge representation in published signatures, not physical abundance in patient microbiomes.
2. **Do Not Conflate Density with Co-selection:** Positive correlations across traits in raw or PC1-residualized profiles do not establish biological co-selection or shared operon architecture.
3. **No Causal Inferences:** The historical *moderately halophilic* contrast (nominal density-adjusted $p = 6.47 \times 10^{-9}$, reconstructed BH $q = 1.94 \times 10^{-8}$ in `mt-script-doc-table-bh-adjusted`, six printed traits, $q \le 0.05$) describes annotation representation among taxa reported increased versus decreased. It does not establish patient-level osmolarity changes, a verified disease-versus-healthy contrast, or causal dysbiosis.
