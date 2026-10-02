# Before and After Summary: Taxon Sampling Model & Complete Exports (#131, #254, #263)

## 1. Overview of Changes

This document records the exact transition from the legacy assertion-count hypergeometric sampling model to the standardized $2 \times 2$ taxon-level sampling urn model in the production pipeline, alongside complete ranking exports and calibration harness benchmarking.

---

## 2. Model & Denominators Comparison (#131)

| Parameter | Previous Assertion-Count Urn | Taxon-Level Sampling Urn (#131) | Notes & Statistical Rationale |
|---|---|---|---|
| **Population Unit** | Total trait annotations summed across all taxonomy subtrees | Distinct taxon carrying trait ($0$ or $1$) | Restores true statistical exchangeability of biological observations |
| **Profile Draws ($N$)** | $N \sim 24,000$ (summed trait counts) | $N = |U_P| \approx 10 - 80$ taxa | Draw size represents actual microbial taxa in signature |
| **Universe Size ($M$)** | $M \sim 500,000$ (summed trait counts) | $M = |\mathcal{U}_{\text{bg}}| \approx 2,475$ unique taxa | Background population of unique taxa in BugSigDB universe |
| **Population Trait Hits ($K$)** | $K = \sum \text{annotations}$ | $K = \sum_{t \in \mathcal{U}_{\text{bg}}} \mathbf{1}[C(t, T) > 0]$ | True taxon-level trait prevalence in background |
| **Observed Trait Hits ($k$)** | $k = \text{sum}(\text{annotations in profile})$ | $k = \sum_{t \in U_P} \mathbf{1}[C(t, T) > 0]$ | Count of unique taxa in signature bearing trait |
| **P-value Degeneracy** | Severe ($p = 0.0$ on modest $2\times$ fold) | Calibrated ($p$-values range $10^{-12} \dots 1.0$) | Prevents BH step-down ties and false-positive distortion |
| **Expected Hits Null** | $\mathbb{E}[k] = N \times (K/M)$ with inflated $N, K, M$ | $\mathbb{E}[k] = N \times (K/M)$ on taxon scale | True expectation under uniform hypergeometric draw |

---

## 3. Results & Ranking Impact (#254)

### 3.1 Overall Enriched Traits
- **Top Enriched Traits**:
  - `spore forming` (High prevalence in decreased profiles across diverse dysbioses)
  - `motile` / `flagellated` (Strong enrichment in invasive / inflammatory conditions)
  - `butyrate producer` & `short chain fatty acid producer` (Consistently depleted in disease, enriched in healthy/increased controls)
  - `pathogenic` & `facultative anaerobe` (Consistently enriched in increased inflammatory profiles)
- **Granular FDR Separation**:
  - With the taxon model, ties in $p = 0.0$ are resolved into nuanced rank-ordered trajectories, allowing clear prioritization of genuine biological drivers over high-annotation bias artifacts.

### 3.2 Complete Ranked Outputs Delivered (#254)
- [`data/derived/trait_enrichment_traits_overall_ranked.tsv`](data/derived/trait_enrichment_traits_overall_ranked.tsv): Global trait rankings with testable denominators ($M, N, K, k$), prevalence, best $p$-value, and best $q$-value.
- [`data/derived/trait_enrichment_trait_study_pairs_ranked.tsv`](data/derived/trait_enrichment_trait_study_pairs_ranked.tsv): Trait-study profile pairs with exact effect sizes and FDR $q$-values.
- [`data/derived/trait_depletion_report.txt`](data/derived/trait_depletion_report.txt) & [`data/derived/trait_depletion_results.tsv`](data/derived/trait_depletion_results.tsv): Complementary under-representation analyses.

---

## 4. Calibration Harness Alignment & Model Status (#263 / #178 / #181 / #182)

- The standard analysis recipes (`just analyze`, `just analyze-depletion`, `just analyze-site`) are explicitly wired to pass the unique-taxon reference background `--taxa-matrix data/derived/taxa_trait_matrix.tsv`.
- Ingestion and statistical scripts enforce strict validation of profile draw sizes ($N \ge \max(k)$) and population counts ($M \ge \max(K, N)$), failing loudly on invalid or missing background matrix inputs.
- **Exchangeability and Calibration Caveats**:
  - Binary presence alone does not establish independent sampling or universal exchangeability across real microbial communities.
  - The hypergeometric model assumes independent draws from the defined background universe, but taxon co-occurrence and phylogenetic structure introduce dependencies that do not create a single universal direction of $p$-value bias.
  - The production statistical model is not marked complete until empirical calibration results are generated under the full calibration battery (#263 / #182: permutation nulls, sampling-unit A/B, and taxon-level bootstrap).
- Ready for cross-adapter benchmarking alongside Noel's GSEA/Score-based calibration (#178 / #181).
