# Taxa Coverage & Profile Trait Mapping Analysis (BugSigDB vs. KG-Microbe)

This investigation quantifies the coverage of BugSigDB taxa and study profiles across KG-Microbe METPO traits, stratified by taxonomic rank and profile mapping status.

---

## 1. Executive Summary across All Three Analytical Tiers

| Metric | Unique Taxa (`taxa.tsv`) | Literature Signature Mentions (`taxa.tsv`) | Study Trait Profiles (`study_trait_profiles_wide.tsv`) |
| :--- | :---: | :---: | :---: |
| **Total Cohort** | **5,829** unique taxa | **59,486** signature mentions | **7,329** study profiles `(study, exp, dir)` |
| **Mapped, Rich** | **3,241** (55.60%) [$\ge 3$ traits] | **52,455** (88.18%) [$\ge 3$ traits] | **6,477** (88.38%) [$\ge 20$ traits] |
| **Mapped, Sparse** | **333** (5.71%) [$1\text{–}2$ traits] | **1,407** (2.37%) [$1\text{–}2$ traits] | **625** (8.53%) [$1\text{–}19$ traits] |
| **Total Trait-Mapped ($\ge 1$ trait)** | **3,574** (**61.31%**) | **53,862** (**90.55%**) | **7,102** (**96.90%**) |
| **Zero-Trait (0 KG-Microbe traits)** | **2,255** (**38.69%**) | **5,624** (**9.45%**) | **227** (**3.10%**) |

### Confirmation on Zero-Trait Profiles in the Pipeline
- Out of **7,329** total BugSigDB study profiles `(study_id, experiment_id, inc_dec)` in `data/derived/study_trait_profiles_wide.tsv`:
  - **7,102 profiles (96.90%)** map to $\ge 1$ KG-Microbe trait (and contain $\ge 1$ trait-mapped taxon).
  - **227 profiles (3.10%)** have **0 mapped traits** (every taxon reported in those 227 profiles has 0 traits in `taxa_trait_matrix.tsv`).
- **Downstream pipeline handling:**
  - **Guild & compositional shift analyses (`src/mtp/trait_guilds.py` `drop_zero_total_pairs` and `src/mtp/trait_composition.py` `compute_log_ratio_shifts`)**: Explicitly filter out zero-trait profiles prior to statistical testing.
  - **Hypergeometric trait enrichment (`scripts/trait_enrichment.py`)**: Zero-trait profiles yield 0 observed counts for all traits and never produce positive enrichments, though `n_taxa` (`taxa_count`) includes all matrix taxa in the signature.
  - **Clustered study-profile heatmaps (`scripts/create_study_trait_profile_heatmap.py`)**: Defaults to `--min-traits 0` unless `--min-traits 1` is passed on the CLI to exclude the 227 zero-trait profiles.

---

## 2. Taxonomic Rank Breakdown: Unique Taxa, Literature Mentions, and Study Profiles

### 2A. Unique Taxa & Literature Occurrences by Rank (`taxa_coverage_by_rank.tsv`)

| Rank | Total Unique Taxa | Zero-Trait Taxa (%) | Sparse 1–2 (%) | Rich $\ge 3$ (%) | Total Trait-Bearing (%) | Total Occurrences | Trait-Bearing Occurrences (%) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Species** | 3,846 | 1,764 (**45.87%**) | 265 (6.89%) | 1,817 (47.24%) | 2,082 (**54.13%**) | 19,113 | 14,908 (**78.00%**) |
| **Genus** | 1,339 | 321 (23.97%) | 50 (3.73%) | 968 (72.29%) | 1,018 (**76.03%**) | 25,388 | 24,499 (**96.50%**) |
| **Family** | 328 | 87 (26.52%) | 5 (1.52%) | 236 (71.95%) | 241 (**73.48%**) | 7,538 | 7,277 (**96.54%**) |
| **Order** | 164 | 46 (28.05%) | 4 (2.44%) | 114 (69.51%) | 118 (**71.95%**) | 2,836 | 2,719 (**95.87%**) |
| **Class** | 78 | 16 (20.51%) | 1 (1.28%) | 61 (78.21%) | 62 (**79.49%**) | 1,930 | 1,886 (**97.72%**) |
| **Phylum** | 72 | 20 (27.78%) | 8 (11.11%) | 44 (61.11%) | 52 (**72.22%**) | 2,672 | 2,568 (**96.11%**) |
| **Superkingdom** | 2 | 1 (50.00%) | 0 (0.00%) | 1 (50.00%) | 1 (**50.00%**) | 9 | 5 (55.56%) |
| **Total / Pooled** | **5,829** | **2,255 (38.69%)** | **333 (5.71%)** | **3,241 (55.60%)** | **3,574 (61.31%)** | **59,486** | **53,862 (90.55%)** |

### 2B. Study Profile Distribution by Taxonomic Rank (`profile_coverage_by_rank.tsv`)

Each study profile (`N = 7,329` total) can contain taxa from one or multiple taxonomic ranks. For profiles containing at least one taxon of a given rank, the table below shows the breakdown by trait-mapping status of the rank's taxa within that profile:

| Rank | Profiles with Rank (% of 7,329) | All Rank Taxa Mapped (%) | Mixed: Mapped + Zero-Trait (%) | Total Mapped at Rank ($\ge 1$ Mapped) | Only Zero-Trait Taxa at Rank (%) | Median Mapped Traits / Profile (IQR) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Species** | 3,332 (45.46%) | 1,934 (58.04%) | 925 (27.76%) | **2,859 (85.80%)** | 473 (**14.20%**) | 109 (43–290) |
| **Genus** | 5,270 (71.91%) | 4,813 (91.33%) | 353 (6.70%) | **5,166 (98.03%)** | 104 (1.97%) | 164 (70–353) |
| **Family** | 2,684 (36.62%) | 2,498 (93.07%) | 117 (4.36%) | **2,615 (97.43%)** | 69 (2.57%) | 280 (135–601) |
| **Order** | 1,440 (19.65%) | 1,362 (94.58%) | 20 (1.39%) | **1,382 (95.97%)** | 58 (4.03%) | 396 (190–762) |
| **Class** | 1,149 (15.68%) | 1,124 (97.82%) | 0 (0.00%) | **1,124 (97.82%)** | 25 (2.18%) | 479 (229–862) |
| **Phylum** | 1,596 (21.78%) | 1,510 (94.61%) | 40 (2.51%) | **1,550 (97.12%)** | 46 (2.88%) | 636 (250–930) |
| **Superkingdom** | 9 (0.12%) | 5 (55.56%) | 0 (0.00%) | **5 (55.56%)** | 4 (44.44%) | 619 (346–866) |

---

## 3. Histogram 1 & Histogram 2 Key Metrics

### Histogram 1: Frequency of Mapped Traits & Mapped Taxa Across Profiles (`N = 7,102` Mapped Profiles; 227 Zero-Trait Profiles Excluded)
- **Panel A (Frequency of Mapped Traits per Profile):**
  - **Median:** **158 unique mapped traits** per profile (IQR: **60–380**, Mean: **288.7**).
  - **Distribution across 7,102 mapped profiles:**
    - `1–19 traits`: 625 profiles (8.8%)
    - `20–49 traits`: 886 profiles (12.5%)
    - `50–99 traits`: 1,035 profiles (14.6%)
    - `100–199 traits`: 1,755 profiles (24.7%) — modal bin
    - `200–399 traits`: 1,174 profiles (16.5%)
    - `400–799 traits`: 763 profiles (10.7%)
    - `800–1,199 traits`: 672 profiles (9.5%)
    - `>=1,200 traits`: 192 profiles (2.7%)
- **Panel B (Mapped Taxa per Profile, Stacked by Taxonomic Rank):**
  - **Median:** **4 trait-mapped taxa** per profile (IQR: **2–9**, Mean: **7.5**).
  - **Distribution across 7,102 mapped profiles:**
    - `1 mapped taxon`: 1,432 profiles (20.2%)
    - `2 mapped taxa`: 944 profiles (13.3%)
    - `3–4 mapped taxa`: 1,282 profiles (18.1%)
    - `5–7 mapped taxa`: 1,171 profiles (16.5%)
    - `8–11 mapped taxa`: 909 profiles (12.8%)
    - `12–19 mapped taxa`: 771 profiles (10.9%)
    - `20–34 mapped taxa`: 442 profiles (6.2%)
    - `>=35 mapped taxa`: 151 profiles (2.1%)

### Histogram 2: Number of Traits from Taxa Profiles & Frequency of Zero-Trait Taxa Across Profiles
- **Panel A (Number of Traits per Taxon, Stacked by Rank, `N = 5,829` Unique Taxa):**
  - **2,255 taxa (38.7%)** have **0 traits** (1,764 Species [78.2% of zero-trait taxa], 321 Genus [14.2%], 170 Family+ [7.5%]).
  - **3,574 taxa (61.3%)** have $\ge 1$ trait: median across all 5,829 taxa is **14 traits** (IQR: 0–42); median among the 3,574 trait-mapped taxa is **32 traits** (IQR: 18–60).
- **Panel B (Frequency of Zero-Trait Taxa per Profile Across All `N = 7,329` Profiles):**
  - **5,494 profiles (75.0%)** contain **0 zero-trait taxa** (100% of their reported taxa map to KG-Microbe traits).
  - **1,835 profiles (25.0%)** contain $\ge 1$ zero-trait taxon:
    - **1,608 profiles (21.9%)** are **mixed profiles** that also contain mapped taxa and are **retained** in trait analyses.
    - **227 profiles (3.1%)** contain **only zero-trait taxa** and are **excluded** when filtering to trait-mapped profiles (`1 zero-trait taxon`: 143 excluded; `2`: 40 excluded; `3–4`: 27 excluded; `5–9`: 12 excluded; `>=10`: 5 excluded).
- **Panel C (Profile Recurrence of the `N = 2,255` Zero-Trait Taxa, Stacked by Rank):**
  - **1,281 zero-trait taxa (56.8%)** appear in **only 1 profile** (singletons).
  - **457 zero-trait taxa (20.3%)** appear in **2 profiles** (77.1% appear in $\le 2$ profiles).
  - Only **26 zero-trait taxa (1.2%)** recur in $\ge 20$ profiles (diagnosable via `top_zero_trait_bugsigdb_taxa.tsv`, led by `Escherichia/Shigella sp.`, `Oscillospiraceae` secondary TaxID `541000`, and `Candida albicans`).

---

## 4. Generated Artifacts & Manifest

### Figures (`.png` at 300 DPI and vector `.svg`):
1. `taxa_coverage_histogram_rank_stratified.png` / `.svg`:
   - **Panel A**: Unique taxa counts per taxonomic rank stacked by trait richness tier with both counts and mapped percentages annotated above each bar.
   - **Panel B**: 100% normalized taxa coverage proportions per rank with counts and percentages.
   - **Panel C**: Distribution of study profiles (`N = 7,329`) containing each taxonomic rank, stacked by mapping status (All Rank Taxa Mapped, Mixed Mapped + Zero-Trait, Only Zero-Trait Taxa).
2. `taxa_coverage_occurrences_vs_unique.png` / `.svg`:
   - Three-tier comparison across Unique Taxa (`N = 5,829`, 61.3% mapped), Literature Signature Mentions (`N = 59,486`, 90.6% mapped), and Study Trait Profiles (`N = 7,329`, 96.9% mapped).
3. `taxa_profile_distribution_by_rank.png` / `.svg`:
   - **Panel A**: 100% stacked composition of study profiles within each taxonomic rank showing exact profile counts and percentages for All-Mapped, Mixed, and Zero-Trait-Only profiles.
   - **Panel B**: Log-scale boxplot + strip distribution of mapped trait richness across `N = 7,102` mapped study profiles stratified by taxonomic rank present.
4. `histogram1_mapped_traits_and_taxa_per_profile.png` / `.svg`:
   - Filtered strictly to mapped profiles (`N = 7,102`, excluding 227 zero-trait profiles):
   - **Panel A**: Frequency of mapped traits across profiles (counts + percentages + median/IQR callout).
   - **Panel B**: Frequency of trait-mapped taxa per profile stacked by taxonomic rank composition.
5. `histogram2_taxa_trait_counts_and_zero_trait_frequency.png` / `.svg`:
   - **Panel A**: Number of traits from taxa profiles (`N = 5,829` unique taxa) stacked by rank, with the 0-trait spike (`2,255 taxa, 38.7%`) highlighted.
   - **Panel B**: Frequency of zero-trait taxa per profile across all `N = 7,329` study profiles, distinguishing retained mixed profiles (`1,608`) from excluded zero-trait-only profiles (`227`).
   - **Panel C**: Profile recurrence distribution of the `2,255` zero-trait taxa across the study cohort, stacked by taxonomic rank.
6. `taxa_trait_count_distribution.png` / `.svg`:
   - Two-panel summary of taxon-level trait counts (with counts and percentages) and rank-stratified log10 trait count boxplots.

### Data Tables (TSVs):
1. `taxa_coverage_by_rank.tsv`: Rank-stratified summary of unique taxa and signature occurrences.
2. `profile_coverage_by_rank.tsv`: Rank-stratified summary of study profile distributions, mapping status counts/percentages, and profile trait richness quantiles.
3. `profile_trait_mapping_summary.tsv`: Profile-level table for all 7,329 study profiles (`study_id`, `experiment_id`, `inc_dec`, total taxa, mapped taxa, zero-trait taxa, mapped trait count, and per-rank taxon breakdowns).
4. `taxa_kgmicrobe_coverage_details.tsv`: Complete audit table for all 5,829 unique taxa.
5. `top_zero_trait_bugsigdb_taxa.tsv`: Top 100 zero-trait BugSigDB taxa ranked by signature occurrences with root-cause category and remediation path.

---

## 5. How to Reproduce

Execute the standalone generator script:

```bash
uv run python scripts/create_taxa_coverage_histogram.py
```
