# Taxa Coverage & Trait Harmonization Analysis (BugSigDB vs. KG-Microbe)

This investigation implements **Approach A** ("Rank-Stratified Stacked Histogram") and in-depth informatics auditing to quantify the coverage of BugSigDB taxa in KG-Microbe METPO traits.

---

## 1. Executive Summary & Key Findings

| Metric | Total Unique Taxa (BugSigDB) | Total Literature Occurrences (BugSigDB Signatures) |
| :--- | :--- | :--- |
| **Total Cohort** | **5,829** unique taxa | **59,486** signature mentions |
| **Mapped, Rich ($\ge 3$ traits)** | **3,241** (55.60%) | **52,455** (88.18%) |
| **Mapped, Sparse ($1\text{–}2$ traits)** | **333** (5.71%) | **1,407** (2.37%) |
| **Total Trait-Bearing ($\ge 1$ traits)** | **3,574** (61.31%) | **53,862** (**90.55%**) |
| **Mapped, 0 Traits (No KG traits)** | **2,255** (**38.69%**) | **5,624** (**9.45%**) |
| **Unmapped (Raw BugSigDB Signatures)** | < 1% of raw signature rows (60/7,425) | 0% in curated `taxa.tsv` (all 5,829 have NCBI TaxIDs) |

### The "Long-Tail" Coverage Paradox:
1. **At the Catalog Level**: 38.69% of unique taxa in BugSigDB have zero traits in KG-Microbe.
2. **At the Real-World Literature Level**: **90.55% of all BugSigDB signature mentions ARE trait-bearing**! Only 9.45% of published literature findings lack phenotypic traits.
3. **Scientific Explanation**: Published microbiome studies repeatedly identify common, well-characterized gut and clinical microbes (*Bacteroides*, *Prevotella*, *Faecalibacterium*, *Bifidobacterium*, *Clostridium*, *Escherichia*), which have dozens to hundreds of phenotypic annotations in KG-Microbe. The 2,255 zero-trait taxa reside predominantly in the "long tail" of uncultivated, rare, or study-specific metagenomic bins.

---

## 2. Approach A: Rank-Stratified Coverage Breakdown

Taxonomic coverage exhibits a steep drop-off at the **Species** rank, while higher taxonomic ranks maintain robust coverage (>72% to >79% trait-bearing for taxa, and >96% for literature occurrences).

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

---

## 3. Informatics Root-Cause Diagnosis of Zero-Trait Taxa

Systematic classification of the 2,255 zero-trait taxa (and specifically the top 100 most frequent in BugSigDB signatures) revealed 5 distinct causes:

```
┌────────────────────────────────────────────────────────────────────────┐
│             Informatics Causes of Zero-Trait Taxa in BugSigDB          │
├────────────────────────────────┬───────────────────┬──────────────────┤
│ Category                       │ Top Examples      │ Remediation Path │
├────────────────────────────────┼───────────────────┼──────────────────┤
│ 1. NCBI TaxID Desynchronization│ Oscillospiraceae  │ Remap alias to   │
│    (Canonical twin has traits; │ (541000 vs 216572)│ canonical TaxID  │
│     BugSigDB used alias)       │ Lactobacillaceae  │ (Immediate win!) │
│                                │ (81850 vs 33958)  │                  │
│                                │ Rothia            │                  │
│                                │ (508215 vs 32207) │                  │
├────────────────────────────────┼───────────────────┼──────────────────┤
│ 2. Ambiguous Bins / Slashes    │ Escherichia/      │ Resolve to LCA   │
│    (Multi-genus composite)     │ Shigella sp.      │ or composite     │
│                                │ (139 occurrences) │ trait vector     │
├────────────────────────────────┼───────────────────┼──────────────────┤
│ 3. Uncultivated MAGs / Bins    │ Clostridiales sp.,│ Impute / roll up │
│    (~60% of zero-trait taxa)   │ UCG-005, CAG:352, │ higher-rank      │
│                                │ UBA1819           │ family traits    │
├────────────────────────────────┼───────────────────┼──────────────────┤
│ 4. Non-Bacterial / Eukaryotes  │ Ascomycota,       │ Integrate fungal │
│    (BacDive lacks micro-       │ Candida albicans, │ trait ontologies │
│     eukaryotic phenotypes)     │ Malassezia        │ (e.g. FungiTraits)│
├────────────────────────────────┼───────────────────┼──────────────────┤
│ 5. Candidatus Lineages         │ Ca. Gastranaero-  │ Metabolic model  │
│    (Candidate uncultured taxa) │ philales, Ca.     │ inference from   │
│                                │ Altimarinota      │ GTDB / GEMs      │
└────────────────────────────────┴───────────────────┴──────────────────┘
```

---

## 4. Generated Artifacts & Manifest

### Figures:
1. `taxa_coverage_histogram_rank_stratified.png` / `.svg`:
   - Two-panel Approach A figure.
   - **Panel A**: Absolute taxa counts per taxonomic rank stacked by trait richness tier (Rich $\ge 3$, Sparse 1–2, Zero traits). Total $N$ labeled above bars.
   - **Panel B**: 100% normalized proportions per rank with inline percentage labels, parity reference line at 50%, and annotations highlighting the species drop-off (45.9% zero) vs. genus/family stability (~75% trait-bearing).
2. `taxa_coverage_occurrences_vs_unique.png` / `.svg`:
   - Two-bar comparative chart directly contrasting Unique Taxa (61.3% trait-bearing) vs. Literature Signature Occurrences (90.55% trait-bearing), illustrating the real-world coverage of microbiome literature.
3. `taxa_trait_count_distribution.png` / `.svg`:
   - **Panel A**: Binned trait count distribution per taxon showing the prominent spike at zero traits (2,255 taxa) alongside the continuous distribution of trait counts for trait-bearing taxa ($1\text{–}1,196$).
   - **Panel B**: Log-scale boxplot of trait counts across taxonomic ranks showing the trait rollup effect in the KG taxonomy.

### Data Tables (TSVs):
1. `taxa_coverage_by_rank.tsv`:
   - Quantitative summary table by rank with counts, percentages, and occurrence weights.
2. `taxa_kgmicrobe_coverage_details.tsv`:
   - Complete audit table for all 5,829 unique taxa with taxon ID, name, rank, trait count, occurrence count, study count, and coverage category.
3. `top_zero_trait_bugsigdb_taxa.tsv`:
   - Top 100 zero-trait BugSigDB taxa ranked by signature occurrences, enriched with diagnosed informatics category and concrete remediation recommendation.

---

## 5. Strategic Recommendations for BugSigDB & KG-Microbe Harmonization

1. **Immediate Quick Win (NCBI Secondary TaxID Remapping)**:
   - Remap `NCBITaxon:541000 Oscillospiraceae` $\rightarrow$ `NCBITaxon:216572` (instantly gains 162 traits across 82 signature occurrences).
   - Remap `NCBITaxon:81850 Lactobacillaceae` $\rightarrow$ `NCBITaxon:33958` (instantly gains 166 traits across 19 signature occurrences).
   - Remap `NCBITaxon:508215 Rothia` $\rightarrow$ `NCBITaxon:32207` (instantly gains 87 traits across 41 signature occurrences).
   - *Total immediate gain*: 142 literature signature occurrences re-connected to rich trait profiles without altering experimental data.

2. **Lowest Common Ancestor (LCA) Mapping for Slashes**:
   - `NCBITaxon:1940338 Escherichia/Shigella sp.` (139 occurrences across 87 studies) cannot map cleanly to a single species in BacDive.
   - Map to family `Enterobacteriaceae` or compute the union of traits shared by *Escherichia* and *Shigella*.

3. **Parental Trait Imputation for Uncultured Metagenomic Taxa**:
   - For taxa like `uncultured Oscillospiraceae bacterium` or `Ruminococcaceae bacterium UCG-005`, propagate essential high-confidence family-level traits (e.g., Gram-stain, oxygen tolerance, endospore formation).

4. **Fungal Trait Integration**:
   - Integrate fungal phenotypic trait databases (e.g. FungiTraits) into KG-Microbe to provide coverage for clinically relevant mycobiome taxa (*Candida albicans*, *Malassezia restricta*, *Ascomycota*).

---

## 6. How to Reproduce

Execute the standalone generator script:

```bash
uv run python scripts/create_taxa_coverage_histogram.py
```
