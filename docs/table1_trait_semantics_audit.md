# Table 1 Trait Semantics, Ontology Overlap, and Source Evidence Audit

**Issue Reference:** [#189](https://github.com/CultureBotAI/bsdbng/issues/189)  
**Leads:** Ilyas Abdolcader, Jed Dongjin Kim-Ozaeta  
**Related Issues:** [#178](https://github.com/CultureBotAI/bsdbng/issues/178), [#179](https://github.com/CultureBotAI/bsdbng/issues/179), [#182](https://github.com/CultureBotAI/bsdbng/issues/182), [#187](https://github.com/CultureBotAI/bsdbng/issues/187)  
**Authoritative Machine-Readable Audit Table:** Generated to `data/derived/table1_trait_semantic_audit.tsv`

---

## 1. Executive Summary & Context

Discussions in [#187](https://github.com/CultureBotAI/bsdbng/issues/187) proposed several analytical adjustments to the six canonical microbial traits highlighted in Table 1:
1. Merging *coccus-shaped* (`METPO:1000668`) and *sphere-shaped* (`METPO:1000683`) into a single morphological feature.
2. Interpreting *arsenate reduction* (`METPO:2000103|pathways:arsenate_reduction`) as *arsC* cytoplasmic detoxification and linking it to heavy-metal exposure or plasmid-borne co-selection.
3. Reading *moderately halophilic* (`METPO:1000623`) as general *halotolerance*.

Biological interpretations cannot be assigned based on display labels. This audit traces each Table 1 feature back to its underlying Knowledge Graph (KG-Microbe) edges, source databases, ontology definitions, and BugSigDB profile distributions before modifying feature naming, groupings, or mechanistic narratives.

---

## 2. Pinned KG-Microbe Provenance

All evidence in this audit is evaluated against the pinned KG-Microbe snapshot recorded in `data/kg-microbe/provenance.json`. Checksums and file metadata match the pinned build exactly:

| File Name | MD5 Checksum | Size (bytes) | Google Drive File ID | Fetch Timestamp |
| :--- | :--- | :--- | :--- | :--- |
| `merged-kg_edges.tsv` | `555ecf5ba5d89288cdc44e368375374f` | 1,034,256,091 | `1iIVd-FEleeRIv7CPS0l1F732p3J9JO5d` | 2026-04-21T00:21:25+00:00 |
| `merged-kg_nodes.tsv` | `dd76b4656f9341db25a7c4afa96cd58e` | 209,545,218 | `1lDVjE67E82LMUdNtqWGvyczrHvNu3BdY` | 2026-04-21T00:21:08+00:00 |

*Provenance rule: Pinned KG inputs must be recovered via the standard provenance workflow; substitute snapshots must be declared as separate scenarios.*

---

## 3. Machine-Readable Feature Audit Summary

The machine-readable table captures all required audit attributes across the 6 features:

| Trait Key | Feature Name | Ontology ID & Label | Primary KG Source (NCBITaxon Edges, Every Scope) | BacDive Strain Edges (Read Under `phenotype`, Not `off`) | Direct KG Taxa | BugSigDB Matrix Taxa (`phenotype`, of 5,829) | Profile Support (`phenotype`, $N=5,629$) | Evidence Status | Decision |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `METPO:1000608` | Facultatively aerobic | `METPO:1000608` (facultatively aerobic) | `infores:madin_etal` (5,496) | 53 | 5,496 | 920 (15.8%) | 4,412 (78.4%) | measured + propagated | **keep_separate** |
| `METPO:2000103\|GO:0030245` | Cellulose catabolic process | `GO:0030245` (cellulose catabolic process) | `infores:madin_etal` (186) | 0 | 186 | 189 (3.2%) | 2,922 (51.9%) | measured + propagated | **keep_separate** |
| `METPO:1000668` | Coccus shaped | `METPO:1000668` (coccus shaped) | `infores:madin_etal` (5,614) | 1,476 | 5,614 | 815 (14.0%) | 4,386 (77.9%) | measured + propagated | **keep_separate** |
| `METPO:1000683` | Sphere shaped | `METPO:1000683` (sphere shaped) | `infores:bactotraits` (381) | 11 | 381 | 307 (5.3%) | 3,490 (62.0%) | measured + propagated | **keep_separate** |
| `METPO:2000103\|pathways:arsenate_reduction` | Arsenate reduction | `pathways:arsenate_reduction` (arsenate_reduction) | `infores:madin_etal` (51) | 0 | 51 | 61 (1.0%) | 1,739 (30.9%) | unresolved mechanism + propagated | **unresolved / keep_separate** |
| `METPO:1000623` | Moderately halophilic | `METPO:1000623` (moderately halophilic) | `infores:madin_etal` (426) | 4 | 426 | 351 (6.0%) | 3,364 (59.8%) | measured (discrete) + unresolved (quant) + propagated | **keep_separate** |

**Source columns ([#774](https://github.com/CultureBotAI/microbial-trait-profiling/issues/774), [#788](https://github.com/CultureBotAI/microbial-trait-profiling/issues/788)).**
- Edges are counted by `primary_knowledge_source` in the pinned edge file.
- **NCBITaxon edges.** The matrix builder (`src/mtp/taxa_trait_matrix.py`) reads them in every strain rollup scope. Each subject carries one edge per trait here, so this column equals *Direct KG Taxa*.
- **BacDive strain edges.** Every `infores:bacdive` edge for these traits sits on a BacDive strain node (`kgmicrobe.strain:`). Each of these strains links to an NCBITaxon node through `biolink:subclass_of`. An earlier version of this table listed these edges as a primary source; they now have their own column.
- **When strain edges count.** They reach the matrix under the primary `phenotype` scope, which counts them at the strains' NCBITaxon nodes. Since #788 it is the default of both the core workflow and the legacy matrix builder. The NCBITaxon nodes carrying each trait then become 5,521, 186, 6,027, 388, 51 and 429. Under the `off` sensitivity scope they are not read.
- **Which corpus each column comes from ([#807](https://github.com/CultureBotAI/microbial-trait-profiling/issues/807)).** Every count was recounted for the final rerun at code commit 511c891:
  - *Primary KG Source*, *BacDive Strain Edges*, *Direct KG Taxa* and the node counts above come from the pinned KG edge file (§2), which the rerun did not change. They were recounted on that file, the node counts with the matrix builder's own edge reader.
  - *BugSigDB Matrix Taxa* counts, among the 5,829 reported taxa of the taxa-trait matrix that the final rerun built under `phenotype`, those that bear the trait themselves or through a descendant; the percentage is of the 5,829. That matrix is not committed. Its sha256 `45df0fe7…` is the pin in `data/analysis/cellulose_decr/cellulose_decr_provenance.json`, and five of the six counts are also its `matrix_taxon_count` in the committed `data/analysis/cellulose_decr/trait_convergence_benchmark.tsv` (arsenate reduction has too few decreased feces taxa to enter that benchmark). Under `off` the six counts are 916, 189, 721, 305, 61 and 350.
  - *Profile Support* counts, among the 5,629 human 16S/WMS profiles with no Tier 1+2, pairing or disease-orientation filter, those with a nonzero count in the committed `data/derived/study_trait_profiles_wide.tsv.gz`, which the final rerun rebuilt under `phenotype` (decompressed sha256 `43b2e230…`). The `off` table it replaced (`519a095b…`) gave 4,411, 2,922, 4,338, 3,490, 1,739 and 3,364.

---

## 4. Deep-Dive Findings by Acceptance Criteria

### 4.1 Coccus vs. Sphere Morphology: Equivalence, Provenance, and Overlap

#### Ontology Structure
- Both `METPO:1000668` (*coccus shaped*) and `METPO:1000683` (*sphere shaped*) are direct subclasses (`biolink:subclass_of`, `rdfs:subClassOf`) of `METPO:1000666` (*cell shape*).
- They are **sibling classes** within METPO, not parent and child.
- `METPO:1000668` defines: *"A cell shape in which an organism has a spherical or nearly spherical morphology, with roughly equal dimensions in all directions."* Synonyms: `coccus`, `coccus-shaped`.
- `METPO:1000683` contains an unpopulated narrative definition in KG nodes. Synonyms: `S_sphere`, `sphere-shaped`.

#### Source Provenance Divergence
The existence of two sibling terms in METPO is an artifact of upstream multi-source extract-transform-load (ETL) harmonization:
- `METPO:1000668` originates from **Madin et al. (2020)** (`infores:madin_etal`), where the primary column `cell_shape` was annotated with `"coccus"`.
- `METPO:1000683` originates from **Bactotraits / BacDive** (`infores:bactotraits`), which utilizes the shape code `S_sphere`.
- BacDive strain records (`infores:bacdive`) also assert both terms (1,476 coccus, 11 sphere), on strain nodes. The primary `phenotype` scope counts them at the strains' NCBITaxon nodes, adding 413 coccus and 7 sphere nodes ([#788](https://github.com/CultureBotAI/microbial-trait-profiling/issues/788)); the `off` matrix does not read them ([#774](https://github.com/CultureBotAI/microbial-trait-profiling/issues/774)). The evidence of the NCBITaxon edges themselves is Madin et al. for coccus and BactoTraits for sphere.

#### Quantified Overlap & Statistical Dependence
1. **Direct KG Taxa Assertions:**
   - Coccus direct taxa: **5,614**
   - Sphere direct taxa: **381**
   - Intersection: **266 taxa** (69.8% of Sphere direct taxa are in Coccus; only 4.7% of Coccus direct taxa are in Sphere).
   - Jaccard similarity: **0.0464** (4.6%).
2. **BugSigDB Matrix Taxa (after taxonomic hierarchy rollup; the final rerun's `phenotype` matrix, §3):**
   - Coccus taxa in matrix ($K_{\text{coccus}}$): **815**
   - Sphere taxa in matrix ($K_{\text{sphere}}$): **307**
   - Intersection: **296 taxa** (96.4% of Sphere matrix taxa are in Coccus; 36.3% of Coccus matrix taxa are in Sphere).
   - Jaccard similarity: **0.3584** (35.8%).
3. **Profile-Level Nonzero Support & Correlation ($N=5,629$ human 16S/WMS profiles, with no Tier 1+2, pairing or disease-orientation filter, of the committed `phenotype` table, §3):**
   - Coccus nonzero profiles: **4,386** (77.9%)
   - Sphere nonzero profiles: **3,490** (62.0%)
   - Pearson correlation across the raw profile counts of all 5,629 profiles: $r = 0.807$ ($p < 10^{-300}$).
   - Spearman rank correlation, same counts: $\rho = 0.915$ ($p < 10^{-300}$).

#### Explicit Grouping Decision: Keep Separate
- **Decision:** **Do NOT merge `METPO:1000668` and `METPO:1000683` in the primary analysis.**
- **Rationale:** While biologically both refer to spherical cell morphology, merging them by naive string collapse or summing profile counts without deduplication would distort the feature space:
  - If summed, the 296 shared BugSigDB taxa would be double-counted in profile vectors, artificially inflating the weight of dual-annotated taxa.
  - Sibling terms in METPO do not carry transitive equivalence. Merging requires constructing a deduplicated union universe ($K_{\text{union}} = 826$ taxa) and formally versioning a new hypothesis family under [#182](https://github.com/CultureBotAI/bsdbng/issues/182).
  - Preserving them separately honors historical hypotheses while explicitly documenting the multi-source annotation artifact.

---

### 4.2 Arsenate Reduction: Cytoplasmic Detoxification vs. Respiratory Reduction

#### Pathway Definitions in Microbial Physiology
Microbial reduction of arsenate ($As^{5+}$) to arsenite ($As^{3+}$) proceeds via two distinct, biologically non-interchangeable mechanisms:
1. **Cytoplasmic Arsenate Detoxification (*ars* operon, *arsC*):**
   - Widespread across bacteria and archaea.
   - Intracellular reduction of toxic arsenate to arsenite mediated by ArsC (glutaredoxin- or thioredoxin-dependent arsenate reductase), followed by active extrusion via ArsB or Acr3 permeases.
   - Does **not** conserve energy or support anaerobic respiration.
2. **Dissimilatory Arsenate Respiration (*arr* operon, *arrAB*):**
   - Periplasmic heterodimeric respiratory complex (*arrA* catalytic subunit, *arrB* electron transfer).
   - Couples the oxidation of organic electron donors (lactate, acetate) or $H_2$ to the reduction of $As^{5+}$ as a terminal electron acceptor in anaerobic respiration, generating a proton-motive force for ATP synthesis.

#### KG-Microbe Source Evidence
- Feature Key: `METPO:2000103|pathways:arsenate_reduction`
- KG-Microbe Object CURIE: `pathways:arsenate_reduction`
- Primary Source: `infores:madin_etal` (51 direct taxon assertions).
- Original Source Column: `pathways.csv` $\to$ `pathway == "arsenate_reduction"`.
- **Evidence Coverage & Unresolved Fraction:**
  - Direct KG edges with mechanistic/gene attribution: **0 / 51 (0.0%)**.
  - **Unresolved fraction: 100%**.
- **Prohibited Inferences:**
  - Do **not** classify this trait as *arsC* detoxification.
  - Do **not** classify this trait as *arrAB* dissimilatory respiration.
  - Do **not** infer specific gene carriage, plasmid carriage, mobile genetic element transfer, patient arsenic exposure, or heavy-metal/antimicrobial co-selection from this taxon-level pathway label.
  - **Decision:** Classified as **unresolved mechanism**; retained under literal label "capable of: arsenate_reduction".

---

### 4.3 Moderately Halophilic: Salt Requirement vs. Tolerance

#### Definition and Ontology Structure
- `METPO:1000623` (*moderately halophilic*) is a subclass of `METPO:1000629` (*halophily preference*).
- Formal METPO definition: *"A halophily preference where growth and proliferation requires high levels of sodium chloride, usually above or about 0.2 M."*
- Synonyms: `moderate-halophilic`, `moderately halophilic`.

#### Physiological Distinction: Halophily vs. Halotolerance
- **Halophiles** exhibit an absolute physiological **requirement** for sodium chloride for structural stability and growth (typically 0.2 M to 3 M NaCl).
- **Halotolerant** microorganisms do not require NaCl for growth, but possess osmoadaptive mechanisms enabling survival and proliferation across elevated osmotic gradients (e.g., *Staphylococcus aureus*).
- METPO defines an explicit **requirement** ($>0.2$ M NaCl, $\approx 1.2\%$ w/v).

#### Missing Quantitative Parameters
- Direct KG edges: 426 taxa from `infores:madin_etal` (qualitative column `range_salinity == "moderate-halophilic"`). Four BacDive strain records (`infores:bacdive`) also assert the term, on strain nodes. The primary `phenotype` scope counts them at the strains' NCBITaxon nodes, adding 3 nodes ([#788](https://github.com/CultureBotAI/microbial-trait-profiling/issues/788)); the `off` matrix does not read them ([#774](https://github.com/CultureBotAI/microbial-trait-profiling/issues/774)).
- Quantitative continuous parameters in KG edges: **None**. The KG records zero values for salt concentration optimum, min/max ranges, ionic composition, or specific growth media.
- **Decision:** Retain original label **"moderately halophilic"**. Reading this trait as "halotolerant" is unsupported by METPO semantics. Missing quantitative culture parameters are classified as **unresolved**.

---

### 4.4 Oxygen Phenotype and Cellulose Degradation Boundaries

#### Facultatively Aerobic (`METPO:1000608`)
- Defined as: *"An oxygen preference in which growth can occur without oxygen but is capable of aerobic growth."*
- **Biological Precision:** Oxygen tolerance class only establishes that growth can occur in the presence or absence of $O_2$. It does **not** establish anaerobic respiration using alternative electron acceptors such as nitrate ($NO_3^-$) or tetrathionate ($S_4O_6^{2-}$).
- As established by [Rojas-Tapias et al. (2022)](https://www.nature.com/articles/s41564-022-01224-7), anaerobic nitrate or tetrathionate respiration depends on specific enzymatic complexes (*nar*, *nap*, *ttr*) and host mucosal inflammatory conditions generating oxidants. Conflating facultative aerobiosis with mucosal nitrate respiration is an unsupported mechanistic leap.

#### Cellulose Catabolic Process (`METPO:2000103|GO:0030245`)
- Defined as: *"The chemical reactions and pathways resulting in the breakdown of cellulose, a linear beta1-4 glucan..."*
- **Biological Precision:**
  1. A curated pathway capability annotation does **not** prove demonstrated in vivo degradation of dietary cellulose in human cohorts.
  2. Cellulose catabolism must **not** be equated with butyrate or short-chain fatty acid (SCFA) production. While some cellulolytic organisms in the rumen or human colon produce SCFAs, primary cellulolytic microbes often ferment cellulose into acetate, succinate, lactate, ethanol, $H_2$, and $CO_2$ (e.g., *Fibrobacter succinogenes* produces succinate and acetate, not butyrate; *Ruminococcus champanellensis* produces acetate and succinate).

---

## 5. Decisions, Approved Wording, and Unsupported Interpretations

| Feature Key | Grouping Decision | Approved Phrasing | Unsupported Interpretations |
| :--- | :--- | :--- | :--- |
| `METPO:1000608` | **Keep separate** | "Taxa annotated as facultatively aerobic (capable of growth in both the presence and absence of oxygen)." | Do not infer anaerobic nitrate/tetrathionate respiration, *nar*/*nap*/*ttr* expression, or host mucosal inflammation. |
| `METPO:2000103\|GO:0030245` | **Keep separate** | "Taxa annotated with pathway capability for cellulose catabolic processes." | Do not infer active in vivo cellulose breakdown in patients. Do not equate cellulose catabolism with butyrate or SCFA production. |
| `METPO:1000668` | **Keep separate** | "Taxa annotated with coccus-shaped cell morphology (primarily derived from Madin et al. isolate compilation)." | Do not merge with sphere-shaped by naive string replacement. Do not claim biological non-equivalence without noting source database divergence. |
| `METPO:1000683` | **Keep separate** | "Taxa annotated with sphere-shaped cell morphology (derived from BactoTraits records and, under the primary `phenotype` scope, BacDive strain records)." | Do not treat as biologically distinct morphology from coccus. Do not sum profile counts without deduplicating shared taxa. |
| `METPO:2000103\|pathways:arsenate_reduction` | **Unresolved mechanism / Keep separate label** | "Taxa annotated with pathway capability for arsenate reduction (derived from Madin et al. compilation; mechanism unresolved between dissimilatory respiration and cytoplasmic detoxification)." | Do not classify as *arsC* detoxification. Do not classify as *arrAB* respiration. Do not infer gene presence, plasmid carriage, arsenic exposure, or co-selection. |
| `METPO:1000623` | **Keep separate label / Unresolved quantitative conditions** | "Taxa annotated as moderately halophilic (defined in METPO as requiring high levels of sodium chloride, typically $\ge 0.2$ M NaCl)." | Do not read as "halotolerant" (halophily requires salt; halotolerance does not). Do not infer continuous salt optima, ranges, or media. |

---

## 6. Downstream Handoff Contracts

### 6.1 Handoff to Issue #178 (Integrated Evidence Table Contract)
Issue [#178](https://github.com/CultureBotAI/bsdbng/issues/178) defines the shared result contract for multi-arm enrichment integration. This audit supplies the following metadata requirements for Table 1 features:
- **Raw trait identifiers:** Must preserve qualified keys (`METPO:1000608`, `METPO:2000103|GO:0030245`, `METPO:1000668`, `METPO:1000683`, `METPO:2000103|pathways:arsenate_reduction`, `METPO:1000623`).
- **Input fingerprinting:** Pinned MD5 checksums for nodes (`dd76b4656f9341db25a7c4afa96cd58e`) and edges (`555ecf5ba5d89288cdc44e368375374f`).
- **Feature universe denominators ($K$):** Coccus ($K=815$) and Sphere ($K=307$) in the `phenotype` matrix (721 and 305 under `off`) must be evaluated as separate testing units unless a deduplicated union universe is formally constructed.
- **Evidence status encoding:** Carry `measured`, `propagated`, and `unresolved` flags into result tables to distinguish direct isolate assertions from taxonomic rollups.

### 6.2 Handoff to Issue #179 (Project-Wide Biological Interpretation)
Issue [#179](https://github.com/CultureBotAI/bsdbng/issues/179) governs biological interpretation across the manuscript, figures, and automated reports. The boundaries established here provide canonical examples:
- *Evidence chain:* Source literature $\to$ curated isolate traits $\to$ NCBI taxonomy rollup $\to$ BugSigDB signature overlap $\to$ hypergeometric statistic.
- *Defensible claims:* A significant shift indicates a higher proportion of taxa in a differential abundance signature possessing the recorded trait annotation.
- *Prohibited claims:* Statistical enrichment in a differential signature does not demonstrate in vivo metabolic fluxes, gene expression levels, environmental chemical exposure, or host causation.
