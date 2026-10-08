# Table 1 Trait Semantics, Ontology Overlap, and Source Evidence Audit

**Issue Reference:** [#189](https://github.com/CultureBotAI/bsdbng/issues/189)  
**Leads:** Ilyas Abdolcader, Jed Dongjin Kim-Ozaeta  
**Related Issues:** [#178](https://github.com/CultureBotAI/bsdbng/issues/178), [#179](https://github.com/CultureBotAI/bsdbng/issues/179), [#182](https://github.com/CultureBotAI/bsdbng/issues/182), [#187](https://github.com/CultureBotAI/bsdbng/issues/187)  
**Machine-Readable Audit Table:** [`data/derived/table1_trait_semantic_audit.tsv`](file:///c:/Users/ilabd/bsdbng/data/derived/table1_trait_semantic_audit.tsv)

---

## 1. Executive Summary & Context

Discussions in [#187](https://github.com/CultureBotAI/bsdbng/issues/187) proposed several analytical adjustments to the six canonical microbial traits highlighted in Table 1:
1. Merging *coccus-shaped* (`METPO:1000668`) and *sphere-shaped* (`METPO:1000683`) into a single morphological feature.
2. Interpreting *arsenate reduction* (`METPO:2000103|pathways:arsenate_reduction`) as *arsC* cytoplasmic detoxification and linking it to heavy-metal exposure or plasmid-borne co-selection.
3. Reading *moderately halophilic* (`METPO:1000623`) as general *halotolerance*.

Biological interpretations cannot be assigned based on display labels. This audit traces each Table 1 feature back to its underlying Knowledge Graph (KG-Microbe) edges, source databases, ontology definitions, and BugSigDB profile distributions before modifying feature naming, groupings, or mechanistic narratives.

---

## 2. Pinned KG-Microbe Provenance

All evidence in this audit is evaluated against the pinned KG-Microbe snapshot recorded in [`data/kg-microbe/provenance.json`](file:///c:/Users/ilabd/bsdbng/data/kg-microbe/provenance.json). Checksums and file metadata match the pinned build exactly:

| File Name | MD5 Checksum | Size (bytes) | Google Drive File ID | Fetch Timestamp |
| :--- | :--- | :--- | :--- | :--- |
| `merged-kg_edges.tsv` | `555ecf5ba5d89288cdc44e368375374f` | 1,034,256,091 | `1iIVd-FEleeRIv7CPS0l1F732p3J9JO5d` | 2026-08-30T18:06:23+00:00 |
| `merged-kg_nodes.tsv` | `dd76b4656f9341db25a7c4afa96cd58e` | 209,545,218 | `1lDVjE67E82LMUdNtqWGvyczrHvNu3BdY` | 2026-08-30T18:05:46+00:00 |

*Provenance rule: Pinned KG inputs must be recovered via the standard provenance workflow; substitute snapshots must be declared as separate scenarios.*

---

## 3. Machine-Readable Feature Audit Summary

The machine-readable table [`data/derived/table1_trait_semantic_audit.tsv`](file:///c:/Users/ilabd/bsdbng/data/derived/table1_trait_semantic_audit.tsv) captures all required audit attributes across the 6 features:

| Trait Key | Feature Name | Ontology ID & Label | Primary KG Source | Direct KG Taxa | BugSigDB Matrix Taxa | Profile Support ($N=5,829$) | Evidence Status | Decision |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `METPO:1000608` | Facultatively aerobic | `METPO:1000608` (facultatively aerobic) | `infores:madin_etal` (5,496); `infores:bacdive` (53) | 5,496 | 861 (16.6%) | 4,491 (77.0%) | measured + propagated | **keep_separate** |
| `METPO:2000103\|GO:0030245` | Cellulose catabolic process | `GO:0030245` (cellulose catabolic process) | `infores:madin_etal` (186) | 186 | 164 (3.2%) | 2,963 (50.8%) | measured + propagated | **keep_separate** |
| `METPO:1000668` | Coccus shaped | `METPO:1000668` (coccus shaped) | `infores:madin_etal` (5,614); `infores:bacdive` (1,476) | 5,614 | 675 (13.0%) | 4,409 (75.6%) | measured + propagated | **keep_separate** |
| `METPO:1000683` | Sphere shaped | `METPO:1000683` (sphere shaped) | `infores:bactotraits` (381); `infores:bacdive` (11) | 381 | 283 (5.5%) | 3,541 (60.8%) | measured + propagated | **keep_separate** |
| `METPO:2000103\|pathways:arsenate_reduction` | Arsenate reduction | `pathways:arsenate_reduction` (arsenate_reduction) | `infores:madin_etal` (51) | 51 | 58 (1.1%) | 1,770 (30.4%) | unresolved mechanism + propagated | **unresolved / keep_separate** |
| `METPO:1000623` | Moderately halophilic | `METPO:1000623` (moderately halophilic) | `infores:madin_etal` (426); `infores:bacdive` (4) | 426 | 337 (6.5%) | 3,424 (58.7%) | measured (discrete) + unresolved (quant) + propagated | **keep_separate** |

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

#### Quantified Overlap & Statistical Dependence
1. **Direct KG Taxa Assertions:**
   - Coccus direct taxa: **5,614**
   - Sphere direct taxa: **381**
   - Intersection: **266 taxa** (69.8% of Sphere direct taxa are in Coccus; only 4.7% of Coccus direct taxa are in Sphere).
   - Jaccard similarity: **0.0464** (4.6%).
2. **BugSigDB Matrix Taxa (after taxonomic hierarchy rollup):**
   - Coccus BugSigDB taxa: **675**
   - Sphere BugSigDB taxa: **283**
   - Intersection: **262 taxa** (92.6% of Sphere taxa in BugSigDB are shared with Coccus; 38.8% of Coccus taxa are shared with Sphere).
   - Coccus-only: 413 taxa; Sphere-only: 21 taxa.
   - Jaccard similarity: **0.3764** (37.6%).
3. **BugSigDB Study Profiles ($N = 5,829$):**
   - Coccus non-zero profiles: **4,409** (75.64%, mean count = 3.15, max = 52).
   - Sphere non-zero profiles: **3,541** (60.75%, mean count = 1.92, max = 34).
   - Numerical correlation: **Pearson $r = 0.9110$**, **Spearman $\rho = 0.8468$**.

#### Preservation of Historical Hypotheses & Collapse Decision
- **Hypothesis A (Independent Features / Source Boundaries):** The terms represent distinct curation streams with different isolate coverage. In BugSigDB, 21 taxa carry Sphere annotations without a corresponding Coccus annotation from Madin et al. Merging them by simple renaming would obscure their distinct provenance and misattribute source evidence.
- **Hypothesis B (Semantic Synonym Artifact):** Biologically, "coccus" and "spherical cell" describe the identical geometric morphology. Their near-perfect profile-level co-variation ($r = 0.911$) confirms they capture the same underlying biological signal in differential abundance studies.
- **Decision:** **Keep separate** in primary Table 1 outputs to preserve exact fidelity to KG-Microbe source records.
- **Prerequisites for Future Union/Collapse:** If a future release harmonizes morphology:
  1. **Unique-assertion deduplication:** The union of direct isolate assertions must be computed *before* taxonomic rollup to avoid double-counting taxa present in both Madin et al. and Bactotraits.
  2. **Feature universe denominator:** The total background universe of morphological features ($K$) must be explicitly recalculated.
  3. **Testing family versioning:** The collapsed feature must be evaluated in a separately versioned testing family (e.g., `morphology_harmonized_v1`) subject to [#182](https://github.com/CultureBotAI/bsdbng/issues/182) before rerunning statistical analyses.

---

### 4.2 Arsenate Reduction: Source Tracing and Mechanistic Boundaries

#### Source Provenance
- `METPO:2000103|pathways:arsenate_reduction` asserts capability (`METPO:2000103`) for `pathways:arsenate_reduction` under relation `RO:0002215` (`capable of`).
- Source: **100% of the 51 direct KG edges** originate from `infores:madin_etal` (derived from FAPROTAX / Bergey's Manual metabolic compilation).

#### Mechanistic Breakdown: Detoxification (*arsC*) vs. Dissimilatory Respiration (*arrAB*)
The Rosalind discussion ([#187](https://github.com/CultureBotAI/bsdbng/issues/187)) suggested interpreting this feature as *arsC*-mediated arsenate detoxification. This interpretation is unsupported by the underlying evidence:
1. **Biological Divergence:** Microbial arsenate reduction proceeds via two fundamentally distinct physiological mechanisms:
   - **Cytoplasmic Detoxification:** Mediated by the ubiquitous small enzyme **ArsC** (e.g., characterized in *Staphylococcus aureus* / *E. coli* by [Zegers et al., 2001](https://www.nature.com/articles/nsb1001-843)). ArsC reduces intracellular $As(V)$ to $As(III)$ for subsequent extrusion by ArsB or Acr3. It operates aerobically or anaerobically and does not conserve energy.
   - **Dissimilatory Anaerobic Respiration:** Mediated by the membrane-bound heterodimeric complex **ArrAB** (respiratory arsenate reductase). ArrAB uses arsenate as an alternative terminal electron acceptor in anaerobic respiration, generating a proton-motive force for ATP synthesis.
2. **Taxonomic Inspection in KG-Microbe:** The 51 direct taxa asserted in KG-Microbe include classic obligate/facultative dissimilatory arsenate respirers (*Chrysiogenes arsenatis* DSM 11915, *Alkalihalobacillus macyae*, *Halarsenatibacter silvermanii*, *Mesobacillus selenatarsenatis*, *Seleniivibrio woodruffii*, *Pyrobaculum ferrireducens*), as well as general environmental isolates carrying *arsC* (*Kocuria rosea*, *Arthrobacter agilis*) and dual-pathway organisms (*Herminiimonas arsenicoxydans*).
3. **Evidence Coverage & Unresolved Fraction:**
   - Direct KG edges with mechanistic/gene attribution: **0 / 51 (0.0%)**.
   - **Unresolved fraction: 100%**.
4. **Prohibited Inferences:**
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
- Direct KG edges: 426 taxa from `infores:madin_etal` (qualitative column `range_salinity == "moderate-halophilic"`), 4 strains from `infores:bacdive`.
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
| `METPO:1000683` | **Keep separate** | "Taxa annotated with sphere-shaped cell morphology (derived from Bactotraits and BacDive morphology records)." | Do not treat as biologically distinct morphology from coccus. Do not sum profile counts without deduplicating shared taxa. |
| `METPO:2000103\|pathways:arsenate_reduction` | **Unresolved mechanism / Keep separate label** | "Taxa annotated with pathway capability for arsenate reduction (derived from Madin et al. compilation; mechanism unresolved between dissimilatory respiration and cytoplasmic detoxification)." | Do not classify as *arsC* detoxification. Do not classify as *arrAB* respiration. Do not infer gene presence, plasmid carriage, arsenic exposure, or co-selection. |
| `METPO:1000623` | **Keep separate label / Unresolved quantitative conditions** | "Taxa annotated as moderately halophilic (defined in METPO as requiring high levels of sodium chloride, typically $\ge 0.2$ M NaCl)." | Do not read as "halotolerant" (halophily requires salt; halotolerance does not). Do not infer continuous salt optima, ranges, or media. |

---

## 6. Downstream Handoff Contracts

### 6.1 Handoff to Issue #178 (Integrated Evidence Table Contract)
Issue [#178](https://github.com/CultureBotAI/bsdbng/issues/178) defines the shared result contract for multi-arm enrichment integration. This audit supplies the following metadata requirements for Table 1 features:
- **Raw trait identifiers:** Must preserve qualified keys (`METPO:1000608`, `METPO:2000103|GO:0030245`, `METPO:1000668`, `METPO:1000683`, `METPO:2000103|pathways:arsenate_reduction`, `METPO:1000623`).
- **Input fingerprinting:** Pinned MD5 checksums for nodes (`dd76b4656f9341db25a7c4afa96cd58e`) and edges (`555ecf5ba5d89288cdc44e368375374f`).
- **Feature universe denominators ($K$):** Coccus ($K=675$) and Sphere ($K=283$) must be evaluated as separate testing units unless a deduplicated union universe is formally constructed.
- **Evidence status encoding:** Carry `measured`, `propagated`, and `unresolved` flags into result tables to distinguish direct isolate assertions from taxonomic rollups.

### 6.2 Handoff to Issue #179 (Project-Wide Biological Interpretation)
Issue [#179](https://github.com/CultureBotAI/bsdbng/issues/179) governs biological interpretation across the manuscript, figures, and automated reports. The boundaries established here provide canonical examples:
- *Evidence chain:* Source literature $\to$ curated isolate traits $\to$ NCBI taxonomy rollup $\to$ BugSigDB signature overlap $\to$ hypergeometric statistic.
- *Defensible claims:* A significant shift indicates a higher proportion of taxa in a differential abundance signature possessing the recorded trait annotation.
- *Prohibited claims:* Statistical enrichment in a differential signature does not demonstrate in vivo metabolic fluxes, gene expression levels, environmental chemical exposure, or host causation.
