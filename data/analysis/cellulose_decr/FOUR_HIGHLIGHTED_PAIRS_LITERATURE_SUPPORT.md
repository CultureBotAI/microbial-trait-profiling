# Four highlighted condition–trait pairs: corpus evidence and literature

This file covers the four exploratory leads in `EXPLORATORY_LEADS`
(`src/mtp/deck_visualizations.py`), all in human feces:

1. Parkinson's disease (PD), `METPO:2000103|GO:0030245` cellulose catabolic process, decreased
2. HIV infection, cellulose catabolic process, decreased
3. Colorectal cancer (CRC), cellulose catabolic process, decreased
4. PD, `METPO:1000623` moderately halophilic, increased

It was rewritten on 2026-10-17 (#552). The previous version gave wrong PubMed
identifiers for 13 of its citations, attributed findings to papers that do not
report them, and presented corpus studies as independent support. Every
identifier below was checked against PubMed, every description of a paper
against its abstract or full text, and every taxon against
[`four_highlighted_pairs_study_evidence.tsv`](four_highlighted_pairs_study_evidence.tsv).
The corrections are listed at the end. The abstract check also found one curation
problem: McHardy et al. 2013 is curated as feces but sampled rectal mucosa (#566).

## Reading this file

- **Corpus studies are not independent support.** A paper curated in this
  BugSigDB vintage (`2026-04-24_00:41_UTC`) contributes to the result itself; its
  agreement is not replication. Each citation below says whether it is in the
  corpus.
- **Trait names are literal.** A taxon counted here bears a KG-Microbe
  annotation because it or a descendant has the annotation edge. None of the 17
  taxa in all three cellulose sets has its own edge; the genera inherit it from
  annotated species or strains (#565). "Cellulose catabolic process"
  does not mean cellulose was degraded, and "moderately halophilic" does not mean
  salt or osmotic tolerance in the sampled strain
  ([`docs/table1_trait_semantics_audit.md`](../../../docs/table1_trait_semantics_audit.md),
  #189). Mechanisms mentioned in the cited papers are hypotheses for this
  analysis, not findings of it.
- Counts use each experiment's own condition and body site, the same rule as the
  rest of this directory, so they agree with
  [`cellulose_decreased_taxa_catalog.tsv`](cellulose_decreased_taxa_catalog.tsv).

## 1. Corpus evidence for each pair

"Studies" are feces studies in the condition whose signature in the stated
direction reports at least one annotated taxon; the decreased / increased counts
are studies reporting each taxon in the stated and the opposite direction
([`four_highlighted_pairs_taxon_drivers.tsv`](four_highlighted_pairs_taxon_drivers.tsv)).

| Pair | Studies | Annotated taxa | Most reported taxa (stated / opposite direction studies) |
| --- | --- | --- | --- |
| 1. PD, cellulose, decreased | 40 | 38 | Lachnospiraceae 20 / 7, *Roseburia* 20 / 4, *Blautia* 13 / 2, Bacillota 8 / 5, Clostridiaceae 8 / 2, Oscillospiraceae 8 / 13, *Clostridium* 7 / 3 |
| 2. HIV, cellulose, decreased | 16 | 24 | *Bacteroides* 8 / 0, *Clostridium* 5 / 2, Bacteroidaceae 5 / 0, Lachnospiraceae 5 / 3, Clostridia 4 / 0, Eubacteriales 4 / 1, *Ruminococcus* 4 / 0 |
| 3. CRC, cellulose, decreased | 22 | 34 | Lachnospiraceae 9 / 5, Clostridia 7 / 1, *Clostridium* 6 / 3, *Blautia* 6 / 1, *Lachnospira* 6 / 1, *Roseburia* 6 / 0 |
| 4. PD, moderately halophilic, increased | 43 | 79 | *Lactobacillus* 12 / 5, Enterobacteriaceae 9 / 2, Pseudomonadota 9 / 1, Bifidobacteriaceae 9 / 0, *Streptococcus* 9 / 5, *Alistipes* 8 / 0 |

For the three cellulose pairs, the paired direction tests in
[`README.md`](README.md#4-direction-within-each-condition) give, at the study
level among disease-versus-healthy contrasts:

- PD: 24 of 34 untied studies down (p = 0.024).
- HIV: 10 of 12 (p = 0.039) with the exact label. Counting two studies filed as
  "Human immunodeficiency virus" gives 10 of 14 (p = 0.18), so the HIV result
  is fragile (#561, #580).
- CRC: 14 of 17 (p = 0.013).

Taxon overlap between the three conditions is not detectably larger than among
random sets of feces studies (U6 p = 0.19;
[`README.md`](README.md#2-is-the-overlap-unusual)). No direction test for pair 4
is run in this directory.

## 2. Corpus studies the previous version cited

Taxa each study's signature reports in the stated direction, from the evidence
table. These studies are part of the result, not support for it.

**PD, cellulose, decreased**

| Study | Citation | Annotated taxa reported decreased |
| --- | --- | --- |
| `bsdb:26179554` | Keshavarzian et al. 2015, *Mov Disord*, PMID 26179554 | Bacillota, *Blautia*, Lachnospiraceae, *Roseburia* |
| `bsdb:27591074` | Unger et al. 2016, *Parkinsonism Relat Disord*, PMID 27591074 | Bacteroidia only |
| `bsdb:28195358` | Hill-Burns et al. 2017, *Mov Disord*, PMID 28195358 | *Blautia*, Lachnospiraceae, *Roseburia* |
| `bsdb:28843021` | Heintz-Buschart et al. 2018, *Mov Disord*, PMID 28843021 | Bacillota, *Bacteroides*, Eubacteriales, Lachnospiraceae, Oscillospiraceae, Pseudomonadota, *Roseburia*, *Ruminococcus* |
| `bsdb:31221587` | Aho et al. 2019, *EBioMedicine*, PMID 31221587 | ten taxa, among them *Roseburia*, *Blautia*, *Clostridium* and Lachnospiraceae |
| `bsdb:32566740` | Wallen et al. 2020, *NPJ Parkinsons Dis*, PMID 32566740 | *Blautia*, *Lachnospira*, Lachnospiraceae, *Roseburia* |
| `bsdb:32934117` | Vascellari et al. 2020, *mSystems*, PMID 32934117 | ten taxa, among them *Bacteroides*, *Blautia*, *Butyrivibrio*, *Lachnospira* and *Roseburia* |
| `bsdb:36376318` | Wallen et al. 2022, *Nature Communications*, PMID 36376318 | *Anaerostipes hadrus*, *Roseburia faecis* |
| `bsdb:38773112` | Nishiwaki et al. 2024, *NPJ Parkinsons Dis*, PMID 38773112 | *Blautia*, *Clostridium*, *Lacrimispora*, *Mediterraneibacter*, *Roseburia*, *Roseburia faecis* |

Wallen et al. 2020 and 2022 come from the same investigators. The 2022 paper
reports that, except for 11 subjects enrolled twice, its 724 subjects are
independent of the 840 in the earlier 16S cohorts (N = 333 and 507). Nishiwaki et
al. 2024 meta-analysed its own Japanese dataset (94 PD, 73 controls) with five
published shotgun datasets, one of them from the USA.

**HIV infection, cellulose, decreased**

| Study | Citation | Annotated taxa reported decreased |
| --- | --- | --- |
| `bsdb:24451087` | McHardy et al. 2013, *Microbiome*, PMID 24451087; rectal mucosa sampled with a sponge, curated as feces (#566) | eight taxa, among them Lachnospiraceae, *Roseburia*, *Ruminococcus* and *Lachnospira*; the only HIV source of *Eubacterium* |
| `bsdb:27477587` | Ling et al. 2016, *Sci Rep*, PMID 27477587 | *Bacteroides* and its family, order, class and phylum |
| `bsdb:29411528` | Zhou et al. 2018, *J Cell Mol Med*, PMID 29411528 | nine taxa, among them *Bacteroides*, Lachnospiraceae and *Lachnoclostridium* |
| `bsdb:30518941` | Kaur et al. 2018, *Sci Rep*, PMID 30518941 | eight taxa, among them *Bacteroides*, Lachnospiraceae and *Ruminococcus* |
| `bsdb:33102244` | Zhou et al. 2020, *Front Cell Infect Microbiol*, PMID 33102244 | *Anaerostipes*, *Bacteroides*, *Blautia*, *Butyrivibrio*, *Clostridium*, *Lachnospira*, *Roseburia*, *Ruminococcus* |
| `bsdb:34378948` | Ishizaka et al. 2021, *Microbiol Spectr*, PMID 34378948 | six taxa, among them *Anaerostipes*, *Bacteroides* and Oscillospiraceae |
| `bsdb:36604616` | Chandiwana et al. 2023, *BMC Microbiol*, PMID 36604616 | *Bacteroides*, Bacteroidaceae, *Clostridium*, Oscillospiraceae, *Ruminococcus* |
| `bsdb:39880958` | Maghini et al. 2025, *Nature*, PMID 39880958 | eight taxa, among them Lachnospiraceae, Oscillospiraceae and *Clostridium* |

Zhou et al. 2020 compared groups by HIV status and sexual orientation. One of its
experiments compares MSM with non-MSM and is classified `unresolved` in the
contrast manifest. Another compares HIV– with HIV+ individuals
(`bsdb:33102244-1`). Since classifier v4 it is `disease_vs_healthy` with
screened controls (Tier 2), so the disease-versus-healthy tests count it (#569).
The HIV study-level result is 10 of 10 (p = 0.0020) when the contrasts whose cases
are named as treated are set aside, or with Tier 1–2 controls only (#580).


Vujkovic-Cvijin et al. 2020 (*Nature Communications*, PMID 32415070,
`bsdb:32415070`) is also a corpus study for pair 2; it reports Lachnospiraceae.
It compared antiretroviral-treated people with HIV with seronegative controls
matched for age, body-mass index, sex and sexual practice. It still reported
Lachnospiraceae and Ruminococcaceae depletion.

**CRC, cellulose, decreased**

| Study | Citation | Annotated taxa reported decreased |
| --- | --- | --- |
| `bsdb:21850056` | Wang et al. 2012, *ISME J*, PMID 21850056 | *Bacteroides*, *Eubacterium*, *Roseburia* |
| `bsdb:24316595` | Ahn et al. 2013, *J Natl Cancer Inst*, PMID 24316595 | Clostridia, Eubacteriales, Lachnospiraceae |
| `bsdb:25104642` | Zackular et al. 2014, *Cancer Prev Res*, PMID 25104642 | *Bacteroides*, *Blautia*, *Clostridium*, Eubacteriales, Lachnospiraceae, *Ruminococcus* |
| `bsdb:30065719` | Deng et al. 2018, *Front Microbiol*, PMID 30065719 | many taxa, among them *Anaerostipes*, *Lachnospira*, *Roseburia* and *Ruminococcus flavefaciens* |
| `bsdb:31281534` | Yang et al. 2019, *Theranostics*, PMID 31281534 | eight taxa, among them *Blautia*, *Lachnospira* and Lachnospiraceae |
| `bsdb:33947803` | Jiang et al. 2021, *mSystems*, PMID 33947803 | *Anaerostipes hadrus*, *Bacteroides cellulosilyticus*, *Ruminococcus bromii*, *[Eubacterium] siraeum* |
| `bsdb:34799562` | Yang et al. 2021, *Nature Communications*, PMID 34799562 | eight taxa, among them *Roseburia*, *Lachnospira*, *Blautia* and *Ruminococcus champanellensis* |

*Ruminococcus champanellensis* is reported decreased by this one CRC study only.

**PD, moderately halophilic, increased**

| Study | Citation | Annotated taxa reported increased |
| --- | --- | --- |
| `bsdb:27591074` | Unger et al. 2016, *Parkinsonism Relat Disord*, PMID 27591074 | Enterobacteriaceae |
| `bsdb:28195358` | Hill-Burns et al. 2017, *Mov Disord*, PMID 28195358 | Bifidobacteriaceae, Lactobacillaceae, *Lactobacillus* |
| `bsdb:28449715` | Bedarf et al. 2017, *Genome Med*, PMID 28449715 | Bacillota, Clostridia |
| `bsdb:29501802` | Qian et al. 2018, *Brain Behav Immun*, PMID 29501802 | many taxa, among them *Alistipes*, *Klebsiella* and *Acinetobacter* |
| `bsdb:30576008` | Barichella et al. 2019, *Mov Disord*, PMID 30576008 | many taxa, among them Enterobacteriaceae, *Escherichia*, *Lactobacillus* and Bifidobacteriaceae |
| `bsdb:31248424` | Lin et al. 2019, *J Neuroinflammation*, PMID 31248424 | many taxa, among them *Enterococcus*, *Lactobacillus* and Rikenellaceae |
| `bsdb:32934117` | Vascellari et al. 2020, *mSystems*, PMID 32934117 | ten taxa, among them Bifidobacteriaceae, *Escherichia* and *Streptococcus* |

Several taxa here, such as Pseudomonadota (211 annotated taxa beneath it), carry
the annotation only through rollup from descendants.

The full list of supporting studies, with every taxon in both directions, is in
[`four_highlighted_pairs_study_evidence.tsv`](four_highlighted_pairs_study_evidence.tsv).

## 3. Literature outside the corpus

None of these papers is in this BugSigDB vintage. They provide context and
hypotheses; none validates the annotation results above.

| Paper | What it reports | Relation to this analysis |
| --- | --- | --- |
| Romano et al. 2021, *NPJ Parkinsons Dis* 7:27, PMID 33692356 | Re-analysis of ten 16S case–control datasets (nine cohorts, six countries; 1,211 samples after filtering, 681 PD and 530 controls). Most consistent changes: Lactobacillus, Akkermansia and Bifidobacterium enriched; Lachnospiraceae and Faecalibacterium depleted. | Eight of its ten datasets are corpus studies (Hill-Burns 2017, Hopfner 2017, Keshavarzian 2015, Pietrucci 2019, Qian 2018, Aho 2019, Heintz-Buschart 2018, Weis 2019); only Scheperjans 2015 and Petrov 2017 are not. It is largely the same evidence, re-analysed. |
| Shen et al. 2021, *Front Aging Neurosci* 13:636545, PMID 33643026 | Meta-analysis of 15 PD case–control studies: Prevotellaceae, Faecalibacterium and Lachnospiraceae lower; Bifidobacteriaceae, Ruminococcaceae, Verrucomicrobiaceae and Christensenellaceae higher in PD. | Overlap with the corpus not assessed. Ruminococcaceae are reported higher in PD. |
| Becker et al. 2022, *Genomics Proteomics Bioinformatics* 20:274, PMID 34839011 | Open-label trial (RESISTA-PD) of resistant starch for 8 weeks: fecal butyrate rose and calprotectin fell in PD patients; alpha and beta diversity, including SCFA-producing bacteria, stayed stable. | Not evidence about cellulose-annotated taxa. |
| Lozupone et al. 2013, *Cell Host Microbe* 14:329, PMID 24034618 | HIV infection associated with characteristic fecal community changes, not restored by antiretroviral therapy, with higher diversity and resemblance to Prevotella-rich communities. | Context for HIV. |
| Vujkovic-Cvijin et al. 2013, *Sci Transl Med* 5:193ra91, PMID 23843452 | A mucosal-adherent community enriched in Proteobacteria and depleted of Bacteroidia, associated with mucosal immune disruption and kynurenine-pathway activity. | Mucosal samples, not feces. |
| Noguera-Julian et al. 2016, *EBioMedicine* 5:135, PMID 27077120 | In Barcelona (n = 156) and Stockholm (n = 84), the Prevotella-rich pattern tracked MSM status independently of HIV. After stratifying for sexual orientation there was no solid evidence of HIV-specific dysbiosis; HIV remained associated with lower richness. | A caveat for pair 2: few HIV studies here separate HIV from sexual behaviour. Vujkovic-Cvijin et al. 2020 (section 2) is one that does, and it still reports Lachnospiraceae depletion. |
| Wirbel et al. 2019, *Nat Med* 25:679, PMID 30936547 | Meta-analysis of CRC fecal metagenomes: protein and mucin catabolism genes enriched, carbohydrate degradation genes depleted. | Context for CRC. Its companion, Thomas et al. 2019 (PMID 30936548), is a corpus study that reports no annotated taxa for pair 3. |
| Chassard et al. 2012, *Int J Syst Evol Microbiol* 62:138, PMID 21357460 | *Ruminococcus champanellensis*, a strictly anaerobic bacterium isolated from human feces that degrades microcrystalline cellulose. | Describes one species in the CRC set; one study reports it decreased. |
| Moraïs et al. 2024, *Science* 383:eadj9223, PMID 38484069 | Cellulosome-producing ruminococcal species are common in ancient, hunter-gatherer and rural populations and rare in industrialized ones. | Context only; it does not compare disease and control. |
| Maini Rekdal et al. 2019, *Science* 364:eaau6323, PMID 31196984; van Kessel et al. 2019, *Nature Communications* 10:310, PMID 30659181 | Gut bacterial tyrosine decarboxylases, including that of *Enterococcus faecalis*, convert levodopa to dopamine and lower its availability. | Levodopa treatment is a possible confounder of PD signatures; neither paper concerns halophily. |
| Savica et al. 2009, *Neurology* 73:1752, PMID 19933976; Cersosimo et al. 2013, *J Neurol* 260:1332, PMID 23263478 | Constipation is more common before PD onset (odds ratio 2.48) and precedes motor symptoms. | Context; stool transit and consistency are not measured in these signatures. |

Flemer et al. 2017 (*Gut*, PMID 26992426), Burns et al. 2015 (*Genome Med*, PMID
26170900) and Vandeputte et al. 2016 (*Gut*, PMID 26069274) are corpus studies
but report no annotated taxa for these pairs. The curated experiments of the
first two are mucosal or tissue samples, although Flemer et al. also sampled
feces. Vandeputte et al. studied stool consistency in 53 healthy women.

## 4. Corrections to the previous version

| Previously cited as | What the cited identifier actually is | Correct record |
| --- | --- | --- |
| Romano et al. 2021, PMID 33640941 | Heiman et al. 2021, *Br J Surg* (breast cancer surgery) | PMID 33692356 |
| Wallen et al. 2022, PMID 36460672 | Chen et al. 2022, *Nature Communications* (aerosol remote sensing) | PMID 36376318; a corpus study supporting pairs 1 and 4 |
| Becker et al. 2023, *NPJ Parkinsons Dis*, PMID 36446995, DOI 10.1038/s41531-022-00430-8 | Bawa et al. 2023 (pharmacist survey); the DOI is Wang et al. 2023 (C/EBPβ/AEP) | Becker et al. 2022, *Genomics Proteomics Bioinformatics*, PMID 34839011 |
| Noguera-Julian et al. 2016, PMID 27077129 | Haahr et al. 2016 (intracavernous stem cells) | PMID 27077120 |
| Serrano-Villar et al. 2017, *EBioMedicine*, PMID 27846282, DOI 10.1016/j.ebiom.2016.11.011 | Yu et al. 2016 (stroke case mix); the DOI is Goheen et al. 2016 (malaria) | not identified; removed |
| Chassard et al. 2012, *Environ Microbiol*, PMID 22176653, DOI 10.1111/j.1462-2920.2011.02678.x | Goldstein et al. 2012 (Fontan circulation); the DOI is Villa et al. 2012 (*E. coli* antifouling) | *Int J Syst Evol Microbiol*, PMID 21357460 |
| Maini Rekdal et al. 2019, PMID 31197014 | Gu et al. 2019 (CO2 electroreduction) | PMID 31196984 |
| Flemer et al. 2017, PMID 27993829, DOI 10.1136/gutjnl-2015-309897 | Westhovens et al. 2017 (filgotinib); the DOI is Blasco-Baque et al. 2017 (periodontitis) | PMID 26992426; a corpus study |
| Burns et al. 2015, PMID 26597098 | Kleiter et al. 2016 (neuromyelitis optica) | PMID 26170900; a corpus study |
| Shen et al. 2021, *Aging*, PMID 33495415 | Umeda-Kameyama et al. 2021 (Alzheimer's facial screening) | *Front Aging Neurosci*, PMID 33643026 |
| Savica et al. 2009, *Arch Neurol*, PMID 19364936 | Brugman et al. 2009 (spastic paraparesis) | *Neurology*, PMID 19933976 |
| Cersosimo et al. 2013, PMID 23180178 | Langford et al. 2013 (THC/CBD spray) | PMID 23263478 |
| Weis et al. 2019, PMID 31391112 | Pereira et al. 2019 (incisional hernia) | *NPJ Parkinsons Dis*, PMID 31815177; a corpus study |
| Qian et al. 2018, *Parkinsonism Relat Disord*, DOI 10.1016/j.parkreldis.2018.02.026 | the DOI is Haehner et al. 2018 (substantia nigra imaging) | *Brain Behav Immun*, PMID 29501802 (PMID was right) |
| Lin et al. 2019, *Parkinsonism Relat Disord*, DOI 10.1016/j.parkreldis.2019.06.001 | the DOI is Mehanna et al. 2019 (young-onset PD) | *J Neuroinflammation*, PMID 31248424 (PMID was right) |
| Scheperjans et al. 2015 as corpus study `bsdb:25476529` | not in this BugSigDB vintage | external |

Claims that the cited papers do not make:

- Wallen et al. 2022: depletion of cellulose, xylan and pectin degradation
  pathways. The full text does not mention cellulose or xylan.
- Romano et al. 2021: a *Roseburia*, *Blautia*, *Lachnospira* and
  Ruminococcaceae "plant-fiber guild" depleted regardless of region. The paper
  reports Lachnospiraceae and Faecalibacterium depletion; Shen et al. 2021 report
  Ruminococcaceae higher in PD.
- Noguera-Julian et al. 2016: HIV depleting Lachnospiraceae, Ruminococcaceae and
  Clostridia independently of MSM status. The paper found no solid evidence of
  HIV-specific dysbiosis after stratifying for sexual orientation.
- Becker et al. 2022: increased Lachnospiraceae after resistant starch. The paper
  reports stable diversity, including SCFA producers.
- Bedarf et al. 2017: depletion of Lachnospiraceae carrying carbohydrate-active
  enzymes. The abstract reports more *Akkermansia muciniphila* and less
  *Prevotella copri* and *Eubacterium biforme*.
- Lozupone et al. 2013: depletion of *Bacteroides* polysaccharide specialists such
  as *B. cellulosilyticus*. The paper reports Prevotella-rich communities and
  higher diversity.
- Size and scope: Wallen et al. 2020 analysed N = 333 and 507, not 515 and 333;
  Nishiwaki et al. 2024 combined six datasets, not a six-country study of more
  than 2,200; Vandeputte et al. 2016 studied 53 women, not 1,126 people.
- Per-study taxa: Unger et al. 2016 reports only Bacteroidia among annotated
  taxa decreased in PD feces, and Jiang et al. 2021 does not report *Roseburia
  faecis*.

The previous version also translated "moderately halophilic" into osmotic, bile
and salt tolerance driven by slow transit and levodopa, and described the
cellulose-annotated taxa as a "fiber-degrading guild" whose loss lowers
short-chain fatty acids. Those are hypotheses that these annotations cannot
test, and they have been removed.
