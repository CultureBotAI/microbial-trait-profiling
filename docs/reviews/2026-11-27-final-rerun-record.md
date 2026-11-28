# Final rerun record: P 358485d / P2 and C2 559ab0b, regenerated in #2397

This document records the accepted final rerun. It names the runs, the commits and inputs they used, each run's population and strain-rollup scope, what was verified, and which tracked outputs this PR (#2398, issue #2397) regenerates from that run. It replaces [the 2026-10-23 final rerun record](2026-10-23-final-rerun-record.md) and [the 2026-10-25 rebuild record](2026-10-25-rebuild-of-F-record.md) as the current record. Both stay as historical records of run F (`511c891`).

How to cite numbers from this record:
- **Core-bound results** must come from one of the verified runs below.
- **Tracked outputs** this PR regenerates are built on the branch from the same pinned inputs. They are not core-verified. Cite them by file, together with the verification status their own provenance records.
- **Population and scope:** whichever source a number comes from, state its population and strain-rollup scope with it.
- **What the numbers measure:** every result here measures how taxa reported in BugSigDB signatures line up with KG-Microbe trait annotations. None of them measures physiological activity, abundance or metabolite flux in the host (see `methods/biological_interpretation.md`).

## Code and inputs

The final run has three layers, each with its own commit. Every layer's records report a clean tree.

| Layer | What it produced | Commit | Dirty |
|---|---|---|---|
| P | preflight; the four cores; the ilyas arm | `358485dd62b804254babe4a9b6d27f0b371205df` | `git_dirty: false` for core, arms and comparison in every unit record [S1, S2] |
| P2, attempt `canary04` | jed and noel arms; comparisons; native influence runs | producer `559ab0b8d094c544f6b56cb50e1142c38fd00d48` (core and ilyas stay at 358485dd) | `git_dirty: false` for every producer under `producers.*` and every native run under `native_evidence.*.execution_after` in each unit's receipt [S3] |
| C2 (`c2-canary04d` + `c2-canary04e`) | additive sidecars; evidence archives; canonical-phenotype metadata | consumer `559ab0b8…`, P commit `358485dd…` | the 04e worktree is clean at `559ab0b8…` [S4] |

The regeneration branch `2397-post-run-regeneration-regen-plan-v5` starts at `origin/main` `83c90ae5a94aba4ca5e2f0887a8e21b811e0f0a3` [S5, S6]. It has since merged main five times: 351b10b0, 25a7f9b8, 315032f3, 98868b49 (main `2c39a7cb`) and 830b23a6 (main `a302aab5`, PR #2452) [S6]. Against current main, the only file under `src/` or `scripts/` that the branch changes is `scripts/adhoc/issue1316_validate_committed.py`, whose `VISUAL_JSON` list e6e14733 and ac41537a edit [S6, S7]. The branch also changes 17 test files, `.pre-commit-config.yaml` and `.secrets.baseline` [S6].

**Environment.**
- **Final run:** P's core records `lock_sha256` `5d4b47085d3266cf…` [S2]. The native runs record Python 3.11.5 and scipy 1.17.1 [S8].
- **Regeneration branch:** main pins pandas 3.1.1 and Python >=3.11.5; P2 and C2 ran pandas 3.0.6 [S5]. So outputs rebuilt on the branch are compared with P2 and C2 by content, not by bytes (REGEN_PLAN_v5 A2) [S5].
- **Branch lock:** `uv.lock` is unchanged from 83c90ae5. Its sha256 is `e54339c6213df7df…` [S6].

**KG-Microbe source.** The edges and nodes are pinned by md5 and sha256 [S9]:

| File | md5 | sha256 |
|---|---|---|
| `merged-kg_edges.tsv` | `555ecf5ba5d89288cdc44e368375374f` | `d13812a436bd0b1c…` |
| `merged-kg_nodes.tsv` | `dd76b4656f9341db25a7c4afa96cd58e` | `8b4cda73ca7e9d02…` |

- **Release:** `data/kg-microbe/provenance.json` records release `20260120` [S10].
- **P preflight:** the `inputs` receipt passed, and so did its `copied-inputs` receipt for the P worktree copies [S9].
- **Branch:** the inputs were cloned from the P worktree before any recipe ran. Their sha256 equal the pins [S11, S6]. `just verify-pinned-kg` passed on the branch ("match the md5s in table1_taxon_influence_protocol.json") [S12].
- **Source mode:** all four P core `run.json` files record the run parameter `source_mode` = `pinned` (source `override`) [S2].

**BugSigDB export.** `full_dump.csv` has sha256 `af84bfca7bfe7691…` and md5 `00df5a86c684eac1…` [S9].

**Corpus.** These tables were built at P preflight and rebuilt on the branch in Phase 1. Every hash below matches between the P worktree and the branch worktree [S6, S13]:

| File | sha256 | Contents |
|---|---|---|
| `data/derived/experiments.tsv` | `16271237a70f61c8195d9f2acdd16f71e88f65859160f95bb734626bc1328e06` | 4,435 rows, 12 columns [S13] |
| `data/derived/study_trait_profiles_wide.tsv` (phenotype rollup, the canonical table) | `dd20b920308d1ab9dd9ff58681e17de2526d2e70cc2cbce0b55e924e0a661ba2` | 5,629 profiles from 3,372 experiments in 1,047 studies; 2,635 columns: 3 key columns + 2,632 trait columns [S14] |
| the off-rollup profile table (`P-358485d/off-derived/`) | `6fdafa21dd8fedacaa757fdf030a6a2d3af0536ba37dbd8aac014b1af7f51977` | 5,629 profiles, 2,609 columns: 3 + 2,606 trait columns [S14] |
| `data/derived/taxa_trait_matrix_columns.tsv` | `db8cdc96…1774cae` | |
| `data/derived/taxa_trait_matrix.tsv` | `6638849b…bf389af` | |
| `data/derived/studies.tsv` | `410b60f4…6e044e69c` | |
| `data/derived/bugsigdb_taxa.tsv` | `e7fa982b…db271c59` | |

- **Paired subset:** the canonical corpus has 2,257 matched pairs from 932 studies. Of its 5,629 profiles, 2,948 come from signatures reported increased and 2,681 from signatures reported decreased [S4].
- **Compared with run F:** F's phenotype table, `43b2e230…`, had 2,637 columns (2,634 traits). The current table has 2 fewer trait columns. The default trait-key policy `facultative_union_v1` merges METPO:1000605 (facultatively anaerobic) and METPO:1000612 (facultative oxygen preference) into METPO:1000608 [S6 (0bbac143, 60cd4086)]. The merged trait is labelled "facultative (aerobe or anaerobe)" in the regenerated outputs that read the policy label [S6 (a01466da)]; some regenerated outputs still carry "facultatively aerobic" for METPO:1000608, for example `table1_traits[0].trait_name` in `composition_calibration_summary.json` (#2437, #2439) [S26].
- **Committed copy:** the tracked `data/derived/study_trait_profiles_wide.tsv.gz` decompresses to `dd20b920…` [S6].

## Verified core-bound runs

Each unit ran these steps:
1. P: `just rerun-core`;
2. P: the ilyas arm;
3. P2 (`559ab0b8`): the jed and noel arms;
4. P2: the comparison;
5. P2: per-unit verification.

The P2 canary was tier12-phenotype. Its approval (2026-11-12) records a completion check of 122 PASS, 0 FAIL and 9 NOT CHECKED, plus a cell-by-cell comparison with canary03 [S15].

| Unit | Contrast population | Rollup | Count profiles | Experiments | Studies | Retained traits |
|---|---|---|---|---|---|---|
| Tier 1+2 × off | `disease_vs_healthy_tier1_2_subject_controls` | off | 1,870 | 935 | 548 | 316 |
| Tier 1+2 × phenotype (P2 canary) | `disease_vs_healthy_tier1_2_subject_controls` | phenotype | 1,870 | 935 | 548 | 342 |
| Full corpus × phenotype (**primary**) | `all_paired` | phenotype | 7,111 | 4,299 | 1,212 | 342 |
| Full corpus × off | `all_paired` | off | 7,111 | 4,299 | 1,212 | 316 |

<!-- counts: P runs/<unit>/verified(.amended).json profiles.count_profiles; retained traits: analysis_bases[].retained_feature_count (P and P2 agree) [S1, S3] -->

Notes on the counts:
- **What is counted:** the profiles, experiments and studies in each core's count-profile product. One profile merges one experiment's signatures in one direction.
- **ORA profiles:** in the full units, the ORA profiles are the canonical human 16S/WMS population: 5,629 rows, 3,372 experiments and 1,047 studies [S1].
- **Trait counts:** the retained trait counts are 316 (off) and 342 (phenotype). Run F's were 317 and 344 [S1, old record].
- **Shared traits:** the paired comparison finds 342 phenotype traits = 316 shared + 26 phenotype-only [S16].
- **Primary unit:** full-phenotype is primary (owner decision 2, REGEN_PLAN_v5 A10). Its population is the canonical corpus that Phases 1–4 rebuild. Tier1+2 × phenotype is kept as an archived subset/sensitivity result [S5]. The off units are declared sensitivity analyses (`scope_role` off = "sensitivity") [S16].

**Full-unit receipt amendment.** P's first unit receipt for full-phenotype and full-off has `status: failed` with `ValueError: Count and ORA profile keys differ`. Both are superseded by `verified.amended.json` (`status: passed`, `p-unit-amendment-v1`).
- **Reason recorded:** the frozen receipt wrongly required identical count and ORA keys for `all_paired`. The amended rule requires the ORA keys to equal the scope-matched 5,629-key preflight population and to be a subset of the count keys.
- **Fresh check:** a fresh `just rerun-verify-comparison` exited 0 [S1].

The P2 receipts for all four units are `status: passed` [S3]. Their scope is `mixed_producer_bindings_and_saved_native_execution_evidence_no_analysis_replay`, so they bind saved execution evidence and don't replay the analysis.

| Unit | Core run (P) | Arm runs ilyas (P), jed, noel (P2) | Comparison (P2) |
|---|---|---|---|
| Tier 1+2 × off | `run-9b73c70503134d2db0bc665b1de84dfd` | `arm-a760cc15043c4e958fe62b5651dfe30f`, `arm-14b3a66fd23c4ad994d9613fd791efe4`, `arm-3ce9483a0c374cc29428107588fcb2b0` | `comparison-disease_vs_healthy_tier1_2_subject_controls-6f465964c99b` |
| Tier 1+2 × phenotype | `run-45e0c71a8fe14187b887a50644283250` | `arm-926bb3573b5c4fa683739f8ad08f6868`, `arm-144ef50414cf422e9904f4da0574a272`, `arm-51abc5214d6749e58d2f7d2784f7dc69` | `comparison-disease_vs_healthy_tier1_2_subject_controls-77b3763d5e1e` |
| Full corpus × phenotype | `run-c06d9783cd8a41c9bc904cfd5abf55a9` | `arm-b609c5cf5228411e9d267059373a3a8e`, `arm-e8407c48405d48a48d05d7316c039b5c`, `arm-4bb9b738c3f245e18fd5fdc028e0f78f` | `comparison-all_paired-67c4d46656d3` |
| Full corpus × off | `run-7b38cd1eba1042ccac719e3e8457f078` | `arm-d89756a799ef43fca8abf79272b7cfc2`, `arm-85664a053f0448cbb7968475fffd1a4d`, `arm-7bb516c31ced483f9b64fec9eeb9895f` | `comparison-all_paired-a4de8800ecfb` |

<!-- P's own jed/noel arms and comparisons (e.g. tier12-off comparison-…-951834d64ec4) are superseded by P2's; IDs above are from P2 runs/<unit>/verified.json [S3] -->

**Calibration verdicts.** All four units report the same two verdicts, in P and in P2 [S1, S3]:
- `noel_table1_six_trait_sign_test_v1`: **partially calibrated.**
  - The null control passed. `power_passed` is null: power is reported, but no target-power criterion was declared.
  - The negative control, `screened_fraction_sign_test`, fails family error, as the protocol declares it will.
  - The verdict records that the six Table 1 traits were selected on this corpus before testing, so their p-values are exploratory (#363).
- `jed_log_ratio_vs_transformed_t_calibration`: **different estimator, uncalibrated.**
  - The harness evaluated `jed_composition_arcsine_sqrt_ttest_calibration_v1`, not the production estimator `jed_table1_paired_log_ratio_sign_test_v1`.

**Native influence runs.** All four units have a completed native run: `completed=True`, `code_dirty=False`, protocol `190-native-v2` [S17]. The full-phenotype run is published in Phase 2d (see below). The tier12-phenotype run is not published into the repository [S5].

### Additive analyses (C2) and where the evidence is

C2 ran at consumer `559ab0b8` against these cores. Its sidecars are separately named and not committed. They are listed here by unit:
- **tier12-off and tier12-phenotype:** delta-z-derived, density-trait-list, influence-extension, label-support (+ trait list), matrix-overlap-coverage, multiplicity, selection-scan, three-test.
- **full-phenotype and full-off:** the same, without influence-extension, matrix-overlap-coverage and selection-scan.
- **paired:** rollup-comparison and selection-scan-scopes.
- **canonical-phenotype:** metadata and row4 [S18].

Lead acceptance, stage by stage:
- **tier12-off:** 9/9 sidecar receipts [S19].
- **paired:** gate passed; phenotype is `primary`, off is `declared_sensitivity` [S16].
- **full-phenotype:** 6/6 receipts. Of its 5,629 input profiles, 2,257 paired experiments (4,514 profiles) from 932 studies form the pairs [S20].
- **full-off:** 6/6 receipts, with the same population [S21].
- **canonical-phenotype:** metadata 5,629 profiles and 2,257 pairs; row4 has 4,814 contrast cells, 0 blank [S4].

This record doesn't re-derive the additive values. Read them from the archived sidecars. Every p-value they report is nominal.

**Acceptance.** C2 was accepted on 2026-11-23. The final result is c2-canary04d stages 1–8 plus c2-canary04e canonical [S4, S22]. c2-canary04d's canonical stage stopped at verify-metadata because of a blank first header (#2396). It was repaired and rerun as c2-canary04e. The 8 metadata outputs are byte-identical between the two tokens [S4].

**Evidence location.** The run directories are not committed. They are kept outside the repository and outside temporary storage:
- **P (`mtp-runs/P-358485d/`):** `preflight/`, `runs/<unit>/` (core, arms, comparison, `verified.json`, plus `verified.amended.json` for the full units), and `wt/` [S1, S13].
- **P2 (`mtp-runs/P2-559ab0b-canary04/`):** `runs/<unit>/verified.json` and `CANARY_APPROVED.json` [S3, S15].
- **C2:** `mtp-runs/P-358485d/{additive,archives,logs}/consumer-559ab0b-c2-canary04{d,e}/`. These are the evidence of record [S23].
- **Archive components:** all are `status: complete` [S24]:

| Component | Entries | Token |
|---|---|---|
| common | 1,115 | 04d |
| tier12-off | 2,609 | 04d |
| tier12-phenotype | 1,882 | 04d |
| paired | 2,410 | 04d |
| full-phenotype | 2,534 | 04d |
| full-off | 2,528 | 04d |
| canonical-phenotype | 1,670 | 04e |

  There are also native archives for each of the four units.
- **Read-only index:** `mtp-runs/C2-final-559ab0b-assembled-20261123/` holds 15,131 APFS clones (04d 11,492, 04e 3,639) and `ASSEMBLY_MANIFEST.tsv`; its verification found 0 hash mismatches. It is an index, not a C2 run [S23, S4].
- **How the two tokens link:** the coordination notes outside the repository, `c2-canary04e-live-20261122/TWO_TOKEN_RECORD.md` [S4].

## Outputs regenerated in this PR

The regeneration ran on the branch in a clean worktree, following REGEN_PLAN_v5 (= v4 plus the binding amendments). Its population is the canonical corpus: 5,629 human 16S/WMS profiles under the phenotype rollup. That is a different denominator from the cores' count profiles above, and the two must not be mixed. The Tier 1+2 taxa-coverage subset is the exception and is labelled as such.

At 830b23a6 the PR changes 543 tracked files against main: 523 modified, 19 added, 1 deleted [S6]; commit D adds this record and edits the Phase 5a files listed below. Of the 505 tracked files under `trait-visualization/trait-profile-investigation/`, 457 change [S6].

| Phase | Commit | Recipe(s) | Tracked files | Notes |
|---|---|---|---|---|
| 1 | 0bbac143 | `verify-pinned-kg`, `verify-trait-terms`, `study-trait-profiles`, `-refresh-gz`, `-check`, `-ora`, `analyze`, `analyze-depletion`, `contrast-orientation`, `trait-direction-screen`, `table1-compositional-stats`, `table1-log-ratio-sensitivity`, `calibration-report`; all ok [S12] | 1 (`study_trait_profiles_wide.tsv.gz`) | Phase 1 gate: the plan's gate covers 4 tables (profiles, matrix columns, studies, matrix); all 6 derived tables above are byte-identical to P's. The global tables go to gitignored paths. `screened_fraction_sign_test` fails family error (0.0832), as `methods/calibration.md` expects [S6]. |
| 2a | 37e25fa6 | `trait-heatmaps` | 17 under `trait-visualization/` | Content parity with C2 04e canonical metadata: profile 5,629, trait 60, paired 2,257 rows match, including the `*_sha256` columns [S6, S25]. Only the profile and paired metadata headers differ from C2 (`profile_id` / `pair_id` against blank); the trait metadata header is blank in both (#2396). |
| 2b | a01466da | `calibrate` (10 replicates, seed 42), `calibrate-composition` (100, seed 182), `audit-trait-density`, `audit-trait-coupling`, `trait-profile-pooled-stats` | 14 | Composition calibration: 932 studies, 2,257 pairs, 342 traits, 100/100 replicates [S26]. ORA calibration: status SUCCESS, 10/10 replicates, 4,514 profiles, 342 testable traits [S26]. Coupling report stamped `37e25fa6` [S26]. |
| 2c | 682aba35 | `cellulose-decr`, then `--repin` | 6 of 19 | Repinned to dd20b920, 6638849b, db8cdc96 and the profile metadata c9871a2a [S6]. Provenance `code_commit` `a01466da` [S26]. Experiment sets: PD 38, HIV 24, CRC 34; core 17, union 47 [S27]. |
| 2d | fa5d1240, da35d9e5 | `python -m mtp.native_influence_publication` from `P2-559ab0b-canary04/runs/full-phenotype/arms/noel/table1_taxon_influence` | 10 (`methods/table1-native-results/`) | Protocol `190-native-v2`, code `559ab0b8`, `code_dirty: false`. Profile input `da0189d4…` (the full-phenotype count-profile table). Contribution ledger: 4,514 profiles, 932 studies, 342 denominator columns, phenotype, 0 mismatches. 26 scenarios, 4 scopes, 1,248 composition and 1,248 ORA comparisons, 208 ORA families [S8]. da35d9e5 corrects the BH-family sentence in `README.md` by hand (owner choice F1 (i)); the frozen #1635 guard is unchanged [S6]. |
| 3 | 112e2cff, 4b3fef1f | `taxa-coverage --strain-rollup phenotype`; `taxa-coverage-548` | 14; 19 added | Canonical: 5,629 profiles [S28]. New Tier 1+2 subset directory `taxa-coverage-tier1-2-subject-controls/`: 1,870 of 5,629 selected, `disease_vs_healthy_tier1_2_subject_controls` [S28]. Median traits per mapped profile 156 → 154 [S6]. |
| 3 | d38e4bad | `trait-profile-investigation`, `oral-six-trait-heatmap`, the six `pd-*` recipes, at 315032f3; all exit 0 [S29] | 417 | aerobe-enrichment-v2 selects METPO:1000602 and METPO:1000604 (#2413) [S6]. |
| 3 | b57c82d3, 9f8d248e, 34fe02c3 | clean-tree re-stamps of the oral heatmap and the two taxa-coverage runs | 2, 1, 1 | The outputs reproduced; only provenance changes. Oral `source_artifact_commit` `d38e4bad`; taxa-coverage `code_commit` `b57c82d3`; the subset `9f8d248e` [S26, S28]. |
| 4 | 033472c4 | `integrated-evidence` (without `--influence-status`); `validation-evidence --report` | 6 (`data/derived/integrated_evidence/`) | `code_commit` `34fe02c3`; 4,514 kept profiles; taxon-influence status `pending #190: influence evidence not supplied` [S26]. No native run carries the canonical `dd20b920` digest that the influence status requires (v5 A11 F2) [S5]. The ignored `integrated_evidence.tsv` is not committed. The validation-evidence Markdown is kept in the phase notes [S6]. |
| test dispositions | e6e14733 (#2432), 689b01ed and 60cd4086 (#2433), 4b542d1e (#2434), 0e29ea54 (#749), 046bc47a | test and doc edits | 1 deleted + tests | e6e14733 retires `stratified_stats_provenance.json` (#2432, option 1). 689b01ed, 60cd4086 and 046bc47a re-pin tests to the regenerated values; 4b542d1e makes two tests hermetic. #749: the integration test is rebound to the full-off table `6fdafa21` (traits 317 → 316; profile pairs 50,086 → 49,770); its canary passes on 6fdafa21, skips when unset, and fails on dd20b920 [S6, S30]. |
| 5a | c27c93dd | prose reconciliation (inventories A and B) | 25 | About 480 changed values reviewed. The seven texts the review found stale are fixed (#2448). The native README is corrected 317 → 342 traits (#2449) [S6, S31]. |
| 5a | f959f81b, 40b69ea9 | detect-secrets exclusion and baseline refresh | 2 | The 122 new findings are sha256 or commit hashes; no new exclusions besides the taxa-coverage subset folder [S6]. |
| 5a | ac41537a, 830b23a6 | #1316 helper mapping; merge of main `a302aab5` | 2 | ac41537a maps the Tier 1+2 taxa-coverage provenance; PR #2452 (#2451) mapped the ten files already on main [S6, S33]. |
| 5a | this PR's data/doc commit D | this record; `just calibration-screen`; owner dispositions for four dated snapshots | 5 | `methods/calibration_screen.md`: Result replaced with a rerun at 830b23a6 (seed 182, 100 replicates, `phenotype`, m = 342, 923 studies; BY q <= 0.05 observed 27, no null call), the `0cc6c66` values kept as one dated sentence. `methods/table1_log_ratio_validation.md`, `docs/table1_trait_semantics_audit.md` and `methods/guild_condition_body_site_heatmap.md`: kept as dated history with a dated current-run note; the guild-maps section is to be rewritten in Phase 6a [S37]. |

**Checks on the branch.**
- **`just check` at 830b23a6:** 10,360 passed, 2 skipped, 9 deselected; detect-secrets and pip-audit passed [S32].
- **#1316 `compare` at 830b23a6 against baseline 83c90ae5:** status `NO_NEW_FAILURES_OR_SERIALIZATION_CHANGES`, as the plan requires before D [S33, S34]. All 211 record-schema pairs validate under both schemas and models and serialize equally; no constraint newly fails; no file is unclassified. The first run, at 40b69ea9, returned `REVIEW_REQUIRED` because 11 files had no reviewed mapping (#2451); ac41537a and PR #2452 mapped them.

## Not regenerated

- **Phase 6 (gitignored), built later on clean commit C.** Phase 5b's commit C sets `RERUN_COMMIT` in `scripts/generate_tufte_plates.py` and `COMMIT` in `scripts/build_slide_deck.py` to this PR's data/doc commit D. Phase 6 then runs on C [S34]:
  - 6a: `cellulose-decr-figures`, `guild-maps`, `table1-trait-maps`, `deck-figures`;
  - 6b: `figure-review-plates --overwrite`;
  - 6c: `figure-review-deck --overwrite`, which verifies cited files against `git ls-tree` at D.

  Don't cite plates or deck figures from main until then. D and C are to be tagged `regen-v5-D` and `regen-v5-C` and merged with a merge commit, so both stay reachable. Pushing the tags needs its own owner approval [S5, S35].
- **NCASM Tufte plates (`ncasm-tufte-plates`):** optional, not run unless the owner asks. They are gated on #2284 [S5].
- **`cellulose-decr-tufte`:** stays audit-only (retired). The plates are rebuilt from the new data as a follow-up after D (#1670) [S5].
- **The tier12-phenotype native influence run:** archived; not published into `methods/table1-native-results/` [S5].
- **Integrated-evidence influence status:** not supplied (see Phase 4). The influence row stays pending #190 [S26].
- **`taxa-oxygen-switches/`:** regenerated in d38e4bad, but held. Its construct counts the merged METPO:1000608 column, which now also counts taxa annotated METPO:1000605 or METPO:1000612, on the aerobic side, and calls it "facultatively aerobic". Its prose was not touched in Phase 5a. The PR stays in draft until #2437 is settled [S6, S35, S36].
- **The trait label "facultatively aerobic":** not renamed in prose (#2436, #2437, #2439). Nothing cites Tufte plate 4 [S6].
- **`trait-visualization/trait-profile-investigation/stratified_stats_provenance.json`:** deleted, not regenerated (#2432 option 1). A producer that keeps a record after a full render is a follow-up [S6].
- **Producer fixes found during the regeneration:** filed, and left to main (the branch changes no producer):
  - #2449: the native README "317-trait" literal;
  - #2450: the "959-study" docstrings;
  - #2435: the ontology-audit `_317`/`_301` labels;
  - #2396: the blank metadata header [S6, S36].
- **`methods/external_cohorts_sun2024.tsv`:** unchanged [S6].
- **Output-affecting PRs #2395, #2391 and #2388:** all three are open and unmerged as of 2026-11-27, so none is in this branch. No record shows an explicit owner hold (REGEN_PLAN_v5 A11 F8) [S36].

## Known issues

- **#2437 (open):** taxa-oxygen-switches counts the merged facultative METPO:1000608 column on the aerobic side and labels it "facultatively aerobic". It gates the PR leaving draft [S35, S36].
- **#2439 (open):** Tufte plates 2 and 7 and the methods plate name METPO:1000608 "Facultatively aerobic" in fixed text. Plates 1 and 4 show the merged label "Facultative (aerobe or anaerobe)" [S36].
- **#2440 (open):** Tufte plate 4 draws a missing condition × trait row as not significant, or skips it, without an error [S36].
- **#2444 (open):** plate 1 cannot check the trait-key policy of its Scale 1–4 tables (pooled, feces, Table 1 sign) [S36].
- **#2445 (open):** the plate 1 subtitle and callout nest parentheses around the merged facultative label [S36].
- **#2449 (open):** `native_influence_publication` hard-codes "317-trait" in the table1-native-results README; corrected by hand on this branch [S36, S31].
- **#2450 (open):** the `validation_evidence` docstrings hard-code the 959-study discovery corpus, which is 932 after #2397 [S36].
- **#2436 (closed), fixed by #2438 (merged):** the Tufte plates looked up METPO:1000608 by its raw label, so under `facultative_union_v1` the facultative column dropped silently. #2438 finds Table 1 traits by trait id in plates 1 and 4, and was merged into the branch at 98868b49 [S6, S36].
- **#749 (open):** the integration-marked trait-screen test was stale. It is rebound to the full-off table here (0e29ea54) [S6, S36].
- **#2432 (open):** the full render deletes `stratified_stats_provenance.json`, while a test and `docs/workflow.md` pinned the committed copy. Option 1 was applied here (e6e14733) [S6, S36].
- **#2433 (open):** five tests pinned values of local production data that this regeneration replaces. They are re-pinned in 689b01ed and 60cd4086 [S6, S36].
- **#2434 (open):** tests that look hermetic read the local `data/derived` and fail once it holds regenerated outputs. Two are fixed in 4b542d1e; the third named in the issue (Tufte plate 2 "/341") was reclassified as a #2433 case and re-pinned in 60cd4086 [S6, S36].
- **#2451 (closed), fixed by PR #2452 (merged):** the #1316 compare helper had no mapping for 11 tracked files [S36].
- **#1363 (open):** whether to retire the 522 tracked generated outputs. For this regeneration they are still committed (owner decision 1), because Phase 6c verifies cited files against git trees at D [S5, S36].
- **Also open:**
  - #2448: the seven stale texts, fixed in c27c93dd;
  - #2435: the ontology-audit labels;
  - #2396: the blank metadata header;
  - #2284: NCASM plate 21 [S36].

## Sources

Paths are outside this repository unless they start with W. K is the directory that holds the repository checkout and the run directories; N = K/mtp-paper-rerun/notes/final-run-downstream-20261113 (coordination notes); W = the regeneration branch worktree.

- S1: `K/mtp-runs/P-358485d/runs/{tier12-off,tier12-phenotype,full-phenotype,full-off}/verified.json` and `{full-phenotype,full-off}/verified.amended.json`
- S2: `K/mtp-runs/P-358485d/runs/full-phenotype/core/run-c06d9783cd8a41c9bc904cfd5abf55a9/run.json` (`lock_sha256`, `git_dirty`); `K/mtp-runs/P-358485d/runs/<unit>/core/run-*/run.json` for all four units (`parameters[].source_mode`)
- S3: `K/mtp-runs/P2-559ab0b-canary04/runs/<unit>/verified.json` (+ `DONE`)
- S4: `K/mtp-paper-rerun/notes/c2-canary04e-live-20261122/TWO_TOKEN_RECORD.md`
- S5: `N/REGEN_PLAN_v5.md` (A1–A12); `N/PHASE0_RECHECK_20261124.md`
- S6: `git log origin/main..HEAD` in W (commit messages 0bbac143 … 830b23a6), `git diff --name-status origin/main HEAD`, and sha256 of W's `data/derived/*`, `data/raw/full_dump.csv`, `data/kg-microbe/*` and `uv.lock`, computed 2026-11-27
- S7: W `scripts/adhoc/issue1316_validate_committed.py` (diff vs origin/main)
- S8: W `methods/table1-native-results/{run_manifest,contribution_manifest,publication_manifest}.json`; `N/phase2-20261124/2d_publish.log`
- S9: `K/mtp-runs/P-358485d/preflight/{inputs,copied-inputs}.json`
- S10: W `data/kg-microbe/provenance.json`
- S11: `K/mtp-paper-rerun/STATUS.md`, row 2026-11-24 21:25 (inputs cloned from P-358485d/wt)
- S12: `N/phase1-20261124/{01_verify-pinned-kg.log,phase1a.out,phase1b.out}`
- S13: `K/mtp-runs/P-358485d/preflight/experiments.json`; sha256 of `K/mtp-runs/P-358485d/wt/data/derived/*`, computed 2026-11-27
- S14: `K/mtp-runs/P-358485d/preflight/profiles-{phenotype,off}.json`
- S15: `K/mtp-runs/P2-559ab0b-canary04/CANARY_APPROVED.json`
- S16: `K/mtp-runs/C2-final-559ab0b-assembled-20261123/records/LEAD_REVIEW_cont03_paired_GATE.md`
- S17: `N/PHASE0_RECHECK_20261124.md` ("Facts for the owner decisions"); REGEN_PLAN_v5 A4
- S18: `ls K/mtp-runs/C2-final-559ab0b-assembled-20261123/additive/*`
- S19: `…/records/LEAD_REVIEW_cont01_additive-tier12-off.md`
- S20: `…/records/LEAD_REVIEW_cont05_additive-full-phenotype.md`
- S21: `…/records/LEAD_REVIEW_cont07_additive-full-off.md`
- S22: `…/records/LEAD_REVIEW_FINAL_ACCEPTANCE.md`; STATUS.md row 2026-11-23 12:10
- S23: `K/mtp-runs/C2-final-559ab0b-assembled-20261123/README.md`
- S24: `K/mtp-runs/C2-final-559ab0b-assembled-20261123/archives/*/archive.json` (`status`, `attempt`, `files` count)
- S25: `N/phase2-20261124/parity_2a.out`
- S26: W `data/derived/calibration/composition/composition_calibration_summary.json`, `data/derived/ora_calibration/calibration_summary.json`, `data/derived/trait_coupling_report.txt`, `data/analysis/cellulose_decr/cellulose_decr_provenance.json`, `data/derived/integrated_evidence/integrated_evidence_provenance.json`, `six-shared-traits/oral/oral_six_trait_provenance.json`
- S27: `N/phase2-20261124/2c_cellulose-decr_repin.log`
- S28: W `trait-visualization/trait-profile-investigation/{taxa-coverage,taxa-coverage-tier1-2-subject-controls}/provenance.json`
- S29: `N/phase3-20261124/run_phase3_315032f3.summary.txt`
- S30: `N/phase5a-20261126/issue749_canary.txt`
- S31: `N/phase5a-20261126/{issue_prose_review.md,issue_native_317.md,issue_ve_959.md}`
- S32: `N/phase5a-20261126/issue1316/just_check_merge_a302aab5.log`
- S33: `N/phase5a-20261126/issue1316/{compare_40b69ea9.json,compare_830b23a6.json}`
- S34: `N/REGEN_PLAN_v4.md` (Phases 5a, 5b, 6)
- S35: PR #2398 (`gh pr view 2398`, head 830b23a6, draft)
- S36: `gh issue view` for #2437, #2439, #2440, #2444, #2445, #2449, #2450, #2436, #2438, #749, #2432, #2433, #2434, #2435, #1363, #2396, #2448, #2413, #2284, #2451, #2395, #2391, #2388 (titles and states, 2026-11-27)
- S37: W `data/derived/calibration/screen/screen_null_{summary.tsv,provenance.json}` (gitignored; `code_commit` 830b23a6); W `data/derived/table1_composition/trait_log_ratio_stats.tsv`; W `trait-visualization/trait-profile-investigation/six-shared-traits/six_shared_trait_summary.tsv`; W `data/derived/integrated_evidence/integrated_evidence_{provenance.json,method_results.tsv}`
