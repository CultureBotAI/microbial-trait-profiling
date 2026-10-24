# Final rerun record: commit 511c891 (2026-10-23)

This records the final rerun before the NCASM talk: the runs made, the commit and inputs they used, their population and strain-rollup scope, what was verified, and which tracked outputs this PR regenerates. A number cited as a core-bound result must come from one of the verified runs below. The tracked outputs this PR regenerates are built from hashed legacy inputs and aren't core-verified, as their own provenance says. Cite them by file, together with the verification status their provenance records. Whichever source a number comes from, state its population and scope with it.

## Code and inputs

- **Commit:** F = `511c891229b659e207cfe0f4281a97c3fc039a80`, main after the pre-rerun merges, ending with #891.
  - Every core, arm and comparison record below reports `git_commit` F and `git_dirty: false`.
  - The downstream pass's provenance records differ. The taxa-coverage record gives the commit as `511c891…+dirty`: the producer's dirty check sees the pass's own regenerated outputs in the working tree, and its recorded script hash is F's. The site-stats record carries no commit, only producer and input hashes (#1017).
- **Environment:** the locked environment at F; `uv.lock` is unchanged.
  - CI's pip-audit flags four advisories in filelock 3.25.2 and setuptools 84.0.0 that were published after F (#999).
  - The lock bump (#1000) merges right after this PR.
- **KG-Microbe source:** pinned (`source_mode: pinned`); `just verify-pinned-kg` passed.
- **Corpus:**
  - `data/derived/experiments.tsv` has sha256 `16271237a70f61c8…`;
  - the downstream pass's profile table, `data/derived/study_trait_profiles_wide.tsv`, has sha256 `43b2e2302b0edf48…` and 2,637 columns. It holds the 5,629 human 16S/WMS profiles under the phenotype rollup. Each core builds its own count profiles instead.

## Verified core-bound runs

Each of the four units ran these steps, and every step exited 0. The last, full corpus × off, finished at 18:49 PDT.
1. `just rerun-core`;
2. the three arms (ilyas, jed, noel) in parallel;
3. `just rerun-verify` on each arm;
4. `just rerun-compare`;
5. `just rerun-verify-comparison`.

Beyond the exit codes, the side effects were checked: the core's profile, experiment and study counts; the trait count and rollup that each analysis basis records; and non-empty outputs.

| Unit | Contrast population | Rollup | Profiles | Experiments | Studies | Traits |
|---|---|---|---|---|---|---|
| Tier 1+2 × off (the talk's primary unit; canary) | `disease_vs_healthy_tier1_2_subject_controls` | off | 1,870 | 935 | 548 | 317 |
| Tier 1+2 × phenotype | `disease_vs_healthy_tier1_2_subject_controls` | phenotype | 1,870 | 935 | 548 | 344 |
| Full corpus × phenotype | `all_paired` | phenotype | 7,111 | 4,299 | 1,212 | 344 |
| Full corpus × off | `all_paired` | off | 7,111 | 4,299 | 1,212 | 317 |

Notes on the counts:
- **What is counted:** the profiles, experiments and studies represented in each core's count-profile product. A profile merges one experiment's signatures in one direction, increased or decreased.
- **Tier 1+2:** every one of the 935 experiments has both profiles, so they are 935 paired contrasts.
- **Full corpus:** the 7,111 profiles come from 4,299 experiments, and some experiments have a signature in one direction only. So 4,299 is not a count of paired contrasts. The full experiments table is larger, with 4,435 rows in 1,245 studies.
- **Hosts:** the full-corpus core counts span all host species (`count_host_species: all`, the default). Over-representation analysis uses *Homo sapiens* (`ora_host_species`).

| Unit | Core run | Arm runs (ilyas, jed, noel) | Comparison |
|---|---|---|---|
| Tier 1+2 × off | `run-09f23bbe33cb42988413bec7caae8623` | `arm-edf15610afce44e3a6e593b6d541c90d`, `arm-114909a09cee4a7095a78afc1aad219b`, `arm-2afce27ed61c4814b737984dd5b93aa8` | `comparison-disease_vs_healthy_tier1_2_subject_controls-f68bdcc0405d` |
| Tier 1+2 × phenotype | `run-ec26436944744c4cbb7f653c16a19287` | `arm-16d11bd4c4b2464daf07bfe6c75ab783`, `arm-1e5709177a3c4dd2940b709176267748`, `arm-b3a766ad095147c0b6c837411096cd3d` | `comparison-disease_vs_healthy_tier1_2_subject_controls-4e1ff7b17579` |
| Full corpus × phenotype | `run-bae7645a0d1243b6aa1da31228f5f014` | `arm-3b3cf87b84464be38712b28bc5ca5c6f`, `arm-27edba337abb4c60917b896bd2ca92b1`, `arm-7ffddf6556d04acfbe43709da6a0f624` | `comparison-all_paired-09666a6843f5` |
| Full corpus × off | `run-a7cf841adbc34221a5673b75f6e2ec3f` | `arm-c458ef1659354ad88d66deb458c91d87`, `arm-8f6a921d76b646c9a59e1bf9026a1715`, `arm-c01764d0e9184cfab3970ba1a112a2e2` | `comparison-all_paired-79968d18fa2e` |

**Calibration verdicts.** All four units report the same two verdicts:
- `noel_table1_six_trait_sign_test_v1`: partially calibrated. The null-control criterion passed. Power was reported (its `power_reported` criterion passed), but no target-power criterion was declared, so the comparison does not certify power. The verdict also records two things: the six Table 1 traits were selected on this corpus before testing, so their p-values are exploratory; and the select-then-test negative control fails, as the protocol declares it will.
- `jed_log_ratio_vs_transformed_t_calibration`: different estimator, uncalibrated.

**Reproducing these values.** `scripts/adhoc/extract_final_rerun_record_values.py --runs-root <runs> --additive-root <Tier 1+2 × off additive dir>` reads the already-verified run directories and prints, as JSON, the run IDs, core counts, scopes, verdict outcomes and the scan, lineage and label-support figures quoted here. It extracts selected values; it isn't an independent verifier. The experiments-table totals and the detailed calibration criteria come from the run records directly.

**Run directories.** They take 4.8 to 6.5 GB each and aren't committed; they're kept outside the repository. To keep them, make a relocatable copy with `just evidence-archive` (see [Retaining execution evidence](../evidence-archives.md)).

## Additive analyses on Tier 1+2 × off

Tier 1+2 × off is the talk's primary unit, and the scan protocol declares off as its primary scope. The repository's primary scope for tracked outputs is phenotype (#788).

These analyses ran at F against core run `run-09f23bbe33cb42988413bec7caae8623`:
- **Bound to the run:** the scan and the lineage extension, through their core-consumption records.
- **Matched by hash only:** the label-support audit reads that core's input files, matched by hash, but isn't bound to the run itself.

Their outputs are separately named sidecars and aren't committed here. Every p-value below is nominal.

### Condition-triple selection scan (#861; for #835)
- **Setup:** protocol `cellulose-decr-selection-scan-v1`, 100,000 U6 study-shuffle permutations, seed 20261022; 34 eligible conditions (at least 3 stool studies), 5,984 triples and 10 traits.
- **How triples are ranked:** each rank below is the triple's rank among the 5,984 by its U6 p-value for the named statistic, not by the observed statistic itself. For example, the rule-chosen triple has the smaller observed distance (0.538 against 0.542), yet it ranks lower because its p-value is larger.
- **PD + HIV + colorectal cancer** (cellulose catabolism): rank 1,368 of 5,984 by U6 p-value for mean pairwise Jaccard distance.
  - Observed sharing is 45.8%, against 43.5% expected under U6 (p = 0.35).
  - By U6 p-value for the three-way core, it ranks 3,446 (p = 0.82).
- **The triple chosen by rule** (the three conditions with the most stool studies: PD 36, obesity 25, colorectal cancer 24): rank 3,873, p = 0.73.
- **The lowest-p triple** (chronic kidney disease + gastric cancer + HIV, for cellulose catabolism): p = 0.00024, the lowest U6 p for mean pairwise Jaccard distance among all 59,820 tested trait-triples in the scan. Gastric cancer has 3 stool studies.

### Lineage extension and study leave-one-out (#891)
- **Status:** protocol `771-extension-v1`, an exploratory sensitivity that leaves the declared `190-native-v2` panel unchanged. Reconciliation mismatches: 0.
- **Effect removed:** the share of the study-mean shared-fraction cellulose effect removed by dropping each lineage (pooled; matched in parentheses):

| Lineage removed | Effect removed |
|---|---|
| Lachnospiraceae | 57.1% (54.1%) |
| Oscillospiraceae | 25.3% (30.1%) |
| Both | 91.2% (90.6%) |
| *Ruminococcus* | 15.9% (19.0%) |
| *Acetivibrio* | 12.2% (12.2%) |
| *Ruminiclostridium* | 10.3% (10.3%) |
| Clostridiaceae | −8.8% (−9.3%) |

### Label-support audit (#877; for #773)
- **Run:** on the Tier 1+2 × off inputs; all 1,870 profiles reconciled.
- **Verification status:** `hashes_and_counts_reconciled_not_core_execution_verified`.
- **Madin table:** not supplied, so Madin support is reported as unknown. Whether to pin the Madin table is still the owner's decision (#773).

## Outputs regenerated in this PR

A downstream pass ran at F with the phenotype rollup, the repository's primary scope for tracked outputs. Its population is the canonical analysis corpus: 5,629 human 16S and shotgun (WMS) profiles, with no Tier 1+2, pairing or disease-orientation filter. This is a different denominator from the cores' counts above, and the two must not be mixed. The steps, in order:
1. ingest and the study TSVs;
2. `just study-trait-profiles`, then `just study-trait-profiles-refresh-gz`;
3. `just trait-heatmaps`;
4. the pooled stats, then `just trait-profile-site-stats`;
5. `just contrast-orientation`;
6. `just cellulose-decr --repin`, then `just cellulose-decr-figures`;
7. `just taxa-coverage --strain-rollup phenotype`.

Every listed step exited 0. The profile check before the gz refresh exits 1 by design, because the committed gz still held the off table. `just integrated-evidence` also ran (exit 0), and its outputs were discarded (see "Not regenerated"). No file this PR adds or changes contains a local path. Elsewhere in the repository, the committed integrated_evidence files still do (#1001).

What changed:
- **Profile table:** `data/derived/study_trait_profiles_wide.tsv.gz` now holds the phenotype profile table; it held the off table before. Its `.strain_rollup.json` record arrives with #936.
- **Taxa coverage:** the outputs are rebuilt from the canonical 5,629 profiles (#820).
- **New files:**
  - `data/analysis/cellulose_decr/cellulose_decreased_taxa_catalog_tier1_2.tsv`, a new output of `cellulose-decr`;
  - the provenance records that `trait-profile-site-stats` and `taxa-coverage` write at F.

## Not regenerated

- **`data/derived/integrated_evidence/`:** stays at its committed version.
  - Its legacy inputs weren't regenerated in this pass: the calibration status, the composition sensitivity, the Arm 2 enrichment and depletion tables, and the full trait screen.
  - Without them the producer recomputes only the arm method results. The calibration and sensitivity statuses come out pending, and the enrichment, aggregated-support and guild-cell rows come out unsupplied (#1003).
  - It also writes absolute paths (#1001).
  - The committed integrated_evidence no longer matches the inputs this PR replaces. Its provenance names the off profile table (`519a095b…`, 317 trait columns) and Arm 1 pooled stats `3f89fca8…`. The committed profile `.gz` now decompresses to `43b2e230…`, and both committed pooled copies are now `46dfce6b…`. Regenerate it with #1003 before citing it next to the phenotype outputs.
- **Eight per-body-site artifacts:** `just trait-profile-site-stats` refreshed seven artifacts at F. It leaves eight from an earlier run, because its pooled/site statistics mode doesn't write them:
  - the caecum, colon, mouth, subgingival dental plaque and vagina statistics;
  - the condition and body-site summaries;
  - the infographic PNG.

  `trait-visualization/trait-profile-investigation/stratified_stats_provenance.json` lists both groups, with hashes.
- **Other tracked investigation outputs:** recipes this pass didn't run leave their tracked outputs at their committed versions. Examples are `just trait-profile-investigation` and the six-shared-trait recipes. This PR changes 29 of the 487 tracked files under `trait-visualization/trait-profile-investigation/`: 27 regenerated outputs and 2 refreshed notes.
- **Outputs still built from the off profile table or matrix** (`519a095b…` / `f6de0576…`, mostly at 8b280a3, #766). Don't combine them with the phenotype outputs above:
  - the composition calibration summary;
  - `data/derived/trait_coupling_report.txt` and its companion TSVs;
  - the density-adjustment redundancy summary;
  - the ORA calibration summary;
  - the PD-cellulose evidence and manifests. Their input guards will refuse the new phenotype tables until those outputs are regenerated (#1016).
- **The Tufte plates and deck figures:** the talk's figures come from the slide pipeline and cite the runs above. Don't rebuild the plates or the deck from main until #1015 lands. `scripts/build_slide_deck.py` and `scripts/generate_tufte_plates.py` read the outputs this PR changes but hard-code pre-rerun numbers, so a rebuild would draw phenotype data under off-era labels.

## Known issues

- **#999:** pip-audit is the only failing check, by owner exception; #1000 follows this PR.
- **#1001, #1003:** integrated evidence (see above).
- **#821:** the taxa-coverage recipe isn't in the regeneration list.
- **#940:** readers that combine the profile table and the matrix don't check their strain-rollup records.
- **#1015:** the deck and plate builders hard-code pre-rerun coverage numbers and a "317 traits" caption.
- **#1016:** tracked outputs still built from the off profile table.
