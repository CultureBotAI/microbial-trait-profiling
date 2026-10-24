# Final rerun record: commit 511c891 (2026-10-23)

This records the final rerun before the NCASM talk: the runs made, the commit and inputs they used, their population and strain-rollup scope, what was verified, and which tracked outputs this PR regenerates. A number cited from this rerun must come from one of the verified runs below. Its population and scope must be stated with it.

## Code and inputs

- **Commit:** F = `511c891229b659e207cfe0f4281a97c3fc039a80`, main after the pre-rerun merges, ending with #891. Every run record below reports `git_commit` F and `git_dirty: false`.
- **Environment:** the locked environment at F; `uv.lock` is unchanged.
  - CI's pip-audit flags four advisories in filelock 3.25.2 and setuptools 84.0.0 that were published after F (#999).
  - The lock bump (#1000) merges right after this PR.
- **KG-Microbe source:** pinned (`source_mode: pinned`); `just verify-pinned-kg` passed.
- **Corpus:**
  - `data/derived/experiments.tsv` has sha256 `16271237a70f61c8…`;
  - the profile table built at F, `data/derived/study_trait_profiles_wide.tsv`, has sha256 `43b2e2302b0edf48…` and 2,637 columns.

## Verified core-bound runs

Each unit ran these steps, and every step exited 0:
1. `just rerun-core`;
2. the three arms (ilyas, jed, noel) in parallel;
3. `just rerun-verify` on each arm;
4. `just rerun-compare`;
5. `just rerun-verify-comparison`.

Beyond the exit codes, the side effects were checked: the core's profile, contrast and study counts; the trait count and rollup that each analysis basis records; and non-empty outputs.

| Unit | Contrast population | Rollup | Profiles | Contrasts | Studies | Traits |
|---|---|---|---|---|---|---|
| Tier 1+2 × off (primary; canary) | `disease_vs_healthy_tier1_2_subject_controls` | off | 1,870 | 935 | 548 | 317 |
| Tier 1+2 × phenotype | `disease_vs_healthy_tier1_2_subject_controls` | phenotype | 1,870 | 935 | 548 | 344 |
| Full corpus × phenotype | `all_paired` | phenotype | 7,111 | 4,299 | 1,212 | 344 |
| Full corpus × off (verification in progress) | `all_paired` | off | pending | pending | pending | pending |

Notes on the counts:
- The full-corpus core counts span all host species (`count_host_species: all`, the default).
- Over-representation analysis uses *Homo sapiens* (`ora_host_species`).

| Unit | Core run | Arm runs (ilyas, jed, noel) | Comparison |
|---|---|---|---|
| Tier 1+2 × off | `run-09f23bbe33cb42988413bec7caae8623` | `arm-edf15610afce44e3a6e593b6d541c90d`, `arm-114909a09cee4a7095a78afc1aad219b`, `arm-2afce27ed61c4814b737984dd5b93aa8` | `comparison-disease_vs_healthy_tier1_2_subject_controls-f68bdcc0405d` |
| Tier 1+2 × phenotype | `run-ec26436944744c4cbb7f653c16a19287` | `arm-16d11bd4c4b2464daf07bfe6c75ab783`, `arm-1e5709177a3c4dd2940b709176267748`, `arm-b3a766ad095147c0b6c837411096cd3d` | `comparison-disease_vs_healthy_tier1_2_subject_controls-4e1ff7b17579` |
| Full corpus × phenotype | `run-bae7645a0d1243b6aa1da31228f5f014` | `arm-3b3cf87b84464be38712b28bc5ca5c6f`, `arm-27edba337abb4c60917b896bd2ca92b1`, `arm-7ffddf6556d04acfbe43709da6a0f624` | `comparison-all_paired-09666a6843f5` |
| Full corpus × off | `run-a7cf841adbc34221a5673b75f6e2ec3f` | pending | pending |

**Calibration verdicts.** The verified units report the same two verdicts:
- `noel_table1_six_trait_sign_test_v1`: partially calibrated. The null control passed; power was not assessed.
- `jed_log_ratio_vs_transformed_t_calibration`: different estimator, uncalibrated.

**Run directories.** They take 4.8 to 6.5 GB each and aren't committed; they're kept outside the repository. To keep them, make a relocatable copy with `just evidence-archive` (see [Retaining execution evidence](../evidence-archives.md)).

## Additive analyses on the primary unit (Tier 1+2 × off)

These ran at F on core run `run-09f23bbe33cb42988413bec7caae8623`. Their outputs are separately named sidecars and aren't committed here. Every p-value below is nominal.

### Condition-triple selection scan (#861; for #835)
- **Setup:** protocol `cellulose-decr-selection-scan-v1`, 100,000 U6 study-shuffle permutations, seed 20261022; 34 eligible conditions (at least 3 stool studies), 5,984 triples and 10 traits.
- **PD + HIV + colorectal cancer** (cellulose catabolism): rank 1,368 of 5,984 by mean pairwise Jaccard distance.
  - Observed sharing is 45.8%, against 43.5% expected under U6 (p = 0.35).
  - By three-way core, it ranks 3,446 (p = 0.82).
- **The triple chosen by rule** (the three conditions with the most stool studies: PD 36, obesity 25, colorectal cancer 24): rank 3,873, p = 0.73.
- **The scan's top triple** (chronic kidney disease + gastric cancer + HIV): p = 0.00024. Gastric cancer has 3 stool studies.

### Lineage extension and study leave-one-out (#891)
- **Status:** protocol `771-extension-v1`, an exploratory sensitivity that leaves the declared `190-native-v2` panel unchanged. Reconciliation mismatches: 0.
- **Effect removed:** the share of the study-mean shared-fraction cellulose effect removed by dropping each lineage (pooled; matched in parentheses):

| Lineage removed | Effect removed |
|---|---|
| Lachnospiraceae | 57.1% (54.1%) |
| Oscillospiraceae | 25.3% (30.1%) |
| Both | 91.2% (90.6%) |
| *Ruminococcus* | 15.9% |
| *Acetivibrio* | 12.2% |
| *Ruminiclostridium* | 10.3% |
| Clostridiaceae | −8.8% (−9.3%) |

### Label-support audit (#877; for #773)
- **Run:** on the Tier 1+2 × off inputs; all 1,870 profiles reconciled.
- **Verification status:** `hashes_and_counts_reconciled_not_core_execution_verified`.
- **Madin table:** not supplied, so Madin support is reported as unknown. Whether to pin the Madin table is still the owner's decision (#773).

## Outputs regenerated in this PR

A downstream pass ran at F on the full corpus with the phenotype rollup, the primary scope for tracked outputs. The steps, in order:
1. ingest and the study TSVs;
2. `just study-trait-profiles`, then `just study-trait-profiles-refresh-gz`;
3. `just trait-heatmaps`;
4. the pooled stats, then `just trait-profile-site-stats`;
5. `just contrast-orientation`;
6. `just cellulose-decr --repin`, then `just cellulose-decr-figures`;
7. `just taxa-coverage --strain-rollup phenotype`.

Every step exited 0. No committed file contains a local path.

What changed:
- **Profile table:** `data/derived/study_trait_profiles_wide.tsv.gz` now holds the phenotype profile table; it held the off table before. Its `.strain_rollup.json` record arrives with #936.
- **Taxa coverage:** the outputs are rebuilt from the canonical 5,629 profiles (#820).
- **New files:**
  - `data/analysis/cellulose_decr/cellulose_decreased_taxa_catalog_tier1_2.tsv`, a new output of `cellulose-decr`;
  - the provenance records that `trait-profile-site-stats` and `taxa-coverage` write at F.

## Not regenerated

- **`data/derived/integrated_evidence/`:** stays at its committed version.
  - Its legacy inputs weren't regenerated in this pass: the calibration status, the composition sensitivity, the Arm 2 enrichment and depletion tables, and the full trait screen.
  - Without them the producer writes every row as pending (#1003).
  - It also writes absolute paths (#1001).
- **The Tufte plates and deck figures:** the talk's figures come from the slide pipeline and cite the runs above.

## Known issues

- **#999:** pip-audit is the only failing check, by owner exception; #1000 follows this PR.
- **#1001, #1003:** integrated evidence (see above).
- **#821:** the taxa-coverage recipe isn't in the regeneration list.
- **#940:** readers that combine the profile table and the matrix don't check their strain-rollup records.
