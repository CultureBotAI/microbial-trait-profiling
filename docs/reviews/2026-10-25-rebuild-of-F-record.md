# Rebuild of run F at commit 511c891 (2026-10-25)

This record supersedes the evidence locations in [the 2026-10-23 final rerun record](2026-10-23-final-rerun-record.md) (#1040). That record's values, scopes and verdicts stand. Its run directories do not: they lived only in temporary storage, which a restart emptied on 2026-10-24, and no archive had been made. This record documents the rebuild of all four units at the same commit, from the same pinned inputs, and the verified, relocatable archives that now hold them.

## What changed and what did not

- **Same:** the commit, the pinned inputs, the locked environment, the recipes, the populations and scopes, every count, both calibration verdicts, and every value checked below.
- **New by construction:** every core, arm and comparison run ID.
- **New:** four verified evidence archives (`just evidence-archive`, then `just evidence-verify`), kept in a durable location outside the repository and outside temporary storage. This record names that location symbolically as `<F-rebuild root>`.

No value in the 2026-10-23 record differs from the rebuild.

## Code and inputs

- **Commit:** F = `511c891229b659e207cfe0f4281a97c3fc039a80`, a detached worktree. Every core, arm and comparison record reports `git_commit` F and `git_dirty: false`.
- **Environment:** `uv sync --locked` at F; F's `uv.lock`, unchanged.
- **Pinned inputs:**
  - the BugSigDB v1.3.1 export, `full_dump.csv`, sha256 `af84bfca7bfe7691…`;
  - KG-Microbe `merged-kg_edges.tsv`, sha256 `d13812a436bd0b1c…`;
  - KG-Microbe `merged-kg_nodes.tsv`, sha256 `8b4cda73ca7e9d02…`.

  `uv run python -m mtp.contribution_inputs` (the pinned-KG check) passed.

### Canary 0: the corpus rebuilt from the pinned dump

Ingest and `mtp.study_tsv` at F reproduce the original tables byte for byte:
- `data/derived/experiments.tsv`: sha256 `16271237a70f61c8…`, the same as the 2026-10-23 record;
- the phenotype profile table, `study_trait_profiles_wide.tsv`: sha256 `43b2e2302b0edf48…`, 2,637 columns and 5,629 rows, the table #1004 tracked as `data/derived/study_trait_profiles_wide.tsv.gz`;
- the off profile table: sha256 `519a095b09e179dc…`, the table tracked before #1004.

## Rebuilt core-bound runs

Each unit ran the same steps as on 2026-10-23, and every step exited 0:
1. `just rerun-core`;
2. the three arms in parallel;
3. `just rerun-verify` on each arm;
4. `just rerun-compare`;
5. `just rerun-verify-comparison`.

Tier 1+2 × off ran first, as the canary; the other units ran only after its counts and values matched.

| Unit | Rollup | Profiles | Experiments | Studies | Traits |
|---|---|---|---|---|---|
| Tier 1+2 × off | off | 1,870 | 935 | 548 | 317 |
| Tier 1+2 × phenotype | phenotype | 1,870 | 935 | 548 | 344 |
| Full corpus × phenotype | phenotype | 7,111 | 4,299 | 1,212 | 344 |
| Full corpus × off | off | 7,111 | 4,299 | 1,212 | 317 |

These match the 2026-10-23 record exactly. **Calibration verdicts**, the same for all four units, as before:
- `noel_table1_six_trait_sign_test_v1`: partially calibrated;
- `jed_log_ratio_vs_transformed_t_calibration`: different estimator, uncalibrated.

| Unit | Core run | Arm runs (ilyas, jed, noel) | Comparison |
|---|---|---|---|
| Tier 1+2 × off | `run-a2d30f9ca1dd45aeb98085f39badc4bf` | `arm-28cbb0883ea14628b12af3b9ec742f8d`, `arm-65728f8a910f4b1189c1cf31db848004`, `arm-f53d3a2215e5445cb106d0e0bdb44c0b` | `comparison-disease_vs_healthy_tier1_2_subject_controls-3d40e23048d0` |
| Tier 1+2 × phenotype | `run-88a731a42ddf47fd95ac2277726a25d9` | `arm-5d06ccb3dd534d4ea336f539c7cc9c1c`, `arm-830ec567172a45d485ab27bb75908832`, `arm-a782b2dfa4b54520b28580d81e2c2107` | `comparison-disease_vs_healthy_tier1_2_subject_controls-2419355910b9` |
| Full corpus × phenotype | `run-2a04494912af4ed8a7ee15e1bf1f10a3` | `arm-a83d9e648b604ee4b57e28392974decc`, `arm-fbf8be7e0120471184b355e5a94f3787`, `arm-d2658609c8594147b51fc6323b4c11a2` | `comparison-all_paired-45100d862968` |
| Full corpus × off | `run-31b57ed4ba364ea3a9f4004be36a3edd` | `arm-477cd4a46ab049adaf546450e48944c0`, `arm-96f555bc60e544199c98adc4ad52ebc4`, `arm-46f8797ee5b548ad86b4927fc561e538` | `comparison-all_paired-c8c75206f020` |

## Additive analyses

They ran at F against the rebuilt Tier 1+2 cores. Every p-value is nominal.

### Tier 1+2 × off: every value in the 2026-10-23 record reproduced

`scripts/adhoc/extract_final_rerun_record_values.py`, run on the rebuilt directories, prints the same values the 2026-10-23 record quotes:

- **Selection scan** (protocol v1; 100,000 permutations; seed 20261022; 34 eligible conditions; 5,984 triples; 10 traits):
  - PD + HIV + colorectal cancer, for cellulose catabolism, ranks 1,368, p = 0.3515 (mean distance); for the three-way core it ranks 3,446, p = 0.822.
  - The rule-chosen triple ranks 3,873, p = 0.730.
  - The lowest-p triple (chronic kidney disease + gastric cancer + HIV) has p = 0.00024, among 59,820 tested trait-triples.
- **Lineage extension**, pooled (matched): Lachnospiraceae 57.1% (54.1%), Oscillospiraceae 25.3% (30.1%), both 91.2% (90.6%), *Ruminococcus* 15.9% (19.0%), *Acetivibrio* 12.2%, *Ruminiclostridium* 10.3%, Clostridiaceae −8.8% (−9.3%). There are 0 reconciliation mismatches.
- **Label support:** 1,870 profiles and 19,388 distinct occurrences, 0 unmapped; status `hashes_and_counts_reconciled_not_core_execution_verified`; Madin not supplied.

### Tier 1+2 × phenotype

These outputs are new. The 2026-10-23 record had none for this unit.
- **Selection scan:** the lead row is the same as off: rank 1,368, p = 0.3515; core rank 3,446, p = 0.822. A phenotype scan of the same recipe on the original Tier 1+2 × phenotype core was seen on 2026-10-24, outside the recorded evidence; protocol v2 (#1050) records this.
  - The rebuilt scan's condition, trait, triple and distribution tables are byte-identical to that run's.
  - Its summary differs only in the core run ID.
- **Lineage extension**, pooled (matched): Lachnospiraceae 56.3% (53.4%), Oscillospiraceae 28.0% (32.9%), both 91.0% (90.5%); 0 reconciliation mismatches.
- **Label support:** 1,870 profiles and 19,388 distinct occurrences, 0 unmapped; Madin not supplied.

Since the 2026-10-23 record, the selection scan's declared primary scope has changed from off to phenotype, through protocol v2 (#1050). That amendment was made after both scopes' results had been seen. Both scopes are reported side by side.

## Evidence archives

Each unit's archive retains:
- the core run, with its raw and KG inputs;
- every consumer record and its outputs: the arms' and comparison's consumer records, and, for the Tier 1+2 units, the scan and the lineage extension;
- the arm, comparison, label-support and log directories, retained explicitly.

`just evidence-verify` re-checked the retained bytes, the core contracts and every consumer chain for each archive.

| Unit | Archive | `archive.json` sha256 | Consumer records |
|---|---|---|---|
| Tier 1+2 × off | `archive:89988fe5-1814-432a-bb7a-41e72a713e79` | `424248ada34ec76b…` | 12 |
| Tier 1+2 × phenotype | `archive:53c358dc-0c4c-48cb-b9b0-927f542f0834` | `16a0804caf056760…` | 12 |
| Full corpus × phenotype | `archive:b9e0fa87-45f2-4973-94d5-71fdc7f4e17f` | `78d7f77fb05cfd84…` | 10 |
| Full corpus × off | `archive:7eb790c0-dea2-4aa2-8c20-dd3f52d53be1` | `fe86d338107950b0…` | 10 |

To use a retained input, resolve it through the archive (`just evidence-resolve <archive.json> <original path>`); see [Retaining execution evidence](../evidence-archives.md). That page also explains why run outputs a record cites must not stay in temporary storage. Since #1049, the rerun and archive recipes refuse such output roots; the sidecars are tracked in #1087.

## Notes for anyone citing these runs

- **Escherichia leave-out for arsenate reduction, Tier 1+2:** two outputs report it, with different transforms.
  - `native_composition_comparison.tsv` (the sign of the study-mean log1p z-contrast) gives 204:143 → 160:160 under off.
  - `effects.tsv.gz`, representation `shared_fraction` (raw fraction difference), gives 204:143 → 161:159 under off and 203:144 → 162:158 under phenotype.

  Name the file and representation whenever you cite it.
- **"Full corpus":** this names three populations.
  - The full-corpus cores above: 7,111 profiles from 1,212 studies.
  - The full experiments table: 4,435 rows from 1,245 studies.
  - The canonical human 16S/WMS profile table: 5,629 profiles from 1,047 studies.

  State which one a number uses.
- **The density adjustment's `relative_change_pct`** doesn't measure attenuation for traits with a negative raw mean (#1084). Coccus and sphere have positive raw means, so their values stand.
