# Manuscript audit

**Audited:** 2026-09-04 · **Last updated:** 2026-09-08 (204→261 correction family applied to the manuscript)
**Manuscript:** *Microbiome trait profiles across health and disease* (working title), targeting mSystems.

This file maps the manuscript onto this repository: where the manuscript lives, which commits and
data releases its numbers are pinned to, which of its published numbers are known to need
correction, and which open bsdbng issues bear on which manuscript claims. It exists because the
manuscript is maintained *outside* this repo, so contributors here otherwise cannot see which
claims their branches and issues affect. Update it when the manuscript's pins move or when an
issue listed here is closed.

## Where the manuscript lives

- The manuscript is maintained in a local rrwrite manuscript directory
  (`repo-research-writer/manuscript/bsdbng_v1`, a git repo with **no remote**), whose **source of
  truth is a shared Google Doc**; local `sections/*.md` and `full_manuscript.md` are derived
  mirrors. Its living backlog is that repo's `NEXT_TASKS.md`.
- Main text is essentially at the 5,000-word limit (last measured 4,995/5,000 on 2026-08-19; the
  2026-09-08 edit pass net-added ~4 words, so ≈4,999 — re-measure before any further additions).
  Ten figures,
  shape/alt-text verified 2026-08-19. Title page, funding, ethics and conflicts sections were
  still absent at the 2026-08-19 state review.

## Data provenance pins

Both analysis arms are now computed on the **same corpus**: BugSigDBExports **v1.3.1**
(Zenodo record 19740674, DOI 10.5281/zenodo.19740674, a reviewed-only export — see #130),
**human samples only**, paired experiments only.

| item | pin |
|---|---|
| Analysis A (composition / paired contrast) | `clustered-heatmap-generation` (CultureBotAI/bsdbng); manuscript pins the v1.3.1 regeneration `3afeecd`. Safe despite later commits: no analysis artifact has changed since, and the tip `c519f52` (2026-09-02) is documentation-only |
| Analysis B (ORA / enrichment) | `enrichment-analysis-heatmap-cluster` (CultureBotAI/bsdbng) @ `c75744a` (2026-08-31 v1.3.1 rerun; also the tip). Note: authored `Zoo Bot <bot@example.com>`, landed without an issue or PR — see #155 |
| KG-Microbe inputs | unchanged since April: `merged-kg_edges.tsv` md5 `555ecf5ba5d89288cdc44e368375374f`, `merged-kg_nodes.tsv` md5 `dd76b4656f9341db25a7c4afa96cd58e` |
| Corpus | 5,829 profiles in the wide matrix, of which 4,670 (= 2 × 2,335 pairs) enter the paired analyses; 4,517 testable; 315 traits tested by A, 267 significant traits exported by B |

The two committed artifacts every number below reproduces from:
`trait-visualization/trait-profile-investigation/trait_profile_paired_contrast_trait_stats.tsv`
(arm A @ `3afeecd`: columns `mean_increased_minus_decreased`, `q_value`, `t_statistic`) and
`data/derived/trait_enrichment_traits_overall_ranked.tsv` (arm B @ `c75744a`).

Consistency was verified independently on 2026-09-01: B's committed
`study_trait_profiles_wide.tsv` row count matches A's regeneration, B's 4,670 paired profiles
equal exactly 2 × A's 2,335 pairs, and the ranked table staged in the manuscript repo is
byte-identical (SHA-256) to the artifact at `c75744a`.

## Headline result and its verification status

**Table 1: six convergent traits** — significant at q ≤ 0.001 in both framings, direction-concordant
in all six (independently recomputed 2026-09-01, all counts match the manuscript table):

| trait | comp q | comp Δ | ORA inc/dec hits | ORA mean fold |
|---|---|---|---|---|
| facultatively aerobic (METPO:1000608) | 2.8e-10 | +0.194 | 531 / 329 | 2.68 |
| cellulose catabolic (METPO:2000103\|GO:0030245) | 6.9e-07 | −0.162 | 275 / 385 | 10.98 |
| coccus shaped (METPO:1000668) | 1.1e-06 | +0.161 | 498 / 366 | 3.20 |
| sphere shaped (METPO:1000683) | 1.9e-05 | +0.144 | 306 / 209 | 3.51 |
| arsenate reduction (METPO:2000103\|pathways:arsenate_reduction) | 1.9e-05 | +0.151 | 178 / 89 | 5.97 |
| moderately halophilic (METPO:1000623) | 3.1e-05 | +0.138 | 157 / 105 | 2.66 |

Composition-only traits: **zero** — every trait the composition arm calls significant is also
called by the ORA arm, so the convergent core's size is set by the composition arm's sensitivity,
not by disagreement between framings. The smallest significant standardized effect it observed is
dz = t/√n = 5.008/√2335 = 0.1036 (METPO:1000623); every enrichment-only trait it measured falls
below that value.

## Corrections applied to the manuscript text on 2026-09-08

A composite-parent keying bug in the manuscript repo's category scripts (fixed there 2026-09-03,
merged 2026-09-04 as `20ce6db`) silently dropped 57 `METPO:2000103|pathways:*`/`GO:*` traits from
the categorisation. The corrected values were **applied to the Google Doc and derived md on
2026-09-08** — 14 verified replace-all edits (11 + 3 adversarial-review follow-ups, every one
with occurrences-changed = expected = 1); edit plans, apply reports and pre-edit backups are
committed in the manuscript repo, merge `9f5f930`. The committed backups are pre-edit by
construction, so the second backup independently evidences the 11 first-pass edits in the doc,
while the 3 follow-ups are evidenced by their apply report and the md mirrors only:

| quantity | manuscript said (pre-2026-09-08) | manuscript now says |
|---|---|---|
| enrichment-only traits | 204 (203 measured by A) | **261 (260 measured by A)** |
| below A's detection floor | 100% of 203 measured (the Abstract over-claimed it for all 204) | 100% — scoped to the **260 A measured** (one trait, non-spore forming METPO:1000872, is untested by A) |
| enrichment-only median / max standardized effect (dz) | 0.022 / 0.080 | 0.022 / **0.088** |
| pairs needed for the median enrichment-only effect | ~53,000 (23×) | **~52,000 (about 22×)** |
| median enriched profiles, ORA-only vs convergent | 203 vs 588 | **176** vs 588 (convergent median is 587.5; 588 is its rounding, unchanged by the pass) |
| reliability of A's sub-threshold dz | 0.57 | **0.606** |
| disattenuated Spearman vs unique_profiles / mean_fold | +0.24 / −0.20 | **+0.250 / −0.155** |

Convergent membership and every Table 1 **data** row were unchanged by the correction (the two
Table 1 edits touched only the caption and the ORA-only section-label count; the
top-10-by-profile-support block was recomputed under the corrected keying before editing —
identical membership). The q ≤ 0.05 sweep values (6→37 convergent, 261→230 ORA-only, 0→11
composition-only) are not in the manuscript text; the old 6→47 was buggy script output that
appeared in the text for ~40 minutes on 2026-09-02 (added `72b6c7b`, removed `6d922b4`) — the
2026-09-04 version of this audit carried it as a stale row.

Spearman method (so the row reproduces from the two artifacts above): correlate arm A's |t|
(equivalently dz = |t|/√2335) against arm B's `unique_profiles` and `mean_fold` over the 260
measured ORA-only traits — raw +0.194 / −0.121 as the script prints them (unrounded +0.1946) —
then divide by the attenuation factor
√reliability = √0.606 ≈ 0.778, where reliability is the estimated true-to-observed variance ratio
of A's per-trait dz (sampling variance taken as 1/n).

## Open manuscript-blocking tasks (tracked in the manuscript repo's NEXT_TASKS.md)

1. ~~Apply the 204 → 261 correction family (above) to the Google Doc and derived md.~~ ✅ Done
   2026-09-08.
2. Corpus numbers in the text still say 12,108 profiles / 4,784 pairs. The 12,108 was a
   wide-matrix profile count, so the correction chain is: **12,108 → 5,829** profiles in the
   matrix, of which **4,670** (= 2 × **2,335** pairs, replacing 4,784) enter the paired analyses,
   **4,517** testable. Substituting 4,670 for 12,108 directly would compare a paired subset to a
   matrix total.
3. Figures 2B/4A/4B are one ORA run behind — regenerated from `45c0593`, not `c75744a`.
4. Table 1 provenance strings reference the pre-rename output paths (see #155: the `_p0.001`
   suffixes were dropped at `c75744a`).
5. Title page, funding, ethics, conflicts; cover letter; figure paste from the staging doc
   (needs a human with Doc access).

## Which open bsdbng issues bear on which manuscript claims

**Directly threaten Table 1 / headline numbers**

- **#133 — composition pairs treated as independent; no effect-size floor.** The q-values behind
  Table 1 membership. Note the issue's numbers predate v1.3.1 (it says 4,784 pairs from 1,704
  studies; the current corpus is 2,335 pairs from 959 studies, with 48 traits at q ≤ 0.05, not
  69). The core defect is corpus-independent, and its measured consequence — study-averaging
  inflated p-values by 2–5 orders of magnitude on the old corpus; not yet re-measured at
  n = 2,335 — matters because three of the six convergent traits sit at q between 1.9e-05 and
  3.1e-05, close enough to the 0.001 cutoff that membership could change. Highest-value open
  statistical issue for the paper.
- **#131 — ORA urn draws trait assertions, not taxa (~100× effective-n inflation).** All ORA
  p/q-values are anticonservative; the 267-trait significant set and the 261 enrichment-only
  count would shrink under a taxon-level urn. Direction concordance (the load-bearing evidence
  for the convergent six) is more robust than the q-values, but the ORA-only block is not.
- **#135 — no golden-value tests in either arm.** Every number above is produced by hand-rolled,
  duplicated, numerically untested statistics code.

**Affect specific sections/figures**

- **#132** — enrichment/depletion mirror artifact and structural-zero depletion: bears on the
  inc/dec hit columns and on any use of the depletion arm (100,999 "depletions" at `c75744a` are
  still dominated by structural zeros).
- **#134** — the condition-trait-module heatmap sums increased and decreased profiles together;
  the "condition × body-site trait modules" interpretation rests on it.
- **#128** — the autism case study pools two labels with opposite directional skew.
- **#155** — `c75744a` regressions: a duplicate Zenodo downloader with no checksum verification
  (main already has the hardened, SHA-pinned one from PR #150), deleted ORA methodology summary
  (Methods' citation target), weakened tests, renamed outputs (breaks the manuscript's provenance
  strings).
- **#142 / #144 / #147** — metadata joins and new directional figures that would extend the
  figure set.
- **#118 / #119 / #123** — per-body-site stratification and rollup-statistic caveats for the ORA
  arm's per-trait claims.
- **#127** — `plot_enrichment_histograms.py` counts experiments where it labels profiles,
  understating trait frequency in Figure 4's histograms.

**Resolved in substance by the v1.3.1 reruns / PR #150 but still open on the tracker**

- **#121** (MediaDive in arm B's trait space): resolved by a code-level exclusion
  (`scripts/trait_enrichment.py`, `trait_depletion.py`, the exporter, histogram and summary
  scripts all filter `mediadive.medium:`), and `trait_enrichment_results.tsv` carries 0 MediaDive
  rows at `c75744a`. The wide matrix itself still contains 2,290 `mediadive.medium:*` columns —
  the exclusion is applied downstream, not upstream.
- **#117** (top-10-per-profile truncation): reruns use `Max-Traits-Per-Row-Report: 0`.
- **#122** (arms disagree on unpaired experiments): both arms' *paired analyses* now drop
  unpaired experiments; arm A's other committed figures still run over all 5,829 profiles,
  unpaired included.
- **#129** (non-human samples): both arms now restrict to human host — B at ingest
  (`ingest(host_species="Homo sapiens")` since `c75744a`), A per-script from `full_dump.csv`
  (`HUMAN_HOST_SPECIES` in four figure scripts). Same net corpus, different mechanism.
- **#130** (vintage decision): both arms pin BugSigDBExports v1.3.1.
- **#151 / #152** (pin by SHA; verify downloaded bytes): implemented on main by PR #150.

These become closeable once the `enrichment-analysis-heatmap-cluster` branch lands through the
normal PR flow — which is currently blocked by the #155 regressions — plus, for #151/#152, a
check that nothing still routes around main's hardened downloader. Closing them now would let
the tracker get ahead of `main`.

**Statistical caveats recorded in the manuscript repo, not yet issues here**

- ~45% of trait columns are exact duplicates of another (both arms retain duplicates after the
  MediaDive exclusion: B's 267 traits cover 261 distinct vectors, A's 315 cover 305), so BH
  corrects over fewer independent hypotheses than it assumes and ranked lists can show one
  measurement as several rows. Table 1's six traits are all duplicate-group size 1.
- The ranked export's `rank` column is prevalence over *testable* profiles; its top rows are
  traits testable in as few as 20 profiles. Any manuscript use of that table should rank by
  absolute profile support instead (Table 1's ORA-only block already does).
