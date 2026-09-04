# Manuscript audit

**Audited:** 2026-09-04
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
- Main text is essentially at the 5,000-word limit (last measured 4,995/5,000). Ten figures,
  shape/alt-text verified 2026-08-19. Title page, funding, ethics and conflicts sections were
  still absent at the 2026-08-19 state review.

## Data provenance pins

Both analysis arms are now computed on the **same corpus**: BugSigDBExports **v1.3.1**
(Zenodo record 19740674, DOI 10.5281/zenodo.19740674, a reviewed-only export — see #130),
**human samples only**, paired experiments only.

| item | pin |
|---|---|
| Analysis A (composition / paired contrast) | branch `clustered-heatmap-generation`, manuscript pins the v1.3.1 regeneration `3afeecd` (tip `c519f52` as of 2026-09-04) |
| Analysis B (ORA / enrichment) | branch `enrichment-analysis-heatmap-cluster` @ `c75744a` (2026-08-31 v1.3.1 rerun; also the tip) |
| KG-Microbe inputs | unchanged since April: `merged-kg_edges.tsv` md5 `555ecf5ba5d89288cdc44e368375374f`, `merged-kg_nodes.tsv` md5 `dd76b4656f9341db25a7c4afa96cd58e` |
| Corpus | 5,829 profiles in the wide matrix; 4,670 paired profiles (= 2 × 2,335 pairs); 4,517 testable; 315 traits tested by A, 267 significant traits exported by B |

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
called by the ORA arm, so the convergent core's size is set by the composition arm's power
(minimum detectable effect |d| ≥ 0.1036 at n = 2,335 pairs), not by disagreement between framings.

## Known corrections not yet in the manuscript text

A composite-parent keying bug in the manuscript repo's category scripts (fixed there 2026-09-03 on
branch `fix/convergence-join-and-ora-only-count`, **not yet merged**) silently dropped 57
`METPO:2000103|pathways:*`/`GO:*` traits from the categorisation. The Google Doc and derived md
still carry the old values:

| quantity | manuscript says | corrected |
|---|---|---|
| enrichment-only traits | 204 (203 measured by A) | **261 (260 measured by A)** |
| below A's detection floor | 100% | 100% (claim survives) |
| enrichment-only median / max \|d\| | 0.022 / 0.080 | 0.022 / **0.088** |
| median enriched profiles, ORA-only vs convergent | 203 vs 588 | **176** vs 588 |
| disattenuated Spearman vs unique_profiles / mean_fold | +0.24 / −0.20 | **+0.25 / −0.155** |
| q ≤ 0.05 sweep | convergent 6→47 | **6→37**, ORA-only 261→230, comp-only 0→11 |

Convergent membership and every Table 1 row are unchanged by the correction.

## Open manuscript-blocking tasks (tracked in the manuscript repo's NEXT_TASKS.md)

1. Apply the 204 → 261 correction family (above) to the Google Doc and derived md.
2. Corpus numbers in the text still say 12,108 profiles / 4,784 pairs; v1.3.1 values are
   4,670 / 2,335 (and 4,517 testable).
3. Figures 2B/4A/4B are one ORA run behind — regenerated from `45c0593`, not `c75744a`.
4. Table 1 provenance strings reference the pre-rename output paths (see #155: the `_p0.001`
   suffixes were dropped at `c75744a`).
5. Title page, funding, ethics, conflicts; cover letter; figure paste from the staging doc
   (needs a human with Doc access).

## Which open bsdbng issues bear on which manuscript claims

**Directly threaten Table 1 / headline numbers**

- **#133 — composition pairs treated as independent; no effect-size floor.** The q-values behind
  Table 1 membership. A study-clustered rerun previously inflated p by 2–5 orders of magnitude;
  three of the six convergent traits sit at q between 1.9e-05 and 3.1e-05, close enough to the
  0.001 cutoff that membership could change. Highest-value open statistical issue for the paper.
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
- **#155** — `c75744a` regressions: duplicate unverified Zenodo downloader (vs #151/#152 on
  main), deleted ORA methodology summary (Methods' citation target), weakened tests, renamed
  outputs (breaks the manuscript's provenance strings).
- **#142 / #144 / #147** — metadata joins and new directional figures that would extend the
  figure set.
- **#118 / #119 / #123 / #127** — per-body-site stratification and rollup-statistic caveats for
  the ORA arm's per-trait claims.

**Resolved in substance by the v1.3.1 reruns but still open on the tracker**

- **#121** (MediaDive in arm B's trait space): 0 of 267 ranked traits are MediaDive at `c75744a`.
- **#117** (top-10-per-profile truncation): reruns use `Max-Traits-Per-Row-Report: 0`.
- **#122** (arms disagree on unpaired experiments): both arms now drop unpaired experiments.
- **#129** (non-human samples): both arms now ingest human-only.
- **#130** (vintage decision): both arms pin BugSigDBExports v1.3.1.

Closing these (with a pointer to `c75744a`/`3afeecd`) would make the tracker reflect reality.

**Statistical caveats recorded in the manuscript repo, not yet issues here**

- ~45% of trait columns are exact duplicates of another (both arms retain duplicates after the
  MediaDive exclusion: B's 267 traits cover 261 distinct vectors, A's 315 cover 305), so BH
  corrects over fewer independent hypotheses than it assumes and ranked lists can show one
  measurement as several rows. Table 1's six traits are all duplicate-group size 1.
- The ranked export's `rank` column is prevalence over *testable* profiles; its top rows are
  traits testable in as few as 20 profiles. Any manuscript use of that table should rank by
  absolute profile support instead (Table 1's ORA-only block already does).
