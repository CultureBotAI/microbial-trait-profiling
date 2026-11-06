# Literature-review prompts

Prompts shared by Jed Dongjin Kim-Ozaeta on 2026-10-14 for the PD-cellulose figures. Each prompt keeps its original wording, with a header naming the figure it is about. Curly quotation marks were normalized to straight ones, and Markdown headings were added to prompt 4.

| File | Figure | Signal(s) |
|---|---|---|
| [01_pd_halophilic_arsenate_reduction.md](01_pd_halophilic_arsenate_reduction.md) | `pd_condition_taxon_trait_network.png`, Panel B | Parkinson's disease: moderately halophilic; arsenate reduction |
| [02_covid19_cellulose_catabolism.md](02_covid19_cellulose_catabolism.md) | same, Panel B | COVID-19: cellulose catabolic process |
| [03_crc_cellulose_catabolism.md](03_crc_cellulose_catabolism.md) | same, Panel B | Colorectal cancer: cellulose catabolic process |
| [04_pd_six_traits_influence_heatmap_audit.md](04_pd_six_traits_influence_heatmap_audit.md) | `pd_exploratory_taxa_six_traits_influence_heatmap.png` | Citation audit of every displayed taxon–trait cell |

## Status

Not run yet. The prompts quote the figures and tables as they stood on branch `clustered-heatmap-generation` at commit `205257c`. The later snapshot summarized below incorporated the canonical profile build (#482) and each profile's own experiment condition and body site (#555), so the "Local result(s)" in prompts 1-3 were already superseded. So was the `q = 0.164` that prompts 2 and 3 repeat in the paragraph they ask for, in the last deliverable bullet ([`02_covid19_cellulose_catabolism.md:53`](02_covid19_cellulose_catabolism.md) and [`03_crc_cellulose_catabolism.md:54`](03_crc_cellulose_catabolism.md)); both values in that later snapshot were `q = 0.0893`. In that snapshot, no condition-trait result passed `q <= 0.05`; prompt 1's statement that moderately halophilic was the only FDR-significant result no longer held. The table is historical, not a statement about current `main` or a later verified rerun.

The original prompts remain quoted records. For reuse, their raw p-values are
nominal. Their "global q" means BH across finite Wilcoxon p-values for supported
fecal condition–trait rows, not just the displayed rows and not the separate
taxon screen (`mt-script-pd-trait-wilcoxon` in the
[family registry](../../methods/multiple_testing_families.yaml)). The existing
decision tiers are q ≤ 0.05 and exploratory 0.05 < q ≤ 0.20 within that family;
BH within this selected family does not correct its selection. Neither the
original prompts nor the table below records the executed family size, so they
are incomplete sources for a fresh significance claim under the
[multiple-testing policy](../../methods/multiple_testing_policy.md).
Before running prompts 1–3, obtain a verified, population-matched
`pd_condition_taxon_trait_network_condition_trait_tests.tsv`, report its executed
family size and selection rule, and replace the local-result sections and the
deliverable q quotations in a new working prompt. Do not silently update these
preserved quotes or reuse the historical significance claims as current results.

| Prompt | Condition × trait | Paired studies | Mean contrast | Median contrast | Positive/negative/zero studies | Nominal p | BH q (condition–trait family) |
|---|---|---:|---:|---:|---|---:|---:|
| 1 | Parkinson's disease × moderately halophilic | 47 | +0.4971 | +0.5205 | 33/8/6 | 0.000195 | 0.0795 |
| 1 | Parkinson's disease × arsenate reduction | 47 | +0.5199 | +0.1788 | 25/6/16 | 0.000698 | 0.0893 |
| 2 | COVID-19 × cellulose catabolic process | 13 | −1.2975 | −1.3061 | 0/12/1 | 0.000488 | 0.0893 |
| 3 | Colorectal cancer × cellulose catabolic process | 28 | −0.8509 | −0.8242 | 3/19/6 | 0.000876 | 0.0893 |

All four rows in this historical snapshot are in that family's exploratory tier (`0.05 < q <= 0.20`), not FDR-significant. Issue #494 tracks readiness and the run plan.

Prompt 4 links the `pd_exploratory_taxa_six_traits_influence_heatmap.png` of commit `205257c`, which gave every experiment its study's majority label. The revision documented under #584 used each profile's own labels and all six traits, each marked with its observed direction in the six-trait lollipop. That revision showed 30 taxa (10 species, 12 genera and 8 families), 25 of them matched to the PD literature reference, so its displayed taxa differed from the linked version. Audit the intended run's figure and its matching `pd_exploratory_taxa_six_traits_reference_table.tsv`; neither historical description supplies the membership of a later run.

`METPO:2000103` in the prompts' identifiers is the KG-Microbe predicate; the traits are `GO:0030245` (cellulose catabolic process) and `pathways:arsenate_reduction`.
