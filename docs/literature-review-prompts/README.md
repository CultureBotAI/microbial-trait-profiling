# Literature-review prompts

Prompts shared by Jed Dongjin Kim-Ozaeta on 2026-10-14 for the PD-cellulose figures. Each prompt keeps its original wording, with a header naming the figure it is about. Curly quotation marks were normalized to straight ones, and Markdown headings were added to prompt 4.

| File | Figure | Signal(s) |
|---|---|---|
| [01_pd_halophilic_arsenate_reduction.md](01_pd_halophilic_arsenate_reduction.md) | `pd_condition_taxon_trait_network.png`, Panel B | Parkinson's disease: moderately halophilic; arsenate reduction |
| [02_covid19_cellulose_catabolism.md](02_covid19_cellulose_catabolism.md) | same, Panel B | COVID-19: cellulose catabolic process |
| [03_crc_cellulose_catabolism.md](03_crc_cellulose_catabolism.md) | same, Panel B | Colorectal cancer: cellulose catabolic process |
| [04_pd_six_traits_influence_heatmap_audit.md](04_pd_six_traits_influence_heatmap_audit.md) | `pd_exploratory_taxa_six_traits_influence_heatmap.png` | Citation audit of every displayed taxon–trait cell |

## Status

Not run yet. The prompts quote the figures and tables as they stood on branch `clustered-heatmap-generation` at commit `205257c`. Those files are now on `main` and have since been regenerated from the canonical profile build (#482) and with each profile's own experiment condition and body site (#555), so the "Local result(s)" in prompts 1-3 are out of date. So is the `q = 0.164` that prompts 2 and 3 repeat in the paragraph they ask for, in the last deliverable bullet ([`02_covid19_cellulose_catabolism.md:53`](02_covid19_cellulose_catabolism.md) and [`03_crc_cellulose_catabolism.md:54`](03_crc_cellulose_catabolism.md)); both values are now `q = 0.0893`. On `main`, no condition-trait result passes `q <= 0.05`; prompt 1's statement that moderately halophilic is the only FDR-significant result no longer holds. Before running prompts 1-3, replace their local results, and the q in those two deliverable bullets, with the current rows of `pd_condition_taxon_trait_network_condition_trait_tests.tsv`:

| Prompt | Condition × trait | Paired studies | Mean contrast | Median contrast | Positive/negative/zero studies | Raw p | Global q |
|---|---|---:|---:|---:|---|---:|---:|
| 1 | Parkinson's disease × moderately halophilic | 47 | +0.4971 | +0.5205 | 33/8/6 | 0.000195 | 0.0795 |
| 1 | Parkinson's disease × arsenate reduction | 47 | +0.5199 | +0.1788 | 25/6/16 | 0.000698 | 0.0893 |
| 2 | COVID-19 × cellulose catabolic process | 13 | −1.2975 | −1.3061 | 0/12/1 | 0.000488 | 0.0893 |
| 3 | Colorectal cancer × cellulose catabolic process | 28 | −0.8509 | −0.8242 | 3/19/6 | 0.000876 | 0.0893 |

All four are exploratory (`0.05 < q <= 0.20`), not FDR-significant. Issue #494 tracks readiness and the run plan.

Prompt 4 audits `pd_exploratory_taxa_six_traits_influence_heatmap.png`, which is stale: it still gives every experiment its study's majority label, and its column headers mark coccus shaped and sphere shaped as increased, whereas the current six-trait lollipop puts both on the negative side. `just pd-trait-taxa-heatmap` stops rather than rebuild it until #584 settles how the analysis picks its traits. Do not run prompt 4 until the heatmap has been regenerated after #584.

`METPO:2000103` in the prompts' identifiers is the KG-Microbe predicate; the traits are `GO:0030245` (cellulose catabolic process) and `pathways:arsenate_reduction`.
