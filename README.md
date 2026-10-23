# LAL-D individual participant data review: data and code

Data and R code for *Lysosomal acid lipase deficiency (LAL-D). It is necessary to change to multimodal therapy: A systematic review and long-term follow-up of Mexican probands* (PROSPERO CRD42024505308).

| File | Content |
|---|---|
| `data/lald_ipd.csv` | 463 participants from 234 published reports and one institutional case series |
| `data/lald_aggregate.csv` | 426 group-level results from 30 reports |
| `data/lald_reports.csv` | the 318 reports retrieved and the institutional series, with PubMed identifier and DOI |
| `data/data_dictionary.csv` | every variable of the three files |
| `01_prepare.R` | shared functions and derived variables |
| `02_estimates.R` | the estimates reported in the text |
| `03_tables.R` | Tables 1, 2 and S4–S8 |
| `04_figures.R` | Figures 1–3 |

R 4.3 or later with `logistf`, `metafor`, `mblm`, `sandwich`, `survival` and `boot`:

```
Rscript run_all.R
```

Tables are written to `results/` and figures to `figures/`. Laboratory values are in mg/dL and U/L.

Code: MIT (`LICENSE`). Data: CC BY 4.0 (`LICENSE-DATA`).
