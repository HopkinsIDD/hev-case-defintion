# hev-case-defintion

Code and datasets for each analyses are as follows:
<img width="4065" height="1365" alt="GitRepo summary table_aim2" src="https://github.com/user-attachments/assets/7e2bc459-e163-4945-8fb9-46a938a2a641" />

## Running the code

All scripts load data with `here()`, so run them from a checkout of this repository (open the folder as an RStudio project, or knit each `.Rmd` from the repo root). The analytic datasets are in the repo root:

| File | Site | Symptom encoding |
|------|------|------------------|
| `dat_ssd_adjusted.rds` | South Sudan (Bentiu) | `"Yes"`/`"No"` factors |
| `dat_ssd_adjusted_children.rds` | South Sudan, children | `"Yes"`/`"No"` factors |
| `dat_bd_adjusted.rds` | Bangladesh | `"0"`/`"1"` factors |

The exhaustive-combination scripts (`*_combos_*.Rmd`, `BD_SSD_joint_*.Rmd`) write intermediate case-definition batches to `batches/` (git-ignored) and read them back from there.

## Random seeds

| Script | Seed | Used for |
|--------|------|----------|
| `SSD_ensembles.Rmd` | `1234` | SuperLearner LOOCV and bootstrap AUC CIs |
| `SSD_ensembles_children.Rmd` | `1234` | SuperLearner LOOCV and bootstrap AUC CIs (seed added Nov 2026; earlier runs were unseeded) |
| `BD_ensembles.Rmd` | `12345` | SuperLearner LOOCV and bootstrap AUC CIs |
| `BD_SSD_joint_atleast1.Rmd` | `123` | Site-stratified bootstrap CIs (1,000 resamples) |
| `BD_SSD_joint_atleast2.Rmd` | `123` | Site-stratified bootstrap CIs (1,000 resamples) |

The seed is set once, before the first stochastic step, so the chunks must be run in order to reproduce results exactly. The other combination scripts are deterministic.
