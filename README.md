# hev-case-defintion

Code and minimal datasets for *Beyond Acute Jaundice: challenges refining case
definitions for suspected hepatitis E in South Sudan and Bangladesh*.

## Layout

```
R/
  data.R       load the datasets, symptom lists, 0/1 recoding, combined-data weights
  metrics.R    sensitivity/specificity/PPV/NPV with Wilson CIs; international case definitions
  combos.R     evaluate every symptom combination; selection rule; weighted bootstrap
  ensemble.R   SuperLearner leave-one-out AUCs and ROC plots
01_descriptives_and_international_case_defs.Rmd   Table 3 (partial), Table 4, Table S4
02_ensemble_models.Rmd                            Figure 1, Figure S2, Figure S4 and their AUCs
03_alternative_case_definitions.Rmd               Table 5, Table S5
04_combined_populations.Rmd                       Table S6, Figure S5, Discussion NPV
dat_*.rds                                         analytic datasets
```

Open the folder as an RStudio project (or knit from the repo root); all paths go
through `here()`. Notebooks 01, 03 and 04 each knit in under a minute. Notebook 02
refits the ensemble once per patient and takes several hours sequentially; see the
note at the top of that file for running it in parallel.

## Data

| File | Population | Symptom encoding |
|------|------------|------------------|
| `dat_ssd_adjusted.rds` | South Sudan (Bentiu), aged 14 and over, n = 421 | `"Yes"`/`"No"` factors |
| `dat_ssd_adjusted_children.rds` | South Sudan, under 14, n = 394 | `"Yes"`/`"No"` factors |
| `dat_bd_adjusted.rds` | Bangladesh, aged 14 and over, n = 1,625 | `"0"`/`"1"` factors |

All patients had jaundice onset within 30 days; `infect_bin` = anti-HEV IgM detected
by ELISA. Use `load_ssd(binary = TRUE)` / `load_bd(binary = TRUE)` to get 0/1
symptom columns. Never compare the raw factors with `== 1`: for the South Sudan
data that is always `FALSE`.

Not in these datasets, so not reproducible from this repository: IgG, the jaundice
signs (yellow skin/eyes, dark urine, pale stools), hospitalisation, death,
co-infection results (Table S1), and Bangladesh symptoms on the day of admission
(Figure S3).

## Methods details that the code fixes

- **Confidence intervals** for single-population sensitivity, specificity, PPV and NPV
  are Wilson score intervals (`prop.test(correct = FALSE)`).
- **Alternative case definitions** are selected among informative definitions
  (sensitivity + specificity > 1) whose sensitivity and specificity both exceed a
  threshold: 0.6 for South Sudan (both rules), 0.5 for Bangladesh, the combined data,
  and children "at least two", and 0.4 for children "at least one". The highest
  specificity is chosen, then the highest sensitivity, then the fewest symptoms.
- **Combined-population CIs** use a bootstrap stratified by population (1,000
  resamples), so the 50/50 weights stay fixed in every resample.

## Random seeds

| Analysis | Seed |
|----------|------|
| Ensemble models, South Sudan (aged 14 and over and under 14) | `1234` |
| Ensemble models, Bangladesh | `12345` |
| Combined-population bootstrap | `123` |

Each seed is set at the top of its chunk, so chunks must be run in order to reproduce
results exactly. Notebooks 01 and 03 are deterministic.
